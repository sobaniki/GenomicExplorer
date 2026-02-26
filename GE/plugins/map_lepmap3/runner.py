#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from __future__ import annotations

"""Lep-MAP3 linkage map scaffold runner.

Design goals (matching the project philosophy):
  * Keep GUI runnable even if Lep-MAP3 isn't installed yet.
  * When Lep-MAP3 is installed (jar + Java), run a standard pipeline:
      ParentCall2 -> Filtering2 -> SeparateChromosomes2 -> OrderMarkers2
    while capturing step logs and intermediate files.
  * Always write the usual artifacts so the ResultsPane can load something:
      map_markers.tsv, map_lengths.tsv, map_lengths.png, artifacts.json

This runner intentionally avoids overfitting to a single Lep-MAP3 invocation style.
The user can pass extra arguments per step, and can disable steps.
"""

import argparse
import json
import shutil
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd


def _log_write(log_path: Path, lines: List[str]) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _ensure_parent(p: Path) -> None:
    p.parent.mkdir(parents=True, exist_ok=True)


def _run_step(
    *,
    java_exe: str,
    jar_path: Path,
    mem_gb: int,
    module: str,
    args_kv: Dict[str, str],
    extra_args: str,
    out_stdout: Path,
    out_stderr: Path,
    cwd: Path,
) -> Tuple[int, List[str]]:
    """Run a Lep-MAP3 module and write stdout/stderr to files."""
    _ensure_parent(out_stdout)
    _ensure_parent(out_stderr)

    cmd: List[str] = [
        java_exe,
        f"-Xmx{int(mem_gb)}g",
        "-cp",
        str(jar_path),
        module,
    ]
    # Lep-MAP3 uses key=value style
    for k, v in args_kv.items():
        if v is None:
            continue
        v = str(v).strip()
        if v == "":
            continue
        cmd.append(f"{k}={v}")

    if extra_args:
        # allow users to pass raw tokens like "removeNonInformative=1" or "lodLimit=10"
        # split by whitespace (simple, consistent with other scaffolds)
        cmd.extend([t for t in extra_args.split() if t.strip()])

    log_lines = [f"[lepmap3] RUN: {' '.join(cmd)}"]
    with out_stdout.open("w", encoding="utf-8", errors="ignore") as fo, out_stderr.open(
        "w", encoding="utf-8", errors="ignore"
    ) as fe:
        try:
            proc = subprocess.run(cmd, cwd=str(cwd), stdout=fo, stderr=fe, check=False)
            return int(proc.returncode), log_lines
        except Exception as e:
            fe.write(str(e) + "\n")
            return 1, log_lines + [f"[lepmap3] EXCEPTION: {e}"]


def _placeholder_outputs(out_dir: Path, msg: str) -> None:
    """Write placeholder outputs so GUI remains functional."""
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "error_message.txt").write_text(msg + "\n", encoding="utf-8")

    # Minimal TSVs
    pd.DataFrame(
        [{"chr": "NA", "marker": "NA", "pos_cM": 0.0, "note": "placeholder"}]
    ).to_csv(out_dir / "map_markers.tsv", sep="\t", index=False)
    pd.DataFrame(
        [{"chr": "NA", "n_markers": 0, "length_cM": 0.0, "note": "placeholder"}]
    ).to_csv(out_dir / "map_lengths.tsv", sep="\t", index=False)

    # Simple plot placeholder (matplotlib only if available; otherwise skip)
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        df = pd.read_csv(out_dir / "map_lengths.tsv", sep="\t")
        plt.figure()
        plt.bar(range(len(df)), df["length_cM"])
        plt.xlabel("chr")
        plt.ylabel("length_cM")
        plt.title("Linkage map length (placeholder)")
        plt.tight_layout()
        plt.savefig(out_dir / "map_lengths.png", dpi=150)
        plt.close()
    except Exception:
        # ignore
        pass


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--params", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    params_path = Path(args.params)
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    log_lines: List[str] = []
    log_lines.append("[map_lepmap3] start")
    log_lines.append(f"[map_lepmap3] params_path={params_path}")
    log_lines.append(f"[map_lepmap3] out_dir={out_dir}")

    params = json.loads(params_path.read_text(encoding="utf-8"))

    jar = Path(str(params.get("lepmap3_jar", "") or "")).expanduser()
    java_exe = str(params.get("java", "java") or "java")
    mem_gb = int(float(params.get("java_mem_gb", 8) or 8))
    threads = int(float(params.get("threads", 4) or 4))
    input_path = Path(str(params.get("input_path", "") or "")).expanduser()
    vcf = str(params.get("vcf_path", ""))
    prefix = str(params.get("output_prefix", "lepmap3") or "lepmap3")

    # pipeline toggles
    do_parent = str(params.get("do_parentcall2", "TRUE")).upper().startswith("T")
    do_filter = str(params.get("do_filtering2", "TRUE")).upper().startswith("T")
    do_sep = str(params.get("do_separatechromosomes2", "TRUE")).upper().startswith("T")
    do_order = str(params.get("do_ordermarkers2", "TRUE")).upper().startswith("T")

    # per-step extra args (raw tokens)
    extra_parent = str(params.get("parentcall2_extra", "") or "")
    extra_filter = str(params.get("filtering2_extra", "") or "")
    extra_sep = str(params.get("separatechromosomes2_extra", "") or "")
    extra_order = str(params.get("ordermarkers2_extra", "") or "")

    # Standard args
    # NOTE: Lep-MAP3 uses data=<file>. Many modules also accept "threads".
    #base_kv = {"threads": str(threads)}
    base_kv = {"numThreads": str(threads)}

    work_dir = out_dir / "lepmap3_work"
    work_dir.mkdir(parents=True, exist_ok=True)

    if not jar.exists():
        msg = f"Lep-MAP3 jar not found: {jar}. Set 'lepmap3_jar' in GUI and ensure Java is available."
        log_lines.append("[map_lepmap3] " + msg)
        _log_write(out_dir / "run.log", log_lines)
        _placeholder_outputs(out_dir, msg)
        _write_artifacts(out_dir, default_table="map_markers.tsv", default_plot="map_lengths.png")
        return 0

    if not input_path.exists():
        msg = f"Input file not found: {input_path}. Provide Lep-MAP3 input data file."
        log_lines.append("[map_lepmap3] " + msg)
        _log_write(out_dir / "run.log", log_lines)
        _placeholder_outputs(out_dir, msg)
        _write_artifacts(out_dir, default_table="map_markers.tsv", default_plot="map_lengths.png")
        return 0

    # ensure Java exists in PATH if java_exe is not absolute
    if Path(java_exe).name == java_exe and shutil.which(java_exe) is None:
        msg = f"Java executable not found in PATH: {java_exe}. Install Java or set 'java' path."
        log_lines.append("[map_lepmap3] " + msg)
        _log_write(out_dir / "run.log", log_lines)
        _placeholder_outputs(out_dir, msg)
        _write_artifacts(out_dir, default_table="map_markers.tsv", default_plot="map_lengths.png")
        return 0

    # Pipeline files
    cur_data = input_path
    parent_out = work_dir / f"{prefix}.parentcall2.txt"
    filter_out = work_dir / f"{prefix}.filtering2.txt"
    sep_out = work_dir / f"{prefix}.separatechromosomes2.txt"
    order_out = work_dir / f"{prefix}.ordermarkers2.txt"
    
    print(vcf)

    # (1) ParentCall2
    if do_parent:
        rc, step_log = _run_step(
            java_exe=java_exe,
            jar_path=jar,
            mem_gb=mem_gb,
            module="ParentCall2",
            #args_kv={**base_kv, "data": str(cur_data)},
            #args_kv={"data": str(cur_data), "vcfFile": str(vcf), "removeNonInformative": int(1)},
            args_kv={"data": str(cur_data), "vcfFile": str(vcf)},
            extra_args=extra_parent,
            out_stdout=parent_out,
            out_stderr=work_dir / f"{prefix}.parentcall2.stderr.txt",
            cwd=work_dir,
        )
        log_lines.extend(step_log)
        log_lines.append(f"[lepmap3] ParentCall2 rc={rc}")
        if rc != 0:
            log_lines.append("[lepmap3] ParentCall2 failed; continuing with input as-is")
        else:
            cur_data = parent_out

    # (2) Filtering2
    if do_filter:
        rc, step_log = _run_step(
            java_exe=java_exe,
            jar_path=jar,
            mem_gb=mem_gb,
            module="Filtering2",
            #args_kv={**base_kv, "data": str(cur_data)},
            #args_kv={"data": str(cur_data), "dataTolerance": float(0.001)},
            args_kv={"data": str(cur_data)},
            extra_args=extra_filter,
            out_stdout=filter_out,
            out_stderr=work_dir / f"{prefix}.filtering2.stderr.txt",
            cwd=work_dir,
        )
        log_lines.extend(step_log)
        log_lines.append(f"[lepmap3] Filtering2 rc={rc}")
        if rc != 0:
            log_lines.append("[lepmap3] Filtering2 failed; continuing with previous data")
        else:
            cur_data = filter_out

    # (3) SeparateChromosomes2
    map_file: Optional[Path] = None
    if do_sep:
        rc, step_log = _run_step(
            java_exe=java_exe,
            jar_path=jar,
            mem_gb=mem_gb,
            module="SeparateChromosomes2",
            #args_kv={**base_kv, "data": str(cur_data), "lodLimit": int(5), "theta": float(0.05), "distortionLod": int(1)},
            args_kv={**base_kv, "data": str(cur_data)},
            extra_args=extra_sep,
            out_stdout=sep_out,
            out_stderr=work_dir / f"{prefix}.separatechromosomes2.stderr.txt",
            cwd=work_dir,
        )
        log_lines.extend(step_log)
        log_lines.append(f"[lepmap3] SeparateChromosomes2 rc={rc}")
        if rc != 0:
            log_lines.append("[lepmap3] SeparateChromosomes2 failed; OrderMarkers2 may be skipped")
        else:
            map_file = sep_out

    # (4) OrderMarkers2
    if do_order and map_file is not None and map_file.exists():
        rc, step_log = _run_step(
            java_exe=java_exe,
            jar_path=jar,
            mem_gb=mem_gb,
            module="OrderMarkers2",
            args_kv={**base_kv, "data": str(cur_data), "map": str(map_file), "sexAveraged": int(1)},
            extra_args=extra_order,
            out_stdout=order_out,
            out_stderr=work_dir / f"{prefix}.ordermarkers2.stderr.txt",
            cwd=work_dir,
        )
        log_lines.extend(step_log)
        log_lines.append(f"[lepmap3] OrderMarkers2 rc={rc}")
    elif do_order:
        log_lines.append("[lepmap3] OrderMarkers2 skipped (no map file from SeparateChromosomes2)")

    # Copy key files to out_dir (stable names)
    produced_files: List[str] = []
    for fp in [parent_out, filter_out, sep_out, order_out]:
        if fp.exists() and fp.stat().st_size > 0:
            dest = out_dir / fp.name
            shutil.copyfile(fp, dest)
            produced_files.append(dest.name)
    # Also include stderr per step if exists
    for fp in work_dir.glob(f"{prefix}.*.stderr.txt"):
        if fp.exists() and fp.stat().st_size > 0:
            dest = out_dir / fp.name
            shutil.copyfile(fp, dest)
            produced_files.append(dest.name)

    _log_write(out_dir / "run.log", log_lines)

    # Minimal outputs: placeholders for now (parsing OrderMarkers2 can be implemented later)
    msg = "Lep-MAP3 pipeline executed (best-effort). Parse/convert outputs later if needed."
    if not produced_files:
        msg = "Lep-MAP3 did not produce any outputs (best-effort scaffold). Check run.log and stderr files."
    _placeholder_outputs(out_dir, msg)

    _write_artifacts(
        out_dir,
        default_table="map_markers.tsv",
        default_plot="map_lengths.png",
        extra_files=sorted(set(produced_files + ["run.log", "error_message.txt"]))
    )
    return 0


def _write_artifacts(
    out_dir: Path,
    default_table: str,
    default_plot: str,
    extra_files: Optional[List[str]] = None,
) -> None:
    art = {
        "table": default_table,
        "plot": default_plot,
        "tables": ["map_markers.tsv", "map_lengths.tsv"],
        "plots": ["map_lengths.png"],
        "files": extra_files or [],
        "note": "Lep-MAP3 scaffold: provide jar + java, then tune step args; parsing outputs can be added later.",
    }
    (out_dir / "artifacts.json").write_text(json.dumps(art, indent=2), encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
