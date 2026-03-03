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


def _parse_ordermarkers2(path: Path) -> pd.DataFrame:
    """Parse OrderMarkers2 stdout and return a marker table.

    The OrderMarkers2 output is a plain text file with blocks like:
      #*** LG = 1 likelihood = ...
      #marker_number\tmale_position\tfemale_position ...
      498\t0.0\t0.0 ...

    We parse only the first 3 columns and keep linkage group id.
    """
    rows = []
    cur_lg: str | None = None
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line:
                continue
            line = line.rstrip("\n")
            if not line:
                continue

            if line.startswith("#*** LG"):
                # e.g. "#*** LG = 1 likelihood = -2472.9"
                try:
                    after = line.split("=", 1)[1].strip()
                    cur_lg = after.split()[0].strip()
                except Exception:
                    cur_lg = None
                continue

            if line.startswith("#"):
                continue

            parts = line.split("\t")
            if len(parts) < 3:
                # some Lep-MAP3 outputs may be space separated
                parts = line.split()
            if len(parts) < 3:
                continue

            try:
                mnum = int(parts[0])
                mpos = float(parts[1])
                fpos = float(parts[2])
            except Exception:
                continue

            rows.append({"lg": str(cur_lg or "NA"), "marker_number": mnum, "male_pos": mpos, "female_pos": fpos})

    df = pd.DataFrame(rows)
    if df.empty:
        return df

    # sexAveraged=1 is commonly used; still be robust.
    df["pos_cM"] = (df["male_pos"].astype(float) + df["female_pos"].astype(float)) / 2.0
    return df


def _iter_vcf_records(vcf_path: Path):
    """Yield (idx_1based, chrom, pos, id) for each variant line in a VCF/VCF.GZ."""
    import gzip

    opener = gzip.open if vcf_path.suffix.lower() == ".gz" else open
    with opener(vcf_path, "rt", encoding="utf-8", errors="ignore") as f:
        idx = 0
        for line in f:
            if not line or line.startswith("#"):
                continue
            idx += 1
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            chrom = parts[0]
            pos = parts[1]
            vid = parts[2]
            yield idx, chrom, pos, vid


def _map_marker_numbers_to_vcf_ids(vcf_path: Path, marker_numbers: list[int]) -> dict[int, str]:
    """Best-effort mapping from Lep-MAP3 marker_number to VCF ID.

    IMPORTANT: This assumes marker_number corresponds to the 1-based order of variants
    in the provided VCF (excluding header lines). This is true for common Lep-MAP3
    workflows that ingest a VCF via ParentCall2.

    If the assumption is violated (e.g., marker numbers are renumbered after filtering),
    the mapping will be incomplete or incorrect. We therefore:
      * fall back to using marker_number as marker id when mapping is missing
      * write a note file to explain the assumption
    """
    want = sorted({int(x) for x in marker_numbers if int(x) > 0})
    if not want:
        return {}

    want_set = set(want)
    max_idx = want[-1]
    out: dict[int, str] = {}
    seen: dict[str, int] = {}

    for idx, chrom, pos, vid in _iter_vcf_records(vcf_path):
        if idx > max_idx and len(out) == len(want_set):
            break
        if idx not in want_set:
            continue
        name = (vid or "").strip()
        if name in {"", "."}:
            name = f"{chrom}:{pos}"

        # make unique
        n = seen.get(name, 0)
        if n > 0:
            name2 = f"{name}__dup{n+1}"
        else:
            name2 = name
        seen[name] = n + 1

        out[int(idx)] = name2
        if len(out) == len(want_set):
            break
    return out


def _write_map_outputs(
    out_dir: Path,
    *,
    order_out: Path,
    vcf_path: Optional[Path] = None,
) -> Tuple[bool, List[str]]:
    """Parse OrderMarkers2 output and write standardized map outputs.

    Returns (ok, notes)
    """
    notes: List[str] = []

    if not order_out.exists() or order_out.stat().st_size <= 0:
        return False, ["OrderMarkers2 output not found or empty."]

    df = _parse_ordermarkers2(order_out)
    if df.empty:
        return False, ["Failed to parse OrderMarkers2 output."]

    # map marker_number -> marker id (best-effort)
    id_map: dict[int, str] = {}
    if vcf_path is not None and vcf_path.exists():
        try:
            id_map = _map_marker_numbers_to_vcf_ids(vcf_path, df["marker_number"].astype(int).tolist())
            notes.append(
                "marker_number→VCF ID mapping was attempted assuming marker_number is the 1-based variant index in the VCF (non-header lines)."
            )
            notes.append(f"VCF mapping hits: {len(id_map)}/{df.shape[0]} markers")
        except Exception as e:
            notes.append(f"VCF mapping failed: {e}")

    df["marker"] = df["marker_number"].map(lambda x: id_map.get(int(x), str(int(x))))

    # Sort and write
    df = df.sort_values(["lg", "pos_cM", "marker"], kind="mergesort")
    map_markers = df[["lg", "marker", "pos_cM", "marker_number", "male_pos", "female_pos"]].copy()
    map_markers.to_csv(out_dir / "map_markers.tsv", sep="\t", index=False)

    # mppR / common format: marker, chr, pos (cM)
    marker_map = df[["marker", "lg", "pos_cM"]].copy()
    marker_map = marker_map.rename(columns={"lg": "chr", "pos_cM": "pos"})
    marker_map.to_csv(out_dir / "marker_map.tsv", sep="\t", index=False)

    # LG length table
    grp = df.groupby("lg", as_index=False).agg(n_markers=("marker", "size"), length_cM=("pos_cM", "max"))
    grp = grp.sort_values(["lg"], kind="mergesort")
    grp.to_csv(out_dir / "map_lengths.tsv", sep="\t", index=False)

    # basic plot
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        plt.figure()
        plt.bar(range(len(grp)), grp["length_cM"].astype(float))
        plt.xticks(range(len(grp)), grp["lg"].astype(str), rotation=90)
        plt.xlabel("LG")
        plt.ylabel("length (cM)")
        plt.title("Linkage map length")
        plt.tight_layout()
        plt.savefig(out_dir / "map_lengths.png", dpi=150)
        plt.close()
    except Exception as e:
        notes.append(f"Failed to write map_lengths.png: {e}")

    # human-readable notes
    (out_dir / "map_note.txt").write_text("\n".join(notes) + "\n", encoding="utf-8")
    return True, notes


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
            #args_kv={**base_kv, "data": str(cur_data), "map": str(map_file), "sexAveraged": int(1)},
            args_kv={**base_kv, "data": str(cur_data), "map": str(map_file)},
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

    # Build standardized outputs for GUI + downstream tools
    order_copy = out_dir / order_out.name
    vcf_path = Path(str(params.get("vcf_path", "") or "")).expanduser()
    ok, notes = _write_map_outputs(out_dir, order_out=order_copy if order_copy.exists() else order_out, vcf_path=vcf_path if vcf_path.exists() else None)

    if ok:
        # keep a short status file (GUI doesn't rely on it, but users appreciate it)
        (out_dir / "status.txt").write_text("Lep-MAP3 map parsed successfully.\n", encoding="utf-8")
        produced_files.extend(["map_markers.tsv", "map_lengths.tsv", "marker_map.tsv", "map_lengths.png", "map_note.txt", "status.txt"])
    else:
        msg = "Lep-MAP3 pipeline executed, but map parsing failed.\n" + "\n".join(notes)
        _placeholder_outputs(out_dir, msg)
        produced_files.append("error_message.txt")

    _write_artifacts(
        out_dir,
        default_table="map_markers.tsv",
        default_plot="map_lengths.png",
        extra_files=sorted(set(produced_files + ["run.log"]))
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
        "tables": ["map_markers.tsv", "map_lengths.tsv", "marker_map.tsv"],
        "plots": ["map_lengths.png"],
        "files": extra_files or [],
        "note": "Lep-MAP3 scaffold: best-effort pipeline (ParentCall2/Filtering2/SeparateChromosomes2/OrderMarkers2) with map_markers.tsv + marker_map.tsv outputs.",
    }
    (out_dir / "artifacts.json").write_text(json.dumps(art, indent=2), encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
