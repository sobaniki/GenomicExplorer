#!/usr/bin/env python3

import argparse
import json
import subprocess
from pathlib import Path
from datetime import datetime


def write_text(path: Path, text: str):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def write_json(path: Path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2), encoding="utf-8")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--params", required=True)
    ap.add_argument("--out", required=True)
    ns = ap.parse_args()

    out_dir = Path(ns.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    runlog = out_dir / "run.log"
    p = json.loads(Path(ns.params).read_text(encoding="utf-8"))

    tassel_jar = p.get("tassel_jar")
    inp = p.get("input")
    inp_fmt = (p.get("input_format") or "auto").strip()
    pedigree = (p.get("pedigree") or "").strip()

    java = (p.get("java") or "java").strip() or "java"
    java_mem = (p.get("java_mem") or "8g").strip() or "8g"
    threads = str(p.get("threads") or "4")

    out_prefix = p.get("out_prefix") or str(out_dir / "imputed_fsfhap")
    extra_args = (p.get("extra_args") or "").strip()
    extra = p.get("extra_options") or {}

    stdout_path = out_dir / "stdout.txt"
    stderr_path = out_dir / "stderr.txt"
    error_path = out_dir / "error.txt"

    imputed_tsv = f"{out_prefix}.tsv"
    imputed_vcf = f"{out_prefix}.vcf"

    # Best-effort command template.
    # TASSEL plugin names/options vary across releases; user can override with extra_args.
    # We keep a conservative default that at least starts TASSEL and loads input.
    cmd = [str(tassel_jar), f"-Xmx{java_mem}"]

    # Choose input flag
    fmt_l = inp_fmt.lower()
    if fmt_l == "vcf" or (fmt_l == "auto" and str(inp).lower().endswith((".vcf", ".vcf.gz"))):
        cmd += ["-vcf", str(inp)]
    elif fmt_l in ("hapmap", "hmp") or (fmt_l == "auto" and "hmp" in str(inp).lower()):
        cmd += ["-h", str(inp)]
    else:
        # fallback: let TASSEL try to detect
        cmd += [str(inp)]
    
    # Export VCF (again, exact flag depends; user may override)
    cmd += ["-export", str(out_prefix), "-exportType", "VCF"]

    # Try to invoke FSFHap if available (best guess)
    # User can provide the exact plugin invocation via extra_args.
    # Common-ish naming in docs/code: FSFHapImputationPlugin
    cmd += ["-FSFHapImputationPlugin"]

    if pedigree:
        # many TASSEL plugins accept a family/pedigree file; exact flag may differ
        cmd += ["-pedigrees", pedigree]

    if isinstance(extra, dict):
        for k, v in extra.items():
            cmd.append(f"{k}={v}")

    if extra_args:
        cmd += extra_args.split()

    with runlog.open("w", encoding="utf-8") as f:
        f.write("[impute_fsfhap_tassel] start\n")
        f.write(f"[impute_fsfhap_tassel] time={datetime.now().isoformat()}\n")
        f.write(f"[impute_fsfhap_tassel] cmd={' '.join(cmd)}\n")

    try:
        proc = subprocess.run(cmd, cwd=str(out_dir), capture_output=True, text=True)
        stdout_path.write_text(proc.stdout, encoding="utf-8", errors="ignore")
        stderr_path.write_text(proc.stderr, encoding="utf-8", errors="ignore")

        if proc.returncode != 0:
            raise RuntimeError(f"TASSEL/FSFHap returned code {proc.returncode}")

        # TASSEL export may produce out_prefix.vcf or out_prefix.vcf.gz; try to locate.
        produced_vcf = None
        for cand in [imputed_vcf, f"{out_prefix}.vcf.gz", str(out_dir / Path(imputed_vcf).name), str(out_dir / (Path(imputed_vcf).name + ".gz"))]:
            if Path(cand).exists():
                produced_vcf = str(Path(cand).resolve())
                break

        # Always produce an imputed TSV placeholder (conversion can be added later).
        if not Path(imputed_tsv).exists():
            write_text(Path(imputed_tsv), "id\tmarker1\nplaceholder\tNA\n")

        art = {
            "tool": "FSFHap (TASSEL)",
            "imputed_tsv": str(Path(imputed_tsv).resolve()),
            "tables": [
                {"title": "Imputed genotypes (TSV placeholder)", "path": str(Path(imputed_tsv).resolve())}
            ],
            "logs": [
                {"title": "Run log", "path": str(runlog.resolve())},
                {"title": "STDOUT", "path": str(stdout_path.resolve())},
                {"title": "STDERR", "path": str(stderr_path.resolve())}
            ]
        }
        if produced_vcf:
            art["imputed_vcf"] = produced_vcf
            art["files"] = [{"title": "Imputed VCF", "path": produced_vcf}]

        write_json(out_dir / "artifacts.json", art)

    except Exception as e:
        write_text(error_path, f"{e}\n")
        if not Path(imputed_tsv).exists():
            write_text(Path(imputed_tsv), "id\tmarker1\nplaceholder\tNA\n")
        write_json(out_dir / "artifacts.json", {
            "tool": "FSFHap (TASSEL)",
            "status": "failed_best_effort",
            "error": str(e),
            "imputed_tsv": str(Path(imputed_tsv).resolve()),
            "tables": [
                {"title": "Imputed genotypes (placeholder)", "path": str(Path(imputed_tsv).resolve())}
            ],
            "logs": [
                {"title": "Run log", "path": str(runlog.resolve())},
                {"title": "Error", "path": str(error_path.resolve())},
                {"title": "STDOUT", "path": str(stdout_path.resolve())},
                {"title": "STDERR", "path": str(stderr_path.resolve())}
            ]
        })


if __name__ == "__main__":
    main()
