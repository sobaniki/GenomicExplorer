#!/usr/bin/env python3

import argparse
import json
import os
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

    jar = p.get("jar")
    inp = p.get("input")
    java = p.get("java") or "java"
    java_mem = p.get("java_mem") or "4g"
    threads = str(p.get("threads") or "4")
    out_prefix = p.get("out_prefix") or str(out_dir / "imputed")
    extra_args = p.get("extra_args") or ""
    extra = p.get("extra_options") or {}

    # We keep this plugin "best-effort" because NOISYmputer packaging varies (jar/zip/scripts).
    # The user can adjust the command through extra_args/extra_options.
    cmd = [java, f"-Xmx{java_mem}", "-jar", str(jar)]

    # Common pattern: accept input and output prefix (may differ by distribution)
    # Users can override by passing full extra_args.
    cmd += [str(inp), str(out_prefix), f"threads={threads}"]

    if isinstance(extra, dict):
        # allow key=value flags
        for k, v in extra.items():
            cmd.append(f"{k}={v}")

    if extra_args:
        # split naively (user can quote if needed)
        cmd += extra_args.split()

    with runlog.open("w", encoding="utf-8") as f:
        f.write("[impute_noisyimputer] start\n")
        f.write(f"[impute_noisyimputer] time={datetime.now().isoformat()}\n")
        f.write(f"[impute_noisyimputer] cmd={' '.join(cmd)}\n")

    stdout_path = out_dir / "stdout.txt"
    stderr_path = out_dir / "stderr.txt"
    error_path = out_dir / "error.txt"

    imputed_tsv = f"{out_prefix}.tsv"

    try:
        # run
        proc = subprocess.run(cmd, cwd=str(out_dir), capture_output=True, text=True)
        stdout_path.write_text(proc.stdout, encoding="utf-8", errors="ignore")
        stderr_path.write_text(proc.stderr, encoding="utf-8", errors="ignore")

        if proc.returncode != 0:
            raise RuntimeError(f"NOISYmputer returned code {proc.returncode}")

        # If tool didn't create TSV, create placeholder so GUI doesn't break
        if not Path(imputed_tsv).exists():
            write_text(Path(imputed_tsv), "id\tmarker1\nplaceholder\tNA\n")

        write_json(out_dir / "artifacts.json", {
            "tool": "NOISYmputer",
            "imputed_tsv": str(Path(imputed_tsv).resolve()),
            "tables": [
                {"title": "Imputed genotypes (TSV)", "path": str(Path(imputed_tsv).resolve())}
            ],
            "logs": [
                {"title": "Run log", "path": str(runlog.resolve())},
                {"title": "STDOUT", "path": str(stdout_path.resolve())},
                {"title": "STDERR", "path": str(stderr_path.resolve())}
            ]
        })

    except Exception as e:
        # Always produce placeholders
        write_text(error_path, f"{e}\n")
        if not Path(imputed_tsv).exists():
            write_text(Path(imputed_tsv), "id\tmarker1\nplaceholder\tNA\n")

        write_json(out_dir / "artifacts.json", {
            "tool": "NOISYmputer",
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
