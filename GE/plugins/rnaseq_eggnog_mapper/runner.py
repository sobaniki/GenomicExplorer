#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""EggNOG-mapper runner (Phase B minimal).

Inputs:
  - protein FASTA
  - EggNOG data_dir

Runs:
  emapper.py (local) -> <prefix>.emapper.annotations

Outputs:
  - emapper annotations (as produced by eggnog-mapper)
  - gene2go.tsv (2 columns: gene_id, go_id)
  - artifacts.json

Notes:
  - This plugin does NOT build OrgDb; pass gene2go.tsv to rnaseq_build_orgdb (Phase A).
  - gene_id is extracted as the first token of the FASTA header (up to first whitespace).
"""

from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Tuple


GO_RE = re.compile(r"GO:\d+")


def _ensure_emapper(emapper_path: str) -> str:
    exe = shutil.which(emapper_path) if Path(emapper_path).name == emapper_path else emapper_path
    if not exe:
        raise RuntimeError(
            "emapper.py not found in PATH. Install eggnog-mapper (e.g., conda install -c bioconda eggnog-mapper) "
            "or set params.emapper_path to the full path."
        )
    return exe

def _ensure_download_script() -> str:
    exe = shutil.which("download_eggnog_data.py")
    if not exe:
        raise RuntimeError(
              "download_eggnog_data.py not found in PATH. "
              "Your eggnog-mapper installation may be incomplete or too old."
        )
    return exe

  
def _has_required_db_files(data_dir: Path) -> Tuple[bool, List[str]]:
    missing: List[str] = []
    eggnog_db = data_dir / "eggnog.db"
    if not eggnog_db.exists():
      missing.append("eggnog.db")

    taxa_hits = list(data_dir.glob("eggnog.taxa*"))
    if not taxa_hits:
      missing.append("eggnog.taxa*")
  
    dmnd = data_dir / "eggnog_proteins.dmnd"
    if not dmnd.exists():
      missing.append("eggnog_proteins.dmnd")
    
    return (len(missing) == 0), missing
    
      
def _auto_download_db(data_dir: Path, download_diamond: bool = True) -> None:
    data_dir.mkdir(parents=True, exist_ok=True)
    dl = _ensure_download_script()
    cmd = [dl, "-y", "--data_dir", str(data_dir)]
    if download_diamond:
      cmd.insert(1, "-D")  # keep close to common usage: download_eggnog_data.py -y -D --data_dir <dir>
    
    proc = subprocess.run(cmd, cwd=str(data_dir), capture_output=True, text=True)
    (data_dir / "download_eggnog_data.stdout.txt").write_text(proc.stdout or "", encoding="utf-8", errors="ignore")
    (data_dir / "download_eggnog_data.stderr.txt").write_text(proc.stderr or "", encoding="utf-8", errors="ignore")
    if proc.returncode != 0:
      raise RuntimeError(
      f"Auto-download failed (code={proc.returncode}). "
      f"See download_eggnog_data.stderr.txt in {data_dir}.\n"
      f"cmd={' '.join(cmd)}"
    )

def _read_params(path: Path) -> Dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _pick_header_id(qseqid: str) -> str:
    # EggNOG-mapper writes query as FASTA header; safest is first whitespace-delimited token.
    token = (qseqid or "").strip().split()[0]
    if not token:
        return ""
    return token


def _find_gos_column(cols: List[str]) -> Optional[int]:
    lower = [c.lower() for c in cols]
    for key in ("gos", "go", "go_terms", "go term", "go terms"):
        if key in lower:
            return lower.index(key)
    return None


def _parse_annotations_to_gene2go(ann_path: Path, out_tsv: Path) -> Tuple[int, int]:
    """Return (n_queries, n_pairs)."""
    n_queries = 0
    n_pairs = 0
    pairs: List[Tuple[str, str]] = []

    with ann_path.open("r", encoding="utf-8", errors="ignore") as f:
        header: Optional[List[str]] = None
        go_col: Optional[int] = None
        q_col: Optional[int] = None

        for line in f:
            raw = line.rstrip("\n")
            if not raw.strip():
                continue
            # eggNOG-mapper frequently prefixes header/comments with '#'
            if header is None:
                if raw.startswith("##") or raw.startswith("!"):
                    continue
                if raw.startswith("#"):
                    raw_hdr = raw.lstrip("#")
                    parts = raw_hdr.split("\t")
                else:
                    parts = raw.split("\t")
                header = parts
                lower = [c.lower() for c in header]
                if "query" in lower:
                    q_col = lower.index("query")
                elif "#query" in lower:
                    q_col = lower.index("#query")
                else:
                    q_col = 0
                go_col = _find_gos_column(header)
                if go_col is None:
                    raise RuntimeError(
                        "Could not locate 'GOs' column in annotations header. "
                        "Please ensure eggnog-mapper outputs a header (TSV) and contains a 'GOs' column."
                    )
                continue

            # after header
            if raw.startswith("#") or raw.startswith("!"):
                continue
            parts = raw.split("\t")

            if q_col is None or go_col is None:
                raise RuntimeError("Internal parsing error: q_col/go_col not set")
            if len(parts) <= max(q_col, go_col):
                continue
            q = _pick_header_id(parts[q_col])
            if not q:
                continue
            n_queries += 1
            gos_raw = parts[go_col].strip()
            if not gos_raw or gos_raw in {"-", "NA", "NaN"}:
                continue
            gos = GO_RE.findall(gos_raw)
            if not gos:
                continue
            for go in sorted(set(gos)):
                pairs.append((q, go))
                n_pairs += 1

    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with out_tsv.open("w", encoding="utf-8") as w:
        w.write("gene_id\tgo_id\n")
        for g, go in pairs:
            w.write(f"{g}\t{go}\n")
    return n_queries, n_pairs


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--params", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    params_path = Path(args.params).resolve()
    out_dir = Path(args.out).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    params = _read_params(params_path)
    protein_fasta = Path(str(params.get("protein_fasta", ""))).expanduser().resolve()
    data_dir = Path(str(params.get("data_dir", ""))).expanduser().resolve()
    output_prefix = str(params.get("output_prefix") or "eggnog").strip() or "eggnog"
    cpu = int(params.get("cpu") or 1)
    tax_scope = str(params.get("tax_scope") or "").strip()
    emapper_path = str(params.get("emapper_path") or "emapper.py").strip() or "emapper.py"
    extra_args = params.get("extra_args")
    auto_download = bool(params.get("auto_download", True))

    if not protein_fasta.exists():
        raise RuntimeError(f"protein_fasta not found: {protein_fasta}")
    if not data_dir.exists() or not data_dir.is_dir():
        if auto_download:
           data_dir.mkdir(parents=True, exist_ok=True)
        else:
           raise RuntimeError(f"data_dir (EggNOG database dir) not found: {data_dir}")

    ok, missing = _has_required_db_files(data_dir)
    if (not ok) and auto_download:
        # Try to fetch minimal DB set + DIAMOND DB in one shot.
        _auto_download_db(data_dir=data_dir, download_diamond=True)
        ok2, missing2 = _has_required_db_files(data_dir)
        if not ok2:
            raise RuntimeError(
                "Auto-download finished but required DB files are still missing: "
                + ", ".join(missing2)
                + f". Please inspect {data_dir} and download logs."
            )
    elif (not ok) and (not auto_download):
        raise RuntimeError(
            "EggNOG DB files missing in data_dir: " + ", ".join(missing)
            + ". Enable auto_download or run download_eggnog_data.py manually."
        )

    emapper_exe = _ensure_emapper(emapper_path)

    cmd: List[str] = [
        emapper_exe,
        "-i", str(protein_fasta),
        "--itype", "proteins",
        "-o", output_prefix,
        "--output_dir", str(out_dir),
        "--data_dir", str(data_dir),
        "--cpu", str(cpu),
    ]
    if tax_scope:
        cmd += ["--tax_scope", tax_scope]

    if isinstance(extra_args, str) and extra_args.strip():
        cmd += extra_args.strip().split()
    elif isinstance(extra_args, list):
        cmd += [str(x) for x in extra_args if str(x).strip()]

    proc = subprocess.run(cmd, cwd=str(out_dir), capture_output=True, text=True)
    (out_dir / "emapper.stdout.txt").write_text(proc.stdout or "", encoding="utf-8", errors="ignore")
    (out_dir / "emapper.stderr.txt").write_text(proc.stderr or "", encoding="utf-8", errors="ignore")
    if proc.returncode != 0:
        raise RuntimeError(
            f"eggnog-mapper failed (code={proc.returncode}). See emapper.stderr.txt for details.\n"
            f"cmd={' '.join(cmd)}"
        )

    candidates = [
        out_dir / f"{output_prefix}.emapper.annotations",
        out_dir / f"{output_prefix}.emapper.annotations.tsv",
    ]
    ann_path = None
    for p in candidates:
        if p.exists() and p.stat().st_size > 0:
            ann_path = p
            break
    if ann_path is None:
        gl = list(out_dir.glob(f"{output_prefix}.emapper.*annotations*"))
        gl = [p for p in gl if p.is_file() and p.stat().st_size > 0]
        if gl:
            ann_path = sorted(gl)[0]
    if ann_path is None:
        raise RuntimeError(
            "Could not find eggnog-mapper annotations file in out_dir. "
            "Expected '<prefix>.emapper.annotations'."
        )

    gene2go_path = out_dir / "gene2go.tsv"
    n_queries, n_pairs = _parse_annotations_to_gene2go(ann_path, gene2go_path)
    if gene2go_path.stat().st_size == 0 or n_pairs == 0:
        raise RuntimeError(
            "gene2go.tsv is empty (no GO terms extracted). "
            "Check the input FASTA and eggnog-mapper settings, and verify the 'GOs' column exists."
        )

    artifacts = {
        "plugin": "rnaseq_eggnog_mapper",
        "protein_fasta": str(protein_fasta),
        "data_dir": str(data_dir),
        "output_prefix": output_prefix,
        "annotations": str(ann_path),
        "gene2go_tsv": str(gene2go_path),
        "n_queries": int(n_queries),
        "n_gene2go_pairs": int(n_pairs),
        "auto_download": bool(auto_download),
        "next_step": "Use gene2go.tsv in 'GO DB builder (Phase A)' to build an OrgDb package.",
        }
    (out_dir / "artifacts.json").write_text(json.dumps(artifacts, ensure_ascii=False, indent=2), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
