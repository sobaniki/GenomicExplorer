#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
from pathlib import Path
from typing import Iterable, List, Optional, Tuple
import hashlib

import pandas as pd


# ---------------------------
# FASTA helpers (no Biopython)
# ---------------------------

def iter_fasta(path: Path) -> Iterable[Tuple[str, str]]:
    """Yield (header, sequence) from a FASTA file."""
    header = None
    seq_chunks: List[str] = []
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks)
                header = line[1:].strip()
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())
        if header is not None:
            yield header, "".join(seq_chunks)


def write_fasta(records: Iterable[Tuple[str, str]], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as w:
        for h, s in records:
            w.write(f">{h}\n")
            for i in range(0, len(s), 60):
                w.write(s[i:i+60] + "\n")


def normalize_id(token: str) -> str:
    """Normalize an identifier: strip whitespace and trailing version like '.1'."""
    token = (token or "").strip()
    if not token:
        return ""
    token = token.split()[0]
    if "." in token:
        left, right = token.rsplit(".", 1)
        if right.isdigit():
            token = left
    return token


def header_id_candidates(header: str) -> List[str]:
    """Get possible ID candidates from FASTA header."""
    h = (header or "").strip()
    if not h:
        return []
    first = h.split()[0]
    cands = [first]
    if "|" in first:
        cands.extend([p for p in first.split("|") if p])

    # add normalized variants
    out: List[str] = []
    seen = set()
    for c in cands:
        for v in (c, normalize_id(c)):
            v = (v or "").strip()
            if not v or v in seen:
                continue
            seen.add(v)
            out.append(v)
    return out


# ---------------------------
# Gene list parsing (optional)
# ---------------------------

def load_gene_ids(path: Path) -> List[str]:
    """Load gene IDs from TSV/CSV. Uses gene/id-like column if present, else first column."""
    try:
        df = pd.read_csv(path, sep="\t", dtype=str)
    except Exception:
        try:
            df = pd.read_csv(path, sep=",", dtype=str)
        except Exception:
            ids = []
            for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
                if not line.strip() or line.startswith("#"):
                    continue
                ids.append(line.strip().split()[0])
            return sorted(set(ids))

    if df.shape[1] == 0:
        return []

    cols = [c.lower() for c in df.columns]
    pick = None
    for key in ["gene_id", "gene", "id", "geneid", "protein_id", "protein", "locus"]:
        if key in cols:
            pick = df.columns[cols.index(key)]
            break
    if pick is None:
        pick = df.columns[0]

    ids = [str(x).strip() for x in df[pick].dropna().tolist()]
    ids = [x for x in ids if x and x.lower() not in {"nan", "none"}]
    return sorted(set(ids))


def build_query_from_gene_list(gene_list: Path, proteome_fasta: Path, out_fasta: Path, log_lines: List[str]) -> Path:
    gene_ids = load_gene_ids(gene_list)
    want = set([normalize_id(x) for x in gene_ids if x])
    if not want:
        raise ValueError(f"No gene IDs found in: {gene_list}")

    found = []
    miss = set(want)

    for header, seq in iter_fasta(proteome_fasta):
        cands = header_id_candidates(header)
        hit = None
        for c in cands:
            if normalize_id(c) in want:
                hit = normalize_id(c)
                break
        if hit is not None:
            found.append((header, seq))
            miss.discard(hit)

    write_fasta(found, out_fasta)

    log_lines.append(f"[query_build] gene_ids: {len(want)}")
    log_lines.append(f"[query_build] found sequences: {len(found)}")
    if miss:
        miss_path = out_fasta.parent / "missing_gene_ids.txt"
        miss_path.write_text("\n".join(sorted(miss)) + "\n", encoding="utf-8")
        log_lines.append(f"[query_build] missing ids written: {miss_path.name} ({len(miss)})")
    return out_fasta


# ---------------------------
# BLAST DB cache
# ---------------------------

def sha1_text(s: str) -> str:
    return hashlib.sha1(s.encode("utf-8")).hexdigest()


def ensure_makeblastdb() -> str:
    exe = shutil.which("makeblastdb")
    if not exe:
        raise RuntimeError("makeblastdb not found in PATH. Install BLAST+ (e.g., conda install -c bioconda blast).")
    return exe


def ensure_blastp() -> str:
    exe = shutil.which("blastp")
    if not exe:
        raise RuntimeError("blastp not found in PATH. Install BLAST+ (e.g., conda install -c bioconda blast).")
    return exe


def compute_cache_key(db_fasta: Path, mode: str = "mtime_size") -> str:
    db_fasta = db_fasta.resolve()
    st = db_fasta.stat()
    if mode == "md5":
        h = hashlib.md5()
        with db_fasta.open("rb") as f:
            for chunk in iter(lambda: f.read(1024 * 1024), b""):
                h.update(chunk)
        sig = f"{db_fasta}|{h.hexdigest()}"
    else:
        sig = f"{db_fasta}|{int(st.st_mtime)}|{st.st_size}"
    return sha1_text(sig)[:16]


def has_blast_db(prefix: Path) -> bool:
    return (prefix.with_suffix(".pin").exists()
            and prefix.with_suffix(".phr").exists()
            and prefix.with_suffix(".psq").exists())


def prepare_db_prefix(db_fasta: Path,
                      cache_dir: Path,
                      cache_mode: str,
                      force: bool,
                      log_lines: List[str]) -> Path:
    cache_dir.mkdir(parents=True, exist_ok=True)
    key = compute_cache_key(db_fasta, mode=cache_mode)
    dbdir = cache_dir / key
    dbdir.mkdir(parents=True, exist_ok=True)
    prefix = dbdir / "db"

    meta = {"db_fasta": str(db_fasta.resolve()), "cache_mode": cache_mode, "key": key}
    (dbdir / "meta.json").write_text(json.dumps(meta, ensure_ascii=False, indent=2), encoding="utf-8")

    if (not force) and has_blast_db(prefix):
        log_lines.append(f"[db] reuse cached blastdb: {prefix}")
        return prefix

    makeblastdb = ensure_makeblastdb()
    cmd = [makeblastdb, "-in", str(db_fasta), "-dbtype", "prot", "-out", str(prefix)]
    log_lines.append("[db] build blastdb: " + " ".join(cmd))
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    log_lines.append(proc.stdout or "")
    log_lines.append(proc.stderr or "")
    if proc.returncode != 0 or not has_blast_db(prefix):
        raise RuntimeError(f"makeblastdb failed (returncode={proc.returncode}). See blast.log.txt")
    return prefix


# ---------------------------
# BLAST run + summarize
# ---------------------------

def run_blastp(query_fasta: Path,
               db_prefix: Path,
               out_tsv: Path,
               task: str,
               evalue: float,
               max_target_seqs: int,
               num_threads: int,
               extra_args: List[str],
               log_lines: List[str]) -> None:
    blastp = ensure_blastp()
    outfmt = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle"
    cmd = [
        blastp,
        "-query", str(query_fasta),
        "-db", str(db_prefix),
        "-task", task,
        "-evalue", str(evalue),
        "-max_target_seqs", str(max_target_seqs),
        "-num_threads", str(num_threads),
        "-outfmt", outfmt,
        "-out", str(out_tsv),
    ] + (extra_args or [])
    log_lines.append("[blastp] run: " + " ".join(cmd))
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    log_lines.append(proc.stdout or "")
    log_lines.append(proc.stderr or "")
    if proc.returncode != 0:
        raise RuntimeError(f"blastp failed (returncode={proc.returncode}). See blast.log.txt")


def best_hit_per_query(blast_tsv: Path, out_best: Path) -> None:
    cols = ["qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","stitle"]
    if not blast_tsv.exists() or blast_tsv.stat().st_size == 0:
        out_best.write_text("", encoding="utf-8")
        return
    df = pd.read_csv(blast_tsv, sep="\t", header=None, names=cols, dtype=str)

    df["bitscore_num"] = pd.to_numeric(df["bitscore"], errors="coerce")
    df["evalue_num"] = pd.to_numeric(df["evalue"], errors="coerce")
    df = df.sort_values(["qseqid","bitscore_num","evalue_num"], ascending=[True, False, True])

    best = df.groupby("qseqid", as_index=False).head(1).drop(columns=["bitscore_num","evalue_num"])
    best.to_csv(out_best, sep="\t", index=False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--params", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    params = json.loads(Path(args.params).read_text(encoding="utf-8"))
    out_dir = Path(args.out).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    # Inputs
    query_fasta = params.get("query_fasta")      # optional if gene_list+proteome provided
    gene_list = params.get("gene_list")          # optional
    proteome_fasta = params.get("proteome_fasta")  # optional

    db_fasta = params.get("db_fasta")            # required unless db_prefix provided
    db_prefix_in = params.get("db_prefix")       # optional (preformatted BLAST DB prefix)

    task = str(params.get("task", "blastp"))
    evalue = float(params.get("evalue", 1e-5))
    max_target_seqs = int(params.get("max_target_seqs", 5))
    num_threads = int(params.get("num_threads", 4))

    # DB options (A/B)
    cache_dir = Path(params.get("cache_dir") or (out_dir.parent / "cache" / "blastdb")).resolve()
    cache_mode = str(params.get("cache_mode") or "mtime_size")
    force_makeblastdb = str(params.get("force_makeblastdb") or "FALSE").upper() == "TRUE"

    extra_args = params.get("extra_args") or []
    if isinstance(extra_args, str):
        extra_args = extra_args.split()

    log_lines: List[str] = []

    # Build / choose query (C)
    if query_fasta and str(query_fasta).strip():
        query_path = Path(query_fasta).resolve()
    else:
        if gene_list and proteome_fasta:
            query_path = out_dir / "query_from_gene_list.faa"
            build_query_from_gene_list(Path(gene_list), Path(proteome_fasta), query_path, log_lines)
        else:
            raise ValueError("Provide either query_fasta, or (gene_list + proteome_fasta).")

    if not query_path.exists():
        raise FileNotFoundError(f"Query FASTA not found: {query_path}")

    # Prepare DB prefix (A/B)
    if db_prefix_in and str(db_prefix_in).strip():
        db_prefix = Path(db_prefix_in).resolve()
        if not has_blast_db(db_prefix):
            raise FileNotFoundError(f"db_prefix missing .pin/.phr/.psq: {db_prefix}")
        log_lines.append(f"[db] use preformatted db_prefix: {db_prefix}")
    else:
        if not db_fasta:
            raise ValueError("Provide db_fasta (Protein FASTA) or db_prefix (preformatted BLAST DB).")
        db_prefix = prepare_db_prefix(Path(db_fasta).resolve(), cache_dir, cache_mode, force_makeblastdb, log_lines)

    # Run BLASTP
    out_tsv = out_dir / "blast_results.tsv"
    run_blastp(query_path, db_prefix, out_tsv, task, evalue, max_target_seqs, num_threads, extra_args, log_lines)

    # Best hit per query
    out_best = out_dir / "best_hit_per_query.tsv"
    best_hit_per_query(out_tsv, out_best)

    # Write log + artifacts
    (out_dir / "blast.log.txt").write_text("\n".join([l for l in log_lines if l is not None]) + "\n", encoding="utf-8")

    artifacts = {
        "blast_results_tsv": str(out_tsv),
        "best_hit_per_query_tsv": str(out_best),
        "log": str(out_dir / "blast.log.txt"),
        "query_used": str(query_path),
        "db_prefix_used": str(db_prefix),
        "cache_dir": str(cache_dir),
    }
    (out_dir / "artifacts.json").write_text(json.dumps(artifacts, ensure_ascii=False, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
