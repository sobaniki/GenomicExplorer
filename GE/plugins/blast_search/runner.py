#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""General BLAST runner for GenomicExplorer.

Supports:
  - Programs: blastn, blastp, blastx, tblastn, tblastx
  - Target:
      * Local DB prefix (preformatted)
      * Local FASTA as DB (auto makeblastdb with caching)
      * NCBI Remote (blast+ -remote)

Outputs:
  - blast_results.tsv (outfmt 6)
  - best_hit_per_query.tsv
  - blast.log.txt
  - artifacts.json
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Iterable, List, Optional, Tuple

import pandas as pd


# ---------------------------
# NCBI taxonomy helper (remote preset)
# ---------------------------

def ncbi_taxid_from_scientific_name(name: str, timeout_sec: int = 20) -> Optional[str]:
    """Resolve a scientific name to a Taxonomy ID (taxid) using NCBI E-utilities.

    Returns the first taxid string if found, else None.
    """
    name = (name or "").strip()
    if not name:
        return None

    # Prefer strict match on Scientific Name.
    term = f'"{name}"[Scientific Name]'
    q = urllib.parse.urlencode({"db": "taxonomy", "term": term, "retmode": "xml"})
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?{q}"
    req = urllib.request.Request(url, headers={"User-Agent": "GenomicExplorer/BLAST-remote"})
    try:
        with urllib.request.urlopen(req, timeout=timeout_sec) as r:
            xml = r.read().decode("utf-8", errors="ignore")
    except Exception:
        return None

    try:
        root = ET.fromstring(xml)
        ids = [e.text.strip() for e in root.findall(".//IdList/Id") if e.text and e.text.strip()]
        return ids[0] if ids else None
    except Exception:
        return None


# ---------------------------
# FASTA helpers (no Biopython)
# ---------------------------

def iter_fasta(path: Path) -> Iterable[Tuple[str, str]]:
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
                w.write(s[i:i + 60] + "\n")


def normalize_id(token: str) -> str:
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
    h = (header or "").strip()
    if not h:
        return []
    first = h.split()[0]
    cands = [first]
    if "|" in first:
        cands.extend([p for p in first.split("|") if p])

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


def guess_fasta_is_nucleotide(path: Path, max_records: int = 3) -> Optional[bool]:
    """Very small heuristic to warn users; returns None if unknown."""
    try:
        checked = 0
        for _h, seq in iter_fasta(path):
            s = (seq or "").upper().replace("-", "")
            if not s:
                continue
            checked += 1
            # fraction of A/C/G/T/N
            good = sum((ch in {"A", "C", "G", "T", "N"}) for ch in s)
            frac = good / max(1, len(s))
            if frac >= 0.90:
                return True
            if frac <= 0.60:
                return False
            if checked >= max_records:
                break
    except Exception:
        return None
    return None


# ---------------------------
# Gene list parsing (optional)
# ---------------------------

def load_gene_ids(path: Path) -> List[str]:
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

    found: List[Tuple[str, str]] = []
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
# DB cache + tool checks
# ---------------------------

def sha1_text(s: str) -> str:
    return hashlib.sha1(s.encode("utf-8")).hexdigest()


def ensure_exe(name: str) -> str:
    exe = shutil.which(name)
    if not exe:
        raise RuntimeError(f"{name} not found in PATH. Install BLAST+ (e.g., conda install -c bioconda blast).")
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


def has_blast_db(prefix: Path, dbtype: str) -> bool:
    if dbtype == "prot":
        return (prefix.with_suffix(".pin").exists()
                and prefix.with_suffix(".phr").exists()
                and prefix.with_suffix(".psq").exists())
    # nucl
    return (prefix.with_suffix(".nin").exists()
            and prefix.with_suffix(".nhr").exists()
            and prefix.with_suffix(".nsq").exists())


def program_to_dbtype(program: str) -> str:
    p = (program or "").strip().lower()
    # program dictates subject database type
    #   - blastn, tblastn, tblastx => nucl
    #   - blastp, blastx => prot
    if p in {"blastp", "blastx"}:
        return "prot"
    if p in {"blastn", "tblastn", "tblastx"}:
        return "nucl"
    raise ValueError(f"Unsupported program: {program}")


def prepare_db_prefix(db_fasta: Path,
                      cache_dir: Path,
                      cache_mode: str,
                      force: bool,
                      dbtype: str,
                      log_lines: List[str]) -> Path:
    cache_dir.mkdir(parents=True, exist_ok=True)
    key = compute_cache_key(db_fasta, mode=cache_mode)
    dbdir = cache_dir / f"{dbtype}_{key}"
    dbdir.mkdir(parents=True, exist_ok=True)
    prefix = dbdir / "db"

    meta = {"db_fasta": str(db_fasta.resolve()), "cache_mode": cache_mode, "key": key, "dbtype": dbtype}
    (dbdir / "meta.json").write_text(json.dumps(meta, ensure_ascii=False, indent=2), encoding="utf-8")

    if (not force) and has_blast_db(prefix, dbtype=dbtype):
        log_lines.append(f"[db] reuse cached blastdb: {prefix} (dbtype={dbtype})")
        return prefix

    makeblastdb = ensure_exe("makeblastdb")
    cmd = [makeblastdb, "-in", str(db_fasta), "-dbtype", dbtype, "-out", str(prefix)]
    log_lines.append("[db] build blastdb: " + " ".join(cmd))
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    log_lines.append(proc.stdout or "")
    log_lines.append(proc.stderr or "")
    if proc.returncode != 0 or not has_blast_db(prefix, dbtype=dbtype):
        raise RuntimeError(f"makeblastdb failed (returncode={proc.returncode}). See blast.log.txt")
    return prefix


# ---------------------------
# BLAST run + summarize
# ---------------------------

def run_blast(program: str,
              query_fasta: Path,
              out_path: Path,
              outfmt: str,
              mode: str,
              db_prefix: Optional[Path],
              db_remote: Optional[str],
              task: Optional[str],
              evalue: float,
              max_target_seqs: int,
              num_threads: int,
              entrez_query: Optional[str],
              extra_args: List[str],
              log_lines: List[str]) -> None:
    program = program.strip().lower()
    exe = ensure_exe(program)

    cmd: List[str] = [
        exe,
        "-query", str(query_fasta),
        "-outfmt", outfmt,
        "-out", str(out_path),
        "-evalue", str(evalue),
        "-max_target_seqs", str(max_target_seqs),
    ]

    # Optional task only for blastp/blastn (blastx/tblast* don't accept -task in many builds)
    if task and program in {"blastp", "blastn"}:
        cmd += ["-task", str(task)]

    if mode == "remote_ncbi":
        if not db_remote:
            raise ValueError("remote_ncbi requires db_remote")
        cmd += ["-db", str(db_remote), "-remote"]
        if entrez_query and str(entrez_query).strip():
            cmd += ["-entrez_query", str(entrez_query).strip()]
        # threads are not applicable to remote
    else:
        if not db_prefix:
            raise ValueError("local mode requires db_prefix")
        cmd += ["-db", str(db_prefix)]
        if num_threads and int(num_threads) > 0:
            cmd += ["-num_threads", str(int(num_threads))]

    cmd += (extra_args or [])

    log_lines.append(f"[{program}] run: " + " ".join(cmd))
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    log_lines.append(proc.stdout or "")
    log_lines.append(proc.stderr or "")
    if proc.returncode != 0:
        raise RuntimeError(f"{program} failed (returncode={proc.returncode}). See blast.log.txt")


def build_html_report(program: str,
                      archive_path: Path,
                      out_html: Path,
                      log_lines: List[str]) -> None:
    """Format a BLAST archive (outfmt 11) into a readable HTML report."""
    blast_formatter = ensure_exe("blast_formatter")
    cmd = [blast_formatter, "-archive", str(archive_path), "-out", str(out_html), "-html"]
    log_lines.append(f"[blast_formatter] run: " + " ".join(cmd))
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    log_lines.append(proc.stdout or "")
    log_lines.append(proc.stderr or "")
    if proc.returncode != 0:
        raise RuntimeError(f"blast_formatter failed (returncode={proc.returncode}). See blast.log.txt")


def best_hit_per_query(blast_tsv: Path, out_best: Path) -> None:
    cols = [
        "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore", "stitle"
    ]
    if not blast_tsv.exists() or blast_tsv.stat().st_size == 0:
        out_best.write_text("", encoding="utf-8")
        return
    df = pd.read_csv(blast_tsv, sep="\t", header=None, names=cols, dtype=str)

    df["bitscore_num"] = pd.to_numeric(df["bitscore"], errors="coerce")
    df["evalue_num"] = pd.to_numeric(df["evalue"], errors="coerce")
    df = df.sort_values(["qseqid", "bitscore_num", "evalue_num"], ascending=[True, False, True])

    best = df.groupby("qseqid", as_index=False).head(1).drop(columns=["bitscore_num", "evalue_num"])
    best.to_csv(out_best, sep="\t", index=False)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--params", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    params = json.loads(Path(args.params).read_text(encoding="utf-8"))
    out_dir = Path(args.out).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    # -----------------
    # Core parameters
    # -----------------
    program = str(params.get("program", "blastp")).strip().lower()
    mode = str(params.get("mode", "local_db")).strip().lower()  # local_db | local_fasta_as_db | remote_ncbi

    # Query inputs
    query_fasta = params.get("query_fasta")  # optional if gene_list+proteome provided
    gene_list = params.get("gene_list")
    proteome_fasta = params.get("proteome_fasta")

    # DB inputs
    db_fasta = params.get("db_fasta")
    db_prefix_in = params.get("db_prefix")
    db_remote = params.get("db_remote")
    entrez_query = params.get("entrez_query")
    organism = params.get("organism")  # optional scientific name (remote preset)

    # Options
    task = params.get("task")  # optional
    evalue = float(params.get("evalue", 1e-5))
    max_target_seqs = int(params.get("max_target_seqs", 5))
    num_threads = int(params.get("num_threads", 4))

    # Optional: create an HTML report (pairwise alignments) using blast_formatter.
    # This will run an additional BLAST pass producing an archive (outfmt 11), then
    # convert it to HTML. Kept optional because it can be slow for very large queries.
    make_html_report = str(params.get("make_html_report") or "FALSE").upper() == "TRUE"

    # DB cache (local_fasta_as_db only)
    cache_dir = Path(params.get("cache_dir") or (out_dir.parent / "cache" / "blastdb")).resolve()
    cache_mode = str(params.get("cache_mode") or "mtime_size")
    force_makeblastdb = str(params.get("force_makeblastdb") or "FALSE").upper() == "TRUE"

    extra_args = params.get("extra_args") or []
    if isinstance(extra_args, str):
        extra_args = extra_args.split()

    log_lines: List[str] = []
    log_lines.append(f"[params] program={program} mode={mode}")

    if organism and str(organism).strip():
        log_lines.append(f"[params] organism={str(organism).strip()}")

    # Validate supported programs
    if program not in {"blastn", "blastp", "blastx", "tblastn", "tblastx"}:
        raise ValueError(f"Unsupported program: {program}")

    # -----------------
    # Build / choose query
    # -----------------
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

    # small warning heuristic
    q_is_nucl = guess_fasta_is_nucleotide(query_path)
    if q_is_nucl is not None:
        log_lines.append(f"[query] heuristic nucleotide={q_is_nucl} ({query_path.name})")

    # -----------------
    # Prepare DB (local) or remote args
    # -----------------
    db_prefix: Optional[Path] = None
    dbtype = program_to_dbtype(program)
    log_lines.append(f"[db] expected dbtype={dbtype} (from program={program})")

    if mode == "remote_ncbi":
        # remote requires only db_remote
        if not db_remote or not str(db_remote).strip():
            raise ValueError("mode=remote_ncbi requires db_remote (e.g., nr, nt, refseq_protein)")
        # make sure exe exists now (better error)
        ensure_exe(program)

        # Preset: if organism (scientific name) is provided and entrez_query is empty,
        # resolve taxid via NCBI taxonomy and build txid... [ORGN].
        if (not (entrez_query and str(entrez_query).strip())) and (organism and str(organism).strip()):
            taxid = ncbi_taxid_from_scientific_name(str(organism).strip())
            if taxid:
                entrez_query = f"txid{taxid}[ORGN]"
                log_lines.append(f"[remote] resolved organism -> {entrez_query}")
            else:
                log_lines.append("[remote][WARN] failed to resolve organism to taxid; running without organism filter")
    else:
        # local_db (prefix) or local_fasta_as_db
        if db_prefix_in and str(db_prefix_in).strip():
            db_prefix = Path(db_prefix_in).resolve()
            if not has_blast_db(db_prefix, dbtype=dbtype):
                suffix = ".pin/.phr/.psq" if dbtype == "prot" else ".nin/.nhr/.nsq"
                raise FileNotFoundError(f"db_prefix missing {suffix}: {db_prefix}")
            log_lines.append(f"[db] use preformatted db_prefix: {db_prefix}")
        else:
            if not db_fasta:
                raise ValueError("Provide db_fasta or db_prefix for local mode")
            if mode not in {"local_fasta_as_db", "local_db"}:
                raise ValueError(f"Unknown mode: {mode}")
            # local_db without prefix behaves same as local_fasta_as_db
            db_prefix = prepare_db_prefix(Path(db_fasta).resolve(), cache_dir, cache_mode, force_makeblastdb, dbtype, log_lines)

    # -----------------
    # Run
    # -----------------
    out_tsv = out_dir / "blast_results.tsv"
    outfmt_tsv = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle"
    run_blast(
        program=program,
        query_fasta=query_path,
        out_path=out_tsv,
        outfmt=outfmt_tsv,
        mode=mode,
        db_prefix=db_prefix,
        db_remote=str(db_remote).strip() if db_remote else None,
        task=str(task).strip() if task else None,
        evalue=evalue,
        max_target_seqs=max_target_seqs,
        num_threads=num_threads,
        entrez_query=str(entrez_query).strip() if entrez_query else None,
        extra_args=extra_args,
        log_lines=log_lines,
    )

    out_best = out_dir / "best_hit_per_query.tsv"
    best_hit_per_query(out_tsv, out_best)

    # Optional HTML report (pairwise alignments)
    out_archive = None
    out_html = None
    if make_html_report:
        try:
            out_archive = out_dir / "blast_archive.asn"
            out_html = out_dir / "blast_report.html"
            # Archive format recommended by BLAST+ for blast_formatter.
            run_blast(
                program=program,
                query_fasta=query_path,
                out_path=out_archive,
                outfmt="11",
                mode=mode,
                db_prefix=db_prefix,
                db_remote=str(db_remote).strip() if db_remote else None,
                task=str(task).strip() if task else None,
                evalue=evalue,
                max_target_seqs=max_target_seqs,
                num_threads=num_threads,
                entrez_query=str(entrez_query).strip() if entrez_query else None,
                extra_args=extra_args,
                log_lines=log_lines,
            )
            blast_formatter = ensure_exe("blast_formatter")
            cmd = [blast_formatter, "-archive", str(out_archive), "-out", str(out_html), "-html"]
            log_lines.append("[blast_formatter] run: " + " ".join(cmd))
            proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
            log_lines.append(proc.stdout or "")
            log_lines.append(proc.stderr or "")
            if proc.returncode != 0 or (not out_html.exists()):
                log_lines.append(f"[WARN] blast_formatter failed (returncode={proc.returncode}); report not generated")
                out_html = None
        except Exception as e:
            log_lines.append(f"[WARN] failed to generate HTML report: {e}")
            out_archive = None
            out_html = None

    (out_dir / "blast.log.txt").write_text("\n".join([l for l in log_lines if l is not None]) + "\n", encoding="utf-8")

    artifacts = {
        "blast_results_tsv": str(out_tsv),
        "best_hit_per_query_tsv": str(out_best),
        "blast_report_html": str(out_html) if out_html else "",
        "blast_archive_asn": str(out_archive) if out_archive else "",
        "log": str(out_dir / "blast.log.txt"),
        "query_used": str(query_path),
        "program": program,
        "mode": mode,
        "db_prefix_used": str(db_prefix) if db_prefix else "",
        "db_remote": str(db_remote or ""),
        "entrez_query": str(entrez_query or ""),
        "organism": str(organism or ""),
        "cache_dir": str(cache_dir),
        "cache_mode": cache_mode,
        "make_html_report": bool(make_html_report),
    }
    (out_dir / "artifacts.json").write_text(json.dumps(artifacts, ensure_ascii=False, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
