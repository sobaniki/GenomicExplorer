#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""annot_gene_pred_abinitio

Gene prediction around GWAS/QTL peaks.

Outputs (written under --out):
  - predicted.gff3
  - predicted_near_peaks.tsv

The GUI (Integrate → GenePred) supplies:
  genome_fasta, peaks_tsv, flank_bp (+/-), tool, threads,
  plus tool-specific settings for AUGUSTUS/SNAP.

This runner additionally supports BRAKER3 (main path) via "tool=braker3".
BRAKER3 is invoked as an external command (default: braker.pl). This runner does
NOT install/configure BRAKER3; it only orchestrates inputs/outputs.

Design note:
  To keep runtime practical for peak-centric exploration, we extract small
  FASTA sequences around each peak and run predictors on those regions.
  Coordinates are then mapped back to the original genome.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple


# -------------------------
# Utilities
# -------------------------


def _sniff_delimiter(path: Path) -> str:
    # best-effort sniff for tsv/csv
    with path.open("r", encoding="utf-8", errors="ignore", newline="") as f:
        head = "".join([f.readline() for _ in range(5)])
    if "\t" in head:
        return "\t"
    if "," in head:
        return ","
    return "\t"


def _pick_col(cols: List[str], candidates: List[str]) -> Optional[str]:
    low = {c.lower(): c for c in cols}
    for k in candidates:
        if k.lower() in low:
            return low[k.lower()]
    return None


def _safe_int(x: Any) -> Optional[int]:
    try:
        if x is None:
            return None
        s = str(x).strip()
        if not s:
            return None
        return int(float(s))
    except Exception:
        return None


def _safe_float(x: Any) -> Optional[float]:
    try:
        if x is None:
            return None
        s = str(x).strip()
        if not s:
            return None
        return float(s)
    except Exception:
        return None


def _read_rows(path: Path) -> Tuple[List[str], List[Dict[str, Any]]]:
    delim = _sniff_delimiter(path)
    with path.open("r", encoding="utf-8", errors="ignore", newline="") as f:
        reader = csv.DictReader(f, delimiter=delim)
        cols = reader.fieldnames or []
        rows = [r for r in reader]
    return cols, rows


def _parse_attrs(attr: str) -> Dict[str, str]:
    """Parse GFF3 (k=v;...) or GTF (k "v"; ...) attributes."""
    attr = (attr or "").strip()
    d: Dict[str, str] = {}
    if not attr:
        return d
    parts = [p.strip() for p in attr.strip(";").split(";") if p.strip()]
    for p in parts:
        if "=" in p:
            k, v = p.split("=", 1)
            d[k.strip()] = v.strip()
        else:
            sp = p.split(None, 1)
            if len(sp) == 2:
                k, v = sp[0].strip(), sp[1].strip().strip('"')
                d[k] = v
    return d


def _split_extra_args(s: str) -> List[str]:
    s = (s or "").strip()
    if not s:
        return []
    return shlex.split(s)


def _which(cmd: str) -> Optional[str]:
    try:
        return shutil.which(cmd)
    except Exception:
        return None


# -------------------------
# Peaks
# -------------------------


CHR_CANDIDATES = ["chr", "chrom", "chromosome", "seqid", "scaffold", "contig", "LG"]
POS_CANDIDATES = ["pos", "position", "bp", "bp_pos", "bp_position", "peak", "peak_pos"]
PVAL_CANDIDATES = ["p", "pvalue", "p_value", "pval", "p_val", "p.value"]
LOD_CANDIDATES = ["lod", "lod_score", "lodscore", "lod_score"]


@dataclass
class Peak:
    index: int
    chr: str
    pos: int
    score: Optional[float] = None  # pvalue or lod
    score_type: str = ""  # 'pvalue' or 'lod'


def read_peaks(path: Path, max_peaks: int = 50) -> List[Peak]:
    """Best-effort peak extraction.

    - If the file looks like a full GWAS results table, we take the top `max_peaks`
      by p-value (ascending).
    - If it looks like a QTL scan/peaks table, we take the top `max_peaks` by LOD
      (descending).
    - Otherwise, we take the first `max_peaks` rows with chr/pos.
    """
    try:
        import pandas as pd  # type: ignore

        df = pd.read_csv(path, sep=_sniff_delimiter(path))
        cols = list(df.columns)
        chr_col = _pick_col(cols, CHR_CANDIDATES)
        pos_col = _pick_col(cols, POS_CANDIDATES)
        p_col = _pick_col(cols, PVAL_CANDIDATES)
        lod_col = _pick_col(cols, LOD_CANDIDATES)
        if not chr_col or not pos_col:
            return []
        sub = df[[chr_col, pos_col] + ([p_col] if p_col else []) + ([lod_col] if lod_col else [])].copy()
        sub.rename(columns={chr_col: "chr", pos_col: "pos"}, inplace=True)
        sub["pos"] = pd.to_numeric(sub["pos"], errors="coerce")
        sub = sub.dropna(subset=["chr", "pos"])
        sub["pos"] = sub["pos"].astype(int)

        score_type = ""
        if p_col and p_col in sub.columns:
            sub["pvalue"] = pd.to_numeric(sub[p_col], errors="coerce")
            score_type = "pvalue"
            sub = sub.dropna(subset=["pvalue"]).sort_values("pvalue", ascending=True)
        elif lod_col and lod_col in sub.columns:
            sub["lod"] = pd.to_numeric(sub[lod_col], errors="coerce")
            score_type = "lod"
            sub = sub.dropna(subset=["lod"]).sort_values("lod", ascending=False)

        sub = sub.head(int(max_peaks))
        peaks: List[Peak] = []
        for i, r in enumerate(sub.itertuples(index=False), start=1):
            chr_ = str(getattr(r, "chr")).strip()
            pos = int(getattr(r, "pos"))
            score = None
            if score_type == "pvalue" and hasattr(r, "pvalue"):
                score = float(getattr(r, "pvalue"))
            if score_type == "lod" and hasattr(r, "lod"):
                score = float(getattr(r, "lod"))
            peaks.append(Peak(index=i, chr=chr_, pos=pos, score=score, score_type=score_type))
        return peaks
    except Exception:
        cols, rows = _read_rows(path)
        chr_col = _pick_col(cols, CHR_CANDIDATES)
        pos_col = _pick_col(cols, POS_CANDIDATES)
        p_col = _pick_col(cols, PVAL_CANDIDATES)
        lod_col = _pick_col(cols, LOD_CANDIDATES)
        if not chr_col or not pos_col:
            return []

        tmp: List[Tuple[float, Peak]] = []
        peaks: List[Peak] = []
        for r in rows:
            chr_ = str(r.get(chr_col, "")).strip()
            pos = _safe_int(r.get(pos_col))
            if not chr_ or pos is None:
                continue
            if p_col and r.get(p_col) is not None:
                p = _safe_float(r.get(p_col))
                if p is None:
                    continue
                tmp.append((p, Peak(index=0, chr=chr_, pos=pos, score=p, score_type="pvalue")))
            elif lod_col and r.get(lod_col) is not None:
                lod = _safe_float(r.get(lod_col))
                if lod is None:
                    continue
                # sort descending => use negative as key
                tmp.append((-lod, Peak(index=0, chr=chr_, pos=pos, score=lod, score_type="lod")))
            else:
                peaks.append(Peak(index=0, chr=chr_, pos=pos, score=None, score_type=""))

        if tmp:
            tmp.sort(key=lambda t: t[0])
            peaks = [pk for _, pk in tmp[: int(max_peaks)]]
        peaks = peaks[: int(max_peaks)]
        for i, pk in enumerate(peaks, start=1):
            pk.index = i
        return peaks


# -------------------------
# FASTA extraction
# -------------------------


def _open_text_maybe_gz(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="ignore")
    return path.open("r", encoding="utf-8", errors="ignore")


def load_fasta_dict(genome_fasta: Path) -> Dict[str, str]:
    """Fallback FASTA loader (loads whole genome; use only if samtools is absent)."""
    seqs: Dict[str, List[str]] = {}
    name: Optional[str] = None
    with _open_text_maybe_gz(genome_fasta) as f:
        for line in f:
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:].strip().split()[0]
                seqs[name] = []
            else:
                if name is None:
                    continue
                seqs[name].append(line.strip())
    return {k: "".join(v).upper() for k, v in seqs.items()}


def ensure_fai(genome_fasta: Path, samtools_bin: str = "samtools") -> None:
    fai = Path(str(genome_fasta) + ".fai")
    if fai.exists():
        return
    subprocess.run([samtools_bin, "faidx", str(genome_fasta)], check=False, capture_output=True, text=True)


def fetch_seq(genome_fasta: Path, chrom: str, start: int, end: int, samtools_bin: str = "samtools", fasta_cache: Optional[Dict[str, str]] = None) -> str:
    """Fetch sequence (1-based inclusive)."""
    start = max(1, int(start))
    end = max(start, int(end))
    # Prefer samtools faidx if available
    if _which(samtools_bin):
        ensure_fai(genome_fasta, samtools_bin=samtools_bin)
        region = f"{chrom}:{start}-{end}"
        proc = subprocess.run([samtools_bin, "faidx", str(genome_fasta), region], capture_output=True, text=True, check=False)
        if proc.returncode == 0 and proc.stdout:
            lines = [ln.strip() for ln in proc.stdout.splitlines() if ln.strip() and not ln.startswith(">")]
            return "".join(lines).upper()
        # fall through

    # Fallback: load into memory
    if fasta_cache is None:
        fasta_cache = load_fasta_dict(genome_fasta)
    seq = fasta_cache.get(chrom)
    if not seq:
        # try stripping chr prefix
        if chrom.lower().startswith("chr"):
            seq = fasta_cache.get(chrom[3:])
    if not seq:
        raise RuntimeError(f"Chromosome/contig not found in genome FASTA: {chrom}")
    # python string slicing is 0-based
    return seq[start - 1 : end]


def wrap_fasta(seq: str, width: int = 60) -> str:
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width))


def make_region_id(pk: Peak, start: int, end: int) -> str:
    # Safe header for downstream tools
    chr_clean = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(pk.chr))
    return f"peak{pk.index:04d}__{chr_clean}__{int(start)}__{int(end)}"


# -------------------------
# GFF parsing for masking / outputs
# -------------------------


@dataclass
class Interval:
    start: int
    end: int


def load_known_gene_intervals(gff_path: Path) -> Dict[str, List[Interval]]:
    """Parse a GFF/GTF and collect gene-like intervals per chromosome."""
    known: Dict[str, List[Interval]] = {}
    with _open_text_maybe_gz(gff_path) as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            seqid, _source, ftype, start, end, _score, _strand, _phase, _attr = parts[:9]
            ftype_l = ftype.lower()
            if ftype_l not in ("gene", "mrna", "transcript", "cds", "exon"):
                continue
            s = _safe_int(start)
            e = _safe_int(end)
            if s is None or e is None:
                continue
            if s > e:
                s, e = e, s
            known.setdefault(str(seqid), []).append(Interval(start=s, end=e))
    # sort
    for k in list(known.keys()):
        known[k].sort(key=lambda iv: (iv.start, iv.end))
    return known


def _overlaps(a_start: int, a_end: int, b_start: int, b_end: int) -> bool:
    return not (a_end < b_start or b_end < a_start)


def interval_overlaps_any(iv: Interval, sorted_intervals: List[Interval]) -> bool:
    """Check overlap against a sorted list (linear scan with early exit).

    For peak-scale use cases this is sufficient. If users provide genome-wide GFFs,
    this still performs fine because we only query a small number of predicted genes.
    """
    for g in sorted_intervals:
        if g.start > iv.end:
            return False
        if _overlaps(iv.start, iv.end, g.start, g.end):
            return True
    return False


# -------------------------
# Predictors
# -------------------------


def run_augustus(
    fasta_in: Path,
    out_gff3: Path,
    augustus_bin: str,
    species: str,
    extra: str,
) -> None:
    cmd = [augustus_bin]
    if species:
        cmd.append(f"--species={species}")
    # Force GFF3 output
    cmd += ["--gff3=on", "--UTR=off"]
    cmd += _split_extra_args(extra)
    cmd.append(str(fasta_in))
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    out_gff3.write_text(proc.stdout or "", encoding="utf-8", errors="ignore")
    (out_gff3.parent / "augustus.stderr.txt").write_text(proc.stderr or "", encoding="utf-8", errors="ignore")
    if proc.returncode != 0:
        raise RuntimeError(f"AUGUSTUS failed (code={proc.returncode}). See augustus.stderr.txt")


def run_snap(
    fasta_in: Path,
    out_gff3: Path,
    snap_bin: str,
    hmm: str,
    extra: str,
) -> None:
    if not hmm:
        raise RuntimeError("SNAP requires snap_hmm")
    cmd = [snap_bin, str(hmm), str(fasta_in), "-gff"]
    cmd += _split_extra_args(extra)
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    out_gff3.write_text(proc.stdout or "", encoding="utf-8", errors="ignore")
    (out_gff3.parent / "snap.stderr.txt").write_text(proc.stderr or "", encoding="utf-8", errors="ignore")
    if proc.returncode != 0:
        raise RuntimeError(f"SNAP failed (code={proc.returncode}). See snap.stderr.txt")


def run_braker3(
    fasta_in: Path,
    out_dir: Path,
    braker_bin: str,
    threads: int,
    species: str,
    prot_faa: str,
    rnaseq_bam: str,
    extra: str,
) -> Path:
    """Run BRAKER3 and return a path to the best-guess GFF/GTF output."""
    out_dir.mkdir(parents=True, exist_ok=True)

    # BRAKER3 generally needs evidence; accept either proteins or BAM.
    if not prot_faa and not rnaseq_bam:
        raise RuntimeError("BRAKER3 requires at least one evidence input: braker_proteins_faa or braker_rnaseq_bam")

    # Prefer full path if provided; otherwise rely on PATH.
    cmd = [braker_bin]
    cmd.append(f"--genome={str(fasta_in)}")
    cmd.append(f"--workingdir={str(out_dir)}")
    if threads:
        cmd.append(f"--cores={int(threads)}")
    if species:
        cmd.append(f"--species={species}")
    if prot_faa:
        cmd.append(f"--prot_seq={prot_faa}")
    if rnaseq_bam:
        cmd.append(f"--bam={rnaseq_bam}")
    cmd += _split_extra_args(extra)

    proc = subprocess.run(cmd, cwd=str(out_dir), capture_output=True, text=True, check=False)
    (out_dir / "braker.stdout.txt").write_text(proc.stdout or "", encoding="utf-8", errors="ignore")
    (out_dir / "braker.stderr.txt").write_text(proc.stderr or "", encoding="utf-8", errors="ignore")
    if proc.returncode != 0:
        raise RuntimeError(f"BRAKER3 failed (code={proc.returncode}). See braker.stderr.txt")

    # Best-effort output discovery
    candidates = [
        out_dir / "braker.gff3",
        out_dir / "braker.gtf",
        out_dir / "braker.gff",
        out_dir / "augustus.hints.gff3",
        out_dir / "augustus.hints.gtf",
        out_dir / "augustus.gff3",
        out_dir / "augustus.gtf",
    ]
    for p in candidates:
        if p.exists() and p.stat().st_size > 0:
            return p

    # search broadly
    for p in out_dir.rglob("*.gff3"):
        if p.stat().st_size > 0:
            return p
    for p in out_dir.rglob("*.gtf"):
        if p.stat().st_size > 0:
            return p
    raise RuntimeError("BRAKER3 finished but no GFF/GTF output was found in workingdir")


def rewrite_gff_source(in_path: Path, out_path: Path, source_name: str) -> None:
    """Rewrite the 2nd column (source) for non-comment lines."""
    out_lines: List[str] = []
    with in_path.open("r", encoding="utf-8", errors="ignore") as f:
        for ln in f:
            if not ln or ln.startswith("#"):
                out_lines.append(ln.rstrip("\n"))
                continue
            parts = ln.rstrip("\n").split("\t")
            if len(parts) >= 2:
                parts[1] = source_name
                out_lines.append("\t".join(parts))
            else:
                out_lines.append(ln.rstrip("\n"))
    out_path.write_text("\n".join(out_lines) + "\n", encoding="utf-8", errors="ignore")


# -------------------------
# Main workflow
# -------------------------


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--params", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    params = json.loads(Path(args.params).read_text(encoding="utf-8"))
    out_dir = Path(args.out).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    genome_fasta = Path(str(params.get("genome_fasta", "")).strip()).expanduser()
    peaks_tsv = Path(str(params.get("peaks_tsv", "")).strip()).expanduser()
    flank_bp = int(float(params.get("flank_bp", 50000) or 50000))
    tool = str(params.get("tool", "braker3")).strip().lower()
    threads = int(float(params.get("threads", 4) or 4))

    # Existing GUI fields (AUGUSTUS/SNAP)
    augustus_bin = str(params.get("augustus_bin", "augustus")).strip() or "augustus"
    augustus_species = str(params.get("augustus_species", "")).strip()
    augustus_extra = str(params.get("augustus_extra", "")).strip()
    snap_bin = str(params.get("snap_bin", "snap")).strip() or "snap"
    snap_hmm = str(params.get("snap_hmm", "")).strip()
    snap_extra = str(params.get("snap_extra", "")).strip()

    # BRAKER3 fields (new)
    braker_bin = str(params.get("braker_bin", "braker.pl")).strip() or "braker.pl"
    braker_species = str(params.get("braker_species", "braker3_peak"))
    braker_prot = str(params.get("braker_proteins_faa", "")).strip()
    braker_bam = str(params.get("braker_rnaseq_bam", "")).strip()
    braker_extra = str(params.get("braker_extra", "")).strip()

    mask_gff = str(params.get("mask_gff", "")).strip()
    mask_mode = str(params.get("mask_mode", "none")).strip().lower()
    samtools_bin = str(params.get("samtools_bin", "samtools")).strip() or "samtools"
    max_peaks = int(float(params.get("max_peaks", 50) or 50))

    if not genome_fasta.exists():
        raise SystemExit(f"genome_fasta not found: {genome_fasta}")
    if not peaks_tsv.exists():
        raise SystemExit(f"peaks_tsv not found: {peaks_tsv}")

    peaks = read_peaks(peaks_tsv, max_peaks=max_peaks)
    if not peaks:
        raise SystemExit("No peaks could be read from peaks_tsv (need chr/pos columns)")

    # Known genes for masking/filtering
    known: Dict[str, List[Interval]] = {}
    if mask_gff and Path(mask_gff).exists() and mask_mode in ("mask_known_genes", "drop_overlapping"):
        known = load_known_gene_intervals(Path(mask_gff))

    regions_tsv = out_dir / "regions.tsv"
    targets_fa = out_dir / "targets.fa"
    targets_masked_fa = out_dir / "targets.masked.fa"

    # If samtools not available, fallback cache
    fasta_cache: Optional[Dict[str, str]] = None
    if not _which(samtools_bin):
        fasta_cache = load_fasta_dict(genome_fasta)

    # Build region FASTA
    region_rows: List[Dict[str, Any]] = []
    fa_lines: List[str] = []
    fa_lines_masked: List[str] = []

    for pk in peaks:
        r_start = max(1, pk.pos - flank_bp)
        r_end = pk.pos + flank_bp
        rid = make_region_id(pk, r_start, r_end)
        seq = fetch_seq(genome_fasta, pk.chr, r_start, r_end, samtools_bin=samtools_bin, fasta_cache=fasta_cache)

        # mask known genes if requested
        seq_masked = seq
        if mask_mode == "mask_known_genes" and known:
            ivs = known.get(pk.chr, [])
            if ivs:
                arr = list(seq_masked)
                for iv in ivs:
                    if not _overlaps(r_start, r_end, iv.start, iv.end):
                        continue
                    ms = max(r_start, iv.start)
                    me = min(r_end, iv.end)
                    # convert to 0-based within region
                    a = ms - r_start
                    b = me - r_start
                    for i in range(a, b + 1):
                        if 0 <= i < len(arr):
                            arr[i] = "N"
                seq_masked = "".join(arr)

        fa_lines.append(f">{rid}")
        fa_lines.append(wrap_fasta(seq))
        fa_lines_masked.append(f">{rid}")
        fa_lines_masked.append(wrap_fasta(seq_masked))

        region_rows.append(
            {
                "peak_index": pk.index,
                "peak_chr": pk.chr,
                "peak_pos": pk.pos,
                "region_id": rid,
                "region_chr": pk.chr,
                "region_start": r_start,
                "region_end": r_end,
            }
        )

    regions_tsv.write_text(
        "\t".join(region_rows[0].keys()) + "\n" + "\n".join("\t".join(str(r[k]) for k in region_rows[0].keys()) for r in region_rows) + "\n",
        encoding="utf-8",
    )
    targets_fa.write_text("\n".join(fa_lines) + "\n", encoding="utf-8")
    targets_masked_fa.write_text("\n".join(fa_lines_masked) + "\n", encoding="utf-8")

    fasta_for_pred = targets_masked_fa if (mask_mode == "mask_known_genes" and known) else targets_fa

    # Run predictor(s)
    predicted_gff3 = out_dir / "predicted.gff3"

    if tool == "braker3":
        braker_work = out_dir / "braker3_work"
        out_path = run_braker3(
            fasta_in=fasta_for_pred,
            out_dir=braker_work,
            braker_bin=braker_bin,
            threads=threads,
            species=braker_species,
            prot_faa=braker_prot,
            rnaseq_bam=braker_bam,
            extra=braker_extra,
        )
        # copy to predicted.gff3
        shutil.copyfile(out_path, predicted_gff3)

    elif tool == "augustus":
        tmp = out_dir / "augustus.raw.gff3"
        run_augustus(fasta_for_pred, tmp, augustus_bin=augustus_bin, species=augustus_species, extra=augustus_extra)
        rewrite_gff_source(tmp, predicted_gff3, "AUGUSTUS")

    elif tool == "snap":
        tmp = out_dir / "snap.raw.gff3"
        run_snap(fasta_for_pred, tmp, snap_bin=snap_bin, hmm=snap_hmm, extra=snap_extra)
        rewrite_gff_source(tmp, predicted_gff3, "SNAP")

    elif tool == "both":
        aug = out_dir / "augustus.raw.gff3"
        snp = out_dir / "snap.raw.gff3"
        run_augustus(fasta_for_pred, aug, augustus_bin=augustus_bin, species=augustus_species, extra=augustus_extra)
        run_snap(fasta_for_pred, snp, snap_bin=snap_bin, hmm=snap_hmm, extra=snap_extra)
        aug2 = out_dir / "augustus.gff3"
        snp2 = out_dir / "snap.gff3"
        rewrite_gff_source(aug, aug2, "AUGUSTUS")
        rewrite_gff_source(snp, snp2, "SNAP")
        predicted_gff3.write_text("", encoding="utf-8")
        with predicted_gff3.open("a", encoding="utf-8") as w:
            w.write(aug2.read_text(encoding="utf-8", errors="ignore"))
            if not aug2.read_text(encoding="utf-8", errors="ignore").endswith("\n"):
                w.write("\n")
            w.write(snp2.read_text(encoding="utf-8", errors="ignore"))

    else:
        raise SystemExit(f"Unsupported tool: {tool}")

    # Build map: region_id -> (peak_chr, peak_pos, region_start)
    region_map: Dict[str, Dict[str, Any]] = {r["region_id"]: r for r in region_rows}

    # Parse predicted.gff3 and emit predicted_near_peaks.tsv
    out_tsv = out_dir / "predicted_near_peaks.tsv"
    rows_out: List[Dict[str, Any]] = []
    pred_i = 0
    with predicted_gff3.open("r", encoding="utf-8", errors="ignore") as f:
        for ln in f:
            if not ln or ln.startswith("#"):
                continue
            parts = ln.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            seqid, source, ftype, start_s, end_s, score_s, strand, _phase, attr = parts[:9]
            ftype_l = ftype.lower()

            # prioritize gene-like features
            if ftype_l not in ("gene", "mrna", "transcript"):
                continue
            if seqid not in region_map:
                continue
            s = _safe_int(start_s)
            e = _safe_int(end_s)
            if s is None or e is None:
                continue
            if s > e:
                s, e = e, s
            score = _safe_float(score_s)
            attrs = _parse_attrs(attr)
            pred_id = (
                attrs.get("ID")
                or attrs.get("gene_id")
                or attrs.get("transcript_id")
                or attrs.get("Name")
                or ""
            )
            if not pred_id:
                pred_i += 1
                pred_id = f"pred_{pred_i:05d}"

            info = region_map[seqid]
            peak_index = int(info["peak_index"])
            peak_chr = str(info["peak_chr"])
            peak_pos = int(info["peak_pos"])
            offset = int(info["region_start"])

            # map back to genome coordinates
            abs_start = offset + s - 1
            abs_end = offset + e - 1
            # distance to peak point
            if abs_start <= peak_pos <= abs_end:
                dist = 0
            else:
                dist = min(abs(peak_pos - abs_start), abs(peak_pos - abs_end))

            # optional filter: drop overlapping known genes
            if mask_mode == "drop_overlapping" and known:
                ivs = known.get(peak_chr, [])
                if ivs and interval_overlaps_any(Interval(abs_start, abs_end), ivs):
                    continue

            rows_out.append(
                {
                    "peak_index": peak_index,
                    "peak_chr": peak_chr,
                    "peak_pos": peak_pos,
                    "region_id": str(seqid),
                    "tool": str(tool),
                    "pred_id": str(pred_id),
                    "pred_start": abs_start,
                    "pred_end": abs_end,
                    "strand": strand,
                    "score": "" if score is None else score,
                    "dist_to_peak": dist,
                    "source": source,
                    "feature": ftype,
                }
            )

    # Write TSV
    header = [
        "peak_index",
        "peak_chr",
        "peak_pos",
        "region_id",
        "tool",
        "pred_id",
        "pred_start",
        "pred_end",
        "strand",
        "score",
        "dist_to_peak",
        "source",
        "feature",
    ]
    with out_tsv.open("w", encoding="utf-8", newline="") as w:
        w.write("\t".join(header) + "\n")
        for r in rows_out:
            w.write("\t".join(str(r.get(k, "")) for k in header) + "\n")

    # Human-readable summary
    summary = out_dir / "summary.txt"
    summary.write_text(
        "\n".join(
            [
                f"tool={tool}",
                f"genome_fasta={genome_fasta}",
                f"peaks_tsv={peaks_tsv}",
                f"flank_bp={flank_bp}",
                f"peaks_used={len(peaks)} (max_peaks={max_peaks})",
                f"known_gff={mask_gff if mask_gff else ''}",
                f"mask_mode={mask_mode}",
                f"targets_fa={targets_fa.name}",
                f"predicted_gff3={predicted_gff3.name}",
                f"predicted_near_peaks_tsv={out_tsv.name} rows={len(rows_out)}",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
