#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""integrate_prioritize

One-shot pipeline to prioritize candidate genes around GWAS/QTL peaks.

Inputs (via --params JSON)
-------------------------
- gwas_results: [path ...] (optional)
- qtl_peaks:    [path ...] (optional)
- mode_gwas: 'top_n'|'p_threshold'
- gwas_top_n, gwas_p_threshold
- mode_qtl: 'top_n'|'lod_threshold'
- qtl_top_n, qtl_lod_threshold
- window_bp: int (± window around each peak)
- chr_normalize: bool (strip 'chr' prefix)

Gene sources
------------
- genome_gff: path to GFF3/GTF (recommended)
- predicted_gff: path to predicted genes GFF3 (optional)

Evidence (optional)
-------------------
- deg_tsv_list: [path ...]
- deg_padj_max: float (default 0.05)
- deg_abs_logfc_min: float (default 0.0)
- blast_tsv_list: [path ...] (best_hit_per_query.tsv or blast_results.tsv)
- blast_evalue_max: float (default 1e-5; optional filter)
- blast_bitscore_min: float (default 0; optional filter)
- eggnog_annotations_list: [path ...] (*.annotations from eggNOG-mapper)
- go_map_file: path (optional; GAF or TSV gene_id→GO)
- go_keywords: string (newline/comma separated; GO:xxxx or free keywords)

Scoring
-------
- max_genes_per_locus: int (0=no limit)
- dist_decay_bp: int (default 50000)
- w_peak, w_dist, w_deg: floats

Outputs (under --out)
---------------------
- loci.tsv
- locus_peaks.tsv
- peaks_merged.tsv
- candidates_ranked.tsv
- meta.json

Design:
- minimal deps (no pandas required)
- robust column inference for chr/pos/pvalue/lod/gene_id
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple


# -------------------------
# Utilities
# -------------------------

def _sniff_delim(path: Path) -> str:
    if path.suffix.lower() == ".csv":
        return ","
    try:
        with path.open("r", encoding="utf-8", errors="ignore") as f:
            head = f.readline()
        if head.count("\t") >= head.count(","):
            return "\t"
        return ","
    except Exception:
        return "\t"


def _lower_map(cols: Sequence[str]) -> Dict[str, str]:
    d: Dict[str, str] = {}
    for c in cols:
        lc = str(c).strip().lower()
        if lc and lc not in d:
            d[lc] = str(c)
    return d


def _pick_col(cols: Sequence[str], candidates: Sequence[str]) -> Optional[str]:
    lm = _lower_map(cols)
    for cand in candidates:
        key = cand.lower()
        if key in lm:
            return lm[key]
    return None


def _safe_float(x: Any) -> Optional[float]:
    if x is None:
        return None
    if isinstance(x, (int, float)):
        try:
            v = float(x)
        except Exception:
            return None
        if math.isnan(v) or math.isinf(v):
            return None
        return v
    s = str(x).strip()
    if not s or s.lower() in {"na", "nan", "none", "."}:
        return None
    try:
        v = float(s.replace(",", ""))
        if math.isnan(v) or math.isinf(v):
            return None
        return v
    except Exception:
        return None


def _safe_int(x: Any) -> Optional[int]:
    v = _safe_float(x)
    if v is None:
        return None
    try:
        return int(round(v))
    except Exception:
        return None


def _strip_chr_prefix(x: str) -> str:
    s = str(x).strip()
    if s.lower().startswith("chr"):
        return s[3:].strip()
    return s


def _dataset_label(path: Path) -> str:
    name = path.name
    for suf in [".tsv", ".csv", ".txt"]:
        if name.lower().endswith(suf):
            name = name[: -len(suf)]
            break
    return name


def read_rows(path: Path) -> Tuple[List[str], List[Dict[str, Any]]]:
    path = Path(path)
    delim = _sniff_delim(path)
    with path.open("r", encoding="utf-8", errors="ignore", newline="") as f:
        reader = csv.DictReader(f, delimiter=delim)
        cols = reader.fieldnames or []
        rows = [r for r in reader]
    return [str(c) for c in cols], rows


def write_tsv(path: Path, header: List[str], rows: List[List[Any]]):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(header)
        for r in rows:
            w.writerow(["" if x is None else x for x in r])


# -------------------------
# Peak picking + loci
# -------------------------

CHR_CANDIDATES = [
    "chr", "chrom", "chromosome", "seqid", "scaffold", "contig", "chromosome_id", "peak_chr", "CHROM",
]
POS_CANDIDATES = [
    "pos", "position", "bp", "bp_pos", "snp_pos", "physical_pos", "peak_pos", "start", "POS",
]
PVAL_CANDIDATES = [
    "p", "pvalue", "p_value", "p-val", "p.val", "pval", "P", "PVAL", "p_wald", "p_wald_dom", "p_wald_add", "p_wald_gen",
]
LOD_CANDIDATES = [
    "lod", "LOD", "lod_score", "lodscore", "lod.max", "lod_max", "maxlod", "max_lod",
]

# QTL peak tables often include a marker column that can be mapped to physical bp by marker_map.tsv
MARKER_CANDIDATES = [
    "marker", "snp", "snp_id", "marker_id", "id", "rs", "locus", "peak_marker", "closest_marker", "best_marker",
    "markername", "marker_name",
]


def _chr_key(x: str) -> tuple:
    """Natural-ish chromosome sort key."""
    s = str(x).strip()
    s2 = _strip_chr_prefix(s)
    if s2.isdigit():
        return (0, int(s2), s)
    up = s2.upper()
    if up in {"X", "Y", "W", "Z", "MT", "M"}:
        order = {"X": 1001, "Y": 1002, "W": 1003, "Z": 1004, "MT": 1005, "M": 1005}
        return (1, order.get(up, 2000), s)
    return (2, s2, s)


def _auto_find_marker_map(qtl_path: Path) -> Optional[Path]:
    """Try to locate marker_map.tsv near a QTL peaks file (same dir or parents)."""
    try:
        p = Path(qtl_path)
        for base in [p.parent, p.parent.parent, p.parent.parent.parent]:
            cand = base / "marker_map.tsv"
            if cand.exists():
                return cand
    except Exception:
        return None
    return None


def load_marker_map(marker_map_tsv: Path, normalize_chr: bool) -> Dict[str, tuple]:
    """Load marker_map.tsv (marker/chr/pos) -> {marker: (chr, pos_bp)}."""
    cols, rows = read_rows(marker_map_tsv)
    m_col = _pick_col(cols, MARKER_CANDIDATES) or _pick_col(cols, ["marker"])
    chr_col = _pick_col(cols, CHR_CANDIDATES)
    pos_col = _pick_col(cols, POS_CANDIDATES)
    if not m_col or not chr_col or not pos_col:
        return {}
    mp: Dict[str, tuple] = {}
    for r in rows:
        mid = str(r.get(m_col, '')).strip()
        if not mid:
            continue
        chr_ = str(r.get(chr_col, '')).strip()
        if normalize_chr:
            chr_ = _strip_chr_prefix(chr_)
        pos = _safe_int(r.get(pos_col))
        if not chr_ or pos is None:
            continue
        mp[mid] = (chr_, int(pos))
    return mp


def _pos_to_bp(pos_val: float, *, assume_cm_when_small: bool = True) -> int:
    """Convert a numeric position to bp.

    If the value looks like cM/Mb (small), apply 1e6 scaling; otherwise treat as bp.
    """
    try:
        v = float(pos_val)
    except Exception:
        return 0
    if assume_cm_when_small and abs(v) < 5e5:
        return int(round(v * 1_000_000.0))
    return int(round(v))


@dataclass
class Peak:
    source_type: str   # gwas/qtl
    dataset: str
    raw_file: str
    chr: str
    pos: int
    start: int
    end: int
    pvalue: Optional[float] = None
    lod: Optional[float] = None


@dataclass
class Locus:
    locus_id: str
    chr: str
    start: int
    end: int
    rep_pos: int
    rep_source: str
    rep_pvalue: Optional[float]
    rep_lod: Optional[float]


def _iter_gwas_peaks(path: Path, mode: str, top_n: int, p_threshold: float, window_bp: int, normalize_chr: bool) -> List[Peak]:
    ds = _dataset_label(path)
    cols, rows = read_rows(path)
    chr_col = _pick_col(cols, CHR_CANDIDATES)
    pos_col = _pick_col(cols, POS_CANDIDATES)
    p_col = _pick_col(cols, PVAL_CANDIDATES)
    if not chr_col or not pos_col or not p_col:
        return []

    recs: List[Tuple[str, int, float]] = []
    for r in rows:
        c = r.get(chr_col)
        p = _safe_int(r.get(pos_col))
        pv = _safe_float(r.get(p_col))
        if c is None or p is None or pv is None:
            continue
        c2 = _strip_chr_prefix(c) if normalize_chr else str(c).strip()
        if not c2:
            continue
        recs.append((c2, int(p), float(pv)))

    if not recs:
        return []

    picked: List[Tuple[str, int, float]]
    if mode == "p_threshold":
        picked = [t for t in recs if t[2] <= p_threshold]
        picked.sort(key=lambda x: x[2])
        picked = picked[:top_n] if top_n > 0 else picked
    else:
        picked = sorted(recs, key=lambda x: x[2])[:top_n]

    out: List[Peak] = []
    for c2, pos, pv in picked:
        out.append(Peak(
            source_type="gwas",
            dataset=ds,
            raw_file=str(path),
            chr=c2,
            pos=int(pos),
            start=max(0, int(pos) - int(window_bp)),
            end=int(pos) + int(window_bp),
            pvalue=float(pv),
            lod=None,
        ))
    return out


def _iter_qtl_peaks(
    path: Path,
    mode: str,
    top_n: int,
    lod_threshold: float,
    window_bp: int,
    normalize_chr: bool,
    marker_map: Optional[Dict[str, tuple]] = None,
) -> List[Peak]:
    ds = _dataset_label(path)
    cols, rows = read_rows(path)
    chr_col = _pick_col(cols, CHR_CANDIDATES)
    pos_col = _pick_col(cols, POS_CANDIDATES)
    lod_col = _pick_col(cols, LOD_CANDIDATES)
    if not chr_col or not pos_col or not lod_col:
        return []

    m_col = _pick_col(cols, MARKER_CANDIDATES)
    recs: List[Tuple[str, float, float, str]] = []
    for r in rows:
        c = r.get(chr_col)
        p = _safe_float(r.get(pos_col))
        lod = _safe_float(r.get(lod_col))
        if c is None or p is None or lod is None:
            continue
        c2 = _strip_chr_prefix(c) if normalize_chr else str(c).strip()
        if not c2:
            continue
        mk = str(r.get(m_col, '')).strip() if m_col else ''
        recs.append((c2, float(p), float(lod), mk))

    if not recs:
        return []

    picked: List[Tuple[str, float, float, str]]
    if mode == "lod_threshold":
        picked = [t for t in recs if t[2] >= lod_threshold]
        picked.sort(key=lambda x: x[2], reverse=True)
        picked = picked[:top_n] if top_n > 0 else picked
    else:
        picked = sorted(recs, key=lambda x: x[2], reverse=True)[:top_n]

    out: List[Peak] = []
    for c2, pos_raw, lod, mk in picked:
        # Convert QTL coordinate to bp.
        chr_use = c2
        pos_bp: int
        if marker_map is not None and mk and mk in marker_map:
            try:
                chr_use, pos_bp = marker_map[mk]
            except Exception:
                pos_bp = _pos_to_bp(pos_raw)
        else:
            pos_bp = _pos_to_bp(pos_raw)
        out.append(Peak(
            source_type="qtl",
            dataset=ds,
            raw_file=str(path),
            chr=str(chr_use),
            pos=int(pos_bp),
            start=max(0, int(pos_bp) - int(window_bp)),
            end=int(pos_bp) + int(window_bp),
            pvalue=None,
            lod=float(lod),
        ))
    return out


def build_loci(peaks: List[Peak]) -> Tuple[List[Locus], List[Peak]]:
    """Cluster peak windows into loci per chromosome."""
    if not peaks:
        return [], []

    peaks_sorted = sorted(peaks, key=lambda x: (x.chr, x.start, x.end, x.pos))
    loci: List[Locus] = []
    locus_peaks: List[Peak] = []

    cur_chr = None
    cur_start = None
    cur_end = None
    cur_peaks: List[Peak] = []

    def _flush():
        nonlocal loci, locus_peaks, cur_chr, cur_start, cur_end, cur_peaks
        if not cur_peaks or cur_chr is None or cur_start is None or cur_end is None:
            cur_peaks = []
            return
        locus_id = f"Locus{len(loci)+1:04d}"

        # Representative peak: prefer smallest pvalue (GWAS), else highest LOD (QTL)
        rep = None
        # Any GWAS peaks?
        gwas = [p for p in cur_peaks if p.pvalue is not None]
        if gwas:
            rep = sorted(gwas, key=lambda p: p.pvalue if p.pvalue is not None else 1.0)[0]
        else:
            rep = sorted(cur_peaks, key=lambda p: p.lod if p.lod is not None else -1e9, reverse=True)[0]

        loci.append(Locus(
            locus_id=locus_id,
            chr=cur_chr,
            start=int(cur_start),
            end=int(cur_end),
            rep_pos=int(rep.pos),
            rep_source=str(rep.source_type),
            rep_pvalue=rep.pvalue,
            rep_lod=rep.lod,
        ))

        for p in cur_peaks:
            locus_peaks.append(p)

        cur_peaks = []

    for p in peaks_sorted:
        if cur_chr is None:
            cur_chr = p.chr
            cur_start = p.start
            cur_end = p.end
            cur_peaks = [p]
            continue
        if p.chr != cur_chr or p.start > (cur_end if cur_end is not None else p.start):
            _flush()
            cur_chr = p.chr
            cur_start = p.start
            cur_end = p.end
            cur_peaks = [p]
        else:
            cur_end = max(int(cur_end), int(p.end))
            cur_start = min(int(cur_start), int(p.start))
            cur_peaks.append(p)

    _flush()

    return loci, locus_peaks


# -------------------------
# GFF/GTF parser
# -------------------------

@dataclass
class Gene:
    gene_id: str
    chr: str
    start: int
    end: int
    strand: str
    source: str  # gff/pred


def _parse_gff_attrs(attr: str) -> Dict[str, str]:
    d: Dict[str, str] = {}
    if not attr:
        return d
    # GTF style: key "value";
    if '"' in attr and "=" not in attr:
        # naive parse
        parts = [p.strip() for p in attr.strip().split(";") if p.strip()]
        for p in parts:
            if " " in p:
                k, v = p.split(" ", 1)
                k = k.strip()
                v = v.strip().strip('"')
                if k and v:
                    d[k] = v
        return d

    # GFF3 style key=value;key2=value2
    parts = [p.strip() for p in attr.strip().split(";") if p.strip()]
    for p in parts:
        if "=" in p:
            k, v = p.split("=", 1)
            k = k.strip(); v = v.strip()
            if k and v:
                d[k] = v
    return d


def read_genes_from_gff(path: Path, source: str, normalize_chr: bool) -> List[Gene]:
    genes: List[Gene] = []
    if not path or not str(path).strip():
        return genes
    path = Path(path)
    if not path.exists():
        return genes

    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            line = line.rstrip("\n")
            parts = line.split("\t")
            if len(parts) < 9:
                continue
            seqid, _src, ftype, start, end, _score, strand, _phase, attrs = parts[:9]
            if ftype.lower() != "gene":
                # allow common gene-like features
                if ftype.lower() not in {"pseudogene", "ncrna_gene", "mirna_gene"}:
                    continue
            s = _safe_int(start)
            e = _safe_int(end)
            if s is None or e is None:
                continue
            chr_ = _strip_chr_prefix(seqid) if normalize_chr else str(seqid).strip()
            if not chr_:
                continue
            ad = _parse_gff_attrs(attrs)
            gid = ad.get("Name") or ad.get("ID") or ad.get("gene_id") or ad.get("gene") or ad.get("locus_tag")
            if not gid:
                # try Parent for some GFFs
                gid = ad.get("Parent") or ""
            gid = str(gid).strip()
            if not gid:
                continue
            genes.append(Gene(
                gene_id=gid,
                chr=chr_,
                start=min(int(s), int(e)),
                end=max(int(s), int(e)),
                strand=str(strand).strip() if strand else ".",
                source=source,
            ))

    return genes


# -------------------------
# DEG evidence
# -------------------------

DEG_GENE_CANDIDATES = ["gene_id", "gene", "Gene", "GeneID", "id", "ID", "symbol", "locus", "locus_tag"]
DEG_LOGFC_CANDIDATES = ["logfc", "logFC", "log2fc", "log2FoldChange", "log2FC", "lfc"]
DEG_PADJ_CANDIDATES = ["padj", "FDR", "fdr", "adj.P.Val", "adj_p_val", "qvalue", "q_value", "qval"]


def load_deg_tables(paths: List[Path]) -> Dict[str, Dict[str, float]]:
    """Return gene_id -> {padj, logFC} aggregated across multiple DEG tables (best padj)."""
    out: Dict[str, Dict[str, float]] = {}
    for p in paths:
        if not p.exists():
            continue
        cols, rows = read_rows(p)
        gene_col = _pick_col(cols, DEG_GENE_CANDIDATES)
        padj_col = _pick_col(cols, [c.lower() for c in DEG_PADJ_CANDIDATES] + DEG_PADJ_CANDIDATES)
        lfc_col = _pick_col(cols, [c.lower() for c in DEG_LOGFC_CANDIDATES] + DEG_LOGFC_CANDIDATES)
        if not gene_col:
            continue
        for r in rows:
            gid = str(r.get(gene_col) or "").strip()
            if not gid:
                continue
            padj = _safe_float(r.get(padj_col)) if padj_col else None
            lfc = _safe_float(r.get(lfc_col)) if lfc_col else None
            if padj is None and lfc is None:
                continue
            prev = out.get(gid)
            if prev is None:
                out[gid] = {}
                prev = out[gid]
            # keep best (smallest) padj, and corresponding logFC if present
            if padj is not None:
                if ("padj" not in prev) or (padj < float(prev.get("padj", 1e9))):
                    prev["padj"] = float(padj)
                    if lfc is not None:
                        prev["logFC"] = float(lfc)
                else:
                    # keep strongest abs logFC if padj ties are absent
                    if lfc is not None and ("logFC" not in prev or abs(lfc) > abs(float(prev.get("logFC", 0.0)))):
                        prev["logFC"] = float(lfc)
            else:
                if lfc is not None and ("logFC" not in prev or abs(lfc) > abs(float(prev.get("logFC", 0.0)))):
                    prev["logFC"] = float(lfc)
    return out




# -------------------------
# BLAST evidence (best hit per query)
# -------------------------

BLAST_QUERY_CANDIDATES = ["gene_id", "qseqid", "query", "query_id", "qacc", "qaccver", "qid"]
BLAST_SUBJECT_CANDIDATES = ["sseqid", "subject", "subject_id", "sacc", "saccver", "hit"]
BLAST_EVALUE_CANDIDATES = ["evalue", "e-value", "e_value", "Evalue", "E-value"]
BLAST_BITSCORE_CANDIDATES = ["bitscore", "bit_score", "score", "Score"]


def _split_multi(s: str) -> List[str]:
    if s is None:
        return []
    txt = str(s).strip()
    if not txt:
        return []
    parts = re.split(r"[\s,;|]+", txt)
    return [p for p in (x.strip() for x in parts) if p]


def load_blast_tables(paths: List[Path]) -> Dict[str, Dict[str, Any]]:
    """Return gene_id(query) -> {blast_evalue, blast_bitscore, blast_best_hit}."""
    out: Dict[str, Dict[str, Any]] = {}
    for fp in paths:
        if not fp.exists():
            continue
        cols, rows = read_rows(fp)
        q_col = _pick_col(cols, BLAST_QUERY_CANDIDATES)
        s_col = _pick_col(cols, BLAST_SUBJECT_CANDIDATES)
        e_col = _pick_col(cols, BLAST_EVALUE_CANDIDATES)
        b_col = _pick_col(cols, BLAST_BITSCORE_CANDIDATES)
        if not q_col:
            continue

        for r in rows:
            q = str(r.get(q_col) or "").strip()
            if not q:
                continue
            e = _safe_float(r.get(e_col)) if e_col else None
            bs = _safe_float(r.get(b_col)) if b_col else None
            hit = str(r.get(s_col) or "").strip() if s_col else ""
            if e is None and bs is None and not hit:
                continue

            cur = out.get(q)
            # prefer smallest evalue; if missing, prefer largest bitscore
            if cur is None:
                out[q] = {"blast_evalue": e, "blast_bitscore": bs, "blast_best_hit": hit, "blast_source": str(fp)}
            else:
                ce = cur.get("blast_evalue")
                cbs = cur.get("blast_bitscore")
                better = False
                if e is not None and (ce is None or float(e) < float(ce)):
                    better = True
                elif e is None and ce is None and bs is not None and (cbs is None or float(bs) > float(cbs)):
                    better = True
                if better:
                    out[q] = {"blast_evalue": e, "blast_bitscore": bs, "blast_best_hit": hit, "blast_source": str(fp)}
                else:
                    # keep best hit string if empty
                    if (not cur.get("blast_best_hit")) and hit:
                        cur["blast_best_hit"] = hit
                    # keep max bitscore
                    if bs is not None and (cbs is None or float(bs) > float(cbs)):
                        cur["blast_bitscore"] = float(bs)
    return out


# -------------------------
# eggNOG-mapper evidence (*.annotations)
# -------------------------

def read_eggnog_annotations(path: Path) -> Tuple[List[str], List[Dict[str, Any]]]:
    """Read eggNOG-mapper .annotations output with comment lines.

    The header line is typically '#query\tseed_ortholog\t...'.
    """
    path = Path(path)
    if not path.exists():
        return [], []

    header: List[str] = []
    rows: List[Dict[str, Any]] = []
    with path.open('r', encoding='utf-8', errors='ignore', newline='') as f:
        for line in f:
            if not line:
                continue
            if line.startswith('#query') or line.startswith('query'):
                header = [h.lstrip('#').strip() for h in line.rstrip('\n').split('\t')]
                break
            if line.startswith('#'):
                continue
            # Some files may omit header; fall back to default header from first data row length
            parts = line.rstrip('\n').split('\t')
            if parts and not header:
                header = [f'col{i+1}' for i in range(len(parts))]
                rows.append(dict(zip(header, parts)))
                break

        if header:
            for line in f:
                if not line or line.startswith('#'):
                    continue
                parts = line.rstrip('\n').split('\t')
                if not parts:
                    continue
                if len(parts) < len(header):
                    parts = parts + [''] * (len(header) - len(parts))
                rows.append({header[i]: parts[i] for i in range(len(header))})

    return header, rows


def load_eggnog_tables(paths: List[Path]) -> Dict[str, Dict[str, Any]]:
    """Return gene_id(query) -> eggNOG fields (desc, GOs, Preferred_name, OGs...)."""
    out: Dict[str, Dict[str, Any]] = {}
    for fp in paths:
        if not fp.exists():
            continue
        cols, rows = read_eggnog_annotations(fp)
        if not cols or not rows:
            continue
        # normalize: allow both '#query' and 'query'
        q_col = None
        for c in cols:
            if c.strip().lower() in {'query', '#query'}:
                q_col = c
                break
        if q_col is None:
            q_col = cols[0]

        # pick common informative cols
        def pick(name_candidates):
            return _pick_col(cols, name_candidates)

        desc_col = pick(['Description', 'description', 'desc', 'product'])
        go_col = pick(['GOs', 'GO', 'go', 'go_terms', 'go_term'])
        og_col = pick(['eggNOG_OGs', 'eggNOG_OG', 'ogs', 'OGs'])
        pref_col = pick(['Preferred_name', 'preferred_name', 'pref_name', 'name'])
        kegg_col = pick(['KEGG_ko', 'kegg_ko', 'ko'])
        pfam_col = pick(['PFAMs', 'pfams', 'PFAM'])

        for r in rows:
            q = str(r.get(q_col) or '').strip()
            if not q:
                continue
            cur = out.get(q)
            if cur is None:
                out[q] = {
                    'eggnog_desc': str(r.get(desc_col) or '').strip() if desc_col else '',
                    'eggnog_go': str(r.get(go_col) or '').strip() if go_col else '',
                    'eggnog_ogs': str(r.get(og_col) or '').strip() if og_col else '',
                    'eggnog_pref_name': str(r.get(pref_col) or '').strip() if pref_col else '',
                    'eggnog_kegg_ko': str(r.get(kegg_col) or '').strip() if kegg_col else '',
                    'eggnog_pfam': str(r.get(pfam_col) or '').strip() if pfam_col else '',
                    'eggnog_source': str(fp),
                }
            else:
                # prefer filling empty fields
                for k, col in [('eggnog_desc', desc_col), ('eggnog_go', go_col), ('eggnog_ogs', og_col), ('eggnog_pref_name', pref_col), ('eggnog_kegg_ko', kegg_col), ('eggnog_pfam', pfam_col)]:
                    if (not cur.get(k)) and col:
                        cur[k] = str(r.get(col) or '').strip()
    return out


# -------------------------
# GO mapping (optional)
# -------------------------

GO_GENE_CANDIDATES = ['gene_id', 'gene', 'id', 'ID', 'locus_tag', 'symbol']
GO_TERM_CANDIDATES = ['go', 'go_id', 'GO', 'GO_ID', 'go_term', 'go_terms', 'GOs']


def load_go_map(path: Optional[Path]) -> Dict[str, List[str]]:
    """Load gene -> GO IDs from a GAF or TSV/CSV-like file.

    - GAF: uses DB Object ID (col2) and Symbol (col3) as keys; GO ID is col5.
    - TSV/CSV: tries to infer gene_id + go columns; GO can be multi-valued.
    """
    if path is None:
        return {}
    path = Path(path)
    if not path.exists():
        return {}

    out: Dict[str, set] = {}
    if path.suffix.lower() in {'.gaf', '.gaf.gz'} or path.name.lower().endswith('.gaf'):
        # minimal GAF (ignore gzip here)
        if str(path).endswith('.gz'):
            import gzip
            opener = lambda: gzip.open(path, 'rt', encoding='utf-8', errors='ignore')
        else:
            opener = lambda: path.open('r', encoding='utf-8', errors='ignore')
        with opener() as f:
            for line in f:
                if not line or line.startswith('!'):
                    continue
                parts = line.rstrip('\n').split('\t')
                if len(parts) < 5:
                    continue
                obj_id = parts[1].strip()
                sym = parts[2].strip()
                go = parts[4].strip()
                if not go or not go.upper().startswith('GO:'):
                    continue
                for key in [obj_id, sym]:
                    if key:
                        out.setdefault(key, set()).add(go)
        return {k: sorted(v) for k, v in out.items()}

    cols, rows = read_rows(path)
    gcol = _pick_col(cols, GO_GENE_CANDIDATES)
    tcol = _pick_col(cols, GO_TERM_CANDIDATES)
    if not gcol or not tcol:
        return {}
    for r in rows:
        gid = str(r.get(gcol) or '').strip()
        if not gid:
            continue
        val = str(r.get(tcol) or '').strip()
        if not val:
            continue
        for go in _split_multi(val):
            if go.upper().startswith('GO:'):
                out.setdefault(gid, set()).add(go.upper())
    return {k: sorted(v) for k, v in out.items()}


def _extract_go_ids(s: str) -> set:
    ids = set()
    if not s:
        return ids
    for token in _split_multi(s):
        if token.upper().startswith('GO:') and len(token) >= 7:
            ids.add(token.upper())
    return ids


def _parse_keywords(s: str) -> List[str]:
    if not s:
        return []
    parts = re.split(r"[\n\r,;]+", str(s))
    out = []
    for p in parts:
        t = p.strip()
        if not t:
            continue
        out.append(t)
    return out

# -------------------------
# Candidate selection + scoring
# -------------------------


def _distance_point_to_interval(pos: int, start: int, end: int) -> int:
    if start <= pos <= end:
        return 0
    if pos < start:
        return start - pos
    return pos - end


def _peak_strength(p: Peak) -> Optional[float]:
    if p.pvalue is not None and p.pvalue > 0:
        try:
            return -math.log10(float(p.pvalue))
        except Exception:
            return None
    if p.lod is not None:
        return float(p.lod)
    return None


def _norm01(x: float, vmax: float) -> float:
    if vmax <= 0:
        return 0.0
    return max(0.0, min(1.0, x / vmax))


def prioritize(
    loci: List[Locus],
    locus_peaks: List[Peak],
    genes: List[Gene],
    deg_map: Dict[str, Dict[str, float]],
    blast_map: Dict[str, Dict[str, Any]],
    blast_evalue_max: float,
    blast_bitscore_min: float,
    eggnog_map: Dict[str, Dict[str, Any]],
    go_map: Dict[str, List[str]],
    go_keywords: List[str],
    max_genes_per_locus: int,
    dist_decay_bp: int,
    deg_padj_max: float,
    deg_abs_logfc_min: float,
    w_peak: float,
    w_dist: float,
    w_deg: float,
    w_blast: float,
    w_eggnog: float,
    w_go: float,
) -> Tuple[List[Dict[str, Any]], Dict[str, int]]:

    genes_by_chr: Dict[str, List[Gene]] = {}
    for g in genes:
        genes_by_chr.setdefault(g.chr, []).append(g)
    for c in genes_by_chr:
        genes_by_chr[c].sort(key=lambda g: (g.start, g.end))

    # strength normalization
    strengths = []
    for p in locus_peaks:
        s = _peak_strength(p)
        if s is not None:
            strengths.append(s)
    vmax_peak = max(strengths) if strengths else 0.0

    # deg normalization
    deg_strengths = []
    for v in deg_map.values():
        padj = v.get("padj")
        if padj is not None and padj > 0:
            deg_strengths.append(-math.log10(float(padj)))
    vmax_deg = max(deg_strengths) if deg_strengths else 0.0

    # BLAST normalization
    blast_strengths = []
    if blast_map:
        for v in blast_map.values():
            ev = v.get('blast_evalue')
            if ev is None:
                continue
            try:
                evf = float(ev)
            except Exception:
                continue
            if evf <= 0:
                evf = 1e-300
            blast_strengths.append(-math.log10(evf))
    vmax_blast = max(blast_strengths) if blast_strengths else 0.0


    # peaks by locus_id
    peaks_by_locus: Dict[str, List[Peak]] = {}
    # locus_peaks list aligns with loci build order? Not tagged. We'll reconstruct by interval membership.
    # We assume that locus_peaks were built from build_loci and the locus order is preserved; but peaks don't store locus_id.
    # So we reassign by checking overlap to locus interval.
    for p in locus_peaks:
        # find matching locus (same chr and within start-end)
        for L in loci:
            if L.chr != p.chr:
                continue
            if L.start <= p.pos <= L.end:
                peaks_by_locus.setdefault(L.locus_id, []).append(p)
                break

    cand_rows: List[Dict[str, Any]] = []
    counts = {"n_loci": len(loci), "n_genes_total": len(genes), "n_candidates": 0}

    for L in loci:
        chr_genes = genes_by_chr.get(L.chr, [])
        # collect genes overlapping locus interval
        cand: List[Tuple[int, Gene]] = []
        for g in chr_genes:
            if g.end < L.start:
                continue
            if g.start > L.end:
                break
            dist = _distance_point_to_interval(L.rep_pos, g.start, g.end)
            cand.append((dist, g))

        cand.sort(key=lambda t: (t[0], t[1].start))
        if max_genes_per_locus and max_genes_per_locus > 0:
            cand = cand[: int(max_genes_per_locus)]

        # per locus rep peak strength
        rep_strength = None
        if L.rep_pvalue is not None and L.rep_pvalue > 0:
            rep_strength = -math.log10(float(L.rep_pvalue))
        elif L.rep_lod is not None:
            rep_strength = float(L.rep_lod)
        rep_peak_norm = _norm01(rep_strength or 0.0, vmax_peak)

        for dist, g in cand:
            dist_score = math.exp(-float(dist) / float(max(1, int(dist_decay_bp))))
            score = w_peak * rep_peak_norm + w_dist * dist_score

            # DEG evidence
            deg_padj = None
            deg_lfc = None
            deg_hit = 0
            deg_score = 0.0
            if deg_map and g.gene_id in deg_map:
                dm = deg_map[g.gene_id]
                deg_padj = dm.get("padj")
                deg_lfc = dm.get("logFC")
                ok = True
                if deg_padj is not None and deg_padj_max is not None:
                    ok = ok and (float(deg_padj) <= float(deg_padj_max))
                if deg_lfc is not None and deg_abs_logfc_min is not None:
                    ok = ok and (abs(float(deg_lfc)) >= float(deg_abs_logfc_min))
                if ok:
                    deg_hit = 1
                    if deg_padj is not None and deg_padj > 0:
                        deg_score = _norm01(-math.log10(float(deg_padj)), vmax_deg)
                    else:
                        deg_score = 0.0

            # BLAST evidence
            blast_hit = 0
            blast_evalue = None
            blast_bitscore = None
            blast_best_hit = ''
            blast_score = 0.0
            if blast_map and g.gene_id in blast_map:
                bm = blast_map[g.gene_id]
                blast_evalue = bm.get('blast_evalue')
                blast_bitscore = bm.get('blast_bitscore')
                blast_best_hit = str(bm.get('blast_best_hit') or '')
                ok = True
                if blast_evalue is not None and blast_evalue_max is not None and float(blast_evalue_max) > 0:
                    try:
                        ok = ok and (float(blast_evalue) <= float(blast_evalue_max))
                    except Exception:
                        pass
                if blast_bitscore is not None and blast_bitscore_min is not None and float(blast_bitscore_min) > 0:
                    try:
                        ok = ok and (float(blast_bitscore) >= float(blast_bitscore_min))
                    except Exception:
                        pass
                if ok and blast_evalue is not None:
                    blast_hit = 1
                    try:
                        ev = float(blast_evalue)
                        if ev <= 0:
                            ev = 1e-300
                        blast_score = _norm01(-math.log10(ev), vmax_blast)
                    except Exception:
                        blast_score = 0.0
                    score += w_blast * blast_score

            # eggNOG evidence
            eggnog_hit = 0
            eggnog_desc = ''
            eggnog_go = ''
            eggnog_ogs = ''
            eggnog_pref_name = ''
            eggnog_kegg_ko = ''
            eggnog_pfam = ''
            if eggnog_map and g.gene_id in eggnog_map:
                em = eggnog_map[g.gene_id]
                eggnog_hit = 1
                eggnog_desc = str(em.get('eggnog_desc') or '')
                eggnog_go = str(em.get('eggnog_go') or '')
                eggnog_ogs = str(em.get('eggnog_ogs') or '')
                eggnog_pref_name = str(em.get('eggnog_pref_name') or '')
                eggnog_kegg_ko = str(em.get('eggnog_kegg_ko') or '')
                eggnog_pfam = str(em.get('eggnog_pfam') or '')
                score += w_eggnog * 1.0

            # GO keyword evidence (search in eggNOG Description and GO IDs)
            go_hit = 0
            go_matches: List[str] = []
            go_score = 0.0
            if go_keywords:
                go_ids = set()
                go_ids |= _extract_go_ids(eggnog_go)
                if go_map and g.gene_id in go_map:
                    go_ids |= set([x.upper() for x in (go_map.get(g.gene_id) or [])])
                text = (eggnog_desc or '')
                text_l = text.lower()
                for kw in go_keywords:
                    k = str(kw).strip()
                    if not k:
                        continue
                    if re.match(r'^GO:\d+$', k.strip(), flags=re.I):
                        if k.upper() in go_ids:
                            go_matches.append(k)
                    else:
                        if k.lower() in text_l:
                            go_matches.append(k)
                if go_matches:
                    go_hit = 1
                    go_score = 1.0
                    score += w_go * go_score

                    score += w_deg * deg_score

            cand_rows.append({
                "locus_id": L.locus_id,
                "chr": L.chr,
                "locus_start": L.start,
                "locus_end": L.end,
                "peak_pos": L.rep_pos,
                "peak_source": L.rep_source,
                "peak_pvalue": L.rep_pvalue,
                "peak_lod": L.rep_lod,
                "peak_strength": rep_strength,
                "gene_id": g.gene_id,
                "gene_start": g.start,
                "gene_end": g.end,
                "strand": g.strand,
                "gene_source": g.source,
                "distance": dist,
                "dist_score": dist_score,
                "peak_score_norm": rep_peak_norm,
                "deg_hit": deg_hit,
                "deg_padj": deg_padj,
                "deg_logFC": deg_lfc,
                "deg_score_norm": deg_score,
                "blast_hit": blast_hit,
                "blast_evalue": blast_evalue,
                "blast_bitscore": blast_bitscore,
                "blast_best_hit": blast_best_hit,
                "blast_score_norm": blast_score,
                "eggnog_hit": eggnog_hit,
                "eggnog_desc": eggnog_desc,
                "eggnog_go": eggnog_go,
                "eggnog_ogs": eggnog_ogs,
                "eggnog_pref_name": eggnog_pref_name,
                "eggnog_kegg_ko": eggnog_kegg_ko,
                "eggnog_pfam": eggnog_pfam,
                "go_hit": go_hit,
                "go_matches": ";".join(go_matches) if go_matches else "",
                "go_score": go_score,
                "score": score,
            })

    counts["n_candidates"] = len(cand_rows)

    # gene-wise best score (one line per locus-gene is also helpful, but ranked table should be by gene overall)
    # We'll keep both: candidates_ranked.tsv is sorted by score and includes locus_id.

    # compute rank
    cand_rows.sort(key=lambda d: (float(d.get("score") or 0.0)), reverse=True)
    for i, d in enumerate(cand_rows, start=1):
        d["rank"] = i

    return cand_rows, counts


# -------------------------
# Main
# -------------------------


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--params", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    params = json.loads(Path(args.params).read_text(encoding="utf-8"))
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    gwas_results = [Path(p) for p in (params.get("gwas_results") or []) if str(p).strip()]
    qtl_peaks = [Path(p) for p in (params.get("qtl_peaks") or []) if str(p).strip()]

    mode_gwas = str(params.get("mode_gwas") or "top_n").strip()
    gwas_top_n = int(params.get("gwas_top_n") or 200)
    gwas_p_threshold = _safe_float(params.get("gwas_p_threshold"))
    if gwas_p_threshold is None:
        gwas_p_threshold = 1e-6

    mode_qtl = str(params.get("mode_qtl") or "top_n").strip()
    qtl_top_n = int(params.get("qtl_top_n") or 50)
    qtl_lod_threshold = _safe_float(params.get("qtl_lod_threshold"))
    if qtl_lod_threshold is None:
        qtl_lod_threshold = 3.0

    window_bp = int(params.get("window_bp") or 50000)
    normalize_chr = bool(params.get("chr_normalize", True))

    # Optional: marker_map.tsv (marker/chr/pos_bp) to convert QTL coordinates (cM/Mb) to bp.
    marker_map_path: Optional[Path] = None
    try:
        mm = str(params.get("marker_map_tsv") or "").strip()
        if mm:
            p = Path(mm)
            if p.exists():
                marker_map_path = p
    except Exception:
        marker_map_path = None

    if marker_map_path is None and qtl_peaks:
        # Auto-detect near the first QTL peaks file
        marker_map_path = _auto_find_marker_map(qtl_peaks[0])

    marker_map: Optional[Dict[str, tuple]] = None
    if marker_map_path is not None and marker_map_path.exists():
        try:
            marker_map = load_marker_map(marker_map_path, normalize_chr=normalize_chr)
        except Exception:
            marker_map = None

    genome_gff = Path(str(params.get("genome_gff") or "").strip()) if str(params.get("genome_gff") or "").strip() else None
    predicted_gff = Path(str(params.get("predicted_gff") or "").strip()) if str(params.get("predicted_gff") or "").strip() else None

    deg_list = [Path(p) for p in (params.get("deg_tsv_list") or []) if str(p).strip()]
    deg_padj_max = _safe_float(params.get("deg_padj_max"))
    if deg_padj_max is None:
        deg_padj_max = 0.05
    deg_abs_logfc_min = _safe_float(params.get("deg_abs_logfc_min"))
    if deg_abs_logfc_min is None:
        deg_abs_logfc_min = 0.0

    max_genes_per_locus = int(params.get("max_genes_per_locus") or 50)
    dist_decay_bp = int(params.get("dist_decay_bp") or 50000)

    w_peak = float(_safe_float(params.get("w_peak")) or 1.0)
    w_dist = float(_safe_float(params.get("w_dist")) or 1.0)
    w_deg = float(_safe_float(params.get("w_deg")) or 1.0)

    # BLAST / eggNOG / GO evidence
    blast_list = [Path(p) for p in (params.get('blast_tsv_list') or []) if str(p).strip()]
    blast_evalue_max = float(_safe_float(params.get('blast_evalue_max')) or 1e-5)
    blast_bitscore_min = float(_safe_float(params.get('blast_bitscore_min')) or 0.0)
    w_blast = float(_safe_float(params.get('w_blast')) or 1.0)

    eggnog_list = [Path(p) for p in (params.get('eggnog_annotations_list') or []) if str(p).strip()]
    w_eggnog = float(_safe_float(params.get('w_eggnog')) or 0.5)

    go_map_file = Path(str(params.get('go_map_file') or '').strip()) if str(params.get('go_map_file') or '').strip() else None
    go_keywords = _parse_keywords(str(params.get('go_keywords') or '').strip())
    w_go = float(_safe_float(params.get('w_go')) or 0.5)

    # --- peaks ---
    peaks: List[Peak] = []
    for fp in gwas_results:
        if fp.exists():
            peaks.extend(_iter_gwas_peaks(fp, mode_gwas, gwas_top_n, float(gwas_p_threshold), window_bp, normalize_chr))
    for fp in qtl_peaks:
        if fp.exists():
            peaks.extend(_iter_qtl_peaks(fp, mode_qtl, qtl_top_n, float(qtl_lod_threshold), window_bp, normalize_chr, marker_map=marker_map))

    # --- loci ---
    loci, locus_peaks = build_loci(peaks)

    # peaks_merged.tsv (rep peaks)
    peaks_merged_rows: List[List[Any]] = []
    for L in loci:
        peaks_merged_rows.append([L.locus_id, L.chr, L.rep_pos, L.rep_source, L.rep_pvalue, L.rep_lod])
    write_tsv(out_dir / "peaks_merged.tsv", ["locus_id", "chr", "pos", "source", "pvalue", "lod"], peaks_merged_rows)

    # locus_peaks.tsv
    lp_rows: List[List[Any]] = []
    for p in locus_peaks:
        # reassign locus_id by overlap
        lid = ""
        for L in loci:
            if L.chr == p.chr and L.start <= p.pos <= L.end:
                lid = L.locus_id
                break
        lp_rows.append([lid, p.chr, p.pos, p.source_type, p.dataset, p.pvalue, p.lod, p.raw_file])
    write_tsv(out_dir / "locus_peaks.tsv", ["locus_id", "chr", "pos", "source_type", "dataset", "pvalue", "lod", "raw_file"], lp_rows)

    # loci.tsv
    loci_rows: List[List[Any]] = []
    for L in loci:
        # count peaks
        n_peaks = sum(1 for r in lp_rows if r[0] == L.locus_id)
        loci_rows.append([L.locus_id, L.chr, L.start, L.end, L.rep_pos, L.rep_source, L.rep_pvalue, L.rep_lod, n_peaks])
    write_tsv(out_dir / "loci.tsv", ["locus_id", "chr", "start", "end", "rep_pos", "rep_source", "rep_pvalue", "rep_lod", "n_peaks"], loci_rows)

    # --- genes ---
    genes: List[Gene] = []
    if genome_gff is not None:
        genes.extend(read_genes_from_gff(genome_gff, source="gff", normalize_chr=normalize_chr))
    if predicted_gff is not None:
        genes.extend(read_genes_from_gff(predicted_gff, source="pred", normalize_chr=normalize_chr))

    # --- evidence ---
    deg_map = load_deg_tables(deg_list) if deg_list else {}
    blast_map = load_blast_tables(blast_list) if blast_list else {}
    eggnog_map = load_eggnog_tables(eggnog_list) if eggnog_list else {}
    go_map = load_go_map(go_map_file) if go_map_file else {}

    # --- prioritize ---
    cand_rows, counts = prioritize(
        loci=loci,
        locus_peaks=locus_peaks,
        genes=genes,
        deg_map=deg_map,
        blast_map=blast_map,
        blast_evalue_max=float(blast_evalue_max),
        blast_bitscore_min=float(blast_bitscore_min),
        eggnog_map=eggnog_map,
        go_map=go_map,
        go_keywords=go_keywords,
        max_genes_per_locus=max_genes_per_locus,
        dist_decay_bp=dist_decay_bp,
        deg_padj_max=float(deg_padj_max),
        deg_abs_logfc_min=float(deg_abs_logfc_min),
        w_peak=float(w_peak),
        w_dist=float(w_dist),
        w_deg=float(w_deg),
        w_blast=float(w_blast),
        w_eggnog=float(w_eggnog),
        w_go=float(w_go),
    )

    # candidates_ranked.tsv
    header = [
        "rank", "score",
        "gene_id", "gene_source", "chr", "gene_start", "gene_end", "strand",
        "locus_id", "locus_start", "locus_end",
        "peak_pos", "peak_source", "peak_pvalue", "peak_lod", "peak_strength", "peak_score_norm",
        "distance", "dist_score",
        "deg_hit", "deg_padj", "deg_logFC", "deg_score_norm",
        "blast_hit", "blast_best_hit", "blast_evalue", "blast_bitscore", "blast_score_norm",
        "eggnog_hit", "eggnog_desc", "eggnog_go", "eggnog_ogs", "eggnog_pref_name", "eggnog_kegg_ko", "eggnog_pfam",
        "go_hit", "go_matches", "go_score",
    ]
    rows = []
    for d in cand_rows:
        rows.append([
            d.get("rank"), d.get("score"),
            d.get("gene_id"), d.get("gene_source"), d.get("chr"), d.get("gene_start"), d.get("gene_end"), d.get("strand"),
            d.get("locus_id"), d.get("locus_start"), d.get("locus_end"),
            d.get("peak_pos"), d.get("peak_source"), d.get("peak_pvalue"), d.get("peak_lod"), d.get("peak_strength"), d.get("peak_score_norm"),
            d.get("distance"), d.get("dist_score"),
            d.get("deg_hit"), d.get("deg_padj"), d.get("deg_logFC"), d.get("deg_score_norm"),
            d.get("blast_hit"), d.get("blast_best_hit"), d.get("blast_evalue"), d.get("blast_bitscore"), d.get("blast_score_norm"),
            d.get("eggnog_hit"), d.get("eggnog_desc"), d.get("eggnog_go"), d.get("eggnog_ogs"), d.get("eggnog_pref_name"), d.get("eggnog_kegg_ko"), d.get("eggnog_pfam"),
            d.get("go_hit"), d.get("go_matches"), d.get("go_score"),
        ])
    write_tsv(out_dir / "candidates_ranked.tsv", header, rows)

    meta = {
        "params": params,
        "counts": counts,
        "n_peaks": len(peaks),
        "n_loci": len(loci),
        "n_genes": len(genes),
        "n_deg_genes": len(deg_map),
        "n_blast_genes": len(blast_map),
        "n_eggnog_genes": len(eggnog_map),
        "n_go_map_genes": len(go_map),
        "outputs": {
            "loci_tsv": "loci.tsv",
            "locus_peaks_tsv": "locus_peaks.tsv",
            "peaks_merged_tsv": "peaks_merged.tsv",
            "candidates_ranked_tsv": "candidates_ranked.tsv",
        },
    }
    (out_dir / "meta.json").write_text(json.dumps(meta, ensure_ascii=False, indent=2), encoding="utf-8")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
