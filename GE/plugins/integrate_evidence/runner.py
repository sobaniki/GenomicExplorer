#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""integrate_evidence (Phase 1)

Goal
----
Create locus-wise evidence tables by joining:
  - integrate_loci outputs (loci.tsv, locus_peaks.tsv)
  - Annotate/GenePred/Candidates gene list TSV
  - optional DEG TSVs (multiple)
  - optional BLAST TSVs (multiple)

This phase prioritizes robustness and minimal assumptions.
If a join key is unavailable, the script falls back to best-effort heuristics.

Outputs (written under --out)
--------------------------------
  - locus_evidence.tsv          : 1 row per locus (aggregated summary)
  - locus_gene_evidence.tsv     : 1 row per (locus, gene)
  - evidence_sources.tsv        : record of loaded files + guessed schemas
  - meta.json                   : params + basic counts

Expected columns (best effort)
------------------------------
 loci.tsv:
   locus_id, chr, start, end, rep_pos, rep_source, rep_pvalue, rep_lod
 locus_peaks.tsv:
   locus_id, chr, pos, source_type (gwas/qtl), pvalue, lod
 annotate genes TSV:
   Ideally contains either:
     (A) locus_id + gene_id
     (B) chr + pos (+/- peak_chr/peak_pos) + gene_id
     (C) gene_id + gene_start/gene_end + chr
 DEG TSV:
   gene_id + (logFC, padj/FDR)
 BLAST TSV:
   gene_id or qseqid/query (+ evalue/bitscore)

"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple


def _sniff_delim(path: Path) -> str:
    if path.suffix.lower() == ".csv":
        return ","
    # try to sniff from first line
    try:
        head = path.read_text(encoding="utf-8", errors="ignore").splitlines()[:5]
        sample = "\n".join(head)
        dialect = csv.Sniffer().sniff(sample, delimiters=["\t", ",", ";", " "])
        return dialect.delimiter
    except Exception:
        return "\t"


def read_table(path: Path) -> List[Dict[str, str]]:
    """Read TSV/CSV as list of dicts. Uses pandas if available, else csv module."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(str(path))

    # Try pandas
    try:
        import pandas as pd  # type: ignore

        sep = _sniff_delim(path)
        df = pd.read_csv(path, sep=sep, dtype=str, comment=None)
        # normalize column names
        df.columns = [str(c).strip() for c in df.columns]
        rows = df.fillna("").to_dict(orient="records")
        return [{str(k): str(v) for k, v in r.items()} for r in rows]
    except Exception:
        # csv fallback
        sep = _sniff_delim(path)
        with path.open("r", encoding="utf-8", errors="ignore", newline="") as f:
            reader = csv.DictReader(f, delimiter=sep)
            out = []
            for row in reader:
                if row is None:
                    continue
                out.append({(k or "").strip(): (v or "").strip() for k, v in row.items()})
        return out


def write_tsv(path: Path, header: List[str], rows: List[List[Any]]):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(header)
        for r in rows:
            w.writerow(["" if x is None else x for x in r])


def _get_any(d: Dict[str, str], keys: List[str]) -> str:
    for k in keys:
        if k in d and str(d[k]).strip() != "":
            return str(d[k]).strip()
    return ""


def _to_float(x: str) -> Optional[float]:
    if x is None:
        return None
    s = str(x).strip()
    if s == "" or s.lower() in {"na", "nan", "none"}:
        return None
    try:
        return float(s)
    except Exception:
        # handle scientific notation in string form
        try:
            return float(s.replace(",", ""))
        except Exception:
            return None


def _best_blast_of_rows(rows: List[Dict[str, str]]) -> Tuple[Optional[float], str, str]:
    """Return (best_evalue, best_hit, best_bitscore) from blast rows."""
    best_e = None
    best_hit = ""
    best_bs = ""
    for r in rows:
        e = _to_float(_get_any(r, ["evalue", "e-value", "Evalue", "E-value"]))
        bs = _get_any(r, ["bitscore", "bit_score", "score", "Score"])  # keep string
        hit = _get_any(r, ["sseqid", "subject", "subject_id", "hit", "sacc"]) or _get_any(r, ["subject_id"]) 
        if e is None:
            continue
        if best_e is None or e < best_e:
            best_e = e
            best_hit = hit
            best_bs = bs
    return best_e, best_hit, best_bs


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--params", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    params_p = Path(args.params)
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    params = json.loads(params_p.read_text(encoding="utf-8"))

    loci_tsv = Path(params.get("loci_tsv", ""))
    locus_peaks_tsv = Path(params.get("locus_peaks_tsv", ""))
    annotate_genes_tsv = Path(params.get("annotate_genes_tsv", ""))

    deg_list = [Path(x) for x in (params.get("deg_tsv_list") or []) if str(x).strip()]
    blast_list = [Path(x) for x in (params.get("blast_tsv_list") or []) if str(x).strip()]

    window_bp = int(params.get("window_bp", 0) or 0)

    # --- load core tables ---
    loci_rows = read_table(loci_tsv) if str(loci_tsv).strip() else []
    peaks_rows = read_table(locus_peaks_tsv) if str(locus_peaks_tsv).strip() else []
    anno_rows = read_table(annotate_genes_tsv) if str(annotate_genes_tsv).strip() else []

    # loci map
    loci: Dict[str, Dict[str, Any]] = {}
    for r in loci_rows:
        locus_id = _get_any(r, ["locus_id", "Locus", "locus"]) or ""
        if not locus_id:
            continue
        chr_ = _get_any(r, ["chr", "chrom", "chromosome", "seqnames"])
        start = _to_float(_get_any(r, ["start", "locus_start", "from"]))
        end = _to_float(_get_any(r, ["end", "locus_end", "to"]))
        rep_pos = _to_float(_get_any(r, ["rep_pos", "pos", "peak_pos"]))
        rep_source = _get_any(r, ["rep_source", "source", "rep_type"]) 
        rep_p = _to_float(_get_any(r, ["rep_pvalue", "pvalue", "p", "min_p"]))
        rep_lod = _to_float(_get_any(r, ["rep_lod", "lod", "max_lod"]))
        loci[locus_id] = {
            "locus_id": locus_id,
            "chr": chr_,
            "start": int(start) if start is not None else "",
            "end": int(end) if end is not None else "",
            "rep_pos": int(rep_pos) if rep_pos is not None else "",
            "rep_source": rep_source,
            "rep_pvalue": rep_p,
            "rep_lod": rep_lod,
        }

    # Peaks by locus for (chr,pos) matching
    peak_index: Dict[Tuple[str, str], str] = {}  # (chr,pos)->locus_id
    peaks_by_locus: Dict[str, List[Dict[str, Any]]] = {}
    for r in peaks_rows:
        locus_id = _get_any(r, ["locus_id", "locus"]) or ""
        chr_ = _get_any(r, ["chr", "chrom", "chromosome"])
        pos = _get_any(r, ["pos", "position", "bp"])
        if locus_id and chr_ and pos:
            peak_index[(chr_, pos)] = locus_id
        peaks_by_locus.setdefault(locus_id or "", []).append(r)

    # --- parse annotate genes: assign to loci ---
    gene_records: List[Dict[str, Any]] = []
    for r in anno_rows:
        gene_id = _get_any(r, ["gene_id", "gene", "ID", "id", "GeneID"])
        if not gene_id:
            # sometimes "query" field for candidates
            gene_id = _get_any(r, ["query", "qseqid"]) 
        if not gene_id:
            continue

        locus_id = _get_any(r, ["locus_id", "Locus", "locus"])  # best

        # try peak match
        if not locus_id:
            chr_ = _get_any(r, ["peak_chr", "chr", "chrom", "chromosome"])
            pos = _get_any(r, ["peak_pos", "pos", "position", "bp"])
            if chr_ and pos and (chr_, pos) in peak_index:
                locus_id = peak_index[(chr_, pos)]

        # try coordinate overlap
        if not locus_id:
            chr_ = _get_any(r, ["chr", "chrom", "chromosome", "seqnames"])
            gs = _to_float(_get_any(r, ["gene_start", "start", "tx_start"]))
            ge = _to_float(_get_any(r, ["gene_end", "end", "tx_end"]))
            if chr_ and gs is not None and ge is not None and loci:
                for lid, L in loci.items():
                    if str(L.get("chr", "")) != chr_:
                        continue
                    ls = L.get("start")
                    le = L.get("end")
                    if isinstance(ls, int) and isinstance(le, int):
                        if not (ge < ls or gs > le):
                            locus_id = lid
                            break
                    elif window_bp and L.get("rep_pos") and isinstance(L.get("rep_pos"), int):
                        rp = int(L["rep_pos"])
                        if abs(((gs + ge) / 2.0) - rp) <= window_bp:
                            locus_id = lid
                            break

        rec = {
            "locus_id": locus_id,
            "gene_id": gene_id,
            "anno_desc": _get_any(r, ["desc", "description", "product", "annotation", "Note"]) ,
            "raw": r,
        }
        # gene coordinates if present
        rec["gene_chr"] = _get_any(r, ["chr", "chrom", "chromosome", "seqnames"]) 
        rec["gene_start"] = _to_float(_get_any(r, ["gene_start", "start"]))
        rec["gene_end"] = _to_float(_get_any(r, ["gene_end", "end"]))
        gene_records.append(rec)

    # index genes by gene_id
    genes_by_id: Dict[str, List[Dict[str, Any]]] = {}
    for g in gene_records:
        genes_by_id.setdefault(g["gene_id"], []).append(g)

    # --- DEG merge ---
    deg_by_gene: Dict[str, Dict[str, Any]] = {}
    for fp in deg_list:
        if not fp.exists():
            continue
        rows = read_table(fp)
        for r in rows:
            gid = _get_any(r, ["gene_id", "gene", "id", "ID", "GeneID"])
            if not gid:
                continue
            logfc = _to_float(_get_any(r, ["logFC", "log2FoldChange", "log2FC", "LFC"]))
            fdr = _to_float(_get_any(r, ["padj", "FDR", "qvalue", "adj.P.Val"]))
            p = _to_float(_get_any(r, ["pvalue", "p", "P.Value"]))
            # keep best (min FDR) across multiple DEG files
            cur = deg_by_gene.get(gid)
            if cur is None or (fdr is not None and (cur.get("deg_fdr") is None or fdr < cur.get("deg_fdr"))):
                deg_by_gene[gid] = {
                    "deg_logfc": logfc,
                    "deg_fdr": fdr,
                    "deg_pvalue": p,
                    "deg_source": str(fp),
                }

    # --- BLAST merge ---
    blast_by_gene: Dict[str, Dict[str, Any]] = {}
    for fp in blast_list:
        if not fp.exists():
            continue
        rows = read_table(fp)
        # group by query/gene
        tmp: Dict[str, List[Dict[str, str]]] = {}
        for r in rows:
            q = _get_any(r, ["gene_id", "qseqid", "query", "query_id", "qacc"]) 
            if not q:
                continue
            tmp.setdefault(q, []).append(r)
        for q, rs in tmp.items():
            best_e, best_hit, best_bs = _best_blast_of_rows(rs)
            cur = blast_by_gene.get(q)
            if cur is None or (best_e is not None and (cur.get("blast_evalue") is None or best_e < cur.get("blast_evalue"))):
                blast_by_gene[q] = {
                    "blast_evalue": best_e,
                    "blast_best_hit": best_hit,
                    "blast_bitscore": best_bs,
                    "blast_source": str(fp),
                }

    # --- build locus_gene_evidence ---
    lg_rows: List[List[Any]] = []
    locus_stats: Dict[str, Dict[str, Any]] = {}

    for g in gene_records:
        lid = g.get("locus_id") or ""
        gid = g.get("gene_id")
        if not gid:
            continue
        L = loci.get(lid, {})
        rp = L.get("rep_pos")
        dist = ""
        if rp and g.get("gene_start") is not None and g.get("gene_end") is not None:
            mid = (float(g["gene_start"]) + float(g["gene_end"])) / 2.0
            dist = int(abs(mid - float(rp)))

        deg = deg_by_gene.get(gid, {})
        bl = blast_by_gene.get(gid, {})

        lg_rows.append([
            lid,
            gid,
            dist,
            deg.get("deg_logfc", ""),
            deg.get("deg_fdr", ""),
            bl.get("blast_best_hit", ""),
            bl.get("blast_evalue", ""),
            bl.get("blast_bitscore", ""),
            g.get("anno_desc", ""),
        ])

        st = locus_stats.setdefault(lid, {"n_genes": 0, "n_deg": 0, "n_blast": 0, "best_fdr": None, "best_e": None, "top_gene": ""})
        st["n_genes"] += 1
        if deg.get("deg_fdr") is not None:
            st["n_deg"] += 1
            if st["best_fdr"] is None or deg.get("deg_fdr") < st["best_fdr"]:
                st["best_fdr"] = deg.get("deg_fdr")
                st["top_gene"] = gid
        if bl.get("blast_evalue") is not None:
            st["n_blast"] += 1
            if st["best_e"] is None or bl.get("blast_evalue") < st["best_e"]:
                st["best_e"] = bl.get("blast_evalue")

    # add loci that had no genes
    for lid in loci.keys():
        locus_stats.setdefault(lid, {"n_genes": 0, "n_deg": 0, "n_blast": 0, "best_fdr": None, "best_e": None, "top_gene": ""})

    # --- build locus_evidence ---
    le_rows: List[List[Any]] = []
    for lid, L in loci.items():
        st = locus_stats.get(lid, {})
        best_p = L.get("rep_pvalue")
        best_lod = L.get("rep_lod")

        # simple evidence_score (optional): higher is better
        score = 0.0
        if isinstance(best_p, (int, float)) and best_p and best_p > 0:
            score += -math.log10(float(best_p))
        if isinstance(best_lod, (int, float)):
            score += float(best_lod)
        bfdr = st.get("best_fdr")
        if isinstance(bfdr, (int, float)) and bfdr and bfdr > 0:
            score += -math.log10(float(bfdr))
        be = st.get("best_e")
        if isinstance(be, (int, float)) and be and be > 0:
            score += -math.log10(float(be)) * 0.25

        le_rows.append([
            lid,
            L.get("chr", ""),
            L.get("start", ""),
            L.get("end", ""),
            L.get("rep_pos", ""),
            L.get("rep_source", ""),
            best_p if best_p is not None else "",
            best_lod if best_lod is not None else "",
            st.get("n_genes", 0),
            st.get("n_deg", 0),
            st.get("n_blast", 0),
            st.get("top_gene", ""),
            st.get("best_fdr") if st.get("best_fdr") is not None else "",
            st.get("best_e") if st.get("best_e") is not None else "",
            round(score, 6),
        ])

    # sort by score desc
    le_rows.sort(key=lambda r: (r[-1] if isinstance(r[-1], (int, float)) else float(r[-1] or 0)), reverse=True)

    # --- sources ---
    sources: List[List[Any]] = []
    def _add_source(kind: str, fp: Path, n: int, note: str):
        sources.append([kind, str(fp), n, note])

    _add_source("loci", loci_tsv, len(loci_rows), "")
    _add_source("locus_peaks", locus_peaks_tsv, len(peaks_rows), "")
    _add_source("annotate_genes", annotate_genes_tsv, len(anno_rows), "")
    for fp in deg_list:
        _add_source("deg", fp, 0, "loaded")
    for fp in blast_list:
        _add_source("blast", fp, 0, "loaded")

    # --- write outputs ---
    write_tsv(out_dir / "locus_gene_evidence.tsv",
              ["locus_id", "gene_id", "distance_to_rep", "deg_logFC", "deg_FDR", "blast_best_hit", "blast_evalue", "blast_bitscore", "anno_desc"],
              lg_rows)

    write_tsv(out_dir / "locus_evidence.tsv",
              ["locus_id", "chr", "start", "end", "rep_pos", "rep_source", "rep_pvalue", "rep_lod",
               "n_genes", "n_deg", "n_blast", "top_gene", "best_deg_FDR", "best_blast_evalue", "evidence_score"],
              le_rows)

    write_tsv(out_dir / "evidence_sources.tsv",
              ["kind", "path", "n_rows", "note"],
              sources)

    meta = {
        "params": params,
        "counts": {
            "n_loci": len(loci),
            "n_locus_peaks": len(peaks_rows),
            "n_annotate_rows": len(anno_rows),
            "n_genes": len(gene_records),
            "n_deg_genes": len(deg_by_gene),
            "n_blast_genes": len(blast_by_gene),
        },
    }
    (out_dir / "meta.json").write_text(json.dumps(meta, indent=2, ensure_ascii=False), encoding="utf-8")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
