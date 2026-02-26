#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 30 10:46:03 2026

@author: soba
"""

#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd

def load_chr_map(path: Optional[str]) -> Dict[str, str]:
    if not path:
        return {}
    p = Path(path)
    if not p.exists():
        return {}
    df = pd.read_csv(p, sep="\t")
    if df.shape[1] < 2:
        return {}
    # prefer columns named from/to, otherwise first 2 cols
    if "from" in df.columns and "to" in df.columns:
        a = df["from"].astype(str).tolist()
        b = df["to"].astype(str).tolist()
    else:
        a = df.iloc[:, 0].astype(str).tolist()
        b = df.iloc[:, 1].astype(str).tolist()
    return {x.strip(): y.strip() for x, y in zip(a, b) if str(x).strip()}


def parse_attrs(attr: str) -> Dict[str, str]:
    """
    Parse GFF3 attributes: key=value;key2=value2
    Parse GTF attributes: key "value"; key2 "value2";
    """
    attr = (attr or "").strip()
    d: Dict[str, str] = {}
    if not attr:
        return d

    parts = [p.strip() for p in attr.strip(";").split(";") if p.strip()]
    for p in parts:
        if "=" in p:  # GFF3
            k, v = p.split("=", 1)
            d[k.strip()] = v.strip()
        else:  # maybe GTF: key "value"
            # split by first space
            sp = p.split(None, 1)
            if len(sp) == 2:
                k, v = sp[0].strip(), sp[1].strip()
                v = v.strip().strip('"').strip()
                d[k] = v
    return d


def choose_gene_id(attrs: Dict[str, str], keys: List[str]) -> str:
    for k in keys:
        if k in attrs and attrs[k]:
            return attrs[k]
    # fallback: try any common keys
    for k in ["ID", "gene_id", "Name", "locus_tag", "Parent"]:
        if k in attrs and attrs[k]:
            return attrs[k]
    return ""


def norm_chr(x: str, chr_map: Dict[str, str] | None = None) -> str:
    x = str(x).strip()
    if not x:
        return x
    lx = x.lower()
    if lx.startswith("chr"):
        x = x[3:]
    x = x.strip()
    if chr_map:
        return chr_map.get(x, x)
    return x



def read_gff_genes(
    gff_path: Path,
    feature_types: List[str],
    gene_id_keys: List[str],
    normalize_chr: bool = True,
    chr_map: Dict[str, str] | None = None,
) -> pd.DataFrame:
    rows = []
    with gff_path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            seqid, source, ftype, start, end, score, strand, phase, attr = parts[:9]
            if feature_types and ftype not in feature_types:
                continue
            try:
                st = int(float(start))
                en = int(float(end))
            except Exception:
                continue
            attrs = parse_attrs(attr)
            gid = choose_gene_id(attrs, gene_id_keys)
            name = attrs.get("Name", "") or attrs.get("gene_name", "") or attrs.get("product", "")
            #chr_ = norm_chr(seqid) if normalize_chr else str(seqid)
            
            chr_ = norm_chr(seqid, chr_map) if normalize_chr else str(seqid)


            rows.append(
                {
                    "chr": chr_,
                    "start": min(st, en),
                    "end": max(st, en),
                    "strand": strand,
                    "feature": ftype,
                    "gene_id": gid,
                    "gene_name": name,
                    "attrs": attr,
                }
            )
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    df = df.sort_values(["chr", "start", "end"], kind="mergesort")
    return df


def interval_distance(pos: int, start: int, end: int) -> Tuple[int, bool]:
    """distance to interval; 0 if inside; overlap flag"""
    if start <= pos <= end:
        return 0, True
    if pos < start:
        return start - pos, False
    return pos - end, False


def annotate_peaks(
    peaks: pd.DataFrame,
    genes: pd.DataFrame,
    chr_col: str,
    pos_col: str,
    flank_bp: int,
    normalize_chr: bool = True,
    closest_n: int = 0,
    chr_map: Dict[str, str] | None = None,
) -> pd.DataFrame:
    out_rows = []

    if genes.empty:
        return pd.DataFrame()

    # index genes by chr for speed
    g_by_chr = {c: sub for c, sub in genes.groupby("chr", sort=False)}

    for i, row in peaks.iterrows():
        chr_raw = row.get(chr_col, "")
        pos_raw = row.get(pos_col, None)

        if pd.isna(pos_raw):
            continue
        try:
            pos = int(float(pos_raw))
        except Exception:
            continue

        #chr_ = norm_chr(chr_raw) if normalize_chr else str(chr_raw)
        chr_ = norm_chr(chr_raw, chr_map) if normalize_chr else str(chr_raw)

        
        sub = g_by_chr.get(chr_)
        if sub is None or sub.empty:
            continue

        win_start = pos - flank_bp
        win_end = pos + flank_bp

        # genes within window (overlap or within flank)
        inwin = sub[(sub["end"] >= win_start) & (sub["start"] <= win_end)].copy()

        if inwin.empty and closest_n and closest_n > 0:
            # take closest genes on same chr
            tmp = sub.copy()
            tmp["distance"] = tmp.apply(lambda r: interval_distance(pos, int(r["start"]), int(r["end"]))[0], axis=1)
            tmp = tmp.sort_values("distance").head(int(closest_n))
            inwin = tmp

        for _, g in inwin.iterrows():
            dist, ov = interval_distance(pos, int(g["start"]), int(g["end"]))
            out_rows.append(
                {
                    "peak_index": int(i),
                    "peak_chr": chr_,
                    "peak_pos": pos,
                    "window_bp": flank_bp,
                    "gene_id": g["gene_id"],
                    "gene_name": g["gene_name"],
                    "gene_chr": g["chr"],
                    "gene_start": int(g["start"]),
                    "gene_end": int(g["end"]),
                    "strand": g["strand"],
                    "feature": g["feature"],
                    "distance_bp": int(dist),
                    "overlap": bool(ov),
                }
            )

    df = pd.DataFrame(out_rows)
    if not df.empty:
        df = df.sort_values(["peak_chr", "peak_pos", "distance_bp", "gene_start"], kind="mergesort")
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--params", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    out_dir = Path(args.out).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    p = json.loads(Path(args.params).read_text(encoding="utf-8"))
    chr_map_tsv = p.get("chr_map_tsv", None)
    
    chr_map = load_chr_map(chr_map_tsv)

    gff_path = Path(p["gff_path"]).expanduser()
    peaks_tsv = Path(p["peaks_tsv"]).expanduser()

    flank_bp = int(p.get("flank_bp", 50000))
    feature_types = list(p.get("feature_types", ["gene"]))
    gene_id_keys = list(p.get("gene_id_keys", ["ID", "gene_id", "locus_tag", "Name"]))
    normalize_chr = bool(p.get("normalize_chr", True))
    closest_n = int(p.get("closest_n", 0))

    chr_col = str(p.get("chr_col", "chr"))
    pos_col = str(p.get("pos_col", "pos"))

    peaks = pd.read_csv(peaks_tsv, sep="\t")
    if chr_col not in peaks.columns or pos_col not in peaks.columns:
        raise SystemExit(f"peaks_tsv must contain columns: {chr_col}, {pos_col}")

    genes = read_gff_genes(gff_path, feature_types, gene_id_keys, normalize_chr=normalize_chr, chr_map=chr_map)

    per_peak = annotate_peaks(
        peaks=peaks,
        genes=genes,
        chr_col=chr_col,
        pos_col=pos_col,
        flank_bp=flank_bp,
        normalize_chr=normalize_chr,
        closest_n=closest_n,
        chr_map=chr_map,
    )

    per_peak_path = out_dir / "per_peak_genes.tsv"
    uniq_path = out_dir / "genes_unique.tsv"

    if per_peak.empty:
        per_peak.to_csv(per_peak_path, sep="\t", index=False)
        pd.DataFrame(columns=["gene_id", "gene_name", "gene_chr", "gene_start", "gene_end", "strand", "feature"]).to_csv(
            uniq_path, sep="\t", index=False
        )
    else:
        per_peak.to_csv(per_peak_path, sep="\t", index=False)

        uniq = (
            per_peak[["gene_id", "gene_name", "gene_chr", "gene_start", "gene_end", "strand", "feature"]]
            .drop_duplicates()
            .sort_values(["gene_chr", "gene_start", "gene_end"], kind="mergesort")
        )
        uniq.to_csv(uniq_path, sep="\t", index=False)

    artifacts = {
        "per_peak_genes_tsv": str(per_peak_path),
        "genes_unique_tsv": str(uniq_path),
        "n_genes_total_in_gff": int(len(genes)) if genes is not None else 0,
        "n_rows_per_peak": int(len(per_peak)) if per_peak is not None else 0,
    }
    (out_dir / "artifacts.json").write_text(json.dumps(artifacts, indent=2, ensure_ascii=False), encoding="utf-8")


if __name__ == "__main__":
    main()
