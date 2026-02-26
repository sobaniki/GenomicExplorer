#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""integrate_loci

Build loci by clustering +/- window_bp intervals around peak positions derived from:
  - multiple GWAS result tables (tsv/csv)
  - multiple QTL peaks tables (tsv/csv)

Outputs (under --out):
  - loci.tsv
  - locus_peaks.tsv
  - peaks_merged.tsv  (representative peak per locus; chr,pos columns)

Design goals:
  - minimal dependencies (works without pandas; uses pandas if available)
  - robust column-name inference for chr/pos/pvalue/lod
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple


try:
    import pandas as _pd  # type: ignore

    _HAVE_PANDAS = True
except Exception:
    _pd = None
    _HAVE_PANDAS = False


CHR_CANDIDATES = [
    "chr",
    "chrom",
    "chromosome",
    "seqid",
    "scaffold",
    "contig",
    "chromosome_id",
    "peak_chr",
    "CHROM",
]
POS_CANDIDATES = [
    "pos",
    "position",
    "bp",
    "bp_pos",
    "snp_pos",
    "physical_pos",
    "peak_pos",
    "start",
    "POS",
]
PVAL_CANDIDATES = [
    "p",
    "pvalue",
    "p_value",
    "p-val",
    "p.val",
    "pval",
    "P",
    "PVAL",
    "p_wald",
    "p_wald_dom",
    "p_wald_add",
    "p_wald_gen",
]
LOD_CANDIDATES = [
    "lod",
    "LOD",
    "lod_score",
    "lodscore",
    "lod.max",
    "lod_max",
    "maxlod",
    "max_lod",
]


def _strip_chr_prefix(x: str) -> str:
    s = str(x).strip()
    if not s:
        return s
    if s.lower().startswith("chr"):
        s = s[3:]
    return s.strip()


def _sniff_delimiter(path: Path) -> str:
    # Heuristic: .csv -> ',', else try sniff from first line
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
    # map lower->original (first occurrence wins)
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
    # fuzzy contains (e.g., "chromosome" in "Chromosome") already covered by lower map
    return None


def _safe_float(x: Any) -> Optional[float]:
    if x is None:
        return None
    if isinstance(x, (int, float)):
        if isinstance(x, float) and (math.isnan(x) or math.isinf(x)):
            return None
        return float(x)
    s = str(x).strip()
    if not s or s.lower() in {"na", "nan", "none", "."}:
        return None
    try:
        v = float(s)
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


def _dataset_label(path: Path) -> str:
    # keep basename without double extensions
    name = path.name
    for suf in [".tsv", ".csv", ".txt"]:
        if name.lower().endswith(suf):
            name = name[: -len(suf)]
            break
    return name


@dataclass
class Peak:
    source_type: str  # gwas/qtl
    dataset: str
    raw_file: str
    chr: str
    pos: int
    start: int
    end: int
    pvalue: Optional[float] = None
    lod: Optional[float] = None


def _read_table_rows(path: Path) -> Tuple[List[str], List[Dict[str, Any]]]:
    """Read tsv/csv into list of dict rows. Minimal deps."""
    delim = _sniff_delimiter(path)
    with path.open("r", encoding="utf-8", errors="ignore", newline="") as f:
        reader = csv.DictReader(f, delimiter=delim)
        cols = reader.fieldnames or []
        rows = [r for r in reader]
    return cols, rows


def _read_table_df(path: Path):
    delim = _sniff_delimiter(path)
    if delim == "\t":
        return _pd.read_csv(path, sep="\t")
    return _pd.read_csv(path)


def _iter_gwas_peaks(
    path: Path,
    mode: str,
    top_n: int,
    p_threshold: float,
    window_bp: int,
    normalize_chr: bool,
) -> List[Peak]:
    ds = _dataset_label(path)
    raw = str(path)

    if _HAVE_PANDAS:
        df = _read_table_df(path)
        cols = list(df.columns)
        chr_col = _pick_col(cols, CHR_CANDIDATES)
        pos_col = _pick_col(cols, POS_CANDIDATES)
        p_col = _pick_col(cols, PVAL_CANDIDATES)
        if not chr_col or not pos_col or not p_col:
            return []
        sub = df[[chr_col, pos_col, p_col]].copy()
        sub.columns = ["chr", "pos", "pvalue"]
        # numeric
        sub["pos"] = _pd.to_numeric(sub["pos"], errors="coerce")
        sub["pvalue"] = _pd.to_numeric(sub["pvalue"], errors="coerce")
        sub = sub.dropna(subset=["chr", "pos", "pvalue"])
        if normalize_chr:
            sub["chr"] = sub["chr"].astype(str).map(_strip_chr_prefix)
        sub["pos"] = sub["pos"].astype(int)
        # filter
        if mode == "threshold":
            sub = sub[sub["pvalue"] <= float(p_threshold)]
        else:
            sub = sub.sort_values("pvalue", ascending=True).head(int(top_n))
        peaks: List[Peak] = []
        for _, r in sub.iterrows():
            chr_ = str(r["chr"]).strip()
            pos = int(r["pos"])
            pval = float(r["pvalue"])
            peaks.append(
                Peak(
                    source_type="gwas",
                    dataset=ds,
                    raw_file=raw,
                    chr=chr_,
                    pos=pos,
                    start=pos - int(window_bp),
                    end=pos + int(window_bp),
                    pvalue=pval,
                    lod=None,
                )
            )
        return peaks

    # no pandas
    cols, rows = _read_table_rows(path)
    chr_col = _pick_col(cols, CHR_CANDIDATES)
    pos_col = _pick_col(cols, POS_CANDIDATES)
    p_col = _pick_col(cols, PVAL_CANDIDATES)
    if not chr_col or not pos_col or not p_col:
        return []

    tmp: List[Tuple[float, Peak]] = []
    for r in rows:
        chr_raw = r.get(chr_col)
        pos_raw = r.get(pos_col)
        p_raw = r.get(p_col)
        pos = _safe_int(pos_raw)
        pval = _safe_float(p_raw)
        if pos is None or pval is None:
            continue
        chr_ = _strip_chr_prefix(chr_raw) if normalize_chr else str(chr_raw).strip()
        if not chr_:
            continue
        pk = Peak(
            source_type="gwas",
            dataset=ds,
            raw_file=raw,
            chr=chr_,
            pos=pos,
            start=pos - int(window_bp),
            end=pos + int(window_bp),
            pvalue=pval,
            lod=None,
        )
        tmp.append((pval, pk))

    if mode == "threshold":
        return [pk for p, pk in tmp if p <= float(p_threshold)]
    tmp.sort(key=lambda t: t[0])
    return [pk for _, pk in tmp[: int(top_n)]]


def _iter_qtl_peaks(
    path: Path,
    mode: str,
    top_n: int,
    lod_threshold: float,
    window_bp: int,
    normalize_chr: bool,
) -> List[Peak]:
    ds = _dataset_label(path)
    raw = str(path)

    if _HAVE_PANDAS:
        df = _read_table_df(path)
        cols = list(df.columns)
        chr_col = _pick_col(cols, CHR_CANDIDATES)
        pos_col = _pick_col(cols, POS_CANDIDATES)
        lod_col = _pick_col(cols, LOD_CANDIDATES)
        if not chr_col or not pos_col or not lod_col:
            return []
        sub = df[[chr_col, pos_col, lod_col]].copy()
        sub.columns = ["chr", "pos", "lod"]
        sub["pos"] = _pd.to_numeric(sub["pos"], errors="coerce")
        sub["lod"] = _pd.to_numeric(sub["lod"], errors="coerce")
        sub = sub.dropna(subset=["chr", "pos", "lod"])
        if normalize_chr:
            sub["chr"] = sub["chr"].astype(str).map(_strip_chr_prefix)
        sub["pos"] = sub["pos"].astype(int)
        if mode == "threshold":
            sub = sub[sub["lod"] >= float(lod_threshold)]
        else:
            sub = sub.sort_values("lod", ascending=False).head(int(top_n))
        peaks: List[Peak] = []
        for _, r in sub.iterrows():
            chr_ = str(r["chr"]).strip()
            pos = int(r["pos"])
            lod = float(r["lod"])
            peaks.append(
                Peak(
                    source_type="qtl",
                    dataset=ds,
                    raw_file=raw,
                    chr=chr_,
                    pos=pos,
                    start=pos - int(window_bp),
                    end=pos + int(window_bp),
                    pvalue=None,
                    lod=lod,
                )
            )
        return peaks

    cols, rows = _read_table_rows(path)
    chr_col = _pick_col(cols, CHR_CANDIDATES)
    pos_col = _pick_col(cols, POS_CANDIDATES)
    lod_col = _pick_col(cols, LOD_CANDIDATES)
    if not chr_col or not pos_col or not lod_col:
        return []

    tmp: List[Tuple[float, Peak]] = []
    for r in rows:
        chr_raw = r.get(chr_col)
        pos_raw = r.get(pos_col)
        lod_raw = r.get(lod_col)
        pos = _safe_int(pos_raw)
        lod = _safe_float(lod_raw)
        if pos is None or lod is None:
            continue
        chr_ = _strip_chr_prefix(chr_raw) if normalize_chr else str(chr_raw).strip()
        if not chr_:
            continue
        pk = Peak(
            source_type="qtl",
            dataset=ds,
            raw_file=raw,
            chr=chr_,
            pos=pos,
            start=pos - int(window_bp),
            end=pos + int(window_bp),
            pvalue=None,
            lod=lod,
        )
        tmp.append((lod, pk))

    if mode == "threshold":
        return [pk for lod, pk in tmp if lod >= float(lod_threshold)]
    tmp.sort(key=lambda t: t[0], reverse=True)
    return [pk for _, pk in tmp[: int(top_n)]]


def _cluster_intervals(peaks: List[Peak]) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]], List[Dict[str, Any]]]:
    """Return loci_rows, locus_peaks_rows, merged_rows."""
    # group by chr
    by_chr: Dict[str, List[Peak]] = {}
    for p in peaks:
        by_chr.setdefault(p.chr, []).append(p)

    loci_rows: List[Dict[str, Any]] = []
    locus_peaks_rows: List[Dict[str, Any]] = []
    merged_rows: List[Dict[str, Any]] = []

    locus_idx = 0
    for chr_ in sorted(by_chr.keys(), key=lambda x: (str(x))):
        arr = by_chr[chr_]
        arr.sort(key=lambda p: (p.start, p.end, p.pos))

        cur: List[Peak] = []
        cur_start: Optional[int] = None
        cur_end: Optional[int] = None

        def flush_cluster(cluster: List[Peak], start: int, end: int):
            nonlocal locus_idx
            locus_idx += 1
            locus_id = f"L{locus_idx}"

            # representative peak
            gwas = [p for p in cluster if p.source_type == "gwas" and p.pvalue is not None]
            qtl = [p for p in cluster if p.source_type == "qtl" and p.lod is not None]
            rep: Peak
            rep_source: str
            rep_p: Optional[float] = None
            rep_lod: Optional[float] = None
            if gwas:
                rep = min(gwas, key=lambda p: float(p.pvalue) if p.pvalue is not None else float("inf"))
                rep_source = "gwas"
                rep_p = rep.pvalue
            else:
                rep = max(qtl, key=lambda p: float(p.lod) if p.lod is not None else float("-inf"))
                rep_source = "qtl"
                rep_lod = rep.lod

            # locus summary
            loci_rows.append(
                {
                    "locus_id": locus_id,
                    "chr": chr_,
                    "start": int(start),
                    "end": int(end),
                    "rep_pos": int(rep.pos),
                    "rep_source": rep_source,
                    "rep_pvalue": rep_p if rep_p is not None else "",
                    "rep_lod": rep_lod if rep_lod is not None else "",
                    "n_peaks": len(cluster),
                    "n_gwas": len(gwas),
                    "n_qtl": len(qtl),
                }
            )

            # per-peak rows
            for p in cluster:
                locus_peaks_rows.append(
                    {
                        "locus_id": locus_id,
                        "source_type": p.source_type,
                        "dataset_label": p.dataset,
                        "chr": p.chr,
                        "pos": p.pos,
                        "start": p.start,
                        "end": p.end,
                        "pvalue": p.pvalue if p.pvalue is not None else "",
                        "lod": p.lod if p.lod is not None else "",
                        "raw_file": p.raw_file,
                    }
                )

            # merged representative peak row
            merged_rows.append(
                {
                    "locus_id": locus_id,
                    "chr": chr_,
                    "pos": int(rep.pos),
                    "start": int(start),
                    "end": int(end),
                    "rep_source": rep_source,
                    "rep_pvalue": rep_p if rep_p is not None else "",
                    "rep_lod": rep_lod if rep_lod is not None else "",
                }
            )

        for p in arr:
            if cur_start is None:
                cur = [p]
                cur_start = p.start
                cur_end = p.end
                continue
            assert cur_end is not None
            if p.start <= cur_end:
                cur.append(p)
                if p.end > cur_end:
                    cur_end = p.end
            else:
                flush_cluster(cur, int(cur_start), int(cur_end))
                cur = [p]
                cur_start = p.start
                cur_end = p.end

        if cur and cur_start is not None and cur_end is not None:
            flush_cluster(cur, int(cur_start), int(cur_end))

    return loci_rows, locus_peaks_rows, merged_rows


def _write_tsv(path: Path, rows: List[Dict[str, Any]], columns: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=columns, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in columns})


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--params", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args(argv)

    params_path = Path(args.params)
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    params = json.loads(params_path.read_text(encoding="utf-8"))

    gwas_files = [Path(p) for p in (params.get("gwas_results") or []) if str(p).strip()]
    qtl_files = [Path(p) for p in (params.get("qtl_peaks") or []) if str(p).strip()]

    gwas_mode = str(params.get("mode_gwas") or "top_n").strip().lower()
    qtl_mode = str(params.get("mode_qtl") or "top_n").strip().lower()
    gwas_top_n = int(params.get("gwas_top_n") or 200)
    qtl_top_n = int(params.get("qtl_top_n") or 50)
    gwas_p_threshold = float(params.get("gwas_p_threshold") or 1e-6)
    qtl_lod_threshold = float(params.get("qtl_lod_threshold") or 3.0)
    window_bp = int(params.get("window_bp") or 100000)
    normalize_chr = bool(params.get("chr_normalize", True))

    peaks: List[Peak] = []
    for p in gwas_files:
        if not p.exists():
            continue
        peaks.extend(
            _iter_gwas_peaks(
                p,
                mode=gwas_mode,
                top_n=gwas_top_n,
                p_threshold=gwas_p_threshold,
                window_bp=window_bp,
                normalize_chr=normalize_chr,
            )
        )

    for p in qtl_files:
        if not p.exists():
            continue
        peaks.extend(
            _iter_qtl_peaks(
                p,
                mode=qtl_mode,
                top_n=qtl_top_n,
                lod_threshold=qtl_lod_threshold,
                window_bp=window_bp,
                normalize_chr=normalize_chr,
            )
        )

    # nothing to do
    if not peaks:
        # still write empty outputs (helps GUI)
        _write_tsv(out_dir / "loci.tsv", [], ["locus_id", "chr", "start", "end", "rep_pos", "rep_source", "rep_pvalue", "rep_lod", "n_peaks", "n_gwas", "n_qtl"])
        _write_tsv(out_dir / "locus_peaks.tsv", [], ["locus_id", "source_type", "dataset_label", "chr", "pos", "start", "end", "pvalue", "lod", "raw_file"])
        _write_tsv(out_dir / "peaks_merged.tsv", [], ["locus_id", "chr", "pos", "start", "end", "rep_source", "rep_pvalue", "rep_lod"])
        (out_dir / "meta.json").write_text(json.dumps({"n_input_peaks": 0, "n_loci": 0}, indent=2), encoding="utf-8")
        return 0

    loci_rows, locus_peaks_rows, merged_rows = _cluster_intervals(peaks)

    _write_tsv(
        out_dir / "loci.tsv",
        loci_rows,
        ["locus_id", "chr", "start", "end", "rep_pos", "rep_source", "rep_pvalue", "rep_lod", "n_peaks", "n_gwas", "n_qtl"],
    )
    _write_tsv(
        out_dir / "locus_peaks.tsv",
        locus_peaks_rows,
        ["locus_id", "source_type", "dataset_label", "chr", "pos", "start", "end", "pvalue", "lod", "raw_file"],
    )
    _write_tsv(
        out_dir / "peaks_merged.tsv",
        merged_rows,
        ["locus_id", "chr", "pos", "start", "end", "rep_source", "rep_pvalue", "rep_lod"],
    )

    (out_dir / "meta.json").write_text(
        json.dumps(
            {
                "n_input_peaks": len(peaks),
                "n_loci": len(loci_rows),
                "have_pandas": _HAVE_PANDAS,
                "window_bp": window_bp,
                "gwas_mode": gwas_mode,
                "qtl_mode": qtl_mode,
            },
            indent=2,
        ),
        encoding="utf-8",
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
