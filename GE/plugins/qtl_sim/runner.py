#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Optional, Tuple, List, Dict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Optional: p-value (if scipy installed)
try:
    from scipy.stats import f as f_dist  # type: ignore
    HAS_SCIPY = True
except Exception:
    HAS_SCIPY = False


def log(msg: str, fp):
    fp.write(msg.rstrip() + "\n")
    fp.flush()


def read_tsv(path: str | Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str)


def to_float_series(s: pd.Series) -> pd.Series:
    return pd.to_numeric(s, errors="coerce")


def align_by_id(
    geno: pd.DataFrame, pheno: pd.DataFrame, cov: Optional[pd.DataFrame] = None
) -> Tuple[pd.DataFrame, pd.DataFrame, Optional[pd.DataFrame]]:
    geno["id"] = geno["id"].astype(str)
    pheno["id"] = pheno["id"].astype(str)

    common = set(geno["id"]).intersection(set(pheno["id"]))
    if cov is not None:
        cov["id"] = cov["id"].astype(str)
        common = common.intersection(set(cov["id"]))

    common = sorted(common)
    geno2 = geno.set_index("id").loc[common].reset_index()
    pheno2 = pheno.set_index("id").loc[common].reset_index()
    cov2 = None
    if cov is not None:
        cov2 = cov.set_index("id").loc[common].reset_index()
    return geno2, pheno2, cov2


def build_design_matrix(n: int, cov: Optional[pd.DataFrame]) -> np.ndarray:
    # Intercept + covariates (if any)
    X = np.ones((n, 1), dtype=float)
    if cov is None:
        return X

    cov_cols = [c for c in cov.columns if c != "id"]
    if not cov_cols:
        return X

    cov_num = cov[cov_cols].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    X = np.concatenate([X, cov_num], axis=1)
    return X


def fit_rss(y: np.ndarray, X: np.ndarray) -> Tuple[float, int]:
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    rss = float(np.sum(resid * resid))
    df_resid = max(0, X.shape[0] - X.shape[1])
    return rss, df_resid


def lod_from_rss(rss0: float, rss1: float, n: int) -> float:
    if rss1 <= 0 or rss0 <= 0:
        return np.nan
    if rss1 > rss0:
        return 0.0
    return (n / 2.0) * np.log10(rss0 / rss1)


def pvalue_f_test(rss0: float, rss1: float, df0: int, df1: int) -> float:
    if not HAS_SCIPY:
        return np.nan
    df_num = df0 - df1
    df_den = df1
    if df_num <= 0 or df_den <= 0:
        return np.nan
    if rss1 <= 0:
        return np.nan
    F = ((rss0 - rss1) / df_num) / (rss1 / df_den)
    if F < 0:
        return 1.0
    return float(f_dist.sf(F, df_num, df_den))


def prepare_map(markers: List[str], map_path: Optional[str]) -> Tuple[pd.DataFrame, Dict[str, str], Dict[str, float]]:
    if map_path:
        mdf = read_tsv(map_path)
        need = {"marker", "chr", "pos"}
        if not need.issubset(set(mdf.columns)):
            raise SystemExit("marker_map_tsv must have columns: marker, chr, pos")
        mdf = mdf[mdf["marker"].isin(markers)].copy()
        mdf["marker"] = mdf["marker"].astype(str)
        mdf["chr"] = mdf["chr"].astype(str)
        mdf["pos"] = pd.to_numeric(mdf["pos"], errors="coerce")
    else:
        mdf = pd.DataFrame({
            "marker": markers,
            "chr": ["1"] * len(markers),
            "pos": np.arange(len(markers), dtype=float)
        })

    chr_map = dict(zip(mdf["marker"], mdf["chr"]))
    pos_map = dict(zip(mdf["marker"], mdf["pos"].astype(float)))
    return mdf, chr_map, pos_map


def scan_markers(
    *,
    geno: pd.DataFrame,
    y_all: np.ndarray,
    X_base: np.ndarray,
    markers: List[str],
    chr_map: Dict[str, str],
    pos_map: Dict[str, float],
    cofactors: Optional[List[str]] = None,
    window: float = 0.0,
) -> pd.DataFrame:
    """
    Per-marker regression:
      Null: y ~ X_base + cofactors
      Alt : y ~ X_base + cofactors + g_m
    window: exclude cofactors on same chr within +/- window (pos distance).
    """
    cofactors = cofactors or []

    # Preconvert all cofactor genotype columns to numeric arrays (lazy)
    g_cache: Dict[str, np.ndarray] = {}

    def get_g(marker: str) -> np.ndarray:
        if marker not in g_cache:
            g_cache[marker] = pd.to_numeric(geno[marker], errors="coerce").to_numpy(dtype=float)
        return g_cache[marker]

    rows = []
    n_total = len(geno)

    for m in markers:
        g_m = get_g(m)

        # choose cofactors for this marker (window exclusion)
        cof_use: List[str] = []
        if cofactors:
            chr_m = chr_map.get(m, "NA")
            pos_m = pos_map.get(m, np.nan)
            for c in cofactors:
                if c == m:
                    continue
                # exclude if same chr and within window (if pos known)
                if window > 0 and chr_map.get(c, "NA") == chr_m:
                    pos_c = pos_map.get(c, np.nan)
                    if np.isfinite(pos_m) and np.isfinite(pos_c):
                        if abs(pos_c - pos_m) <= window:
                            continue
                cof_use.append(c)

        # build X0 and X1 with per-sample masking
        # base finite mask
        ok = np.isfinite(y_all) & np.isfinite(g_m)
        ok &= np.all(np.isfinite(X_base), axis=1)

        # add cofactors (and also require they are finite)
        cof_cols = []
        for c in cof_use:
            if c not in geno.columns:
                continue
            g_c = get_g(c)
            ok &= np.isfinite(g_c)
            cof_cols.append(g_c)

        n = int(np.sum(ok))
        if n < 10:
            rows.append((m, chr_map.get(m, "NA"), pos_map.get(m, np.nan), np.nan, np.nan, n, len(cof_cols)))
            continue

        y = y_all[ok]
        X0 = X_base[ok, :]
        if cof_cols:
            C = np.stack([g[ok] for g in cof_cols], axis=1)
            X0 = np.concatenate([X0, C], axis=1)

        X1 = np.concatenate([X0, g_m[ok][:, None]], axis=1)

        rss0, df0 = fit_rss(y, X0)
        rss1, df1 = fit_rss(y, X1)

        lod = lod_from_rss(rss0, rss1, n)
        pval = pvalue_f_test(rss0, rss1, df0, df1)

        rows.append((m, chr_map.get(m, "NA"), pos_map.get(m, np.nan), lod, pval, n, len(cof_cols)))

    scan = pd.DataFrame(rows, columns=["marker", "chr", "pos", "lod", "pvalue", "n", "n_cofactors"])
    scan["pos"] = pd.to_numeric(scan["pos"], errors="coerce")
    scan = scan.sort_values(["chr", "pos"], kind="mergesort")
    return scan


def pick_cofactors_from_sim(scan: pd.DataFrame, k: int) -> List[str]:
    if k <= 0:
        return []
    s = scan.dropna(subset=["lod"]).copy()
    s = s[np.isfinite(s["lod"])]
    s = s.sort_values("lod", ascending=False)
    return s["marker"].astype(str).head(k).tolist()


def make_peaks(scan: pd.DataFrame) -> pd.DataFrame:
    peaks = []
    scan2 = scan.dropna(subset=["lod"]).copy()
    scan2 = scan2[np.isfinite(scan2["lod"])]

    for c, sub in scan2.groupby("chr", sort=False):
        sub2 = sub.sort_values("lod", ascending=False)
        if len(sub2) > 0:
            peaks.append(sub2.iloc[0])

    peaks_df = pd.DataFrame(peaks) if peaks else pd.DataFrame(columns=scan.columns)

    top1 = scan2.sort_values("lod", ascending=False).head(1)
    if len(top1) == 1:
        top1 = top1.assign(chr="ALL_TOP")
        peaks_df = pd.concat([peaks_df, top1], ignore_index=True)

    return peaks_df


def plot_lod(scan: pd.DataFrame, trait: str, png_path: Path, title_prefix: str):
    plot_df = scan.dropna(subset=["lod", "pos"]).copy()
    if len(plot_df) == 0:
        plt.figure()
        plt.title(f"{title_prefix} (trait={trait})")
        plt.savefig(png_path, dpi=150)
        plt.close()
        return

    def chr_key(x: str):
        try:
            return (0, int(float(x)))
        except Exception:
            return (1, x)

    chrs = sorted(plot_df["chr"].astype(str).unique().tolist(), key=chr_key)

    offset = 0.0
    xs = []
    ys = []
    ticks = []
    ticklabs = []

    for c in chrs:
        sub = plot_df[plot_df["chr"].astype(str) == str(c)].sort_values("pos")
        if len(sub) == 0:
            continue
        x = sub["pos"].to_numpy(dtype=float) + offset
        y = sub["lod"].to_numpy(dtype=float)
        xs.append(x)
        ys.append(y)
        ticks.append((x.min() + x.max()) / 2.0)
        ticklabs.append(str(c))
        offset = x.max() + 1.0

    plt.figure()
    x_all = np.concatenate(xs)
    y_all = np.concatenate(ys)
    plt.plot(x_all, y_all)
    plt.xticks(ticks, ticklabs, rotation=0)
    plt.xlabel("Chromosome")
    plt.ylabel("LOD")
    plt.title(f"{title_prefix} (trait={trait})")
    plt.tight_layout()
    plt.savefig(png_path, dpi=150)
    plt.close()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--params", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    out_dir = Path(args.out).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    plot_dir = out_dir / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)

    log_path = out_dir / "run.log"
    with open(log_path, "w", encoding="utf-8") as fp:
        log("[qtl_sim/cim] start", fp)
        log(f"[qtl_sim/cim] params={args.params}", fp)
        log(f"[qtl_sim/cim] out={out_dir}", fp)

        params = json.loads(Path(args.params).read_text(encoding="utf-8"))
        geno_path = params.get("genotype_tsv")
        pheno_path = params.get("phenotype_tsv")
        map_path = params.get("marker_map_tsv", None)
        cov_path = params.get("covariates_tsv", None)
        trait = params.get("trait")

        mode = str(params.get("mode", "sim")).lower()  # sim or cim
        cofactors_k = int(params.get("cofactors_k", 5))
        window = float(params.get("window", 10.0))

        if mode not in {"sim", "cim"}:
            raise SystemExit("mode must be 'sim' or 'cim'")

        if not geno_path or not pheno_path or not trait:
            raise SystemExit("params must contain genotype_tsv, phenotype_tsv, trait")

        log(f"[qtl_sim/cim] genotype_tsv={geno_path}", fp)
        log(f"[qtl_sim/cim] phenotype_tsv={pheno_path}", fp)
        log(f"[qtl_sim/cim] marker_map_tsv={map_path}", fp)
        log(f"[qtl_sim/cim] covariates_tsv={cov_path}", fp)
        log(f"[qtl_sim/cim] trait={trait}", fp)
        log(f"[qtl_sim/cim] mode={mode}", fp)
        log(f"[qtl_sim/cim] cofactors_k={cofactors_k}", fp)
        log(f"[qtl_sim/cim] window={window}", fp)

        geno = read_tsv(geno_path)
        pheno = read_tsv(pheno_path)

        if "id" not in geno.columns:
            raise SystemExit("genotype_tsv must have column: id")
        if "id" not in pheno.columns:
            raise SystemExit("phenotype_tsv must have column: id")
        if trait not in pheno.columns:
            raise SystemExit(f"trait not found in phenotype_tsv: {trait}")

        cov = None
        if cov_path:
            cov = read_tsv(cov_path)
            if "id" not in cov.columns:
                raise SystemExit("covariates_tsv must have column: id")

        geno, pheno, cov = align_by_id(geno, pheno, cov)
        n_total = len(geno)
        log(f"[qtl_sim/cim] N aligned samples = {n_total}", fp)
        if n_total < 5:
            raise SystemExit("Too few overlapping samples after aligning by id")

        y_all = to_float_series(pheno[trait]).to_numpy(dtype=float)

        markers = [c for c in geno.columns if c != "id"]
        if not markers:
            raise SystemExit("No marker columns found in genotype_tsv (besides id)")
        log(f"[qtl_sim/cim] markers={len(markers)}", fp)

        mdf, chr_map, pos_map = prepare_map(markers, map_path)

        X_base = build_design_matrix(n_total, cov)

        # 1) SIM always available (and used for cofactors selection in CIM)
        log("[qtl_sim/cim] running SIM...", fp)
        sim_scan = scan_markers(
            geno=geno, y_all=y_all, X_base=X_base,
            markers=markers, chr_map=chr_map, pos_map=pos_map,
            cofactors=None, window=0.0
        )

        # If mode=sim, output SIM
        if mode == "sim":
            scan = sim_scan
            peaks_df = make_peaks(scan)
            scan.to_csv(out_dir / "scan.tsv", sep="\t", index=False)
            peaks_df.to_csv(out_dir / "peaks.tsv", sep="\t", index=False)

            plot_lod(scan, trait, plot_dir / "lod.png", "QTL scan (SIM)")

            # alias
            peaks_df.to_csv(out_dir / "results.tsv", sep="\t", index=False)

            log("[qtl_sim/cim] wrote scan.tsv / peaks.tsv / plots/lod.png", fp)
            log("[qtl_sim/cim] done (SIM)", fp)
            return

        # 2) CIM
        log("[qtl_sim/cim] selecting cofactors from SIM...", fp)
        cofactors = pick_cofactors_from_sim(sim_scan, cofactors_k)
        log(f"[qtl_sim/cim] cofactors={cofactors}", fp)

        # save sim_scan for reference
        sim_scan.to_csv(out_dir / "sim_scan.tsv", sep="\t", index=False)

        log("[qtl_sim/cim] running CIM...", fp)
        cim_scan = scan_markers(
            geno=geno, y_all=y_all, X_base=X_base,
            markers=markers, chr_map=chr_map, pos_map=pos_map,
            cofactors=cofactors, window=window
        )

        peaks_df = make_peaks(cim_scan)

        cim_scan.to_csv(out_dir / "scan.tsv", sep="\t", index=False)
        peaks_df.to_csv(out_dir / "peaks.tsv", sep="\t", index=False)

        plot_lod(cim_scan, trait, plot_dir / "lod.png", "QTL scan (CIM)")

        # alias
        peaks_df.to_csv(out_dir / "results.tsv", sep="\t", index=False)

        log("[qtl_sim/cim] wrote scan.tsv / peaks.tsv / plots/lod.png (+ sim_scan.tsv)", fp)
        log("[qtl_sim/cim] done (CIM)", fp)


if __name__ == "__main__":
    main()

