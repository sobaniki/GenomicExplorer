#!/usr/bin/env python3

import argparse
import json
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

try:
    import plotly.graph_objects as go  # type: ignore
    import plotly.io as pio  # type: ignore
    from plotly.subplots import make_subplots  # type: ignore
except Exception:
    go = None
    pio = None
    make_subplots = None


def _sort_chr_labels(chrs: List[str]) -> List[str]:
    def _key(c: str):
        m = re.search(r"\d+", str(c))
        if m:
            return (0, int(m.group(0)), str(c))
        return (1, 10**9, str(c))

    return sorted([str(x) for x in chrs], key=_key)


def _standardize_lod_df(df: pd.DataFrame) -> pd.DataFrame:
    d = df.copy()
    d.columns = [str(c).strip() for c in d.columns]

    #if "LOD" in d.columns and "lod" not in d.columns:
    #    d.rename(columns={"LOD": "lod"}, inplace=True)
    
    if "lod" not in d.columns:
        for cand in ["LOD", "main"]:
            if cand in d.columns:
                d.rename(columns={cand: "pos"}, inplace=True)
                break

    if "pos" not in d.columns:
        for cand in ["cM", "cm", "position", "Position", "pos_cM", "pos_cm"]:
            if cand in d.columns:
                d.rename(columns={cand: "pos"}, inplace=True)
                break
    if "chr" not in d.columns:
        for cand in ["chrom", "chromosome", "Chr", "CHR"]:
            if cand in d.columns:
                d.rename(columns={cand: "chr"}, inplace=True)
                break
    if "lod" not in d.columns:
        for c in d.columns:
            if c.lower() in ("marker", "chr", "pos"):
                continue
            if pd.api.types.is_numeric_dtype(d[c]):
                d.rename(columns={c: "lod"}, inplace=True)
                break

    need = {"chr", "pos", "lod"}
    if not need.issubset(set(d.columns)):
        raise ValueError(f"LOD profile missing required columns {need}. Found: {list(d.columns)}")

    d = d[["chr", "pos", "lod"]].copy()
    d["chr"] = d["chr"].astype(str)
    d["pos"] = pd.to_numeric(d["pos"], errors="coerce")
    d["lod"] = pd.to_numeric(d["lod"], errors="coerce")
    d = d.dropna(subset=["pos", "lod"])
    return d


def _genome_offsets(lod_df: pd.DataFrame) -> Tuple[pd.DataFrame, Dict[str, float], List[str]]:
    d = lod_df.copy()
    chrs = _sort_chr_labels(d["chr"].dropna().astype(str).unique().tolist())
    offsets: Dict[str, float] = {}
    cur = 0.0
    gap = 5.0
    for c in chrs:
        sub = d[d["chr"].astype(str) == c]
        if sub.empty:
            offsets[c] = cur
            continue
        mn = float(np.nanmin(sub["pos"].values))
        mx = float(np.nanmax(sub["pos"].values))
        offsets[c] = cur - mn
        cur += (mx - mn) + gap
    d["x"] = d.apply(lambda r: float(r["pos"]) + offsets.get(str(r["chr"]), 0.0), axis=1)
    return d, offsets, chrs


def _infer_ci_cols(peaks_df: pd.DataFrame) -> Tuple[Optional[str], Optional[str]]:
    cols = list(peaks_df.columns)
    for l, r in [("ci.lo", "ci.hi"), ("ci_lo", "ci_hi"), ("ci_l", "ci_r"), ("ci_left", "ci_right"), ("l", "r")]:
        if l in cols and r in cols:
            return l, r
    for l in cols:
        if re.search(r"ci.*(lo|left|lower)", l, flags=re.I):
            for r in cols:
                if re.search(r"ci.*(hi|right|upper)", r, flags=re.I):
                    return l, r
    return None, None


def _read_peaks(peaks_tsv: Path) -> Optional[pd.DataFrame]:
    if not peaks_tsv or not peaks_tsv.exists():
        return None
    try:
        df = pd.read_csv(peaks_tsv, sep="\t")
    except Exception:
        return None

    if "chr" not in df.columns:
        for cand in ["Chr", "chrom", "chromosome", "CHR"]:
            if cand in df.columns:
                df.rename(columns={cand: "chr"}, inplace=True)
                break
    if "pos" not in df.columns:
        for cand in ["cM", "cm", "position", "pos_cM", "pos_cm"]:
            if cand in df.columns:
                df.rename(columns={cand: "pos"}, inplace=True)
                break
    if "lod" not in df.columns:
        if "LOD" in df.columns:
            df.rename(columns={"LOD": "lod"}, inplace=True)

    if not {"chr", "pos", "lod"}.issubset(set(df.columns)):
        return None

    df["chr"] = df["chr"].astype(str)
    df["pos"] = pd.to_numeric(df["pos"], errors="coerce")
    df["lod"] = pd.to_numeric(df["lod"], errors="coerce")
    df = df.dropna(subset=["chr", "pos", "lod"])
    return df


def _read_thresholds(perm_thresholds_tsv: Optional[Path]) -> List[Tuple[float, str]]:
    out: List[Tuple[float, str]] = []
    if not perm_thresholds_tsv or not perm_thresholds_tsv.exists():
        return out
    try:
        df = pd.read_csv(perm_thresholds_tsv, sep="\t")
    except Exception:
        return out
    if "threshold" not in df.columns:
        return out
    for _, row in df.iterrows():
        try:
            t = float(row.get("threshold"))
        except Exception:
            continue
        if not np.isfinite(t):
            continue
        try:
            a = float(row.get("alpha", np.nan))
        except Exception:
            a = np.nan
        lab = f"thr (alpha={a:g})" if np.isfinite(a) else "threshold"
        out.append((t, lab))
    return out


def build_qtl_plotly_html(
    lod_profile_tsv: Path,
    peaks_tsv: Optional[Path],
    perm_thresholds_tsv: Optional[Path],
    out_html: Path,
    *,
    title: str,
    pos_unit: str = "cM",
    zoom_window: float = 20.0,
) -> bool:
    if go is None or pio is None or make_subplots is None:
        return False

    lod_raw = pd.read_csv(lod_profile_tsv, sep="\t")
    lod_df = _standardize_lod_df(lod_raw)
    d, offsets, chrs = _genome_offsets(lod_df)

    peaks_df = _read_peaks(peaks_tsv) if peaks_tsv else None
    top_peak = None
    ci_l = ci_r = None
    if peaks_df is not None and not peaks_df.empty:
        top_peak = peaks_df.sort_values("lod", ascending=False).iloc[0]
        ci_l, ci_r = _infer_ci_cols(peaks_df)

    thr_lines = _read_thresholds(perm_thresholds_tsv)

    fig = make_subplots(
        rows=2,
        cols=1,
        shared_xaxes=False,
        vertical_spacing=0.12,
        row_heights=[0.62, 0.38],
        subplot_titles=("Genome-wide scan", "Peak region"),
    )

    fig.add_trace(
        go.Scatter(
            x=d["x"],
            y=d["lod"],
            mode="lines",
            name="LOD",
            hovertemplate=(
                "chr=%{customdata[0]}<br>pos=%{customdata[1]:.3f} "
                + pos_unit
                + "<br>LOD=%{y:.3f}<extra></extra>"
            ),
            customdata=np.stack([d["chr"].values, d["pos"].values], axis=1),
        ),
        row=1,
        col=1,
    )

    # chromosome boundaries / alternating background
    for i, c in enumerate(chrs):
        sub = d[d["chr"].astype(str) == c]
        if sub.empty:
            continue
        xmin = float(np.nanmin(sub["x"].values))
        xmax = float(np.nanmax(sub["x"].values))
        if i % 2 == 0:
            fig.add_vrect(x0=xmin, x1=xmax, fillcolor="rgba(0,0,0,0.03)", line_width=0, row=1, col=1)
        fig.add_vline(x=xmax, line_width=1, line_dash="dot", line_color="rgba(0,0,0,0.25)", row=1, col=1)

    tick_x: List[float] = []
    tick_t: List[str] = []
    for c in chrs:
        sub = d[d["chr"].astype(str) == c]
        if sub.empty:
            continue
        xmin = float(np.nanmin(sub["x"].values))
        xmax = float(np.nanmax(sub["x"].values))
        tick_x.append((xmin + xmax) / 2.0)
        tick_t.append(str(c))

    fig.update_xaxes(row=1, col=1, tickmode="array", tickvals=tick_x, ticktext=tick_t, title_text="Chromosome")
    fig.update_yaxes(row=1, col=1, title_text="LOD")

    if peaks_df is not None and not peaks_df.empty:
        xpk, ypk, hcd = [], [], []
        for _, r in peaks_df.iterrows():
            c = str(r["chr"])
            pos = float(r["pos"])
            lodv = float(r["lod"])
            xpk.append(pos + offsets.get(c, 0.0))
            ypk.append(lodv)
            hcd.append([c, pos])
        fig.add_trace(
            go.Scatter(
                x=xpk,
                y=ypk,
                mode="markers",
                name="Peaks",
                marker=dict(size=8),
                hovertemplate=(
                    "PEAK<br>chr=%{customdata[0]}<br>pos=%{customdata[1]:.3f} "
                    + pos_unit
                    + "<br>LOD=%{y:.3f}<extra></extra>"
                ),
                customdata=np.array(hcd),
            ),
            row=1,
            col=1,
        )

    for t, lab in thr_lines:
        fig.add_hline(y=t, line_dash="dash", line_width=2, annotation_text=lab, annotation_position="top left", row=1, col=1)
        fig.add_hline(y=t, line_dash="dash", line_width=2, row=2, col=1)

    # zoom panel
    if top_peak is not None:
        c0 = str(top_peak["chr"])
        p0 = float(top_peak["pos"])

        zmin, zmax = (p0 - zoom_window / 2.0), (p0 + zoom_window / 2.0)
        if peaks_df is not None and ci_l and ci_r and ci_l in peaks_df.columns and ci_r in peaks_df.columns:
            try:
                l0 = float(top_peak[ci_l])
                r0 = float(top_peak[ci_r])
                if np.isfinite(l0) and np.isfinite(r0):
                    zmin, zmax = min(l0, r0), max(l0, r0)
            except Exception:
                pass

        sub = lod_df[(lod_df["chr"].astype(str) == c0) & (lod_df["pos"] >= zmin) & (lod_df["pos"] <= zmax)].copy()
        if sub.empty:
            sub = lod_df[lod_df["chr"].astype(str) == c0].copy()

        fig.add_trace(
            go.Scatter(
                x=sub["pos"],
                y=sub["lod"],
                mode="lines",
                name=f"LOD (chr {c0})",
                showlegend=False,
                hovertemplate=(
                    "chr=%{customdata[0]}<br>pos=%{x:.3f} " + pos_unit + "<br>LOD=%{y:.3f}<extra></extra>"
                ),
                customdata=np.stack([sub["chr"].values], axis=1),
            ),
            row=2,
            col=1,
        )

        fig.add_vrect(x0=zmin, x1=zmax, fillcolor="rgba(0,0,0,0.05)", line_width=0, row=2, col=1)
        fig.add_trace(
            go.Scatter(
                x=[p0],
                y=[float(top_peak["lod"])],
                mode="markers+text",
                text=[f"Peak: chr{c0}:{p0:.3f} {pos_unit}"],
                textposition="top center",
                marker=dict(size=10),
                showlegend=False,
            ),
            row=2,
            col=1,
        )

        fig.update_xaxes(row=2, col=1, title_text=f"Position ({pos_unit})")
        fig.update_yaxes(row=2, col=1, title_text="LOD")
    else:
        fig.update_xaxes(row=2, col=1, title_text=f"Position ({pos_unit})")
        fig.update_yaxes(row=2, col=1, title_text="LOD")

    fig.update_layout(
        title=title,
        hovermode="closest",
        margin=dict(l=60, r=25, t=70, b=60),
        height=900,
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="left", x=0.0),
    )

    out_html.parent.mkdir(parents=True, exist_ok=True)
    html = pio.to_html(fig, include_plotlyjs="inline", full_html=True)
    out_html.write_text(html, encoding="utf-8")
    return True


def _detect_files(qtl_out_dir: Path, trait: Optional[str]) -> Tuple[str, Path, Optional[Path], Optional[Path]]:
    """Return (trait_used, lod_profile, peaks, perm_thresholds)."""
    # r/qtl: single-file convention
    lod_single = qtl_out_dir / "lod_profile.tsv"
    if lod_single.exists():
        return (
            trait or "trait",
            lod_single,
            (qtl_out_dir / "peaks.tsv") if (qtl_out_dir / "peaks.tsv").exists() else None,
            (qtl_out_dir / "perm_thresholds.tsv") if (qtl_out_dir / "perm_thresholds.tsv").exists() else None,
        )

    # qtl2: per trait convention
    candidates = sorted(qtl_out_dir.glob("lod_profile_*.tsv"))
    if not candidates:
        raise FileNotFoundError("No lod_profile.tsv or lod_profile_*.tsv found")

    if trait:
        lod = qtl_out_dir / f"lod_profile_{trait}.tsv"
        if lod.exists():
            t = trait
        else:
            # try exact match by stripping spaces
            t = trait.strip()
            lod = qtl_out_dir / f"lod_profile_{t}.tsv"
            if not lod.exists():
                # fallback to first
                lod = candidates[0]
                t = lod.stem.replace("lod_profile_", "", 1)
    else:
        lod = candidates[0]
        t = lod.stem.replace("lod_profile_", "", 1)

    peaks = qtl_out_dir / f"peaks_{t}.tsv"
    perm = qtl_out_dir / f"perm_thresholds_{t}.tsv"
    return (t, lod, peaks if peaks.exists() else None, perm if perm.exists() else None)


def _try_export_images(fig, out_base: Path, formats: List[str], width: int, height: int) -> Dict[str, str]:
    """Export static images via kaleido if available. Returns mapping format->path."""
    if pio is None:
        return {}
    out: Dict[str, str] = {}
    for fmt in formats:
        fmt2 = fmt.lower().strip().lstrip(".")
        if fmt2 not in ("png", "svg", "pdf"):
            continue
        out_path = out_base.with_suffix("." + fmt2)
        try:
            pio.write_image(fig, str(out_path), format=fmt2, width=width, height=height, scale=1)
            out[fmt2] = str(out_path)
        except Exception:
            # kaleido not installed or export failed
            continue
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--params", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    params = json.loads(Path(args.params).read_text(encoding="utf-8"))
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    qtl_out_dir = Path(params.get("qtl_out_dir", "")).expanduser()
    if not qtl_out_dir.exists():
        raise FileNotFoundError(f"qtl_out_dir not found: {qtl_out_dir}")

    trait = params.get("trait")
    title = str(params.get("title") or "QTL LOD profile")
    pos_unit = str(params.get("pos_unit") or "cM")
    zoom_window = float(params.get("zoom_window") or 20.0)

    formats = params.get("export_formats") or []
    if isinstance(formats, str):
        formats = [s.strip() for s in formats.split(",") if s.strip()]
    width = int(params.get("export_width") or 1400)
    height = int(params.get("export_height") or 900)

    trait_used, lod_tsv, peaks_tsv, perm_tsv = _detect_files(qtl_out_dir, trait)

    plots_dir = out_dir / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)
    html_path = plots_dir / f"lod_plotly_{trait_used}.html"

    # Build figure twice: once for html via helper, and again for static exports if requested.
    ok = build_qtl_plotly_html(
        lod_profile_tsv=lod_tsv,
        peaks_tsv=peaks_tsv,
        perm_thresholds_tsv=perm_tsv,
        out_html=html_path,
        title=title,
        pos_unit=pos_unit,
        zoom_window=zoom_window,
    )
    if not ok:
        raise RuntimeError("Plotly not available (install plotly and optional kaleido)")

    artifacts = {
        "plot_html": str(html_path),
        "trait": trait_used,
        "lod_profile": str(lod_tsv),
        "peaks": str(peaks_tsv) if peaks_tsv else "",
        "perm_thresholds": str(perm_tsv) if perm_tsv else "",
    }

    # If exports requested, reconstruct figure to export images (no JS embedding).
    if formats and go is not None and pio is not None and make_subplots is not None:
        # Reuse internal builder by reading html? Instead, rebuild a fig quickly.
        lod_raw = pd.read_csv(lod_tsv, sep="\t")
        lod_df = _standardize_lod_df(lod_raw)
        d, offsets, chrs = _genome_offsets(lod_df)
        peaks_df = _read_peaks(peaks_tsv) if peaks_tsv else None
        top_peak = None
        ci_l = ci_r = None
        if peaks_df is not None and not peaks_df.empty:
            top_peak = peaks_df.sort_values("lod", ascending=False).iloc[0]
            ci_l, ci_r = _infer_ci_cols(peaks_df)
        thr_lines = _read_thresholds(perm_tsv)

        fig = make_subplots(rows=2, cols=1, shared_xaxes=False, vertical_spacing=0.12, row_heights=[0.62, 0.38], subplot_titles=("Genome-wide scan", "Peak region"))
        fig.add_trace(go.Scatter(x=d["x"], y=d["lod"], mode="lines", name="LOD", customdata=np.stack([d["chr"].values, d["pos"].values], axis=1)), row=1, col=1)
        for i, c in enumerate(chrs):
            sub = d[d["chr"].astype(str) == c]
            if sub.empty:
                continue
            xmin = float(np.nanmin(sub["x"].values))
            xmax = float(np.nanmax(sub["x"].values))
            if i % 2 == 0:
                fig.add_vrect(x0=xmin, x1=xmax, fillcolor="rgba(0,0,0,0.03)", line_width=0, row=1, col=1)
            fig.add_vline(x=xmax, line_width=1, line_dash="dot", line_color="rgba(0,0,0,0.25)", row=1, col=1)
        tick_x, tick_t = [], []
        for c in chrs:
            sub = d[d["chr"].astype(str) == c]
            if sub.empty:
                continue
            xmin = float(np.nanmin(sub["x"].values))
            xmax = float(np.nanmax(sub["x"].values))
            tick_x.append((xmin + xmax) / 2.0)
            tick_t.append(str(c))
        fig.update_xaxes(row=1, col=1, tickmode="array", tickvals=tick_x, ticktext=tick_t, title_text="Chromosome")
        fig.update_yaxes(row=1, col=1, title_text="LOD")
        if peaks_df is not None and not peaks_df.empty:
            xpk, ypk = [], []
            for _, r in peaks_df.iterrows():
                c = str(r["chr"])
                pos = float(r["pos"])
                lodv = float(r["lod"])
                xpk.append(pos + offsets.get(c, 0.0))
                ypk.append(lodv)
            fig.add_trace(go.Scatter(x=xpk, y=ypk, mode="markers", name="Peaks", marker=dict(size=8)), row=1, col=1)
        for t, _lab in thr_lines:
            fig.add_hline(y=t, line_dash="dash", line_width=2, row=1, col=1)
            fig.add_hline(y=t, line_dash="dash", line_width=2, row=2, col=1)
        if top_peak is not None:
            c0 = str(top_peak["chr"])
            p0 = float(top_peak["pos"])
            zmin, zmax = (p0 - zoom_window / 2.0), (p0 + zoom_window / 2.0)
            if peaks_df is not None and ci_l and ci_r and ci_l in peaks_df.columns and ci_r in peaks_df.columns:
                try:
                    l0 = float(top_peak[ci_l])
                    r0 = float(top_peak[ci_r])
                    if np.isfinite(l0) and np.isfinite(r0):
                        zmin, zmax = min(l0, r0), max(l0, r0)
                except Exception:
                    pass
            sub = lod_df[(lod_df["chr"].astype(str) == c0) & (lod_df["pos"] >= zmin) & (lod_df["pos"] <= zmax)].copy()
            if sub.empty:
                sub = lod_df[lod_df["chr"].astype(str) == c0].copy()
            fig.add_trace(go.Scatter(x=sub["pos"], y=sub["lod"], mode="lines", showlegend=False), row=2, col=1)
            fig.add_vrect(x0=zmin, x1=zmax, fillcolor="rgba(0,0,0,0.05)", line_width=0, row=2, col=1)
            fig.add_trace(go.Scatter(x=[p0], y=[float(top_peak["lod"])], mode="markers", showlegend=False), row=2, col=1)
        fig.update_xaxes(row=2, col=1, title_text=f"Position ({pos_unit})")
        fig.update_yaxes(row=2, col=1, title_text="LOD")
        fig.update_layout(title=title, margin=dict(l=60, r=25, t=70, b=60), height=900)

        exports = _try_export_images(fig, plots_dir / f"lod_plotly_{trait_used}", formats, width, height)
        artifacts["exports"] = exports

    (out_dir / "artifacts.json").write_text(json.dumps(artifacts, indent=2), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
