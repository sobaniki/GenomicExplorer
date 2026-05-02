#!/usr/bin/env python3

"""Preprocess: Map converter (bp↔cM estimation / format conversion)

In practice, teams rarely share the exact same marker set when building linkage maps.
Still, it is often useful to *approximate* genetic positions (cM) for a new marker set
using a reference map that contains both physical positions (bp) and genetic positions (cM).

This plugin implements a simple and robust "Marey map" style conversion:
  - Fit per-chromosome bp→cM relationships using reference (bp,cM) anchor points
  - Enforce monotonicity (optional; recommended)
  - Predict cM for target markers by interpolation/extrapolation

It also supports basic map format conversion.

Inputs (params.json)
--------------------
out_dir_override: str (optional)

input_kind: str  # auto|marker_map|plink
input_path: str  # marker_map.tsv/csv/xlsx OR plink prefix OR .bim

reference_map: str (optional)  # file containing chr,bp,cM (marker optional)
reference_sheet: str (optional; xlsx)
ref_auto_cols: bool (default True)
ref_chr_col: str (optional)
ref_bp_col: str (optional)
ref_cm_col: str (optional)
ref_marker_col: str (optional)

fit_method: str  # piecewise|linear_chr|linear_global
enforce_monotone: bool (default True)
default_cM_per_Mb: float (default 1.0)  # fallback slope

write_plink_bim_cm: bool (default False)  # if input_kind is plink

Outputs
-------
marker_map_bp.tsv     : marker, chr, pos (bp)
marker_map_cm.tsv     : marker, chr, pos (cM; when reference is provided)
marker_map.tsv        : marker, chr, pos (pos=cM when estimated, else bp)
marey_fit_report.tsv  : per-chrom fit summary
artifacts.json
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple


def _read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _write_json(path: Path, obj: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2), encoding="utf-8")


def _write_tsv(path: Path, rows: Iterable[Iterable], header: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(header)
        for r in rows:
            w.writerow(list(r))


def _as_bool(x, default: bool = False) -> bool:
    if x is None:
        return default
    if isinstance(x, bool):
        return x
    s = str(x).strip().lower()
    if s in {"1", "true", "t", "yes", "y"}:
        return True
    if s in {"0", "false", "f", "no", "n"}:
        return False
    return default


def _to_float(x) -> Optional[float]:
    if x is None:
        return None
    if isinstance(x, (int, float)):
        if isinstance(x, float) and (math.isnan(x) or math.isinf(x)):
            return None
        return float(x)
    s = str(x).strip()
    if s == "" or s == "." or s.lower() in {"na", "nan"}:
        return None
    try:
        v = float(s)
        if math.isnan(v) or math.isinf(v):
            return None
        return v
    except Exception:
        return None


def _read_table(path: Path, sheet: str = "") -> Tuple[List[str], List[List[str]]]:
    suf = path.suffix.lower()
    if suf in {".xlsx", ".xls"}:
        try:
            import pandas as pd  # type: ignore
        except Exception as e:
            raise RuntimeError("Reading .xlsx requires pandas + openpyxl") from e
        df = pd.read_excel(path, sheet_name=(sheet if sheet else 0))
        df.columns = [str(c) for c in df.columns]
        header = list(df.columns)
        rows = df.astype(object).where(df.notna(), None).values.tolist()
        return header, [["" if v is None else str(v) for v in row] for row in rows]

    # TSV/CSV
    # Try to sniff delimiter from the first KB
    with path.open("r", encoding="utf-8", errors="ignore", newline="") as f:
        sample = f.read(2048)
        f.seek(0)
        delim = "\t" if sample.count("\t") >= sample.count(",") else ","
        rdr = csv.reader(f, delimiter=delim)
        rows = list(rdr)
    if not rows:
        return [], []
    header = [str(x).strip() for x in rows[0]]
    body = rows[1:]
    return header, body


@dataclass
class MapRecord:
    marker: str
    chr: str
    bp: float
    cm: Optional[float] = None
    extra: Tuple[str, ...] = ()


def _read_plink_bim(bim_path: Path) -> List[MapRecord]:
    recs: List[MapRecord] = []
    with bim_path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            chr_ = parts[0]
            marker = parts[1]
            cm = _to_float(parts[2])
            bp = _to_float(parts[3])
            if bp is None:
                continue
            extra = tuple(parts[4:])
            recs.append(MapRecord(marker=marker, chr=str(chr_), bp=float(bp), cm=cm, extra=extra))
    return recs


def _guess_cols(header: List[str]) -> Dict[str, str]:
    hlow = [h.strip().lower() for h in header]

    def pick(cands: List[str]) -> str:
        for c in cands:
            if c in hlow:
                return header[hlow.index(c)]
        return ""

    chr_col = pick(["chr", "chrom", "chromosome", "scaffold", "lg"])
    bp_col = pick(["bp", "pos", "position", "physical", "physical_pos", "physpos"])
    cm_col = pick(["cm", "c_m", "c.m", "genpos", "genetic", "genetic_pos", "geneticpos"])
    marker_col = pick(["marker", "id", "snp", "rs", "name"])
    return {"chr": chr_col, "bp": bp_col, "cm": cm_col, "marker": marker_col}


def _read_map_generic(
    path: Path,
    *,
    sheet: str = "",
    auto_cols: bool = True,
    chr_col: str = "",
    bp_col: str = "",
    cm_col: str = "",
    marker_col: str = "",
) -> List[MapRecord]:
    header, rows = _read_table(path, sheet=sheet)
    if not header:
        raise RuntimeError(f"Empty map file: {path}")

    if auto_cols:
        guessed = _guess_cols(header)
        chr_col = chr_col or guessed.get("chr", "")
        bp_col = bp_col or guessed.get("bp", "")
        cm_col = cm_col or guessed.get("cm", "")
        marker_col = marker_col or guessed.get("marker", "")

    # fallback: GE marker_map.tsv => marker, chr, pos
    if not marker_col:
        marker_col = header[0]
    if not chr_col and len(header) >= 2:
        chr_col = header[1]
    if not bp_col and len(header) >= 3:
        bp_col = header[2]

    idx = {h: i for i, h in enumerate(header)}
    if chr_col not in idx or bp_col not in idx:
        raise RuntimeError(
            f"Failed to detect columns in {path}. header={header[:12]}... "
            f"Need chr_col and bp_col (or pos)."
        )

    chr_i = idx[chr_col]
    bp_i = idx[bp_col]
    cm_i = idx.get(cm_col, -1) if cm_col else -1
    mk_i = idx.get(marker_col, 0) if marker_col else 0

    recs: List[MapRecord] = []
    for r in rows:
        if not r:
            continue
        if len(r) < len(header):
            r = r + [""] * (len(header) - len(r))
        marker = str(r[mk_i]).strip() if mk_i < len(r) else ""
        chr_ = str(r[chr_i]).strip() if chr_i < len(r) else ""
        bp = _to_float(r[bp_i] if bp_i < len(r) else None)
        cm = _to_float(r[cm_i] if (cm_i >= 0 and cm_i < len(r)) else None)
        if not marker:
            marker = f"m{len(recs)+1}"
        if not chr_ or bp is None:
            continue
        recs.append(MapRecord(marker=marker, chr=chr_, bp=float(bp), cm=cm))
    return recs


def _pava_isotonic(y: List[float]) -> List[float]:
    """Pool Adjacent Violators Algorithm for isotonic regression (non-decreasing)."""
    n = len(y)
    if n == 0:
        return []
    blocks: List[Tuple[int, int, float, float]] = []  # s, e, mean, weight
    for i, v in enumerate(y):
        blocks.append((i, i, float(v), 1.0))
        while len(blocks) >= 2 and blocks[-2][2] > blocks[-1][2]:
            s1, e1, m1, w1 = blocks[-2]
            s2, e2, m2, w2 = blocks[-1]
            ww = w1 + w2
            mm = (m1 * w1 + m2 * w2) / ww if ww > 0 else (m1 + m2) / 2
            blocks = blocks[:-2]
            blocks.append((s1, e2, mm, ww))
    out = [0.0] * n
    for s, e, m, _w in blocks:
        for i in range(s, e + 1):
            out[i] = float(m)
    return out


@dataclass
class FitModel:
    method: str
    enforce_monotone: bool
    default_cM_per_Mb: float
    per_chr: Dict[str, Dict[str, List[float] | float]]
    global_slope: float


def _build_fit_model(
    ref: List[MapRecord],
    *,
    method: str,
    enforce_monotone: bool,
    default_cM_per_Mb: float,
) -> FitModel:
    by_chr: Dict[str, List[Tuple[float, float]]] = {}
    for r in ref:
        if r.cm is None:
            continue
        by_chr.setdefault(str(r.chr), []).append((float(r.bp), float(r.cm)))

    # global slope fallback
    global_slope = default_cM_per_Mb / 1e6
    all_pairs = [(bp, cm) for pairs in by_chr.values() for (bp, cm) in pairs]
    if len(all_pairs) >= 2:
        bps = [p[0] for p in all_pairs]
        cms = [p[1] for p in all_pairs]
        bp0, bp1 = min(bps), max(bps)
        cm0, cm1 = min(cms), max(cms)
        if bp1 > bp0:
            global_slope = (cm1 - cm0) / (bp1 - bp0)

    per_chr_model: Dict[str, Dict[str, List[float] | float]] = {}
    for chr_, pairs in by_chr.items():
        if len(pairs) < 2:
            continue
        pairs.sort(key=lambda x: x[0])
        bp = [p[0] for p in pairs]
        cm = [p[1] for p in pairs]
        # collapse duplicate bp
        bp2: List[float] = []
        cm2: List[float] = []
        i = 0
        while i < len(bp):
            j = i
            s = 0.0
            c = 0
            while j < len(bp) and bp[j] == bp[i]:
                s += cm[j]
                c += 1
                j += 1
            bp2.append(bp[i])
            cm2.append(s / c if c else cm[i])
            i = j
        bp, cm = bp2, cm2
        if enforce_monotone and len(cm) >= 2:
            cm = _pava_isotonic(cm)

        if method == "linear_chr":
            slope = global_slope
            if bp[-1] > bp[0]:
                slope = (cm[-1] - cm[0]) / (bp[-1] - bp[0])
            intercept = cm[0] - slope * bp[0]
            per_chr_model[str(chr_)] = {"slope": float(slope), "intercept": float(intercept)}
        else:
            per_chr_model[str(chr_)] = {"bp": bp, "cm": cm}

    return FitModel(
        method=method,
        enforce_monotone=enforce_monotone,
        default_cM_per_Mb=default_cM_per_Mb,
        per_chr=per_chr_model,
        global_slope=float(global_slope),
    )


def _predict_cm(model: FitModel, chr_: str, bp: float) -> float:
    chr_ = str(chr_)
    if model.method == "linear_global":
        return float(model.global_slope) * float(bp)

    m = model.per_chr.get(chr_)
    if m is None:
        return float(model.global_slope) * float(bp)

    if model.method == "linear_chr":
        slope = float(m.get("slope", model.global_slope))
        intercept = float(m.get("intercept", 0.0))
        return slope * float(bp) + intercept

    bps = list(m.get("bp", []))
    cms = list(m.get("cm", []))
    if len(bps) < 2:
        return float(model.global_slope) * float(bp)

    # extrapolate
    if bp <= bps[0]:
        b0, b1 = bps[0], bps[1]
        c0, c1 = cms[0], cms[1]
        slope = (c1 - c0) / (b1 - b0) if b1 != b0 else model.global_slope
        return c0 + slope * (bp - b0)
    if bp >= bps[-1]:
        b0, b1 = bps[-2], bps[-1]
        c0, c1 = cms[-2], cms[-1]
        slope = (c1 - c0) / (b1 - b0) if b1 != b0 else model.global_slope
        return c1 + slope * (bp - b1)

    # binary search interval
    lo, hi = 0, len(bps) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if bp < bps[mid]:
            hi = mid
        else:
            lo = mid
    b0, b1 = bps[lo], bps[hi]
    c0, c1 = cms[lo], cms[hi]
    if b1 == b0:
        return float(c0)
    t = (bp - b0) / (b1 - b0)
    return c0 + t * (c1 - c0)


def _write_outputs(
    out_dir: Path,
    base_out: Path,
    input_recs: List[MapRecord],
    cm_recs: List[Tuple[str, str, float]],
    report_rows: List[List],
    plink_bim_out: Optional[Path] = None,
):
    bp_rows = [(r.marker, r.chr, int(round(r.bp))) for r in input_recs]
    cm_rows = [(m, c, f"{cm:.6f}") for (m, c, cm) in cm_recs]
    use_cm = len(cm_recs) == len(input_recs) and len(cm_recs) > 0
    #main_rows = cm_rows if use_cm else bp_rows
    main_rows = []
    for cyc1 in range(0, len(bp_rows)):
      main_rows.append([bp_rows[cyc1][0], bp_rows[cyc1][1], bp_rows[cyc1][2], cm_rows[cyc1][2]])

    for root in {out_dir, base_out}:
        root.mkdir(parents=True, exist_ok=True)
        _write_tsv(root / "marker_map_bp.tsv", bp_rows, ["marker", "chr", "pos"])
        if use_cm:
            _write_tsv(root / "marker_map_cm.tsv", cm_rows, ["marker", "chr", "pos"])
        #_write_tsv(root / "marker_map.tsv", main_rows, ["marker", "chr", "pos"])
        _write_tsv(root / "marker_map.tsv", main_rows, ["marker", "chr", "pos", "cM"])
        _write_tsv(
            root / "marey_fit_report.tsv",
            report_rows,
            ["chr", "n_anchors", "bp_min", "bp_max", "cm_min", "cm_max", "method", "monotone"],
        )

        if plink_bim_out and plink_bim_out.exists():
            try:
                (root / plink_bim_out.name).write_text(
                    plink_bim_out.read_text(encoding="utf-8", errors="ignore"),
                    encoding="utf-8",
                )
            except Exception:
                pass

        artifacts = {
            "marker_map": str((root / "marker_map.tsv").resolve()),
            "marker_map_bp": str((root / "marker_map_bp.tsv").resolve()),
            "marey_fit_report": str((root / "marey_fit_report.tsv").resolve()),
        }
        if use_cm:
            artifacts["marker_map_cm"] = str((root / "marker_map_cm.tsv").resolve())
        if plink_bim_out and (root / plink_bim_out.name).exists():
            artifacts["plink_bim_cm"] = str((root / plink_bim_out.name).resolve())
        _write_json(root / "artifacts.json", artifacts)


def _make_unique_run_dir(parent: Path, prefix: str) -> Path:
    """Return a unique subfolder path under parent (not created here)."""
    parent = parent.expanduser().resolve()
    ts = time.strftime("%Y%m%d_%H%M%S")
    base = parent / f"{prefix}_{ts}"
    if not base.exists():
        return base
    for i in range(1, 1000):
        cand = parent / f"{prefix}_{ts}_{i}"
        if not cand.exists():
            return cand
    return parent / f"{prefix}_{ts}_{int(time.time())}"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--params", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    params = _read_json(Path(args.params))
    out_dir = Path(args.out).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    base_override_raw = str(params.get("out_dir_override") or "").strip()
    base_override = Path(base_override_raw).expanduser().resolve() if base_override_raw else out_dir
    base_override.mkdir(parents=True, exist_ok=True)

    # IMPORTANT UI/UX:
    # Map converter is commonly used inside the same output_folder as other Preprocess tools.
    # To avoid mixing files directly under the chosen folder, publish into a dedicated run subfolder.
    published_out = base_override
    if base_override_raw and base_override != out_dir:
        published_out = _make_unique_run_dir(base_override, "map_converter")
        published_out.mkdir(parents=True, exist_ok=True)

    input_kind = str(params.get("input_kind") or "auto").strip().lower()
    input_path = Path(str(params.get("input_path") or "").strip()).expanduser()
    if not str(input_path):
        raise SystemExit("input_path is required")

    # determine input records
    bim_path: Optional[Path] = None
    if input_kind in {"plink", "auto"}:
        if input_path.suffix.lower() == ".bim" and input_path.exists():
            bim_path = input_path
        elif (Path(str(input_path) + ".bim")).exists():
            bim_path = Path(str(input_path) + ".bim")

    if input_kind == "plink" or (input_kind == "auto" and bim_path is not None):
        if bim_path is None or not bim_path.exists():
            raise SystemExit("PLINK input requires .bim (or prefix with .bim)")
        input_recs = _read_plink_bim(bim_path)
        input_kind = "plink"
    else:
        if not input_path.exists():
            raise SystemExit(f"input_path not found: {input_path}")
        input_recs = _read_map_generic(input_path, sheet=str(params.get("input_sheet") or ""), auto_cols=True)
        input_kind = "marker_map"

    if not input_recs:
        raise SystemExit("No records found in input map")

    # reference map
    ref_path_raw = str(params.get("reference_map") or "").strip()
    ref_recs: List[MapRecord] = []
    if ref_path_raw:
        ref_path = Path(ref_path_raw).expanduser()
        if not ref_path.exists():
            raise SystemExit(f"reference_map not found: {ref_path}")
        ref_recs = _read_map_generic(
            ref_path,
            sheet=str(params.get("reference_sheet") or ""),
            auto_cols=_as_bool(params.get("ref_auto_cols"), True),
            chr_col=str(params.get("ref_chr_col") or ""),
            bp_col=str(params.get("ref_bp_col") or ""),
            cm_col=str(params.get("ref_cm_col") or ""),
            marker_col=str(params.get("ref_marker_col") or ""),
        )

    fit_method = str(params.get("fit_method") or "piecewise").strip().lower()
    if fit_method not in {"piecewise", "linear_chr", "linear_global"}:
        fit_method = "piecewise"

    enforce_monotone = _as_bool(params.get("enforce_monotone"), True)
    default_cM_per_Mb = float(params.get("default_cM_per_Mb") or 1.0)

    cm_recs: List[Tuple[str, str, float]] = []
    report_rows: List[List] = []
    plink_bim_out: Optional[Path] = None

    if ref_recs:
        model = _build_fit_model(
            ref_recs,
            method=fit_method,
            enforce_monotone=enforce_monotone,
            default_cM_per_Mb=default_cM_per_Mb,
        )

        for chr_, m in sorted(model.per_chr.items(), key=lambda x: str(x[0])):
            if fit_method == "linear_chr":
                report_rows.append([chr_, "(linear)", "", "", "", "", fit_method, str(enforce_monotone)])
            else:
                bps = list(m.get("bp", []))
                cms = list(m.get("cm", []))
                if len(bps) >= 2:
                    report_rows.append([
                        chr_, len(bps), int(round(min(bps))), int(round(max(bps))),
                        f"{min(cms):.6f}", f"{max(cms):.6f}", fit_method, str(enforce_monotone)
                    ])

        for r in input_recs:
            cm = _predict_cm(model, r.chr, r.bp)
            cm_recs.append((r.marker, r.chr, float(cm)))

        if _as_bool(params.get("write_plink_bim_cm"), False) and input_kind == "plink" and bim_path is not None:
            plink_bim_out = out_dir / "map_cm.bim"
            with plink_bim_out.open("w", encoding="utf-8", newline="") as f:
                for r, (_, _, cm) in zip(input_recs, cm_recs):
                    a = list(r.extra) if r.extra else ["A", "G"]
                    if len(a) < 2:
                        a = (a + ["A", "G"])[:2]
                    f.write(f"{r.chr}\t{r.marker}\t{cm:.6f}\t{int(round(r.bp))}\t{a[0]}\t{a[1]}\n")
    else:
        report_rows.append(["ALL", 0, "", "", "", "", "(no_reference)", ""])
        for r in input_recs:
            if r.cm is not None:
                cm_recs.append((r.marker, r.chr, float(r.cm)))

    _write_outputs(out_dir, published_out, input_recs, cm_recs, report_rows, plink_bim_out=plink_bim_out)

    # Enrich artifacts in the temporary out_dir so GUI can show the published path.
    try:
        art_tmp = out_dir / "artifacts.json"
        art = {}
        if art_tmp.exists():
            art = json.loads(art_tmp.read_text(encoding="utf-8", errors="ignore"))
        art["published_dir"] = str(published_out.resolve())
        art["published_marker_map"] = str((published_out / "marker_map.tsv").resolve())
        art["published_marker_map_bp"] = str((published_out / "marker_map_bp.tsv").resolve())
        if (published_out / "marker_map_cm.tsv").exists():
            art["published_marker_map_cm"] = str((published_out / "marker_map_cm.tsv").resolve())
        art["published_marey_fit_report"] = str((published_out / "marey_fit_report.tsv").resolve())
        _write_json(art_tmp, art)
    except Exception:
        pass

    print(f"[prep_map_converter] wrote outputs to: {published_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
