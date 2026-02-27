#!/usr/bin/env python3

"""QC Filtering plugin using PLINK2.

Inputs (params.json)
--------------------
plink_prefix: str
    PLINK prefix (either bfile prefix or pfile prefix). Suffixes .bed/.bim/.fam/.pgen/.pvar/.psam are accepted.
phenotype_tsv: str (optional)
    phenotype file; when provided, non-numeric columns are ignored and samples with phenotype missingness
    greater than max_missing_pheno are removed.

Filters
-------
max_missing_geno: float
min_call_rate: float
min_maf: float

Optional
--------
ld_prune: bool
ld_window: int (variant count)
ld_step: int (variant count)
ld_r2: float

segdist: bool
cross_type: one of f2, bc, dh, riself, risib
seg_p: float (keep if p >= seg_p)
het_max: float (max heterozygosity allowed when segdist enabled; always applied for dh/ri*)

Outputs (written under out_dir_override or --out)
------------------------------------------------
plink/filtered.{bed,bim,fam}
filter_summary.tsv
removed_samples.tsv
removed_markers.tsv
artifacts.json

Notes
-----
- Uses PLINK2 binary available on PATH by default.
- Always outputs PLINK1 binary set (bed/bim/fam) for downstream compatibility.
"""

import argparse
import json
import math
import os
import re
import shutil
import subprocess
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


def _read_text(path: Path) -> str:
    try:
        return path.read_text(encoding="utf-8", errors="ignore")
    except Exception:
        return ""


def _write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text or "", encoding="utf-8", errors="ignore")


def _run(cmd: Sequence[str], cwd: Optional[Path], stdout_path: Path, stderr_path: Path) -> None:
    stdout_path.parent.mkdir(parents=True, exist_ok=True)
    stderr_path.parent.mkdir(parents=True, exist_ok=True)
    proc = subprocess.run(
        list(cmd),
        cwd=str(cwd) if cwd else None,
        capture_output=True,
        text=True,
        check=False,
    )
    _write_text(stdout_path, proc.stdout or "")
    _write_text(stderr_path, proc.stderr or "")
    if proc.returncode != 0:
        tail = (proc.stderr or proc.stdout or "")[-2000:]
        raise RuntimeError(f"Command failed (code={proc.returncode}): {' '.join(cmd)}\n--- tail ---\n{tail}")


def _resolve_plink_prefix(prefix_like: str) -> Tuple[str, str]:
    """Return (mode, prefix) where mode is '--bfile' or '--pfile'."""
    p = Path(prefix_like).expanduser()
    sfx = p.suffix.lower()
    if sfx in {".bed", ".bim", ".fam", ".pgen", ".pvar", ".psam"}:
        p = p.with_suffix("")

    # prefer pfile when both exist
    if (Path(str(p) + ".pgen")).exists() and (Path(str(p) + ".pvar")).exists() and (Path(str(p) + ".psam")).exists():
        return "--pfile", str(p)
    if (Path(str(p) + ".bed")).exists() and (Path(str(p) + ".bim")).exists() and (Path(str(p) + ".fam")).exists():
        return "--bfile", str(p)

    # fallback: allow prefix without checking all companions
    if (Path(str(p) + ".pgen")).exists():
        return "--pfile", str(p)
    return "--bfile", str(p)


def _read_samples_from_prefix(mode: str, prefix: str) -> List[Tuple[str, str]]:
    """Return list of (FID, IID) in file order."""
    p = Path(prefix)
    if mode == "--bfile":
        fam = Path(str(p) + ".fam")
        rows: List[Tuple[str, str]] = []
        if fam.exists():
            for line in fam.read_text(encoding="utf-8", errors="ignore").splitlines():
                if not line.strip():
                    continue
                toks = line.split()
                if len(toks) >= 2:
                    rows.append((toks[0], toks[1]))
        return rows

    psam = Path(str(p) + ".psam")
    rows = []
    if psam.exists():
        lines = psam.read_text(encoding="utf-8", errors="ignore").splitlines()
        # header begins with '#'
        for ln in lines:
            if not ln.strip() or ln.startswith("#"):
                continue
            toks = ln.split()
            if len(toks) >= 2:
                rows.append((toks[0], toks[1]))
    return rows


def _read_table_auto(path: Path):
    import pandas as pd

    # try tab first
    for sep in ["\t", ","]:
        try:
            df = pd.read_csv(path, sep=sep, dtype=str)
            if df.shape[1] >= 2:
                return df
        except Exception:
            continue
    # last resort
    return pd.read_csv(path, sep=None, engine="python", dtype=str)


def _detect_numeric_trait_cols(df) -> List[str]:
    import pandas as pd

    cols = [str(c) for c in df.columns]
    if len(cols) <= 1:
        return []

    out: List[str] = []
    # First col is treated as ID
    for c in cols[1:]:
        s = df[c]
        # treat empty strings as NA
        s2 = s.replace({"": None, "NA": None, "NaN": None, "nan": None})
        non_empty = s2.notna()
        if int(non_empty.sum()) < 5:
            continue
        num = pd.to_numeric(s2, errors="coerce")
        ok = num.notna() & non_empty
        # accept if >=95% of non-empty entries are numeric
        denom = int(non_empty.sum())
        if denom == 0:
            continue
        if float(ok.sum()) / float(denom) < 0.95:
            continue
        # need variation
        if num.dropna().nunique() < 2:
            continue
        out.append(c)
    return out


def _phenotype_keep_list(
    phenotype_tsv: Path,
    samples_fid_iid: List[Tuple[str, str]],
    max_missing_pheno: float,
    out_dir: Path,
) -> Tuple[Optional[Path], Optional[Path], Dict[str, float]]:
    """Return (keep_file_path, filtered_pheno_path, stats)."""
    if not phenotype_tsv or not phenotype_tsv.exists():
        return None, None, {"n_pheno_samples": 0, "n_numeric_traits": 0}

    import pandas as pd

    df = _read_table_auto(phenotype_tsv)
    df.columns = [str(c).strip() for c in df.columns]
    if df.shape[1] < 1:
        return None, None, {"n_pheno_samples": 0, "n_numeric_traits": 0}

    id_col = df.columns[0]
    df[id_col] = df[id_col].astype(str)

    trait_cols = _detect_numeric_trait_cols(df)
    stats: Dict[str, float] = {
        "n_pheno_samples": float(df.shape[0]),
        "n_numeric_traits": float(len(trait_cols)),
    }

    if not trait_cols:
        # no numeric traits => do not filter by phenotype
        return None, None, stats

    sub = df[[id_col] + trait_cols].copy()
    # coerce
    for c in trait_cols:
        sub[c] = pd.to_numeric(sub[c].replace({"": None, "NA": None, "NaN": None, "nan": None}), errors="coerce")

    miss = sub[trait_cols].isna().mean(axis=1)
    keep_mask = miss <= float(max_missing_pheno)
    keep_ids = set(sub.loc[keep_mask, id_col].astype(str).tolist())

    # match to genotype IIDs/FIDs
    iid_to_fid: Dict[str, str] = {iid: fid for fid, iid in samples_fid_iid}
    fid_set = {fid for fid, _ in samples_fid_iid}

    keep_rows: List[Tuple[str, str]] = []
    n_matched = 0
    for sid in keep_ids:
        if sid in iid_to_fid:
            keep_rows.append((iid_to_fid[sid], sid))
            n_matched += 1
        elif sid in fid_set:
            # if phenotype uses FID, keep all rows with that FID
            for fid, iid in samples_fid_iid:
                if fid == sid:
                    keep_rows.append((fid, iid))
                    n_matched += 1

    stats["n_pheno_keep_ids"] = float(len(keep_ids))
    stats["n_pheno_keep_matched"] = float(n_matched)

    if not keep_rows:
        return None, None, stats

    keep_path = out_dir / "_tmp" / "keep_samples.txt"
    keep_path.parent.mkdir(parents=True, exist_ok=True)
    keep_path.write_text("\n".join([f"{fid}\t{iid}" for fid, iid in keep_rows]) + "\n", encoding="utf-8")

    # write filtered phenotype (subset by IID match)
    kept_iids = {iid for _, iid in keep_rows}
    df_f = df[df[id_col].astype(str).isin(kept_iids)].copy()
    ph_out = out_dir / "phenotype_filtered.tsv"
    df_f.to_csv(ph_out, sep="\t", index=False)

    return keep_path, ph_out, stats


def _load_bim(prefix: Path) -> List[str]:
    bim = Path(str(prefix) + ".bim")
    ids: List[str] = []
    if bim.exists():
        for ln in bim.read_text(encoding="utf-8", errors="ignore").splitlines():
            if not ln.strip():
                continue
            toks = ln.split()
            if len(toks) >= 2:
                ids.append(toks[1])
    return ids


def _load_fam(prefix: Path) -> List[Tuple[str, str]]:
    fam = Path(str(prefix) + ".fam")
    out: List[Tuple[str, str]] = []
    if fam.exists():
        for ln in fam.read_text(encoding="utf-8", errors="ignore").splitlines():
            if not ln.strip():
                continue
            toks = ln.split()
            if len(toks) >= 2:
                out.append((toks[0], toks[1]))
    return out


def _chi2_pvalue(stat: float, df: int) -> float:
    """Survival function for chi-square (df=1 or df=2)."""
    if not math.isfinite(stat) or stat < 0:
        return 1.0
    if df == 1:
        # p = erfc(sqrt(x/2))
        return float(math.erfc(math.sqrt(stat / 2.0)))
    if df == 2:
        return float(math.exp(-stat / 2.0))
    # fallback (conservative)
    return float(math.exp(-stat / 2.0))


def _segdist_filter(
    plink2_bin: str,
    prefix: Path,
    out_dir: Path,
    cross_type: str,
    p_keep: float,
    het_max: float,
) -> Optional[Path]:
    """Return snplist to exclude (variants failing)."""
    tmp = out_dir / "_tmp"
    tmp.mkdir(parents=True, exist_ok=True)

    raw_prefix = tmp / "geno_export"
    # export additive dosages to .raw
    _run(
        [plink2_bin, "--bfile", str(prefix), "--allow-extra-chr", "--export", "A", "--out", str(raw_prefix)],
        cwd=None,
        stdout_path=tmp / "plink_export_stdout.txt",
        stderr_path=tmp / "plink_export_stderr.txt",
    )

    raw_path = Path(str(raw_prefix) + ".raw")
    if not raw_path.exists():
        return None

    import pandas as pd

    df = pd.read_csv(raw_path, sep="\s+", dtype=str)
    # marker columns start after 6 leading cols in PLINK raw
    lead = ["FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"]
    mcols = [c for c in df.columns if c not in lead]
    if not mcols:
        return None

    # convert to numeric
    G = df[mcols].replace({"NA": None, "-9": None, "": None}).apply(pd.to_numeric, errors="coerce")

    rows = []
    fail = []

    ct = (cross_type or "").strip().lower()

    for c in mcols:
        g = G[c]
        non = g.dropna()
        n = int(non.shape[0])
        if n < 10:
            # too few
            rows.append((c, n, 1.0, 0.0, 0.0))
            continue

        c0 = int((non == 0).sum())
        c1 = int((non == 1).sum())
        c2 = int((non == 2).sum())
        het = c1 / float(n) if n > 0 else 0.0

        # default: no filtering
        p = 1.0

        if ct == "f2":
            # expected 1:2:1
            exp = [0.25 * n, 0.5 * n, 0.25 * n]
            obs = [c0, c1, c2]
            stat = 0.0
            for o, e in zip(obs, exp):
                if e > 0:
                    stat += (o - e) ** 2 / e
            p = _chi2_pvalue(stat, df=2)

        elif ct in {"bc", "dh", "riself", "risib"}:
            # two-class tests on dominant expected 1:1 between two main genotypes
            # determine which genotypes are the two main ones
            counts = {0: c0, 1: c1, 2: c2}
            # drop the smallest class
            top = sorted(counts.items(), key=lambda kv: kv[1], reverse=True)
            gA, a = top[0]
            gB, b = top[1]

            # heterozygosity guard for DH/RIL
            if ct in {"dh", "riself", "risib"}:
                if het > float(het_max):
                    p = 0.0
                else:
                    # prefer homozygotes 0/2 when present
                    a = c0
                    b = c2
                    if a + b < max(10, int(0.6 * n)):
                        # fallback to top-2
                        a = top[0][1]
                        b = top[1][1]
                    n2 = a + b
                    if n2 >= 10:
                        exp = [0.5 * n2, 0.5 * n2]
                        stat = (a - exp[0]) ** 2 / exp[0] + (b - exp[1]) ** 2 / exp[1]
                        p = _chi2_pvalue(stat, df=1)
                    else:
                        p = 1.0
            else:
                # BC: use top-2
                n2 = a + b
                if n2 >= 10:
                    exp = [0.5 * n2, 0.5 * n2]
                    stat = (a - exp[0]) ** 2 / exp[0] + (b - exp[1]) ** 2 / exp[1]
                    p = _chi2_pvalue(stat, df=1)
                else:
                    p = 1.0

        rows.append((c, n, p, het, float(c1)))
        if float(p) < float(p_keep):
            fail.append(c)
        # additional heterozygosity filter for dh/ri
        if ct in {"dh", "riself", "risib"} and het > float(het_max):
            if c not in fail:
                fail.append(c)

    # write report
    rep = out_dir / "segdist.tsv"
    rep.parent.mkdir(parents=True, exist_ok=True)
    rep.write_text("marker\tn_nonmissing\tp_value\thet_rate\thet_count\n", encoding="utf-8")
    with rep.open("a", encoding="utf-8") as f:
        for r in rows:
            f.write(f"{r[0]}\t{r[1]}\t{r[2]:.6g}\t{r[3]:.6g}\t{r[4]:.0f}\n")

    if not fail:
        return None

    snplist = out_dir / "_tmp" / "segdist_fail.snplist"
    snplist.write_text("\n".join(fail) + "\n", encoding="utf-8")
    return snplist


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--params", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    params = json.loads(Path(args.params).read_text(encoding="utf-8"))

    out_dir = Path(params.get("out_dir_override") or args.out).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    plink2_bin = str(params.get("plink2_bin") or "plink2").strip() or "plink2"

    mode, plink_prefix = _resolve_plink_prefix(str(params.get("plink_prefix") or "").strip())
    if not plink_prefix:
        raise ValueError("plink_prefix is required")

    # thresholds
    max_missing_geno = float(params.get("max_missing_geno", 0.2))
    max_missing_pheno = float(params.get("max_missing_pheno", 0.5))
    min_maf = float(params.get("min_maf", 0.05))
    min_call_rate = float(params.get("min_call_rate", 0.9))

    ld_prune = bool(params.get("ld_prune", False))
    ld_window = int(params.get("ld_window", 100))
    ld_step = int(params.get("ld_step", 10))
    ld_r2 = float(params.get("ld_r2", 0.2))

    segdist = bool(params.get("segdist", False))
    cross_type = str(params.get("cross_type") or "f2")
    seg_p = float(params.get("seg_p", 0.001))
    het_max = float(params.get("het_max", 0.2))

    out_name = str(params.get("out_name") or "filtered")
    out_name = re.sub(r"[^A-Za-z0-9_.-]+", "_", out_name)

    tmp = out_dir / "_tmp"
    tmp.mkdir(parents=True, exist_ok=True)

    # phenotype keep file
    samples = _read_samples_from_prefix(mode, plink_prefix)
    phenotype_path = str(params.get("phenotype_tsv") or "").strip()
    keep_file = None
    filtered_pheno_path = None
    ph_stats: Dict[str, float] = {}
    if phenotype_path:
        keep_file, filtered_pheno_path, ph_stats = _phenotype_keep_list(Path(phenotype_path), samples, max_missing_pheno, out_dir)

    # Stage 1: basic QC filters
    plink_dir = out_dir / "plink"
    plink_dir.mkdir(parents=True, exist_ok=True)

    stage1 = plink_dir / f"{out_name}_stage1"
    stage1.parent.mkdir(parents=True, exist_ok=True)

    cmd = [plink2_bin, mode, plink_prefix, "--allow-extra-chr"]
    # sample missingness
    if max_missing_geno is not None and max_missing_geno > 0:
        cmd += ["--mind", f"{max_missing_geno:.6g}"]
    # marker call rate
    marker_max_miss = max(0.0, min(0.999, 1.0 - float(min_call_rate)))
    cmd += ["--geno", f"{marker_max_miss:.6g}"]
    # maf
    if min_maf is not None and min_maf > 0:
        cmd += ["--maf", f"{min_maf:.6g}"]
    # phenotype keep
    if keep_file is not None and keep_file.exists():
        cmd += ["--keep", str(keep_file)]

    cmd += ["--make-bed", "--out", str(stage1)]

    _write_text(out_dir / "cmd_stage1.txt", " ".join(cmd))
    _run(cmd, cwd=None, stdout_path=tmp / "stage1_stdout.txt", stderr_path=tmp / "stage1_stderr.txt")

    current = stage1

    # LD pruning
    if ld_prune:
        prune_prefix = tmp / "prune"
        cmd_prune = [plink2_bin, "--bfile", str(current), "--allow-extra-chr", "--indep-pairwise", str(int(ld_window)), str(int(ld_step)), f"{ld_r2:.6g}", "--out", str(prune_prefix)]
        _write_text(out_dir / "cmd_ld_prune.txt", " ".join(cmd_prune))
        _run(cmd_prune, cwd=None, stdout_path=tmp / "ld_stdout.txt", stderr_path=tmp / "ld_stderr.txt")

        prune_in = Path(str(prune_prefix) + ".prune.in")
        if prune_in.exists():
            stage2 = plink_dir / f"{out_name}_stage2_ld"
            cmd_extract = [plink2_bin, "--bfile", str(current), "--allow-extra-chr", "--extract", str(prune_in), "--make-bed", "--out", str(stage2)]
            _write_text(out_dir / "cmd_ld_apply.txt", " ".join(cmd_extract))
            _run(cmd_extract, cwd=None, stdout_path=tmp / "ld_apply_stdout.txt", stderr_path=tmp / "ld_apply_stderr.txt")
            current = stage2

    # Segregation distortion
    if segdist:
        snplist = _segdist_filter(plink2_bin, current, out_dir, cross_type, seg_p, het_max)
        if snplist is not None and snplist.exists():
            stage3 = plink_dir / f"{out_name}_stage3_seg"
            cmd_seg = [plink2_bin, "--bfile", str(current), "--allow-extra-chr", "--exclude", str(snplist), "--make-bed", "--out", str(stage3)]
            _write_text(out_dir / "cmd_seg_apply.txt", " ".join(cmd_seg))
            _run(cmd_seg, cwd=None, stdout_path=tmp / "seg_apply_stdout.txt", stderr_path=tmp / "seg_apply_stderr.txt")
            current = stage3

    # Final prefix
    final_prefix = plink_dir / out_name
    cmd_final = [plink2_bin, "--bfile", str(current), "--allow-extra-chr", "--make-bed", "--out", str(final_prefix)]
    _write_text(out_dir / "cmd_final.txt", " ".join(cmd_final))
    _run(cmd_final, cwd=None, stdout_path=tmp / "final_stdout.txt", stderr_path=tmp / "final_stderr.txt")

    # Removed samples/markers (best-effort)
    before_fam = _load_fam(Path(plink_prefix)) if mode == "--bfile" else _load_fam(Path(str(stage1)))
    after_fam = _load_fam(final_prefix)
    before_iids = {iid for _, iid in before_fam} if before_fam else set(iid for _, iid in samples)
    after_iids = {iid for _, iid in after_fam}
    removed_samples = sorted(list(before_iids - after_iids))

    before_bim = _load_bim(Path(plink_prefix)) if mode == "--bfile" else _load_bim(Path(str(stage1)))
    after_bim = _load_bim(final_prefix)
    removed_markers = sorted(list(set(before_bim) - set(after_bim))) if before_bim and after_bim else []

    (out_dir / "removed_samples.tsv").write_text("sample_id\n" + "\n".join(removed_samples) + "\n", encoding="utf-8")
    (out_dir / "removed_markers.tsv").write_text("marker_id\n" + "\n".join(removed_markers) + "\n", encoding="utf-8")

    # Summary
    n_before_samples = len(before_iids) if before_iids else len(samples)
    n_after_samples = len(after_iids)
    n_before_markers = len(before_bim) if before_bim else 0
    n_after_markers = len(after_bim) if after_bim else 0

    summary = [
        ("n_samples_before", n_before_samples),
        ("n_samples_after", n_after_samples),
        ("n_markers_before", n_before_markers),
        ("n_markers_after", n_after_markers),
        ("max_missing_geno", max_missing_geno),
        ("max_missing_pheno", max_missing_pheno),
        ("min_maf", min_maf),
        ("min_call_rate", min_call_rate),
        ("ld_prune", int(ld_prune)),
        ("ld_window", ld_window),
        ("ld_step", ld_step),
        ("ld_r2", ld_r2),
        ("segdist", int(segdist)),
        ("cross_type", cross_type),
        ("seg_p_keep", seg_p),
        ("het_max", het_max),
        ("plink2_bin", plink2_bin),
        ("input_mode", mode),
        ("input_prefix", plink_prefix),
        ("output_prefix", str(final_prefix)),
    ]
    for k, v in ph_stats.items():
        summary.append((k, v))

    out_summary = out_dir / "filter_summary.tsv"
    out_summary.write_text("key\tvalue\n", encoding="utf-8")
    with out_summary.open("a", encoding="utf-8") as f:
        for k, v in summary:
            f.write(f"{k}\t{v}\n")

    artifacts = {
        "filtered_plink_prefix": str(final_prefix),
        "filtered_phenotype_tsv": str(filtered_pheno_path) if filtered_pheno_path else "",
        "filter_summary": str(out_summary),
        "removed_samples": str(out_dir / "removed_samples.tsv"),
        "removed_markers": str(out_dir / "removed_markers.tsv"),
    }
    _write_text(out_dir / "artifacts.json", json.dumps(artifacts, ensure_ascii=False, indent=2))

    # convenience pointers
    _write_text(out_dir / "filtered_plink_prefix.txt", str(final_prefix))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
