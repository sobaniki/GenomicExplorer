#!/usr/bin/env python3
"""Preprocess: Merge (phenotype + genotype)

Goal
----
Merge multiple *already matched* phenotype/genotype outputs (typically produced
by Split/Merge/Extract + Extract). This plugin does **not** attempt to match or
rename IDs; it assumes inputs are already internally consistent (per dataset).

However, when merging multiple datasets, the following issues can arise:
- Genotype files may have different marker sets (union vs intersection)
- Sample IDs may overlap across genotype files (duplicate sample names)
- Phenotype files may have different trait columns
- Trait names may collide across files
- Sample IDs may overlap across phenotype files

This plugin provides conservative, explicit policies for these cases and
produces detailed reports.

Inputs (params.json)
--------------------
out_dir_override: str (optional)

Phenotype
- phenotype_files: list[str] (optional)
- phenotype_sheet: str (optional; applied to all xlsx)
- phenotype_id_col: str (optional; default 'id' or first column)
- phenotype_merge_mode: str (auto|hstack|vstack) default auto
- phenotype_join: str (inner|outer|left) default outer (for hstack)
- phenotype_dup_policy: str (mean_numeric|first|keep_all) default mean_numeric
- trait_collision: str (prefix|suffix|error) default prefix

Genotype
- genotype_mode: str (none|vcf|plink) default none
- genotype_inputs: list[str] (vcf paths if mode=vcf; PLINK prefixes if mode=plink)
- bcftools_bin: str (default bcftools)
- plink2_bin: str (default plink2)
- genotype_merge_strategy: str (auto|merge_samples|concat_variants) default auto
- dup_sample_policy: str (rename|drop_later|error) default rename
- output_genotype: str (vcf|plink|both) default vcf
- drop_duplicate_variants: bool default False

Final sync
- final_sample_policy: str (intersection|none) default intersection

Outputs
-------
phenotype_merged.tsv
phenotype_merge_report.tsv

If genotype_mode != none:
- genotype_merged.vcf.gz (+.tbi)
- genotype_merge_report.tsv
- logs/*

If final_sample_policy=intersection and both phenotype+genotype exist:
- phenotype_final.tsv
- genotype_final.vcf.gz (+.tbi)
- final_sync_report.tsv

Notes
-----
For PLINK inputs, this plugin converts each dataset to a bgzipped VCF via
PLINK2 and performs the merge in VCF space using bcftools. This avoids
implementing two parallel merge stacks.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple


# ----------------- small utilities -----------------

def _write_text(p: Path, s: str) -> None:
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(s, encoding="utf-8")


def _read_text(p: Path) -> str:
    return p.read_text(encoding="utf-8", errors="ignore")


def _run_cmd(cmd: List[str], *, cwd: Optional[Path] = None, stdout_path: Optional[Path] = None, stderr_path: Optional[Path] = None) -> int:
    if stdout_path is not None:
        stdout_path.parent.mkdir(parents=True, exist_ok=True)
    if stderr_path is not None:
        stderr_path.parent.mkdir(parents=True, exist_ok=True)

    out_f = stdout_path.open("w", encoding="utf-8") if stdout_path else open(os.devnull, "w", encoding="utf-8")
    err_f = stderr_path.open("w", encoding="utf-8") if stderr_path else open(os.devnull, "w", encoding="utf-8")
    try:
        r = subprocess.run(cmd, cwd=str(cwd) if cwd else None, stdout=out_f, stderr=err_f)
        return int(r.returncode)
    finally:
        try:
            out_f.close()
        except Exception:
            pass
        try:
            err_f.close()
        except Exception:
            pass


def _safe_label(path_or_prefix: str) -> str:
    """Create a short label used for renaming columns/samples."""
    s = Path(path_or_prefix).name
    # For PLINK prefix, remove common extensions
    for ext in (".bed", ".bim", ".fam", ".pgen", ".pvar", ".psam", ".vcf", ".gz"):
        if s.endswith(ext):
            s = s[: -len(ext)]
    s = re.sub(r"[^A-Za-z0-9_.-]+", "_", s)
    return s or "file"


# ----------------- phenotype merge -----------------

def _read_table(path: Path, sheet: str = ""):
    import pandas as pd

    suf = path.suffix.lower()
    if suf in (".xlsx", ".xls"):
        df = pd.read_excel(path, sheet_name=sheet if sheet else 0, engine="openpyxl")
    else:
        # autodetect separator by extension
        if suf in (".csv",):
            df = pd.read_csv(path, dtype=object)
        else:
            df = pd.read_csv(path, sep="\t", dtype=object)
    return df


def _ensure_id_col(df, id_col: str) -> str:
    cols = list(df.columns)
    if id_col and id_col in cols:
        return id_col
    # default: 'id' if exists else first col
    if "id" in cols:
        return "id"
    return cols[0]


def _reduce_duplicates(df, id_col: str, policy: str):
    import pandas as pd

    df = df.copy()
    df = df[df[id_col].notna()]
    df[id_col] = df[id_col].astype(str)

    # count duplicates
    dup_n = int(df.duplicated(subset=[id_col]).sum())
    if dup_n == 0 or policy == "keep_all":
        return df, dup_n

    # convert possible numeric columns
    # Keep id_col as str; try numeric conversion for others
    other_cols = [c for c in df.columns if c != id_col]
    num_df = df[other_cols].apply(pd.to_numeric, errors="coerce")

    if policy == "first":
        out = df.drop_duplicates(subset=[id_col], keep="first")
        return out, dup_n

    if policy == "mean_numeric":
        # numeric mean per id; for non-numeric keep first non-null
        # Build aggregation dict
        agg = {}
        for c in other_cols:
            if num_df[c].notna().any():
                agg[c] = "mean"
            else:
                agg[c] = lambda x: x.dropna().astype(str).iloc[0] if len(x.dropna()) else ""
        g = df.groupby(id_col, as_index=False).agg(agg)
        # Ensure id first
        cols = [id_col] + [c for c in g.columns if c != id_col]
        g = g[cols]
        return g, dup_n

    # fallback
    out = df.drop_duplicates(subset=[id_col], keep="first")
    return out, dup_n


def _phenotype_merge(phenotype_files: List[str], out_dir: Path, *, sheet: str, id_col: str, merge_mode: str, join: str,
                     dup_policy: str, trait_collision: str) -> Tuple[Optional[Path], Optional[Path]]:
    import pandas as pd

    if not phenotype_files:
        return None, None

    rows = []
    dfs = []
    for fp in phenotype_files:
        p = Path(fp)
        df = _read_table(p, sheet=sheet)
        idc = _ensure_id_col(df, id_col)
        df, dup_n = _reduce_duplicates(df, idc, dup_policy)
        df = df.copy()
        df[idc] = df[idc].astype(str)
        label = _safe_label(fp)
        rows.append({
            "file": str(p),
            "label": label,
            "n_rows": int(df.shape[0]),
            "n_cols": int(df.shape[1]),
            "id_col": idc,
            "n_ids": int(df[idc].nunique()),
            "dup_rows_removed": int(dup_n) if dup_policy != "keep_all" else 0,
        })
        dfs.append((label, idc, df))

    report_path = out_dir / "phenotype_merge_report.tsv"
    pd.DataFrame(rows).to_csv(report_path, sep="\t", index=False)

    if merge_mode == "auto":
        # Heuristic: if sample overlaps are high, do hstack; else vstack
        # Compute avg Jaccard overlap
        sets = [set(df[idc].astype(str).tolist()) for _, idc, df in dfs]
        if len(sets) <= 1:
            mode = "vstack"
        else:
            inter = set.intersection(*sets)
            union = set.union(*sets)
            j = (len(inter) / len(union)) if union else 0.0
            mode = "hstack" if j >= 0.5 else "vstack"
        merge_mode_eff = mode
    else:
        merge_mode_eff = merge_mode

    if merge_mode_eff == "vstack":
        # Align columns (union); preserve id col as 'id'
        normed = []
        all_cols = set()
        for label, idc, df in dfs:
            df2 = df.rename(columns={idc: "id"})
            all_cols |= set(df2.columns)
            normed.append((label, df2))
        all_cols = ["id"] + [c for c in sorted(all_cols) if c != "id"]
        df_outs = []
        for label, df2 in normed:
            for c in all_cols:
                if c not in df2.columns:
                    df2[c] = ""
            df2 = df2[all_cols]
            df2["_source"] = label
            df_outs.append(df2)
        out = pd.concat(df_outs, axis=0, ignore_index=True)

        # If duplicates in final: reduce unless keep_all
        if dup_policy != "keep_all":
            out, _ = _reduce_duplicates(out, "id", dup_policy)

    else:
        # hstack join by id
        base_label, base_idc, base_df = dfs[0]
        out = base_df.rename(columns={base_idc: "id"}).copy()
        out["id"] = out["id"].astype(str)
        used_cols = set(out.columns)

        collisions = []

        for label, idc, df in dfs[1:]:
            df2 = df.rename(columns={idc: "id"}).copy()
            df2["id"] = df2["id"].astype(str)
            cols = [c for c in df2.columns if c != "id"]
            rename_map = {}
            for c in cols:
                if c in used_cols:
                    collisions.append({"trait": c, "label": label})
                    if trait_collision == "error":
                        raise RuntimeError(f"Trait column collision: '{c}' already exists; set trait_collision=prefix/suffix")
                    if trait_collision == "prefix":
                        rename_map[c] = f"{label}__{c}"
                    else:
                        rename_map[c] = f"{c}__{label}"
            if rename_map:
                df2 = df2.rename(columns=rename_map)
            # join
            how = "outer" if join == "outer" else ("inner" if join == "inner" else "left")
            out = out.merge(df2, on="id", how=how)
            used_cols |= set(out.columns)

        # write collision report
        if collisions:
            pd.DataFrame(collisions).to_csv(out_dir / "phenotype_trait_collisions.tsv", sep="\t", index=False)

        # optionally reduce duplicates caused by 'keep_all' earlier
        if dup_policy != "keep_all":
            out, _ = _reduce_duplicates(out, "id", dup_policy)

    ph_out = out_dir / "phenotype_merged.tsv"
    out.to_csv(ph_out, sep="\t", index=False)

    _write_text(out_dir / "phenotype_merge_mode.txt", merge_mode_eff + "\n")

    return ph_out, report_path


# ----------------- genotype merge -----------------

def _bcftools_list_samples(bcftools: str, vcf: Path, log_dir: Path) -> List[str]:
    tmp_out = log_dir / f"samples_{vcf.name}.txt"
    tmp_err = log_dir / f"samples_{vcf.name}.stderr.txt"
    rc = _run_cmd([bcftools, "query", "-l", str(vcf)], stdout_path=tmp_out, stderr_path=tmp_err)
    if rc != 0:
        raise RuntimeError(f"bcftools query -l failed for {vcf}. See {tmp_err}")
    s = _read_text(tmp_out).splitlines()
    return [x.strip() for x in s if x.strip()]


def _ensure_bgz_and_index(bcftools: str, tabix: str, vcf_in: Path, out_dir: Path, log_dir: Path) -> Path:
    # If already .vcf.gz, assume bgz; otherwise bgzip via bcftools view -Oz
    v = vcf_in
    if v.suffixes[-2:] == [".vcf", ".gz"] or v.name.endswith(".vcf.gz"):
        vgz = v
    else:
        vgz = out_dir / (v.name + ".vcf.gz")
        rc = _run_cmd([bcftools, "view", "-Oz", "-o", str(vgz), str(v)], stdout_path=log_dir / "bcftools_view_stdout.txt", stderr_path=log_dir / "bcftools_view_stderr.txt")
        if rc != 0:
            raise RuntimeError(f"bcftools view -Oz failed for {v}. See logs")
    # index
    rc = _run_cmd([tabix, "-f", "-p", "vcf", str(vgz)], stdout_path=log_dir / "tabix_stdout.txt", stderr_path=log_dir / "tabix_stderr.txt")
    if rc != 0:
        # try bcftools index
        rc2 = _run_cmd([bcftools, "index", "-f", "-t", str(vgz)], stdout_path=log_dir / "bcftools_index_stdout.txt", stderr_path=log_dir / "bcftools_index_stderr.txt")
        if rc2 != 0:
            raise RuntimeError(f"Failed to index VCF {vgz}. See logs")
    return vgz


def _plink_to_vcf_bgz(plink2: str, prefix: str, out_dir: Path, label: str, log_dir: Path) -> Path:
    outp = out_dir / f"{label}.plink_export"
    cmd = [plink2, "--bfile", prefix, "--export", "vcf", "bgz", "--out", str(outp)]
    rc = _run_cmd(cmd, stdout_path=log_dir / f"plink2_export_{label}_stdout.txt", stderr_path=log_dir / f"plink2_export_{label}_stderr.txt")
    if rc != 0:
        raise RuntimeError(f"plink2 export vcf failed for {prefix}. See logs")
    vcf = outp.with_suffix(".vcf.gz")
    if not vcf.exists():
        # some plink2 versions output .vcf.gz directly without .vcf.gz suffix logic
        # keep robust: search
        cand = list(out_dir.glob(f"{label}.plink_export*.vcf.gz"))
        if cand:
            vcf = cand[0]
        else:
            raise RuntimeError("plink2 export did not produce .vcf.gz")
    return vcf


def _rename_vcf_samples_if_needed(bcftools: str, vcf: Path, label: str, existing: set, policy: str, out_dir: Path, log_dir: Path) -> Tuple[Path, List[str]]:
    samples = _bcftools_list_samples(bcftools, vcf, log_dir)
    overlap = [s for s in samples if s in existing]
    if not overlap:
        existing.update(samples)
        return vcf, samples

    if policy == "error":
        raise RuntimeError(f"Duplicate sample IDs across genotype inputs: {overlap[:10]} ...")

    if policy == "drop_later":
        # Keep only non-overlapping samples from this file
        keep = [s for s in samples if s not in existing]
        keep_path = out_dir / f"keep_no_dups_{label}.txt"
        _write_text(keep_path, "\n".join(keep) + "\n")
        out_vcf = out_dir / f"{vcf.stem}.nodups.vcf.gz"
        rc = _run_cmd([bcftools, "view", "-S", str(keep_path), "-Oz", "-o", str(out_vcf), str(vcf)], stdout_path=log_dir / f"bcftools_dropdups_{label}_stdout.txt", stderr_path=log_dir / f"bcftools_dropdups_{label}_stderr.txt")
        if rc != 0:
            raise RuntimeError(f"bcftools view -S failed while dropping duplicates for {vcf}. See logs")
        existing.update(keep)
        return out_vcf, keep

    # rename policy
    new_names = [f"{label}__{s}" if s in existing else s for s in samples]
    names_file = out_dir / f"reheader_names_{label}.txt"
    _write_text(names_file, "\n".join(new_names) + "\n")
    out_vcf = out_dir / f"{vcf.stem}.renamed.vcf.gz"
    rc = _run_cmd([bcftools, "reheader", "-s", str(names_file), "-o", str(out_vcf), str(vcf)], stdout_path=log_dir / f"bcftools_reheader_{label}_stdout.txt", stderr_path=log_dir / f"bcftools_reheader_{label}_stderr.txt")
    if rc != 0:
        raise RuntimeError(f"bcftools reheader failed for {vcf}. See logs")
    # update existing with new names
    existing.update(new_names)
    return out_vcf, new_names


def _genotype_merge(genotype_mode: str, genotype_inputs: List[str], out_dir: Path, *, bcftools: str, plink2: str,
                    merge_strategy: str, dup_sample_policy: str, output_genotype: str, drop_dup_variants: bool,
                    phenotype_ids: Optional[set], final_policy: str, out_prefix_base: str) -> Tuple[Optional[Path], Optional[Path]]:
    import pandas as pd

    if genotype_mode == "none" or not genotype_inputs:
        return None, None

    log_dir = out_dir / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)

    # locate tabix (optional)
    tabix = shutil.which("tabix") or "tabix"

    # Convert inputs to bgz VCF + index
    vcf_dir = out_dir / "_vcf_inputs"
    vcf_dir.mkdir(parents=True, exist_ok=True)

    vcf_paths: List[Path] = []
    in_rows = []

    for i, inp in enumerate(genotype_inputs):
        label = _safe_label(inp)
        inp_s = str(inp)
        pth = Path(inp_s)
        is_vcf = pth.exists() and (pth.name.lower().endswith('.vcf') or pth.name.lower().endswith('.vcf.gz') or pth.name.lower().endswith('.bcf') or pth.name.lower().endswith('.bcf.gz'))
        if genotype_mode == 'plink' or (genotype_mode == 'auto' and not is_vcf):
            # inp is PLINK prefix (or a bed/pgen path stored as prefix)
            vcf = _plink_to_vcf_bgz(plink2, inp_s, vcf_dir, label, log_dir)
        else:
            vcf = Path(inp_s)
        vcf = _ensure_bgz_and_index(bcftools, tabix, vcf, vcf_dir, log_dir)
        samp = _bcftools_list_samples(bcftools, vcf, log_dir)
        in_rows.append({
            "input": inp,
            "label": label,
            "vcf": str(vcf),
            "n_samples": len(samp),
        })
        vcf_paths.append(vcf)

    pd.DataFrame(in_rows).to_csv(out_dir / "genotype_inputs.tsv", sep="\t", index=False)

    # Decide strategy
    sample_lists = [ _bcftools_list_samples(bcftools, v, log_dir) for v in vcf_paths ]
    all_same_samples = True
    if sample_lists:
        ref = sample_lists[0]
        for s in sample_lists[1:]:
            if s != ref:
                all_same_samples = False
                break

    if merge_strategy == "auto":
        strat = "concat_variants" if all_same_samples else "merge_samples"
    else:
        strat = merge_strategy

    _write_text(out_dir / "genotype_merge_strategy.txt", strat + "\n")

    # Handle duplicate samples if merging by samples
    vcf_paths2: List[Path] = []
    existing = set()
    renamed_sample_lists = []

    for v, inp in zip(vcf_paths, genotype_inputs):
        label = _safe_label(inp)
        if strat == "merge_samples":
            v2, s2 = _rename_vcf_samples_if_needed(bcftools, v, label, existing, dup_sample_policy, vcf_dir, log_dir)
            vcf_paths2.append(v2)
            renamed_sample_lists.append(s2)
        else:
            # concat requires identical sample sets; still track
            s = _bcftools_list_samples(bcftools, v, log_dir)
            existing.update(s)
            vcf_paths2.append(v)
            renamed_sample_lists.append(s)

    # Merge/concat
    merged_vcf = out_dir / f"{out_prefix_base}.vcf.gz"
    if strat == "concat_variants":
        cmd = [bcftools, "concat", "-a", "-Oz", "-o", str(merged_vcf)] + [str(v) for v in vcf_paths2]
        rc = _run_cmd(cmd, stdout_path=log_dir / "bcftools_concat_stdout.txt", stderr_path=log_dir / "bcftools_concat_stderr.txt")
        if rc != 0:
            raise RuntimeError(f"bcftools concat failed. See logs/bcftools_concat_stderr.txt")
        # Sort
        sorted_vcf = out_dir / f"{out_prefix_base}.sorted.vcf.gz"
        rc = _run_cmd([bcftools, "sort", "-Oz", "-o", str(sorted_vcf), str(merged_vcf)], stdout_path=log_dir / "bcftools_sort_stdout.txt", stderr_path=log_dir / "bcftools_sort_stderr.txt")
        if rc == 0:
            merged_vcf.unlink(missing_ok=True)
            merged_vcf = sorted_vcf
        # drop duplicate variants (optional)
        if drop_dup_variants:
            dedup_vcf = out_dir / f"{out_prefix_base}.dedup.vcf.gz"
            rc = _run_cmd([bcftools, "norm", "-d", "all", "-Oz", "-o", str(dedup_vcf), str(merged_vcf)], stdout_path=log_dir / "bcftools_norm_stdout.txt", stderr_path=log_dir / "bcftools_norm_stderr.txt")
            if rc == 0:
                merged_vcf = dedup_vcf
    else:
        cmd = [bcftools, "merge", "-m", "none", "-Oz", "-o", str(merged_vcf)] + [str(v) for v in vcf_paths2]
        rc = _run_cmd(cmd, stdout_path=log_dir / "bcftools_merge_stdout.txt", stderr_path=log_dir / "bcftools_merge_stderr.txt")
        if rc != 0:
            raise RuntimeError("bcftools merge failed. See logs/bcftools_merge_stderr.txt")

    # Index merged
    rc = _run_cmd([tabix, "-f", "-p", "vcf", str(merged_vcf)], stdout_path=log_dir / "tabix_merged_stdout.txt", stderr_path=log_dir / "tabix_merged_stderr.txt")
    if rc != 0:
        rc2 = _run_cmd([bcftools, "index", "-f", "-t", str(merged_vcf)], stdout_path=log_dir / "bcftools_index_merged_stdout.txt", stderr_path=log_dir / "bcftools_index_merged_stderr.txt")
        if rc2 != 0:
            raise RuntimeError("Failed to index merged VCF. See logs")

    # Report
    merged_samples = _bcftools_list_samples(bcftools, merged_vcf, log_dir)
    rep = {
        "genotype_mode": genotype_mode,
        "n_inputs": len(genotype_inputs),
        "strategy": strat,
        "dup_sample_policy": dup_sample_policy,
        "n_samples_merged": len(merged_samples),
    }
    pd.DataFrame([rep]).to_csv(out_dir / "genotype_merge_report.tsv", sep="\t", index=False)

    # Convert to PLINK if requested
    plink_out = None
    if output_genotype in ("plink", "both"):
        outp = out_dir / out_prefix_base
        cmd = [plink2, "--vcf", str(merged_vcf), "--make-bed", "--out", str(outp)]
        rc = _run_cmd(cmd, stdout_path=log_dir / "plink2_makebed_stdout.txt", stderr_path=log_dir / "plink2_makebed_stderr.txt")
        if rc != 0:
            raise RuntimeError("plink2 --vcf --make-bed failed. See logs/plink2_makebed_stderr.txt")
        plink_out = outp.with_suffix(".bed")

    # Final sync to phenotype
    final_vcf = None
    if final_policy == "intersection" and phenotype_ids is not None:
        # intersect
        inter = [s for s in merged_samples if s in phenotype_ids]
        _write_text(out_dir / "final_samples_intersection.txt", "\n".join(inter) + "\n")
        final_vcf = out_dir / f"{out_prefix_base}.final.vcf.gz"
        rc = _run_cmd([bcftools, "view", "-S", str(out_dir / "final_samples_intersection.txt"), "-Oz", "-o", str(final_vcf), str(merged_vcf)], stdout_path=log_dir / "bcftools_view_final_stdout.txt", stderr_path=log_dir / "bcftools_view_final_stderr.txt")
        if rc != 0:
            raise RuntimeError("bcftools view -S (final intersection) failed. See logs")
        rc = _run_cmd([tabix, "-f", "-p", "vcf", str(final_vcf)], stdout_path=log_dir / "tabix_final_stdout.txt", stderr_path=log_dir / "tabix_final_stderr.txt")
        if rc != 0:
            _run_cmd([bcftools, "index", "-f", "-t", str(final_vcf)], stdout_path=log_dir / "bcftools_index_final_stdout.txt", stderr_path=log_dir / "bcftools_index_final_stderr.txt")

        # replace merged_vcf for downstream reference

    return merged_vcf, final_vcf


# ----------------- main -----------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--params", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    params = json.loads(Path(args.params).read_text(encoding="utf-8"))
    out_dir = Path(args.out).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    # Allow GUI to override out_dir
    out_override = (params.get("out_dir_override") or "").strip()
    if out_override:
        out_dir = Path(out_override).expanduser().resolve()
        out_dir.mkdir(parents=True, exist_ok=True)

    # phenotype
    ph_files = params.get("phenotype_files") or []
    sheet = (params.get("phenotype_sheet") or "").strip()
    id_col = (params.get("phenotype_id_col") or "").strip()
    merge_mode = (params.get("phenotype_merge_mode") or "auto").strip().lower()
    join = (params.get("phenotype_join") or "outer").strip().lower()
    dup_pol = (params.get("phenotype_dup_policy") or "mean_numeric").strip().lower()
    trait_col = (params.get("trait_collision") or "prefix").strip().lower()

    ph_path, _ph_rep = _phenotype_merge(ph_files, out_dir, sheet=sheet, id_col=id_col, merge_mode=merge_mode, join=join,
                                        dup_policy=dup_pol, trait_collision=trait_col)

    phenotype_ids = None
    if ph_path and ph_path.exists():
        try:
            import pandas as pd
            df = pd.read_csv(ph_path, sep="\t", dtype=str)
            if "id" in df.columns:
                phenotype_ids = set(df["id"].astype(str).tolist())
            else:
                phenotype_ids = set(df.iloc[:, 0].astype(str).tolist())
        except Exception:
            phenotype_ids = None

    # genotype
    merge_genotype = bool(params.get('merge_genotype', False))
    geno_inputs = params.get('genotype_inputs') or []
    geno_mode = (params.get('genotype_mode') or 'none').strip().lower()
    if merge_genotype and geno_inputs:
        geno_mode = 'auto' if geno_mode in ('', 'none') else geno_mode
    if not merge_genotype:
        geno_mode = 'none'
        geno_inputs = []
    bcftools = (params.get("bcftools_bin") or "bcftools").strip() or "bcftools"
    plink2 = (params.get("plink2_bin") or "plink2").strip() or "plink2"
    strat = (params.get("genotype_merge_strategy") or "auto").strip().lower()
    dup_samp_pol = (params.get("dup_sample_policy") or "rename").strip().lower()
    out_geno = (params.get('output_genotype') or ('plink' if geno_mode != 'none' else 'vcf')).strip().lower()
    out_prefix_base = re.sub(r'[^A-Za-z0-9_.-]+', '_', (params.get('genotype_out_prefix') or 'genotype_merged'))
    drop_dup_var = bool(params.get("drop_duplicate_variants", False))

    final_pol = (params.get("final_sample_policy") or "intersection").strip().lower()

    merged_vcf, final_vcf = _genotype_merge(geno_mode, geno_inputs, out_dir, bcftools=bcftools, plink2=plink2,
                                            merge_strategy=strat, dup_sample_policy=dup_samp_pol, output_genotype=out_geno,
                                            drop_dup_variants=drop_dup_var, phenotype_ids=phenotype_ids, final_policy=final_pol, out_prefix_base=out_prefix_base)

    # final phenotype sync if requested
    if final_pol == "intersection" and ph_path and merged_vcf:
        try:
            import pandas as pd
            # Determine genotype samples
            log_dir = out_dir / "logs"
            g_samples = _bcftools_list_samples(bcftools, final_vcf if final_vcf else merged_vcf, log_dir)
            gset = set(g_samples)

            df = pd.read_csv(ph_path, sep="\t", dtype=str)
            idc = "id" if "id" in df.columns else df.columns[0]
            df2 = df[df[idc].astype(str).isin(gset)].copy()
            df2.to_csv(out_dir / "phenotype_final.tsv", sep="\t", index=False)

            rep = {
                "phenotype_ids": int(len(set(df[idc].astype(str).tolist()))),
                "genotype_samples": int(len(gset)),
                "intersection": int(df2[idc].nunique()),
            }
            pd.DataFrame([rep]).to_csv(out_dir / "final_sync_report.tsv", sep="\t", index=False)
        except Exception as e:
            _write_text(out_dir / "final_sync_error.txt", str(e) + "\n")

    # minimal sentinel
    _write_text(out_dir / "DONE.txt", "OK\n")


if __name__ == "__main__":
    main()
