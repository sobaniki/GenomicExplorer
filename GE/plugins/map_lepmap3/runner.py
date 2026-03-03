#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from __future__ import annotations

"""Lep-MAP3 linkage map scaffold runner.

Design goals (matching the project philosophy):
  * Keep GUI runnable even if Lep-MAP3 isn't installed yet.
  * When Lep-MAP3 is installed (jar + Java), run a standard pipeline:
      ParentCall2 -> Filtering2 -> SeparateChromosomes2 -> OrderMarkers2
    while capturing step logs and intermediate files.
  * Always write the usual artifacts so the ResultsPane can load something:
      map_markers.tsv, map_lengths.tsv, map_lengths.png, artifacts.json

This runner intentionally avoids overfitting to a single Lep-MAP3 invocation style.
The user can pass extra arguments per step, and can disable steps.
"""

import argparse
import json
import re
import shutil
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd
import numpy as np


def _norm_suffixes(p: Path) -> str:
    """Return a normalized extension string for format inference.

    IMPORTANT:
      Users often have filenames like "geno.merge.vcf" or "geno.final.tsv".
      Path.suffixes would return [".merge", ".vcf"], which would *not* match
      ".vcf" if we simply join suffixes. Therefore we infer by *ending*.

    Examples:
      a.vcf.gz        -> ".vcf.gz"
      geno.merge.vcf  -> ".vcf"
      geno.final.tsv  -> ".tsv"
      x.pvar.zst      -> ".pvar.zst"
    """
    n = p.name.lower()
    known = [
        ".vcf.gz",
        ".vcf",
        ".bed",
        ".bim",
        ".fam",
        ".pgen",
        ".pvar.zst",
        ".psam.zst",
        ".pvar",
        ".psam",
        ".tsv",
        ".csv",
        ".txt",
    ]
    for ext in known:
        if n.endswith(ext):
            return ext
    return p.suffix.lower()


def _ensure_vcf_from_genotype(
    *,
    geno_path: Path,
    work_dir: Path,
    out_prefix: str,
    threads: int = 1,
) -> Tuple[Optional[Path], List[str]]:
    """Return a VCF path to be used by Lep-MAP3.

    Lep-MAP3 reads genotypes from `vcfFile` in ParentCall2, so for non-VCF inputs we
    convert to a temporary VCF in `work_dir`.

    Supported inputs:
      - VCF / VCF.GZ: used as-is
      - PLINK 1 (.bed/.bim/.fam): plink2 --bfile ... --export vcf
      - PLINK 2 (.pgen/.pvar/.psam): plink2 --pfile ... --export vcf
      - TSV/CSV genotype matrix: converted to a minimal VCF (GT only)
    """
    notes: List[str] = []
    if not geno_path or not geno_path.exists():
        return None, [f"Genotype file not found: {geno_path}"]

    suf = _norm_suffixes(geno_path)

    # VCF
    if suf in {".vcf", ".vcf.gz"}:
        return geno_path, [f"Using VCF input: {geno_path}"]

    # PLINK1
    if suf in {".bed", ".bim", ".fam"}:
        prefix = geno_path.with_suffix("")
        bed, bim, fam = prefix.with_suffix(".bed"), prefix.with_suffix(".bim"), prefix.with_suffix(".fam")
        missing = [str(x) for x in (bed, bim, fam) if not x.exists()]
        if missing:
            return None, [
                "PLINK input requires .bed/.bim/.fam with the same prefix.",
                f"Missing: {', '.join(missing)}",
            ]
        vcf_out_prefix = work_dir / f"{out_prefix}.from_plink1"
        vcf_path = vcf_out_prefix.with_suffix(".vcf")
        ok, n2 = _plink2_export_vcf(
            plink2_exe="plink2",
            kind="bfile",
            prefix=prefix,
            out_prefix=vcf_out_prefix,
            threads=threads,
        )
        notes.extend(n2)
        if not ok or not vcf_path.exists():
            return None, notes + ["PLINK1 -> VCF conversion failed."]
        return vcf_path, notes

    # PLINK2
    if suf in {".pgen", ".pvar", ".psam", ".pvar.zst", ".psam.zst"}:
        # Accept selecting any component file; derive prefix.
        # For .pvar.zst, Path.with_suffix drops only last suffix, so do manual.
        if suf.endswith(".pvar.zst"):
            prefix = Path(str(geno_path)[: -len(".pvar.zst")])
        elif suf.endswith(".psam.zst"):
            prefix = Path(str(geno_path)[: -len(".psam.zst")])
        else:
            prefix = geno_path.with_suffix("")

        pgen = prefix.with_suffix(".pgen")
        pvar = prefix.with_suffix(".pvar")
        psam = prefix.with_suffix(".psam")
        # Some installs use compressed .pvar.zst/.psam.zst
        if not pvar.exists() and prefix.with_suffix(".pvar.zst").exists():
            pvar = prefix.with_suffix(".pvar.zst")
        if not psam.exists() and prefix.with_suffix(".psam.zst").exists():
            psam = prefix.with_suffix(".psam.zst")

        missing = [str(x) for x in (pgen, pvar, psam) if not x.exists()]
        if missing:
            return None, [
                "PLINK2 input requires .pgen/.pvar/.psam with the same prefix.",
                f"Missing: {', '.join(missing)}",
            ]
        vcf_out_prefix = work_dir / f"{out_prefix}.from_plink2"
        vcf_path = vcf_out_prefix.with_suffix(".vcf")
        ok, n2 = _plink2_export_vcf(
            plink2_exe="plink2",
            kind="pfile",
            prefix=prefix,
            out_prefix=vcf_out_prefix,
            threads=threads,
        )
        notes.extend(n2)
        if not ok or not vcf_path.exists():
            return None, notes + ["PLINK2 -> VCF conversion failed."]
        return vcf_path, notes

    # TSV/CSV genotype matrix
    if suf in {".tsv", ".csv", ".txt"}:
        vcf_path = work_dir / f"{out_prefix}.from_table.vcf"

        marker_map_path = None
        # If a marker_map.tsv (marker/chr/pos) exists next to the genotype table, use it
        # to populate CHROM/POS in the generated VCF (avoids dummy sequential POS).
        cand = geno_path.parent / "marker_map.tsv"
        if cand.exists():
            marker_map_path = cand

        ok, n2 = _table_to_minimal_vcf(geno_path=geno_path, out_vcf=vcf_path, marker_map_path=marker_map_path)
        notes.extend(n2)
        if not ok or not vcf_path.exists():
            return None, notes + ["Table -> VCF conversion failed."]
        return vcf_path, notes
    return None, [
        f"Unsupported genotype file extension: {geno_path.name}",
        "Supported: .vcf/.vcf.gz, PLINK (.bed/.bim/.fam), PLINK2 (.pgen/.pvar/.psam), .tsv/.csv.",
    ]


def _plink2_export_vcf(
    *,
    plink2_exe: str,
    kind: str,
    prefix: Path,
    out_prefix: Path,
    threads: int = 1,
) -> Tuple[bool, List[str]]:
    notes: List[str] = []
    if shutil.which(plink2_exe) is None:
        return False, [
            f"plink2 not found in PATH: {plink2_exe}",
            "Install plink2 (conda/micromamba) or ensure it is on PATH.",
        ]

    cmd = [plink2_exe]
    if kind == "bfile":
        cmd += ["--bfile", str(prefix)]
    else:
        cmd += ["--pfile", str(prefix)]

    # Export to VCF (hard calls). Keep it simple and widely compatible.
    cmd += ["--export", "vcf", "--out", str(out_prefix), "--threads", str(max(1, int(threads)))]
    notes.append("[plink2] " + " ".join(cmd))
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True)
    except Exception as e:
        return False, notes + [f"Failed to run plink2: {e}"]

    if proc.stdout:
        notes.append(proc.stdout.strip()[:2000])
    if proc.stderr:
        notes.append(proc.stderr.strip()[:2000])
    if proc.returncode != 0:
        return False, notes + [f"plink2 exited with code {proc.returncode}"]
    return True, notes


def _gt_from_value(v) -> str:
    if v is None:
        return "./."
    s = str(v).strip()
    if s == "" or s in {".", "./.", "NA", "NaN", "nan", "None", "-9"}:
        return "./."
    # pass-through GT like 0/1, 1|1, etc.
    if "/" in s or "|" in s:
        # basic validation; otherwise keep as missing
        if any(t in s for t in ["0", "1", "."]):
            return s.replace("|", "/")
        return "./."

    up = s.upper()
    # common allele encodings
    if up in {"AA", "A"}:
        return "0/0"
    if up in {"AB", "BA", "H", "HET", "HETER", "HETERO", "HETEROZYGOUS"}:
        return "0/1"
    if up in {"BB", "B"}:
        return "1/1"

    # numeric dosage 0/1/2
    try:
        x = float(s)
        if x != x:  # NaN
            return "./."
        xi = int(round(x))
        if xi == 0:
            return "0/0"
        if xi == 1:
            return "0/1"
        if xi == 2:
            return "1/1"
        return "./."
    except Exception:
        return "./."



def _load_marker_map(marker_map_path: Path) -> Tuple[Dict[str, Tuple[str, int, str]], List[str]]:
    """Load marker_map.tsv (marker/chr/pos) and return mapping for VCF CHROM/POS.

    Returns:
      mapping: marker -> (chrom, pos_int, pos_raw_str)
    """
    notes: List[str] = []
    if marker_map_path is None or (not Path(marker_map_path).exists()):
        return {}, notes
    try:
        df = _read_table(Path(marker_map_path))
    except Exception as e:
        return {}, [f"Failed to read marker_map: {marker_map_path} ({e})"]
    if df.shape[0] == 0 or df.shape[1] < 2:
        return {}, [f"marker_map is empty or invalid: {marker_map_path}"]

    mcol = _infer_col(df, ["marker", "id", "snp", "locus", "variant"], 0)
    ccol = _infer_col(df, ["chr", "chrom", "chromosome"], 1)
    pcol = _infer_col(df, ["pos", "position", "bp", "cm", "cM"], 2)

    pos_num = pd.to_numeric(df[pcol], errors="coerce")
    scale = 1
    if pos_num.notna().any():
        nonna = pos_num.dropna().astype(float).values
        has_frac = bool((abs(nonna - np.round(nonna)) > 1e-6).any())
        # If positions are fractional, likely cM. Scale to preserve ordering and avoid too many duplicates.
        if has_frac:
            scale = 1_000_000

    mapping: Dict[str, Tuple[str, int, str]] = {}
    used: Dict[str, set] = {}
    for _, row in df.iterrows():
        mk = str(row[mcol]).strip()
        if mk == "" or mk.lower() in {"nan", "na", "none"}:
            continue

        chrom_raw = row[ccol]
        chrom = str(chrom_raw).strip()
        # normalize "1.0" -> "1"
        if re.fullmatch(r"[0-9]+\.0", chrom):
            try:
                chrom = str(int(float(chrom)))
            except Exception:
                pass

        pos_raw = row[pcol]
        pos = pd.to_numeric(pos_raw, errors="coerce")
        if pd.isna(pos):
            continue
        pos_i = int(round(float(pos) * scale))
        if pos_i < 1:
            pos_i = 1

        s = used.setdefault(chrom, set())
        while pos_i in s:
            pos_i += 1
        s.add(pos_i)
        mapping[mk] = (chrom, pos_i, str(pos_raw))

    notes.append(f"Using marker_map for CHROM/POS: {marker_map_path} (rows={df.shape[0]}, mapped={len(mapping)}, scale={scale})")
    return mapping, notes


def _table_to_minimal_vcf(*, geno_path: Path, out_vcf: Path, marker_map_path: Optional[Path] = None) -> Tuple[bool, List[str]]:
    """Convert a genotype matrix TSV/CSV into a minimal VCF with GT only.

    Supported layouts (auto-detected):
      A) sample rows: first column is sample id, remaining columns are markers
      B) marker rows: first column is marker id, remaining columns are samples
    """
    notes: List[str] = []
    try:
        df = _read_table(geno_path)
    except Exception as e:
        return False, [f"Failed to read genotype table: {e}"]
    if df.shape[0] == 0 or df.shape[1] < 2:
        return False, ["Genotype table is empty or has too few columns."]

    cols = [str(c) for c in df.columns.tolist()]
    first = cols[0].lower().strip()

    id_like = {"id", "iid", "sample", "samples", "ind", "individual", "taxa", "genotype", "line"}

    # Heuristic: if the first column header looks like a sample-id column, assume sample rows.
    # Otherwise assume marker rows.
    sample_rows = first in id_like

    # If there is no meaningful header (e.g. pandas auto "Unnamed: 0"), peek at the first cell.
    # Best-effort only; when unsure, we keep marker-rows.
    if not sample_rows:
        try:
            first_cell = str(df.iloc[0, 0]).lower().strip()
            if first_cell in id_like:
                sample_rows = True
        except Exception:
            pass

    out_vcf.parent.mkdir(parents=True, exist_ok=True)

    # Load optional marker_map.tsv to populate CHROM/POS (avoid dummy sequential POS)
    marker_map: Dict[str, Tuple[str, int, str]] = {}
    if marker_map_path is not None and Path(marker_map_path).exists():
        marker_map, mm_notes = _load_marker_map(Path(marker_map_path))
        notes.extend(mm_notes)

    n_markers: int = 0

    # Prepare header + body
    with open(out_vcf, "wt", encoding="utf-8") as w:
        w.write("##fileformat=VCFv4.2\n")
        w.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>\n")
        if marker_map:
            w.write("##INFO=<ID=GE_POS,Number=1,Type=String,Description=Original position from marker_map.tsv>\n")

        # Track occupied (CHROM, POS) from marker_map so that fallback sequential positions do not collide.
        used_pos: Dict[str, set] = {}
        if marker_map:
            for (ch, pi, _) in marker_map.values():
                used_pos.setdefault(str(ch), set()).add(int(pi))
        fallback_next: Dict[str, int] = {"1": 1}

        if sample_rows:
            id_col = df.columns[0]
            samples = [str(x) for x in df[id_col].tolist()]
            marker_cols = [c for c in df.columns.tolist()[1:]]
            w.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n")

            for m in marker_cols:
                mname = str(m)
                gts = [_gt_from_value(v) for v in df[m].tolist()]
                if marker_map and mname in marker_map:
                    chrom, pos_i, pos_raw = marker_map[mname]
                    info = f"GE_POS={pos_raw}"
                else:
                    chrom = "1"
                    pos_i = int(fallback_next.get(chrom, 1))
                    s = used_pos.setdefault(chrom, set())
                    while pos_i in s:
                        pos_i += 1
                    s.add(pos_i)
                    fallback_next[chrom] = pos_i + 1
                    info = "."
                w.write(f"{chrom}\t{pos_i}\t{mname}\tA\tC\t.\tPASS\t{info}\tGT\t" + "\t".join(gts) + "\n")
            n_markers = len(marker_cols)

        else:
            marker_col = df.columns[0]
            samples = [str(c) for c in df.columns.tolist()[1:]]
            w.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n")

            for _, row in df.iterrows():
                mname = str(row[marker_col])
                vals = [row[c] for c in df.columns.tolist()[1:]]
                gts = [_gt_from_value(v) for v in vals]
                if marker_map and mname in marker_map:
                    chrom, pos_i, pos_raw = marker_map[mname]
                    info = f"GE_POS={pos_raw}"
                else:
                    chrom = "1"
                    pos_i = int(fallback_next.get(chrom, 1))
                    s = used_pos.setdefault(chrom, set())
                    while pos_i in s:
                        pos_i += 1
                    s.add(pos_i)
                    fallback_next[chrom] = pos_i + 1
                    info = "."
                w.write(f"{chrom}\t{pos_i}\t{mname}\tA\tC\t.\tPASS\t{info}\tGT\t" + "\t".join(gts) + "\n")
            n_markers = int(df.shape[0])

    notes.append(f"Converted table to VCF: {out_vcf} (markers={int(n_markers)})")
    return True, notes


def _detect_sep(path: Path) -> str:
    """Best-effort delimiter detection for small TSV/CSV files."""
    try:
        line = path.read_text(encoding="utf-8", errors="ignore").splitlines()[0]
    except Exception:
        return "\t"
    if "\t" in line:
        return "\t"
    if "," in line:
        return ","
    return None  # whitespace


def _read_table(path: Path) -> pd.DataFrame:
    sep = _detect_sep(path)
    if sep is None:
        return pd.read_csv(path, sep=r"\s+", engine="python")
    return pd.read_csv(path, sep=sep)


def _read_vcf_samples(vcf_path: Path) -> list[str]:
    """Read sample names from a VCF/VCF.GZ header (#CHROM line)."""
    import gzip
    opener = gzip.open if vcf_path.suffix.lower() == ".gz" else open
    with opener(vcf_path, "rt", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("#CHROM"):
                parts = line.rstrip("\n").split("\t")
                if len(parts) > 9:
                    return parts[9:]
                return []
    return []


def _infer_col(df: pd.DataFrame, candidates: list[str], fallback_idx: int) -> str:
    cols = [c.strip() for c in df.columns.tolist()]
    lower_map = {c.lower(): c for c in cols}
    for cand in candidates:
        if cand.lower() in lower_map:
            return lower_map[cand.lower()]
    # fallback by index
    if 0 <= fallback_idx < len(cols):
        return cols[fallback_idx]
    return cols[0]


def _build_pedigree_from_cross_files(
    *,
    out_path: Path,
    vcf_path: Path,
    cross_ind_path: Path,
    par_per_cross_path: Path,
    pedigree_model: str = "founder_parents",
) -> tuple[bool, list[str]]:
    """Create Lep-MAP3 pedigree file (6 lines, n+2 tab cols) from cross_ind + par_per_cross.

    ParentCall2 can take genotypes from `vcfFile=` and pedigree from `data=` (this file).
    In this mode, the pedigree file may contain "dummy" parents (no genotype data needed).

    IMPORTANT: Lep-MAP3 expects each family to have parent individuals in the pedigree.
    Therefore we explicitly add (parent1, parent2) columns for each family (and keep
    offspring columns after them). This matches the examples used by Lep-MAP3 authors.

    `pedigree_model`:
      - founder_parents: offspring have father=parent1, mother=parent2
      - dummy_f1_self: offspring have father=F1_<family>, mother=F1_<family> (F1 is dummy)

    Notes:
      - Shared parents across families can be utilised by ParentCall2 with halfSibs=1.
    """
    notes: list[str] = []

    # Robust guard: GUI may pass "." (current directory) when no file is selected.
    # Treat non-files as missing and let caller fall back to manual pedigree.
    if (not cross_ind_path) or (not Path(cross_ind_path).is_file()):
        return False, [f"cross_ind.tsv is not a file: {cross_ind_path}"]
    if (not par_per_cross_path) or (not Path(par_per_cross_path).is_file()):
        return False, [f"par_per_cross.tsv is not a file: {par_per_cross_path}"]

    samples = _read_vcf_samples(vcf_path)
    if not samples:
        return False, [f"Failed to read VCF samples from: {vcf_path}"]

    df_ci = _read_table(cross_ind_path)
    df_pp = _read_table(par_per_cross_path)

    # infer columns
    ci_id = _infer_col(df_ci, ["id", "sample", "ind", "individual"], 0)
    ci_cross = _infer_col(df_ci, ["cross", "family", "pop", "population"], 1)

    pp_cross = _infer_col(df_pp, ["cross", "family", "pop", "population"], 0)
    pp_p1 = _infer_col(df_pp, ["parent1", "p1", "father", "male", "parent_a", "par1"], 1)
    pp_p2 = _infer_col(df_pp, ["parent2", "p2", "mother", "female", "parent_b", "par2"], 2)

    id2cross: dict[str, str] = {}
    for _, r in df_ci.iterrows():
        sid = str(r[ci_id]).strip()
        if sid == "" or sid.lower() == "nan":
            continue
        fam = str(r[ci_cross]).strip()
        if fam == "" or fam.lower() == "nan":
            continue
        id2cross[sid] = fam

    cross2par: dict[str, tuple[str, str]] = {}
    for _, r in df_pp.iterrows():
        fam = str(r[pp_cross]).strip()
        if fam == "" or fam.lower() == "nan":
            continue
        p1 = str(r[pp_p1]).strip()
        p2 = str(r[pp_p2]).strip()
        # keep blanks as "0" for now; we will replace per-family with dummy IDs if needed
        if p1 == "" or p1.lower() == "nan":
            p1 = "0"
        if p2 == "" or p2.lower() == "nan":
            p2 = "0"
        cross2par[fam] = (p1, p2)

    if not id2cross:
        notes.append("cross_ind.tsv parsed, but no (id -> cross) rows were found. Falling back to FAM1 for all VCF samples.")

    # assign samples to families (preserve VCF order within family)
    fam2offspring: dict[str, list[str]] = {}
    for sid in samples:
        fam = id2cross.get(sid, "FAM1")
        fam2offspring.setdefault(fam, []).append(sid)

    # ensure every family has a parent definition (or dummy later)
    for fam in list(fam2offspring.keys()):
        if fam not in cross2par:
            cross2par[fam] = ("0", "0")
            notes.append(f"Family {fam} not found in par_per_cross.tsv; using dummy parents for this family.")

    # shared parent detection (for halfSibs suggestion)
    parent_counts: dict[str, int] = {}
    for fam, (p1, p2) in cross2par.items():
        for p in (p1, p2):
            if p and p != "0":
                parent_counts[p] = parent_counts.get(p, 0) + 1
    shared_parents = [p for p, n in parent_counts.items() if n > 1]
    if shared_parents:
        notes.append(
            "Shared parent(s) across families detected (consider ParentCall2 halfSibs=1): "
            + (", ".join(shared_parents[:8]) + ("..." if len(shared_parents) > 8 else ""))
        )

    fam_line = ["CHR", "POS"]
    ind_line = ["CHR", "POS"]
    father_line = ["CHR", "POS"]
    mother_line = ["CHR", "POS"]
    sex_line = ["CHR", "POS"]
    pheno_line = ["CHR", "POS"]

    # build columns grouped by family: [parent1, parent2, offspring...]
    # Note: parent IDs can be dummy and do not need genotype data in VCF.
    for fam in sorted(fam2offspring.keys()):
        p1, p2 = cross2par.get(fam, ("0", "0"))

        if pedigree_model == "dummy_f1_self":
            # one dummy parent (F1) selfed; this is a heuristic for RIL-like designs
            f1 = f"F1_{fam}"
            # add F1 as the only parent column
            fam_line.append(fam)
            ind_line.append(f1)
            father_line.append("0")
            mother_line.append("0")
            sex_line.append("0")
            pheno_line.append("0")

            dad_for_off = f1
            mom_for_off = f1
        else:
            # founder parents model
            if not p1 or p1 == "0":
                p1 = f"DAD_{fam}"
                notes.append(f"Family {fam}: parent1 missing -> using dummy parent {p1}")
            if not p2 or p2 == "0":
                p2 = f"MOM_{fam}"
                notes.append(f"Family {fam}: parent2 missing -> using dummy parent {p2}")

            # add parents as columns in this family (even if repeated across families)
            fam_line.append(fam)
            ind_line.append(p1)
            father_line.append("0")
            mother_line.append("0")
            sex_line.append("1")
            pheno_line.append("0")

            fam_line.append(fam)
            ind_line.append(p2)
            father_line.append("0")
            mother_line.append("0")
            sex_line.append("2")
            pheno_line.append("0")

            dad_for_off = p1
            mom_for_off = p2

        # offspring columns
        for sid in fam2offspring[fam]:
            fam_line.append(fam)
            ind_line.append(sid)
            father_line.append(dad_for_off)
            mother_line.append(mom_for_off)
            sex_line.append("0")
            pheno_line.append("0")

    # Write file (tab-separated)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as f:
        f.write("\t".join(fam_line) + "\n")
        f.write("\t".join(ind_line) + "\n")
        f.write("\t".join(father_line) + "\n")
        f.write("\t".join(mother_line) + "\n")
        f.write("\t".join(sex_line) + "\n")
        f.write("\t".join(pheno_line) + "\n")

    notes.append(f"Pedigree columns written: {len(ind_line)-2} (including parents/dummies)")
    return True, notes


def _parse_ordermarkers2(path: Path) -> pd.DataFrame:
    """Parse OrderMarkers2 stdout and return a marker table.

    The OrderMarkers2 output is a plain text file with blocks like:
      #*** LG = 1 likelihood = ...
      #marker_number\tmale_position\tfemale_position ...
      498\t0.0\t0.0 ...

    We parse only the first 3 columns and keep linkage group id.
    """
    rows = []
    cur_lg: str | None = None
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line:
                continue
            line = line.rstrip("\n")
            if not line:
                continue

            if line.startswith("#*** LG"):
                # e.g. "#*** LG = 1 likelihood = -2472.9"
                try:
                    after = line.split("=", 1)[1].strip()
                    cur_lg = after.split()[0].strip()
                except Exception:
                    cur_lg = None
                continue

            if line.startswith("#"):
                continue

            parts = line.split("\t")
            if len(parts) < 3:
                # some Lep-MAP3 outputs may be space separated
                parts = line.split()
            if len(parts) < 3:
                continue

            try:
                mnum = int(parts[0])
                mpos = float(parts[1])
                fpos = float(parts[2])
            except Exception:
                continue

            rows.append({"lg": str(cur_lg or "NA"), "marker_number": mnum, "male_pos": mpos, "female_pos": fpos})

    df = pd.DataFrame(rows)
    if df.empty:
        return df

    # sexAveraged=1 is commonly used; still be robust.
    df["pos_cM"] = (df["male_pos"].astype(float) + df["female_pos"].astype(float)) / 2.0
    return df


def _iter_vcf_records(vcf_path: Path):
    """Yield (idx_1based, chrom, pos, id) for each variant line in a VCF/VCF.GZ."""
    import gzip

    opener = gzip.open if vcf_path.suffix.lower() == ".gz" else open
    with opener(vcf_path, "rt", encoding="utf-8", errors="ignore") as f:
        idx = 0
        for line in f:
            if not line or line.startswith("#"):
                continue
            idx += 1
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            chrom = parts[0]
            pos = parts[1]
            vid = parts[2]
            yield idx, chrom, pos, vid


def _map_marker_numbers_to_vcf_ids(vcf_path: Path, marker_numbers: list[int]) -> dict[int, str]:
    """Best-effort mapping from Lep-MAP3 marker_number to VCF ID.

    IMPORTANT: This assumes marker_number corresponds to the 1-based order of variants
    in the provided VCF (excluding header lines). This is true for common Lep-MAP3
    workflows that ingest a VCF via ParentCall2.

    If the assumption is violated (e.g., marker numbers are renumbered after filtering),
    the mapping will be incomplete or incorrect. We therefore:
      * fall back to using marker_number as marker id when mapping is missing
      * write a note file to explain the assumption
    """
    want = sorted({int(x) for x in marker_numbers if int(x) > 0})
    if not want:
        return {}

    want_set = set(want)
    max_idx = want[-1]
    out: dict[int, str] = {}
    seen: dict[str, int] = {}

    for idx, chrom, pos, vid in _iter_vcf_records(vcf_path):
        if idx > max_idx and len(out) == len(want_set):
            break
        if idx not in want_set:
            continue
        name = (vid or "").strip()
        if name in {"", "."}:
            name = f"{chrom}:{pos}"

        # make unique
        n = seen.get(name, 0)
        if n > 0:
            name2 = f"{name}__dup{n+1}"
        else:
            name2 = name
        seen[name] = n + 1

        out[int(idx)] = name2
        if len(out) == len(want_set):
            break
    return out


def _write_map_outputs(
    out_dir: Path,
    *,
    order_out: Path,
    vcf_path: Optional[Path] = None,
) -> Tuple[bool, List[str]]:
    """Parse OrderMarkers2 output and write standardized map outputs.

    Returns (ok, notes)
    """
    notes: List[str] = []

    if not order_out.exists() or order_out.stat().st_size <= 0:
        return False, ["OrderMarkers2 output not found or empty."]

    df = _parse_ordermarkers2(order_out)
    if df.empty:
        return False, ["Failed to parse OrderMarkers2 output."]

    # map marker_number -> marker id (best-effort)
    id_map: dict[int, str] = {}
    if vcf_path is not None and vcf_path.exists():
        try:
            id_map = _map_marker_numbers_to_vcf_ids(vcf_path, df["marker_number"].astype(int).tolist())
            notes.append(
                "marker_number→VCF ID mapping was attempted assuming marker_number is the 1-based variant index in the VCF (non-header lines)."
            )
            notes.append(f"VCF mapping hits: {len(id_map)}/{df.shape[0]} markers")
        except Exception as e:
            notes.append(f"VCF mapping failed: {e}")

    df["marker"] = df["marker_number"].map(lambda x: id_map.get(int(x), str(int(x))))

    # Sort and write
    df = df.sort_values(["lg", "pos_cM", "marker"], kind="mergesort")
    map_markers = df[["lg", "marker", "pos_cM", "marker_number", "male_pos", "female_pos"]].copy()
    map_markers.to_csv(out_dir / "map_markers.tsv", sep="\t", index=False)

    # mppR / common format: marker, chr, pos (cM)
    marker_map = df[["marker", "lg", "pos_cM"]].copy()
    marker_map = marker_map.rename(columns={"lg": "chr", "pos_cM": "pos"})
    marker_map.to_csv(out_dir / "marker_map.tsv", sep="\t", index=False)

    # LG length table
    grp = df.groupby("lg", as_index=False).agg(n_markers=("marker", "size"), length_cM=("pos_cM", "max"))
    grp = grp.sort_values(["lg"], kind="mergesort")
    grp.to_csv(out_dir / "map_lengths.tsv", sep="\t", index=False)

    # basic plot
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        plt.figure()
        plt.bar(range(len(grp)), grp["length_cM"].astype(float))
        plt.xticks(range(len(grp)), grp["lg"].astype(str), rotation=90)
        plt.xlabel("LG")
        plt.ylabel("length (cM)")
        plt.title("Linkage map length")
        plt.tight_layout()
        plt.savefig(out_dir / "map_lengths.png", dpi=150)
        plt.close()
    except Exception as e:
        notes.append(f"Failed to write map_lengths.png: {e}")

    # human-readable notes
    (out_dir / "map_note.txt").write_text("\n".join(notes) + "\n", encoding="utf-8")
    return True, notes


def _log_write(log_path: Path, lines: List[str]) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _ensure_parent(p: Path) -> None:
    p.parent.mkdir(parents=True, exist_ok=True)


def _run_step(
    *,
    java_exe: str,
    jar_path: Path,
    mem_gb: int,
    module: str,
    args_kv: Dict[str, str],
    extra_args: str,
    out_stdout: Path,
    out_stderr: Path,
    cwd: Path,
) -> Tuple[int, List[str]]:
    """Run a Lep-MAP3 module and write stdout/stderr to files."""
    _ensure_parent(out_stdout)
    _ensure_parent(out_stderr)

    cmd: List[str] = [
        java_exe,
        f"-Xmx{int(mem_gb)}g",
        "-cp",
        str(jar_path),
        module,
    ]
    # Lep-MAP3 uses key=value style
    for k, v in args_kv.items():
        if v is None:
            continue
        v = str(v).strip()
        if v == "":
            continue
        cmd.append(f"{k}={v}")

    if extra_args:
        # allow users to pass raw tokens like "removeNonInformative=1" or "lodLimit=10"
        # split by whitespace (simple, consistent with other scaffolds)
        cmd.extend([t for t in extra_args.split() if t.strip()])

    log_lines = [f"[lepmap3] RUN: {' '.join(cmd)}"]
    with out_stdout.open("w", encoding="utf-8", errors="ignore") as fo, out_stderr.open(
        "w", encoding="utf-8", errors="ignore"
    ) as fe:
        try:
            proc = subprocess.run(cmd, cwd=str(cwd), stdout=fo, stderr=fe, check=False)
            return int(proc.returncode), log_lines
        except Exception as e:
            fe.write(str(e) + "\n")
            return 1, log_lines + [f"[lepmap3] EXCEPTION: {e}"]


def _placeholder_outputs(out_dir: Path, msg: str) -> None:
    """Write placeholder outputs so GUI remains functional."""
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "error_message.txt").write_text(msg + "\n", encoding="utf-8")

    # Minimal TSVs
    pd.DataFrame(
        [{"chr": "NA", "marker": "NA", "pos_cM": 0.0, "note": "placeholder"}]
    ).to_csv(out_dir / "map_markers.tsv", sep="\t", index=False)
    pd.DataFrame(
        [{"chr": "NA", "n_markers": 0, "length_cM": 0.0, "note": "placeholder"}]
    ).to_csv(out_dir / "map_lengths.tsv", sep="\t", index=False)

    # Simple plot placeholder (matplotlib only if available; otherwise skip)
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        df = pd.read_csv(out_dir / "map_lengths.tsv", sep="\t")
        plt.figure()
        plt.bar(range(len(df)), df["length_cM"])
        plt.xlabel("chr")
        plt.ylabel("length_cM")
        plt.title("Linkage map length (placeholder)")
        plt.tight_layout()
        plt.savefig(out_dir / "map_lengths.png", dpi=150)
        plt.close()
    except Exception:
        # ignore
        pass


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--params", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    params_path = Path(args.params)
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    log_lines: List[str] = []
    log_lines.append("[map_lepmap3] start")
    log_lines.append(f"[map_lepmap3] params_path={params_path}")
    log_lines.append(f"[map_lepmap3] out_dir={out_dir}")

    params = json.loads(params_path.read_text(encoding="utf-8"))

    jar = Path(str(params.get("lepmap3_jar", "") or "")).expanduser()
    java_exe = str(params.get("java", "java") or "java")
    mem_gb = int(float(params.get("java_mem_gb", 8) or 8))
    threads = int(float(params.get("threads", 4) or 4))
    input_path = Path(str(params.get("input_path", "") or "")).expanduser()
    geno = str(params.get("genotype_path", "") or params.get("vcf_path", "") or "")
    geno_path = Path(geno).expanduser() if geno else Path("")
    cross_ind_path = Path(str(params.get("cross_ind_path", "") or "")).expanduser()
    par_per_cross_path = Path(str(params.get("par_per_cross_path", "") or "")).expanduser()
    pedigree_model = str(params.get("pedigree_model", "founder_parents") or "founder_parents")
    halfsibs_mode = str(params.get("parentcall2_halfsibs", "AUTO") or "AUTO")
    prefix = str(params.get("output_prefix", "lepmap3") or "lepmap3")

    # pipeline toggles
    do_parent = str(params.get("do_parentcall2", "TRUE")).upper().startswith("T")
    do_filter = str(params.get("do_filtering2", "TRUE")).upper().startswith("T")
    do_sep = str(params.get("do_separatechromosomes2", "TRUE")).upper().startswith("T")
    do_order = str(params.get("do_ordermarkers2", "TRUE")).upper().startswith("T")

    # per-step extra args (raw tokens)
    extra_parent = str(params.get("parentcall2_extra", "") or "")
    extra_filter = str(params.get("filtering2_extra", "") or "")
    extra_sep = str(params.get("separatechromosomes2_extra", "") or "")
    extra_order = str(params.get("ordermarkers2_extra", "") or "")

    # Standard args
    # NOTE: Lep-MAP3 uses data=<file>. Many modules also accept "threads".
    #base_kv = {"threads": str(threads)}
    base_kv = {"numThreads": str(threads)}

    work_dir = out_dir / "lepmap3_work"
    work_dir.mkdir(parents=True, exist_ok=True)

    if not jar.exists():
        msg = f"Lep-MAP3 jar not found: {jar}. Set 'lepmap3_jar' in GUI and ensure Java is available."
        log_lines.append("[map_lepmap3] " + msg)
        _log_write(out_dir / "run.log", log_lines)
        _placeholder_outputs(out_dir, msg)
        _write_artifacts(out_dir, default_table="map_markers.tsv", default_plot="map_lengths.png")
        return 0

    # Genotypes are required (Lep-MAP3 reads genotypes from vcfFile in ParentCall2).
    # We accept VCF / PLINK / TSV-CSV, and convert to a temporary VCF when needed.
    vcf_path, conv_notes = _ensure_vcf_from_genotype(
        geno_path=geno_path,
        work_dir=work_dir,
        out_prefix=prefix,
        threads=threads,
    )
    for n in conv_notes:
        log_lines.append("[map_lepmap3] " + n)

    if not vcf_path or not vcf_path.exists():
        msg = f"Genotype file not usable: {geno_path}. Provide VCF/PLINK/TSV (conversion creates VCF for Lep-MAP3)."
        log_lines.append("[map_lepmap3] " + msg)
        _log_write(out_dir / "run.log", log_lines)
        _placeholder_outputs(out_dir, msg)
        _write_artifacts(out_dir, default_table="map_markers.tsv", default_plot="map_lengths.png")
        return 0

    # Auto pedigree from NAM-style files (preferred): cross_ind.tsv + par_per_cross.tsv
    ped_notes: List[str] = []
    if cross_ind_path.is_file() and par_per_cross_path.is_file():
        auto_ped = work_dir / f"{prefix}.pedigree.auto.txt"
        ok_ped, ped_notes = _build_pedigree_from_cross_files(
            out_path=auto_ped,
            vcf_path=vcf_path,
            cross_ind_path=cross_ind_path,
            par_per_cross_path=par_per_cross_path,
            pedigree_model=pedigree_model,
        )
        if not ok_ped:
            msg = "Failed to auto-generate pedigree from cross_ind/par_per_cross:\n" + "\n".join(ped_notes)
            log_lines.append("[map_lepmap3] " + msg)
            _log_write(out_dir / "run.log", log_lines)
            _placeholder_outputs(out_dir, msg)
            _write_artifacts(out_dir, default_table="map_markers.tsv", default_plot="map_lengths.png")
            return 0
        input_path = auto_ped
        log_lines.append(f"[map_lepmap3] auto pedigree generated: {auto_ped}")
        for n in ped_notes:
            log_lines.append("[map_lepmap3] " + n)
    else:
        if not input_path.exists():
            msg = f"Input pedigree file not found: {input_path}. Provide either (cross_ind.tsv + par_per_cross.tsv) or a manual pedigree file."
            log_lines.append("[map_lepmap3] " + msg)
            _log_write(out_dir / "run.log", log_lines)
            _placeholder_outputs(out_dir, msg)
            _write_artifacts(out_dir, default_table="map_markers.tsv", default_plot="map_lengths.png")
            return 0

    # ParentCall2: halfSibs auto/override
    try:
        hs = None
        if str(halfsibs_mode).strip().upper().startswith("A") and par_per_cross_path.is_file():
            df_pp = _read_table(par_per_cross_path)
            pp_cross = _infer_col(df_pp, ["cross", "family", "pop", "population"], 0)
            pp_p1 = _infer_col(df_pp, ["parent1", "p1", "father", "male", "parent_a", "par1"], 1)
            pp_p2 = _infer_col(df_pp, ["parent2", "p2", "mother", "female", "parent_b", "par2"], 2)
            parent_counts: dict[str, int] = {}
            for _, r in df_pp.iterrows():
                for p in [str(r[pp_p1]).strip(), str(r[pp_p2]).strip()]:
                    if p and p != "0" and p.lower() != "nan":
                        parent_counts[p] = parent_counts.get(p, 0) + 1
            hs = 1 if any(n > 1 for n in parent_counts.values()) else 0
        elif str(halfsibs_mode).strip() in {"0", "1"}:
            hs = int(str(halfsibs_mode).strip())

        if hs == 1 and "halfsibs=" not in extra_parent.lower():
            extra_parent = (extra_parent + " halfSibs=1").strip()
            log_lines.append("[map_lepmap3] ParentCall2: appended halfSibs=1")
    except Exception as e:
        log_lines.append(f"[map_lepmap3] halfSibs auto failed: {e}")

    # ensure Java exists in PATH if java_exe is not absolute
    if Path(java_exe).name == java_exe and shutil.which(java_exe) is None:
        msg = f"Java executable not found in PATH: {java_exe}. Install Java or set 'java' path."
        log_lines.append("[map_lepmap3] " + msg)
        _log_write(out_dir / "run.log", log_lines)
        _placeholder_outputs(out_dir, msg)
        _write_artifacts(out_dir, default_table="map_markers.tsv", default_plot="map_lengths.png")
        return 0

    # Pipeline files
    cur_data = input_path
    parent_out = work_dir / f"{prefix}.parentcall2.txt"
    filter_out = work_dir / f"{prefix}.filtering2.txt"
    sep_out = work_dir / f"{prefix}.separatechromosomes2.txt"
    order_out = work_dir / f"{prefix}.ordermarkers2.txt"

    vcf_file = str(vcf_path)

    # (1) ParentCall2
    if do_parent:
        rc, step_log = _run_step(
            java_exe=java_exe,
            jar_path=jar,
            mem_gb=mem_gb,
            module="ParentCall2",
            #args_kv={**base_kv, "data": str(cur_data)},
            #args_kv={"data": str(cur_data), "vcfFile": str(vcf), "removeNonInformative": int(1)},
            args_kv={"data": str(cur_data), "vcfFile": vcf_file},
            extra_args=extra_parent,
            out_stdout=parent_out,
            out_stderr=work_dir / f"{prefix}.parentcall2.stderr.txt",
            cwd=work_dir,
        )
        log_lines.extend(step_log)
        log_lines.append(f"[lepmap3] ParentCall2 rc={rc}")
        if rc != 0:
            log_lines.append("[lepmap3] ParentCall2 failed; continuing with input as-is")
        else:
            cur_data = parent_out

    # (2) Filtering2
    if do_filter:
        rc, step_log = _run_step(
            java_exe=java_exe,
            jar_path=jar,
            mem_gb=mem_gb,
            module="Filtering2",
            #args_kv={**base_kv, "data": str(cur_data)},
            #args_kv={"data": str(cur_data), "dataTolerance": float(0.001)},
            args_kv={"data": str(cur_data)},
            extra_args=extra_filter,
            out_stdout=filter_out,
            out_stderr=work_dir / f"{prefix}.filtering2.stderr.txt",
            cwd=work_dir,
        )
        log_lines.extend(step_log)
        log_lines.append(f"[lepmap3] Filtering2 rc={rc}")
        if rc != 0:
            log_lines.append("[lepmap3] Filtering2 failed; continuing with previous data")
        else:
            cur_data = filter_out

    # (3) SeparateChromosomes2
    map_file: Optional[Path] = None
    if do_sep:
        rc, step_log = _run_step(
            java_exe=java_exe,
            jar_path=jar,
            mem_gb=mem_gb,
            module="SeparateChromosomes2",
            #args_kv={**base_kv, "data": str(cur_data), "lodLimit": int(5), "theta": float(0.05), "distortionLod": int(1)},
            args_kv={**base_kv, "data": str(cur_data)},
            extra_args=extra_sep,
            out_stdout=sep_out,
            out_stderr=work_dir / f"{prefix}.separatechromosomes2.stderr.txt",
            cwd=work_dir,
        )
        log_lines.extend(step_log)
        log_lines.append(f"[lepmap3] SeparateChromosomes2 rc={rc}")
        if rc != 0:
            log_lines.append("[lepmap3] SeparateChromosomes2 failed; OrderMarkers2 may be skipped")
        else:
            map_file = sep_out

    # (4) OrderMarkers2
    if do_order and map_file is not None and map_file.exists():
        rc, step_log = _run_step(
            java_exe=java_exe,
            jar_path=jar,
            mem_gb=mem_gb,
            module="OrderMarkers2",
            #args_kv={**base_kv, "data": str(cur_data), "map": str(map_file), "sexAveraged": int(1)},
            args_kv={**base_kv, "data": str(cur_data), "map": str(map_file)},
            extra_args=extra_order,
            out_stdout=order_out,
            out_stderr=work_dir / f"{prefix}.ordermarkers2.stderr.txt",
            cwd=work_dir,
        )
        log_lines.extend(step_log)
        log_lines.append(f"[lepmap3] OrderMarkers2 rc={rc}")
    elif do_order:
        log_lines.append("[lepmap3] OrderMarkers2 skipped (no map file from SeparateChromosomes2)")

    # Copy key files to out_dir (stable names)
    produced_files: List[str] = []
    for fp in [parent_out, filter_out, sep_out, order_out]:
        if fp.exists() and fp.stat().st_size > 0:
            dest = out_dir / fp.name
            shutil.copyfile(fp, dest)
            produced_files.append(dest.name)
    # Also include stderr per step if exists
    for fp in work_dir.glob(f"{prefix}.*.stderr.txt"):
        if fp.exists() and fp.stat().st_size > 0:
            dest = out_dir / fp.name
            shutil.copyfile(fp, dest)
            produced_files.append(dest.name)

    _log_write(out_dir / "run.log", log_lines)

    # Build standardized outputs for GUI + downstream tools
    order_copy = out_dir / order_out.name
    ok, notes = _write_map_outputs(
        out_dir,
        order_out=order_copy if order_copy.exists() else order_out,
        vcf_path=vcf_path if vcf_path and vcf_path.exists() else None,
    )

    if ok:
        # keep a short status file (GUI doesn't rely on it, but users appreciate it)
        (out_dir / "status.txt").write_text("Lep-MAP3 map parsed successfully.\n", encoding="utf-8")
        produced_files.extend(["map_markers.tsv", "map_lengths.tsv", "marker_map.tsv", "map_lengths.png", "map_note.txt", "status.txt"])
    else:
        msg = "Lep-MAP3 pipeline executed, but map parsing failed.\n" + "\n".join(notes)
        _placeholder_outputs(out_dir, msg)
        produced_files.append("error_message.txt")

    _write_artifacts(
        out_dir,
        default_table="map_markers.tsv",
        default_plot="map_lengths.png",
        extra_files=sorted(set(produced_files + ["run.log"]))
    )
    return 0


def _write_artifacts(
    out_dir: Path,
    default_table: str,
    default_plot: str,
    extra_files: Optional[List[str]] = None,
) -> None:
    art = {
        "table": default_table,
        "plot": default_plot,
        "tables": ["map_markers.tsv", "map_lengths.tsv", "marker_map.tsv"],
        "plots": ["map_lengths.png"],
        "files": extra_files or [],
        "note": "Lep-MAP3 scaffold: best-effort pipeline (ParentCall2/Filtering2/SeparateChromosomes2/OrderMarkers2) with map_markers.tsv + marker_map.tsv outputs.",
    }
    (out_dir / "artifacts.json").write_text(json.dumps(art, indent=2), encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
