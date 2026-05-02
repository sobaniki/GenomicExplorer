#!/usr/bin/env python3

"""Preprocess: Split/Merge/Extract (ID matching & subset generation)

This plugin helps users align complex genotype sample IDs with phenotype spreadsheets
and optionally extract matching samples from VCF/PLINK.

Typical use case (maizegenetics integrated datasets):
- Genotype sample IDs: e.g. "PI233313:250007467"
- Phenotype IDs:      e.g. "PI233313" (or other naming variants)

The plugin supports:
- Key transforms (before/after delimiter, regex)
- String normalization
- Optional alias table
- Phenotype row filtering (e.g., Panel == NAM)
- Duplicate resolution
- Outputs: match reports, phenotype.tsv for downstream GE modules
- Optional extraction: bcftools view -S, plink2 --keep

Inputs (params.json)
--------------------
out_dir_override: str (optional)
output_prefix: str (optional)  # if GUI passes, used only for reporting

# genotype IDs
id_table_tsv: str  # required (TSV with sample ids; can be a VCF header-derived table)
id_col: str (optional, default: 'id')
fid_col: str (optional)  # for PLINK keep files; when empty, use '0'

# phenotype
phenotype_path: str  # required (xlsx/csv/tsv)
phenotype_sheet: str (optional; only for xlsx)
phenotype_id_col: str (optional; default: first column)

# panel-aware phenotype ID (for mixed datasets, e.g., maize NAM)
panel_aware: bool (default False)
panel_col: str (optional; default: '')
panel_values: str (optional; comma-separated; e.g. 'NAM')
panel_id_col: str (optional; e.g. 'Z_Num')

# transforms
geno_transform: str  # none|before_delim|after_delim|regex
pheno_transform: str # none|before_delim|after_delim|regex
delim: str (optional, default ':')
regex: str (optional; used when *_transform == regex)

# normalization
strip: bool (default True)
collapse_spaces: bool (default True)
case: str (keep|lower|upper)
remove_parentheses: bool (default False)
remove_chars_regex: str (optional)  # regex to remove characters

# optional filters
pheno_filter_col: str (optional)
pheno_filter_values: str (optional; comma-separated)

# alias
alias_tsv: str (optional)  # first two columns: source, target

# duplicate resolution
keep_genotype: str  # keep_all|first
phenotype_reduce: str  # keep_all|mean_numeric|first

# outputs
write_phenotype_tsv: bool (default True)
phenotype_out_name: str (default 'phenotype.tsv')
write_sample_list: bool (default True)
sample_list_name: str (default 'selected_samples.txt')
max_pairs_tsv: int (default 200000)  # cap to avoid gigantic TSVs

# extraction (optional)
extract_mode: str  # none|bcftools|plink2
vcf_path: str (optional)
bcftools_bin: str (optional, default 'bcftools')
vcf_out: str (optional)

plink_prefix: str (optional)
plink2_bin: str (optional, default 'plink2')
plink_out_prefix: str (optional)

Outputs
-------
match_summary.tsv
matched_keys.tsv
matched_pairs.tsv (capped)
unmatched_genotype.tsv
unmatched_phenotype.tsv
duplicate_keys.tsv
phenotype.tsv (optional; first col is full genotype ID)
selected_samples.txt (optional; full genotype IDs)
artifacts.json
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import shlex
from dataclasses import dataclass
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


def _as_bool(x, default: bool = False) -> bool:
    if x is None:
        return default
    if isinstance(x, bool):
        return x
    s = str(x).strip().lower()
    if s in {"1", "true", "t", "yes", "y", "on"}:
        return True
    if s in {"0", "false", "f", "no", "n", "off"}:
        return False
    return default


def _run(cmd: Sequence[str], cwd: Optional[Path], stdout_path: Path, stderr_path: Path) -> None:
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


def _run_capture(cmd: Sequence[str], cwd: Optional[Path] = None) -> Tuple[int, str, str]:
    """Run a command and capture stdout/stderr (no exception)."""
    proc = subprocess.run(
        list(cmd),
        cwd=str(cwd) if cwd else None,
        capture_output=True,
        text=True,
        check=False,
    )
    return int(proc.returncode), (proc.stdout or ""), (proc.stderr or "")


def _run_shell(cmdline: str, cwd: Optional[Path], stdout_path: Path, stderr_path: Path) -> None:
    """Run a shell pipeline via bash -lc (captures stdout/stderr)."""
    proc = subprocess.run(
        ["bash", "-lc", cmdline],
        cwd=str(cwd) if cwd else None,
        capture_output=True,
        text=True,
        check=False,
    )
    _write_text(stdout_path, proc.stdout or "")
    _write_text(stderr_path, proc.stderr or "")
    if proc.returncode != 0:
        tail = (proc.stderr or proc.stdout or "")[-2000:]
        raise RuntimeError(f"Command failed (code={proc.returncode}): {cmdline}\n--- tail ---\n{tail}")


def _vcf_has_records(vcf_path: Path, max_bytes: int = 2_000_000) -> bool:
    """Cheap check: does a VCF/VCF.GZ contain at least one non-header record?"""
    try:
        if str(vcf_path).endswith(".gz"):
            import gzip

            with gzip.open(vcf_path, "rt", encoding="utf-8", errors="ignore") as f:
                read_bytes = 0
                for ln in f:
                    read_bytes += len(ln.encode("utf-8", errors="ignore"))
                    if not ln.startswith("#"):
                        return True
                    if read_bytes >= max_bytes:
                        break
        else:
            with open(vcf_path, "rt", encoding="utf-8", errors="ignore") as f:
                read_bytes = 0
                for ln in f:
                    read_bytes += len(ln.encode("utf-8", errors="ignore"))
                    if not ln.startswith("#"):
                        return True
                    if read_bytes >= max_bytes:
                        break
    except Exception:
        return False
    return False


def _read_table_auto(path: Path):
    import pandas as pd

    # Try common separators first for speed and determinism
    for sep in ["\t", ","]:
        try:
            df = pd.read_csv(path, sep=sep, dtype=str)
            if df.shape[1] >= 1:
                return df
        except Exception:
            continue
    # Fallback: let pandas infer
    return pd.read_csv(path, sep=None, engine="python", dtype=str)


def _read_phenotype(path: Path, sheet: Optional[str]):
    import pandas as pd

    if path.suffix.lower() in {".xlsx", ".xls"}:
        if sheet and str(sheet).strip():
            return pd.read_excel(path, sheet_name=str(sheet).strip(), dtype=str)
        return pd.read_excel(path, sheet_name=0, dtype=str)
    return _read_table_auto(path)


def _parse_csv_list(s: str) -> List[str]:
    if not s:
        return []
    out: List[str] = []
    for tok in str(s).split(","):
        t = tok.strip()
        if t:
            out.append(t)
    return out


@dataclass
class NormOpt:
    strip: bool = True
    collapse_spaces: bool = True
    case: str = "keep"  # keep|lower|upper
    remove_parentheses: bool = False
    remove_chars_regex: str = ""


def _normalize(s: str, opt: NormOpt) -> str:
    if s is None:
        return ""
    x = str(s)
    if opt.strip:
        x = x.strip()
    if opt.collapse_spaces:
        x = re.sub(r"\s+", " ", x)
    if opt.remove_parentheses:
        x = re.sub(r"\([^)]*\)", "", x)
        if opt.strip:
            x = x.strip()
        if opt.collapse_spaces:
            x = re.sub(r"\s+", " ", x)
    if opt.remove_chars_regex:
        try:
            x = re.sub(opt.remove_chars_regex, "", x)
        except re.error:
            # ignore invalid regex
            pass
    cm = (opt.case or "keep").strip().lower()
    if cm == "lower":
        x = x.lower()
    elif cm == "upper":
        x = x.upper()
    return x


def _apply_transform(s: str, mode: str, delim: str, rx: str) -> str:
    mode = (mode or "none").strip().lower()
    x = "" if s is None else str(s)
    if mode == "none":
        return x
    if mode == "before_delim":
        if delim and delim in x:
            return x.split(delim, 1)[0]
        return x
    if mode == "after_delim":
        if delim and delim in x:
            return x.split(delim, 1)[1]
        return x
    if mode == "regex":
        if not rx:
            return x
        try:
            m = re.search(rx, x)
        except re.error:
            return x
        if not m:
            return x
        # prefer named group 'id' if present
        if "id" in m.groupdict() and m.group("id") is not None:
            return str(m.group("id"))
        # otherwise group1
        if m.groups():
            return str(m.group(1))
        return x
    return x


def _load_alias_map(path: Path) -> Dict[str, str]:
    if not path or not path.exists():
        return {}
    try:
        import pandas as pd

        df = _read_table_auto(path)
        if df.shape[1] < 2:
            return {}
        src = df.iloc[:, 0].astype(str)
        tgt = df.iloc[:, 1].astype(str)
        mp: Dict[str, str] = {}
        for a, b in zip(src.tolist(), tgt.tolist()):
            a2 = str(a).strip()
            b2 = str(b).strip()
            if a2 and b2:
                mp[a2] = b2
        return mp
    except Exception:
        # fallback: plain text
        mp: Dict[str, str] = {}
        for ln in _read_text(path).splitlines():
            if not ln.strip() or ln.lstrip().startswith("#"):
                continue
            toks = re.split(r"\t|,", ln.strip())
            if len(toks) >= 2:
                mp[toks[0].strip()] = toks[1].strip()
        return mp


def _safe_col(df, name: str) -> Optional[str]:
    if df is None or df.shape[1] < 1:
        return None
    if name and str(name).strip() and str(name).strip() in df.columns:
        return str(name).strip()
    return None


def _first_col(df) -> str:
    return str(df.columns[0])


def _to_numeric_series(s):
    import pandas as pd

    s2 = s.replace({"": None, "NA": None, "NaN": None, "nan": None, ".": None, "-": None})
    return pd.to_numeric(s2, errors="coerce")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--params", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    params_path = Path(args.params)
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    p = json.loads(_read_text(params_path) or "{}")

    out_override = str(p.get("out_dir_override") or "").strip()
    root = Path(out_override).expanduser().resolve() if out_override else out_dir
    root.mkdir(parents=True, exist_ok=True)

    # --- read inputs ---
    id_table_tsv_str = str(p.get("id_table_tsv") or "").strip()
    id_table_tsv = Path(id_table_tsv_str).expanduser() if id_table_tsv_str else None
    phenotype_path = Path(str(p.get("phenotype_path") or "").strip()).expanduser()

    if not phenotype_path.exists():
        raise SystemExit(f"phenotype_path not found: {phenotype_path}")

    geno_id_col_req = str(p.get("id_col") or "id").strip() or "id"
    geno_fid_col_req = str(p.get("fid_col") or "").strip()

    ph_sheet = str(p.get("phenotype_sheet") or "").strip()
    ph_id_col_req = str(p.get("phenotype_id_col") or "").strip()

    # transforms
    geno_transform = str(p.get("geno_transform") or "before_delim").strip() or "before_delim"
    pheno_transform = str(p.get("pheno_transform") or "none").strip() or "none"
    delim = str(p.get("delim") or ":")
    rx = str(p.get("regex") or "").strip()

    # normalization
    opt = NormOpt(
        strip=_as_bool(p.get("strip"), True),
        collapse_spaces=_as_bool(p.get("collapse_spaces"), True),
        case=str(p.get("case") or "keep"),
        remove_parentheses=_as_bool(p.get("remove_parentheses"), False),
        remove_chars_regex=str(p.get("remove_chars_regex") or "").strip(),
    )

    alias_tsv = str(p.get("alias_tsv") or "").strip()
    alias_map = _load_alias_map(Path(alias_tsv).expanduser()) if alias_tsv else {}

    # filters
    ph_filter_col = str(p.get("pheno_filter_col") or "").strip()
    ph_filter_values = _parse_csv_list(str(p.get("pheno_filter_values") or "").strip())

    # panel-aware phenotype ID (e.g., NAM: use Z_Num; others: use Family_Inbred_Name)
    panel_aware = _as_bool(p.get("panel_aware"), False)
    panel_col = str(p.get("panel_col") or "").strip()
    panel_values = _parse_csv_list(str(p.get("panel_values") or "").strip())
    panel_id_col = str(p.get("panel_id_col") or "").strip()

    keep_genotype = str(p.get("keep_genotype") or "first").strip().lower()
    phenotype_reduce = str(p.get("phenotype_reduce") or "mean_numeric").strip().lower()

    write_pheno = _as_bool(p.get("write_phenotype_tsv"), True)
    pheno_out_name = str(p.get("phenotype_out_name") or "phenotype.tsv").strip() or "phenotype.tsv"
    write_list = _as_bool(p.get("write_sample_list"), True)
    list_name = str(p.get("sample_list_name") or "selected_samples.txt").strip() or "selected_samples.txt"
    max_pairs = int(p.get("max_pairs_tsv") or 200000)
    # extraction
    # New simplified UI: extract_genotype + genotype_file + genotype_out_prefix (engine defaults to PLINK2)
    extract_genotype = _as_bool(p.get('extract_genotype') or p.get('do_extract_genotype'), False)
    genotype_file = str(p.get('genotype_file') or p.get('genotype_path') or '').strip()
    genotype_out_prefix = str(p.get('genotype_out_prefix') or p.get('genotype_prefix') or '').strip()

    extract_mode = str(p.get("extract_mode") or "none").strip().lower()
    vcf_path = str(p.get("vcf_path") or "").strip()
    bcftools_bin = str(p.get("bcftools_bin") or "bcftools").strip() or "bcftools"
    vcf_out = str(p.get("vcf_out") or "").strip()

    plink_prefix = str(p.get("plink_prefix") or "").strip()
    plink2_bin = str(p.get("plink2_bin") or "plink2").strip() or "plink2"
    plink_out_prefix = str(p.get("plink_out_prefix") or "").strip()

    if extract_genotype:
        # Override old extract_mode settings
        extract_mode = 'plink2_auto'
        if not genotype_out_prefix:
            genotype_out_prefix = 'genotype_subset'

    # --- load genotype id table ---
    import pandas as pd

    def _gid_from_plink_fam(fam_path: Path) -> "pd.DataFrame":
        rows = []
        for ln in _read_text(fam_path).splitlines():
            if not ln.strip():
                continue
            parts = ln.strip().split()
            if len(parts) < 2:
                continue
            rows.append((parts[0], parts[1]))
        return pd.DataFrame(rows, columns=["fid", "id"]) if rows else pd.DataFrame(columns=["fid", "id"])

    def _gid_from_genotype(geno_path: str) -> "pd.DataFrame":
        gp = Path(geno_path).expanduser()
        nm = gp.name.lower()

        # PLINK prefix/.bed/.pgen
        if gp.suffix.lower() in ('.bed', '.bim', '.fam', '.pgen', '.pvar', '.psam'):
            gp = gp.with_suffix('')
        fam = gp.with_suffix('.fam')
        psam = gp.with_suffix('.psam')
        if fam.exists():
            return _gid_from_plink_fam(fam)
        if psam.exists():
            # plink2 psam: first column is '#IID' (or IID). FID may be missing.
            try:
                df = _read_table_auto(psam)
                df.columns = [str(c).strip() for c in df.columns]
                iid_col = '#IID' if '#IID' in df.columns else ('IID' if 'IID' in df.columns else _first_col(df))
                iids = df[iid_col].astype(str).map(lambda z: str(z).strip())
                iids = [x for x in iids.tolist() if x]
                return pd.DataFrame({'fid': ['0']*len(iids), 'id': iids})
            except Exception:
                return pd.DataFrame(columns=['fid','id'])

        # VCF/BCF: ask PLINK2 to write a .fam
        if gp.exists() and (nm.endswith('.vcf') or nm.endswith('.vcf.gz') or nm.endswith('.bcf') or nm.endswith('.bcf.gz')):
            tmp_dir = root / "_tmp_id_from_genotype"
            tmp_dir.mkdir(parents=True, exist_ok=True)
            outp = tmp_dir / "ids"
            stdout_p = tmp_dir / "plink2_makejustfam_stdout.txt"
            stderr_p = tmp_dir / "plink2_makejustfam_stderr.txt"

            cmd = [plink2_bin, '--allow-extra-chr']
            if nm.endswith('.bcf') or nm.endswith('.bcf.gz'):
                cmd += ['--bcf', str(gp)]
            else:
                cmd += ['--vcf', str(gp)]
            cmd += ['--make-just-fam', '--out', str(outp)]
            try:
                _run(cmd, cwd=tmp_dir, stdout_path=stdout_p, stderr_path=stderr_p)
                fam2 = Path(str(outp) + '.fam')
                if fam2.exists():
                    return _gid_from_plink_fam(fam2)
            except Exception:
                pass

            # Fallback: bcftools query -l (may fail if header is non-standard)
            try:
                rc, out, err = _run_capture([bcftools_bin, 'query', '-l', str(gp)])
                _write_text(tmp_dir / 'bcftools_query_l_stdout.txt', out)
                _write_text(tmp_dir / 'bcftools_query_l_stderr.txt', err)
                if rc == 0:
                    ids = [ln.strip() for ln in out.splitlines() if ln.strip()]
                    return pd.DataFrame({'fid': ['0']*len(ids), 'id': ids})
            except Exception:
                pass

        return pd.DataFrame(columns=['fid','id'])

    # Load genotype IDs from the provided ID table, or derive from genotype_file when extracting genotype.
    gid_df = None
    if id_table_tsv and id_table_tsv.exists():
        gid_df = _read_table_auto(id_table_tsv)
    else:
        if extract_genotype and genotype_file:
            gid_df = _gid_from_genotype(genotype_file)
            if gid_df is None or gid_df.shape[0] == 0:
                raise SystemExit(
                    "id_table_tsv is missing, and failed to derive sample IDs from genotype_file. "
                    "Please provide genotype IDs table or a readable genotype file."
                )
            # If we synthesized a table, prefer its conventional column names.
            if geno_id_col_req not in {'id', 'iid'}:
                geno_id_col_req = 'id'
            if not geno_fid_col_req:
                geno_fid_col_req = 'fid'
        else:
            raise SystemExit(f"id_table_tsv not found: {id_table_tsv_str}")
    gid_df.columns = [str(c).strip() for c in gid_df.columns]

    iid_col = _safe_col(gid_df, geno_id_col_req) or _first_col(gid_df)
    fid_col = _safe_col(gid_df, geno_fid_col_req) if geno_fid_col_req else None

    full_ids = gid_df[iid_col].astype(str)
    full_ids = full_ids[full_ids.notna()].map(lambda z: str(z).strip())
    full_ids = full_ids[full_ids != ""]

    if full_ids.empty:
        raise SystemExit("No genotype IDs found in id_table_tsv")

    # Keep original order; drop duplicates (full IDs)
    seen = set()
    ordered_full_ids: List[str] = []
    for x in full_ids.tolist():
        if x not in seen:
            seen.add(x)
            ordered_full_ids.append(x)

    # FID list (optional)
    fid_list: List[str] = []
    if fid_col:
        fid_s = gid_df[fid_col].astype(str)
        fid_s = fid_s[full_ids.index].fillna("0")
        fid_list = [str(v).strip() if str(v).strip() else "0" for v in fid_s.tolist()]
    else:
        fid_list = ["0"] * len(ordered_full_ids)

    # compute genotype keys
    geno_keys: List[str] = []
    for x in ordered_full_ids:
        k = _apply_transform(x, geno_transform, delim, rx)
        k = _normalize(k, opt)
        if alias_map and k in alias_map:
            k = alias_map[k]
        geno_keys.append(k)

    # mapping key -> list of genotype (fid,iid)
    geno_map: Dict[str, List[Tuple[str, str]]] = {}
    for fid, iid, k in zip(fid_list, ordered_full_ids, geno_keys):
        if k == "":
            continue
        geno_map.setdefault(k, []).append((fid, iid))

    # Apply genotype duplicate reduction
    if keep_genotype == "first":
        geno_map = {k: [v[0]] for k, v in geno_map.items()}

    # --- load phenotype ---
    ph_df = _read_phenotype(phenotype_path, ph_sheet)
    ph_df.columns = [str(c).strip() for c in ph_df.columns]

    # optional phenotype filter
    if ph_filter_col and ph_filter_col in ph_df.columns and ph_filter_values:
        ph_df = ph_df[ph_df[ph_filter_col].astype(str).isin(ph_filter_values)].copy()

    # determine phenotype id column
    ph_id_col = _safe_col(ph_df, ph_id_col_req)
    if not ph_id_col:
        # heuristic for common maize sheet
        if "Family_Inbred_Name" in ph_df.columns:
            ph_id_col = "Family_Inbred_Name"
        else:
            ph_id_col = _first_col(ph_df)

    # choose phenotype id raw values (optionally panel-aware)
    if panel_aware and panel_col and panel_id_col and panel_values:
        if panel_col not in ph_df.columns:
            raise SystemExit(f"panel_aware requested but panel_col not found in phenotype: {panel_col}")
        if panel_id_col not in ph_df.columns:
            raise SystemExit(f"panel_aware requested but panel_id_col not found in phenotype: {panel_id_col}")
        pv = set([str(x) for x in panel_values])
        panel_series = ph_df[panel_col].astype(str).fillna("")
        use_panel = panel_series.isin(pv)
        ph_ids_raw = ph_df[ph_id_col].astype(str).fillna("")
        ph_ids_alt = ph_df[panel_id_col].astype(str).fillna("")
        ph_ids_used = ph_ids_raw.where(~use_panel, ph_ids_alt)
        ph_df = ph_df.copy()
        ph_df["__pheno_id_raw"] = ph_ids_used
        ph_df["__pheno_id_source"] = [panel_id_col if bool(u) else ph_id_col for u in use_panel.tolist()]
        ph_ids_raw = ph_ids_used
    else:
        ph_df = ph_df.copy()
        ph_df["__pheno_id_raw"] = ph_df[ph_id_col].astype(str).fillna("")
        ph_df["__pheno_id_source"] = ph_id_col
        ph_ids_raw = ph_df["__pheno_id_raw"].astype(str).fillna("")

    ph_keys: List[str] = []
    for x in ph_ids_raw.tolist():
        k = _apply_transform(x, pheno_transform, delim, rx)
        k = _normalize(k, opt)
        if alias_map and k in alias_map:
            k = alias_map[k]
        ph_keys.append(k)

    ph_df["__pheno_key"] = ph_keys

    # drop empty keys
    ph_df2 = ph_df[ph_df["__pheno_key"].astype(str).str.len() > 0].copy()

    # phenotype reduction
    if phenotype_reduce == "first":
        ph_red = ph_df2.drop_duplicates(subset=["__pheno_key"], keep="first").copy()
    elif phenotype_reduce == "mean_numeric":
        # group by key; take first for non-numeric, mean for numeric columns
        key = ph_df2["__pheno_key"]
        num_cols = []
        for c in ph_df2.columns:
            if c in {"__pheno_key"} or str(c).startswith("__"):
                continue
            if c == ph_id_col:
                continue
            # numeric detection: if >=80% convertible
            s = ph_df2[c]
            num = _to_numeric_series(s)
            denom = int(s.replace({"": None}).notna().sum())
            if denom >= 10 and float(num.notna().sum()) / float(max(1, denom)) >= 0.8:
                num_cols.append(c)
        # build aggregated
        import pandas as pd

        agg_rows = []
        for k, sub in ph_df2.groupby("__pheno_key", sort=False):
            row = {}
            # keep representative raw id (first) for reporting
            row["__pheno_id_raw"] = str(sub["__pheno_id_raw"].iloc[0])
            row["__pheno_id_source"] = str(sub["__pheno_id_source"].iloc[0])
            row["__pheno_key"] = k
            # keep filter col if present
            for c in ph_df2.columns:
                if c in {ph_id_col, "__pheno_key"} or str(c).startswith("__"):
                    continue
                if c in num_cols:
                    row[c] = float(_to_numeric_series(sub[c]).mean(skipna=True)) if _to_numeric_series(sub[c]).notna().any() else ""
                else:
                    # take first non-empty
                    vals = [str(v) for v in sub[c].tolist()]
                    v2 = ""
                    for v in vals:
                        if v and v.strip() and v.strip().lower() not in {"na", "nan", ".", "-"}:
                            v2 = v
                            break
                    row[c] = v2
            agg_rows.append(row)
        ph_red = pd.DataFrame(agg_rows)
    else:
        ph_red = ph_df2.copy()

    # build phenotype key set
    ph_key_set = set([str(k) for k in ph_red["__pheno_key"].astype(str).tolist()])
    geno_key_set = set([str(k) for k in geno_map.keys()])

    matched_keys = sorted(list(geno_key_set & ph_key_set))
    unmatched_geno = sorted(list(geno_key_set - ph_key_set))
    unmatched_pheno = sorted(list(ph_key_set - geno_key_set))

    # duplicate key info
    dup_geno_keys = sorted([k for k, v in geno_map.items() if len(v) > 1])

    # phenotype duplicates from original
    ph_counts = ph_df2["__pheno_key"].value_counts()
    dup_ph_keys = sorted([k for k, n in ph_counts.items() if int(n) > 1])

    # Build matched pairs (key, fid, iid, pheno_row)
    # Use reduced phenotype table for values
    ph_red_index = {str(r["__pheno_key"]): i for i, r in ph_red.iterrows()}

    pairs: List[Dict[str, str]] = []
    selected_samples: List[Tuple[str, str]] = []

    for k in matched_keys:
        g_list = geno_map.get(k, [])
        if not g_list:
            continue
        selected_samples.extend(g_list)
        # create a row per selected genotype sample (capped for report)
        if len(pairs) < max_pairs:
            ph_i = ph_red_index.get(k)
            ph_raw = ""
            if ph_i is not None:
                try:
                    ph_raw = str(ph_red.loc[ph_i, "__pheno_id_raw"]) if "__pheno_id_raw" in ph_red.columns else ""
                except Exception:
                    ph_raw = ""
            for fid, iid in g_list:
                if len(pairs) >= max_pairs:
                    break
                pairs.append({
                    "key": k,
                    "geno_fid": str(fid),
                    "geno_iid": str(iid),
                    "pheno_id_raw": ph_raw,
                })

    # --- write outputs ---
    root.mkdir(parents=True, exist_ok=True)

    # summary
    n_geno_full = len(ordered_full_ids)
    n_geno_keys = len(geno_key_set)
    n_ph_rows = int(ph_df.shape[0])
    n_ph_red = int(ph_red.shape[0])

    summary_rows = [
        ("genotype_full_ids", n_geno_full),
        ("genotype_unique_keys", n_geno_keys),
        ("phenotype_rows", n_ph_rows),
        ("phenotype_rows_after_filter", int(ph_df2.shape[0])),
        ("phenotype_unique_keys", n_ph_red),
        ("matched_keys", len(matched_keys)),
        ("matched_rate_vs_pheno_keys", float(len(matched_keys)) / float(max(1, n_ph_red))),
        ("unmatched_genotype_keys", len(unmatched_geno)),
        ("unmatched_phenotype_keys", len(unmatched_pheno)),
        ("duplicate_genotype_keys", len(dup_geno_keys)),
        ("duplicate_phenotype_keys", len(dup_ph_keys)),
        ("keep_genotype", keep_genotype),
        ("phenotype_reduce", phenotype_reduce),
        ("geno_transform", geno_transform),
        ("pheno_transform", pheno_transform),
        ("delim", delim),
        ("phenotype_id_col", ph_id_col),
        ("genotype_id_col", iid_col),
        ("panel_aware", panel_aware),
        ("panel_col", panel_col),
        ("panel_values", ",".join(panel_values)),
        ("panel_id_col", panel_id_col),
    ]

    import pandas as pd

    df_sum = pd.DataFrame(summary_rows, columns=["metric", "value"])
    sum_path = root / "match_summary.tsv"
    df_sum.to_csv(sum_path, sep="\t", index=False)

    # matched keys
    pd.DataFrame({"key": matched_keys}).to_csv(root / "matched_keys.tsv", sep="\t", index=False)

    # unmatched lists
    pd.DataFrame({"key": unmatched_geno}).to_csv(root / "unmatched_genotype.tsv", sep="\t", index=False)
    pd.DataFrame({"key": unmatched_pheno}).to_csv(root / "unmatched_phenotype.tsv", sep="\t", index=False)

    # duplicates
    dup_rows = []
    for k in dup_geno_keys:
        for fid, iid in geno_map.get(k, []):
            dup_rows.append({"key": k, "geno_fid": fid, "geno_iid": iid})
    for k in dup_ph_keys:
        dup_rows.append({"key": k, "pheno_rows": int(ph_counts.get(k, 0))})
    if dup_rows:
        pd.DataFrame(dup_rows).to_csv(root / "duplicate_keys.tsv", sep="\t", index=False)
    else:
        pd.DataFrame({"note": ["no duplicate keys detected"]}).to_csv(root / "duplicate_keys.tsv", sep="\t", index=False)

    # matched pairs (capped)
    pd.DataFrame(pairs).to_csv(root / "matched_pairs.tsv", sep="\t", index=False)

    # sample list (full genotype IDs)
    sample_list_path = root / list_name
    if write_list:
        _write_text(sample_list_path, "\n".join([iid for _fid, iid in selected_samples]) + "\n")

    # phenotype.tsv for downstream GE modules
    phenotype_out_path = root / pheno_out_name
    if write_pheno:
        # expand reduced phenotype to full genotype IDs
        out_rows = []
        cols = [c for c in ph_red.columns if c != "__pheno_key" and not str(c).startswith("__")]
        # ensure id is first
        trait_cols = [c for c in cols if c != ph_id_col]

        for k in matched_keys:
            ph_i = ph_red_index.get(k)
            if ph_i is None:
                continue
            ph_row = ph_red.loc[ph_i]
            for fid, iid in geno_map.get(k, []):
                row = {"id": str(iid)}
                for c in trait_cols:
                    row[c] = ph_row.get(c, "")
                out_rows.append(row)

        out_df = pd.DataFrame(out_rows)
        if out_df.shape[0] == 0:
            out_df = pd.DataFrame({"id": []})
        out_df.to_csv(phenotype_out_path, sep="\t", index=False)

    # artifacts.json
    artifacts = {
        "match_summary_tsv": str(sum_path),
        "matched_keys_tsv": str(root / "matched_keys.tsv"),
        "matched_pairs_tsv": str(root / "matched_pairs.tsv"),
        "phenotype_tsv": str(phenotype_out_path) if write_pheno else "",
        "sample_list": str(sample_list_path) if write_list else "",
        "n_matched_keys": len(matched_keys),
        "n_selected_samples": len(selected_samples),
        "extract_mode": extract_mode,
        "subset_vcf": "",
        "subset_plink_prefix": "",
    }

    # --- optional extraction ---
    try:

        if extract_mode == "plink2_auto":
            if not genotype_file:
                raise RuntimeError("extract_genotype is enabled but genotype_file is empty")

            in_path = Path(genotype_file).expanduser()
            # determine input type by extension / existence
            name = in_path.name.lower()
            is_vcf = (name.endswith('.vcf') or name.endswith('.vcf.gz') or name.endswith('.bcf') or name.endswith('.bcf.gz')) and in_path.exists()

            # output prefix lives under root
            safe_pref = re.sub(r"[^A-Za-z0-9_.-]+", "_", genotype_out_prefix)
            out_pref = (root / safe_pref)
            out_pref.parent.mkdir(parents=True, exist_ok=True)

            keep_path = root / "keep.fid_iid.tsv"
            lines = []
            for fid, iid in selected_samples:
                fid2 = fid if fid and str(fid).strip() else "0"
                lines.append(f"{fid2}	{iid}")
            _write_text(keep_path, "\n".join(lines) + "\n")

            if is_vcf:
                # PLINK2 is often more tolerant than bcftools for public VCF headers.
                cmd = [
                    plink2_bin,
                    "--vcf", str(in_path),
                    "--keep", str(keep_path),
                    "--export", "vcf", "bgz",
                    "--out", str(out_pref),
                ]
                _run(cmd, cwd=None, stdout_path=root / "plink2_extract_stdout.txt", stderr_path=root / "plink2_extract_stderr.txt")
                out_v = Path(str(out_pref) + ".vcf.gz")
                if out_v.exists():
                    artifacts["subset_vcf"] = str(out_v)
                    # optional index
                    try:
                        _run([bcftools_bin, "index", "-t", str(out_v)], cwd=None,
                             stdout_path=root / "bcftools_index_stdout.txt", stderr_path=root / "bcftools_index_stderr.txt")
                    except Exception:
                        pass
                else:
                    raise RuntimeError("plink2 --export vcf bgz did not produce .vcf.gz")

            else:
                # Treat as PLINK prefix (accept .bed/.bim/.fam or .pgen/.pvar/.psam)
                pref = in_path
                if pref.suffix.lower() in {".bed", ".bim", ".fam", ".pgen", ".pvar", ".psam"}:
                    pref = pref.with_suffix("")
                cmd = [plink2_bin]
                if (Path(str(pref) + ".pgen")).exists():
                    cmd += ["--pfile", str(pref)]
                else:
                    cmd += ["--bfile", str(pref)]
                cmd += ["--keep", str(keep_path), "--make-bed", "--out", str(out_pref)]
                _run(cmd, cwd=None, stdout_path=root / "plink2_extract_stdout.txt", stderr_path=root / "plink2_extract_stderr.txt")
                if Path(str(out_pref) + ".bed").exists():
                    artifacts["subset_plink_prefix"] = str(out_pref)
                else:
                    raise RuntimeError("plink2 --make-bed did not produce .bed")

        elif extract_mode == "bcftools":

            if not vcf_path:
                raise RuntimeError("extract_mode=bcftools requires vcf_path")
            vcf_in = Path(vcf_path).expanduser()
            if not vcf_in.exists():
                raise RuntimeError(f"VCF not found: {vcf_in}")

            # decide output format by extension
            if vcf_out:
                out_v = Path(vcf_out).expanduser()
            else:
                out_v = root / "subset.vcf.gz"
            out_v.parent.mkdir(parents=True, exist_ok=True)

            requested = []
            seen2 = set()
            for _fid, iid in selected_samples:
                s = str(iid).strip()
                if s and s not in seen2:
                    seen2.add(s)
                    requested.append(s)

            if not requested:
                raise RuntimeError("No selected samples (requested list is empty). Check matching settings.")

            # Get sample list from VCF header (fast; header-only)
            rc, so, se = _run_capture([bcftools_bin, "query", "-l", str(vcf_in)])
            _write_text(root / "bcftools_query_l_stdout.txt", so)
            _write_text(root / "bcftools_query_l_stderr.txt", se)
            if rc != 0:
                tail = (se or so or "")[-2000:]
                raise RuntimeError(f"bcftools query -l failed (code={rc}). tail: {tail}")

            vcf_samples = [ln.strip() for ln in so.splitlines() if ln.strip()]
            _write_text(root / "vcf_samples.txt", "\n".join(vcf_samples) + "\n")
            vcf_set = set(vcf_samples)

            present = [s for s in requested if s in vcf_set]
            missing = [s for s in requested if s not in vcf_set]

            import pandas as pd

            pd.DataFrame([
                {"metric": "requested_samples", "value": len(requested)},
                {"metric": "present_in_vcf", "value": len(present)},
                {"metric": "missing_from_vcf", "value": len(missing)},
                {"metric": "vcf_total_samples", "value": len(vcf_samples)},
            ]).to_csv(root / "bcftools_sample_check.tsv", sep="\t", index=False)

            # Save small previews
            _write_text(root / "keep_samples_requested_head.txt", "\n".join(requested[:50]) + "\n")
            _write_text(root / "keep_samples_present_head.txt", "\n".join(present[:50]) + "\n")
            _write_text(root / "keep_samples_missing_head.txt", "\n".join(missing[:50]) + "\n")
            _write_text(root / "vcf_samples_head.txt", "\n".join(vcf_samples[:50]) + "\n")

            if not present:
                raise RuntimeError(
                    "0 requested samples were found in the VCF header. "
                    "Likely causes: (1) ids.tsv was generated from a different VCF, "
                    "(2) sample IDs differ (spaces/case), or (3) you matched on wrong phenotype ID (e.g., NAM needs Z_Num). "
                    "See bcftools_sample_check.tsv and vcf_samples_head.txt."
                )

            # bcftools view -S expects sample names exactly as in VCF (full IDs)
            keep_path = root / "keep_samples.txt"
            _write_text(keep_path, "\n".join(present) + "\n")

            # output format
            out_flag = "-Oz"
            if out_v.suffix.lower() == ".bcf":
                out_flag = "-Ob"
            elif not (out_v.name.lower().endswith(".vcf.gz") or out_v.suffix.lower() == ".gz"):
                out_flag = "-Ov"

            # ---- robust VCF extraction ----
            # Some public maize VCFs contain malformed '##SAMPLE=<...>' header meta-lines or incomplete tag declarations,
            # which can cause bcftools view -S to abort even when sample IDs match. In that case we retry after stripping
            # problematic meta-lines, and finally fall back to plink2 export.
            engine_used = ""
            ok = False
            err_notes: List[str] = []

            try:
                _run(
                    [
                        bcftools_bin,
                        "view",
                        "--force-samples",
                        "-S",
                        str(keep_path),
                        out_flag,
                        "-o",
                        str(out_v),
                        str(vcf_in),
                    ],
                    cwd=None,
                    stdout_path=root / "bcftools_stdout.txt",
                    stderr_path=root / "bcftools_stderr.txt",
                )
                ok = True
                engine_used = "bcftools"
            except Exception as e1:
                err_notes.append(str(e1))
                se1 = _read_text(root / "bcftools_stderr.txt")
                # retry with sanitized header (drop ##SAMPLE structured lines)
                if ("bcf_hdr_parse_line" in se1) or ("Could not parse the header line" in se1) or ("Undefined tags in the header" in se1):
                    cat = "zcat" if str(vcf_in).endswith(".gz") else "cat"
                    # NOTE: use grep -v to drop malformed SAMPLE meta-lines; keep the rest unchanged.
                    cmdline = (
                        f"{cat} {shlex.quote(str(vcf_in))} "
                        f"| grep -v '^##SAMPLE=<' "
                        f"| {shlex.quote(bcftools_bin)} view --force-samples -S {shlex.quote(str(keep_path))} "
                        f"{out_flag} -o {shlex.quote(str(out_v))} -"
                    )
                    try:
                        _run_shell(
                            cmdline,
                            cwd=None,
                            stdout_path=root / "bcftools_sanitize_stdout.txt",
                            stderr_path=root / "bcftools_sanitize_stderr.txt",
                        )
                        ok = True
                        engine_used = "bcftools_sanitized"
                    except Exception as e2:
                        err_notes.append(str(e2))

            if not ok:
                # final fallback: plink2 (often more tolerant of non-standard VCF meta-lines)
                keep2 = root / "keep_plink2.tsv"
                _write_text(keep2, "\n".join([f"0\t{s}" for s in present]) + "\n")
                out_pref2 = (root / "subset_plink2" / "subset")
                out_pref2.parent.mkdir(parents=True, exist_ok=True)

                # try plink2 export vcf bgz
                cmd = [
                    plink2_bin,
                    "--vcf", str(vcf_in),
                    "--keep", str(keep2),
                    "--export", "vcf", "bgz",
                    "--out", str(out_pref2),
                ]
                try:
                    _run(cmd, cwd=None, stdout_path=root / "plink2_export_stdout.txt", stderr_path=root / "plink2_export_stderr.txt")
                    # plink2 writes <out>.vcf.gz
                    out_v2 = Path(str(out_pref2) + ".vcf.gz")
                    if out_v2.exists():
                        out_v = out_v2
                        ok = True
                        engine_used = "plink2_export_vcf"
                except Exception as e3:
                    err_notes.append(str(e3))
            # record which engine succeeded (if any)
            _write_text(root / "vcf_extract_engine.txt", engine_used + "\n")
            if not ok:
                raise RuntimeError(
                    "VCF extraction failed via bcftools (raw) and bcftools (sanitized header), and plink2 fallback.\n"
                    + "\n---\n".join(err_notes[-3:])
                )
            
            # index (best-effort)
            try:
                if out_flag == "-Oz":
                    _run(
                        [bcftools_bin, "index", "-t", str(out_v)],
                        cwd=None,
                        stdout_path=root / "bcftools_index_stdout.txt",
                        stderr_path=root / "bcftools_index_stderr.txt",
                    )
            except Exception:
                pass

            # sanity check: at least one variant record
            if not _vcf_has_records(out_v):
                raise RuntimeError(
                    "Subset VCF contains no variant records (header-only). "
                    "If you expected records, double-check the VCF path and sample IDs. "
                    "See bcftools_stdout/stderr and bcftools_sample_check.tsv."
                )

            artifacts["subset_vcf"] = str(out_v)

        elif extract_mode == "plink2":
            if not plink_prefix:
                raise RuntimeError("extract_mode=plink2 requires plink_prefix")
            pref = Path(plink_prefix).expanduser()
            # accept suffix
            if pref.suffix.lower() in {".bed", ".bim", ".fam", ".pgen", ".pvar", ".psam"}:
                pref = pref.with_suffix("")

            out_pref = Path(plink_out_prefix).expanduser() if plink_out_prefix else (root / "subset_plink" / "subset")
            out_pref.parent.mkdir(parents=True, exist_ok=True)

            keep_path = root / "keep.fid_iid.tsv"
            lines = []
            for fid, iid in selected_samples:
                fid2 = fid if fid and str(fid).strip() else "0"
                lines.append(f"{fid2}\t{iid}")
            _write_text(keep_path, "\n".join(lines) + "\n")

            # Use plink2; allow both bfile and pfile by trying pfile first
            cmd = [plink2_bin]
            if (Path(str(pref) + ".pgen")).exists():
                cmd += ["--pfile", str(pref)]
            else:
                cmd += ["--bfile", str(pref)]
            cmd += ["--keep", str(keep_path), "--make-bed", "--out", str(out_pref)]

            _run(cmd, cwd=None, stdout_path=root / "plink2_stdout.txt", stderr_path=root / "plink2_stderr.txt")
            artifacts["subset_plink_prefix"] = str(out_pref)

    except Exception as e:
        _write_text(root / "extract_error.txt", str(e))
        artifacts["extract_error"] = str(e)

    _write_text(root / "artifacts.json", json.dumps(artifacts, ensure_ascii=False, indent=2))

    print("[prep_split_merge_extract] done")
    print(json.dumps({
        "out": str(root),
        "matched_keys": len(matched_keys),
        "selected_samples": len(selected_samples),
    }, ensure_ascii=False))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
