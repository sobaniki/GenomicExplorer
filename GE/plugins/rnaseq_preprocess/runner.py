#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from collections import Counter
from pathlib import Path

import pandas as pd

csv.field_size_limit(sys.maxsize)


def as_bool(x, default=False):
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


def read_table_auto(path: Path, sheet_name: str = "") -> pd.DataFrame:
    suf = path.suffix.lower()
    if suf in {".xlsx", ".xls"}:
        return pd.read_excel(path, sheet_name=(sheet_name or 0), dtype=object)
    for sep in [r"\s+", "\t", ","]:
        try:
            df = pd.read_csv(path, sep=sep, engine="python", dtype=object)
            if df.shape[1] > 1:
                return df
        except Exception:
            pass
    return pd.read_csv(path, sep=None, engine="python", dtype=object)


def is_numberish_series(s: pd.Series, min_frac: float = 0.8) -> bool:
    x = pd.to_numeric(s, errors="coerce")
    frac = float(x.notna().mean()) if len(x) else 0.0
    return frac >= min_frac


# -------------------------
# ID normalization / parsing
# -------------------------

RE_ZM = re.compile(r"(Zm\d{5}[A-Za-z]{1,3}\d{6})")
RE_GRM = re.compile(r"(GRMZM\d+G\d{6})")


def _strip_refgen_suffix(s: str) -> str:
    # e.g. Zm00001d027231_T007.RefGen_V4 -> Zm00001d027231_T007
    return re.sub(r"\.RefGen_V\d+$", "", s)


def _strip_version_suffix(s: str) -> str:
    # e.g. ABC.2 -> ABC
    return re.sub(r"\.\d+$", "", s)


def _strip_isoform_suffix(s: str) -> str:
    # e.g. Zm..._T007 / _P001 / .t1 / .p1
    s2 = re.sub(r"_(?:T|P)\d{3}$", "", s, flags=re.I)
    s2 = re.sub(r"\.t\d+$", "", s2, flags=re.I)
    s2 = re.sub(r"\.p\d+$", "", s2, flags=re.I)
    return s2


def extract_id_tokens(s: str, strip_isoform=True, strip_version=True):
    """Return candidate tokens in a *priority order*.

    We keep multiple variants because inputs often mix:
      - gene vs transcript vs protein
      - old IDs (GRMZM2G...) vs new IDs (Zm...)
      - RefGen suffix / version suffix / isoform suffix
    """
    s = str(s or "").strip()
    if not s:
        return []

    out: list[str] = []

    def add(v: str):
        v = str(v or "").strip()
        if v and v not in out:
            out.append(v)

    add(s)

    for m in RE_ZM.finditer(s):
        add(m.group(1))
    for m in RE_GRM.finditer(s):
        add(m.group(1))

    t = _strip_refgen_suffix(s)
    add(t)

    # add isoform-stripped variant (gene-level) if requested
    if strip_isoform:
        add(_strip_isoform_suffix(t))

    # add version-stripped variant (e.g., .1/.2)
    if strip_version:
        add(_strip_version_suffix(t))
        if strip_isoform:
            add(_strip_version_suffix(_strip_isoform_suffix(t)))

    # common separators
    add(s.split()[0])

    return [x for x in out if x]


def parse_attrs(attr_str: str) -> dict[str, str]:
    d: dict[str, str] = {}
    for item in (attr_str or "").split(";"):
        if not item:
            continue
        if "=" not in item:
            continue
        k, v = item.split("=", 1)
        k = k.strip()
        v = v.strip()
        if k and v:
            d[k] = v
    return d


def parse_gff3_gene_index(
    path: Path | None,
    gene_field: str = "Name",
    strip_isoform: bool = True,
    strip_version: bool = True,
):
    """Parse GFF3 and build mappings to *gene-level* IDs.

    Returns:
      gene_names: set[str]
      gff_any_to_gene: dict[str,str]  # gene tokens and transcript tokens -> gene_name
      gene_to_rep_transcript: dict[str,str]  # gene_name -> representative transcript (if tag exists)

    Representative transcript:
      - If mRNA has attribute longest=1 (Phytozome RefGen style), use it.
    """
    gene_names: set[str] = set()
    gene_token_to_gene: dict[str, str] = {}
    gff_any_to_gene: dict[str, str] = {}
    gene_to_rep_tr: dict[str, str] = {}

    if not path or not path.exists():
        return gene_names, gff_any_to_gene, gene_to_rep_tr

    gene_field = (gene_field or "Name").strip() or "Name"

    # pass1: gene lines -> gene_name + tokens
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for ln in f:
            if ln.startswith("#"):
                continue
            parts = ln.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            ftype = parts[2]
            if ftype != "gene":
                continue
            attrs = parse_attrs(parts[8])
            gname_raw = attrs.get(gene_field) or attrs.get("Name") or attrs.get("ID") or ""
            gid_raw = attrs.get("ID") or ""
            gname = extract_id_tokens(gname_raw, strip_isoform=False, strip_version=strip_version)[0] if gname_raw else ""
            if not gname:
                gname = extract_id_tokens(gid_raw, strip_isoform=False, strip_version=strip_version)[0] if gid_raw else ""
            gname = _strip_refgen_suffix(gname)
            if strip_version:
                gname = _strip_version_suffix(gname)
            if not gname:
                continue
            gene_names.add(gname)

            # map multiple tokens (Name/ID variants) to gene_name
            for tok in extract_id_tokens(gname_raw, strip_isoform=strip_isoform, strip_version=strip_version):
                gene_token_to_gene[tok] = gname
            for tok in extract_id_tokens(gid_raw, strip_isoform=strip_isoform, strip_version=strip_version):
                gene_token_to_gene[tok] = gname

    # pass2: transcript lines -> gene
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for ln in f:
            if ln.startswith("#"):
                continue
            parts = ln.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            ftype = parts[2]
            if ftype not in {"mRNA", "transcript"}:
                continue
            attrs = parse_attrs(parts[8])
            tid_raw = attrs.get("Name") or attrs.get("ID") or ""
            pid_raw = attrs.get("Parent") or ""
            if not tid_raw or not pid_raw:
                continue

            parent_gene = None
            for ptok in extract_id_tokens(pid_raw, strip_isoform=strip_isoform, strip_version=strip_version):
                if ptok in gene_token_to_gene:
                    parent_gene = gene_token_to_gene[ptok]
                    break
            if not parent_gene:
                continue

            # transcript token(s) -> gene
            for tok in extract_id_tokens(tid_raw, strip_isoform=strip_isoform, strip_version=strip_version):
                gff_any_to_gene[tok] = parent_gene

            # rep transcript by tag
            if str(attrs.get("longest", "")).strip() == "1":
                # prefer Name-like (no RefGen suffix) for reporting
                tr_id = extract_id_tokens(tid_raw, strip_isoform=False, strip_version=strip_version)[0]
                tr_id = _strip_refgen_suffix(tr_id)
                if strip_version:
                    tr_id = _strip_version_suffix(tr_id)
                gene_to_rep_tr[parent_gene] = tr_id

    # gene tokens -> gene
    gff_any_to_gene.update(gene_token_to_gene)
    return gene_names, gff_any_to_gene, gene_to_rep_tr


def parse_fasta_id_tokens(
    path: Path | None,
    strip_isoform: bool = True,
    strip_version: bool = True,
):
    ids = set()
    if not path or not path.exists():
        return ids
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for ln in f:
            if not ln.startswith(">"):
                continue
            header = ln[1:].strip()
            for part in header.split():
                if "=" in part:
                    _, v = part.split("=", 1)
                    for tok in extract_id_tokens(v, strip_isoform=strip_isoform, strip_version=strip_version):
                        ids.add(tok)
                else:
                    for tok in extract_id_tokens(part, strip_isoform=strip_isoform, strip_version=strip_version):
                        ids.add(tok)
    return ids


def _infer_maizedb_target_tag_from_gff(gff_any_to_gene: dict[str, str]) -> str | None:
    """Infer which AGPv tag matches the provided GFF3 gene namespace (best-effort).

    For maize B73:
      - AGPv4 gene IDs typically look like: Zm00001d027231
      - AGPv5 gene IDs typically look like: Zm00001eb000020
    """
    if not gff_any_to_gene:
        return None
    genes = list(set(gff_any_to_gene.values()))
    if not genes:
        return None
    # sample up to 5000 genes
    genes = genes[:5000]
    n4 = sum(1 for g in genes if re.fullmatch(r"Zm\d{5}d\d{6}", str(g)))
    n5 = sum(1 for g in genes if re.fullmatch(r"Zm\d{5}eb\d{6}", str(g)))
    if n4 >= max(10, n5 * 2):
        return "AGPv4"
    if n5 >= max(10, n4 * 2):
        return "AGPv5"
    return None


def build_alias_map(
    path: Path | None,
    target_mode: str,
    gff_any_to_gene: dict[str, str],
    valid_fa: set[str],
    strip_isoform: bool = True,
    strip_version: bool = True,
    alias_format: str = "auto",
    alias_target_assembly: str = "auto (from GFF3)",
    custom_target_regex: str = "",
):
    """Build alias mapping.

    Supported formats:
      1) Generic alias TSV: each row may contain multiple IDs. We map *every* token found in the row
         to a chosen target (prefer GFF3 gene when target_mode is gff_gene_name/auto).

      2) MaizeDB genes_to_alias_ids.tsv-like format (special):
         Rows are grouped by a *primary* ID (usually AGPv5-like, e.g. Zm00001eb...).
         Different assemblies/versions appear on different rows (AGPv3/AGPv4/AGPv5 ...).
         In this case we **group by primary** and map *all* aliases in that group to the chosen target
         (prefer the alias that belongs to the GFF3 namespace).

    Returns:
      mp: dict[str,str] mapping any token -> target
      rows: number of non-empty lines processed
    """
    mp: dict[str, str] = {}
    rows = 0
    if not path or not path.exists():
        return mp, rows

    # --- Detect MaizeDB special format (best-effort) ---
    def _detect_maizedb_format(p: Path, n: int = 200) -> bool:
        tot = 0
        hit = 0
        with p.open("r", encoding="utf-8", errors="ignore") as f:
            for ln in f:
                if not ln.strip():
                    continue
                row = ln.rstrip("\n").split("\t")
                if len(row) < 4:
                    continue
                tot += 1
                if re.search(r"AGPv\d", str(row[3])):
                    hit += 1
                if tot >= n:
                    break
        return tot >= 10 and hit >= max(5, int(tot * 0.2))

    alias_format_norm = str(alias_format or "auto").strip().lower()
    force_maizedb = "maizedb" in alias_format_norm
    force_generic = "generic" in alias_format_norm and not force_maizedb

    is_maizedb = False
    if force_maizedb:
        is_maizedb = True
    elif force_generic:
        is_maizedb = False
    else:
        is_maizedb = _detect_maizedb_format(path)

    # desired target assembly tag for MaizeDB tables
    desired_tag = None
    if is_maizedb:
        ta = str(alias_target_assembly or "").strip()
        if ta in {"AGPv3", "AGPv4", "AGPv5"}:
            desired_tag = ta
        else:
            desired_tag = _infer_maizedb_target_tag_from_gff(gff_any_to_gene)

    # Optional custom regex for target gene namespace (useful for non-maize crops / custom tables)
    target_re = None
    if custom_target_regex:
        try:
            target_re = re.compile(custom_target_regex)
        except Exception:
            target_re = None

    # --- Helper: pick a token list but keep only gene-like tokens for maize ---
    def _looks_like_gene_id(tok: str) -> bool:
        if not tok:
            return False
        t = str(tok)
        return t.startswith("Zm") or t.startswith("GRMZM") or bool(RE_ZM.search(t) or RE_GRM.search(t))

    # --- MaizeDB grouped parsing ---
    if is_maizedb:
        groups: dict[str, list[tuple[list[str], str | None]]] = {}
        # store (tokens_in_row, tag)
        ambi: list[dict[str, str]] = []

        def _tag(meta: str) -> str | None:
            m = re.search(r"(AGPv\d+)", str(meta))
            if m:
                return m.group(1)
            # fallback heuristics
            if "RefGen_V4" in str(meta):
                return "AGPv4"
            if "RefGen_V5" in str(meta):
                return "AGPv5"
            return None

        with path.open("r", encoding="utf-8", errors="ignore") as f:
            for ln in f:
                if not ln.strip():
                    continue
                row = ln.rstrip("\n").split("\t")
                if len(row) < 3:
                    continue
                rows += 1
                primary_raw = row[0]
                alias_raw = row[2] if len(row) >= 3 else ""
                meta = row[3] if len(row) >= 4 else ""

                primary_toks = extract_id_tokens(primary_raw, strip_isoform=strip_isoform, strip_version=strip_version)
                if not primary_toks:
                    continue
                primary = primary_toks[0]

                # collect tokens from primary + alias + (optionally) column1
                toks: list[str] = []
                for cell in [primary_raw, alias_raw, row[1] if len(row) >= 2 else ""]:
                    toks.extend(extract_id_tokens(cell, strip_isoform=strip_isoform, strip_version=strip_version))
                # keep only gene-like tokens (avoid mapping AGPv meta strings / plain B73 etc)
                toks = [t for t in toks if _looks_like_gene_id(t)]

                groups.setdefault(primary, []).append((toks, _tag(meta)))

        def _pick_target_without_gff(primary: str, entries: list[tuple[list[str], str | None]]) -> str | None:
            """When GFF3 cannot decide, choose a stable target ID within the primary group.

            Priority:
              1) custom target regex match
              2) desired_tag (AGPv3/4/5)
              3) primary itself
              4) first gene-like token
            """
            if target_re is not None:
                for toks, _tg in entries:
                    for t in toks:
                        if target_re.search(str(t)):
                            return t

            if desired_tag == "AGPv4":
                for toks, _tg in entries:
                    for t in toks:
                        if re.fullmatch(r"Zm\d{5}d\d{6}", str(t)):
                            return t
            if desired_tag == "AGPv3":
                for toks, _tg in entries:
                    for t in toks:
                        if str(t).startswith("GRMZM"):
                            return t
            if desired_tag == "AGPv5":
                # primary itself is AGPv5-like for MaizeDB tables
                if primary:
                    return primary
                for toks, _tg in entries:
                    for t in toks:
                        if re.fullmatch(r"Zm\d{5}eb\d{6}", str(t)):
                            return t

            if primary:
                return primary
            for toks, _tg in entries:
                if toks:
                    return toks[0]
            return None

        # choose target per primary and map all tokens in that group to it
        for primary, entries in groups.items():
            # candidates in GFF3 namespace
            cand_targets: list[str] = []

            # 1) prefer desired_tag if inferred
            if desired_tag:
                for toks, tg in entries:
                    if tg != desired_tag:
                        continue
                    for t in toks:
                        if t in gff_any_to_gene:
                            cand_targets.append(gff_any_to_gene[t])
                            break

            # 2) otherwise any token that maps to a GFF3 gene
            if not cand_targets:
                for toks, _tg in entries:
                    for t in toks:
                        if t in gff_any_to_gene:
                            cand_targets.append(gff_any_to_gene[t])
                            break

            # 3) fallback: if primary itself maps to GFF3 gene
            if not cand_targets and primary in gff_any_to_gene:
                cand_targets.append(gff_any_to_gene[primary])

            if not cand_targets:
                # No GFF3 guidance -> choose within-group canonical target
                fallback_target = _pick_target_without_gff(primary, entries) or primary
                # map all tokens to fallback_target
                for toks, _tg in entries:
                    for t in toks:
                        if t and t not in mp:
                            mp[t] = fallback_target
                if primary and primary not in mp:
                    mp[primary] = fallback_target
                continue

            # finalize target
            uniq = sorted(set(cand_targets))
            target = uniq[0]
            if len(uniq) > 1:
                ambi.append({"primary": primary, "targets": ",".join(uniq)})

            # map everything in this primary group to the target gene
            for toks, _tg in entries:
                for t in toks:
                    if t and t not in mp:
                        mp[t] = target
            if primary and primary not in mp:
                mp[primary] = target

        # stash ambiguity into a reserved key (written later in main)
        # (caller will ignore unknown keys; we only use it internally)
        if ambi:
            mp["__AMBIGUOUS__"] = json.dumps(ambi, ensure_ascii=False)

        return mp, rows

    # --- Generic alias table (row-wise) ---
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for ln in f:
            row = ln.rstrip("\n").split("\t")
            if not row or all(not c.strip() for c in row):
                continue
            rows += 1
            candidates: list[str] = []
            for cell in row:
                candidates.extend(extract_id_tokens(cell, strip_isoform=strip_isoform, strip_version=strip_version))

            target = None

            if target_mode in {"auto", "gff_gene_name"}:
                for c in candidates:
                    if c in gff_any_to_gene:
                        target = gff_any_to_gene[c]
                        break

            if target is None and target_mode == "auto":
                for c in candidates:
                    if c in valid_fa:
                        target = c
                        break

            if target is None:
                target_tokens = extract_id_tokens(row[0], strip_isoform=strip_isoform, strip_version=strip_version)
                target = target_tokens[0] if target_tokens else str(row[0]).strip()

            for c in candidates:
                if c and c not in mp:
                    mp[c] = target

    return mp, rows

def resolve_gene_id(
    raw,
    alias_map: dict[str, str],
    gff_any_to_gene: dict[str, str],
    valid_fa: set[str],
    target_mode: str,
    strip_isoform: bool = True,
    strip_version: bool = True,
):
    toks = extract_id_tokens(raw, strip_isoform=strip_isoform, strip_version=strip_version)

    
# 1) alias mapping (with short chaining to handle indirectly-linked version tables)
    for t in toks:
        if t in alias_map:
            cur = alias_map[t]
            if target_mode in {"auto", "gff_gene_name"}:
                seen = set()
                for _ in range(4):
                    if cur in seen:
                        break
                    seen.add(cur)
                    if cur in gff_any_to_gene:
                        return gff_any_to_gene[cur], t, "alias->gff_gene"
                    if cur in alias_map:
                        cur = alias_map[cur]
                        continue
                    break
            return cur, t, "alias"

# 2) GFF3 mapping (gene-level)
    if target_mode in {"auto", "gff_gene_name"}:
        for t in toks:
            if t in gff_any_to_gene:
                return gff_any_to_gene[t], t, "gff_gene"

    # 3) FASTA token match (only as a last resort in auto)
    if target_mode == "auto":
        for t in toks:
            if t in valid_fa:
                return t, t, "fasta"

    # 4) normalized/original
    if toks:
        return toks[0], toks[0], "normalized"
    return str(raw), str(raw), "original"


# -------------------------
# table orientation detection
# -------------------------


def detect_sample_rows(df: pd.DataFrame, metadata_cols_hint: int = 0) -> int:
    if metadata_cols_hint and metadata_cols_hint > 0:
        return min(metadata_cols_hint, df.shape[1] - 1)
    sample = df.head(min(len(df), 50))
    numeric = [is_numberish_series(sample.iloc[:, j]) for j in range(df.shape[1])]
    first_numeric = next((j for j, v in enumerate(numeric) if v), None)
    if first_numeric is None:
        return max(1, min(10, df.shape[1] - 1))
    return max(1, first_numeric)


# -------------------------
# representative protein FASTA
# -------------------------


def _infer_gene_from_seqid(seqid: str) -> str:
    s = _strip_refgen_suffix(seqid)
    s = re.sub(r"\.\w+$", "", s)  # drop extension-like
    s = _strip_isoform_suffix(s)
    return s


def _extract_kv_from_header(header: str) -> dict[str, str]:
    d = {}
    for tok in header.split():
        if "=" in tok:
            k, v = tok.split("=", 1)
            d[k.strip()] = v.strip()
    return d


def _score_first_transcript(transcript: str) -> tuple[int, int]:
    # lower is better
    m = re.search(r"_T(\d{3})$", transcript)
    if m:
        return (0, int(m.group(1)))
    return (1, 10**9)


def export_representative_protein_fasta(
    in_fasta: Path,
    out_fasta: Path,
    out_map_tsv: Path,
    isoform_mode: str,
    header_mode: str,
    gene_to_rep_transcript: dict[str, str],
    strip_version: bool,
):
    """Export one representative protein per gene.

    isoform_mode:
      - gff_longest: prefer transcript tagged longest=1 in GFF3
      - longest: longest protein sequence per gene
      - first: smallest transcript index per gene

    header_mode:
      - gene
      - gene|transcript
      - gene_meta (gene plus transcript/protein fields)
    """

    isoform_mode = (isoform_mode or "gff_longest").strip() or "gff_longest"
    header_mode = (header_mode or "gene").strip() or "gene"

    # If GFF-based is requested but we have no mapping, fall back
    if isoform_mode == "gff_longest" and not gene_to_rep_transcript:
        isoform_mode = "longest"

    def normalize_id(s: str) -> str:
        s2 = _strip_refgen_suffix(str(s or "").strip())
        if strip_version:
            s2 = _strip_version_suffix(s2)
        return s2

    # Parser: stream records
    def iter_fasta_records(fp: Path):
        header = None
        seq_lines = []
        with fp.open("r", encoding="utf-8", errors="ignore") as f:
            for ln in f:
                if ln.startswith(">"):
                    if header is not None:
                        yield header, "".join(seq_lines)
                    header = ln[1:].strip()
                    seq_lines = []
                else:
                    seq_lines.append(ln.strip())
            if header is not None:
                yield header, "".join(seq_lines)

    written = set()
    map_rows = []

    if isoform_mode == "gff_longest":
        # write on-the-fly (memory friendly)
        want = {g: normalize_id(t) for g, t in gene_to_rep_transcript.items()}
        with out_fasta.open("w", encoding="utf-8") as out:
            for header, seq in iter_fasta_records(in_fasta):
                if not seq:
                    continue
                seqid = header.split()[0]
                kv = _extract_kv_from_header(header)
                gene = normalize_id(kv.get("locus") or kv.get("gene") or _infer_gene_from_seqid(seqid))
                tr = normalize_id(kv.get("transcript") or "")
                prot = normalize_id(seqid)

                if gene in written:
                    continue
                if gene not in want:
                    continue
                if tr and tr != want[gene]:
                    continue

                if header_mode == "gene|transcript" and tr:
                    new_id = f"{gene}|{tr}"
                    new_header = new_id
                elif header_mode == "gene_meta":
                    new_header = f"{gene} transcript={tr or '-'} protein={prot}"
                else:
                    new_header = gene

                out.write(f">{new_header}\n")
                for i in range(0, len(seq), 60):
                    out.write(seq[i : i + 60] + "\n")

                written.add(gene)
                map_rows.append(
                    {
                        "gene_id": gene,
                        "transcript_id": tr or "",
                        "protein_id": prot,
                        "orig_header": header,
                        "mode": "gff_longest",
                    }
                )

        # genes in want but not written are missing
        missing = sorted([g for g in want.keys() if g not in written])
        return map_rows, missing

    # else: keep best per gene
    best = {}

    for header, seq in iter_fasta_records(in_fasta):
        if not seq:
            continue
        seqid = header.split()[0]
        kv = _extract_kv_from_header(header)
        gene = normalize_id(kv.get("locus") or kv.get("gene") or _infer_gene_from_seqid(seqid))
        tr = normalize_id(kv.get("transcript") or "")
        prot = normalize_id(seqid)

        if isoform_mode == "first":
            score = (_score_first_transcript(tr), -len(seq))
        else:  # longest
            score = (len(seq),)

        if gene not in best or score > best[gene][0]:
            best[gene] = (score, header, seq, tr, prot)

    with out_fasta.open("w", encoding="utf-8") as out:
        for gene in sorted(best.keys()):
            _, header, seq, tr, prot = best[gene]
            if header_mode == "gene|transcript" and tr:
                new_header = f"{gene}|{tr}"
            elif header_mode == "gene_meta":
                new_header = f"{gene} transcript={tr or '-'} protein={prot}"
            else:
                new_header = gene
            out.write(f">{new_header}\n")
            for i in range(0, len(seq), 60):
                out.write(seq[i : i + 60] + "\n")

            map_rows.append(
                {
                    "gene_id": gene,
                    "transcript_id": tr or "",
                    "protein_id": prot,
                    "orig_header": header,
                    "mode": isoform_mode,
                }
            )

    # no missing list for non-gff modes
    return map_rows, []


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--params", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    params = json.loads(Path(args.params).read_text(encoding="utf-8"))

    input_path = Path(str(params.get("input_path", "")).strip())
    if not input_path.exists():
        raise SystemExit("input_path not found")

    sheet_name = str(params.get("sheet_name", "") or "")
    orientation = str(params.get("orientation", "auto") or "auto")
    metadata_cols_hint = int(params.get("metadata_cols", 0) or 0)
    sample_id_col = str(params.get("sample_id_col", "") or "")

    # ID mapping controls
    target_mode = str(params.get("target_gene_id", "gff_gene_name") or "gff_gene_name")
    gff_gene_field = str(params.get("gff_gene_field", "Name") or "Name")
    strip_isoform = as_bool(params.get("strip_isoform"), True)
    strip_version = as_bool(params.get("strip_version"), True)
    sum_duplicate_genes = as_bool(params.get("sum_duplicate_genes"), True)
    force_integer = as_bool(params.get("force_integer"), True)

    alias_path = Path(str(params.get("alias_tsv", "")).strip()) if str(params.get("alias_tsv", "")).strip() else None
    alias_format = str(params.get("alias_format", "auto") or "auto")
    alias_target_assembly = str(params.get("alias_target_assembly", "auto (from GFF3)") or "auto (from GFF3)")
    custom_target_regex = str(params.get("custom_target_regex", "") or "").strip()
    gff_path = Path(str(params.get("gff3_path", "")).strip()) if str(params.get("gff3_path", "")).strip() else None
    fa_path = Path(str(params.get("fasta_path", "")).strip()) if str(params.get("fasta_path", "")).strip() else None

    # representative FASTA export controls
    export_rep_fasta = as_bool(params.get("export_rep_fasta"), False)
    isoform_mode = str(params.get("isoform_mode", "gff_longest") or "gff_longest")
    fasta_header_mode = str(params.get("fasta_header_mode", "gene") or "gene")

    df = read_table_auto(input_path, sheet_name=sheet_name)
    df = df.copy()
    df.columns = [str(c) for c in df.columns]

    # Build GFF3 index (gene-level)
    gene_names, gff_any_to_gene, gene_to_rep_tr = parse_gff3_gene_index(
        gff_path,
        gene_field=gff_gene_field,
        strip_isoform=strip_isoform,
        strip_version=strip_version,
    )

    valid_fa = parse_fasta_id_tokens(fa_path, strip_isoform=strip_isoform, strip_version=strip_version)
    alias_map, alias_rows = build_alias_map(
        alias_path,
        target_mode=target_mode,
        gff_any_to_gene=gff_any_to_gene,
        valid_fa=valid_fa,
        strip_isoform=strip_isoform,
        strip_version=strip_version,
        alias_format=alias_format,
        alias_target_assembly=alias_target_assembly,
        custom_target_regex=custom_target_regex,
    )

    alias_ambiguity = []
    if "__AMBIGUOUS__" in alias_map:
        try:
            alias_ambiguity = json.loads(alias_map.pop("__AMBIGUOUS__"))
        except Exception:
            alias_ambiguity = []


    # Orientation
    if orientation == "auto":
        meta_cols = detect_sample_rows(df, metadata_cols_hint=metadata_cols_hint)
        orientation = "samples_rows"
    else:
        meta_cols = detect_sample_rows(df, metadata_cols_hint=metadata_cols_hint)

    if orientation == "samples_rows":
        meta_cols = max(1, min(meta_cols, df.shape[1] - 1))
        meta_df = df.iloc[:, :meta_cols].copy()
        expr_df = df.iloc[:, meta_cols:].copy()

        if sample_id_col and sample_id_col in meta_df.columns:
            sample_series = meta_df[sample_id_col].astype(str)
        else:
            sample_id_col = str(meta_df.columns[0])
            sample_series = meta_df.iloc[:, 0].astype(str)

        sample_ids = [str(x).strip() for x in sample_series.tolist()]

        # make unique but keep originals in design
        seen = Counter()
        unique_sample_ids = []
        for s in sample_ids:
            seen[s] += 1
            unique_sample_ids.append(s if seen[s] == 1 else f"{s}_{seen[s]}")

        expr_num = expr_df.apply(pd.to_numeric, errors="coerce").fillna(0.0)
        counts = expr_num.T
        counts.columns = unique_sample_ids
        design_out = meta_df.copy()
        design_out.insert(0, "sample", unique_sample_ids)
    else:
        # genes_rows: first column is gene_id, remaining columns are samples (best effort)
        gene_col = str(df.columns[0])
        counts = df.set_index(gene_col).apply(pd.to_numeric, errors="coerce").fillna(0.0)
        unique_sample_ids = [str(c) for c in counts.columns]
        sample_id_col = "sample"
        design_out = pd.DataFrame({"sample": unique_sample_ids})

    gene_map_rows = []
    resolved_ids = []
    source_counter = Counter()
    unmatched = []

    for raw in counts.index.astype(str):
        resolved, matched_token, source = resolve_gene_id(
            raw,
            alias_map=alias_map,
            gff_any_to_gene=gff_any_to_gene,
            valid_fa=valid_fa,
            target_mode=target_mode,
            strip_isoform=strip_isoform,
            strip_version=strip_version,
        )
        source_counter[source] += 1
        resolved_ids.append(resolved)
        gene_map_rows.append(
            {
                "raw_gene_id": raw,
                "resolved_gene_id": resolved,
                "matched_token": matched_token,
                "source": source,
            }
        )
        if source in {"normalized", "original", "fasta"}:
            unmatched.append(raw)

    counts.index = resolved_ids

    if sum_duplicate_genes:
        counts = counts.groupby(level=0).sum()
    if force_integer:
        counts = counts.round().astype(int)

    counts_out = counts.reset_index().rename(columns={"index": "gene_id"})

    counts_path = out_dir / "counts.tsv"
    design_path = out_dir / "design.tsv"
    gene_map_path = out_dir / "gene_id_map.tsv"
    unmatched_path = out_dir / "unmatched_gene_ids.tsv"
    summary_path = out_dir / "summary.tsv"
    artifacts_path = out_dir / "artifacts.json"

    counts_out.to_csv(counts_path, sep="\t", index=False)
    design_out.to_csv(design_path, sep="\t", index=False)
    pd.DataFrame(gene_map_rows).to_csv(gene_map_path, sep="\t", index=False)
    pd.DataFrame({"raw_gene_id": unmatched}).to_csv(unmatched_path, sep="\t", index=False)


    # If MaizeDB alias table caused ambiguous primary->target mappings, write them out
    alias_ambig_path = out_dir / "alias_ambiguity.tsv"
    if alias_ambiguity:
        try:
            pd.DataFrame(alias_ambiguity).to_csv(alias_ambig_path, sep="\t", index=False)
        except Exception:
            pass

    # Optional: representative protein FASTA export
    rep_fa_path = None
    rep_map_path = None
    rep_missing_path = None
    rep_missing = []
    rep_written_count = 0

    if export_rep_fasta and fa_path and fa_path.exists():
        rep_fa_path = out_dir / "proteins_rep.fa"
        rep_map_path = out_dir / "proteins_rep_map.tsv"
        rep_missing_path = out_dir / "proteins_rep_missing_genes.tsv"
        map_rows, rep_missing = export_representative_protein_fasta(
            in_fasta=fa_path,
            out_fasta=rep_fa_path,
            out_map_tsv=rep_map_path,
            isoform_mode=isoform_mode,
            header_mode=fasta_header_mode,
            gene_to_rep_transcript=gene_to_rep_tr,
            strip_version=strip_version,
        )
        rep_written_count = len(map_rows)
        pd.DataFrame(map_rows).to_csv(rep_map_path, sep="\t", index=False)
        pd.DataFrame({"gene_id": rep_missing}).to_csv(rep_missing_path, sep="\t", index=False)

    summary = pd.DataFrame(
        [
            {"metric": "input_rows", "value": int(df.shape[0])},
            {"metric": "input_cols", "value": int(df.shape[1])},
            {"metric": "orientation_used", "value": orientation},
            {"metric": "metadata_cols", "value": int(meta_cols)},
            {"metric": "samples", "value": int(len(unique_sample_ids))},
            {"metric": "genes_output", "value": int(counts.shape[0])},
            {"metric": "alias_rows", "value": int(alias_rows)},
            {"metric": "alias_format", "value": alias_format},
            {"metric": "alias_target_assembly", "value": alias_target_assembly},
            {"metric": "custom_target_regex", "value": custom_target_regex or ""},
            {"metric": "alias_ambiguity_groups", "value": int(len(alias_ambiguity))},
            {"metric": "gff_genes_indexed", "value": int(len(gene_names))},
            {"metric": "gff_rep_transcripts", "value": int(len(gene_to_rep_tr))},
            {"metric": "resolved_alias_to_gff_gene", "value": int(source_counter.get("alias->gff_gene", 0))},
            {"metric": "resolved_by_alias", "value": int(source_counter.get("alias", 0))},
            {"metric": "resolved_by_gff_gene", "value": int(source_counter.get("gff_gene", 0))},
            {"metric": "resolved_by_fasta", "value": int(source_counter.get("fasta", 0))},
            {"metric": "resolved_by_normalized", "value": int(source_counter.get("normalized", 0))},
            {"metric": "resolved_by_original", "value": int(source_counter.get("original", 0))},
            {"metric": "design_sample_col_used", "value": sample_id_col},
            {"metric": "target_gene_id", "value": target_mode},
            {"metric": "gff_gene_field", "value": gff_gene_field},
            {"metric": "export_rep_fasta", "value": bool(export_rep_fasta)},
            {"metric": "isoform_mode", "value": isoform_mode},
            {"metric": "fasta_header_mode", "value": fasta_header_mode},
            {"metric": "proteins_rep_written", "value": int(rep_written_count)},
            {"metric": "proteins_rep_missing_genes", "value": int(len(rep_missing))},
        ]
    )
    summary.to_csv(summary_path, sep="\t", index=False)

    artifacts = {
        "counts": str(counts_path),
        "design": str(design_path),
        "gene_id_map": str(gene_map_path),
        "unmatched_gene_ids": str(unmatched_path),
        "summary": str(summary_path),
    }
    if rep_fa_path is not None:
        artifacts["proteins_rep_fasta"] = str(rep_fa_path)
    if rep_map_path is not None:
        artifacts["proteins_rep_map"] = str(rep_map_path)
    if rep_missing_path is not None:
        artifacts["proteins_rep_missing_genes"] = str(rep_missing_path)

    artifacts_path.write_text(json.dumps(artifacts, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(artifacts, ensure_ascii=False))


if __name__ == "__main__":
    main()
