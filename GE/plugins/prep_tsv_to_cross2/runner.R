#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(jsonlite)
  library(data.table)
})

args <- commandArgs(trailingOnly=TRUE)
get_arg <- function(flag, default=NULL) {
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) return(args[i+1])
  default
}
params_path <- get_arg("--params")
out_dir <- get_arg("--out")
if (is.null(params_path) || is.null(out_dir)) stop("Usage: --params params.json --out out_dir")
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

p <- fromJSON(params_path)

gt_path <- p$genotype_tsv
if (is.null(gt_path) || !file.exists(gt_path)) stop("genotype_tsv not found")
mm_path <- if (!is.null(p$marker_map_tsv) && nzchar(p$marker_map_tsv) && file.exists(p$marker_map_tsv)) p$marker_map_tsv else NULL
ph_path <- if (!is.null(p$phenotype_tsv) && nzchar(p$phenotype_tsv) && file.exists(p$phenotype_tsv)) p$phenotype_tsv else NULL

crosstype <- "f2"
if (!is.null(p$cross_type) && nzchar(p$cross_type)) crosstype <- as.character(p$cross_type)
if (!is.null(p$crosstype) && nzchar(p$crosstype)) crosstype <- as.character(p$crosstype)

base_out <- out_dir
if (!is.null(p$out_dir_override) && nzchar(p$out_dir_override)) {
  base_out <- as.character(p$out_dir_override)
}
base_out <- normalizePath(base_out, winslash="/", mustWork=FALSE)
dir.create(base_out, recursive=TRUE, showWarnings=FALSE)

cat("[prep_tsv_to_cross2] genotype_tsv=", gt_path, "\n")
cat("[prep_tsv_to_cross2] marker_map_tsv=", ifelse(is.null(mm_path), "(none)", mm_path), "\n")
cat("[prep_tsv_to_cross2] phenotype_tsv=", ifelse(is.null(ph_path), "(none)", ph_path), "\n")
cat("[prep_tsv_to_cross2] crosstype=", crosstype, "\n")
cat("[prep_tsv_to_cross2] out_dir=", base_out, "\n")

# ---- read genotype.tsv (wide) ----
df <- fread(gt_path, sep="\t", header=TRUE, data.table=FALSE, check.names=FALSE)
if (ncol(df) < 2) stop("genotype.tsv must have at least 2 columns (id + markers)")
ids <- as.character(df[[1]])
markers <- colnames(df)[-1]
G <- df[, -1, drop=FALSE]
# Drop accidental marker column named 'id' (can appear as a duplicated column in some exports)
drop_id <- which(tolower(markers) == 'id')
if (length(drop_id) > 0) {
  markers <- markers[-drop_id]
  G <- G[, -drop_id, drop=FALSE]
}

# ---- normalize to qtl2 codes (diploid: 1/2/3) ----
# Accept 0/1/2 (dosage), 1/2/3 (qtl2-like), A/H/B, AA/AB/BB, etc.
#
# Heuristic for whether '0' should be treated as missing (common in some qtl exports).
# Default: if we ever see a 3 (or higher) in numeric genotypes, we treat 0 as missing;
# otherwise, we treat 0 as a valid genotype (dosage=0). You can override with params: zero_as_missing.
detect_zero_as_missing <- function(Gdf) {
  nr <- nrow(Gdf); nc <- ncol(Gdf)
  if (is.null(nr) || is.null(nc) || nr == 0 || nc == 0) return(FALSE)
  r_idx <- seq_len(min(nr, 200))
  c_idx <- seq_len(min(nc, 200))
  sub <- as.matrix(Gdf[r_idx, c_idx, drop=FALSE])
  z <- as.vector(sub)
  z <- z[!is.na(z)]
  if (length(z) == 0) return(FALSE)
  z <- trimws(as.character(z))
  z <- z[z != '' & z != '.' & toupper(z) != 'NA' & z != '-9']
  if (length(z) == 0) return(FALSE)
  iv <- suppressWarnings(as.integer(z))
  iv <- iv[!is.na(iv)]
  if (length(iv) == 0) return(FALSE)
  max(iv) >= 3L
}
zero_as_missing <- if (!is.null(p$zero_as_missing)) as.logical(p$zero_as_missing) else detect_zero_as_missing(G)
cat('[prep_tsv_to_cross2] zero_as_missing=', zero_as_missing, '\n')

normalize_cell <- function(v) {
  if (is.na(v) || v=="" || v=="." || v=="NA" || v=="NaN" || v=="-9") return(NA_integer_)
  vv <- as.character(v)
  vv <- trimws(vv)
  if (vv=="") return(NA_integer_)

  # ABH-style
  if (vv %in% c("A","a")) return(1L)
  if (vv %in% c("H","h")) return(2L)
  if (vv %in% c("B","b")) return(3L)

  # genotype strings
  if (vv %in% c("AA","aa")) return(1L)
  if (vv %in% c("AB","ab","BA","ba")) return(2L)
  if (vv %in% c("BB","bb")) return(3L)

  # numeric
  iv <- suppressWarnings(as.integer(vv))
  if (!is.na(iv)) {
    if (isTRUE(zero_as_missing) && iv == 0L) return(NA_integer_)
    if (!isTRUE(zero_as_missing) && iv %in% c(0L,1L,2L)) return(iv + 1L)
    if (iv %in% c(1L,2L,3L,4L,5L)) return(iv)
  }
  return(NA_integer_)
}

# Apply conversion column-wise (avoid huge apply on data.frame)
Gq <- matrix(NA_integer_, nrow=nrow(G), ncol=ncol(G))
colnames(Gq) <- markers
rownames(Gq) <- ids
for (j in seq_along(markers)) {
  v <- G[[j]]
  Gq[, j] <- vapply(v, normalize_cell, integer(1))
}

geno_csv <- file.path(base_out, "geno.csv")
# write with id column first
out_dt <- as.data.table(Gq)
out_dt[, id := ids]
setcolorder(out_dt, c("id", markers))
fwrite(out_dt, geno_csv, sep = ",", quote = F, na = "NA")

# ---- phenotype ----
pheno_csv <- file.path(base_out, "pheno.csv")
if (!is.null(ph_path)) {
  ph <- fread(ph_path, sep="\t", header=TRUE, data.table=FALSE, check.names=FALSE)
  if (ncol(ph) < 1) stop("phenotype.tsv has no columns")
  # ensure id column exists
  if (!("id" %in% colnames(ph))) {
    colnames(ph)[1] <- "id"
  }
  # keep only rows in ids order if possible
  if ("id" %in% colnames(ph)) {
    idx <- match(ids, as.character(ph[["id"]]))
    if (any(!is.na(idx))) ph <- ph[idx, , drop=FALSE]
  }
  fwrite(as.data.table(ph), pheno_csv, sep = ",", quote = F, na = "NA")
} else {
  fwrite(data.table(id=ids, pheno=seq_along(ids)), pheno_csv, sep = ",", quote = F, na = "NA")
}

# ---- marker map (pmap/gmap) ----
pmap_csv <- file.path(base_out, 'pmap.csv')
gmap_csv <- file.path(base_out, 'gmap.csv')

parse_chrpos <- function(x) {
  x2 <- as.character(x)
  chr <- rep(NA_character_, length(x2))
  pos <- rep(NA_real_, length(x2))
  # VCF-style 'chr:pos'
  m <- regexec('^([^:]+):([0-9]+)$', x2)
  r <- regmatches(x2, m)
  ok <- lengths(r) == 3
  if (any(ok)) {
    chr[ok] <- vapply(r[ok], function(z) z[2], character(1))
    pos[ok] <- suppressWarnings(as.numeric(vapply(r[ok], function(z) z[3], character(1))))
  }
  # Also allow 'chr_pos' or 'chr-pos' when chr part is numeric-ish
  m2 <- regexec('^([^_\\-]+)[_\\-]([0-9]+)$', x2)
  r2 <- regmatches(x2, m2)
  ok2 <- lengths(r2) == 3
  if (any(ok2)) {
    idx2 <- which(ok2)
    chr2 <- vapply(r2[ok2], function(z) z[2], character(1))
    pos2 <- suppressWarnings(as.numeric(vapply(r2[ok2], function(z) z[3], character(1))))
    chr2s <- sub('^(chr|Chr|CHR)', '', chr2)
    keep <- grepl('^[0-9]+$', chr2s)
    if (any(keep)) {
      chr[idx2[keep]] <- chr2[keep]
      pos[idx2[keep]] <- pos2[keep]
    }
  }
  list(chr=chr, pos=pos)
}
normalize_chr <- function(x) {
  # Spec: do not auto-normalize/convert chromosome labels.
  as.character(x)
}

make_map <- function(mm, use_pos_col) {
  # Spec: marker_map.tsv is fixed as: 1st=marker, 2nd=chr, 3rd=pos.
  # (use_pos_col is kept for compatibility but ignored.)
  if (is.null(mm) || nrow(mm) == 0) {
    out <- data.table(marker=markers, chr='1', pos=seq_along(markers))
    parsed <- parse_chrpos(markers)
    ok <- !is.na(parsed$chr) & !is.na(parsed$pos)
    if (any(ok)) {
      out$chr[ok] <- normalize_chr(parsed$chr[ok])
      out$pos[ok] <- parsed$pos[ok]
    }
    return(out)
  }
  mm <- as.data.table(mm)
  if (ncol(mm) < 3) {
    out <- data.table(marker=markers, chr='1', pos=seq_along(markers))
    return(out)
  }

  mm_marker <- trimws(as.character(mm[[1]]))
  mm_chr <- as.character(mm[[2]])
  pos_str <- gsub(",", "", trimws(as.character(mm[[3]])))
  mm_pos <- suppressWarnings(as.numeric(pos_str))
  if (all(is.na(mm_pos))) mm_pos <- seq_along(mm_marker)
  missp <- is.na(mm_pos)
  if (any(missp)) mm_pos[missp] <- seq_len(sum(missp))
  mm_chr[is.na(mm_chr) | mm_chr==''] <- '1'
  mdt <- data.table(marker=mm_marker, chr=mm_chr, pos=as.numeric(mm_pos))
  # align to genotype marker order
  idx <- match(markers, mdt$marker)
  ok <- !is.na(idx)
  out <- data.table(marker=markers, chr="1", pos=seq_along(markers))
  out$chr[ok] <- mdt$chr[idx[ok]]
  out$pos[ok] <- mdt$pos[idx[ok]]
  out$chr <- normalize_chr(out$chr)
  # fill missing chr/pos by parsing marker names when possible
  parsed <- parse_chrpos(out$marker)
  miss_chr <- is.na(out$chr) | out$chr==''
  if (any(miss_chr & !is.na(parsed$chr))) out$chr[miss_chr] <- normalize_chr(parsed$chr[miss_chr])
  miss_pos <- is.na(out$pos)
  if (any(miss_pos & !is.na(parsed$pos))) out$pos[miss_pos] <- parsed$pos[miss_pos]
  out$chr[is.na(out$chr) | out$chr==''] <- '1'
  out$pos[is.na(out$pos)] <- seq_along(out$pos)[is.na(out$pos)]
  out
}

if (!is.null(mm_path)) {
  mm <- fread(mm_path, sep="\t", header=TRUE, data.table=FALSE, check.names=FALSE)
  # Prefer explicit bp/cM columns when available
  pmap_dt <- make_map(mm, c("pos_bp","bp","POS","pos","position"))
  gmap_dt <- make_map(mm, c("pos_cM","cM","pos","position"))
} else {
  pmap_dt <- make_map(NULL, character(0))
  gmap_dt <- make_map(NULL, character(0))
}

# Optional convenience scaling for gmap: if pos looks bp-scale, convert to Mb (~cM)
gmap_max <- suppressWarnings(max(gmap_dt$pos, na.rm=TRUE))
if (is.finite(gmap_max) && gmap_max > 1e4) {
  cat('[prep_tsv_to_cross2] INFO: gmap pos seems bp-scale (max=', gmap_max, '); applying pos/1e6 (1Mb≈1cM)\n')
  gmap_dt$pos <- gmap_dt$pos / 1e6
}

fwrite(pmap_dt, pmap_csv, sep = ",", quote = F, na = "NA")
fwrite(gmap_dt, gmap_csv, sep = ",", quote = F, na = "NA")

# ---- cross2.yaml ----
cross2_yaml <- file.path(base_out, "cross2.yaml")
# minimal yaml
yaml_lines <- c(
  paste0("crosstype: ", crosstype),
  paste0("geno: ", basename(geno_csv)),
  paste0("pheno: ", basename(pheno_csv)),
  paste0("gmap: ", basename(gmap_csv)),
  paste0("pmap: ", basename(pmap_csv)),
  "genotypes:",
  "  1: 1",
  "  2: 2",
  "  3: 3",
  "  4: 4",
  "  5: 5",
  "alleles:",
  "- A",
  "- B",
  "x_chr: FALSE"
)
writeLines(yaml_lines, cross2_yaml)

meta <- list(
  cross2_yaml = cross2_yaml,
  geno_csv = geno_csv,
  pheno_csv = pheno_csv,
  pmap_csv = pmap_csv,
  gmap_csv = gmap_csv,
  crosstype = crosstype,
  genotype_tsv = gt_path,
  marker_map_tsv = ifelse(is.null(mm_path), "", mm_path),
  phenotype_tsv = ifelse(is.null(ph_path), "", ph_path),
  out_dir = base_out
)
writeLines(toJSON(meta, auto_unbox=TRUE, pretty=TRUE), con=file.path(out_dir, "artifacts.json"))
