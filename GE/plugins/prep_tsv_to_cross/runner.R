#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(jsonlite)
  library(data.table)
  library(qtl)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) return(args[i + 1])
  default
}

params_path <- get_arg("--params")
out_dir <- get_arg("--out")
if (is.null(params_path) || is.null(out_dir)) stop("Usage: --params params.json --out out_dir")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_dir, "run.log")
sink(log_file, split = TRUE)

cat("[prep_tsv_to_cross] start\n")
cat("[prep_tsv_to_cross] params_path=", params_path, "\n")
cat("[prep_tsv_to_cross] out_dir=", out_dir, "\n")

p <- fromJSON(params_path)

geno_tsv <- p$genotype_tsv
phe_tsv  <- if (!is.null(p$phenotype_tsv)) p$phenotype_tsv else NULL
marker_map_tsv <- if (!is.null(p$marker_map_tsv)) p$marker_map_tsv else NULL

cross_type <- if (!is.null(p$cross_type) && nchar(p$cross_type) > 0) p$cross_type else "f2"
geno_code  <- if (!is.null(p$geno_code) && nchar(p$geno_code) > 0) p$geno_code else "auto"

na_strings <- c("NA", ".", "-", "", "nan", "NaN")
if (!is.null(p$na_strings) && length(p$na_strings) > 0) {
  na_strings <- unique(c(na_strings, as.character(p$na_strings)))
}

if (is.null(geno_tsv) || !file.exists(geno_tsv)) stop("genotype_tsv is required and must exist")
cat("[prep_tsv_to_cross] geno_tsv=", geno_tsv, "\n")
cat("[prep_tsv_to_cross] phe_tsv=", phe_tsv, "\n")
cat("[prep_tsv_to_cross] marker_map_tsv=", marker_map_tsv, "\n")
cat("[prep_tsv_to_cross] cross_type=", cross_type, " geno_code=", geno_code, "\n")

# Read genotype.tsv while preserving the *exact* header marker names.
# In some environments, parsers may sanitize/modify column names; we force the header from the first line.
hdr <- tryCatch(readLines(geno_tsv, n=1, warn=FALSE), error=function(e) NULL)
hdr_fields <- NULL
if (!is.null(hdr) && length(hdr) == 1L) {
  hdr_fields <- strsplit(hdr, "\t", fixed=TRUE)[[1]]
}
gt <- fread(geno_tsv, sep="\t", header=TRUE, data.table=FALSE, fill=TRUE, quote="")
if (!is.null(hdr_fields) && length(hdr_fields) == ncol(gt)) {
  colnames(gt) <- hdr_fields
}
#if (!("id" %in% colnames(gt))) stop("genotype_tsv must have 'id' column")
#gt$id <- as.character(gt$id)
gt$id <- as.character(gt[, 1])

#markers <- setdiff(colnames(gt), "id")
markers <- colnames(gt)[2:ncol(gt)]
# Drop accidental marker column named 'id' (can appear as a duplicated column in some exports)
markers <- markers[tolower(markers) != 'id']
if (length(markers) == 0) stop("No marker columns found in genotype_tsv")

# phenotype
ph <- NULL
if (!is.null(phe_tsv) && file.exists(phe_tsv)) {
  ph <- fread(phe_tsv, sep="\t", header=TRUE, data.table=FALSE)
  #if (!("id" %in% colnames(ph))) stop("phenotype_tsv must have 'id' column")
  #ph$id <- as.character(ph$id)
  # merge to keep genotype order
  ph <- ph[match(gt$id, ph[, 1]), , drop=FALSE]
} else {
  # NOTE: R/qtl requires at least one phenotype column with at least one
  # non-missing value. We provide a simple index phenotype by default.
  ph <- data.frame(id=gt$id, dummy=seq_along(gt$id))
}

# Extract genotype matrix (as character first)
Xc <- as.matrix(gt[, markers, drop=FALSE])
Xc <- apply(Xc, 2, as.character)

# Normalize NAs
is_na_token <- function(z) {
  z2 <- trimws(z)
  z2 %in% na_strings
}

# Determine encoding
guess_numeric012 <- function(col) {
  z <- trimws(col)
  z[is_na_token(z)] <- NA
  z_num <- suppressWarnings(as.numeric(z))
  ok <- !is.na(z_num)
  if (!any(ok)) return(FALSE)
  all(unique(z_num[ok]) %in% c(0,1,2))
}

guess_AHB <- function(col) {
  z <- toupper(trimws(col))
  z[is_na_token(z)] <- NA
  ok <- !is.na(z)
  if (!any(ok)) return(FALSE)
  all(unique(z[ok]) %in% c("A","H","B"))
}

guess_AAABBB <- function(col) {
  z <- toupper(trimws(col))
  z[is_na_token(z)] <- NA
  ok <- !is.na(z)
  if (!any(ok)) return(FALSE)
  all(unique(z[ok]) %in% c("AA","AB","BB"))
}

if (geno_code == "auto") {
  # use first non-empty column to guess
  idx <- which(colSums(!is.na(Xc) & !(trimws(Xc) %in% na_strings)) > 0)
  if (length(idx) == 0) stop("All genotypes are missing")
  col0 <- Xc[, idx[1]]
  if (guess_numeric012(col0)) geno_code <- "012"
  else if (guess_AHB(col0)) geno_code <- "AHB"
  else if (guess_AAABBB(col0)) geno_code <- "AAABBB"
  else geno_code <- "AAABBB"  # fallback: treat as genotype strings
  cat("[prep_tsv_to_cross] auto-detected geno_code=", geno_code, "\n")
}

# Convert to qtl genotype strings (AA/AB/BB)
to_AAABBB <- function(col) {
  z <- trimws(col)
  z[is_na_token(z)] <- NA
  if (geno_code == "012") {
    zn <- suppressWarnings(as.numeric(z))
    out <- rep(NA_character_, length(zn))
    out[!is.na(zn) & zn==0] <- "AA"
    out[!is.na(zn) & zn==1] <- "AB"
    out[!is.na(zn) & zn==2] <- "BB"
    return(out)
  }
  if (geno_code == "AHB") {
    z2 <- toupper(z)
    out <- rep(NA_character_, length(z2))
    out[!is.na(z2) & z2=="A"] <- "AA"
    out[!is.na(z2) & z2=="H"] <- "AB"
    out[!is.na(z2) & z2=="B"] <- "BB"
    return(out)
  }
  # already AA/AB/BB (or similar)
  z2 <- toupper(z)
  out <- rep(NA_character_, length(z2))
  out[!is.na(z2) & z2 %in% c("AA","AB","BB")] <- z2[!is.na(z2) & z2 %in% c("AA","AB","BB")]
  return(out)
}

Xg <- apply(Xc, 2, to_AAABBB)


# ----- Fix3: build cross directly (no read.cross) -----

# Map AA/AB/BB strings to numeric codes (R/qtl conventions): 1=AA, 2=AB, 3=BB
Xcode <- apply(Xg, 2, function(v) {
  out <- rep(NA_integer_, length(v))
  out[v=="AA"] <- 1L
  out[v=="AB"] <- 2L
  out[v=="BB"] <- 3L
  out
})
storage.mode(Xcode) <- "integer"
rownames(Xcode) <- gt$id
colnames(Xcode) <- markers

ct <- tolower(cross_type)
if (ct %in% c("bc","bc1","bc2")) {
  n_bb <- sum(Xcode == 3L, na.rm=TRUE)
  if (n_bb > 0) cat("[prep_tsv_to_cross] warning: found BB genotypes in backcross; set to NA: ", n_bb, "\n")
  Xcode[Xcode == 3L] <- NA_integer_
}
if (ct %in% c("riself","risib","ril")) {
  n_ab <- sum(Xcode == 2L, na.rm=TRUE)
  if (n_ab > 0) cat("[prep_tsv_to_cross] warning: found AB genotypes in RIL; set to NA: ", n_ab, "\n")
  Xcode[Xcode == 2L] <- NA_integer_
}

# Build genotype list; default puts all markers into chr "1". If marker_map_tsv is provided
# with columns marker/chr/pos, we split markers by chr and set map positions to pos.
parse_chrpos <- function(x) {
  x2 <- as.character(x)
  chr <- rep(NA_character_, length(x2))
  pos <- rep(NA_real_, length(x2))
  m <- regexec('^([^:]+):([0-9]+)$', x2)
  r <- regmatches(x2, m)
  ok <- lengths(r) == 3
  if (any(ok)) {
    chr[ok] <- vapply(r[ok], function(z) z[2], character(1))
    pos[ok] <- suppressWarnings(as.numeric(vapply(r[ok], function(z) z[3], character(1))))
  }
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

# Build genotype list; if marker_map_tsv is provided with columns marker/chr/pos,
# we split markers by chr and set map positions to pos. If marker_map_tsv is absent,
# we infer chr/pos from marker ids when possible (e.g. VCF 'chr:pos').
build_geno_list <- function(mat, marker_map_tsv=NULL) {
  mk <- colnames(mat)

  # 1) marker_map.tsv
  if (!is.null(marker_map_tsv) && file.exists(marker_map_tsv)) {
    # Spec: marker_map.tsv is fixed as: 1st=marker, 2nd=chr, 3rd=pos.
    # marker/chr are kept as-is (string). pos is numeric; if non-numeric then sequential.
    mp0 <- tryCatch(fread(marker_map_tsv, sep='\t', header=TRUE, data.table=FALSE, fill=TRUE, quote=""), error=function(e) NULL)
    if (!is.null(mp0) && ncol(mp0) >= 3) {
      mp <- mp0[, 1:3, drop=FALSE]
      colnames(mp) <- c('marker','chr','pos')

      # robust match with optional trimming of spaces
      mk_trim <- trimws(as.character(mk))
      if (any(duplicated(mk_trim))) {
        stop("Duplicate marker names after trimming spaces in genotype.tsv header; cannot match marker_map.tsv")
      }
      mk_dict <- setNames(mk, mk_trim)

      mp$marker <- trimws(as.character(mp$marker))
      mp$chr <- as.character(mp$chr)

      # pos: numeric parse; if non-numeric then NA (later replaced)
      pos_str <- gsub(",", "", trimws(as.character(mp$pos)))
      mp$pos <- suppressWarnings(as.numeric(pos_str))

      # keep only markers that exist in genotype (by trimmed match)
      in_gt <- mp$marker %in% names(mk_dict)
      # If nothing matches, try the common R 'make.names' transformation (e.g. "1" -> "X1").
      # We do NOT change the genotype marker names; this is only for matching.
      if (!any(in_gt)) {
        mp_mk2 <- make.names(mp$marker, unique=FALSE)
        in_gt2 <- mp_mk2 %in% names(mk_dict)
        if (any(in_gt2)) {
          mp$marker[in_gt2] <- mp_mk2[in_gt2]
          in_gt <- in_gt2
        }
      }
      mp <- mp[in_gt, , drop=FALSE]
      if (nrow(mp) > 0) {
        # Replace marker names with the exact genotype header names
        mp$marker <- unname(mk_dict[mp$marker])

        # chr: if missing/blank, fallback to '1' (only when absent)
        mp$chr[is.na(mp$chr) | mp$chr==''] <- '1'
        mp$chr <- normalize_chr(mp$chr)

        # pos: if entirely non-numeric, use sequential per chr; otherwise fill only missing
        if (all(is.na(mp$pos))) {
          mp$pos <- ave(seq_len(nrow(mp)), mp$chr, FUN=function(x) seq_along(x))
        } else {
          for (cc in unique(mp$chr)) {
            ii <- which(mp$chr == cc)
            miss <- is.na(mp$pos[ii])
            if (any(miss)) {
              base <- suppressWarnings(max(mp$pos[ii], na.rm=TRUE))
              if (!is.finite(base)) base <- 0
              mp$pos[ii][miss] <- base + seq_len(sum(miss))
            }
          }
        }

        # Optional convenience scaling: if positions look like bp, convert to Mb (= approx cM)
        mp_max <- suppressWarnings(max(mp$pos, na.rm=TRUE))
        if (is.finite(mp_max) && mp_max > 1e4) {
          cat("[prep_tsv_to_cross] INFO: marker_map pos seems bp-scale (max=", mp_max, "); applying pos/1e6 (1Mb≈1cM)\n")
          mp$pos <- mp$pos / 1e6
        }

        g <- list()
        chr_order <- unique(mp$chr)
        for (cc in chr_order) {
          sub <- mp[mp$chr == cc, , drop=FALSE]
          o <- order(sub$pos, na.last=TRUE)
          mks <- sub$marker[o]
          submat <- mat[, mks, drop=FALSE]
          g[[as.character(cc)]] <- list(
            data = submat,
            map = setNames(as.numeric(sub$pos[o]), mks),
            alleles = c('A','B')
          )
        }
        rest <- setdiff(mk, mp$marker)
        if (length(rest) > 0) {
          submat <- mat[, rest, drop=FALSE]
          g[['un']] <- list(data=submat, map=setNames(seq_len(ncol(submat)), colnames(submat)), alleles=c('A','B'))
        }
        return(g)
      }
    }
  }

  # 2) infer from marker ids
  parsed <- parse_chrpos(mk)
  ok <- !is.na(parsed$chr) & !is.na(parsed$pos)
  if (any(ok)) {
    mp <- data.frame(marker=mk[ok], chr=normalize_chr(parsed$chr[ok]), pos=parsed$pos[ok], stringsAsFactors=FALSE)
    g <- list()
    for (cc in unique(mp$chr)) {
      sub <- mp[mp$chr == cc, , drop=FALSE]
      o <- order(sub$pos)
      mks <- sub$marker[o]
      submat <- mat[, mks, drop=FALSE]
      g[[as.character(cc)]] <- list(data=submat, map=setNames(sub$pos[o], mks), alleles=c('A','B'))
    }
    rest <- setdiff(mk, mp$marker)
    if (length(rest) > 0) {
      submat <- mat[, rest, drop=FALSE]
      g[['un']] <- list(data=submat, map=setNames(seq_len(ncol(submat)), colnames(submat)), alleles=c('A','B'))
    }
    return(g)
  }

  # 3) default
  list('1' = list(data = mat, map = setNames(seq_len(ncol(mat)), colnames(mat)), alleles = c('A','B')))
}

# pheno: rownames must match individual ids
ph_df <- as.data.frame(ph)
rownames(ph_df) <- ph_df$id
ph_df$id <- NULL

# ASMap::mstmap expects a unique genotype identifier column in cross$pheno.
# The default is "Genotype" (see mstmap.cross(id=...)).
if (!("Genotype" %in% colnames(ph_df))) {
  ph_df$Genotype <- rownames(ph_df)
}

cross <- list(
  pheno = ph_df,
  geno = build_geno_list(Xcode, marker_map_tsv)
)
# IMPORTANT: R/qtl/ASMap expect the cross-type class to be the first element,
# e.g. c("bc","cross"), not c("cross","bc"). ASMap::mstmap.cross branches
# on class(object)[1] and will fail if "cross" is first.
class(cross) <- c(cross_type, "cross")

# Summaries from Xcode (dose 0/1/2)
dose <- Xcode
# 1/2/3 -> 0/1/2
suppressWarnings({
  dose[dose == 1L] <- 0L
  dose[dose == 2L] <- 1L
  dose[dose == 3L] <- 2L
})

miss_by_ind <- rowMeans(is.na(dose))
sample_dt <- data.table(id = rownames(dose), missing_rate = miss_by_ind)
fwrite(sample_dt, file.path(out_dir, "sample_summary.tsv"), sep="\t")

miss_by_mar <- colMeans(is.na(dose))
pA <- colMeans(dose == 0L, na.rm=TRUE) + 0.5 * colMeans(dose == 1L, na.rm=TRUE)
maf <- pmin(pA, 1 - pA)
mar_dt <- data.table(
  marker = colnames(dose),
  missing_rate = miss_by_mar,
  maf = maf,
  n_nonmissing = colSums(!is.na(dose))
)
setorder(mar_dt, missing_rate, -maf)
fwrite(mar_dt, file.path(out_dir, "marker_summary.tsv"), sep="\t")

# Save cross
out_rds <- file.path(out_dir, "cross.rds")
saveRDS(cross, out_rds)
cat("[prep_tsv_to_cross] saved cross=", out_rds, "\n")

art <- list(
  cross_rds = out_rds,
  sample_summary_tsv = file.path(out_dir, "sample_summary.tsv"),
  marker_summary_tsv = file.path(out_dir, "marker_summary.tsv")
)
write(toJSON(art, auto_unbox=TRUE, pretty=TRUE), file=file.path(out_dir, "artifacts.json"))

cat("[prep_tsv_to_cross] done\n")
sink()
