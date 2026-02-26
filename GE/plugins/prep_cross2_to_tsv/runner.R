#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(jsonlite)
  library(data.table)
  library(qtl2)
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

cross2_path <- p$cross2_path
if (is.null(cross2_path) || !file.exists(cross2_path)) stop("cross2_path not found")

base_out <- out_dir
if (!is.null(p$out_dir_override) && nzchar(p$out_dir_override)) base_out <- as.character(p$out_dir_override)
base_out <- normalizePath(base_out, winslash="/", mustWork=FALSE)
dir.create(base_out, recursive=TRUE, showWarnings=FALSE)

cat("[prep_cross2_to_tsv] cross2_path=", cross2_path, "\n")
cat("[prep_cross2_to_tsv] out_dir=", base_out, "\n")

# Read cross2
cross2 <- NULL
if (grepl("\\.rds$", cross2_path, ignore.case=TRUE)) {
  cross2 <- readRDS(cross2_path)
} else {
  cross2 <- qtl2::read_cross2(cross2_path)
}

if (is.null(cross2$geno) || length(cross2$geno) == 0) stop("cross2 has no geno")

# Combine geno matrices across chromosomes
mats <- lapply(cross2$geno, function(m) {
  if (is.data.frame(m)) m <- as.matrix(m)
  if (!is.matrix(m)) m <- as.matrix(m)
  m
})

# Ensure consistent rownames
ids <- rownames(mats[[1]])
if (is.null(ids)) ids <- as.character(seq_len(nrow(mats[[1]])))

# cbind
X <- do.call(cbind, mats)
rownames(X) <- ids

# Missing might be 0 in some exports
X[X == 0] <- NA

u <- sort(unique(as.integer(na.omit(as.vector(X)))))
maxv <- if (length(u)==0) NA_integer_ else max(u)

X_dose <- X
if (!is.na(maxv)) {
  if (maxv <= 3L) {
    # 1/2/3 -> 0/1/2
    if (all(u %in% c(1L,2L,3L))) {
      X_dose <- X - 1L
    } else if (all(u %in% c(1L,2L))) {
      X_dose <- (X - 1L) * 2L
    } else {
      X_dose <- X
    }
  } else {
    # Polyploid / unknown coding: keep numeric as-is (with missing already NA)
    X_dose <- X
  }
}

markers <- colnames(X_dose)
if (is.null(markers)) markers <- paste0("m", seq_len(ncol(X_dose)))

# Write genotype.tsv
out_geno <- file.path(base_out, "genotype.tsv")
DT <- as.data.table(X_dose)
setnames(DT, markers)
DT[, id := ids]
setcolorder(DT, c("id", markers))
fwrite(DT, out_geno, sep="\t", quote=FALSE, na="NA")

# marker_map.tsv (SPEC: 1st=marker, 2nd=chr, 3rd=pos)
out_map <- file.path(base_out, "marker_map.tsv")

# In qtl2, pmap/gmap are often *named numeric vectors* (names are markers, values are positions).
# The previous implementation treated them as tables and lost marker names, causing chr=1 / pos=seq fallback.
extract_pos <- function(map_list, chr, markers_chr) {
  if (is.null(map_list) || length(map_list) == 0) return(rep(NA_real_, length(markers_chr)))
  chr <- as.character(chr)
  if (!(chr %in% names(map_list))) return(rep(NA_real_, length(markers_chr)))
  x <- map_list[[chr]]

  # Named atomic vector (most common)
  if (is.atomic(x) && !is.null(names(x))) {
    nm <- names(x)
    vv <- suppressWarnings(as.numeric(x))
    idx <- match(markers_chr, nm)
    out <- rep(NA_real_, length(markers_chr))
    ok <- !is.na(idx)
    out[ok] <- vv[idx[ok]]
    return(out)
  }

  # Data frame / matrix map
  if (is.data.frame(x) || is.matrix(x)) {
    df <- as.data.frame(x, stringsAsFactors=FALSE)
    mk <- NULL
    if ("marker" %in% names(df)) {
      mk <- df$marker
    } else if (!is.null(rownames(df)) && any(nzchar(rownames(df)))) {
      mk <- rownames(df)
    } else if (ncol(df) >= 1) {
      mk <- df[[1]]
    } else {
      mk <- rep(NA_character_, nrow(df))
    }

    pos <- NULL
    if ("pos" %in% names(df)) pos <- df$pos
    else if ("position" %in% names(df)) pos <- df$position
    else if (ncol(df) >= 2) pos <- df[[2]]
    else pos <- rep(NA, nrow(df))

    mk <- as.character(mk)
    pos_num <- suppressWarnings(as.numeric(pos))
    idx <- match(markers_chr, mk)
    out <- rep(NA_real_, length(markers_chr))
    ok <- !is.na(idx)
    out[ok] <- pos_num[idx[ok]]
    return(out)
  }

  # Unnamed vector but same length
  if (is.atomic(x) && length(x) == length(markers_chr)) {
    return(suppressWarnings(as.numeric(x)))
  }

  rep(NA_real_, length(markers_chr))
}

# Build mapping marker -> (chr,pos) directly from cross2$geno chromosome structure
rows <- list()
for (chr in names(cross2$geno)) {
  m <- cross2$geno[[chr]]
  if (is.data.frame(m)) m <- as.matrix(m)
  if (!is.matrix(m)) m <- as.matrix(m)
  mks <- colnames(m)
  if (is.null(mks)) mks <- paste0(as.character(chr), "_m", seq_len(ncol(m)))

  pos <- extract_pos(cross2$pmap, chr, mks)
  if (all(is.na(pos))) pos <- extract_pos(cross2$gmap, chr, mks)

  # If positions are missing/non-numeric, assign sequential (only then)
  if (all(is.na(pos))) {
    pos <- seq_along(mks)
  } else {
    miss <- is.na(pos)
    if (any(miss)) {
      basev <- suppressWarnings(max(pos, na.rm=TRUE))
      if (!is.finite(basev)) basev <- 0
      pos[miss] <- basev + seq_len(sum(miss))
    }
  }

  rows[[as.character(chr)]] <- data.table(
    marker = as.character(mks),
    chr    = as.character(chr),
    pos    = as.numeric(pos)
  )
}

mm <- rbindlist(rows, use.names=TRUE, fill=TRUE)

# Align to genotype marker order
mm_key <- mm[!duplicated(marker)]
idx <- match(markers, mm_key$marker)

# If no match (rare), try R make.names only for matching (do not change marker strings)
if (any(is.na(idx))) {
  mk2 <- make.names(markers, unique=FALSE)
  mm2 <- make.names(mm_key$marker, unique=FALSE)
  idx2 <- match(mk2, mm2)
  idx[is.na(idx)] <- idx2[is.na(idx)]
}

out_mm <- data.table(marker=markers, chr=NA_character_, pos=as.numeric(NA))
ok <- !is.na(idx)
if (any(ok)) {
  out_mm$chr[ok] <- as.character(mm_key$chr[idx[ok]])
  out_mm$pos[ok] <- as.numeric(mm_key$pos[idx[ok]])
}

# Fallback only for missing entries (should be rare)
miss <- is.na(out_mm$chr) | out_mm$chr == ""
if (any(miss)) out_mm$chr[miss] <- "1"

badpos <- is.na(out_mm$pos)
if (any(badpos)) out_mm$pos[badpos] <- seq_along(out_mm$pos)[badpos]

fwrite(out_mm, out_map, sep="\t", quote=FALSE, na="NA")

# phenotype.tsv
out_ph <- file.path(base_out, "phenotype.tsv")
ph <- cross2$pheno
if (is.null(ph) || nrow(ph) == 0) {
  ph_dt <- data.table(id=ids, pheno=seq_along(ids))
} else {
  ph_dt <- as.data.table(ph)
  if (!("id" %in% names(ph_dt))) {
    rid <- rownames(ph)
    if (!is.null(rid)) {
      ph_dt[, id := as.character(rid)]
    } else {
      ph_dt[, id := ids]
    }
  }
  setcolorder(ph_dt, c("id", setdiff(names(ph_dt), "id")))
  idx <- match(ids, as.character(ph_dt$id))
  if (any(!is.na(idx))) ph_dt <- ph_dt[idx, , drop=FALSE]
}

fwrite(ph_dt, out_ph, sep="\t", quote=FALSE, na="NA")

meta <- list(
  genotype_tsv = out_geno,
  marker_map_tsv = out_map,
  phenotype_tsv = out_ph,
  source_cross2_path = cross2_path,
  out_dir = base_out
)
writeLines(toJSON(meta, auto_unbox=TRUE, pretty=TRUE), con=file.path(out_dir, "artifacts.json"))
