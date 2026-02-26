#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(jsonlite)
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
mm_path <- if (!is.null(p$marker_map_tsv)) p$marker_map_tsv else NULL
out_vcf <- if (!is.null(p$out_vcf) && nzchar(p$out_vcf)) p$out_vcf else file.path(out_dir, "out.vcf")

cat("[prep_tsv_to_vcf] genotype_tsv=", gt_path, "\n")
cat("[prep_tsv_to_vcf] marker_map_tsv=", ifelse(is.null(mm_path), "(none)", mm_path), "\n")
cat("[prep_tsv_to_vcf] out_vcf=", out_vcf, "\n")

# Read genotype.tsv (wide)
df <- read.table(gt_path, header=TRUE, sep="\t", check.names=FALSE, stringsAsFactors=FALSE, comment.char="", quote="")
if (ncol(df) < 2) stop("genotype.tsv must have at least 2 columns (id + markers)")
ids <- df[[1]]
markers <- colnames(df)[-1]
G <- df[, -1, drop=FALSE]

# Read marker_map.tsv (optional)
mm <- NULL
if (!is.null(mm_path) && file.exists(mm_path)) {
  mm <- read.table(mm_path, header=TRUE, sep="\t", check.names=FALSE, stringsAsFactors=FALSE, comment.char="", quote="")
}

get_col <- function(x, keys) {
  for (k in keys) {
    if (!is.null(x) && k %in% names(x)) return(x[[k]])
  }
  NULL
}

# Spec: marker_map.tsv fixed as 1st=marker, 2nd=chr, 3rd=pos.
mm_marker <- NULL
mm_chr <- NULL
mm_pos <- NULL
if (!is.null(mm) && ncol(mm) >= 3) {
  mm_marker <- mm[[1]]
  mm_chr <- mm[[2]]
  mm_pos <- mm[[3]]
}
mm_ref <- get_col(mm, c("ref","REF","A1"))
mm_alt <- get_col(mm, c("alt","ALT","A2"))

# Build per-marker metadata aligned to genotype columns
chr <- rep("1", length(markers))
pos <- seq_len(length(markers))
ref <- rep("A", length(markers))
alt <- rep("C", length(markers))

if (!is.null(mm) && !is.null(mm_marker)) {
  idx <- match(markers, trimws(as.character(mm_marker)))
  ok <- !is.na(idx)
  if (!is.null(mm_chr)) chr[ok] <- as.character(mm_chr[idx[ok]])
  if (!is.null(mm_pos)) {
    pos_str <- gsub(",", "", trimws(as.character(mm_pos[idx[ok]])))
    pos_num <- suppressWarnings(as.numeric(pos_str))

    # VCF POS must be integer. If the input looks like genetic position (cM/Mb; small values and/or decimals),
    # convert by assuming 1 cM ≈ 1 Mbp => POS_bp ≈ POS_cM * 1e6.
    has <- !is.na(pos_num)
    if (any(has)) {
      frac <- abs(pos_num - round(pos_num)) > 1e-9
      mx <- suppressWarnings(max(pos_num[has], na.rm=TRUE))

      if (any(frac) || (is.finite(mx) && mx < 1e5)) {
        pos_int <- round(pos_num * 1e6)
      } else {
        pos_int <- round(pos_num)
      }
      pos[ok][has] <- as.integer(pos_int[has])
    }
  }
  if (!is.null(mm_ref)) ref[ok] <- as.character(mm_ref[idx[ok]])
  if (!is.null(mm_alt)) alt[ok] <- as.character(mm_alt[idx[ok]])
}

# Normalize missing/empty alleles
ref[is.na(ref) | ref==""] <- "A"
alt[is.na(alt) | alt==""] <- "C"
chr[is.na(chr) | chr==""] <- "1"
pos[is.na(pos)] <- seq_len(sum(is.na(pos)))

# Open output connection (supports .gz)
con <- NULL
if (grepl("\\.gz$", out_vcf, ignore.case=TRUE)) {
  con <- gzfile(out_vcf, "wt")
} else {
  con <- file(out_vcf, "wt")
}

on.exit({ try(close(con), silent=TRUE) }, add=TRUE)

# Write header
hdr <- c(
  "##fileformat=VCFv4.2",
  "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
  paste0("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t", paste(ids, collapse="\t"))
)
writeLines(hdr, con=con)

# Helper: convert genotype value to GT string
as_gt <- function(v) {
  if (is.na(v) || v=="" || v=="." || v=="NA" || v=="NaN" || v=="-9") return("./.")
  # already GT-like
  if (grepl("[\\/\\|]", v)) {
    if (v=="." || v=="./." ) return("./.")
    return(v)
  }
  # numeric 0/1/2
  vv <- suppressWarnings(as.integer(v))
  if (is.na(vv)) return("./.")
  if (vv==0) return("0/0")
  if (vv==1) return("0/1")
  if (vv==2) return("1/1")
  return("./.")
}

# Write records (stream)
for (j in seq_along(markers)) {
  m <- markers[j]
  gcol <- G[[j]]
  gts <- vapply(as.character(gcol), as_gt, character(1))
  line <- paste(chr[j], pos[j], m, ref[j], alt[j], ".", "PASS", ".", "GT", paste(gts, collapse="\t"), sep="\t")
  writeLines(line, con=con)
}

meta <- list(
  vcf_path = out_vcf,
  genotype_tsv = gt_path,
  marker_map_tsv = ifelse(is.null(mm_path), "", mm_path)
)
writeLines(toJSON(meta, auto_unbox=TRUE, pretty=TRUE),
           con=file.path(out_dir, "artifacts.json"))
