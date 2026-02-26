#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(jsonlite)
  library(gaston)
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
vcf_path <- p$vcf_path
if (is.null(vcf_path) || !file.exists(vcf_path)) stop("vcf_path not found")

tsv_dir <- NULL
if (!is.null(p$out_dir_override) && nzchar(p$out_dir_override)) {
  tsv_dir <- p$out_dir_override
} else {
  tsv_dir <- file.path(out_dir, "tsv_from_vcf")
}
dir.create(tsv_dir, recursive=TRUE, showWarnings=FALSE)

cat("[prep_vcf_to_tsv] vcf_path=", vcf_path, "\n")
cat("[prep_vcf_to_tsv] tsv_dir=", tsv_dir, "\n")

bm <- read.vcf(vcf_path, convert.chr = FALSE)

# Save bedmatrix
rds_path <- file.path(out_dir, "bedmatrix.rds")
saveRDS(bm, rds_path)
cat("[prep_vcf_to_tsv] wrote:", rds_path, "\n")

# Extract matrix
G <- tryCatch(as.matrix(bm), error=function(e) NULL)
if (is.null(G)) stop("Failed to convert bed.matrix to matrix")

# IDs
ids <- NULL
ids <- tryCatch(bm@ped$id, error=function(e) NULL)
if (is.null(ids) || length(ids) != nrow(G)) {
  ids <- rownames(G)
}
if (is.null(ids) || length(ids) != nrow(G)) {
  ids <- paste0("ind", seq_len(nrow(G)))
}

# Marker names
mnames <- NULL
mnames <- tryCatch(bm@snps$id, error=function(e) NULL)
if (is.null(mnames) || length(mnames) != ncol(G)) {
  mnames <- colnames(G)
}
if (is.null(mnames) || length(mnames) != ncol(G)) {
  mnames <- paste0("m", seq_len(ncol(G)))
}

# Ensure marker names are unique and consistent between genotype.tsv and marker_map.tsv
mnames_raw <- as.character(mnames)
mnames <- make.unique(mnames_raw)

# Helper: parse chr:pos from marker ids when available (e.g. VCF without explicit map)
parse_chrpos <- function(x) {
  x2 <- as.character(x)
  m <- regexec('^([^:]+):([0-9]+)$', x2)
  r <- regmatches(x2, m)
  chr <- rep(NA_character_, length(x2))
  pos <- rep(NA_real_, length(x2))
  ok <- lengths(r) == 3
  if (any(ok)) {
    chr[ok] <- vapply(r[ok], function(z) z[2], character(1))
    pos[ok] <- suppressWarnings(as.numeric(vapply(r[ok], function(z) z[3], character(1))))
  }
  list(chr=chr, pos=pos)
}

normalize_chr <- function(x) {
  # Spec: do not auto-normalize/convert chromosome labels.
  as.character(x)
}

# Write genotype.tsv (wide)
geno_path <- file.path(tsv_dir, "genotype.tsv")
df <- data.frame(id = ids, stringsAsFactors = FALSE)
# Convert NA to empty; keep numeric
G2 <- G
# ensure numeric (may be integer)
# leave as is; write.table will write NA as blank by default when na=""
df <- cbind(df, as.data.frame(G2, check.names = FALSE))
colnames(df)[-1] <- mnames
write.table(df, file=geno_path, sep="\t", quote=FALSE, row.names=FALSE, na="")
cat("[prep_vcf_to_tsv] wrote:", geno_path, "\n")

# marker_map.tsv (SPEC: 1st=marker, 2nd=chr, 3rd=pos)
mm_path <- file.path(tsv_dir, "marker_map.tsv")
mm <- NULL
mm <- tryCatch(bm@snps, error=function(e) NULL)
if (!is.null(mm) && all(c('id','chr','pos') %in% names(mm)) && length(mm$id) == length(mnames)) {
  chr <- normalize_chr(as.character(mm$chr))
  pos <- suppressWarnings(as.numeric(as.character(mm$pos)))
  # If pos is non-numeric, fallback to sequential indices
  if (all(is.na(pos))) pos <- seq_len(length(mnames))
  miss_pos <- is.na(pos)
  if (any(miss_pos)) pos[miss_pos] <- seq_len(sum(miss_pos))
  chr[is.na(chr) | chr==''] <- '1'
  out_mm <- data.frame(
    marker = as.character(mnames),
    chr = chr,
    pos = pos,
    stringsAsFactors = FALSE
  )

  # Keep raw IDs / alleles in a sidecar file (optional; not used by converters)
  vi_path <- file.path(tsv_dir, "variant_info.tsv")
  vi <- data.frame(marker=as.character(mnames), marker_raw=as.character(mm$id), stringsAsFactors=FALSE)
  if ('A1' %in% names(mm)) vi$ref <- as.character(mm$A1)
  if ('A2' %in% names(mm)) vi$alt <- as.character(mm$A2)
  write.table(vi, file=vi_path, sep="\t", quote=FALSE, row.names=FALSE)
} else {
  out_mm <- data.frame(
    marker = as.character(mnames),
    chr = rep('1', length(mnames)),
    pos = seq_len(length(mnames)),
    stringsAsFactors = FALSE
  )
  parsed <- parse_chrpos(as.character(mnames_raw))
  ok <- !is.na(parsed$chr) & !is.na(parsed$pos)
  if (any(ok)) {
    out_mm$chr[ok] <- normalize_chr(parsed$chr[ok])
    out_mm$pos[ok] <- parsed$pos[ok]
  }
}
write.table(out_mm, file=mm_path, sep="\t", quote=FALSE, row.names=FALSE)
cat("[prep_vcf_to_tsv] wrote:", mm_path, "\n")

# pedigree.tsv (if available)
ped_path <- file.path(tsv_dir, "pedigree.tsv")
ped <- NULL
ped <- tryCatch(bm@ped, error=function(e) NULL)
if (!is.null(ped)) {
  write.table(ped, file=ped_path, sep="\t", quote=FALSE, row.names=FALSE)
  cat("[prep_vcf_to_tsv] wrote:", ped_path, "\n")
}

meta <- list(
  bedmatrix_rds = rds_path,
  genotype_tsv  = geno_path,
  marker_map_tsv = mm_path,
  pedigree_tsv = if (file.exists(ped_path)) ped_path else "",
  tsv_dir = tsv_dir
)
writeLines(toJSON(meta, auto_unbox=TRUE, pretty=TRUE),
           con=file.path(out_dir, "artifacts.json"))
