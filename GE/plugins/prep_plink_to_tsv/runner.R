#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(jsonlite)
  library(data.table)
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

plink_prefix <- p$plink_prefix
if (is.null(plink_prefix) || !nzchar(plink_prefix)) stop("plink_prefix is required")
plink_prefix <- as.character(plink_prefix)

# Allow passing .bed path too
if (grepl("\\.bed$", plink_prefix, ignore.case=TRUE)) {
  plink_prefix <- sub("\\.bed$", "", plink_prefix, ignore.case=TRUE)
}

bed <- paste0(plink_prefix, ".bed")
bim <- paste0(plink_prefix, ".bim")
fam <- paste0(plink_prefix, ".fam")
if (!file.exists(bed) || !file.exists(bim) || !file.exists(fam)) {
  stop("PLINK files not found (.bed/.bim/.fam): ", plink_prefix)
}

out_override <- if (!is.null(p$out_dir_override)) as.character(p$out_dir_override) else ""
base_out <- if (nzchar(out_override)) out_override else file.path(out_dir, "tsv_from_plink")
dir.create(base_out, recursive=TRUE, showWarnings=FALSE)

cat("[prep_plink_to_tsv] plink_prefix=", plink_prefix, "\n")
cat("[prep_plink_to_tsv] base_out=", base_out, "\n")

# Read as bed.matrix
bm <- read.bed.matrix(bed)

# Save bed.matrix for interoperability
rds_path <- file.path(base_out, "bedmatrix.rds")
saveRDS(bm, rds_path)

ids <- as.character(bm@ped$id)
markers_raw <- as.character(bm@snps$id)
markers <- make.unique(markers_raw)
# Genotype matrix
X <- NULL
X <- tryCatch({
  as.matrix(bm)
}, error=function(e) {
  NULL
})
if (is.null(X)) {
  # Fallback: try to use internal conversion if available
  if ("as.matrix.bed.matrix" %in% getNamespaceExports("gaston")) {
    X <- gaston::as.matrix.bed.matrix(bm)
  } else {
    stop("Failed to convert bed.matrix to numeric matrix (as.matrix failed).")
  }
}

# Ensure numeric and use NA for missing
X <- suppressWarnings(matrix(as.numeric(X), nrow=nrow(X), ncol=ncol(X), dimnames=dimnames(X)))

# Write marker map (SPEC: 1st=marker, 2nd=chr, 3rd=pos)
marker_map_tsv <- file.path(base_out, "marker_map.tsv")
chr <- as.character(bm@snps$chr)
pos <- suppressWarnings(as.numeric(as.character(bm@snps$pos)))
chr[is.na(chr) | chr==''] <- '1'
if (all(is.na(pos))) pos <- seq_len(length(markers))
miss_pos <- is.na(pos)
if (any(miss_pos)) pos[miss_pos] <- seq_len(sum(miss_pos))
mm <- data.table(
  marker = as.character(markers),
  chr    = chr,
  pos    = pos
)
fwrite(mm, marker_map_tsv, sep="\t", na = "NA", quote = FALSE)

# Optional sidecar: raw marker IDs
variant_info_tsv <- file.path(base_out, "variant_info.tsv")
try({
  vi <- data.table(marker=as.character(markers), marker_raw=as.character(markers_raw))
  fwrite(vi, variant_info_tsv, sep="\t", na="NA", quote=FALSE)
}, silent=TRUE)

# Write pedigree (as available)
ped_tsv <- file.path(base_out, "pedigree.tsv")
try({
  ped <- as.data.table(bm@ped)
  fwrite(ped, ped_tsv, sep="\t", na = "NA", quote = F)
}, silent=TRUE)

# Write genotype wide TSV (id + markers)
genotype_tsv <- file.path(base_out, "genotype.tsv")
# Use data.table for efficiency
DT <- data.table(id = as.character(ids))
# Convert to data.table column-wise to avoid losing colnames
Xdt <- as.data.table(X)
setnames(Xdt, as.character(markers))
DT <- cbind(DT, Xdt)
fwrite(DT, genotype_tsv, sep="\t", na="NA", quote = F)

meta <- list(
  plink_prefix  = plink_prefix,
  bedmatrix_rds = rds_path,
  genotype_tsv  = genotype_tsv,
  marker_map_tsv = marker_map_tsv,
  pedigree_tsv = ped_tsv
)
writeLines(toJSON(meta, auto_unbox=TRUE, pretty=TRUE), con=file.path(out_dir, "artifacts.json"))

cat("[prep_plink_to_tsv] wrote:", genotype_tsv, "\n")
cat("[prep_plink_to_tsv] wrote:", marker_map_tsv, "\n")
cat("[prep_plink_to_tsv] wrote:", rds_path, "\n")
