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
out_prefix <- if (!is.null(p$out_prefix)) p$out_prefix else file.path(out_dir, "plink", "data")
dir.create(dirname(out_prefix), recursive=TRUE, showWarnings=FALSE)

cat("[prep_vcf_to_plink] vcf_path=", vcf_path, "\n")
cat("[prep_vcf_to_plink] out_prefix=", out_prefix, "\n")

# VCF¨bed.matrixigastonj
bm <- read.vcf(vcf_path, convert.chr = F)

# •K‚¸RDS
rds_path <- file.path(out_dir, "bedmatrix.rds")
saveRDS(bm, rds_path)
cat("[prep_vcf_to_plink] wrote:", rds_path, "\n")

plink_ok <- FALSE
try({
  if ("write.bed.matrix" %in% getNamespaceExports("gaston")) {
    gaston::write.bed.matrix(bm, out_prefix)
    plink_ok <- TRUE
  }
}, silent=TRUE)

if (!plink_ok) {
  cat("[prep_vcf_to_plink] WARN: could not write PLINK files. bedmatrix.rds is available.\n")
} else {
  cat("[prep_vcf_to_plink] wrote PLINK prefix:", out_prefix, "\n")
}

meta <- list(
  bedmatrix_rds = rds_path,
  plink_prefix  = out_prefix,
  plink_written = plink_ok
)
writeLines(toJSON(meta, auto_unbox=TRUE, pretty=TRUE),
           con=file.path(out_dir, "artifacts.json"))
