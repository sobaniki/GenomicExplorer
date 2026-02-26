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

geno_tsv <- p$genotype_tsv
map_tsv  <- if (!is.null(p$marker_map_tsv)) p$marker_map_tsv else NULL
out_prefix <- if (!is.null(p$out_prefix)) p$out_prefix else file.path(out_dir, "plink", "data")

dir.create(dirname(out_prefix), recursive=TRUE, showWarnings=FALSE)

cat("[prep_tsv_to_plink] geno_tsv=", geno_tsv, "\n")
cat("[prep_tsv_to_plink] map_tsv=", map_tsv, "\n")
cat("[prep_tsv_to_plink] out_prefix=", out_prefix, "\n")

gt <- fread(geno_tsv, sep="\t", header=TRUE, data.table=FALSE)
#stopifnot("id" %in% colnames(gt))
ids <- as.character(gt[, 1])

markers <- colnames(gt)[2:ncol(gt)]
if (length(markers) == 0) stop("No marker columns in genotype_tsv")

X <- as.matrix(gt[, markers, drop=FALSE])
mode_num <- function(z){
  z <- suppressWarnings(as.numeric(z))
  z[is.nan(z)] <- NA
  z
}
X <- apply(X, 2, mode_num)

# map（あれば利用）
chr <- rep(1, length(markers))
pos <- seq_along(markers)
if (!is.null(map_tsv) && file.exists(map_tsv)) {
  mp <- fread(map_tsv, sep="\t", header=TRUE, data.table=FALSE)
  #stopifnot(all(c("marker","chr","pos") %in% colnames(mp)))
  #mp$marker <- as.character(mp$marker)
  #m2 <- match(markers, mp$marker)
  m2 <- match(markers, mp[, 1])
  ok <- !is.na(m2)
  #chr[ok] <- as.character(mp$chr[m2[ok]])
  #pos[ok] <- as.numeric(mp$pos[m2[ok]])
  chr[ok] <- as.character(mp[m2[ok], 2])
  pos[ok] <- as.numeric(mp[m2[ok], 3])
}

# bed.matrix を作る（gaston）
# 注：as.bed.matrix の引数は環境で差がある可能性があるので、必要ならここだけ修正してください。
bm <- as.bed.matrix(X)
bm@snps$id <- markers
bm@snps$chr <- chr
bm@snps$pos <- pos
bm@ped$id <- ids

# 必ずRDSで保存（最悪これが成果物）
rds_path <- file.path(out_dir, "bedmatrix.rds")
saveRDS(bm, rds_path)
cat("[prep_tsv_to_plink] wrote:", rds_path, "\n")

# PLINK出力（可能なら）
plink_ok <- FALSE
try({
  if ("write.bed.matrix" %in% getNamespaceExports("gaston")) {
    gaston::write.bed.matrix(bm, out_prefix)
    plink_ok <- TRUE
  }
}, silent=TRUE)

# もし上がダメなら、あなたがここを環境に合わせて1行直す想定
if (!plink_ok) {
  cat("[prep_tsv_to_plink] WARN: could not write PLINK files (write.bed.matrix not found / failed).\n")
  cat("[prep_tsv_to_plink]       bedmatrix.rds is available; you can adjust the PLINK-writer call here.\n")
} else {
  cat("[prep_tsv_to_plink] wrote PLINK prefix:", out_prefix, "\n")
}

# 成果物メタ（GUI側で拾いやすく）
meta <- list(
  bedmatrix_rds = rds_path,
  plink_prefix  = out_prefix,
  plink_written = plink_ok
)
writeLines(toJSON(meta, auto_unbox=TRUE, pretty=TRUE),
           con=file.path(out_dir, "artifacts.json"))

