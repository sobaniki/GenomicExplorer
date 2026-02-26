#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

suppressPackageStartupMessages({
  library(jsonlite)
  library(data.table)
})

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 4 || args[1] != "--params" || args[3] != "--out") {
  cat("Usage: runner.R --params params.json --out out_dir\n")
  quit(status=1)
}
params_path <- args[2]
out_dir <- args[4]
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
plot_dir <- file.path(out_dir, "plots")
dir.create(plot_dir, recursive=TRUE, showWarnings=FALSE)

p <- fromJSON(params_path)

# -----------------------
# Helper: standardize covariates
# -----------------------
read_covariates <- function(path, id_col="id") {
  if (is.null(path) || !file.exists(path)) return(NULL)
  cv <- fread(path, sep="\t", header=TRUE)
  if (!(id_col %in% names(cv))) {
    stop("covariates_tsv must contain an 'id' column")
  }
  cv[[id_col]] <- as.character(cv[[id_col]])
  m <- as.matrix(cv[, setdiff(names(cv), id_col), with=FALSE])
  rownames(m) <- cv[[id_col]]
  storage.mode(m) <- "numeric"
  return(m)
}

# -----------------------
# Load phenotype
# -----------------------
pheno_tsv <- p$phenotype_tsv
trait <- p$trait
if (is.null(pheno_tsv) || !file.exists(pheno_tsv)) stop("phenotype_tsv not found")
ph <- fread(pheno_tsv, sep="\t", header=TRUE)
if (!("id" %in% names(ph))) stop("phenotype_tsv must have 'id' column")
if (!(trait %in% names(ph))) stop(paste0("trait not found: ", trait))
ph$id <- as.character(ph$id)
y <- ph[[trait]]
names(y) <- ph$id

# covariates (optional)
cov_path <- p$covariates_tsv
COV <- read_covariates(cov_path)

# -----------------------
# Load genotype from PLINK (.bed)
# -----------------------
# plink_prefix <- p$plink_prefix
# if (is.null(plink_prefix)) stop("plink_prefix missing")
# bed_path <- paste0(plink_prefix, ".bed")
# if (!file.exists(bed_path)) stop(paste0("PLINK .bed not found: ", bed_path))
# 
# suppressPackageStartupMessages({
#   library(gaston)
# })
# bm <- read.bed.matrix(plink_prefix)
genotype_tsv <- if (!is.null(p$genotype_tsv)) p$genotype_tsv else NULL

bm <- NULL
if (grepl(".bed$", genotype_tsv)) {
  cat("[GWAS] loading PLINK via read.bed.matrix...\n")
  bm <- read.bed.matrix(genotype_tsv)
} else if (grepl(".vcf(|.gz)$", genotype_tsv, perl = T)) {
  cat("[GWAS] loading VCF via read.vcf...\n")
  bm <- read.vcf(genotype_tsv, convert.chr = F)
} else {
  cat("[GWAS] loading TSV -> as.bed.matrix...\n")
  gt <- fread(genotype_tsv, sep="\t", header=TRUE, data.table=FALSE)
  ids <- as.character(gt[, 1])
  markers <- colnames(gt[, 2:ncol(gt)])
  if (length(markers) == 0) stop("No marker columns in genotype_tsv")
  
  X <- as.matrix(gt[, markers, drop=FALSE])
  X <- apply(X, 2, function(z) suppressWarnings(as.numeric(z)))
  X[is.nan(X)] <- NA
  
  chr <- rep(1, length(markers))
  pos <- seq_along(markers)
  
  if (!is.null(marker_map_tsv) && file.exists(marker_map_tsv)) {
    mp <- fread(marker_map_tsv, sep="\t", header=TRUE, data.table=FALSE)
    m2 <- match(markers, as.character(mp[ ,1]))
    ok <- !is.na(m2)
    chr[ok] <- as.character(mp[ok, 2])
    pos[ok] <- as.numeric(mp[ok, 3])
  } else {
    cat("[GWAS] marker_map_tsv not provided; using dummy chr/pos\n")
  }
  
  # NOTE: as.bed.matrix signatures may vary by version.
  # If your environment needs different args, adjust here.
  bm <- as.bed.matrix(X)
  bm@snps$id <- markers
  bm@snps$chr <- chr
  bm@snps$pos <- pos
  bm@ped$id <- ids
}

# match individuals by id
iid <- as.character(bm@ped$id)
common <- intersect(names(y), iid)
if (length(common) < 5) stop("Too few common individuals between phenotype and genotype")

y2 <- y[common]
bm2 <- bm[which(iid %in% common), ]

# reorder bm2 rows to match y2
ord <- match(names(y2), as.character(bm2@ped$id))
bm2 <- bm2[ord, ]

# -----------------------
# Basic SNP filtering (MAF / missing)
# -----------------------
maf <- if (!is.null(p$maf)) as.numeric(p$maf) else 0.01
missing_max <- if (!is.null(p$missing_max)) as.numeric(p$missing_max) else 0.20

g <- as.matrix(bm2)  # n x m, values 0/1/2 with NA
# missing per SNP
miss <- colMeans(is.na(g))
keep_miss <- miss <= missing_max
# MAF per SNP (ignore missing)
af <- colMeans(g, na.rm=TRUE) / 2
maf_v <- pmin(af, 1-af)
keep_maf <- maf_v >= maf
keep <- keep_miss & keep_maf

if (sum(keep) < 10) {
  warning("Too few markers after filtering; creating empty results.tsv")
  fwrite(data.frame(), file.path(out_dir, "results.tsv"), sep="\t")
  quit(status=0)
}

g <- g[, keep, drop=FALSE]
snps <- bm2@snps[keep, ]

# -----------------------
# Impute missing (MLMM does not accept NAs)
# -----------------------
if (any(is.na(g))) {
  cat("[gwas_mlmm] imputing missing genotypes by SNP mean\n")
  mu <- colMeans(g, na.rm=TRUE)
  for (j in seq_len(ncol(g))) {
    idx_na <- which(is.na(g[, j]))
    if (length(idx_na) > 0) g[idx_na, j] <- mu[j]
  }
}

K_paths <- if (!is.null(p$K)) p$K else NULL  # may be string or vector
K_user <- NULL
if (!is.null(K_paths)) {
  # jsonlite may parse a single string or an array; normalize to character vector
  if (is.list(K_paths)) K_paths <- unlist(K_paths)
  K_paths <- as.character(K_paths)
  K_paths <- K_paths[nzchar(K_paths)]
  if (length(K_paths) > 0) {
    for (kp in K_paths) {
      if (!file.exists(kp)) stop(paste0("K file not found: ", kp))
    }
    if (length(K_paths) == 1) {
      cat("[gwas_gaston] loading K from:", K_paths[1], "\n")
      if (grepl("*.rds$", K_paths[1])) {
        K_user <- readRDS(K_paths[1])
      } else {
        K_user <- as.matrix(data.frame(data.table::fread(K_paths[1], header = T), row.names = 1))
      }
    } else {
      cat("[gwas_gaston] loading multiple K matrices (list) n=", length(K_paths), "\n")
      K_user <- lapply(K_paths, readRDS)
    }
  }
}

# -----------------------
# Kinship matrix (GRM)
# -----------------------
if (!is.null(K_user)) {
  Kmat <- K_user
} else {
  cat("[GWAS] computing GRM from genotype (K = NULL)\n")
  Kmat <- GRM(bm2, autosome.only = F)
}
Kmat <- Kmat[common, common]
rownames(Kmat) <- colnames(Kmat) <- common

# COV to matrix aligned with individuals
cofs <- NULL
if (!is.null(COV)) {
  cofs <- COV[common, , drop=FALSE]
}

# -----------------------
# Run MLMM
# -----------------------
max_steps <- if (!is.null(p$mlmm_max_steps)) as.integer(p$mlmm_max_steps) else 10

# Ensure names/rownames
names(y2) <- common
rownames(g) <- common
colnames(g) <- as.character(snps$id)

if (!(is.null(cofs))) {
  mlmm_out <- mlmm::mlmm(Y = y2,
                         X = g,
                         K = Kmat,
                         maxsteps = max_steps,
                         nbchunks = 2)
} else {
  mlmm_out <- mlmm::mlmm_cof(Y = y2,
                             X = g,
                             K = Kmat,
                             maxsteps = max_steps,
                             cof = cofs,
                             nbchunks = 2)
}
saveRDS(mlmm_out, file.path(out_dir, "mlmm_result.rds"))

res <- data.frame(
  marker = as.character(snps$id),
  chr    = snps$chr,
  pos    = snps$pos,
  pvalue = mlmm_out$opt_mbonf$out$pval,
  #thresh = mlmm_out$bonf_thresh,
  #cof = mlmm_out$opt_mbonf$cof,
  stringsAsFactors = FALSE
)
res_path <- file.path(out_dir, "results.tsv")
fwrite(res, res_path, sep="\t")
cat("[gwas_mlmm] wrote:", res_path, "\n")

# -----------------------
# Manhattan plot (static PNG; interactive HTML is generated in GUI if Plotly is available)
# -----------------------
png_path <- file.path(plot_dir, "manhattan.png")
png(png_path, width=1200, height=450)
dfp <- res
dfp$chr <- as.character(dfp$chr)
dfp$pos <- suppressWarnings(as.numeric(dfp$pos))
dfp$mlogp <- -log10(dfp$pvalue)
dfp <- dfp[order(as.numeric(gsub("[^0-9]", "", dfp$chr)), dfp$pos), ]
dfp$x <- seq_len(nrow(dfp))
plot(dfp$x, dfp$mlogp, pch=20, cex=0.4, xlab="Markers", ylab="-log10(p)", main=paste0("MLMM - ", trait))


# QQ plot (simple)
qq_png <- file.path(plot_dir, "qq.png")
pp <- suppressWarnings(as.numeric(res$pvalue))
pp <- pp[is.finite(pp) & pp > 0 & pp <= 1]
if (length(pp) >= 10) {
  pp <- sort(pp)
  nqq <- length(pp)
  exp <- -log10((1:nqq) / (nqq + 1))
  obs <- -log10(pp)
  # flip to conventional axis direction (0 -> high)
  exp <- rev(exp)
  obs <- rev(obs)
  # lambda GC (df=1)
  lam <- NA_real_
  chisq <- qchisq(1 - pp, df=1)
  if (length(chisq) > 0) {
    lam <- median(chisq, na.rm=TRUE) / qchisq(0.5, df=1)
  }

  png(qq_png, width=700, height=700)
  plot(exp, obs, pch=20, cex=0.6,
       xlab="Expected -log10(p)", ylab="Observed -log10(p)",
       main=paste0("QQ plot (trait=", trait, ")"))
  abline(0, 1, col="red", lwd=2)
  if (is.finite(lam)) {
    mtext(paste0("lambdaGC=", sprintf("%.3f", lam)), side=3, adj=1, line=-1)
  }
  dev.off()

  qq_html <- file.path(plot_dir, "qq.html")
  if (!file.exists(qq_html)) {
    writeLines("<html><body><p>Interactive QQ plot will be generated by GUI (Plotly).</p></body></html>", qq_html)
  }
  cat("[GWAS] wrote:", qq_png, "\n")
}

dev.off()

# placeholder for HTML path so GUI can detect it (created later)
html_path <- file.path(plot_dir, "manhattan.html")
if (!file.exists(html_path)) {
  writeLines("<html><body><p>Interactive Manhattan plot will be generated by GUI (Plotly).</p></body></html>", html_path)
}

cat("[gwas_mlmm] wrote:", png_path, "\n")
