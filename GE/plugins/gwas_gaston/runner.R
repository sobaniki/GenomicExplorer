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
plot_dir <- file.path(out_dir, "plots")
dir.create(plot_dir, recursive=TRUE, showWarnings=FALSE)

log_file <- file.path(out_dir, "run.log")
sink(log_file, split=TRUE)

cat("[GWAS] start\n")
cat("[GWAS] params_path=", params_path, "\n")
cat("[GWAS] out_dir=", out_dir, "\n")

p <- fromJSON(params_path)

# -----------------------
# Params (common)
# -----------------------
method <- p$method
phenotype_tsv <- p$phenotype_tsv
trait <- p$trait
maf <- if (!is.null(p$maf)) as.numeric(p$maf) else 0.05
missing_max <- if (!is.null(p$missing_max)) as.numeric(p$missing_max) else 0.1
use_lmm <- if (!is.null(p$use_lmm)) as.logical(p$use_lmm) else TRUE
cov_tsv <- if (!is.null(p$covariates_tsv)) p$covariates_tsv else NULL
if (!is.null(cov_tsv) && nchar(cov_tsv) == 0) cov_tsv <- NULL

stopifnot(!is.null(phenotype_tsv), file.exists(phenotype_tsv))
stopifnot(!is.null(trait), nchar(trait) > 0)

cat("[GWAS] phenotype_tsv=", phenotype_tsv, "\n")
cat("[GWAS] trait=", trait, "\n")
cat("[GWAS] maf=", maf, " missing_max=", missing_max, " use_lmm=", use_lmm, "\n")
cat("[GWAS] covariates_tsv=", cov_tsv, "\n")

# -----------------------
# Gaston-specific options
# -----------------------
if (method == "LM/LMM") {
  test <- if (!is.null(p$test)) as.character(p$test) else "wald"
  pc_n <- if (!is.null(p$p)) as.integer(p$p) else 0
  
  if (!test %in% c("score", "wald", "lrt")) {
    cat("[gwas_gaston] WARN: unknown test=", test, " -> fallback to wald\n")
    test <- "wald"
  }
  if (is.na(pc_n) || pc_n < 0) pc_n <- 0
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
# Input switch (plink/vcf/tsv)
# -----------------------
#plink_prefix <- if (!is.null(p$plink_prefix)) p$plink_prefix else NULL
#vcf_path <- if (!is.null(p$vcf_path)) p$vcf_path else NULL
genotype_tsv <- if (!is.null(p$genotype_tsv)) p$genotype_tsv else NULL
marker_map_tsv <- if (!is.null(p$marker_map_tsv)) p$marker_map_tsv else NULL

#cat("[GWAS] plink_prefix=", plink_prefix, "\n")
#cat("[GWAS] vcf_path=", vcf_path, "\n")
cat("[GWAS] genotype_file=", genotype_tsv, "\n")
cat("[GWAS] marker_map_tsv=", marker_map_tsv, "\n")

# has_plink <- !is.null(plink_prefix) && file.exists(paste0(plink_prefix, ".bed"))
# has_vcf   <- !is.null(vcf_path) && file.exists(vcf_path)
# has_tsv   <- !is.null(genotype_tsv) && file.exists(genotype_tsv)
# 
# if (!has_plink && !has_vcf && !has_tsv) {
#   stop("Must provide one of: plink_prefix(.bed), vcf_path, genotype_tsv")
# }
# if (sum(c(has_plink, has_vcf, has_tsv)) > 1) {
#   cat("[GWAS] WARN: multiple genotype inputs provided; priority plink > vcf > tsv\n")
# }

# -----------------------
# Load genotype as bed.matrix
# -----------------------
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
cat("[GWAS] genotype loaded. n=", nrow(bm@ped), " m=", nrow(bm@snps), "\n")

# -----------------------
# Load phenotype and align IDs
# -----------------------
ph <- fread(phenotype_tsv, sep="\t", header=TRUE, data.table=FALSE)
stopifnot(trait %in% colnames(ph))

gid <- as.character(bm@ped$id)
ord <- match(gid, ph[, 1])

keep <- !is.na(ord)
if (sum(keep) < 5) stop("Too few overlapping samples between genotype and phenotype (by id)")

# subset genotype & phenotype to overlap
bm <- bm[keep, ]
ph2 <- ph[ord[keep], , drop=FALSE]

y <- suppressWarnings(as.numeric(ph2[[trait]]))
keep_y <- !is.na(y)
if (sum(keep_y) < 5) stop("Too few non-missing phenotype values")

bm <- bm[keep_y, ]
ph2 <- ph2[keep_y, , drop=FALSE]
y <- y[keep_y]

cat("[GWAS] aligned samples n=", length(y), "\n")

# -----------------------
# Optional covariates
# -----------------------
Xcov <- NULL
if (!is.null(cov_tsv) && file.exists(cov_tsv)) {
  cv <- fread(cov_tsv, sep="\t", header=TRUE, data.table=FALSE)
  
  ord2 <- match(as.character(bm@ped$id), cv[, 1])
  keep2 <- !is.na(ord2)
  if (sum(keep2) != nrow(bm@ped)) {
    cat("[GWAS] WARN: covariates missing for some samples; dropping those samples\n")
    bm <- bm[keep2, ]
    ph2 <- ph2[keep2, , drop=FALSE]
    y <- y[keep2]
    ord2 <- ord2[keep2]
  }
  
  cv2 <- cv[ord2, , drop=FALSE]
  cov_cols <- colnames(cv2)[2:ncol(cv2)]
  if (length(cov_cols) > 0) {
    Xcov <- as.matrix(cv2[, cov_cols, drop=FALSE])
    Xcov <- apply(Xcov, 2, function(z) suppressWarnings(as.numeric(z)))
    Xcov[is.nan(Xcov)] <- NA
    if (any(!is.finite(Xcov))) {
      cat("[GWAS] WARN: NA/NaN in covariates; dropping rows with any NA\n")
      okc <- apply(Xcov, 1, function(r) all(is.finite(r)))
      bm <- bm[okc, ]
      ph2 <- ph2[okc, , drop=FALSE]
      y <- y[okc]
      Xcov <- Xcov[okc, , drop=FALSE]
    }
    cat("[GWAS] covariates p=", ncol(Xcov), "\n")
  }
}

X <- NULL
if (is.null(Xcov)) {
  X <- matrix(1, nrow(bm))
} else {
  X <- cbind(rep(1, nrow(bm)), Xcov)
}

# -----------------------
# QC filters using gaston fields
# -----------------------
# callrate & maf are stored at bm@snps after loading
miss_rate <- 1 - bm@snps$callrate
maf_vec <- bm@snps$maf

keep_snp <- rep(TRUE, nrow(bm@snps))
if (!is.null(missing_max) && is.finite(missing_max)) keep_snp <- keep_snp & (miss_rate <= missing_max)
if (!is.null(maf) && is.finite(maf)) keep_snp <- keep_snp & (maf_vec >= maf)

cat("[GWAS] QC keep SNPs=", sum(keep_snp), "/", length(keep_snp), "\n")
bm2 <- select.snps(bm, which(keep_snp))

# -----------------------
# GWAS
# -----------------------
cat("[GWAS] association...\n")

if (method != "FarmCPU") {
  if (!is.null(K_user)) {
    Kmat <- K_user
    eigenK <- eigen(Kmat, symmetric=TRUE)
  } else {
    cat("[GWAS] computing GRM from genotype (K = NULL)\n")
    Kmat <- GRM(bm2, autosome.only = F)
    eigenK <- eigen(Kmat, symmetric=TRUE)
  }
}

if (method == "LM/LMM") {
  binary <- length(unique(y)) == 2
  if (binary == F) {
    response <- "quantitative"
  } else {
    response <- "binary"
  }
  
  # Run association
  if (use_lmm) {
    if (binary == T) {
      if (test != "score") {
        cat("[GWAS] WARN: test=", test, " requested with LMM for a binary trait; forcing score\n")
        test <- "score"
      }
      # score test uses K (or list of K). If user didn't provide, we computed GRM above.
      ans <- association.test(bm2, 
                              Y = y, 
                              X = X, 
                              eigenK = eigenK,
                              K = Kmat, 
                              method = "lmm",
                              response=response, 
                              test = test,
                              p = pc_n)
    } else {
      if (test == "score") {
        cat("[GWAS] WARN: test=", test, " requested with LMM for a quantitative trait; forcing wald\n")
        test <- "wald"
      }
      ans <- association.test(bm2, 
                              Y = y, 
                              X = X, 
                              eigenK = eigenK, 
                              method = "lmm",
                              response = response, 
                              test = test,
                              p = pc_n)
    }
  } else {
    if (test != "wald") {
      cat("[GWAS] WARN: test=", test, " requested with LM; forcing wald\n")
      test <- "wald"
    }
    ans <- association.test(bm2, 
                            Y = y, 
                            X = X, 
                            method = "lm",
                            response = response, 
                            test = test,
                            p = pc_n,
                            eigenK = eigenK)
  }
  saveRDS(ans, file.path(out_dir, "gaston_result.rds"))
  
  # -----------------------
  # Write results
  # -----------------------
  if (test == "wald") {
    res <- data.frame(
      marker = bm2@snps$id,
      chr    = bm2@snps$chr,
      pos    = bm2@snps$pos,
      beta   = ans$beta,
      sd     = ans$sd,
      pvalue = ans$p
    )
  } else if (test == "lrt") {
      res <- data.frame(
        marker = bm2@snps$id,
        chr    = bm2@snps$chr,
        pos    = bm2@snps$pos,
        LRT = ans$LRT,
        pvalue = ans$p
      )
  } else {
    res <- data.frame(
      marker = bm2@snps$id,
      chr    = bm2@snps$chr,
      pos    = bm2@snps$pos,
      score = ans$score,
      pvalue = ans$p
    )
  }
} else if (method == "MLMM") {
  max_steps <- if (!is.null(p$mlmm_max_steps)) as.integer(p$mlmm_max_steps) else 10
  
  if (ncol(X) == 1) {
    mlmm_out <- mlmm::mlmm(Y = y,
                           X = as.matrix(bm2),
                           K = Kmat,
                           maxsteps = max_steps,
                           nbchunks = 2)
  } else {
    mlmm_out <- mlmm::mlmm_cof(Y = y,
                               X = as.matrix(bm2),
                               K = Kmat,
                               maxsteps = max_steps,
                               cof = X,
                               nbchunks = 2)
  }
  saveRDS(mlmm_out, file.path(out_dir, "mlmm_result.rds"))
  
  res <- data.frame(
    marker = as.character(bm2@snps$id),
    chr    = bm2@snps$chr,
    pos    = bm2@snps$pos,
    pvalue = mlmm_out$opt_mbonf$out$pval,
    stringsAsFactors = FALSE
  )
} else {
  BK <- bigmemory::as.big.matrix(as.matrix(bm2), type = "double")
  map_df <- data.frame(snp = bm2@snps$id,
                       chr = as.numeric(as.factor(bm2@snps$chr)),
                       pos = bm2@snps$pos)
  
  gwas_df <- FarmCPUpp::farmcpu(Y = data.frame(id = bm2@ped$id,
                                               pheno = y),
                                GD = BK,
                                GM = map_df)
  
  saveRDS(gwas_df, file.path(out_dir, "farmcpu_result.rds"))
  
  res <- data.frame(
    marker = as.character(bm2@snps$id),
    chr    = bm2@snps$chr,
    pos    = bm2@snps$pos,
    pvalue = gwas_df[[1]]$GWAS$p.value,
    stringsAsFactors = FALSE
  )
}
res_path <- file.path(out_dir, "results.tsv")
fwrite(res, res_path, sep="\t")

# Add BH-FDR q-values (best-effort)
if (!("qvalue" %in% names(res))) {
  res$qvalue <- p.adjust(res$pvalue, method="BH")
}

fwrite(res, res_path, sep="\t")
cat("[GWAS] wrote:", res_path, "\n")

# Manhattan plot (simple)
png_path <- file.path(plot_dir, "manhattan.png")
png(png_path, width=1200, height=450)
# create x coordinate by chr blocks
dfp <- res
dfp$chr <- as.character(dfp$chr)
dfp$pos <- suppressWarnings(as.numeric(dfp$pos))
dfp$mlogp <- -log10(dfp$pvalue)

# order chromosomes numerically when possible
chr_key <- function(x) {
  suppressWarnings(v <- as.numeric(x))
  ifelse(is.na(v), 1e9, v)
}
chrs <- unique(dfp$chr[order(chr_key(dfp$chr))])

offset <- 0
dfp$x <- NA_real_
ticks <- c()
ticklabs <- c()
for (c in chrs) {
  idx <- which(dfp$chr == c)
  sub <- dfp[idx, ]
  sub <- sub[order(sub$pos), ]
  if (all(is.na(sub$pos))) {
    dfp$x[idx] <- seq_along(idx) + offset
  } else {
    dfp$x[idx][order(sub$pos)] <- sub$pos[order(sub$pos)] + offset
  }
  ticks <- c(ticks, mean(dfp$x[idx], na.rm=TRUE))
  ticklabs <- c(ticklabs, c)
  offset <- max(dfp$x[idx], na.rm=TRUE) + 1
}

plot(dfp$x, dfp$mlogp, pch=20, cex=0.6,
     xlab="Chromosome", ylab="-log10(p)",
     main=paste0("GWAS (trait=", trait, ")"))
axis(1, at=ticks, labels=ticklabs)


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

cat("[GWAS] wrote:", png_path, "\n")
cat("[GWAS] done\n")
sink()
