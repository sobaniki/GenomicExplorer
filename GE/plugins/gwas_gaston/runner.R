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
# trait(s) may be empty => auto-detect numeric traits
traits_raw <- NULL
if (!is.null(p$traits) && nzchar(p$traits)) traits_raw <- p$traits
if (is.null(traits_raw) && !is.null(p$trait) && nzchar(p$trait)) traits_raw <- p$trait
traits <- character(0)
if (!is.null(traits_raw) && nzchar(traits_raw)) {
  traits <- trimws(unlist(strsplit(as.character(traits_raw), '[,;]+', perl=TRUE)))
  traits <- traits[nzchar(traits)]
}
maf <- if (!is.null(p$maf)) as.numeric(p$maf) else 0.05
missing_max <- if (!is.null(p$missing_max)) as.numeric(p$missing_max) else 0.1
use_lmm <- if (!is.null(p$use_lmm)) as.logical(p$use_lmm) else TRUE
cov_tsv <- if (!is.null(p$covariates_tsv)) p$covariates_tsv else NULL
if (!is.null(cov_tsv) && nchar(cov_tsv) == 0) cov_tsv <- NULL

stopifnot(!is.null(phenotype_tsv), file.exists(phenotype_tsv))

cat("[GWAS] phenotype_tsv=", phenotype_tsv, "\n")
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
# Load phenotype and align IDs (multi-trait capable)
# -----------------------
ph <- fread(phenotype_tsv, sep="\t", header=TRUE, data.table=FALSE)

# Helper: detect numeric traits (exclude id col and non-numeric columns)
detect_numeric_traits <- function(ph, id_col=1L, min_nonmiss=5L, max_new_na_frac=0.05) {
  cols <- colnames(ph)
  if (length(cols) <= id_col) return(character(0))
  out <- character(0)
  for (c in cols[-id_col]) {
    v <- ph[[c]]
    if (is.numeric(v) || is.integer(v)) {
      vv <- suppressWarnings(as.numeric(v))
      ok <- is.finite(vv)
      if (sum(ok) >= min_nonmiss && length(unique(vv[ok])) >= 2) out <- c(out, c)
    } else {
      v_chr <- as.character(v)
      suppressWarnings(v_num <- as.numeric(v_chr))
      n_orig <- sum(!is.na(v_chr) & nzchar(v_chr))
      n_new_na <- sum(is.na(v_num) & !is.na(v_chr) & nzchar(v_chr))
      ok <- is.finite(v_num)
      if (sum(ok) >= min_nonmiss && (n_new_na / max(1, n_orig)) <= max_new_na_frac && length(unique(v_num[ok])) >= 2) {
        out <- c(out, c)
      }
    }
  }
  out
}

if (length(traits) == 0) {
  traits <- detect_numeric_traits(ph, id_col=1L)
  cat('[GWAS] auto-detected traits n=', length(traits), '\n')
}
if (length(traits) == 0) stop('No numeric traits found in phenotype_tsv (excluding first id column). Please specify trait(s).')

# validate specified traits
missing_traits <- setdiff(traits, colnames(ph))
if (length(missing_traits) > 0) {
  cat('[GWAS] available phenotype columns:\n')
  cat(paste(colnames(ph), collapse=', '), '\n')
  stop(paste0('trait(s) not found: ', paste(missing_traits, collapse=', ')))
}

# Align genotype and phenotype IDs once (by first column)
gid <- as.character(bm@ped$id)
ord <- match(gid, ph[, 1])
keep <- !is.na(ord)
if (sum(keep) < 5) stop('Too few overlapping samples between genotype and phenotype (by id)')

bm_base <- bm[keep, ]
ph_base <- ph[ord[keep], , drop=FALSE]

# Optional covariates aligned to bm_base
Xcov_base <- NULL
if (!is.null(cov_tsv) && file.exists(cov_tsv)) {
  cv <- fread(cov_tsv, sep='\t', header=TRUE, data.table=FALSE)
  ord2 <- match(as.character(bm_base@ped$id), cv[, 1])
  keep2 <- !is.na(ord2)
  if (sum(keep2) != nrow(bm_base@ped)) {
    cat('[GWAS] WARN: covariates missing for some samples; dropping those samples\n')
    bm_base <- bm_base[keep2, ]
    ph_base <- ph_base[keep2, , drop=FALSE]
    ord2 <- ord2[keep2]
  }
  cv2 <- cv[ord2, , drop=FALSE]
  cov_cols <- colnames(cv2)[2:ncol(cv2)]
  if (length(cov_cols) > 0) {
    Xcov_base <- as.matrix(cv2[, cov_cols, drop=FALSE])
    Xcov_base <- apply(Xcov_base, 2, function(z) suppressWarnings(as.numeric(z)))
    Xcov_base[is.nan(Xcov_base)] <- NA
    if (any(!is.finite(Xcov_base))) {
      cat('[GWAS] WARN: NA/NaN in covariates; dropping rows with any NA\n')
      okc <- apply(Xcov_base, 1, function(r) all(is.finite(r)))
      bm_base <- bm_base[okc, ]
      ph_base <- ph_base[okc, , drop=FALSE]
      Xcov_base <- Xcov_base[okc, , drop=FALSE]
    }
    cat('[GWAS] covariates p=', ncol(Xcov_base), '\n')
  }
}

safe_name <- function(x) {
  s <- gsub('[^A-Za-z0-9._-]+', '_', as.character(x))
  if (nchar(s) == 0) s <- 'trait'
  if (nchar(s) > 64) s <- substr(s, 1, 64)
  s
}

subset_K_to_ids <- function(K, ids) {
  if (is.null(K)) return(NULL)
  if (is.list(K)) {
    return(lapply(K, function(k) subset_K_to_ids(k, ids)))
  }
  if (is.matrix(K) || is.data.frame(K)) {
    K <- as.matrix(K)
    rn <- rownames(K)
    cn <- colnames(K)
    if (!is.null(rn) && !is.null(cn) && all(ids %in% rn) && all(ids %in% cn)) {
      K <- K[ids, ids, drop=FALSE]
      return(K)
    }
  }
  # fallback: return as-is
  K
}

run_one_trait <- function(trait, out_dir_trait) {
  dir.create(out_dir_trait, recursive=TRUE, showWarnings=FALSE)
  plot_dir_trait <- file.path(out_dir_trait, 'plots')
  dir.create(plot_dir_trait, recursive=TRUE, showWarnings=FALSE)

  cat('[GWAS] trait=', trait, ' out=', out_dir_trait, '\n')

  # phenotype vector
  y <- suppressWarnings(as.numeric(ph_base[[trait]]))
  keep_y <- !is.na(y)
  if (sum(keep_y) < 5) {
    cat('[GWAS] WARN: too few non-missing phenotype values for trait=', trait, '\n')
    # create empty artifacts
    fwrite(data.frame(), file.path(out_dir_trait, 'results.tsv'), sep='\t')
    png(file.path(plot_dir_trait, 'manhattan.png'), width=1200, height=450)
    plot.new(); text(0.5, 0.5, paste0('No data for trait: ', trait))
    dev.off()
    return(list(trait=trait, out_dir=out_dir_trait, n=0L, m=0L, ok=FALSE))
  }

  bm_t <- bm_base[keep_y, ]
  y_t <- y[keep_y]

  X <- NULL
  if (is.null(Xcov_base)) {
    X <- matrix(1, nrow(bm_t))
  } else {
    Xcov_t <- Xcov_base[keep_y, , drop=FALSE]
    X <- cbind(rep(1, nrow(bm_t)), Xcov_t)
  }

  # QC filters using gaston fields
  miss_rate <- 1 - bm_t@snps$callrate
  maf_vec <- bm_t@snps$maf

  keep_snp <- rep(TRUE, nrow(bm_t@snps))
  if (!is.null(missing_max) && is.finite(missing_max)) keep_snp <- keep_snp & (miss_rate <= missing_max)
  if (!is.null(maf) && is.finite(maf)) keep_snp <- keep_snp & (maf_vec >= maf)

  cat('[GWAS] QC keep SNPs=', sum(keep_snp), '/', length(keep_snp), '\n')
  bm2 <- select.snps(bm_t, which(keep_snp))

  # GWAS
  if (method != 'FarmCPU') {
    if (!is.null(K_user)) {
      Kmat <- subset_K_to_ids(K_user, as.character(bm2@ped$id))
      eigenK <- eigen(Kmat, symmetric=TRUE)
    } else {
      cat('[GWAS] computing GRM from genotype (K = NULL)\n')
      Kmat <- GRM(bm2, autosome.only = F)
      eigenK <- eigen(Kmat, symmetric=TRUE)
    }
  }

  if (method == 'LM/LMM') {
    test <- if (!is.null(p$test)) as.character(p$test) else 'wald'
    pc_n <- if (!is.null(p$p)) as.integer(p$p) else 0
    if (!test %in% c('score','wald','lrt')) test <- 'wald'
    if (is.na(pc_n) || pc_n < 0) pc_n <- 0

    binary <- length(unique(y_t)) == 2
    response <- if (binary) 'binary' else 'quantitative'

    if (use_lmm) {
      if (binary) {
        if (test != 'score') {
          cat('[GWAS] WARN: forcing score test for binary LMM\n')
          test <- 'score'
        }
        ans <- association.test(bm2, Y=y_t, X=X, eigenK=eigenK, K=Kmat, method='lmm', response=response, test=test, p=pc_n)
      } else {
        if (test == 'score') {
          cat('[GWAS] WARN: forcing wald test for quantitative LMM\n')
          test <- 'wald'
        }
        ans <- association.test(bm2, Y=y_t, X=X, eigenK=eigenK, method='lmm', response=response, test=test, p=pc_n)
      }
    } else {
      if (test != 'wald') {
        cat('[GWAS] WARN: forcing wald test for LM\n')
        test <- 'wald'
      }
      ans <- association.test(bm2, Y=y_t, X=X, method='lm', response=response, test=test, p=pc_n, eigenK=eigenK)
    }
    saveRDS(ans, file.path(out_dir_trait, 'gaston_result.rds'))

    if (test == 'wald') {
      res <- data.frame(marker=bm2@snps$id, chr=bm2@snps$chr, pos=bm2@snps$pos, beta=ans$beta, sd=ans$sd, pvalue=ans$p)
    } else if (test == 'lrt') {
      res <- data.frame(marker=bm2@snps$id, chr=bm2@snps$chr, pos=bm2@snps$pos, LRT=ans$LRT, pvalue=ans$p)
    } else {
      res <- data.frame(marker=bm2@snps$id, chr=bm2@snps$chr, pos=bm2@snps$pos, score=ans$score, pvalue=ans$p)
    }

  } else if (method == 'MLMM') {
    if (!requireNamespace('mlmm', quietly=TRUE)) stop('Package mlmm is required for method=MLMM')
    max_steps <- if (!is.null(p$mlmm_max_steps)) as.integer(p$mlmm_max_steps) else 10
    if (ncol(X) == 1) {
      mlmm_out <- mlmm::mlmm(Y=y_t, X=as.matrix(bm2), K=Kmat, maxsteps=max_steps, nbchunks=2)
    } else {
      mlmm_out <- mlmm::mlmm_cof(Y=y_t, X=as.matrix(bm2), K=Kmat, maxsteps=max_steps, cof=X, nbchunks=2)
    }
    saveRDS(mlmm_out, file.path(out_dir_trait, 'mlmm_result.rds'))
    res <- data.frame(marker=as.character(bm2@snps$id), chr=bm2@snps$chr, pos=bm2@snps$pos,
                      pvalue=mlmm_out$opt_mbonf$out$pval, stringsAsFactors=FALSE)

  } else {
    if (!requireNamespace('FarmCPUpp', quietly=TRUE)) stop('Package FarmCPUpp is required for method=FarmCPU')
    if (!requireNamespace('bigmemory', quietly=TRUE)) stop('Package bigmemory is required for method=FarmCPU')
    BK <- bigmemory::as.big.matrix(as.matrix(bm2), type='double')
    map_df <- data.frame(snp=bm2@snps$id, chr=as.numeric(as.factor(bm2@snps$chr)), pos=bm2@snps$pos)
    gwas_df <- FarmCPUpp::farmcpu(Y=data.frame(id=bm2@ped$id, pheno=y_t), GD=BK, GM=map_df)
    saveRDS(gwas_df, file.path(out_dir_trait, 'farmcpu_result.rds'))
    res <- data.frame(marker=as.character(bm2@snps$id), chr=bm2@snps$chr, pos=bm2@snps$pos,
                      pvalue=gwas_df[[1]]$GWAS$p.value, stringsAsFactors=FALSE)
  }

  # Add BH-FDR q-values
  if (!('qvalue' %in% names(res)) && ('pvalue' %in% names(res))) {
    res$qvalue <- p.adjust(res$pvalue, method='BH')
  }

  res_path <- file.path(out_dir_trait, 'results.tsv')
  fwrite(res, res_path, sep='\t')

  # Manhattan plot
  png_path <- file.path(plot_dir_trait, 'manhattan.png')
  png(png_path, width=1200, height=450)
  dfp <- res
  dfp$chr <- as.character(dfp$chr)
  dfp$pos <- suppressWarnings(as.numeric(dfp$pos))
  dfp$mlogp <- -log10(dfp$pvalue)

  chr_key <- function(x) {
    suppressWarnings(v <- as.numeric(x))
    ifelse(is.na(v), 1e9, v)
  }
  chrs <- unique(dfp$chr[order(chr_key(dfp$chr))])

  offset <- 0
  dfp$x <- NA_real_
  ticks <- c(); ticklabs <- c()
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

  plot(dfp$x, dfp$mlogp, pch=20, cex=0.6, xlab='Chromosome', ylab='-log10(p)', main=paste0('GWAS (trait=', trait, ')'))
  axis(1, at=ticks, labels=ticklabs)
  dev.off()

  # QQ plot
  qq_png <- file.path(plot_dir_trait, 'qq.png')
  pp <- suppressWarnings(as.numeric(res$pvalue))
  pp <- pp[is.finite(pp) & pp > 0 & pp <= 1]
  if (length(pp) >= 10) {
    pp <- sort(pp)
    nqq <- length(pp)
    exp <- -log10((1:nqq)/(nqq+1))
    obs <- -log10(pp)
    exp <- rev(exp); obs <- rev(obs)
    lam <- NA_real_
    chisq <- qchisq(1-pp, df=1)
    if (length(chisq) > 0) lam <- median(chisq, na.rm=TRUE)/qchisq(0.5, df=1)
    png(qq_png, width=700, height=700)
    plot(exp, obs, pch=20, cex=0.6, xlab='Expected -log10(p)', ylab='Observed -log10(p)', main=paste0('QQ plot (trait=', trait, ')'))
    abline(0,1,col='red', lwd=2)
    if (is.finite(lam)) mtext(paste0('lambdaGC=', sprintf('%.3f', lam)), side=3, adj=1, line=-1)
    dev.off()

    qq_html <- file.path(plot_dir_trait, 'qq.html')
    if (!file.exists(qq_html)) writeLines('<html><body><p>Interactive QQ plot will be generated by GUI (Plotly).</p></body></html>', qq_html)
  }

  # placeholder for HTML paths
  html_path <- file.path(plot_dir_trait, 'manhattan.html')
  if (!file.exists(html_path)) writeLines('<html><body><p>Interactive Manhattan plot will be generated by GUI (Plotly).</p></body></html>', html_path)

  return(list(trait=trait, out_dir=out_dir_trait, n=nrow(bm2@ped), m=nrow(bm2@snps), ok=TRUE))
}

# Run traits
index <- data.frame(trait=character(0), out_subdir=character(0), n=integer(0), m=integer(0), ok=logical(0), stringsAsFactors=FALSE)

for (i in seq_along(traits)) {
  t <- traits[i]
  out_dir_trait <- out_dir
  out_sub <- '.'
  if (length(traits) > 1 && i > 1) {
    out_sub <- file.path('traits', safe_name(t))
    out_dir_trait <- file.path(out_dir, out_sub)
  }
  resi <- run_one_trait(t, out_dir_trait)
  index <- rbind(index, data.frame(trait=t, out_subdir=out_sub, n=as.integer(resi$n), m=as.integer(resi$m), ok=as.logical(resi$ok), stringsAsFactors=FALSE))
}

if (nrow(index) > 1) {
  fwrite(index, file.path(out_dir, 'traits_index.tsv'), sep='\t')
  writeLines(toJSON(list(traits=index$trait, default_trait=index$trait[1], index='traits_index.tsv'), auto_unbox=TRUE, pretty=TRUE),
             con=file.path(out_dir, 'traits.json'))
}

cat('[GWAS] done\n')
sink()
