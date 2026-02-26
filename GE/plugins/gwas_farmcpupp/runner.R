#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(jsonlite)
  library(data.table)
  library(gaston)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) return(args[i + 1])
  default
}

params_path <- get_arg("--params")
out_dir <- get_arg("--out")
if (is.null(params_path) || is.null(out_dir)) {
  stop("Usage: runner.R --params params.json --out out_dir")
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
plot_dir <- file.path(out_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_dir, "run.log")
sink(log_file, split = TRUE)

write_error_artifacts <- function(msg) {
  cat("[gwas_farmcpupp] ERROR:", msg, "\n")
  fwrite(data.frame(
    snp = character(0), chr = character(0), pos = integer(0),
    pval = numeric(0), beta = numeric(0), se = numeric(0), method = character(0)
  ), file.path(out_dir, "results.tsv"), sep = "\t")
  png(file.path(plot_dir, "manhattan.png"), width = 1200, height = 600)
  plot.new(); text(0.5, 0.5, paste("FarmCPUpp failed:\n", msg))
  dev.off()
  writeLines(msg, file.path(out_dir, "error_message.txt"))
  art <- list(
    status = "error",
    message = msg,
    table = "results.tsv",
    plot = "plots/manhattan.png",
    method = "FarmCPUpp"
  )
  writeLines(toJSON(art, auto_unbox = TRUE, pretty = TRUE), file.path(out_dir, "artifacts.json"))
}

make_manhattan <- function(df, out_png) {
  df <- df[is.finite(df$pval) & df$pval > 0, , drop = FALSE]
  png(out_png, width = 1400, height = 700)
  if (nrow(df) == 0) {
    plot.new(); text(0.5, 0.5, "no results")
    dev.off();
    return(invisible(NULL))
  }
  x <- seq_len(nrow(df))
  y <- -log10(df$pval)
  plot(x, y, pch = 20, cex = 0.6, main = "Manhattan (FarmCPUpp)", xlab = "marker", ylab = "-log10(p)")
  dev.off()
}


cat("[gwas_farmcpupp] start\n")
cat("[gwas_farmcpupp] params_path=", params_path, "\n")
cat("[gwas_farmcpupp] out_dir=", out_dir, "\n")

p <- fromJSON(params_path)

# -----------------------
# Params (common)
# -----------------------
phenotype_tsv <- p$phenotype_tsv
trait <- p$trait
maf <- if (!is.null(p$maf)) as.numeric(p$maf) else 0.05
missing_max <- if (!is.null(p$missing_max)) as.numeric(p$missing_max) else 0.1
cov_tsv <- if (!is.null(p$covariates_tsv)) p$covariates_tsv else NULL
if (!is.null(cov_tsv) && nchar(cov_tsv) == 0) cov_tsv <- NULL

stopifnot(!is.null(phenotype_tsv), file.exists(phenotype_tsv))
stopifnot(!is.null(trait), nchar(trait) > 0)

cat("[FarmCPU] phenotype_tsv=", phenotype_tsv, "\n")
cat("[FarmCPU] trait=", trait, "\n")
cat("[FarmCPU] maf=", maf, " missing_max=", missing_max, " use_lmm=", use_lmm, "\n")
cat("[FarmCPU] covariates_tsv=", cov_tsv, "\n")

# -----------------------
# Input switch
# -----------------------
genotype_tsv <- if (!is.null(p$genotype_tsv)) p$genotype_tsv else NULL
marker_map_tsv <- if (!is.null(p$marker_map_tsv)) p$marker_map_tsv else NULL

cat("[FarmCPU] genotype_file=", genotype_tsv, "\n")
cat("[FarmCPU] marker_map_tsv=", marker_map_tsv, "\n")

# -----------------------
# Load genotype as bed.matrix
# -----------------------
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

cat("[FarmCPU] genotype loaded. n=", nrow(bm@ped), " m=", nrow(bm@snps), "\n")

# -----------------------
# Load phenotype and align IDs
# -----------------------
ph <- fread(phenotype_tsv, sep="\t", header=TRUE, data.table=FALSE)
stopifnot(trait %in% colnames(ph))

# bm sample IDs: bm@ped$id (gaston?d?l)
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

cat("[FarmCPU] aligned samples n=", length(y), "\n")

# -----------------------
# Optional covariates
# -----------------------
Xcov <- NULL
if (!is.null(cov_tsv) && file.exists(cov_tsv)) {
  cv <- fread(cov_tsv, sep="\t", header=TRUE, data.table=FALSE)
  stopifnot("id" %in% colnames(cv))
  cv$id <- as.character(cv$id)
  
  ord2 <- match(as.character(bm@ped$id), cv$id)
  keep2 <- !is.na(ord2)
  if (sum(keep2) != nrow(bm@ped)) {
    cat("[FarmCPU] WARN: covariates missing for some samples; dropping those samples\n")
    bm <- bm[keep2, ]
    ph2 <- ph2[keep2, , drop=FALSE]
    y <- y[keep2]
    ord2 <- ord2[keep2]
  }
  
  cv2 <- cv[ord2, , drop=FALSE]
  cov_cols <- setdiff(colnames(cv2), "id")
  if (length(cov_cols) > 0) {
    Xcov <- as.matrix(cv2[, cov_cols, drop=FALSE])
    Xcov <- apply(Xcov, 2, function(z) suppressWarnings(as.numeric(z)))
    Xcov[is.nan(Xcov)] <- NA
    if (any(!is.finite(Xcov))) {
      cat("[FarmCPU] WARN: NA/NaN in covariates; dropping rows with any NA\n")
      okc <- apply(Xcov, 1, function(r) all(is.finite(r)))
      bm <- bm[okc, ]
      ph2 <- ph2[okc, , drop=FALSE]
      y <- y[okc]
      Xcov <- Xcov[okc, , drop=FALSE]
    }
    cat("[FarmCPU] covariates p=", ncol(Xcov), "\n")
  }
}
X <- NULL
if (is.null(Xcov)) {
  #X <- matrix(1, nrow(bm))
  X <- NULL
} else {
  #X <- cbind(rep(1, nrow(bm)), Xcov)
  X <- Xcov
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

cat("[FarmCPU] QC keep SNPs=", sum(keep_snp), "/", length(keep_snp), "\n")
bm2 <- select.snps(bm, which(keep_snp))

BK <- bigmemory::as.big.matrix(as.matrix(bm2), type = "double")
map_df <- data.frame(snp = bm2@snps$id,
                     chr = bm2@snps$chr,
                     pos = bm2@snps$pos)
  
gwas_df <- FarmCPUpp::farmcpu(Y = data.frame(id = bm2@ped$id,
                                             pheno = y2),
                              GD = BK,
                              GM = map_df,
                              CV = X,
                              GP = NULL,
                              method.sub = as.character(p$farmcpu_method),
                              method.sub.final = as.character(p$farmcpu_method),
                              method.bin = as.character(p$farmcpu_method_bin),
                              bin.size = c(5e+05, 5e+06, 5e+07), 
                              bin.selection = seq(10, 100, 10),
                              memo = NULL, 
                              Prior = NULL, 
                              ncores.glm = 1, 
                              maxLoop = as.numeric(p$farmcpu_maxloop),
                              converge = 1, 
                              iteration.output = F, 
                              p.threshold = as.numeric(p$farmcpu_p_threshold),
                              ncores.reml = 1, 
                              threshold = 0.7)

fwrite(gwas_df, file.path(out_dir, "results.tsv"), sep = "\t", quote = F)
make_manhattan(gwas_df, file.path(plot_dir, "manhattan.png"))
# placeholder for HTML path so GUI can detect it (created later)
html_path <- file.path(plot_dir, "manhattan.html")
if (!file.exists(html_path)) {
  writeLines("<html><body><p>Interactive Manhattan plot will be generated by GUI (Plotly).</p></body></html>", html_path)
}

# QQ plot (simple)
qq_png <- file.path(plot_dir, "qq.png")
pp <- NA
if ("pval" %in% names(gwas_df)) {
  pp <- suppressWarnings(as.numeric(gwas_df$pval))
} else if ("pvalue" %in% names(gwas_df)) {
  pp <- suppressWarnings(as.numeric(gwas_df$pvalue))
} else if ("P.value" %in% names(gwas_df)) {
  pp <- suppressWarnings(as.numeric(gwas_df$P.value))
}
pp <- pp[is.finite(pp) & pp > 0 & pp <= 1]
if (length(pp) >= 10) {
  pp <- sort(pp)
  nqq <- length(pp)
  exp <- -log10((1:nqq) / (nqq + 1))
  obs <- -log10(pp)
  exp <- rev(exp)
  obs <- rev(obs)
  
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

art <- list(
  status = status,
  message = message,
  table = "results.tsv",
  plot = "plots/manhattan.png",
  method = "FarmCPUpp"
)
writeLines(toJSON(art, auto_unbox = TRUE, pretty = TRUE), file.path(out_dir, "artifacts.json"))

# Always exit 0 so the GUI can proceed; errors are recorded in outputs.
quit(status = 0)
