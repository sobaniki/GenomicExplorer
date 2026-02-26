#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(jsonlite)
  library(data.table)
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
  stop("Usage: --params params.json --out out_dir")
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(out_dir, "run.log")
sink(log_file, split = TRUE)

cat("[impute_phenotypes] start\n")
cat("[impute_phenotypes] params_path=", params_path, "\n")
cat("[impute_phenotypes] out_dir=", out_dir, "\n")

p <- fromJSON(params_path)

# ---------------- inputs ----------------
phenotype_tsv <- if (!is.null(p$phenotype_tsv)) p$phenotype_tsv else NULL
kinship_tsv   <- if (!is.null(p$kinship_tsv)) p$kinship_tsv else NULL

method <- if (!is.null(p$method) && nchar(p$method) > 0) tolower(p$method) else "mean"

out_prefix <- if (!is.null(p$out_prefix) && nchar(p$out_prefix) > 0) p$out_prefix else file.path(out_dir, "imputed")

na_strings <- c("NA", "", ".", "-", "nan", "NaN", "NAN")
if (!is.null(p$na_strings) && length(p$na_strings) > 0) {
  na_strings <- unique(c(na_strings, as.character(p$na_strings)))
}

cat("[impute_phenotypes] method=", method, "\n")
cat("[impute_phenotypes] phenotype_tsv=", phenotype_tsv, "\n")
cat("[impute_phenotypes] kinship_tsv=", kinship_tsv, "\n")

if (is.null(phenotype_tsv) || !file.exists(phenotype_tsv)) {
  stop("phenotype_tsv is required and must exist")
}

# ---------------- helpers ----------------
write_artifacts <- function(lst) {
  writeLines(toJSON(lst, auto_unbox = TRUE, pretty = TRUE), file.path(out_dir, "artifacts.json"))
}

calc_missing_rate <- function(X) {
  mean(is.na(X))
}

load_pheno_tsv <- function(path) {
  dt <- fread(path, sep = "\t", header = TRUE, data.table = FALSE, na.strings = na_strings)
  #if (!("id" %in% colnames(dt))) stop("phenotype_tsv must have an 'id' column")
  #ids <- as.character(dt$id)
  ids <- dt[, 1]
  #traits <- setdiff(colnames(dt), "id")
  traits <- colnames(dt)[2:ncol(dt)]
  if (length(traits) == 0) stop("No trait columns in phenotype_tsv")
  Y <- as.matrix(dt[, traits, drop = FALSE])
  suppressWarnings(storage.mode(Y) <- "numeric")

  # reconcile row mismatch if any
  if (length(ids) != nrow(Y)) {
    cat(sprintf("[WARN] phenotype_tsv rows mismatch: length(id)=%d but nrow(Y)=%d. ", length(ids), nrow(Y)))
    if (length(ids) > nrow(Y)) {
      cat("Truncating ids to nrow(Y).\n")
      ids <- ids[seq_len(nrow(Y))]
    } else {
      cat("Padding ids to nrow(Y).\n")
      ids <- c(ids, paste0("sample_", seq_len(nrow(Y) - length(ids)) + length(ids)))
    }
  }

  rownames(Y) <- ids
  list(ids = ids, traits = traits, Y = Y)
}

save_pheno_tsv <- function(ids, traits, Y, path) {
  dt <- data.table(id = ids)
  dt <- cbind(dt, as.data.table(Y))
  setnames(dt, c("id", traits))
  fwrite(dt, path, sep = "\t", na = "NA")
}

load_kinship_tsv <- function(path) {
  # expected: wide square matrix with first col 'id' (row ids) and column headers are sample ids
  Kdt <- fread(path, sep = "\t", header = TRUE, data.table = FALSE, na.strings = na_strings)
  #if (!("id" %in% colnames(Kdt))) stop("kinship_tsv must have an 'id' column as the first column")
  #row_ids <- as.character(Kdt$id)
  row_ids <- as.character(Kdt[, 1])
  #col_ids <- setdiff(colnames(Kdt), "id")
  col_ids <- colnames(Kdt)[2:ncol(Kdt)]
  if (length(col_ids) == 0) stop("kinship_tsv has no sample columns")
  K <- as.matrix(Kdt[, col_ids, drop = FALSE])
  suppressWarnings(storage.mode(K) <- "numeric")
  rownames(K) <- row_ids
  colnames(K) <- col_ids
  K
}

make_trait_summary <- function(traits, Y_before, Y_after) {
  miss_before <- colMeans(is.na(Y_before))
  miss_after  <- colMeans(is.na(Y_after))
  mu_before <- apply(Y_before, 2, function(x) mean(x, na.rm = TRUE))
  sd_before <- apply(Y_before, 2, function(x) sd(x, na.rm = TRUE))
  mu_after <- apply(Y_after, 2, function(x) mean(x, na.rm = TRUE))
  sd_after <- apply(Y_after, 2, function(x) sd(x, na.rm = TRUE))
  data.frame(
    trait = traits,
    missing_before = miss_before,
    missing_after  = miss_after,
    mean_before = mu_before,
    sd_before = sd_before,
    mean_after  = mu_after,
    sd_after  = sd_after,
    stringsAsFactors = FALSE
  )
}

make_sample_summary <- function(ids, Y_before, Y_after) {
  miss_before <- rowMeans(is.na(Y_before))
  miss_after  <- rowMeans(is.na(Y_after))
  data.frame(
    id = ids,
    missing_before = miss_before,
    missing_after = miss_after,
    stringsAsFactors = FALSE
  )
}

# ---------------- main ----------------
dat <- load_pheno_tsv(phenotype_tsv)
ids <- dat$ids
traits <- dat$traits
Y0 <- dat$Y

cat("[impute_phenotypes] n_samples=", nrow(Y0), " n_traits=", ncol(Y0), "\n")
cat("[impute_phenotypes] missing_rate(before)=", calc_missing_rate(Y0), "\n")

Yimp <- Y0
extra_artifacts <- list()

if (method %in% c("mean", "median")) {
  if (method == "mean") {
    m <- colMeans(Y0, na.rm = TRUE)
    for (j in seq_along(m)) Yimp[is.na(Yimp[, j]), j] <- m[j]
  } else {
    med <- apply(Y0, 2, function(x) median(x, na.rm = TRUE))
    for (j in seq_along(med)) Yimp[is.na(Yimp[, j]), j] <- med[j]
  }

} else if (method %in% c("missforest", "rf")) {
  if (!requireNamespace("missForest", quietly = TRUE)) {
    stop("R package 'missForest' is not installed. Install it to use method=missForest.")
  }
  # parameters
  maxiter <- if (!is.null(p$rf_maxiter)) as.integer(p$rf_maxiter) else 10L
  ntree   <- if (!is.null(p$rf_ntree)) as.integer(p$rf_ntree) else 100L
  seed    <- if (!is.null(p$rf_seed)) as.integer(p$rf_seed) else 1L

  set.seed(seed)
  df <- as.data.frame(Y0)
  out <- missForest::missForest(
    xmis = df,
    maxiter = maxiter,
    ntree = ntree,
    verbose = TRUE
  )
  Yimp <- as.matrix(out$ximp)
  rownames(Yimp) <- ids
  colnames(Yimp) <- traits
  extra_artifacts$missForest_OOBerror <- out$OOBerror

} else if (method == "mice") {
  if (!requireNamespace("mice", quietly = TRUE)) {
    stop("R package 'mice' is not installed. Install it to use method=mice.")
  }
  m <- if (!is.null(p$mice_m)) as.integer(p$mice_m) else 5L
  maxit <- if (!is.null(p$mice_maxit)) as.integer(p$mice_maxit) else 10L
  seed <- if (!is.null(p$mice_seed)) as.integer(p$mice_seed) else 1L
  printFlag <- if (!is.null(p$mice_printFlag)) isTRUE(p$mice_printFlag) else FALSE

  df <- as.data.frame(Y0)
  set.seed(seed)
  imp <- mice::mice(df, m = m, maxit = maxit, printFlag = printFlag)
  # return one completed dataset (the 1st); users can re-run with different seeds/m
  dfc <- mice::complete(imp, action = 1)
  Yimp <- as.matrix(dfc)
  rownames(Yimp) <- ids
  colnames(Yimp) <- traits
  extra_artifacts$mice_m <- m
  extra_artifacts$mice_maxit <- maxit

} else if (method == "phenix") {
  if (is.null(kinship_tsv) || !file.exists(kinship_tsv)) {
    stop("method=phenix requires kinship_tsv (a square matrix with id column)")
  }
  if (!requireNamespace("phenix", quietly = TRUE)) {
    stop("R package 'phenix' is not installed. Install it to use method=phenix.")
  }

  K <- load_kinship_tsv(kinship_tsv)

  # Align IDs between Y and K
  common <- intersect(ids, intersect(rownames(K), colnames(K)))
  if (length(common) < 3) {
    stop("Too few overlapping sample IDs between phenotype_tsv and kinship_tsv")
  }
  if (length(common) < length(ids)) {
    cat(sprintf("[WARN] phenix: using %d/%d samples that overlap between phenotype and kinship\n", length(common), length(ids)))
  }
  # subset and order
  Y_sub <- Y0[common, , drop = FALSE]
  K_sub <- K[common, common]

  quantnorm <- if (!is.null(p$phenix_quantnorm)) isTRUE(p$phenix_quantnorm) else FALSE
  scale_cols <- if (!is.null(p$phenix_scale)) isTRUE(p$phenix_scale) else TRUE
  trim <- if (!is.null(p$phenix_trim)) isTRUE(p$phenix_trim) else FALSE
  trim_sds <- if (!is.null(p$phenix_trim_sds)) as.numeric(p$phenix_trim_sds) else 4
  seed <- if (!is.null(p$phenix_seed)) as.integer(p$phenix_seed) else 8473L
  maxit <- if (!is.null(p$phenix_maxit)) as.integer(p$phenix_maxit) else 1000L
  reltol <- if (!is.null(p$phenix_reltol)) as.numeric(p$phenix_reltol) else 1e-8
  tau <- if (!is.null(p$phenix_tau)) as.numeric(p$phenix_tau) else 0

  set.seed(seed)
  out <- phenix::phenix(
    Y = Y_sub,
    K = K_sub,
    test = FALSE,
    seed = seed,
    quantnorm = quantnorm,
    scale = scale_cols,
    trim = trim,
    trim.sds = trim_sds,
    maxit = maxit,
    reltol = reltol,
    tau = tau
  )

  # out$imp is imputed phenotype matrix
  Yimp_sub <- out$imp
  rownames(Yimp_sub) <- common
  colnames(Yimp_sub) <- traits

  # Put back into full matrix (if some ids not used, keep as NA there)
  Yimp[,] <- NA_real_
  Yimp[common, ] <- Yimp_sub

  # save extra outputs
  if (!is.null(out$h2)) {
    h2 <- data.frame(trait = traits, h2 = as.numeric(out$h2), stringsAsFactors = FALSE)
    fwrite(h2, file.path(out_dir, "phenix_h2.tsv"), sep = "\t")
    extra_artifacts$phenix_h2_tsv <- file.path(out_dir, "phenix_h2.tsv")
  }
  if (!is.null(out$U)) {
    U <- out$U
    rownames(U) <- common
    colnames(U) <- traits
    # write as TSV with id
    dtU <- data.table(id = common)
    dtU <- cbind(dtU, as.data.table(U))
    setnames(dtU, c("id", traits))
    fwrite(dtU, file.path(out_dir, "phenix_breeding_values_U.tsv"), sep = "\t")
    extra_artifacts$phenix_U_tsv <- file.path(out_dir, "phenix_breeding_values_U.tsv")
  }

  extra_artifacts$phenix_used_samples <- length(common)
  extra_artifacts$phenix_seed <- seed

} else if (method == "softimpute") {
  if (!requireNamespace("softImpute", quietly = TRUE)) {
    stop("R package 'softImpute' is not installed. Install it to use method=softImpute.")
  }
  # softImpute expects matrix; center/scale optional
  rank_max <- if (!is.null(p$si_rank_max)) as.integer(p$si_rank_max) else 50L
  lambda <- if (!is.null(p$si_lambda)) as.numeric(p$si_lambda) else 30
  maxit <- if (!is.null(p$si_maxit)) as.integer(p$si_maxit) else 100L
  type <- if (!is.null(p$si_type)) as.character(p$si_type) else "svd"

  fit <- softImpute::softImpute(Y0, rank.max = rank_max, lambda = lambda, maxit = maxit, type = type,
                                thresh = 1e-05,
                                trace.it = F,
                                warm.start = NULL,
                                final.svd = T)
  Yimp <- softImpute::complete(Y0, fit)
  rownames(Yimp) <- ids
  colnames(Yimp) <- traits
  extra_artifacts$soft_rank_max <- rank_max
  extra_artifacts$soft_lambda <- lambda
  extra_artifacts$soft_type <- type
  extra_artifacts$soft_maxit <- maxit

} else {
  stop(sprintf("Unknown method '%s'. Supported: mean, median, missForest, mice, phenix, softImpute", method))
}

# write outputs
#out_tsv <- paste0(out_prefix, ".tsv")
out_tsv <- file.path(out_dir, "pheno.imp.tsv")
save_pheno_tsv(ids, traits, Yimp, out_tsv)

# summaries
trait_sum <- make_trait_summary(traits, Y0, Yimp)
sample_sum <- make_sample_summary(ids, Y0, Yimp)
fwrite(trait_sum, file.path(out_dir, "trait_summary.tsv"), sep = "\t")
fwrite(sample_sum, file.path(out_dir, "sample_summary.tsv"), sep = "\t")

# imputed mask
mask <- is.na(Y0) & !is.na(Yimp)
mask_dt <- data.table(id = ids)
mask_dt <- cbind(mask_dt, as.data.table(apply(mask, 2, as.integer)))
setnames(mask_dt, c("id", traits))
fwrite(mask_dt, file.path(out_dir, "imputed_mask.tsv"), sep = "\t")

write_artifacts(c(list(
  method = method,
  phenotype_tsv_in = phenotype_tsv,
  kinship_tsv_in = kinship_tsv,
  imputed_tsv = out_tsv,
  n_samples = nrow(Y0),
  n_traits = ncol(Y0),
  missing_before = calc_missing_rate(Y0),
  missing_after = calc_missing_rate(Yimp),
  trait_summary_tsv = file.path(out_dir, "trait_summary.tsv"),
  sample_summary_tsv = file.path(out_dir, "sample_summary.tsv"),
  imputed_mask_tsv = file.path(out_dir, "imputed_mask.tsv")
), extra_artifacts))

cat("[impute_phenotypes] missing_rate(after)=", calc_missing_rate(Yimp), "\n")
cat("[impute_phenotypes] output=", out_tsv, "\n")
cat("[impute_phenotypes] done\n")

