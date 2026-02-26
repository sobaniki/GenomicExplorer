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
if (is.null(params_path) || is.null(out_dir)) stop("Usage: --params params.json --out out_dir")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(out_dir, "run.log")
sink(log_file, split = TRUE)

cat("[gp_glmnet] start\n")
cat("[gp_glmnet] params_path=", params_path, "\n")
cat("[gp_glmnet] out_dir=", out_dir, "\n")

p <- fromJSON(params_path)

geno_tsv <- p$genotype_tsv
pheno_tsv <- p$phenotype_tsv

# Trait selection:
# - Default: ALL traits in phenotype.tsv excluding the first ID column.
# - Optional spec: indices (1-based excluding ID) like "1", "1-5", "2,3,5".
# - You may also specify trait names matching the header.
traits_spec <- ""
if (!is.null(p$traits_spec) && nchar(as.character(p$traits_spec)) > 0) traits_spec <- as.character(p$traits_spec)
if (!is.null(p$traits) && nchar(as.character(p$traits)) > 0) traits_spec <- as.character(p$traits)
if (!is.null(p$trait) && nchar(as.character(p$trait)) > 0) traits_spec <- as.character(p$trait)

stopifnot(!is.null(geno_tsv), file.exists(geno_tsv))
stopifnot(!is.null(pheno_tsv), file.exists(pheno_tsv))

mode <- if (!is.null(p$mode) && nchar(p$mode) > 0) as.character(p$mode) else "fit"  # fit / loo / kfold
k_folds <- if (!is.null(p$k_folds)) as.integer(p$k_folds) else 5
seed <- if (!is.null(p$seed)) as.integer(p$seed) else 1

center <- if (!is.null(p$center_markers)) as.logical(p$center_markers) else TRUE
scale <- if (!is.null(p$scale_markers)) as.logical(p$scale_markers) else TRUE
impute <- if (!is.null(p$impute_missing)) as.character(p$impute_missing) else "mean"  # mean / median

penalty <- NULL
if (!is.null(p$penalty) && nchar(as.character(p$penalty)) > 0) penalty <- as.character(p$penalty)
if (is.null(penalty) && !is.null(p$model) && nchar(as.character(p$model)) > 0) penalty <- as.character(p$model)

alpha <- NA_real_
if (!is.null(p$alpha)) {
  suppressWarnings(alpha <- as.numeric(as.character(p$alpha)))
}
if (is.na(alpha)) {
  if (!is.null(penalty) && tolower(penalty) %in% c("ridge", "lasso")) {
    alpha <- if (tolower(penalty) == "ridge") 0 else 1
  } else {
    # default: ridge
    alpha <- 0
    penalty <- "ridge"
  }
}
if (is.null(penalty) || !nzchar(penalty)) {
  penalty <- if (alpha == 0) "ridge" else "lasso"
}

lambda_mode <- if (!is.null(p$lambda_mode) && nchar(as.character(p$lambda_mode)) > 0) as.character(p$lambda_mode) else "cv.min"
lambda_value <- NA_real_
if (!is.null(p$lambda_value) && nchar(as.character(p$lambda_value)) > 0) {
  suppressWarnings(lambda_value <- as.numeric(as.character(p$lambda_value)))
}
if (!is.null(p$lambda) && nchar(as.character(p$lambda)) > 0) {
  suppressWarnings(lambda_value <- as.numeric(as.character(p$lambda)))
}

nlambda <- if (!is.null(p$nlambda)) as.integer(p$nlambda) else 100L
lambda_min_ratio <- NA_real_
if (!is.null(p$lambda_min_ratio) && nchar(as.character(p$lambda_min_ratio)) > 0) {
  suppressWarnings(lambda_min_ratio <- as.numeric(as.character(p$lambda_min_ratio)))
}
if (!is.finite(lambda_min_ratio) || is.na(lambda_min_ratio) || lambda_min_ratio <= 0) {
  lambda_min_ratio <- 1e-4
}

nfolds_fit <- if (!is.null(p$nfolds_fit)) as.integer(p$nfolds_fit) else 10L
top_coef_n <- if (!is.null(p$top_coef_n)) as.integer(p$top_coef_n) else 100L

cat("[gp_glmnet] genotype_tsv=", geno_tsv, "\n")
cat("[gp_glmnet] phenotype_tsv=", pheno_tsv, "\n")
cat("[gp_glmnet] traits_spec=", traits_spec, "\n")
cat("[gp_glmnet] mode=", mode, " k_folds=", k_folds, " seed=", seed, "\n")
cat("[gp_glmnet] penalty=", penalty, " alpha=", alpha, "\n")
cat("[gp_glmnet] lambda_mode=", lambda_mode, " lambda_value=", lambda_value, " nlambda=", nlambda, " lambda_min_ratio=", lambda_min_ratio, "\n")
cat("[gp_glmnet] center=", center, " scale=", scale, " impute=", impute, "\n")

if (!(mode %in% c("fit", "loo", "kfold"))) stop("mode must be one of: fit, loo, kfold")
if (mode == "kfold" && (is.na(k_folds) || k_folds < 2)) stop("k_folds must be >= 2")
if (!(lambda_mode %in% c("cv.min", "cv.1se", "fixed"))) stop("lambda_mode must be one of: cv.min, cv.1se, fixed")
if (lambda_mode == "fixed" && !is.finite(lambda_value)) stop("lambda_mode=fixed requires lambda_value")

rmse <- function(a, b) {
  ok <- is.finite(a) & is.finite(b)
  if (!any(ok)) return(NA_real_)
  sqrt(mean((a[ok] - b[ok])^2))
}
mae <- function(a, b) {
  ok <- is.finite(a) & is.finite(b)
  if (!any(ok)) return(NA_real_)
  mean(abs(a[ok] - b[ok]))
}
pearson <- function(a, b) {
  ok <- is.finite(a) & is.finite(b)
  if (sum(ok) < 3) return(NA_real_)
  suppressWarnings(cor(a[ok], b[ok], method = "pearson"))
}

impute_vec <- function(v, method = "mean") {
  if (!anyNA(v)) return(v)
  if (method == "median") {
    mu <- stats::median(v, na.rm = TRUE)
  } else {
    mu <- mean(v, na.rm = TRUE)
  }
  if (!is.finite(mu)) mu <- 0
  v[is.na(v)] <- mu
  v
}

parse_trait_spec <- function(spec, trait_cols) {
  if (is.null(spec)) spec <- ""
  spec <- trimws(as.character(spec))
  if (nchar(spec) == 0) return(seq_along(trait_cols))
  spec2 <- gsub("\\s+", "", spec)
  toks <- unlist(strsplit(spec2, ","))
  idxs <- integer(0)
  for (tk in toks) {
    if (!nchar(tk)) next
    if (grepl("^[0-9]+-[0-9]+$", tk)) {
      ab <- as.integer(unlist(strsplit(tk, "-")))
      a <- ab[1]; b <- ab[2]
      if (is.na(a) || is.na(b)) stop(paste0("Invalid trait range: ", tk))
      if (a > b) { tmp <- a; a <- b; b <- tmp }
      idxs <- c(idxs, seq.int(a, b))
    } else if (grepl("^[0-9]+$", tk)) {
      idxs <- c(idxs, as.integer(tk))
    } else {
      j <- match(tk, trait_cols)
      if (is.na(j)) stop(paste0("Trait name not found in phenotype_tsv header: ", tk))
      idxs <- c(idxs, j)
    }
  }
  idxs <- idxs[!is.na(idxs)]
  if (!length(idxs)) stop("No valid traits parsed from trait(s) spec.")
  idxs <- idxs[!duplicated(idxs)]
  if (any(idxs < 1 | idxs > length(trait_cols))) {
    stop(paste0("Trait index out of range. Allowed: 1..", length(trait_cols)))
  }
  idxs
}

set.seed(seed)

# Load genotype
Gt <- fread(geno_tsv, sep = "\t", header = TRUE, data.table = FALSE)
if (ncol(Gt) < 2) stop("genotype_tsv must have at least 2 columns: <ID> + markers")
ids_g <- as.character(Gt[[1]])
marker_cols <- colnames(Gt)[-1]
if (length(marker_cols) == 0) stop("No marker columns in genotype_tsv")
X <- as.matrix(Gt[, marker_cols, drop = FALSE])
X <- apply(X, 2, function(z) suppressWarnings(as.numeric(z)))
X[is.nan(X)] <- NA

# Load phenotype
Pt <- fread(pheno_tsv, sep = "\t", header = TRUE, data.table = FALSE)
if (ncol(Pt) < 2) stop("phenotype_tsv must have at least 2 columns: <ID> + trait(s)")
ids_p <- as.character(Pt[[1]])
trait_cols <- colnames(Pt)[-1]

trait_idxs <- parse_trait_spec(traits_spec, trait_cols)
traits <- trait_cols[trait_idxs]

# Align
m <- match(ids_g, ids_p)
ok <- !is.na(m)
if (!any(ok)) stop("No overlapping IDs between genotype and phenotype")
X <- X[ok, , drop = FALSE]
ids <- ids_g[ok]

cat("[gp_glmnet] n_samples=", length(ids), " n_markers=", ncol(X), "\n")

# Optional covariates (id + covariate columns)
cov_path <- p$covariate_tsv
if (!is.null(cov_path) && nzchar(cov_path) && file.exists(cov_path)) {
  cat(sprintf("[gp_glmnet] covariate_tsv=%s\n", cov_path))
  C <- fread(cov_path, sep = "\t", header = TRUE, data.table = FALSE)
  if (ncol(C) < 2) stop("covariate_tsv must have at least 2 columns: <ID> + covariates")
  cov_cols <- colnames(C)[-1]
  if (length(cov_cols) > 0) {
    ids_c <- as.character(C[[1]])
    mm <- match(ids, ids_c)
    Z <- as.matrix(C[mm, cov_cols, drop = FALSE])
    Z <- apply(Z, 2, function(z) suppressWarnings(as.numeric(z)))
    Z[is.nan(Z)] <- NA
    for (j in seq_len(ncol(Z))) Z[, j] <- impute_vec(Z[, j], method = impute)
    colnames(Z) <- paste0("cov_", cov_cols)
    X <- cbind(Z, X)
  }
}

# Impute missing markers
for (j in seq_len(ncol(X))) {
  X[, j] <- impute_vec(X[, j], method = impute)
}

# Optional center/scale
if (center || scale) {
  X <- as.matrix(scale(X, center = center, scale = scale))
}
X <- t(na.omit(t(X)))

if (!requireNamespace("glmnet", quietly = TRUE)) {
  stop("R package 'glmnet' is required but not installed. Please run: install.packages('glmnet')")
}

choose_lambda_s <- function(cvfit) {
  if (lambda_mode == "fixed") return(lambda_value)
  if (lambda_mode == "cv.1se") return(cvfit$lambda.1se)
  cvfit$lambda.min
}

lambda_index <- function(lambda_seq, s) {
  if (!is.finite(s) || is.na(s)) return(NA_integer_)
  # exact match is best, otherwise nearest
  j <- which(lambda_seq == s)
  if (length(j) > 0) return(j[1])
  which.min(abs(lambda_seq - s))
}

fit_and_predict_trait <- function(y) {
  n <- length(y)
  pred <- rep(NA_real_, n)
  fold_id <- rep(0L, n)
  fold_metrics <- NULL

  obs <- which(is.finite(y))
  if (length(obs) < 3) {
    cat("[gp_glmnet] WARN: too few observed phenotypes. Using mean predictor. n_obs=", length(obs), "\n")
    mu <- mean(y[obs], na.rm = TRUE)
    if (!is.finite(mu)) mu <- 0
    pred[] <- mu
    return(list(
      pred = pred,
      fold_id = fold_id,
      fold_metrics = NULL,
      cor = pearson(y, pred),
      rmse = rmse(y, pred),
      mae = mae(y, pred),
      n = n,
      n_obs = length(obs),
      lambda = NA_real_,
      lambda_mode = lambda_mode,
      coef = NULL
    ))
  }

  Xo <- X[obs, , drop = FALSE]
  yo <- as.numeric(y[obs])

  # Fit (and possibly CV) on observed data
  cvfit <- NULL
  full_fit <- NULL
  s <- NA_real_

  if (mode == "fit") {
    if (lambda_mode == "fixed") {
      full_fit <- glmnet::glmnet(
        x = Xo, y = yo, alpha = alpha,
        lambda = lambda_value,
        intercept = TRUE, standardize = FALSE
      )
      s <- lambda_value
    } else {
      nf <- min(max(3L, nfolds_fit), length(obs))
      cvfit <- glmnet::cv.glmnet(
        x = Xo, y = yo, alpha = alpha,
        nfolds = nf,
        keep = FALSE,
        intercept = TRUE, standardize = FALSE,
        nlambda = nlambda,
        lambda.min.ratio = lambda_min_ratio
      )
      s <- choose_lambda_s(cvfit)
      full_fit <- cvfit$glmnet.fit
    }
    pred <- as.numeric(predict(full_fit, newx = X, s = s))
  } else if (mode == "kfold") {
    # CV predictions for observed samples; missing-y samples get predictions from the full fit.
    perm <- sample.int(length(obs))
    foldid_o <- rep(seq_len(k_folds), length.out = length(obs))
    foldid_o[perm] <- foldid_o
    if (lambda_mode == "fixed") {
      cvfit <- glmnet::cv.glmnet(
        x = Xo, y = yo, alpha = alpha,
        foldid = foldid_o,
        keep = TRUE,
        intercept = TRUE, standardize = FALSE,
        lambda = lambda_value
      )
    } else {
      cvfit <- glmnet::cv.glmnet(
        x = Xo, y = yo, alpha = alpha,
        foldid = foldid_o,
        keep = TRUE,
        intercept = TRUE, standardize = FALSE,
        nlambda = nlambda,
        lambda.min.ratio = lambda_min_ratio
      )
    }
    s <- choose_lambda_s(cvfit)
    j <- lambda_index(cvfit$lambda, s)
    pred_o <- as.numeric(cvfit$fit.preval[, j])
    pred[obs] <- pred_o
    fold_id[obs] <- as.integer(foldid_o)
    # fill missing-y individuals using full fit
    miss <- which(!is.finite(y))
    if (length(miss) > 0) {
      pred[miss] <- as.numeric(predict(cvfit$glmnet.fit, newx = X[miss, , drop = FALSE], s = s))
    }
    # fold metrics
    fold_metrics <- data.frame(
      fold = integer(0), n_test = integer(0), cor = numeric(0), rmse = numeric(0), mae = numeric(0),
      stringsAsFactors = FALSE
    )
    for (k in seq_len(k_folds)) {
      te_o <- which(foldid_o == k)
      yy <- yo[te_o]
      pp <- pred_o[te_o]
      fold_metrics <- rbind(
        fold_metrics,
        data.frame(
          fold = k,
          n_test = length(te_o),
          cor = pearson(yy, pp),
          rmse = rmse(yy, pp),
          mae = mae(yy, pp),
          stringsAsFactors = FALSE
        )
      )
    }
  } else if (mode == "loo") {
    # LOO CV for observed samples using foldid=1..n_obs.
    foldid_o <- seq_len(length(obs))
    if (lambda_mode == "fixed") {
      cvfit <- glmnet::cv.glmnet(
        x = Xo, y = yo, alpha = alpha,
        foldid = foldid_o,
        keep = TRUE,
        intercept = TRUE, standardize = FALSE,
        lambda = lambda_value
      )
    } else {
      cvfit <- glmnet::cv.glmnet(
        x = Xo, y = yo, alpha = alpha,
        foldid = foldid_o,
        keep = TRUE,
        intercept = TRUE, standardize = FALSE,
        nlambda = nlambda,
        lambda.min.ratio = lambda_min_ratio
      )
    }
    s <- choose_lambda_s(cvfit)
    j <- lambda_index(cvfit$lambda, s)
    pred_o <- as.numeric(cvfit$fit.preval[, j])
    pred[obs] <- pred_o
    fold_id[obs] <- as.integer(foldid_o)
    miss <- which(!is.finite(y))
    if (length(miss) > 0) {
      pred[miss] <- as.numeric(predict(cvfit$glmnet.fit, newx = X[miss, , drop = FALSE], s = s))
    }
  }

  # coefficients (from full fit)
  coef_df <- NULL
  try({
    fit_obj <- if (!is.null(cvfit)) cvfit$glmnet.fit else full_fit
    cc <- glmnet::coef.glmnet(fit_obj, s = s)
    # sparse matrix -> named numeric
    vv <- as.numeric(cc)
    nn <- rownames(cc)
    if (!is.null(nn) && length(nn) == length(vv)) {
      df0 <- data.frame(term = nn, coef = vv, stringsAsFactors = FALSE)
      df0$abs <- abs(df0$coef)
      df0 <- df0[df0$term != "(Intercept)", , drop = FALSE]
      df0 <- df0[order(-df0$abs), , drop = FALSE]
      if (is.finite(top_coef_n) && top_coef_n > 0 && nrow(df0) > top_coef_n) df0 <- df0[seq_len(top_coef_n), , drop = FALSE]
      coef_df <- df0[, c("term", "coef", "abs"), drop = FALSE]
    }
  }, silent = TRUE)

  list(
    pred = pred,
    fold_id = fold_id,
    fold_metrics = fold_metrics,
    cor = pearson(y, pred),
    rmse = rmse(y, pred),
    mae = mae(y, pred),
    n = n,
    n_obs = length(obs),
    lambda = s,
    lambda_mode = lambda_mode,
    coef = coef_df
  )
}

summary_rows <- list()
fold_rows <- list()
artifact_tables <- character(0)
plots <- character(0)

for (ii in seq_along(traits)) {
  tr_name <- traits[ii]
  tr_idx <- trait_idxs[ii]
  cat("[gp_glmnet] trait_idx=", tr_idx, " trait=", tr_name, "\n")

  y_raw <- suppressWarnings(as.numeric(Pt[[tr_name]]))
  y <- y_raw[m[ok]]

  res <- fit_and_predict_trait(y)

  pred_df <- data.frame(
    id = ids,
    trait_idx = tr_idx,
    trait = tr_name,
    y = y,
    yhat = res$pred,
    fold = res$fold_id,
    stringsAsFactors = FALSE
  )

  out_name_idx <- paste0("predictions_", tr_idx, ".tsv")
  fwrite(pred_df, file.path(out_dir, out_name_idx), sep = "\t")
  artifact_tables <- c(artifact_tables, out_name_idx)
  if (length(traits) == 1) {
    fwrite(pred_df, file.path(out_dir, "predictions.tsv"), sep = "\t")
    artifact_tables <- c(artifact_tables, "predictions.tsv")
  }

  # coefficients (top)
  if (!is.null(res$coef) && nrow(res$coef) > 0) {
    coef_name <- paste0("coef_top_", tr_idx, ".tsv")
    fwrite(res$coef, file.path(out_dir, coef_name), sep = "\t")
    artifact_tables <- c(artifact_tables, coef_name)
    if (length(traits) == 1) {
      fwrite(res$coef, file.path(out_dir, "coef_top.tsv"), sep = "\t")
      artifact_tables <- c(artifact_tables, "coef_top.tsv")
    }
  }

  # summary row
  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    trait_idx = tr_idx,
    trait = tr_name,
    n = res$n,
    n_obs = res$n_obs,
    mode = mode,
    backend = "glmnet",
    penalty = penalty,
    alpha = alpha,
    lambda_mode = res$lambda_mode,
    lambda = res$lambda,
    cor = res$cor,
    rmse = res$rmse,
    mae = res$mae,
    stringsAsFactors = FALSE
  )

  if (!is.null(res$fold_metrics)) {
    fm <- res$fold_metrics
    fm$trait_idx <- tr_idx
    fm$trait <- tr_name
    fold_rows[[length(fold_rows) + 1]] <- fm
  }
}

summary_df <- rbindlist(summary_rows, fill = TRUE)
fwrite(summary_df, file.path(out_dir, "summary.tsv"), sep = "\t")
artifact_tables <- unique(c("summary.tsv", artifact_tables))

if (length(fold_rows) > 0) {
  fold_df <- rbindlist(fold_rows, fill = TRUE)
  fwrite(fold_df, file.path(out_dir, "fold_metrics.tsv"), sep = "\t")
  artifact_tables <- unique(c(artifact_tables, "fold_metrics.tsv"))
}

# simple barplot across traits (cor)
try({
  png(file.path(out_dir, "trait_cor.png"), width = 900, height = 450)
  par(mar=c(7,4,2,1))
  barplot(summary_df$cor, names.arg = paste0(summary_df$trait_idx, ":", summary_df$trait), las = 2)
  abline(h=0)
  dev.off()
  plots <- c(plots, "trait_cor.png")
}, silent = TRUE)

mean_or_na <- function(v) {
  v <- as.numeric(v)
  if (!length(v) || all(!is.finite(v))) return(NA_real_)
  mean(v[is.finite(v)])
}

overall <- list(
  mode = mode,
  model = paste0("glmnet_", tolower(penalty)),
  penalty = penalty,
  alpha = alpha,
  lambda_mode = lambda_mode,
  n = length(ids),
  n_markers = ncol(X),
  n_traits = length(traits),
  traits = as.list(traits),
  cor = mean_or_na(summary_df$cor),
  rmse = mean_or_na(summary_df$rmse),
  mae = mean_or_na(summary_df$mae)
)

writeLines(toJSON(overall, auto_unbox = TRUE, pretty = TRUE), file.path(out_dir, "metrics.json"))

artifacts <- list(
  plugin_id = "gp_glmnet",
  tables = artifact_tables,
  plots = unique(plots),
  default_table = if (length(traits) > 1) "summary.tsv" else "predictions.tsv",
  default_plot = if ("trait_cor.png" %in% plots) "trait_cor.png" else "",
  outputs = list(
    metrics_json = "metrics.json"
  )
)
writeLines(toJSON(artifacts, auto_unbox = TRUE, pretty = TRUE), file.path(out_dir, "artifacts.json"))

cat("[gp_glmnet] done\n")

sink()
