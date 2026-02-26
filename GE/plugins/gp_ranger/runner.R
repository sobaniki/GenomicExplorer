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

cat("[gp_ranger] start\n")
cat("[gp_ranger] params_path=", params_path, "\n")
cat("[gp_ranger] out_dir=", out_dir, "\n")

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

mode <- if (!is.null(p$mode) && nchar(p$mode) > 0) p$mode else "fit"  # fit / loo / kfold
k_folds <- if (!is.null(p$k_folds)) as.integer(p$k_folds) else 5
seed <- if (!is.null(p$seed)) as.integer(p$seed) else 1

num_trees <- if (!is.null(p$num_trees)) as.integer(p$num_trees) else 500
mtry <- if (!is.null(p$mtry) && nchar(as.character(p$mtry)) > 0) as.integer(p$mtry) else NA_integer_
min_node_size <- if (!is.null(p$min_node_size)) as.integer(p$min_node_size) else 5
importance <- if (!is.null(p$importance) && nchar(p$importance) > 0) p$importance else "none"  # none / impurity / permutation

center <- if (!is.null(p$center_markers)) as.logical(p$center_markers) else FALSE
scale <- if (!is.null(p$scale_markers)) as.logical(p$scale_markers) else FALSE
impute <- if (!is.null(p$impute_missing)) p$impute_missing else "mean"

cat("[gp_ranger] genotype_tsv=", geno_tsv, "\n")
cat("[gp_ranger] phenotype_tsv=", pheno_tsv, "\n")
cat("[gp_ranger] traits_spec=", traits_spec, "
")
cat("[gp_ranger] mode=", mode, " k_folds=", k_folds, " seed=", seed, "\n")
cat("[gp_ranger] num_trees=", num_trees, " mtry=", mtry, " min_node_size=", min_node_size, " importance=", importance, "\n")

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

set.seed(seed)

# Load genotype
Gt <- fread(geno_tsv, sep = "	", header = TRUE, data.table = FALSE)
if (ncol(Gt) < 2) stop("genotype_tsv must have at least 2 columns: <ID> + markers")
id_col <- colnames(Gt)[1]
ids_g <- as.character(Gt[[1]])
marker_cols <- colnames(Gt)[-1]
if (length(marker_cols) == 0) stop("No marker columns in genotype_tsv")
X <- as.matrix(Gt[, marker_cols, drop = FALSE])
X <- apply(X, 2, function(z) suppressWarnings(as.numeric(z)))
X[is.nan(X)] <- NA

# Load phenotype
Pt <- fread(pheno_tsv, sep = "	", header = TRUE, data.table = FALSE)
if (ncol(Pt) < 2) stop("phenotype_tsv must have at least 2 columns: <ID> + trait(s)")
id_col_p <- colnames(Pt)[1]
ids_p <- as.character(Pt[[1]])
trait_cols <- colnames(Pt)[-1]

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

trait_idxs <- parse_trait_spec(traits_spec, trait_cols)
traits <- trait_cols[trait_idxs]

m <- match(ids_g, ids_p)
ok <- !is.na(m)
if (!any(ok)) stop("No overlapping IDs between genotype and phenotype")
X <- X[ok, , drop = FALSE]
ids <- ids_g[ok]

# Impute missing markers (mean)
if (impute == "mean") {
  for (j in seq_len(ncol(X))) {
    v <- X[, j]
    if (anyNA(v)) {
      mu <- mean(v, na.rm = TRUE)
      if (!is.finite(mu)) mu <- 0
      v[is.na(v)] <- mu
      X[, j] <- v
    }
  }
}

# Optional center/scale
if (center || scale) {
  X <- as.matrix(scale(X, center = center, scale = scale))
}

# Optional covariates (id + covariate columns)
cov_path <- p$covariate_tsv
if (!is.null(cov_path) && nzchar(cov_path) && file.exists(cov_path)) {
  cat(sprintf("[gp_ranger] covariate_tsv=%s\n", cov_path))
  C <- fread(cov_path, sep = "\t", header = TRUE, data.table = FALSE)
  if (ncol(C) < 2) stop("covariate_tsv must have at least 2 columns: <ID> + covariates")
  id_col_c <- colnames(C)[1]
  cov_cols <- colnames(C)[-1]
  if (length(cov_cols) > 0) {
    ids_c <- as.character(C[[1]])
    mm <- match(ids, ids_c)
    Z <- as.matrix(C[mm, cov_cols, drop = FALSE])
    Z <- apply(Z, 2, function(z) suppressWarnings(as.numeric(z)))
    Z[is.nan(Z)] <- NA
    # impute mean for covariates
    for (j in seq_len(ncol(Z))) {
      v <- Z[, j]
      if (anyNA(v)) {
        mu <- mean(v, na.rm = TRUE)
        if (!is.finite(mu)) mu <- 0
        v[is.na(v)] <- mu
        Z[, j] <- v
      }
    }
    colnames(Z) <- paste0("cov_", cov_cols)
    X <- cbind(Z, X)
  }
}


# -----------------------
# Run per-trait
# -----------------------

n_total <- length(ids)

use_ranger <- requireNamespace("ranger", quietly = TRUE)
if (!use_ranger) {
  cat("[gp_ranger] WARN: ranger package not available. Using mean predictor as placeholder.\n")
}

fit_predict <- function(y, train_idx, test_idx) {
  if (!use_ranger) {
    mu <- mean(y[train_idx], na.rm = TRUE)
    if (!is.finite(mu)) mu <- 0
    return(rep(mu, length(test_idx)))
  }
  # drop missing y in training
  tr_ok <- train_idx[is.finite(y[train_idx])]
  if (length(tr_ok) < 2) {
    mu <- mean(y[train_idx], na.rm = TRUE)
    if (!is.finite(mu)) mu <- 0
    return(rep(mu, length(test_idx)))
  }
  df_tr <- data.frame(y = y[tr_ok], X[tr_ok, , drop = FALSE])

  # auto mtry if not provided
  mtry_use <- mtry
  if (is.na(mtry_use) || !is.finite(mtry_use) || mtry_use <= 0) {
    mtry_use <- max(1L, floor(sqrt(ncol(X))))
  }
  rf <- ranger::ranger(
    dependent.variable.name = "y",
    data = df_tr,
    num.trees = num_trees,
    mtry = mtry_use,
    min.node.size = min_node_size,
    importance = importance,
    seed = seed,
    respect.unordered.factors = TRUE,
    write.forest = TRUE
  )
  df_te <- data.frame(y = rep(NA_real_, length(test_idx)), X[test_idx, , drop = FALSE])
  pr <- predict(rf, data = df_te)$predictions
  as.numeric(pr)
}

run_one_trait <- function(y) {
  n <- length(y)
  pred <- rep(NA_real_, n)
  fold_id <- rep(0L, n)

  obs <- which(is.finite(y))
  if (!(mode %in% c("fit", "loo", "kfold"))) stop("mode must be one of: fit, loo, kfold")

  fold_metrics <- data.frame(
    fold = integer(0),
    n_test = integer(0),
    cor = numeric(0),
    rmse = numeric(0),
    mae = numeric(0),
    stringsAsFactors = FALSE
  )

  if (mode == "fit") {
    # train on observed, predict all
    if (length(obs) >= 2) {
      pred <- fit_predict(y, obs, seq_len(n))
    } else {
      mu <- mean(y, na.rm = TRUE); if (!is.finite(mu)) mu <- 0
      pred <- rep(mu, n)
    }
    fold_id <- rep(0L, n)
  } else if (mode == "loo") {
    cat("[gp_ranger] running LOO-CV ...\n")
    for (ii in seq_along(obs)) {
      i <- obs[ii]
      tr <- setdiff(obs, i)
      pred[i] <- fit_predict(y, tr, i)
      fold_id[i] <- i
      if (ii %% 25 == 0) cat("  ...", ii, "/", length(obs), "\n")
    }
    # predict missing-y rows (if any) from model trained on all observed
    miss <- setdiff(seq_len(n), obs)
    if (length(miss) > 0 && length(obs) >= 2) {
      pred[miss] <- fit_predict(y, obs, miss)
      fold_id[miss] <- 0L
    }
  } else if (mode == "kfold") {
    if (is.na(k_folds) || k_folds < 2) stop("k_folds must be >= 2")
    cat("[gp_ranger] running k-fold CV ...\n")
    perm <- sample(obs)
    folds <- split(perm, rep(seq_len(k_folds), length.out = length(perm)))
    for (k in seq_len(k_folds)) {
      te <- sort(folds[[k]])
      tr <- setdiff(obs, te)
      pred[te] <- fit_predict(y, tr, te)
      fold_id[te] <- k
      fold_metrics <- rbind(
        fold_metrics,
        data.frame(
          fold = k,
          n_test = length(te),
          cor = pearson(y[te], pred[te]),
          rmse = rmse(y[te], pred[te]),
          mae = mae(y[te], pred[te]),
          stringsAsFactors = FALSE
        )
      )
    }
    # predict missing-y rows (if any) from model trained on all observed
    miss <- setdiff(seq_len(n), obs)
    if (length(miss) > 0 && length(obs) >= 2) {
      pred[miss] <- fit_predict(y, obs, miss)
      fold_id[miss] <- 0L
    }
  }

  list(
    pred = pred,
    fold_id = fold_id,
    fold_metrics = fold_metrics,
    cor = pearson(y, pred),
    rmse = rmse(y, pred),
    mae = mae(y, pred),
    n = n,
    n_obs = length(obs)
  )
}

summary_rows <- list()
fold_rows <- list()
imp_rows <- list()
pred_files <- character(0)
wrote_scatter <- FALSE

for (ii in seq_along(traits)) {
  tr_name <- traits[ii]
  tr_idx <- trait_idxs[ii]
  cat("[gp_ranger] trait_idx=", tr_idx, " trait=", tr_name, "\n")

  y_raw <- suppressWarnings(as.numeric(Pt[[tr_name]]))
  y <- y_raw[m[ok]]

  res <- run_one_trait(y)

  pred_df <- data.frame(
    id = ids,
    trait_idx = tr_idx,
    trait = tr_name,
    y_true = y,
    y_pred = res$pred,
    fold = res$fold_id,
    stringsAsFactors = FALSE
  )

  out_name_idx <- paste0("predictions_", tr_idx, ".tsv")
  fwrite(pred_df, file.path(out_dir, out_name_idx), sep = "	")
  pred_files <- c(pred_files, out_name_idx)

  if (length(traits) == 1) {
    fwrite(pred_df, file.path(out_dir, "predictions.tsv"), sep = "	")
    pred_files <- c(pred_files, "predictions.tsv")
  }

  # fold metrics
  if (nrow(res$fold_metrics) > 0) {
    fm <- res$fold_metrics
    fm$trait_idx <- tr_idx
    fm$trait <- tr_name
    fold_rows[[length(fold_rows) + 1]] <- fm
  }

  # summary row
  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    trait_idx = tr_idx,
    trait = tr_name,
    n = res$n,
    n_obs = res$n_obs,
    mode = mode,
    backend = "ranger",
    cor = res$cor,
    rmse = res$rmse,
    mae = res$mae,
    stringsAsFactors = FALSE
  )

  # variable importance per trait (fit on all observed)
  if (use_ranger && importance != "none") {
    obs <- which(is.finite(y))
    if (length(obs) >= 2) {
      df_all <- data.frame(y = y[obs], X[obs, , drop = FALSE])
      mtry_use <- mtry
      if (is.na(mtry_use) || !is.finite(mtry_use) || mtry_use <= 0) {
        mtry_use <- max(1L, floor(sqrt(ncol(X))))
      }
      rf_all <- ranger::ranger(
        dependent.variable.name = "y",
        data = df_all,
        num.trees = num_trees,
        mtry = mtry_use,
        min.node.size = min_node_size,
        importance = importance,
        seed = seed
      )
      imp <- ranger::importance(rf_all)
      imp_df <- data.frame(
        trait_idx = tr_idx,
        trait = tr_name,
        marker = names(imp),
        importance = as.numeric(imp),
        stringsAsFactors = FALSE
      )
      imp_rows[[length(imp_rows) + 1]] <- imp_df
    }
  }

  # scatter plot (first trait only)
  if (!wrote_scatter) {
    png(file.path(out_dir, "pred_vs_true.png"), width = 900, height = 700)
    plot(pred_df$y_true, pred_df$y_pred, xlab = "Observed", ylab = "Predicted", main = paste0("ranger GP (trait ", tr_idx, ")"))
    abline(0, 1)
    dev.off()
    wrote_scatter <- TRUE
  }
}

summary_df <- rbindlist(summary_rows, fill = TRUE)
fwrite(summary_df, file.path(out_dir, "summary.tsv"), sep = "\t")

if (length(fold_rows) > 0) {
  fold_df <- rbindlist(fold_rows, fill = TRUE)
  fwrite(fold_df, file.path(out_dir, "fold_metrics.tsv"), sep = "\t")
}

if (length(imp_rows) > 0) {
  imp_all <- rbindlist(imp_rows, fill = TRUE)
  imp_all <- imp_all[order(imp_all$trait_idx, -imp_all$importance), , drop = FALSE]
  fwrite(imp_all, file.path(out_dir, "variable_importance.tsv"), sep = "\t")
}

# trait correlation barplot
try({
  png(file.path(out_dir, "trait_cor.png"), width = 900, height = 450)
  par(mar=c(7,4,2,1))
  barplot(summary_df$cor, names.arg = paste0(summary_df$trait_idx, ":", summary_df$trait), las = 2)
  abline(h=0)
  dev.off()
}, silent = TRUE)

mean_or_na <- function(v) {
  v <- as.numeric(v)
  if (!length(v) || all(!is.finite(v))) return(NA_real_)
  mean(v[is.finite(v)])
}

m_all <- list(
  mode = mode,
  model_type = "ranger",
  n = n_total,
  p = ncol(X),
  n_traits = length(traits),
  traits = as.list(traits),
  cor = mean_or_na(summary_df$cor),
  rmse = mean_or_na(summary_df$rmse),
  mae = mean_or_na(summary_df$mae),
  k_folds = if (mode == "kfold") k_folds else NA,
  num_trees = num_trees,
  mtry = if (is.na(mtry)) NA else mtry,
  min_node_size = min_node_size,
  importance = importance,
  used_ranger = use_ranger
)
writeLines(toJSON(m_all, auto_unbox = TRUE, pretty = TRUE), file.path(out_dir, "metrics.json"))

# Artifacts
tables_list <- list(
  list(name = "summary.tsv", path = "summary.tsv")
)
# predictions files
for (nm in pred_files) {
  tables_list[[length(tables_list) + 1]] <- list(name = nm, path = nm)
}
if (file.exists(file.path(out_dir, "fold_metrics.tsv"))) {
  tables_list[[length(tables_list) + 1]] <- list(name = "fold_metrics.tsv", path = "fold_metrics.tsv")
}
if (file.exists(file.path(out_dir, "variable_importance.tsv"))) {
  tables_list[[length(tables_list) + 1]] <- list(name = "variable_importance.tsv", path = "variable_importance.tsv")
}

plots_list <- list(
  list(name = "trait_cor.png", path = "trait_cor.png")
)
if (file.exists(file.path(out_dir, "pred_vs_true.png"))) {
  plots_list[[length(plots_list) + 1]] <- list(name = "pred_vs_true.png", path = "pred_vs_true.png")
}

art <- list(
  tables = tables_list,
  plots = plots_list,
  default_table = if (length(traits) > 1) "summary.tsv" else "predictions.tsv",
  default_plot = "trait_cor.png"
)
writeLines(toJSON(art, auto_unbox = TRUE, pretty = TRUE), file.path(out_dir, "artifacts.json"))

cat("[gp_ranger] done\n")
sink()
