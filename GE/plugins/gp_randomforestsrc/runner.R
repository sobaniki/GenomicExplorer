#!/usr/bin/env Rscript

# randomForestSRC runner (multi-trait)
# - Inputs: genotype_tsv (id + markers), phenotype_tsv (id + traits)
# - Optional: covariate_tsv (id + covariates)
# - Trait selection: indices / ranges / names (comma-separated)
# - Modes: fit / loo / kfold
# - Output: predictions_<trait_idx>.tsv (+ predictions.tsv if single),
#           summary.tsv, fold_metrics.tsv, variable_importance.tsv (if enabled),
#           trait_cor.png, pred_vs_true.png, metrics.json, artifacts.json

suppressPackageStartupMessages({
  library(jsonlite)
  library(data.table)
  library(randomForestSRC)
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

cat("[gp_randomforestsrc] start\n")
cat("[gp_randomforestsrc] params_path=", params_path, "\n", sep = "")
cat("[gp_randomforestsrc] out_dir=", out_dir, "\n", sep = "")

p <- fromJSON(params_path, simplifyVector = FALSE)

geno_tsv <- p$genotype_tsv
pheno_tsv <- p$phenotype_tsv

traits <- ""
if (!is.null(p$traits) && nchar(as.character(p$traits)) > 0) traits <- as.character(p$traits)

mode <- if (!is.null(p$mode) && nchar(as.character(p$mode)) > 0) as.character(p$mode) else "fit"
k_folds <- if (!is.null(p$k_folds)) as.integer(p$k_folds) else 5L
seed <- if (!is.null(p$seed)) as.integer(p$seed) else 1L

ntree <- if (!is.null(p$ntree)) as.integer(p$ntree) else 500L
mtry <- if (!is.null(p$mtry) && nchar(as.character(p$mtry)) > 0) suppressWarnings(as.integer(p$mtry)) else NA_integer_
nodesize <- if (!is.null(p$nodesize)) as.integer(p$nodesize) else 5L
importance <- if (!is.null(p$importance) && nchar(as.character(p$importance)) > 0) as.character(p$importance) else "none"
block_size <- if (!is.null(p$block_size) && nchar(as.character(p$block_size)) > 0) suppressWarnings(as.integer(p$block_size)) else NA_integer_

center <- if (!is.null(p$center_markers)) as.logical(p$center_markers) else FALSE
scale <- if (!is.null(p$scale_markers)) as.logical(p$scale_markers) else FALSE
impute <- if (!is.null(p$impute_missing)) as.character(p$impute_missing) else "mean"

if (is.null(geno_tsv) || !file.exists(geno_tsv)) stop("genotype_tsv not found")
if (is.null(pheno_tsv) || !file.exists(pheno_tsv)) stop("phenotype_tsv not found")

cat("[gp_randomforestsrc] genotype_tsv=", geno_tsv, "\n", sep = "")
cat("[gp_randomforestsrc] phenotype_tsv=", pheno_tsv, "\n", sep = "")
cat("[gp_randomforestsrc] traits=", traits, "\n", sep = "")
cat("[gp_randomforestsrc] mode=", mode, " k_folds=", k_folds, " seed=", seed, "\n", sep = "")
cat("[gp_randomforestsrc] ntree=", ntree, " mtry=", mtry, " nodesize=", nodesize, " importance=", importance, " block_size=", block_size, "\n", sep = "")
cat("[gp_randomforestsrc] center=", center, " scale=", scale, " impute=", impute, "\n", sep = "")

if (mode == "kfold" && (is.na(k_folds) || k_folds < 2)) stop("k_folds must be >= 2")

# ------------ Helpers ------------
parse_trait_spec <- function(spec, trait_cols) {
  if (is.null(spec)) spec <- ""
  spec <- trimws(as.character(spec))
  if (nchar(spec) == 0) return(seq_along(trait_cols))
  spec2 <- gsub("[[:space:]]+", "", spec)
  toks <- unlist(strsplit(spec2, ","))
  idxs <- integer(0)
  for (tk in toks) {
    if (!nchar(tk)) next
    if (grepl("^[0-9]+-[0-9]+$", tk)) {
      ab <- suppressWarnings(as.integer(unlist(strsplit(tk, "-"))))
      a <- ab[1]; b <- ab[2]
      if (is.na(a) || is.na(b)) stop(paste0("Invalid trait range: ", tk))
      if (a > b) { tmp <- a; a <- b; b <- tmp }
      idxs <- c(idxs, seq.int(a, b))
    } else if (grepl("^[0-9]+$", tk)) {
      idxs <- c(idxs, as.integer(tk))
    } else {
      j <- match(tk, trait_cols)
      if (is.na(j)) stop(paste0("Trait name not found: ", tk))
      idxs <- c(idxs, j)
    }
  }
  idxs <- idxs[!is.na(idxs)]
  if (!length(idxs)) stop("No valid traits parsed from trait(s) spec.")
  idxs <- idxs[!duplicated(idxs)]
  if (any(idxs < 1 | idxs > length(trait_cols))) stop(paste0("Trait index out of range: 1..", length(trait_cols)))
  idxs
}

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

# Load tables
geno <- data.frame(fread(geno_tsv,
                         header = T),
                   row.names = 1)
if (ncol(geno) < 2) {
  stop("genotype_tsv must have at least 2 columns: <ID> + markers")
} 
ids_geno <- rownames(geno)


pheno <- data.frame(fread(pheno_tsv, 
                          header = T),
                    row.names = 1)
if (ncol(pheno) < 2) {
  stop("phenotype_tsv must have at least 2 columns: <ID> + trait(s)")
}
ids_pheno <- rownames(pheno)
trait_names <- colnames(pheno)

trait_idxs <- parse_trait_spec(traits, trait_names)
trait_target <- trait_names[trait_idxs]

ids_match <- intersect(ids_geno,
                       ids_pheno)
if (length(ids_match) < 1) {
  stop("No overlapping IDs between genotype and phenotype")
}

geno <- geno[ids_match, , drop = F]
pheno <- pheno[ids_match, trait_target, drop = F]

cat("[gp_gp_randomforestsrc] n_samples=", length(ids_match), " n_markers=", ncol(geno), " n_traits=", ncol(pheno), "\n")

for (j in seq_len(ncol(geno))) {
  v <- geno[, j]
  if (anyNA(v)) {
    if (impute == "mean") {
      mu <- mean(v, na.rm = TRUE)
    } else {
      mu <- median(v, na.rm = TRUE)
    }
    if (!is.finite(mu)) mu <- 0
    v[is.na(v)] <- mu
    geno[, j] <- v
  }
}

# Optional center/scale
if (center || scale) {
  geno <- as.matrix(scale(geno, center = center, scale = scale))
}
geno <- t(na.omit(t(geno)))

# ------------ Optional covariates ------------
pred_orig_names <- colnames(geno)
cov_path <- p$covariate_tsv
if (!is.null(cov_path) && nzchar(as.character(cov_path)) && file.exists(as.character(cov_path))) {
  cat("[gp_randomforestsrc] covariate_tsv=", as.character(cov_path), "\n", sep = "")
  C <- fread(as.character(cov_path), sep = "\t", header = TRUE, data.table = FALSE)
  if (ncol(C) < 2) stop("covariate_tsv must have: <ID> + covariates")
  ids_c <- as.character(C[[1]])
  cov_cols <- colnames(C)[-1]
  if (length(cov_cols) > 0) {
    mm <- match(ids, ids_c)
    Z <- as.matrix(C[mm, cov_cols, drop = FALSE])
    Z <- apply(Z, 2, function(z) suppressWarnings(as.numeric(z)))
    Z[is.nan(Z)] <- NA_real_
    # impute covariates
    for (j in seq_len(ncol(Z))) {
      v <- Z[, j]
      if (anyNA(v)) {
        mu <- mean(v, na.rm = TRUE); if (!is.finite(mu)) mu <- 0
        v[is.na(v)] <- mu
        Z[, j] <- v
      }
    }
    cov_names <- paste0("cov_", cov_cols)
    colnames(Z) <- cov_names
    geno <- cbind(Z, geno)
    pred_orig_names <- c(cov_names, pred_orig_names)
  }
}

set.seed(seed)

# ------------ Safe column names for formula ------------
n <- nrow(geno)
pX <- ncol(geno)
kY <- ncol(pheno)

pred_safe <- paste0("x", seq_len(pX))
y_safe <- paste0("y", seq_len(kY))

Xdf <- as.data.frame(geno)
colnames(Xdf) <- pred_safe

Ydf <- as.data.frame(pheno)
colnames(Ydf) <- y_safe

# training requires complete responses
trainable <- complete.cases(Ydf)

DF <- cbind(Ydf, Xdf)

# mappings
pred_map <- data.frame(predictor = pred_safe, marker = pred_orig_names, stringsAsFactors = FALSE)
trait_map <- data.frame(y = y_safe, trait_idx = trait_idxs, trait = traits, stringsAsFactors = FALSE)

# ------------ Model / Predict helpers ------------
make_formula <- function() {
  as.formula(paste0("Multivar(", paste(y_safe, collapse = ","), ") ~ ."))
}

fit_model <- function(train_idx) {
  df_tr <- DF[train_idx, , drop = FALSE]
  args <- list(formula = make_formula(),
               data = df_tr,
               ntree = ntree,
               nodesize = nodesize,
               importance = importance,
               seed = seed)
  if (is.finite(mtry) && !is.na(mtry) && mtry > 0) args$mtry <- mtry
  if (is.finite(block_size) && !is.na(block_size) && block_size > 0) args$block.size <- block_size
  do.call(randomForestSRC::rfsrc, args)
}

# ------------ Run ------------
pred <- matrix(NA, nrow = n, ncol = kY)
colnames(pred) <- trait_target

fold_id <- rep(NA_integer_, n)

if (mode == "fit") {
  tr_idx <- which(trainable)

  fit <- fit_model(tr_idx)
  obj <- predict.rfsrc(fit, DF)
  pred[1:nrow(pred), ] <- get.mv.predicted(obj)
  #saveRDS(pred, file = "/media/soba/Noc4/GenomicExplorer/P0213/sim/ex2.rds")
  
  fold_id[] <- 0L
} else if (mode == "kfold") {
  perm <- sample.int(n)
  folds <- split(perm, rep(seq_len(k_folds), length.out = n))
  for (k in seq_len(k_folds)) {
    test_idx <- folds[[k]]
    train_idx <- setdiff(seq_len(n), test_idx)
    train_idx <- train_idx[trainable[train_idx]]
    fit <- fit_model(train_idx)
    pred[test_idx, ] <- predict_model(fit, test_idx)
    fold_id[test_idx] <- k
    cat("  fold ", k, "/", k_folds, "\n", sep = "")
  }
} else if (mode == "loo") {
  for (i in seq_len(n)) {
    test_idx <- i
    train_idx <- setdiff(seq_len(n), test_idx)
    train_idx <- train_idx[trainable[train_idx]]
    fit <- fit_model(train_idx)
    pred[i, ] <- as.numeric(predict_model(fit, test_idx))
    fold_id[i] <- i
    if (i %% 10 == 0) cat("  loo ", i, "/", n, "\n", sep = "")
  }
}

# ------------ Write outputs ------------
summary_rows <- list()
fold_rows <- list()
pred_files <- character(0)
wrote_scatter <- FALSE

for (ii in seq_along(traits)) {
  tr <- trait_target[ii]
  tr_idx <- trait_idxs[ii]
  y_true <- as.numeric(pheno[, ii])
  y_pred <- as.numeric(pred[, ii])

  DT <- data.frame(
    id = ids_match,
    trait_idx = tr_idx,
    trait = tr,
    y_true = y_true,
    y_pred = y_pred,
    fold = fold_id,
    stringsAsFactors = FALSE
  )

  out_tsv <- paste0("predictions_", tr_idx, ".tsv")
  fwrite(DT, file.path(out_dir, out_tsv), sep = "\t")
  pred_files <- c(pred_files, out_tsv)
  if (length(traits) == 1) {
    fwrite(DT, file.path(out_dir, "predictions.tsv"), sep = "\t")
    pred_files <- c(pred_files, "predictions.tsv")
  }

  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    trait_idx = tr_idx,
    trait = tr,
    n = nrow(DT),
    n_obs = sum(is.finite(y_true)),
    mode = mode,
    backend = "randomForestSRC",
    cor = pearson(y_true, y_pred),
    rmse = rmse(y_true, y_pred),
    mae = mae(y_true, y_pred),
    stringsAsFactors = FALSE
  )

  folds_u <- sort(unique(fold_id[!is.na(fold_id)]))
  for (fk in folds_u) {
    idx <- which(fold_id == fk)
    if (!length(idx)) next
    yt <- y_true[idx]; yp <- y_pred[idx]
    fold_rows[[length(fold_rows) + 1]] <- data.frame(
      trait_idx = tr_idx,
      trait = tr,
      fold = fk,
      n_test = length(idx),
      cor = pearson(yt, yp),
      rmse = rmse(yt, yp),
      mae = mae(yt, yp),
      stringsAsFactors = FALSE
    )
  }

  if (!wrote_scatter) {
    try({
      png(file.path(out_dir, "pred_vs_true.png"), width = 900, height = 700)
      plot(DT$y_true, DT$y_pred, xlab = "Observed", ylab = "Predicted",
           main = paste0("randomForestSRC GP (trait ", tr_idx, ": ", tr, ")"))
      abline(0, 1)
      dev.off()
    }, silent = TRUE)
    wrote_scatter <- TRUE
  }
}

summary_df <- rbindlist(summary_rows, fill = TRUE)
fwrite(summary_df, file.path(out_dir, "summary.tsv"), sep = "\t")

fold_df <- rbindlist(fold_rows, fill = TRUE)
if (nrow(fold_df) == 0) {
  fold_df <- data.frame(trait_idx = integer(0), trait = character(0), fold = integer(0), n_test = integer(0),
                        cor = numeric(0), rmse = numeric(0), mae = numeric(0))
}
fwrite(fold_df, file.path(out_dir, "fold_metrics.tsv"), sep = "\t")

# ------------ Variable importance (fit on all complete cases) ------------
if (!is.null(importance) && importance != "none") {
  tr_idx_cc <- which(trainable)
  if (length(tr_idx_cc) >= 2) {
    cat("[gp_randomforestsrc] computing variable importance via vimp() on complete cases\n")
    fit_all <- fit_model(tr_idx_cc)
    vi_obj <- try(randomForestSRC::vimp(fit_all), silent = TRUE)
    if (!inherits(vi_obj, "try-error") && !is.null(vi_obj$importance)) {
      vi <- vi_obj$importance
      out_rows <- list()

      if (is.matrix(vi)) {
        rn <- rownames(vi); if (is.null(rn)) rn <- pred_safe
        cn <- colnames(vi)
        if (is.null(cn)) cn <- y_safe
        for (jj in seq_along(cn)) {
          yj <- cn[jj]
          tm <- trait_map[trait_map$y == yj, , drop = FALSE]
          if (nrow(tm) == 0) next
          vj <- as.numeric(vi[, jj])
          dfj <- data.frame(
            trait_idx = tm$trait_idx,
            trait = tm$trait,
            predictor = rn,
            importance = vj,
            stringsAsFactors = FALSE
          )
          out_rows[[length(out_rows) + 1]] <- dfj
        }
        vi_df <- rbindlist(out_rows, fill = TRUE)
      } else {
        nm <- names(vi); if (is.null(nm)) nm <- pred_safe
        vi_df <- data.frame(
          trait_idx = 0,
          trait = "ALL",
          predictor = nm,
          importance = as.numeric(vi),
          stringsAsFactors = FALSE
        )
      }

      vi_df <- merge(vi_df, pred_map, by.x = "predictor", by.y = "predictor", all.x = TRUE, sort = FALSE)
      if (!("marker" %in% colnames(vi_df))) vi_df$marker <- vi_df$predictor
      vi_df <- vi_df[, c("trait_idx", "trait", "marker", "importance")]
      vi_df <- vi_df[order(vi_df$trait_idx, -vi_df$importance), , drop = FALSE]
      fwrite(vi_df, file.path(out_dir, "variable_importance.tsv"), sep = "\t")
    } else {
      cat("[gp_randomforestsrc] WARN: vimp() did not return importance.\n")
    }
  } else {
    cat("[gp_randomforestsrc] WARN: not enough complete cases for vimp().\n")
  }
}

# ------------ Plot: trait-wise cor ------------
try({
  png(file.path(out_dir, "trait_cor.png"), width = 900, height = 450)
  par(mar = c(7, 4, 2, 1))
  barplot(summary_df$cor, names.arg = summary_df$trait, las = 2)
  abline(h = 0)
  dev.off()
}, silent = TRUE)

mean_or_na <- function(v) {
  v <- as.numeric(v)
  if (!length(v) || all(!is.finite(v))) return(NA_real_)
  mean(v[is.finite(v)])
}

metrics <- list(
  mode = mode,
  model_type = "randomForestSRC",
  n = n,
  p = pX,
  n_traits = length(traits),
  traits = traits,
  cor = mean_or_na(summary_df$cor),
  rmse = mean_or_na(summary_df$rmse),
  mae = mean_or_na(summary_df$mae),
  k_folds = if (mode == "kfold") k_folds else NA,
  ntree = ntree,
  mtry = if (is.na(mtry)) NA else mtry,
  nodesize = nodesize,
  importance = importance,
  block_size = if (is.na(block_size)) NA else block_size,
  n_complete_cases = sum(trainable)
)
writeLines(toJSON(metrics, auto_unbox = TRUE, pretty = TRUE), file.path(out_dir, "metrics.json"))

# Artifacts
tables_list <- list(list(name = "summary.tsv", path = "summary.tsv"))
for (nm in pred_files) {
  tables_list[[length(tables_list) + 1]] <- list(name = nm, path = nm)
}
tables_list[[length(tables_list) + 1]] <- list(name = "fold_metrics.tsv", path = "fold_metrics.tsv")
if (file.exists(file.path(out_dir, "variable_importance.tsv"))) {
  tables_list[[length(tables_list) + 1]] <- list(name = "variable_importance.tsv", path = "variable_importance.tsv")
}

plots_list <- list(list(name = "trait_cor.png", path = "trait_cor.png"))
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

cat("[gp_randomforestsrc] done\n")
sink()
