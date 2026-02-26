#!/usr/bin/env Rscript

# BGLR::Multitrait runner (multi-trait genomic prediction)
# - Inputs: genotype_tsv (id + markers), phenotype_tsv (id + traits)
# - traits spec: indices (1-based excluding ID), ranges, comma list, or trait names
# - Modes: fit / loo / kfold
# - ETA: multi-kernel (BRR / SpikeSlab / RKHS) + optional JSON extra terms, matching single-trait BGLR UI
# - resCov: user-provided R code (evaluated) passed to BGLR::Multitrait(resCov=...)
# - Output: predictions_<trait_idx>.tsv (+ predictions.tsv when single trait), summary.tsv,
#           fold_metrics.tsv, metrics.json, trait_cor.png, artifacts.json

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

cat("[gp_bglr_multitrait] start\n")
cat("[gp_bglr_multitrait] params_path=", params_path, "\n", sep="")
cat("[gp_bglr_multitrait] out_dir=", out_dir, "\n", sep="")

p <- fromJSON(params_path, simplifyVector = FALSE)

# ------------ Inputs ------------
geno_tsv <- p$genotype_tsv
pheno_tsv <- p$phenotype_tsv
cov_tsv  <- if (!is.null(p$covariate_tsv) && nchar(as.character(p$covariate_tsv)) > 0) as.character(p$covariate_tsv) else ""

traits_spec <- ""
if (!is.null(p$traits_spec) && nchar(as.character(p$traits_spec)) > 0) traits_spec <- as.character(p$traits_spec)
if (!is.null(p$traits) && nchar(as.character(p$traits)) > 0) traits_spec <- as.character(p$traits)
if (!is.null(p$trait) && nchar(as.character(p$trait)) > 0) traits_spec <- as.character(p$trait)

mode <- if (!is.null(p$mode) && nchar(as.character(p$mode)) > 0) as.character(p$mode) else "fit"
k_folds <- if (!is.null(p$k_folds)) as.integer(p$k_folds) else 5L
seed <- if (!is.null(p$seed)) as.integer(p$seed) else 1L

nIter <- if (!is.null(p$nIter)) as.integer(p$nIter) else 2500L
burnIn <- if (!is.null(p$burnIn)) as.integer(p$burnIn) else 800L
thin <- if (!is.null(p$thin)) as.integer(p$thin) else 5L

center <- if (!is.null(p$center_markers)) as.logical(p$center_markers) else TRUE
scale <- if (!is.null(p$scale_markers)) as.logical(p$scale_markers) else TRUE
impute <- if (!is.null(p$impute_missing)) as.character(p$impute_missing) else "mean"

# Kernel/ETA params (same names as single-trait gp_bglr_brr)
get_chr <- function(x, default="") {
  if (is.null(x)) return(default)
  z <- as.character(x)
  if (is.na(z)) return(default)
  z
}

kernel1_model <- get_chr(p$kernel1_model, "BRR")
kernel1_probIn <- get_chr(p$kernel1_probIn, "")
kernel1_rkhs_source <- get_chr(p$kernel1_rkhs_source, "markers")
kernel1_rkhs_kernel <- get_chr(p$kernel1_rkhs_kernel, "linear")
kernel1_rkhs_h <- get_chr(p$kernel1_rkhs_h, "")
kernel1_rkhs_matrix_kind <- get_chr(p$kernel1_rkhs_matrix_kind, "kernel")
kernel1_rkhs_matrix_tsv <- get_chr(p$kernel1_rkhs_matrix_tsv, "")

kernel2_model <- get_chr(p$kernel2_model, "None")
kernel2_probIn <- get_chr(p$kernel2_probIn, "")
kernel2_rkhs_source <- get_chr(p$kernel2_rkhs_source, "markers")
kernel2_rkhs_kernel <- get_chr(p$kernel2_rkhs_kernel, "linear")
kernel2_rkhs_h <- get_chr(p$kernel2_rkhs_h, "")
kernel2_rkhs_matrix_kind <- get_chr(p$kernel2_rkhs_matrix_kind, "kernel")
kernel2_rkhs_matrix_tsv <- get_chr(p$kernel2_rkhs_matrix_tsv, "")

eta_json <- get_chr(p$eta_json, "")
resCov_code <- get_chr(p$resCov, "")

if (is.null(geno_tsv) || !file.exists(geno_tsv)) stop("genotype_tsv not found")
if (is.null(pheno_tsv) || !file.exists(pheno_tsv)) stop("phenotype_tsv not found")
if (nchar(cov_tsv) > 0 && !file.exists(cov_tsv)) stop("covariate_tsv not found")

cat("[gp_bglr_multitrait] genotype_tsv=", geno_tsv, "\n", sep="")
cat("[gp_bglr_multitrait] phenotype_tsv=", pheno_tsv, "\n", sep="")
cat("[gp_bglr_multitrait] covariate_tsv=", cov_tsv, "\n", sep="")
cat("[gp_bglr_multitrait] traits_spec=", traits_spec, "\n", sep="")
cat("[gp_bglr_multitrait] mode=", mode, " k_folds=", k_folds, " seed=", seed, "\n", sep="")
cat("[gp_bglr_multitrait] nIter=", nIter, " burnIn=", burnIn, " thin=", thin, "\n", sep="")
cat("[gp_bglr_multitrait] center=", center, " scale=", scale, " impute=", impute, "\n", sep="")
cat("[gp_bglr_multitrait] kernel1_model=", kernel1_model, " kernel2_model=", kernel2_model, "\n", sep="")

if (!(mode %in% c("fit", "loo", "kfold"))) stop("mode must be one of: fit, loo, kfold")
if (mode == "kfold" && (is.na(k_folds) || k_folds < 2)) stop("k_folds must be >= 2")

# ------------ Helpers ------------
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

as_num <- function(x, default = NA_real_) {
  if (is.null(x)) return(default)
  v <- suppressWarnings(as.numeric(as.character(x)))
  if (!is.finite(v)) default else v
}

normalize_model <- function(m) {
  m <- trimws(as.character(m))
  if (!nchar(m)) return("None")
  # Keep case as BGLR expects
  m
}

read_square_matrix <- function(path, ids) {
  M <- as.matrix(fread(path, sep = "\t", header = TRUE, data.table = FALSE, check.names = FALSE))
  if (nrow(M) == 0 || ncol(M) == 0) stop(paste0("Empty matrix: ", path))

  # If first column is IDs, use it as rownames
  if (!is.null(colnames(M)) && (colnames(M)[1] %in% c("id", "ID", "sample", "Sample"))) {
    rn <- as.character(M[, 1])
    M <- M[, -1, drop = FALSE]
    rownames(M) <- rn
  }

  # Ensure square with row/col names
  if (is.null(rownames(M)) || is.null(colnames(M))) {
    stop("Matrix TSV must have rownames and colnames (or first column as id).")
  }

  common <- intersect(ids, rownames(M))
  common <- intersect(common, colnames(M))
  if (length(common) < length(ids)) {
    missing <- setdiff(ids, common)
    stop(paste0("Matrix is missing IDs: ", paste(head(missing, 10), collapse=","), ifelse(length(missing) > 10, "...", "")))
  }
  M2 <- M[ids, ids, drop = FALSE]
  as.matrix(M2)
}

build_kernel_from_markers <- function(X, kernel = "linear", h = NA_real_) {
  kernel <- tolower(trimws(as.character(kernel)))
  if (kernel == "linear") {
    K <- tcrossprod(X) / max(1, ncol(X))
    return(K)
  }
  # gaussian: exp(-D^2/h)
  G <- tcrossprod(X)
  d2 <- outer(diag(G), diag(G), "+") - 2 * G
  d2[d2 < 0] <- 0
  if (!is.finite(h) || h <= 0) {
    # heuristic: median of distances
    vals <- as.numeric(d2[upper.tri(d2)])
    h <- median(vals[is.finite(vals)], na.rm = TRUE)
    if (!is.finite(h) || h <= 0) h <- 1
  }
  exp(-d2 / h)
}

# ------------ Load tables ------------
GT <- fread(geno_tsv, sep = "\t", header = TRUE, data.table = FALSE)
if (ncol(GT) < 2) stop("genotype_tsv must have: <ID> + markers")
ids_g <- as.character(GT[[1]])
marker_cols <- colnames(GT)[-1]
X <- as.matrix(GT[, marker_cols, drop = FALSE])
X <- apply(X, 2, function(z) suppressWarnings(as.numeric(z)))
X[is.nan(X)] <- NA_real_

PT <- fread(pheno_tsv, sep = "\t", header = TRUE, data.table = FALSE)
if (ncol(PT) < 2) stop("phenotype_tsv must have: <ID> + trait(s)")
ids_p <- as.character(PT[[1]])
trait_cols <- colnames(PT)[-1]

trait_idxs <- parse_trait_spec(traits_spec, trait_cols)
traits <- trait_cols[trait_idxs]

m <- match(ids_g, ids_p)
ok <- !is.na(m)
if (!any(ok)) stop("No overlapping IDs between genotype and phenotype")

ids <- ids_g[ok]
X <- X[ok, , drop = FALSE]
Y <- sapply(traits, function(tr) suppressWarnings(as.numeric(PT[[tr]][m[ok]])))
Y <- as.matrix(Y)
colnames(Y) <- traits

cat("[gp_bglr_multitrait] n_samples=", nrow(X), " n_markers=", ncol(X), " n_traits=", ncol(Y), "\n", sep="")

# ------------ Preprocess markers ------------
if (impute %in% c("mean", "median")) {
  for (j in seq_len(ncol(X))) {
    v <- X[, j]
    if (anyNA(v)) {
      mu <- if (impute == "median") median(v, na.rm = TRUE) else mean(v, na.rm = TRUE)
      if (!is.finite(mu)) mu <- 0
      v[is.na(v)] <- mu
      X[, j] <- v
    }
  }
}
if (center || scale) {
  X <- as.matrix(scale(X, center = center, scale = scale))
}
X <- t(na.omit(t(X)))

# ------------ Optional covariates ------------
Z <- NULL
if (nchar(cov_tsv) > 0) {
  CT <- fread(cov_tsv, sep = "\t", header = TRUE, data.table = FALSE)
  if (ncol(CT) < 2) stop("covariate_tsv must have: <ID> + covariates")
  ids_c <- as.character(CT[[1]])
  m2 <- match(ids, ids_c)
  if (any(is.na(m2))) {
    miss <- ids[is.na(m2)]
    stop(paste0("covariate_tsv missing IDs: ", paste(head(miss, 10), collapse=","), ifelse(length(miss) > 10, "...", "")))
  }
  Z <- as.matrix(CT[m2, -1, drop = FALSE])
  Z <- apply(Z, 2, function(z) suppressWarnings(as.numeric(z)))
  Z[is.nan(Z)] <- NA_real_
  # simple impute
  for (j in seq_len(ncol(Z))) {
    v <- Z[, j]
    if (anyNA(v)) {
      mu <- mean(v, na.rm = TRUE)
      if (!is.finite(mu)) mu <- 0
      v[is.na(v)] <- mu
      Z[, j] <- v
    }
  }
  cat("[gp_bglr_multitrait] covariates n=", ncol(Z), "\n", sep="")
}

set.seed(seed)

# ------------ Build ETA list ------------
SUPPORTED <- c("BRR", "SpikeSlab", "RKHS", "FIXED")

term_from_prefix <- function(prefix) {
  model <- normalize_model(get_chr(p[[paste0(prefix, "_model")]], ifelse(prefix=="kernel1", "BRR", "None")))
  probIn <- get_chr(p[[paste0(prefix, "_probIn")]], "")
  rkhs_source <- normalize_model(get_chr(p[[paste0(prefix, "_rkhs_source")]], "markers"))
  rkhs_kernel <- normalize_model(get_chr(p[[paste0(prefix, "_rkhs_kernel")]], "linear"))
  rkhs_h <- get_chr(p[[paste0(prefix, "_rkhs_h")]], "")
  rkhs_matrix_kind <- normalize_model(get_chr(p[[paste0(prefix, "_rkhs_matrix_kind")]], "kernel"))
  rkhs_matrix_tsv <- get_chr(p[[paste0(prefix, "_rkhs_matrix_tsv")]], "")

  list(
    model = model,
    probIn = probIn,
    rkhs_source = rkhs_source,
    rkhs_kernel = rkhs_kernel,
    rkhs_h = rkhs_h,
    rkhs_matrix_kind = rkhs_matrix_kind,
    rkhs_matrix_tsv = rkhs_matrix_tsv
  )
}

build_term <- function(t) {
  m <- t$model
  if (m == "None") return(NULL)
  if (!(m %in% SUPPORTED)) stop(paste0("Unsupported model: ", m, " (supported: ", paste(SUPPORTED, collapse=", "), ")"))

  if (m == "RKHS") {
    source <- tolower(trimws(as.character(t$rkhs_source)))
    kern <- tolower(trimws(as.character(t$rkhs_kernel)))
    h <- as_num(t$rkhs_h, NA_real_)
    if (source == "matrix") {
      if (!nchar(t$rkhs_matrix_tsv) || !file.exists(t$rkhs_matrix_tsv)) stop("RKHS matrix_tsv not found")
      M <- read_square_matrix(t$rkhs_matrix_tsv, ids)
      kind <- tolower(trimws(as.character(t$rkhs_matrix_kind)))
      if (kind == "distance") {
        # distance -> kernel
        D2 <- M
        D2[D2 < 0] <- 0
        if (!is.finite(h) || h <= 0) {
          vals <- as.numeric(D2[upper.tri(D2)])
          h <- median(vals[is.finite(vals)], na.rm = TRUE)
          if (!is.finite(h) || h <= 0) h <- 1
        }
        K <- exp(-D2 / h)
      } else {
        K <- M
      }
      return(list(K = K, model = "RKHS"))
    } else {
      K <- build_kernel_from_markers(X, kernel = kern, h = h)
      return(list(K = K, model = "RKHS"))
    }
  }

  # Marker regression
  term <- list(X = X, model = m)
  if (m == "SpikeSlab") {
    pr <- as_num(t$probIn, NA_real_)
    if (is.finite(pr) && pr > 0 && pr < 1) term$probIn <- pr
  }
  term
}

ETA <- list()

# Kernel 1
k1 <- term_from_prefix("kernel1")
eta1 <- build_term(k1)
if (!is.null(eta1)) ETA[[length(ETA) + 1]] <- eta1

# Kernel 2
k2 <- term_from_prefix("kernel2")
eta2 <- build_term(k2)
if (!is.null(eta2)) ETA[[length(ETA) + 1]] <- eta2

# Covariates as FIXED
if (!is.null(Z)) {
  ETA[[length(ETA) + 1]] <- list(X = Z, model = "FIXED")
}

# Extra ETA JSON (optional)
if (nchar(eta_json) > 0 && file.exists(eta_json)) {
  cat("[gp_bglr_multitrait] reading extra ETA JSON: ", eta_json, "\n", sep="")
  extra <- fromJSON(eta_json, simplifyVector = FALSE)
  if (!is.list(extra)) stop("eta_json must be a JSON array")
  for (obj in extra) {
    if (is.null(obj$model)) stop("Each ETA term in JSON must include 'model'")
    m <- normalize_model(obj$model)
    if (!(m %in% SUPPORTED)) stop(paste0("Unsupported model in eta_json: ", m))
    if (m == "RKHS") {
      if (!is.null(obj$K_tsv) && nchar(as.character(obj$K_tsv)) > 0) {
        K <- read_square_matrix(as.character(obj$K_tsv), ids)
        ETA[[length(ETA) + 1]] <- list(K = K, model = "RKHS")
      } else {
        kern <- if (!is.null(obj$kernel)) as.character(obj$kernel) else "linear"
        h <- as_num(obj$h, NA_real_)
        K <- build_kernel_from_markers(X, kernel = kern, h = h)
        ETA[[length(ETA) + 1]] <- list(K = K, model = "RKHS")
      }
    } else if (m == "FIXED") {
      if (is.null(obj$X_tsv)) stop("FIXED term requires X_tsv")
      XT <- fread(as.character(obj$X_tsv), sep = "\t", header = TRUE, data.table = FALSE)
      if (ncol(XT) < 2) stop("X_tsv must have <ID> + columns")
      ids_x <- as.character(XT[[1]])
      m3 <- match(ids, ids_x)
      if (any(is.na(m3))) stop("X_tsv is missing some IDs")
      XX <- as.matrix(XT[m3, -1, drop = FALSE])
      XX <- apply(XX, 2, function(z) suppressWarnings(as.numeric(z)))
      XX[is.nan(XX)] <- NA_real_
      for (j in seq_len(ncol(XX))) {
        v <- XX[, j]
        if (anyNA(v)) {
          mu <- mean(v, na.rm = TRUE)
          if (!is.finite(mu)) mu <- 0
          v[is.na(v)] <- mu
          XX[, j] <- v
        }
      }
      ETA[[length(ETA) + 1]] <- list(X = XX, model = "FIXED")
    } else {
      # BRR / SpikeSlab marker regression
      term <- list(X = X, model = m)
      if (m == "SpikeSlab" && !is.null(obj$probIn)) {
        pr <- as_num(obj$probIn, NA_real_)
        if (is.finite(pr) && pr > 0 && pr < 1) term$probIn <- pr
      }
      ETA[[length(ETA) + 1]] <- term
    }
  }
}

if (!length(ETA)) stop("ETA is empty (kernel1_model=None and no other terms)")

# ------------ Backend selection ------------
has_bglr <- requireNamespace("BGLR", quietly = TRUE)
backend <- "fallback(BGLR-per-trait)"

extract_yhat <- function(res) {
  if (!is.null(res$yHat)) return(as.matrix(res$yHat))
  if (!is.null(res$YHat)) return(as.matrix(res$YHat))
  if (!is.null(res$yh)) return(as.matrix(res$yh))
  if (!is.null(res$fitted.values)) return(as.matrix(res$fitted.values))
  stop("Could not extract predicted values from Multitrait result.")
}

make_resCov <- function(Y_in) {
  if (!nchar(trimws(resCov_code))) return(NULL)
  env <- new.env(parent = baseenv())
  env$Y <- Y_in
  env$n_traits <- ncol(Y_in)
  env$V <- diag(env$n_traits)
  S <- tryCatch(cov(Y_in, use = "pairwise.complete.obs"), error = function(e) NULL)
  if (is.null(S) || any(!is.finite(S))) S <- diag(env$n_traits)
  env$S <- S
  env$S0 <- S
  tryCatch(
    eval(parse(text = resCov_code), envir = env),
    error = function(e) stop(paste0("Failed to evaluate resCov code: ", e$message))
  )
}

fit_multitrait <- function(Y_in) {
  if (!has_bglr) stop("BGLR package is not installed.")
  if (!("Multitrait" %in% getNamespaceExports("BGLR"))) stop("BGLR::Multitrait is not available in this BGLR version.")

  rc <- NULL
  if (nchar(trimws(resCov_code))) rc <- make_resCov(Y_in)

  # try with resCov if provided; fallback to default if resCov is incompatible
  if (!is.null(rc)) {
    res <- tryCatch(
      BGLR::Multitrait(y = Y_in, ETA = ETA, nIter = nIter, burnIn = burnIn, thin = thin, resCov = rc, verbose = FALSE),
      error = function(e) e
    )
    if (!inherits(res, "error")) {
      backend <<- "BGLR::Multitrait"
      return(res$ETAHat)
    }
    cat("[gp_bglr_multitrait] warning: Multitrait with resCov failed: ", res$message, "\n", sep="")
  }

  res2 <- BGLR::Multitrait(y = Y_in, ETA = ETA, nIter = nIter, burnIn = burnIn, thin = thin, verbose = FALSE)
  backend <<- "BGLR::Multitrait"
  res2$ETAHat
}

fit_bglr_per_trait <- function(y_in) {
  if (!has_bglr) {
    mu <- mean(y_in[is.finite(y_in)], na.rm = TRUE)
    if (!is.finite(mu)) mu <- 0
    return(rep(mu, length(y_in)))
  }
  fm <- BGLR::BGLR(y = y_in, ETA = ETA, nIter = nIter, burnIn = burnIn, thin = thin, verbose = FALSE)
  fm$yHat
}

# ------------ Fit / CV ------------
pred <- matrix(NA_real_, nrow = nrow(Y), ncol = ncol(Y))
colnames(pred) <- colnames(Y)
fold_id <- rep(NA_integer_, nrow(Y))
use_fallback <- FALSE

try({
  if (mode == "fit") {
    pred <- fit_multitrait(Y)
    fold_id[] <- 0L
  } else if (mode == "kfold") {
    n <- nrow(Y)
    perm <- sample.int(n)
    folds <- split(perm, rep(seq_len(k_folds), length.out = n))
    for (k in seq_len(k_folds)) {
      test_idx <- folds[[k]]
      Y2 <- Y
      Y2[test_idx, ] <- NA_real_
      ph <- fit_multitrait(Y2)
      pred[test_idx, ] <- ph[test_idx, , drop = FALSE]
      fold_id[test_idx] <- k
      cat("  fold ", k, "/", k_folds, "\n", sep="")
    }
  } else if (mode == "loo") {
    n <- nrow(Y)
    for (i in seq_len(n)) {
      Y2 <- Y
      Y2[i, ] <- NA_real_
      ph <- fit_multitrait(Y2)
      pred[i, ] <- ph[i, ]
      fold_id[i] <- i
      if (i %% 10 == 0) cat("  loo ", i, "/", n, "\n", sep="")
    }
  }
}, silent = FALSE)

if (all(is.na(pred))) {
  cat("[gp_bglr_multitrait] Multitrait failed; using per-trait fallback\n")
  use_fallback <- TRUE
}

if (use_fallback) {
  backend <- "fallback(BGLR-per-trait)"
  if (mode == "fit") {
    for (j in seq_len(ncol(Y))) pred[, j] <- fit_bglr_per_trait(Y[, j])
    fold_id[] <- 0L
  } else if (mode == "kfold") {
    n <- nrow(Y)
    perm <- sample.int(n)
    folds <- split(perm, rep(seq_len(k_folds), length.out = n))
    for (k in seq_len(k_folds)) {
      test_idx <- folds[[k]]
      for (j in seq_len(ncol(Y))) {
        y <- Y[, j]
        y2 <- y
        y2[test_idx] <- NA_real_
        ph <- fit_bglr_per_trait(y2)
        pred[test_idx, j] <- ph[test_idx]
      }
      fold_id[test_idx] <- k
      cat("  [fallback] fold ", k, "/", k_folds, "\n", sep="")
    }
  } else if (mode == "loo") {
    n <- nrow(Y)
    for (i in seq_len(n)) {
      for (j in seq_len(ncol(Y))) {
        y <- Y[, j]
        y2 <- y
        y2[i] <- NA_real_
        ph <- fit_bglr_per_trait(y2)
        pred[i, j] <- ph[i]
      }
      fold_id[i] <- i
      if (i %% 10 == 0) cat("  [fallback] loo ", i, "/", n, "\n", sep="")
    }
  }
}

# ------------ Write outputs ------------
summary_rows <- list()
fold_rows <- list()

for (ii in seq_along(traits)) {
  tr <- traits[ii]
  tr_idx <- trait_idxs[ii]
  y_true <- as.numeric(Y[, ii])
  y_pred <- as.numeric(pred[, ii])

  DT <- data.frame(
    id = ids,
    trait_idx = tr_idx,
    trait = tr,
    y_true = y_true,
    y_pred = y_pred,
    fold = fold_id,
    stringsAsFactors = FALSE
  )

  out_tsv <- file.path(out_dir, paste0("predictions_", tr_idx, ".tsv"))
  fwrite(DT, out_tsv, sep = "\t")
  if (length(traits) == 1) fwrite(DT, file.path(out_dir, "predictions.tsv"), sep = "\t")

  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    trait_idx = tr_idx,
    trait = tr,
    n = nrow(DT),
    n_obs = sum(is.finite(y_true)),
    mode = mode,
    backend = backend,
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
}

summary_df <- rbindlist(summary_rows, fill = TRUE)
fwrite(summary_df, file.path(out_dir, "summary.tsv"), sep = "\t")

fold_df <- rbindlist(fold_rows, fill = TRUE)
if (nrow(fold_df) == 0) {
  fold_df <- data.frame(trait_idx = integer(0), trait = character(0), fold = integer(0), n_test = integer(0),
                        cor = numeric(0), rmse = numeric(0), mae = numeric(0))
}
fwrite(fold_df, file.path(out_dir, "fold_metrics.tsv"), sep = "\t")

metrics <- list(
  mode = mode,
  n = nrow(Y),
  n_markers = ncol(X),
  n_traits = length(traits),
  traits = traits,
  trait_idxs = trait_idxs,
  backend = backend,
  seed = seed,
  nIter = nIter,
  burnIn = burnIn,
  thin = thin,
  kernel1_model = kernel1_model,
  kernel2_model = kernel2_model,
  resCov_code = resCov_code,
  eta_json = eta_json,
  has_covariates = !is.null(Z)
)
writeLines(toJSON(metrics, auto_unbox = TRUE, pretty = TRUE), file.path(out_dir, "metrics.json"))

# Plot: trait-wise cor
try({
  png(file.path(out_dir, "trait_cor.png"), width = 900, height = 450)
  par(mar = c(7, 4, 2, 1))
  barplot(summary_df$cor, names.arg = summary_df$trait, las = 2)
  abline(h = 0)
  dev.off()
}, silent = TRUE)

art <- list(
  tables = c("summary.tsv", "fold_metrics.tsv"),
  plots = c("trait_cor.png"),
  default_table = "summary.tsv",
  default_plot = "trait_cor.png"
)
writeLines(toJSON(art, auto_unbox = TRUE, pretty = TRUE), file.path(out_dir, "artifacts.json"))

cat("[gp_bglr_multitrait] done\n")
sink()
