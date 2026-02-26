#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(jsonlite)
  library(data.table)
  library(BGLR)
  library(RAINBOWR)
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

cat("[gp_bglr_brr] start\n")
cat("[gp_bglr_brr] params_path=", params_path, "\n")
cat("[gp_bglr_brr] out_dir=", out_dir, "\n")

p <- fromJSON(params_path, simplifyVector = FALSE)

# -----------------------
# Inputs
# -----------------------
geno_tsv <- p$genotype_tsv
pheno_tsv <- p$phenotype_tsv

# Trait selection:
# - If phenotype.tsv has multiple trait columns (excluding the first ID column), default is ALL traits.
# - Optionally specify trait indices (1-based, excluding ID) like: "1", "1-5", "2,3,5".
# - You may also mix in trait names (must match the header).
traits_spec <- ""
if (!is.null(p$traits_spec) && nchar(as.character(p$traits_spec)) > 0) traits_spec <- as.character(p$traits_spec)
if (!is.null(p$traits) && nchar(as.character(p$traits)) > 0) traits_spec <- as.character(p$traits)
if (!is.null(p$trait) && nchar(as.character(p$trait)) > 0) traits_spec <- as.character(p$trait)

stopifnot(!is.null(geno_tsv), file.exists(geno_tsv))
stopifnot(!is.null(pheno_tsv), file.exists(pheno_tsv))
cat("[gp_bglr_brr] genotype_tsv=", geno_tsv, "\n")
cat("[gp_bglr_brr] phenotype_tsv=", pheno_tsv, "\n")
#cat("[gp_bglr_brr] trait=", trait, " id_col=", id_col, "\n")

cov_tsv <- NULL
if (!is.null(p$covariate_tsv) && nchar(p$covariate_tsv) > 0) cov_tsv <- p$covariate_tsv
if (is.null(cov_tsv) && !is.null(p$covariate_file) && nchar(p$covariate_file) > 0) cov_tsv <- p$covariate_file
if (is.null(cov_tsv) && !is.null(p$covariates_tsv) && nchar(p$covariates_tsv) > 0) cov_tsv <- p$covariates_tsv
if (!is.null(cov_tsv)) {
  cov_tsv <- as.character(cov_tsv)
  if (!file.exists(cov_tsv)) stop(paste0("covariate_tsv not found: ", cov_tsv))
  cat("[gp_bglr_brr] covariate_tsv=", cov_tsv, "\n")
}

# -----------------------
# Options
# -----------------------
mode <- if (!is.null(p$mode) && nchar(p$mode) > 0) p$mode else "fit"  # fit / loo / kfold
k_folds <- if (!is.null(p$k_folds)) as.integer(p$k_folds) else 5
seed <- if (!is.null(p$seed)) as.integer(p$seed) else 1

nIter <- if (!is.null(p$nIter)) as.integer(p$nIter) else 2000
burnIn <- if (!is.null(p$burnIn)) as.integer(p$burnIn) else 500
thin <- if (!is.null(p$thin)) as.integer(p$thin) else 1
#verbose <- if (!is.null(p$verbose)) as.logical(p$verbose) else TRUE

center <- if (!is.null(p$center_markers)) as.logical(p$center_markers) else TRUE
scale <- if (!is.null(p$scale_markers)) as.logical(p$scale_markers) else TRUE
impute <- if (!is.null(p$impute_missing)) p$impute_missing else "mean"

cat("[gp_bglr_brr] mode=", mode, " k_folds=", k_folds, " seed=", seed, "\n")
cat("[gp_bglr_brr] nIter=", nIter, " burnIn=", burnIn, " thin=", thin, "\n")
cat("[gp_bglr_brr] center=", center, " scale=", scale, " impute=", impute, "\n")

if (!(mode %in% c("fit", "loo", "kfold"))) stop("mode must be one of: fit, loo, kfold")
if (mode == "kfold" && (is.na(k_folds) || k_folds < 2)) stop("k_folds must be >= 2")

# -----------------------
# Helpers
# -----------------------
as_num <- function(x) {
  if (is.null(x)) return(NA_real_)
  suppressWarnings(as.numeric(as.character(x)))
}

normalize_model <- function(m) {
  if (is.null(m) || nchar(as.character(m)) == 0) return(NULL)
  m <- as.character(m)
  if (m == "BayesCpi") return("BayesC")
  m
}

read_square_matrix <- function(path, ids, kind = "kernel", h = NA_real_) {
  cat("[gp_bglr_brr] reading matrix file: ", path, " (kind=", kind, ")\n", sep = "")
  DT <- fread(path, sep = "\t", header = TRUE, data.table = FALSE)
  if (nrow(DT) == 0) stop(paste0("Matrix file has 0 rows: ", path))

  row_ids <- NULL
  # Detect rownames in first column
  c1 <- DT[[1]]
  if (!is.numeric(c1)) {
    c1c <- as.character(c1)
    if (sum(c1c %in% ids) >= max(3L, floor(0.1 * length(ids)))) {
      row_ids <- c1c
      DT <- DT[, -1, drop = FALSE]
    }
  }

  mat <- as.matrix(DT)
  mat <- apply(mat, 2, function(z) suppressWarnings(as.numeric(z)))
  mat[is.nan(mat)] <- NA

  col_ids <- colnames(DT)
  if (!is.null(col_ids)) col_ids <- as.character(col_ids)

  # Align / subset
  if (!is.null(row_ids) && !is.null(col_ids) && all(ids %in% row_ids) && all(ids %in% col_ids)) {
    rr <- match(ids, row_ids)
    cc <- match(ids, col_ids)
    mat <- mat[rr, cc, drop = FALSE]
  } else if (nrow(mat) == length(ids) && ncol(mat) == length(ids)) {
    # assume already in the same order as ids
    # nothing
  } else {
    stop(paste0(
      "Matrix dimensions / labels do not match ids. ",
      "nrow=", nrow(mat), " ncol=", ncol(mat), " n_ids=", length(ids),
      ". Provide row/col names (ids) or a square matrix in the same order as genotype/phenotype alignment."
    ))
  }

  if (anyNA(mat)) stop("Matrix contains NA after numeric conversion / alignment.")

  # Symmetrize (best effort)
  mat <- (mat + t(mat)) / 2

  if (kind %in% c("distance", "squared_distance", "dist")) {
    D2 <- mat
    D2[D2 < 0] <- 0
    if (!is.finite(h) || is.na(h) || h <= 0) {
      v <- D2[upper.tri(D2, diag = FALSE)]
      v <- v[is.finite(v) & v > 0]
      if (length(v) == 0) {
        h <- 1
      } else {
        h <- stats::median(v)
      }
      cat("[gp_bglr_brr] rkhs_h auto-set to median(distance)=", h, "\n")
    }
    K <- exp(-D2 / h)
    return(list(K = K, h = h))
  }

  list(K = mat, h = h)
}

# -----------------------
# Load tables (genotype / phenotype)
# -----------------------
gt <- fread(geno_tsv, sep = "	", header = TRUE, data.table = FALSE)
if (ncol(gt) < 2) stop("genotype_tsv must have at least 2 columns: <ID> + markers")
id_col <- colnames(gt)[1]
ids_g <- as.character(gt[[1]])
marker_cols <- colnames(gt)[-1]
if (length(marker_cols) == 0) stop("No marker columns in genotype_tsv")

X <- as.matrix(gt[, marker_cols, drop = FALSE])
X <- apply(X, 2, function(z) suppressWarnings(as.numeric(z)))
X[is.nan(X)] <- NA

pt <- fread(pheno_tsv, sep = "	", header = TRUE, data.table = FALSE)
if (ncol(pt) < 2) stop("phenotype_tsv must have at least 2 columns: <ID> + trait(s)")
id_col_p <- colnames(pt)[1]
ids_p <- as.character(pt[[1]])
trait_cols <- colnames(pt)[-1]

# parse trait spec (indices excluding ID; 1-based) or trait names
parse_trait_spec <- function(spec, trait_cols) {
  if (is.null(spec)) spec <- ""
  spec <- trimws(as.character(spec))
  if (nchar(spec) == 0) {
    return(seq_along(trait_cols))
  }
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

# Align by id
m <- match(ids_g, ids_p)
ok <- !is.na(m)
if (!any(ok)) stop("No overlapping IDs between genotype and phenotype")

X <- X[ok, , drop = FALSE]
ids <- ids_g[ok]

cat("[gp_bglr_brr] n_samples=", length(ids), " n_markers=", ncol(X), "\n")

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
} else {
  for (j in seq_len(ncol(X))) {
    v <- X[, j]
    if (anyNA(v)) {
      mu <- median(v, na.rm = TRUE)
      if (!is.finite(mu)) mu <- 0
      v[is.na(v)] <- mu
      X[, j] <- v
    }
  }
}

X_ori <- X
# Center/scale markers
if (center || scale) {
  sc <- scale(X, center = center, scale = scale)
  X <- as.matrix(sc)
}
# Drop markers that became NA (e.g. zero variance after scaling)
X <- t(na.omit(t(X)))

# -----------------------
# Covariates (FIXED)
# -----------------------
Z <- NULL
cov_cols <- character(0)
if (!is.null(cov_tsv)) {
  ct <- fread(cov_tsv, sep = "\t", header = TRUE, data.table = FALSE)
  if (!(id_col %in% colnames(ct))) stop(paste0("covariate_tsv must contain id column: ", id_col))

  ids_c <- as.character(ct[[id_col]])
  cov_cols <- setdiff(colnames(ct), id_col)
  if (length(cov_cols) == 0) {
    cat("[gp_bglr_brr] WARN: covariate_tsv has no covariate columns (only id_col). Ignoring.\n")
  } else {
    mmc <- match(ids, ids_c)
    if (all(is.na(mmc))) {
      cat("[gp_bglr_brr] WARN: no overlapping IDs between covariate and genotype/phenotype. Ignoring covariates.\n")
    } else {
      Z <- as.matrix(ct[mmc, cov_cols, drop = FALSE])
      Z <- apply(Z, 2, function(z) suppressWarnings(as.numeric(z)))
      Z[is.nan(Z)] <- NA
      # mean-impute covariates
      for (j in seq_len(ncol(Z))) {
        v <- Z[, j]
        if (anyNA(v)) {
          mu <- mean(v, na.rm = TRUE)
          if (!is.finite(mu)) mu <- 0
          v[is.na(v)] <- mu
          Z[, j] <- v
        }
      }
      cat("[gp_bglr_brr] covariates added: n_cov=", ncol(Z), "\n")
    }
  }
}

# -----------------------
# ETA terms (multi-kernel)
# -----------------------
# New params (preferred): kernel1_*, kernel2_*, eta_json
# Backward compatible: model_type/probIn/rkhs_kernel/rkhs_h

get_term_from_prefix <- function(prefix) {
  m <- normalize_model(paste0(""))
  # Access via p[[..]]
  m <- normalize_model(p[[paste0(prefix, "_model")]])
  if (is.null(m) || nchar(m) == 0) return(NULL)
  if (m == "None") return(NULL)

  term <- list(
    model = m,
    probIn = as_num(p[[paste0(prefix, "_probIn")]]),
    rkhs_source = if (!is.null(p[[paste0(prefix, "_rkhs_source")]]) && nchar(p[[paste0(prefix, "_rkhs_source")]]) > 0) as.character(p[[paste0(prefix, "_rkhs_source")]]) else "markers",
    rkhs_kernel = if (!is.null(p[[paste0(prefix, "_rkhs_kernel")]]) && nchar(p[[paste0(prefix, "_rkhs_kernel")]]) > 0) as.character(p[[paste0(prefix, "_rkhs_kernel")]]) else "addNOIA",
    rkhs_h = as_num(p[[paste0(prefix, "_rkhs_h")]]),
    rkhs_matrix_tsv = if (!is.null(p[[paste0(prefix, "_rkhs_matrix_tsv")]]) && nchar(p[[paste0(prefix, "_rkhs_matrix_tsv")]]) > 0) as.character(p[[paste0(prefix, "_rkhs_matrix_tsv")]]) else "",
    rkhs_matrix_kind = if (!is.null(p[[paste0(prefix, "_rkhs_matrix_kind")]]) && nchar(p[[paste0(prefix, "_rkhs_matrix_kind")]]) > 0) as.character(p[[paste0(prefix, "_rkhs_matrix_kind")]]) else "kernel"
  )
  term
}

terms <- list()

k1 <- get_term_from_prefix("kernel1")
if (!is.null(k1)) {
  terms[[length(terms) + 1]] <- k1
} else {
  # backward compatible single-kernel
  mt <- normalize_model(p$model_type)
  if (is.null(mt)) mt <- "BRR"
  term <- list(
    model = mt,
    probIn = as_num(p$probIn),
    rkhs_source = "markers",
    rkhs_kernel = if (!is.null(p$rkhs_kernel) && nchar(p$rkhs_kernel) > 0) as.character(p$rkhs_kernel) else "addNOIA",
    rkhs_h = as_num(p$rkhs_h),
    rkhs_matrix_tsv = "",
    rkhs_matrix_kind = "kernel"
  )
  terms[[length(terms) + 1]] <- term
}

k2 <- get_term_from_prefix("kernel2")
if (!is.null(k2)) terms[[length(terms) + 1]] <- k2

eta_json <- NULL
if (!is.null(p$eta_json) && nchar(p$eta_json) > 0) eta_json <- as.character(p$eta_json)
if (!is.null(p$eta_terms_json) && nchar(p$eta_terms_json) > 0) eta_json <- as.character(p$eta_terms_json)
if (!is.null(eta_json)) {
  if (!file.exists(eta_json)) stop(paste0("eta_json not found: ", eta_json))
  cat("[gp_bglr_brr] extra ETA json=", eta_json, "\n")
  extra_terms <- fromJSON(eta_json, simplifyVector = FALSE)
  # allow single object
  if (!is.null(extra_terms$model) && is.null(extra_terms[[1]])) {
    extra_terms <- list(extra_terms)
  }
  if (!is.list(extra_terms)) stop("eta_json must be a JSON array (list) of ETA term objects")
  for (tt in extra_terms) {
    if (is.null(tt$model)) stop("eta_json term missing 'model'")
    t2 <- list(
      model = normalize_model(tt$model),
      probIn = as_num(tt$probIn),
      rkhs_source = if (!is.null(tt$rkhs_source) && nchar(tt$rkhs_source) > 0) as.character(tt$rkhs_source) else if (!is.null(tt$rkhs_matrix_tsv) && nchar(tt$rkhs_matrix_tsv) > 0) "matrix" else "markers",
      rkhs_kernel = if (!is.null(tt$rkhs_kernel) && nchar(tt$rkhs_kernel) > 0) as.character(tt$rkhs_kernel) else "addNOIA",
      rkhs_h = as_num(tt$rkhs_h),
      rkhs_matrix_tsv = if (!is.null(tt$rkhs_matrix_tsv) && nchar(tt$rkhs_matrix_tsv) > 0) as.character(tt$rkhs_matrix_tsv) else "",
      rkhs_matrix_kind = if (!is.null(tt$rkhs_matrix_kind) && nchar(tt$rkhs_matrix_kind) > 0) as.character(tt$rkhs_matrix_kind) else "kernel"
    )
    terms[[length(terms) + 1]] <- t2
  }
}

# Validate models
supported <- c("BRR", "BayesA", "BayesB", "BayesC", "BL", "RKHS")
for (i in seq_along(terms)) {
  if (!(terms[[i]]$model %in% supported)) {
    stop(paste0("Unsupported model for term ", i, ": ", terms[[i]]$model, ". Supported: ", paste(supported, collapse = ", ")))
  }
}

cat("[gp_bglr_brr] ETA terms: ", length(terms), " (covariates=", if (!is.null(Z)) "YES" else "NO", ")\n", sep = "")
for (i in seq_along(terms)) {
  tt <- terms[[i]]
  cat("  term", i, ": model=", tt$model,
      " probIn=", tt$probIn,
      " rkhs_source=", tt$rkhs_source,
      " rkhs_kernel=", tt$rkhs_kernel,
      " rkhs_h=", tt$rkhs_h,
      " rkhs_matrix_tsv=", tt$rkhs_matrix_tsv,
      " rkhs_matrix_kind=", tt$rkhs_matrix_kind,
      "\n", sep="")
}

# -----------------------
# Prepare ETA list (precompute K per RKHS term)
# -----------------------
ETA <- list()
eta_meta <- list()

if (!is.null(Z)) {
  ETA[[length(ETA) + 1]] <- list(X = Z, model = "FIXED")
  eta_meta[[length(eta_meta) + 1]] <- list(model = "FIXED", n_cov = ncol(Z), columns = cov_cols)
}

for (i in seq_along(terms)) {
  tt <- terms[[i]]
  if (tt$model == "RKHS") {
    # matrix from file (kernel or distance) OR build from markers
    use_file <- FALSE
    if (!is.null(tt$rkhs_matrix_tsv) && nchar(tt$rkhs_matrix_tsv) > 0) use_file <- TRUE
    if (tt$rkhs_source %in% c("matrix", "file", "matrix_file")) use_file <- TRUE

    if (use_file) {
      path <- tt$rkhs_matrix_tsv
      if (is.null(path) || nchar(path) == 0) stop(paste0("RKHS term ", i, " requires rkhs_matrix_tsv when rkhs_source=matrix"))
      if (!file.exists(path)) stop(paste0("rkhs_matrix_tsv not found: ", path))
      kk <- read_square_matrix(path, ids, kind = tt$rkhs_matrix_kind, h = tt$rkhs_h)
      ETA[[length(ETA) + 1]] <- list(K = kk$K, model = "RKHS")
      eta_meta[[length(eta_meta) + 1]] <- list(model = "RKHS", source = "matrix", matrix = basename(path), matrix_kind = tt$rkhs_matrix_kind, rkhs_h = kk$h)
      # Update h in term (for output)
      terms[[i]]$rkhs_h <- kk$h
    } else {
      #kk <- build_kernel_from_markers(X, kernel = tt$rkhs_kernel, h = tt$rkhs_h)
      kk <- RAINBOWR::calcGRM(as.matrix(X_ori), 
                              methodGRM = tt$rkhs_kernel)
      ETA[[length(ETA) + 1]] <- list(K = kk, model = "RKHS")
      eta_meta[[length(eta_meta) + 1]] <- list(model = "RKHS", source = "markers", rkhs_kernel = tt$rkhs_kernel)
      #terms[[i]]$rkhs_h <- kk$h
    }
  } else {
    e1 <- list(X = X, model = tt$model)
    if ((tt$model %in% c("BayesB", "BayesC")) && is.finite(tt$probIn) && !is.na(tt$probIn) && tt$probIn > 0 && tt$probIn < 1) {
      e1$probIn <- tt$probIn
    }
    ETA[[length(ETA) + 1]] <- e1
    eta_meta[[length(eta_meta) + 1]] <- list(model = tt$model, source = "markers", probIn = if ((tt$model %in% c("BayesB", "BayesC"))) tt$probIn else NULL)
  }
}

# -----------------------
# Metrics helper
# -----------------------
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

# -----------------------
# Core runner wrapper
# -----------------------
fit_once <- function(y_in) {
  fm <- BGLR(
    y = y_in,
    ETA = ETA,
    nIter = nIter,
    burnIn = burnIn,
    thin = thin,
    verbose = F
  )
  fm
}


# -----------------------
# Run per-trait
# -----------------------

run_one_trait <- function(y) {
  pred <- rep(NA_real_, length(y))
  fold_id <- rep(NA_integer_, length(y))
  fold_metrics <- NULL

  if (mode == "fit") {
    fm <- fit_once(y)
    pred <- fm$yHat
    fold_id <- rep(0L, length(y))
  } else if (mode == "loo") {
    cat("[gp_bglr_brr] running LOO-CV ...\n")
    for (i in seq_along(y)) {
      y2 <- y
      y2[i] <- NA
      fm <- fit_once(y2)
      pred[i] <- fm$yHat[i]
      fold_id[i] <- i
      if (i %% 10 == 0) cat("  ...", i, "/", length(y), "\n")
    }
  } else if (mode == "kfold") {
    cat("[gp_bglr_brr] running k-fold CV ...\n")
    n <- length(y)
    perm <- sample.int(n)
    folds <- split(perm, rep(seq_len(k_folds), length.out = n))
    fold_metrics <- data.frame(
      fold = integer(0),
      n_test = integer(0),
      cor = numeric(0),
      rmse = numeric(0),
      mae = numeric(0),
      stringsAsFactors = FALSE
    )
    for (k in seq_len(k_folds)) {
      te <- folds[[k]]
      y2 <- y
      y2[te] <- NA
      fm <- fit_once(y2)
      pred[te] <- fm$yHat[te]
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
  } else {
    stop("mode must be one of: fit, loo, kfold")
  }

  list(
    pred = pred,
    fold_id = fold_id,
    fold_metrics = fold_metrics,
    cor = pearson(y, pred),
    rmse = rmse(y, pred),
    mae = mae(y, pred),
    n = length(y),
    n_obs = sum(is.finite(y))
  )
}

summary_rows <- list()
fold_rows <- list()
artifact_tables <- character(0)

for (ii in seq_along(traits)) {
  tr_name <- traits[ii]
  tr_idx <- trait_idxs[ii]
  cat("[gp_bglr_brr] trait_idx=", tr_idx, " trait=", tr_name, "\n")
  y_raw <- suppressWarnings(as.numeric(pt[[tr_name]]))
  y <- y_raw[m[ok]]

  res <- run_one_trait(y)

  # write predictions
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
  fwrite(pred_df, file.path(out_dir, out_name_idx), sep = "	")
  artifact_tables <- c(artifact_tables, out_name_idx)

  if (length(traits) == 1) {
    fwrite(pred_df, file.path(out_dir, "predictions.tsv"), sep = "	")
    artifact_tables <- c(artifact_tables, "predictions.tsv")
  }

  # summary row
  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    trait_idx = tr_idx,
    trait = tr_name,
    n = res$n,
    n_obs = res$n_obs,
    mode = mode,
    backend = "BGLR",
    cor = res$cor,
    rmse = res$rmse,
    mae = res$mae,
    stringsAsFactors = FALSE
  )

  # fold metrics (kfold)
  if (!is.null(res$fold_metrics)) {
    fm <- res$fold_metrics
    fm$trait_idx <- tr_idx
    fm$trait <- tr_name
    fold_rows[[length(fold_rows) + 1]] <- fm
  }
}

summary_df <- rbindlist(summary_rows, fill = TRUE)
summary_tsv <- file.path(out_dir, "summary.tsv")
fwrite(summary_df, summary_tsv, sep = "\t")
artifact_tables <- unique(c("summary.tsv", artifact_tables))

fold_metrics_tsv <- NULL
if (length(fold_rows) > 0) {
  fold_df <- rbindlist(fold_rows, fill = TRUE)
  fold_metrics_tsv <- file.path(out_dir, "fold_metrics.tsv")
  fwrite(fold_df, fold_metrics_tsv, sep = "\t")
  artifact_tables <- unique(c(artifact_tables, "fold_metrics.tsv"))
}

# simple barplot across traits
try({
  png(file.path(out_dir, "trait_cor.png"), width = 900, height = 450)
  par(mar=c(7,4,2,1))
  barplot(summary_df$cor, names.arg = paste0(summary_df$trait_idx, ":", summary_df$trait), las = 2)
  abline(h=0)
  dev.off()
}, silent = TRUE)

# summary metrics.json (mean across traits)
mean_or_na <- function(v) {
  v <- as.numeric(v)
  if (!length(v) || all(!is.finite(v))) return(NA_real_)
  mean(v[is.finite(v)])
}

overall <- list(
  mode = mode,
  model = "BGLR",
  n = length(ids),
  n_markers = ncol(X),
  n_traits = length(traits),
  traits = as.list(traits),
  cor = mean_or_na(summary_df$cor),
  rmse = mean_or_na(summary_df$rmse),
  mae = mean_or_na(summary_df$mae),
  eta = eta_meta
)

metrics_json <- file.path(out_dir, "metrics.json")
writeLines(toJSON(overall, auto_unbox = TRUE, pretty = TRUE), metrics_json)

# Artifacts
artifacts <- list(
  plugin_id = "gp_bglr_brr",
  tables = artifact_tables,
  plots = c("trait_cor.png"),
  default_table = if (length(traits) > 1) "summary.tsv" else "predictions.tsv",
  default_plot = "trait_cor.png",
  outputs = list(
    metrics_json = basename(metrics_json)
  )
)
writeLines(toJSON(artifacts, auto_unbox = TRUE, pretty = TRUE), file.path(out_dir, "artifacts.json"))

cat("[gp_bglr_brr] done\n")

sink()

