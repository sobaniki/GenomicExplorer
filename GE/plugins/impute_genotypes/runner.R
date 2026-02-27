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

cat("[impute_genotypes] start\n")
cat("[impute_genotypes] params_path=", params_path, "\n")
cat("[impute_genotypes] out_dir=", out_dir, "\n")

p <- fromJSON(params_path)

# Inputs
input_mode <- if (!is.null(p$input_mode) && nchar(p$input_mode) > 0) tolower(p$input_mode) else "tsv"

geno_tsv <- if (!is.null(p$genotype_tsv)) p$genotype_tsv else NULL
vcf_path <- if (!is.null(p$vcf_path)) p$vcf_path else NULL
marker_map_tsv <- if (!is.null(p$marker_map_tsv)) p$marker_map_tsv else NULL

method <- if (!is.null(p$method) && nchar(p$method) > 0) tolower(p$method) else "mean"

# Output prefix (optional)
out_prefix <- if (!is.null(p$out_prefix) && nchar(p$out_prefix) > 0) p$out_prefix else file.path(out_dir, "imputed")

# Common options
round_to_012 <- if (!is.null(p$round_to_012)) isTRUE(p$round_to_012) else TRUE

# Quick accuracy check (mask evaluation)
quick_accuracy_check <- if (!is.null(p$quick_accuracy_check)) isTRUE(p$quick_accuracy_check) else FALSE
mask_rate <- if (!is.null(p$mask_rate)) as.numeric(p$mask_rate) else 0.01
mask_seed <- if (!is.null(p$mask_seed)) as.integer(p$mask_seed) else 1L

na_strings <- c("NA", ".", "-", "nan", "NaN", "NAN", "")
if (!is.null(p$na_strings) && length(p$na_strings) > 0) {
  na_strings <- unique(c(na_strings, as.character(p$na_strings)))
}

cat("[impute_genotypes] input_mode=", input_mode, " method=", method, "\n")
cat("[impute_genotypes] out_prefix=", out_prefix, " round_to_012=", round_to_012, "\n")

# --- helpers ---
write_artifacts <- function(lst) {
  # jsonlite cannot serialize some S3 objects (e.g., class 'table') by default.
  # Sanitize recursively to keep artifacts robust.
  sanitize_for_json <- function(x) {
    if (inherits(x, "table")) return(unclass(x))
    if (inherits(x, "ftable")) return(unclass(as.matrix(x)))
    if (is.factor(x)) return(as.character(x))
    if (inherits(x, "Date") || inherits(x, "POSIXct") || inherits(x, "POSIXt")) return(as.character(x))
    if (is.data.frame(x)) {
      out <- lapply(x, sanitize_for_json)
      out <- as.data.frame(out, stringsAsFactors = FALSE, check.names = FALSE)
      return(out)
    }
    if (is.list(x)) return(lapply(x, sanitize_for_json))
    return(x)
  }
  lst2 <- sanitize_for_json(lst)
  writeLines(toJSON(lst2, auto_unbox = TRUE, pretty = TRUE), file.path(out_dir, "artifacts.json"))
}

calc_missing_rate <- function(X) {
  mean(is.na(X))
}

clip01 <- function(x, lo, hi) {
  x[x < lo] <- lo
  x[x > hi] <- hi
  x
}

# Convert numeric 0/1/2 to GT string
num_to_gt <- function(x) {
  if (is.na(x)) return("./.")
  if (x <= 0.5) return("0/0")
  if (x < 1.5) return("0/1")
  return("1/1")
}

# Convert GT string (0/0, 0/1, 1/1, phased allowed) to 0/1/2
parse_gt_to_dosage <- function(gt) {
  if (is.na(gt)) return(NA_real_)
  if (gt == "." || gt == "./." || gt == ".|.") return(NA_real_)
  # take only GT (before :)
  gt2 <- strsplit(gt, ":", fixed = TRUE)[[1]][1]
  if (gt2 == "./." || gt2 == ".|.") return(NA_real_)
  gt2 <- gsub("\\|", "/", gt2)
  a <- strsplit(gt2, "/", fixed = TRUE)[[1]]
  if (length(a) != 2) return(NA_real_)
  if (a[1] == "." || a[2] == ".") return(NA_real_)
  as.numeric(a[1]) + as.numeric(a[2])
}

# Load genotype TSV into (ids, X numeric)
load_geno_tsv <- function(path) {
  gt <- fread(path, sep = "\t", header = TRUE, data.table = FALSE, na.strings = na_strings)
  #if (!("id" %in% colnames(gt))) stop("genotype_tsv must have 'id' column")
  #ids <- as.character(gt$id)
  #markers <- setdiff(colnames(gt), "id")
  #if (length(markers) == 0) stop("No marker columns in genotype_tsv")
  ids <- as.character(gt[, 1])
  markers <- colnames(gt)[2:ncol(gt)]
  #X <- as.matrix(gt[, markers, drop = FALSE])
  X <- as.matrix(gt[, 2:ncol(gt), drop = FALSE])
  suppressWarnings(storage.mode(X) <- "numeric")

  # Safety: ensure ids length matches matrix rows.
  # In some malformed TSVs, fread() may keep the 'id' column length while the
  # marker matrix is shorter (e.g., due to embedded tabs/newlines). Beagle
  # requires exact sample count match, so reconcile here with a loud message.
  if (length(ids) != nrow(X)) {
    cat(sprintf("[WARN] genotype_tsv rows mismatch: length(id)=%d but nrow(X)=%d. ", length(ids), nrow(X)))
    if (length(ids) > nrow(X)) {
      cat("Truncating ids to nrow(X).\n")
      ids <- ids[seq_len(nrow(X))]
    } else {
      cat("Padding ids to nrow(X).\n")
      ids <- c(ids, paste0("sample_", seq_len(nrow(X) - length(ids)) + length(ids)))
    }
  }
  list(ids = ids, markers = markers, X = X)
}

# Save genotype matrix to TSV (wide)
save_geno_tsv <- function(ids, markers, X, path) {
  dt <- data.table(id = ids)
  dt <- cbind(dt, as.data.table(X))
  setnames(dt, c("id", markers))
  #dt <- data.frame(id = ids,
  #                 as.data.frame(X))
  fwrite(dt, path, sep = "\t", na = "NA", quote = F)
}

# Summaries
make_marker_summary <- function(markers, X_before, X_after) {
  miss_before <- colMeans(is.na(X_before))
  miss_after  <- colMeans(is.na(X_after))

  calc_maf <- function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0L) return(NA_real_)
    pA <- mean(x) / 2
    pA <- min(max(pA, 0.0), 1.0)
    min(pA, 1.0 - pA)
  }

  maf_before <- apply(X_before, 2, calc_maf)
  maf_after  <- apply(X_after,  2, calc_maf)

  data.frame(marker = markers,
             missing_before = miss_before,
             missing_after = miss_after,
             maf_before = maf_before,
             maf_after = maf_after,
             stringsAsFactors = FALSE)
}

make_sample_summary <- function(ids, X_before, X_after) {
  miss_before <- rowMeans(is.na(X_before))
  miss_after  <- rowMeans(is.na(X_after))
  data.frame(id = ids,
             missing_before = miss_before,
             missing_after = miss_after,
             stringsAsFactors = FALSE)
}

impute_matrix <- function(X_in, method, p) {
  Ximp <- X_in
  if (method == "mean") {
    m <- colMeans(X_in, na.rm = TRUE)
    # If a column is entirely NA, colMeans returns NaN; fall back to 1 (heterozygous)
    m[is.nan(m)] <- 1
    for (j in seq_along(m)) {
      Ximp[is.na(Ximp[, j]), j] <- m[j]
    }
  } else if (method == "em") {
    suppressPackageStartupMessages(library(rrBLUP))
    # rrBLUP expects {-1,0,1}
    X_rr <- X_in - 1
    res <- rrBLUP::A.mat(
      X = X_rr,
      min.MAF = 0,
      max.missing = 1,
      impute.method = "EM",
      tol = if (!is.null(p$em_tol)) as.numeric(p$em_tol) else 0.02,
      n.core = if (!is.null(p$em_n_core)) as.integer(p$em_n_core) else 1,
      shrink = F,
      return.imputed = TRUE
    )
    X_rr_imp <- res$imputed
    if (isTRUE(!is.null(p$em_clip) && isTRUE(p$em_clip))) {
      X_rr_imp <- clip01(X_rr_imp, -1, 1)
    }
    Ximp <- X_rr_imp + 1
  } else if (method == "rf") {
    suppressPackageStartupMessages(library(missRanger))

    rf_max_markers <- if (!is.null(p$rf_max_markers)) as.integer(p$rf_max_markers) else 5000L
    rf_allow_large <- if (!is.null(p$rf_allow_large)) isTRUE(p$rf_allow_large) else FALSE

    if (ncol(X_in) > rf_max_markers && !rf_allow_large) {
      stop(sprintf("RF imputation disabled: n_markers=%d exceeds rf_max_markers=%d. Set rf_allow_large=TRUE to force.",
                   ncol(X_in), rf_max_markers))
    }

    df <- as.data.frame(X_in)
    num_trees <- if (!is.null(p$rf_num_trees)) as.integer(p$rf_num_trees) else 100L
    maxiter   <- if (!is.null(p$rf_maxiter)) as.integer(p$rf_maxiter) else 5L
    seed      <- if (!is.null(p$rf_seed)) as.integer(p$rf_seed) else 1L
    pmm_k     <- if (!is.null(p$rf_pmm_k)) as.integer(p$rf_pmm_k) else 3L

    set.seed(seed)
    df_imp <- missRanger::missRanger(
      data = df,
      num.trees = num_trees,
      maxiter = maxiter,
      pmm.k = pmm_k,
      verbose = 1
    )
    Ximp <- as.matrix(df_imp)
    storage.mode(Ximp) <- "numeric"
  }
  Ximp
}

round_to_012_matrix <- function(X) {
  Xr <- round(X)
  Xr[Xr < 0] <- 0
  Xr[Xr > 2] <- 2
  Xr
}

quick_accuracy_eval <- function(X0, method, p, mask_rate, mask_seed) {
  obs_idx <- which(!is.na(X0) & is.finite(X0))
  n_obs <- length(obs_idx)
  if (n_obs < 100L) {
    return(list(status = "skipped", reason = "too_few_observed", n_observed = n_obs))
  }

  n_mask <- as.integer(round(mask_rate * n_obs))
  # Cap masked cells to keep this check "quick".
  max_mask <- 200000L
  if (n_mask > max_mask) n_mask <- max_mask
  if (n_mask < 100L) {
    return(list(status = "skipped", reason = "n_mask_too_small", n_observed = n_obs, n_mask = n_mask, mask_rate = mask_rate))
  }

  set.seed(mask_seed)
  idx_mask <- sample(obs_idx, n_mask)

  X_train <- X0
  X_train[idx_mask] <- NA

  X_imp_cv <- impute_matrix(X_train, method, p)

  y_true <- X0[idx_mask]
  y_pred <- X_imp_cv[idx_mask]

  # classification-style metrics (0/1/2)
  y_true_r <- round_to_012_matrix(y_true)
  y_pred_r <- round_to_012_matrix(y_pred)

  acc <- mean(y_true_r == y_pred_r, na.rm = TRUE)
  rmse <- sqrt(mean((y_pred - y_true)^2, na.rm = TRUE))
  mae  <- mean(abs(y_pred - y_true), na.rm = TRUE)

  acc_by_class <- c(NA_real_, NA_real_, NA_real_)
  names(acc_by_class) <- c("acc_0", "acc_1", "acc_2")
  for (k in 0:2) {
    sel <- which(y_true_r == k)
    if (length(sel) > 0) {
      acc_by_class[k + 1] <- mean(y_pred_r[sel] == k, na.rm = TRUE)
    }
  }

  # confusion table (3x3)
  # NOTE: keep as matrix to be JSON-serializable
  conf <- as.matrix(table(factor(y_true_r, levels = 0:2), factor(y_pred_r, levels = 0:2)))

  list(
    status = "ok",
    method = method,
    mask_rate = mask_rate,
    mask_seed = mask_seed,
    n_observed = n_obs,
    n_mask = n_mask,
    acc_overall = acc,
    rmse = rmse,
    mae = mae,
    acc_0 = acc_by_class[["acc_0"]],
    acc_1 = acc_by_class[["acc_1"]],
    acc_2 = acc_by_class[["acc_2"]],
    confusion = conf
  )
}

# ---------------- main ----------------

if (method %in% c("mean", "em", "rf")) {
  if (is.null(geno_tsv) || !file.exists(geno_tsv)) stop("genotype_tsv is required for method mean/em/rf")
  cat("[impute_genotypes] genotype_tsv=", geno_tsv, "\n")
  dat <- load_geno_tsv(geno_tsv)
  
  ids <- dat$ids
  markers <- dat$markers
  X0 <- dat$X

  cat("[impute_genotypes] n_samples=", nrow(X0), " n_markers=", ncol(X0), "\n")
  cat("[impute_genotypes] missing_rate(before)=", calc_missing_rate(X0), "\n")

  Ximp <- impute_matrix(X0, method, p)

  # Optional: quick accuracy check by masking observed cells
  qc <- NULL
  if (isTRUE(quick_accuracy_check)) {
    cat(sprintf("[impute_genotypes] quick_accuracy_check=TRUE (mask_rate=%.4f seed=%d)\n", mask_rate, mask_seed))
    qc <- tryCatch(
      quick_accuracy_eval(X0, method, p, mask_rate, mask_seed),
      error = function(e) list(status = "error", reason = as.character(e))
    )

    # Write summary tables (best-effort)
    if (!is.null(qc) && is.list(qc) && identical(qc$status, "ok")) {
      acc_df <- data.frame(
        method = qc$method,
        mask_rate = qc$mask_rate,
        mask_seed = qc$mask_seed,
        n_observed = qc$n_observed,
        n_mask = qc$n_mask,
        acc_overall = qc$acc_overall,
        acc_0 = qc$acc_0,
        acc_1 = qc$acc_1,
        acc_2 = qc$acc_2,
        rmse = qc$rmse,
        mae = qc$mae,
        stringsAsFactors = FALSE
      )
      fwrite(acc_df, file.path(out_dir, "accuracy_summary.tsv"), sep = "\t", quote = F, na = "NA")

      # Confusion matrix (rows=true, cols=pred)
      conf <- qc$confusion
      conf_df <- data.frame(true = rownames(conf), conf, check.names = FALSE)
      fwrite(conf_df, file.path(out_dir, "accuracy_confusion.tsv"), sep = "\t", quote = F, na = "NA")
    } else if (!is.null(qc) && is.list(qc)) {
      # Record why it was skipped/failed
      acc_df <- data.frame(
        status = as.character(qc$status),
        reason = as.character(qc$reason),
        mask_rate = mask_rate,
        mask_seed = mask_seed,
        stringsAsFactors = FALSE
      )
      fwrite(acc_df, file.path(out_dir, "accuracy_summary.tsv"), sep = "\t", quote = F, na = "NA")
    }
  }

  # outputs
  #out_tsv_cont <- paste0(out_prefix, ".tsv")
  out_tsv_cont <- file.path(out_dir, "geno_imp.tsv")
  save_geno_tsv(ids, markers, Ximp, out_tsv_cont)

  out_tsv_round <- NULL
  Xround <- NULL
  if (round_to_012) {
    Xround <- round_to_012_matrix(Ximp)
    #out_tsv_round <- paste0(out_prefix, "_rounded.tsv")
    out_tsv_round <- file.path(out_dir, "geno_imp_rounded.tsv")
    save_geno_tsv(ids, markers, Xround, out_tsv_round)
  }

  # summaries
  marker_sum <- make_marker_summary(markers, X0, if (round_to_012) Xround else Ximp)
  sample_sum <- make_sample_summary(ids, X0, if (round_to_012) Xround else Ximp)
  fwrite(marker_sum, file.path(out_dir, "marker_summary.tsv"), sep = "\t", quote = F, na = "NA")
  fwrite(sample_sum, file.path(out_dir, "sample_summary.tsv"), sep = "\t", quote = F, na = "NA")

  # Keep artifacts.json lightweight and robust: confusion matrix is written as TSV;
  # store only scalar metrics + file path in artifacts.
  qc_art <- qc
  if (!is.null(qc_art) && is.list(qc_art) && !is.null(qc_art$confusion)) {
    qc_art$confusion <- NULL
    qc_art$confusion_tsv <- file.path(out_dir, "accuracy_confusion.tsv")
  }

  write_artifacts(list(
    method = method,
    input_mode = "tsv",
    genotype_tsv_in = geno_tsv,
    imputed_tsv = out_tsv_cont,
    imputed_tsv_rounded = out_tsv_round,
    n_samples = nrow(X0),
    n_markers = ncol(X0),
    missing_before = calc_missing_rate(X0),
    missing_after = calc_missing_rate(if (round_to_012) Xround else Ximp),
    quick_accuracy_check = qc_art
  ))

  cat("[impute_genotypes] done\n")

} else if (method == "beagle") {
  # Beagle requires VCF input
  qc <- NULL
  if (isTRUE(quick_accuracy_check)) {
    qc <- list(
      status = "skipped",
      reason = "quick_accuracy_check_not_supported_for_beagle",
      mask_rate = mask_rate,
      mask_seed = mask_seed
    )
    fwrite(data.frame(status = qc$status, reason = qc$reason, mask_rate = mask_rate, mask_seed = mask_seed, stringsAsFactors = FALSE),
           file.path(out_dir, "accuracy_summary.tsv"), sep = "\t", quote = F, na = "NA")
    cat("[WARN] quick_accuracy_check requested but not supported for beagle. Skipping.\n")
  }
  beagle_jar <- if (!is.null(p$beagle_jar)) p$beagle_jar else NULL
  if (is.null(beagle_jar) || !file.exists(beagle_jar)) {
    stop("beagle_jar is required (path to beagle*.jar)")
  }

  nthreads <- if (!is.null(p$beagle_nthreads)) as.integer(p$beagle_nthreads) else 4L
  java_mem <- if (!is.null(p$beagle_java_mem)) as.character(p$beagle_java_mem) else "4g"
  
  #Beagle options
  burnin <- p$burnin
  iterations <- p$iterations
  phase_states <- p$phase_states
  ne <- p$ne
  em <- p$em
  window <- p$window
  overlap <- p$overlap
  seed <- p$seed

  input_vcf <- NULL

  if (!is.null(vcf_path) && file.exists(vcf_path)) {
    input_vcf <- vcf_path
    cat("[impute_genotypes] vcf_path=", vcf_path, "\n")
  } else {
    if (is.null(geno_tsv) || !file.exists(geno_tsv)) stop("For beagle: set vcf_path or genotype_tsv")
    if (is.null(marker_map_tsv) || !file.exists(marker_map_tsv)) {
      stop("For beagle with TSV input: marker_map_tsv with columns marker, chr, pos is required")
    }

    cat("[impute_genotypes] building VCF from TSV\n")
    dat <- load_geno_tsv(geno_tsv)
    ids <- dat$ids
    markers <- dat$markers
    X0 <- dat$X

    # VCF sample IDs must not contain whitespace and should be unique.
    ids_raw <- ids
    ids <- gsub("[[:space:]]+", "_", ids)
    ids <- make.unique(ids)
    if (any(ids != ids_raw)) {
      cat("[WARN] Sample IDs contained whitespace/duplicates. Writing sanitized IDs to VCF and saving mapping file.\n")
      idmap <- data.frame(original_id = ids_raw, vcf_id = ids, stringsAsFactors = FALSE)
      fwrite(idmap, file.path(out_dir, "sample_id_map.tsv"), sep = "\t", quote = F, na = "NA")
    }

    mm <- fread(marker_map_tsv, sep = "\t", header = TRUE, data.table = FALSE)
    if (!("marker" %in% colnames(mm))) stop("marker_map_tsv must have 'marker' column")
    if (!("chr" %in% colnames(mm))) stop("marker_map_tsv must have 'chr' column")
    if (!("pos" %in% colnames(mm))) stop("marker_map_tsv must have 'pos' column")

    mm$marker <- as.character(mm$marker)
    mm$chr <- as.character(mm$chr)
    mm$pos <- as.integer(mm$pos)

    idx <- match(markers, mm$marker)
    if (any(is.na(idx))) {
      missing_markers <- markers[is.na(idx)]
      stop(sprintf("marker_map_tsv is missing %d markers (e.g. %s)",
                   length(missing_markers), paste(head(missing_markers, 5), collapse = ",")))
    }

    mm2 <- mm[idx, , drop = FALSE]
    # default REF/ALT if not present
    ref <- if ("ref" %in% colnames(mm2)) as.character(mm2$ref) else rep("A", nrow(mm2))
    alt <- if ("alt" %in% colnames(mm2)) as.character(mm2$alt) else rep("C", nrow(mm2))

    vcf_out <- file.path(out_dir, "input_for_beagle.vcf")
    con <- file(vcf_out, open = "wt")
    on.exit({if (!is.null(con)) try(close(con), silent=TRUE)}, add = TRUE)

    writeLines("##fileformat=VCFv4.3", con)
    writeLines("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">", con)

    # contig lines (unique chr)
    for (cc in unique(mm2$chr)) {
      writeLines(sprintf("##contig=<ID=%s>", cc), con)
    }

    header <- c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT", ids)
    writeLines(paste(header, collapse = "\t"), con)

    for (j in seq_along(markers)) {
      gts <- vapply(X0[, j], num_to_gt, character(1))
      # Beagle is strict: number of sample fields in each data line must match header.
      if (length(gts) != length(ids)) {
        cat(sprintf("[WARN] VCF field count mismatch at marker=%s: length(gts)=%d but length(ids)=%d. ",
                    markers[j], length(gts), length(ids)))
        if (length(gts) < length(ids)) {
          cat("Padding missing genotypes with ./..\n")
          gts <- c(gts, rep("./.", length(ids) - length(gts)))
        } else {
          cat("Truncating genotypes to match header.\n")
          gts <- gts[seq_len(length(ids))]
        }
      }
      line <- c(mm2$chr[j], as.character(mm2$pos[j]), markers[j], ref[j], alt[j], ".", "PASS", ".", "GT", gts)
      writeLines(paste(line, collapse = "\t"), con)
    }

    input_vcf <- vcf_out

    # IMPORTANT: close VCF file before running Beagle (avoid partial write being read)
    try(flush(con), silent=TRUE)
    close(con)
    con <- NULL

    cat("[impute_genotypes] VCF written:", input_vcf, "\n")
  }

  out_prefix_beagle <- file.path(out_dir, "beagle_out")
  cmd <- sprintf('java -Xmx%s -jar "%s" gt="%s" out="%s" nthreads=%d burnin=%d iterations=%d phase-states=%d ne=%d em="%s" window=%d overlap=%d seed=%d',
                 java_mem, beagle_jar, input_vcf, out_prefix_beagle, nthreads, burnin, iterations, phase_states, ne, em, window, overlap, seed)
  cat("[impute_genotypes] running:", cmd, "\n")
  status <- system(cmd)
  if (status != 0) stop(sprintf("Beagle failed with status=%d", status))

  out_vcf_gz <- paste0(out_prefix_beagle, ".vcf.gz")
  if (!file.exists(out_vcf_gz)) {
    # sometimes output is .vcf
    out_vcf <- paste0(out_prefix_beagle, ".vcf")
    if (file.exists(out_vcf)) {
      out_vcf_gz <- out_vcf
    } else {
      stop("Beagle output VCF not found")
    }
  }

  # Convert Beagle VCF to TSV (dosage 0/1/2 from GT)
  cat("[impute_genotypes] parsing Beagle output to TSV\n")

  # read VCF line-by-line to avoid huge memory overhead
  con <- if (grepl("\\.gz$", out_vcf_gz)) gzfile(out_vcf_gz, open = "rt") else file(out_vcf_gz, open = "rt")
  on.exit({if (!is.null(con)) try(close(con), silent=TRUE)}, add = TRUE)

  ids <- NULL
  markers <- character(0)
  rows <- list()
  header_done <- FALSE

  while (TRUE) {
    ln <- readLines(con, n = 1)
    if (length(ln) == 0) break
    if (startsWith(ln, "##")) next
    if (startsWith(ln, "#CHROM")) {
      parts <- strsplit(ln, "\t", fixed = TRUE)[[1]]
      ids <- parts[10:length(parts)]
      header_done <- TRUE
      next
    }
    if (!header_done) next

    parts <- strsplit(ln, "\t", fixed = TRUE)[[1]]
    mid <- parts[3]
    markers <- c(markers, mid)

    sample_fields <- parts[10:length(parts)]
    dos <- vapply(sample_fields, parse_gt_to_dosage, numeric(1))
    rows[[length(rows) + 1L]] <- dos
  }

  if (is.null(ids)) stop("Failed to parse VCF header (no sample IDs)")

  Ximp <- do.call(rbind, rows)
  colnames(Ximp) <- ids

  # Transpose to sample x marker
  Ximp2 <- t(Ximp)
  ids2 <- rownames(Ximp2)

  # output TSV: id + markers
  out_tsv_cont <- paste0(out_prefix, ".tsv")
  save_geno_tsv(ids2, markers, Ximp2, out_tsv_cont)

  out_tsv_round <- NULL
  Xround <- NULL
  if (round_to_012) {
    Xround <- round(Ximp2)
    Xround[Xround < 0] <- 0
    Xround[Xround > 2] <- 2
    out_tsv_round <- paste0(out_prefix, "_rounded.tsv")
    save_geno_tsv(ids2, markers, Xround, out_tsv_round)
  }

  marker_missing_before <- if (exists("X0")) colMeans(is.na(X0)) else rep(NA_real_, length(markers))
  sample_missing_before <- if (exists("X0")) rowMeans(is.na(X0)) else rep(NA_real_, length(ids2))

  Xafter <- if (round_to_012 && !is.null(Xround)) Xround else Ximp2
  marker_missing_after <- colMeans(is.na(Xafter))
  sample_missing_after <- rowMeans(is.na(Xafter))

  calc_maf <- function(v) {
    v <- v[!is.na(v)]
    if (length(v) == 0L) return(NA_real_)
    p <- mean(v) / 2.0
    p <- min(max(p, 0.0), 1.0)
    min(p, 1.0 - p)
  }
  maf_after <- apply(Xafter, 2, calc_maf)
  maf_before <- if (exists("X0")) apply(X0, 2, calc_maf) else rep(NA_real_, length(markers))

  marker_sum <- data.frame(
    marker = markers,
    missing_before = marker_missing_before,
    missing_after  = marker_missing_after,
    maf_before     = maf_before,
    maf_after      = maf_after,
    stringsAsFactors = FALSE
  )
  sample_sum <- data.frame(
    id = ids2,
    missing_before = sample_missing_before,
    missing_after  = sample_missing_after,
    stringsAsFactors = FALSE
  )

  fwrite(marker_sum, file.path(out_dir, "marker_summary.tsv"), sep = "\t", quote = F, na = "NA")
  fwrite(sample_sum, file.path(out_dir, "sample_summary.tsv"), sep = "\t", quote = F, na = "NA")

  write_artifacts(list(
    method = method,
    input_mode = if (!is.null(vcf_path) && file.exists(vcf_path)) "vcf" else "tsv",
    genotype_tsv_in = geno_tsv,
    vcf_in = input_vcf,
    beagle_jar = beagle_jar,
    beagle_out_vcf = out_vcf_gz,
    imputed_tsv = out_tsv_cont,
    imputed_tsv_rounded = out_tsv_round,
    nthreads = nthreads,
    quick_accuracy_check = qc
  ))

  cat("[impute_genotypes] done\n")

} else {
  stop(sprintf("Unknown method: %s", method))
}

sink()
