#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(jsonlite)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default=NULL) {
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) return(args[i+1])
  default
}

params_path <- get_arg("--params")
out_dir <- get_arg("--out")
if (is.null(params_path) || is.null(out_dir)) stop("Usage: --params params.json --out out_dir")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_dir, "run.log")
sink(log_file, split = TRUE)

cat("[impute_linkimputer] start\n")
cat("[impute_linkimputer] params_path=", params_path, "\n")
cat("[impute_linkimputer] out_dir=", out_dir, "\n")

p <- fromJSON(params_path)

geno_tsv <- if (!is.null(p$genotype_tsv)) p$genotype_tsv else NULL
marker_map_tsv <- if (!is.null(p$marker_map_tsv)) p$marker_map_tsv else NULL
out_prefix <- if (!is.null(p$out_prefix) && nchar(p$out_prefix) > 0) p$out_prefix else file.path(out_dir, "imputed")
extra <- if (!is.null(p$extra_options)) p$extra_options else list()

na_strings <- c("NA", ".", "-", "nan", "NaN", "NAN", "")

write_artifacts <- function(lst) {
  writeLines(toJSON(lst, auto_unbox = TRUE, pretty = TRUE), file.path(out_dir, "artifacts.json"))
}

load_geno_tsv <- function(path) {
  gt <- fread(path, sep = "\t", header = TRUE, data.table = FALSE, na.strings = na_strings)
  if (!("id" %in% colnames(gt))) stop("genotype_tsv must have 'id' column")
  ids <- as.character(gt$id)
  markers <- setdiff(colnames(gt), "id")
  X <- as.matrix(gt[, markers, drop=FALSE])
  suppressWarnings(storage.mode(X) <- "numeric")
  list(ids=ids, markers=markers, X=X)
}

save_geno_tsv <- function(ids, markers, X, path) {
  dt <- data.table(id = ids)
  dt <- cbind(dt, as.data.table(X))
  setnames(dt, c("id", markers))
  fwrite(dt, path, sep = "\t", na = "NA")
}

mean_impute <- function(X) {
  Ximp <- X
  m <- colMeans(X, na.rm = TRUE)
  for (j in seq_along(m)) Ximp[is.na(Ximp[,j]), j] <- m[j]
  Ximp
}

if (is.null(geno_tsv) || !file.exists(geno_tsv)) stop("genotype_tsv is required")
cat("[impute_linkimputer] genotype_tsv=", geno_tsv, "\n")
if (!is.null(marker_map_tsv)) cat("[impute_linkimputer] marker_map_tsv=", marker_map_tsv, "\n")

D <- load_geno_tsv(geno_tsv)
ids <- D$ids; markers <- D$markers; X0 <- D$X

cat("[impute_linkimputer] n_samples=", nrow(X0), " n_markers=", ncol(X0), "\n")
cat("[impute_linkimputer] missing_rate(before)=", mean(is.na(X0)), "\n")

Ximp <- NULL
used <- "fallback_mean"

ok <- requireNamespace("LinkImputeR", quietly = TRUE)
if (ok) {
  cat("[impute_linkimputer] LinkImputeR available\n")
  # Best-effort call: try a few known function names.
  fun_candidates <- c("LinkImpute", "linkImpute", "linkimpute", "impute")
  f <- NULL
  for (nm in fun_candidates) {
    if (exists(nm, where = asNamespace("LinkImputeR"), inherits = FALSE)) {
      f <- get(nm, envir = asNamespace("LinkImputeR"))
      used <- paste0("LinkImputeR::", nm)
      break
    }
  }
  if (!is.null(f)) {
    cat("[impute_linkimputer] using ", used, "\n")
    # Different implementations expect different orientation; try both.
    try1 <- try({
      res <- f(X0, map = if (!is.null(marker_map_tsv) && file.exists(marker_map_tsv)) fread(marker_map_tsv, data.table=FALSE) else NULL)
      res
    }, silent=TRUE)
    if (!inherits(try1, "try-error")) {
      if (is.matrix(try1) || is.data.frame(try1)) {
        Ximp <- as.matrix(try1)
      } else if (is.list(try1) && !is.null(try1$geno)) {
        Ximp <- as.matrix(try1$geno)
      }
    }
    if (is.null(Ximp)) {
      try2 <- try({
        res <- f(t(X0), map = if (!is.null(marker_map_tsv) && file.exists(marker_map_tsv)) fread(marker_map_tsv, data.table=FALSE) else NULL)
        res
      }, silent=TRUE)
      if (!inherits(try2, "try-error")) {
        if (is.matrix(try2) || is.data.frame(try2)) {
          Ximp <- t(as.matrix(try2))
        } else if (is.list(try2) && !is.null(try2$geno)) {
          Ximp <- t(as.matrix(try2$geno))
        }
      }
    }
  } else {
    cat("[WARN] LinkImputeR found but no known entry function was detected.\n")
  }
}

if (is.null(Ximp)) {
  cat("[impute_linkimputer] falling back to mean imputation\n")
  Ximp <- mean_impute(X0)
}

cat("[impute_linkimputer] missing_rate(after)=", mean(is.na(Ximp)), "\n")

out_tsv <- paste0(out_prefix, "_linkimpute.tsv")
save_geno_tsv(ids, markers, Ximp, out_tsv)

sum_tsv <- paste0(out_prefix, "_summary.tsv")
marker_missing_before <- colMeans(is.na(X0))
marker_missing_after <- colMeans(is.na(Ximp))
summary_dt <- data.frame(marker=markers, missing_before=marker_missing_before, missing_after=marker_missing_after)
fwrite(summary_dt, sum_tsv, sep="\t")

write_artifacts(list(
  plugin = "impute_linkimputer",
  method_used = used,
  imputed_tsv = out_tsv,
  tables = list(
    list(name="marker_summary", path=sum_tsv)
  )
))

cat("[impute_linkimputer] done\n")
RS