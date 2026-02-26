#!/usr/bin/env Rscript

# QTL mapping via R/qtlbim (Bayesian interval mapping)
#
# Runner contract (GenomicExplorer):
#   Rscript runner.R --params <params.json> --out <out_dir>
#
# Expected params (best-effort; extra keys are ignored):
#   cross_rds   : path to R/qtl cross object saved via saveRDS()
#   trait       : phenotype column name or 1-based index (as string/int)
#   step        : pseudomarker step (cM)
#   error_prob  : genotyping error probability
#   map_function: kosambi|haldane
#   model       : add | add_dom | add_epi
#   n_iter      : number of samples saved
#   burn_in     : number of burn-in iterations
#   thin        : thinning interval
#   seed        : RNG seed
#   covar_tsv   : optional TSV, first column is id, remaining columns covariates
#   scale_covar : TRUE/FALSE (z-score covariates)
#   prior_nqtl  : expected number of main-effect QTL (main.nqtl)
#   extra_options: dict passed to qb.mcmc (e.g., {"genoupdate": true})
#
# Outputs (for GUI):
#   mcmc_summary.tsv
#   bim_profile.tsv
#   peaks.tsv
#   bim_plot.png
#   artifacts.json

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

`%||%` <- function(a, b) {
  if (is.null(a) || length(a) == 0) b else a
}

as_int_safe <- function(x, default = NA_integer_) {
  suppressWarnings({
    v <- as.integer(x)
    if (length(v) == 0 || is.na(v)) default else v
  })
}

as_num_safe <- function(x, default = NA_real_) {
  suppressWarnings({
    v <- as.numeric(x)
    if (length(v) == 0 || is.na(v)) default else v
  })
}

write_tsv <- function(dt, path) {
  if (is.null(dt)) dt <- data.table::data.table()
  dt <- as.data.table(dt)
  # fwrite() warns and exits when there are zero columns; create an empty file instead
  if (ncol(dt) == 0) {
    file.create(path)
    return(invisible(TRUE))
  }
  fwrite(dt, path, sep = "	")
}

write_artifacts <- function(out_dir, artifacts) {
  writeLines(toJSON(artifacts, auto_unbox = TRUE, pretty = TRUE),
             con = file.path(out_dir, "artifacts.json"))
}

params_path <- get_arg("--params")
out_dir <- get_arg("--out")
if (is.null(params_path) || is.null(out_dir)) {
  stop("Usage: runner.R --params params.json --out out_dir")
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- logging ----
log_file <- file.path(out_dir, "run.log")
sink(log_file, split = TRUE)
cat("[qtl_bim] start\n")
cat("[qtl_bim] params_path=", params_path, "\n")
cat("[qtl_bim] out_dir=", out_dir, "\n")



# Robust params reader: prefers JSON, but tolerates TSV (1-row) or a single-cell JSON
read_params <- function(path) {
  # try json first
  x <- tryCatch(jsonlite::fromJSON(path), error = function(e) NULL)
  if (!is.null(x)) return(x)

  # try TSV
  dt <- tryCatch(data.table::fread(path), error = function(e) NULL)
  if (is.null(dt)) stop("Failed to read params: ", path)

  # single-cell JSON
  if (nrow(dt) >= 1 && ncol(dt) == 1) {
    v <- as.character(dt[[1]][1])
    if (grepl("^\\s*\\{", v)) {
      x2 <- tryCatch(jsonlite::fromJSON(v), error = function(e) NULL)
      if (!is.null(x2)) return(x2)
    }
  }

  # assume first row key-value table
  as.list(dt[1])
}

# show a helpful traceback in stderr.txt when R errors
options(error = function() {
  cat("[qtl_bim] ERROR\n")
  traceback(2)
  quit(status = 1)
})
# ---- parse params ----
p <- read_params(params_path)
cat("[qtl_bim] params keys=", paste(names(p), collapse=", "), "\n", sep="")

cross_rds <- as.character(p$cross_rds %||% "")
trait_in <- p$trait %||% ""

step <- as_num_safe(p$step %||% 1.0, 1.0)
if (!is.finite(step) || step <= 0) step <- 1.0

error_prob <- as_num_safe(p$error_prob %||% 0.001, 0.001)
if (!is.finite(error_prob) || error_prob <= 0 || error_prob >= 0.5) error_prob <- 0.001

map_function <- tolower(as.character(p$map_function %||% "kosambi"))
if (!(map_function %in% c("kosambi", "haldane", "c-f", "morgan"))) map_function <- "kosambi"

model <- tolower(as.character(p$model %||% "add"))
if (!(model %in% c("add", "add_epi"))) model <- "add"

n_iter <- as_int_safe(p$n_iter %||% 3000L, 3000L)
burn_in <- as_int_safe(p$burn_in %||% 1000L, 1000L)
thin <- as_int_safe(p$thin %||% 20L, 20L)
seed <- as_int_safe(p$seed %||% 1L, 1L)

covar_tsv <- as.character(p$covar_tsv %||% "")
scale_covar <- isTRUE(p$scale_covar %||% TRUE)

prior_nqtl <- as_int_safe(p$prior_nqtl %||% 3L, 3L)
if (is.na(prior_nqtl) || prior_nqtl < 0) prior_nqtl <- 3L

#extra_options <- p$extra_options
#if (!is.null(extra_options) && !is.list(extra_options)) extra_options <- NULL

# ---- dependency check ----
has_qtl <- requireNamespace("qtl", quietly = TRUE)
has_qtlbim <- requireNamespace("qtlbim", quietly = TRUE)

safe_scalar <- function(x) {
  if (is.null(x) || length(x) == 0) return("")
  as.character(x)[1]
}

summary_keys <- c(
  "status", "message",
  "cross_rds", "trait", "model",
  "step", "error_prob", "map_function",
  "n_iter", "burn_in", "thin", "seed",
  "prior_nqtl",
  "covar_tsv", "scale_covar",
  "has_qtl", "has_qtlbim", "r_version"
)

summary_values <- c(
  "init", "",
  safe_scalar(cross_rds), safe_scalar(trait_in), safe_scalar(model),
  safe_scalar(step), safe_scalar(error_prob), safe_scalar(map_function),
  safe_scalar(n_iter), safe_scalar(burn_in), safe_scalar(thin), safe_scalar(seed),
  safe_scalar(prior_nqtl),
  safe_scalar(covar_tsv), safe_scalar(scale_covar),
  safe_scalar(has_qtl), safe_scalar(has_qtlbim), safe_scalar(R.version.string)
)

summary_dt <- data.table::setDT(data.frame(
  key = c(
    "status", "message",
    "cross_rds", "trait", "model",
    "step", "error_prob", "map_function",
    "n_iter", "burn_in", "thin", "seed",
    "prior_nqtl",
    "covar_tsv", "scale_covar",
    "has_qtl", "has_qtlbim", "r_version"
  ),
  value = c(
    "init", "",
    cross_rds, as.character(trait_in), model,
    step, error_prob, map_function,
    n_iter, burn_in, thin, seed,
    prior_nqtl,
    covar_tsv, scale_covar,
    has_qtl, has_qtlbim, R.version.string
  ),
  stringsAsFactors = FALSE
))

summary_path <- file.path(out_dir, "mcmc_summary.tsv")
profile_path <- file.path(out_dir, "bim_profile.tsv")
peaks_path <- file.path(out_dir, "peaks.tsv")
plot_path <- file.path(out_dir, "bim_plot.png")

# Always create placeholder outputs so GUI doesn't break
write_tsv(summary_dt, summary_path)
write_tsv(data.table::data.table(), profile_path)
write_tsv(data.table::data.table(), peaks_path)

if (!has_qtl || !has_qtlbim) {
  msg <- "Required R packages not installed: need 'qtl' and 'qtlbim'."
  summary_dt[key == "status", value := "error"]
  summary_dt[key == "message", value := msg]
  write_tsv(summary_dt, summary_path)

  write_artifacts(out_dir, list(
    plugin = "qtl_bim",
    main_table = "mcmc_summary.tsv",
    tables = list(
      list(name = "MCMC summary", path = "mcmc_summary.tsv")
    ),
    plots = list()
  ))
  cat("[qtl_bim] ", msg, "\n", sep = "")
  cat("[qtl_bim] done (deps missing)\n")
  quit(status = 0)
}

suppressPackageStartupMessages({
  library(qtl)
  library(qtlbim)
})

if (!nzchar(cross_rds) || !file.exists(cross_rds)) {
  stop("cross_rds not found")
}

cat("[qtl_bim] reading cross\n")
cross <- readRDS(cross_rds)

# ---- trait resolution ----
trait <- trait_in
if (is.numeric(trait_in)) {
  trait <- as.integer(trait_in)
} else {
  # allow numeric-as-string
  tnum <- suppressWarnings(as.integer(trait_in))
  if (!is.na(tnum) && nzchar(as.character(trait_in))) {
    trait <- tnum
  } else {
    trait <- as.character(trait_in)
  }
}

# ---- covariates: merge into cross$pheno and reference by name ----
fixcov_names <- NULL
if (nzchar(covar_tsv) && file.exists(covar_tsv)) {
  cat("[qtl_bim] reading covariates_tsv=", covar_tsv, "\n", sep = "")
  dtc <- tryCatch(fread(covar_tsv), error = function(e) NULL)
  if (!is.null(dtc) && ncol(dtc) >= 2) {
    id <- as.character(dtc[[1]])
    dtc[[1]] <- NULL
    # numeric conversion
    for (j in seq_len(ncol(dtc))) {
      dtc[[j]] <- suppressWarnings(as.numeric(dtc[[j]]))
    }
    # keep numeric columns
    keep <- names(dtc)[vapply(dtc, function(x) is.numeric(x) && any(is.finite(x)), logical(1))]
    if (length(keep) > 0) {
      dtc <- dtc[, ..keep]
      if (scale_covar) {
        dtc <- as.data.table(scale(as.matrix(dtc)))
        setnames(dtc, keep)
      }
      # align to individuals
      ids <- rownames(cross$pheno)
      m <- match(ids, id)
      dtc2 <- dtc[m, , drop = FALSE]
      # prefix names to avoid collisions
      cov_names <- paste0("cov_", keep)
      setnames(dtc2, cov_names)
      cross$pheno <- cbind(cross$pheno, as.data.frame(dtc2))
      fixcov_names <- cov_names
      cat("[qtl_bim] merged covariates: ", paste(fixcov_names, collapse = ","), "\n", sep = "")
    } else {
      cat("[qtl_bim] WARN: no numeric covariate columns found\n")
    }
  } else {
    cat("[qtl_bim] WARN: covariates.tsv read failed or has <2 columns\n")
  }
}

# ---- qtlbim limitation: X chromosome ----
# GUI currently sends "extra_options" (dict). Allow drop_x either top-level or inside extra_options.
# drop_x <- isTRUE((p$drop_x %||% (if (is.list(extra_options)) extra_options$drop_x else NULL)) %||% TRUE)
# if (drop_x) {
#   # subset out chr X if present
#   chr_names <- tryCatch(qtl::chrnames(cross), error = function(e) NULL)
#   if (!is.null(chr_names)) {
#     # common representations: "X" or "x" or "chrX"
#     is_x <- tolower(chr_names) %in% c("x", "chrx")
#     if (any(is_x)) {
#       keep_chr <- chr_names[!is_x]
#       cat("[qtl_bim] dropping X chromosome: ", paste(chr_names[is_x], collapse = ","), "\n", sep = "")
#       cross <- subset(cross, chr = keep_chr)
#     }
#   }
# }

# ---- qb.genoprob ----
cat("[qtl_bim] qb.genoprob step=", step, " error_prob=", error_prob, " map=", map_function, "\n", sep = "")
cross2 <- qb.genoprob(cross = cross, 
                      map.function = map_function,
                      step = step, 
                      tolerance = 1e-6,
                      stepwidth = "variable",
                      error.prob = error_prob)

# ---- qb.data / qb.model ----
trait_dist <- as.character(p$trait_dist %||% "normal")
if (!(trait_dist %in% c("normal", "binary", "ordinal"))) trait_dist <- "normal"

qbData <- if (is.null(fixcov_names)) {
  qb.data(cross = cross2, 
          pheno.col = trait, 
          trait = trait_dist,
          
          censor = NULL, 
          fixcov = c(0), 
          rancov = c(0), 
          boxcox = F,
          standardize = F)
} else {
  qb.data(cross = cross2, 
          pheno.col = trait, 
          trait = trait_dist, 
          fixcov = fixcov_names,
          
          censor = NULL, 
          rancov = c(0), 
          boxcox = F,
          standardize = F)
}

epi <- (model == "add_epi")
qbModel <- qb.model(
  cross = cross2,
  epistasis = epi,
  main.nqtl = prior_nqtl,
  mean.nqtl = as_int_safe(p$mean_nqtl %||% (prior_nqtl + 3L), prior_nqtl + 3L),
  max.nqtl = NULL,
  interval = NULL,
  chr.nqtl = NULL,
  intcov = c(0), 
  depen = F,
  prop = c(0.5, 0.1, 0.05), 
  contrast = T
)

# ---- qb.mcmc ----
set.seed(seed)
cat("[qtl_bim] qb.mcmc n.iter=", n_iter, " n.burnin=", burn_in, " n.thin=", thin, " seed=", seed, "\n", sep = "")

# NOTE: Do NOT use do.call() for qb.mcmc().
# qb.mcmc() internally uses deparse(substitute(cross)); with do.call the first arg is a
# fully-evaluated cross object, so deparse() returns multiple strings and triggers
# 'if (is.transient.cross) ... the condition has length > 1'.

# Optional qb.mcmc formal args we allow to override via extra_options
genoupdate <- TRUE
# if (!is.null(extra_options) && is.list(extra_options) && !is.null(extra_options$genoupdate)) {
#   genoupdate <- isTRUE(extra_options$genoupdate)
# }

qbObj <- qb.mcmc(
  cross = cross2,
  data = qbData,
  model = qbModel,
  mydir = out_dir,
  n.iter = n_iter,
  n.thin = thin,
  n.burnin = burn_in,
  genoupdate = genoupdate,
  seed = seed,
  verbose = T
)
saveRDS(qbObj, file.path(out_dir, "qb_object.rds"))

# ---- qb.scanone ----
type_scan <- as.character(p$type_scan %||% "BF")
#if (!(type_scan %in% c("LPD", "BF", "2logBF", "heritability", "detection"))) type_scan <- "BF"

scan <- c("heritability", 
          "LPD", 
          "LR", 
          "deviance", 
          "detection", 
          "variance",
          "estimate", 
          "cellmean", 
          "count", 
          "log10", 
          "posterior", 
          "logposterior", 
          "BF", 
          "2logBF",
          "nqtl")

cat("[qtl_bim] qb.scanone type.scan=", type_scan, "\n", sep = "")
bim_prof <- qb.scanone(qbObject = qbObj, 
                       epistasis = epi,
                       #scan = scan,
                       type.scan = type_scan, 
                       #covar = ,
                       #adjust.covar = NA,
                       #chr = NULL,
                       sum.scan = "yes",
                       min.iter = 1,
                       aggregate = T,
                       smooth = 3,
                       weight = "sqrt",
                       #split.chr = ,
                       center.type = "scan",
                       half = F,
                       verbose = T)
write_tsv(as.data.table(bim_prof), profile_path)

# summary peaks (max per chromosome; may yield multiple rows per chr if linked loci)
bim_peaks <- tryCatch({
  as.data.table(summary(bim_prof))
}, error = function(e) {
  cat("[qtl_bim] WARN: summary(qb.scanone) failed: ", conditionMessage(e), "\n", sep = "")
  data.table::setDT(data.frame())
})

if (nrow(bim_peaks) > 0) {
  bim_peaks[, trait := as.character(trait_in)]
}
write_tsv(bim_peaks, peaks_path)

# ---- plot ----
tryCatch({
  grDevices::png(plot_path, width = 1200, height = 700)
  plot(bim_prof)
  grDevices::dev.off()
}, error = function(e) {
  cat("[qtl_bim] WARN: plot failed: ", conditionMessage(e), "\n", sep = "")
  try(grDevices::dev.off(), silent = TRUE)
})

# ---- finalize summary + artifacts ----
sm_lines <- capture.output(summary(qbObj))

summary_dt[key == "status", value := "ok"]
summary_dt[key == "message", value := "qtlbim run completed"]
summary_dt <- rbind(
  summary_dt,
  data.table::setDT(data.frame(key = sprintf("qb_summary_line_%03d", seq_along(sm_lines)), value = sm_lines)),
  fill = TRUE
)
write_tsv(summary_dt, summary_path)

write_artifacts(out_dir, list(
  plugin = "qtl_bim",
  default_trait = as.character(trait_in),
  traits = list(as.character(trait_in)),
  main_table = "peaks.tsv",
  trait_files = setNames(list(list(
    peaks = "peaks.tsv",
    lod_profile = "bim_profile.tsv",
    plot_lod = "bim_plot.png"
  )), as.character(trait_in)),
  tables = list(
    list(name = "Peaks", path = "peaks.tsv"),
    list(name = "Profile", path = "bim_profile.tsv"),
    list(name = "MCMC summary", path = "mcmc_summary.tsv")
  ),
  plots = list(
    list(name = "BIM scan", path = "bim_plot.png")
  )
))

cat("[qtl_bim] done\n")
