#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("jsonlite is required")
  if (!requireNamespace("data.table", quietly = TRUE)) stop("data.table is required")
  if (!requireNamespace("GWASpoly", quietly = TRUE)) stop("GWASpoly is required")
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
plot_dir <- file.path(out_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_dir, "run.log")
sink(log_file, split = TRUE)

cat("[poly_gwaspoly_scan] start\n")
cat("params_path=", params_path, "\n")
cat("out_dir=", out_dir, "\n")

p <- tryCatch(jsonlite::fromJSON(params_path), error = function(e) list())

# -----------------------
# Params
# -----------------------
ploidy <- if (!is.null(p$ploidy)) as.integer(p$ploidy) else 4L

dosage_tsv <- if (!is.null(p$dosage_tsv)) as.character(p$dosage_tsv) else ""
marker_info_tsv <- if (!is.null(p$marker_info_tsv)) as.character(p$marker_info_tsv) else ""
pheno_tsv <- if (!is.null(p$pheno_tsv)) as.character(p$pheno_tsv) else ""
trait <- if (!is.null(p$trait)) as.character(p$trait) else "trait1"

models <- if (!is.null(p$models)) p$models else c("additive")
if (is.list(models)) models <- unlist(models)
models <- as.character(models)
models <- models[nzchar(models)]
if (length(models) == 0) models <- c("additive")

fixed_effects <- if (!is.null(p$fixed_effects)) p$fixed_effects else character(0)
if (is.list(fixed_effects)) fixed_effects <- unlist(fixed_effects)
fixed_effects <- as.character(fixed_effects)
fixed_effects <- fixed_effects[nzchar(fixed_effects)]

loco <- if (!is.null(p$loco)) as.logical(p$loco) else TRUE
n_pcs <- if (!is.null(p$n_pcs)) as.integer(p$n_pcs) else 0L
maf <- if (!is.null(p$maf)) as.numeric(p$maf) else NULL
geno_freq <- if (!is.null(p$geno_freq)) as.numeric(p$geno_freq) else NULL

threshold_method <- if (!is.null(p$threshold_method)) as.character(p$threshold_method) else "M.eff"
threshold_level <- if (!is.null(p$threshold_level)) as.numeric(p$threshold_level) else 0.05
n_permute <- if (!is.null(p$n_permute)) as.integer(p$n_permute) else 1000L
bp_window <- if (!is.null(p$bp_window)) as.numeric(p$bp_window) else NA_real_

ncores <- if (!is.null(p$ncores)) as.integer(p$ncores) else 1L
if (is.na(ncores) || ncores < 1) ncores <- 1L

stopifnot(nzchar(dosage_tsv), file.exists(dosage_tsv))
stopifnot(nzchar(marker_info_tsv), file.exists(marker_info_tsv))
stopifnot(nzchar(pheno_tsv), file.exists(pheno_tsv))
stopifnot(nzchar(trait))

cat("ploidy=", ploidy, " trait=", trait, "\n")
cat("dosage_tsv=", dosage_tsv, "\n")
cat("marker_info_tsv=", marker_info_tsv, "\n")
cat("pheno_tsv=", pheno_tsv, "\n")
cat("models=", paste(models, collapse=","), "\n")
cat("fixed_effects=", paste(fixed_effects, collapse=","), "\n")
cat("loco=", loco, " n_pcs=", n_pcs, " maf=", maf, " geno_freq=", geno_freq, "\n")
cat("threshold_method=", threshold_method, " level=", threshold_level, " n_permute=", n_permute, "\n")
cat("bp_window=", bp_window, " ncores=", ncores, "\n")

have_gwaspoly <- requireNamespace("GWASpoly", quietly = TRUE)
if (!have_gwaspoly) {
  msg <- "GWASpoly is not installed. Please install it in your GenomicExplorer R environment."
  writeLines(msg, con = file.path(out_dir, "error_message.txt"))
  stop(msg)
}

# Optional parallel backend: GWASpoly has n.core, but some systems need foreach/doParallel.
if (ncores > 1) {
  if (requireNamespace("doParallel", quietly = TRUE) && requireNamespace("foreach", quietly = TRUE)) {
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    on.exit({
      try(parallel::stopCluster(cl), silent=TRUE)
    }, add = TRUE)
  } else {
    cat("[WARN] doParallel/foreach not installed; continuing without explicit parallel backend\n")
  }
}

# -----------------------
# Prepare GWASpoly object
# -----------------------
# Expected formats:
#  - dosage.tsv: first column sample id, remaining columns marker ids (dosage 0..ploidy)
#  - marker_info.tsv: marker, chr, pos (bp)
#  - phenotype.tsv: first column sample id, includes trait and optional fixed_effects

# GWASpoly can read files directly.
# Note: GWASpoly expects specific column names in map; we standardize to marker/chrom/pos if possible.
map_dt <- data.table::fread(marker_info_tsv)
# heuristics for column names
colnames(map_dt) <- sub("^#", "", colnames(map_dt))

# Ensure marker id column is named 'Marker' for GWASpoly
if (!"Marker" %in% names(map_dt)) {
  cand <- intersect(names(map_dt), c("marker", "MarkerName", "SNP", "id", "ID"))
  if (length(cand) >= 1) data.table::setnames(map_dt, cand[1], "Marker")
}
# Ensure chromosome column is named 'Chrom'
if (!"Chrom" %in% names(map_dt)) {
  cand <- intersect(names(map_dt), c("chr", "Chr", "chrom", "chromosome", "Chromosome"))
  if (length(cand) >= 1) data.table::setnames(map_dt, cand[1], "Chrom")
}
# Ensure position column is named 'Position'
if (!"Position" %in% names(map_dt)) {
  cand <- intersect(names(map_dt), c("pos", "Pos", "bp", "BP", "position"))
  if (length(cand) >= 1) data.table::setnames(map_dt, cand[1], "Position")
}

map_tmp <- file.path(out_dir, "_map_for_gwaspoly.tsv")
data.table::fwrite(map_dt, map_tmp, sep = "\t", quote = F, na = "NA")

dosage_dt <- data.frame(data.table::fread(dosage_tsv, 
                                          header = T),
                        row.names = 1)
dosage_dt_t <- t(dosage_dt)
colnames(dosage_dt_t) <- rownames(dosage_dt)

dosage_map_dt <- cbind(map_dt[, 1:3],
                       dosage_dt_t)
dosage_map_tmp <- file.path(out_dir, "_dosage_map_for_gwaspoly.tsv")
data.table::fwrite(dosage_map_dt, dosage_map_tmp, sep = "\t", quote = F, na = "NA")

# phenotype + dosage: we let GWASpoly read directly
# If fixed effects are provided, they should exist in pheno file.

# read GWASpoly object
obj <- GWASpoly::read.GWASpoly(
  ploidy     = ploidy,
  pheno.file = pheno_tsv,
  #geno.file  = dosage_tsv,
  geno.file  = dosage_map_tmp,
  format = "numeric",
  n.traits = 3,
  delim = "\t"
  #map.file   = map_tmp,
)

# Filters (optional)
if (!is.null(maf)) obj@params$MAF <- maf
if (!is.null(geno_freq)) obj@params$geno.freq <- geno_freq

# Parameters
fixed <- GWASpoly::set.params(fixed  = if (length(fixed_effects) > 0) fixed_effects else NULL,
                              fixed.type = NULL,
                              n.PC   = n_pcs,
                              MAF = maf,
                              geno.freq = geno_freq,
                              P3D = T)

# Kinship
obj <- GWASpoly::set.K(data = obj, 
                       K = NULL,
                       n.core = 1,
                       LOCO = loco)

# GWAS
ans <- GWASpoly::GWASpoly(data = obj,
                          models = models,
                          traits = NULL,
                          params = fixed,
                          n.core = ncores,
                          quiet = F)

# -----------------------
# Collect results
# -----------------------
res_list <- list()
for (m in models) {
  tbl <- data.frame(ans@map[, 1:3],
                    effects = ans@effects,
                    pval = ans@scores)
  dt <- data.table::as.data.table(tbl)
  # Standardize columns to chr/pos/pvalue
  # Guess p-value column
  pcol <- intersect(names(dt), c("P.value", "p.value", "p", "P", "pval", "pvalue", "Pvalue"))
  if (length(pcol) == 0) pcol <- names(dt)[grepl("pval", names(dt), ignore.case = TRUE)][1]
  if (is.na(pcol) || !pcol %in% names(dt)) pcol <- names(dt)[ncol(dt)]

  # chrom/pos columns
  ccol <- intersect(names(dt), c("Chrom", "CHR", "chr"))
  if (length(ccol) == 0) ccol <- "Chrom"
  ppos <- intersect(names(dt), c("Position", "POS", "pos"))
  if (length(ppos) == 0) ppos <- "Position"

  out <- data.table::data.table(
    model = m,
    chr = as.character(dt[[ccol[1]]]),
    pos = as.numeric(dt[[ppos[1]]]),
    pvalue = 10 ^ -(as.numeric(dt[[pcol[1]]]))
  )
  res_list[[m]] <- out
}

if (length(res_list) == 0) {
  msg <- "No results were produced by GWASpoly. Check input formats and GWASpoly version."
  writeLines(msg, con = file.path(out_dir, "error_message.txt"))
  stop(msg)
}

res_long <- data.table::rbindlist(res_list, use.names = TRUE, fill = TRUE)
# default table = first model
first_model <- names(res_list)[1]
res_def <- res_long[model == first_model]

# Write outputs
results_tsv <- file.path(out_dir, "results.tsv")
results_long_tsv <- file.path(out_dir, "results_long.tsv")
data.table::fwrite(res_def, results_tsv, sep = "\t", quote = F, na = "NA")
data.table::fwrite(res_long, results_long_tsv, sep = "\t", quote = F, na = "NA")

# threshold and QTL candidates (best effort)
thr_tsv <- file.path(out_dir, "threshold.tsv")
qtl_tsv <- file.path(out_dir, "qtl.tsv")

thr_dt <- data.table::data.table(method = threshold_method, level = threshold_level, model = first_model)
try({
  if (threshold_method == "permute") {
    # set.threshold may exist; if so, compute and overwrite level
    if ("set.threshold" %in% getNamespaceExports("GWASpoly")) {
      obj2 <- GWASpoly::set.threshold(obj, method = "permute", n.perm = n_permute, alpha = threshold_level)
      thr_dt$threshold <- as.numeric(obj2@params$threshold)
    }
  }
}, silent = TRUE)

data.table::fwrite(thr_dt, thr_tsv, sep = "\t", quote = F, na = "NA")

# get.QTL if available
try({
  if ("get.QTL" %in% getNamespaceExports("GWASpoly")) {
    bw <- if (is.na(bp_window)) NULL else bp_window
    q <- GWASpoly::get.QTL(obj, model = first_model, method = threshold_method, alpha = threshold_level, bp.window = bw)
    qdt <- data.table::as.data.table(q)
    data.table::fwrite(qdt, qtl_tsv, sep = "\t", quote = F, na = "NA")
  }
}, silent = TRUE)

# Static plots (best effort)
try({
  if ("manhattan.plot" %in% getNamespaceExports("GWASpoly")) {
    png(file.path(plot_dir, "manhattan.png"), width = 1600, height = 900)
    GWASpoly::manhattan.plot(obj, model = first_model)
    dev.off()
  }
}, silent = TRUE)
try({
  if ("qq.plot" %in% getNamespaceExports("GWASpoly")) {
    png(file.path(plot_dir, "qq.png"), width = 1200, height = 900)
    GWASpoly::qq.plot(obj, model = first_model)
    dev.off()
  }
}, silent = TRUE)

# Artifacts manifest (optional)
meta <- list(
  module = "poly_gwaspoly_scan",
  ploidy = ploidy,
  trait = trait,
  models = models,
  default_model = first_model,
  results_tsv = "results.tsv",
  results_long_tsv = "results_long.tsv",
  threshold_tsv = "threshold.tsv",
  qtl_tsv = if (file.exists(qtl_tsv)) "qtl.tsv" else NULL,
  manhattan_png = if (file.exists(file.path(plot_dir, "manhattan.png"))) file.path("plots", "manhattan.png") else NULL,
  qq_png = if (file.exists(file.path(plot_dir, "qq.png"))) file.path("plots", "qq.png") else NULL,
  default_table = "results.tsv"
)
writeLines(jsonlite::toJSON(meta, auto_unbox = TRUE, pretty = TRUE), con = file.path(out_dir, "artifacts.json"))

cat("[poly_gwaspoly_scan] done\n")

sink(NULL)
