#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(jsonlite)
  library(GBScleanR)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) return(args[i + 1])
  default
}

params_path <- get_arg("--params")
out_dir <- get_arg("--out")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

logf <- file.path(out_dir, "run.log")
log <- function(...) {
  msg <- paste0("[impute_gbscleanr] ", paste(..., collapse = " "))
  cat(msg, "\n")
  cat(msg, "\n", file = logf, append = TRUE)
}

errf <- file.path(out_dir, "error.txt")
#sink(file.path(out_dir, "stdout.txt"))
#sink(file.path(out_dir, "stderr.txt"), type = "message")

on.exit({
  sink(type = "message")
  sink()
}, add = TRUE)

log("start")
log("params_path=", params_path)
log("out_dir=", out_dir)

params <- tryCatch({
  jsonlite::fromJSON(params_path)
}, error = function(e) {
  writeLines(paste("Failed to read params:", e$message), errf)
  stop(e)
})

vcf_path <- params$vcf_path
if (is.null(vcf_path) || !file.exists(vcf_path)) {
  writeLines("vcf_path is required and must exist (.vcf/.vcf.gz)", errf)
  stop("vcf_path is required")
}

out_prefix <- params$out_prefix
if (is.null(out_prefix) || !nzchar(out_prefix)) {
  out_prefix <- file.path(out_dir, "gbscleanr")
}

parse_num <- function(x, default) {
  if (is.null(x) || !nzchar(as.character(x))) return(default)
  suppressWarnings({
    v <- as.numeric(x)
    if (is.na(v)) default else v
  })
}
parse_int <- function(x, default) {
  if (is.null(x) || !nzchar(as.character(x))) return(default)
  suppressWarnings({
    v <- as.integer(x)
    if (is.na(v)) default else v
  })
}
parse_bool <- function(x, default = FALSE) {
  if (is.null(x)) return(default)
  if (is.logical(x)) return(x)
  s <- tolower(as.character(x))
  if (s %in% c("true", "t", "1", "yes", "y")) return(TRUE)
  if (s %in% c("false", "f", "0", "no", "n")) return(FALSE)
  default
}

parents <- params$parents
parents <- as.integer(unlist(strsplit(parents, ",")))

cross_type <- as.character(params$cross_type)
log(cross_type)

recomb_rate <- parse_num(params$recomb_rate, 0.04)
error_rate <- parse_num(params$error_rate, 0.0025)
call_threshold <- parse_num(params$call_threshold, 0.9)
het_parent <- parse_bool(params$het_parent, FALSE)
optim <- parse_bool(params$optim, TRUE)
iter <- parse_int(params$iter, 2)
n_threads <- parse_int(params$n_threads, 1)
dummy_reads <- parse_int(params$dummy_reads, 5)

log("inputs:", "vcf_path=", vcf_path)
log("out_prefix=", out_prefix)
log("opts:", "recomb_rate=", recomb_rate, "error_rate=", error_rate,
    "call_threshold=", call_threshold, "het_parent=", het_parent,
    "optim=", optim, "iter=", iter, "n_threads=", n_threads,
    "dummy_reads=", dummy_reads)

# Try to run GBScleanR (Bioconductor). If not installed, fail with a helpful message.
if (!requireNamespace("GBScleanR", quietly = TRUE)) {
  msg <- paste0(
    "GBScleanR is not installed in this R environment.\n",
    "Install (Bioconductor):\n",
    "  if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager')\n",
    "  BiocManager::install('GBScleanR')\n",
    "Then rerun.\n"
  )
  writeLines(msg, errf)
  stop(msg)
}

gds_fn <- paste0(out_prefix, ".gds")

log("converting VCF -> GDS:", gds_fn)
tryCatch({
  GBScleanR::gbsrVCF2GDS(vcf_path, gds_fn, force = T)
}, error = function(e) {
  writeLines(paste("gbsrVCF2GDS failed:", e$message), errf)
  stop(e)
})

log("loading GDS")
gds <- tryCatch({
  GBScleanR::loadGDS(gds_fn)
}, error = function(e) {
  writeLines(paste("loadGDS failed:", e$message), errf)
  stop(e)
})

log("countRead")
gds <- tryCatch({
  GBScleanR::countRead(gds)
}, error = function(e) {
  writeLines(paste("countRead failed:", e$message), errf)
  try(GBScleanR::closeGDS(gds), silent = TRUE)
  stop(e)
})

#gds <- setParents(gds, parents = )
gds <- initScheme(gds, mating = cbind(parents))
gds <- addScheme(gds, crosstype = cross_type)

log("estGeno")
gds2 <- tryCatch({
  GBScleanR::estGeno(
    gds,
    recomb_rate = recomb_rate,
    error_rate = error_rate,
    call_threshold = call_threshold,
    het_parent = het_parent,
    optim = optim,
    iter = iter,
    n_threads = n_threads,
    dummy_reads = dummy_reads
  )
}, error = function(e) {
  writeLines(paste("estGeno failed:", e$message), errf)
  try(GBScleanR::closeGDS(gds), silent = TRUE)
  stop(e)
})

log("export corrected genotype matrix (reference allele counts)")
geno <- tryCatch({
  GBScleanR::getGenotype(gds2, node = "cor", parents = FALSE, valid = TRUE)
}, error = function(e) {
  writeLines(paste("getGenotype failed:", e$message), errf)
  try(GBScleanR::closeGDS(gds2), silent = TRUE)
  stop(e)
})

mar <- tryCatch(GBScleanR::getMarID(gds2), error = function(e) NULL)
sam <- tryCatch(GBScleanR::getSamID(gds2), error = function(e) NULL)

if (is.null(mar)) mar <- paste0("m", seq_len(nrow(geno)))
if (is.null(sam)) sam <- paste0("s", seq_len(ncol(geno)))

out_tsv <- paste0(out_prefix, "_imputed.tsv")

df <- data.frame(marker = mar, geno, check.names = FALSE)
colnames(df) <- c("marker", sam)

utils::write.table(df, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
#write.table(geno, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

log("closing GDS")
try(GBScleanR::closeGDS(gds2), silent = TRUE)

artifacts <- list(
  imputed_tsv = out_tsv,
  gds = gds_fn,
  method = "GBScleanR",
  notes = "imputed_tsv is a marker x sample table (0/1/2 = #reference alleles) exported from node='cor'"
)
writeLines(jsonlite::toJSON(artifacts, auto_unbox = TRUE, pretty = TRUE), file.path(out_dir, "artifacts.json"))

log("done")
