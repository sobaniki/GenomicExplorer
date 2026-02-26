#!/usr/bin/env Rscript

suppressWarnings(suppressMessages({
  if (!requireNamespace("jsonlite", quietly=TRUE)) stop("jsonlite is required")
}))

args <- commandArgs(trailingOnly=TRUE)
get_arg <- function(flag, default=NULL){
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) return(args[[i+1]])
  default
}
params_path <- get_arg("--params")
out_dir <- get_arg("--out")
if (is.null(params_path) || is.null(out_dir)) {
  stop("Usage: Rscript runner.R --params params.json --out out_dir")
}
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

logf <- file.path(out_dir, "run.log")
log <- function(...){
  msg <- paste0("[", Sys.time(), "] ", paste0(..., collapse=""))
  cat(msg, "\n")
  cat(msg, "\n", file=logf, append=TRUE)
}
write_error <- function(msg){
  writeLines(as.character(msg), con=file.path(out_dir, "error_message.txt"), useBytes=TRUE)
}
write_artifacts <- function(meta){
  p <- file.path(out_dir, "artifacts.json")
  writeLines(jsonlite::toJSON(meta, auto_unbox=TRUE, pretty=TRUE), con=p, useBytes=TRUE)
}
safe_png <- function(path, fun){
  tryCatch({
    png(path, width=1200, height=700)
    on.exit(dev.off(), add=TRUE)
    fun()
  }, error=function(e){
    try(dev.off(), silent=TRUE)
  })
}

params <- tryCatch(jsonlite::fromJSON(params_path), error=function(e) list())
log("[poly_qc_counts] start")
log("[poly_qc_counts] params_path=", params_path)
log("[poly_qc_counts] out_dir=", out_dir)

counts_rds <- if(!is.null(params$counts_rds)) as.character(params$counts_rds) else ""
if (!nzchar(counts_rds) || !file.exists(counts_rds)) {
  msg <- paste0("counts_rds not found: ", counts_rds)
  write_error(msg)
  write_artifacts(list(module="poly_qc_counts", error=msg))
  quit(status=1)
}

# thresholds
marker_max_missing <- if(!is.null(params$marker_max_missing)) as.numeric(params$marker_max_missing) else 0.5
sample_max_missing <- if(!is.null(params$sample_max_missing)) as.numeric(params$sample_max_missing) else 0.5
marker_min_mean_depth <- if(!is.null(params$marker_min_mean_depth)) as.numeric(params$marker_min_mean_depth) else 5
sample_min_mean_depth <- if(!is.null(params$sample_min_mean_depth)) as.numeric(params$sample_min_mean_depth) else 5
min_ref_frac <- if(!is.null(params$min_ref_frac)) as.numeric(params$min_ref_frac) else 0.02
max_ref_frac <- if(!is.null(params$max_ref_frac)) as.numeric(params$max_ref_frac) else 0.98

# load
pc <- readRDS(counts_rds)
ref <- as.matrix(pc$ref)
total <- as.matrix(pc$total)
marker <- pc$marker
sample <- pc$sample
meta_in <- pc$meta

if (is.null(marker)) marker <- data.frame(marker=rownames(total), stringsAsFactors=FALSE)
if (is.null(sample)) sample <- data.frame(ind=colnames(total), stringsAsFactors=FALSE)

# compute stats
marker_missing <- rowMeans(is.na(total))
marker_mean_depth <- rowMeans(total, na.rm=TRUE)
marker_ref_frac <- rowMeans(ref / total, na.rm=TRUE)

sample_missing <- colMeans(is.na(total))
sample_mean_depth <- colMeans(total, na.rm=TRUE)
sample_ref_frac <- colMeans(ref / total, na.rm=TRUE)

# filter markers
keep_marker <- rep(TRUE, nrow(total))
keep_marker <- keep_marker & (marker_missing <= marker_max_missing)
keep_marker <- keep_marker & (marker_mean_depth >= marker_min_mean_depth)
keep_marker <- keep_marker & (is.na(marker_ref_frac) | (marker_ref_frac >= min_ref_frac & marker_ref_frac <= max_ref_frac))

# filter samples
keep_sample <- rep(TRUE, ncol(total))
keep_sample <- keep_sample & (sample_missing <= sample_max_missing)
keep_sample <- keep_sample & (sample_mean_depth >= sample_min_mean_depth)

# apply
ref_f <- ref[keep_marker, keep_sample, drop=FALSE]
total_f <- total[keep_marker, keep_sample, drop=FALSE]
marker_f <- marker[match(rownames(total_f), marker$marker, nomatch=0), , drop=FALSE]
if (nrow(marker_f) != nrow(total_f)) marker_f <- data.frame(marker=rownames(total_f), stringsAsFactors=FALSE)
sample_f <- sample[match(colnames(total_f), sample$ind, nomatch=0), , drop=FALSE]
if (nrow(sample_f) != ncol(total_f)) sample_f <- data.frame(ind=colnames(total_f), stringsAsFactors=FALSE)

# summary report
report <- list(
  before=list(n_markers=nrow(total), n_samples=ncol(total)),
  after=list(n_markers=nrow(total_f), n_samples=ncol(total_f)),
  removed=list(
    markers=sum(!keep_marker),
    samples=sum(!keep_sample)
  ),
  thresholds=list(
    marker_max_missing=marker_max_missing,
    sample_max_missing=sample_max_missing,
    marker_min_mean_depth=marker_min_mean_depth,
    sample_min_mean_depth=sample_min_mean_depth,
    min_ref_frac=min_ref_frac,
    max_ref_frac=max_ref_frac
  )
)
writeLines(jsonlite::toJSON(report, auto_unbox=TRUE, pretty=TRUE), con=file.path(out_dir, "qc_report.json"))

rm_markers <- data.frame(marker=rownames(total)[!keep_marker], stringsAsFactors=FALSE)
rm_samples <- data.frame(ind=colnames(total)[!keep_sample], stringsAsFactors=FALSE)
utils::write.table(rm_markers, file.path(out_dir, "removed_markers.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
utils::write.table(rm_samples, file.path(out_dir, "removed_samples.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

# plots before/after
safe_png(file.path(out_dir, "qc_counts_depth.png"), function(){
  par(mfrow=c(2,2))
  hist(marker_mean_depth[is.finite(marker_mean_depth)], main="Marker mean depth (before)", xlab="mean DP", breaks=50)
  hist(rowMeans(total_f, na.rm=TRUE)[is.finite(rowMeans(total_f, na.rm=TRUE))], main="Marker mean depth (after)", xlab="mean DP", breaks=50)
  hist(sample_mean_depth[is.finite(sample_mean_depth)], main="Sample mean depth (before)", xlab="mean DP", breaks=50)
  hist(colMeans(total_f, na.rm=TRUE)[is.finite(colMeans(total_f, na.rm=TRUE))], main="Sample mean depth (after)", xlab="mean DP", breaks=50)
})
safe_png(file.path(out_dir, "qc_counts_missingness.png"), function(){
  par(mfrow=c(2,2))
  hist(marker_missing, main="Marker missingness (before)", xlab="missing rate", breaks=50)
  hist(rowMeans(is.na(total_f)), main="Marker missingness (after)", xlab="missing rate", breaks=50)
  hist(sample_missing, main="Sample missingness (before)", xlab="missing rate", breaks=50)
  hist(colMeans(is.na(total_f)), main="Sample missingness (after)", xlab="missing rate", breaks=50)
})

# save filtered rds
pc_out <- list(
  ref=ref_f,
  total=total_f,
  marker=marker_f,
  sample=sample_f,
  meta=c(meta_in, list(qc=report, qc_created_at=as.character(Sys.time())))
)
saveRDS(pc_out, file.path(out_dir, "poly_counts_filtered.rds"))

meta <- list(
  module="poly_qc_counts",
  counts_rds_in=basename(counts_rds),
  counts_rds_out="poly_counts_filtered.rds",
  qc_report_json="qc_report.json",
  removed_markers_tsv="removed_markers.tsv",
  removed_samples_tsv="removed_samples.tsv",
  plot_depth="qc_counts_depth.png",
  plot_missingness="qc_counts_missingness.png",
  default_table="qc_report.json",
  default_plot="qc_counts_missingness.png"
)
write_artifacts(meta)
log("[poly_qc_counts] done")
