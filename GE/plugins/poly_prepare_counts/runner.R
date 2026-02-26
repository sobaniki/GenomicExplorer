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

# ---------- helpers ----------
have_pkg <- function(pkg) requireNamespace(pkg, quietly=TRUE)

read_tsv_matrix <- function(path, id_col_guess=c("marker","id","ID","variant","rsid","snp","SNP"), prefer_dt=TRUE){
  if (!file.exists(path)) stop(paste0("File not found: ", path))
  if (prefer_dt && have_pkg("data.table")) {
    dt <- data.table::fread(path, data.table=FALSE, check.names=FALSE)
  } else {
    dt <- utils::read.table(path, header=TRUE, sep="\t", quote="", comment.char="", check.names=FALSE, stringsAsFactors=FALSE)
  }
  if (nrow(dt) < 1) stop(paste0("Empty TSV: ", path))
  # Determine marker column
  id_col <- NULL
  for (nm in id_col_guess) {
    if (nm %in% colnames(dt)) { id_col <- nm; break }
  }
  if (is.null(id_col)) {
    # assume first column is marker id
    id_col <- colnames(dt)[1]
  }
  ids <- dt[[id_col]]
  dt[[id_col]] <- NULL
  mat <- as.matrix(dt)
  storage.mode(mat) <- "numeric"
  rownames(mat) <- make.unique(as.character(ids))
  return(mat)
}

# parse AD string like "10,3" or "10,3,0" -> c(ref=10, total=sum)
parse_ad <- function(x){
  if (is.na(x) || !nzchar(x)) return(c(NA_real_, NA_real_))
  sp <- strsplit(x, ",", fixed=TRUE)[[1]]
  v <- suppressWarnings(as.numeric(sp))
  if (length(v) < 1 || all(is.na(v))) return(c(NA_real_, NA_real_))
  ref <- v[1]
  tot <- sum(v, na.rm=TRUE)
  c(ref, tot)
}

# ---------- main ----------
params <- tryCatch(jsonlite::fromJSON(params_path), error=function(e) list())
log("[poly_prepare_counts] start")
log("[poly_prepare_counts] params_path=", params_path)
log("[poly_prepare_counts] out_dir=", out_dir)

ploidy <- if(!is.null(params$ploidy)) as.integer(params$ploidy) else 4L
input_mode <- if(!is.null(params$input_mode)) as.character(params$input_mode) else "vcf"

# common metadata
meta <- list(
  module="poly_prepare_counts",
  ploidy=ploidy,
  input_mode=input_mode,
  outputs=list()
)

# inputs
vcf_path <- if(!is.null(params$vcf_path)) as.character(params$vcf_path) else ""
ad_field <- if(!is.null(params$ad_field)) as.character(params$ad_field) else "AD"
dp_field <- if(!is.null(params$dp_field)) as.character(params$dp_field) else "DP"
use_dp_for_total <- if(!is.null(params$use_dp_for_total)) as.logical(params$use_dp_for_total) else FALSE
drop_multiallelic <- if(!is.null(params$drop_multiallelic)) as.logical(params$drop_multiallelic) else TRUE

ref_tsv <- if(!is.null(params$ref_count_tsv)) as.character(params$ref_count_tsv) else ""
total_tsv <- if(!is.null(params$total_count_tsv)) as.character(params$total_count_tsv) else ""
counts_rds <- if(!is.null(params$counts_rds)) as.character(params$counts_rds) else ""
orientation <- if(!is.null(params$orientation)) as.character(params$orientation) else "markers_rows"  # markers_rows or samples_rows
write_matrix_tsv <- if(!is.null(params$write_matrix_tsv)) as.logical(params$write_matrix_tsv) else FALSE

# outputs paths
counts_rds_out <- file.path(out_dir, "poly_counts.rds")
marker_tsv_out <- file.path(out_dir, "marker_info.tsv")
sample_tsv_out <- file.path(out_dir, "sample_info.tsv")
summary_marker_out <- file.path(out_dir, "counts_marker_summary.tsv")
summary_sample_out <- file.path(out_dir, "counts_sample_summary.tsv")
plot_depth_out <- file.path(out_dir, "counts_depth.png")
plot_missing_out <- file.path(out_dir, "counts_missingness.png")

ref_mat <- NULL
total_mat <- NULL
marker_info <- NULL
sample_info <- NULL

# ---- load/build PolyCounts ----
tryCatch({
  if (input_mode == "counts_rds") {
    if (!nzchar(counts_rds) || !file.exists(counts_rds)) stop("counts_rds not found")
    pc <- readRDS(counts_rds)
    if (is.null(pc$ref) || is.null(pc$total)) stop("counts_rds must contain list(ref=..., total=...)")
    ref_mat <- as.matrix(pc$ref)
    total_mat <- as.matrix(pc$total)
    marker_info <- pc$marker
    sample_info <- pc$sample
    if (is.null(marker_info)) marker_info <- data.frame(marker=rownames(ref_mat), stringsAsFactors=FALSE)
    if (is.null(sample_info)) sample_info <- data.frame(ind=colnames(ref_mat), stringsAsFactors=FALSE)
    meta$input_mode_detail <- "counts_rds"
  } else if (input_mode == "counts_tsv") {
    if (!nzchar(ref_tsv) || !file.exists(ref_tsv)) stop("ref_count_tsv not found")
    if (!nzchar(total_tsv) || !file.exists(total_tsv)) stop("total_count_tsv not found")
    ref_mat <- read_tsv_matrix(ref_tsv)
    total_mat <- read_tsv_matrix(total_tsv)
    if (orientation == "samples_rows") {
      ref_mat <- t(ref_mat)
      total_mat <- t(total_mat)
    }
    # harmonize ids
    common_markers <- intersect(rownames(ref_mat), rownames(total_mat))
    common_samples <- intersect(colnames(ref_mat), colnames(total_mat))
    ref_mat <- ref_mat[common_markers, common_samples, drop=FALSE]
    total_mat <- total_mat[common_markers, common_samples, drop=FALSE]
    marker_info <- data.frame(marker=common_markers, stringsAsFactors=FALSE)
    sample_info <- data.frame(ind=common_samples, stringsAsFactors=FALSE)
    meta$input_mode_detail <- "counts_tsv"
  } else if (input_mode == "vcf") {
    if (!nzchar(vcf_path) || !file.exists(vcf_path)) stop("vcf_path not found")
    if (!have_pkg("vcfR")) stop("R package 'vcfR' is required for VCF input. Install it (install.packages('vcfR')).")
    vcf <- vcfR::read.vcfR(vcf_path, verbose=FALSE)
    fix <- vcf@fix
    # optionally drop multiallelic
    if (drop_multiallelic && any(grepl(",", fix[, "ALT"], fixed=TRUE))) {
      keep <- !grepl(",", fix[, "ALT"], fixed=TRUE)
      vcf <- vcf[keep, ]
      fix <- vcf@fix
      log("[poly_prepare_counts] dropped multiallelic sites: ", sum(!keep))
    }
    # marker ids
    marker_id <- fix[, "ID"]
    marker_id[is.na(marker_id) | marker_id == "."] <- paste0(fix[, "CHROM"], ":", fix[, "POS"])
    marker_id <- make.unique(marker_id)
    # sample names
    samp <- colnames(vcf@gt)[-1]
    # extract AD and DP
    ad <- vcfR::extract.gt(vcf, element=ad_field, as.numeric=FALSE)
    dp <- NULL
    if (dp_field %in% vcfR::extract.gt(vcf, element="FORMAT", as.numeric=FALSE)[,1]) {
      # ignore; FORMAT varies per row; better attempt extract.gt
    }
    dp <- tryCatch(vcfR::extract.gt(vcf, element=dp_field, as.numeric=TRUE), error=function(e) NULL)
    # parse AD to ref and total
    nvar <- nrow(ad)
    nsamp <- ncol(ad)
    ref_mat <- matrix(NA_real_, nrow=nvar, ncol=nsamp)
    total_mat <- matrix(NA_real_, nrow=nvar, ncol=nsamp)
    for (i in seq_len(nvar)) {
      # vectorized parse for row
      parsed <- t(vapply(ad[i, ], parse_ad, FUN.VALUE=c(NA_real_, NA_real_)))
      ref_mat[i, ] <- parsed[,1]
      total_mat[i, ] <- parsed[,2]
    }
    if (use_dp_for_total && !is.null(dp)) {
      # dp may include NAs; keep AD-based where DP missing
      total_mat <- ifelse(!is.na(dp), dp, total_mat)
    }
    rownames(ref_mat) <- marker_id
    rownames(total_mat) <- marker_id
    colnames(ref_mat) <- samp
    colnames(total_mat) <- samp
    marker_info <- data.frame(marker=marker_id, chr=fix[, "CHROM"], pos=as.integer(fix[, "POS"]), ref=fix[, "REF"], alt=fix[, "ALT"], stringsAsFactors=FALSE)
    sample_info <- data.frame(ind=samp, stringsAsFactors=FALSE)
    meta$input_mode_detail <- "vcfR"
  } else {
    stop(paste0("Unknown input_mode: ", input_mode))
  }

  # basic normalization: ensure total >= ref, set invalid to NA
  bad <- which(!is.na(ref_mat) & !is.na(total_mat) & (total_mat < ref_mat), arr.ind=TRUE)
  if (nrow(bad) > 0) {
    log("[poly_prepare_counts] found total<ref in ", nrow(bad), " cells; set to NA")
    ref_mat[bad] <- NA_real_
    total_mat[bad] <- NA_real_
  }
  total_mat[total_mat <= 0] <- NA_real_
  ref_mat[is.na(total_mat)] <- NA_real_

  # summaries
  marker_missing <- rowMeans(is.na(total_mat), na.rm=FALSE)
  marker_mean_dp <- rowMeans(total_mat, na.rm=TRUE)
  marker_ref_frac <- rowMeans(ref_mat / total_mat, na.rm=TRUE)

  sample_missing <- colMeans(is.na(total_mat), na.rm=FALSE)
  sample_mean_dp <- colMeans(total_mat, na.rm=TRUE)
  sample_ref_frac <- colMeans(ref_mat / total_mat, na.rm=TRUE)

  msum <- data.frame(
    marker=rownames(total_mat),
    missing_rate=marker_missing,
    mean_depth=marker_mean_dp,
    mean_ref_frac=marker_ref_frac,
    stringsAsFactors=FALSE
  )
  ssum <- data.frame(
    ind=colnames(total_mat),
    missing_rate=sample_missing,
    mean_depth=sample_mean_dp,
    mean_ref_frac=sample_ref_frac,
    stringsAsFactors=FALSE
  )
  utils::write.table(msum, summary_marker_out, sep="\t", quote=FALSE, row.names=FALSE)
  utils::write.table(ssum, summary_sample_out, sep="\t", quote=FALSE, row.names=FALSE)

  # plots
  safe_png(plot_depth_out, function(){
    par(mfrow=c(1,2))
    hist(marker_mean_dp[is.finite(marker_mean_dp)], main="Marker mean depth", xlab="mean DP", breaks=50)
    hist(sample_mean_dp[is.finite(sample_mean_dp)], main="Sample mean depth", xlab="mean DP", breaks=50)
  })
  safe_png(plot_missing_out, function(){
    par(mfrow=c(1,2))
    hist(marker_missing, main="Marker missingness", xlab="missing rate", breaks=50)
    hist(sample_missing, main="Sample missingness", xlab="missing rate", breaks=50)
  })

  # write marker/sample info
  if (is.null(marker_info)) marker_info <- data.frame(marker=rownames(total_mat), stringsAsFactors=FALSE)
  if (is.null(sample_info)) sample_info <- data.frame(ind=colnames(total_mat), stringsAsFactors=FALSE)
  utils::write.table(marker_info, marker_tsv_out, sep="\t", quote=FALSE, row.names=FALSE)
  utils::write.table(sample_info, sample_tsv_out, sep="\t", quote=FALSE, row.names=FALSE)

  # optional matrix TSVs (can be huge)
  if (isTRUE(write_matrix_tsv)) {
    # write as marker rows with marker column
    ref_out <- file.path(out_dir, "ref_count.tsv")
    tot_out <- file.path(out_dir, "total_count.tsv")
    ref_df <- data.frame(marker=rownames(ref_mat), ref_mat, check.names=FALSE)
    tot_df <- data.frame(marker=rownames(total_mat), total_mat, check.names=FALSE)
    utils::write.table(ref_df, ref_out, sep="\t", quote=FALSE, row.names=FALSE)
    utils::write.table(tot_df, tot_out, sep="\t", quote=FALSE, row.names=FALSE)
    meta$outputs$ref_count_tsv <- basename(ref_out)
    meta$outputs$total_count_tsv <- basename(tot_out)
  }

  # save PolyCounts RDS
  poly_counts <- list(
    ref=ref_mat,
    total=total_mat,
    marker=marker_info,
    sample=sample_info,
    meta=list(
      ploidy=ploidy,
      input_mode=input_mode,
      source=list(vcf=vcf_path, ref_tsv=ref_tsv, total_tsv=total_tsv, counts_rds=counts_rds),
      created_at=as.character(Sys.time())
    )
  )
  saveRDS(poly_counts, counts_rds_out)

  meta$outputs$counts_rds <- "poly_counts.rds"
  meta$outputs$marker_info_tsv <- "marker_info.tsv"
  meta$outputs$sample_info_tsv <- "sample_info.tsv"
  meta$outputs$counts_marker_summary_tsv <- "counts_marker_summary.tsv"
  meta$outputs$counts_sample_summary_tsv <- "counts_sample_summary.tsv"
  meta$outputs$counts_depth_png <- "counts_depth.png"
  meta$outputs$counts_missingness_png <- "counts_missingness.png"
  meta$default_table <- "counts_marker_summary.tsv"
  meta$default_plot <- "counts_missingness.png"

  write_artifacts(meta)
  log("[poly_prepare_counts] wrote poly_counts.rds and summaries")
  log("[poly_prepare_counts] done")
}, error=function(e){
  write_error(conditionMessage(e))
  meta$error <- conditionMessage(e)
  meta$default_table <- NULL
  meta$default_plot <- NULL
  write_artifacts(meta)
  log("[poly_prepare_counts] ERROR: ", conditionMessage(e))
  quit(status=1)
})
