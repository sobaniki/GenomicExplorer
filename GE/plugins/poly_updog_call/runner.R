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
if (is.null(params_path) || is.null(out_dir)) stop("Usage: Rscript runner.R --params params.json --out out_dir")
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
log("[poly_updog_call] start")
log("[poly_updog_call] params_path=", params_path)
log("[poly_updog_call] out_dir=", out_dir)

ploidy <- if(!is.null(params$ploidy)) as.integer(params$ploidy) else 4L
counts_rds <- if(!is.null(params$counts_rds)) as.character(params$counts_rds) else ""
marker_file <- if(!is.null(params$marker_file)) as.character(params$marker_file) else ""
method <- if(!is.null(params$method)) as.character(params$method) else "multidog"
ncores <- if(!is.null(params$ncores)) as.integer(params$ncores) else 1L
subset_markers <- if(!is.null(params$subset_markers)) as.character(params$subset_markers) else NULL

# --- NEW: guaranteed posterior probabilities ---
#output_prob <- if(!is.null(params$output_prob)) as.logical(params$output_prob) else TRUE
output_prob <- if(!is.null(params$output_prob)) as.logical(params$output_prob) else FALSE
call_prob_threshold <- if(!is.null(params$call_prob_threshold)) as.numeric(params$call_prob_threshold) else 0.9
prob_format <- if(!is.null(params$prob_format)) as.character(params$prob_format) else "long"  # long|wide|none

# parent/progeny options (optional; for f1/s1 models)
parent1_id <- if(!is.null(params$parent1_id)) as.character(params$parent1_id) else ""
parent2_id <- if(!is.null(params$parent2_id)) as.character(params$parent2_id) else ""
progeny_ids_raw <- if(!is.null(params$progeny_ids)) as.character(params$progeny_ids) else ""
progeny_file <- if(!is.null(params$progeny_file)) as.character(params$progeny_file) else ""

parse_id_list <- function(x){
  x <- as.character(x)
  x <- gsub("[\r\n\t]", " ", x)
  parts <- unlist(strsplit(x, "[,; ]+"))
  parts <- trimws(parts)
  parts[nzchar(parts)]
}
extra_options <- params$extra_options
if (is.null(extra_options)) extra_options <- list()

if (!nzchar(counts_rds) || !file.exists(counts_rds)) {
  msg <- paste0("counts_rds not found: ", counts_rds)
  write_error(msg)
  write_artifacts(list(module="poly_updog_call", error=msg))
  quit(status=1)
}

have_updog <- requireNamespace("updog", quietly=TRUE)
if (!have_updog) {
  msg <- "R package 'updog' is not installed. Install from CRAN/Bioconductor as appropriate."
  write_error(msg)
  write_artifacts(list(module="poly_updog_call", error=msg))
  quit(status=1)
}

pc <- readRDS(counts_rds)
ref <- as.matrix(pc$ref)
total <- as.matrix(pc$total)
marker_info <- pc$marker
sample_info <- pc$sample

if (!is.null(subset_markers) && length(subset_markers) > 0) {
  keep <- intersect(rownames(total), subset_markers)
  ref <- ref[keep, , drop=FALSE]
  total <- total[keep, , drop=FALSE]
  if (!is.null(marker_info) && "marker" %in% colnames(marker_info)) {
    marker_info <- marker_info[marker_info$marker %in% keep, , drop=FALSE]
  }
  log("[poly_updog_call] subset markers: ", nrow(total))
}

# Ensure numeric and consistent
storage.mode(ref) <- "numeric"
storage.mode(total) <- "numeric"
total[total <= 0] <- NA_real_
ref[is.na(total)] <- NA_real_

markers <- rownames(total)
samples <- colnames(total)

# Resolve parents/progeny
parents <- character(0)
if (nzchar(parent1_id)) parents <- c(parents, parent1_id)
if (nzchar(parent2_id)) parents <- c(parents, parent2_id)
parents <- unique(parents)

missing_parents <- setdiff(parents, samples)
if (length(missing_parents) > 0) {
  log("[poly_updog_call] WARNING: parent IDs not found in counts_rds columns: ", paste(missing_parents, collapse=","))
}

progeny_samples <- NULL
if (nzchar(progeny_file) && file.exists(progeny_file)) {
  ids <- readLines(progeny_file, warn=FALSE)
  ids <- trimws(ids)
  ids <- ids[nzchar(ids)]
  progeny_samples <- intersect(ids, samples)
  log("[poly_updog_call] progeny from file: n=", length(progeny_samples))
} else if (nzchar(progeny_ids_raw)) {
  ids <- parse_id_list(progeny_ids_raw)
  progeny_samples <- intersect(ids, samples)
  log("[poly_updog_call] progeny from list: n=", length(progeny_samples))
} else if (length(parents) > 0) {
  progeny_samples <- setdiff(samples, parents)
  log("[poly_updog_call] progeny=all non-parents: n=", length(progeny_samples))
} else {
  progeny_samples <- samples
  log("[poly_updog_call] progeny=all samples: n=", length(progeny_samples))
}

if (length(progeny_samples) == 0) {
  log("[poly_updog_call] WARNING: no progeny samples selected; falling back to all samples")
  progeny_samples <- samples
}

# Parents present in counts (kept for downstream linkage mapping).
parents_present <- parents[parents %in% samples]
sample_order_all <- unique(c(parents_present, progeny_samples))


# Parent count vectors (per marker)
p1ref <- p1size <- p2ref <- p2size <- NULL
if (nzchar(parent1_id) && parent1_id %in% samples) {
  p1ref <- ref[, parent1_id]
  p1size <- total[, parent1_id]
}
if (nzchar(parent2_id) && parent2_id %in% samples) {
  p2ref <- ref[, parent2_id]
  p2size <- total[, parent2_id]
}

# Subset to progeny matrices for calling (parents are provided separately via p1*/p2*)
ref_call <- ref[, progeny_samples, drop=FALSE]
total_call <- total[, progeny_samples, drop=FALSE]

# Keep sample_info in sync (if available)
if (!is.null(sample_info)) {
  if ("sample" %in% colnames(sample_info)) {
    sample_info <- sample_info[sample_info$sample %in% sample_order_all, , drop=FALSE]
  } else if ("ind" %in% colnames(sample_info)) {
    sample_info <- sample_info[sample_info$ind %in% sample_order_all, , drop=FALSE]
  }
}

# helper: extract genotype from multidog result
extract_geno <- function(res){
  if (is.null(res)) return(NULL)
  for (nm in c("geno", "genotype", "call", "calls", "z", "dosage")) {
    if (!is.null(res[[nm]])) return(res[[nm]])
  }
  if (!is.null(res$fullout) && !is.null(res$fullout$geno)) return(res$fullout$geno)
  NULL
}

# helper: extract posterior probabilities from flexdog result
extract_postprob <- function(fd, ploidy){
  cand <- NULL
  for (nm in c("postprob", "postprobs", "posterior", "posteriors", "prob", "probs")) {
    if (!is.null(fd[[nm]])) { cand <- fd[[nm]]; break }
  }
  if (is.null(cand)) return(NULL)
  m <- as.matrix(cand)
  if (ncol(m) != (ploidy + 1)) return(NULL)
  if (is.null(colnames(m))) colnames(m) <- as.character(0:ploidy)
  suppressWarnings({
    if (all(is.finite(as.numeric(colnames(m))))) {
      ord <- order(as.numeric(colnames(m)))
      m <- m[, ord, drop=FALSE]
    }
  })
  m
}

# enforce probability output: use flexdog for guaranteed posterior probs
if (isTRUE(output_prob) && method != "flexdog") {
  log("[poly_updog_call] output_prob=TRUE -> force method=flexdog (guaranteed posterior probs)")
  method <- "flexdog"
}

markers_n <- length(markers)
samples_n <- length(progeny_samples)
dosage_levels <- as.character(0:ploidy)

samples_all_n <- length(sample_order_all)


# dosage_mat <- matrix(NA_integer_, nrow=markers_n, ncol=samples_n,
#                      dimnames=list(markers, progeny_samples))
fit_summary <- data.frame(marker=markers, converged=NA, stringsAsFactors=FALSE)

# probability array (marker x sample x (P+1))
# prob_arr <- array(NA_real_,
#                   dim=c(markers_n, samples_n, ploidy+1),
#                   dimnames=list(markers, progeny_samples, dosage_levels))
# 
# # Full containers (parents + progeny). Parents are appended before progeny.
# prob_arr_all <- array(NA_real_,
#                       dim=c(markers_n, samples_all_n, ploidy+1),
#                       dimnames=list(markers, sample_order_all, dosage_levels))
# dosage_mat_all <- matrix(NA_integer_, nrow=markers_n, ncol=samples_all_n,
#                          dimnames=list(markers, sample_order_all))

failed <- data.frame(marker_id=character(), reason=character(), stringsAsFactors=FALSE)

log("[poly_updog_call] using updog::multidog")

mout <- updog::multidog(refmat = ref,
                        sizemat = total,
                        ploidy = ploidy,
                        model = as.character(params$model),
                        nc = params$ncores,
                        p1_id = params$parent1_id,
                        p2_id = params$parent2_id,
                        seq = params$seq,
                        bias = params$bias,
                        od = params$od,
                        update_seq = params$update_seq,
                        update_bias = params$update_bias,
                        update_od = params$update_od)

dosage_mat_all <- updog::format_multidog(x = mout,
                                         varname = "geno")
if (!is.null(mout$snpdf$p1geno)) {
  p1geno <- mout$snpdf$p1geno
  p2geno <- mout$snpdf$p2geno
  dosage_mat_all <- cbind(p1geno,
                          p2geno,
                          dosage_mat_all)
  colnames(dosage_mat_all)[1:2] <- c(params$parent1_id,
                                     params$parent2_id)
} else if (!is.null(mout$snpdf$pgeno)) {
  pgeno <- mout$snpdf$pgeno
  dosage_mat_all <- cbind(pgeno,
                          dosage_mat_all)
  colnames(dosage_mat_all)[1] <- params$parent1_id
}
#dosage_mat <- dosage_mat_all[, progeny_samples]

pr <- colnames(mout$inddf)[grep("^Pr_",
                                colnames(mout$inddf))]
prob_arr <- c()
for (cyc1 in 1:length(pr)) {
  obj <- updog::format_multidog(x = mout,
                                varname = pr[cyc1])
  prob_arr <- c(prob_arr,
                list(obj))
}
prob_arr <- as.array(prob_arr)
#prob_arr_all <- updog::format_multidog(x = mout,
#                                       varname = "maxpostprob")
#prob_arr <- prob_arr_all[, progeny_samples]

fit_summary$converged <- TRUE

# derived matrices (posterior-based summaries)
dvec <- 0:ploidy

# # progeny-only (legacy)
# dosage_pp <- apply(prob_arr, c(1,2), function(p){
#   if (all(is.na(p))) return(NA_real_)
#   sum(p * dvec, na.rm=TRUE)
# })
# dosage_maxprob <- apply(prob_arr, c(1,2), function(p){
#   if (all(is.na(p))) return(NA_real_)
#   max(p, na.rm=TRUE)
# })
dosage_maxprob <- updog::format_multidog(x = mout,
                                         varname = "maxpostprob")
# dosage_entropy <- apply(prob_arr, c(1,2), function(p){
#   if (all(is.na(p))) return(NA_real_)
#   p <- p[is.finite(p) & p > 0]
#   -sum(p * log(p))
# })

# all samples (parents + progeny)
# dosage_pp_all <- apply(prob_arr_all, c(1,2), function(p){
#   if (all(is.na(p))) return(NA_real_)
#   sum(p * dvec, na.rm=TRUE)
# })
# dosage_maxprob_all <- apply(prob_arr_all, c(1,2), function(p){
#   if (all(is.na(p))) return(NA_real_)
#   max(p, na.rm=TRUE)
# })
# dosage_entropy_all <- apply(prob_arr_all, c(1,2), function(p){
#   if (all(is.na(p))) return(NA_real_)
#   p <- p[is.finite(p) & p > 0]
#   -sum(p * log(p))
# })

# write dosage outputs
# utils::write.table(
#   data.frame(ind=progeny_samples, t(dosage_mat), check.names=FALSE),
#   file.path(out_dir, "dosage.tsv"),
#   sep="\t", quote=FALSE, row.names=FALSE
# )
if (marker_file != "") {
  marker <- data.table::fread(marker_file,
                              header = T,
                              data.table = F)
  marker <- marker[marker[, 1] %in% rownames(dosage_mat_all), ]
}

utils::write.table(
  data.frame(ind=sample_order_all, t(dosage_mat_all), check.names=FALSE),
  file.path(out_dir, "dosage_all.tsv"),
  sep="\t", quote=FALSE, row.names=FALSE
)
if (params$model %in% c("f1", "f1pp") & marker_file != "") {
  utils::write.table(
    data.frame(marker = rownames(dosage_mat_all),
               p1 = dosage_mat_all[, params$parent1_id],
               p2 = dosage_mat_all[, params$parent2_id],
               chr = marker[, 2],
               pos = marker[, 3],
               dosage_mat_all[, progeny_samples]),
    file.path(out_dir, "dosage_all.mappoly.csv"),
    sep=",", quote=FALSE, row.names=FALSE
  )
} else if (params$model %in% c("s1", "s1pp") & marker_file != "") {
  utils::write.table(
    data.frame(marker = rownames(dosage_mat_all),
               p1 = dosage_mat_all[, params$parent1_id],
               p2 = dosage_mat_all[, params$parent1_id],
               chr = marker[, 2],
               pos = marker[, 3],
               dosage_mat_all[, progeny_samples]),
    file.path(out_dir, "dosage_all.mappoly.csv"),
    sep=",", quote=FALSE, row.names=FALSE
  )
}
utils::write.table(
  data.frame(sample_id=progeny_samples),
  file.path(out_dir, "progeny_ids.tsv"),
  sep="\t", quote=FALSE, row.names=FALSE
)

# QC tables (based on all samples)
#called_all <- dosage_maxprob_all >= call_prob_threshold
marker_stats <- data.frame(
  marker_id = markers,
  #call_rate = rowMeans(called_all, na.rm=TRUE),
  mean_maxprob = rowMeans(dosage_maxprob, na.rm=TRUE),
  #mean_entropy = rowMeans(dosage_entropy_all, na.rm=TRUE),
  n_nonmissing = rowSums(!is.na(dosage_mat_all)),
  stringsAsFactors=FALSE
)
sample_stats <- data.frame(
  sample_id = sample_order_all,
  #call_rate = colMeans(called_all, na.rm=TRUE),
  #mean_maxprob = colMeans(dosage_maxprob, na.rm=TRUE),
  #mean_entropy = colMeans(dosage_entropy_all, na.rm=TRUE),
  n_nonmissing = colSums(!is.na(dosage_mat_all)),
  stringsAsFactors=FALSE
)

utils::write.table(marker_stats, file.path(out_dir, "marker_stats.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
utils::write.table(sample_stats, file.path(out_dir, "sample_stats.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

# legacy marker_callrate.tsv (keep for backward compatibility)
# call_rate_legacy <- 1 - rowMeans(is.na(dosage_mat))
# summ <- data.frame(marker=markers, call_rate=call_rate_legacy, stringsAsFactors=FALSE)
# if (!is.null(marker_info) && "marker" %in% colnames(marker_info)) {
#   summ <- merge(marker_info, summ, by="marker", all.y=TRUE, sort=FALSE)
# }
# utils::write.table(summ, file.path(out_dir, "marker_callrate.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

# plots
safe_png(file.path(out_dir, "dosage_hist.png"), function(){
  x <- as.numeric(dosage_mat_all)
  x <- x[is.finite(x)]
  hist(x, breaks=seq(-0.5, ploidy+0.5, by=1), main="Dosage histogram", xlab="dosage")
})
safe_png(file.path(out_dir, "maxprob_hist.png"), function(){
  x <- as.numeric(dosage_maxprob)
  x <- x[is.finite(x)]
  hist(x, breaks=30, main=paste0("Max posterior prob (threshold=", call_prob_threshold, ")"), xlab="maxprob")
  abline(v=call_prob_threshold, lty=2)
})

# prob TSV output (optional)
if (isTRUE(output_prob) && prob_format == "long") {
  gz <- gzfile(file.path(out_dir, "prob_long.tsv.gz"), "wt")
  on.exit(close(gz), add=TRUE)
  writeLines("marker_id\tsample_id\td\tp", gz)
  for (mi in seq_along(markers)) {
    m <- markers[mi]
    pp <- prob_arr[mi,,]
    for (si in seq_along(progeny_samples)) {
      s <- progeny_samples[si]
      pv <- pp[si,]
      if (all(is.na(pv))) next
      for (di in seq_along(dosage_levels)) {
        writeLines(paste(m, s, dosage_levels[di], format(pv[di], scientific=FALSE), sep="\t"), gz)
      }
    }
  }
  
  # also write parents+progeny (useful for linkage mapping)
  if (length(parents_present) > 0) {
    gz2 <- gzfile(file.path(out_dir, "prob_long_all.tsv.gz"), "wt")
    on.exit(close(gz2), add=TRUE)
    writeLines("marker_id\tsample_id\td\tp", gz2)
    for (mi in seq_along(markers)) {
      m <- markers[mi]
      pp <- prob_arr_all[mi,,]
      for (si in seq_along(sample_order_all)) {
        s <- sample_order_all[si]
        pv <- pp[si,]
        if (all(is.na(pv))) next
        for (di in seq_along(dosage_levels)) {
          writeLines(paste(m, s, dosage_levels[di], format(pv[di], scientific=FALSE), sep="\t"), gz2)
        }
      }
    }
  }
} else if (isTRUE(output_prob) && prob_format == "wide") {
  # wide format (marker_id, sample_id, p0..pP)
  gz <- gzfile(file.path(out_dir, "prob_wide.tsv.gz"), "wt")
  on.exit(close(gz), add=TRUE)
  hdr <- paste(c("marker_id","sample_id", paste0("p", dosage_levels)), collapse="\t")
  writeLines(hdr, gz)
  for (mi in seq_along(markers)) {
    m <- markers[mi]
    pp <- prob_arr[mi,,]
    for (si in seq_along(progeny_samples)) {
      s <- progeny_samples[si]
      pv <- pp[si,]
      if (all(is.na(pv))) next
      writeLines(paste(c(m, s, format(pv, scientific=FALSE)), collapse="\t"), gz)
    }
  }
  
  # also write parents+progeny (useful for linkage mapping)
  if (length(parents_present) > 0) {
    gz2 <- gzfile(file.path(out_dir, "prob_wide_all.tsv.gz"), "wt")
    on.exit(close(gz2), add=TRUE)
    hdr <- paste(c("marker_id","sample_id", paste0("p", dosage_levels)), collapse="\t")
    writeLines(hdr, gz2)
    for (mi in seq_along(markers)) {
      m <- markers[mi]
      pp <- prob_arr_all[mi,,]
      for (si in seq_along(sample_order_all)) {
        s <- sample_order_all[si]
        pv <- pp[si,]
        if (all(is.na(pv))) next
        writeLines(paste(c(m, s, format(pv, scientific=FALSE)), collapse="\t"), gz2)
      }
    }
  }
}

# save legacy dosage.rds (keep)
poly_dosage <- list(
  dosage=dosage_mat_all,
  marker=marker_info,
  sample=sample_info,
  meta=list(
    ploidy=ploidy,
    method=method,
    counts_rds=counts_rds,
    parent1_id=parent1_id,
    parent2_id=parent2_id,
    n_progeny=length(progeny_samples),
    progeny_samples=progeny_samples,
    extra_options=extra_options,
    output_prob=output_prob,
    call_prob_threshold=call_prob_threshold,
    prob_format=prob_format,
    created_at=as.character(Sys.time())
  )
)
saveRDS(poly_dosage, file.path(out_dir, "dosage.rds"))

# build geno.rds (common schema)
# markers
markers_df <- NULL
if (!is.null(marker_info)) {
  mi <- marker_info
  if (!("marker_id" %in% colnames(mi))) {
    if ("marker" %in% colnames(mi)) mi$marker_id <- as.character(mi$marker)
    else mi$marker_id <- markers
  }
  markers_df <- mi
} else {
  markers_df <- data.frame(marker_id=markers, stringsAsFactors=FALSE)
}

# samples
samples_df <- NULL
if (!is.null(sample_info)) {
  si <- sample_info
  if (!("sample_id" %in% colnames(si))) {
    if ("ind" %in% colnames(si)) si$sample_id <- as.character(si$ind)
    else if ("sample" %in% colnames(si)) si$sample_id <- as.character(si$sample)
    else si$sample_id <- sample_order_all
  }
  samples_df <- si
} else {
  samples_df <- data.frame(sample_id=sample_order_all, stringsAsFactors=FALSE)
}

# Ensure order and add role labels
samples_df <- samples_df[samples_df$sample_id %in% sample_order_all, , drop=FALSE]
samples_df <- samples_df[match(sample_order_all, samples_df$sample_id), , drop=FALSE]
samples_df$role <- ifelse(samples_df$sample_id == parent1_id, "parent1",
                          ifelse(samples_df$sample_id == parent2_id, "parent2", "progeny"))

# ensure matrices dimnames match marker_id/sample_id order
# (dosage_mat already dimnames=markers/progeny_samples)
geno <- list(
  version = "1.0.0",
  created_at = as.character(Sys.time()),
  ploidy = as.integer(ploidy),
  allele_coding = "ref_alt",
  caller = list(method="updog", method_version=NA_character_, params=params, seed=NA_integer_, notes="prob guaranteed via flexdog"),
  samples = samples_df,
  markers = markers_df,
  dosage = dosage_mat_all,
  prob = if (isTRUE(output_prob)) prob_arr_all else array(numeric(0), dim=c(0,0,0)),
  marker_stats = marker_stats,
  sample_stats = sample_stats,
  #dosage_pp = dosage_pp_all,
  dosage_maxprob = dosage_maxprob,
  #dosage_entropy = dosage_entropy_all,
  missing = list(dosage_value=NA_integer_, prob_value=NA_real_, definition="depth0_or_filtered_or_failed"),
  run_log = list(
    n_markers_input = markers_n,
    n_markers_output = markers_n,
    n_samples = samples_all_n,
    failed_markers = failed
  ),
  extra = list(
    counts_ref = NULL,
    counts_total = NULL,
    caller_internal = NULL,
    parents_present = parents_present,
    progeny_samples = progeny_samples
  )
)
saveRDS(geno, file.path(out_dir, "geno.rds"))

meta <- list(
  module="poly_updog_call",
  ploidy=ploidy,
  method=method,
  output_prob=output_prob,
  call_prob_threshold=call_prob_threshold,
  prob_format=prob_format,
  counts_rds_in=basename(counts_rds),
  parent1_id=parent1_id,
  parent2_id=parent2_id,
  n_progeny=length(progeny_samples),
  progeny_file=if (nzchar(progeny_file)) basename(progeny_file) else "",
  dosage_tsv="dosage.tsv",
  dosage_all_tsv="dosage_all.tsv",
  progeny_ids_tsv="progeny_ids.tsv",
  dosage_rds="dosage.rds",
  geno_rds="geno.rds",
  marker_stats_tsv="marker_stats.tsv",
  sample_stats_tsv="sample_stats.tsv",
  marker_callrate_tsv="marker_callrate.tsv",
  plot="dosage_hist.png",
  plot_maxprob="maxprob_hist.png",
  prob_long_tsv_gz=if (isTRUE(output_prob) && prob_format=="long") "prob_long.tsv.gz" else "",
  prob_long_all_tsv_gz=if (isTRUE(output_prob) && prob_format=="long" && length(parents_present)>0) "prob_long_all.tsv.gz" else "",
  prob_wide_tsv_gz=if (isTRUE(output_prob) && prob_format=="wide") "prob_wide.tsv.gz" else "",
  prob_wide_all_tsv_gz=if (isTRUE(output_prob) && prob_format=="wide" && length(parents_present)>0) "prob_wide_all.tsv.gz" else "",
  default_table="dosage.tsv",
  default_plot="dosage_hist.png"
)
write_artifacts(meta)
log("[poly_updog_call] done")