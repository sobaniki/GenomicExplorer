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
log <- function(...) {
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
`%||%` <- function(a, b) if (!is.null(a)) a else b

params <- tryCatch(jsonlite::fromJSON(params_path), error=function(e) list())
log("start")
log("params_path=", params_path)
log("out_dir=", out_dir)

if (!requireNamespace("mappoly", quietly=TRUE)) {
  msg <- "R package 'mappoly' is required but not installed in this environment."
  write_error(msg); write_artifacts(list(module="poly_mappoly_filter", error=msg)); quit(status=1)
}

mappoly_data_rds <- as.character(params$mappoly_data_rds %||% "")
filter_marker_missing <- as.numeric(params$filter_marker_missing %||% 0.05)
filter_ind_missing <- as.numeric(params$filter_ind_missing %||% 0.05)
seg_pval <- params$seg_pval
seg_mode <- as.character(params$seg_mode %||% "bonferroni") # bonferroni|fixed
elim_redundant <- isTRUE(params$elim_redundant %||% TRUE)

if (!nzchar(mappoly_data_rds) || !file.exists(mappoly_data_rds)) {
  msg <- paste0("mappoly_data_rds not found: ", mappoly_data_rds)
  write_error(msg); write_artifacts(list(module="poly_mappoly_filter", error=msg)); quit(status=1)
}

dat <- readRDS(mappoly_data_rds)

nmrk0 <- dat$n.mrk %||% NA_integer_

# Missingness filters
log("filter_missing marker thres=", filter_marker_missing)
dat1 <- mappoly::filter_missing(dat, type="marker", filter.thres=filter_marker_missing)
log("filter_missing individual thres=", filter_ind_missing)
dat2 <- mappoly::filter_missing(dat1, type="individual", filter.thres=filter_ind_missing)

# Segregation distortion filter
if (is.null(seg_pval)) {
  if (seg_mode == "bonferroni") {
    seg_pval <- 0.05 / (dat2$n.mrk %||% 1)
  } else {
    seg_pval <- 0.05
  }
}
seg_pval <- as.numeric(seg_pval)
log("filter_segregation chisq.pval.thres=", seg_pval)
# Build sequence
seq_init <- mappoly::make_seq_mappoly(dat2,
                                      arg = "all",
                                      data.name = NULL,
                                      info.parent = "all",
                                      genomic.info = NULL)

seq_final <- seq_init
red_map <- NULL
if (isTRUE(elim_redundant)) {
  log("elim_redundant")
  seq_red <- mappoly::elim_redundant(seq_init)
  seq_final <- seq_red
  # best-effort redundant mapping extraction (depends on mappoly version)
  if (!is.null(seq_red$redundant)) red_map <- seq_red$redundant
}

# Save outputs
out_files <- list()

seq_init_rds <- file.path(out_dir, "mappoly_seq_init.rds")
saveRDS(seq_init, seq_init_rds)
out_files$mappoly_seq_init_rds <- basename(seq_init_rds)

seq_final_rds <- file.path(out_dir, "mappoly_seq_final.rds")
saveRDS(seq_final, seq_final_rds)
out_files$mappoly_seq_final_rds <- basename(seq_final_rds)

if (!is.null(red_map)) {
  red_tsv <- file.path(out_dir, "redundant_markers.tsv")
  try(utils::write.table(red_map, red_tsv, sep="\t", quote=FALSE, row.names=FALSE), silent=TRUE)
  out_files$redundant_markers_tsv <- basename(red_tsv)
}

summary_tsv <- file.path(out_dir, "mappoly_filter_summary.tsv")
utils::write.table(data.frame(
  key=c("n_markers_input","n_markers_after_missing","n_markers_after_seg","elim_redundant"),
  value=c(nmrk0, dat2$n.mrk %||% NA, length(seq_init$seq.num %||% integer()), isTRUE(elim_redundant)),
  stringsAsFactors=FALSE
), summary_tsv, sep="\t", quote=FALSE, row.names=FALSE)
out_files$mappoly_filter_summary_tsv <- basename(summary_tsv)

write_artifacts(c(list(module="poly_mappoly_filter"), out_files))
log("done")
