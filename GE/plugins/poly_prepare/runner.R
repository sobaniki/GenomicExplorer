
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
placeholder_plot <- function(path, title="plot"){
  png(path, width=1200, height=700)
  plot.new()
  text(0.5, 0.6, title, cex=2)
  text(0.5, 0.4, "placeholder", cex=1.4)
  dev.off()
}
params <- tryCatch(jsonlite::fromJSON(params_path), error=function(e) list())
log("start")
log("params_path=", params_path)
log("out_dir=", out_dir)


ploidy <- if(!is.null(params$ploidy)) as.integer(params$ploidy) else 4L
input_mode <- if(!is.null(params$input_mode)) as.character(params$input_mode) else "dosage_tsv"
dosage_in <- if(!is.null(params$dosage_tsv)) as.character(params$dosage_tsv) else ""
marker_in <- if(!is.null(params$marker_info_tsv)) as.character(params$marker_info_tsv) else ""
sample_in <- if(!is.null(params$sample_info_tsv)) as.character(params$sample_info_tsv) else ""

dosage_out <- file.path(out_dir, "dosage.tsv")
if (nzchar(dosage_in) && file.exists(dosage_in)) {
  file.copy(dosage_in, dosage_out, overwrite=TRUE)
} else {
  set.seed(1)
  df <- data.frame(ind=paste0("id",1:10))
  for (m in 1:20) df[[paste0("m",m)]] <- sample(c(0:ploidy, NA), 10, replace=TRUE)
  utils::write.table(df, dosage_out, sep="\t", quote=FALSE, row.names=FALSE)
}

marker_out <- file.path(out_dir, "marker_info.tsv")
if (nzchar(marker_in) && file.exists(marker_in)) {
  file.copy(marker_in, marker_out, overwrite=TRUE)
} else {
  mi <- data.frame(marker=paste0("m",1:20), chr=rep(1:4, length.out=20), pos=seq(1,20)*1e6)
  utils::write.table(mi, marker_out, sep="\t", quote=FALSE, row.names=FALSE)
}

sample_out <- file.path(out_dir, "sample_info.tsv")
if (nzchar(sample_in) && file.exists(sample_in)) {
  file.copy(sample_in, sample_out, overwrite=TRUE)
} else {
  si <- data.frame(ind=paste0("id",1:10))
  utils::write.table(si, sample_out, sep="\t", quote=FALSE, row.names=FALSE)
}

qc_png <- file.path(out_dir, "qc_missingness.png")
placeholder_plot(qc_png, "Missingness QC (placeholder)")

meta <- list(
  module="poly_prepare",
  ploidy=ploidy,
  input_mode=input_mode,
  dosage_tsv="dosage.tsv",
  marker_info_tsv="marker_info.tsv",
  sample_info_tsv="sample_info.tsv",
  plot="qc_missingness.png",
  default_table="dosage.tsv"
)
write_artifacts(meta)
log("done")
