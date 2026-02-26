
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
read_counts <- if(!is.null(params$read_counts_tsv)) as.character(params$read_counts_tsv) else ""
have_polyrad <- requireNamespace("polyRAD", quietly=TRUE)

dosage_out <- file.path(out_dir, "dosage.tsv")
note <- NULL
if (!have_polyrad) note <- "polyRAD not installed; wrote placeholder dosage."
if (have_polyrad && nzchar(read_counts) && file.exists(read_counts)) {
  note <- "polyRAD installed; runner scaffold (implement calling here)."
}
set.seed(1)
df <- data.frame(ind=paste0("id",1:10))
for (m in 1:20) df[[paste0("m",m)]] <- sample(c(0:ploidy, NA), 10, replace=TRUE)
utils::write.table(df, dosage_out, sep="\t", quote=FALSE, row.names=FALSE)

qc_png <- file.path(out_dir, "posterior_plot.png")
placeholder_plot(qc_png, "Posterior genotypes (placeholder)")

meta <- list(
  module="poly_polyrad_geno",
  ploidy=ploidy,
  read_counts_tsv=basename(read_counts),
  note=note,
  dosage_tsv="dosage.tsv",
  plot="posterior_plot.png",
  default_table="dosage.tsv"
)
write_artifacts(meta)
log("done")
