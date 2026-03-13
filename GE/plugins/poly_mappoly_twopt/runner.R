suppressWarnings(suppressMessages({
  library(mappoly)
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
  write_error(msg); write_artifacts(list(module="poly_mappoly_twopt", error=msg)); quit(status=1)
}

seq_rds <- as.character(params$seq_rds %||% "")
ncpus <- as.integer(params$ncpus %||% 1)
verbose <- isTRUE(params$verbose %||% TRUE)
expected_groups <- as.integer(params$expected_groups %||% 1)

if (!nzchar(seq_rds) || !file.exists(seq_rds)) {
  msg <- paste0("seq_rds not found: ", seq_rds)
  write_error(msg); write_artifacts(list(module="poly_mappoly_twopt", error=msg)); quit(status=1)
}

seq_init <- readRDS(seq_rds)
log("est_pairwise_rf2 ncpus=", ncpus)

mappoly_data_rds <- as.character(params$mappoly_data_rds %||% "")
dat2 <- readRDS(mappoly_data_rds)

tpt <- mappoly::est_pairwise_rf2(input.seq = seq_init$unique.seq,
                                 ncpus = ncpus, 
                                 mrk.pairs = NULL,
                                 verbose = verbose,
                                 tol = .Machine$double.eps ^ 0.25)
rfmat <- mappoly::rf_list_to_matrix(input.twopt = tpt,
                                    thresh.LOD.ph = 0,
                                    thresh.LOD.rf = 0,
                                    thresh.rf = 0.5,
                                    ncpus = ncpus,
                                    shared.alleles = F,
                                    verbose = verbose)

# save
out_files <- list()

tpt_rds <- file.path(out_dir, "mappoly_twopt.rds")
saveRDS(tpt, tpt_rds)
out_files$mappoly_twopt_rds <- basename(tpt_rds)

rf_rds <- file.path(out_dir, "mappoly_rf_matrix.rds")
saveRDS(rfmat, rf_rds)
out_files$mappoly_rf_matrix_rds <- basename(rf_rds)

# optional simple plot
png(file.path(out_dir, "rf_matrix.png"), width=1200, height=1000)
try(plot(rfmat, fact=2), silent=TRUE)
dev.off()
out_files$rf_matrix_png <- "rf_matrix.png"

summary_tsv <- file.path(out_dir, "mappoly_twopt_summary.tsv")
infer_n_pairs <- function(tpt){
  # Different MAPpoly versions store this differently.
  if (is.null(tpt)) return(NA_integer_)
  if (!is.null(tpt$pairwise) && is.data.frame(tpt$pairwise)) return(nrow(tpt$pairwise))
  if (!is.null(tpt$mrk.pairs) && is.matrix(tpt$mrk.pairs)) return(nrow(tpt$mrk.pairs) %/% 2)
  if (is.list(tpt)) {
    # often it's a list of length N with per-pair entries
    if (!is.null(tpt[[1]]) && !is.list(tpt[[1]])) {
      # unlikely
      return(length(tpt))
    }
    if (!is.null(names(tpt))) return(length(tpt))
  }
  NA_integer_
}
utils::write.table(data.frame(
  key=c("n_markers","n_pairs"),
  value=c(length(seq_init$seq.num %||% integer()), infer_n_pairs(tpt)),
  stringsAsFactors=FALSE
), summary_tsv, sep="\t", quote=FALSE, row.names=FALSE)
out_files$mappoly_twopt_summary_tsv <- basename(summary_tsv)

write_artifacts(c(list(module="poly_mappoly_twopt"), out_files))
log("done")
