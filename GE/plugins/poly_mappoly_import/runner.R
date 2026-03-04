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
  msg <- "R package 'mappoly' is required but not installed in this environment. Please install it in the plugin R environment."
  write_error(msg); write_artifacts(list(module="poly_mappoly_import", error=msg)); quit(status=1)
}

geno_rds <- as.character(params$geno_rds %||% "")
ploidy <- as.numeric(params$ploidy %||% 4)
parent1_id <- as.character(params$parent1_id %||% "")
parent2_id <- as.character(params$parent2_id %||% "")
use_prob <- isTRUE(params$use_prob %||% TRUE) # default TRUE (we guarantee prob)

if (!nzchar(geno_rds) || !file.exists(geno_rds)) {
  msg <- paste0("geno_rds not found: ", geno_rds)
  write_error(msg); write_artifacts(list(module="poly_mappoly_import", error=msg)); quit(status=1)
}

out_files <- list()
if (grepl("*.csv$", geno_rds)) {
  if (isTRUE(use_prob)) {
    dat <- mappoly::read_geno_prob(file.in = geno_rds,
                                   prob.thres = 0.95,
                                   filter.non.conforming = F,
                                   elim.redundant = F,
                                   verbose = T)
  } else {
    dat <- mappoly::read_geno_csv(file.in = geno_rds,
                                  ploidy = ploidy,
                                  filter.non.conforming = F,
                                  elim.redundant = F,
                                  verbose = T)
  }
  dat_rds <- file.path(out_dir, "mappoly_data.rds")
  saveRDS(dat, dat_rds)
  out_files$mappoly_data_rds <- basename(dat_rds)
  
  markers <- dat$n.mrk
  samples <- dat$n.ind
  offspring <- dat$n.ind - 2
  
  mappoly_csv <- file.path(out_dir, "mappoly_dosage.csv")
  out_files$mappoly_dosage_csv <- basename(mappoly_csv)
} else if (grepl("*.vcf\\(\\|.tar.gz\\)$", geno_rds, perl = F)) {
  dat <- mappoly::read_vcf(file.in = geno_rds,
                           parent.1,
                           parent.2,
                           ploidy = ploidy,
                           filter.non.conforming = T,
                           thresh.line = 0.05,
                           min.gt.depth = 0,
                           min.av.depth = 0,
                           max.missing = 1,
                           elim.redundant = T,
                           verbose = T,
                           read.geno.prob = F,
                           prob.thres = 0.95)
  
  dat_rds <- file.path(out_dir, "mappoly_data.rds")
  saveRDS(dat, dat_rds)
  out_files$mappoly_data_rds <- basename(dat_rds)
  
  markers <- dat$n.mrk
  samples <- dat$n.ind
  offspring <- dat$n.ind - 2
  
  mappoly_csv <- file.path(out_dir, "mappoly_dosage.csv")
  out_files$mappoly_dosage_csv <- basename(mappoly_csv)
}

# summary
summary_tsv <- file.path(out_dir, "mappoly_import_summary.tsv")
utils::write.table(data.frame(
  key=c("ploidy","n_markers","n_individuals","parent1_id","parent2_id","n_offspring"),
  value=c(ploidy, length(markers), length(samples), parent1_id, parent2_id, length(offspring)),
  stringsAsFactors=FALSE
), summary_tsv, sep="\t", quote=FALSE, row.names=FALSE)
out_files$mappoly_import_summary_tsv <- basename(summary_tsv)

write_artifacts(c(list(module="poly_mappoly_import"), out_files))
log("done")
