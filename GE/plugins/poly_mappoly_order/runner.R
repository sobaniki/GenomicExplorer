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
  write_error(msg); write_artifacts(list(module="poly_mappoly_order", error=msg)); quit(status=1)
}

rf_matrix_rds <- as.character(params$rf_matrix_rds %||% "")
group_rds <- as.character(params$group_rds %||% "")
lg_id <- params$lg_id %||% "1"
order_method <- as.character(params$order_method %||% "mds") # mds|genome
use_genomic_info <- isTRUE(params$use_genomic_info %||% TRUE)

if (!nzchar(rf_matrix_rds) || !file.exists(rf_matrix_rds)) {
  msg <- paste0("rf_matrix_rds not found: ", rf_matrix_rds)
  write_error(msg); write_artifacts(list(module="poly_mappoly_order", error=msg)); quit(status=1)
}
if (!nzchar(group_rds) || !file.exists(group_rds)) {
  msg <- paste0("group_rds not found: ", group_rds)
  write_error(msg); write_artifacts(list(module="poly_mappoly_order", error=msg)); quit(status=1)
}

m <- readRDS(rf_matrix_rds)
g <- readRDS(group_rds)

lg_ids <- NULL
if (is.character(lg_id) && lg_id == "all") {
  # best-effort infer number of groups
  ng <- g$n.groups %||% g$ngroups %||% NA_integer_
  if (is.na(ng)) {
    lg_ids <- 1:12
    log("[WARN] cannot infer number of groups; defaulting lg_ids=1:12")
  } else {
    lg_ids <- seq_len(as.integer(ng))
  }
} else {
  lg_ids <- as.integer(lg_id)
}

#dat2 <- readRDS("/tmp/poly_mappoly_import_47ysk20g/out/mappoly_data.rds")
mappoly_data_rds <- as.character(params$mappoly_data_rds %||% "")
dat2 <- readRDS(mappoly_data_rds)

results <- list()
for (i in lg_ids) {
  log("ordering LG", i)

  s <- mappoly::make_seq_mappoly(g, i, genomic.info = if (isTRUE(use_genomic_info)) i else NULL)
 
  m1 <- mappoly::make_mat_mappoly(m, s)

  if (order_method == "genome") {
    go <- mappoly::get_genomic_order(s)
    s_ord <- mappoly::make_seq_mappoly(go)
    mds_o <- NULL
  } else {
    mds_o <- mappoly::mds_mappoly(input.mat = m1)
    s_ord <- mappoly::make_seq_mappoly(mds_o)
    go <- tryCatch(mappoly::get_genomic_order(s), error=function(e) NULL)
  }

  results[[as.character(i)]] <- list(
    lg=i,
    seq_base=s,
    mat=m1,
    order_method=order_method,
    mds=mds_o,
    seq_ordered=s_ord,
    genomic_order=go
  )

  # plot diagnostic
  png(file.path(out_dir, paste0("lg", i, "_order.png")), width=1200, height=1000)
  try(plot(m1, ord=s_ord, fact=5), silent=TRUE)
  dev.off()
}

order_rds <- file.path(out_dir, "mappoly_order.rds")
saveRDS(results, order_rds)

# summary
sum_df <- data.frame(
  lg=lg_ids,
  n_markers=sapply(results, function(x) length(x$seq_ordered$seq.num %||% integer())),
  order_method=order_method,
  stringsAsFactors=FALSE
)
summary_tsv <- file.path(out_dir, "mappoly_order_summary.tsv")
utils::write.table(sum_df, summary_tsv, sep="\t", quote=FALSE, row.names=FALSE)

write_artifacts(list(
  module="poly_mappoly_order",
  mappoly_order_rds=basename(order_rds),
  mappoly_order_summary_tsv=basename(summary_tsv)
))
log("done")
