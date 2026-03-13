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
  write_error(msg); write_artifacts(list(module="poly_mappoly_group", error=msg)); quit(status=1)
}

rf_matrix_rds <- as.character(params$rf_matrix_rds %||% "")
expected_groups <- params$expected_groups
comp_mat <- isTRUE(params$comp_mat %||% TRUE)
inter <- isTRUE(params$inter %||% TRUE)

if (!nzchar(rf_matrix_rds) || !file.exists(rf_matrix_rds)) {
  msg <- paste0("rf_matrix_rds not found: ", rf_matrix_rds)
  write_error(msg); write_artifacts(list(module="poly_mappoly_group", error=msg)); quit(status=1)
}

m <- readRDS(rf_matrix_rds)

if (!is.null(expected_groups)) expected_groups <- as.integer(expected_groups)
log("group_mappoly expected_groups=", ifelse(is.null(expected_groups), "NULL", expected_groups), ", comp.mat=", comp_mat)

mappoly_data_rds <- as.character(params$mappoly_data_rds %||% "")
dat2 <- readRDS(mappoly_data_rds)

g <- mappoly::group_mappoly(input.mat = m, 
                            expected.groups = expected_groups, 
                            inter = inter, 
                            comp.mat = comp_mat,
                            LODweight = F,
                            verbose = T)

out_files <- list()

g_rds <- file.path(out_dir, "mappoly_group.rds")
saveRDS(g, g_rds)
out_files$mappoly_group_rds <- basename(g_rds)

# summary table (best-effort)
summary_tsv <- file.path(out_dir, "mappoly_group_summary.tsv")
try({
  tab <- as.data.frame(g$groups)
  utils::write.table(tab, summary_tsv, sep="\t", quote=FALSE, row.names=FALSE)
}, silent=TRUE)
out_files$mappoly_group_summary_tsv <- basename(summary_tsv)

png(file.path(out_dir, "group_plot.png"), width=1200, height=900)
try(plot(g), silent=TRUE)
dev.off()
out_files$group_plot_png <- "group_plot.png"

write_artifacts(c(list(module="poly_mappoly_group"), out_files))
log("done")
