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
  write_error(msg); write_artifacts(list(module="poly_mappoly_export_map", error=msg)); quit(status=1)
}

# NOTE: do not hard-code any upstream temp paths here; export should depend only on
# the provided map_list_rds and parameters.
dat2 <- readRDS("/tmp/poly_mappoly_import_47ysk20g/out/mappoly_data.rds")
map_list_rds <- as.character(params$map_list_rds %||% "")
use_updated <- isTRUE(params$use_updated %||% TRUE)
export_csv <- isTRUE(params$export_csv %||% TRUE)
export_qtlpoly <- isTRUE(params$export_qtlpoly %||% FALSE)
step <- as.numeric(params$step %||% 1)
error <- as.numeric(params$error %||% 0.05)

if (!nzchar(map_list_rds) || !file.exists(map_list_rds)) {
  msg <- paste0("map_list_rds not found: ", map_list_rds)
  write_error(msg); write_artifacts(list(module="poly_mappoly_export_map", error=msg)); quit(status=1)
}

maps <- readRDS(map_list_rds)

out_files <- list()

# Export each LG map to CSV using export_map_list
if (isTRUE(export_csv)) {
  dir.create(file.path(out_dir, "maps_csv"), showWarnings=FALSE)
  for (nm in names(maps)) {
    obj <- maps[[nm]]
    mp <- if (isTRUE(use_updated) && !is.null(obj$map_updated)) obj$map_updated else obj$map
    csvf <- file.path(out_dir, "maps_csv", paste0("LG", nm, ".csv"))
    try(mappoly::export_map_list(mp, file = csvf), silent=TRUE)
  }
  out_files$maps_csv_dir <- "maps_csv"
}

# Optional: export to qtlpoly object (per LG)
if (isTRUE(export_qtlpoly)) {
  qtl_list <- list()
  for (nm in names(maps)) {
    obj <- maps[[nm]]
    mp <- if (isTRUE(use_updated) && !is.null(obj$map_updated)) obj$map_updated else obj$map
    g1 <- mappoly::calc_genoprob_error(mp, step = step, error = error)
    qtl_list[[nm]] <- mappoly::export_qtlpoly(g1)
  }
  qtl_rds <- file.path(out_dir, "qtlpoly_export.rds")
  saveRDS(qtl_list, qtl_rds)
  out_files$qtlpoly_export_rds <- basename(qtl_rds)
}

maps2 <- c()
for (cyc1 in 1:length(maps)) {
  maps2 <- c(maps2,
             list(maps[[cyc1]]$map))
}
g1_all <- lapply(maps2, mappoly::calc_genoprob)
qtlall_rds <- file.path(out_dir, "qtlpolyall_export.rds")
saveRDS(g1_all, qtlall_rds)

# Summary list of exported files
summary_tsv <- file.path(out_dir, "export_summary.tsv")
utils::write.table(data.frame(
  key=names(out_files), value=unlist(out_files), stringsAsFactors=FALSE
), summary_tsv, sep="\t", quote=FALSE, row.names=FALSE)
out_files$export_summary_tsv <- basename(summary_tsv)

write_artifacts(c(list(module="poly_mappoly_export_map"), out_files))
log("done")
