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
#dat2 <- readRDS("/tmp/poly_mappoly_import_47ysk20g/out/mappoly_data.rds")
mappoly_data_rds <- as.character(params$mappoly_data_rds %||% "")
dat2 <- readRDS(mappoly_data_rds)

map_list_rds <- as.character(params$map_list_rds %||% "")
use_updated <- isTRUE(params$use_updated %||% TRUE)
export_csv <- isTRUE(params$export_csv %||% TRUE)
export_qtlpoly <- isTRUE(params$export_qtlpoly %||% FALSE)
step <- as.numeric(params$step %||% 0)
error <- as.numeric(params$error %||% 0.01)

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
  
  map_files <- list.files(path = file.path(out_dir, "maps_csv"),
                          pattern = "*.csv",
                          full.names = T)
  map_markers <- c()
  map_lengths <- c()
  for (cyc1 in 1:length(map_files)) {
    obj <- data.table::fread(map_files[[cyc1]],
                             header = T,
                             data.table = F)
    map_markers <- rbind(map_markers,
                         obj)
    
    map_lengths <- rbind(map_lengths,
                         data.frame(chr = obj[1, 3],
                                    n_markers = nrow(obj),
                                    length_cM = max(obj[, 7])
))
  }
  map_markers_out <- data.frame(chr = map_markers[, 3],
                                marker = map_markers[, 1],
                                cM = map_markers[, 7])
  data.table::fwrite(as.data.frame(map_markers_out),
                     file.path(out_dir, "map_markers.tsv"),
                     quote = F,
                     sep = "\t",
                     na = "NA",
                     row.names = F,
                     col.names = T)
  data.table::fwrite(as.data.frame(map_lengths),
                     file.path(out_dir, "map_lengths.tsv"),
                     quote = F,
                     sep = "\t",
                     na = "NA",
                     row.names = F,
                     col.names = T)
}

# Optional: export to qtlpoly object (per LG)
if (isTRUE(export_qtlpoly)) {
  qtl_list <- list()
  for (nm in names(maps)) {
    obj <- maps[[nm]]
    mp <- if (isTRUE(use_updated) && !is.null(obj$map_updated)) obj$map_updated else obj$map
    g1 <- mappoly::calc_genoprob_error(mp, 
                                       step = step,
                                       phase.config = "best",
                                       error = error,
                                       th.prob = 0.95,
                                       restricted = TRUE,
                                       verbose = TRUE)
    qtl_list[[nm]] <- mappoly::export_qtlpoly(g1)
  }
  qtl_rds <- file.path(out_dir, "qtlpoly_export.rds")
  saveRDS(qtl_list, qtl_rds)
  out_files$qtlpoly_export_rds <- basename(qtl_rds)
}

maps2 <- c()
for (cyc1 in 1:length(maps)) {
  obj <- maps[[cyc1]]
  mp <- if (isTRUE(use_updated) && !is.null(obj$map_updated)) obj$map_updated else obj$map
  maps2 <- c(maps2,
             list(mp))
}
#g1_all <- lapply(maps2, mappoly::calc_genoprob)
g1_all <- lapply(maps2, 
                 mappoly::calc_genoprob_error,
                 step = step,
                 phase.config = "best",
                 error = error,
                 th.prob = 0.95,
                 restricted = T,
                 verbose = F)
qtlall_rds <- file.path(out_dir, "qtlpoly_export_complete.rds")
saveRDS(g1_all, qtlall_rds)

# # Optional: export to qtlpoly object (per LG)
# if (isTRUE(export_qtlpoly)) {
#   qtl_list3 <- list()
#   for (nm in names(maps)) {
#     obj <- maps[[nm]]
#     mp <- if (isTRUE(use_updated) && !is.null(obj$map_updated)) obj$map_updated else obj$map
#     g3 <- mappoly::calc_genoprob_error(mp, 
#                                        step = step,
#                                        phase.config = "best",
#                                        error = error,
#                                        th.prob = 0.95,
#                                        restricted = TRUE,
#                                        verbose = TRUE)
#     qtl_list[[nm]] <- mappoly::export_qtlpoly(g3)
#     #qtl_list[[nm]] <- g3
#   }
#   qtl_rds3 <- file.path(out_dir, "qtlpoly_export3.rds")
#   saveRDS(qtl_list3, qtl_rds3)
#   out_files$qtlpoly_export_rds3 <- basename(qtl_rds3)
# }

# Summary list of exported files
summary_tsv <- file.path(out_dir, "export_summary.tsv")
utils::write.table(data.frame(
  key=names(out_files), value=unlist(out_files), stringsAsFactors=FALSE
), summary_tsv, sep="\t", quote=FALSE, row.names=FALSE)
out_files$export_summary_tsv <- basename(summary_tsv)

write_artifacts(c(list(module="poly_mappoly_export_map"), out_files))
log("done")
