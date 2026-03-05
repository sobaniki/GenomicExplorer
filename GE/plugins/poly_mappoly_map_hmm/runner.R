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
  write_error(msg); write_artifacts(list(module="poly_mappoly_map_hmm", error=msg)); quit(status=1)
}

#dat2 <- readRDS("/tmp/poly_mappoly_import_47ysk20g/out/mappoly_data.rds")
mappoly_data_rds <- as.character(params$mappoly_data_rds %||% "")
dat2 <- readRDS(mappoly_data_rds)

order_rds <- as.character(params$order_rds %||% "")
ncpus <- as.integer(params$ncpus %||% 1)
start_set <- as.integer(params$start_set %||% 4)
thres_twopt <- as.numeric(params$thres_twopt %||% 5)
thres_hmm <- as.numeric(params$thres_hmm %||% 50)
extend_tail <- as.integer(params$extend_tail %||% 50)
info_tail <- isTRUE(params$info_tail %||% TRUE)
sub_map_size_diff_limit <- as.integer(params$sub_map_size_diff_limit %||% 5)
phase_number_limit <- as.integer(params$phase_number_limit %||% 20)
tol <- as.numeric(params$tol %||% 1e-2)
tol_final <- as.numeric(params$tol_final %||% 1e-4)
update_global_error <- isTRUE(params$update_global_error %||% FALSE)
global_error <- as.numeric(params$global_error %||% 0.05)

if (!nzchar(order_rds) || !file.exists(order_rds)) {
  msg <- paste0("order_rds not found: ", order_rds)
  write_error(msg); write_artifacts(list(module="poly_mappoly_map_hmm", error=msg)); quit(status=1)
}

ord <- readRDS(order_rds)

maps <- list()
summary <- data.frame(lg=integer(), n_markers=integer(), map_length=numeric(), updated=numeric(), stringsAsFactors=FALSE)

for (nm in names(ord)) {
  o <- ord[[nm]]
  lg <- o$lg
  s <- o$seq_ordered
  log("HMM LG", lg, ": markers=", length(s$seq.num %||% integer()))

  # detailed pairwise for this LG
  tpt1 <- mappoly::est_pairwise_rf(input.seq = s,
                                   count.cache = NULL,
                                   count.matrix = NULL,
                                   ncpus = ncpus,
                                   mrk.pairs = NULL,
                                   n.batches = 1L,
                                   est.type = "disc",
                                   verbose = TRUE,
                                   memory.warning = TRUE,
                                   parallelization.type = "PSOCK",
                                   tol = .Machine$double.eps ^ 0.25,
                                   ll = FALSE)

  mp <- mappoly::est_rf_hmm_sequential(
    input.seq = s,
    start.set = start_set,
    thres.twopt = thres_twopt,
    thres.hmm = thres_hmm,
    extend.tail = extend_tail,
    info.tail = info_tail,
    twopt = tpt1,
    sub.map.size.diff.limit = sub_map_size_diff_limit,
    phase.number.limit = phase_number_limit,
    reestimate.single.ph.configuration = TRUE,
    tol = tol,
    tol.final = tol_final,
    verbose = TRUE,
    detailed.verbose = FALSE,
    high.prec = FALSE
  )

  mp_up <- NULL
  if (isTRUE(update_global_error)) {
    mp_up <- mappoly::est_full_hmm_with_global_error(input.map = mp, 
                                                     error = global_error,
                                                     tol = 0.001,
                                                     restricted = TRUE,
                                                     th.prob = 0.95,
                                                     verbose = TRUE)
  }

  maps[[nm]] <- list(map=mp, map_updated=mp_up, twopt=tpt1)

  # best-effort map length extraction (depends on MAPpoly version)
  infer_map_length <- function(map){
    if (is.null(map)) return(NA_real_)
    for (nm in c("pos", "position", "cM", "cm")) {
      if (!is.null(map[[nm]]) && is.numeric(map[[nm]])) {
        return(suppressWarnings(max(map[[nm]], na.rm=TRUE)))
      }
    }
    # sometimes stored in a data.frame returned by extract_map()
    ex <- tryCatch(mappoly::extract_map(map), error=function(e) NULL)
    if (!is.null(ex)) {
      for (nm in c("pos", "position", "cM", "cm")) if (nm %in% colnames(ex)) {
        return(suppressWarnings(max(as.numeric(ex[[nm]]), na.rm=TRUE)))
      }
    }
    NA_real_
  }
  len0 <- infer_map_length(mp)
  len1 <- if (isTRUE(update_global_error) && !is.null(mp_up)) infer_map_length(mp_up) else NA_real_
  summary <- rbind(summary, data.frame(lg=lg, n_markers=length(s$seq.num %||% integer()), map_length=len0, updated_length=len1, stringsAsFactors=FALSE))

  # plot
  png(file.path(out_dir, paste0("lg", lg, "_map.png")), width=1400, height=900)
  try(plot(mp, mrk.names=FALSE, cex=0.6), silent=TRUE)
  dev.off()
  if (!is.null(mp_up)) {
    png(file.path(out_dir, paste0("lg", lg, "_map_updated.png")), width=1400, height=900)
    try(plot(mp_up, mrk.names=FALSE, cex=0.6), silent=TRUE)
    dev.off()
  }
}

map_rds <- file.path(out_dir, "mappoly_map_list.rds")
saveRDS(maps, map_rds)
summary_tsv <- file.path(out_dir, "mappoly_map_summary.tsv")
utils::write.table(summary, summary_tsv, sep="\t", quote=FALSE, row.names=FALSE)

write_artifacts(list(
  module="poly_mappoly_map_hmm",
  mappoly_map_list_rds=basename(map_rds),
  mappoly_map_summary_tsv=basename(summary_tsv)
))
log("done")
