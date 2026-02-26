#!/usr/bin/env Rscript
suppressWarnings(suppressMessages({
  library(jsonlite)
  library(onemap)
}))

args <- commandArgs(trailingOnly=TRUE)
get_arg <- function(flag, default=NULL) {
  if (!(flag %in% args)) return(default)
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) return(default)
  args[[i+1]]
}

params_path <- get_arg("--params")
out_dir <- get_arg("--out")
if (is.null(params_path) || is.null(out_dir)) {
  stop("Usage: runner.R --params params.json --out out_dir")
}
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

log_file <- file.path(out_dir, "run.log")
log <- function(...) {
  msg <- paste0("[map_onemap] ", paste0(..., collapse=""))
  cat(msg, "\n")
  cat(msg, "\n", file=log_file, append=TRUE)
}

write_artifacts <- function(primary_table="map_markers.tsv", primary_plot="map_lengths.png", extra=list()) {
  art <- list(
    primary_table=primary_table,
    primary_plot=primary_plot,
    tables=unique(c(primary_table, "map_lengths.tsv", "groups.tsv", "marker_order.tsv", "onemap_summary.tsv")),
    plots=unique(c(primary_plot)),
    extra=extra
  )
  writeLines(toJSON(art, auto_unbox=TRUE, pretty=TRUE), con=file.path(out_dir, "artifacts.json"))
}

write_placeholder <- function(errmsg) {
  log("WARN: writing placeholder outputs. reason=", errmsg)
  writeLines(errmsg, con=file.path(out_dir, "error_message.txt"))
  # minimal table so GUI doesn't break
  df <- data.frame(group=character(), marker=character(), pos_cM=numeric(), order=integer(), stringsAsFactors=FALSE)
  fwrite_ok <- FALSE
  try({
    data.table::fwrite(df, file.path(out_dir, "map_markers.tsv"), sep="\t")
    fwrite_ok <- TRUE
  }, silent=TRUE)
  if (!fwrite_ok) {
    write.table(df, file.path(out_dir, "map_markers.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
  }
  df2 <- data.frame(group=character(), length_cM=numeric(), n_markers=integer(), stringsAsFactors=FALSE)
  try({
    if (requireNamespace("data.table", quietly=TRUE)) data.table::fwrite(df2, file.path(out_dir, "map_lengths.tsv"), sep="\t")
    else write.table(df2, file.path(out_dir, "map_lengths.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
  }, silent=TRUE)

  # create an empty placeholder plot
  png(file.path(out_dir, "map_lengths.png"), width=800, height=400)
  plot.new()
  title("OneMap: no plot (runner placeholder)")
  text(0.5, 0.5, errmsg, cex=0.9)
  dev.off()

  # other optional tables
  for (nm in c("groups.tsv","marker_order.tsv","onemap_summary.tsv")) {
    if (!file.exists(file.path(out_dir, nm))) {
      writeLines("", con=file.path(out_dir, nm))
    }
  }
  write_artifacts()
}

log("start")
log("params_path=", params_path)
log("out_dir=", out_dir)

params <- fromJSON(params_path)

# Inputs
input_mode <- if (!is.null(params$input_mode) && nzchar(params$input_mode)) params$input_mode else "rds"
onemap_rds <- if (!is.null(params$onemap_rds)) params$onemap_rds else ""
input_file <- if (!is.null(params$input_file)) params$input_file else ""
cross_type <- if (!is.null(params$cross_type) && nzchar(params$cross_type)) params$cross_type else "outcross"
parent1 <- params$parent1
parent2 <- params$parent2
#sep <- if (!is.null(params$sep) && nzchar(params$sep)) params$sep else ","
#na_str <- if (!is.null(params$na_str) && nzchar(params$na_str)) params$na_str else "NA"

# Options
do_group <- if (!is.null(params$do_group)) as.logical(toupper(params$do_group)) else TRUE
lod <- if (!is.null(params$lod) && nzchar(params$lod)) as.numeric(params$lod) else 3
max_rf <- if (!is.null(params$max_rf) && nzchar(params$max_rf)) as.numeric(params$max_rf) else 0.25

do_order <- if (!is.null(params$do_order)) as.logical(toupper(params$do_order)) else TRUE
order_search <- if (!is.null(params$order_search)) as.character(params$order_search) else "twopt"
order_twopt_alg <- if (!is.null(params$order_twopt_alg)) as.character(params$order_twopt_alg) else "rec"
do_mds <- if (!is.null(params$do_mds)) as.logical(toupper(params$do_mds)) else FALSE
#order_method <- if (!is.null(params$order_method) && nzchar(params$order_method)) params$order_method else "twopt"
subset_group <- if (!is.null(params$subset_group) && nzchar(params$subset_group)) params$subset_group else ""

#do_ripple <- if (!is.null(params$do_ripple)) as.logical(toupper(params$do_ripple)=="TRUE") else FALSE
#ripple_ws <- if (!is.null(params$ripple_ws) && nzchar(params$ripple_ws)) as.integer(params$ripple_ws) else 5
#ripple_method <- if (!is.null(params$ripple_method) && nzchar(params$ripple_method)) params$ripple_method else "twopt"

map_function <- if (!is.null(params$map_function) && nzchar(params$map_function)) params$map_function else "kosambi"
#do_save_rds <- if (!is.null(params$save_rds)) as.logical(toupper(params$save_rds)=="TRUE") else TRUE

#extra <- params$extra_options
#if (is.null(extra)) extra <- list()

# Ensure data.table for fast write if available
if (requireNamespace("data.table", quietly=TRUE)) {
  `%fwrite%` <- function(x, f) data.table::fwrite(x, f, sep="\t", quote = F)
} else {
  `%fwrite%` <- function(x, f) write.table(x, f, sep="\t", row.names=FALSE, quote=FALSE)
}

log("loading input: mode=", input_mode)

om <- NULL
if (tolower(input_mode) == "rds") {
  if (!nzchar(onemap_rds) || !file.exists(onemap_rds)) {
    write_placeholder("onemap_rds not set / not found")
    quit(status=0)
  }
  om <- tryCatch(readRDS(onemap_rds), error=function(e) e)
  if (inherits(om, "error")) {
    write_placeholder(paste0("Failed to read onemap_rds: ", om$message))
    quit(status=0)
  }
} else if (grepl(".vcf(|.tar.gz)$", input_file, perl = T)) {
  om <- onemap::onemap_read_vcfR(vcf = input_file,
                                 cross = cross_type,
                                 parent1 = parent1,
                                 parent2 = parent2)
} else {
  om <- onemap::read_onemap(inputfile = input_file,
                            dir = NULL,
                            verbose = T)
}

# Save original object for debugging
# if (do_save_rds) {
#   try(saveRDS(om, file=file.path(out_dir, "onemap_input.rds")), silent=TRUE)
# }

summary_lines <- c(
  paste0("input_mode\t", input_mode),
  paste0("cross_type\t", cross_type),
  paste0("do_group\t", do_group),
  paste0("lod\t", lod),
  paste0("max_rf\t", max_rf),
  paste0("do_order\t", do_order),
  #paste0("order_method\t", order_method),
  #paste0("do_ripple\t", do_ripple),
  #paste0("ripple_ws\t", ripple_ws),
  paste0("map_function\t", map_function)
)
writeLines(summary_lines, con=file.path(out_dir, "onemap_summary.tsv"))

bins <- find_bins(input.obj = om, 
                  exact = F)
bins_out <- create_data_bins(input.obj = om, 
                             bins = bins)
LOD_sug <- suggest_lod(bins_out)
twopts_f2 <- rf_2pts(input.obj = bins_out,
                     LOD = lod,
                     max.rf = max_rf,
                     verbose = T,
                     rm_mks = T)
mark_all_f2 <- make_seq(input.obj = twopts_f2, 
                        arg = "all",
                        phase = NULL,
                        data.name = NULL,
                        twopt = NULL)

# Grouping
if (do_group) {
  log("grouping: LOD=", lod, " max_rf=", max_rf)
  
  grp <- onemap::group(input.seq = mark_all_f2, 
                       LOD = LOD_sug, 
                       max.rf = max_rf,
                       verbose = T)
  
  gi <- grp$groups
  groups_tbl <- as.data.frame(table(gi))
  colnames(groups_tbl) <- c("group", "n_markers")
  groups_tbl %fwrite% file.path(out_dir, "groups.tsv")
} else {
  seq_obj <- mark_all_f2
}

# Choose group(s)
group_ids <- NULL
if (nzchar(subset_group)) {
  # allow comma list
  group_ids <- strsplit(subset_group, ",")[[1]]
  group_ids <- trimws(group_ids)
} else if (do_group) {
  group_ids <- 1:grp$n.groups
} else {
  group_ids <- unique(om$CHROM)
}

set_map_fun(type = map_function)

all_markers_out <- list()
map_lengths <- list()
marker_order_out <- list()
map_list <- c()
for (gid in group_ids) {
  if (do_order) {
    seq_obj <- onemap::make_seq(grp, gid)
    if (do_mds) {
      seq_obj <- mds_onemap(
        input.seq = seq_obj,
        out.file = NULL,
        p = NULL,
        ispc = T,
        displaytext = F,
        weightfn = "lod2",
        mapfn = map_function,
        ndim = 2,
        rm_unlinked = T,
        size = NULL,
        overlap = NULL,
        phase_cores = 1,
        tol = 1e-05,
        hmm = T,
        parallelization.type = "PSOCK")
    }
    ordered <- onemap::order_seq(input.seq = seq_obj,
                                 n.init = 5,
                                 subset.search = order_search,
                                 subset.n.try = 30,
                                 subset.THRES = 3,
                                 twopt.alg = order_twopt_alg,
                                 THRES = 3,
                                 touchdown = F,
                                 tol = 0.1,
                                 rm_unlinked = F,
                                 verbose = T)
    seq_ord <- onemap::make_seq(ordered, "force")
    
    mp <- onemap::map(input.seq = seq_ord,
                      tol = 1e-04,
                      verbose = FALSE,
                      rm_unlinked = FALSE,
                      phase_cores = 1,
                      parallelization.type = "PSOCK",
                      global_error = NULL,
                      genotypes_errors = NULL,
                      genotypes_probs = NULL)
  } else {
    mp <- onemap::map(input.seq = seq_obj,
                      tol = 1e-04,
                      verbose = FALSE,
                      rm_unlinked = FALSE,
                      phase_cores = 1,
                      parallelization.type = "PSOCK",
                      global_error = NULL,
                      genotypes_errors = NULL,
                      genotypes_probs = NULL)
  }
  
  map_list <- c(map_list,
                list(mp))
  
  pos <- c(0, cumsum(mp$seq.rf * 100))
  mk <- colnames(mp$data.name$geno)[mp$seq.num]
  
  tbl <- data.frame(group = gid, 
                    marker = mk, 
                    pos_cM = as.numeric(pos), 
                    order = mp$seq.num, 
                    stringsAsFactors = F)
  len_cM <- max(pos)
  nmk <- length(mp$seq.num)
        
  all_markers_out[[gid]] <- tbl
  marker_order_out[[gid]] <- tbl[, c("group", "marker", "order")]
  map_lengths[[gid]] <- data.frame(group = gid, 
                                   length_cM = len_cM, 
                                   n_markers = nmk, 
                                   stringsAsFactors = F)
}
saveRDS(map_list, file=file.path(out_dir, "map_group.rds"))
onemap::write_map(map.list = map_list,
                  file.out = file.path(out_dir, "map_group.onemap.map"))

map_markers <- do.call(rbind, all_markers_out)
map_lengths_df <- do.call(rbind, map_lengths)
marker_order_df <- do.call(rbind, marker_order_out)

map_markers %fwrite% file.path(out_dir, "map_markers.tsv")
map_lengths_df %fwrite% file.path(out_dir, "map_lengths.tsv")
marker_order_df %fwrite% file.path(out_dir, "marker_order.tsv")

# Plot map lengths
png(file.path(out_dir, "map_lengths.png"), width=800, height=500)
par(mar=c(10,5,2,1))
labs <- map_lengths_df$group
vals <- map_lengths_df$length_cM
if (all(is.na(vals))) {
  plot.new(); title("Map lengths (cM)"); text(0.5,0.5,"lengths not available", cex=1.1)
} else {
  barplot(vals, names.arg=labs, las=2, ylab="cM", main="OneMap: linkage group lengths")
}
dev.off()

# Save final map objects if any
try(saveRDS(list(input=om, grouped=grp), file=file.path(out_dir, "onemap_objects.rds")), silent=TRUE)

write_artifacts(primary_table="map_markers.tsv", primary_plot="map_lengths.png",
                extra=list(groups=group_ids, do_group=do_group, do_order=do_order))

log("done")
