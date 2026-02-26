#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(jsonlite)
  library(data.table)
  library(qtl)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) return(args[i + 1])
  default
}

params_path <- get_arg("--params")
out_dir <- get_arg("--out")
if (is.null(params_path) || is.null(out_dir)) {
  stop("Usage: --params params.json --out out_dir")
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(out_dir, "run.log")
sink(log_file, split = TRUE)

cat("[map_rqtl] start\n")
cat("[map_rqtl] params_path=", params_path, "\n")
cat("[map_rqtl] out_dir=", out_dir, "\n")

p <- fromJSON(params_path)

cross_rds <- p$cross_rds
if (is.null(cross_rds) || !file.exists(cross_rds)) {
  stop("cross_rds is required and must exist (an RDS of an R/qtl 'cross' object)")
}

map_function <- "kosambi"
if (!is.null(p$map_function) && nchar(p$map_function) > 0) map_function <- tolower(p$map_function)

error_prob <- 0.001
if (!is.null(p$error_prob)) suppressWarnings({ error_prob <- as.numeric(p$error_prob) })

jittermap_flag <- "TRUE"
if (!is.null(p$jittermap) && nchar(p$jittermap) > 0) jittermap_flag <- toupper(p$jittermap)

do_jittermap <- identical(jittermap_flag, "TRUE")

# Advanced (best-effort)
form_lg_flag <- "FALSE"
if (!is.null(p$form_linkage_groups) && nchar(p$form_linkage_groups) > 0) form_lg_flag <- toupper(p$form_linkage_groups)
do_form_lg <- identical(form_lg_flag, "TRUE")

max_rf <- 0.25
if (!is.null(p$max_rf)) suppressWarnings({ max_rf <- as.numeric(p$max_rf) })
min_lod <- 3
if (!is.null(p$min_lod)) suppressWarnings({ min_lod <- as.numeric(p$min_lod) })

order_flag <- "FALSE"
if (!is.null(p$order_markers) && nchar(p$order_markers) > 0) order_flag <- toupper(p$order_markers)
do_order <- identical(order_flag, "TRUE")
order_window <- 7
if (!is.null(p$order_window)) suppressWarnings({ order_window <- as.integer(p$order_window) })
#order_method <- "likelihood"
#if (!is.null(p$order_method) && nchar(p$order_method) > 0) order_method <- tolower(p$order_method)

ripple_flag <- "FALSE"
if (!is.null(p$ripple) && nchar(p$ripple) > 0) ripple_flag <- toupper(p$ripple)
do_ripple <- identical(ripple_flag, "TRUE")
ripple_window <- 7
if (!is.null(p$ripple_window)) suppressWarnings({ ripple_window <- as.integer(p$ripple_window) })
ripple_method <- "likelihood"
if (!is.null(p$ripple_method) && nchar(p$ripple_method) > 0) ripple_method <- tolower(p$ripple_method)

chr_sel <- NULL
# if (!is.null(p$chromosomes) && nchar(trimws(p$chromosomes)) > 0) {
#   chr_sel <- trimws(unlist(strsplit(p$chromosomes, ",")))
#   if (length(chr_sel) == 0) chr_sel <- NULL
# }

cat("[map_rqtl] cross_rds=", cross_rds, "\n")
cat("[map_rqtl] map_function=", map_function, " error_prob=", error_prob, " jittermap=", do_jittermap, "\n")
cat("[map_rqtl] formLinkageGroups=", do_form_lg, " max.rf=", max_rf, " min.lod=", min_lod, "\n")
cat("[map_rqtl] orderMarkers=", do_order, " window=", order_window, "\n")
cat("[map_rqtl] ripple=", do_ripple, " window=", ripple_window, " method=", ripple_method, "\n")
#if (!is.null(chr_sel)) cat("[map_rqtl] chromosomes selected:", paste(chr_sel, collapse = ","), "\n")

cross <- readRDS(cross_rds)
if (!inherits(cross, "cross")) {
  stop("The provided RDS is not an object of class 'cross' (R/qtl).")
}

g <- pull.geno(cross)  # individuals x markers

mono_markers <- colnames(g)[apply(g, 2, function(x) {
  ux <- unique(x[!is.na(x)])
  length(ux) <= 1
})]

length(mono_markers)
head(mono_markers)

cross <- drop.markers(cross, mono_markers)

## Best-effort cleanup: jitter identical marker positions to avoid singularities
if (do_jittermap) {
  cross <- jittermap(cross,
                     amount = 1e-6)
}

## (Optional) Re-estimate linkage groups
if (do_form_lg) {
  cross <- est.rf(cross = cross,
                  maxit = 10000,
                  tol = 1e-6)
  cross <- formLinkageGroups(cross = cross, 
                             max.rf = max_rf, 
                             min.lod = min_lod,
                             reorgMarkers = T,
                             verbose = T)
}

## Determine chromosomes to process
chrs <- qtl::chrnames(cross)
if (is.null(chrs) || length(chrs) == 0) {
  chrs <- names(qtl::pull.map(cross))
}
if (!is.null(chr_sel)) {
  chrs <- intersect(chrs, chr_sel)
}

## (Optional) Marker ordering and/or ripple evaluation (best-effort)
marker_order_dt <- NULL
if (do_order) {
  for (chr in chrs) {
    cross <- orderMarkers(cross = cross, 
                          chr = chr, 
                          window = order_window,
                          use.ripple = T,
                          error.prob = error_prob,
                          map.function = map_function,
                          maxit = 4000,
                          tol = 1e-4,
                          sex.sp = T,
                          verbose = F)
  }
}

## Save current marker order (after optional ordering/linkage grouping)
m <- qtl::pull.map(cross)
ord_rows <- list()
for (chr in names(m)) {
  ord_rows[[length(ord_rows) + 1]] <- data.table(chr = chr, 
                                                 order_index = seq_along(m[[chr]]), 
                                                 marker = names(m[[chr]]))
}
marker_order_dt <- rbindlist(ord_rows, use.names = TRUE, fill = TRUE)
fwrite(marker_order_dt, file.path(out_dir, "marker_order.tsv"), sep = "\t")

ripple_rows <- list()
if (do_ripple) {
  for (chr in chrs) {
    cat("[map_rqtl] ripple chr=", chr, "\n")
    rr <- ripple(cross = cross, 
                 chr = chr, 
                 window = ripple_window, 
                 method = ripple_method,
                 error.prob = error_prob, 
                 map.function = map_function,
                 maxit = 4000,
                 tol = 1e-6,
                 sex.sp = T,
                 verbose = T,
                 n.cluster = 1)
    if (!is.null(rr)) {
      dt <- as.data.table(rr)
      dt[, chr := chr]
      ripple_rows[[length(ripple_rows) + 1]] <- dt
    }
  }
  if (length(ripple_rows) > 0) {
    ripple_dt <- rbindlist(ripple_rows, use.names = TRUE, fill = TRUE)
    fwrite(ripple_dt, file.path(out_dir, "ripple_results.tsv"), sep = "\t")
  }
}

# Estimate map given current marker order
cat("[map_rqtl] running qtl::est.map\n")
map <- est.map(cross = cross,
               #chr,
               error.prob = error_prob, 
               map.function = map_function,
               m = 0,
               p = 0,
               maxit = 10000,
               tol = 1e-6,
               sex.sp = T,
               verbose = T,
               omit.noninformative = T,
               #offset,
               n.cluster = 1)

# Replace map and export
out_cross <- replace.map(cross, map)
saveRDS(out_cross, "/media/soba/Noc4/GenomicExplorer/P0213/ex1.rds")
# Export marker map table
map_list <- qtl::pull.map(out_cross)
rows <- list()
for (chr in names(map_list)) {
  v <- map_list[[chr]]
  if (length(v) == 0) next
  dt <- data.table(chr = chr, marker = names(v), cM = as.numeric(v))
  rows[[length(rows) + 1]] <- dt
}
map_dt <- rbindlist(rows, use.names = TRUE, fill = TRUE)
setorder(map_dt, chr, cM)

map_tsv <- file.path(out_dir, "map_markers.tsv")
fwrite(map_dt, map_tsv, sep = "\t")

# Length summary
len_dt <- map_dt[, .(n_markers = .N, length_cM = max(cM, na.rm = TRUE)), by = chr]
setorder(len_dt, -length_cM)
len_tsv <- file.path(out_dir, "map_lengths.tsv")
fwrite(len_dt, len_tsv, sep = "\t")

# Plot lengths
png_path <- file.path(out_dir, "map_lengths.png")
png(png_path, width = 1100, height = 600)
op <- par(mar = c(8, 5, 3, 2))
barplot(height = len_dt$length_cM,
        names.arg = len_dt$chr,
        las = 2,
        ylab = "Length (cM)",
        main = "Linkage group lengths (r/qtl est.map)")
par(op)
dev.off()

# Save mapped cross
saveRDS(out_cross, file.path(out_dir, "linkage_map_cross.rds"))

# Artifacts (match map_asmap output keys where possible)
art <- list(
  plugin = "map_rqtl",
  algorithm = "est.map",
  inputs = list(cross_rds = cross_rds),
  params = list(
    map_function = map_function,
    error_prob = error_prob,
    jittermap = do_jittermap,
    form_linkage_groups = do_form_lg,
    max_rf = max_rf,
    min_lod = min_lod,
    order_markers = do_order,
    order_window = order_window,
    #order_method = order_method,
    ripple = do_ripple,
    ripple_window = ripple_window,
    ripple_method = ripple_method
    #chromosomes = if (is.null(chr_sel)) "" else paste(chr_sel, collapse = ",")
  ),
  outputs = list(
    map_markers_tsv = "map_markers.tsv",
    map_lengths_tsv = "map_lengths.tsv",
    plot_png = "map_lengths.png",
    map_cross_rds = "linkage_map_cross.rds",
    marker_order_tsv = if (file.exists(file.path(out_dir, "marker_order.tsv"))) "marker_order.tsv" else "",
    ripple_results_tsv = if (file.exists(file.path(out_dir, "ripple_results.tsv"))) "ripple_results.tsv" else ""
  )
)
writeLines(toJSON(art, auto_unbox = TRUE, pretty = TRUE), file.path(out_dir, "artifacts.json"))

cat("[map_rqtl] done\n")
