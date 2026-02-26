#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(jsonlite)
  library(data.table)
  library(qtl)
  library(ASMap)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) return(args[i + 1])
  default
}

params_path <- get_arg("--params")
out_dir <- get_arg("--out")
if (is.null(params_path) || is.null(out_dir)) stop("Usage: --params params.json --out out_dir")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(out_dir, "run.log")
sink(log_file, split = TRUE)

cat("[map_asmap] start\n")
cat("[map_asmap] params_path=", params_path, "\n")
cat("[map_asmap] out_dir=", out_dir, "\n")

p <- fromJSON(params_path)

cross_rds <- p$cross_rds
if (is.null(cross_rds) || !file.exists(cross_rds)) {
  stop("cross_rds is required and must exist (an RDS of an R/qtl 'cross' object)")
}

# Options (best-effort; we will only pass args supported by the installed ASMap)
bychr <- as.logical(p$bychr)
anchor <- as.logical(p$anchor)

dist_fun <- if (!is.null(p$dist_fun) && nchar(p$dist_fun) > 0) p$dist_fun else "kosambi"
objective_fun <- if (!is.null(p$objective_fun) && nchar(p$objective_fun) > 0) p$objective_fun else "COUNT"

p_value <- 1e-6
if (!is.null(p$p_value)) suppressWarnings({ p_value <- as.numeric(p$p_value) })
missing_threshold <- NA_real_
if (!is.null(p$missing_threshold)) suppressWarnings({ missing_threshold <- as.numeric(p$missing_threshold) })
noMap_dist <- NA_real_
if (!is.null(p$noMap_dist)) suppressWarnings({ noMap_dist <- as.numeric(p$noMap_dist) })
noMap_size <- NA_integer_
if (!is.null(p$noMap_size)) suppressWarnings({ noMap_size <- as.integer(p$noMap_size) })

cat("[map_asmap] cross_rds=", cross_rds, "\n")
cat("[map_asmap] dist_fun=", dist_fun, " objective_fun=", objective_fun, "\n")
cat("[map_asmap] p_value=", p_value, " missing_threshold=", missing_threshold, " noMap_dist=", noMap_dist, " noMap_size=", noMap_size, "\n")

cross <- readRDS(cross_rds)
if (!inherits(cross, "cross")) {
  stop("The provided RDS is not an object of class 'cross' (R/qtl).")
}

# IMPORTANT: ASMap::mstmap.cross branches on class(object)[1]. Some upstream
# constructors may create objects with class c("cross", "bc") etc. Ensure the
# cross-type class is first, e.g. c("bc","cross").
cls <- class(cross)
if (length(cls) >= 2 && cls[1] == "cross") {
  # Move "cross" to the end while preserving order of other classes.
  cls2 <- c(cls[cls != "cross"], "cross")
  class(cross) <- cls2
  cat("[map_asmap] reordered class: ", paste(cls, collapse=","), " -> ", paste(cls2, collapse=","), "\n")
}
if (class(cross)[1] == "f2") {
  cross <- convert2bcsft(cross, F.gen = 2)
}

g <- pull.geno(cross)  # individuals x markers

mono_markers <- colnames(g)[apply(g, 2, function(x) {
  ux <- unique(x[!is.na(x)])
  length(ux) <= 1
})]

length(mono_markers)
head(mono_markers)

cross <- drop.markers(cross, mono_markers)


# Ensure ASMap can find the genotype identifier column. By default mstmap.cross
# looks for a column named "Genotype" in cross$pheno.
if (!is.null(cross$pheno)) {
  if (!("Genotype" %in% colnames(cross$pheno))) {
    cross$pheno$Genotype <- rownames(cross$pheno)
    cat("[map_asmap] added pheno$Genotype from rownames\n")
  }
}

# Run mapping
f <- ASMap::mstmap
#fn_formals <- names(formals(f))
#first_arg <- fn_formals[1]

call_args <- list()
#call_args[[first_arg]] <- cross
call_args[["object"]] <- cross
# common options (names vary by ASMap versions; we filter by formals)
call_args[["dist.fun"]] <- dist_fun
call_args[["objective.fun"]] <- objective_fun
call_args[["p.value"]] <- p_value

call_args[["bychr"]] <- bychr
call_args[["anchor"]] <- anchor

if (!is.na(missing_threshold)) {
  call_args[["miss.thresh"]] <- missing_threshold
}
if (!is.na(noMap_dist)) {
  call_args[["noMap.dist"]] <- noMap_dist
}
if (!is.na(noMap_size)) {
  call_args[["noMap.size"]] <- noMap_size
}

# Quiet mode if available
call_args[["trace"]] <- FALSE
call_args[["verbose"]] <- FALSE

# Keep only args supported by this installed ASMap::mstmap
#call_args <- call_args[names(call_args) %in% fn_formals]

cat("[map_asmap] calling mstmap with args:\n")
print(call_args)

out_cross <- do.call(f, call_args)

# Export map table
map_list <- qtl::pull.map(out_cross)
# map_list is a named list of numeric vectors (positions) with marker names
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

# Lengths summary
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
        main = "Linkage group lengths")
par(op)
dev.off()

# Save RDS
saveRDS(out_cross, file.path(out_dir, "linkage_map_cross.rds"))

# Artifacts
art <- list(
  plugin = "map_asmap",
  inputs = list(cross_rds = cross_rds),
  outputs = list(
    map_markers_tsv = "map_markers.tsv",
    map_lengths_tsv = "map_lengths.tsv",
    plot_png = "map_lengths.png",
    map_cross_rds = "linkage_map_cross.rds"
  )
)
writeLines(toJSON(art, auto_unbox = TRUE, pretty = TRUE), file.path(out_dir, "artifacts.json"))

cat("[map_asmap] done\n")
