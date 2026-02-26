#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(jsonlite)
  library(data.table)
  library(qtl)
})

args <- commandArgs(trailingOnly=TRUE)
get_arg <- function(flag, default=NULL) {
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) return(args[i+1])
  default
}
params_path <- get_arg("--params")
out_dir <- get_arg("--out")
if (is.null(params_path) || is.null(out_dir)) stop("Usage: --params params.json --out out_dir")
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

p <- fromJSON(params_path)

cross_rds <- p$cross_rds
if (is.null(cross_rds) || !file.exists(cross_rds)) stop("cross_rds not found")

base_out <- out_dir
if (!is.null(p$out_dir_override) && nzchar(p$out_dir_override)) base_out <- as.character(p$out_dir_override)
base_out <- normalizePath(base_out, winslash="/", mustWork=FALSE)
dir.create(base_out, recursive=TRUE, showWarnings=FALSE)

cat("[prep_rqtl_to_tsv] cross_rds=", cross_rds, "\n")
cat("[prep_rqtl_to_tsv] out_dir=", base_out, "\n")

cross <- readRDS(cross_rds)

# genotype matrix (individuals x markers)
X <- qtl::pull.geno(cross)
# In R/qtl, 0 often indicates missing
X[X == 0] <- NA

# Convert genotype codes to dosage-like 0/1/2 when possible
u <- sort(unique(as.integer(na.omit(as.vector(X)))))
maxv <- if (length(u)==0) NA_integer_ else max(u)
minv <- if (length(u)==0) NA_integer_ else min(u)

X_dose <- X
if (!is.na(maxv)) {
  if (all(u %in% c(1L,2L,3L))) {
    X_dose <- X - 1L
  } else if (all(u %in% c(1L,2L))) {
    # BC-style (no heterozygotes): map 1->0, 2->2
    X_dose <- (X - 1L) * 2L
  } else if (all(u %in% c(1L,3L))) {
    X_dose <- X
    X_dose[X == 1L] <- 0L
    X_dose[X == 3L] <- 2L
  } else if (all(u %in% c(0L,1L,2L))) {
    X_dose <- X
  } else {
    # Fallback: keep numeric as-is (with missing already NA)
    X_dose <- X
  }
}

ids <- NULL
# Prefer phenotype rownames / id-like columns when available
if (!is.null(cross$pheno) && nrow(cross$pheno) == nrow(X_dose)) {
  rid <- rownames(cross$pheno)
  if (!is.null(rid) && length(rid) == nrow(X_dose) && all(nzchar(as.character(rid)))) {
    ids <- as.character(rid)
  } else {
    ph <- cross$pheno
    cand_cols <- c('id','ID','sample','Sample','ind','Ind','individual','Individual','genotype','Genotype','name','Name')
    hit <- cand_cols[cand_cols %in% colnames(ph)]
    if (length(hit) > 0) {
      v <- as.character(ph[[hit[1]]])
      if (length(v) == nrow(X_dose) && any(nzchar(v))) ids <- v
    }
  }
}
if (is.null(ids)) {
  ridX <- rownames(X_dose)
  if (!is.null(ridX) && length(ridX) == nrow(X_dose) && all(nzchar(as.character(ridX)))) ids <- as.character(ridX)
}
if (is.null(ids)) ids <- as.character(seq_len(nrow(X_dose)))
markers <- colnames(X_dose)

# Write genotype.tsv
out_geno <- file.path(base_out, "genotype.tsv")
DT <- as.data.table(X_dose)
DT[, id := ids]
setcolorder(DT, c("id", markers))
fwrite(DT, out_geno, sep="\t", quote=FALSE, na="NA")

# marker_map.tsv
out_map <- file.path(base_out, "marker_map.tsv")
mp <- qtl::pull.map(cross)
map_dt_list <- list()
for (chr in names(mp)) {
  v <- mp[[chr]]
  if (is.null(names(v))) {
    mnames <- paste0("m", seq_along(v))
  } else {
    mnames <- names(v)
  }
  map_dt_list[[chr]] <- data.table(marker=mnames, chr=as.character(chr), pos=as.numeric(v))
}
mm <- rbindlist(map_dt_list, use.names=TRUE, fill=TRUE)
# align to genotype marker order
if (!is.null(markers) && length(markers) > 0 && nrow(mm) > 0) {
  idx <- match(markers, mm$marker)
  ok <- !is.na(idx)
  out_mm <- data.table(marker=markers, chr="1", pos=seq_along(markers))
  out_mm$chr[ok] <- mm$chr[idx[ok]]
  out_mm$pos[ok] <- mm$pos[idx[ok]]
  out_mm$chr[is.na(out_mm$chr) | out_mm$chr==""] <- "1"
  out_mm$pos[is.na(out_mm$pos)] <- seq_along(out_mm$pos)[is.na(out_mm$pos)]
  mm <- out_mm
}
fwrite(mm, out_map, sep="\t", quote=FALSE, na="NA")

# phenotype.tsv
out_ph <- file.path(base_out, "phenotype.tsv")
ph <- cross$pheno
if (is.null(ph) || nrow(ph) == 0) {
  ph_dt <- data.table(id=ids, pheno=seq_along(ids))
} else {
  ph_dt <- as.data.table(ph)
  ph_dt[, id := as.character(ids)]
  setcolorder(ph_dt, c("id", setdiff(names(ph_dt), "id")))
  # order by ids if possible
  if ("id" %in% names(ph_dt)) {
    idx <- match(ids, as.character(ph_dt$id))
    if (any(!is.na(idx))) ph_dt <- ph_dt[idx, , drop=FALSE]
  }
}
fwrite(ph_dt, out_ph, sep="\t", quote=FALSE, na="NA")

meta <- list(
  genotype_tsv = out_geno,
  marker_map_tsv = out_map,
  phenotype_tsv = out_ph,
  source_cross_rds = cross_rds,
  out_dir = base_out
)
writeLines(toJSON(meta, auto_unbox=TRUE, pretty=TRUE), con=file.path(out_dir, "artifacts.json"))
