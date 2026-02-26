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
if (is.null(params_path) || is.null(out_dir)) stop("Usage: --params params.json --out out_dir")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_dir, "run.log")
sink(log_file, split = TRUE)

cat("[prep_vcf_to_cross] start\n")
p <- fromJSON(params_path)

vcf_path <- p$vcf_path
phe_tsv  <- if (!is.null(p$phenotype_tsv)) p$phenotype_tsv else NULL
cross_type <- if (!is.null(p$cross_type) && nchar(p$cross_type) > 0) p$cross_type else "f2"

missing_max <- 1.0
if (!is.null(p$missing_max)) suppressWarnings({ missing_max <- as.numeric(p$missing_max) })
min_maf <- 0.0
if (!is.null(p$min_maf)) suppressWarnings({ min_maf <- as.numeric(p$min_maf) })

if (is.null(vcf_path) || !file.exists(vcf_path)) stop("vcf_path is required and must exist")

cat("[prep_vcf_to_cross] vcf_path=", vcf_path, "\n")
cat("[prep_vcf_to_cross] cross_type=", cross_type, " missing_max=", missing_max, " min_maf=", min_maf, "\n")

# read VCF with data.table, skipping meta lines
dt <- fread(vcf_path, sep="\t", header=TRUE, skip="#CHROM", data.table=FALSE, showProgress=FALSE)

# Normalize column name "#CHROM"
if ("#CHROM" %in% colnames(dt)) colnames(dt)[colnames(dt)=="#CHROM"] <- "CHROM"

fixed <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT")
if (!all(fixed %in% colnames(dt))) stop("VCF is missing required columns; expected CHROM POS ID REF ALT ...")

sample_cols <- setdiff(colnames(dt), fixed)
if (length(sample_cols) == 0) stop("No samples found in VCF")

# marker names
marker <- dt$ID
marker[is.na(marker) | marker=="." | marker==""] <- paste0(dt$CHROM, ":", dt$POS)
marker <- make.unique(marker)

# biallelic check
is_biallelic <- !grepl(",", dt$ALT, fixed=TRUE)
dt <- dt[is_biallelic, , drop=FALSE]
marker <- marker[is_biallelic]

cat("[prep_vcf_to_cross] sites after biallelic filter=", nrow(dt), "\n")

# Extract GT field (first subfield before ':')
extract_gt <- function(x) sub(":.*$", "", x)
gt_mat <- sapply(sample_cols, function(cn) extract_gt(dt[[cn]]), simplify="matrix")
gt_mat <- t(gt_mat)  # samples x markers
rownames(gt_mat) <- sample_cols
colnames(gt_mat) <- marker

# Convert GT to dosage 0/1/2 (diploid biallelic); others -> NA
gt_to_dose <- function(z) {
  z <- gsub("\\|", "/", z)
  z[z == "." | z == "./." | z == ".|."] <- NA
  out <- rep(NA_real_, length(z))
  out[!is.na(z) & z %in% c("0/0")] <- 0
  out[!is.na(z) & z %in% c("0/1","1/0")] <- 1
  out[!is.na(z) & z %in% c("1/1")] <- 2
  out
}

dose_mat <- apply(gt_mat, 2, gt_to_dose)  # markers columns
dose_mat <- as.matrix(dose_mat)           # samples x markers

# marker filtering
miss_rate <- colMeans(is.na(dose_mat))
pA <- colMeans(dose_mat == 0, na.rm=TRUE) + 0.5 * colMeans(dose_mat == 1, na.rm=TRUE)
maf <- pmin(pA, 1 - pA)
keep <- (miss_rate <= missing_max) & (maf >= min_maf)
cat("[prep_vcf_to_cross] markers kept=", sum(keep), " / ", length(keep), "\n")

dose_mat <- dose_mat[, keep, drop=FALSE]
miss_rate <- miss_rate[keep]
maf <- maf[keep]
marker_kept <- colnames(dose_mat)

# phenotype
ph <- NULL
if (!is.null(phe_tsv) && file.exists(phe_tsv)) {
  ph <- fread(phe_tsv, sep="\t", header=TRUE, data.table=FALSE)
  if (!("id" %in% colnames(ph))) stop("phenotype_tsv must have 'id' column")
  ph$id <- as.character(ph$id)
  # match to sample order
  ph <- ph[match(rownames(dose_mat), ph$id), , drop=FALSE]
} else {
  # R/qtl requires at least one phenotype column with a non-missing value.
  ph <- data.frame(id=rownames(dose_mat), dummy=seq_len(nrow(dose_mat)))
}

# Convert to AA/AB/BB strings
to_gt_str <- function(v) {
  out <- rep(NA_character_, length(v))
  out[!is.na(v) & v==0] <- "AA"
  out[!is.na(v) & v==1] <- "AB"
  out[!is.na(v) & v==2] <- "BB"
  out
}
Xg <- apply(dose_mat, 2, to_gt_str)


# -------------------------------
# Fix3: Build cross object directly (no read.cross)
# -------------------------------

to_code <- function(v) {
  out <- rep(NA_integer_, length(v))
  out[v=="AA"] <- 1L
  out[v=="AB"] <- 2L
  out[v=="BB"] <- 3L
  out
}
Xcode <- apply(Xg, 2, to_code)
storage.mode(Xcode) <- "integer"
rownames(Xcode) <- rownames(dose_mat)
colnames(Xcode) <- marker_kept

ct <- tolower(cross_type)
if (ct %in% c("bc","bc1","bc2")) {
  n_bb <- sum(Xcode == 3L, na.rm=TRUE)
  if (n_bb > 0) cat("[prep_vcf_to_cross] warning: found BB in backcross; set to NA: ", n_bb, "\n")
  Xcode[Xcode == 3L] <- NA_integer_
}
if (ct %in% c("riself","risib","ril")) {
  n_ab <- sum(Xcode == 2L, na.rm=TRUE)
  if (n_ab > 0) cat("[prep_vcf_to_cross] warning: found AB in RIL; set to NA: ", n_ab, "\n")
  Xcode[Xcode == 2L] <- NA_integer_
}

# Build geno_list using VCF CHROM/POS (spec: keep chr labels as-is; pos numeric)
chr_kept <- as.character(dt$CHROM[keep])
pos_kept <- suppressWarnings(as.numeric(dt$POS[keep]))
pos_kept[is.na(pos_kept)] <- seq_len(sum(is.na(pos_kept)))

# Optional convenience scaling: if positions look like bp, convert to Mb (~cM)
mx <- suppressWarnings(max(pos_kept, na.rm=TRUE))
if (is.finite(mx) && mx > 1e4) {
  cat('[prep_vcf_to_cross] INFO: POS seems bp-scale (max=', mx, '); applying pos/1e6 (1Mb≈1cM)\n')
  pos_kept <- pos_kept / 1e6
}

geno_list <- list()
chr_order <- unique(chr_kept)
for (cc in chr_order) {
  idx <- which(chr_kept == cc)
  submat <- Xcode[, idx, drop=FALSE]
  mks <- colnames(submat)
  mp <- pos_kept[idx]
  o <- order(mp, na.last=TRUE)
  geno_list[[as.character(cc)]] <- list(
    data = submat[, o, drop=FALSE],
    map = setNames(as.numeric(mp[o]), mks[o]),
    alleles = c('A','B')
  )
}

ph_df <- as.data.frame(ph)
rownames(ph_df) <- ph_df$id
ph_df$id <- NULL

# ASMap::mstmap expects a unique genotype identifier column in cross$pheno.
if (!("Genotype" %in% colnames(ph_df))) {
  ph_df$Genotype <- rownames(ph_df)
}

cross <- list(pheno = ph_df, geno = geno_list)
# IMPORTANT: R/qtl/ASMap expect the cross-type class to be the first element
# (e.g. c("bc","cross")). ASMap::mstmap.cross branches on class(object)[1]
# and will fail if "cross" is first.
class(cross) <- c(cross_type, "cross")

# Summaries
miss_by_ind <- rowMeans(is.na(Xcode))
sample_dt <- data.table(id = rownames(Xcode), missing_rate = miss_by_ind)
fwrite(sample_dt, file.path(out_dir, "sample_summary.tsv"), sep="\t")

mar_dt <- data.table(marker = marker_kept, missing_rate = miss_rate, maf = maf,
                     n_nonmissing = colSums(!is.na(dose_mat)))
setorder(mar_dt, missing_rate, -maf)
fwrite(mar_dt, file.path(out_dir, "marker_summary.tsv"), sep="\t")

out_rds <- file.path(out_dir, "cross.rds")
saveRDS(cross, out_rds)
cat("[prep_vcf_to_cross] saved cross=", out_rds, "\n")

art <- list(
  cross_rds = out_rds,
  sample_summary_tsv = file.path(out_dir, "sample_summary.tsv"),
  marker_summary_tsv = file.path(out_dir, "marker_summary.tsv")
)
write(toJSON(art, auto_unbox=TRUE, pretty=TRUE), file=file.path(out_dir, "artifacts.json"))

cat("[prep_vcf_to_cross] done\n")
sink()

