#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(jsonlite)
  library(data.table)
  library(qtl2)
})

# -------- utils --------
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) return(args[i + 1])
  default
}

call_supported <- function(fun, arglist) {
  # call function with only supported named args (helps qtl2 version differences)
  fml <- tryCatch(names(formals(fun)), error=function(e) NULL)
  if (is.null(fml)) return(do.call(fun, arglist))
  if ("..." %in% fml) return(do.call(fun, arglist))
  keep <- names(arglist) %in% fml
  do.call(fun, arglist[keep])
}

calc_kinship_compat <- function(genoprobs, kind=c("loco", "overall")) {
  kind <- match.arg(kind)
  f <- get("calc_kinship", envir=asNamespace("qtl2"))

  # Try common named-argument patterns across qtl2 versions
  tries <- list(
    list(probs=genoprobs, type=kind),
    list(genoprobs=genoprobs, type=kind),
    list(probs=genoprobs, kind=kind),
    list(genoprobs=genoprobs, kind=kind),
    list(probs=genoprobs, method=kind),
    list(genoprobs=genoprobs, method=kind),
    list(probs=genoprobs, what=kind),
    list(genoprobs=genoprobs, what=kind)
  )
  for (a in tries) {
    res <- tryCatch(call_supported(f, a), error=function(e) NULL)
    if (!is.null(res)) return(res)
  }

  # Fallback: positional second argument
  res <- tryCatch(f(genoprobs, kind), error=function(e) NULL)
  if (!is.null(res)) return(res)

  stop("calc_kinship failed (kind=", kind, ")")
}

read_covar_tsv <- function(path, ids) {
  # returns numeric matrix with rownames=ids order; can return NULL
  if (is.null(path) || !nzchar(path) || !file.exists(path)) return(NULL)
  dt <- tryCatch(fread(path), error=function(e) NULL)
  if (is.null(dt) || ncol(dt) < 2) return(NULL)

  # Detect id column
  id_col <- NULL
  if ("id" %in% names(dt)) id_col <- "id"
  if (is.null(id_col) && "sample" %in% names(dt)) id_col <- "sample"
  if (is.null(id_col) && "IID" %in% names(dt)) id_col <- "IID"
  if (is.null(id_col)) id_col <- names(dt)[1]

  idv <- as.character(dt[[id_col]])
  dt[[id_col]] <- NULL

  # Keep only numeric covars
  num_cols <- names(dt)[vapply(dt, is.numeric, logical(1))]
  if (length(num_cols) == 0) return(NULL)
  dt <- dt[, ..num_cols]
  mat <- as.matrix(dt)
  rownames(mat) <- idv

  # Reorder + subset to ids (pheno rownames)
  mat <- mat[ids, , drop=FALSE]
  return(mat)
}

write_artifacts <- function(out_dir, artifacts) {
  writeLines(toJSON(artifacts, auto_unbox=TRUE, pretty=TRUE),
             con=file.path(out_dir, "artifacts.json"))
}

# -------- main --------
params_path <- get_arg("--params")
out_dir <- get_arg("--out")
if (is.null(params_path) || is.null(out_dir)) stop("Usage: --params params.json --out out_dir")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_dir, "run.log")
sink(log_file, split = TRUE)

cat("[qtl2_scan1] start\n")
cat("[qtl2_scan1] params_path=", params_path, "\n")
cat("[qtl2_scan1] out_dir=", out_dir, "\n")

p <- fromJSON(params_path)

cross2_rds <- if (!is.null(p$cross2_rds)) p$cross2_rds else NULL
if (is.null(cross2_rds) || !file.exists(cross2_rds)) stop("cross2_rds not found")

# trait(s): allow single string (trait) or comma-separated (traits)
traits_raw <- NULL
if (!is.null(p$traits) && nzchar(p$traits)) traits_raw <- p$traits
if (is.null(traits_raw) && !is.null(p$trait) && nzchar(p$trait)) traits_raw <- p$trait
if (is.null(traits_raw)) stop("trait(s) is required")
traits <- trimws(unlist(strsplit(traits_raw, "[,;]+", perl=TRUE)))
traits <- traits[nzchar(traits)]
if (length(traits) == 0) stop("trait(s) is required")

n_perm <- if (!is.null(p$n_perm)) as.integer(p$n_perm) else 100L
n_perm <- max(0L, n_perm)

alpha <- if (!is.null(p$alpha)) as.numeric(p$alpha) else 0.05
if (is.na(alpha) || alpha <= 0 || alpha >= 1) alpha <- 0.05

step <- if (!is.null(p$step)) as.numeric(p$step) else 1.0
if (is.na(step) || step <= 0) step <- 1.0

stepwidth <- p$stepwidth

error_prob <- if (!is.null(p$error_prob)) as.numeric(p$error_prob) else 0.0001
if (is.na(error_prob) || error_prob <= 0 || error_prob >= 0.5) error_prob <- 0.0001

# use_kinship <- if (!is.null(p$use_kinship)) as.logical(p$use_kinship) else FALSE
# if (is.na(use_kinship)) use_kinship <- FALSE
# kinship_type <- if (!is.null(p$kinship_type) && nzchar(p$kinship_type)) tolower(p$kinship_type) else "loco"
# if (!(kinship_type %in% c("loco", "overall"))) kinship_type <- "loco"
kinship_type <- p$kinship_type
trait_model <- p$model
map_function <- p$map_function
effect <- p$effect

method <- if (!is.null(p$method) && nzchar(p$method)) as.character(p$method) else NULL
n_cores <- if (!is.null(p$n_cores)) as.integer(p$n_cores) else 1L
n_cores <- max(1L, n_cores)

peak_drop <- if (!is.null(p$peak_drop)) as.numeric(p$peak_drop) else 1.5
if (is.na(peak_drop) || peak_drop <= 0) peak_drop <- 1.5

save_rds <- if (!is.null(p$save_rds)) as.logical(p$save_rds) else FALSE
if (is.na(save_rds)) save_rds <- FALSE

# Optional: coefficient / effect extraction (best-effort; never fails the run)
do_coef <- if (!is.null(p$do_coef)) as.logical(p$do_coef) else FALSE
if (is.na(do_coef)) do_coef <- FALSE
coef_max_peaks <- if (!is.null(p$coef_max_peaks)) as.integer(p$coef_max_peaks) else 3L
coef_max_peaks <- max(1L, min(50L, coef_max_peaks))

covar_tsv <- if (!is.null(p$covar_tsv) && nzchar(p$covar_tsv)) p$covar_tsv else NULL
scale_covar <- if (!is.null(p$scale_covar)) as.logical(p$scale_covar) else TRUE
if (is.na(scale_covar)) scale_covar <- TRUE

seed <- if (!is.null(p$seed)) as.integer(p$seed) else 1L
set.seed(seed)

cat("[qtl2_scan1] reading cross2\n")
if (grepl("\\.rds$", cross2_rds, perl = T)) {
  cross2 <- readRDS(cross2_rds)
} else {
  cross2 <- read_cross2(cross2_rds)
}
if (is.null(cross2$pheno)) stop("cross2$pheno is missing")

# validate traits
avail_traits <- colnames(cross2$pheno)
missing_traits <- setdiff(traits, avail_traits)
if (length(missing_traits) > 0) {
  cat("[qtl2_scan1] available phenotypes:\n")
  cat(paste(avail_traits, collapse=", "), "\n")
  stop(paste0("trait not found: ", paste(missing_traits, collapse=", ")))
}

# genotype probabilities
cat("[qtl2_scan1] genotype probs: step=", step, " error_prob=", error_prob, "\n")

map <- cross2$gmap

# Insert pseudomarkers (grid) if available
# if (exists("insert_pseudomarkers", where=asNamespace("qtl2"), inherits=FALSE)) {
#   f <- get("insert_pseudomarkers", envir=asNamespace("qtl2"))
#   cross2 <- tryCatch({
#     call_supported(f, list(map=map, step=step, stepwidth=stepwidth))
#   }, error=function(e) {
#     cat("[qtl2_scan1] WARN: insert_pseudomarkers failed: ", conditionMessage(e), "\n")
#     cross2
#   })
# }
pmap <- insert_pseudomarkers(map = map, 
                             step = step, 
                             off_end = 0,
                             stepwidth = stepwidth,
                             pseudomarker_map = NULL,
                             tol = 0.01,
                             cores = n_cores)


pr <- qtl2::calc_genoprob(cross = cross2, 
                          map = pmap,
                          error_prob = error_prob, 
                          map_function = map_function,
                          lowmem = FALSE,
                          quiet = TRUE,
                          cores = n_cores)

# Build map table (chr/marker/pos) for joining later
map_df_list <- lapply(names(map), function(chr) {
  data.frame(chr=chr, marker=names(map[[chr]]), pos=as.numeric(map[[chr]]), stringsAsFactors=FALSE)
})
map_df <- do.call(rbind, map_df_list)

# kinship
kin <- NULL
if (kinship_type != "none") {
  cat("[qtl2_scan1] calc_kinship type=", kinship_type, "\n")
  kin <- qtl2::calc_kinship(probs = pr,
                            type = kinship_type,
                            omit_x = F,
                            use_allele_probs = T,
                            quiet = F,
                            cores = n_cores)
}

# For each trait, run scan + write outputs
all_peaks_paths <- c()
all_plot_paths <- c()
all_lod_paths <- c()
all_perm_paths <- c()
all_coef_paths <- c()
all_effplot_paths <- c()

# store per-trait file mapping for GUI
trait_files <- list()

for (trait in traits) {
  cat("[qtl2_scan1] trait=", trait, "\n")

  ph <- cross2$pheno[, trait, drop=FALSE]
  # Ensure rownames exist
  if (is.null(rownames(ph))) {
    # try to set from ind IDs if present
    if (!is.null(cross2$geno) && length(cross2$geno) > 0) {
      rn <- rownames(cross2$geno[[1]])
      if (!is.null(rn) && nrow(ph) == length(rn)) rownames(ph) <- rn
    }
  }

  ids <- rownames(ph)
  addcovar <- NULL
  if (!is.null(ids)) {
    addcovar <- read_covar_tsv(covar_tsv, ids)
    if (!is.null(addcovar) && scale_covar) {
      addcovar <- scale(addcovar)
      addcovar <- as.matrix(addcovar)
    }
  }

  # scan1
  cat("[qtl2_scan1] scan1\n")
  out <- qtl2::scan1(genoprobs = pr,
                     pheno = ph,
                     kinship = kin,
                     addcovar = addcovar,
                     Xcovar = NULL,
                     intcovar = NULL,
                     weights = NULL,
                     reml = T,
                     model = trait_model,
                     hsq = NULL,
                     cores = n_cores)
  if (save_rds) {
    saveRDS(out, file=file.path(out_dir, paste0("scan1_", trait, ".rds")))
  }

  # permutation + threshold
  thr <- NA_real_
  perm_path <- file.path(out_dir, paste0("perm_thresholds_", trait, ".tsv"))
  perm_tbl <- data.frame(trait=trait, alpha=alpha, threshold=NA_real_)
  if (n_perm > 0) {
    cat("[qtl2_scan1] scan1perm: n_perm=", n_perm, " cores=", n_cores, "\n")
    perms <- qtl2::scan1perm(genoprobs = pr,
                             pheno = ph,
                             kinship = kin,
                             addcovar = addcovar,
                             Xcovar = NULL,
                             intcovar = NULL,
                             weights = NULL,
                             reml = TRUE,
                             model = trait_model,
                             n_perm = n_perm,
                             perm_Xsp = F,
                             perm_strata = NULL,
                             chr_lengths = NULL,
                             cores = n_cores)
    
    thr_vec <- summary(perms, alpha=alpha)
    thr <- as.numeric(thr_vec[1])
    perm_tbl$threshold <- thr
    fwrite(as.data.frame(thr_vec), file=perm_path, sep="\t")
  } else {
    fwrite(perm_tbl, file=perm_path, sep="\t")
  }

  # peaks
  peaks <- NULL
  if (!is.na(thr)) {
    peaks <- find_peaks(out, map, threshold=thr, drop=peak_drop)
  } else {
    peaks <- find_peaks(out, map, drop=peak_drop)
  }
  if (is.null(peaks) || nrow(peaks) == 0) {
    peaks <- data.frame(chr=character(), pos=numeric(), lod=numeric(), stringsAsFactors = FALSE)
  }
  peaks_path <- file.path(out_dir, paste0("peaks_", trait, ".tsv"))
  fwrite(peaks, file=peaks_path, sep="\t")

  # optional: coefficients / effect plot
  coef_path <- ""
  effplot_path <- ""
  if (effect != "none") {
    cat("[qtl2_scan1] calculate QTL effects \n")
    if (effect == "coef") {
      cf <- qtl2::scan1coef(genoprobs = pr,
                            pheno = ph,
                            kinship = kin,
                            addcovar = addcovar,
                            nullcovar = NULL,
                            intcovar = NULL,
                            weights = NULL,
                            contrasts = NULL,
                            model = trait_model,
                            zerosum = T,
                            se = F,
                            hsq = NULL,
                            reml = T)
    } else {
      cf <- qtl2::scan1blup(genoprobs = pr,
                            pheno = ph,
                            kinship = kin,
                            addcovar = addcovar,
                            nullcovar = NULL,
                            contrasts = NULL,
                            se = F,
                            reml = T,
                            tol = 0.000000000001,
                            cores = n_cores,
                            quiet = T)
    }
    
    # Convert to data.frame with marker names
    cf_mat <- NULL
    if (is.matrix(cf)) cf_mat <- cf
    if (is.array(cf) && length(dim(cf)) == 2) cf_mat <- as.matrix(cf)
    if (is.null(cf_mat)) stop("scan1coef returned non-2D object")
    cf_df <- as.data.frame(cf_mat)
    cf_df$marker <- rownames(cf_mat)
    if (is.null(cf_df$marker)) cf_df$marker <- rownames(as.data.frame(out))
    
    # pick top peaks
    pk <- peaks
    if (nrow(pk) > 0 && ("lod" %in% names(pk))) {
      pk <- pk[order(-pk$lod), , drop=FALSE]
    }
    pk <- head(pk, coef_max_peaks)
    
    # map peak pos -> nearest marker name
    peak_markers <- character(0)
    if (nrow(pk) > 0 && all(c("chr","pos") %in% names(pk))) {
      for (i in seq_len(nrow(pk))) {
        chr_i <- as.character(pk$chr[i])
        pos_i <- as.numeric(pk$pos[i])
        if (!chr_i %in% names(map)) next
        v <- as.numeric(map[[chr_i]])
        nm <- names(map[[chr_i]])
        if (length(v) == 0) next
        j <- which.min(abs(v - pos_i))
        peak_markers <- c(peak_markers, nm[j])
      }
    }
    peak_markers <- unique(peak_markers)
    if (length(peak_markers) == 0) {
      # fallback: take the top marker from scan1 output
      scan_df2 <- as.data.frame(out)
      if (ncol(scan_df2) == 1) colnames(scan_df2) <- "lod"
      scan_df2$marker <- rownames(scan_df2)
      scan_df2 <- scan_df2[order(-scan_df2[[1]]), , drop=FALSE]
      peak_markers <- head(scan_df2$marker, 1)
    }
    
    #sub_df <- cf_df[cf_df$marker %in% peak_markers, , drop=FALSE]
    #if (nrow(sub_df) == 0) sub_df <- cf_df[1, , drop=FALSE]
    sub_df <- data.frame(marker = rownames(cf),
                         cf)
    peak_marker <- which.max(sub_df$ac1)
    
    coef_path <- file.path(out_dir, paste0("scan1coef_peaks_", trait, ".tsv"))
    fwrite(sub_df, file=coef_path, sep="\t", quote = F)
    
    # effect plot for the first selected marker
    effplot_path <- file.path(out_dir, paste0("effect_plot_", trait, ".png"))
    png(effplot_path, width=1100, height=700)
    m0 <- sub_df$marker[1]
    y <- as.numeric(sub_df[1, setdiff(names(sub_df), "marker")])
    labs <- setdiff(names(sub_df), "marker")
    # m0 <- sub_df$marker[peak_marker]
    # if (effect == "coef") {
    #   y <- as.numeric(sub_df[peak_marker, 2:4])
    #   labs <- colnames(sub_df)[2:4]
    # } else {
    #   y <- as.numeric(sub_df[peak_marker, 2:3])
    #   labs <- colnames(sub_df)[2:3]
    # }
    # stop(labs)
    
    op <- par(mar=c(9,5,4,2))
    barplot(y, names.arg=labs, las=2, main=paste0("Allele effects (scan1coef): ", trait, " @ ", m0), ylab="effect")
    par(op)
    dev.off()
  }

  # LOD profile (join with map_df via marker names)
  scan_df <- as.data.frame(out)
  if (ncol(scan_df) == 1) colnames(scan_df) <- "lod"
  scan_df$marker <- rownames(scan_df)
  
  lod_profile <- merge(map_df, scan_df, by="marker", all.x=FALSE, all.y=TRUE)
  #lod_profile <- na.omit(lod_profile)
  lod_path <- file.path(out_dir, paste0("lod_profile_", trait, ".tsv"))
  fwrite(lod_profile, file=lod_path, sep="\t", quote = F)

  # plot
  png_path <- file.path(out_dir, paste0("lod_plot_", trait, ".png"))
  cat("[qtl2_scan1] writing plot: ", png_path, "\n")
  png(png_path, width=1200, height=800)
  try({
    plot(out, map, main=paste0("qtl2 scan1: ", trait))
    if (!is.na(thr)) abline(h=thr, lty=2)
  }, silent=TRUE)
  dev.off()

  all_peaks_paths <- c(all_peaks_paths, basename(peaks_path))
  all_plot_paths <- c(all_plot_paths, basename(png_path))
  all_lod_paths <- c(all_lod_paths, basename(lod_path))
  all_perm_paths <- c(all_perm_paths, basename(perm_path))

  if (nzchar(coef_path)) all_coef_paths <- c(all_coef_paths, basename(coef_path))
  if (nzchar(effplot_path)) all_effplot_paths <- c(all_effplot_paths, basename(effplot_path))

  trait_files[[trait]] <- list(
    peaks = basename(peaks_path),
    plot_lod = basename(png_path),
    lod_profile = basename(lod_path),
    perm_threshold = basename(perm_path),
    scan1coef_peaks = if (nzchar(coef_path)) basename(coef_path) else "",
    plot_effect = if (nzchar(effplot_path)) basename(effplot_path) else ""
  )
}

# GUI artifacts
# Use the first trait as default for single-view panes, but keep lists for later UI extension.
default_trait <- traits[1]
art <- list(
  traits = traits,
  default_trait = default_trait,
  peaks = all_peaks_paths,
  plots = all_plot_paths,
  lod_profiles = all_lod_paths,
  perm_thresholds = all_perm_paths,
  scan1coef_peaks = all_coef_paths,
  effect_plots = all_effplot_paths,
  trait_files = trait_files,
  table = paste0("peaks_", default_trait, ".tsv"),
  plot = paste0("lod_plot_", default_trait, ".png"),
  lod_profile = paste0("lod_profile_", default_trait, ".tsv"),
  perm_threshold = paste0("perm_thresholds_", default_trait, ".tsv")
)
write_artifacts(out_dir, art)

cat("[qtl2_scan1] done\n")
sink()
