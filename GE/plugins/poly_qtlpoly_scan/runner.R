suppressWarnings(suppressMessages({
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("jsonlite is required")
  library(qtlpoly)
}))

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) return(args[[i + 1]])
  default
}

params_path <- get_arg("--params")
out_dir <- get_arg("--out")
if (is.null(params_path) || is.null(out_dir)) {
  stop("Usage: Rscript runner.R --params params.json --out out_dir")
}
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

logf <- file.path(out_dir, "run.log")
log <- function(...) {
  msg <- paste0("[", Sys.time(), "] ", paste0(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = logf, append = TRUE)
}
write_error <- function(msg) {
  writeLines(as.character(msg), con = file.path(out_dir, "error_message.txt"), useBytes = TRUE)
}
write_artifacts <- function(meta) {
  p <- file.path(out_dir, "artifacts.json")
  writeLines(jsonlite::toJSON(meta, auto_unbox = TRUE, pretty = TRUE), con = p, useBytes = TRUE)
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

read_tsv <- function(path) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  } else {
    data.table::fread(path, sep = "\t", data.table = FALSE, check.names = FALSE)
  }
}

as_num <- function(x, default = NA_real_) {
  if (is.null(x)) return(default)
  suppressWarnings(as.numeric(x))
}

parse_traits <- function(x) {
  if (is.null(x)) return(character())
  if (is.character(x) && length(x) == 1L) {
    x <- trimws(unlist(strsplit(x, ",")))
  }
  x <- as.character(x)
  x <- x[nzchar(x)]
  unique(x)
}

get_ind_names <- function(g) {
  if (!is.null(g$ind.names)) return(as.character(g$ind.names))
  if (!is.null(g$ind.names2)) return(as.character(g$ind.names2))
  p <- g$probs %||% NULL
  if (!is.null(p) && length(dim(p)) == 3L) {
    dn <- dimnames(p)
    if (!is.null(dn) && length(dn) >= 3L && !is.null(dn[[3]])) return(as.character(dn[[3]]))
  }
  NULL
}

infer_map_cols <- function(df) {
  nms <- names(df)
  pick <- function(cands) {
    hit <- cands[cands %in% nms]
    if (length(hit) > 0) hit[[1]] else NA_character_
  }
  list(
    chr = pick(c("chr", "LG", "linkage_group", "group")),
    cm  = pick(c("pos_cM", "pos_cm", "cM", "cm", "pos_cm", "pos_cM", "pos")),
    bp  = pick(c("pos_bp", "bp", "pos_bp", "pos_physical", "pos_bp", "physical_pos"))
  )
}

interp_bp <- function(map_df, lg, cm_vec) {
  if (nrow(map_df) == 0) return(rep(NA_real_, length(cm_vec)))
  cols <- infer_map_cols(map_df)
  if (is.na(cols$chr) || is.na(cols$cm) || is.na(cols$bp)) return(rep(NA_real_, length(cm_vec)))
  sub <- map_df[map_df[[cols$chr]] == lg, , drop = FALSE]
  if (nrow(sub) < 2) return(rep(NA_real_, length(cm_vec)))
  x <- suppressWarnings(as.numeric(sub[[cols$cm]]))
  y <- suppressWarnings(as.numeric(sub[[cols$bp]]))
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  if (length(x) < 2) return(rep(NA_real_, length(cm_vec)))
  o <- order(x)
  x <- x[o]; y <- y[o]
  approx(x = x, y = y, xout = cm_vec, rule = 2)$y
}

get_trait_result <- function(remim_mod, pheno_col_idx, trait_label) {
  res <- remim_mod$results[[as.character(pheno_col_idx)]] %||% NULL
  if (is.null(res)) {
    # Fallback: try by name
    if (!is.null(names(remim_mod$results)) && trait_label %in% names(remim_mod$results)) {
      res <- remim_mod$results[[trait_label]]
    }
  }
  res
}

extract_profile_lop <- function(remim_mod, pheno_col_idx, trait_label) {
  res <- get_trait_result(remim_mod, pheno_col_idx, trait_label)
  if (is.null(res)) return(NULL)
  pval <- res$pval %||% NULL
  if (is.null(pval)) return(NULL)
  pv <- suppressWarnings(as.numeric(pval))
  pv <- pmax(pv, .Machine$double.xmin)
  lop <- -log10(pv)
  lop[!is.finite(lop)] <- NA_real_
  lop
}

extract_qtls_table <- function(remim_mod, data_obj, pheno_col_idx, trait_label) {
  res <- get_trait_result(remim_mod, pheno_col_idx, trait_label)
  if (is.null(res)) return(NULL)
  qtls <- res$qtls %||% NULL
  if (is.null(qtls) || nrow(qtls) == 0) return(data.frame())
  lower <- res$lower %||% NULL
  upper <- res$upper %||% NULL
  df <- as.data.frame(qtls, stringsAsFactors = FALSE)
  df$trait <- trait_label
  
  # Normalize column names
  # Expected: LG, Pos, Nmrk, Mrk, Score, Pval
  if (!"LG" %in% names(df) && "Chr" %in% names(df)) df$LG <- df$Chr
  if (!"Pos" %in% names(df) && "pos" %in% names(df)) df$Pos <- df$pos
  if (!"Pval" %in% names(df) && "pval" %in% names(df)) df$Pval <- df$pval
  
  df$Pval_num <- suppressWarnings(as.numeric(df$Pval))
  df$LOP <- -log10(df$Pval_num)
  df$chr <- df$LG
  df$pos_cM <- suppressWarnings(as.numeric(df$Pos))
  df$pos_unit <- "cM"
  
  # Compatibility columns for the GUI's generic QTL Plotly panel
  # (build_qtl_scan_figure expects chr/pos/lod and optional CI columns).
  df$pos <- df$pos_cM
  df$lod <- df$LOP
  
  if (!is.null(lower) && nrow(lower) == nrow(df)) {
    lo <- as.data.frame(lower, stringsAsFactors = FALSE)
    if (!"Pos_lower" %in% names(lo) && "Pos" %in% names(lo)) lo$Pos_lower <- lo$Pos
    df$pos_lower_cM <- suppressWarnings(as.numeric(lo$Pos_lower))
    df$Mrk_lower <- lo$Mrk_lower %||% NA_character_
    df$Score_lower <- lo$Score_lower %||% NA_real_
    df$Pval_lower <- lo$Pval_lower %||% NA_character_
  } else {
    df$pos_lower_cM <- NA_real_
  }
  if (!is.null(upper) && nrow(upper) == nrow(df)) {
    up <- as.data.frame(upper, stringsAsFactors = FALSE)
    if (!"Pos_upper" %in% names(up) && "Pos" %in% names(up)) up$Pos_upper <- up$Pos
    df$pos_upper_cM <- suppressWarnings(as.numeric(up$Pos_upper))
    df$Mrk_upper <- up$Mrk_upper %||% NA_character_
    df$Score_upper <- up$Score_upper %||% NA_real_
    df$Pval_upper <- up$Pval_upper %||% NA_character_
  } else {
    df$pos_upper_cM <- NA_real_
  }
  
  # CI columns expected by the GUI helper
  df$ci_lower <- df$pos_lower_cM
  df$ci_upper <- df$pos_upper_cM
  df
}

log("start")
log("params_path=", params_path)
log("out_dir=", out_dir)

params <- tryCatch(jsonlite::fromJSON(params_path), error = function(e) list())

if (!requireNamespace("qtlpoly", quietly = TRUE)) {
  msg <- "R package 'qtlpoly' is required but not installed in this environment. Try: install.packages('qtlpoly')"
  write_error(msg)
  write_artifacts(list(module = "poly_qtlpoly_scan", error = msg))
  quit(status = 1)
}

qtlpoly_export_rds <- as.character(params$qtlpoly_export_rds %||% params$geno_prob_rds %||% "")
pheno_tsv <- as.character(params$pheno_tsv %||% "")
map_markers_tsv <- as.character(params$map_markers_tsv %||% "")

if (!nzchar(qtlpoly_export_rds) || !file.exists(qtlpoly_export_rds)) {
  msg <- paste0("qtlpoly_export_rds not found: ", qtlpoly_export_rds)
  write_error(msg)
  write_artifacts(list(module = "poly_qtlpoly_scan", error = msg))
  quit(status = 1)
}

# if (!nzchar(pheno_tsv) || !file.exists(pheno_tsv)) {
#   msg <- paste0("pheno_tsv not found: ", pheno_tsv)
#   write_error(msg)
#   write_artifacts(list(module = "poly_qtlpoly_scan", error = msg))
#   quit(status = 1)
# }

ploidy <- as.integer(params$ploidy %||% 4L)
step_cM <- as_num(params$step_cM %||% params$step %||% 1, 1)
w_size <- as.integer(params$w_size_cM %||% params$w.size %||% 15L)
d_sint <- as_num(params$d_sint %||% 1.5, 1.5)
n_cores <- as.integer(params$n_cores %||% params$n.clusters %||% 1L)

threshold_mode <- as.character(params$threshold_mode %||% "pointwise_manual")
sig_fwd <- as_num(params$sig_fwd %||% 0.01, 0.01)
sig_bwd <- as_num(params$sig_bwd %||% 1e-4, 1e-4)
alpha_fwd <- as_num(params$alpha_fwd %||% 0.2, 0.2)
alpha_bwd <- as_num(params$alpha_bwd %||% 0.05, 0.05)
n_sim <- as.integer(params$n_sim %||% 200L)
seed <- as.integer(params$seed %||% 1L)

plot_grid <- isTRUE(params$plot_grid %||% FALSE)
plot_supint <- isTRUE(params$plot_supint %||% FALSE)
fit_model <- isTRUE(params$fit_model %||% FALSE)
do_effects <- isTRUE(params$effects %||% FALSE)

id_col <- as.character(params$id_col %||% "id")
trait_list <- parse_traits(params$trait %||% params$traits %||% character())
if (length(trait_list) == 0) trait_list <- character()

log("ploidy=", ploidy, ", step_cM=", step_cM, ", w_size=", w_size, ", d_sint=", d_sint, ", n_cores=", n_cores)
log("threshold_mode=", threshold_mode)

# Load genotype probabilities exported for qtlpoly
geno_prob <- readRDS(qtlpoly_export_rds)
if (!is.list(geno_prob) || length(geno_prob) == 0) {
  msg <- "Invalid qtlpoly_export_rds: expected a non-empty list. (Tip: use mappoly::export_qtlpoly output)"
  write_error(msg)
  write_artifacts(list(module = "poly_qtlpoly_scan", error = msg))
  quit(status = 1)
}
#stop("!")
# Load phenotype and align to individuals
if (nchar(pheno_tsv) < 1) {
  ph <- geno_prob$pheno
} else {
  #ph <- read_tsv(pheno_tsv)
  ph <- data.frame(data.table::fread(pheno_tsv, header = T),
                   row.names = 1)
  if (nrow(ph) == 0) {
    msg <- "Empty phenotype table"
    write_error(msg)
    write_artifacts(list(module = "poly_qtlpoly_scan", error = msg))
    quit(status = 1)
  }
}

# if (id_col %in% colnames(ph)) {
#   rn <- as.character(ph[[id_col]])
#   ph[[id_col]] <- NULL
# } else {
#   rn <- as.character(ph[[1]])
#   ph[[1]] <- NULL
# }
# rownames(ph) <- rn

if (ncol(ph) == 0) {
  msg <- "Phenotype table has no trait columns after removing the ID column"
  write_error(msg)
  write_artifacts(list(module = "poly_qtlpoly_scan", error = msg))
  quit(status = 1)
}

ind_names <- NULL
if (nchar(pheno_tsv) < 1) {
  ind_names <- geno_prob$ind.names
} else {
  for (g in geno_prob) {
    ind_names <- get_ind_names(g)
    if (!is.null(ind_names)) break
  }
}

if (!is.null(ind_names)) {
  common <- intersect(ind_names, rownames(ph))
  if (length(common) == 0) {
    msg <- "No overlapping individual IDs between genotype probabilities and phenotype table"
    write_error(msg)
    write_artifacts(list(module = "poly_qtlpoly_scan", error = msg))
    quit(status = 1)
  }
  ph <- ph[common, , drop = FALSE]
  # reorder to match genoprob object
  ph <- ph[ind_names[ind_names %in% rownames(ph)], , drop = FALSE]
}

# Decide which traits to analyze
if (length(trait_list) > 0) {
  missing_traits <- setdiff(trait_list, colnames(ph))
  if (length(missing_traits) > 0) {
    msg <- paste0("Trait(s) not found in phenotype table: ", paste(missing_traits, collapse = ", "))
    write_error(msg)
    write_artifacts(list(module = "poly_qtlpoly_scan", error = msg))
    quit(status = 1)
  }
}

pheno_col_idx <- if (length(trait_list) == 0) seq_len(ncol(ph)) else match(trait_list, colnames(ph))

# Build qtlpoly data object
log("Building qtlpoly.data")
if (nchar(pheno_tsv) < 1) {
  data_obj <- geno_prob
} else {
  data_obj <- qtlpoly::read_data(ploidy = ploidy, 
                                 geno.prob = geno_prob, 
                                 pheno = ph, 
                                 step = step_cM)
}
# Optional score.null computation
score_null <- NULL
if (threshold_mode %in% c("resampling_genomewide", "resampling_pointwise")) {
  log("Computing score.null via simulate_qtl + null_model: n_sim=", n_sim, ", seed=", seed)
  set.seed(seed)
  sim <- qtlpoly::simulate_qtl(data = data_obj, 
                               mu = 0, 
                               h2.qtl = NULL, 
                               var.error = 1,
                               n.sim = n_sim, 
                               missing = TRUE, 
                               seed = seed)
  score_null <- qtlpoly::null_model(data = sim$results)
  # score_null <- qtlpoly::null_model(data = data_obj,
  #                                   offset.data = NULL,
  #                                   pheno.col = pheno_col_idx,
  #                                   n.clusters = n_cores,
  #                                   plot = "null",
  #                                   verbose = T)
}

# If requested, derive pointwise thresholds from score.null (min p-value quantiles)
if (threshold_mode == "resampling_pointwise") {
  # score_null contains simulated p-values per trait. For multi-trait runs,
  # we use the most stringent (smallest) threshold across selected traits.
  sig_fwd_vec <- c(); sig_bwd_vec <- c()
  for (k in seq_along(pheno_col_idx)) {
    pvals <- score_null$results[[k]]$pval %||% NULL
    if (is.null(pvals)) next
    min_p <- apply(pvals, 2, min, na.rm = TRUE)
    sig_fwd_vec <- c(sig_fwd_vec, as.numeric(stats::quantile(min_p, probs = alpha_fwd, na.rm = TRUE)))
    sig_bwd_vec <- c(sig_bwd_vec, as.numeric(stats::quantile(min_p, probs = alpha_bwd, na.rm = TRUE)))
  }
  if (length(sig_fwd_vec) == 0 || length(sig_bwd_vec) == 0) {
    msg <- "score.null did not contain expected pval matrices"
    write_error(msg)
    write_artifacts(list(module = "poly_qtlpoly_scan", error = msg))
    quit(status = 1)
  }
  sig_fwd <- min(sig_fwd_vec, na.rm = TRUE)
  sig_bwd <- min(sig_bwd_vec, na.rm = TRUE)
  log("Derived pointwise thresholds from resampling (min across traits): sig_fwd=", sig_fwd, ", sig_bwd=", sig_bwd)
  score_null <- NULL
}

use_score_null <- (threshold_mode == "resampling_genomewide")
if (use_score_null) {
  sig_fwd <- alpha_fwd
  sig_bwd <- alpha_bwd
}

# Run REMIM
log("Running remim")
remim_mod <- qtlpoly::remim(
  data = data_obj,
  pheno.col = pheno_col_idx,
  w.size = w_size,
  sig.fwd = sig_fwd,
  sig.bwd = sig_bwd,
  score.null = if (use_score_null) score_null else NULL,
  d.sint = d_sint,
  polygenes = F,
  n.clusters = n_cores,
  n.rounds = Inf,
  plot = "remim",
  verbose = T
)

saveRDS(remim_mod, file.path(out_dir, "remim_model.rds"))

# Export peaks table
trait_labels <- colnames(ph)[pheno_col_idx]
all_peaks <- list()
for (k in seq_along(pheno_col_idx)) {
  col_i <- pheno_col_idx[[k]]
  tr <- trait_labels[[k]]
  df <- extract_qtls_table(remim_mod, data_obj, col_i, tr)
  if (is.null(df)) next
  all_peaks[[length(all_peaks) + 1L]] <- df
}
peaks <- if (length(all_peaks) == 0) data.frame() else do.call(rbind, all_peaks)

meta <- list(
  module = "poly_qtlpoly_scan",
  ploidy = ploidy,
  step_cM = step_cM,
  traits = trait_labels,
  threshold_mode = threshold_mode,
  sig_fwd = sig_fwd,
  sig_bwd = sig_bwd,
  peaks_tsv = "peaks.tsv",
  lod_profile_tsv = "lod_profile.tsv",
  thresholds_tsv = "thresholds.tsv",
  plot = "lod_plot.png",
  default_table = "peaks.tsv",
  remim_model_rds = "remim_model.rds"
)

if (file.exists(file.path(out_dir, "fitted_model.rds"))) meta$fitted_model_rds <- "fitted_model.rds"
if (file.exists(file.path(out_dir, "qtl_effects.tsv"))) meta$qtl_effects_tsv <- "qtl_effects.tsv"

write_artifacts(meta)

# Optional interpolation to bp, if user provides a map table with both cM and bp
map_df <- NULL
if (nzchar(map_markers_tsv) && file.exists(map_markers_tsv)) {
  map_df <- tryCatch(read_tsv(map_markers_tsv), error = function(e) NULL)
}
if (!is.null(map_df) && nrow(peaks) > 0) {
  cols <- infer_map_cols(map_df)
  if (!is.na(cols$chr) && !is.na(cols$cm) && !is.na(cols$bp)) {
    peaks$pos_bp <- NA_real_
    peaks$pos_lower_bp <- NA_real_
    peaks$pos_upper_bp <- NA_real_
    for (lg in unique(peaks$chr)) {
      idx <- which(peaks$chr == lg)
      peaks$pos_bp[idx] <- interp_bp(map_df, lg, peaks$pos_cM[idx])
      peaks$pos_lower_bp[idx] <- interp_bp(map_df, lg, peaks$pos_lower_cM[idx])
      peaks$pos_upper_bp[idx] <- interp_bp(map_df, lg, peaks$pos_upper_cM[idx])
    }
    peaks$pos_bp_unit <- "bp"
  }
}

peaks_out <- file.path(out_dir, "peaks.tsv")
utils::write.table(peaks, peaks_out, sep = "\t", quote = FALSE, row.names = FALSE)

# Export thresholds (as -log10(p) lines for Plotly reuse)
thr_out <- file.path(out_dir, "thresholds.tsv")
thr_df <- data.frame(
  alpha = c(sig_fwd, sig_bwd),
  threshold = c(-log10(sig_fwd), -log10(sig_bwd)),
  stringsAsFactors = FALSE
)
utils::write.table(thr_df, thr_out, sep = "\t", quote = FALSE, row.names = FALSE)

# Make plots
lod_png <- file.path(out_dir, "lod_plot.png")
sint_png <- file.path(out_dir, "sint_plot.png")
qtl_png <- file.path(out_dir, "qtl_plot.png")

log("Plotting profile")
gp <- qtlpoly::plot_profile(data = data_obj, model = remim_mod, 
                            pheno.col = pheno_col_idx,
                            grid = isTRUE(plot_grid), sup.int = isTRUE(plot_supint))

# Export profile points for Plotly (best-effort via ggplot build)
# This enables the GUI to reuse the QTL Plotly panel (genome-wide + peak zoom).
# Output schema: chr, pos, lod (where lod == -log10(p)).
if (requireNamespace("ggplot2", quietly = TRUE)) {
  prof_out <- file.path(out_dir, "lod_profile.tsv")
  prof_ok <- FALSE
  try({
    b <- ggplot2::ggplot_build(gp)
    # pick the layer that most likely corresponds to the profile curve.
    # Some plots include helper layers (e.g. y=0 baselines or shaded ribbons)
    # that can have as many points as the true curve; we prefer layers with
    # varying positive y-values, then many unique x positions.
    best_i <- NA_integer_
    best_score <- -Inf
    for (i in seq_along(b$data)) {
      di <- b$data[[i]]
      if (is.data.frame(di) && all(c("x", "y", "PANEL") %in% names(di))) {
        xx <- suppressWarnings(as.numeric(di$x))
        yy <- suppressWarnings(as.numeric(di$y))
        xx <- xx[is.finite(xx)]
        yy <- yy[is.finite(yy)]
        if (length(xx) < 2 || length(yy) < 2) next
        uniq_x <- length(unique(round(xx, 8)))
        y_sd <- suppressWarnings(stats::sd(yy))
        y_max <- suppressWarnings(max(yy, na.rm = TRUE))
        varying <- is.finite(y_sd) && y_sd > 0
        positive <- is.finite(y_max) && y_max > 0
        score <- (if (varying) 1e9 else 0) + (if (positive) 1e8 else 0) + uniq_x * 1e3 + nrow(di)
        if (is.finite(score) && score > best_score) {
          best_score <- score
          best_i <- i
        }
      }
    }
    if (!is.na(best_i)) {
      di <- b$data[[best_i]]
      lay <- b$layout$layout
      # Try to find the facet variable that corresponds to linkage group.
      skip_cols <- c("PANEL", "ROW", "COL", "SCALE_X", "SCALE_Y", "AXIS_X", "AXIS_Y")
      facet_cols <- setdiff(names(lay), skip_cols)
      lg_col <- NULL
      if ("LG" %in% facet_cols) {
        lg_col <- "LG"
      } else if ("chr" %in% facet_cols) {
        lg_col <- "chr"
      } else if ("Chr" %in% facet_cols) {
        lg_col <- "Chr"
      } else if (length(facet_cols) >= 1) {
        lg_col <- facet_cols[[1]]
      }
      if (!is.null(lg_col) && lg_col %in% names(lay)) {
        map_panel <- lay[, c("PANEL", lg_col), drop = FALSE]
        names(map_panel) <- c("PANEL", "chr")
        dd <- merge(di, map_panel, by = "PANEL", all.x = TRUE)
      } else {
        dd <- di
        dd$chr <- NA_character_
      }
      out <- data.frame(
        chr = as.character(dd$chr),
        pos = suppressWarnings(as.numeric(dd$x)),
        lod = suppressWarnings(as.numeric(dd$y)),
        stringsAsFactors = FALSE
      )
      
      # Prefer the actual REMIM p-value profile when available. This avoids
      # accidentally exporting a zero-valued helper layer from ggplot.
      lop_true <- extract_profile_lop(remim_mod, pheno_col_idx[[1]], trait_labels[[1]])
      if (!is.null(lop_true)) {
        if (length(lop_true) == nrow(out)) {
          out$lod <- lop_true
        } else if (length(lop_true) > 0 && abs(length(lop_true) - nrow(out)) <= 2L) {
          n_keep <- min(length(lop_true), nrow(out))
          out <- out[seq_len(n_keep), , drop = FALSE]
          out$lod <- lop_true[seq_len(n_keep)]
        }
      }
      
      out <- out[is.finite(out$pos) & is.finite(out$lod), , drop = FALSE]
      if (nrow(out) > 0) {
        utils::write.table(out, prof_out, sep = "	", quote = FALSE, row.names = FALSE)
        prof_ok <- TRUE
      }
    }
  }, silent = TRUE)
  if (!prof_ok) {
    # Ensure the file exists to avoid GUI confusion.
    utils::write.table(data.frame(chr=character(), pos=numeric(), lod=numeric()), prof_out, sep = "\t", quote = FALSE, row.names = FALSE)
  }
}

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  # ggsave requires ggplot2; fallback to base png capture
  png(lod_png, width = 1400, height = 800)
  print(gp)
  dev.off()
} else {
  ggplot2::ggsave(filename = lod_png, plot = gp, width = 12, height = 6, dpi = 150)
}

# Support intervals plot (best-effort; may be empty when no QTL)
try({
  log("Plotting support intervals")
  p2 <- qtlpoly::plot_sint(data = data_obj, model = remim_mod)
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    png(sint_png, width = 1400, height = 800)
    print(p2)
    dev.off()
  } else {
    ggplot2::ggsave(filename = sint_png, plot = p2, width = 12, height = 6, dpi = 150)
  }
}, silent = TRUE)

if (file.exists(sint_png)) meta$sint_plot <- basename(sint_png)
if (file.exists(qtl_png)) meta$qtl_plot <- basename(qtl_png)

fitted_mod <- NULL
if (isTRUE(fit_model) && nrow(peaks) > 0) {
  log("Fitting REML model")
  fitted_mod <- qtlpoly::fit_model(data = data_obj, model = remim_mod)
  saveRDS(fitted_mod, file.path(out_dir, "fitted_model.rds"))
  log("Plotting QTL heritability/significance")
  p3 <- qtlpoly::plot_qtl(data = data_obj, model = remim_mod, fitted = fitted_mod)
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    png(qtl_png, width = 1400, height = 800)
    print(p3)
    dev.off()
  } else {
    ggplot2::ggsave(filename = qtl_png, plot = p3, width = 12, height = 6, dpi = 150)
  }
}

if (isTRUE(do_effects) && !is.null(fitted_mod)) {
  log("Estimating allele effects")
  eff <- qtlpoly::qtl_effects(data = data_obj, model = remim_mod, fitted = fitted_mod)
  saveRDS(eff, file.path(out_dir, "qtl_effects.rds"))
  # Export a lightweight summary table if available
  if (!is.null(eff$effects)) {
    try({
      ef <- eff$effects
      utils::write.table(ef, file.path(out_dir, "qtl_effects.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    }, silent = TRUE)
  }
}

log("done")
