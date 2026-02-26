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

cat("[QTL] start\n")
cat("[QTL] params_path=", params_path, "\n")
cat("[QTL] out_dir=", out_dir, "\n")

p <- fromJSON(params_path)

cross_rds <- if (!is.null(p$cross_rds)) p$cross_rds else NULL
if (is.null(cross_rds) || !file.exists(cross_rds)) stop("cross_rds not found")

trait <- if (!is.null(p$trait) && nchar(p$trait) > 0) p$trait else NULL
#if (is.null(trait)) stop("trait is required")

analysis_mode <- if (!is.null(p$analysis_mode) && nchar(p$analysis_mode) > 0) tolower(p$analysis_mode) else "scanone"
if (!(analysis_mode %in% c("scanone", "cim", "mqm", "stepwise"))) analysis_mode <- "scanone"

n_perm <- if (!is.null(p$n_perm)) as.integer(p$n_perm) else 100L
n_perm <- max(0L, n_perm)

alpha <- if (!is.null(p$alpha)) as.numeric(p$alpha) else 0.05
if (is.na(alpha) || alpha <= 0 || alpha >= 1) alpha <- 0.05

step <- if (!is.null(p$step)) as.numeric(p$step) else 1.0
if (is.na(step) || step <= 0) step <- 1.0

error_prob <- if (!is.null(p$error_prob)) as.numeric(p$error_prob) else 0.001
if (is.na(error_prob) || error_prob <= 0 || error_prob >= 0.5) error_prob <- 0.001

map_function <- if (!is.null(p$map_function) && nchar(p$map_function) > 0) tolower(p$map_function) else "kosambi"
if (!(map_function %in% c("kosambi", "haldane"))) map_function <- "kosambi"

## shared window/n_marcovar (CIM/MQM)
n_marcovar <- if (!is.null(p$n_marcovar)) as.integer(p$n_marcovar) else 10L
n_marcovar <- max(0L, n_marcovar)
window <- if (!is.null(p$window)) as.numeric(p$window) else 10.0
if (is.na(window) || window < 0) window <- 10.0

## cofactor selection
cofactor_mode <- if (!is.null(p$cofactor_mode) && nchar(p$cofactor_mode) > 0) tolower(p$cofactor_mode) else "auto"
if (!(cofactor_mode %in% c("auto", "top", "file"))) cofactor_mode <- "auto"
cofactor_file <- if (!is.null(p$cofactor_file) && nchar(p$cofactor_file) > 0) p$cofactor_file else NULL
cofactor_lod_threshold <- if (!is.null(p$cofactor_lod_threshold)) as.numeric(p$cofactor_lod_threshold) else 0.0
if (is.na(cofactor_lod_threshold)) cofactor_lod_threshold <- 0.0
cofactor_min_dist <- if (!is.null(p$cofactor_min_dist)) as.numeric(p$cofactor_min_dist) else 10.0
if (is.na(cofactor_min_dist) || cofactor_min_dist < 0) cofactor_min_dist <- 10.0

## stability
stabilize <- if (!is.null(p$stabilize)) as.logical(p$stabilize) else TRUE
if (is.na(stabilize)) stabilize <- TRUE
stabilize_max_iter <- if (!is.null(p$stabilize_max_iter)) as.integer(p$stabilize_max_iter) else 5L
stabilize_max_iter <- max(1L, min(10L, stabilize_max_iter))

## stepwise params
stepwise_max_qtl <- if (!is.null(p$stepwise_max_qtl)) as.integer(p$stepwise_max_qtl) else 5L
stepwise_max_qtl <- max(1L, min(20L, stepwise_max_qtl))
stepwise_penalty_mode <- if (!is.null(p$stepwise_penalty_mode) && nchar(p$stepwise_penalty_mode) > 0) tolower(p$stepwise_penalty_mode) else "from_perm"
if (!(stepwise_penalty_mode %in% c("from_perm", "manual"))) stepwise_penalty_mode <- "from_perm"
stepwise_penalty_main <- if (!is.null(p$stepwise_penalty_main)) as.numeric(p$stepwise_penalty_main) else NA_real_

set.seed(if (!is.null(p$seed)) as.integer(p$seed) else 1L)

with_warn_capture <- function(expr) {
  warns <- character(0)
  val <- withCallingHandlers(
    expr,
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(val = val, warns = warns)
}

has_overspec_warn <- function(warns) {
  if (length(warns) == 0) return(FALSE)
  any(grepl("addcovar appears to be over-specified", warns, fixed = TRUE))
}

build_marker_map <- function(cross) {
  mp <- qtl::pull.map(cross)
  out <- list()
  for (chr in names(mp)) {
    v <- mp[[chr]]
    if (length(v) == 0) next
    out[[chr]] <- data.frame(chr = chr, marker = names(v), pos = as.numeric(v), stringsAsFactors = FALSE)
  }
  if (length(out) == 0) return(data.frame(chr=character(), marker=character(), pos=numeric()))
  do.call(rbind, out)
}

marker_lod_scan <- function(cross, trait) {
  ph <- qtl::pull.pheno(cross)
  y_all <- ph[[trait]]
  if (is.null(y_all)) stop("trait not found")
  G <- qtl::pull.geno(cross)
  if (is.null(G) || ncol(G) == 0) stop("No marker genotypes found in cross")
  map_df <- build_marker_map(cross)
  markers <- colnames(G)
  lods <- rep(NA_real_, length(markers))
  for (j in seq_along(markers)) {
    g <- G[, j]
    ok <- !is.na(g) & !is.na(y_all)
    if (sum(ok) < 10) next
    y <- y_all[ok]
    g2 <- g[ok]
    if (length(unique(g2)) < 2) next
    rss0 <- sum((y - mean(y))^2)
    fit <- stats::lm(y ~ factor(g2))
    rss1 <- sum(stats::residuals(fit)^2)
    if (is.finite(rss0) && is.finite(rss1) && rss1 > 0) {
      n <- length(y)
      lods[j] <- (n/2.0) * log10(rss0 / rss1)
    }
  }
  df <- data.frame(marker = markers, lod = lods, stringsAsFactors = FALSE)
  if (nrow(map_df) > 0) {
    df <- merge(df, map_df, by = "marker", all.x = TRUE, sort = FALSE)
  } else {
    df$chr <- NA_character_
    df$pos <- NA_real_
  }
  df
}

# Robust marker existence check (mqmsetcofactors() can be strict about marker naming)
marker_exists_in_cross <- function(cross, marker) {
  if (is.null(marker) || is.na(marker) || nchar(marker) == 0) return(FALSE)
  if (is.null(cross$geno) || length(cross$geno) == 0) return(FALSE)
  any(vapply(cross$geno, function(g) {
    if (is.null(g) || is.null(g$data)) return(FALSE)
    marker %in% colnames(g$data)
  }, logical(1)))
}

safe_mqmsetcofactors <- function(cross, markers) {
  markers <- unique(as.character(markers))
  markers <- markers[is.finite(nchar(markers)) & nchar(markers) > 0]
  markers <- markers[vapply(markers, function(m) marker_exists_in_cross(cross, m), logical(1))]
  if (length(markers) == 0) return(list(cross = cross, good = character(0), dropped = character(0)))

  good <- character(0)
  dropped <- character(0)

  # Try all at once first
  res_all <- try(mqmsetcofactors(cross, cofactors = markers), silent = TRUE)
  if (!inherits(res_all, "try-error")) {
    return(list(cross = res_all, good = markers, dropped = character(0)))
  }

  # Fallback: add cofactors one by one, keeping those that work
  cur <- cross
  for (m in markers) {
    res <- try(mqmsetcofactors(cur, cofactors = c(good, m)), silent = TRUE)
    if (inherits(res, "try-error")) {
      dropped <- c(dropped, m)
    } else {
      cur <- res
      good <- c(good, m)
    }
  }
  list(cross = cur, good = good, dropped = dropped)
}

select_top_cofactors <- function(df, nmax, min_dist, lod_thr) {
  df2 <- df[is.finite(df$lod) & !is.na(df$chr) & !is.na(df$pos), , drop=FALSE]
  df2 <- df2[df2$lod >= lod_thr, , drop=FALSE]
  df2 <- df2[order(df2$lod, decreasing = TRUE), , drop=FALSE]
  picked <- character(0)
  picked_rows <- list()
  for (i in seq_len(nrow(df2))) {
    m <- df2$marker[i]
    chr <- as.character(df2$chr[i])
    pos <- as.numeric(df2$pos[i])
    ok <- TRUE
    if (length(picked) > 0) {
      prev <- df2[df2$marker %in% picked, , drop=FALSE]
      prev_chr <- prev[as.character(prev$chr) == chr, , drop=FALSE]
      if (nrow(prev_chr) > 0) {
        if (any(abs(prev_chr$pos - pos) < min_dist)) ok <- FALSE
      }
    }
    if (!ok) next
    picked <- c(picked, m)
    picked_rows[[length(picked_rows)+1]] <- df2[i, , drop=FALSE]
    if (length(picked) >= nmax) break
  }
  if (length(picked_rows) == 0) {
    return(list(markers = character(0), df = df2[0, , drop=FALSE]))
  }
  outdf <- do.call(rbind, picked_rows)
  outdf$rank <- seq_len(nrow(outdf))
  list(markers = outdf$marker, df = outdf)
}

cap_n_marcovar <- function(n_marcovar, cross) {
  nind <- try(qtl::nind(cross), silent = TRUE)
  if (inherits(nind, "try-error") || is.null(nind) || is.na(nind)) return(n_marcovar)
  # conservative cap: keep some degrees of freedom
  cap <- max(0, as.integer(nind) - 10L)
  min(n_marcovar, cap)
}

cat("[QTL] loading cross: ", cross_rds, "\n")
cross <- readRDS(cross_rds)

if (is.null(trait)) {
  trait <- names(cross$pheno)[2]
}

if (inherits(cross, "cross")) {
  cl <- class(cross)
  if (length(cl) >= 2 && cl[1] == "cross") class(cross) <- c(cl[-1], "cross")
}

# Remove no variant markers
g <- pull.geno(cross)  # individuals x markers
mono_markers <- colnames(g)[apply(g, 2, function(x) {
  ux <- unique(x[!is.na(x)])
  length(ux) <= 1
})]
cross <- drop.markers(cross, mono_markers)

# Remove duplicated markers
dup <- findDupMarkers(cross, exact=TRUE)
to_drop <- unlist(lapply(dup, function(z) z[-1]))
cross <- drop.markers(cross, to_drop)

ph <- qtl::pull.pheno(cross)
if (!(trait %in% colnames(ph))) {
  stop(sprintf("trait '%s' not found in cross$pheno. Available: %s", trait, paste(colnames(ph), collapse=", ")))
}

cat("[QTL] calc.genoprob step=", step, " error_prob=", error_prob, " map=", map_function, "\n")
cross <- qtl::calc.genoprob(cross, 
                            step = step, 
                            off.end = 0,
                            error.prob = error_prob, 
                            map.function = map_function,
                            stepwidth = "fixed")

cofactors_selected <- character(0)
cofactors_df <- NULL

build_cofactors <- function(cross, trait, cofactor_mode, n_marcovar) {
  cof_selected <- character(0)
  cof_df <- data.frame()
  if (cofactor_mode == "file") {
    if (is.null(cofactor_file) || !file.exists(cofactor_file)) stop("cofactor_mode='file' but cofactor_file not found")
    cf <- try(data.table::fread(cofactor_file, header = TRUE, sep = "\t"), silent = TRUE)
    if (inherits(cf, "try-error")) {
      lines <- readLines(cofactor_file, warn = FALSE)
      lines <- lines[nchar(trimws(lines)) > 0]
      cof_selected <- unique(trimws(lines))
      cof_df <- data.frame(marker = cof_selected, stringsAsFactors = FALSE)
    } else {
      cn <- tolower(names(cf))
      if ("marker" %in% cn) {
        cof_selected <- unique(as.character(cf[[ which(cn=="marker")[1] ]]))
      } else {
        cof_selected <- unique(as.character(cf[[1]]))
      }
      cof_selected <- cof_selected[nchar(cof_selected) > 0]
      cof_df <- data.frame(marker = cof_selected, stringsAsFactors = FALSE)
    }
  } else if (cofactor_mode == "top") {
    cat("[QTL] selecting cofactors by marker-wise LOD: n=", n_marcovar,
        " lod_thr=", cofactor_lod_threshold, " min_dist=", cofactor_min_dist, "\n")
    df_lod <- marker_lod_scan(cross, trait)
    sel <- select_top_cofactors(df_lod, nmax = n_marcovar, min_dist = cofactor_min_dist, lod_thr = cofactor_lod_threshold)
    cof_selected <- sel$markers
    cof_df <- sel$df
  } else {
    cof_selected <- character(0)
    cof_df <- data.frame()
  }
  list(markers = cof_selected, df = cof_df)
}

run_cim_stable <- function(cross, trait, method, window, n_marcovar, do_perm = FALSE, n_perm = 0L) {
  cim_fm <- names(formals(qtl::cim))
  n_eff <- cap_n_marcovar(n_marcovar, cross)
  if (n_eff != n_marcovar) {
    cat("[QTL] stabilize: cap n_marcovar ", n_marcovar, " -> ", n_eff, "\n")
    n_marcovar <- n_eff
  }

  # base args
  args0 <- list(cross, pheno.col = trait, method = method, window = window)
  if ("n.marcovar" %in% cim_fm) args0$n.marcovar <- n_marcovar

  # if user provided cofactors, try to pass them
  cof_sel <- cofactors_selected
  if (length(cof_sel) > 0) {
    allm <- colnames(qtl::pull.geno(cross))
    cof_sel <- intersect(cof_sel, allm)
  }
  if (length(cof_sel) > 0) {
    if ("cofactors" %in% cim_fm) {
      args0$cofactors <- cof_sel
    } else if ("selected.markers" %in% cim_fm) {
      args0$selected.markers <- cof_sel
    } else if ("addcovar" %in% cim_fm) {
      args0$addcovar <- qtl::pull.geno(cross)[, cof_sel, drop = FALSE]
    } else if ("covar" %in% cim_fm) {
      args0$covar <- qtl::pull.geno(cross)[, cof_sel, drop = FALSE]
    } else {
      cat("[QTL] WARN: cim() does not expose cofactors args in this qtl version; using internal selection.\n")
    }
  }

  # iterative stabilization on over-spec warning
  it <- 0L
  warns_all <- character(0)
  res <- NULL
  while (TRUE) {
    it <- it + 1L
    a <- args0
    if (do_perm) a$n.perm <- as.integer(n_perm)

    out <- with_warn_capture(do.call(qtl::cim, a))
    warns_all <- unique(c(warns_all, out$warns))
    res <- out$val

    if (!(stabilize && has_overspec_warn(out$warns) && it < stabilize_max_iter)) break

    # reduce n.marcovar if possible
    if ("n.marcovar" %in% cim_fm && !is.null(a$n.marcovar) && a$n.marcovar > 0) {
      new_n <- max(0L, as.integer(floor(a$n.marcovar / 2)))
      cat("[QTL] stabilize: over-specified addcovar warning; retry with n_marcovar ", a$n.marcovar, " -> ", new_n, "\n")
      args0$n.marcovar <- new_n
    } else {
      break
    }
  }

  n_used <- if ("n.marcovar" %in% cim_fm && !is.null(args0$n.marcovar)) as.integer(args0$n.marcovar) else NA_integer_
  list(lod = res, warns = warns_all, n_marcovar_used = n_used)
}

lod <- NULL
perm <- NULL
thr <- NA_real_
extra_files <- list()

if (analysis_mode == "scanone") {
  cat("[QTL] scanone method=hk\n")
  lod <- qtl::scanone(cross, 
                      pheno.col = trait, 
                      method = p$sim_method,
                      model = p$sim_model,
                      addcovar = NULL, 
                      intcovar = NULL, 
                      weights = NULL,
                      #use = "all.obs",
                      upper = F,
                      ties.random = F, 
                      start = NULL,
                      maxit = 4000,
                      tol = 1e-4
                      #batchsize=250,
                      #n.cluster=1,
                      #ind.noqtl
                      )
  if (n_perm > 0) {
    cat("[QTL] permutation (scanone) n_perm=", n_perm, "\n")
    perm <- qtl::scanone(cross, 
                         pheno.col = trait, 
                         method = p$sim_method, 
                         n.perm = n_perm,
                         model = p$sim_model,
                         addcovar = NULL, 
                         intcovar = NULL, 
                         weights = NULL,
                         #use = "all.obs",
                         upper = F,
                         ties.random = F, 
                         start = NULL,
                         maxit = 4000,
                         tol = 1e-4,
                         perm.Xsp = F, 
                         perm.strata = NULL, 
                         verbose = T,
                         #batchsize=250,
                         n.cluster = 1
                         #ind.noqtl
                         )
    thr <- as.numeric(summary(perm, 
                              alpha = alpha))
    cat("[QTL] threshold(alpha=", alpha, ")=", thr, "\n")
  }

} else if (analysis_mode == "cim") {
  cat("[QTL] CIM n_marcovar=", n_marcovar, " window=", window, " cofactor_mode=", cofactor_mode, " stabilize=", stabilize, "\n")

  # Build cofactors (if requested)
  sel <- build_cofactors(cross, trait, cofactor_mode, n_marcovar)
  cofactors_selected <- sel$markers
  cofactors_df <- sel$df

  out1 <- run_cim_stable(cross, 
                         trait, 
                         method = p$cim_method, 
                         window = window, 
                         n_marcovar = n_marcovar, 
                         do_perm = FALSE)
  lod <- out1$lod
  n_marcovar_used <- out1$n_marcovar_used

  if (n_perm > 0) {
    cat("[QTL] permutation (CIM) n_perm=", n_perm, "\n")
    outp <- run_cim_stable(cross, 
                           trait, 
                           method = p$cim_method, 
                           window = window, 
                           n_marcovar = if (is.na(n_marcovar_used)) n_marcovar else n_marcovar_used,
                           do_perm = TRUE, 
                           n_perm = n_perm)
    perm <- outp$lod
    thr <- as.numeric(summary(perm, alpha = alpha))
    cat("[QTL] threshold(alpha=", alpha, ")=", thr, "\n")
  }

} else if (analysis_mode == "mqm") {
  cat("[QTL] MQM (mqmscan) requested\n")

  # # We will still plot a LOD profile from mqmscan.
  # # Cofactors are strongly recommended; we reuse the same selection logic.
  # sel <- build_cofactors(cross, trait, if (cofactor_mode == "auto") "top" else cofactor_mode, n_marcovar)
  # cofactors_selected <- sel$markers
  # cofactors_df <- sel$df
  # 
  # # Validate marker names (mqmsetcofactors() is strict and will error on any unknown marker)
  # all_markers <- qtl::markernames(cross)
  # if (length(cofactors_selected) > 0) {
  #   cofactors_selected <- unique(cofactors_selected[!is.na(cofactors_selected) & nchar(cofactors_selected) > 0])
  #   dropped <- setdiff(cofactors_selected, all_markers)
  #   if (length(dropped) > 0) {
  #     cat("[QTL] MQM: dropping non-existent cofactor markers: ", paste(dropped, collapse=","), "\n")
  #   }
  #   cofactors_selected <- intersect(cofactors_selected, all_markers)
  # }
  # 
  # if (length(cofactors_selected) > 0) {
  #   # mqmsetcofactors() will error if even one marker name is not recognized.
  #   # Therefore we attempt all-at-once, and if that fails we add cofactors one-by-one.
  #   res <- safe_mqmsetcofactors(cross, 
  #                               cofactors_selected)
  #   cross <- res$cross
  #   if (length(res$good) > 0) {
  #     cat("[QTL] MQM cofactors set: ", length(res$good), "\n")
  #   } else {
  #     cat("[QTL] WARN: MQM could not set any cofactors; running MQM without cofactors\n")
  #   }
  #   if (length(res$dropped) > 0) {
  #     cat("[QTL] MQM dropped cofactors: ", paste(res$dropped, collapse=","), "\n")
  #   }
  #   # Keep only those that actually got set
  #   cofactors_selected <- res$good
  #   if (nrow(cofactors_df) > 0 && length(cofactors_selected) > 0) {
  #     cofactors_df <- cofactors_df[cofactors_df$marker %in% cofactors_selected, , drop=FALSE]
  #   }
  # } else {
  #   cat("[QTL] MQM running without explicit cofactors\n")
  # }
  # 
  cof <- mqmautocofactors(cross = cross,
                          num = as.numeric(p$mqm_cof),
                          distance = as.numeric(p$mqm_dis),
                          dominance = as.logical(p$mqm_dom),
                          plot = F,
                          verbose = T)
  cofactors_df <- data.frame(marker = cof)
  
  if (as.logical(p$mqm_dom) == T) {
    dom <- "dominance"
  } else {
    dom <- "additive"
  }
  lod <- mqmscan(cross = cross, 
                 cofactors = cof,
                 pheno.col = trait,
                 model = , 
                 forceML = F,
                 cofactor.significance = 0.02,
                 em.iter = 1000,
                 window.size = p$mqm_win, 
                 step.size = p$mqm_step,
                 logtransform = F, 
                 estimate.map = F,
                 plot = F, 
                 verbose = T, 
                 outputmarkers = T,
                 multicore = T, 
                 batchsize = 10, 
                 n.clusters = 1, 
                 test.normality = F,
                 off.end = 0)
  colnames(lod)[1:3] <- c("chr", "pos", "lod")
  
  if (n_perm > 0) {
    cat("[QTL] MQM permutation n_perm=", n_perm, "\n")
    perm <- mqmpermutation(cross, 
                           scanfunction = mqmscan, 
                           pheno.col = trait,
                           multicore = T,
                           n.perm = n_perm,
                           file = "MQM_output.txt",
                           n.cluster = 1, 
                           method = "permutation",
                           cofactors = cof, 
                           plot = F, 
                           verbose = T)
    perm <- mqmprocesspermutation(mqmpermutationresult = perm)
    thr <- as.numeric(summary(perm, 
                              alpha = alpha))
    cat("[QTL] threshold(alpha=", alpha, ")=", thr, "\n")
  } else if (n_perm > 0) {
    cat("[QTL] WARN: skipping permutation.\n")
  }

} else if (analysis_mode == "stepwise") {
  cat("[QTL] Stepwise (additive) via stepwiseqtl\n")

  # Use scanone LOD curve as the profile for plotting and for penalty (optional)
  lod0 <- qtl::scanone(cross, pheno.col = trait, method = "hk")
  lod <- lod0

  thr_main <- NA_real_
  if (n_perm > 0) {
    cat("[QTL] permutation (scanone hk) for stepwise penalty n_perm=", n_perm, "\n")
    perm0 <- qtl::scanone(cross, pheno.col = trait, method = "hk", n.perm = n_perm)
    perm <- perm0
    thr_main <- as.numeric(summary(perm0, alpha = alpha))
    cat("[QTL] scanone perm threshold(alpha=", alpha, ")=", thr_main, "\n")
  }

  penalty_main <- thr_main
  if (stepwise_penalty_mode == "manual") {
    if (!is.finite(stepwise_penalty_main)) stop("stepwise_penalty_mode='manual' but stepwise_penalty_main is not numeric")
    penalty_main <- as.numeric(stepwise_penalty_main)
  }
  if (!is.finite(penalty_main)) {
    # safe default if no permutation
    penalty_main <- 3.0
    cat("[QTL] WARN: penalty_main not available; using default 3.0\n")
  }
  cat("[QTL] stepwise penalties(main)=", penalty_main, " max_qtl=", stepwise_max_qtl, "\n")

  # stepwiseqtl args vary across versions; keep minimal and additive-only
  sw_fm <- names(formals(qtl::stepwiseqtl))
  sw_args <- list(cross, pheno.col = trait, method = "hk", max.qtl = stepwise_max_qtl, penalties = c(penalty_main))
  if ("additive.only" %in% sw_fm) sw_args$additive.only <- TRUE
  if ("keeptrace" %in% sw_fm) sw_args$keeptrace <- TRUE

  sw <- do.call(qtl::stepwiseqtl, sw_args)
  saveRDS(sw, file.path(out_dir, "qtl_model.rds"))

  # --- summarize model robustly (avoid as.data.frame(summary(sw)) which is version-dependent)
  # stepwiseqtl() return type can vary across qtl versions; it can be a list, a qtl object, or even NULL.
  qtl_obj <- NULL
  if (is.list(sw) && ("qtl" %in% names(sw)) && !is.null(sw[["qtl"]])) {
    qtl_obj <- sw[["qtl"]]
  } else if (is.list(sw) && ("best.qtl" %in% names(sw)) && !is.null(sw[["best.qtl"]])) {
    qtl_obj <- sw[["best.qtl"]]
  } else if (inherits(sw, "qtl")) {
    qtl_obj <- sw
  } else {
    cat("[QTL] WARN: stepwiseqtl returned unexpected type; treating as no QTL selected\n")
  }

  qtl_model_df <- data.frame()
  if (!is.null(qtl_obj) && !is.null(qtl_obj$chr) && length(qtl_obj$chr) > 0) {
    qtl_model_df <- data.frame(
      qtl = paste0("Q", seq_along(qtl_obj$chr)),
      chr = as.character(qtl_obj$chr),
      pos = as.numeric(qtl_obj$pos),
      stringsAsFactors = FALSE
    )
  }
  qtl_model_tsv <- file.path(out_dir, "qtl_model.tsv")
  data.table::fwrite(qtl_model_df, qtl_model_tsv, sep = "\t")
  extra_files$qtl_model_tsv <- qtl_model_tsv

  # fitqtl for reporting
  fit_tsv <- file.path(out_dir, "fitqtl.tsv")
  if (is.null(qtl_obj) || is.null(qtl_obj$chr) || length(qtl_obj$chr) == 0) {
    data.table::fwrite(data.frame(note = "no qtl selected"), fit_tsv, sep = "\t")
    extra_files$fitqtl_tsv <- fit_tsv
  } else {
    fit <- try(qtl::fitqtl(cross, pheno.col = trait, qtl = qtl_obj, method = "hk", get.ests = TRUE, dropone = TRUE), silent = TRUE)
    if (!inherits(fit, "try-error")) {
      if (!is.null(fit$ests)) {
        fit_df <- as.data.frame(fit$ests)
      } else {
        fit_df <- data.frame(lod = fit$lod, n.qtl = fit$n.qtl, df = fit$df, stringsAsFactors = FALSE)
      }
      data.table::fwrite(fit_df, fit_tsv, sep = "\t")
    } else {
      data.table::fwrite(data.frame(error = as.character(fit)), fit_tsv, sep = "\t")
    }
    extra_files$fitqtl_tsv <- fit_tsv
  }

  # peaks.tsv should show selected QTL
  peaks_df <- qtl_model_df
  if (nrow(peaks_df) == 0) {
    # fallback: report max of scanone
    lod_df0 <- as.data.frame(lod0)
    i <- which.max(lod_df0$lod)
    peaks_df <- lod_df0[i, , drop = FALSE]
    peaks_df$p.value <- NA_real_
  }
  peaks_tsv <- file.path(out_dir, "peaks.tsv")
  data.table::fwrite(peaks_df, peaks_tsv, sep = "\t")
  # store threshold as "thr" for plot line
  thr <- thr_main

  # also store a trace if present
  if (!is.null(attr(sw, "trace"))) {
    tr <- attr(sw, "trace")
    trace_tsv <- file.path(out_dir, "stepwise_trace.tsv")
    # trace can be a list with variable-length rows depending on qtl count; coerce safely
    tr_df <- try(as.data.frame(tr), silent = TRUE)
    if (inherits(tr_df, "try-error")) {
      # fallback: store as text
      writeLines(capture.output(print(tr)), trace_tsv)
    } else {
      data.table::fwrite(tr_df, trace_tsv, sep = "\t")
    }
    extra_files$stepwise_trace_tsv <- trace_tsv
  }
}

## Save cofactors (CIM/MQM)
cofactors_tsv <- file.path(out_dir, "cofactors.tsv")
if (analysis_mode %in% c("cim", "mqm")) {
  if (is.null(cofactors_df)) cofactors_df <- data.frame()
  if (nrow(cofactors_df) == 0) cofactors_df <- data.frame(marker = character(0), stringsAsFactors = FALSE)
  cofactors_df$cofactor_mode <- rep(cofactor_mode, nrow(cofactors_df))
  cofactors_df$n_marcovar <- rep(n_marcovar, nrow(cofactors_df))
  cofactors_df$window <- rep(window, nrow(cofactors_df))
  data.table::fwrite(cofactors_df, cofactors_tsv, sep = "\t")
}

## Save LOD profile
lod_df <- as.data.frame(lod)
lod_tsv <- file.path(out_dir, "lod_profile.tsv")
data.table::fwrite(lod_df, lod_tsv, sep = "\t")

## Save permutation threshold
thr_tsv <- file.path(out_dir, "perm_thresholds.tsv")
data.table::fwrite(data.frame(alpha = alpha, threshold = thr, n_perm = n_perm), thr_tsv, sep = "\t")

## Peaks table for non-stepwise modes
if (analysis_mode %in% c("scanone", "cim", "mqm")) {
  peaks_tsv <- file.path(out_dir, "peaks.tsv")
  peaks_df <- data.frame()
  if (!is.null(perm)) {
    s <- try(summary(lod, perms = perm, alpha = alpha, pvalues = TRUE), silent = TRUE)
    if (!inherits(s, "try-error") && !is.null(s) && nrow(as.data.frame(s)) > 0) {
      peaks_df <- as.data.frame(s)
    }
  } else {
    pk <- try(qtl::find.pks(lod, cutoff = max(lod_df$lod, na.rm=TRUE) - 1e-9), silent = TRUE)
    if (!inherits(pk, "try-error") && !is.null(pk) && nrow(pk) > 0) {
      peaks_df <- pk
    }
  }
  if (nrow(peaks_df) == 0) {
    i <- which.max(lod_df$lod)
    peaks_df <- lod_df[i, , drop = FALSE]
    peaks_df$p.value <- NA_real_
  }
  data.table::fwrite(peaks_df, peaks_tsv, sep = "\t", quote = F)
}

## Plot
plot_png <- file.path(out_dir, "lod_plot.png")
png(plot_png, width = 1200, height = 700)
main_title <- switch(
  analysis_mode,
  scanone = sprintf("scanone (hk): %s", trait),
  cim = sprintf("CIM (hk): %s", trait),
  mqm = sprintf("MQM: %s", trait),
  stepwise = sprintf("Stepwise (additive): %s", trait),
  sprintf("QTL: %s", trait)
)
plot(lod, main = main_title)
if (!is.na(thr)) {
  abline(h = thr, lty = 2)
  legend("topright", legend = sprintf("thr (alpha=%.3f)=%.3f", alpha, thr), bty = "n")
}
dev.off()

artifacts <- list(
  lod_profile_tsv = lod_tsv,
  cofactors_tsv = if (analysis_mode %in% c("cim", "mqm")) cofactors_tsv else NULL,
  peaks_tsv = file.path(out_dir, "peaks.tsv"),
  perm_thresholds_tsv = thr_tsv,
  lod_plot_png = plot_png,
  method = paste0(analysis_mode, "_hk"),
  n_perm = n_perm,
  alpha = alpha,
  step = step,
  error_prob = error_prob,
  map_function = map_function,
  trait = trait,
  analysis_mode = analysis_mode,
  n_marcovar = n_marcovar,
  window = window,
  cofactor_mode = if (analysis_mode %in% c("cim", "mqm")) cofactor_mode else NULL,
  cofactor_file = if (analysis_mode %in% c("cim", "mqm")) cofactor_file else NULL,
  cofactor_lod_threshold = if (analysis_mode %in% c("cim", "mqm")) cofactor_lod_threshold else NULL,
  cofactor_min_dist = if (analysis_mode %in% c("cim", "mqm")) cofactor_min_dist else NULL,
  stabilize = stabilize,
  stabilize_max_iter = stabilize_max_iter,
  stepwise_max_qtl = if (analysis_mode == "stepwise") stepwise_max_qtl else NULL,
  stepwise_penalty_mode = if (analysis_mode == "stepwise") stepwise_penalty_mode else NULL,
  stepwise_penalty_main = if (analysis_mode == "stepwise") stepwise_penalty_main else NULL
)

for (nm in names(extra_files)) artifacts[[nm]] <- extra_files[[nm]]

write(toJSON(artifacts, auto_unbox = TRUE, pretty = TRUE), file = file.path(out_dir, "artifacts.json"))
cat("[QTL] done\n")
