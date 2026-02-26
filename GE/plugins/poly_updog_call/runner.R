#!/usr/bin/env Rscript

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
log <- function(...){
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
safe_png <- function(path, fun){
  tryCatch({
    png(path, width=1200, height=700)
    on.exit(dev.off(), add=TRUE)
    fun()
  }, error=function(e){
    try(dev.off(), silent=TRUE)
  })
}

params <- tryCatch(jsonlite::fromJSON(params_path), error=function(e) list())
log("[poly_updog_call] start")
log("[poly_updog_call] params_path=", params_path)
log("[poly_updog_call] out_dir=", out_dir)

ploidy <- if(!is.null(params$ploidy)) as.integer(params$ploidy) else 4L
counts_rds <- if(!is.null(params$counts_rds)) as.character(params$counts_rds) else ""
method <- if(!is.null(params$method)) as.character(params$method) else "multidog"
ncores <- if(!is.null(params$ncores)) as.integer(params$ncores) else 1L
subset_markers <- if(!is.null(params$subset_markers)) as.character(params$subset_markers) else NULL

# --- NEW: guaranteed posterior probabilities ---
output_prob <- if(!is.null(params$output_prob)) as.logical(params$output_prob) else TRUE
call_prob_threshold <- if(!is.null(params$call_prob_threshold)) as.numeric(params$call_prob_threshold) else 0.9
prob_format <- if(!is.null(params$prob_format)) as.character(params$prob_format) else "long"  # long|wide|none

# parent/progeny options (optional; for f1/s1 models)
parent1_id <- if(!is.null(params$parent1_id)) as.character(params$parent1_id) else ""
parent2_id <- if(!is.null(params$parent2_id)) as.character(params$parent2_id) else ""
progeny_ids_raw <- if(!is.null(params$progeny_ids)) as.character(params$progeny_ids) else ""
progeny_file <- if(!is.null(params$progeny_file)) as.character(params$progeny_file) else ""

parse_id_list <- function(x){
  x <- as.character(x)
  x <- gsub("[\r\n\t]", " ", x)
  parts <- unlist(strsplit(x, "[,; ]+"))
  parts <- trimws(parts)
  parts[nzchar(parts)]
}
extra_options <- params$extra_options
if (is.null(extra_options)) extra_options <- list()

if (!nzchar(counts_rds) || !file.exists(counts_rds)) {
  msg <- paste0("counts_rds not found: ", counts_rds)
  write_error(msg)
  write_artifacts(list(module="poly_updog_call", error=msg))
  quit(status=1)
}

have_updog <- requireNamespace("updog", quietly=TRUE)
if (!have_updog) {
  msg <- "R package 'updog' is not installed. Install from CRAN/Bioconductor as appropriate."
  write_error(msg)
  write_artifacts(list(module="poly_updog_call", error=msg))
  quit(status=1)
}

pc <- readRDS(counts_rds)
ref <- as.matrix(pc$ref)
total <- as.matrix(pc$total)
marker_info <- pc$marker
sample_info <- pc$sample

if (!is.null(subset_markers) && length(subset_markers) > 0) {
  keep <- intersect(rownames(total), subset_markers)
  ref <- ref[keep, , drop=FALSE]
  total <- total[keep, , drop=FALSE]
  if (!is.null(marker_info) && "marker" %in% colnames(marker_info)) {
    marker_info <- marker_info[marker_info$marker %in% keep, , drop=FALSE]
  }
  log("[poly_updog_call] subset markers: ", nrow(total))
}

# Ensure numeric and consistent
storage.mode(ref) <- "numeric"
storage.mode(total) <- "numeric"
total[total <= 0] <- NA_real_
ref[is.na(total)] <- NA_real_

markers <- rownames(total)
samples <- colnames(total)

# Resolve parents/progeny
parents <- character(0)
if (nzchar(parent1_id)) parents <- c(parents, parent1_id)
if (nzchar(parent2_id)) parents <- c(parents, parent2_id)
parents <- unique(parents)

missing_parents <- setdiff(parents, samples)
if (length(missing_parents) > 0) {
  log("[poly_updog_call] WARNING: parent IDs not found in counts_rds columns: ", paste(missing_parents, collapse=","))
}

progeny_samples <- NULL
if (nzchar(progeny_file) && file.exists(progeny_file)) {
  ids <- readLines(progeny_file, warn=FALSE)
  ids <- trimws(ids)
  ids <- ids[nzchar(ids)]
  progeny_samples <- intersect(ids, samples)
  log("[poly_updog_call] progeny from file: n=", length(progeny_samples))
} else if (nzchar(progeny_ids_raw)) {
  ids <- parse_id_list(progeny_ids_raw)
  progeny_samples <- intersect(ids, samples)
  log("[poly_updog_call] progeny from list: n=", length(progeny_samples))
} else if (length(parents) > 0) {
  progeny_samples <- setdiff(samples, parents)
  log("[poly_updog_call] progeny=all non-parents: n=", length(progeny_samples))
} else {
  progeny_samples <- samples
  log("[poly_updog_call] progeny=all samples: n=", length(progeny_samples))
}

if (length(progeny_samples) == 0) {
  log("[poly_updog_call] WARNING: no progeny samples selected; falling back to all samples")
  progeny_samples <- samples
}

# Parents present in counts (kept for downstream linkage mapping).
parents_present <- parents[parents %in% samples]
sample_order_all <- unique(c(parents_present, progeny_samples))


# Parent count vectors (per marker)
p1ref <- p1size <- p2ref <- p2size <- NULL
if (nzchar(parent1_id) && parent1_id %in% samples) {
  p1ref <- ref[, parent1_id]
  p1size <- total[, parent1_id]
}
if (nzchar(parent2_id) && parent2_id %in% samples) {
  p2ref <- ref[, parent2_id]
  p2size <- total[, parent2_id]
}

# Subset to progeny matrices for calling (parents are provided separately via p1*/p2*)
ref_call <- ref[, progeny_samples, drop=FALSE]
total_call <- total[, progeny_samples, drop=FALSE]

# Keep sample_info in sync (if available)
if (!is.null(sample_info)) {
  if ("sample" %in% colnames(sample_info)) {
    sample_info <- sample_info[sample_info$sample %in% sample_order_all, , drop=FALSE]
  } else if ("ind" %in% colnames(sample_info)) {
    sample_info <- sample_info[sample_info$ind %in% sample_order_all, , drop=FALSE]
  }
}

# helper: extract genotype from multidog result
extract_geno <- function(res){
  if (is.null(res)) return(NULL)
  for (nm in c("geno", "genotype", "call", "calls", "z", "dosage")) {
    if (!is.null(res[[nm]])) return(res[[nm]])
  }
  if (!is.null(res$fullout) && !is.null(res$fullout$geno)) return(res$fullout$geno)
  NULL
}

# helper: extract posterior probabilities from flexdog result
extract_postprob <- function(fd, ploidy){
  cand <- NULL
  for (nm in c("postprob", "postprobs", "posterior", "posteriors", "prob", "probs")) {
    if (!is.null(fd[[nm]])) { cand <- fd[[nm]]; break }
  }
  if (is.null(cand)) return(NULL)
  m <- as.matrix(cand)
  if (ncol(m) != (ploidy + 1)) return(NULL)
  if (is.null(colnames(m))) colnames(m) <- as.character(0:ploidy)
  suppressWarnings({
    if (all(is.finite(as.numeric(colnames(m))))) {
      ord <- order(as.numeric(colnames(m)))
      m <- m[, ord, drop=FALSE]
    }
  })
  m
}

# enforce probability output: use flexdog for guaranteed posterior probs
if (isTRUE(output_prob) && method != "flexdog") {
  log("[poly_updog_call] output_prob=TRUE -> force method=flexdog (guaranteed posterior probs)")
  method <- "flexdog"
}

markers_n <- length(markers)
samples_n <- length(progeny_samples)
dosage_levels <- as.character(0:ploidy)

samples_all_n <- length(sample_order_all)

# main containers
# NOTE: we keep dosage_mat as integer matrix marker x sample
# For output_prob=TRUE, dosage_mat is derived from prob_arr (MAP)
dosage_mat <- matrix(NA_integer_, nrow=markers_n, ncol=samples_n,
                     dimnames=list(markers, progeny_samples))
fit_summary <- data.frame(marker=markers, converged=NA, stringsAsFactors=FALSE)

# probability array (marker x sample x (P+1))
prob_arr <- array(NA_real_,
                  dim=c(markers_n, samples_n, ploidy+1),
                  dimnames=list(markers, progeny_samples, dosage_levels))

# Full containers (parents + progeny). Parents are appended before progeny.
prob_arr_all <- array(NA_real_,
                      dim=c(markers_n, samples_all_n, ploidy+1),
                      dimnames=list(markers, sample_order_all, dosage_levels))
dosage_mat_all <- matrix(NA_integer_, nrow=markers_n, ncol=samples_all_n,
                         dimnames=list(markers, sample_order_all))

failed <- data.frame(marker_id=character(), reason=character(), stringsAsFactors=FALSE)

tryCatch({
  if (method == "multidog" && "multidog" %in% getNamespaceExports("updog")) {
    log("[poly_updog_call] using updog::multidog (dosage only)")
    args_md <- c(list(refmat=ref_call, sizemat=total_call, ploidy=ploidy), extra_options)
    if (!is.null(p1ref) && !is.null(p1size)) { args_md$p1ref <- p1ref; args_md$p1size <- p1size }
    if (!is.null(p2ref) && !is.null(p2size)) { args_md$p2ref <- p2ref; args_md$p2size <- p2size }
    res <- do.call(updog::multidog, args_md)
    g <- extract_geno(res)
    if (is.null(g)) {
      if (is.matrix(res) || is.data.frame(res)) g <- res
    }
    if (is.null(g)) stop("Could not extract genotype from multidog")
    g <- as.matrix(g)
    if (nrow(g) == samples_n && ncol(g) == markers_n) g <- t(g)
    if (!all(dim(g) == dim(dosage_mat))) stop("genotype matrix dimension mismatch")
    dosage_mat[,] <- as.integer(round(g))
    dosage_mat_all[, progeny_samples] <- dosage_mat
    fit_summary$converged <- TRUE
    # prob_arr remains NA (and will be empty output if output_prob=FALSE)
  }

  if (method == "flexdog") {
    log("[poly_updog_call] using per-marker updog::flexdog")
    if (!("flexdog" %in% getNamespaceExports("updog"))) stop("updog::flexdog not available in this updog version")

    idx <- seq_along(markers)
    run_one <- function(i){
      m <- markers[i]
      r <- ref_call[m, ]
      s <- total_call[m, ]
      ok <- which(!is.na(r) & !is.na(s))
      if (length(ok) < 4) {
        return(list(marker=m, converged=FALSE,
                    prob=matrix(NA_real_, nrow=samples_n, ncol=ploidy+1,
                                dimnames=list(progeny_samples, dosage_levels)),
                    p1prob=rep(NA_real_, ploidy+1),
                    p2prob=rep(NA_real_, ploidy+1),
                    err="too_few_nonmissing"))
      }
      out <- tryCatch({
        args_fd <- c(list(refvec=r, sizevec=s, ploidy=ploidy), extra_options)
        if (!is.null(p1ref) && !is.null(p1size)) { args_fd$p1ref <- p1ref[i]; args_fd$p1size <- p1size[i] }
        if (!is.null(p2ref) && !is.null(p2size)) { args_fd$p2ref <- p2ref[i]; args_fd$p2size <- p2size[i] }
        fd <- do.call(updog::flexdog, args_fd)

        # Extract parameter estimates to enable parent posterior calls with fixed params.
        get1 <- function(obj, nm){ if (!is.null(obj[[nm]])) return(obj[[nm]]); if (!is.null(obj$par) && !is.null(obj$par[[nm]])) return(obj$par[[nm]]); NULL }
        bias_hat <- get1(fd, "bias"); if (is.null(bias_hat) && !is.null(extra_options$bias)) bias_hat <- extra_options$bias
        seq_hat  <- get1(fd, "seq");  if (is.null(seq_hat)  && !is.null(extra_options$seq))  seq_hat  <- extra_options$seq
        od_hat   <- get1(fd, "od");   if (is.null(od_hat)   && !is.null(extra_options$od))   od_hat   <- extra_options$od

        parent_postprob <- function(rp, sp){
          if (is.null(rp) || is.null(sp) || is.na(rp) || is.na(sp) || !is.finite(rp) || !is.finite(sp) || sp <= 0) {
            return(rep(NA_real_, ploidy+1))
          }
          args_p <- c(list(refvec=c(rp), sizevec=c(sp), ploidy=ploidy), extra_options)
          # override with fixed params and disable parameter updates (single sample)
          args_p$bias <- as.numeric(bias_hat)
          args_p$seq  <- as.numeric(seq_hat)
          args_p$od   <- as.numeric(od_hat)
          args_p$update_bias <- FALSE
          args_p$update_seq  <- FALSE
          args_p$update_od   <- FALSE
          # parents should not be used as parents when we are calling the parents themselves
          args_p$p1ref <- NULL; args_p$p1size <- NULL; args_p$p2ref <- NULL; args_p$p2size <- NULL
          fd2 <- do.call(updog::flexdog, args_p)
          pp2 <- extract_postprob(fd2, ploidy)
          if (is.null(pp2)) {
            g2 <- fd2$geno
            if (is.null(g2) && !is.null(fd2$genotype)) g2 <- fd2$genotype
            g2 <- as.integer(round(g2))
            outp <- rep(0, ploidy+1)
            if (!is.na(g2)) outp[g2+1] <- 1
            return(as.numeric(outp))
          }
          if (nrow(pp2) != 1) pp2 <- pp2[1, , drop=FALSE]
          as.numeric(pp2[1, ])
        }

        pp <- extract_postprob(fd, ploidy)
        if (is.null(pp)) {
          # fallback to degenerate probability at MAP genotype
          g <- fd$geno
          if (is.null(g) && !is.null(fd$genotype)) g <- fd$genotype
          g <- as.integer(round(g))
          pp <- matrix(0, nrow=samples_n, ncol=ploidy+1,
                       dimnames=list(progeny_samples, dosage_levels))
          for (j in seq_along(g)) if (!is.na(g[j])) pp[j, as.character(g[j])] <- 1
        } else {
          # ensure row order matches progeny_samples if rownames exist
          if (!is.null(rownames(pp))) {
            common <- intersect(progeny_samples, rownames(pp))
            pp2 <- matrix(NA_real_, nrow=samples_n, ncol=ploidy+1,
                          dimnames=list(progeny_samples, dosage_levels))
            pp2[common, ] <- pp[common, , drop=FALSE]
            pp <- pp2
          } else {
            # assume order matches input
            rownames(pp) <- progeny_samples
          }
        }
        p1pp <- if (!is.null(p1ref) && !is.null(p1size) && nzchar(parent1_id) && parent1_id %in% samples) parent_postprob(p1ref[i], p1size[i]) else rep(NA_real_, ploidy+1)
        p2pp <- if (!is.null(p2ref) && !is.null(p2size) && nzchar(parent2_id) && parent2_id %in% samples) parent_postprob(p2ref[i], p2size[i]) else rep(NA_real_, ploidy+1)
        list(marker=m, converged=TRUE, prob=pp, p1prob=p1pp, p2prob=p2pp, err=NA_character_)
      }, error=function(e){
        list(marker=m, converged=FALSE,
             prob=matrix(NA_real_, nrow=samples_n, ncol=ploidy+1,
                         dimnames=list(progeny_samples, dosage_levels)),
             p1prob=rep(NA_real_, ploidy+1),
             p2prob=rep(NA_real_, ploidy+1),
             err=conditionMessage(e))
      })
      out
    }

    res_list <- if (ncores > 1 && .Platform$OS.type != "windows") {
      parallel::mclapply(idx, run_one, mc.cores=ncores)
    } else {
      lapply(idx, run_one)
    }

    for (k in seq_along(res_list)) {
      rr <- res_list[[k]]
      prob_arr[rr$marker, , ] <- rr$prob
      # Fill full arrays (parents + progeny)
      prob_arr_all[rr$marker, progeny_samples, ] <- rr$prob
      if (nzchar(parent1_id) && parent1_id %in% sample_order_all) {
        prob_arr_all[rr$marker, parent1_id, ] <- rr$p1prob
      }
      if (nzchar(parent2_id) && parent2_id %in% sample_order_all) {
        prob_arr_all[rr$marker, parent2_id, ] <- rr$p2prob
      }
      fit_summary$converged[fit_summary$marker == rr$marker] <- isTRUE(rr$converged)
      if (!isTRUE(rr$converged)) {
        failed <- rbind(failed, data.frame(marker_id=rr$marker, reason=rr$err, stringsAsFactors=FALSE))
      }
    }

    # derive MAP dosage from prob
    dosage_mat <- apply(prob_arr, c(1,2), function(p){
      if (all(is.na(p))) return(NA_integer_)
      which.max(p) - 1L
    })
    storage.mode(dosage_mat) <- "integer"
    dimnames(dosage_mat) <- list(markers, progeny_samples)

    # fill dosage for all samples: progeny MAP, parents from their posterior
    dosage_mat_all[, progeny_samples] <- dosage_mat
    if (nzchar(parent1_id) && parent1_id %in% sample_order_all) {
      d1 <- apply(prob_arr_all[, parent1_id, , drop=FALSE], 1, function(p){ if (all(is.na(p))) NA_integer_ else which.max(p) - 1L })
      dosage_mat_all[, parent1_id] <- as.integer(d1)
    }
    if (nzchar(parent2_id) && parent2_id %in% sample_order_all) {
      d2 <- apply(prob_arr_all[, parent2_id, , drop=FALSE], 1, function(p){ if (all(is.na(p))) NA_integer_ else which.max(p) - 1L })
      dosage_mat_all[, parent2_id] <- as.integer(d2)
    }
  }

  # derived matrices (posterior-based summaries)
  dvec <- 0:ploidy

  # progeny-only (legacy)
  dosage_pp <- apply(prob_arr, c(1,2), function(p){
    if (all(is.na(p))) return(NA_real_)
    sum(p * dvec, na.rm=TRUE)
  })
  dosage_maxprob <- apply(prob_arr, c(1,2), function(p){
    if (all(is.na(p))) return(NA_real_)
    max(p, na.rm=TRUE)
  })
  dosage_entropy <- apply(prob_arr, c(1,2), function(p){
    if (all(is.na(p))) return(NA_real_)
    p <- p[is.finite(p) & p > 0]
    -sum(p * log(p))
  })

  # all samples (parents + progeny)
  dosage_pp_all <- apply(prob_arr_all, c(1,2), function(p){
    if (all(is.na(p))) return(NA_real_)
    sum(p * dvec, na.rm=TRUE)
  })
  dosage_maxprob_all <- apply(prob_arr_all, c(1,2), function(p){
    if (all(is.na(p))) return(NA_real_)
    max(p, na.rm=TRUE)
  })
  dosage_entropy_all <- apply(prob_arr_all, c(1,2), function(p){
    if (all(is.na(p))) return(NA_real_)
    p <- p[is.finite(p) & p > 0]
    -sum(p * log(p))
  })

  # write dosage outputs
  utils::write.table(
    data.frame(ind=progeny_samples, t(dosage_mat), check.names=FALSE),
    file.path(out_dir, "dosage.tsv"),
    sep="\t", quote=FALSE, row.names=FALSE
  )
  utils::write.table(
    data.frame(ind=sample_order_all, t(dosage_mat_all), check.names=FALSE),
    file.path(out_dir, "dosage_all.tsv"),
    sep="\t", quote=FALSE, row.names=FALSE
  )
  utils::write.table(
    data.frame(sample_id=progeny_samples),
    file.path(out_dir, "progeny_ids.tsv"),
    sep="\t", quote=FALSE, row.names=FALSE
  )

  # QC tables (based on all samples)
  called_all <- dosage_maxprob_all >= call_prob_threshold
  marker_stats <- data.frame(
    marker_id = markers,
    call_rate = rowMeans(called_all, na.rm=TRUE),
    mean_maxprob = rowMeans(dosage_maxprob_all, na.rm=TRUE),
    mean_entropy = rowMeans(dosage_entropy_all, na.rm=TRUE),
    n_nonmissing = rowSums(!is.na(dosage_mat_all)),
    stringsAsFactors=FALSE
  )
  sample_stats <- data.frame(
    sample_id = sample_order_all,
    call_rate = colMeans(called_all, na.rm=TRUE),
    mean_maxprob = colMeans(dosage_maxprob_all, na.rm=TRUE),
    mean_entropy = colMeans(dosage_entropy_all, na.rm=TRUE),
    n_nonmissing = colSums(!is.na(dosage_mat_all)),
    stringsAsFactors=FALSE
  )

  utils::write.table(marker_stats, file.path(out_dir, "marker_stats.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
  utils::write.table(sample_stats, file.path(out_dir, "sample_stats.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

  # legacy marker_callrate.tsv (keep for backward compatibility)
  call_rate_legacy <- 1 - rowMeans(is.na(dosage_mat))
  summ <- data.frame(marker=markers, call_rate=call_rate_legacy, stringsAsFactors=FALSE)
  if (!is.null(marker_info) && "marker" %in% colnames(marker_info)) {
    summ <- merge(marker_info, summ, by="marker", all.y=TRUE, sort=FALSE)
  }
  utils::write.table(summ, file.path(out_dir, "marker_callrate.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

  # plots
  safe_png(file.path(out_dir, "dosage_hist.png"), function(){
    x <- as.numeric(dosage_mat)
    x <- x[is.finite(x)]
    hist(x, breaks=seq(-0.5, ploidy+0.5, by=1), main="Dosage histogram", xlab="dosage")
  })
  safe_png(file.path(out_dir, "maxprob_hist.png"), function(){
    x <- as.numeric(dosage_maxprob)
    x <- x[is.finite(x)]
    hist(x, breaks=30, main=paste0("Max posterior prob (threshold=", call_prob_threshold, ")"), xlab="maxprob")
    abline(v=call_prob_threshold, lty=2)
  })

  # prob TSV output (optional)
  if (isTRUE(output_prob) && prob_format == "long") {
    gz <- gzfile(file.path(out_dir, "prob_long.tsv.gz"), "wt")
    on.exit(close(gz), add=TRUE)
    writeLines("marker_id\tsample_id\td\tp", gz)
    for (mi in seq_along(markers)) {
      m <- markers[mi]
      pp <- prob_arr[mi,,]
      for (si in seq_along(progeny_samples)) {
        s <- progeny_samples[si]
        pv <- pp[si,]
        if (all(is.na(pv))) next
        for (di in seq_along(dosage_levels)) {
          writeLines(paste(m, s, dosage_levels[di], format(pv[di], scientific=FALSE), sep="\t"), gz)
        }
      }
    }

    # also write parents+progeny (useful for linkage mapping)
    if (length(parents_present) > 0) {
      gz2 <- gzfile(file.path(out_dir, "prob_long_all.tsv.gz"), "wt")
      on.exit(close(gz2), add=TRUE)
      writeLines("marker_id\tsample_id\td\tp", gz2)
      for (mi in seq_along(markers)) {
        m <- markers[mi]
        pp <- prob_arr_all[mi,,]
        for (si in seq_along(sample_order_all)) {
          s <- sample_order_all[si]
          pv <- pp[si,]
          if (all(is.na(pv))) next
          for (di in seq_along(dosage_levels)) {
            writeLines(paste(m, s, dosage_levels[di], format(pv[di], scientific=FALSE), sep="\t"), gz2)
          }
        }
      }
    }
  } else if (isTRUE(output_prob) && prob_format == "wide") {
    # wide format (marker_id, sample_id, p0..pP)
    gz <- gzfile(file.path(out_dir, "prob_wide.tsv.gz"), "wt")
    on.exit(close(gz), add=TRUE)
    hdr <- paste(c("marker_id","sample_id", paste0("p", dosage_levels)), collapse="\t")
    writeLines(hdr, gz)
    for (mi in seq_along(markers)) {
      m <- markers[mi]
      pp <- prob_arr[mi,,]
      for (si in seq_along(progeny_samples)) {
        s <- progeny_samples[si]
        pv <- pp[si,]
        if (all(is.na(pv))) next
        writeLines(paste(c(m, s, format(pv, scientific=FALSE)), collapse="\t"), gz)
      }
    }

    # also write parents+progeny (useful for linkage mapping)
    if (length(parents_present) > 0) {
      gz2 <- gzfile(file.path(out_dir, "prob_wide_all.tsv.gz"), "wt")
      on.exit(close(gz2), add=TRUE)
      hdr <- paste(c("marker_id","sample_id", paste0("p", dosage_levels)), collapse="\t")
      writeLines(hdr, gz2)
      for (mi in seq_along(markers)) {
        m <- markers[mi]
        pp <- prob_arr_all[mi,,]
        for (si in seq_along(sample_order_all)) {
          s <- sample_order_all[si]
          pv <- pp[si,]
          if (all(is.na(pv))) next
          writeLines(paste(c(m, s, format(pv, scientific=FALSE)), collapse="\t"), gz2)
        }
      }
    }
  }

  # save legacy dosage.rds (keep)
  poly_dosage <- list(
    dosage=dosage_mat,
    marker=marker_info,
    sample=sample_info,
    meta=list(
      ploidy=ploidy,
      method=method,
      counts_rds=counts_rds,
      parent1_id=parent1_id,
      parent2_id=parent2_id,
      n_progeny=length(progeny_samples),
      progeny_samples=progeny_samples,
      extra_options=extra_options,
      output_prob=output_prob,
      call_prob_threshold=call_prob_threshold,
      prob_format=prob_format,
      created_at=as.character(Sys.time())
    )
  )
  saveRDS(poly_dosage, file.path(out_dir, "dosage.rds"))

  # build geno.rds (common schema)
  # markers
  markers_df <- NULL
  if (!is.null(marker_info)) {
    mi <- marker_info
    if (!("marker_id" %in% colnames(mi))) {
      if ("marker" %in% colnames(mi)) mi$marker_id <- as.character(mi$marker)
      else mi$marker_id <- markers
    }
    markers_df <- mi
  } else {
    markers_df <- data.frame(marker_id=markers, stringsAsFactors=FALSE)
  }

  # samples
  samples_df <- NULL
  if (!is.null(sample_info)) {
    si <- sample_info
    if (!("sample_id" %in% colnames(si))) {
      if ("ind" %in% colnames(si)) si$sample_id <- as.character(si$ind)
      else if ("sample" %in% colnames(si)) si$sample_id <- as.character(si$sample)
      else si$sample_id <- sample_order_all
    }
    samples_df <- si
  } else {
    samples_df <- data.frame(sample_id=sample_order_all, stringsAsFactors=FALSE)
  }

  # Ensure order and add role labels
  samples_df <- samples_df[samples_df$sample_id %in% sample_order_all, , drop=FALSE]
  samples_df <- samples_df[match(sample_order_all, samples_df$sample_id), , drop=FALSE]
  samples_df$role <- ifelse(samples_df$sample_id == parent1_id, "parent1",
                            ifelse(samples_df$sample_id == parent2_id, "parent2", "progeny"))

  # ensure matrices dimnames match marker_id/sample_id order
  # (dosage_mat already dimnames=markers/progeny_samples)
  geno <- list(
    version = "1.0.0",
    created_at = as.character(Sys.time()),
    ploidy = as.integer(ploidy),
    allele_coding = "ref_alt",
    caller = list(method="updog", method_version=NA_character_, params=params, seed=NA_integer_, notes="prob guaranteed via flexdog"),
    samples = samples_df,
    markers = markers_df,
    dosage = dosage_mat_all,
    prob = if (isTRUE(output_prob)) prob_arr_all else array(numeric(0), dim=c(0,0,0)),
    marker_stats = marker_stats,
    sample_stats = sample_stats,
    dosage_pp = dosage_pp_all,
    dosage_maxprob = dosage_maxprob_all,
    dosage_entropy = dosage_entropy_all,
    missing = list(dosage_value=NA_integer_, prob_value=NA_real_, definition="depth0_or_filtered_or_failed"),
    run_log = list(
      n_markers_input = markers_n,
      n_markers_output = markers_n,
      n_samples = samples_all_n,
      failed_markers = failed
    ),
    extra = list(
      counts_ref = NULL,
      counts_total = NULL,
      caller_internal = NULL,
      parents_present = parents_present,
      progeny_samples = progeny_samples
    )
  )
  saveRDS(geno, file.path(out_dir, "geno.rds"))

  meta <- list(
    module="poly_updog_call",
    ploidy=ploidy,
    method=method,
    output_prob=output_prob,
    call_prob_threshold=call_prob_threshold,
    prob_format=prob_format,
    counts_rds_in=basename(counts_rds),
    parent1_id=parent1_id,
    parent2_id=parent2_id,
    n_progeny=length(progeny_samples),
    progeny_file=if (nzchar(progeny_file)) basename(progeny_file) else "",
    dosage_tsv="dosage.tsv",
    dosage_all_tsv="dosage_all.tsv",
    progeny_ids_tsv="progeny_ids.tsv",
    dosage_rds="dosage.rds",
    geno_rds="geno.rds",
    marker_stats_tsv="marker_stats.tsv",
    sample_stats_tsv="sample_stats.tsv",
    marker_callrate_tsv="marker_callrate.tsv",
    plot="dosage_hist.png",
    plot_maxprob="maxprob_hist.png",
    prob_long_tsv_gz=if (isTRUE(output_prob) && prob_format=="long") "prob_long.tsv.gz" else "",
    prob_long_all_tsv_gz=if (isTRUE(output_prob) && prob_format=="long" && length(parents_present)>0) "prob_long_all.tsv.gz" else "",
    prob_wide_tsv_gz=if (isTRUE(output_prob) && prob_format=="wide") "prob_wide.tsv.gz" else "",
    prob_wide_all_tsv_gz=if (isTRUE(output_prob) && prob_format=="wide" && length(parents_present)>0) "prob_wide_all.tsv.gz" else "",
    default_table="dosage.tsv",
    default_plot="dosage_hist.png"
  )
  write_artifacts(meta)

  log("[poly_updog_call] done")
}, error=function(e){
  write_error(conditionMessage(e))
  write_artifacts(list(module="poly_updog_call", error=conditionMessage(e)))
  log("[poly_updog_call] ERROR: ", conditionMessage(e))
  quit(status=1)
})
