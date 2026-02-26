#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(jsonlite)
  library(data.table)
  library(TCC)
  library(Biostrings)
})

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 4 || args[1] != "--params" || args[3] != "--out") {
  cat("Usage: runner.R --params params.json --out out_dir\n")
  quit(status=2)
}
params_path <- args[2]
out_dir <- args[4]
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
log_file <- file.path(out_dir, "run.log")
sink(log_file, split=TRUE)

cat("[simulate_rnaseq] start\n")
cat("[simulate_rnaseq] params_path=", params_path, "\n", sep="")
cat("[simulate_rnaseq] out_dir=", out_dir, "\n", sep="")

p <- fromJSON(params_path)
`%||%` <- function(a, b) if (!is.null(a) && length(a)>0 && !(is.character(a) && a=="")) a else b

seed <- as.integer(p$seed %||% 1)
set.seed(seed)

if ((!(is.null(p$genename_file)) || p$genename_file != "")) genename_file <- p$genename_file
n_genes <- as.integer(p$n_genes %||% 20000)
n_use_org_at <- isTRUE(p$use_org_at_tair_db %||% FALSE)
n_org_at_keytype <- as.character(p$org_at_keytype %||% "TAIR")
if (is.na(n_org_at_keytype) || n_org_at_keytype == "") n_org_at_keytype <- "TAIR"
n_group <- as.integer(p$n_samples_per_group %||% 3)
#conditions <- p$conditions %||% c("A","B")
conditions <- p$conditions
if (is.list(conditions)) conditions <- unlist(conditions)
conditions <- as.character(conditions)
if (length(conditions) < 2) conditions <- c("A","B")

# Optional batch
n_batches <- as.integer(p$n_batches %||% 1)

# DE settings
prop_de <- as.numeric(p$prop_de %||% 0.1)
prop_de <- max(0, min(1, prop_de))
logfc_sd <- as.numeric(p$logfc_sd %||% 1.0)
logfc_mean <- as.numeric(p$logfc_mean %||% 0.0)

# Count model (NB)
mean_count <- as.numeric(p$mean_count %||% 50)
# NB size: larger -> less dispersion
nb_size <- as.numeric(p$nb_size %||% 10)
if (nb_size <= 0) nb_size <- 10

# Library size variation
libsize_mean <- as.numeric(p$libsize_mean %||% 1.0)  # multiplicative
libsize_sd <- as.numeric(p$libsize_sd %||% 0.25)

# Use TCC if available
use_tcc <- isTRUE(p$use_tcc %||% FALSE)
has_tcc <- FALSE
if (use_tcc) {
  has_tcc <- suppressWarnings(requireNamespace("TCC", quietly=TRUE))
  if (!has_tcc) cat("[simulate_rnaseq] WARN: TCC not available; falling back to internal NB generator\n")
}

# Optional: real Arabidopsis gene IDs from org.At.tair.db
gene_ids <- NULL
has_org_at <- FALSE
if (n_use_org_at) {
  has_org_at <- suppressWarnings(requireNamespace("org.At.tair.db", quietly=TRUE)) &&
                suppressWarnings(requireNamespace("AnnotationDbi", quietly=TRUE))
  if (!has_org_at) {
    cat("[simulate_rnaseq] WARN: use_org_at_tair_db=TRUE but org.At.tair.db/AnnotationDbi not available; falling back to synthetic gene IDs\n")
    n_use_org_at <- FALSE
  } else {
    cat("[simulate_rnaseq] Using org.At.tair.db gene IDs (keytype=", n_org_at_keytype, ")\n", sep="")
    # keys() order can be arbitrary; sort for reproducibility across runs
    gene_ids <- tryCatch({
      sort(AnnotationDbi::keys(org.At.tair.db::org.At.tair.db, keytype=n_org_at_keytype))
    }, error=function(e) {
      cat("[simulate_rnaseq] WARN: keys() failed for keytype=", n_org_at_keytype, ": ", e$message, "\n", sep="")
      NULL
    })
    if (is.null(gene_ids) || length(gene_ids) < 50) {
      cat("[simulate_rnaseq] WARN: Could not obtain sufficient gene IDs from org.At.tair.db; falling back to synthetic gene IDs\n")
      n_use_org_at <- FALSE
      gene_ids <- NULL
    } else {
      n_genes <- length(gene_ids)
      cat("[simulate_rnaseq] n_genes overridden to ", n_genes, " from org.At.tair.db\n", sep="")
    }
  }
}

gene_length <- NULL
#if (!(n_use_org_at) & (!(is.null(genename_file)) || genename_file != "")) {
if (nchar(genename_file) > 0) {
  if (grepl(".(fa|fasta)$", genename_file, perl = T)) {
    fasta <- readDNAStringSet(genename_file)
    gene_ids <- names(fasta)
    gene_length <- fasta@ranges@width
  } else {
    gene_ids <- fread(genename_file, header = F, data.table = F)[, 1]
    gene_length <- rep(1000, length(gene_ids))
  }
  gene_ids <- gsub(" .*$", "", gene_ids)
  n_genes <- length(gene_ids)
}

# Output folder layout
root <- file.path(out_dir, "dataset")
rnaseq_dir <- file.path(root, "rnaseq")
dir.create(rnaseq_dir, recursive=TRUE, showWarnings=FALSE)

# ---- generate design
samples <- unlist(lapply(seq_along(conditions), function(i) {
  paste0("S", i, "_", sprintf("%03d", seq_len(n_group)))
}))
cond_vec <- rep(conditions, each=n_group)

batch_vec <- rep("batch1", length(samples))
if (n_batches > 1) {
  # simple cyclic assignment
  batch_vec <- paste0("batch", ((seq_along(samples)-1) %% n_batches) + 1)
}

design_dt <- data.table(sample=samples, condition=cond_vec, batch=batch_vec)

# ---- generate counts
# We try to use TCC::simulateReadCounts if it exists; otherwise fall back.
counts <- NULL
truth_de <- NULL

if (use_tcc && has_tcc) {
  sim_fun <- NULL
  if ("simulateReadCounts" %in% getNamespaceExports("TCC")) sim_fun <- getExportedValue("TCC", "simulateReadCounts")
  if (!is.null(sim_fun)) {
    cat("[simulate_rnaseq] Using TCC::simulateReadCounts\n")

    # IMPORTANT:
    # Many TCC builds use the built-in dataset `arab` (from the NBPSeq package)
    # to estimate parameters for the NB simulation. If `arab` is missing,
    # simulateReadCounts() can error with: "object 'arab' not found".
    # To make GUI execution robust across environments, we proactively load it.
    if (!suppressWarnings(requireNamespace("NBPSeq", quietly=TRUE))) {
      cat("[simulate_rnaseq] WARN: NBPSeq not available; cannot load dataset 'arab'. Falling back to internal NB generator.\n")
      sim_fun <- NULL
    } else {
      ok_arab <- tryCatch({
        data("arab", package="NBPSeq", envir=.GlobalEnv)
        exists("arab", envir=.GlobalEnv)
      }, error=function(e) FALSE)
      if (!isTRUE(ok_arab)) {
        cat("[simulate_rnaseq] WARN: Could not load dataset 'arab' from NBPSeq; falling back to internal NB generator.\n")
        sim_fun <- NULL
      }
    }

    if (is.null(sim_fun)) {
      cat("[simulate_rnaseq] WARN: TCC cannot run because required dataset 'arab' is unavailable; falling back to internal NB generator.\n")
    } else {

    # TCC's simulateReadCounts signature differs by version.
    # - Newer docs (as provided by the user): simulateReadCounts(Ngene, PDEG, DEG.assign, DEG.foldchange, replicates, group, fc.matrix)
    # - Older (some forks): simulateReadCounts(numGene, numDEG, numSample, DEGfold, randomSeed, ...)
    # We inspect formals() and call the matching signature.
    fmls <- tryCatch(names(formals(sim_fun)), error=function(e) character(0))

    k <- length(conditions)
    reps <- rep(n_group, k)
    deg_assign <- rep(1/k, k)
    # fold-change per group; must be positive
    deg_fc <- pmax(0.05, 2^(rnorm(k, mean=logfc_mean, sd=logfc_sd)))

    tcc_obj <- NULL
    if ("Ngene" %in% fmls && "PDEG" %in% fmls) {
      # --- Documented Bioconductor-style interface
      tcc_obj <- tryCatch({
        sim_fun(
          Ngene = n_genes,
          PDEG = prop_de,
          DEG.assign = deg_assign,
          DEG.foldchange = deg_fc,
          replicates = reps
        )
      }, error=function(e) {
        cat("[simulate_rnaseq] WARN: TCC::simulateReadCounts (Ngene/PDEG) failed: ", e$message, "\n", sep="")
        NULL
      })
    } else if ("numGene" %in% fmls || "numDEG" %in% fmls) {
      # --- Legacy / alternate interface (best-effort)
      tcc_obj <- tryCatch({
        sim_fun(
          numGene = n_genes,
          numDEG  = as.integer(round(n_genes * prop_de)),
          numSample = reps,
          DEGfold = median(deg_fc),
          randomSeed = seed
        )
      }, error=function(e) {
        cat("[simulate_rnaseq] WARN: TCC::simulateReadCounts (legacy) failed: ", e$message, "\n", sep="")
        NULL
      })
    } else {
      cat("[simulate_rnaseq] WARN: Unrecognized simulateReadCounts signature; falling back\n")
    }

    if (!is.null(tcc_obj) && !is.null(tcc_obj$count)) {
      counts <- as.matrix(tcc_obj$count)
      # Some TCC versions may return a data.frame-like object that becomes
      # a character matrix after as.matrix(). Coerce to numeric safely.
      if (!is.numeric(counts) && !is.integer(counts)) {
        suppressWarnings(storage.mode(counts) <- "numeric")
      }
      # Replace any NA (written as blanks by fwrite) with 0 for stable test datasets.
      na_n <- sum(is.na(counts))
      if (na_n > 0) {
        cat("[simulate_rnaseq] WARN: counts contains ", na_n, " NA values; replacing with 0\n", sep="")
        counts[is.na(counts)] <- 0
      }
      
      ### Adjust to libsize_mean
      sample_mean <- mean(apply(counts, 2, sum))
      lib_size <- (sample_mean / libsize_mean) 
      counts <- counts / lib_size
      
      storage.mode(counts) <- "integer"
      # TCC returns genes x samples (group order). We generate samples in the same order.
      if (ncol(counts) == length(samples)) colnames(counts) <- samples
      if (is.null(rownames(counts))) {
        rownames(counts) <- paste0("gene", sprintf("%05d", seq_len(nrow(counts))))
      }

      # truth from TCC (preferred)
      # trueDEG: 0 non-DEG, 1 DEG up in group1, 2 up in group2, ...
      trueDEG <- NULL
      if (!is.null(tcc_obj$simulation) && !is.null(tcc_obj$simulation$trueDEG)) trueDEG <- tcc_obj$simulation$trueDEG
      if (!is.null(trueDEG) && length(trueDEG) == nrow(counts)) {
        truth_de <- data.table(
          gene = rownames(counts),
          is_de = as.integer(trueDEG != 0),
          de_group = as.integer(trueDEG),
          log2fc = 0
        )
        # encode log2 fold-change of the "up" group (as provided to simulateReadCounts)
        for (g in seq_len(k)) {
          truth_de[de_group == g, log2fc := log2(deg_fc[g])]
        }
      } else {
        # fallback truth if TCC doesn't provide it
        n_de <- as.integer(round(n_genes * prop_de))
        de_idx <- if (n_de>0) sample(seq_len(n_genes), n_de) else integer(0)
        truth_de <- data.table(gene=rownames(counts), is_de=0L, de_group=0L, log2fc=0)
        if (n_de>0) {
          grp <- sample(seq_len(k), n_de, replace=TRUE)
          truth_de[de_idx, `:=`(is_de=1L, de_group=grp, log2fc=log2(deg_fc[grp]))]
        }
      }
    }

    } # end else(sim_fun)
  } else {
    cat("[simulate_rnaseq] WARN: TCC does not export simulateReadCounts; falling back\n")
  }
}

if (is.null(counts)) {
  cat("[simulate_rnaseq] Using internal NB generator\n")
  n_samp <- length(samples)
  gene_ids_local <- if (!is.null(gene_ids)) gene_ids else paste0("gene", sprintf("%05d", seq_len(n_genes)))

  # baseline mean per gene: log-normal around mean_count
  base_mu <- rlnorm(n_genes, meanlog=log(mean_count + 1e-9), sdlog=1.0)

  # choose DE genes per non-reference condition vs first condition
  n_de <- as.integer(round(n_genes * prop_de))
  de_idx <- if (n_de>0) sample(seq_len(n_genes), n_de) else integer(0)

  # logFC for DE genes (relative to condition1)
  logfc <- if (n_de>0) rnorm(n_de, mean=logfc_mean, sd=logfc_sd) else numeric(0)

  # library size factor per sample
  lib_f <- rlnorm(n_samp, meanlog=log(libsize_mean), sdlog=libsize_sd)

  # batch factor (multiplicative)
  batch_levels <- unique(batch_vec)
  batch_f_map <- setNames(rlnorm(length(batch_levels), meanlog=0, sdlog=0.1), batch_levels)
  batch_f <- batch_f_map[batch_vec]

  counts <- matrix(0L, nrow=n_genes, ncol=n_samp)
  rownames(counts) <- gene_ids_local
  colnames(counts) <- samples

  # reference condition = conditions[1]
  ref <- conditions[1]
  for (j in seq_len(n_samp)) {
    cond <- cond_vec[j]
    mu <- base_mu * lib_f[j] * batch_f[j]
    if (cond != ref && n_de>0) {
      # apply DE to selected genes for any non-reference condition
      mu[de_idx] <- mu[de_idx] * exp(logfc)
    }
    # NB sampling: Var = mu + mu^2/size
    counts[, j] <- rnbinom(n_genes, mu=mu, size=nb_size)
  }

  truth_de <- data.table(gene=gene_ids_local, is_de=0L, logFC=0)
  if (n_de>0) truth_de[de_idx, `:=`(is_de=1L, logFC=logfc)]
}

# Final sanity: ensure integer matrix and no NA
if (!is.null(counts)) {
  if (!is.numeric(counts) && !is.integer(counts)) {
    suppressWarnings(storage.mode(counts) <- "numeric")
  }
  na_n <- sum(is.na(counts))
  if (na_n > 0) {
    cat("[simulate_rnaseq] WARN: counts contains ", na_n, " NA values at final check; replacing with 0\n", sep="")
    counts[is.na(counts)] <- 0
  }
  storage.mode(counts) <- "integer"
}

# If counts were generated by TCC path, apply org.At gene IDs if requested
#if (!is.null(counts) && n_use_org_at && !is.null(gene_ids)) {
if (!is.null(counts) & !is.null(gene_ids)) {
  if (nrow(counts) == length(gene_ids)) {
    rownames(counts) <- gene_ids
    if (!is.null(truth_de) && nrow(truth_de) == length(gene_ids)) {
      truth_de[, gene := gene_ids]
    }
  } else {
    cat("[simulate_rnaseq] WARN: TCC-generated nrow(counts) != length(gene_ids); keeping synthetic gene IDs\n")
  }
}

# Write outputs
counts_path <- file.path(rnaseq_dir, "counts.tsv")
design_path <- file.path(rnaseq_dir, "design.tsv")
truth_path <- file.path(rnaseq_dir, "truth_de.tsv")
genelist_path <- file.path(rnaseq_dir, "gene_list.txt")
genelength_path <- file.path(rnaseq_dir, "gene_length.tsv")

# counts as genes x samples with first column gene
counts_dt <- as.data.table(counts, keep.rownames="gene")
fwrite(counts_dt, counts_path, sep="\t", quote = F, na = "NA")
fwrite(design_dt, design_path, sep="\t", quote = F, na = "NA")
fwrite(truth_de, truth_path, sep="\t", quote = F, na = "NA")
fwrite(as.data.frame(gene_ids), genelist_path, sep="\t", quote = F, na = "NA", row.names = F, col.names = F)
if (is.null(gene_length)) gene_length <- rep(1000, length(gene_ids))
fwrite(data.frame(gene = gene_ids, length = gene_length), genelength_path, sep="\t", quote = F, na = "NA", row.names = F, col.names = T)

manifest <- list(
  generator="simulate_rnaseq",
  version="0.2.0",
  seed=seed,
  n_genes=n_genes,
  use_org_at_tair_db=as.logical(n_use_org_at),
  org_at_keytype=if (n_use_org_at) n_org_at_keytype else NULL,
  conditions=conditions,
  n_samples_per_group=n_group,
  n_batches=n_batches,
  prop_de=prop_de,
  logfc_mean=logfc_mean,
  logfc_sd=logfc_sd,
  mean_count=mean_count,
  nb_size=nb_size,
  libsize_mean=libsize_mean,
  libsize_sd=libsize_sd,
  use_tcc=as.logical(use_tcc && has_tcc),
  outputs=list(
    counts=counts_path,
    design=design_path,
    truth_de=truth_path
  )
)
writeLines(toJSON(manifest, auto_unbox=TRUE, pretty=TRUE), file.path(rnaseq_dir, "manifest.json"))

cat("[simulate_rnaseq] wrote\n")
cat("  ", counts_path, "\n")
cat("  ", design_path, "\n")
cat("  ", truth_path, "\n")

cat("[simulate_rnaseq] done\n")
