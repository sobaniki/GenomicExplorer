#!/usr/bin/env Rscript
# rnaseq DEG + optional GOstats
suppressWarnings(suppressMessages({
  library(jsonlite)
  library(data.table)
}))

fail <- function(msg, code=1) {
  cat("[rnaseq_deg][ERROR] ", msg, "\n", sep="")
  quit(status=code)
}

as_bool <- function(x) {
  if (is.null(x)) return(FALSE)
  if (is.logical(x)) return(isTRUE(x))
  if (is.character(x)) return(tolower(x) %in% c("1","true","t","yes","y"))
  if (is.numeric(x)) return(x != 0)
  FALSE
}

read_table_auto <- function(path) {
  # fread handles tsv/csv; fall back to read.csv if needed
  tryCatch({
    data.table::fread(path, data.table=FALSE)
  }, error=function(e) {
    read.csv(path, check.names=FALSE, stringsAsFactors=FALSE)
  })
}

save_plot_png <- function(path, expr) {
  png(path, width=1000, height=700, res=120)
  on.exit(dev.off(), add=TRUE)
  expr
}

make_qc_plots <- function(out_dir, counts_mat, expr_mat, des_df, condition_col, method_label) {
  if (is.null(expr_mat) || ncol(expr_mat) < 2) return(invisible(NULL))
  # library size
  lib_png <- file.path(out_dir, "qc_libsize.png")
  save_plot_png(lib_png, {
    lib <- colSums(counts_mat)
    barplot(lib, las=2, main=paste0("Library size (", method_label, ")"), ylab="sum(counts)")
  })

  # PCA (on samples)
  pca_png <- file.path(out_dir, "qc_pca.png")
  save_plot_png(pca_png, {
    X <- t(expr_mat)
    X <- X[, colSums(is.na(X)) == 0, drop=FALSE]
    X[is.na(X)] <- 0
    pr <- prcomp(X, scale.=TRUE)
    x <- pr$x[,1]; y <- pr$x[,2]
    grp <- NULL
    if (!is.null(des_df) && condition_col %in% colnames(des_df)) {
      grp <- as.factor(des_df[[condition_col]])
    } else {
      grp <- factor(rep("sample", length(x)))
    }
    plot(x, y, pch=19, xlab="PC1", ylab="PC2", main=paste0("PCA (", method_label, ")"))
    text(x, y, labels=rownames(pr$x), pos=3, cex=0.7)
    legend("topright", legend=levels(grp), pch=19, bty="n")
  })

  # Sample distance heatmap
  dist_png <- file.path(out_dir, "qc_sample_distance.png")
  save_plot_png(dist_png, {
    d <- dist(t(expr_mat))
    m <- as.matrix(d)
    heatmap(m, symm=TRUE, main=paste0("Sample distance (", method_label, ")"))
  })

  invisible(list(qc_libsize=lib_png, qc_pca=pca_png, qc_sample_distance=dist_png))
}



args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 4 || args[1] != "--params" || args[3] != "--out") {
  cat("Usage: runner.R --params params.json --out out_dir\n")
  quit(status=2)
}

params_path <- args[2]
out_dir <- args[4]
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

cat("[rnaseq_deg] start\n")
cat("[rnaseq_deg] params_path=", params_path, "\n", sep="")
cat("[rnaseq_deg] out_dir=", out_dir, "\n", sep="")

params <- tryCatch(jsonlite::fromJSON(params_path), error=function(e) fail(paste("Failed to parse params.json:", e$message)))
`%||%` <- function(a, b) { if (is.null(a) || length(a)==0 || (is.character(a)&&a=="")) b else a }

## Accept both legacy keys (counts_path/design_path/condition_A) and
## GUI keys (counts_matrix/design_table/cond_a/cond_b) for robustness.
counts_path <- params$counts_path %||% params$counts_matrix %||% params$counts
design_path <- params$design_path %||% params$design_table %||% params$design
method <- tolower(params$method %||% "edger")
sample_col <- params$sample_col %||% "sample"
condition_col <- params$condition_col %||% "condition"
cond_A <- params$condition_A %||% params$cond_a %||% params$condA
cond_B <- params$condition_B %||% params$cond_b %||% params$condB
fdr_cutoff <- as.numeric(params$fdr_cutoff %||% params$fdr %||% 0.05)
lfc_cutoff <- as.numeric(params$lfc_cutoff %||% params$lfc %||% 0.0)
do_go <- as_bool(params$do_go %||% (params$go_mode == "gostats"))
orgdb_pkg <- params$orgdb_package %||% params$orgdb %||% "org.Hs.eg.db"
keytype <- params$keytype %||% "ENTREZID"
ontology <- toupper(params$ontology %||% "BP")


# Advanced design/contrast (optional)
design_formula <- params$design_formula %||% params$formula
contrast_str <- params$contrast %||% params$contrast_string
reference_level <- params$reference_level %||% params$ref_level
advanced <- as_bool(params$advanced %||% FALSE)
subset_to_conditions <- as_bool(params$subset_to_conditions %||% TRUE)

# Step3: probe-only mode (emit design/contrast info and exit)
probe_only <- as_bool(params$probe_only %||% params$dry_run %||% FALSE)

# Step3: safer selection helpers
coef_name <- params$coef_name %||% params$coefName
deseq2_results_name <- params$deseq2_results_name %||% params$deseq2_name

# method-specific (Step1)
edger_test <- toupper(params$edger_test %||% params$test_type %||% "QLF")  # QLF or LRT
edger_robust <- as_bool(params$edger_robust %||% params$robust_dispersion %||% FALSE)

deseq2_fitType <- tolower(params$deseq2_fitType %||% params$fitType %||% "parametric")
deseq2_test <- toupper(params$deseq2_test %||% params$deseq2_test_type %||% "WALD")  # WALD or LRT
deseq2_independentFiltering <- as_bool(params$independentFiltering %||% TRUE)
deseq2_cooksCutoff <- params$cooksCutoff %||% TRUE

limma_quality_weights <- as_bool(params$limma_quality_weights %||% FALSE)

# Step2: additional method options
edger_trend_method <- tolower(params$edger_trend_method %||% params$trend_method %||% "")
edger_use_treat <- as_bool(params$edger_use_treat %||% params$treat_lfc %||% FALSE)
edger_treat_lfc <- as.numeric(params$edger_treat_lfc %||% params$treat_lfc_threshold %||% params$treat_threshold %||% 0.0)

deseq2_use_shrink <- as_bool(params$deseq2_use_shrink %||% params$lfc_shrink %||% FALSE)
deseq2_shrink_type <- tolower(params$deseq2_shrink_type %||% params$shrink_method %||% "apeglm")
deseq2_qc_transform <- tolower(params$deseq2_qc_transform %||% params$qc_transform %||% "vst")

limma_use_treat <- as_bool(params$limma_use_treat %||% params$treat_lfc %||% FALSE)
limma_treat_lfc <- as.numeric(params$limma_treat_lfc %||% params$treat_lfc_threshold %||% params$treat_threshold %||% 0.0)



if (is.null(counts_path) || !file.exists(counts_path)) fail("counts_path not found")
if (is.null(design_path) || !file.exists(design_path)) fail("design_path not found")
# Load counts
cnt <- read_table_auto(counts_path)
if (ncol(cnt) < 3) fail("Counts matrix must have gene_id column + >=2 sample columns")

gene_id_col <- params$gene_id_col %||% colnames(cnt)[1]
if (!(gene_id_col %in% colnames(cnt))) {
  gene_id_col <- colnames(cnt)[1]
}
genes <- as.character(cnt[[gene_id_col]])
cnt[[gene_id_col]] <- NULL

# Ensure numeric counts
for (j in seq_len(ncol(cnt))) {
  cnt[[j]] <- suppressWarnings(as.integer(cnt[[j]]))
}
cnt[is.na(cnt)] <- 0L

# Load design
des <- read_table_auto(design_path)
if (!(sample_col %in% colnames(des))) fail(paste0("design table lacks sample_col=", sample_col))
if (!(condition_col %in% colnames(des))) fail(paste0("design table lacks condition_col=", condition_col))
des[[sample_col]] <- as.character(des[[sample_col]])
des[[condition_col]] <- as.character(des[[condition_col]])

# Determine modes (legacy two-group vs advanced formula/contrast)
#use_simple_cmp <- !(is.null(cond_A) || is.null(cond_B) || cond_A=="" || cond_B=="")
use_simple_cmp <- F
if (!is.null(design_formula) && is.character(design_formula) && nchar(design_formula) > 0) {
  # ensure it starts with "~"
  # Use POSIX character class to avoid invalid escape sequences like "\\s" on some R builds
  if (!grepl("^[[:space:]]*~", design_formula)) design_formula <- paste0("~", design_formula)
} else {
  design_formula <- NULL
}
if (!is.null(contrast_str) && is.character(contrast_str) && nchar(contrast_str) > 0) {
  contrast_str <- contrast_str
} else {
  contrast_str <- NULL
}
if (!is.null(reference_level) && is.character(reference_level) && nchar(reference_level) > 0) {
  reference_level <- reference_level
} else {
  reference_level <- NULL
}

# auto-enable advanced if formula/contrast is provided
# if (isTRUE(advanced) || !is.null(design_formula) || !is.null(contrast_str)) {
#   advanced <- TRUE
# } else {
#   advanced <- FALSE
# }
advanced <- TRUE

# if (!advanced && !use_simple_cmp) {
#   fail("condition_A and condition_B must be set (two-group comparison) in Simple mode")
# }

# Optionally subset to two conditions (useful even in Advanced mode)
des_use <- des
if (use_simple_cmp && subset_to_conditions) {
  des_use <- des[des[[condition_col]] %in% c(cond_A, cond_B), , drop=FALSE]
}
if (nrow(des_use) < 2) fail("Not enough samples after filtering design table")

# Align samples between counts and design
sample_names <- colnames(cnt)
keep_samples <- intersect(sample_names, as.character(des_use[[sample_col]]))
if (length(keep_samples) < 2) fail("No overlapping samples between counts columns and design sample IDs")
# reorder by design
des_use <- des_use[match(keep_samples, as.character(des_use[[sample_col]])), , drop=FALSE]
cnt2 <- cnt[, keep_samples, drop=FALSE]

counts_mat <- as.matrix(cnt2)
rownames(counts_mat) <- make.unique(genes)
mode(counts_mat) <- "integer"

# Make design columns usable for formulas
for (cc in colnames(des_use)) {
  if (is.character(des_use[[cc]])) des_use[[cc]] <- factor(des_use[[cc]])
}

# Ensure condition factor and reference level (if provided)
if (condition_col %in% colnames(des_use)) {
  if (!is.factor(des_use[[condition_col]])) des_use[[condition_col]] <- factor(as.character(des_use[[condition_col]]))
  if (is.null(reference_level) && use_simple_cmp) reference_level <- cond_B
  if (!is.null(reference_level) && reference_level %in% levels(des_use[[condition_col]])) {
    des_use[[condition_col]] <- stats::relevel(des_use[[condition_col]], ref=reference_level)
  }
}

# Simple-mode factor for two-group comparisons (kept for backward compatibility)
cond <- NULL
if (use_simple_cmp && condition_col %in% colnames(des_use)) {
  cond <- factor(as.character(des_use[[condition_col]]), levels=c(cond_B, cond_A)) # coef2 = A vs B
}

# Step3: Probe mode (emit model matrix / levels to help build contrasts safely)
if (probe_only) {
  df_use <- design_formula
  if (is.null(df_use) || !nzchar(df_use)) {
    df_use <- paste0("~ ", condition_col)
  }
  if (!grepl("^[[:space:]]*~", df_use)) df_use <- paste0("~", df_use)
  mm <- model.matrix(as.formula(df_use), data=des_use)
  info <- list(
    n_samples = ncol(counts_mat),
    samples = colnames(counts_mat),
    design_formula = df_use,
    model_matrix_cols = colnames(mm),
    condition_col = condition_col,
    condition_levels = if (condition_col %in% colnames(des_use) && is.factor(des_use[[condition_col]])) levels(des_use[[condition_col]]) else NULL
  )
  # Optional: try to compute DESeq2 results names (may be slow for huge data; best-effort)
  if (requireNamespace("DESeq2", quietly=TRUE)) {
    try({
      suppressWarnings(suppressMessages(library(DESeq2)))
      dds0 <- DESeqDataSetFromMatrix(countData=counts_mat, colData=des_use, design=as.formula(df_use))
      dds0 <- DESeq(dds0, test="Wald", quiet=TRUE)
      info$deseq2_results_names <- resultsNames(dds0)
    }, silent=TRUE)
  }
  jsonlite::write_json(info, file.path(out_dir, "design_info.json"), auto_unbox=TRUE, pretty=TRUE)
  cat("[rnaseq_deg] probe_only: wrote design_info.json\n")
  quit(status=0)
}

deg_out <- file.path(out_dir, "deg.tsv")
ma_png <- file.path(out_dir, "ma_plot.png")
vol_png <- file.path(out_dir, "volcano_plot.png")
norm_counts_out <- file.path(out_dir, "normalized_counts.tsv")

res_df <- NULL
qc_expr <- NULL
qc_method <- NULL
if (method == "edger") {
  suppressWarnings(suppressMessages({
    if (!requireNamespace("edgeR", quietly=TRUE)) fail("edgeR not installed")
    library(edgeR)
  }))

  y <- DGEList(counts=counts_mat)

  # Build design matrix (Advanced) or two-group design (Simple)
  design <- NULL
  if (advanced && !is.null(design_formula)) {
    design <- model.matrix(as.formula(design_formula), data=des_use)
    keep <- filterByExpr(y, design=design)
  } else {
    if (is.null(cond)) fail("edgeR: condition_A/B required in Simple mode (or provide design_formula)")
    y$samples$group <- cond
    design <- model.matrix(~cond)
    keep <- filterByExpr(y, group=cond)
  }

  y <- y[keep, , keep.lib.sizes=FALSE]
  y <- calcNormFactors(y)

  # Dispersion estimation
  if (!is.null(edger_trend_method) && nzchar(edger_trend_method)) {
    y <- estimateDisp(y, design, robust=edger_robust, trend.method=edger_trend_method)
  } else {
    y <- estimateDisp(y, design, robust=edger_robust)
  }

  # Determine test target
  coef_idx <- NULL
  if (!is.null(params$coef) && !is.null(suppressWarnings(as.integer(params$coef)))) {
    coef_idx <- as.integer(params$coef)
  }
  # Step3: allow selecting coefficient by name (safer than numeric index)
  if (is.null(coef_idx) && !is.null(coef_name) && is.character(coef_name) && nzchar(coef_name)) {
    cn <- trimws(coef_name)
    if (cn %in% colnames(design)) {
      coef_idx <- which(colnames(design) == cn)[1]
    }
  }
  if (is.null(coef_idx) && use_simple_cmp && !is.null(cond)) {
    # default for ~cond is coef=2 (A vs B if ref=B)
    coef_idx <- 2L
  }
  if (is.null(coef_idx) && use_simple_cmp && !is.null(cond_A) && condition_col %in% colnames(des_use)) {
    # try to match coefficient name for cond_A (treatment coding)
    guess <- make.names(paste0(condition_col, cond_A))
    if (guess %in% colnames(design)) coef_idx <- which(colnames(design) == guess)[1]
  }
  if (is.null(coef_idx) | is.na(coef_idx)) coef_idx <- 2L

  # Contrast vector if provided
  contrast_vec <- NULL
  if (!is.null(contrast_str)) {
    if (!requireNamespace("limma", quietly=TRUE)) {
      fail("edgeR: 'contrast' requires limma (makeContrasts). Please install limma or specify coef.")
    }
    contrast_mat <- limma::makeContrasts(contrasts=contrast_str, levels=colnames(design))
    contrast_vec <- as.numeric(contrast_mat)
    names(contrast_vec) <- rownames(contrast_mat)
  }

  # Fit + test
  if (edger_test == "LRT") {
    fit <- glmFit(y, design)
    if (edger_use_treat && edger_treat_lfc > 0) {
      if (!is.null(contrast_vec)) {
        ttst <- glmTreat(fit, contrast=contrast_vec, lfc=edger_treat_lfc)
      } else {
        ttst <- glmTreat(fit, coef=coef_idx, lfc=edger_treat_lfc)
      }
    } else {
      if (!is.null(contrast_vec)) {
        ttst <- glmLRT(fit, contrast=contrast_vec)
      } else {
        ttst <- glmLRT(fit, coef=coef_idx)
      }
    }
  } else {
    fit <- glmQLFit(y, design, robust=edger_robust)
    if (edger_use_treat && edger_treat_lfc > 0) {
      if (!is.null(contrast_vec)) {
        ttst <- glmTreat(fit, contrast=contrast_vec, lfc=edger_treat_lfc)
      } else {
        ttst <- glmTreat(fit, coef=coef_idx, lfc=edger_treat_lfc)
      }
    } else {
      if (!is.null(contrast_vec)) {
        #ttst <- glmQLFTest(fit)
        ttst <- glmQLFTest(fit, contrast=contrast_vec)
      } else {
        #ttst <- glmQLFTest(fit)
        ttst <- glmQLFTest(fit, coef=coef_idx)
      }
    }
  }
  tt <- topTags(ttst, n=Inf)$table
  tt$gene <- rownames(tt)
  res_df <- tt[, c("gene","logFC","logCPM","PValue","FDR")]
  colnames(res_df) <- c("gene","log2FoldChange","logCPM","pvalue","padj")

  # normalized counts (CPM)
  cpm_mat <- cpm(y, log=FALSE, prior.count=0)
  cpm_df <- data.frame(gene=rownames(cpm_mat), cpm_mat, check.names=FALSE)
  write.table(cpm_df, norm_counts_out, sep="	", quote=FALSE, row.names=FALSE)

  # QC expression matrix (logCPM)
  qc_expr <- cpm(y, log=TRUE, prior.count=1)
  qc_method <- "edgeR"

  # plots
  save_plot_png(ma_png, {
    with(res_df, {
      plot(logCPM, log2FoldChange, pch=16, cex=0.5, xlab="logCPM", ylab="log2FC", main="MA-like plot (edgeR)")
      abline(h=c(-lfc_cutoff, lfc_cutoff), lty=2)
    })
  })
  save_plot_png(vol_png, {
    with(res_df, {
      plot(log2FoldChange, -log10(pmax(pvalue, 1e-300)), pch=16, cex=0.5,
           xlab="log2FC", ylab="-log10(p)", main="Volcano (edgeR)")
      abline(v=c(-lfc_cutoff, lfc_cutoff), lty=2)
      abline(h=-log10(fdr_cutoff), lty=2)
    })
  })
} else if (method == "deseq2" || method == "deseq") {
  suppressWarnings(suppressMessages({
    if (!requireNamespace("DESeq2", quietly=TRUE)) fail("DESeq2 not installed")
    library(DESeq2)
  }))

  # Build colData
  coldata <- des_use
  rownames(coldata) <- as.character(des_use[[sample_col]])
  # ensure order matches counts
  coldata <- coldata[colnames(counts_mat), , drop=FALSE]

  # Simple mode uses a synthetic 'cond' column; Advanced mode uses design_formula
  dds_design <- NULL
  if (advanced && !is.null(design_formula)) {
    dds_design <- as.formula(design_formula)
  } else {
    if (is.null(cond)) fail("DESeq2: condition_A/B required in Simple mode (or provide design_formula)")
    coldata$cond <- cond
    dds_design <- as.formula("~ cond")
  }

  # Ensure characters are factors (DESeq2 prefers factors)
  for (cc in colnames(coldata)) {
    if (is.character(coldata[[cc]])) coldata[[cc]] <- factor(coldata[[cc]])
  }

  dds <- DESeqDataSetFromMatrix(countData=counts_mat, colData=coldata, design=dds_design)

  # filter low counts
  keep <- rowSums(counts(dds)) >= as.numeric(params$min_total_count %||% 10)
  dds <- dds[keep,]

  # Fit
  reduced_formula <- params$reduced_formula %||% params$reduced
  if (deseq2_test == "LRT") {
    if (is.null(reduced_formula) || reduced_formula == "") {
      fail("DESeq2 LRT requires 'reduced_formula' (e.g., ~ batch).")
    }
    if (!grepl("^[[:space:]]*~", reduced_formula)) reduced_formula <- paste0("~", reduced_formula)
    dds <- DESeq(dds, test="LRT", reduced=as.formula(reduced_formula), fitType=deseq2_fitType, quiet=TRUE)
  } else {
    dds <- DESeq(dds, test="Wald", fitType=deseq2_fitType, quiet=TRUE)
  }

  # Decide results extraction
  res <- NULL
  # Step3: if results name is explicitly provided, use it (safest for complex designs)
  if (!is.null(deseq2_results_name) && is.character(deseq2_results_name) && nzchar(deseq2_results_name)) {
    nm0 <- trimws(deseq2_results_name)
    res <- results(dds, name=nm0,
                   independentFiltering=deseq2_independentFiltering, cooksCutoff=deseq2_cooksCutoff)
    # also treat as the preferred coef for shrinkage
    contrast_str <- paste0("name=", nm0)
  } else if (!is.null(contrast_str)) {
    cs <- gsub("[[:space:]]+", "", contrast_str)
    if (grepl("^name=", cs, ignore.case=TRUE)) {
      nm <- sub("^name=", "", cs, ignore.case=TRUE)
      res <- results(dds, name=nm, independentFiltering=deseq2_independentFiltering, cooksCutoff=deseq2_cooksCutoff)
    } else {
      parts <- unlist(strsplit(contrast_str, "[,;]"))
      parts <- parts[nchar(trimws(parts))>0]
      if (length(parts) == 3) {
        res <- results(dds, contrast=c(trimws(parts[1]), trimws(parts[2]), trimws(parts[3])),
                       independentFiltering=deseq2_independentFiltering, cooksCutoff=deseq2_cooksCutoff)
      } else if (use_simple_cmp && condition_col %in% colnames(coldata)) {
        res <- results(dds, contrast=c(condition_col, cond_A, cond_B),
                       independentFiltering=deseq2_independentFiltering, cooksCutoff=deseq2_cooksCutoff)
      } else {
        # fall back to name
        res <- results(dds, name=contrast_str, independentFiltering=deseq2_independentFiltering, cooksCutoff=deseq2_cooksCutoff)
      }
    }
  } else if (use_simple_cmp && condition_col %in% colnames(coldata)) {
    res <- results(dds, contrast=c(condition_col, cond_A, cond_B),
                   independentFiltering=deseq2_independentFiltering, cooksCutoff=deseq2_cooksCutoff)
  } else if (!is.null(cond)) {
    res <- results(dds, contrast=c("cond", cond_A, cond_B),
                   independentFiltering=deseq2_independentFiltering, cooksCutoff=deseq2_cooksCutoff)
  } else {
    # If no contrast is given, take the last coefficient
    rn <- resultsNames(dds)
    if (length(rn) < 2) fail("DESeq2: cannot determine a coefficient. Provide 'contrast' or 'cond_a/cond_b'.")
    res <- results(dds, name=rn[length(rn)], independentFiltering=deseq2_independentFiltering, cooksCutoff=deseq2_cooksCutoff)
  }

  # Optional LFC shrinkage (Step2)
  if (deseq2_use_shrink) {
    shrinked <- NULL
    cs_trim <- NULL
    if (!is.null(contrast_str)) cs_trim <- trimws(contrast_str)
    if (!is.null(cs_trim) && grepl("^name[[:space:]]*=", cs_trim, ignore.case=TRUE)) {
      nm <- sub("^name[[:space:]]*=", "", cs_trim, ignore.case=TRUE)
      shrinked <- tryCatch({ lfcShrink(dds, coef=nm, type=deseq2_shrink_type) }, error=function(e) NULL)
    }
    if (is.null(shrinked) && use_simple_cmp && condition_col %in% colnames(coldata)) {
      shrinked <- tryCatch({ lfcShrink(dds, contrast=c(condition_col, cond_A, cond_B), type=deseq2_shrink_type) }, error=function(e) NULL)
    }
    if (is.null(shrinked) && !is.null(cond)) {
      shrinked <- tryCatch({ lfcShrink(dds, contrast=c("cond", cond_A, cond_B), type=deseq2_shrink_type) }, error=function(e) NULL)
    }
    if (is.null(shrinked)) shrinked <- res
    res <- shrinked
  }

  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  # keep baseMean for MA plot
  keep_cols <- c("gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")
  keep_cols <- keep_cols[keep_cols %in% colnames(res_df)]
  res_df <- res_df[, keep_cols, drop=FALSE]

  # normalized counts
  nc <- counts(dds, normalized=TRUE)
  nc_df <- data.frame(gene=rownames(nc), nc, check.names=FALSE)
  write.table(nc_df, norm_counts_out, sep="	", quote=FALSE, row.names=FALSE)

  # QC expression matrix (Step2: choose transform)
  qc_method <- NULL
  if (deseq2_qc_transform == "rlog") {
    tr <- tryCatch({ DESeq2::rlog(dds, blind=TRUE) }, error=function(e) NULL)
    if (!is.null(tr)) {
      qc_expr <- assay(tr)
      qc_method <- "DESeq2 (rlog)"
    }
  }
  if (is.null(qc_method) && deseq2_qc_transform == "vst") {
    tr <- tryCatch({ DESeq2::vst(dds, blind=TRUE) }, error=function(e) NULL)
    if (!is.null(tr)) {
      qc_expr <- assay(tr)
      qc_method <- "DESeq2 (vst)"
    }
  }
  if (is.null(qc_method)) {
    qc_expr <- log2(nc + 1)
    qc_method <- "DESeq2 (log2 norm)"
  }
  
  # plots
  save_plot_png(ma_png, {
    x <- NULL
    if ("baseMean" %in% colnames(res_df)) {
      x <- log10(res_df$baseMean + 1)
      xl <- "log10(baseMean + 1)"
    } else {
      x <- log10(rowMeans(nc)+1)
      xl <- "log10(mean normalized count + 1)"
    }
    plot(x, res_df$log2FoldChange, pch=16, cex=0.5,
         xlab=xl, ylab="log2FC", main="MA-like plot (DESeq2)")
    abline(h=c(-lfc_cutoff, lfc_cutoff), lty=2)
  })
  save_plot_png(vol_png, {
    plot(res_df$log2FoldChange, -log10(pmax(res_df$pvalue, 1e-300)), pch=16, cex=0.5,
         xlab="log2FC", ylab="-log10(p)", main="Volcano (DESeq2)")
    abline(v=c(-lfc_cutoff, lfc_cutoff), lty=2)
    abline(h=-log10(fdr_cutoff), lty=2)
  })
} else if (method == "limma" || method == "limmavoom" || method == "voom") {
  suppressWarnings(suppressMessages({
    if (!requireNamespace("limma", quietly=TRUE)) fail("limma not installed")
    if (!requireNamespace("edgeR", quietly=TRUE)) fail("edgeR not installed (needed for TMM+voom)")
    library(limma)
    library(edgeR)
  }))

  y <- DGEList(counts=counts_mat)
  y <- calcNormFactors(y)

  # Build design matrix
  design <- NULL
  if (advanced && !is.null(design_formula)) {
    design <- model.matrix(as.formula(design_formula), data=des_use)
  } else {
    if (is.null(cond)) fail("limma-voom: condition_A/B required in Simple mode (or provide design_formula)")
    design <- model.matrix(~cond)
  }

  # voom
  if (limma_quality_weights) {
    v <- voomWithQualityWeights(y, design, plot=FALSE)
  } else {
    v <- voom(y, design, plot=FALSE)
  }

  fit <- lmFit(v, design)

  if (!is.null(contrast_str)) {
    cont <- makeContrasts(contrasts=contrast_str, levels=colnames(design))
    fit <- contrasts.fit(fit, cont)
    fit <- eBayes(fit)
    if (limma_use_treat && limma_treat_lfc > 0) {
      fit <- treat(fit, lfc=limma_treat_lfc)
    }
    tt <- topTable(fit, coef=1, number=Inf, sort.by="P")
  } else {
    coef_idx <- NULL
    if (!is.null(params$coef) && !is.null(suppressWarnings(as.integer(params$coef)))) {
      coef_idx <- as.integer(params$coef)
    }
    # Step3: allow selecting coefficient by name
    if (is.null(coef_idx) && !is.null(coef_name) && is.character(coef_name) && nzchar(coef_name)) {
      cn <- trimws(coef_name)
      if (cn %in% colnames(design)) {
        coef_idx <- which(colnames(design) == cn)[1]
      }
    }
    if (is.null(coef_idx) && use_simple_cmp && !is.null(cond)) coef_idx <- 2L
    if (is.null(coef_idx) && use_simple_cmp && !is.null(cond_A) && condition_col %in% colnames(des_use)) {
      guess <- make.names(paste0(condition_col, cond_A))
      if (guess %in% colnames(design)) coef_idx <- which(colnames(design) == guess)[1]
    }
    if (is.null(coef_idx) | is.na(coef_idx)) coef_idx <- 2L
    fit <- eBayes(fit)
    if (limma_use_treat && limma_treat_lfc > 0) {
      fit <- treat(fit, lfc=limma_treat_lfc)
    }
    tt <- topTable(fit, coef=coef_idx, number=Inf, sort.by="P")
  }

  tt$gene <- rownames(tt)
  res_df <- tt[, c("gene","logFC","AveExpr","P.Value","adj.P.Val")]
  colnames(res_df) <- c("gene","log2FoldChange","logCPM","pvalue","padj")

  # normalized expression (voom logCPM)
  nc_df <- data.frame(gene=rownames(v$E), v$E, check.names=FALSE)
  write.table(nc_df, norm_counts_out, sep="	", quote=FALSE, row.names=FALSE)

  # QC expression matrix
  qc_expr <- v$E
  qc_method <- "limma-voom"

  save_plot_png(ma_png, {
    with(res_df, {
      plot(logCPM, log2FoldChange, pch=16, cex=0.5, xlab="AveExpr", ylab="log2FC", main="MA-like plot (limma-voom)")
      abline(h=c(-lfc_cutoff, lfc_cutoff), lty=2)
    })
  })
  save_plot_png(vol_png, {
    with(res_df, {
      plot(log2FoldChange, -log10(pmax(pvalue, 1e-300)), pch=16, cex=0.5,
           xlab="log2FC", ylab="-log10(p)", main="Volcano (limma-voom)")
      abline(v=c(-lfc_cutoff, lfc_cutoff), lty=2)
      abline(h=-log10(fdr_cutoff), lty=2)
    })
  })
} else if (method == "ebseq") {
  eb_fast <- as.logical(params$eb_fast)
  eb_maxround <- as.numeric(params$eb_maxround)
  
  sizes <- EBSeq::MedianNorm(Data = counts_mat,
                             alternative = F)
  eb <- EBSeq::EBTest(Data = counts_mat, 
                      NgVector = NULL,
                      Conditions = factor(des_use[, condition_col]), 
                      sizeFactors = sizes, 
                      fast = eb_fast,
                      Alpha = NULL,
                      Beta = NULL,
                      Qtrm = 1,
                      QtrmCut = 0,
                      maxround = eb_maxround,
                      step1 = 1e-06,
                      step2 = 0.01,
                      thre = log(2),
                      sthre = 0,
                      filter = 10,
                      stopthre = 1e-4)
  pp <- EBSeq::GetPPMat(EBout = eb)
  
  mA <- eb$Mean[, 1]
  mB <- eb$Mean[, 2]
  l2fc <- log2(mA/mB)
  names(l2fc) <- rownames(pp)
  l2fc2 <- l2fc[rownames(pp)]
  mean_exp <- eb$MeanList[rownames(pp)] + 1
  
  res_df <- data.frame(
    gene = rownames(pp),
    log2FoldChange = l2fc2,
    logCPM = log2(mean_exp),
    pvalue = NA_real_,
    padj = NA_real_,
    pp,
    stringsAsFactors = FALSE
  )
  # Convert posterior to a pseudo-q (best-effort): higher pp = more DE
  res_df$padj <- pmax(1 - res_df$PPDE, 0)
  res_df$pvalue <- res_df$padj
  # normalized counts placeholder
  #nc_df <- data.frame(gene=rownames(counts_mat), counts_mat, check.names=FALSE)
  #write.table(nc_df, norm_counts_out, sep="\t", quote=FALSE, row.names=FALSE)
  write.table(res_df, norm_counts_out, sep="\t", quote=FALSE, row.names=FALSE)
  
  save_plot_png(ma_png, {
    with(res_df, {
      plot(-log(padj), l2fc2, pch=16, cex=0.5, xlab="PPDE", ylab="log2FC", main="MA-like plot (EBSeq)")
      abline(h=c(-lfc_cutoff, lfc_cutoff), lty=2)
    })
  })
  
  save_plot_png(vol_png, {
    with(res_df, {
      plot(l2fc2, -log(padj), pch=16, cex=0.5,
           xlab="log2FC", ylab="PPDE", main="Volcano-like (EBSeq)")
      abline(v=c(-lfc_cutoff, lfc_cutoff), lty=2)
      abline(h=-log10(fdr_cutoff), lty=2)
    })
  })
# } else if (method == "samr" || method == "sam") {
#   ans <- samr::SAMseq(x = counts_mat,
#                       y = des_use[, condition_col],
#                       censoring.status = NULL,
#                       resp.type = "Two class unpaired",
#                       geneid = as.character(1:nrow(counts_mat)),
#                       genenames = rownames(counts_mat),
#                       nperms = 100,
#                       random.seed = NULL,
#                       nresamp = 20,
#                       fdr.output = 0.20)
#   
#   g_all <- data.frame(rownames(counts_mat),
#                       score = NA,
#                       fold_change = NA,
#                       q = NA)
#   rownames(g_all) <- rownames(counts_mat)
#   
#   g_down <- g_up <- NULL
#   if (!(is.null(ans$siggenes.table$genes.lo))) {
#     g_down <- ans$siggenes.table$genes.lo
#     rownames(g_down) <- g_down[, 2]
#     g_all[rownames(g_down), 2:4] <- g_down[, 3:5]
#   }
#   if (!(is.null(ans$siggenes.table$genes.up))) {
#     g_up <- ans$siggenes.table$genes.up
#     rownames(g_up) <- g_up[, 2]
#     g_all[rownames(g_up), 2:4] <- g_up[, 3:5]
#   }
#   
#   res_df <- data.frame(gene = rownames(counts_mat),
#                        g_all,
#                        padj = NA,
#                        pvalue = NA,
#                        log2FoldChange = NA)
#   write.table(res_df, norm_counts_out, sep="\t", quote=FALSE, row.names=FALSE)
  # if (is.null(res_df) || nrow(res_df)==0) {
  #   # fallback: no sig genes; still output ranking by d-score
  #   d <- sam_fit$tt[, "d"]
  #   res_df <- data.frame(
  #     gene=rownames(counts_mat),
  #     log2FoldChange=NA_real_,
  #     logCPM=NA_real_,
  #     pvalue=NA_real_,
  #     padj=NA_real_,
  #     d_score=as.numeric(d),
  #     stringsAsFactors=FALSE
  #   )
  # }
  #nc_df <- data.frame(gene=rownames(x), x, check.names=FALSE)
  #write.table(nc_df, norm_counts_out, sep="\t", quote=FALSE, row.names=FALSE)
  # save_plot_png(vol_png, {
  #   with(res_df, {
  #     plot(log2FoldChange, -log10(pmax(padj, 1e-300)), pch=16, cex=0.5,
  #          xlab="log2FC", ylab="-log10(q)", main="Volcano-like (samr)")
  #     abline(v=c(-lfc_cutoff, lfc_cutoff), lty=2)
  #     abline(h=-log10(fdr_cutoff), lty=2)
  #   })
  # })
  # save_plot_png(ma_png, {
  #   with(res_df, {
  #     plot(seq_along(log2FoldChange), log2FoldChange, pch=16, cex=0.5,
  #          xlab="rank", ylab="log2FC", main="Rank vs log2FC (samr)")
  #     abline(h=c(-lfc_cutoff, lfc_cutoff), lty=2)
  #   })
  # })
} else if (method %in% c("NOISeq", "noiseq", "noiseseq", "noi-seq")) {
  k <- as.numeric(params$k)
  norm <- as.character(params$norm)
  nclust <- as.numeric(params$nclust)
  lc <- as.numeric(params$lc)
  r <- as.numeric(params$r)
  adj <- as.numeric(params$adj)
  a0per <- as.numeric(params$a0per)
  filter <- as.numeric(params$filter)
  
  data_noi <- NOISeq::readData(data = counts_mat,
                               factors = des_use,
                               length = NULL,
                               biotype = NULL,
                               chromosome = NULL,
                               gc = NULL)
  ans <- NOISeq::noiseqbio(input = data_noi,
                           k = k,
                           norm = norm,
                           nclust = nclust,
                           plot = F,
                           factor = condition_col,
                           conditions = NULL,
                           lc = lc,
                           r = r,
                           adj = adj,
                           a0per = a0per,
                           random.seed = 12345,
                           filter = filter,
                           depth = NULL,
                           cv.cutoff = 500,
                           cpm = 1)
  # ans <- noiseq(input = data_noi,
  #               k = 0.5,
  #               norm = "rpkm",
  #               replicates = "biological",
  #               factor = "condition",
  #               conditions = NULL,
  #               pnr = 0.2,
  #               nss = 0.5,
  #               v = 0.02,
  #               lc = 0)
  
  res_df <- ans@results[[1]]
  mean_exp <- apply(res_df[, 1:2], 1, mean, na.rm = T)
  
  nc_df <- data.frame(gene=rownames(res_df), res_df)
  write.table(nc_df, norm_counts_out, sep="\t", quote=FALSE, row.names=FALSE)
  
  res_df <- data.frame(gene = rownames(res_df),
                       baseMean = mean_exp,
                       res_df, 
                       padj = -res_df$prob, 
                       pvalue = 10^(-abs(res_df$theta)), 
                       log2FoldChange = res_df$log2FC)
  res_df <- na.omit(res_df)
  
  save_plot_png(ma_png, {
    with(res_df, {
      plot(res_df$baseMean, res_df$log2FC, pch=16, cex=0.5, xlab="AveExpr", ylab="log2FC", main="MA-like plot")
      abline(h=c(-lfc_cutoff, lfc_cutoff), lty=2)
    })
  })
  save_plot_png(vol_png, {
    with(res_df, {
      plot(res_df$log2FC, res_df$theta, pch=16, cex=0.5,
           xlab="log2FC", ylab="theta", main="Volcano (NOISeq)")
      abline(v=c(-lfc_cutoff, lfc_cutoff), lty=2)
      abline(h=-log10(fdr_cutoff), lty=2)
    })
  })
# } else if (method == "bayseq") {
#   replicates <- des_use[, condition_col]
#   DE <- as.numeric(factor(replicates))
#   groups <- list(NDE = rep(1, length(replicates)),
#                  DE = DE)
#   
#   CD <- new("countData", 
#             data = count_mat, 
#             replicates = replicates, 
#             groups = groups)
#   libsizes(CD) <- getLibsizes(cD = CD)
#   
#   CD <- getPriors.NB(cD = CD, 
#                      samplesize = 10^5, 
#                      samplingSubset = NULL,
#                      equalDispersions = T,
#                      estimation = "QL", 
#                      verbose = T,
#                      zeroML = F,
#                      consensus = F,
#                      cl = NULL)
#   CD <- getLikelihoods.NB(cD = CD, 
#                           pET = "BIC", 
#                           marginalise = FALSE, 
#                           subset = NULL,
#                           priorSubset = NULL, 
#                           bootStraps = 1, 
#                           conv = 1e-4, 
#                           nullData = FALSE,
#                           returnAll = FALSE, 
#                           returnPD = FALSE, 
#                           verbose = TRUE, 
#                           discardSampling = FALSE, 
#                           cl = NULL)
#   
#   mA <- rowMeans(counts_mat[, cond == cond_A, drop=FALSE] + 0.5)
#   mB <- rowMeans(counts_mat[, cond == cond_B, drop=FALSE] + 0.5)
#   l2fc <- log2(mA/mB)
#   res_df <- data.frame(
#     gene=rownames(counts_mat),
#     log2FoldChange=l2fc,
#     logCPM=log2(rowMeans(counts_mat)+1),
#     pvalue=NA_real_,
#     padj=NA_real_,
#     stringsAsFactors=FALSE
#   )
#   nc_df <- data.frame(gene=rownames(counts_mat), counts_mat, check.names=FALSE)
#   write.table(nc_df, norm_counts_out, sep="\t", quote=FALSE, row.names=FALSE)
#   save_plot_png(ma_png, {
#     with(res_df, {
#       plot(logCPM, log2FoldChange, pch=16, cex=0.5, xlab="log2(mean count + 1)", ylab="log2FC", main="MA-like plot (baySeq placeholder)")
#       abline(h=c(-lfc_cutoff, lfc_cutoff), lty=2)
#     })
#   })
#   save_plot_png(vol_png, {
#     with(res_df, {
#       plot(log2FoldChange, rep(0, length(log2FoldChange)), pch=16, cex=0.5,
#            xlab="log2FC", ylab="(no p-values)", main="Volcano placeholder (baySeq)")
#       abline(v=c(-lfc_cutoff, lfc_cutoff), lty=2)
#     })
#   })
} else {
  fail(paste0("Unknown method: ", method))
}

# write DEG table
res_df <- res_df[order(res_df$padj, res_df$pvalue), , drop=FALSE]
write.table(res_df, deg_out, sep="	", quote=FALSE, row.names=FALSE)

# QC plots (PCA / distance / library size)
qc_lib_png <- file.path(out_dir, "qc_libsize.png")
qc_pca_png <- file.path(out_dir, "qc_pca.png")
qc_dist_png <- file.path(out_dir, "qc_sample_distance.png")
if (!is.null(qc_expr)) {
  try(make_qc_plots(out_dir, counts_mat, qc_expr, des_use, condition_col, qc_method %||% method), silent=TRUE)
}

# Backward-compatible plot filenames expected by older GUI code
try({
  if (file.exists(vol_png)) file.copy(vol_png, file.path(out_dir, "volcano.png"), overwrite=TRUE)
}, silent=TRUE)
try({
  if (file.exists(ma_png)) file.copy(ma_png, file.path(out_dir, "ma_plot.png"), overwrite=TRUE)
}, silent=TRUE)

# select DE genes
sel <- res_df[!is.na(res_df$padj) & res_df$padj <= fdr_cutoff & !is.na(res_df$log2FoldChange) & abs(res_df$log2FoldChange) >= lfc_cutoff, , drop=FALSE]
deg_list_out <- file.path(out_dir, "deg_list.txt")
writeLines(sel$gene, deg_list_out)

go_out <- NULL
if (do_go) {
  suppressWarnings(suppressMessages({
    if (!requireNamespace("GOstats", quietly=TRUE)) fail("GOstats not installed")
    if (!requireNamespace("AnnotationDbi", quietly=TRUE)) fail("AnnotationDbi not installed")
    if (!requireNamespace("GO.db", quietly=TRUE)) fail("GO.db not installed")
    library(GOstats)
    library(AnnotationDbi)
    library(GO.db)
  }))
  # Load orgdb package by name
  if (!requireNamespace(orgdb_pkg, quietly=TRUE)) {
    fail(paste0("OrgDb package not installed: ", orgdb_pkg))
  }
  suppressWarnings(suppressMessages(library(orgdb_pkg, character.only=TRUE)))
  orgdb <- get(orgdb_pkg)

  genes_all <- unique(res_df$gene)
  genes_sel <- unique(sel$gene)

  # Map to ENTREZID if needed
  to_entrez <- function(keys) {
    keys <- as.character(keys)
    if (toupper(keytype) == "ENTREZID") return(unique(keys))
    m <- AnnotationDbi::select(orgdb, keys=keys, keytype=keytype, columns=c("ENTREZID"))
    m <- m[!is.na(m$ENTREZID), , drop=FALSE]
    unique(as.character(m$ENTREZID))
  }
  universe_entrez <- to_entrez(genes_all)
  selected_entrez <- to_entrez(genes_sel)

  if (length(selected_entrez) < 5) fail("Too few selected genes for GOstats (after mapping)")
  if (length(universe_entrez) < 50) fail("Too small gene universe for GOstats (after mapping)")

  go_out <- file.path(out_dir, sprintf("go_enrichment_%s.tsv", ontology))

  params_go <- new("GOHyperGParams",
                   geneIds=selected_entrez,
                   universeGeneIds=universe_entrez,
                   annotation=orgdb_pkg,
                   ontology=ontology,
                   pvalueCutoff=1.0,
                   conditional=FALSE,
                   testDirection="over")
  hg <- hyperGTest(params_go)
  sm <- summary(hg)
  write.table(sm, go_out, sep="\t", quote=FALSE, row.names=FALSE)
}

# artifacts
art <- list(
  deg_tsv = deg_out,
  deg_list = deg_list_out,
  ma_plot = ma_png,
  volcano_plot = vol_png,
  normalized_counts = norm_counts_out,
  qc_libsize = qc_lib_png,
  qc_pca = qc_pca_png,
  qc_sample_distance = qc_dist_png,
  go_tsv = go_out
)
writeLines(jsonlite::toJSON(art, auto_unbox=TRUE, pretty=TRUE), file.path(out_dir, "artifacts.json"))

cat("[rnaseq_deg] done\n")
