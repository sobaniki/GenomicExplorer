#!/usr/bin/env Rscript
# rnaseq_cluster: transform -> HVG -> PCA -> UMAP/t-SNE -> clustering

suppressWarnings(suppressMessages({
  library(jsonlite)
  library(data.table)
  library(mclust)
}))

fail <- function(msg, code=1) {
  cat("[rnaseq_cluster][ERROR] ", msg, "\n", sep="")
  quit(status=code)
}

`%||%` <- function(a, b) {
  if (is.null(a) || length(a)==0 || (is.character(a) && a=="")) b else a
}

as_bool <- function(x, default=FALSE) {
  if (is.null(x)) return(default)
  if (is.logical(x)) return(isTRUE(x))
  if (is.numeric(x)) return(x != 0)
  if (is.character(x)) return(tolower(x) %in% c("1","true","t","yes","y"))
  default
}

read_table_auto <- function(path) {
  tryCatch({
    data.table::fread(path, data.table=FALSE, check.names=FALSE)
  }, error=function(e) {
    read.csv(path, check.names=FALSE, stringsAsFactors=FALSE)
  })
}

save_plot_png <- function(path, expr) {
  png(path, width=1000, height=700, res=120)
  on.exit(dev.off(), add=TRUE)
  expr
}

require_or <- function(pkg, alt_msg=NULL) {
  ok <- suppressWarnings(requireNamespace(pkg, quietly=TRUE))
  if (!ok) {
    if (is.null(alt_msg)) alt_msg <- paste0("R package '", pkg, "' is required but not installed.")
    fail(alt_msg)
  }
  TRUE
}

row_vars <- function(mat) {
  if (suppressWarnings(requireNamespace("matrixStats", quietly=TRUE))) {
    return(matrixStats::rowVars(mat, na.rm=TRUE))
  }
  apply(mat, 1, stats::var, na.rm=TRUE)
}

# ---------- args ----------
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 4 || args[1] != "--params" || args[3] != "--out") {
  cat("Usage: runner.R --params params.json --out out_dir\n")
  quit(status=2)
}
params_path <- args[2]
out_dir <- args[4]
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

cat("[rnaseq_cluster] start\n")
cat("[rnaseq_cluster] params_path=", params_path, "\n", sep="")
cat("[rnaseq_cluster] out_dir=", out_dir, "\n", sep="")

params <- tryCatch(jsonlite::fromJSON(params_path), error=function(e) fail(paste("Failed to parse params.json:", e$message)))

# Inputs
counts_path <- params$counts_path %||% params$counts_matrix %||% params$counts
design_path <- params$design_path %||% params$design_table %||% params$design
sample_col <- params$sample_col %||% "sample"
group_col  <- params$group_col  %||% params$condition_col %||% "condition"

if (is.null(counts_path) || !file.exists(counts_path)) fail(paste0("Counts matrix not found: ", counts_path))
if (is.null(design_path) || !file.exists(design_path)) fail(paste0("Design table not found: ", design_path))

# Options
transform <- tolower(params$transform %||% "vst")   # vst / rlog / logcpm / log2norm
n_hvg <- as.integer(params$n_hvg %||% params$hvg %||% 2000)
min_cpm <- as.numeric(params$min_cpm %||% 0)
min_samples <- as.integer(params$min_samples %||% 0)
scale_genes <- as_bool(params$scale_genes %||% params$zscore %||% FALSE)

pca_method <- tolower(params$pca_method %||% "prcomp")  # prcomp / irlba
n_pcs <- as.integer(params$n_pcs %||% params$pca_npcs %||% 50)
cluster_input_pcs <- as.integer(params$cluster_input_pcs %||% min(20, n_pcs))
umap_input_pcs <- as.integer(params$umap_input_pcs %||% min(30, n_pcs))
tsne_input_pcs <- as.integer(params$tsne_input_pcs %||% min(30, n_pcs))

do_umap <- as_bool(params$do_umap %||% TRUE)
do_tsne <- as_bool(params$do_tsne %||% FALSE)

seed <- as.integer(params$seed %||% 1)

# UMAP params
umap_n_neighbors <- as.integer(params$umap_n_neighbors %||% params$n_neighbors %||% 15)
umap_min_dist <- as.numeric(params$umap_min_dist %||% params$min_dist %||% 0.1)
umap_metric <- params$umap_metric %||% "euclidean"

# tSNE params
tsne_perplexity <- as.numeric(params$tsne_perplexity %||% 30)
tsne_theta <- as.numeric(params$tsne_theta %||% 0.5)

# Clustering params
cluster_method <- tolower(params$cluster_method %||% params$cluster %||% "kmeans") # none/kmeans/hclust/mclust
k <- as.integer(params$k %||% params$n_clusters %||% 3)
hclust_method <- params$hclust_method %||% "ward.D2"

set.seed(seed)

# ---------- read counts + design ----------
cnt_df <- read_table_auto(counts_path)
if (ncol(cnt_df) < 3) fail("Counts matrix must have at least: gene_id + 2 samples")
genes <- as.character(cnt_df[[1]])
cnt_df <- cnt_df[,-1, drop=FALSE]
cnt <- as.matrix(cnt_df)
mode(cnt) <- "numeric"
rownames(cnt) <- make.unique(genes)

des <- read_table_auto(design_path)
if (!(sample_col %in% colnames(des))) fail(paste0("design table must contain sample_col: ", sample_col))

# align samples
sample_names <- colnames(cnt)
keep_samples <- intersect(sample_names, as.character(des[[sample_col]]))
if (length(keep_samples) < 2) fail("No overlapping samples between counts columns and design sample IDs")

des_use <- des[match(keep_samples, as.character(des[[sample_col]])), , drop=FALSE]
cnt_use <- cnt[, keep_samples, drop=FALSE]

# coerce group col (if present)
grp <- NULL
if (!is.null(group_col) && group_col %in% colnames(des_use)) {
  grp <- as.factor(des_use[[group_col]])
  names(grp) <- as.character(des_use[[sample_col]])
}

# optional low-expression filter (CPM-based)
if (min_cpm > 0 && min_samples > 0) {
  if (suppressWarnings(requireNamespace("edgeR", quietly=TRUE))) {
    cpm0 <- edgeR::cpm(cnt_use)
    keep_gene <- rowSums(cpm0 >= min_cpm) >= min_samples
    cnt_use <- cnt_use[keep_gene, , drop=FALSE]
  } else {
    cat("[rnaseq_cluster] edgeR not installed; skipping CPM filter\n")
  }
}
if (nrow(cnt_use) < 2) fail("Too few genes after filtering")

# ---------- transform ----------
expr <- NULL
transform_used <- transform

if (transform %in% c("vst","rlog","log2norm")) {
  if (!suppressWarnings(requireNamespace("DESeq2", quietly=TRUE))) {
    cat("[rnaseq_cluster] DESeq2 not installed; falling back to logCPM transform\n")
    transform_used <- "logcpm"
  }
}

if (transform_used == "logcpm") {
  require_or("edgeR", "edgeR is required for logCPM transform (fallback). Please install edgeR or DESeq2.")
  y <- edgeR::DGEList(counts=cnt_use)
  y <- edgeR::calcNormFactors(y)
  expr <- edgeR::cpm(y, log=TRUE, prior.count=1)
} else if (transform_used == "log2norm") {
  require_or("DESeq2")
  coldata <- data.frame(row.names=colnames(cnt_use), dummy=factor(rep("all", ncol(cnt_use))))
  dds <- DESeq2::DESeqDataSetFromMatrix(countData=round(cnt_use), colData=coldata, design=~1)
  dds <- DESeq2::estimateSizeFactors(dds)
  nc <- DESeq2::counts(dds, normalized=TRUE)
  expr <- log2(nc + 1)
} else if (transform_used == "rlog") {
  require_or("DESeq2")
  coldata <- data.frame(row.names=colnames(cnt_use), dummy=factor(rep("all", ncol(cnt_use))))
  dds <- DESeq2::DESeqDataSetFromMatrix(countData=round(cnt_use), colData=coldata, design=~1)
  dds <- DESeq2::estimateSizeFactors(dds)
  rld <- DESeq2::rlog(dds, blind=TRUE)
  expr <- SummarizedExperiment::assay(rld)
} else { # vst
  require_or("DESeq2")
  coldata <- data.frame(row.names=colnames(cnt_use), dummy=factor(rep("all", ncol(cnt_use))))
  dds <- DESeq2::DESeqDataSetFromMatrix(countData=round(cnt_use), colData=coldata, design=~1)
  dds <- DESeq2::estimateSizeFactors(dds)
  vsd <- DESeq2::varianceStabilizingTransformation(dds, blind=TRUE)
  expr <- SummarizedExperiment::assay(vsd)
}

if (is.null(expr) || ncol(expr) < 2) fail("Failed to build transformed expression matrix")

# ---------- HVG ----------
vars <- row_vars(expr)
ord <- order(vars, decreasing=TRUE)
if (is.na(n_hvg) || n_hvg <= 0 || n_hvg >= nrow(expr)) {
  keep_idx <- ord
} else {
  keep_idx <- ord[seq_len(n_hvg)]
}
expr_hvg <- expr[keep_idx, , drop=FALSE]
vars_hvg <- vars[keep_idx]
gene_hvg <- rownames(expr_hvg)

hvg_tbl <- data.frame(gene=gene_hvg, variance=vars_hvg, stringsAsFactors=FALSE)
data.table::fwrite(hvg_tbl, file.path(out_dir, "hvg.tsv"), sep="\t")

# optional z-score per gene
if (scale_genes) {
  expr_hvg <- t(scale(t(expr_hvg)))
  expr_hvg[is.na(expr_hvg)] <- 0
}

# ---------- PCA ----------
X <- t(expr_hvg) # samples x genes

if (n_pcs < 2) n_pcs <- 2
if (cluster_input_pcs < 2) cluster_input_pcs <- 2
if (umap_input_pcs < 2) umap_input_pcs <- 2
if (tsne_input_pcs < 2) tsne_input_pcs <- 2

pca <- NULL
if (pca_method == "irlba" && suppressWarnings(requireNamespace("irlba", quietly=TRUE))) {
  pca <- irlba::prcomp_irlba(X, n=n_pcs, center=TRUE, scale.=TRUE)
} else {
  pca_method <- "prcomp"
  pca <- prcomp(X, center=TRUE, scale.=TRUE, rank.=n_pcs)
}
scores <- as.data.frame(pca$x)
scores$sample <- rownames(scores)

# explained variance
expl <- (pca$sdev^2) / sum(pca$sdev^2)
expl_tbl <- data.frame(PC=paste0("PC", seq_along(expl)), variance_explained=expl, stringsAsFactors=FALSE)
data.table::fwrite(expl_tbl, file.path(out_dir, "pca_variance.tsv"), sep="\t")

# ---------- embeddings ----------
emb <- data.frame(sample=colnames(expr_hvg), stringsAsFactors=FALSE)
if (!is.null(grp)) emb[[group_col]] <- as.character(grp[colnames(expr_hvg)])
for (i in seq_len(min(n_pcs, ncol(pca$x)))) {
  emb[[paste0("PC", i)]] <- pca$x[, i]
}

umap_xy <- NULL
tsne_xy <- NULL

if (do_umap) {
  require_or("uwot", "UMAP requires the 'uwot' package. Install it or disable UMAP.")
  use_pcs <- min(umap_input_pcs, ncol(pca$x))
  um <- uwot::umap(pca$x[, seq_len(use_pcs), drop=FALSE],
                   n_neighbors=umap_n_neighbors,
                   min_dist=umap_min_dist,
                   metric=umap_metric,
                   n_components=2,
                   verbose=FALSE)
  colnames(um) <- c("UMAP1","UMAP2")
  umap_xy <- um
  emb$UMAP1 <- um[,1]
  emb$UMAP2 <- um[,2]
}

if (do_tsne) {
  require_or("Rtsne", "t-SNE requires the 'Rtsne' package. Install it or disable t-SNE.")
  use_pcs <- min(tsne_input_pcs, ncol(pca$x))
  perp <- tsne_perplexity
  max_perp <- floor((nrow(pca$x) - 1) / 3)
  if (is.finite(max_perp) && max_perp >= 5) perp <- min(perp, max_perp)
  ts <- Rtsne::Rtsne(pca$x[, seq_len(use_pcs), drop=FALSE],
                     dims=2, perplexity=perp, theta=tsne_theta,
                     verbose=FALSE, check_duplicates=FALSE)
  tsne_xy <- ts$Y
  colnames(tsne_xy) <- c("TSNE1","TSNE2")
  emb$TSNE1 <- tsne_xy[,1]
  emb$TSNE2 <- tsne_xy[,2]
}

data.table::fwrite(emb, file.path(out_dir, "sample_embedding.tsv"), sep="\t")

# ---------- clustering ----------
clusters <- NULL
cluster_note <- ""

if (cluster_method %in% c("none","")) {
  cluster_method <- "none"
} else if (cluster_method == "kmeans") {
  use_pcs <- min(cluster_input_pcs, ncol(pca$x))
  km <- kmeans(pca$x[, seq_len(use_pcs), drop=FALSE], centers=k, nstart=10)
  clusters <- km$cluster
} else if (cluster_method == "hclust") {
  use_pcs <- min(cluster_input_pcs, ncol(pca$x))
  d <- dist(pca$x[, seq_len(use_pcs), drop=FALSE])
  hc <- hclust(d, method=hclust_method)
  clusters <- cutree(hc, k=k)
  save_plot_png(file.path(out_dir, "dendrogram.png"), {
    plot(hc, main=paste0("Hierarchical clustering (k=", k, ")"))
    rect.hclust(hc, k=k, border="red")
  })
} else if (cluster_method == "mclust") {
  require_or("mclust", "Model-based clustering requires the 'mclust' package.")
  use_pcs <- min(cluster_input_pcs, ncol(pca$x))
  if (is.finite(k) && k >= 1) {
    mc <- mclust::Mclust(pca$x[, seq_len(use_pcs), drop=FALSE], G=k)
    cluster_note <- paste0("mclust model=", mc$modelName)
  } else {
    mc <- mclust::Mclust(pca$x[, seq_len(use_pcs), drop=FALSE])
    cluster_note <- paste0("mclust chose G=", mc$G, " model=", mc$modelName)
    k <- mc$G
  }
  clusters <- mc$classification
} else {
  fail(paste0("Unknown cluster_method: ", cluster_method))
}

# ---------- plots ----------
make_scatter <- function(x, y, xlab, ylab, title, color_by=NULL, out_png) {
  save_plot_png(out_png, {
    colv <- "black"
    legend_items <- NULL
    if (!is.null(color_by)) {
      f <- as.factor(color_by)
      pal <- grDevices::rainbow(length(levels(f)))
      colv <- pal[as.integer(f)]
      legend_items <- list(levels=levels(f), cols=pal)
    }
    plot(x, y, pch=19, col=colv, xlab=xlab, ylab=ylab, main=title)
    text(x, y, labels=names(x), pos=3, cex=0.6)
    if (!is.null(legend_items)) {
      legend("topright", legend=legend_items$levels, col=legend_items$cols, pch=19, bty="n", cex=0.8)
    }
  })
}

color_vec <- NULL
if (!is.null(clusters)) {
  names(clusters) <- colnames(expr_hvg)
  color_vec <- as.factor(clusters[colnames(expr_hvg)])
} else if (!is.null(grp)) {
  color_vec <- grp[colnames(expr_hvg)]
}

make_scatter(pca$x[,1], pca$x[,2],
             "PC1", "PC2",
             paste0("PCA (", transform_used, ", HVG=", nrow(expr_hvg), ", pca=", pca_method, ")"),
             color_by=color_vec,
             out_png=file.path(out_dir, "pca.png"))

if (!is.null(umap_xy)) {
  make_scatter(umap_xy[,1], umap_xy[,2],
               "UMAP1", "UMAP2",
               paste0("UMAP (n_neighbors=", umap_n_neighbors, ", min_dist=", umap_min_dist, ")"),
               color_by=color_vec,
               out_png=file.path(out_dir, "umap.png"))
}

if (!is.null(tsne_xy)) {
  make_scatter(tsne_xy[,1], tsne_xy[,2],
               "tSNE1", "tSNE2",
               paste0("t-SNE (perplexity=", tsne_perplexity, ", theta=", tsne_theta, ")"),
               color_by=color_vec,
               out_png=file.path(out_dir, "tsne.png"))
}

if (!is.null(clusters)) {
  cl_tbl <- data.frame(sample=colnames(expr_hvg),
                       cluster=as.integer(clusters[colnames(expr_hvg)]),
                       method=cluster_method,
                       note=cluster_note,
                       stringsAsFactors=FALSE)
  data.table::fwrite(cl_tbl, file.path(out_dir, "sample_clusters.tsv"), sep="\t")

  sz <- as.data.frame(table(cluster=as.integer(clusters)))
  colnames(sz) <- c("cluster","n")
  data.table::fwrite(sz, file.path(out_dir, "cluster_sizes.tsv"), sep="\t")

  save_plot_png(file.path(out_dir, "cluster_sizes.png"), {
    barplot(sz$n, names.arg=sz$cluster, las=1, xlab="cluster", ylab="n",
            main=paste0("Cluster sizes (", cluster_method, ")"))
  })
} else {
  data.table::fwrite(data.frame(sample=colnames(expr_hvg), cluster=NA_integer_, method="none"),
                     file.path(out_dir, "sample_clusters.tsv"), sep="\t")
}

art <- list(
  plugin="rnaseq_cluster",
  transform=transform_used,
  n_hvg=nrow(expr_hvg),
  n_samples=ncol(expr_hvg),
  pca_method=pca_method,
  n_pcs=n_pcs,
  do_umap=do_umap,
  do_tsne=do_tsne,
  cluster_method=cluster_method,
  k=k,
  files=list(
    sample_embedding="sample_embedding.tsv",
    sample_clusters="sample_clusters.tsv",
    hvg="hvg.tsv",
    pca_variance="pca_variance.tsv",
    pca_plot="pca.png",
    umap_plot=if (!is.null(umap_xy)) "umap.png" else NULL,
    tsne_plot=if (!is.null(tsne_xy)) "tsne.png" else NULL,
    dendrogram=if (file.exists(file.path(out_dir, "dendrogram.png"))) "dendrogram.png" else NULL,
    cluster_sizes_plot=if (file.exists(file.path(out_dir, "cluster_sizes.png"))) "cluster_sizes.png" else NULL
  )
)
writeLines(jsonlite::toJSON(art, auto_unbox=TRUE, pretty=TRUE), file.path(out_dir, "artifacts.json"))

cat("[rnaseq_cluster] done\n")
