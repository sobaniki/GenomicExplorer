#!/usr/bin/env Rscript
# rnaseq_genecluster_modelbased: gene-level model-based clustering
# backends: HTSCluster / MBCluster.Seq / coseq / clusterSeq

suppressWarnings(suppressMessages({
  library(jsonlite)
  library(data.table)
}))

fail <- function(msg, code=1) {
  cat("[rnaseq_genecluster_modelbased][ERROR] ", msg, "\n", sep="")
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

require_or <- function(pkg, alt_msg=NULL) {
  ok <- suppressWarnings(requireNamespace(pkg, quietly=TRUE))
  if (!ok) {
    if (is.null(alt_msg)) alt_msg <- paste0("R package '", pkg, "' is required but not installed.")
    fail(alt_msg)
  }
  TRUE
}

write_session_info <- function(path) {
  tryCatch({
    si <- utils::capture.output(sessionInfo())
    writeLines(si, con=path)
  }, error=function(e) {
    # ignore
  })
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

cat("[rnaseq_genecluster_modelbased] start\n")
cat("[rnaseq_genecluster_modelbased] params_path=", params_path, "\n", sep="")
cat("[rnaseq_genecluster_modelbased] out_dir=", out_dir, "\n", sep="")

params <- tryCatch(jsonlite::fromJSON(params_path), error=function(e) fail(paste("Failed to parse params.json:", e$message)))

# Inputs
counts_path <- params$counts_path %||% params$counts_matrix %||% params$counts
design_path <- params$design_path %||% params$design_table %||% params$design
sample_col <- params$sample_col %||% "sample"
group_col  <- params$group_col  %||% params$condition_col %||% "condition"

if (is.null(counts_path) || !file.exists(counts_path)) fail(paste0("Counts matrix not found: ", counts_path))
if (is.null(design_path) || !file.exists(design_path)) fail(paste0("Design table not found: ", design_path))

# Options
backend <- tolower(params$backend %||% params$method %||% "htscluster")  # htscluster / mbclusterseq / coseq / clusterseq
seed <- as.integer(params$seed %||% 1)

# Gene filter (recommended)
n_hvg <- as.integer(params$n_hvg %||% params$hvg %||% 2000)
min_cpm <- as.numeric(params$min_cpm %||% 0)
min_samples <- as.integer(params$min_samples %||% 0)

# HTSCluster options
gmin <- as.integer(params$gmin %||% 2)
gmax <- as.integer(params$gmax %||% params$g %||% params$k %||% 8)
choose_by <- toupper(params$choose_by %||% "ICL")  # ICL or BIC
force_g <- as.integer(params$force_g %||% 0)       # if >0, use this g
norm_method <- params$norm_method %||% "TMM"       # TMM/TC/UQ/Med/DESeq or numeric vector
alg_type <- toupper(params$alg_type %||% "EM")     # EM/CEM
init_type <- params$init_type %||% "small-em"      # small-em/kmeans
iter <- as.integer(params$iter %||% 1000)
cutoff <- as.numeric(params$cutoff %||% 1e-5)

# MBCluster.Seq options (if installed)
mb_model <- tolower(params$mb_model %||% "nbinom")     # nb / poisson
mb_k <- as.integer(params$mb_k %||% params$k %||% 6)
mb_iter <- as.integer(params$mb_iter %||% 100)

# coseq options (if installed)
coseq_model <- params$coseq_model %||% params$model %||% "Poisson"
coseq_kmin <- as.integer(params$coseq_kmin %||% params$kmin %||% 2)
coseq_kmax <- as.integer(params$coseq_kmax %||% params$kmax %||% 8)
coseq_transform <- params$coseq_transform %||% "log"

# clusterSeq options (if installed)
clusterseq_k <- as.integer(params$clusterseq_k %||% params$k %||% 6)

set.seed(seed)

# ---------- read counts + design ----------
cnt_df <- read_table_auto(counts_path)
if (ncol(cnt_df) < 3) fail("Counts matrix must have at least: gene_id + 2 samples")
genes <- as.character(cnt_df[[1]])
cnt_df <- cnt_df[,-1, drop=FALSE]
cnt <- as.matrix(cnt_df)
mode(cnt) <- "numeric"
colnames(cnt) <- colnames(cnt_df)

des <- read_table_auto(design_path)
if (!(sample_col %in% colnames(des))) fail(paste0("Design table must contain column: ", sample_col))
if (!(group_col %in% colnames(des))) fail(paste0("Design table must contain column: ", group_col))

# Align samples by sample_col order
des[[sample_col]] <- as.character(des[[sample_col]])
sample_names <- colnames(cnt)
idx <- match(sample_names, des[[sample_col]])
if (any(is.na(idx))) {
  missing <- sample_names[is.na(idx)]
  fail(paste0("Some samples in counts are missing in design table (", sample_col, "): ", paste(missing, collapse=", ")))
}
des2 <- des[idx, , drop=FALSE]

conds <- as.factor(des2[[group_col]])
if (nlevels(conds) < 2) fail("Need at least 2 groups/conditions for gene clustering.")

# Optional low-expression filter using CPM-like
libsize <- colSums(cnt)
cpm <- sweep(cnt, 2, libsize/1e6, "/")
keep <- rep(TRUE, nrow(cnt))
if (min_samples > 0 && min_cpm > 0) {
  keep <- rowSums(cpm >= min_cpm, na.rm=TRUE) >= min_samples
}

# HVG selection (variance on logCPM)
logcpm <- log2(cpm + 1)
if (!is.null(n_hvg) && n_hvg > 0 && n_hvg < nrow(logcpm)) {
  v <- row_vars(logcpm[keep, , drop=FALSE])
  ord <- order(v, decreasing=TRUE)
  top_idx <- which(keep)[ord[seq_len(min(n_hvg, length(ord)))]]
  keep2 <- rep(FALSE, nrow(cnt)); keep2[top_idx] <- TRUE
  keep <- keep2
}

cnt_f <- cnt[keep, , drop=FALSE]
genes_f <- genes[keep]

# Write minimal reproducibility artifacts early
write_session_info(file.path(out_dir, "sessionInfo.txt"))

# ---------- run backend clustering ----------
artifacts <- list(
  plugin="rnaseq_genecluster_modelbased",
  backend=backend,
  seed=seed,
  inputs=list(
    counts_path=counts_path,
    design_path=design_path,
    sample_col=sample_col,
    group_col=group_col,
    samples=as.character(sample_names)
  ),
  n_genes_input=nrow(cnt),
  n_genes_used=nrow(cnt_f),
  n_samples=ncol(cnt_f),
  group_col=group_col,
  groups=levels(conds),
  filter=list(n_hvg=n_hvg, min_cpm=min_cpm, min_samples=min_samples),
  selected_g=NA,
  params=list(
    backend=backend,
    seed=seed,
    n_hvg=n_hvg,
    min_cpm=min_cpm,
    min_samples=min_samples,
    # backend-specific knobs (record even if unused)
    htscluster=list(gmin=gmin, gmax=gmax, choose_by=choose_by, force_g=force_g, norm_method=norm_method, alg_type=alg_type, init_type=init_type, iter=iter, cutoff=cutoff),
    mbclusterseq=list(model=mb_model, k=mb_k, iter=mb_iter),
    coseq=list(model=coseq_model, kmin=coseq_kmin, kmax=coseq_kmax, transform=coseq_transform),
    clusterseq=list(k=clusterseq_k)
  )
)

cluster_labels <- NULL
selected_g <- NA
model_scores <- NULL

if (backend %in% c("htscluster","hts")) {
  require_or("HTSCluster", "R package 'HTSCluster' is required for backend=htscluster. Install from CRAN: install.packages('HTSCluster').")
  y <- cnt_f
  storage.mode(y) <- "numeric"

  # Decide whether wrapper or single g
  if (!is.null(force_g) && force_g > 0) {
    gmin2 <- force_g; gmax2 <- force_g
  } else {
    gmin2 <- gmin
    gmax2 <- gmax
    if (gmax2 < gmin2) { tmp <- gmin2; gmin2 <- gmax2; gmax2 <- tmp }
  }

  res <- tryCatch({
    if (gmin2 == gmax2) {
      HTSCluster::PoisMixClus(y, g=gmin2, conds=as.vector(conds), norm=norm_method,
                              init.type=init_type, alg.type=alg_type,
                              cutoff=cutoff, iter=iter, verbose=FALSE)
    } else {
      HTSCluster::PoisMixClusWrapper(y, gmin=gmin2, gmax=gmax2, conds=as.vector(conds), norm=norm_method,
                                     gmin.init.type=init_type, init.runs=1, init.iter=10,
                                     split.init=TRUE, alg.type=alg_type, cutoff=cutoff, iter=iter,
                                     verbose=FALSE)
    }
  }, error=function(e) fail(paste0("HTSCluster failed: ", e$message)))

  if (inherits(res, "HTSClusterWrapper")) {
    lol <- res$logLike.all
    icl <- res$ICL.all
    g_candidates <- seq(gmin2, gmax2)
    model_scores <- data.frame(g=g_candidates, LL=lol, ICL=icl)
    if (choose_by == "ICL") {
      selected_g <- model_scores$g[which.max(model_scores$ICL)]
      cluster_labels <- res$ICL.results$labels
    } else {
      selected_g <- model_scores$g[which.max(model_scores$LL)]
      cluster_labels <- res$BIC.results$labels
    }
  } else if (inherits(res, "HTSCluster")) {
    selected_g <- gmin2
    cluster_labels <- res$labels
    pp <- res$probaPost
  }

  artifacts$backend_details <- list(
    gmin=gmin2, gmax=gmax2, selected_g=selected_g,
    choose_by=choose_by, norm=norm_method, alg=alg_type, init=init_type,
    iter=iter, cutoff=cutoff
  )

} else if (backend %in% c("mbclusterseq","mbcluster.seq","mbcluster")) {
  Count <- cnt_f
  Treatment <- as.character(conds)
  GeneID <- genes_f

  obj <- MBCluster.Seq::RNASeq.Data(Count = cnt_f, 
                                    Normalizer = NULL, 
                                    Treatment = conds, 
                                    GeneID = genes_f)
  nk <- MBCluster.Seq::KmeansPlus.RNASeq(data = obj,
                                         nK = 10,
                                         model = mb_model,
                                         print.steps = F)
  res <- MBCluster.Seq::Cluster.RNASeq(data = obj, 
                                       model = mb_model,
                                       centers = nk$centers,
                                       #centers = NULL,
                                       method = "EM",
                                       iter.max = mb_iter,
                                       TMP = NULL)
  # trr <- MBCluster.Seq::Hybrid.Tree(data = obj,
  #                                   cluster0 = res,
  #                                   model = mb_model)
  cluster_labels <- res$cluster
  selected_g <- max(cluster_labels)
  artifacts$backend_details <- list(engine="MBCluster.Seq", model=mb_model, k=mb_k, iter=mb_iter)
} else if (backend %in% c("coseq")) {
  if (coseq_model == "Poisson") {
    warning("Only Transformation: none can be used with Model: Poisson")
    coseq_transform = "none"
  }
  res <- coseq::coseq(object = cnt_f, 
                      K = seq(coseq_kmin, 
                              coseq_kmax),
                      subset = NULL,
                      model = coseq_model, 
                      transformation = coseq_transform, 
                      normFactors = "TMM",
                      meanFilterCutoff = NULL,
                      modelChoice = "ICL",
                      parallel = F,
                      #BPPARAM = bpparam(),
                      seed = seed)
  cluster_labels <- coseq::clusters(res)
  selected_g <- max(cluster_labels)
  artifacts$backend_details <- list(engine="coseq", model=coseq_model, kmin=coseq_kmin, kmax=coseq_kmax, transform=coseq_transform, selected_k=selected_g)

} else if (backend %in% c("clusterseq","clusterseq")) {
  library(clusterSeq)
  # CD <- new("countData", 
  #           data = cnt_f, 
  #           replicates = des2[[group_col]])
  libsizes <- getLibsizes(
    #cD = CD
    data = cnt_f
    )
  cnt_f[cnt_f == 0] <- 1
  normRT <- log2(t(t(cnt_f / libsizes)) * mean(libsizes))
  kClust <- kCluster(cD = normRT,
                     maxK = 100,
                     matrixFile = NULL,
                     replicates = des2[[group_col]],
                     algorithm = "Lloyd",
                     B = 1000,
                     sdm = 1)
  mkClust <- makeClusters(aM = kClust, 
                          cD = normRT,
                          threshold = 1)
  selected_g <- length(mkClust@listData)
  names(mkClust@listData) <- paste0("Cls",
                                    1:selected_g,
                                    "_")
  cls <- unlist(mkClust@listData)
  cls <- rownames(cnt_f)[cls]
  cluster_labels <- as.numeric(gsub("^Cls",
                                    "",
                                    gsub("_.*$",
                                         "",
                                         names(unlist(mkClust@listData)))))
  names(cluster_labels) <- cls

  artifacts$backend_details <- list(engine="clusterSeq", k=clusterseq_k, selected_k=selected_g)
} else {
  fail(paste0("Unknown backend: ", backend))
}

# ---------- write outputs ----------
gene_clusters <- data.frame(
  gene_id = genes_f,
  cluster = cluster_labels,
  stringsAsFactors = FALSE
)
data.table::fwrite(gene_clusters, file=file.path(out_dir, "gene_clusters.tsv"), sep="\t", quote=FALSE)

# Write per-cluster gene lists for downstream GO/enrichment plugins
dir.create(file.path(out_dir, "cluster_gene_lists"), recursive=TRUE, showWarnings=FALSE)
cluster_levels <- sort(unique(cluster_labels))
cluster_gene_lists <- lapply(cluster_levels, function(k) {
  gg <- genes_f[cluster_labels == k]
  fn <- sprintf("cluster_gene_lists/cluster_%s_genes.txt", k)
  writeLines(gg, con=file.path(out_dir, fn))
  data.frame(cluster=k, n_genes=length(gg), file=fn, stringsAsFactors=FALSE)
})
cluster_gene_lists <- do.call(rbind, cluster_gene_lists)
data.table::fwrite(cluster_gene_lists, file=file.path(out_dir, "cluster_gene_lists.tsv"), sep="\t", quote=FALSE)

data.table::fwrite(data.frame(gene_id=genes_f, stringsAsFactors=FALSE),
                   file=file.path(out_dir, "genes_used.tsv"), sep="\t", quote=FALSE)

if (!is.null(model_scores)) {
  data.table::fwrite(model_scores, file=file.path(out_dir, "model_scores.tsv"), sep="\t", quote=FALSE)
}

cs <- as.data.frame(table(cluster_labels), stringsAsFactors=FALSE)
colnames(cs) <- c("cluster","n_genes")
data.table::fwrite(cs, file=file.path(out_dir, "cluster_sizes.tsv"), sep="\t", quote=FALSE)

# per-condition mean profile (logCPM)
logcpm_f <- log2(sweep(cnt_f, 2, colSums(cnt_f)/1e6, "/") + 1)
conds_chr <- as.character(conds)
conds_levels <- levels(conds)
profile <- do.call(cbind, lapply(conds_levels, function(g) {
  cols <- which(conds_chr == g)
  if (length(cols) == 1) logcpm_f[, cols] else rowMeans(logcpm_f[, cols, drop=FALSE])
}))
colnames(profile) <- conds_levels

cluster_profile <- do.call(rbind, lapply(cluster_levels, function(k) {
  idxg <- which(cluster_labels == k)
  if (length(idxg) == 1) profile[idxg, , drop=FALSE] else colMeans(profile[idxg, , drop=FALSE])
}))
cluster_profile <- as.data.frame(cluster_profile)
cluster_profile$cluster <- cluster_levels
cluster_profile <- cluster_profile[, c("cluster", conds_levels), drop=FALSE]
data.table::fwrite(cluster_profile, file=file.path(out_dir, "cluster_profiles.tsv"), sep="\t", quote=FALSE)

# Representative genes: top by variance across conditions
rep_list <- lapply(cluster_levels, function(k) {
  idxg <- which(cluster_labels == k)
  vv <- apply(profile[idxg, , drop=FALSE], 1, var)
  ord <- order(vv, decreasing=TRUE)
  topn <- min(20, length(ord))
  data.frame(cluster=k, gene_id=genes_f[idxg][ord[seq_len(topn)]], score=vv[ord[seq_len(topn)]], stringsAsFactors=FALSE)
})
rep_df <- do.call(rbind, rep_list)
data.table::fwrite(rep_df, file=file.path(out_dir, "representative_genes.tsv"), sep="\t", quote=FALSE)

# ---------- plots helper ----------
save_plot_png <- function(path, plot_expr) {
  png(path, width=1100, height=750, res=120)
  on.exit(dev.off(), add=TRUE)
  plot_expr
}

# ---------- additional visualizations ----------
save_plot_png(file.path(out_dir, "cluster_condition_heatmap.png"), {
  mat <- as.matrix(cluster_profile[, conds_levels, drop=FALSE])
  rownames(mat) <- paste0("C", cluster_profile$cluster)
  op <- par(mar=c(8,7,3,2))
  heatmap(mat, Rowv=NA, Colv=NA, scale="row", margins=c(8,8), main="Cluster x condition heatmap (row-scaled)")
  par(op)
})

# Representative gene profiles per cluster (top 10 by variance)
dir.create(file.path(out_dir, "rep_gene_profiles"), recursive=TRUE, showWarnings=FALSE)
top_per_cluster <- 10
for (k in cluster_levels) {
  rr <- rep_df[rep_df$cluster == k, , drop=FALSE]
  if (nrow(rr) < 1) next
  rr <- rr[order(rr$score, decreasing=TRUE), , drop=FALSE]
  rr <- rr[seq_len(min(top_per_cluster, nrow(rr))), , drop=FALSE]
  idxg <- match(rr$gene_id, genes_f)
  matg <- profile[idxg, , drop=FALSE]
  fn <- sprintf("rep_gene_profiles/cluster_%s_rep_profiles.png", k)
  save_plot_png(file.path(out_dir, fn), {
    x <- seq_along(conds_levels)
    op <- par(mar=c(7,5,3,1))
    matplot(x, t(matg), type="l", lty=1, xaxt="n", xlab="Condition", ylab="logCPM", main=paste0("Representative genes: C", k))
    axis(1, at=x, labels=conds_levels, las=2)
    legend("topright", legend=rr$gene_id, cex=0.6, bty="n")
    par(op)
  })
}

# UMAP/PCA embedding of genes (colored by cluster)
embed <- NULL
embed_method <- NULL
if (suppressWarnings(requireNamespace("uwot", quietly=TRUE))) {
  embed <- tryCatch({
    uwot::umap(profile, n_neighbors=15, min_dist=0.1, metric="euclidean", ret_model=FALSE, verbose=FALSE)
  }, error=function(e) NULL)
  if (!is.null(embed) && is.matrix(embed) && ncol(embed) >= 2) {
    embed_method <- "UMAP (uwot)"
  } else {
    embed <- NULL
  }
}
if (is.null(embed)) {
  pc <- stats::prcomp(profile, center=TRUE, scale.=TRUE)
  embed <- pc$x[, 1:2, drop=FALSE]
  embed_method <- "PCA"
}

embed_df <- data.frame(gene_id=genes_f, cluster=cluster_labels, dim1=embed[,1], dim2=embed[,2], stringsAsFactors=FALSE)
data.table::fwrite(embed_df, file=file.path(out_dir, "gene_embedding.tsv"), sep="\t", quote=FALSE)

save_plot_png(file.path(out_dir, "gene_embedding.png"), {
  op <- par(mar=c(5,5,3,1))
  plot(embed_df$dim1, embed_df$dim2, pch=16, cex=0.5, col=embed_df$cluster,
       xlab="dim1", ylab="dim2", main=paste0("Gene embedding colored by cluster (", embed_method, ")"))
  par(op)
})

save_plot_png(file.path(out_dir, "cluster_sizes.png"), {
  op <- par(mar=c(7,5,3,1))
  barplot(height=cs$n_genes, names.arg=cs$cluster, las=2, ylab="Number of genes", main="Gene cluster sizes")
  par(op)
})

save_plot_png(file.path(out_dir, "cluster_profiles.png"), {
  mat <- as.matrix(cluster_profile[, conds_levels, drop=FALSE])
  x <- seq_along(conds_levels)
  op <- par(mar=c(7,5,3,1))
  plot(x, mat[1,], type="l", xaxt="n", xlab="Condition", ylab="Mean logCPM", main="Cluster mean profiles")
  axis(1, at=x, labels=conds_levels, las=2)
  if (nrow(mat) > 1) {
    for (i in 2:nrow(mat)) lines(x, mat[i,])
  }
  legend("topright", legend=paste0("C", cluster_profile$cluster), cex=0.8, bty="n")
  par(op)
})

artifacts$selected_g <- selected_g
artifacts$gene_lists <- list(
  per_cluster_tsv="cluster_gene_lists.tsv",
  per_cluster_dir="cluster_gene_lists",
  universe="genes_used.tsv"
)
artifacts$outputs <- list(
  tables=c("gene_clusters.tsv","genes_used.tsv","cluster_sizes.tsv","cluster_profiles.tsv","representative_genes.tsv","cluster_gene_lists.tsv","gene_embedding.tsv"),
  plots=c("cluster_sizes.png","cluster_profiles.png","cluster_condition_heatmap.png","gene_embedding.png"),
  rep_gene_profiles_dir="rep_gene_profiles",
  session_info="sessionInfo.txt"
)
jsonlite::write_json(artifacts, path=file.path(out_dir, "artifacts.json"), auto_unbox=TRUE, pretty=TRUE)

cat("[rnaseq_genecluster_modelbased] done\n")
