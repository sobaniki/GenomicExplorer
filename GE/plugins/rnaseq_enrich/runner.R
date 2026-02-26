#!/usr/bin/env Rscript
# RNA-seq enrichment: clusterProfiler ORA / goseq / topGO
suppressWarnings(suppressMessages({
  library(jsonlite)
  library(data.table)
}))

fail <- function(msg, code=1) {
  cat("[rnaseq_enrich][ERROR] ", msg, "\n", sep="")
  quit(status=code)
}

as_bool <- function(x) {
  if (is.null(x)) return(FALSE)
  if (is.logical(x)) return(isTRUE(x))
  if (is.character(x)) return(tolower(x) %in% c("1","true","t","yes","y"))
  if (is.numeric(x)) return(x != 0)
  FALSE
}

`%||%` <- function(a, b) { if (is.null(a) || length(a)==0 || (is.character(a)&&a=="")) b else a }

read_lines_nonempty <- function(path) {
  if (is.null(path) || path == "" || !file.exists(path)) return(character(0))
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("tsv", "csv")) {
    dt <- tryCatch(
      data.table::fread(path, sep=ifelse(ext == "csv", ",", "\t"), header=TRUE, data.table=FALSE),
      error=function(e) NULL
    )
    if (!is.null(dt) && nrow(dt) > 0 && ncol(dt) >= 1) {
      x <- as.character(dt[[1]])
      x <- trimws(x)
      x <- x[nzchar(x)]
      return(unique(x))
    }
  }
  x <- readLines(path, warn=FALSE, encoding="UTF-8")
  x <- trimws(x)
  x <- x[nzchar(x)]
  unique(x)
}

read_gene_length_tsv <- function(path, gene_col=1L, len_col=2L) {
  if (is.null(path) || path=="" || !file.exists(path)) return(NULL)
  dt <- tryCatch(data.table::fread(path, sep="\t", header=TRUE, data.table=FALSE), error=function(e) NULL)
  if (is.null(dt) || nrow(dt)<1) return(NULL)
  if (is.character(gene_col)) {
    if (!gene_col %in% colnames(dt)) return(NULL)
    g <- dt[[gene_col]]
  } else {
    if (gene_col < 1 || gene_col > ncol(dt)) return(NULL)
    g <- dt[[gene_col]]
  }
  if (is.character(len_col)) {
    if (!len_col %in% colnames(dt)) return(NULL)
    l <- dt[[len_col]]
  } else {
    if (len_col < 1 || len_col > ncol(dt)) return(NULL)
    l <- dt[[len_col]]
  }
  g <- as.character(g)
  l <- suppressWarnings(as.numeric(l))
  ok <- !is.na(g) & nzchar(g) & is.finite(l) & l > 0
  out <- data.frame(gene=g[ok], length=l[ok], stringsAsFactors=FALSE)
  out <- out[!duplicated(out$gene), , drop=FALSE]
  out
}


# ---------- plotting helpers ----------
have_pkg <- function(pkg) requireNamespace(pkg, quietly=TRUE)

safe_png <- function(path, fun){
  tryCatch({
    png(path, width=1200, height=800)
    on.exit(dev.off(), add=TRUE)
    fun()
  }, error=function(e){
    try(dev.off(), silent=TRUE)
  })
}

safe_ggsave <- function(path, plot, width=12, height=8) {
  tryCatch({
    if (!have_pkg('ggplot2')) return(FALSE)
    ggplot2::ggsave(filename=path, plot=plot, width=width, height=height, dpi=150)
    TRUE
  }, error=function(e){
    FALSE
  })
}

parse_ratio <- function(x){
  if (is.null(x) || length(x)==0) return(NA_real_)
  x <- as.character(x)[1]
  sp <- strsplit(x, "/", fixed=TRUE)[[1]]
  if (length(sp) != 2) return(suppressWarnings(as.numeric(x)))
  a <- suppressWarnings(as.numeric(sp[1])); b <- suppressWarnings(as.numeric(sp[2]))
  if (!is.finite(a) || !is.finite(b) || b == 0) return(NA_real_)
  a / b
}

go_term_map <- function(go_ids){
  go_ids <- as.character(go_ids)
  out <- setNames(rep(NA_character_, length(go_ids)), go_ids)
  if (!have_pkg('GO.db') || !have_pkg('AnnotationDbi')) return(out)
  ok <- unique(go_ids[!is.na(go_ids) & nzchar(go_ids)])
  if (length(ok) < 1) return(out)
  m <- tryCatch({
    AnnotationDbi::select(GO.db::GO.db, keys=ok, keytype='GOID', columns=c('TERM'))
  }, error=function(e) NULL)
  if (is.null(m) || nrow(m) < 1) return(out)
  m <- m[!is.na(m$TERM) & nzchar(m$TERM), , drop=FALSE]
  if (nrow(m) < 1) return(out)
  # If duplicated GOID, keep first
  m <- m[!duplicated(m$GOID), , drop=FALSE]
  out[m$GOID] <- as.character(m$TERM)
  out
}

standardize_clusterprofiler <- function(df){
  if (is.null(df) || nrow(df) < 1) return(data.frame())
  term <- if ('Description' %in% colnames(df)) as.character(df$Description) else if ('ID' %in% colnames(df)) as.character(df$ID) else as.character(df[[1]])
  count <- if ('Count' %in% colnames(df)) suppressWarnings(as.numeric(df$Count)) else rep(NA_real_, nrow(df))
  padj <- if ('p.adjust' %in% colnames(df)) suppressWarnings(as.numeric(df[['p.adjust']])) else if ('pvalue' %in% colnames(df)) suppressWarnings(as.numeric(df$pvalue)) else rep(NA_real_, nrow(df))
  gr <- if ('GeneRatio' %in% colnames(df)) vapply(df$GeneRatio, parse_ratio, numeric(1)) else rep(NA_real_, nrow(df))
  br <- if ('BgRatio' %in% colnames(df)) vapply(df$BgRatio, parse_ratio, numeric(1)) else rep(NA_real_, nrow(df))
  # set size (genes annotated to term in universe)
  bg_num <- if ('BgRatio' %in% colnames(df)) vapply(df$BgRatio, function(x){
    x <- as.character(x)[1]; sp <- strsplit(x, "/", fixed=TRUE)[[1]]
    if (length(sp) == 2) suppressWarnings(as.numeric(sp[1])) else NA_real_
  }, numeric(1)) else rep(NA_real_, nrow(df))
  set_size <- if ('setSize' %in% colnames(df)) suppressWarnings(as.numeric(df$setSize)) else bg_num
  rich_factor <- ifelse(is.finite(gr) & is.finite(br) & br > 0, gr / br, NA_real_)

  out <- data.frame(term=term, count=count, gene_ratio=gr, bg_ratio=br, rich_factor=rich_factor, set_size=set_size, padj=padj, stringsAsFactors=FALSE)
  out$mlog10 <- -log10(pmax(out$padj, 1e-300))
  out
}

standardize_goseq <- function(df, total_de, total_uni){
  if (is.null(df) || nrow(df) < 1) return(data.frame())
  go_id <- if ('category' %in% colnames(df)) as.character(df$category) else as.character(df[[1]])
  pval <- if ('over_represented_pvalue' %in% colnames(df)) suppressWarnings(as.numeric(df$over_represented_pvalue)) else rep(NA_real_, nrow(df))
  padj <- if ('over_represented_padj' %in% colnames(df)) suppressWarnings(as.numeric(df$over_represented_padj)) else p.adjust(pval, method='BH')
  cnt <- if ('numDEInCat' %in% colnames(df)) suppressWarnings(as.numeric(df$numDEInCat)) else rep(NA_real_, nrow(df))
  num_in <- if ('numInCat' %in% colnames(df)) suppressWarnings(as.numeric(df$numInCat)) else rep(NA_real_, nrow(df))
  # term names
  tm <- go_term_map(go_id)
  term <- unname(tm[go_id])
  term[is.na(term) | !nzchar(term)] <- go_id[is.na(term) | !nzchar(term)]

  gr <- if (is.finite(total_de) && total_de > 0) cnt / total_de else NA_real_
  bg_ratio <- if (is.finite(total_uni) && total_uni > 0) num_in / total_uni else NA_real_
  # enrichment (DE proportion in category / global DE proportion)
  de_prop_cat <- ifelse(is.finite(num_in) & num_in > 0, cnt / num_in, NA_real_)
  de_prop_bg <- if (is.finite(total_de) && is.finite(total_uni) && total_uni > 0) total_de / total_uni else NA_real_
  rich_factor <- ifelse(is.finite(de_prop_cat) & is.finite(de_prop_bg) & de_prop_bg > 0, de_prop_cat / de_prop_bg, NA_real_)
  out <- data.frame(term=term, count=cnt, gene_ratio=gr, bg_ratio=bg_ratio, rich_factor=rich_factor, set_size=num_in, padj=padj, stringsAsFactors=FALSE)
  out$mlog10 <- -log10(pmax(out$padj, 1e-300))
  out
}

standardize_topgo <- function(df){
  if (is.null(df) || nrow(df) < 1) return(data.frame())
  term <- if ('Term' %in% colnames(df)) as.character(df$Term) else as.character(df[[1]])
  cnt <- if ('Significant' %in% colnames(df)) suppressWarnings(as.numeric(df$Significant)) else rep(NA_real_, nrow(df))
  ann <- if ('Annotated' %in% colnames(df)) suppressWarnings(as.numeric(df$Annotated)) else rep(NA_real_, nrow(df))
  pval <- if ('pvalue' %in% colnames(df)) suppressWarnings(as.numeric(df$pvalue)) else rep(NA_real_, nrow(df))
  padj <- p.adjust(pval, method='BH')
  gr <- ifelse(is.finite(ann) & ann > 0, cnt/ann, NA_real_)
  out <- data.frame(term=term, count=cnt, gene_ratio=gr, bg_ratio=NA_real_, rich_factor=NA_real_, set_size=ann, padj=padj, stringsAsFactors=FALSE)
  out$mlog10 <- -log10(pmax(out$padj, 1e-300))
  out
}

make_dotplot <- function(std, out_png, title='GO enrichment dotplot', top_n=plot_top_n){
  safe_png(out_png, function(){
    if (is.null(std) || nrow(std) < 1) {
      plot.new(); title(main=title); text(0.5, 0.5, 'No enriched terms')
      return()
    }
    s <- std
    s <- s[is.finite(s$padj), , drop=FALSE]
    if (nrow(s) < 1) s <- std
    s <- s[order(s$padj, decreasing=FALSE), , drop=FALSE]
    s <- head(s, top_n)
    # shorten labels
    lab <- as.character(s$term)
    if (is.finite(label_max_chars) && label_max_chars > 0) {
      lab <- ifelse(nchar(lab) > label_max_chars,
                    paste0(substr(lab, 1, max(1, label_max_chars-3)), '...'),
                    lab)
    }

    get_metric <- function(df, key){
      key <- tolower(as.character(key))
      if (key %in% c('gene_ratio','generatio')) return(df$gene_ratio)
      if (key %in% c('count')) return(df$count)
      if (key %in% c('mlog10','-log10','log10','neglog10')) return(df$mlog10)
      if (key %in% c('rich_factor','richfactor')) return(df$rich_factor)
      if (key %in% c('bg_ratio','bgratio')) return(df$bg_ratio)
      if (key %in% c('set_size','setsize','size')) return(df$set_size)
      return(NULL)
    }
    pick_metric <- function(df, key, fallback){
      v <- get_metric(df, key)
      if (is.null(v) || all(!is.finite(v))) v <- get_metric(df, fallback)
      v
    }

    x <- pick_metric(s, plot_x_metric, 'gene_ratio')
    if (is.null(x) || all(!is.finite(x))) x <- s$count
    y <- rev(seq_len(nrow(s)))
    # sizes
    cnt <- pick_metric(s, plot_size_metric, 'count')
    if (is.null(cnt) || all(!is.finite(cnt))) cnt <- s$count
    if (all(!is.finite(cnt))) cnt <- rep(1, nrow(s))
    cnt[!is.finite(cnt)] <- min(cnt[is.finite(cnt)], na.rm=TRUE)
    cexv <- 0.8 + 2.5 * (cnt - min(cnt, na.rm=TRUE)) / max(1e-9, (max(cnt, na.rm=TRUE) - min(cnt, na.rm=TRUE)))

    colv <- pick_metric(s, plot_color_metric, 'mlog10')
    if (is.null(colv) || all(!is.finite(colv))) colv <- s$mlog10
    if (!all(is.finite(colv))) colv <- rep(1, nrow(s))
    z <- (colv - min(colv, na.rm=TRUE)) / max(1e-9, (max(colv, na.rm=TRUE) - min(colv, na.rm=TRUE)))
    cols <- grDevices::colorRampPalette(c('grey70', 'black'))(100)
    colp <- cols[pmax(1, pmin(100, 1 + floor(z*99)))]

    op <- par(mar=c(5, 12, 4, 2) + 0.1)
    on.exit(par(op), add=TRUE)
    plot(x, y, pch=16, cex=cexv, col=colp, yaxt='n', xlab=plot_x_metric, ylab='', main=title)
    axis(2, at=y, labels=rev(lab), las=2, cex.axis=0.7)
  })
}

make_barplot <- function(std, out_png, title='GO enrichment barplot', top_n=plot_top_n){
  safe_png(out_png, function(){
    if (is.null(std) || nrow(std) < 1) {
      plot.new(); title(main=title); text(0.5, 0.5, 'No enriched terms')
      return()
    }
    s <- std
    s <- s[is.finite(s$padj), , drop=FALSE]
    if (nrow(s) < 1) s <- std
    s <- s[order(s$padj, decreasing=FALSE), , drop=FALSE]
    s <- head(s, top_n)
    lab <- as.character(s$term)
    if (is.finite(label_max_chars) && label_max_chars > 0) {
      lab <- ifelse(nchar(lab) > label_max_chars,
                    paste0(substr(lab, 1, max(1, label_max_chars-3)), '...'),
                    lab)
    }
    score <- NULL
    if (!is.null(plot_bar_metric)) {
      if (tolower(plot_bar_metric) %in% c('gene_ratio','generatio')) score <- s$gene_ratio
      if (tolower(plot_bar_metric) %in% c('count')) score <- s$count
      if (tolower(plot_bar_metric) %in% c('rich_factor','richfactor')) score <- s$rich_factor
      if (tolower(plot_bar_metric) %in% c('bg_ratio','bgratio')) score <- s$bg_ratio
      if (tolower(plot_bar_metric) %in% c('set_size','setsize','size')) score <- s$set_size
      if (tolower(plot_bar_metric) %in% c('mlog10','-log10','log10','neglog10')) score <- s$mlog10
    }
    if (is.null(score) || all(!is.finite(score))) score <- s$mlog10
    if (!all(is.finite(score))) score <- rep(0, nrow(s))
    op <- par(mar=c(5, 12, 4, 2) + 0.1)
    on.exit(par(op), add=TRUE)
    barplot(rev(score), names.arg=rev(lab), horiz=TRUE, las=1, cex.names=0.7,
            xlab=plot_bar_metric, main=title)
  })
}

make_dot_bar <- function(std, prefix, ont, out_dir, label){
  dot_png <- file.path(out_dir, sprintf('%s_dotplot_%s.png', prefix, ont))
  bar_png <- file.path(out_dir, sprintf('%s_barplot_%s.png', prefix, ont))
  make_dotplot(std, dot_png, title=paste0(label, ' ', ont))
  make_barplot(std, bar_png, title=paste0(label, ' ', ont))
  list(dot=basename(dot_png), bar=basename(bar_png))
}


# Build gene -> GO mapping from OrgDb
build_gene2go_from_orgdb <- function(keys, orgdb_pkg, keytype) {
  if (!requireNamespace("AnnotationDbi", quietly=TRUE)) fail("AnnotationDbi not installed")
  if (!requireNamespace(orgdb_pkg, quietly=TRUE)) fail(paste0("OrgDb package not installed: ", orgdb_pkg))
  suppressWarnings(suppressMessages(library(orgdb_pkg, character.only=TRUE)))
  orgdb <- get(orgdb_pkg)

  # Prefer GOALL/ONTOLOGYALL if available
  cols <- c("GOALL", "ONTOLOGYALL")
  cols2 <- c("GO", "ONTOLOGY")
  have_all <- all(cols %in% AnnotationDbi::columns(orgdb))
  have_plain <- all(cols2 %in% AnnotationDbi::columns(orgdb))

  if (have_all) {
    m <- AnnotationDbi::select(orgdb, keys=as.character(keys), keytype=keytype, columns=cols)
    m <- m[!is.na(m$GOALL) & nzchar(m$GOALL), , drop=FALSE]
    # Keep only valid ontology strings (BP/MF/CC)
    m$ONTOLOGYALL <- toupper(as.character(m$ONTOLOGYALL))
    gene2go <- split(as.character(m$GOALL), as.character(m[[keytype]]))
    gene2go <- lapply(gene2go, function(v) unique(v[!is.na(v) & nzchar(v)]))
    return(list(gene2go=gene2go, mapping=m, col_go="GOALL", col_ont="ONTOLOGYALL"))
  }

  if (have_plain) {
    m <- AnnotationDbi::select(orgdb, keys=as.character(keys), keytype=keytype, columns=cols2)
    m <- m[!is.na(m$GO) & nzchar(m$GO), , drop=FALSE]
    m$ONTOLOGY <- toupper(as.character(m$ONTOLOGY))
    gene2go <- split(as.character(m$GO), as.character(m[[keytype]]))
    gene2go <- lapply(gene2go, function(v) unique(v[!is.na(v) & nzchar(v)]))
    return(list(gene2go=gene2go, mapping=m, col_go="GO", col_ont="ONTOLOGY"))
  }

  fail(paste0("OrgDb does not provide GO columns (GOALL/ONTOLOGYALL or GO/ONTOLOGY): ", orgdb_pkg))
}

# Parse args
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 4 || args[1] != "--params" || args[3] != "--out") {
  cat("Usage: runner.R --params params.json --out out_dir\n")
  quit(status=2)
}

params_path <- args[2]
out_dir <- args[4]
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

cat("[rnaseq_enrich] start\n")
cat("[rnaseq_enrich] params_path=", params_path, "\n", sep="")
cat("[rnaseq_enrich] out_dir=", out_dir, "\n", sep="")

params <- tryCatch(jsonlite::fromJSON(params_path), error=function(e) fail(paste("Failed to parse params.json:", e$message)))

# Optional extra library paths (e.g., OrgDb built by rnaseq_build_orgdb)
apply_extra_libs <- function(x) {
  if (is.null(x) || length(x) == 0) return()
  libs <- c()
  if (is.character(x)) libs <- as.character(x)
  if (is.list(x)) libs <- unlist(x)
  libs <- unique(trimws(as.character(libs)))
  libs <- libs[nzchar(libs)]

  # split common separators (allow users to paste multiple paths)
  libs2 <- c()
  for (s in libs) {
    if (is.null(s) || length(s) == 0) next
    ss <- trimws(as.character(s)[1])
    if (!nzchar(ss)) next
    parts <- unlist(strsplit(ss, "[;\n,]+", perl=TRUE))
    parts <- trimws(parts)
    parts <- parts[nzchar(parts)]

    # PATH-like ':' separator (mainly for Linux). Avoid breaking Windows drive letters by using this
    # only when the whole string is not a directory but all split parts are.
    if (length(parts) == 1 && !dir.exists(parts[1]) && grepl(":", ss, fixed=TRUE)) {
      parts_colon <- unlist(strsplit(ss, ":", fixed=TRUE))
      parts_colon <- trimws(parts_colon)
      parts_colon <- parts_colon[nzchar(parts_colon)]
      if (length(parts_colon) > 1 && all(dir.exists(parts_colon))) {
        parts <- parts_colon
      }
    }

    libs2 <- c(libs2, parts)
  }

  libs <- unique(libs2)
  if (length(libs) == 0) return()

  # Expand user-friendly inputs:
  # - If user points to an OrgDb build out_dir, also add out_dir/Rlib
  # - If user points to a package directory itself, also add its parent
  expanded <- c()
  for (p in libs) {
    p <- trimws(as.character(p)[1])
    if (!nzchar(p)) next

    # prefer Rlib if present
    if (dir.exists(file.path(p, "Rlib"))) expanded <- c(expanded, file.path(p, "Rlib"))

    expanded <- c(expanded, p)

    # if p looks like a package dir (DESCRIPTION exists), use its parent as lib root
    if (file.exists(file.path(p, "DESCRIPTION"))) expanded <- c(expanded, dirname(p))
  }

  expanded <- unique(expanded)
  expanded <- expanded[dir.exists(expanded)]
  if (length(expanded) > 0) {
    .libPaths(c(expanded, .libPaths()))
    cat("[rnaseq_enrich] .libPaths prepended: ", paste(expanded, collapse=", "), "\n", sep="")
  }
}


apply_extra_libs(params$r_libs %||% params$r_lib %||% params$rlibs %||% params$extra_libs)

engine <- tolower(params$engine %||% "clusterprofiler")
ontology <- toupper(params$ontology %||% "BP")
# Allow "ALL" to run BP/MF/CC
ontologies <- if (ontology == "ALL") c("BP","MF","CC") else c(ontology)

orgdb_pkg <- params$orgdb %||% params$orgdb_package %||% "org.Hs.eg.db"
keytype <- params$keytype %||% "ENTREZID"

p_adjust_method <- params$p_adjust_method %||% "BH"
pvalue_cutoff <- as.numeric(params$pvalue_cutoff %||% 0.05)
qvalue_cutoff <- as.numeric(params$qvalue_cutoff %||% 0.2)

# ---------- plot settings (GUI-controlled) ----------
mode <- tolower(params$mode %||% "enrich")

plot_top_n <- suppressWarnings(as.integer(params$plot_top_n %||% params$top_n %||% 20))
if (!is.finite(plot_top_n) || plot_top_n < 1) plot_top_n <- 20L

label_max_chars <- suppressWarnings(as.integer(params$plot_label_max_chars %||% params$label_max_chars %||% 60))
if (!is.finite(label_max_chars) || label_max_chars < 0) label_max_chars <- 60L

# metric keys:
#   gene_ratio, count, mlog10, rich_factor, bg_ratio, set_size
plot_x_metric <- tolower(as.character(params$plot_x_metric %||% "gene_ratio"))
plot_size_metric <- tolower(as.character(params$plot_size_metric %||% "count"))
plot_color_metric <- tolower(as.character(params$plot_color_metric %||% "mlog10"))
plot_bar_metric <- tolower(as.character(params$plot_bar_metric %||% "mlog10"))

make_cnetplot <- as_bool(params$make_cnetplot %||% params$cnetplot %||% FALSE)
make_emapplot <- as_bool(params$make_emapplot %||% params$emapplot %||% FALSE)

# cluster x GO heatmap mode
cluster_gene_lists_tsv <- params$cluster_gene_lists_tsv %||% params$cluster_gene_lists %||% ""
cluster_heatmap_top_terms <- suppressWarnings(as.integer(params$cluster_heatmap_top_terms %||% params$heatmap_top_terms %||% 30))
if (!is.finite(cluster_heatmap_top_terms) || cluster_heatmap_top_terms < 5) cluster_heatmap_top_terms <- 30L

# Inputs
# gene_list can be explicit list file, or derived from deg_table

gene_list_path <- params$gene_list %||% ""
universe_path <- params$universe_list %||% params$universe %||% ""

deg_table <- params$deg_table %||% ""
gene_col <- params$gene_col %||% "gene"
padj_col <- params$padj_col %||% "padj"
logfc_col <- params$logfc_col %||% "log2FoldChange"
use_deg_filter <- as_bool(params$use_deg_filter %||% FALSE)
fdr_cutoff <- as.numeric(params$fdr_cutoff %||% 0.05)
lfc_cutoff <- as.numeric(params$lfc_cutoff %||% 0.0)

# goseq specific
length_tsv <- params$gene_length_tsv %||% params$gene_length_table %||% params$gene_length %||% ""
length_gene_col <- params$length_gene_col %||% "gene"
length_len_col <- params$length_len_col %||% "length"

# topGO specific
# Keep it simple for now

topgo_algorithm <- params$topgo_algorithm %||% "weight01"
topgo_statistic <- params$topgo_statistic %||% "fisher"

# Resolve gene list (skip when running cluster heatmap mode)
selected_genes <- character(0)
universe_genes <- character(0)

if (mode != "cluster_heatmap") {
if (nzchar(deg_table) && file.exists(deg_table)) {
  dt <- tryCatch(data.table::fread(deg_table, sep="\t", header=TRUE, data.table=FALSE), error=function(e) NULL)
  if (!is.null(dt) && nrow(dt) > 0) {
    if (!(gene_col %in% colnames(dt))) {
      # fallback: first column
      gene_col_eff <- colnames(dt)[1]
    } else {
      gene_col_eff <- gene_col
    }
    g <- as.character(dt[[gene_col_eff]])
    g <- g[!is.na(g) & nzchar(g)]

    if (use_deg_filter) {
      # padj
      if (padj_col %in% colnames(dt)) {
        pv <- suppressWarnings(as.numeric(dt[[padj_col]]))
      } else if ("FDR" %in% colnames(dt)) {
        pv <- suppressWarnings(as.numeric(dt[["FDR"]]))
      } else {
        pv <- rep(NA_real_, nrow(dt))
      }
      # logfc
      if (logfc_col %in% colnames(dt)) {
        lfcv <- suppressWarnings(as.numeric(dt[[logfc_col]]))
      } else if ("logFC" %in% colnames(dt)) {
        lfcv <- suppressWarnings(as.numeric(dt[["logFC"]]))
      } else {
        lfcv <- rep(0, nrow(dt))
      }

      keep <- rep(TRUE, length(g))
      if (is.finite(fdr_cutoff)) keep <- keep & is.finite(pv) & pv <= fdr_cutoff
      if (is.finite(lfc_cutoff) && lfc_cutoff > 0) keep <- keep & is.finite(lfcv) & abs(lfcv) >= lfc_cutoff
      selected_genes <- unique(g[keep])
    } else {
      selected_genes <- unique(g)
    }
  }
}

if (length(selected_genes) < 1) {
  if (!nzchar(gene_list_path) || !file.exists(gene_list_path)) {
    fail("No genes provided. Set gene_list or deg_table.")
  }
  selected_genes <- read_lines_nonempty(gene_list_path)
}

if (length(selected_genes) < 1) fail("Selected gene list is empty")

universe_genes <- read_lines_nonempty(universe_path)
if (length(universe_genes) < 1) {
  cat("[rnaseq_enrich] universe not provided; using selected genes as universe (not recommended)\n")
  universe_genes <- selected_genes
}

# Write resolved inputs for traceability
writeLines(selected_genes, file.path(out_dir, "selected_genes.txt"))
writeLines(universe_genes, file.path(out_dir, "universe_genes.txt"))
}

# Run engines
artifacts <- list()

do_cluster_heatmap <- function() {
  if (!nzchar(cluster_gene_lists_tsv) || !file.exists(cluster_gene_lists_tsv)) {
    fail("cluster_heatmap mode requires cluster_gene_lists_tsv")
  }

  if (!requireNamespace("clusterProfiler", quietly=TRUE)) fail("clusterProfiler not installed")
  if (!requireNamespace("AnnotationDbi", quietly=TRUE)) fail("AnnotationDbi not installed")
  if (!requireNamespace(orgdb_pkg, quietly=TRUE)) fail(paste0("OrgDb package not installed: ", orgdb_pkg))

  suppressWarnings(suppressMessages({
    library(clusterProfiler)
    library(AnnotationDbi)
    library(orgdb_pkg, character.only=TRUE)
  }))
  OrgDb <- get(orgdb_pkg)

  cdt <- tryCatch(data.table::fread(cluster_gene_lists_tsv, sep="\t", header=TRUE, data.table=FALSE), error=function(e) NULL)
  if (is.null(cdt) || nrow(cdt) < 1) fail("cluster_gene_lists_tsv is empty")
  if (!('cluster' %in% colnames(cdt))) cdt$cluster <- as.character(cdt[[1]])
  if (!('file' %in% colnames(cdt))) {
    if ('path' %in% colnames(cdt)) cdt$file <- as.character(cdt[['path']]) else fail("cluster_gene_lists_tsv must contain 'file' column")
  }

  base_dir <- dirname(cluster_gene_lists_tsv)
  gene_lists <- list()
  for (i in seq_len(nrow(cdt))) {
    cl <- as.character(cdt$cluster[i])
    fp <- as.character(cdt$file[i])
    if (!nzchar(cl) || !nzchar(fp)) next
    gg <- read_lines_nonempty(file.path(base_dir, fp))
    if (length(gg) > 0) gene_lists[[cl]] <- gg
  }
  if (length(gene_lists) < 1) fail("No cluster gene lists could be loaded")

  # Universe: prefer provided universe_list, otherwise union of all clusters
  universe_hm <- read_lines_nonempty(universe_path)
  if (length(universe_hm) < 1) {
    universe_hm <- sort(unique(unlist(gene_lists)))
    cat("[rnaseq_enrich] cluster_heatmap: universe_list not provided; using union of clusters\n")
  }

  writeLines(universe_hm, file.path(out_dir, "universe_genes.txt"))

  hm_dir <- file.path(out_dir, "cluster_heatmap")
  dir.create(hm_dir, recursive=TRUE, showWarnings=FALSE)

  for (ont in ontologies) {
    res_by_cluster <- list()
    for (cl in names(gene_lists)) {
      gg <- gene_lists[[cl]]
      if (length(gg) < 5) next
      ego <- tryCatch({
        clusterProfiler::enrichGO(
          gene=gg,
          universe=universe_hm,
          OrgDb=OrgDb,
          keyType=keytype,
          ont=ont,
          pAdjustMethod=p_adjust_method,
          pvalueCutoff=pvalue_cutoff,
          qvalueCutoff=qvalue_cutoff,
          readable=FALSE
        )
      }, error=function(e) NULL)
      if (is.null(ego)) next
      df <- as.data.frame(ego)
      if (is.null(df) || nrow(df) < 1) next
      res_by_cluster[[cl]] <- df
      # save per cluster
      out_tsv <- file.path(hm_dir, sprintf("clusterprofiler_%s_cluster_%s.tsv", ont, cl))
      data.table::fwrite(df, out_tsv, sep="\t")
    }

    clusters <- names(res_by_cluster)
    if (length(clusters) < 1) {
      # still emit empty artifacts
      blank_png <- file.path(hm_dir, sprintf("cluster_go_heatmap_%s.png", ont))
      safe_png(blank_png, function(){ plot.new(); title(main=paste0("Cluster×GO heatmap ", ont)); text(0.5,0.5,"No enriched terms") })
      artifacts[[paste0("cluster_go_heatmap_", ont, "_png")]] <<- file.path("cluster_heatmap", basename(blank_png))
      next
    }

    # pick union of top terms across clusters
    top_ids <- character(0)
    id_to_desc <- list()
    for (cl in clusters) {
      d <- res_by_cluster[[cl]]
      if (!('ID' %in% colnames(d))) next
      d <- d[is.finite(d[['p.adjust']]), , drop=FALSE]
      d <- d[order(d[['p.adjust']]), , drop=FALSE]
      d <- head(d, cluster_heatmap_top_terms)
      ids <- as.character(d[['ID']])
      top_ids <- unique(c(top_ids, ids))
      if ('Description' %in% colnames(d)) {
        for (j in seq_len(nrow(d))) {
          id <- as.character(d[['ID']][j])
          if (!id %in% names(id_to_desc)) id_to_desc[[id]] <- as.character(d[['Description']][j])
        }
      }
    }

    if (length(top_ids) < 1) {
      blank_png <- file.path(hm_dir, sprintf("cluster_go_heatmap_%s.png", ont))
      safe_png(blank_png, function(){ plot.new(); title(main=paste0("Cluster×GO heatmap ", ont)); text(0.5,0.5,"No enriched terms") })
      artifacts[[paste0("cluster_go_heatmap_", ont, "_png")]] <<- file.path("cluster_heatmap", basename(blank_png))
      next
    }

    # matrix: rows=terms, cols=clusters, values=-log10(p.adjust)
    mat <- matrix(0, nrow=length(top_ids), ncol=length(clusters))
    rownames(mat) <- top_ids
    colnames(mat) <- clusters
    for (ci in seq_along(clusters)) {
      cl <- clusters[ci]
      d <- res_by_cluster[[cl]]
      if (is.null(d) || nrow(d) < 1) next
      v <- suppressWarnings(as.numeric(d[['p.adjust']]))
      ids <- as.character(d[['ID']])
      m <- -log10(pmax(v, 1e-300))
      idx <- match(ids, top_ids)
      ok <- is.finite(idx) & is.finite(m)
      mat[idx[ok], ci] <- pmax(mat[idx[ok], ci], m[ok])
    }

    # labels
    desc <- vapply(top_ids, function(id) {
      x <- id_to_desc[[id]]
      if (is.null(x) || !nzchar(x)) id else as.character(x)[1]
    }, character(1))
    if (is.finite(label_max_chars) && label_max_chars > 0) {
      desc <- ifelse(nchar(desc) > label_max_chars,
                     paste0(substr(desc, 1, max(1, label_max_chars-3)), '...'),
                     desc)
    }

    # write tsv
    tsv_path <- file.path(hm_dir, sprintf("cluster_go_heatmap_%s.tsv", ont))
    out_df <- data.frame(ID=top_ids, Description=desc, mat, check.names=FALSE, stringsAsFactors=FALSE)
    data.table::fwrite(out_df, tsv_path, sep="\t")

    # plot heatmap
    png_path <- file.path(hm_dir, sprintf("cluster_go_heatmap_%s.png", ont))
    safe_png(png_path, function(){
      # base heatmap (avoid dependencies)
      op <- par(mar=c(5, 12, 4, 2) + 0.1)
      on.exit(par(op), add=TRUE)
      mat2 <- mat
      rownames(mat2) <- desc
      if (nrow(mat2) == 1 && ncol(mat2) == 1) {
        plot.new(); title(main=paste0("Cluster×GO heatmap ", ont)); text(0.5,0.5, sprintf("%s: %.2f", rownames(mat2)[1], mat2[1,1]))
      } else {
        heatmap(mat2, Rowv=NA, Colv=NA, scale="none", margins=c(6, 12), main=paste0("Cluster×GO heatmap ", ont), xlab="Cluster", ylab="GO term")
      }
    })

    artifacts[[paste0("cluster_go_heatmap_", ont)]] <<- file.path("cluster_heatmap", basename(tsv_path))
    artifacts[[paste0("cluster_go_heatmap_", ont, "_png")]] <<- file.path("cluster_heatmap", basename(png_path))
  }
}

if (mode == "cluster_heatmap") {
  do_cluster_heatmap()
  art_path <- file.path(out_dir, "artifacts.json")
  jsonlite::write_json(list(artifacts=artifacts), art_path, auto_unbox=TRUE, pretty=TRUE)
  cat("[rnaseq_enrich] done (cluster_heatmap)\n")
  quit(status=0)
}

do_clusterprofiler <- function() {
  if (!requireNamespace("clusterProfiler", quietly=TRUE)) fail("clusterProfiler not installed")
  if (!requireNamespace("AnnotationDbi", quietly=TRUE)) fail("AnnotationDbi not installed")
  if (!requireNamespace(orgdb_pkg, quietly=TRUE)) fail(paste0("OrgDb package not installed: ", orgdb_pkg))

  suppressWarnings(suppressMessages({
    library(clusterProfiler)
    library(AnnotationDbi)
    library(orgdb_pkg, character.only=TRUE)
  }))
  OrgDb <- get(orgdb_pkg)

  for (ont in ontologies) {
    ego <- tryCatch({
      clusterProfiler::enrichGO(
        gene=selected_genes,
        universe=universe_genes,
        OrgDb=OrgDb,
        keyType=keytype,
        ont=ont,
        pAdjustMethod=p_adjust_method,
        pvalueCutoff=pvalue_cutoff,
        qvalueCutoff=qvalue_cutoff,
        readable=FALSE
      )
    }, error=function(e) {
      fail(paste0("clusterProfiler enrichGO failed (", ont, "): ", e$message))
    })

    df <- as.data.frame(ego)
    out_tsv <- file.path(out_dir, sprintf("clusterprofiler_enrichGO_%s.tsv", ont))
    data.table::fwrite(df, out_tsv, sep="\t")
    saveRDS(ego, file.path(out_dir, sprintf("clusterprofiler_enrichGO_%s.rds", ont)))
    # plots
    std <- standardize_clusterprofiler(df)
    pp <- make_dot_bar(std, prefix='clusterprofiler', ont=ont, out_dir=out_dir, label='clusterProfiler enrichGO')
    artifacts[[paste0("clusterprofiler_", ont)]] <<- basename(out_tsv)
    artifacts[[paste0("clusterprofiler_", ont, "_dotplot")]] <<- pp$dot
    artifacts[[paste0("clusterprofiler_", ont, "_barplot")]] <<- pp$bar

    # Optional: cnetplot / enrichment map (clusterProfiler only)
    if (make_cnetplot || make_emapplot) {
      if (requireNamespace("enrichplot", quietly=TRUE) && requireNamespace("ggplot2", quietly=TRUE)) {
        suppressWarnings(suppressMessages({
          library(enrichplot)
          library(ggplot2)
        }))

        if (make_cnetplot) {
          p1 <- tryCatch({
            enrichplot::cnetplot(ego, showCategory=plot_top_n, circular=FALSE, colorEdge=TRUE)
          }, error=function(e) NULL)
          if (!is.null(p1)) {
            fn <- file.path(out_dir, sprintf("clusterprofiler_cnetplot_%s.png", ont))
            if (safe_ggsave(fn, p1)) {
              artifacts[[paste0("clusterprofiler_", ont, "_cnetplot")]] <<- basename(fn)
            }
          }
        }

        if (make_emapplot) {
          p2 <- tryCatch({
            ego2 <- enrichplot::pairwise_termsim(ego)
            enrichplot::emapplot(ego2, showCategory=plot_top_n)
          }, error=function(e) NULL)
          if (!is.null(p2)) {
            fn <- file.path(out_dir, sprintf("clusterprofiler_emapplot_%s.png", ont))
            if (safe_ggsave(fn, p2)) {
              artifacts[[paste0("clusterprofiler_", ont, "_emapplot")]] <<- basename(fn)
            }
          }
        }
      } else {
        cat("[rnaseq_enrich] enrichplot/ggplot2 not available; skip cnetplot/emapplot\n")
      }
    }
  }
}

do_goseq <- function() {
  if (!requireNamespace("goseq", quietly=TRUE)) fail("goseq not installed")
  # gene length table is strongly recommended
  gl <- read_gene_length_tsv(length_tsv, gene_col=length_gene_col, len_col=length_len_col)
  if (is.null(gl) || nrow(gl) < 10) {
    fail("goseq requires gene_length_tsv with sufficient rows (gene,length)")
  }

  # Limit to universe genes that have length
  gl <- gl[gl$gene %in% universe_genes, , drop=FALSE]
  if (nrow(gl) < 50) fail("Too few universe genes with lengths for goseq")

  # DE indicator vector (named)
  de <- as.integer(gl$gene %in% selected_genes)
  names(de) <- gl$gene

  suppressWarnings(suppressMessages({
    library(goseq)
  }))

  pwf <- tryCatch({
    goseq::nullp(de, bias.data=gl$length)
  }, error=function(e) {
    fail(paste0("goseq nullp failed: ", e$message))
  })

  # PWF diagnostic plot
  pwf_png <- file.path(out_dir, 'goseq_pwf.png')
  safe_png(pwf_png, function(){
    try(plot(pwf, main='goseq PWF (length bias)', xlab='Gene length', ylab='P(DE)'), silent=TRUE)
  })
  artifacts[['goseq_pwf']] <<- basename(pwf_png)

  # Build gene2go from OrgDb
  m <- build_gene2go_from_orgdb(keys=gl$gene, orgdb_pkg=orgdb_pkg, keytype=keytype)
  map <- m$mapping
  go_col <- m$col_go
  ont_col <- m$col_ont

  for (ont in ontologies) {
    map_ont <- map[map[[ont_col]] == ont, , drop=FALSE]
    gene2cat <- split(as.character(map_ont[[go_col]]), as.character(map_ont[[keytype]]))
    gene2cat <- lapply(gene2cat, function(v) unique(v[!is.na(v) & nzchar(v)]))

    # Ensure all genes appear (even empty) for goseq consistency
    # goseq will ignore those without categories

    res <- tryCatch({
      goseq::goseq(pwf, gene2cat=gene2cat)
    }, error=function(e) {
      fail(paste0("goseq goseq() failed (", ont, "): ", e$message))
    })

    # Adjust p-values
    if (!is.null(res$over_represented_pvalue)) {
      res$over_represented_padj <- p.adjust(res$over_represented_pvalue, method=p_adjust_method)
    }

    out_tsv <- file.path(out_dir, sprintf("goseq_%s.tsv", ont))
    data.table::fwrite(res, out_tsv, sep="\t")
    # plots
    std <- standardize_goseq(res, total_de=sum(de, na.rm=TRUE), total_uni=length(gl$gene))
    pp <- make_dot_bar(std, prefix='goseq', ont=ont, out_dir=out_dir, label='goseq enrichment')
    artifacts[[paste0("goseq_", ont)]] <<- basename(out_tsv)
    artifacts[[paste0("goseq_", ont, "_dotplot")]] <<- pp$dot
    artifacts[[paste0("goseq_", ont, "_barplot")]] <<- pp$bar
  }
}

do_topgo <- function() {
  if (!requireNamespace("topGO", quietly=TRUE)) fail("topGO not installed")

  # Use universe genes; build gene->GO mapping
  m <- build_gene2go_from_orgdb(keys=universe_genes, orgdb_pkg=orgdb_pkg, keytype=keytype)
  gene2go <- m$gene2go

  # All genes factor (0/1)
  geneList <- as.integer(universe_genes %in% selected_genes)
  names(geneList) <- universe_genes

  suppressWarnings(suppressMessages({
    library(topGO)
  }))

  for (ont in ontologies) {
    GOdata <- tryCatch({
      new("topGOdata",
          ontology=ont,
          allGenes=geneList,
          geneSel=function(x) x == 1,
          annot=topGO::annFUN.gene2GO,
          gene2GO=gene2go,
          nodeSize=10)
    }, error=function(e) {
      fail(paste0("topGOdata init failed (", ont, "): ", e$message))
    })

    test <- tryCatch({
      topGO::runTest(GOdata, algorithm=topgo_algorithm, statistic=topgo_statistic)
    }, error=function(e) {
      fail(paste0("topGO runTest failed (", ont, "): ", e$message))
    })

    tab <- tryCatch({
      topGO::GenTable(GOdata, pvalue=test, orderBy="pvalue", ranksOf="pvalue", topNodes=200)
    }, error=function(e) {
      fail(paste0("topGO GenTable failed (", ont, "): ", e$message))
    })

    out_tsv <- file.path(out_dir, sprintf("topgo_%s.tsv", ont))
    data.table::fwrite(tab, out_tsv, sep="\t")
    # plots
    std <- standardize_topgo(tab)
    pp <- make_dot_bar(std, prefix='topgo', ont=ont, out_dir=out_dir, label='topGO enrichment')
    artifacts[[paste0("topgo_", ont)]] <<- basename(out_tsv)
    artifacts[[paste0("topgo_", ont, "_dotplot")]] <<- pp$dot
    artifacts[[paste0("topgo_", ont, "_barplot")]] <<- pp$bar
  }
}

cat("[rnaseq_enrich] engine=", engine, "\n", sep="")
cat("[rnaseq_enrich] ontology=", paste(ontologies, collapse=","), "\n", sep="")

if (engine %in% c("clusterprofiler", "clusterprofiler_ora", "ora")) {
  do_clusterprofiler()
} else if (engine %in% c("goseq")) {
  do_goseq()
} else if (engine %in% c("topgo")) {
  do_topgo()
} else if (engine %in% c("all")) {
  do_clusterprofiler()
  do_goseq()
  do_topgo()
} else {
  fail(paste0("Unknown engine: ", engine))
}

# artifacts.json
art_path <- file.path(out_dir, "artifacts.json")
jsonlite::write_json(list(artifacts=artifacts), art_path, auto_unbox=TRUE, pretty=TRUE)
cat("[rnaseq_enrich] done\n")
