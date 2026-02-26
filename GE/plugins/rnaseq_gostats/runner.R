#!/usr/bin/env Rscript
# Standalone GOstats ORA
suppressWarnings(suppressMessages({
  library(jsonlite)
  library(data.table)
}))

fail <- function(msg, code=1) {
  cat("[rnaseq_gostats][ERROR] ", msg, "\n", sep="")
  quit(status=code)
}

as_bool <- function(x) {
  if (is.null(x)) return(FALSE)
  if (is.logical(x)) return(isTRUE(x))
  if (is.character(x)) return(tolower(x) %in% c("1","true","t","yes","y"))
  if (is.numeric(x)) return(x != 0)
  FALSE
}

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

standardize_gostats <- function(df, sel_n){
  if (is.null(df) || nrow(df) < 1) return(data.frame())
  term <- if ('Term' %in% colnames(df)) as.character(df$Term) else if ('GOTerm' %in% colnames(df)) as.character(df$GOTerm) else as.character(df[[1]])
  pval <- if ('Pvalue' %in% colnames(df)) suppressWarnings(as.numeric(df$Pvalue)) else if ('pvalue' %in% colnames(df)) suppressWarnings(as.numeric(df$pvalue)) else rep(NA_real_, nrow(df))
  padj <- p.adjust(pval, method='BH')
  cnt <- if ('Count' %in% colnames(df)) suppressWarnings(as.numeric(df$Count)) else if ('Significant' %in% colnames(df)) suppressWarnings(as.numeric(df$Significant)) else rep(NA_real_, nrow(df))
  size <- if ('Size' %in% colnames(df)) suppressWarnings(as.numeric(df$Size)) else rep(NA_real_, nrow(df))
  gr <- ifelse(is.finite(size) & size > 0, cnt/size, ifelse(sel_n > 0, cnt/sel_n, NA_real_))
  rich_factor <- gr
  out <- data.frame(term=term, count=cnt, gene_ratio=gr, bg_ratio=NA_real_, rich_factor=rich_factor, set_size=size, padj=padj, stringsAsFactors=FALSE)
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

make_dot_bar_gostats <- function(std, ont, out_dir){
  dot_png <- file.path(out_dir, sprintf('gostats_dotplot_%s.png', ont))
  bar_png <- file.path(out_dir, sprintf('gostats_barplot_%s.png', ont))
  make_dotplot(std, dot_png, title=paste0('GOstats ORA ', ont))
  make_barplot(std, bar_png, title=paste0('GOstats ORA ', ont))
  list(dot=basename(dot_png), bar=basename(bar_png))
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 4 || args[1] != "--params" || args[3] != "--out") {
  cat("Usage: runner.R --params params.json --out out_dir\n")
  quit(status=2)
}

params_path <- args[2]
out_dir <- args[4]
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

cat("[rnaseq_gostats] start\n")
cat("[rnaseq_gostats] params_path=", params_path, "\n", sep="")
cat("[rnaseq_gostats] out_dir=", out_dir, "\n", sep="")

params <- tryCatch(jsonlite::fromJSON(params_path), error=function(e) fail(paste("Failed to parse params.json:", e$message)))
`%||%` <- function(a, b) { if (is.null(a) || length(a)==0 || (is.character(a)&&a=="")) b else a }

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
    cat("[rnaseq_gostats] .libPaths prepended: ", paste(expanded, collapse=", "), "\n", sep="")
  }
}


apply_extra_libs(params$r_libs %||% params$r_lib %||% params$rlibs %||% params$extra_libs)

gene_list_path <- params$gene_list %||% params$gene_list_path %||% params$deg_list
universe_path <- params$universe_list %||% params$universe %||% params$universe_list_path

orgdb_pkg <- params$orgdb_package %||% params$orgdb %||% "org.Hs.eg.db"
keytype <- params$keytype %||% "ENTREZID"
ontology <- toupper(params$ontology %||% "BP")
conditional <- as_bool(params$conditional %||% FALSE)

# Optional: pvalue cutoff for reporting (GOstats summary already contains pvalue/odds)
report_p_cutoff <- as.numeric(params$report_p_cutoff %||% 1.0)

# ---------- plot settings (GUI-controlled) ----------
plot_top_n <- suppressWarnings(as.integer(params$plot_top_n %||% params$top_n %||% 20))
if (!is.finite(plot_top_n) || plot_top_n < 1) plot_top_n <- 20L

label_max_chars <- suppressWarnings(as.integer(params$plot_label_max_chars %||% params$label_max_chars %||% 60))
if (!is.finite(label_max_chars) || label_max_chars < 0) label_max_chars <- 60L

plot_x_metric <- tolower(as.character(params$plot_x_metric %||% "gene_ratio"))
plot_size_metric <- tolower(as.character(params$plot_size_metric %||% "count"))
plot_color_metric <- tolower(as.character(params$plot_color_metric %||% "mlog10"))
plot_bar_metric <- tolower(as.character(params$plot_bar_metric %||% "mlog10"))

if (is.null(gene_list_path) || gene_list_path=="" || !file.exists(gene_list_path)) {
  fail("gene_list file not found")
}

genes_sel <- read_lines_nonempty(gene_list_path)
if (length(genes_sel) < 1) fail("gene_list is empty")

genes_uni <- read_lines_nonempty(universe_path)
# If universe not provided, fall back to selected genes (warn) so pipeline still produces output.
if (length(genes_uni) < 1) {
  cat("[rnaseq_gostats] universe_list not provided; using gene_list as universe (not recommended)\n")
  genes_uni <- genes_sel
}

suppressWarnings(suppressMessages({
  if (!requireNamespace("GOstats", quietly=TRUE)) fail("GOstats not installed")
  if (!requireNamespace("AnnotationDbi", quietly=TRUE)) fail("AnnotationDbi not installed")
  if (!requireNamespace("GO.db", quietly=TRUE)) fail("GO.db not installed")
  library(GOstats)
  library(AnnotationDbi)
  library(GO.db)
}))

# Load OrgDb
if (!requireNamespace(orgdb_pkg, quietly=TRUE)) {
  fail(paste0("OrgDb package not installed: ", orgdb_pkg))
}
suppressWarnings(suppressMessages(library(orgdb_pkg, character.only=TRUE)))
orgdb <- get(orgdb_pkg)

# Map to ENTREZID if needed
map_to_entrez <- function(keys) {
  keys <- as.character(keys)
  if (toupper(keytype) == "ENTREZID") return(unique(keys))
  m <- AnnotationDbi::select(orgdb, keys=keys, keytype=keytype, columns=c("ENTREZID"))
  m <- m[!is.na(m$ENTREZID), , drop=FALSE]
  unique(as.character(m$ENTREZID))
}

sel_entrez <- map_to_entrez(genes_sel)
uni_entrez <- map_to_entrez(genes_uni)

if (length(sel_entrez) < 5) {
  fail("Too few selected genes for GOstats (after mapping)")
}
if (length(uni_entrez) < 50) {
  cat("[rnaseq_gostats] WARNING: universe is small after mapping (", length(uni_entrez), ")\n", sep="")
}

out_tsv <- file.path(out_dir, sprintf("go_enrichment_%s.tsv", ontology))

params_go <- new("GOHyperGParams",
                 geneIds=sel_entrez,
                 universeGeneIds=uni_entrez,
                 annotation=orgdb_pkg,
                 ontology=ontology,
                 pvalueCutoff=1.0,
                 conditional=conditional,
                 testDirection="over")

hg <- hyperGTest(params_go)
sm <- summary(hg)
if (!is.null(sm) && nrow(sm) > 0 && is.finite(report_p_cutoff) && report_p_cutoff < 1.0) {
  sm <- sm[sm$Pvalue <= report_p_cutoff, , drop=FALSE]
}
write.table(sm, out_tsv, sep="\t", quote=FALSE, row.names=FALSE)

# plots
std_plot <- standardize_gostats(sm, sel_n=length(sel_entrez))
pp_plot <- make_dot_bar_gostats(std_plot, ontology, out_dir)

# artifacts
art <- list(
  go_tsv = out_tsv,
  gene_list = gene_list_path,
  universe_list = if (!is.null(universe_path) && universe_path != "") universe_path else NULL,
  ontology = ontology,
  orgdb = orgdb_pkg,
  keytype = keytype
)
# add plot paths (basenames)
if (exists('pp_plot')) {
  art$dotplot_png <- pp_plot$dot
  art$barplot_png <- pp_plot$bar
}
writeLines(jsonlite::toJSON(art, auto_unbox=TRUE, pretty=TRUE), file.path(out_dir, "artifacts.json"))

cat("[rnaseq_gostats] done\n")
