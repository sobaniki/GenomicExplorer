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
log <- function(...) {
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

`%||%` <- function(a, b) if (!is.null(a)) a else b

safe_write_tsv <- function(df, path){
  utils::write.table(df, path, sep="\t", quote=FALSE, row.names=FALSE)
}

params <- tryCatch(jsonlite::fromJSON(params_path), error=function(e) list())
log("start")
log("params_path=", params_path)
log("out_dir=", out_dir)

geno_rds <- as.character(params$geno_rds %||% "")
ploidy_in <- as.integer(params$ploidy %||% NA)
parent1_id <- as.character(params$parent1_id %||% "")
parent2_id <- as.character(params$parent2_id %||% "")
engine <- as.character(params$engine %||% "both")  # mappoly|polymapr|both
write_prob <- isTRUE(params$write_prob %||% TRUE)

if (!nzchar(geno_rds) || !file.exists(geno_rds)) {
  msg <- paste0("geno_rds not found: ", geno_rds)
  write_error(msg)
  write_artifacts(list(module="poly_polymapr_map", error=msg))
  quit(status=1)
}

geno <- readRDS(geno_rds)
if (is.null(geno$ploidy)) stop("geno.rds missing 'ploidy'")
ploidy <- as.integer(geno$ploidy)
if (!is.na(ploidy_in) && ploidy_in != ploidy) {
  log("[WARN] params$ploidy (", ploidy_in, ") != geno$ploidy (", ploidy, "); using geno$ploidy")
}

samples_df <- geno$samples
markers_df <- geno$markers
dosage <- geno$dosage
prob <- geno$prob

if (is.null(samples_df) || is.null(markers_df) || is.null(dosage)) stop("geno.rds missing samples/markers/dosage")

samples <- as.character(samples_df$sample_id %||% samples_df$sample %||% samples_df$ind)
if (is.null(samples) || anyNA(samples)) stop("could not infer sample_id from geno$samples")

markers <- as.character(markers_df$marker_id %||% markers_df$marker)
if (is.null(markers) || anyNA(markers)) stop("could not infer marker_id from geno$markers")

if (!all(rownames(dosage) %in% markers) && !all(markers %in% rownames(dosage))) {
  # best-effort: assume rownames(dosage) are markers
  markers <- rownames(dosage)
}
if (!all(colnames(dosage) %in% samples) && !all(samples %in% colnames(dosage))) {
  samples <- colnames(dosage)
}

# infer parents from roles if not provided
if (!nzchar(parent1_id) && !is.null(samples_df$role)) {
  x <- samples_df$sample_id[samples_df$role == "parent1"]
  if (length(x) >= 1) parent1_id <- as.character(x[[1]])
}
if (!nzchar(parent2_id) && !is.null(samples_df$role)) {
  x <- samples_df$sample_id[samples_df$role == "parent2"]
  if (length(x) >= 1) parent2_id <- as.character(x[[1]])
}

parents_present <- intersect(c(parent1_id, parent2_id), samples)
offspring <- setdiff(samples, parents_present)

if (length(parents_present) < 2) {
  log("[WARN] Could not confirm both parents in geno.rds. parent1_id=", parent1_id, ", parent2_id=", parent2_id)
}

# marker positions (optional)
get_chr <- function(df){
  for (nm in c("chr","chrom","chromosome","seq")) if (nm %in% colnames(df)) return(df[[nm]])
  rep(NA, nrow(df))
}
get_pos <- function(df){
  for (nm in c("pos","position","bp","seqpos")) if (nm %in% colnames(df)) return(df[[nm]])
  rep(NA, nrow(df))
}

chr <- get_chr(markers_df)
pos <- get_pos(markers_df)
if (length(chr) != length(markers)) chr <- rep(NA, length(markers))
if (length(pos) != length(markers)) pos <- rep(NA, length(markers))

summary_tsv <- file.path(out_dir, "linkage_prep_summary.tsv")
pkg_tsv <- file.path(out_dir, "package_status.tsv")

have_mappoly <- requireNamespace("mappoly", quietly=TRUE)
have_polymapr <- requireNamespace("polymapR", quietly=TRUE)
safe_write_tsv(
  data.frame(
    package=c("mappoly","polymapR"),
    installed=c(have_mappoly, have_polymapr),
    stringsAsFactors=FALSE
  ),
  pkg_tsv
)

safe_write_tsv(
  data.frame(
    key=c("ploidy","n_markers","n_samples","n_parents_present","n_offspring"),
    value=c(ploidy, length(markers), length(samples), length(parents_present), length(offspring)),
    stringsAsFactors=FALSE
  ),
  summary_tsv
)

out_files <- list(
  linkage_prep_summary_tsv=basename(summary_tsv),
  package_status_tsv=basename(pkg_tsv)
)

# -----------------------------
# MAPpoly export
# -----------------------------
if (engine %in% c("mappoly","both")) {
  # dosage CSV for mappoly::read_geno_csv
  mappoly_csv <- file.path(out_dir, "mappoly_dosage.csv")
  miss_code <- ploidy + 1L
  get_d <- function(sid){
    if (!nzchar(sid) || !(sid %in% colnames(dosage))) return(rep(NA_integer_, nrow(dosage)))
    as.integer(dosage[, sid])
  }
  dP <- get_d(parent1_id)
  dQ <- get_d(parent2_id)
  dP[is.na(dP)] <- miss_code
  dQ[is.na(dQ)] <- miss_code

  off_mat <- if (length(offspring) > 0) as.matrix(dosage[, offspring, drop=FALSE]) else matrix(nrow=nrow(dosage), ncol=0)
  off_mat <- apply(off_mat, 2, function(x){
    x <- as.integer(x)
    x[is.na(x)] <- miss_code
    x
  })
  if (!is.matrix(off_mat)) off_mat <- matrix(off_mat, nrow=nrow(dosage))

  df <- data.frame(
    marker=markers,
    dosageP1=dP,
    dosageP2=dQ,
    chrom=chr,
    pos=pos,
    stringsAsFactors=FALSE
  )
  if (ncol(off_mat) > 0) {
    colnames(off_mat) <- offspring
    df <- cbind(df, as.data.frame(off_mat, check.names=FALSE))
  }
  utils::write.csv(df, mappoly_csv, quote=FALSE, row.names=FALSE)
  out_files$mappoly_dosage_csv <- basename(mappoly_csv)

  # probability TXT for mappoly::read_geno_prob (optional)
  if (isTRUE(write_prob) && !is.null(prob) && length(dim(prob)) == 3 && dim(prob)[3] == (ploidy+1)) {
    mappoly_prob <- file.path(out_dir, "mappoly_prob.txt")
    con <- file(mappoly_prob, open="wt")
    on.exit(try(close(con), silent=TRUE), add=TRUE)

    writeLines(paste("ploidy", ploidy), con)
    writeLines(paste("n.ind", length(offspring)), con)
    writeLines(paste("n.mrk", length(markers)), con)
    writeLines(paste("mrk.names", paste(markers, collapse=" ")), con)
    writeLines(paste("ind.names", paste(offspring, collapse=" ")), con)
    writeLines(paste("dosageP", paste(dP, collapse=" ")), con)
    writeLines(paste("dosageQ", paste(dQ, collapse=" ")), con)
    writeLines(paste("seq", paste(ifelse(is.na(chr), "NA", chr), collapse=" ")), con)
    writeLines(paste("seqpos", paste(ifelse(is.na(pos), "NA", pos), collapse=" ")), con)
    writeLines("nphen 0", con)
    writeLines("", con)   # line 11 skipped
    writeLines("", con)   # line 12 + nphen skipped

    # probability table
    # marker offspring P0..Pploidy
    # stream writing to avoid huge memory
    dose_levels <- 0:ploidy
    hdr <- paste(c("marker","offspring", paste0("P", dose_levels)), collapse="\t")
    writeLines(hdr, con)

    # align prob dims (markers x samples x (P+1))
    # best-effort reorder
    if (!is.null(dimnames(prob)[[1]])) {
      mi <- match(markers, dimnames(prob)[[1]])
    } else {
      mi <- seq_along(markers)
    }
    if (!is.null(dimnames(prob)[[2]])) {
      si <- match(offspring, dimnames(prob)[[2]])
    } else {
      si <- match(offspring, colnames(dosage))
    }
    for (ii in seq_along(markers)) {
      m <- markers[ii]
      i_m <- mi[ii]
      if (is.na(i_m)) next
      for (jj in seq_along(offspring)) {
        s <- offspring[jj]
        i_s <- si[jj]
        if (is.na(i_s)) next
        pv <- prob[i_m, i_s, ]
        if (all(is.na(pv))) {
          line <- paste(c(m, s, rep("NA", ploidy+1)), collapse="\t")
        } else {
          line <- paste(c(m, s, format(as.numeric(pv), scientific=FALSE)), collapse="\t")
        }
        writeLines(line, con)
      }
    }
    close(con)
    out_files$mappoly_prob_txt <- basename(mappoly_prob)

    # sanity read if package installed
    if (have_mappoly) {
      log("[mappoly] sanity check: read_geno_csv / read_geno_prob")
      try({
        md <- mappoly::read_geno_csv(mappoly_csv, ploidy=ploidy)
        saveRDS(md, file.path(out_dir, "mappoly_data_dosage.rds"))
        out_files$mappoly_data_dosage_rds <- "mappoly_data_dosage.rds"
      }, silent=TRUE)
      try({
        mp <- mappoly::read_geno_prob(mappoly_prob)
        saveRDS(mp, file.path(out_dir, "mappoly_data_prob.rds"))
        out_files$mappoly_data_prob_rds <- "mappoly_data_prob.rds"
      }, silent=TRUE)
    }
  }
}

# -----------------------------
# polymapR export
# -----------------------------
if (engine %in% c("polymapr","both")) {
  polymapr_csv <- file.path(out_dir, "polymapr_prob.csv")

  if (is.null(prob) || length(dim(prob)) != 3 || dim(prob)[3] != (ploidy+1)) {
    log("[WARN] geno$prob missing; polymapR probabilistic export skipped")
  } else {
    dose_levels <- 0:ploidy
    Pcols <- paste0("P", dose_levels)
    # marker name
    marker_name <- ifelse(!is.na(chr) & !is.na(pos), paste0(chr, "_", pos), markers)
    # MAP genotype
    map_dosage <- apply(prob, c(1,2), function(p){ if (all(is.na(p))) NA_integer_ else which.max(p) - 1L })

    # Build long table (can be large; write with streaming)
    con <- file(polymapr_csv, open="wt")
    on.exit(try(close(con), silent=TRUE), add=TRUE)
    hdr <- paste(c("marker","MarkerName","SampleName","ratio", Pcols, "maxgeno","maxP","geno"), collapse=",")
    writeLines(hdr, con)

    # align order
    mnames_prob <- dimnames(prob)[[1]]
    snames_prob <- dimnames(prob)[[2]]
    mi <- if (!is.null(mnames_prob)) match(markers, mnames_prob) else seq_along(markers)
    si <- if (!is.null(snames_prob)) match(samples, snames_prob) else seq_along(samples)

    for (ii in seq_along(markers)) {
      i_m <- mi[ii]
      if (is.na(i_m)) next
      for (jj in seq_along(samples)) {
        i_s <- si[jj]
        if (is.na(i_s)) next
        pv <- prob[i_m, i_s, ]
        if (all(is.na(pv))) next
        mx <- max(pv, na.rm=TRUE)
        gx <- which.max(pv) - 1L
        gmap <- map_dosage[i_m, i_s]
        line <- paste(
          c(ii,
            paste0('"', marker_name[ii], '"'),
            paste0('"', samples[jj], '"'),
            "NA",
            format(as.numeric(pv), scientific=FALSE),
            gx,
            format(mx, scientific=FALSE),
            ifelse(is.na(gmap), "NA", gmap)
          ),
          collapse="," )
        writeLines(line, con)
      }
    }
    close(con)
    out_files$polymapr_prob_csv <- basename(polymapr_csv)
  }
}

meta <- c(
  list(
    module="poly_polymapr_map",
    note="Prepared linkage mapping inputs (MAPpoly / polymapR).",
    geno_rds_in=basename(geno_rds),
    ploidy=ploidy,
    engine=engine,
    parent1_id=parent1_id,
    parent2_id=parent2_id,
    n_samples=length(samples),
    n_markers=length(markers)
  ),
  out_files,
  list(default_table=basename(summary_tsv))
)

write_artifacts(meta)
log("done")
