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

params <- tryCatch(jsonlite::fromJSON(params_path), error=function(e) list())
log("start")
log("params_path=", params_path)
log("out_dir=", out_dir)

if (!requireNamespace("mappoly", quietly=TRUE)) {
  msg <- "R package 'mappoly' is required but not installed in this environment. Please install it in the plugin R environment."
  write_error(msg); write_artifacts(list(module="poly_mappoly_import", error=msg)); quit(status=1)
}

geno_rds <- as.character(params$geno_rds %||% "")
ploidy <- as.numeric(params$ploidy %||% 4)
parent1_id <- as.character(params$parent1_id %||% "")
parent2_id <- as.character(params$parent2_id %||% "")
use_prob <- isTRUE(params$use_prob %||% TRUE) # default TRUE (we guarantee prob)

if (!nzchar(geno_rds) || !file.exists(geno_rds)) {
  msg <- paste0("geno_rds not found: ", geno_rds)
  write_error(msg); write_artifacts(list(module="poly_mappoly_import", error=msg)); quit(status=1)
}

out_files <- list()
if (grepl("*.csv$", geno_rds)) {
  if (isTRUE(use_prob)) {
    dat <- mappoly::read_geno_prob(file.in = geno_rds,
                                   prob.thres = 0.95,
                                   filter.non.conforming = F,
                                   elim.redundant = F,
                                   verbose = T)
    # dat <- mappoly::read_geno(file.in = geno_rds,
    #                           filter.non.conforming = F,
    #                           elim.redundant = F,
    #                           verbose = T)
  } else {
    dat <- mappoly::read_geno_csv(file.in = geno_rds,
                                  ploidy = ploidy,
                                  filter.non.conforming = F,
                                  elim.redundant = F,
                                  verbose = T)
  }
  dat_rds <- file.path(out_dir, "mappoly_data.rds")
  saveRDS(dat, dat_rds)
  out_files$mappoly_data_rds <- basename(dat_rds)
  
  markers <- dat$n.mrk
  samples <- dat$n.ind
  offspring <- dat$n.ind - 2
  
  mappoly_csv <- file.path(out_dir, "mappoly_dosage.csv")
  out_files$mappoly_dosage_csv <- basename(mappoly_csv)
} else if (grepl("*.vcf(|.tar.gz)$", geno_rds, perl = T)) {
  dat <- mappoly::read_vcf(file.in = geno_rds,
                           parent.1,
                           parent.2,
                           ploidy = ploidy,
                           filter.non.conforming = T,
                           thresh.line = 0.05,
                           min.gt.depth = 0,
                           min.av.depth = 0,
                           max.missing = 1,
                           elim.redundant = T,
                           verbose = T,
                           read.geno.prob = F,
                           prob.thres = 0.95)
  
  dat_rds <- file.path(out_dir, "mappoly_data.rds")
  saveRDS(dat, dat_rds)
  out_files$mappoly_data_rds <- basename(dat_rds)
  
  markers <- dat$n.mrk
  samples <- dat$n.ind
  offspring <- dat$n.ind - 2
  
  mappoly_csv <- file.path(out_dir, "mappoly_dosage.csv")
  out_files$mappoly_dosage_csv <- basename(mappoly_csv)
} else {
  geno <- readRDS(geno_rds)
  ploidy <- as.integer(geno$ploidy)
  if (is.na(ploidy)) stop("geno.rds missing ploidy")
  
  samples_df <- geno$samples
  markers_df <- geno$markers
  if (is.null(samples_df) || is.null(markers_df)) stop("geno.rds missing samples/markers")
  
  samples <- as.character(samples_df$sample_id %||% samples_df$sample %||% samples_df$ind)
  markers <- as.character(markers_df$marker_id %||% markers_df$marker)
  if (anyNA(samples) || anyNA(markers)) stop("failed to infer sample_id or marker_id")
  
  # infer parents from role if not provided
  if (!nzchar(parent1_id) && !is.null(samples_df$role)) {
    x <- samples_df$sample_id[samples_df$role == "parent1"]
    if (length(x) >= 1) parent1_id <- as.character(x[[1]])
  }
  if (!nzchar(parent2_id) && !is.null(samples_df$role)) {
    x <- samples_df$sample_id[samples_df$role == "parent2"]
    if (length(x) >= 1) parent2_id <- as.character(x[[1]])
  }
  
  if (!nzchar(parent1_id) || !nzchar(parent2_id)) {
    log("[WARN] parent1_id/parent2_id not provided and not found in samples$role. read_geno_prob requires parental dosages; proceeding but parents may be missing.")
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
  chr <- get_chr(markers_df); pos <- get_pos(markers_df)
  if (length(chr) != length(markers)) chr <- rep(NA, length(markers))
  if (length(pos) != length(markers)) pos <- rep(NA, length(markers))
  
  # Build MAPpoly input file(s)
  dosage <- geno$dosage
  prob <- geno$prob
  
  parents_present <- intersect(c(parent1_id, parent2_id), colnames(dosage))
  offspring <- setdiff(colnames(dosage), parents_present)
  
  get_d <- function(sid){
    if (!nzchar(sid) || !(sid %in% colnames(dosage))) return(rep(NA_integer_, nrow(dosage)))
    as.integer(dosage[, sid])
  }
  dP <- get_d(parent1_id); dQ <- get_d(parent2_id)
  
  out_files <- list()
  
  if (isTRUE(use_prob) && !is.null(prob) && length(dim(prob)) == 3 && dim(prob)[3] == (ploidy+1)) {
    # write mappoly_prob.txt for read_geno_prob
    mappoly_prob <- file.path(out_dir, "mappoly_prob.txt")
    con <- file(mappoly_prob, open="wt")
    on.exit(try(close(con), silent=TRUE), add=TRUE)
  
    writeLines(paste("ploidy", ploidy), con)
    writeLines(paste("n.ind", length(offspring)), con)
    writeLines(paste("n.mrk", length(markers)), con)
    writeLines(paste("mrk.names", paste(markers, collapse=" ")), con)
    writeLines(paste("ind.names", paste(offspring, collapse=" ")), con)
    writeLines(paste("dosageP", paste(ifelse(is.na(dP), ploidy+1L, dP), collapse=" ")), con)
    writeLines(paste("dosageQ", paste(ifelse(is.na(dQ), ploidy+1L, dQ), collapse=" ")), con)
    writeLines(paste("seq", paste(ifelse(is.na(chr), "NA", chr), collapse=" ")), con)
    writeLines(paste("seqpos", paste(ifelse(is.na(pos), "NA", pos), collapse=" ")), con)
    writeLines("nphen 0", con)
    writeLines("", con)
    writeLines("", con)
  
    dose_levels <- 0:ploidy
    hdr <- paste(c("marker","offspring", paste0("P", dose_levels)), collapse="\t")
    writeLines(hdr, con)
  
    # align dimnames
    mi <- if (!is.null(dimnames(prob)[[1]])) match(markers, dimnames(prob)[[1]]) else seq_along(markers)
    si <- if (!is.null(dimnames(prob)[[2]])) match(offspring, dimnames(prob)[[2]]) else match(offspring, colnames(dosage))
    
    for (ii in seq_along(markers)) {
      i_m <- mi[ii]
      if (is.na(i_m)) next
      m <- markers[ii]
      for (jj in seq_along(offspring)) {
        i_s <- si[jj]
        if (is.na(i_s)) next
        s <- offspring[jj]
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
  
    # NOTE: mappoly versions differ in read_geno_prob() signature.
    # In some versions, 'ploidy' is read from the header of the probability file
    # and passing ploidy as an argument raises "unused argument".
    dat <- mappoly::read_geno_prob(file.in = mappoly_prob)
    dat_rds <- file.path(out_dir, "mappoly_data.rds")
    saveRDS(dat, dat_rds)
    out_files$mappoly_data_rds <- basename(dat_rds)
  
  } else {
    # fallback dosage CSV
    mappoly_csv <- file.path(out_dir, "mappoly_dosage.csv")
    miss_code <- ploidy + 1L
    off_mat <- if (length(offspring) > 0) as.matrix(dosage[, offspring, drop=FALSE]) else matrix(nrow=nrow(dosage), ncol=0)
    
    off_mat <- apply(off_mat, 2, function(x){x <- as.integer(x); x[is.na(x)] <- miss_code; x})
    
    if (!is.matrix(off_mat)) off_mat <- matrix(off_mat, nrow=nrow(dosage))
    df <- data.frame(marker=markers, dosageP1=ifelse(is.na(dP), miss_code, dP), dosageP2=ifelse(is.na(dQ), miss_code, dQ), chrom=chr, pos=pos, stringsAsFactors=FALSE)
    if (ncol(off_mat) > 0) {colnames(off_mat) <- offspring; df <- cbind(df, as.data.frame(off_mat, check.names=FALSE))}
    utils::write.csv(df, mappoly_csv, quote=FALSE, row.names=FALSE)
    out_files$mappoly_dosage_csv <- basename(mappoly_csv)
    
    #dat <- mappoly::read_geno_csv(file.in = "/media/soba/Noc4/GenomicExplorer/P0210/sim/simulate_dataset_20260210_184738/dataset/polyploid/truth_dosage_merge_marker.csv", ploidy = ploidy)
    dat_rds <- file.path(out_dir, "mappoly_data.rds")
    saveRDS(dat, dat_rds)
    out_files$mappoly_data_rds <- basename(dat_rds)
  }
}

# summary
summary_tsv <- file.path(out_dir, "mappoly_import_summary.tsv")
utils::write.table(data.frame(
  key=c("ploidy","n_markers","n_individuals","parent1_id","parent2_id","n_offspring"),
  value=c(ploidy, length(markers), length(samples), parent1_id, parent2_id, length(offspring)),
  stringsAsFactors=FALSE
), summary_tsv, sep="\t", quote=FALSE, row.names=FALSE)
out_files$mappoly_import_summary_tsv <- basename(summary_tsv)

write_artifacts(c(list(module="poly_mappoly_import"), out_files))
log("done")
