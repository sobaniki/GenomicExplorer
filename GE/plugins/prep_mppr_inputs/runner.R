#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(jsonlite)
  library(data.table)
})

args <- commandArgs(trailingOnly=TRUE)
get_arg <- function(flag, default=NULL) {
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) return(args[i+1])
  default
}

params_path <- get_arg("--params")
out_dir <- get_arg("--out")
if (is.null(params_path) || is.null(out_dir)) {
  stop("Usage: runner.R --params params.json --out out_dir")
}
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

p <- fromJSON(params_path)
`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && nzchar(as.character(a)[1])) a else b

as_bool <- function(x, default=TRUE) {
  if (is.null(x) || length(x) == 0) return(default)
  if (is.logical(x)) return(isTRUE(x[1]))
  s <- tolower(as.character(x[1]))
  if (s %in% c('1','true','t','yes','y')) return(TRUE)
  if (s %in% c('0','false','f','no','n')) return(FALSE)
  default
}

# ---- params ----
plink_prefix <- as.character(p$plink_prefix %||% "")
phenotype_path <- as.character(p$phenotype_path %||% p$phenotype_tsv %||% "")
phenotype_sheet <- as.character(p$phenotype_sheet %||% "")
phenotype_id_col <- as.character(p$phenotype_id_col %||% "")
cross_col <- as.character(p$cross_col %||% "Family_Inbred_Name")
panel_col <- as.character(p$panel_col %||% "")
panel_values_raw <- as.character(p$panel_values %||% "")
iid_key_delim <- as.character(p$iid_key_delim %||% ":")
plink2_bin <- as.character(p$plink2_bin %||% "plink2")
out_override <- as.character(p$out_dir_override %||% "")

# extras
parent_select_mode <- as.character(p$parent_select_mode %||% 'best_callrate')
genetic_map_path <- as.character(p$genetic_map_path %||% '')
genetic_map_sheet <- as.character(p$genetic_map_sheet %||% '')
map_merge_mode <- as.character(p$map_merge_mode %||% 'marker_then_chr_bp')
interpolate_missing_cm <- as_bool(p$interpolate_missing_cm, TRUE)

# normalize prefix
if (grepl("\\.bed$", plink_prefix, ignore.case=TRUE)) plink_prefix <- sub("\\.bed$", "", plink_prefix, ignore.case=TRUE)
if (grepl("\\.bim$", plink_prefix, ignore.case=TRUE)) plink_prefix <- sub("\\.bim$", "", plink_prefix, ignore.case=TRUE)
if (grepl("\\.fam$", plink_prefix, ignore.case=TRUE)) plink_prefix <- sub("\\.fam$", "", plink_prefix, ignore.case=TRUE)

if (!nzchar(plink_prefix)) stop("plink_prefix is required")
if (!nzchar(phenotype_path)) stop("phenotype_path is required")

bed <- paste0(plink_prefix, ".bed")
bim <- paste0(plink_prefix, ".bim")
fam <- paste0(plink_prefix, ".fam")
if (!file.exists(bed) || !file.exists(bim) || !file.exists(fam)) {
  stop("PLINK files not found (.bed/.bim/.fam): ", plink_prefix)
}

base_out <- if (nzchar(out_override)) out_override else file.path(out_dir, "mppr_inputs")
dir.create(base_out, recursive=TRUE, showWarnings=FALSE)

log_file <- file.path(base_out, "prep.log")
sink(log_file, split=TRUE)

cat("[prep_mppr_inputs] start\n")


run_plink2_safe <- function(args) {
      cat("[prep_mppr_inputs] plink2 ", plink2_bin, " ", paste(args, collapse=" "), "\n", sep="")
      out <- tryCatch({
        tmp_out <- tempfile("plink2_out_")
        tmp_err <- tempfile("plink2_err_")
        rc <- system2(plink2_bin, args=args, stdout=tmp_out, stderr=tmp_err)
        so <- paste(readLines(tmp_out, warn=FALSE), collapse="\n")
        se <- paste(readLines(tmp_err, warn=FALSE), collapse="\n")
        file.remove(tmp_out)
        file.remove(tmp_err)
        list(code=rc, out=so, err=se)
      }, error=function(e) {
        list(code=1, out="", err=e$message)
      })
      if (is.null(out$code) || out$code != 0) {
        tail <- paste0(out$out, "\n", out$err)
        if (nchar(tail) > 2000) tail <- substr(tail, nchar(tail)-2000, nchar(tail))
        stop("plink2 failed: ", tail)
      }
      TRUE
    }

cat("[prep_mppr_inputs] plink_prefix=", plink_prefix, "\n")
cat("[prep_mppr_inputs] phenotype_path=", phenotype_path, "\n")
cat("[prep_mppr_inputs] base_out=", base_out, "\n")

# ---- read phenotype (TSV/CSV or XLSX) ----
read_phenotype_any <- function(path, sheet="") {
  path <- as.character(path)
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("xlsx", "xls")) {
    if (!requireNamespace("readxl", quietly=TRUE)) {
      cat("[prep_mppr_inputs] NOTE: installing readxl ...\n")
      install.packages("readxl", repos="https://cloud.r-project.org", quiet=TRUE)
    }
    sh <- sheet
    if (!nzchar(sh)) sh <- 1
    # Treat common missing value tokens as NA to avoid noisy parsing warnings
    df <- as.data.table(readxl::read_excel(path, sheet=sh, na=c("", "NA", "NaN", ".")))
    return(df)
  }
  # default: fread auto
  dt <- tryCatch({
    fread(path)
  }, error=function(e) {
    fread(path, sep=",", header=TRUE)
  })
  dt
}

ph <- read_phenotype_any(phenotype_path, phenotype_sheet)
if (nrow(ph) < 1) stop("Phenotype table is empty")

# Detect id column if not provided
if (!nzchar(phenotype_id_col) || !(phenotype_id_col %in% names(ph))) {
  for (cand in c("id","ID","sample","Sample","IID","Z_Num","Entry","Entry_Num")) {
    if (cand %in% names(ph)) { phenotype_id_col <- cand; break }
  }
}
if (!nzchar(phenotype_id_col) || !(phenotype_id_col %in% names(ph))) phenotype_id_col <- names(ph)[1]

if (!(cross_col %in% names(ph))) {
  # fallback heuristics
  for (cand in c("Family_Inbred_Name","family","Family","cross","Cross","pop","population")) {
    if (cand %in% names(ph)) { cross_col <- cand; break }
  }
}
if (!(cross_col %in% names(ph))) stop("cross_col not found in phenotype: ", cross_col)

# optional panel filter
panel_values <- character(0)
if (nzchar(panel_values_raw)) {
  parts <- unlist(strsplit(panel_values_raw, "[,;\\s]+"))
  parts <- parts[nzchar(parts)]
  panel_values <- unique(parts)
}
if (nzchar(panel_col) && length(panel_values) > 0 && (panel_col %in% names(ph))) {
  ph <- ph[get(panel_col) %in% panel_values]
}

if (nrow(ph) < 1) stop("No rows after panel filter")

# Keep unique rows by raw phenotype id
ph[[phenotype_id_col]] <- as.character(ph[[phenotype_id_col]])
ph[[cross_col]] <- as.character(ph[[cross_col]])
ph <- ph[!is.na(ph[[phenotype_id_col]]) & nzchar(ph[[phenotype_id_col]]), ]
if (nrow(ph) < 1) stop("Phenotype has no valid IDs")
ph <- ph[!duplicated(ph[[phenotype_id_col]]), ]

# ---- read fam and build key mapping ----
fam_dt <- fread(fam, header=FALSE)
if (ncol(fam_dt) < 2) stop("Invalid .fam")
setnames(fam_dt, c("FID","IID", paste0("V", 3:ncol(fam_dt))))
fam_dt[, IID := as.character(IID)]
fam_dt[, FID := as.character(FID)]

# compute base key from IID
get_base_key <- function(x, delim=":") {
  x <- as.character(x)
  if (!nzchar(delim)) return(x)
  # take substring before first delim
  i <- regexpr(delim, x, fixed=TRUE)
  out <- ifelse(i > 0, substr(x, 1, i-1), x)
  out
}

fam_dt[, base_key := get_base_key(IID, iid_key_delim)]

# Map phenotype IDs to actual IIDs.
# If phenotype IDs already match IIDs, keep as-is.
# Otherwise, try matching by base_key.
ph[, raw_id := as.character(get(phenotype_id_col))]
ph[, iid := raw_id]
setkey(fam_dt, IID)
ph[, has_exact := iid %in% fam_dt$IID]

exact_rate <- mean(ph$has_exact)
cat(sprintf("[prep_mppr_inputs] phenotype exact IID match rate: %.3f\n", exact_rate))

if (exact_rate < 0.5) {
  # attempt base-key mapping
  fam_key_map <- fam_dt[, .(IID, FID), by=base_key]
  setkey(fam_key_map, base_key)
  ph[, base_key := raw_id]
  ph <- fam_key_map[ph, on=.(base_key), nomatch=0]
  # fam_key_map[ph] keeps columns from ph plus IID/FID from map (as IID/FID)
  # rename
  setnames(ph, c("IID","FID"), c("iid","fid"))
  if (!("raw_id" %in% names(ph))) ph[, raw_id := base_key]
} else {
  ph[, fid := fam_dt[.(iid), FID]]
}

# drop rows without iid
ph <- ph[!is.na(iid) & nzchar(iid)]
if (nrow(ph) < 1) stop("No phenotype IDs could be mapped to PLINK IIDs")

offspring_iids <- unique(ph$iid)
cat("[prep_mppr_inputs] offspring n=", length(offspring_iids), "\n")

# ---- build cross_ind ----
# normalize cross strings
normalize_cross <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x
}
ph[, cross := normalize_cross(get(cross_col))]

cross_ind <- ph[, .(id=iid, cross=cross)]

# ---- parse parents per cross and map to parent IIDs ----
parse_cross_parents <- function(cross) {
  s <- gsub("\\s+", "", as.character(cross))
  # allow delimiters: x, X, ×
  parts <- unlist(strsplit(s, "[xX×]"))
  parts <- parts[nzchar(parts)]
  if (length(parts) >= 2) return(parts[1:2])
  character(0)
}

crosses <- unique(cross_ind$cross)
pp_list <- lapply(crosses, function(cr) {
  pr <- parse_cross_parents(cr)
  if (length(pr) != 2) {
    data.table(cross=cr, parent_key1=NA_character_, parent_key2=NA_character_)
  } else {
    data.table(cross=cr, parent_key1=pr[1], parent_key2=pr[2])
  }
})
pp <- rbindlist(pp_list, fill=TRUE)

bad_cross <- pp[is.na(parent_key1) | is.na(parent_key2)]
if (nrow(bad_cross) > 0) {
  fwrite(bad_cross, file.path(base_out, "bad_crosses.tsv"), sep="\t", quote=FALSE)
  stop("Some cross names could not be parsed into two parents. See bad_crosses.tsv")
}

parent_keys <- unique(c(pp$parent_key1, pp$parent_key2))


# map parent key -> parent IID (handle replicate IIDs per parent key)
escape_regex <- function(x) gsub("([\\.\\+\\*\\?\\^\\$\\(\\)\\[\\]\\{\\}\\|\\\\])", "\\\\\\1", x)

get_parent_candidates <- function(key) {
  key <- as.character(key)
  cand <- fam_dt[base_key == key, .(FID, IID)]
  if (nrow(cand) > 0) return(cand)
  idx <- grep(paste0("^", escape_regex(key)), fam_dt$IID, ignore.case=TRUE)
  if (length(idx) > 0) return(fam_dt[idx, .(FID, IID)])
  data.table(FID=character(0), IID=character(0))
}

parent_select_mode <- tolower(trimws(as.character(parent_select_mode)))
if (!(parent_select_mode %in% c("best_callrate","lexicographic","first"))) parent_select_mode <- "best_callrate"

parent_cand <- rbindlist(lapply(parent_keys, function(k) {
  cand <- get_parent_candidates(k)
  if (nrow(cand) < 1) return(data.table(parent_key=k, FID=NA_character_, IID=NA_character_))
  cand[, parent_key := as.character(k)]
  cand
}), fill=TRUE)

# Missing candidates?
missing_keys <- unique(parent_cand[is.na(IID) | !nzchar(IID), parent_key])
if (length(missing_keys) > 0) {
  fwrite(parent_cand, file.path(base_out, "parent_mapping_candidates.tsv"), sep="\t", quote=FALSE, na="NA")
  stop("Some parents could not be mapped to PLINK IIDs. See parent_mapping_candidates.tsv")
}

# Compute per-sample missingness for all candidate parent replicates (optional)
parent_cand[, miss := NA_real_]
if (parent_select_mode == "best_callrate") {
  cand_keep <- unique(parent_cand[, .(FID, IID)])
  keep_path2 <- file.path(base_out, "keep_parent_candidates.fid_iid.tsv")
  fwrite(cand_keep, keep_path2, sep="\t", col.names=FALSE, quote=FALSE)
  miss_pref <- file.path(base_out, "parent_candidates_missing")
  miss_dt <- NULL
  tryCatch({
    run_plink2_safe(c("--bfile", plink_prefix, "--keep", keep_path2, "--missing", "--out", miss_pref))
    smiss <- paste0(miss_pref, ".smiss")
    imiss <- paste0(miss_pref, ".imiss")
    if (file.exists(smiss)) miss_dt <- fread(smiss)
    if (is.null(miss_dt) && file.exists(imiss)) miss_dt <- fread(imiss)
  }, error=function(e) {
    cat("[prep_mppr_inputs] WARN: could not compute parent missingness: ", e$message, "\n")
    miss_dt <<- NULL
  })
  if (!is.null(miss_dt)) {
    fid_col <- grep("^#?fid$", names(miss_dt), ignore.case=TRUE, value=TRUE)[1]
    iid_col <- grep("^iid$", names(miss_dt), ignore.case=TRUE, value=TRUE)[1]
    miss_col <- grep("f_miss", names(miss_dt), ignore.case=TRUE, value=TRUE)[1]
    if (is.na(miss_col) || !nzchar(miss_col)) miss_col <- grep("miss", names(miss_dt), ignore.case=TRUE, value=TRUE)[1]
    if (!is.na(fid_col) && !is.na(iid_col) && !is.na(miss_col)) {
      setnames(miss_dt, c(fid_col, iid_col, miss_col), c("FID", "IID", "F_MISS"))
      miss_dt[, FID := as.character(FID)]
      miss_dt[, IID := as.character(IID)]
      miss_dt[, F_MISS := suppressWarnings(as.numeric(as.character(F_MISS)))]
      miss_dt <- unique(miss_dt, by=c("FID","IID"))
      setkey(miss_dt, FID, IID)
      # IMPORTANT: ensure we join using parent_cand's FID/IID, not miss_dt's columns.
      # Using miss_dt[.(FID,IID)] can accidentally capture miss_dt's own columns in NSE.
      parent_cand[, miss := miss_dt[.SD, on=.(FID, IID), x.F_MISS]]
    }
  }
}

# Choose one replicate per parent_key
pick_one <- function(dt) {
  dt <- dt[!is.na(IID) & nzchar(IID)]
  if (nrow(dt) < 1) return(data.table(parent_iid=NA_character_, parent_fid=NA_character_, miss=NA_real_, candidates=""))
  cand_str <- paste(dt$IID, collapse=",")
  if (parent_select_mode == "best_callrate") {
    # NOTE: dt is typically .SD (locked) when called from data.table.
    # Avoid modifying by reference on locked .SD by working on a copy.
    dt2 <- copy(dt)
    dt2[, miss2 := ifelse(is.na(miss), Inf, miss)]
    setorder(dt2, miss2, IID)
    best <- dt2[1]
    return(data.table(parent_iid=as.character(best$IID), parent_fid=as.character(best$FID), miss=best$miss, candidates=cand_str))
  }
  if (parent_select_mode == "lexicographic") {
    setorder(dt, IID)
    best <- dt[1]
    return(data.table(parent_iid=as.character(best$IID), parent_fid=as.character(best$FID), miss=best$miss, candidates=cand_str))
  }
  # first (fam order)
  best <- dt[1]
  data.table(parent_iid=as.character(best$IID), parent_fid=as.character(best$FID), miss=best$miss, candidates=cand_str)
}

parent_choice <- parent_cand[, pick_one(.SD), by=.(parent_key)]

# Write replicate choice report
fwrite(parent_choice, file.path(base_out, "parent_replicate_choice.tsv"), sep="\t", quote=FALSE, na="NA")

parent_map <- parent_choice[, .(parent_key=as.character(parent_key),
                               parent_iid=as.character(parent_iid),
                               parent_fid=as.character(parent_fid))]

missing_par <- parent_map[is.na(parent_iid) | !nzchar(parent_iid)]
if (nrow(missing_par) > 0) {
  fwrite(parent_map, file.path(base_out, "parent_mapping.tsv"), sep="\t", quote=FALSE, na="NA")
  stop("Some parents could not be mapped to PLINK IIDs. See parent_mapping.tsv")
}

# Build par_per_cross with IIDs
pm <- setNames(parent_map$parent_iid, parent_map$parent_key)
pp_out <- pp[, .(cross=cross, parent1=pm[parent_key1], parent2=pm[parent_key2])]

# ---- extract subset plink (offspring + parents) ----
keep_dt <- unique(rbind(
  ph[, .(FID=fid, IID=iid)],
  parent_map[, .(FID=parent_fid, IID=parent_iid)]
))
keep_path <- file.path(base_out, "keep.fid_iid.tsv")
fwrite(keep_dt, keep_path, sep="\t", col.names=FALSE, quote=FALSE)

subset_pref <- file.path(base_out, "subset")

# run plink2
run_plink2 <- function(args) {
  cat("[prep_mppr_inputs] plink2 ", plink2_bin, " ", paste(args, collapse=" "), "\n", sep="")
  out <- tryCatch({
    tmp_out <- tempfile("plink2_out_")
    tmp_err <- tempfile("plink2_err_")
    rc <- system2(plink2_bin, args=args, stdout=tmp_out, stderr=tmp_err)
    so <- paste(readLines(tmp_out, warn=FALSE), collapse="\n")
    se <- paste(readLines(tmp_err, warn=FALSE), collapse="\n")
    file.remove(tmp_out)
    file.remove(tmp_err)
    list(code=rc, out=so, err=se)
  }, error=function(e) {
    list(code=1, out="", err=e$message)
  })
  if (is.null(out$code) || out$code != 0) {
    tail <- paste0(out$out, "\n", out$err)
    if (nchar(tail) > 2000) tail <- substr(tail, nchar(tail)-2000, nchar(tail))
    stop("plink2 failed: ", tail)
  }
  TRUE
}

run_plink2(c("--bfile", plink_prefix, "--keep", keep_path, "--make-bed", "--out", subset_pref))

# ---- convert subset plink -> genotype.tsv + marker_map.tsv ----
if (!requireNamespace("gaston", quietly=TRUE)) {
  stop("R package 'gaston' is required (used elsewhere in GE). Please install it.")
}
suppressPackageStartupMessages(library(gaston))

bm <- read.bed.matrix(paste0(subset_pref, ".bed"))
ids <- as.character(bm@ped$id)
markers_raw <- as.character(bm@snps$id)
markers <- make.unique(markers_raw)
X <- as.matrix(bm)
X <- suppressWarnings(matrix(as.numeric(X), nrow=nrow(X), ncol=ncol(X), dimnames=dimnames(X)))


# marker map (bp from .bim / bed.matrix; optionally merge external cM map)
chr <- as.character(bm@snps$chr)
bp <- suppressWarnings(as.numeric(as.character(bm@snps$pos)))
chr[is.na(chr) | chr==""] <- "1"
if (all(is.na(bp))) bp <- seq_len(length(markers))

mm_bp <- data.table(marker=as.character(markers), marker_raw=as.character(markers_raw), chr=chr, bp=bp)

marker_map_bp_tsv <- file.path(base_out, "marker_map_bp.tsv")
fwrite(mm_bp[, .(marker, chr, bp)], marker_map_bp_tsv, sep="\t", na="NA", quote=FALSE)

normalize_chr <- function(x) {
  x <- tolower(as.character(x))
  x <- gsub("^chr", "", x)
  x <- gsub("^lg", "", x)
  x <- trimws(x)
  x <- sub("^0+", "", x)
  ifelse(nzchar(x), x, "1")
}

read_table_any <- function(path, sheet="") {
  path <- as.character(path)
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("xlsx", "xls")) {
    if (!requireNamespace("readxl", quietly=TRUE)) {
      cat("[prep_mppr_inputs] NOTE: installing readxl ...\n")
      install.packages("readxl", repos="https://cloud.r-project.org", quiet=TRUE)
    }
    sh <- sheet
    if (!nzchar(sh)) sh <- 1
    df <- as.data.table(readxl::read_excel(path, sheet=sh, na=c("", "NA", "NaN", ".")))
    return(df)
  }
  dt <- tryCatch({
    fread(path)
  }, error=function(e) {
    fread(path, sep=",", header=TRUE)
  })
  dt
}

mm_out <- copy(mm_bp)
mm_out[, chr2 := normalize_chr(chr)]
mm_out[, cm := NA_real_]

map_report <- data.table(step=character(0), n=integer(0), note=character(0))

if (nzchar(genetic_map_path) && file.exists(genetic_map_path)) {
  cat("[prep_mppr_inputs] genetic_map_path=", genetic_map_path, "\n")
  gm <- read_table_any(genetic_map_path, genetic_map_sheet)
  if (nrow(gm) > 0) {
    nms <- names(gm)
    nms_l <- tolower(nms)
    pick_col <- function(cands) {
      hit <- which(nms_l %in% cands)
      if (length(hit) > 0) return(nms[hit[1]])
      for (c in cands) {
        hit2 <- grep(c, nms_l, fixed=TRUE)
        if (length(hit2) > 0) return(nms[hit2[1]])
      }
      NA_character_
    }
    col_marker <- pick_col(c("marker","snp","rs","rsid","id","name"))
    col_chr <- pick_col(c("chr","chrom","chromosome","lg","linkage_group","linkagegroup"))
    col_bp <- pick_col(c("bp","pos","position","physical_pos","phys_pos"))
    col_cm <- pick_col(c("cm","cM","genetic_pos","gpos","genpos"))

    if (is.na(col_cm) && ncol(gm) >= 4) {
      col_chr <- nms[1]; col_marker <- nms[2]; col_cm <- nms[3]; col_bp <- nms[4]
    }

    gm2 <- as.data.table(gm)
    if (!is.na(col_marker)) gm2[, gm_marker := as.character(get(col_marker))] else gm2[, gm_marker := NA_character_]
    if (!is.na(col_chr)) gm2[, gm_chr := normalize_chr(get(col_chr))] else gm2[, gm_chr := NA_character_]
    if (!is.na(col_bp)) gm2[, gm_bp := suppressWarnings(as.numeric(as.character(get(col_bp))))] else gm2[, gm_bp := NA_real_]
    if (!is.na(col_cm)) gm2[, gm_cm := suppressWarnings(as.numeric(as.character(get(col_cm))))] else gm2[, gm_cm := NA_real_]
    gm2 <- gm2[!is.na(gm_cm)]

    if (nrow(gm2) > 0) {
      gm_m <- gm2[!is.na(gm_marker) & nzchar(gm_marker)]
      if (nrow(gm_m) > 0) gm_m <- gm_m[!duplicated(gm_marker)]
      gm_cb <- gm2[!is.na(gm_chr) & !is.na(gm_bp)]
      if (nrow(gm_cb) > 0) gm_cb <- gm_cb[!duplicated(paste0(gm_chr, ":", gm_bp))]

      # merge by marker ids
      if (tolower(map_merge_mode) %in% c("marker_then_chr_bp","by_marker")) {
        if (nrow(gm_m) > 0) {
          setkey(gm_m, gm_marker)
          mm_out[, cm := gm_m[.(marker_raw), gm_cm]]
          n_hit <- sum(!is.na(mm_out$cm))
          map_report <- rbind(map_report, data.table(step="by_marker_raw", n=n_hit, note="match marker_raw"), fill=TRUE)
          if (n_hit < nrow(mm_out)) {
            mm_out[is.na(cm), cm := gm_m[.(marker), gm_cm]]
            n_hit2 <- sum(!is.na(mm_out$cm))
            map_report <- rbind(map_report, data.table(step="by_marker", n=n_hit2, note="match marker"), fill=TRUE)
          }
        }
      }

      # merge by chr+bp
      if (tolower(map_merge_mode) %in% c("marker_then_chr_bp","by_chr_bp")) {
        if (nrow(gm_cb) > 0) {
          setkey(gm_cb, gm_chr, gm_bp)
          mm_out[is.na(cm) & !is.na(bp), cm := gm_cb[.(chr2, bp), gm_cm]]
          n_hit3 <- sum(!is.na(mm_out$cm))
          map_report <- rbind(map_report, data.table(step="by_chr_bp", n=n_hit3, note="match chr+bp"), fill=TRUE)
        }
      }

      # fill missing cM
      if (interpolate_missing_cm && any(is.na(mm_out$cm))) {
        slopes <- c()
        for (cc in unique(mm_out$chr2)) {
          dtc <- mm_out[chr2 == cc & !is.na(cm) & !is.na(bp)]
          if (nrow(dtc) >= 2) slopes <- c(slopes, (max(dtc$cm) - min(dtc$cm)) / (max(dtc$bp) - min(dtc$bp)))
        }
        slope_med <- if (length(slopes) > 0 && is.finite(median(slopes))) median(slopes) else 1e-6
        for (cc in unique(mm_out$chr2)) {
          idx <- which(mm_out$chr2 == cc)
          bpc <- mm_out$bp[idx]
          cmc <- mm_out$cm[idx]
          ok <- which(!is.na(cmc) & !is.na(bpc))
          if (length(ok) >= 2) {
            ord <- order(bpc[ok])
            x <- bpc[ok][ord]; y <- cmc[ok][ord]
            mm_out$cm[idx] <- approx(x, y, xout=bpc, rule=2)$y
          } else if (length(ok) == 1) {
            x0 <- bpc[ok[1]]; y0 <- cmc[ok[1]]
            mm_out$cm[idx] <- y0 + (bpc - x0) * slope_med
          }
        }
        if (any(is.na(mm_out$cm))) {
          mm_out[is.na(cm), cm := (bp - min(bp, na.rm=TRUE)) * 1e-6]
        }
        map_report <- rbind(map_report, data.table(step="interpolate", n=sum(!is.na(mm_out$cm)), note="filled cM"), fill=TRUE)
      }
    }
  }
}

use_cm <- any(!is.na(mm_out$cm)) && mean(!is.na(mm_out$cm)) > 0.1
if (use_cm) {
  mm <- mm_out[, .(marker=marker, chr=chr2, pos=cm)]
  map_report <- rbind(map_report, data.table(step="final", n=nrow(mm), note="marker_map.tsv uses cM"), fill=TRUE)
} else {
  mm <- mm_out[, .(marker=marker, chr=chr2, pos=bp)]
  map_report <- rbind(map_report, data.table(step="final", n=nrow(mm), note="marker_map.tsv uses bp"), fill=TRUE)
}

marker_map_tsv <- file.path(base_out, "marker_map.tsv")
fwrite(mm, marker_map_tsv, sep="\t", na="NA", quote=FALSE)

map_report_tsv <- file.path(base_out, "marker_map_merge_report.tsv")
fwrite(map_report, map_report_tsv, sep="\t", na="NA", quote=FALSE)

# helper to write geno tsv
write_geno_tsv <- function(iid_vec, out_path) {
  keep <- ids %in% iid_vec
  if (!any(keep)) stop("No rows to write for ", out_path)
  DT <- data.table(id = ids[keep])
  Xsub <- X[keep, , drop=FALSE]
  Xdt <- as.data.table(Xsub)
  setnames(Xdt, as.character(markers))
  DT <- cbind(DT, Xdt)
  fwrite(DT, out_path, sep="\t", na="NA", quote=FALSE)
}

# Determine which offspring/parents were actually kept
offspring_kept <- intersect(offspring_iids, ids)
parent_iids <- unique(parent_map$parent_iid)
parents_kept <- intersect(parent_iids, ids)

if (length(offspring_kept) < 10) {
  stop("Too few offspring samples after mapping/extraction: ", length(offspring_kept))
}
if (length(parents_kept) < 2) {
  stop("Too few parents samples after mapping/extraction: ", length(parents_kept))
}

# write geno matrices
geno_off_tsv <- file.path(base_out, "geno_off.tsv")
geno_par_tsv <- file.path(base_out, "geno_par.tsv")
write_geno_tsv(offspring_kept, geno_off_tsv)
write_geno_tsv(parents_kept, geno_par_tsv)

# phenotype_mppr.tsv
# Use original phenotype columns but replace id and cross
ph2 <- copy(ph)
ph2[, id := iid]
ph2[, cross := cross]
# Drop helper cols
for (c in c("raw_id","iid","fid","has_exact","base_key")) if (c %in% names(ph2)) ph2[, (c) := NULL]
# Keep only rows for offspring_kept
ph2 <- ph2[id %in% offspring_kept]
# Ensure id first
cols <- names(ph2)
cols <- c("id", setdiff(cols, "id"))
setcolorder(ph2, cols)
phenotype_mppr_tsv <- file.path(base_out, "phenotype_mppr.tsv")
fwrite(ph2, phenotype_mppr_tsv, sep="\t", na="NA", quote=FALSE)

# cross_ind.tsv + par_per_cross.tsv
cross_ind_tsv <- file.path(base_out, "cross_ind.tsv")
par_per_cross_tsv <- file.path(base_out, "par_per_cross.tsv")
fwrite(cross_ind[id %in% offspring_kept], cross_ind_tsv, sep="\t", na="NA", quote=FALSE)
fwrite(pp_out, par_per_cross_tsv, sep="\t", na="NA", quote=FALSE)

# also save parent mapping for transparency
fwrite(parent_map, file.path(base_out, "parent_mapping.tsv"), sep="\t", na="NA", quote=FALSE)

artifacts <- list(
  plink_prefix = plink_prefix,
  subset_plink_prefix = subset_pref,
  marker_map_tsv = marker_map_tsv,
  marker_map_bp_tsv = marker_map_bp_tsv,
  marker_map_merge_report_tsv = map_report_tsv,
  parent_replicate_choice_tsv = file.path(base_out, "parent_replicate_choice.tsv"),
  geno_off_tsv = geno_off_tsv,
  geno_par_tsv = geno_par_tsv,
  phenotype_mppr_tsv = phenotype_mppr_tsv,
  cross_ind_tsv = cross_ind_tsv,
  par_per_cross_tsv = par_per_cross_tsv,
  parent_mapping_tsv = file.path(base_out, "parent_mapping.tsv")
)
writeLines(toJSON(artifacts, auto_unbox=TRUE, pretty=TRUE), con=file.path(base_out, "artifacts.json"))
writeLines(toJSON(artifacts, auto_unbox=TRUE, pretty=TRUE), con=file.path(out_dir, "artifacts.json"))

cat("[prep_mppr_inputs] DONE\n")
cat("[prep_mppr_inputs] geno_off_tsv=", geno_off_tsv, "\n")
cat("[prep_mppr_inputs] geno_par_tsv=", geno_par_tsv, "\n")
cat("[prep_mppr_inputs] marker_map_tsv=", marker_map_tsv, "\n")
cat("[prep_mppr_inputs] phenotype_mppr_tsv=", phenotype_mppr_tsv, "\n")
cat("[prep_mppr_inputs] cross_ind_tsv=", cross_ind_tsv, "\n")
cat("[prep_mppr_inputs] par_per_cross_tsv=", par_per_cross_tsv, "\n")

