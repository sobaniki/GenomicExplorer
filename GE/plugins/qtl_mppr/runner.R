#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(jsonlite)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) return(args[i + 1])
  default
}

params_path <- get_arg("--params")
out_dir <- get_arg("--out")
if (is.null(params_path) || is.null(out_dir)) {
  stop("Usage: runner.R --params params.json --out out_dir")
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_dir, "run.log")
sink(log_file, split = TRUE)

cat("[qtl_mppr] start\n")
cat("[qtl_mppr] params_path=", params_path, "\n")
cat("[qtl_mppr] out_dir=", out_dir, "\n")

p <- fromJSON(params_path)

require_or_install <- function(pkg) {
  if (requireNamespace(pkg, quietly = TRUE)) return(TRUE)
  cat(sprintf("[qtl_mppr] NOTE: installing missing R package '%s' (CRAN) ...\n", pkg))
  ok <- tryCatch({
    install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
    TRUE
  }, error = function(e) {
    cat(sprintf("[qtl_mppr] ERROR: install.packages('%s') failed: %s\n", pkg, e$message))
    FALSE
  })
  if (!ok || !requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("R package '%s' is required but not installed. Try: install.packages('%s')", pkg, pkg))
  }
  TRUE
}

require_or_install("mppR")
require_or_install("ggplot2")

suppressPackageStartupMessages({
  library(mppR)
  library(ggplot2)
})

write_artifacts <- function(artifacts) {
  writeLines(toJSON(artifacts, auto_unbox = TRUE, pretty = TRUE),
             con = file.path(out_dir, "artifacts.json"))
}

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && nzchar(as.character(a)[1])) a else b

detect_id_col <- function(dt) {
  nms <- names(dt)
  for (cand in c("id", "ID", "sample", "Sample", "IID", "geno", "geno.id", "genotype", "Genotype")) {
    if (cand %in% nms) return(cand)
  }
  nms[1]
}

normalize_chr <- function(x) {
  x <- as.character(x)
  x <- gsub("^chr", "", x, ignore.case = TRUE)
  # allow e.g. "1", "01", "1A" -> extract first number
  n <- suppressWarnings(as.integer(gsub("[^0-9].*$", "", x)))
  ifelse(is.na(n), x, n)
}

normalize_geno_matrix <- function(mat) {
  # mat: data.frame / matrix (may be numeric or character)
  if (is.data.frame(mat)) {
    # keep as data.frame for efficient operations
    df <- mat
  } else {
    df <- as.data.frame(mat, stringsAsFactors = FALSE)
  }

  # fast path: numeric
  if (all(vapply(df, function(col) is.numeric(col) || is.integer(col), logical(1)))) {
    for (j in seq_len(ncol(df))) {
      v <- df[[j]]
      vv <- as.integer(round(v))
      vv[is.na(v)] <- NA_integer_
      vv[vv < 0 | vv > 2] <- NA_integer_
      out <- rep(NA_character_, length(vv))
      out[vv == 0] <- "AA"
      out[vv == 1] <- "AB"
      out[vv == 2] <- "BB"
      df[[j]] <- out
    }
    return(as.matrix(df))
  }

  # character-ish path
  for (j in seq_len(ncol(df))) {
    v <- as.character(df[[j]])
    v[!nzchar(v)] <- NA_character_
    v[v %in% c(".", "./.", ".|.", "NA", "NaN", "nan")] <- NA_character_

    # numeric strings
    is_num <- grepl("^-?\\d+(\\.\\d+)?$", v)
    if (any(is_num, na.rm = TRUE)) {
      vv <- suppressWarnings(as.integer(round(as.numeric(v))))
      vv[is.na(v)] <- NA_integer_
      vv[vv < 0 | vv > 2] <- NA_integer_
      vv[!is_num] <- NA_integer_
      out <- rep(NA_character_, length(v))
      out[vv == 0] <- "AA"
      out[vv == 1] <- "AB"
      out[vv == 2] <- "BB"
      # keep original non-numeric where present
      keep <- !is_num & !is.na(v)
      out[keep] <- toupper(gsub("[\\/|]", "", v[keep]))
      # fix single-letter alleles -> duplicate
      one <- keep & nchar(out) == 1
      out[one] <- paste0(out[one], out[one])
      df[[j]] <- out
      next
    }

    # VCF-style 0/0 etc
    v2 <- v
    v2[v2 %in% c("0/0", "0|0")] <- "AA"
    v2[v2 %in% c("0/1", "1/0", "0|1", "1|0")] <- "AB"
    v2[v2 %in% c("1/1", "1|1")] <- "BB"

    # allele strings like A/G, A|G, AG
    v2 <- toupper(gsub("[\\/|]", "", v2))
    one <- !is.na(v2) & nchar(v2) == 1
    v2[one] <- paste0(v2[one], v2[one])
    df[[j]] <- v2
  }
  as.matrix(df)
}

read_geno_tsv <- function(path) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) stop("geno TSV not found: ", path)
  dt <- fread(path)
  if (ncol(dt) < 2) stop("geno TSV must have >=2 columns: id + markers")
  id_col <- detect_id_col(dt)
  ids <- as.character(dt[[id_col]])
  dt[[id_col]] <- NULL
  df <- as.data.frame(dt)
  rownames(df) <- ids
  normalize_geno_matrix(df)
}

read_pheno_tsv <- function(path, trait) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) stop("pheno TSV not found: ", path)
  dt <- fread(path)
  if (ncol(dt) < 2) stop("pheno TSV must have >=2 columns: id + trait(s)")
  id_col <- detect_id_col(dt)
  ids <- as.character(dt[[id_col]])
  dt[[id_col]] <- NULL

  # optional cross column that should not be treated as phenotype
  for (cand in c("cross", "Cross", "family", "Family", "pop", "population")) {
    if (cand %in% names(dt)) {
      # keep for auto cross.ind if needed
      attr(dt, "_cross_col") <- cand
      break
    }
  }

  if (!(trait %in% names(dt))) {
    # auto pick first numeric column, but avoid common metadata columns (especially in NAM phenotypes)
    meta_cols <- c(
      "Panel", "panel",
      "Family_Num", "Family.num", "FamilyNum", "family_num",
      "Family_Inbred_Name", "family_inbred_name",
      "Entry_Num", "Entry.num", "EntryNum", "entry_num",
      "Z_Num", "Z.num", "ZNum", "z_num"
    )
    # start with numeric/integer columns
    num_cols <- names(dt)[vapply(dt, function(x) is.numeric(x) || is.integer(x), logical(1))]
    # if none, try coercion (tolerate '.', 'NA', etc.)
    if (length(num_cols) == 0) {
      for (c in names(dt)) dt[[c]] <- suppressWarnings(as.numeric(dt[[c]]))
      num_cols <- names(dt)[vapply(dt, function(x) is.numeric(x) || is.integer(x), logical(1))]
    }
    if (length(num_cols) == 0) stop("trait column not found and no numeric columns detected in pheno TSV")
    num_cols2 <- setdiff(num_cols, meta_cols)
    pick <- if (length(num_cols2) > 0) num_cols2[1] else num_cols[1]
    cat("[qtl_mppr] WARN: trait column '", trait, "' not found/empty; using numeric column '", pick, "'\n", sep="")
    trait <- pick
  }

  y <- suppressWarnings(as.numeric(dt[[trait]]))

  # IMPORTANT: some versions of mppR (QC.mppData) choke when pheno is simplified to a vector.
  # To keep it safely 2D, always provide >=2 phenotype columns for SIM/CIM.
  ph <- cbind(y, .dummy_keep2d = rep(0, length(y)))
  colnames(ph)[1] <- trait
  rownames(ph) <- ids

  list(pheno = ph, ids = ids, trait = trait, cross_col = attr(dt, "_cross_col"))
}


read_pheno_tsv_multi <- function(path, traits) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) stop("pheno TSV not found: ", path)
  dt <- fread(path)
  if (ncol(dt) < 2) stop("pheno TSV must have >=2 columns: id + trait(s)")
  id_col <- detect_id_col(dt)
  ids <- as.character(dt[[id_col]])
  dt[[id_col]] <- NULL

  # optional cross column
  for (cand in c("cross", "Cross", "family", "Family", "pop", "population")) {
    if (cand %in% names(dt)) {
      attr(dt, "_cross_col") <- cand
      break
    }
  }

  # traits may come as a single comma/space separated string
  if (length(traits) == 1) {
    traits0 <- as.character(traits)[1]
    parts <- unlist(strsplit(traits0, "[,;\\s]+"))
    parts <- parts[nzchar(parts)]
    traits <- parts
  }
  traits <- unique(as.character(traits))
  traits <- traits[nzchar(traits)]

  if (length(traits) == 0) stop("GE traits are required (comma-separated column names).")

  miss <- setdiff(traits, names(dt))
  if (length(miss) > 0) {
    stop("Trait columns not found in pheno TSV: ", paste(miss, collapse = ", "))
  }

  mat <- sapply(traits, function(tr) suppressWarnings(as.numeric(dt[[tr]])))
  if (is.null(dim(mat))) mat <- matrix(mat, ncol = 1)
  colnames(mat) <- traits
  rownames(mat) <- ids

  list(pheno = mat, ids = ids, traits = traits, cross_col = attr(dt, "_cross_col"))
}

read_map_tsv <- function(path) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) stop("map TSV not found: ", path)
  dt <- fread(path)
  if (ncol(dt) < 3) stop("map TSV must have >=3 columns: marker, chr, pos")
  nms <- names(dt)
  mk_col <- NULL
  for (cand in c("mk.names", "marker", "Marker", "id", "ID", "SNP", "rs")) {
    if (cand %in% nms) { mk_col <- cand; break }
  }
  if (is.null(mk_col)) mk_col <- nms[1]
  chr_col <- NULL
  for (cand in c("chr", "chrom", "chromosome", "Chr", "CHR")) {
    if (cand %in% nms) { chr_col <- cand; break }
  }
  if (is.null(chr_col)) chr_col <- nms[2]
  pos_col <- NULL
  for (cand in c("pos.cM", "pos_cm", "pos", "cM", "cm", "position", "bp", "POS")) {
    if (cand %in% nms) { pos_col <- cand; break }
  }
  if (is.null(pos_col)) pos_col <- nms[3]

  map <- data.frame(
    mk.names = as.character(dt[[mk_col]]),
    chr = normalize_chr(dt[[chr_col]]),
    pos.cM = suppressWarnings(as.numeric(dt[[pos_col]])),
    stringsAsFactors = FALSE
  )
  map <- map[!is.na(map$mk.names) & nzchar(map$mk.names), , drop = FALSE]
  map$chr <- suppressWarnings(as.numeric(map$chr))
  # if chr is still NA, keep as-is string; mppR expects numeric, but we try best-effort
  if (all(is.na(map$chr))) {
    stop("map 'chr' column could not be parsed as numeric. Please use 1,2,3,...")
  }
  map <- map[order(map$chr, map$pos.cM), , drop = FALSE]
  map
}

read_cross_ind <- function(path, ids, pheno_dt_fallback = NULL, cross_col = NULL) {
  if (!is.null(path) && nzchar(path) && file.exists(path)) {
    dt <- fread(path)
    if (ncol(dt) < 2) stop("cross_ind TSV must have >=2 columns: id + cross")
    id_col <- detect_id_col(dt)
    cross_col2 <- setdiff(names(dt), id_col)[1]
    key <- as.character(dt[[id_col]])
    val <- as.character(dt[[cross_col2]])
    m <- setNames(val, key)
    out <- unname(m[ids])
    if (any(is.na(out))) stop("cross_ind is missing some phenotype ids (first missing: ", ids[which(is.na(out))[1]], ")")
    return(out)
  }

  # fallback: try pheno data column (already captured by read_pheno_tsv)
  if (!is.null(pheno_dt_fallback) && !is.null(cross_col) && cross_col %in% names(pheno_dt_fallback)) {
    m <- setNames(as.character(pheno_dt_fallback[[cross_col]]), ids)
    return(unname(m[ids]))
  }

  # last resort: derive from prefix before '_' or first 4 chars
  cat("[qtl_mppr] WARN: cross_ind_tsv not provided; inferring cross.ind from genotype id prefix.\n")
  if (any(grepl("_", ids))) {
    return(sub("_.*$", "", ids))
  }
  substr(ids, 1, 4)
}

read_par_per_cross <- function(path) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) stop("par_per_cross TSV not found: ", path)
  dt <- fread(path)
  if (ncol(dt) < 3) stop("par_per_cross TSV must have >=3 columns: cross, parent1, parent2")
  nms <- names(dt)
  cross_col <- NULL
  for (cand in c("cross", "Cross", "cross.ind", "cross_id")) {
    if (cand %in% nms) { cross_col <- cand; break }
  }
  if (is.null(cross_col)) cross_col <- nms[1]
  # parents: next two columns
  others <- setdiff(nms, cross_col)
  p1 <- others[1]
  p2 <- others[2]
  mat <- cbind(as.character(dt[[cross_col]]), as.character(dt[[p1]]), as.character(dt[[p2]]))
  mat
}

# ---- params ----
geno_off_tsv <- p$geno_off_tsv %||% ""
geno_par_tsv <- p$geno_par_tsv %||% ""
map_tsv <- p$map_tsv %||% ""
pheno_tsv <- p$pheno_tsv %||% ""
trait <- as.character(p$trait %||% "")
cross_ind_tsv <- p$cross_ind_tsv %||% ""
par_per_cross_tsv <- p$par_per_cross_tsv %||% ""


method <- toupper(as.character(p$method %||% "SIM"))
if (!(method %in% c("SIM", "CIM", "GE"))) method <- "SIM"

# Cofactor selection (used by CIM and GE)
cof_threshold <- suppressWarnings(as.numeric(p$cof_threshold %||% ifelse(method == "GE", 4, 3)))
if (is.na(cof_threshold) || cof_threshold <= 0) cof_threshold <- ifelse(method == "GE", 4, 3)

cof_window <- suppressWarnings(as.numeric(p$cof_window %||% 50))
if (is.na(cof_window) || cof_window <= 0) cof_window <- 50

cim_window <- suppressWarnings(as.numeric(p$cim_window %||% 20))
if (is.na(cim_window) || cim_window <= 0) cim_window <- 20

# GE (GxE) settings
traits_ge_raw <- as.character(p$traits_ge %||% "")
env_names_raw <- as.character(p$env_names %||% "")
VCOV <- toupper(as.character(p$VCOV %||% "UN"))
if (!(VCOV %in% c("UN", "CS", "CSE", "CS_CSE"))) VCOV <- "UN"
VCOV_data <- tolower(as.character(p$VCOV_data %||% "unique"))
if (!(VCOV_data %in% c("unique", "minus_cof"))) VCOV_data <- "unique"
ge_sim_only <- as.logical(p$SIM_only %||% FALSE)
if (is.na(ge_sim_only)) ge_sim_only <- FALSE
cof_red <- as.logical(p$cof_red %||% FALSE)
if (is.na(cof_red)) cof_red <- FALSE
cof_pval_sign <- suppressWarnings(as.numeric(p$cof_pval_sign %||% 0.1))
if (is.na(cof_pval_sign) || cof_pval_sign <= 0 || cof_pval_sign >= 1) cof_pval_sign <- 0.1
ref_par <- as.character(p$ref_par %||% "")
maxIter <- as.integer(p$maxIter %||% 100L)
if (is.na(maxIter) || maxIter < 1) maxIter <- 100L
msMaxIter <- as.integer(p$msMaxIter %||% 100L)
if (is.na(msMaxIter) || msMaxIter < 1) msMaxIter <- 100L

split_list <- function(x) {
  if (is.null(x)) return(character(0))
  x <- as.character(x)
  x <- paste(x, collapse = ",")
  parts <- unlist(strsplit(x, "[,;\\s]+"))
  parts <- parts[nzchar(parts)]
  unique(parts)
}

Q_eff <- tolower(as.character(p$Q_eff %||% "cr"))
if (!(Q_eff %in% c("cr", "par", "biall", "anc"))) Q_eff <- "cr"

plot_gen_eff <- as.logical(p$plot_gen_eff %||% FALSE)
if (is.na(plot_gen_eff)) plot_gen_eff <- FALSE

pop_type <- toupper(as.character(p$pop_type %||% "RIL"))
if (!(pop_type %in% c("F", "BC", "BCSFT", "DH", "RIL"))) pop_type <- "RIL"

type_mating <- tolower(as.character(p$type_mating %||% "selfing"))
if (!(type_mating %in% c("selfing", "sib.mat"))) type_mating <- "selfing"

F_gen <- p$F_gen
if (!is.null(F_gen)) F_gen <- as.integer(F_gen)
BC_gen <- p$BC_gen
if (!is.null(BC_gen)) BC_gen <- as.integer(BC_gen)

error_prob <- as.numeric(p$error_prob %||% 1e-4)
if (is.na(error_prob) || error_prob <= 0 || error_prob >= 0.5) error_prob <- 1e-4
map_function <- tolower(as.character(p$map_function %||% "haldane"))
if (!(map_function %in% c("haldane", "kosambi", "c-f", "morgan"))) map_function <- "haldane"

n_cores <- as.integer(p$n_cores %||% 1L)
if (is.na(n_cores) || n_cores < 1) n_cores <- 1L

threshold <- suppressWarnings(as.numeric(p$threshold %||% ifelse(method == "GE", 4, 3)))
if (is.na(threshold) || threshold <= 0) threshold <- ifelse(method == "GE", 4, 3)
window <- suppressWarnings(as.numeric(p$window %||% ifelse(method == "GE", 20, 50)))
if (is.na(window) || window <= 0) window <- ifelse(method == "GE", 20, 50)

do_perm <- as.logical(p$do_perm %||% FALSE)
if (is.na(do_perm)) do_perm <- FALSE
n_perm <- as.integer(p$n_perm %||% 0L)
if (is.na(n_perm) || n_perm < 0) n_perm <- 0L

# Backward/GUI-friendly behavior:
# If user specifies n_perm > 0, we should run permutation even when do_perm is missing.
if (n_perm > 0L && !isTRUE(do_perm)) {
  do_perm <- TRUE
}
q_val <- as.numeric(p$q_val %||% 0.95)
if (is.na(q_val) || q_val <= 0 || q_val >= 1) q_val <- 0.95

seed <- as.integer(p$seed %||% 1L)
if (is.na(seed) || seed < 1) seed <- 1L
set.seed(seed)

pos_unit <- as.character(p$pos_unit %||% "cM")

cat("[qtl_mppr] Q_eff=", Q_eff, " pop_type=", pop_type, " n_cores=", n_cores, "\n")

# ---- read inputs ----
cat("[qtl_mppr] reading genotype matrices ...\n")
geno_off <- read_geno_tsv(geno_off_tsv)
geno_par <- read_geno_tsv(geno_par_tsv)

cat("[qtl_mppr] reading map ...\n")
map <- read_map_tsv(map_tsv)

# align marker order
mk <- map$mk.names
off_cols <- colnames(geno_off)
par_cols <- colnames(geno_par)

common <- intersect(mk, intersect(off_cols, par_cols))
if (length(common) < 5) {
  stop("Too few common markers between map and genotype matrices. map markers=", length(mk),
       " geno_off=", length(off_cols), " geno_par=", length(par_cols),
       " common=", length(common))
}
if (length(common) < length(mk)) {
  cat("[qtl_mppr] WARN: some map markers are missing from genotype matrices. Using common markers only: ",
      length(common), "/", length(mk), "\n", sep="")
}
map <- map[map$mk.names %in% common, , drop=FALSE]
mk <- map$mk.names
geno_off <- geno_off[, mk, drop=FALSE]
geno_par <- geno_par[, mk, drop=FALSE]


cat("[qtl_mppr] reading phenotype ...\n")
if (method == "GE") {
  traits_ge <- split_list(traits_ge_raw %||% trait)
  if (length(traits_ge) < 2) stop("GE method requires 2+ trait columns (comma-separated), e.g. trait_env1,trait_env2")
  ph <- read_pheno_tsv_multi(pheno_tsv, traits_ge)
  pheno <- ph$pheno
  ids <- rownames(pheno)
  trait <- traits_ge
} else {
  ph <- read_pheno_tsv(pheno_tsv, trait)
  pheno <- ph$pheno
  trait <- ph$trait
  ids <- rownames(pheno)
}

cat("[qtl_mppr] reading cross.ind ...\n")
cross.ind <- read_cross_ind(cross_ind_tsv, ids)

cat("[qtl_mppr] reading par.per.cross ...\n")
par.per.cross <- read_par_per_cross(par_per_cross_tsv)

# sanity checks
parents <- unique(c(as.character(par.per.cross[,2]), as.character(par.per.cross[,3])))
missing_par <- setdiff(parents, rownames(geno_par))
if (length(missing_par) > 0) {
  stop("Parents listed in par_per_cross are missing from geno_par matrix: ", paste(missing_par, collapse=", "))
}

# ---- build mppData ----
cat("[qtl_mppr] create.mppData ...\n")
mppData <- create.mppData(
  geno.off = geno_off,
  geno.par = geno_par,
  map = map,
  pheno = pheno,
  cross.ind = cross.ind,
  par.per.cross = par.per.cross
)

cat("[qtl_mppr] QC.mppData ...\n")
mppData <- QC.mppData(mppData)

cat("[qtl_mppr] IBS.mppData ...\n")
mppData <- IBS.mppData(mppData)

cat("[qtl_mppr] IBD.mppData ...\n")
ibd_ok <- TRUE
mppData <- tryCatch({
  args_ibd <- list(
    mppData = mppData,
    het.miss.par = TRUE,
    type = pop_type,
    error.prob = error_prob,
    map.function = map_function
  )
  if (pop_type == "RIL") args_ibd$type.mating <- type_mating
  if (pop_type == "F" && !is.null(F_gen) && !is.na(F_gen)) args_ibd$F.gen <- F_gen
  if (pop_type == "BC" && !is.null(BC_gen) && !is.na(BC_gen)) args_ibd$BC.gen <- BC_gen
  do.call(IBD.mppData, args_ibd)
}, error = function(e) {
  ibd_ok <<- FALSE
  cat("[qtl_mppr] WARN: IBD.mppData failed: ", e$message, "\n", sep="")
  mppData
})

if (Q_eff == "anc") {
  stop("Q_eff='anc' requires parental clustering (par_clu) which is not generated automatically in the CRAN version of mppR. Use Q_eff='cr'/'par'/'biall'.")
}

# ---- scan ----
cat("[qtl_mppr] method=", method, "\n", sep="")

Qprof <- NULL
cofactors_obj <- NULL
cofactors_df <- NULL
env_profiles <- list()
env_names <- character(0)
profile_kind <- ""

if (method == "SIM") {
  cat("[qtl_mppr] mpp_SIM ...\n")
  Qprof <- mpp_SIM(mppData = mppData, trait = trait, Q.eff = Q_eff, plot.gen.eff = plot_gen_eff, n.cores = n_cores)
  profile_kind <- "SIM"
} else if (method == "CIM") {
  cat("[qtl_mppr] CIM: running mpp_SIM for cofactor selection ...\n")
  SIM0 <- mpp_SIM(mppData = mppData, trait = trait, Q.eff = Q_eff, plot.gen.eff = FALSE, n.cores = n_cores)
  cofactors_obj <- QTL_select(Qprof = SIM0, threshold = cof_threshold, window = cof_window, verbose = FALSE)
  cat("[qtl_mppr] mpp_CIM ... window=", cim_window, "\n", sep="")
  Qprof <- mpp_CIM(mppData = mppData, trait = trait, Q.eff = Q_eff, cofactors = cofactors_obj,
                   window = cim_window, plot.gen.eff = plot_gen_eff, n.cores = n_cores)
  profile_kind <- "CIM"
} else if (method == "GE") {
  cat("[qtl_mppr] GE: mppGE_SIM ... VCOV=", VCOV, " nEnv=", length(trait), "\n", sep="")
  env_names <- split_list(env_names_raw)
  if (length(env_names) == 0) env_names <- paste0("Env", seq_along(trait))
  if (length(env_names) != length(trait)) {
    cat("[qtl_mppr] WARN: env_names length != traits length; using auto Env1.. \n")
    env_names <- paste0("Env", seq_along(trait))
  }

  if (!ibd_ok) {
    stop("GE method requires IBD.mppData. Please check population type/settings (RIL/F/BC/DH).")
  }

  SIM_GE <- mppGE_SIM(mppData = mppData, trait = trait, VCOV = VCOV, ref_par = if (nzchar(ref_par)) ref_par else NULL,
                      n.cores = n_cores, maxIter = maxIter, msMaxIter = msMaxIter)

  cofactors_obj <- QTL_select(Qprof = SIM_GE, threshold = cof_threshold, window = cof_window, verbose = FALSE)

  if (ge_sim_only) {
    Qprof <- SIM_GE
    profile_kind <- "GE_SIM"
  } else {
    cat("[qtl_mppr] GE: mppGE_CIM ... VCOV_data=", VCOV_data, " cof_red=", cof_red, " window=", cim_window, "\n", sep="")
    Qprof <- mppGE_CIM(
      mppData = mppData,
      trait = trait,
      VCOV = VCOV,
      VCOV_data = VCOV_data,
      cofactors = cofactors_obj,
      cof_red = cof_red,
      cof_pval_sign = cof_pval_sign,
      window = cim_window,
      ref_par = if (nzchar(ref_par)) ref_par else NULL,
      n.cores = n_cores,
      maxIter = maxIter,
      msMaxIter = msMaxIter
    )
    profile_kind <- "GE_CIM"
  }
} else {
  stop("Unknown method: ", method)
}

qdf <- as.data.frame(Qprof)

# ---- standardize profile: chr, pos, lod ----
mk_col <- NULL
for (cand in c("mk.names", "marker", "mk")) if (cand %in% names(qdf)) { mk_col <- cand; break }
if (is.null(mk_col)) mk_col <- names(qdf)[1]

chr_col <- NULL
for (cand in c("chr", "Chr", "CHR")) if (cand %in% names(qdf)) { chr_col <- cand; break }
if (is.null(chr_col)) chr_col <- names(qdf)[2]

pos_cm_col <- NULL
for (cand in c("pos.cM", "pos_cm", "cM", "cm")) if (cand %in% names(qdf)) { pos_cm_col <- cand; break }
pos_col <- NULL
for (cand in c("pos", "position", "Position")) if (cand %in% names(qdf)) { pos_col <- cand; break }
pos_use <- if (!is.null(pos_cm_col)) pos_cm_col else pos_col
if (is.null(pos_use)) stop("QTL profile is missing a position column (pos.cM/pos)")

lod_col <- NULL
for (cand in c("-log10(p.val)", "log10.p.val", "log10p", "LOD", "lod")) {
  if (cand %in% names(qdf)) { lod_col <- cand; break }
}
if (is.null(lod_col)) {
  numc <- names(qdf)[vapply(qdf, is.numeric, logical(1))]
  if (length(numc) == 0) stop("QTL profile has no numeric columns to use as LOD")
  lod_col <- tail(numc, 1)
}

lod_profile <- data.frame(
  chr = as.character(qdf[[chr_col]]),
  pos = suppressWarnings(as.numeric(qdf[[pos_use]])),
  lod = suppressWarnings(as.numeric(qdf[[lod_col]])),
  stringsAsFactors = FALSE
)
lod_profile <- lod_profile[is.finite(lod_profile$pos) & is.finite(lod_profile$lod), , drop = FALSE]

fwrite(lod_profile, file = file.path(out_dir, "lod_profile.tsv"), sep = "	", quote = FALSE)

# ---- GE: write per-environment profiles (best-effort) ----
if (method == "GE") {
  n_env <- length(trait)
  base_cols <- unique(c(mk_col, chr_col, pos_use, lod_col))
  cand_cols <- setdiff(names(qdf), base_cols)
  num_cols <- cand_cols[vapply(qdf[cand_cols], function(x) is.numeric(x) || is.integer(x), logical(1))]
  pval_cols <- num_cols[vapply(qdf[num_cols], function(x) {
    xx <- suppressWarnings(as.numeric(x))
    xx <- xx[is.finite(xx)]
    if (length(xx) < 20) return(FALSE)
    mn <- min(xx); mx <- max(xx)
    mn >= 0 && mx <= 1
  }, logical(1))]
  pval_cols <- pval_cols[order(match(pval_cols, names(qdf)))]
  pval_cols <- head(pval_cols, n_env)

  if (length(pval_cols) == n_env) {
    for (i in seq_len(n_env)) {
      pv <- suppressWarnings(as.numeric(qdf[[pval_cols[i]]]))
      lod_env <- -log10(pmax(pv, 1e-300))
      df_env <- data.frame(chr = lod_profile$chr, pos = lod_profile$pos, lod = lod_env)
      env <- env_names[i]
      safe_env <- gsub("[^A-Za-z0-9_\\-]+", "_", env)
      fn <- paste0("lod_profile_env_", safe_env, ".tsv")
      fwrite(df_env, file = file.path(out_dir, fn), sep = "	", quote = FALSE)
      env_profiles[[env]] <- fn
    }
  } else {
    cat("[qtl_mppr] NOTE: could not infer per-environment p-value columns from GE profile; only global profile will be saved.\n")
  }
}

# ---- write cofactors (optional) ----
if (!is.null(cofactors_obj)) {
  cf <- tryCatch(as.data.frame(cofactors_obj), error = function(e) NULL)
  if (!is.null(cf) && nrow(cf) > 0) {
    cofactors_df <- cf
    fwrite(cf, file = file.path(out_dir, "cofactors.tsv"), sep = "	", quote = FALSE)
  }
}

# ---- permutation threshold (optional; SIM/CIM only) ----
thr_used <- threshold
perm_tsv <- file.path(out_dir, "perm_thresholds.tsv")
if (method != "GE" && do_perm && n_perm > 0) {
  cat("[qtl_mppr] mpp_perm N=", n_perm, " q.val=", q_val, "\n", sep="")
  Perm <- mpp_perm(mppData = mppData, trait = trait, Q.eff = Q_eff, N = n_perm, q.val = q_val, verbose = FALSE, n.cores = n_cores)
  thr <- tryCatch(as.numeric(Perm$q.val)[1], error = function(e) NA_real_)
  if (is.finite(thr)) thr_used <- thr
  alpha_out <- 1 - q_val
  data.table::fwrite(
    data.frame(alpha = alpha_out, threshold = thr_used, n_perm = n_perm, q_val = q_val),
    perm_tsv,
    sep = "	",
    quote = FALSE
  )
} else {
  data.table::fwrite(
    data.frame(alpha = NA_real_, threshold = thr_used, n_perm = ifelse(method == "GE", NA_integer_, 0L), q_val = ifelse(method == "GE", NA_real_, NA_real_)),
    perm_tsv,
    sep = "	",
    quote = FALSE
  )
}

# ---- peaks ----
cat("[qtl_mppr] QTL_select threshold=", thr_used, " window=", window, "\n", sep="")
QTL <- QTL_select(Qprof = Qprof, threshold = thr_used, window = window, verbose = FALSE)
peaks <- tryCatch(as.data.frame(QTL), error = function(e) data.frame())
if (nrow(peaks) == 0) {
  i <- which.max(lod_profile$lod)
  peaks <- data.frame(chr = lod_profile$chr[i], pos = lod_profile$pos[i], lod = lod_profile$lod[i])
} else {
  if (!("chr" %in% names(peaks))) {
    for (cand in c("Chr", "chrom", "chromosome")) if (cand %in% names(peaks)) peaks$chr <- peaks[[cand]]
  }
  pos_p <- NULL
  for (cand in c("pos.cM", "pos", "cM", "cm")) if (cand %in% names(peaks)) { pos_p <- cand; break }
  if (is.null(pos_p)) pos_p <- names(peaks)[min(4, ncol(peaks))]
  lod_p <- NULL
  for (cand in c("-log10(p.val)", lod_col, "log10p", "LOD", "lod")) if (cand %in% names(peaks)) { lod_p <- cand; break }
  if (is.null(lod_p)) {
    numc <- names(peaks)[vapply(peaks, is.numeric, logical(1))]
    lod_p <- tail(numc, 1)
  }
  peaks$pos <- suppressWarnings(as.numeric(peaks[[pos_p]]))
  peaks$lod <- suppressWarnings(as.numeric(peaks[[lod_p]]))
  peaks$chr <- as.character(peaks$chr)
  keep_cols <- unique(c("chr", "pos", "lod", intersect(names(peaks), c("mk.names", "pos.cM"))))
  peaks <- peaks[, keep_cols, drop = FALSE]
}
fwrite(peaks, file = file.path(out_dir, "peaks.tsv"), sep = "	", quote = FALSE)

# ---- plots ----
cat("[qtl_mppr] plotting ...\n")
plot_png <- file.path(out_dir, "lod_plot.png")

dfp <- lod_profile
dfp$chr_num <- suppressWarnings(as.integer(dfp$chr))
dfp <- dfp[order(dfp$chr_num, dfp$pos), , drop=FALSE]
chrs <- unique(dfp$chr_num)
chr_max <- tapply(dfp$pos, dfp$chr_num, max, na.rm=TRUE)
gap <- max(1, median(chr_max, na.rm=TRUE) * 0.02)
offs <- cumsum(c(0, head(chr_max + gap, -1)))
names(offs) <- as.character(sort(chrs))
dfp$x <- dfp$pos + offs[as.character(dfp$chr_num)]
peaks$x <- peaks$pos + offs[as.character(suppressWarnings(as.integer(peaks$chr)))]

png(plot_png, width=1400, height=800)
op <- par(mar=c(5,5,4,2))
main_lab <- if (method == "GE") {
  paste0("mppR ", profile_kind, " (VCOV=", VCOV, "): ", paste(trait, collapse=","))
} else {
  paste0("mppR ", profile_kind, " (", Q_eff, "): ", trait)
}
plot(dfp$x, dfp$lod, type="l", lwd=2, xlab=paste0("Genome position (", pos_unit, ")"), ylab="-log10(p)",
     main=main_lab)
if (is.finite(thr_used)) abline(h=thr_used, lty=2)
if (nrow(peaks) > 0) points(peaks$x, peaks$lod, pch=19)

ticks <- sapply(sort(chrs), function(c) {
  rng <- dfp$x[dfp$chr_num == c]
  if (length(rng) == 0) return(NA_real_)
  mean(range(rng))
})
ticks <- ticks[is.finite(ticks)]
axis(1, at=ticks, labels=names(ticks))
par(op)
dev.off()

# ---- artifacts ----
art <- list(
  lod_profile_tsv = "lod_profile.tsv",
  peaks_tsv = "peaks.tsv",
  perm_thresholds_tsv = "perm_thresholds.tsv",
  lod_plot_png = "lod_plot.png",
  cofactors_tsv = if (!is.null(cofactors_df)) "cofactors.tsv" else "",
  method = method,
  profile_kind = profile_kind,
  trait = trait,
  env_names = env_names,
  env_profiles = env_profiles,
  Q_eff = Q_eff,
  pop_type = pop_type,
  threshold = thr_used,
  window = window,
  cof_threshold = cof_threshold,
  cof_window = cof_window,
  cim_window = cim_window,
  do_perm = do_perm,
  n_perm = n_perm,
  q_val = q_val,
  seed = seed,
  pos_unit = pos_unit,
  ibd_ok = ibd_ok,
  VCOV = VCOV,
  VCOV_data = VCOV_data
)
write_artifacts(art)

cat("[qtl_mppr] done\n")
sink()
