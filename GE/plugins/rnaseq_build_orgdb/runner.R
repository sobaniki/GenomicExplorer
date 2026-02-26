#!/usr/bin/env Rscript
# Build OrgDb package from user-provided gene-to-GO mapping (Phase A)

suppressWarnings(suppressMessages({
  library(jsonlite)
  library(data.table)
}))

fail <- function(msg, code=1) {
  cat("[rnaseq_build_orgdb][ERROR] ", msg, "\n", sep="")
  quit(status=code)
}

`%||%` <- function(a, b) { if (is.null(a) || length(a)==0 || (is.character(a)&&a=="")) b else a }

as_bool <- function(x) {
  if (is.null(x)) return(FALSE)
  if (is.logical(x)) return(isTRUE(x))
  if (is.character(x)) return(tolower(x) %in% c("1","true","t","yes","y"))
  if (is.numeric(x)) return(x != 0)
  FALSE
}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 4 || args[1] != "--params" || args[3] != "--out") {
  cat("Usage: runner.R --params params.json --out out_dir\n")
  quit(status=2)
}

params_path <- args[2]
out_dir <- args[4]
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

cat("[rnaseq_build_orgdb] start\n")
cat("[rnaseq_build_orgdb] params_path=", params_path, "\n", sep="")
cat("[rnaseq_build_orgdb] out_dir=", out_dir, "\n", sep="")

params <- tryCatch(jsonlite::fromJSON(params_path), error=function(e) fail(paste("Failed to parse params.json:", e$message)))

mode <- tolower(params$mode %||% "gene2go_tsv")
pkg_name <- params$pkg_name %||% params$orgdb_pkg %||% params$orgdb %||% "org.Custom.eg.db"
prefix <- unlist(strsplit(pkg_name, ""))
# Raw taxonomy inputs (may be empty); resolved later against NCBI Taxonomy when possible
user_tax_id <- params$tax_id %||% ""
user_genus <- params$genus %||% ""
user_species <- params$species %||% ""
user_organism <- params$organism %||% params$scientific_name %||% params$species_name %||% ""

gene2go_path <- params$gene2go_tsv %||% params$gene2go %||% ""
gaf_path <- params$gaf %||% params$gaf_path %||% ""
eggnog_ann_path <- params$eggnog_annotations %||% params$eggnog_ann %||% params$annotations %||% ""
gene_info_path <- params$gene_info_tsv %||% params$gene_info %||% ""
protein_fasta <- params$protein_fasta %||% params$fasta %||% ""

# if (mode == "gene2go_tsv") {
#   if (!nzchar(gene2go_path) || !file.exists(gene2go_path)) fail("gene2go_tsv file not found")
# } else if (mode == "gaf") {
#   if (!nzchar(gaf_path) || !file.exists(gaf_path)) fail("gaf file not found")
# } else if (mode %in% c("eggnog_annotations", "eggnog", "emapper_annotations")) {
#   mode <- "eggnog_annotations"
#   if (!nzchar(eggnog_ann_path) || !file.exists(eggnog_ann_path)) fail("eggnog_annotations file not found")
# } else {
#   fail(paste0("Unsupported mode: ", mode))
# }

# required packages
suppressWarnings(suppressMessages({
  if (!requireNamespace("AnnotationForge", quietly=TRUE)) fail("AnnotationForge not installed")
  if (!requireNamespace("AnnotationDbi", quietly=TRUE)) fail("AnnotationDbi not installed")
  if (!requireNamespace("GO.db", quietly=TRUE)) fail("GO.db not installed")
}))

# ---- Resolve taxonomy (tax_id / genus / species) ----
# Priority:
#   1) If both tax_id and genus/species are provided: prefer tax_id
#   2) If tax_id is invalid (no match): fall back to genus/species
#   3) If genus/species are also invalid: use dummy values
# If NCBI lookup fails due to network issues, we keep user-provided values when possible.

dummy_tax_id <- suppressWarnings(as.integer(params$dummy_tax_id %||% 1))
if (is.na(dummy_tax_id) || dummy_tax_id < 1) dummy_tax_id <- 1L
dummy_genus <- params$dummy_genus %||% "Unknown"
dummy_species <- params$dummy_species %||% "sp"
use_ncbi <- as_bool(params$ncbi_taxonomy_lookup %||% TRUE)

norm_str <- function(x) {
  if (is.null(x) || length(x) == 0) return(NA_character_)
  x <- trimws(as.character(x)[1])
  if (!nzchar(x)) return(NA_character_)
  x
}

to_int <- function(x) {
  x <- norm_str(x)
  if (is.na(x)) return(NA_integer_)
  suppressWarnings({
    v <- as.integer(x)
    if (is.na(v)) return(NA_integer_)
    v
  })
}

# Allow inputs like genus="Homo sapiens" or organism="Homo sapiens"
split_sci_name <- function(genus, species, organism) {
  genus <- norm_str(genus)
  species <- norm_str(species)
  organism <- norm_str(organism)

  if (!is.na(genus) && is.na(species) && grepl("\\s", genus)) {
    parts <- strsplit(genus, "\\s+")[[1]]
    if (length(parts) >= 2) {
      genus <- parts[[1]]
      species <- paste(parts[-1], collapse=" ")
    }
  }
  if (is.na(genus) && !is.na(species) && grepl("\\s", species)) {
    parts <- strsplit(species, "\\s+")[[1]]
    if (length(parts) >= 2) {
      genus <- parts[[1]]
      species <- paste(parts[-1], collapse=" ")
    }
  }
  if ((is.na(genus) || is.na(species)) && !is.na(organism) && grepl("\\s", organism)) {
    parts <- strsplit(organism, "\\s+")[[1]]
    if (length(parts) >= 2) {
      if (is.na(genus)) genus <- parts[[1]]
      if (is.na(species)) species <- paste(parts[-1], collapse=" ")
    }
  }
  list(genus=genus, species=species)
}

ncbi_add_query_params <- function(u) {
  tool <- "GenomicExplorer"
  email <- norm_str(params$ncbi_email %||% params$email %||% "")
  api_key <- norm_str(params$ncbi_api_key %||% params$api_key %||% "")
  u <- paste0(u, "&tool=", utils::URLencode(tool, reserved=TRUE))
  if (!is.na(email)) u <- paste0(u, "&email=", utils::URLencode(email, reserved=TRUE))
  if (!is.na(api_key)) u <- paste0(u, "&api_key=", utils::URLencode(api_key, reserved=TRUE))
  u
}

fetch_json <- function(u, timeout_sec=15) {
  old_timeout <- getOption("timeout")
  options(timeout=timeout_sec)
  on.exit(options(timeout=old_timeout), add=TRUE)
  txt <- tryCatch(readLines(u, warn=FALSE, encoding="UTF-8"), error=function(e) NULL)
  if (is.null(txt) || length(txt) == 0) return(list(status="error", data=NULL))
  j <- tryCatch(jsonlite::fromJSON(paste(txt, collapse="
")), error=function(e) NULL)
  if (is.null(j)) return(list(status="error", data=NULL))
  list(status="ok", data=j)
}

ncbi_esearch_taxonomy <- function(term) {
  term <- norm_str(term)
  if (is.na(term)) return(list(status="nomatch", taxid=NA_integer_))
  u <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy&term=",
              utils::URLencode(term, reserved=TRUE), "&retmode=json")
  u <- ncbi_add_query_params(u)
  r <- fetch_json(u)
  if (r$status != "ok") return(list(status="error", taxid=NA_integer_))
  idlist <- r$data$esearchresult$idlist
  if (is.null(idlist) || length(idlist) < 1) return(list(status="nomatch", taxid=NA_integer_))
  tid <- suppressWarnings(as.integer(idlist[[1]]))
  if (is.na(tid)) return(list(status="nomatch", taxid=NA_integer_))
  list(status="ok", taxid=tid)
}

ncbi_taxonomy_summary <- function(taxid) {
  taxid <- suppressWarnings(as.integer(taxid))
  if (is.na(taxid) || taxid < 1) return(list(status="nomatch", scientific=NA_character_))
  u <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=taxonomy&id=", taxid, "&retmode=json")
  u <- ncbi_add_query_params(u)
  r <- fetch_json(u)
  if (r$status != "ok") return(list(status="error", scientific=NA_character_))
  rec <- NULL
  if (!is.null(r$data$result) && !is.null(r$data$result[[as.character(taxid)]])) {
    rec <- r$data$result[[as.character(taxid)]]
  }
  if (is.null(rec)) return(list(status="nomatch", scientific=NA_character_))
  sci <- rec$scientificname %||% rec$scientificName %||% rec$scientific_name %||% ""
  sci <- norm_str(sci)
  if (is.na(sci)) return(list(status="nomatch", scientific=NA_character_))
  list(status="ok", scientific=sci)
}

taxid_to_genus_species <- function(taxid) {
  s <- ncbi_taxonomy_summary(taxid)
  if (s$status != "ok") return(list(status=s$status, taxid=taxid, genus=NA_character_, species=NA_character_))
  parts <- strsplit(s$scientific, "\\s+")[[1]]
  genus <- if (length(parts) >= 1) parts[[1]] else NA_character_
  species <- if (length(parts) >= 2) paste(parts[-1], collapse=" ") else "sp"
  genus <- norm_str(genus)
  species <- norm_str(species)
  list(status="ok", taxid=as.integer(taxid), genus=genus, species=species)
}

name_to_taxid <- function(genus, species) {
  genus <- norm_str(genus); species <- norm_str(species)
  if (is.na(genus) || is.na(species)) return(list(status="nomatch", taxid=NA_integer_))
  sci <- paste(genus, species)
  # Use [Scientific Name] for more precise matching
  term <- paste0(sci, "[Scientific Name]")
  ncbi_esearch_taxonomy(term)
}

resolve_taxonomy <- function(user_tax_id, user_genus, user_species, user_organism,
                             use_ncbi=TRUE, dummy_tax_id=1L, dummy_genus="Unknown", dummy_species="sp") {
  sp <- split_sci_name(user_genus, user_species, user_organism)
  genus <- sp$genus
  species <- sp$species
  taxid <- to_int(user_tax_id)

  # Treat common placeholder values as "not provided"
  if (!is.na(genus) && !is.na(species) &&
      tolower(genus) == tolower(dummy_genus) && tolower(species) == tolower(dummy_species)) {
    genus <- NA_character_
    species <- NA_character_
  }
  if (!is.na(taxid) && taxid == dummy_tax_id && (is.na(genus) || is.na(species))) {
    taxid <- NA_integer_
  }

  has_tax <- !is.na(taxid)
  has_name <- !is.na(genus) && !is.na(species)

  ret_dummy <- function(src) list(tax_id=dummy_tax_id, genus=dummy_genus, species=dummy_species, source=src)

  if (has_tax) {
    if (use_ncbi) {
      r <- taxid_to_genus_species(taxid)
      if (r$status == "ok") {
        return(list(tax_id=r$taxid, genus=r$genus, species=r$species, source="ncbi_taxid"))
      } else if (r$status == "error") {
        return(list(tax_id=taxid,
                    genus=ifelse(has_name, genus, dummy_genus),
                    species=ifelse(has_name, species, dummy_species),
                    source="user_taxid_network_error"))
      }
      # nomatch: fall through
    } else {
      return(list(tax_id=taxid,
                  genus=ifelse(has_name, genus, dummy_genus),
                  species=ifelse(has_name, species, dummy_species),
                  source="user_taxid"))
    }

    # taxid nomatch
    if (has_name) {
      if (use_ncbi) {
        t2 <- name_to_taxid(genus, species)
        if (t2$status == "ok") {
          r2 <- taxid_to_genus_species(t2$taxid)
          if (r2$status == "ok") {
            return(list(tax_id=r2$taxid, genus=r2$genus, species=r2$species, source="ncbi_name_fallback"))
          } else if (r2$status == "error") {
            return(list(tax_id=t2$taxid, genus=genus, species=species, source="user_name_network_error"))
          }
        } else if (t2$status == "error") {
          return(list(tax_id=taxid, genus=genus, species=species, source="user_inputs_network_error"))
        }
        return(ret_dummy("dummy_no_match"))
      } else {
        return(list(tax_id=dummy_tax_id, genus=genus, species=species, source="user_name_fallback_no_lookup"))
      }
    } else {
      return(ret_dummy("dummy_taxid_no_match"))
    }
  } else {
    # no taxid
    if (has_name) {
      if (use_ncbi) {
        t2 <- name_to_taxid(genus, species)
        if (t2$status == "ok") {
          r2 <- taxid_to_genus_species(t2$taxid)
          if (r2$status == "ok") {
            return(list(tax_id=r2$taxid, genus=r2$genus, species=r2$species, source="ncbi_name"))
          } else if (r2$status == "error") {
            return(list(tax_id=t2$taxid, genus=genus, species=species, source="user_name_network_error"))
          }
        } else if (t2$status == "error") {
          return(list(tax_id=dummy_tax_id, genus=genus, species=species, source="user_name_network_error"))
        }
        return(ret_dummy("dummy_name_no_match"))
      } else {
        return(list(tax_id=dummy_tax_id, genus=genus, species=species, source="user_name_no_lookup"))
      }
    } else {
      return(ret_dummy("dummy_no_input"))
    }
  }
}

tx <- resolve_taxonomy(user_tax_id, user_genus, user_species, user_organism,
                       use_ncbi=use_ncbi,
                       dummy_tax_id=dummy_tax_id,
                       dummy_genus=dummy_genus,
                       dummy_species=dummy_species)

tax_id <- suppressWarnings(as.integer(tx$tax_id))
if (is.na(tax_id) || tax_id < 1) tax_id <- dummy_tax_id
genus <- tx$genus %||% dummy_genus
species <- tx$species %||% dummy_species

cat("[rnaseq_build_orgdb] taxonomy resolved: tax_id=", tax_id,
    " genus=", genus, " species=", species,
    " source=", (tx$source %||% "unknown"), "
", sep="")

parse_fasta_ids <- function(path) {
  if (!nzchar(path) || !file.exists(path)) return(character(0))
  x <- readLines(path, warn=FALSE, encoding="UTF-8")
  hdr <- x[startsWith(x, ">")]
  if (length(hdr) == 0) return(character(0))
  # take the first token after '>' (up to whitespace); also split by '|' and keep the first non-empty token
  ids <- sub("^>\\s*", "", hdr)
  ids <- sub("\\s.*$", "", ids)
  ids <- vapply(strsplit(ids, "\\|"), function(v) {
    v <- v[nzchar(v)]
    if (length(v) == 0) return(NA_character_)
    v[[1]]
  }, character(1))
  ids <- ids[!is.na(ids) & nzchar(ids)]
  unique(ids)
}

read_gene_info <- function(path, ids_fallback=character(0)) {
  if (nzchar(path) && file.exists(path)) {
    dt <- tryCatch(data.table::fread(path, sep="\t", header=TRUE, data.table=FALSE), error=function(e) NULL)
    if (!is.null(dt) && nrow(dt) > 0) {
      # Accept common column names; otherwise use first column as gene id
      cn <- colnames(dt)
      gene_col <- if ("gene" %in% cn) "gene" else if ("gene_id" %in% cn) "gene_id" else cn[[1]]
      sym_col <- if ("symbol" %in% cn) "symbol" else if ("SYMBOL" %in% cn) "SYMBOL" else NULL
      name_col <- if ("genename" %in% cn) "genename" else if ("GENENAME" %in% cn) "GENENAME" else NULL
      gid <- as.character(dt[[gene_col]])
      gid <- gid[!is.na(gid) & nzchar(gid)]
      sym <- if (!is.null(sym_col)) as.character(dt[[sym_col]]) else gid
      gnm <- if (!is.null(name_col)) as.character(dt[[name_col]]) else gid
      out <- data.frame(GID=gid, SYMBOL=sym, GENENAME=gnm, stringsAsFactors=FALSE)
      out <- out[!duplicated(out$GID), , drop=FALSE]
      return(out)
    }
  }
  # fallback
  ids <- unique(ids_fallback)
  ids <- ids[!is.na(ids) & nzchar(ids)]
  if (length(ids) == 0) {
    return(data.frame(GID=character(0), SYMBOL=character(0), GENENAME=character(0), stringsAsFactors=FALSE))
  }
  data.frame(GID=ids, SYMBOL=ids, GENENAME=ids, stringsAsFactors=FALSE)
}

normalize_go_id <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x == ""] <- NA_character_
  # ensure prefix
  x <- ifelse(is.na(x), NA_character_, ifelse(grepl("^GO:", x), x, paste0("GO:", x)))
  x
}

infer_ontology_from_go <- function(go_ids) {
  go_ids <- normalize_go_id(go_ids)
  ont <- rep(NA_character_, length(go_ids))
  ok <- !is.na(go_ids) & nzchar(go_ids)
  if (any(ok)) {
    # GO.db::Ontology returns BP/MF/CC
    ont[ok] <- suppressWarnings(as.character(GO.db::Ontology(go_ids[ok])))
  }
  ont <- toupper(ont)
  ont[!(ont %in% c("BP","MF","CC"))] <- NA_character_
  ont
}

read_gene2go_tsv <- function(path) {
  dt <- tryCatch(data.table::fread(path, sep="\t", header=TRUE, data.table=FALSE), error=function(e) NULL)
  if (is.null(dt) || nrow(dt) < 1) {
    # try no-header
    dt <- tryCatch(data.table::fread(path, sep="\t", header=FALSE, data.table=FALSE), error=function(e) NULL)
  }
  if (is.null(dt) || nrow(dt) < 1) fail("Failed to read gene2go TSV")

  # Normalize to at least two columns
  if (ncol(dt) < 2) fail("gene2go TSV must have at least 2 columns: gene_id, go_id")

  gid <- as.character(dt[[1]])
  go <- normalize_go_id(dt[[2]])

  evidence <- rep("IEA", length(gid))

  if (ncol(dt) >= 3) {
    ev <- as.character(dt[[3]])
    ev[is.na(ev) | !nzchar(ev)] <- "IEA"
    evidence <- ev
  }
  out <- data.frame(GID=gid, GO=go, EVIDENCE=evidence, stringsAsFactors=FALSE)
  out <- out[!is.na(out$GID) & nzchar(out$GID) & !is.na(out$GO) & nzchar(out$GO), , drop=FALSE]
  out
}

read_gaf <- function(path) {
  # GAF 2.x: tab-delimited, comments start with '!'
  x <- readLines(path, warn=FALSE, encoding="UTF-8")
  x <- x[!startsWith(x, "!")]
  if (length(x) == 0) fail("GAF contains no data lines")
  # fread on textConnection
  dt <- tryCatch(data.table::fread(text=paste(x, collapse="\n"), sep="\t", header=FALSE, data.table=FALSE, fill=TRUE, quote=""),
                 error=function(e) NULL)
  if (is.null(dt) || nrow(dt) < 1) fail("Failed to parse GAF")
  # columns: 1 DB, 2 DB_Object_ID, 3 DB_Object_Symbol, 4 Qualifier, 5 GO_ID, 6 DB:Reference, 7 Evidence, 9 Aspect
  if (ncol(dt) < 9) {
    # still accept but require at least 7 for evidence, 5 for GO
    if (ncol(dt) < 7) fail("GAF must have at least 7 columns")
  }
  gid <- as.character(dt[[2]])
  go <- normalize_go_id(dt[[5]])
  evidence <- as.character(dt[[7]])
  evidence[is.na(evidence) | !nzchar(evidence)] <- "IEA"
  out <- data.frame(GID=gid, GO=go, EVIDENCE=evidence, stringsAsFactors=FALSE)
  out <- out[!is.na(out$GID) & nzchar(out$GID) & !is.na(out$GO) & nzchar(out$GO), , drop=FALSE]
  out
}

read_eggnog_annotations <- function(path) {
  # eggNOG-mapper output (typically *.emapper.annotations)
  # Header line is often commented: "#query\t...\tGOs\t..."
  x <- readLines(path, warn=FALSE, encoding="UTF-8")
  x <- x[nzchar(trimws(x))]
  if (length(x) == 0) fail("EggNOG annotations file is empty")

  # find header line
  hdr_idx <- which(grepl("^#?query\\b", x, ignore.case=TRUE))
  if (length(hdr_idx) == 0) {
    fail("Could not locate a header line starting with 'query' in eggNOG annotations")
  }
  hdr_line <- x[[hdr_idx[[1]]]]
  hdr_line <- sub("^#", "", hdr_line)
  hdr <- strsplit(hdr_line, "\\t", fixed=FALSE)[[1]]
  hdr_l <- tolower(hdr)

  q_col <- if ("query" %in% hdr_l) which(hdr_l == "query")[[1]] else 1L
  # locate GOs column
  go_keys <- c("gos", "go", "go_terms", "go terms", "go_term")
  go_col <- NA_integer_
  for (k in go_keys) {
    if (k %in% hdr_l) { go_col <- which(hdr_l == k)[[1]]; break }
  }
  if (is.na(go_col)) fail("Could not locate 'GOs' column in eggNOG annotations header")

  # parse data lines (after header), skipping comments
  if (hdr_idx[[1]] >= length(x)) fail("EggNOG annotations contains no data lines after header")
  data_lines <- x[(hdr_idx[[1]] + 1):length(x)]
  data_lines <- data_lines[!startsWith(data_lines, "#")]
  if (length(data_lines) == 0) fail("EggNOG annotations contains no non-comment data lines")

  pick_id <- function(s) {
    s <- trimws(as.character(s))
    s <- sub("\\s.*$", "", s)              # first whitespace token
    # also split by '|' and take first non-empty
    ss <- unlist(strsplit(s, "\\|", fixed=FALSE))
    ss <- ss[nzchar(ss)]
    if (length(ss) == 0) return(NA_character_)
    ss[[1]]
  }

  go_re <- gregexpr("GO:[0-9]+", data_lines, perl=TRUE)

  out_pairs <- list()
  k <- 0L
  for (ln in data_lines) {
    parts <- strsplit(ln, "\\t", fixed=FALSE)[[1]]
    if (length(parts) < max(q_col, go_col)) next
    gid <- pick_id(parts[[q_col]])
    if (is.na(gid) || !nzchar(gid)) next
    gos_raw <- parts[[go_col]]
    if (is.na(gos_raw) || !nzchar(trimws(gos_raw)) || gos_raw %in% c("-", "NA", "NaN")) next

    gos <- regmatches(gos_raw, gregexpr("GO:[0-9]+", gos_raw, perl=TRUE))[[1]]
    if (length(gos) == 0) next
    gos <- unique(normalize_go_id(gos))
    gos <- gos[!is.na(gos) & nzchar(gos)]
    if (length(gos) == 0) next
    k <- k + 1L
    out_pairs[[k]] <- data.frame(GID=rep(gid, length(gos)), GO=gos, EVIDENCE=rep("IEA", length(gos)), stringsAsFactors=FALSE)
  }

  if (length(out_pairs) == 0) {
    fail("No GO terms could be extracted from eggNOG annotations (check the 'GOs' column)")
  }
  out <- data.table::rbindlist(out_pairs, use.names=TRUE, fill=TRUE)
  out <- unique(as.data.frame(out))
  out <- out[!is.na(out$GID) & nzchar(out$GID) & !is.na(out$GO) & nzchar(out$GO), , drop=FALSE]
  out
}

ids_fa <- parse_fasta_ids(protein_fasta)

go_df <- NULL
if (grepl(".(tsv|txt)$", gene2go_path, perl = T)) {
  go_df <- read_gene2go_tsv(gene2go_path)
} else if (grepl(".gaf(|.gz)$", gene2go_path, perl = T)) {
  go_df <- read_gaf(gene2go_path)
} else {
  go_df <- read_eggnog_annotations(gene2go_path)
}

# if (mode == "gene2go_tsv") {
#   go_df <- read_gene2go_tsv(gene2go_path)
# } else if (mode == "gaf") {
#   go_df <- read_gaf(gaf_path)
# } else {
#   go_df <- read_eggnog_annotations(eggnog_ann_path)
# }

if (nrow(go_df) < 1) fail("No valid gene-to-GO records after parsing")

if (length(ids_fa) > 0 && as_bool(params$filter_by_fasta %||% TRUE)) {
  before <- nrow(go_df)
  go_df <- go_df[go_df$GID %in% ids_fa, , drop=FALSE]
  cat("[rnaseq_build_orgdb] filtered by FASTA ids: ", before, " -> ", nrow(go_df), "\n", sep="")
  if (nrow(go_df) < 1) fail("All mappings were filtered out by FASTA IDs")
}

# Deduplicate
go_df <- unique(go_df[, c("GID","GO","EVIDENCE")])

all_ids <- unique(go_df$GID)
if (length(ids_fa) > 0) all_ids <- unique(c(all_ids, ids_fa))

gene_info <- read_gene_info(gene_info_path, ids_fallback=all_ids)
if (nrow(gene_info) < 1) {
  # at minimum, require gene ids present
  gene_info <- data.frame(GID=all_ids, SYMBOL=all_ids, GENENAME=all_ids, stringsAsFactors=FALSE)
}

# Keep gene_info only for IDs present (or FASTA)
gene_info <- gene_info[gene_info$GID %in% all_ids, , drop=FALSE]
gene_info <- gene_info[!duplicated(gene_info$GID), , drop=FALSE]

# Many Bioconductor tools (including GOstats) expect an ENTREZID column.
# For custom organisms without NCBI Gene IDs, we alias ENTREZID to our internal GID.
if (!("ENTREZID" %in% colnames(gene_info))) {
  gene_info$ENTREZID <- gene_info$GID
}
# Order columns with GID first (required by makeOrgPackage)
base_cols <- c("GID","ENTREZID","SYMBOL","GENENAME")
extra_cols <- setdiff(colnames(gene_info), base_cols)
gene_info <- gene_info[, c(base_cols, extra_cols), drop=FALSE]

cat("[rnaseq_build_orgdb] genes=", length(unique(go_df$GID)), " go_records=", nrow(go_df), "\n", sep="")

# Prepare output library
rlib <- file.path(out_dir, "Rlib")
dir.create(rlib, recursive=TRUE, showWarnings=FALSE)

# Make package in out_dir
oldwd <- getwd()
setwd(out_dir)
on.exit(setwd(oldwd), add=TRUE)

maintainer <- params$maintainer %||% "GenomicExplorer <noreply@example.com>"
author <- params$author %||% "GenomicExplorer"
pkg_version <- params$pkg_version %||% "0.1.0"

cat("[rnaseq_build_orgdb] building package: ", pkg_name, "\n", sep="")

# IMPORTANT: pass goTable="go" so AnnotationForge generates GOALL tables and ontology-specific mappings.
# This is required for downstream tools like GOstats.
AnnotationForge::makeOrgPackage(gene_info=gene_info, go=go_df,
                                version=pkg_version,
                                maintainer=maintainer,
                                author=author,
                                outputDir=out_dir,
                                tax_id=tax_id,
                                # genus=genus,
                                # species=species,
                                genus = prefix[1],
                                species = prefix[2:length(prefix)],
                                goTable="go")

pkg_dir <- file.path(out_dir, pkg_name)
if (!dir.exists(pkg_dir)) {
  # some AnnotationForge versions create a directory based on the package name returned
  # try to find a single org.*.eg.db directory
  cands <- list.dirs(out_dir, full.names=TRUE, recursive=FALSE)
  cands <- cands[basename(cands) != "Rlib"]
  cands <- cands[grepl("^org\\..+\\.eg\\.db$", basename(cands))]
  if (length(cands) == 1) {
    pkg_dir <- cands[[1]]
    pkg_name <- basename(pkg_dir)
  }
}
if (!dir.exists(pkg_dir)) fail(paste0("Package directory not found after build: ", pkg_dir))

# Build source tarball then install to out_dir/Rlib (keeps environment clean)
r_bin <- file.path(R.home("bin"), "R")
cat("[rnaseq_build_orgdb] R CMD build\n")
res_build <- system2(r_bin, c("CMD", "build", shQuote(pkg_dir)), stdout=TRUE, stderr=TRUE)
cat(paste(res_build, collapse="\n"), "\n")

tarballs <- list.files(out_dir, pattern=paste0("^", gsub("\\.", "\\\\.", pkg_name), "_.*\\.tar\\.gz$"), full.names=TRUE)
if (length(tarballs) < 1) {
  # fallback: any tar.gz created
  tarballs <- list.files(out_dir, pattern="\\.tar\\.gz$", full.names=TRUE)
}
if (length(tarballs) < 1) fail("Failed to locate built source tar.gz")
tarball <- tarballs[[1]]

cat("[rnaseq_build_orgdb] installing to the R environment", "\n", sep="")
tryCatch({
  install.packages(tarball, repos=NULL, type="source", quiet=F)
}, error=function(e) fail(paste("install.packages to the R environment failed:", e$message)))

cat("[rnaseq_build_orgdb] installing to Rlib=", rlib, "\n", sep="")
tryCatch({
  install.packages(tarball, repos=NULL, type="source", lib=rlib, quiet=TRUE)
}, error=function(e) fail(paste("install.packages failed:", e$message)))

# Summary artifacts
summary_tsv <- file.path(out_dir, "orgdb_build_summary.tsv")
kv <- data.frame(
  key=c("orgdb_pkg", "r_libs", "mode", "gene_records", "unique_genes", "unique_go", "tax_id", "genus", "species", "source_mapping"),
  value=c(
    pkg_name,
    rlib,
    mode,
    nrow(go_df),
    length(unique(go_df$GID)),
    length(unique(go_df$GO)),
    ifelse(is.na(tax_id), "", as.character(tax_id)),
    genus,
    species,
    ifelse(mode=="gene2go_tsv", gene2go_path, ifelse(mode=="gaf", gaf_path, eggnog_ann_path))
  ),
  stringsAsFactors=FALSE
)
write.table(kv, summary_tsv, sep="\t", quote=FALSE, row.names=FALSE)

art <- list(
  orgdb_pkg=pkg_name,
  r_libs=rlib,
  summary_tsv=summary_tsv,
  mode=mode,
  gene_records=nrow(go_df),
  unique_genes=length(unique(go_df$GID)),
  unique_go=length(unique(go_df$GO))
)
writeLines(jsonlite::toJSON(art, auto_unbox=TRUE, pretty=TRUE), file.path(out_dir, "artifacts.json"))

cat("[rnaseq_build_orgdb] done\n")
