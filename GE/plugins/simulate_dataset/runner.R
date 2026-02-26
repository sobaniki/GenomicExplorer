#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(jsonlite)
  library(data.table)
  library(AlphaSimR)
  library(VariantAnnotation)
  library(GenomicRanges)
  library(IRanges)
  library(Biostrings)
  library(S4Vectors)
  library(updog)
  library(Rsamtools)
})

args <- commandArgs(trailingOnly=TRUE)
get_arg <- function(flag, default=NULL) {
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) return(args[i+1])
  default
}

params_path <- get_arg("--params")
out_dir <- get_arg("--out")
if (is.null(params_path) || is.null(out_dir)) stop("Usage: --params params.json --out out_dir")

dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
log_file <- file.path(out_dir, "run.log")
sink(log_file, split=TRUE)

cat("[simulate_dataset] start\n")
cat("[simulate_dataset] params_path=", params_path, "\n")
cat("[simulate_dataset] out_dir=", out_dir, "\n")

p <- fromJSON(params_path)

`%||%` <- function(a, b) if (!is.null(a)) a else b

seed <- as.integer(p$seed %||% 1)
set.seed(seed)

species <- as.character(p$species_template %||% "rice")
pop_type <- as.character(p$population_type %||% "F2")
n_ind <- as.integer(p$n_individuals %||% 200)
n_markers <- as.integer(p$n_markers %||% 2000)
missing_rate <- as.numeric(p$missing_rate %||% 0.0)
geno_error <- as.numeric(p$genotype_error %||% 0.0)

# Population options (used mainly by the internal simulator backend)
ril_selfing <- as.integer(p$ril_selfing %||% 6)
bc_generation <- as.integer(p$bc_generation %||% 1)
bc_recurrent <- as.character(p$bc_recurrent %||% "P1") # P1 or P2
magic_n_founders <- as.integer(p$magic_n_founders %||% 8)
magic_n_families <- as.integer(p$magic_n_families %||% 1)
nam_n_families <- as.integer(p$nam_n_families %||% 5)
nam_ril_selfing <- as.integer(p$nam_ril_selfing %||% 6)

trait_num <- as.numeric(p$trait_number)
geno_mean <- as.numeric(rep(p$geno_mean, trait_num))
geno_var <- as.numeric(rep(p$geno_var, trait_num))
meanDD <- as.numeric(p$meanDD)
varDD <- as.numeric(p$varDD)
relAA <- as.numeric(p$relAA)
pheno_missing <- as.numeric(p$pheno_missing)

#Real chromosomes from FASTA
fasta_file <- NULL
if ((!(is.null(p$fasta_file)) || p$fasta_file != "")) fasta_file <- p$fasta_file

#trait_arch <- as.character(p$trait_architecture %||% "qtl")  # qtl/poly/both

if (p$trait_name == "" | is.null(p$trait_name)) {
  trait_name <- paste0("trait",
                       1:trait_num)
} else if (trait_num == length(strsplit(p$trait_name, ",")[[1]])) {
  trait_name <- p$trait_name
} else {
  trait_name <- paste0(strsplit(p$trait_name, ",")[[1]],
                       1:trait_num)
}

h2 <- as.numeric(rep(p$h2, trait_num))
H2 <- as.numeric(rep(p$H2, trait_num))
n_qtl_per_chr <- as.numeric(p$n_qtl)

seq <- as.numeric(p$seq)
bias <- as.numeric(p$bias)
od <- as.numeric(p$od)

# Genome template (simple, bp-based). We keep cM as a rough proportional proxy.
tmpl <- list(
  rice = list(n_chr=12, chr_len_mb=rep(30, 12)),
  tomato = list(n_chr=12, chr_len_mb=c(95, 55, 65, 66, 63, 50, 67, 61, 68, 65, 54, 61))
)

custom_n_chr <- p$custom_n_chr
custom_chr_len_mb <- p$custom_chr_len_mb

if (species == "custom") {
  n_chr <- as.integer(custom_n_chr %||% 12)
  if (!is.null(custom_chr_len_mb)) {
    if (is.list(custom_chr_len_mb)) custom_chr_len_mb <- unlist(custom_chr_len_mb)
    chr_len_mb <- as.numeric(custom_chr_len_mb)
    if (length(chr_len_mb) != n_chr) {
      stop("custom_chr_len_mb must have length custom_n_chr")
    }
  } else {
    chr_len_mb <- rep(50, n_chr)
  }
} else {
  if (!species %in% names(tmpl)) {
    cat("[simulate_dataset] WARN: unknown species_template=", species, " -> use rice\n")
    species <- "rice"
  }
  n_chr <- tmpl[[species]]$n_chr
  chr_len_mb <- tmpl[[species]]$chr_len_mb
}

chr_names <- sprintf("chr%02d", seq_len(n_chr))
chr_len_bp <- as.integer(round(chr_len_mb * 1e6))

if (nchar(fasta_file) > 0) {
  fa <- Rsamtools::scanFa(fasta_file)
  chr_names <- names(fa)
  chr_len_bp <- fa@ranges@width
  n_chr <- length(chr_names)
  #n_qtl_per_chr <- rep(n_qtl_per_chr, n_chr)
  species <- "custom"
}

cat("[simulate_dataset] species=", species, " n_chr=", n_chr, "\n")
cat("[simulate_dataset] pop_type=", pop_type, " n_ind=", n_ind, " n_markers=", n_markers, "\n")
cat("[simulate_dataset] missing_rate=", missing_rate, " geno_error=", geno_error, "\n")
#cat("[simulate_dataset] trait_arch=", trait_arch, " trait=", trait_name, " h2=", h2, " n_qtl=", n_qtl, "\n")

# -------------------------------
# Marker map
# -------------------------------

# Allocate markers proportional to chr length
prop <- chr_len_bp / sum(chr_len_bp)
markers_per_chr <- pmax(1L, as.integer(round(prop * n_markers)))
diffm <- n_markers - sum(markers_per_chr)
if (diffm != 0) {
  # adjust largest chromosomes
  ord <- order(chr_len_bp, decreasing=TRUE)
  for (k in seq_len(abs(diffm))) {
    i <- ord[((k-1) %% length(ord)) + 1]
    markers_per_chr[i] <- markers_per_chr[i] + sign(diffm)
    if (markers_per_chr[i] < 1) markers_per_chr[i] <- 1
  }
}

marker_id <- character(0)
chr <- character(0)
pos_bp <- integer(0)
pos_cM <- numeric(0)
for (i in seq_len(n_chr)) {
  m <- markers_per_chr[i]
  bp <- sort(sample.int(chr_len_bp[i], m, replace=FALSE))
  # crude cM: assume 1 cM per 1 Mb
  cm <- bp / 1e6
  marker_id <- c(marker_id, sprintf("m%02d_%06d", i, seq_len(m)))
  chr <- c(chr, rep(chr_names[i], m))
  pos_bp <- c(pos_bp, bp)
  pos_cM <- c(pos_cM, cm)
}

map_dt <- data.table(marker_id=marker_id, chr=chr, pos_bp=pos_bp, pos_cM=pos_cM)

# -------------------------------
# Genotype simulation
# -------------------------------

use_alphasimr <- FALSE
if (!is.null(p$use_alphasimr)) use_alphasimr <- as.logical(p$use_alphasimr)

# Export options
export_qtl2 <- TRUE
export_rqtl <- TRUE
if (!is.null(p$export_qtl2)) export_qtl2 <- as.logical(p$export_qtl2)
if (!is.null(p$export_rqtl)) export_rqtl <- as.logical(p$export_rqtl)

ploidy <- as.integer(p$ploidy)
generate_polyploid_depth <- as.logical(p$generate_polyploid_depth %||% T)
depth_mean <- as.numeric(p$depth_mean %||% 20)
depth_nb_size <- as.numeric(p$depth_nb_size %||% 30)  # larger -> less overdispersion
allele_bias_logodds <- as.numeric(p$allele_bias_logodds %||% 0.0)  # log-odds shift (+ favors ALT)
seq_error <- as.numeric(p$seq_error %||% 0.01)

geno_mat <- NULL
truth_qtl <- data.table()

geno_merge <- NULL
mergepop <- NULL

# Create a simple founder population with biallelic markers.
# We use a rough approximation: per chr, markers placed on a genetic map in Morgans.
gen_map <- split(map_dt$pos_cM / 100, map_dt$chr)  # Morgans

# NOTE: AlphaSimR typically needs founder haplotypes; we generate them with runMacs if available.
# To avoid external dependencies, we fallback to internal simulator if runMacs is not available.
if (!exists("runMacs", where=asNamespace("AlphaSimR"), inherits=FALSE)) {
  cat("[simulate_dataset] WARN: AlphaSimR::runMacs not available in this installation; fallback.\n")
  use_alphasimr <- FALSE
} else {
  pop_type_u0 <- toupper(as.character(pop_type))
  n_ind_macs <- max(50L, n_ind)
  inbred_macs <- if (pop_type_u0 %in% c("FOUNDER_OUTBRED","F1_OUTBRED","S1_OUTBRED")) FALSE else TRUE
  
  founder <- AlphaSimR::runMacs(nInd = n_ind_macs, 
                                nChr = n_chr, 
                                segSites = markers_per_chr + n_qtl_per_chr * n_chr, 
                                inbred = inbred_macs, 
                                species = "GENERIC",
                                ploidy = ploidy)
  SP <- SimParam$new(founder)
  

  SP$addTraitADE(
    nQtlPerChr = n_qtl_per_chr,
    mean = geno_mean,
    var = geno_var,
    meanDD = meanDD,
    varDD = varDD,
    relAA = relAA,
    corA = NULL,
    corDD = NULL,
    corAA = NULL,
    useVarA = TRUE,
    gamma = FALSE,
    shape = 1,
    force = FALSE,
    name = trait_name
  )
  
  qtl_info <- SP$traits
  geno_map <- SP$genMap
  
  SP$addSnpChip(nSnp = sapply(gen_map, length))
  
  pop <- newPop(founder, 
                simParam = SP)
  
  pop <- setPheno(pop, 
                  h2 = h2,
                  H2 = H2,
                  simParam = SP)
  pop_type_u <- toupper(as.character(pop_type))
  if (pop_type_u %in% c("FOUNDER_INBRED", "FOUNDER_OUTBRED")) {
    pop1 <- pop
    ped_dt_internal <- data.table(sample_id = pop1@id,
                                  sire = NA,
                                  dam = NA,
                                  generation = "Founder",
                                  family = 0)
  } else if (pop_type_u == "F1_OUTBRED") {
    pop0 <- randCross(pop, 
                      nCrosses= 2, 
                      nProgeny = 1, 
                      simParam = SP)
    pop0 <- setPheno(pop0, 
                     h2 = h2,
                     H2 = H2,
                     simParam = SP)
    pop1 <- randCross(pop0, 
                      nCrosses= 1, 
                      nProgeny = n_ind_macs, 
                      simParam = SP)
    pop1 <- setPheno(pop1, 
                     h2 = h2,
                     H2 = H2,
                     simParam = SP)
    ped_dt_internal <- data.table(sample_id = c(pop1@mother[1], pop1@father[1], pop1@id),
                                  sire = c(NA, NA, pop1@father),
                                  dam = c(NA, NA, pop1@mother),
                                  generation = c("Founder", "Founder", rep("F1", length(pop1@id))),
                                  family = c(0, 0, rep(1, length(pop1@id))))
    mergepop <- mergePops(list(pop0, pop1))
  } else if (pop_type_u == "S1_OUTBRED") {
    pop0 <- randCross(pop, 
                      nCrosses= 1, 
                      nProgeny = 1, 
                      simParam = SP)
    pop0 <- setPheno(pop0, 
                     h2 = h2,
                     H2 = H2,
                     simParam = SP)
    pop1 <- self(pop0,
                 nProgeny = n_ind_macs,
                 keepParents = F, 
                 simParam = SP)
    pop1 <- setPheno(pop1, 
                     h2 = h2,
                     H2 = H2,
                     simParam = SP)
    ped_dt_internal <- data.table(sample_id = c(pop0@mother[1], pop0@father[1], pop0@id, pop1@id),
                                  sire = c(NA, NA, pop0@father, pop1@father),
                                  dam = c(NA, NA, pop0@mother, pop1@mother),
                                  generation = c("Founder", "Founder", "F1", rep("S1", length(pop1@id))),
                                  family = c(0, 0, 0, rep(1, length(pop1@id))))
    mergepop <- mergePops(list(pop0, pop1))
  } else if (pop_type_u == "F2") {
    pop0 <- randCross(pop, 
                      nCrosses= 1, 
                      nProgeny = 1, 
                      simParam = SP)
    pop0 <- setPheno(pop0, 
                     h2 = h2,
                     H2 = H2,
                     simParam = SP)
    pop1 <- self(pop0,
                 nProgeny = n_ind_macs,
                 keepParents = T, 
                 simParam = SP)
    pop1 <- setPheno(pop1, 
                     h2 = h2,
                     H2 = H2,
                     simParam = SP)
    ped_dt_internal <- data.table(sample_id = c(pop0@mother, pop0@father, pop0@id, pop1@id),
                                  sire = c(NA, NA, pop0@father, rep(pop0@id, length(pop1@id))),
                                  dam = c(NA, NA, pop0@mother, rep(pop0@id, length(pop1@id))),
                                  generation = c("Founder", "Founder", "F1", rep("F2", length(pop1@id))),
                                  family = c(0, 0, 1, rep(1, length(pop1@id))))
    mergepop <- mergePops(list(pop, pop1))
  } else if (toupper(pop_type) == "DH") {
    pop0 <- randCross(pop, 
                      nCrosses= 1, 
                      nProgeny = 1, 
                      simParam = SP)
    pop0 <- setPheno(pop0, 
                     h2 = h2,
                     H2 = H2,
                     simParam = SP)
    pop1 <- makeDH(pop0,
                   nDH = n_ind_macs,
                   useFemale = T,
                   keepParents = T, 
                   simParam = SP)
    pop1 <- setPheno(pop1, 
                     h2 = h2,
                     H2 = H2,
                     simParam = SP)
    ped_dt_internal <- data.table(sample_id = c(pop0@mother, pop0@father, pop0@id, pop1@id),
                                  sire = c(NA, NA, pop0@father, rep(pop0@id, length(pop1@id))),
                                  dam = c(NA, NA, pop0@mother, rep(pop0@id, length(pop1@id))),
                                  generation = c("Founder", "Founder", "F1", rep("DH", length(pop1@id))),
                                  family = c(0, 0, 1, rep(1, length(pop1@id))))
    mergepop <- mergePops(list(pop, pop1))
  } else if (toupper(pop_type) == "RIL") {
    n_self <- as.integer(p$ril_selfing %||% 6)
    pop0 <- randCross(pop, 
                      nCrosses= 1, 
                      nProgeny = 1, 
                      simParam = SP)
    pop0 <- setPheno(pop0, 
                     h2 = h2,
                     H2 = H2,
                     simParam = SP)
    pop1 <- self(pop0,
                 nProgeny = n_ind_macs,
                 keepParents = T, 
                 simParam = SP)
    pop1 <- setPheno(pop1, 
                     h2 = h2,
                     H2 = H2,
                     simParam = SP)
    ped_dt_internal <- data.table(sample_id = c(pop0@mother, pop0@father, pop0@id, pop1@id),
                                  sire = c(NA, NA, pop0@father, rep(pop0@id, length(pop1@id))),
                                  dam = c(NA, NA, pop0@mother, rep(pop0@id, length(pop1@id))),
                                  generation = c("Founder", "Founder", "F1", rep("F2", length(pop1@id))),
                                  family = c(0, 0, 1, rep(1, length(pop1@id))))
    for (cyc1 in 3:n_self) {
      pop1 <- self(pop1,
                   nProgeny = 1,
                   keepParents = F, 
                   simParam = SP)
      pop1 <- setPheno(pop1, 
                       h2 = h2,
                       H2 = H2,
                       simParam = SP)
      ped_dt_internal <- rbind(ped_dt_internal,
                               data.table(sample_id = pop1@id,
                                          sire = pop1@father,
                                          dam = pop1@mother,
                                          generation = rep(paste0("F", cyc1), length(pop1@id)),
                                          family = rep(1, length(pop1@id))))
    }
    mergepop <- mergePops(list(pop, pop1))
  } else if (toupper(pop_type) == "BC") {
    pop0 <- randCross(pop, 
                      nCrosses= 1, 
                      nProgeny = 1, 
                      simParam = SP)
    pop0 <- setPheno(pop0, 
                     h2 = h2,
                     H2 = H2,
                     simParam = SP)
    ped_dt_internal <- data.table(sample_id = c(pop0@mother, pop0@father, pop0@id),
                                  sire = c(NA, NA, pop0@father),
                                  dam = c(NA, NA, pop0@mother),
                                  generation = c("Founder", "Founder", "F1"),
                                  family = c(0, 0, 1))
    if (bc_recurrent == "P1") {
      pop1 <- makeCross2(females = pop,
                         males = pop0,
                         crossPlan = as.matrix(data.frame(recurrent = pop0@mother,
                                                          cross = pop0@iid[1])),
                         nProgeny = n_ind_macs, 
                         simParam = SP)
      pop1 <- setPheno(pop1, 
                       h2 = h2,
                       H2 = H2,
                       simParam = SP)
    } else {
      pop1 <- makeCross2(females = pop0,
                         males = pop,
                         crossPlan = as.matrix(data.frame(cross = pop0@iid[1],
                                                          recurrent = pop0@father)),
                         nProgeny = n_ind_macs, 
                         simParam = SP)
      pop1 <- setPheno(pop1, 
                       h2 = h2,
                       H2 = H2,
                       simParam = SP)
    }
    ped_dt_internal <- rbind(ped_dt_internal,
                             data.table(sample_id = pop1@id,
                                        sire = pop1@father,
                                        dam = pop1@mother,
                                        generation = rep("BC1", length(pop1@id)),
                                        family = rep(1, length(pop1@id))))
    if (bc_generation > 1) {
      for (cyc1 in 2:bc_generation) {
        if (bc_recurrent == "P1") {
          pop1 <- makeCross2(females = pop,
                             males = pop1,
                             crossPlan = as.matrix(data.frame(recurrent = pop0@mother,
                                                              cross = pop1@iid[1])),
                             nProgeny = n_ind_macs, 
                             simParam = SP)
          pop1 <- setPheno(pop1, 
                           h2 = h2,
                           H2 = H2,
                           simParam = SP)
        } else {
          pop1 <- makeCross2(females = pop1,
                             males = pop,
                             crossPlan = as.matrix(data.frame(cross = pop1@iid[1],
                                                              recurrent = pop0@father)),
                             nProgeny = n_ind_macs, 
                             simParam = SP)
          pop1 <- setPheno(pop1, 
                           h2 = h2,
                           H2 = H2,
                           simParam = SP)
        }
      }
      ped_dt_internal <- rbind(ped_dt_internal,
                               data.table(sample_id = pop1@id,
                                          sire = pop1@father,
                                          dam = pop1@mother,
                                          generation = rep(paste0("BC", 2), length(pop1@id)),
                                          family = rep(1, length(pop1@id))))
    }
    mergepop <- mergePops(list(pop, pop1))
  } else if (toupper(pop_type) == "MAGIC") {
    pop0 <- randCross(pop, 
                      nCrosses= round(magic_n_founders / 2), 
                      nProgeny = 1, 
                      simParam = SP)
    pop0 <- setPheno(pop0, 
                     h2 = h2,
                     H2 = H2,
                     simParam = SP)
    pop00 <- makeCross(pop0, 
                       crossPlan = as.matrix(data.frame(p1 = pop0@id[1:(magic_n_founders / 4)],
                                                        p2 = pop0@id[(magic_n_founders / 4 + 1):length(pop0@id)])),
                       nProgeny = 1, 
                       simParam = SP)
    pop00 <- setPheno(pop00, 
                      h2 = h2,
                      H2 = H2,
                      simParam = SP)
    ped_dt_internal <- data.table(sample_id = c(pop0@mother, pop0@father, pop0@id, pop00@id),
                                  sire = c(rep(NA, length(pop0@id) * 2), pop0@father, pop00@father),
                                  dam = c(rep(NA, length(pop0@id) * 2), pop0@mother, pop00@mother),
                                  generation = c(rep("Founder", length(pop0@id) * 2), 
                                                 rep("F1", length(pop0@id)),
                                                 rep("M2", length(pop00@id))),
                                  family = c(rep(0, length(pop0@id) * 2), rep(1, length(pop0@id) + length(pop00@id))))
    if (magic_n_founders > 4) {
      for (cyc1 in 1:(magic_n_founders / 8)) {
        pop00 <- makeCross(pop00, 
                           crossPlan = as.matrix(data.frame(p1 = pop00@id[1:(length(pop00@id) / 2)],
                                                            p2 = pop00@id[(length(pop00@id) / 2 + 1):length(pop00@id)])),
                           nProgeny = 1, 
                           simParam = SP)
        pop00 <- setPheno(pop00, 
                          h2 = h2,
                          H2 = H2,
                          simParam = SP)
        ped_dt_internal <- rbind(ped_dt_internal,
                                 data.table(sample_id = pop00@id,
                                            sire = pop00@father,
                                            dam = pop00@mother,
                                            generation = rep(paste0("M", cyc1 + 2), length(pop00@id)),
                                            family = rep(1, length(pop00@id))))
      }
    }
    pop1 <- self(pop00,
                 nProgeny = n_ind_macs,
                 keepParents = F, 
                 simParam = SP)
    pop1 <- setPheno(pop1, 
                     h2 = h2,
                     H2 = H2,
                     simParam = SP)
    ped_dt_internal <- rbind(ped_dt_internal,
                             data.table(sample_id = pop1@id,
                                        sire = pop1@father,
                                        dam = pop1@mother,
                                        generation = rep(paste0("S", 1), length(pop1@id)),
                                        family = rep(1, length(pop1@id))))
    for (cyc1 in 2:6) {
      pop1 <- self(pop1,
                   nProgeny = 1,
                   keepParents = F, 
                   simParam = SP)
      pop1 <- setPheno(pop1, 
                       h2 = h2,
                       H2 = H2,
                       simParam = SP)
      ped_dt_internal <- rbind(ped_dt_internal,
                               data.table(sample_id = pop1@id,
                                          sire = pop1@father,
                                          dam = pop1@mother,
                                          generation = rep(paste0("S", cyc1), length(pop1@id)),
                                          family = rep(1, length(pop1@id))))
    }
  } else if (toupper(pop_type) == "NAM") {
    pop0 <- makeCross(pop,
                      crossPlan = as.matrix(data.frame(p1 = pop@id[1],
                                                       p2 = pop@id[2:(nam_n_families + 1)])),
                      nProgeny = 1, 
                      simParam = SP)
    pop0 <- setPheno(pop0, 
                     h2 = h2,
                     H2 = H2,
                     simParam = SP)
    pop1 <- self(pop0,
                 nProgeny = round(n_ind_macs / nam_n_families),
                 keepParents = F, 
                 simParam = SP)
    pop1 <- setPheno(pop1, 
                     h2 = h2,
                     H2 = H2,
                     simParam = SP)
    ped_dt_internal <- data.table(sample_id = c(unique(pop0@mother), pop0@father, pop0@id, pop1@id),
                                  sire = c(rep(NA, nam_n_families + 1), pop0@father, pop1@father),
                                  dam = c(rep(NA, nam_n_families + 1), pop0@mother, pop1@mother),
                                  generation = c(rep("Founder", nam_n_families + 1), 
                                                 rep("F1", length(pop0@id)),
                                                 rep("F2", length(pop1@id))),
                                  family = c(rep(0, nam_n_families + 1), 
                                             1:nam_n_families,
                                             rep(1:nam_n_families, 
                                                 each = round(n_ind_macs / nam_n_families))))
    for (cyc1 in 3:nam_ril_selfing) {
      pop1 <- self(pop1,
                   nProgeny = 1,
                   keepParents = F, 
                   simParam = SP)
      pop1 <- setPheno(pop1, 
                       h2 = h2,
                       H2 = H2,
                       simParam = SP)
      ped_dt_internal <- rbind(ped_dt_internal,
                               data.table(sample_id = pop1@id,
                                          sire = pop1@father,
                                          dam = pop1@mother,
                                          generation = rep(paste0("F", cyc1), length(pop1@id)),
                                          family = rep(1:nam_n_families, 
                                                       each = round(n_ind_macs / nam_n_families))))
    }
  }
  G <- pullSnpGeno(pop1)
  geno_mat <- G
  pheno_mat <- pop1@pheno
  qtl_mat <- pullQtlGeno(pop1)
  founders_ids <- unique(c(pop1@mother,
                           pop1@father))
  if (!(is.null(mergepop))) {
    mergeG <- pullSnpGeno(mergepop)
    geno_merge <- mergeG
    if (!(toupper(pop_type) %in% c("S1_INBRED", "S1_OUTBRED"))) {
      geno_merge <- geno_merge[rownames(geno_merge) %in% c(pop0@mother, pop0@father, rownames(G)) , ]
    }
    pheno_merge <- mergepop@pheno
    qtl_merge <- pullQtlGeno(mergepop)
  }
}


# Remove internal helper column
if ("marker_index" %in% names(map_dt)) map_dt[, marker_index := NULL]

#stopifnot(!is.null(geno_mat), nrow(geno_mat) == n_ind, ncol(geno_mat) == n_markers)

geno_mask <- geno_mat
# Apply genotype error (flip among 0/1/2 randomly)
if (geno_error > 0) {
  n <- length(geno_mask)
  idx <- which(runif(n) < geno_error / 100)
  if (length(idx) > 0) {
    geno_mask[idx] <- sample(c(0,1,2), length(idx), replace=TRUE)
  }
}

# Apply missingness
if (missing_rate > 0) {
  n <- length(geno_mask)
  idx <- which(runif(n) < missing_rate / 100)
  if (length(idx) > 0) geno_mask[idx] <- NA
}

if (!exists("sample_id")) sample_id <- sprintf("id%04d", seq_len(n_ind))
colnames(geno_mat) <- map_dt$marker_id
rownames(geno_mat) <- sample_id
if (!(is.null(geno_merge))) {
  colnames(geno_merge) <- map_dt$marker_id
}

# -------------------------------
# Trait simulation + truth
# -------------------------------
truth_qtl <- c()
for (cyc1 in 1:length(qtl_info)) {
  qtl_info_each <- qtl_info[[cyc1]]
  x <- 1
  qtl_pos <- c()
  for (cyc2 in 1:length(qtl_info_each@lociPerChr)) {
    loc_each <- qtl_info_each@lociLoc[x:(x + qtl_info_each@lociPerChr[cyc2] - 1)]
    geno_map_each <- geno_map[[cyc2]][loc_each]
    qtl_pos <- c(qtl_pos,
                 geno_map_each)
    
    x <- x + qtl_info_each@lociPerChr[cyc2]
  }
  truth_qtl <- rbind(truth_qtl,
                     data.frame(trait = trait_name[cyc1],
                                qtl_id = paste0(trait_name[cyc1],
                                                "_q",
                                                1:length(qtl_pos)),
                                marker_id = names(qtl_pos),
                                chr = gsub("_.*$",
                                           "",
                                           names(qtl_pos)),
                                pos_bp = round(qtl_pos * 10000000),
                                pos_cM = qtl_pos,
                                effect.additive = qtl_info_each@addEff,
                                effect.dominance = qtl_info_each@domEff))
}

pheno_dt <- data.frame(sample_id = sample_id,
                       pheno_mat)
pheno_mask <- pheno_mat
if (pheno_missing > 0) {
  n <- length(pheno_mask)
  idx <- which(runif(n) < pheno_missing / 100)
  if (length(idx) > 0) {
    pheno_mask[idx] <- NA
  }
}
pheno_mask_dt <- data.frame(sample_id = sample_id,
                            pheno_mask)
# -------------------------------
# Split table for quick CV/testing
# -------------------------------
split <- rep("train", n_ind)
u <- runif(n_ind)
split[u >= 0.70 & u < 0.85] <- "val"
split[u >= 0.85] <- "test"
split_dt <- data.table(sample_id=sample_id, split=split)

# -------------------------------
# Write dataset structure
# -------------------------------
dataset_dir <- file.path(out_dir, "dataset")
dir.create(dataset_dir, recursive=TRUE, showWarnings=FALSE)

dirs <- c(
  "genome", "genotype", "phenotype", "pedigree", "splits",
  "qtl2", "rqtl",
  "rnaseq", "polyploid", "logs"
)
for (d in dirs) dir.create(file.path(dataset_dir, d), recursive=TRUE, showWarnings=FALSE)

fwrite(map_dt, file.path(dataset_dir, "genome", "map.tsv"), sep="\t", na = "NA", quote = F)
if (nrow(truth_qtl) > 0) fwrite(truth_qtl, file.path(dataset_dir, "genome", "truth_qtl.tsv"), sep="\t", na = "NA", quote = F)

# Founders list (useful for pedigree inspection)
if (exists("founders_ids") && length(founders_ids) > 0) {
  fwrite(data.table(founder_id=founders_ids), file.path(dataset_dir, "genome", "founders.tsv"), sep="\t", na = "NA", quote = F)
}

geno_dt <- as.data.table(geno_mat)
geno_dt <- cbind(data.table(sample_id=sample_id), geno_dt)
fwrite(geno_dt, file.path(dataset_dir, "genotype", "geno.tsv"), sep="\t", na="NA", quote = F)

geno_mask_dt <- cbind(data.table(sample_id=sample_id), as.data.frame(geno_mask))
fwrite(geno_mask_dt, file.path(dataset_dir, "genotype", "geno.missing.tsv"), sep="\t", na="NA", quote = F)

Amat <- rrBLUP::A.mat(geno_mat - 1,
                      min.MAF = NULL,
                      max.missing = NULL,
                      impute.method = "mean",
                      tol = 0.02,
                      n.core = 1,
                      shrink = F,
                      return.imputed = F)
fwrite(as.data.frame(Amat), 
       file.path(dataset_dir, "genotype", "Amat.tsv"), sep="\t", na="NA", row.names = T, col.names = T, quote = F)

if (!(is.null(geno_merge))) {
  geno_merge_dt <- as.data.table(geno_merge)
  merge_id <- rownames(geno_merge)
  merge_id[(length(merge_id) - n_ind + 1):length(merge_id)] <- sample_id
  geno_merge_dt <- cbind(data.table(sample_id=merge_id), geno_merge_dt)
  fwrite(geno_merge_dt, file.path(dataset_dir, "genotype", "geno_merge.tsv"), sep="\t", na="NA", quote = F)
}

#PLINK and VCF outputs
#bed <- gaston::as.bed.matrix(x = geno_mat + 1)
bed <- gaston::as.bed.matrix(x = geno_mat)
bed@snps$chr <- chr
bed@snps$id <- marker_id
bed@snps$dist <- pos_cM
bed@snps$pos <- pos_bp
gaston::write.bed.matrix(x = bed,
                         basename = file.path(dataset_dir, "genotype", "geno"))

fwrite(pheno_dt, file.path(dataset_dir, "phenotype", "pheno.tsv"), sep="\t", na = "NA", quote = F)
fwrite(pheno_mask_dt, file.path(dataset_dir, "phenotype", "pheno.missing.tsv"), sep="\t", na = "NA", quote = F)
fwrite(split_dt, file.path(dataset_dir, "splits", "split.tsv"), sep="\t", na = "NA", quote = F)

# Pedigree
if (exists("ped_dt_internal")) {
  ped_dt <- ped_dt_internal
  # ensure sample rows exist
  miss <- setdiff(sample_id, ped_dt$sample_id)
  if (length(miss) > 0) {
    ped_dt <- rbind(ped_dt, data.table(sample_id=miss, sire=NA_character_, dam=NA_character_, generation=pop_type, family=1), fill=TRUE)
  }
} else {
  ped_dt <- data.table(sample_id=sample_id, sire=NA_character_, dam=NA_character_, generation=pop_type, family=1)
}
fwrite(ped_dt, file.path(dataset_dir, "pedigree", "pedigree.tsv"), sep="\t", na="NA", quote = F)

# -------------------------------
# Export formats: qtl2 / r/qtl
# -------------------------------

# Use family covariate (useful for MAGIC/NAM as a fixed effect)
# NOTE:
# Avoid data.table's `..var` scoping for maximum compatibility across versions.
# (Some older data.table builds used in conda/R environments can error with `..sample_id`.)
sample_ids <- sample_id
covar_dt <- unique(ped_dt[sample_id %in% sample_ids, .(id=sample_id, family=as.integer(family))])

if (isTRUE(export_qtl2)) {
  qtl2_dir <- file.path(dataset_dir, "qtl2")
  dir.create(qtl2_dir, recursive=TRUE, showWarnings=FALSE)

  # qtl2 expects id column named "id"
  geno_qtl2 <- as.data.table(geno_mat + 1)
  geno_qtl2 <- cbind(data.table(id=sample_id), geno_qtl2)
  fwrite(geno_qtl2, file.path(qtl2_dir, "geno.csv"), na = "NA", quote = F)

  ph_qtl2 <- copy(pheno_dt)
  setnames(ph_qtl2, "sample_id", "id")
  fwrite(ph_qtl2, file.path(qtl2_dir, "pheno.csv"), na = "NA", quote = F)

  if (nrow(covar_dt) > 0) fwrite(covar_dt, file.path(qtl2_dir, "covar.csv"), quote = F)

  pmap <- map_dt[, .(marker=marker_id, chr=chr, pos=pos_bp)]
  gmap <- map_dt[, .(marker=marker_id, chr=chr, pos=pos_cM)]
  fwrite(pmap, file.path(qtl2_dir, "pmap.csv"), na = "NA", quote = F)
  fwrite(gmap, file.path(qtl2_dir, "gmap.csv"), na = "NA", quote = F)

  # Determine crosstype for qtl2
  pop_u <- toupper(pop_type)
  crosstype <- if (pop_u == "F2") "f2" else if (pop_u == "BC") "bc" else "riself"
  if (pop_u %in% c("MAGIC", "NAM")) {
    cat("[simulate_dataset] NOTE: qtl2 export for ", pop_u, " uses crosstype=riself with covar family.\n")
  }

  yaml_lines <- c(
    paste0("crosstype: ", crosstype),
    "geno:",
    "  file: geno.csv",
    "pheno:",
    "  file: pheno.csv",
    "gmap:",
    "  file: gmap.csv",
    "pmap:",
    "  file: pmap.csv",
    "genotypes:",
    "  1: 1",
    "  2: 2",
    "  3: 3",
    "  4: 4",
    "  5: 5"
  )
  if (file.exists(file.path(qtl2_dir, "covar.csv"))) {
    yaml_lines <- c(yaml_lines, "covar:", "  file: covar.csv")
  }
  writeLines(yaml_lines, file.path(qtl2_dir, "cross2.yaml"))
}
asf <- function(x) {
  as.numeric(factor(x))
}
if (isTRUE(export_rqtl)) {
  rqtl_dir <- file.path(dataset_dir, "rqtl")
  dir.create(rqtl_dir, recursive=TRUE, showWarnings=FALSE)
  
  p1 <- as.vector(t(geno_merge[1, ]))
  p2 <- as.vector(t(geno_merge[2, ]))
  pp <- geno_merge[3:nrow(geno_merge), ]
  p1_2 <- p1 == 2
  pp_2 <- pp[, p1_2]
  # write.csv(pp_2,
  #           "/media/soba/Noc4/GenomicExplorer/P0213/sim/ex3.csv",
            # quote = F)
  pp_2 <- (pp_2 - 1) * (-1) + 1
  pp[, p1_2] <- pp_2
  # write.csv(pp,
  #           "/media/soba/Noc4/GenomicExplorer/P0213/sim/ex4.csv",
  #           quote = F)
  #geno_mat <- pp

  if (!requireNamespace("qtl", quietly=TRUE)) {
    cat("[simulate_dataset] WARN: package 'qtl' not installed; skipping r/qtl cross export.\n")
  } else {
    suppressPackageStartupMessages(library(qtl))

    # Convert 0/1/2 dosage to R/qtl codes: 1=AA, 2=AB, 3=BB
    #Xcode <- geno_mat
    Xcode <- pp
    Xcode2 <- matrix(NA, nrow=nrow(Xcode), ncol=ncol(Xcode))
    # Xcode2[!is.na(Xcode) & Xcode==0] <- 1L
    # Xcode2[!is.na(Xcode) & Xcode==1] <- 2L
    # Xcode2[!is.na(Xcode) & Xcode==2] <- 3L
    Xcode2[!is.na(Xcode) & Xcode==0] <- "AA"
    Xcode2[!is.na(Xcode) & Xcode==1] <- "AB"
    Xcode2[!is.na(Xcode) & Xcode==2] <- "BB"
    #Xcode2[!is.na(Xcode) & Xcode==0] <- "1"
    #Xcode2[!is.na(Xcode) & Xcode==1] <- "2"
    #Xcode2[!is.na(Xcode) & Xcode==2] <- "3"
    rownames(Xcode2) <- sample_id
    colnames(Xcode2) <- map_dt$marker_id

    ct <- tolower(if (toupper(pop_type) == "F2") "f2" else if (toupper(pop_type) == "BC") "bc" else "riself")
    if (toupper(pop_type) %in% c("MAGIC", "NAM")) {
      cat("[simulate_dataset] NOTE: r/qtl export for ", toupper(pop_type), " uses cross type 'riself' and family as covariate.\n")
    }

    # Split markers by chr and build geno list
    idx_by_chr <- split(seq_len(nrow(map_dt)), map_dt$chr)
    geno_list <- list()
    if (toupper(pop_type) == "F2") {
      #genotypes = c("1", "2", "3")
      genotypes = c("AA", "AB", "BB")
      #genotypes = c("A", "H", "B")
    } else {
      #genotypes = c("1", "2")
      genotypes = c("AA", "BB")
      #genotypes = c("A", "B")
    }
    for (cc in names(idx_by_chr)) {
      idx <- idx_by_chr[[cc]]
      mk <- map_dt$marker_id[idx]
      o <- order(map_dt$pos_cM[idx])
      mk <- mk[o]
      submat <- Xcode2[, mk, drop=FALSE]
      geno_list[[as.character(cc)]] <- list(
        data = apply(submat, 2, asf),
        data = submat,
        map = setNames(as.numeric(map_dt$pos_cM[idx][o]), mk),
        alleles = c("A", "B"),
        genotypes = genotypes
      )
    }

    ph <- as.data.frame(pheno_dt)
    ph$family <- covar_dt$family[match(ph$sample_id, covar_dt$id)]
    rownames(ph) <- ph$sample_id

    cross <- list(pheno=ph, geno=geno_list)
    class(cross) <- c(ct, "cross")

    saveRDS(cross, file.path(rqtl_dir, "cross.rds"))
    save(cross, file=file.path(rqtl_dir, "cross.RData"))
    
    # if (ct == "f2") {
    #   cross_asmap <- cross
    #   class(cross_asmap) <- c("bcsft", "cross")
    #   saveRDS(cross_asmap, file.path(rqtl_dir, "cross.asmap.rds"))
    # }
    ### Lep-Map3 pedigree
    if (!(is.null(geno_merge))) {
      lepped <- cbind(rep("CHR", 6),
                      rep("POS", 6),
                      rbind(rep("F", nrow(geno_merge)),
                            rownames(geno_merge),
                            c(0, 0, rep(rownames(geno_merge)[1], nrow(geno_merge) - 2)),
                            c(0, 0, rep(rownames(geno_merge)[2], nrow(geno_merge) - 2)),
                            c(1, 2, rep(0, nrow(geno_merge) - 2)),
                            rep(0, nrow(geno_merge))))
      write.table(lepped,
                  file.path(dataset_dir, "genotype", "lep-map.ped.tsv"),
                  quote = F,
                  sep = "\t",
                  row.names = F,
                  col.names = F)
    }
  }
}

# Polyploid outputs (Updog depth + truth dosage)
dosage_mat <- geno_mat
if (isTRUE(generate_polyploid_depth)) {
  poly_dir <- file.path(dataset_dir, "polyploid")
  dir.create(poly_dir, recursive=TRUE, showWarnings=FALSE)

  # Write truth dosage matrix
  dose_dt <- as.data.table(dosage_mat)
  dose_dt <- cbind(data.table(sample_id=sample_id), dose_dt)
  fwrite(dose_dt, file.path(poly_dir, "truth_dosage.tsv"), sep="\t", na="NA", quote = F)
  if (!(is.null(geno_merge))) {
    dose_merge_dt <- as.data.table(geno_merge)
    merge_id <- rownames(geno_merge)
    merge_id[(length(merge_id) - n_ind + 1):length(merge_id)] <- sample_id
    dose_merge_dt <- cbind(data.table(sample_id=merge_id), dose_merge_dt)
    fwrite(dose_merge_dt, file.path(poly_dir, "truth_dosage_merge.tsv"), sep="\t", na="NA", quote = F)
    
    rownames(geno_merge) <- merge_id
    dose_merge_mk <- data.frame(marker_id = marker_id,
                                t(geno_merge[1:2, ]),
                                chr = chr,
                                pos = pos_bp,
                                t(geno_merge[3:nrow(geno_merge), ]))
    fwrite(dose_merge_mk, file.path(poly_dir, "truth_dosage_merge_marker.csv"), sep=",", na="NA", quote = F)
  }

  total <- c()
  ref <- c()
  for (cyc1 in 1:ncol(dosage_mat)) {
    sizevec <- stats::rpois(n = nrow(dosage_mat), 
                            lambda = depth_mean)
    refvec  <- rflexdog(sizevec = sizevec, 
                        geno = dosage_mat[, cyc1],
                        ploidy = ploidy, 
                        seq = seq,
                        bias = bias, 
                        od = od)
    total <- cbind(total,
                   sizevec)
    ref <- cbind(ref,
                 refvec)
  }
  rownames(ref) <- rownames(total) <- rownames(dosage_mat)
  colnames(ref) <- colnames(total) <- colnames(dosage_mat)
  alt <- total - ref
  
  ref1 <- data.frame(sample_id = rownames(dosage_mat), ref)
  total1 <- data.frame(sample_id = rownames(dosage_mat), total)
  
  ref_depth <- file.path(poly_dir, "ref_depth.tsv")
  tot_depth <- file.path(poly_dir, "total_depth.tsv")
  fwrite(ref1, ref_depth, sep="\t", row.names = F, na = "NA", quote = F)
  fwrite(total1, tot_depth, sep="\t", row.names = F, na = "NA", quote = F)
  
  if (!(is.null(geno_merge))) {
    dosage_mat2 <- geno_merge
    total2 <- c()
    ref2 <- c()
    for (cyc1 in 1:ncol(dosage_mat2)) {
      sizevec <- stats::rpois(n = nrow(dosage_mat2), 
                              lambda = depth_mean)
      refvec  <- rflexdog(sizevec = sizevec, 
                          geno = dosage_mat2[, cyc1],
                          ploidy = ploidy, 
                          seq = seq,
                          bias = bias, 
                          od = od)
      total2 <- cbind(total2,
                     sizevec)
      ref2 <- cbind(ref2,
                   refvec)
    }
    colnames(ref2) <- colnames(total2) <- colnames(dosage_mat2)
    ref3 <- data.frame(sample_id = rownames(dosage_mat2), ref2)
    total3 <- data.frame(sample_id = rownames(dosage_mat2), total2)
    
    ref_merge_depth <- file.path(poly_dir, "ref_merge_depth.tsv")
    tot_merge_depth <- file.path(poly_dir, "total_merge_depth.tsv")
    fwrite(ref3, ref_merge_depth, sep="\t", row.names = F, na = "NA", quote = F)
    fwrite(total3, tot_merge_depth, sep="\t", row.names = F, na = "NA", quote = F)
  }

  # Also export a qtl2-side helper (NOT directly readable by qtl2)
  qtl2_poly_dir <- file.path(dataset_dir, "qtl2_poly")
  dir.create(qtl2_poly_dir, recursive=TRUE, showWarnings=FALSE)
  dose_csv <- file.path(qtl2_poly_dir, "dosage.csv")
  dose_csv_dt <- as.data.table(dosage_mat)
  dose_csv_dt <- cbind(data.table(id=sample_id), dose_csv_dt)
  fwrite(dose_csv_dt, dose_csv, na = "NA", quote = F)
  writeLines(c(
    "This folder contains polyploid helper files.",
    "- dosage.csv: integer dosage (0..P) per marker.",
    "- updog_depth.tsv (in polyploid/): ref/alt depths suitable for Updog.",
    "To use with qtl2-style workflows, estimate genotype probabilities per marker",
    "(e.g., with Updog) and convert them into a genoprobs object before running scans."
  ), file.path(qtl2_poly_dir, "README.txt"))
  
  # タビックス索引も作りたい場合（任意）
  #library(Rsamtools)
  #indexTabix(file.path(dataset_dir, "genotype", "geno.vcf.gz"), format = "vcf")
  
  message("Wrote: ", file.path(dataset_dir, "genotype", "geno.vcf"))
  ################################################################################
  
  dose_csv2 <- file.path(qtl2_poly_dir, "dosage-map.tsv")
  dose_csv_dt2 <- as.data.table(data.frame(map_dt[, 1:3], REF=0, ALT=1, t(dosage_mat)))
  dose_csv_dt2 <- cbind(data.table(dose_csv_dt2))
  fwrite(dose_csv_dt2, dose_csv2, sep = "\t", na = "NA", quote = F)
}

################################################################################
### VCF ###
rr <- GRanges(
  seqnames = chr,
  ranges   = IRanges(start = pos_bp, width = rep(1, length(pos_bp))),
  strand   = "*"
)
names(rr) <- marker_id

ref_dna <- DNAStringSet(rep("A", length(chr)))
alt_dna <- DNAStringSetList("T")

mcols(rr)$REF <- ref_dna
mcols(rr)$ALT <- alt_dna

samples <- rownames(geno_mat)
variant_ids <- marker_id
nvar <- length(variant_ids)
ns  <- length(samples)

GT_mat <- geno_mat
GT_mat[geno_mat == 0] <- "0/0"
GT_mat[geno_mat == 1] <- "0/1"
GT_mat[geno_mat == 2] <- "1/1"
GT_mat[is.na(geno_mat)] <- "./."

# ---------------------------
# 6) VCF Header を作る（FORMAT/INFO 定義）
# ---------------------------
hdr <- VCFHeader()
meta(hdr) <- DataFrameList(fileformat = DataFrame(Value = "VCFv4.3"))
g <- geno(hdr)
g <- rbind(g, DataFrame(Number="1", Type="String", Description="Genotype",
                        row.names="GT"))

DP_mat <- AD_arr <- NULL
if (isTRUE(generate_polyploid_depth)) {
  DP_mat <- total
  storage.mode(DP_mat) <- "integer"
  
  # AD（ref, alt の2値を想定。複数ALTの場合は alt 列数を増やす必要あり）
  AD_arr <- NULL
  AD_ref <- ref
  AD_alt <- alt
  storage.mode(AD_ref) <- "integer"
  storage.mode(AD_alt) <- "integer"
  
  # 3次元配列: variants x samples x alleles(REF,ALT)
  AD_arr <- array(NA_integer_, dim = c(nvar, ns, 2),
                  dimnames = list(variant_ids, samples, c("REF","ALT")))
  AD_arr[, , 1] <- AD_ref
  AD_arr[, , 2] <- AD_alt
  g <- rbind(g, DataFrame(Number="1", Type="Integer", Description="Read depth",
                          row.names="DP"))
  g <- rbind(g, DataFrame(Number="R", Type="Integer",
                          Description="Allelic depths for ref and alt alleles",
                          row.names="AD"))
}
geno(hdr) <- g

# ---------------------------
# 7) VCF
# ---------------------------
# fixed: CHROM, POS, ID, REF, ALT, QUAL, FILTER
fixed <- DataFrame(
  REF    = mcols(rr)$REF,
  ALT    = mcols(rr)$ALT,
  QUAL   = if ("QUAL" %in% names(mcols(rr))) mcols(rr)$QUAL else rep(NA_real_, nvar),
  FILTER = if ("FILTER" %in% names(mcols(rr))) mcols(rr)$FILTER else rep("PASS", nvar)
)

# info: rr の mcols に INFO が入っていればそのまま DataFrame 化
info_keys <- setdiff(names(mcols(rr)), c("REF","ALT","QUAL","FILTER"))
info <- if (length(info_keys) > 0) DataFrame(mcols(rr)[, info_keys, drop=FALSE]) else DataFrame()

# geno: GT/DP/AD
geno <- SimpleList(GT = t(GT_mat))
if (!is.null(DP_mat)) geno$DP <- t(DP_mat)
if (!is.null(AD_arr)) geno$AD <- t(matrix(paste(ref, alt, sep = ","), nrow = nrow(ref), ncol = ncol(ref)))
vcf <- VCF(
  rowRanges = rr,
  colData   = DataFrame(row.names = rownames(geno_mat)),
  fixed     = fixed,
  geno      = geno,
)
header(vcf) <- hdr
writeVcf(vcf, filename = file.path(dataset_dir, "genotype", "geno.vcf"))

if (!(is.null(geno_merge))) {
  ################################################################################
  ### Merge VCF ###
  rr <- GRanges(
    seqnames = chr,
    ranges   = IRanges(start = pos_bp, width = rep(1, length(pos_bp))),
    strand   = "*"
  )
  names(rr) <- marker_id
  
  ref_dna <- DNAStringSet(rep("A", length(chr)))
  alt_dna <- DNAStringSetList("T")
  
  mcols(rr)$REF <- ref_dna
  mcols(rr)$ALT <- alt_dna
  
  samples <- rownames(geno_merge)
  variant_ids <- marker_id
  nvar <- length(variant_ids)
  ns  <- length(samples)
  
  GT_mat <- geno_merge
  GT_mat[geno_merge == 0] <- "0/0"
  GT_mat[geno_merge == 1] <- "0/1"
  GT_mat[geno_merge == 2] <- "1/1"
  GT_mat[is.na(geno_merge)] <- "./."
  
  # ---------------------------
  # 6) VCF Header を作る（FORMAT/INFO 定義）
  # ---------------------------
  hdr <- VCFHeader()
  meta(hdr) <- DataFrameList(fileformat = DataFrame(Value = "VCFv4.3"))
  g <- geno(hdr)
  g <- rbind(g, DataFrame(Number="1", Type="String", Description="Genotype",
                          row.names="GT"))
  
  DP_mat <- AD_arr <- NULL
  if (isTRUE(generate_polyploid_depth)) {
    DP_mat <- total2
    storage.mode(DP_mat) <- "integer"
    
    # AD（ref, alt の2値を想定。複数ALTの場合は alt 列数を増やす必要あり）
    AD_arr <- NULL
    AD_ref <- ref2
    AD_alt <- alt2 <- total2 - ref2
    storage.mode(AD_ref) <- "integer"
    storage.mode(AD_alt) <- "integer"
    
    # 3次元配列: variants x samples x alleles(REF,ALT)
    AD_arr <- array(NA_integer_, dim = c(nvar, ns, 2),
                    dimnames = list(variant_ids, samples, c("REF","ALT")))
    AD_arr[, , 1] <- AD_ref
    AD_arr[, , 2] <- AD_alt
    g <- rbind(g, DataFrame(Number="1", Type="Integer", Description="Read depth",
                            row.names="DP"))
    g <- rbind(g, DataFrame(Number="R", Type="Integer",
                            Description="Allelic depths for ref and alt alleles",
                            row.names="AD"))
  }
  geno(hdr) <- g
  
  # ---------------------------
  # 7) VCF
  # ---------------------------
  # fixed: CHROM, POS, ID, REF, ALT, QUAL, FILTER
  fixed <- DataFrame(
    REF    = mcols(rr)$REF,
    ALT    = mcols(rr)$ALT,
    QUAL   = if ("QUAL" %in% names(mcols(rr))) mcols(rr)$QUAL else rep(NA_real_, nvar),
    FILTER = if ("FILTER" %in% names(mcols(rr))) mcols(rr)$FILTER else rep("PASS", nvar)
  )
  
  # info: rr の mcols に INFO が入っていればそのまま DataFrame 化
  info_keys <- setdiff(names(mcols(rr)), c("REF","ALT","QUAL","FILTER"))
  info <- if (length(info_keys) > 0) DataFrame(mcols(rr)[, info_keys, drop=FALSE]) else DataFrame()
  
  # geno: GT/DP/AD
  geno <- SimpleList(GT = t(GT_mat))
  if (!is.null(DP_mat)) geno$DP <- t(DP_mat)
  if (!is.null(AD_arr)) geno$AD <- t(matrix(paste(ref2, alt2, sep = ","), nrow = nrow(ref2), ncol = ncol(ref2)))
  vcf <- VCF(
    rowRanges = rr,
    colData   = DataFrame(row.names = rownames(geno_merge)),
    fixed     = fixed,
    geno      = geno,
  )
  header(vcf) <- hdr
  writeVcf(vcf, filename = file.path(dataset_dir, "genotype", "geno.merge.vcf"))
  
  ### Lep-Map3 pedigree
  if (!(is.null(geno_merge))) {
    lepped <- cbind(rep("CHR", 6),
                    rep("POS", 6),
                    rbind(rep("F", nrow(geno_merge)),
                          merge_id,
                          c(0, 0, rep(rownames(geno_merge)[1], nrow(geno_merge) - 2)),
                          c(0, 0, rep(rownames(geno_merge)[2], nrow(geno_merge) - 2)),
                          c(1, 2, rep(0, nrow(geno_merge) - 2)),
                          rep(0, nrow(geno_merge))))
    write.table(lepped,
                file.path(dataset_dir, "genotype", "lep-map.ped.tsv"),
                quote = F,
                sep = "\t",
                row.names = F,
                col.names = F)
  }
}

manifest <- list(
  generator = list(id="simulate_dataset", version="0.4.0"),
  seed = seed,
  species_template = species,
  genome = list(n_chr=n_chr, chr_names=chr_names, chr_len_bp=chr_len_bp),
  population = list(
    type=pop_type,
    n_individuals=n_ind,
    n_markers=n_markers,
    options=list(
      ril_selfing=ril_selfing,
      bc_generation=bc_generation,
      bc_recurrent=bc_recurrent,
      magic_n_founders=magic_n_founders,
      magic_n_families=magic_n_families,
      nam_n_families=nam_n_families,
      nam_ril_selfing=nam_ril_selfing
    )
  ),
  genotype = list(missing_rate=missing_rate, genotype_error=geno_error),
  polyploid = list(
    ploidy = ploidy,
    generate_polyploid_depth = generate_polyploid_depth,
    depth_mean = depth_mean,
    depth_nb_size = depth_nb_size,
    allele_bias_logodds = allele_bias_logodds,
    seq_error = seq_error,
    files = list(
      truth_dosage_tsv = if (generate_polyploid_depth && file.exists(file.path(dataset_dir, "polyploid", "truth_dosage.tsv"))) "polyploid/truth_dosage.tsv" else NULL,
      updog_depth_tsv = if (generate_polyploid_depth && file.exists(file.path(dataset_dir, "polyploid", "updog_depth.tsv"))) "polyploid/updog_depth.tsv" else NULL,
      qtl2_poly_dosage_csv = if (generate_polyploid_depth && file.exists(file.path(dataset_dir, "qtl2_poly", "dosage.csv"))) "qtl2_poly/dosage.csv" else NULL
    )
  ),
  trait = list(name=trait_name, 
               #architecture=trait_arch, 
               h2=h2,
               n_qtl_per_chr=n_qtl_per_chr, 
               #qtl_effect_sd=qtl_effect_sd,
               truth_qtl=if (nrow(truth_qtl) > 0) "genome/truth_qtl.tsv" else NULL),
  files = list(
    map_tsv="genome/map.tsv",
    founders_tsv=if (exists("founders_ids") && length(founders_ids) > 0) "genome/founders.tsv" else NULL,
    genotype_tsv="genotype/geno.tsv",
    phenotype_tsv="phenotype/pheno.tsv",
    split_tsv="splits/split.tsv",
    pedigree_tsv="pedigree/pedigree.tsv",
    qtl2_yaml=if (export_qtl2 && file.exists(file.path(dataset_dir, "qtl2", "cross2.yaml"))) "qtl2/cross2.yaml" else NULL,
    rqtl_cross_rds=if (export_rqtl && file.exists(file.path(dataset_dir, "rqtl", "cross.rds"))) "rqtl/cross.rds" else NULL
  )
)
writeLines(toJSON(manifest, auto_unbox=TRUE, pretty=TRUE), file.path(dataset_dir, "manifest.json"))

cat("[simulate_dataset] wrote dataset_dir=", dataset_dir, "\n")

# Optional: copy/export to user-selected directory
export_dir <- p$export_dir
if (!is.null(export_dir) && nchar(export_dir) > 0) {
  export_dir <- normalizePath(export_dir, winslash="/", mustWork=FALSE)
  dir.create(export_dir, recursive=TRUE, showWarnings=FALSE)
  target <- file.path(export_dir, basename(dataset_dir))
  cat("[simulate_dataset] exporting dataset to ", target, "\n")
  # copy recursively
  # If target exists, create unique suffix
  if (dir.exists(target)) {
    target <- file.path(export_dir, paste0(basename(dataset_dir), "_", format(Sys.time(), "%Y%m%d_%H%M%S")))
  }
  ok <- file.copy(dataset_dir, target, recursive=TRUE)
  cat("[simulate_dataset] export_ok=", ok, "\n")
}

cat("[simulate_dataset] done\n")
sink()
