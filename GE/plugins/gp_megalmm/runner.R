#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(jsonlite)
  library(data.table)
  library(MegaLMM)
  library(rrBLUP)
})

# ------------ Helpers ------------
parse_trait_spec <- function(spec, trait_cols) {
  if (is.null(spec)) spec <- ""
  spec <- trimws(as.character(spec))
  if (nchar(spec) == 0) return(seq_along(trait_cols))
  spec2 <- gsub("\\s+", "", spec)
  toks <- unlist(strsplit(spec2, ","))
  idxs <- integer(0)
  for (tk in toks) {
    if (!nchar(tk)) next
    if (grepl("^[0-9]+-[0-9]+$", tk)) {
      ab <- as.integer(unlist(strsplit(tk, "-")))
      a <- ab[1]; b <- ab[2]
      if (is.na(a) || is.na(b)) stop(paste0("Invalid trait range: ", tk))
      if (a > b) { tmp <- a; a <- b; b <- tmp }
      idxs <- c(idxs, seq.int(a, b))
    } else if (grepl("^[0-9]+$", tk)) {
      idxs <- c(idxs, as.integer(tk))
    } else {
      j <- match(tk, trait_cols)
      if (is.na(j)) stop(paste0("Trait name not found: ", tk))
      idxs <- c(idxs, j)
    }
  }
  idxs <- idxs[!is.na(idxs)]
  if (!length(idxs)) stop("No valid traits parsed from trait(s) spec.")
  idxs <- idxs[!duplicated(idxs)]
  if (any(idxs < 1 | idxs > length(trait_cols))) stop(paste0("Trait index out of range: 1..", length(trait_cols)))
  idxs
}

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) return(args[i + 1])
  default
}

params_path <- get_arg("--params")
out_dir <- get_arg("--out")
if (is.null(params_path) || is.null(out_dir)) stop("Usage: --params params.json --out out_dir")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(out_dir, "run.log")
sink(log_file, split = TRUE)

cat("[gp_megalmm] start\n")
cat("[gp_megalmm] params_path=", params_path, "\n")
cat("[gp_megalmm] out_dir=", out_dir, "\n")

p <- fromJSON(params_path)

geno_tsv <- p$genotype_tsv
pheno_tsv <- p$phenotype_tsv

traits <- ""
if (!is.null(p$traits) && nchar(as.character(p$traits)) > 0) {
  traits <- as.character(p$traits)
}

mode <- p$mode
k_folds <- as.integer(p$k_folds)
seed <- as.integer(p$seed)

center <- as.logical(p$center_markers)
scale <- as.logical(p$scale_markers)
impute <- p$impute_missing

K_factors <- as.integer(p$n_factors)
nIter <- as.integer(p$nIter)
burnIn <- as.integer(p$burnIn)
thin <- as.integer(p$thin)

extra <- p$extra_options

cat("[gp_megalmm] genotype_tsv=", geno_tsv, "\n")
cat("[gp_megalmm] phenotype_tsv=", pheno_tsv, "\n")
cat("[gp_megalmm] mode=", mode, " k_folds=", k_folds, " seed=", seed, "\n")
cat("[gp_megalmm] K_factors=", K_factors, " nIter=", nIter, " burnIn=", burnIn, " thin=", thin, "\n")

if (mode == "kfold" && (is.na(k_folds) || k_folds < 2)) stop("k_folds must be >=2")

# Load tables
geno <- data.frame(fread(geno_tsv,
                         header = T),
                   row.names = 1)
if (ncol(geno) < 2) {
  stop("genotype_tsv must have at least 2 columns: <ID> + markers")
} 
ids_geno <- rownames(geno)


pheno <- data.frame(fread(pheno_tsv, 
                          header = T),
                    row.names = 1)
if (ncol(pheno) < 2) {
  stop("phenotype_tsv must have at least 2 columns: <ID> + trait(s)")
}
ids_pheno <- rownames(pheno)
trait_names <- colnames(pheno)

trait_idxs <- parse_trait_spec(traits, trait_names)
trait_target <- trait_names[trait_idxs]

ids_match <- intersect(ids_geno,
                       ids_pheno)
if (length(ids_match) < 1) {
  stop("No overlapping IDs between genotype and phenotype")
}

geno <- geno[ids_match, , drop = F]
pheno <- pheno[ids_match, trait_target, drop = F]

cat("[gp_megalmm] n_samples=", length(ids_match), " n_markers=", ncol(geno), " n_traits=", ncol(pheno), "\n")

for (j in seq_len(ncol(geno))) {
  v <- geno[, j]
  if (anyNA(v)) {
    if (impute == "mean") {
      mu <- mean(v, na.rm = TRUE)
    } else {
      mu <- median(v, na.rm = TRUE)
    }
    if (!is.finite(mu)) mu <- 0
    v[is.na(v)] <- mu
    geno[, j] <- v
  }
}

Kmat <- rrBLUP::A.mat(geno)

# ---------- helpers ----------
rmse <- function(a, b) {
  ok <- is.finite(a) & is.finite(b)
  if (!any(ok)) return(NA_real_)
  sqrt(mean((a[ok] - b[ok])^2))
}
mae <- function(a, b) {
  ok <- is.finite(a) & is.finite(b)
  if (!any(ok)) return(NA_real_)
  mean(abs(a[ok] - b[ok]))
}
pearson <- function(a, b) {
  ok <- is.finite(a) & is.finite(b)
  if (sum(ok) < 3) return(NA_real_)
  suppressWarnings(cor(a[ok], b[ok], method = "pearson"))
}

set.seed(seed)
cat("[gp_megalmm] attempting MegaLMM multi-trait fit ...\n")
# Data preparation
trait_names <- trait_target
if (!(is.null(p$covariate_tsv))) {
  cov <- data.frame(fread(cov,
                          header = T),
                    row.names = 1)
  cov1 <- cov[rownames(Y), , drop = F]
  cov2 <- as.data.frame(na.omit(cov1))
  if (nrow(cov2) == nrow(Y)) {
    cov_names <- colnames(cov2)
    tc_data <- data.frame(Y,
                          cov2)
    colnames(tc_data) <- c(trait_names,
                           cov_names)
  }
} else {
  tc_data <- pheno
}
tc_data1 <- data.frame(ID = rep(rownames(pheno), 
                                length(trait_names)),
                       Pheno = rep(trait_names, 
                               each = nrow(pheno)),
                       Value = as.vector(as.matrix(tc_data[, trait_names])))
if (!(is.null(p$covariate_tsv))) {
  for (cyc1 in 1:length(cov_names)) {
    tc_data1 <- cbind(tc_data1,
                      rep(tc_data[, cov_names[cyc1]], length(trait_names)))
  }
  colnames(tc_data1)[4:ncol(tc_data1)] <- cov_names
  id_cols <- c("ID",
               cov_names)
} else {
  id_cols <- "ID"
}
tc_data1 <- as.data.frame(tc_data1)

data_matrices = create_data_matrices(
  tall_data = tc_data1,
  id_cols = id_cols,
  names_from = 'Pheno',
  values_from = 'Value'
)

run_parameters = MegaLMM_control(
  which_sampler = list(Y = 1, 
                       F = 1),
  run_sampler_times = 1,
  #scale_Y = rep(T, F),
  K = K_factors,
  h2_divisions = 20,
  h2_step_size = NULL,
  drop0_tol = 1e-14,
  K_eigen_tol = 1e-10,
  burn = burnIn,
  thin = thin,
  max_NA_groups = Inf,
  svd_K = TRUE,
  verbose = F,
  save_current_state = TRUE,
  diagonalize_ZtZ_Kinv = TRUE
)

Lambda_prior = list(
  sampler = sample_Lambda_prec_ARD,
  # function that implements the ARD Lambda prior
  # described in Runcie et al 2013 paper.
  #See code to see requirements for this function.
  # other options are:
  # ?sample_Lambda_prec_horseshoe
  # ?sample_Lambda_prec_BayesC
  Lambda_df = 3,
  delta_1   = list(shape = 2, rate = 1),
  delta_2   = list(shape = 3, rate = 1),  # parameters of the gamma distribution giving the expected change in proportion of non-zero loadings in each consecutive factor
  # parameters of the gamma distribution giving the expected change
  # in proportion of non-zero loadings in each consecutive factor
  delta_iterations_factor = 100
  # parameter that affects mixing of the MCMC sampler. This value is generally fine.
)

priors = MegaLMM_priors(
  tot_Y_var = list(V = 0.5, nu = 5),      
  tot_F_var = list(V = 18 / 20, nu = 20),     
  h2_priors_resids_fun = function(h2s, n)  1,  
  h2_priors_factors_fun = function(h2s, n) 1, 
  Lambda_prior = Lambda_prior,
  B2_prior = list(sampler = sample_B2_prec_horseshoe, prop_0 = 0.1),
  cis_effects_prior = list(prec = 1)
)
# ---------- CV / predictions per trait ----------
summary_rows <- list()
artifact_tables <- list()

Y <- data_matrices$Y
pred <- matrix(NA, 
               nrow = nrow(Y),
               ncol = ncol(Y))
rownames(pred) <- rownames(Y)
colnames(pred) <- colnames(Y)
fold_id <- rep(NA, nrow(Y))

fold_ID_matrix = matrix(NA,
                        nrow = nrow(Y),
                        ncol = ncol(Y),
                        dimnames = dimnames(Y))

mega <- function(x) {
  MegaLMM_state = setup_model_MegaLMM(
    Y = x,  
    formula = ~ (1 | ID),  
    extra_regressions = NULL,
    data = data_matrices$data,         
    relmat = list(ID = Kmat), 
    cis_genotypes = NULL,
    Lambda_fixed = NULL,
    run_parameters=run_parameters,
    posteriorSample_params = c("Lambda", "U_F", "F", "delta", "tot_F_prec", "F_h2",
                               "tot_Eta_prec", "resid_h2", "B1", "B2_F", "B2_R", "U_R", "cis_effects",
                               "Lambda_m_eff", "Lambda_pi", "B2_R_pi", "B2_F_pi"),
    posteriorMean_params = c(),
    posteriorFunctions = list(),
    run_ID = "MegaLMM_fit"
  )
  MegaLMM_state = set_priors_MegaLMM(MegaLMM_state,
                                     priors)
  MegaLMM_state = initialize_variables_MegaLMM(MegaLMM_state)
  MegaLMM_state$run_parameters$burn = run_parameters$burn
  MegaLMM_state = initialize_MegaLMM(MegaLMM_state,
                                     verbose = F)
  MegaLMM_state$Posterior$posteriorFunctions = list(
    U = 'U_F %*% Lambda + U_R + X1 %*% B1',
    G = 't(Lambda) %*% diag(F_h2[1,]) %*% Lambda + diag(resid_h2[1,]/tot_Eta_prec[1,])',
    R = 't(Lambda) %*% diag(1-F_h2[1,]) %*% Lambda + diag((1-resid_h2[1,])/tot_Eta_prec[1,])',
    h2 = '(colSums(F_h2[1,]*Lambda^2)+resid_h2[1,]/tot_Eta_prec[1,])/(colSums(Lambda^2)+1/tot_Eta_prec[1,])'
  )
  for(i in 1:5) {
    MegaLMM_state = reorder_factors(MegaLMM_state,
                                    drop_cor_threshold = 0.6)
    MegaLMM_state = clear_Posterior(MegaLMM_state)
    MegaLMM_state = sample_MegaLMM(MegaLMM_state,
                                   nIter)
  }
  MegaLMM_state = clear_Posterior(MegaLMM_state)
  
  for(i in 1:4) {
    MegaLMM_state = sample_MegaLMM(MegaLMM_state,
                                   nIter) 
    MegaLMM_state = save_posterior_chunk(MegaLMM_state)
  }
  U_samples = load_posterior_param(MegaLMM_state,
                                   'U')
  U_hat = get_posterior_mean(U_samples)
  return(U_hat)
}

if (mode == "fit") {
  ans <- mega(Y)
  
  pred[1:nrow(pred), ] <- ans
  
  fold_id <- rep(0L, nrow(Y))
} else if (mode == "loo") {
  for (i in 1:nrow(Y)) {
    Y_test <- Y
    Y_test[i, ] <- NA
    
    ans <- mega(Y_test)
    pred[i, ] <- ans[i, ]
    
    fold_id[i] <- i
    if (i %% 20 == 0) cat(" loo ", i, "/", nrow(Y), "\n")
  }
} else if (mode == "kfold") {
  n <- nrow(Y)
  perm <- sample.int(n)
  folds <- split(perm, rep(seq_len(k_folds), length.out = n))
  for (k in seq_len(k_folds)) {
    test_idx <- folds[[k]]
    
    Y_test <- Y
    Y_test[test_idx, ] <- NA
    
    ans <- mega(Y_test)
    pred[test_idx, ] <- ans[test_idx, ]
    
    fold_id[test_idx] <- k
    cat(" fold ", k, "/", k_folds, "\n")
  }
}

out_tsv <- file.path(out_dir, "predictions_all.tsv")
DT <- data.frame(
  id = rep(ids_match, length(trait_target)),
  trait_idx = rep(trait_idxs, each = nrow(Y)),
  trait = rep(trait_target, each = nrow(Y)),
  y_true = as.vector(as.matrix(Y)),
  y_pred = as.vector(as.matrix(pred)),
  stringsAsFactors = FALSE
)
fwrite(DT, 
       out_tsv, 
       sep = "\t")

corv <- c()
rmsev <- c()
maev <- c()
for (cyc1 in 1:length(trait_target)) {
  corv <- c(corv,
            pearson(DT$y_true[DT$trait == trait_target[cyc1]], DT$y_pred[DT$trait == trait_target[cyc1]]))
  rmsev <- c(rmsev,
             rmse(DT$y_true[DT$trait == trait_target[cyc1]], DT$y_pred[DT$trait == trait_target[cyc1]]))
  maev <- c(maev,
            mae(DT$y_true[DT$trait == trait_target[cyc1]], DT$y_pred[DT$trait == trait_target[cyc1]]))
}

summary_rows[[length(summary_rows) + 1]] <- data.frame(
  trait_idx = trait_idxs,
  trait = trait_target,
  n = nrow(Y),
  #n_obs = sum(is.finite(y)),
  mode = mode,
  backend = "MegaLMM",
  cor = corv,
  rmse = rmsev,
  mae = maev,
  stringsAsFactors = FALSE
)
#artifact_tables[[trait_target]] <- basename(out_tsv)

summary_df <- rbindlist(summary_rows, fill = TRUE)
summary_path <- file.path(out_dir, "summary.tsv")
fwrite(summary_df, summary_path, sep = "\t")

metrics <- list(
  mode = mode,
  n = length(ids_match),
  n_markers = ncol(geno),
  n_traits = length(trait_target),
  traits = trait_target,
  backend = "MegaLMM",
  K_factors = K_factors,
  seed = seed
)
writeLines(toJSON(metrics, auto_unbox = TRUE, pretty = TRUE), file.path(out_dir, "metrics.json"))

for (cyc1 in 1:length(trait_target)) {
  out_tsv <- file.path(out_dir, paste0("predictions_", trait_target[cyc1], ".tsv"))
  DT <- data.frame(
    id = ids_match,
    trait_idx = trait_idxs[cyc1],
    trait = trait_target[cyc1],
    y_true = Y[, cyc1],
    y_pred = pred[, cyc1],
    stringsAsFactors = FALSE
  )
  fwrite(DT, 
         out_tsv, 
         sep = "\t")
}

# Make a simple plot: trait-wise cor barplot (base R)
try({
  png(file.path(out_dir, "trait_cor.png"), width = 900, height = 450)
  par(mar=c(7,4,2,1))
  barplot(summary_df$cor, names.arg = summary_df$trait, las = 2)
  abline(h=0)
  dev.off()
}, silent = TRUE)

art <- list(
  tables = c("summary.tsv", unname(unlist(artifact_tables))),
  plots = c("trait_cor.png"),
  default_table = "summary.tsv",
  default_plot = "trait_cor.png"
)
writeLines(toJSON(art, auto_unbox = TRUE, pretty = TRUE), file.path(out_dir, "artifacts.json"))

cat("[gp_megalmm] done\n")
sink()
