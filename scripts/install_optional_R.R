#!/usr/bin/env Rscript

u <- Sys.getenv("R_LIBS_USER")
#if (nzchar(u)) {
#  dir.create(u, showWarnings = FALSE, recursive = TRUE)
#  .libPaths(c(normalizePath(u, winslash="\\", mustWork=FALSE), .libPaths()))
#}

suppressWarnings(suppressMessages({
  args <- commandArgs(trailingOnly = TRUE)
}))

get_arg <- function(flag, default=NULL) {
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) return(args[i+1])
  default
}

profile <- get_arg("--profile", "full")

cat("[install_optional_R] profile =", profile, "\n")
#cat("[install_optional_R] R_LIBS_USER =", Sys.getenv("R_LIBS_USER"), "\n")
cat("[install_optional_R] .libPaths() =", paste(.libPaths(), collapse=" | "), "\n")

# --- installers ---
install_cran <- function(pkgs) {
  pkgs <- pkgs[!pkgs %in% rownames(installed.packages())]
  if (length(pkgs) == 0) return(invisible(TRUE))
  install.packages(pkgs, repos="https://cloud.r-project.org")
}

install_cran_ver <- function(pkgs, vers) {
  for (cyc1 in 1:length(pkgs)) {
    remotes::install_version(pkgs[cyc1], version = vers[cyc1], repos = "https://cloud.r-project.org")
  }
}

install_bioc <- function(pkgs) {
  if (!requireNamespace("BiocManager", quietly=TRUE)) {
    install.packages("BiocManager", repos="https://cloud.r-project.org")
  }
  pkgs <- pkgs[!pkgs %in% rownames(installed.packages())]
  if (length(pkgs) == 0) return(invisible(TRUE))
  BiocManager::install(pkgs, ask=FALSE, update=FALSE)
}

install_github <- function(repos) {
  if (!requireNamespace("devtools", quietly=TRUE)) {
    install.packages("devtools", repos="https://cloud.r-project.org")
  }
  for (r in repos) {
    pkg <- sub(".*/", "", r)
    if (pkg %in% rownames(installed.packages())) next
    devtools::install_github(r, upgrade="never", dependencies=TRUE)
  }
}

install_github_ver <- function(pkgs, vers) {
  for (cyc1 in 1:length(pkgs)) {
    devtools::install_github(pkgs[cyc1], ref = vers[cyc1], upgrade="never", dependencies=TRUE)
  }
}

# ---- packages inferred from GE plugins (core vs heavy) ----
# core: 比較的入れやすく、機能の穴埋めに効く
cran_core <- c(
  "ASMap", "onemap", "mappoly", "AlphaSimR", "updog", 
  "rrBLUP", "RAINBOWR", "NBPSeq", "missRanger", "softImpute", 
  "samr", "CompQuadForm", "RLRsim", "quadprog", "doParallel", "lobstr", "coda"
)
cran_version <- c("1.0-8", "3.2.4", "0.4.2", "2.1.0", "2.1.6",
                  "4.6.3", "0.1.38", "0.3.1", "2.6.1", "1.4-3",
                  "3.0", "1.4.4", "3.1-9", "1.5-8", "1.0.17", "1.2.0", "0.19-4.1")

bioc_core <- c(
  #"GBScleanR", 
  #"TCC", 
  #"clusterSeq", 
  #"org.At.tair.db", 
  #"baySeq"
)

dev_github <- c(#"hadley/pryr",
		#"fboehm/qtlbim",
		"Gregor-Mendel-Institute/MultLocMixMod",
		"amkusmec/FarmCPUpp",
		"jlboat/PHENIX",
		"deruncie/MegaLMM",
		"jendelman/GWASpoly")
dev_github_ver <- c(#"0.1.6.9000",
                    #"2.1.0.9000",
                    "0.1.1",
                    "1.2.0",
                    "1.0.1",
                    "0.9.5",
                    "2.14")
                    
install.packages("https://github.com/Gregor-Mendel-Institute/mlmm/files/1356516/emma_1.1.2.tar.gz", repos = NULL)

# ---- profile selection ----
if (profile == "core") {
  #install_cran(cran_core)
  #install_bioc(bioc_core)
  install_cran_ver(cran_core, cran_version)
} else {
  #install_cran(cran_core)
  #install_bioc(bioc_core)
  #install_github(dev_github)
  install_cran_ver(cran_core, cran_version)
  #install_github_ver(dev_github, dev_github_ver)
  install_github(dev_github)
}

install.packages("https://cran.r-project.org/src/contrib/Archive/pryr/pryr_0.1.6.tar.gz", repos = NULL)
install.packages("https://cran.r-project.org/src/contrib/Archive/qtlbim/qtlbim_2.0.7.tar.gz", repos = NULL)
install.packages("https://cran.r-project.org/src/contrib/Archive/MBCluster.Seq/MBCluster.Seq_1.0.tar.gz", repos = NULL)
install.packages("https://cran.r-project.org/src/contrib/Archive/qtlpoly/qtlpoly_0.2.2.tar.gz", repos = NULL)
install.packages("https://cran.r-project.org/src/contrib/Archive/impute/impute_1.26.0.tar.gz", repos = NULL)

cat("[install_optional_R] done.\n")
