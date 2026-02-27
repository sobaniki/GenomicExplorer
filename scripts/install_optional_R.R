#!/usr/bin/env Rscript

u <- Sys.getenv("R_LIBS_USER")
if (nzchar(u)) {
  dir.create(u, showWarnings = FALSE, recursive = TRUE)
  .libPaths(c(normalizePath(u, winslash="\\", mustWork=FALSE), .libPaths()))
}

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

# ---- packages inferred from GE plugins (core vs heavy) ----
# core: 比較的入れやすく、機能の穴埋めに効く
cran_core <- c(
  "ASMap", "onemap", "mappoly", "AlphaSimR", "updog", "rrBLUP", "RAINBOWR", "NBPSeq", "missRanger", "softImpute", "samr"
)

bioc_core <- c(
  "GBScleanR", "TCC", "clusterSeq", "org.At.tair.db", "baySeq"
)

dev_github <- c("hadley/pryr",
		"fboehm/qtlbim",
		"Gregor-Mendel-Institute/MultLocMixMod",
		"amkusmec/FarmCPUpp",
		"jlboat/PHENIX",
		"deruncie/MegaLMM",
		"jendelman/GWASpoly")
#For MLMM
install.packages("https://github.com/Gregor-Mendel-Institute/mlmm/files/1356516/emma_1.1.2.tar.gz", repos = NULL)
install.packages("https://cran.r-project.org/src/contrib/Archive/MBCluster.Seq/MBCluster.Seq_1.0.tar.gz", repos = NULL)
install.packages("https://cran.r-project.org/src/contrib/Archive/qtlpoly/qtlpoly_0.2.4.tar.gz", repos = NULL)
install.packages("https://cran.r-project.org/src/contrib/Archive/impute/impute_1.26.0.tar.gz", repos = NULL)

# ---- profile selection ----
if (profile == "core") {
  install_cran(cran_core)
  install_bioc(bioc_core)
} else {
  install_cran(cran_core)
  install_bioc(bioc_core)
  install_github(dev_github)
}

cat("[install_optional_R] done.\n")
