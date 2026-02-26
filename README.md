<img width="1834" height="1035" alt="Screenshot from 2026-02-26 16-01-30" src="https://github.com/user-attachments/assets/12250a31-e7f8-42d4-b4e1-4e06d0438198" />

# GenomicExplorer: a GUI software for statistical genomics

GenomicExplorer combines many R and Python packages and third party tools, and utilizes those for statistical genomics. 
I recommend the use of Micromamba to prepare GenomicExplorer.

## Prerequisites

GenomicExplorer needs Micromamba.
The detailed procedure to install Micromamba is described in the webpage.

Here, I show the procudure using terminal in Linux OS (Ubuntu):

```install_micromamba
curl -fsSL https://micro.mamba.pm/install.sh -o ~/install_micromamba.sh
bash ~/install_micromamba.sh
```

## Download and Install

In Linux, GenomicExplorer can be downloaded and installed using terminal:

```install_genomicexplorer
git clone https://github.com/sobaniki/GenomicExplorer.git
cd GenomicExplorer
bash ./scripts/install_linux.sh
bash ./scripts/install_optional_R.sh
```

## To use interactive plot panel

In GenomicExplorer, some useful interactive plots can be drawn. PySide6 via micromamba or conda is often unstable. Please re-install PySide6 with pip. 

```use_plotly
micromamba activate GenomicExplorer
micromamba remove -y pyside6  || true
pip install --upgrade PySide6 PySide6-Addons
python -c "from PySide6.QtWebEngineWidgets import QWebEngineView; print('QtWebEngine OK')"
```

## Run

In Linux, GenomicExplorer can be started using terminal:

```start_genomicexplorer
bash ./scripts/run_linux.sh
```

## For Windows OS users

The most functions of GenomicExplorer should work in Windows OS.

Please click "install_genomicexplorer.cmd", "install_optional_R.cmd", and "run_genomicexplorer.cmd" to install and start GenomicExplorer.

## File format
In GenomicExplorer, several file format for genotype data can be used.

GenomicExplorer expects VCF (Varinat Call Format), PLINK BED or TSV/CSV formats for many genomics applications.
The genotype TSV/CSV files compose of n (samples as row) x m (markers as column).
The rownames or first column should be sample (individual) names/IDs.
The genotype code should be 0/1/2 (1 is heterozygous) for diploids.

The marker (map) TSV/CSV files include marker names/IDs, chromosome names/IDs, and positions in chromosome (BP or cM) as the first three columns.

The phenotype/covariate TSV/CSV files include sample (individual) names/IDs as the rownames or first column. The following columns are treated as phenotypes/covariates.

The some original formats or files may be necessary for Map/QTL, which can be prepared with GenomicExplorer.

For RNA-seq, count data compose of m (genes as row) x n (samples as column) as TSV/CSV formats. Design data compose of n (samples as row) x m (conditions as column) as TSV/CSV formats.

## How to use GenomicExplorer

GenomicExplorer includes various functions for statistical genomics. These functions are separated as tabs (Preprocess/GP/GWAS/RNA-seq/Map/QTL/Polyploid/Integrate/Sim) in GenomicExplorer.

The popular and reliable packages/tools can be used with GenomicExplorer. Please read those manuals or documents to understand the detailed usage of each packages/tool.

### Preprocess
The Preprocess tab implements File converter, Genotype Imputation and Phenotype Imputation.

### GP (Genomic Prediction)

Both single trait and multi-trait models are implemented.

GenomicExplorer mainly depends on BGLR, a popular R package for GP. In addition, RandomForest, an non-linear model, can also be used with GenomicExplorer.

### GWAS (Genome-Wide Association Study)

Linear Mixed Model (LMM) and other methods are implemented.

### RNA-seq

<img width="916" height="843" alt="Screenshot from 2026-02-26 16-05-41" src="https://github.com/user-attachments/assets/ba66c7a4-cdd6-4203-a0db-2d93d2a70302" />

DEG (Differentially Expressed Genes), sample clustering/dimension reduction, gene clustering model and GO (Gene Ontology) enrichment analyses are implemented.

For non-model organisms, GO DB (DataBase) builder is also implemented.

### Map

Three R packages for linkage map construction can be used in GenomicExplorer.

GenomicExplorer depends on three R packages (ASMap/qtl/onemap) for linkage map construction.

(In future, Lep-MAP3 may be used with GenomicExplorer)

### QTL

QTL (Quantitative Trait Loci) analysis can be performed with GenomicExplorer.

GenomicExplorer depends on three R packages (qtl/qtl2/qtlbim) for QTL analysis.

### Polyploid

For polyploids, GenomicExplorer prepares an independent tab.

### Integrate

Various results derived from GWAS/QTL/RNA-seq/others can be integrated.

### Sim

The Sim tab can generate genotype and RNA-seq count data. This is mainly used for debug.
