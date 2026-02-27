<img width="1834" height="1035" alt="Screenshot from 2026-02-26 16-01-30" src="https://github.com/user-attachments/assets/12250a31-e7f8-42d4-b4e1-4e06d0438198" />

# GenomicExplorer: a GUI software for statistical genomics

GenomicExplorer combines many R and Python packages and third party tools, and utilizes those for statistical genomics. 
I recommend the use of Micromamba to prepare GenomicExplorer.

## Prerequisites

GenomicExplorer needs Micromamba.
The detailed procedure to install Micromamba is described in the webpage.

[Micromamba webpage](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html)

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

The most functions of GenomicExplorer should work in Windows OS. Although GenomicExplorer can work with PowerShell, I prepared the shortcut to install and start GenomicExplorer.

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

The popular and reliable packages/tools can be used with GenomicExplorer. I will prepare the detailed manual of GenomicExplorer for these packages in future. Please read those manuals or documents to understand the detailed usage of each packages/tool.

### Preprocess
The Preprocess tab implements File converter, Genotype Imputation and Phenotype Imputation.
File converter accepts VCF/PLINK BED/TSV/CSV or RDS of r/qtl and qtl2 YAML and convertd into the other formats.

In Genotype Imputation, users can select a method from mean/EM/RF (RandomForest). Users can prepare the Beagle v5.4 software (Browning et al., 2018) and use it for Genotype Imputation.

[rrBLUP](https://cran.r-project.org/web/packages/rrBLUP/index.html)
[missRanger](https://cran.r-project.org/web/packages/missRanger/index.html)  
[Beagle v5.4](https://faculty.washington.edu/browning/beagle/b5_4.html)

[PHENIX](https://github.com/jlboat/PHENIX?tab=readme-ov-file)

### GP (Genomic Prediction)

Both single trait and multi-trait models are implemented.

GenomicExplorer mainly depends on BGLR, a popular R package for GP (Perez and de los Campos, 2014). 
Basic priors of BGLR, such as BRR/BayesA/BayesB/BayesC/RKHS, can be used with single kernel or multiple kernels in GenomicExplorer.

MegaLMM (Runcie et al., 2021) can also be used for multi-trait model.

In addition, RandomForest, an non-linear model, can also be used with GenomicExplorer. 
For single trait, a fast package, ranger, is adopted in GenomicExplorer. For multi-trait model, randomForestSRC is used in GenomicExplorer.

[BGLR github](https://github.com/gdlc/BGLR-R)    
[MegaLMM github](https://github.com/deruncie/MegaLMM/tree/master)    
[ranger](https://cran.r-project.org/web/packages/ranger/index.html)    
[ranfomForestSRC](https://cran.r-project.org/web/packages/randomForestSRC/index.html)

[RAINBOWR](https://github.com/KosukeHamazaki/RAINBOWR)

### GWAS (Genome-Wide Association Study)

Linear Mixed Model (LMM) and other methods are implemented. The R package gaston is used for LMM.

[gaston](https://cran.r-project.org/web/packages/gaston/index.html)  
[MLMM](https://github.com/Gregor-Mendel-Institute/MultLocMixMod)  
[FarmCPU](https://github.com/amkusmec/FarmCPUpp)

### RNA-seq (Transcriptome analysis)

<img width="458" height="422" alt="Screenshot from 2026-02-26 16-05-41" src="https://github.com/user-attachments/assets/ba66c7a4-cdd6-4203-a0db-2d93d2a70302" />

DEG (Differentially Expressed Genes), sample clustering/dimension reduction, gene clustering model and GO (Gene Ontology) enrichment analyses are implemented.

For non-model organisms, GO DB (DataBase) builder is also implemented.

[Bioconductor](https://www.bioconductor.org/)  
[GOstats](https://www.bioconductor.org/packages/release/bioc/html/GOstats.html)  

### Map (Linkage Map Construction)

Three R packages for linkage map construction can be used in GenomicExplorer.

GenomicExplorer depends on three R packages (ASMap/qtl/onemap) for linkage map construction.

(In future, Lep-MAP3 may be used with GenomicExplorer)

[onemap](https://cran.r-project.org/web/packages/onemap/index.html)

### QTL (Quantitative Trait Loci analysis)

QTL (Quantitative Trait Loci) analysis can be performed with GenomicExplorer.

GenomicExplorer depends on three R packages (qtl/qtl2/qtlbim) for QTL analysis.

[qtl2 github](https://github.com/rqtl/qtl2)  
[qtlbim github](https://github.com/fboehm/qtlbim)

### Polyploid (Dosage calling, GWAS, Map and QTL)

For polyploids, GenomicExplorer prepares an independent tab.

[updog](https://github.com/dcgerard/updog)  
[GWASpoly](https://github.com/jendelman/GWASpoly)  
[MAPpoly](https://github.com/mmollina/MAPpoly)  
[QTLpoly](https://github.com/guilherme-pereira/QTLpoly)

### Integrate

Various results derived from GWAS/QTL/RNA-seq/others can be integrated.

### Sim (Simulation of genotype/phenotype and RNA-seq count data)

The Sim tab can generate genotype/phenotype and RNA-seq count data. This is mainly used for debug.
