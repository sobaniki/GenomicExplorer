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

In Genotype Imputation, users can select a method from mean/EM (Endelman, 2011)/RF (Random Forest). Users can prepare the Beagle v5.4 software (Browning et al., 2018) and use it for Genotype Imputation.   
For phenotype imputation, GenomicExplorer supports mean and median imputation, missForest (Stekhoven and Buhlmann, 2012), mice (van Buuren and Groothuis-Oudshoorn, 2011), PHENIX (Dahl et al., 2016) and softImpute (Mazumder et al., 2010).

[rrBLUP](https://cran.r-project.org/web/packages/rrBLUP/index.html)   
[missRanger](https://cran.r-project.org/web/packages/missRanger/index.html)    
[Beagle v5.4](https://faculty.washington.edu/browning/beagle/b5_4.html)   

[PHENIX](https://github.com/jlboat/PHENIX?tab=readme-ov-file)

### GP (Genomic Prediction)

<img width="300" height="300" alt="gp_scatter_2x2_plotly" src="https://github.com/user-attachments/assets/a908ef65-76dd-4fbd-9db6-65f3a01c3e2c" />

Both single trait and multi-trait models are implemented.

GenomicExplorer mainly depends on BGLR, a popular R package for GP (Perez and de los Campos, 2014). 
Basic priors of BGLR, such as BRR/BayesA/BayesB/BayesC/RKHS, can be used with single kernel or multiple kernels in GenomicExplorer.   
Some kernels (e.g., additive relationship matrix) can be calculated using RAINBOWR (Hamazaki and Iwata, 2020).

Using glmnet (Friedman et al., 2010), Ridge and Lasso regressions are also supported.

In addition to BGLR (Pérez-Rodríguez and de los Campos, 2022), MegaLMM (Runcie et al., 2021) can also be used for multi-trait model.

Random Forest, an non-linear model, can also be used with GenomicExplorer. 
For single trait, a fast package, ranger (Wright and Ziegler, 2017), is adopted in GenomicExplorer.   
For multi-trait model, randomForestSRC is used in GenomicExplorer.  

[BGLR github](https://github.com/gdlc/BGLR-R)    
[MegaLMM github](https://github.com/deruncie/MegaLMM/tree/master)    
[ranger](https://cran.r-project.org/web/packages/ranger/index.html)    
[ranfomForestSRC](https://cran.r-project.org/web/packages/randomForestSRC/index.html)

[RAINBOWR](https://github.com/KosukeHamazaki/RAINBOWR)

### GWAS (Genome-Wide Association Study)
<img width="350" height="200" alt="manhattan" src="https://github.com/user-attachments/assets/a7bc9abc-a625-4756-8644-2300ae0241fd" />

Linear Mixed Model (LMM) and other methods are implemented. The R package gaston is used for LMM.   
Multi-locus models, such as MLMM  (Segura et al., 2012) and FarmCPUpp (Kusmec and Schnable, 2018), can be performed in GenomicExplorer.

[gaston](https://cran.r-project.org/web/packages/gaston/index.html)  
[MLMM](https://github.com/Gregor-Mendel-Institute/MultLocMixMod)  
[FarmCPU](https://github.com/amkusmec/FarmCPUpp)

### RNA-seq (Transcriptome analysis)

<img width="229" height="211" alt="Screenshot from 2026-02-26 16-05-41" src="https://github.com/user-attachments/assets/ba66c7a4-cdd6-4203-a0db-2d93d2a70302" />

DEG (Differentially Expressed Genes) (Robinson et al., 2010; Love et al., 2014; Law et al., 2014; Ritchie et al., 2015; Tarazona et al., 2011; Leng et al., 2013), sample clustering/dimension reduction (Scrucca et al., 2023), model-based gene clustering (Rau et al., 2015; Si et al., 2014; Rau and Maugis-Rabusseau, 2018) and GO (Gene Ontology) enrichment analyses are implemented.   
R packages for enrichment analyses are generally available from Bioconductor (Falcon and Gentleman, 2007; Yu et al., 2012; Young et al., 2010).   
For non-model organisms, GO DB (DataBase) builder is also implemented.   

[Bioconductor](https://www.bioconductor.org/)  
[GOstats](https://www.bioconductor.org/packages/release/bioc/html/GOstats.html)  

### Map (Linkage Map Construction)

<img width="350" height="225" alt="map_plotly" src="https://github.com/user-attachments/assets/acc5f135-ac75-4927-a1a1-8d783cb5e183" />

Linkage map construction can be roughly performed in GenomicExplorer.
GenomicExplorer depends on three R packages (ASMap/Rqtl/onemap) for linkage map construction(Taylor and Butler, 2017; Broman et al., 2003; Margarido et al., 2007).   
Lep-MAP3 (Rastas, 2017) can be used with GenomicExplorer if user provide the Java distribution.   

Please note that constructing linkage maps typically requires careful manual work, and GenomicExplorer currently offers limited options for creating highly accurate linkage maps.

[onemap](https://cran.r-project.org/web/packages/onemap/index.html)
[Lep-MAP3](https://sourceforge.net/p/lep-map3/wiki/LM3%20Home/)

### QTL (Quantitative Trait Loci analysis)

QTL (Quantitative Trait Loci) analysis can be performed with GenomicExplorer.

GenomicExplorer depends on three R packages (Rqtl/Rqtl2/qtlbim) for general QTL analysis (Broman et al., 2003, 2019; Yandell et al, 2007).   
For NAM-like populations, the mppR package (Garin et al., 2024) can be used in GenomicExplorer.    

[qtl2 github](https://github.com/rqtl/qtl2)  
[qtlbim github](https://github.com/fboehm/qtlbim)   
[mppR github](https://github.com/vincentgarin/mppR)

### Polyploid (Dosage calling, GWAS, Map and QTL)

For polyploids, GenomicExplorer prepares an independent tab.   
Genotype dosage calling, GWAS, linkage map construction and QTL scan are supported using R packgaes (Gerard et al., 2018; Mollinari and Garcia, 2019; da Silva Pereira et al., 2020; Rosyara et al., 2016).   

[updog](https://github.com/dcgerard/updog)  
[GWASpoly](https://github.com/jendelman/GWASpoly)  
[MAPpoly](https://github.com/mmollina/MAPpoly)  
[QTLpoly](https://github.com/guilherme-pereira/QTLpoly)

### Integrate

Various results derived from GWAS/QTL/RNA-seq/others can be integrated.

### Sim (Simulation of genotype/phenotype and RNA-seq count data)

The Sim tab can generate genotype/phenotype and RNA-seq count data. This is mainly used for debug.   
Bi-parental crosses, such as F2 and RIL, inbred and outbred population and others are supported.    
GenomicExplorer use AlphaSimR (Gaynor et al., 2021) for simulation of genotype/phenotype data.

## Paper
Ishimori, M. (2026) GenomicExplorer: a graphical user interface platform for statistical genetics and genomics. (submitted)

## Citation
Broman, K.W., D.M. Gatti, P. Simecek, N.A. Furlotte, P. Prins, S. Sen, B.S. Yandell and G.A. Churchill (2019) R/qtl2: Software for mapping quantitative trait loci with high-dimensional data and multiparent populations. Genetics 211: 495–502.   

Broman, K.W., H. Wu, S. Sen and G.A. Churchill (2003) R/qtl: QTL mapping in experimental crosses. Bioinformatics 19: 889–890.   

Browning, B.L., Y. Zhou and S.R. Browning (2018) A one-penny imputed genome from next-generation reference panels. Am. J. Hum. Genet. 103: 338–348.   

Cantalapiedra, C.P., A. Hernández-Plaza, I. Letunic, P. Bork and J. Huerta-Cepas (2021) eggNOG-mapper v2: Functional annotation, orthology assignments, and domain prediction at the metagenomic scale. Mol. Biol. Evol. 38: 5825–5829.   

Chang, C.C., C.C. Chow, L.C. Tellier, S. Vattikuti, S.M. Purcell and J.J. Lee (2015) Second-generation PLINK: rising to the challenge of larger and richer datasets. GigaScience 4: 7.   

Dahl, A., V. Iotchkova, A. Baud, Å. Johansson, U. Gyllensten, N. Soranzo, R. Mott, A. Kranis and J. Marchini (2016) A multiple-phenotype imputation method for genetic studies. Nat. Genet. 48: 466–472.   

da Silva Pereira, G., D.C. Gemenet, M. Mollinari, B.A. Olukolu, J.C. Wood, F. Diaz, V. Mosquera, W.J. Gruneberg, A. Khan, C.R. Buell, et al. (2020) Multiple QTL mapping in autopolyploids: a random-effect model approach with application in a hexaploid sweetpotato full-sib population. Genetics 215: 579–595.   

Endelman, J.B. (2011) Ridge regression and other kernels for genomic selection with R package rrBLUP. Plant Genome 4: 250–255.   

Falcon, S., and R. Gentleman (2007) Using GOstats to test gene lists for GO term association. Bioinformatics 23: 257–258.   

Friedman, J., T. Hastie and R. Tibshirani (2010) Regularization paths for generalized linear models via coordinate descent. J. Stat. Softw. 33: 1–22.   

Gabriel, L., T. Brüna, K.J. Hoff, M. Ebel, A. Lomsadze, M. Borodovsky and M. Stanke (2024) BRAKER3: fully automated genome annotation using RNA-seq and protein evidence with GeneMark-ETP, AUGUSTUS and TSEBRA. Genome Res. 34: 769–777.   

Gaynor, R.C., G. Gorjanc and J.M. Hickey (2021) AlphaSimR: an R package for breeding program simulations. G3 (Bethesda) 11: jkaa017.   

Gerard, D., L.F.V. Ferrao, A.A.F. Garcia and M. Stephens (2018) Genotyping polyploids from messy sequencing data. Genetics 210: 789–807.   

Hamazaki, K. and H. Iwata (2020) RAINBOW: Haplotype-based genome-wide association study using a novel SNP-set method. PLoS Comput. Biol. 16: e1007663.   

Korf, I. (2004) Gene finding in novel genomes. BMC Bioinformatics 5: 59.   

Kusmec, A. and P.S. Schnable (2018) FarmCPUpp: Efficient large-scale genomewide association studies. Plant Direct 2: e00053.   

Law, C.W., Y. Chen, W. Shi and G.K. Smyth (2014) voom: precision weights unlock linear model analysis tools for RNA-seq read counts. Genome Biol. 15: R29.   

Leng, N., J.A. Dawson, J.A. Thomson, V. Ruotti, A.I. Rissman, B.M.G. Smits, J.D. Haag, M.N. Gould, R.M. Stewart and C. Kendziorski (2013) EBSeq: an empirical Bayes hierarchical model for inference in RNA-seq experiments. Bioinformatics 29: 1035–1043.   

Love, M.I., W. Huber and S. Anders (2014) Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biol. 15: 550.   

Margarido, G.R.A., A.P. Souza and A.A.F. Garcia (2007) OneMap: software for genetic mapping in outcrossing species. Hereditas 144: 78–79.   

Mazumder, R., T. Hastie and R. Tibshirani (2010) Spectral regularization algorithms for learning large incomplete matrices. J. Mach. Learn. Res. 11: 2287–2322.  

McGinnis, S. and T.L. Madden (2004) BLAST: at the core of a powerful and diverse set of sequence analysis tools. Nucleic Acids Res. 32: W20–W25.   

Mollinari, M. and A.A.F. Garcia (2019) Linkage analysis and haplotype phasing in experimental autopolyploid populations with high ploidy level using hidden Markov models. G3 (Bethesda): 3297–3314.   

Pérez, P. and G. de los Campos (2014) Genome-wide regression and prediction with the BGLR statistical package. Genetics 198: 483–495.   

Pérez-Rodríguez, P. and G. de los Campos (2022) Multitrait Bayesian shrinkage and variable selection models with the BGLR-R package. Genetics 222: iyac112.  

Purcell, S., B. Neale, K. Todd-Brown, L. Thomas, M.A.R. Ferreira, D. Bender, J. Maller, P. Sklar, P.I.W. de Bakker, M.J. Daly and P.C. Sham (2007) PLINK: a tool set for whole-genome association and population-based linkage analyses. Am. J. Hum. Genet. 81: 559–575.   

Rastas, P. (2017) Lep-MAP3: robust linkage mapping even for low-coverage whole-genome sequencing data. Bioinformatics 33: 3726–3732.   

Rau, A. and C. Maugis-Rabusseau (2018) Transformation and model choice for RNA-seq co-expression analysis. Brief. Bioinform. 19: 425–436.   

Rau, A., C. Maugis-Rabusseau, M.-L. Martin-Magniette and G. Celeux (2015) Co-expression analysis of high-throughput transcriptome sequencing data with Poisson mixture models. Bioinformatics 31: 1420–1427.   

Ritchie, M.E., B. Phipson, D. Wu, Y. Hu, C.W. Law, W. Shi and G.K. Smyth (2015) limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Res. 43: e47.   

Robinson, M.D., D.J. McCarthy and G.K. Smyth (2010) edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics 26: 139–140.   

Rosyara, U.R., W.S. De Jong, D.S. Douches and J.B. Endelman (2016) Software for genome-wide association studies in autopolyploids and its application to potato. Plant Genome 9: 1–10.   

Runcie, D.E., J. Qu, H. Cheng and L. Crawford (2021) MegaLMM: Mega-scale linear mixed models for genomic predictions with thousands of traits. Genome Biol. 22: 213.   

Scrucca, L., C. Fraley, T.B. Murphy and A.E. Raftery (2023) Model-based clustering, classification, and density estimation using mclust in R. Chapman and Hall/CRC, New York, p. 268.   

Segura, V., B.J. Vilhjalmsson, A. Platt, A. Korte, Ü. Seren, Q. Long and M. Nordborg (2012) An efficient multi-locus mixed-model approach for genome-wide association studies in structured populations. Nat. Genet. 44: 825–830.   

Si, Y., P. Liu, P. Li and T.P. Brutnell (2014) Model-based clustering for RNA-seq data. Bioinformatics 30: 197–205.   

Sievert, C. (2020) Interactive Web-Based Data Visualization with R, Plotly, and Shiny. Chapman and Hall/CRC, New York, p. 470.   

Stanke, M. and B. Morgenstern (2005) AUGUSTUS: a web server for gene prediction in eukaryotes that allows user-defined constraints. Nucleic Acids Res. 33: W465–W467.   

Stekhoven, D. J. and P. Bühlmann (2012) MissForest—non-parametric missing value imputation for mixed-type data. Bioinformatics 28: 112–118.   

Tarazona, S., F. García-Alcalde, J. Dopazo, A. Ferrer and A. Conesa (2011) Differential expression in RNA-seq: a matter of depth. Genome Res. 21: 2213–2223.  

Taylor, J. and D. Butler (2017) R package ASMap: Efficient genetic linkage map construction and diagnosis. J. Stat. Softw. 79: 1–29.   

van Buuren, S. and K. Groothuis-Oudshoorn (2011) mice: Multivariate imputation by chained equations in R. J. Stat. Softw. 45: 1–67.   

Wright, M.N. and A. Ziegler (2017) ranger: a fast implementation of random forests for high dimensional data in C++ and R. J. Stat. Softw. 77: 1–17.   

Yandell, B.S., T. Mehta, S. Banerjee, D. Shriner, R. Venkataraman, J.Y. Moon, W.W. Neely, H. Wu, R. von Smith and N. Yi (2007) R/qtlbim: QTL with Bayesian interval mapping in experimental crosses. Bioinformatics 23: 641–643.   

Young, M.D., M.J. Wakefield, G.K. Smyth and A. Oshlack (2010) Gene ontology analysis for RNA-seq: accounting for selection bias. Genome Biol. 11: R14.  

Yu, G., L.-G. Wang, Y. Han and Q.-Y. He (2012) clusterProfiler: an R package for comparing biological themes among gene clusters. OMICS 16: 284–287.  


