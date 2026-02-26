# GenomicExplorer: a GUI software for statistical genomics.

GenomicExplorer combines many R and Python packages and third party tools, and utilizes those for statistical genomics. 
I recommend the use of Micromamba to prepare GenomicExplorer.

## **Prerequisites**

GenomicExplorer needs Micromamba.
The detailed procedure to install Micromamba is described in the webpage.

Here, I show the procudure using terminal in Linux OS (Ubuntu):

```install_micromamba
curl -fsSL https://micro.mamba.pm/install.sh -o ~/install_micromamba.sh
bash ~/install_micromamba.sh
```

## **Download and Install**

In Linux, GenomicExplorer can be downloaded and installed using terminal:

```install_genomicexplorer
git clone https://github.com/sobaniki/GenomicExplorer.git
cd GenomicExplorer
bash ./scripts/install_linux.sh
bash ./scripts/install_optional_R.sh
```

## **Run**

In Linux, GenomicExplorer can be started using terminal:

```start_genomicexplorer
bash ./scripts/run_linux.sh
```

## **How to use GenomicExplorer**

GenomicExplorer includes various functions for statistical genomics. These functions are separated as tabs (Preprocess/GP/GWAS/RNA-seq/Map/QTL/Polyploid/Integrate/Sim) in GenomicExplorer.

### Preprocess
The Preprocess tab implements File converter, Genotype Imputation and Phenotype Imputation.

### GP

### GWAS

### RNA-seq

### Map

### QTL

### Polyploid

### Integrate

### Sim
