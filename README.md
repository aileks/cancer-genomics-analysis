
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Cancer Genomics Analysis

<!-- badges: start -->
<!-- badges: end -->

## Overview

Analysis of TCGA’s acute myeloid leukemia data (TCGA-LAML) to explore
mutation patterns, gene expression, and clinical outcomes relationships.

## Data

All data are sourced from the Genomic Data Commons (GDC) through the
TCGA project:

- Gene (RNA-Seq) expression data
- Somatic mutations
- Clinical information on patient outcomes

## Workflow

1.  Data acquisition (TCGAbiolinks): Retrieving data from TCGA using
    TCGAbiolinks
2.  Mutation analysis (maftools): Identifying significant mutations
    using maftools
3.  Expression analysis (DESeq2): Differential expression analysis with
    DESeq2
4.  Survival analysis: Correlating genomic features with patient
    outcomes
5.  Pathway enrichment: Functional interpretation of genomic alterations
6.  Visualization: Creating publication-quality figures

## Installation

### Prerequisites

- R (version 4.1.0 or higher)
- Bioconductor (version 3.14 or higher)

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c(
  "TCGAbiolinks", "maftools", "DESeq2",
  "survival", "survminer", "pathview",
  "ComplexHeatmap", "BSgenome.Hsapiens.UCSC.hg38"
))
```

## Project Structure

- `01_setup_and_data_acquisition.R`: Setup and data download
- `02_mutation_analysis.R`: TODO
- `03_expression_analysis.R`: TODO
- `04_survival_analysis.R`: TODO
- `05_pathway_enrichment.R`: TODO
- `06_visualization.R`: TODO

## Requirements

- Storage: \>=10GB
- RAM: \>=8GB
- Stable internet connection

## Acknowledgments

- [The Cancer Genome Atlas
  (TCGA)](https://www.cancer.gov/ccg/research/genome-sequencing/tcga)
  for providing the data
- [Bioconductor](https://www.bioconductor.org) for the suite of
  bioinformatics tools
