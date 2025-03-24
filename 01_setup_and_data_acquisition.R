if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c(
  "TCGAbiolinks", "maftools", "DESeq2",
  "survival", "survminer", "pathview",
  "ComplexHeatmap", "BSgenome.Hsapiens.UCSC.hg38"
))

library(TCGAbiolinks)

query <- GDCquery(
  project = "TCGA-LAML",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

GDCdownload(query)

data <- GDCprepare(query)

saveRDS(data, "results/laml_expression.rds")
