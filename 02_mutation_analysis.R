library(TCGAbiolinks)
library(maftools)
library(SummarizedExperiment)
library(BSgenome.Hsapiens.UCSC.hg38)

## Create output directories if they do not exist
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

## Get mutation data
query_mutations <- GDCquery(
  project = "TCGA-LAML",
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking",
)
GDCdownload(query_mutations)
laml_maf <- GDCprepare(query_mutations)
saveRDS(laml_maf, "results/laml_mutations_raw.rds")

## Crate MAF (mutation annotation format) objects
laml <- read.maf(maf = laml_maf)
laml_summary <- getSampleSummary(laml)
laml_gene_summary <- getGeneSummary(laml)

## Generate summaries from MAF objects
write.csv(laml_summary, "results/laml_summary.csv", row.names = FALSE)
write.csv(laml_gene_summary, "results/laml_gene_summary.csv", row.names = FALSE)

## Create oncoplot for mutation landscape visualization
pdf("figures/laml_oncoplot_top20.pdf", width = 10, height = 8)
oncoplot(maf = laml, top = 20, removeNonMutated = TRUE)
dev.off()

## Transition and transversion analysis
pdf("figures/laml_titv.pdf", width = 10, height = 8)
titv <- titv(maf = laml, plot = TRUE)
dev.off()

## Lollipop plots for key AML genes
## FLT3 is commonly mutated in AML
pdf("figures/laml_lollipop_FLT3.pdf", width = 10, height = 6)
lollipopPlot(maf = laml, gene = "FLT3", AACol = "HGVSp_Short", showMutationRate = TRUE)
dev.off()

## DNMT3A is also commonly mutated in AML
pdf("figures/laml_lollipop_DNMT3A.pdf", width = 10, height = 6)
lollipopPlot(maf = laml, gene = "DNMT3A", AACol = "HGVSp_Short", showMutationRate = TRUE)
dev.off()

## Identify significantly mutated genes
laml_sig_genes <- oncodrive(maf = laml, AACol = "HGVSp_Short", minMut = 5)
saveRDS(laml_sig_genes, "results/laml_significant_genes.rds")

## Analyze mutation co-occurrence and mutual exclusivity
pdf("figures/laml_cooccurrence.pdf", width = 10, height = 8)
somaticInteractions(maf = laml, top = 25, pvalue = c(0.05, 0.1))
dev.off()

## Compare mutation profiles with clinical data
if (file.exists("results/laml_expression.rds")) {
  # Load previously saved expression data
  expr_data <- readRDS("results/laml_expression.rds")

  # Extract clinical data
  clinical_data <- as.data.frame(colData(expr_data))

  # Identify and fix list columns
  list_cols <- sapply(clinical_data, is.list)
  if (any(list_cols)) {
    for (col in names(which(list_cols))) {
      clinical_data[[col]] <- sapply(
        clinical_data[[col]],
        function(x) paste(unlist(x), collapse = ";")
      )
    }
  }

  write.csv(clinical_data, "results/laml_clinical_data.csv", row.names = FALSE)

  # Compare survival by key gene mutations using
  # mafSurvival function for Kaplan-Meier analysis
  if ("days_to_last_followup" %in% colnames(clinical_data) &
    "vital_status" %in% colnames(clinical_data)) {
    pdf("figures/laml_FLT3_survival.pdf", width = 8, height = 6)
    mafSurvival(
      maf = laml, genes = "FLT3", time = "days_to_last_followup",
      Status = "vital_status", isTCGA = TRUE
    )
    dev.off()
  }
}

## Extract mutation signatures
laml_tnm <- trinucleotideMatrix(
  maf = laml,
  prefix = NULL,
  add = FALSE,
  ref_genome = "BSgenome.Hsapiens.UCSC.hg38"
)

## There is a column with no mutations in the T[T>G]A context
ttga_col <- which(colnames(laml_tnm$nmf_matrix) == "T[T>G]A")
if (length(ttga_col) > 0) {
  # Add a very small pseudo-count just to that column
  laml_tnm$nmf_matrix[, ttga_col] <- 0.001
  print("Added pseudocount to T[T>G]A context")
}

## Extract signatures
laml_signatures <- extractSignatures(
  mat = laml_tnm,
  n = 1,
  plotBestFitRes = FALSE
)

## Plot the signatures
pdf("figures/laml_mutation_signatures.pdf", width = 10, height = 8)
plotSignatures(laml_signatures)
dev.off()

cat("Mutation analysis completed. Results saved in 'results' and 'figures' directories.\n")
