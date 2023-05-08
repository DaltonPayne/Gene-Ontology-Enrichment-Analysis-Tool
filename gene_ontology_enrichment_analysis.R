#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  cat("Usage: Rscript gene_ontology_enrichment_analysis.R <expression_file>\n")
  quit("no", 1)
}

expression_file <- args[1]

library(clusterProfiler)
library(DESeq2)
library(org.Hs.eg.db)

## Load expression data
countData <- read.table(expression_file, header=TRUE, row.names=1)
condition <- factor(rep(c("treated", "control"), each=3))
colData <- data.frame(condition, row.names=colnames(countData))

## DESeq2 analysis
dds <- DESeqDataSetFromMatrix(countData, colData, design=~condition)
dds <- DESeq(dds)

## Get differentially expressed genes
res <- results(dds, alpha=0.05, lfcThreshold=1)
degs <- rownames(res)[which(res$padj < 0.05)]

## Gene Ontology enrichment analysis
ego <- enrichGO(gene     = degs,
               OrgDb    = org.Hs.eg.db,
               keyType  = "ENSEMBL",
               ont      = "BP",
               pAdjustMethod = "BH",
               pvalueCutoff  = 0.05,
               qvalueCutoff  = 0.1)

## Print enriched GO terms
print(ego)
