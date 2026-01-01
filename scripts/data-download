library(TCGAbiolinks)
library(SummarizedExperiment)
library(tidyverse)
library(sesameData)

# download gene expression (rna-seq)

query_rna <- GDCquery(
    project = "TCGA-GBM",
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
)

GDCdownload(query_rna, directory = "data/raw")

data_rna_se <- GDCprepare(query_rna, directory = "data/raw")

# download methylation (HM27) data

query_met <- GDCquery(
    project = "TCGA-GBM",
    data.category = "DNA Methylation",
    platform = "Illumina Human Methylation 27",
    data.type = "Methylation Beta Value"
)

GDCdownload(query_met, directory = "data/raw")

data_met_se <- GDCprepare(query_met, directory = "data/raw")

# save raw data

saveRDS(data_rna_se, file = "data/raw/tcga_gbm_rnaseq_raw.rds")
saveRDS(data_met_se, file = "data/raw/tcga_gbm_methyl_raw.rds")
