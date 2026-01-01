library(tidyverse)
library(SummarizedExperiment)

data_rna <- readRDS("data/raw/tcga_gbm_rnaseq_raw.rds")
data_met <- readRDS("data/raw/tcga_gbm_methyl_raw.rds")

met_matrix <- assay(data_met)

# select only the 2 mgmt probes and check intersection

target_probes <- c("cg12434587", "cg12981137")

mgmt_values <- met_matrix[target_probes, ] |>
    t() |>
    as.data.frame() |>
    rownames_to_column(var = "barcode")

# clean ids and create labels

mgmt_labels <- mgmt_values |>
    mutate(patient_id = substr(barcode, 1, 12)) |>
    rowwise() |>
    mutate(avg_beta = mean(c(cg12434587, cg12981137), na.rm = TRUE)) |>
    ungroup() |>
    mutate(MGMT_status = if_else(avg_beta > 0.21, 1, 0)) |>
    select(patient_id, MGMT_status) |>
    distinct(patient_id, .keep_all = TRUE)

print(paste("labels created. count:", nrow(mgmt_labels)))
table(mgmt_labels$MGMT_status)


# process gene expression

rna_matrix <- assay(data_rna)

# filter for top 2000 variable genes

print("filtering for top 2000 most variable genes...")
gene_variances <- apply(rna_matrix, 1, var)
top_genes <- names(sort(gene_variances, decreasing = TRUE))[1:2000]
rna_matrix_filtered <- rna_matrix[top_genes, ]

# transpose and clean ids

rna_features <- rna_matrix_filtered |>
    t() |>
    as.data.frame() |>
    rownames_to_column(var = "barcode") |>
    mutate(patient_id = substr(barcode, 1, 12)) |>
    select(patient_id, everything(), -barcode) |>
    distinct(patient_id, .keep_all = TRUE)

print(paste("gene features processed. patients:", nrow(rna_features)))


# merge datasets

final_dataset <- rna_features |>
    inner_join(mgmt_labels, by = "patient_id")

print("--- merge complete ---")
print(paste("final patient count:", nrow(final_dataset)))
print(paste("final feature count:", ncol(final_dataset) - 2))

# save final dataset

dir.create("data/processed", showWarnings = FALSE)
write_csv(final_dataset, "data/processed/gbm_final_dataset.csv")
print("saved to data/processed/gbm_final_dataset.csv")
