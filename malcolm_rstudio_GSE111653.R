# Install, update & load packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.18")

library(BiocManager)

BiocManager::install("biomaRt", update = TRUE, force = TRUE)
BiocManager::install("cowplot", update = TRUE)
BiocManager::install("ggplot2", update = TRUE)
BiocManager::install("tximport", update = TRUE)
BiocManager::install("rhdf5", update = TRUE)
BiocManager::install("gridExtra", update = TRUE)
BiocManager::install("devtools", update = TRUE)

library(devtools)
library(rhdf5)

devtools::install_github("pachterlab/sleuth", force = TRUE)

library(biomaRt)
library(ggplot2)
library(tximport)
library(dplyr)
library(gridExtra)

suppressMessages({
  library('cowplot')
  library('sleuth')})

# build ensembl mart
listEnsembl()
listEnsemblArchives()
ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", version = 111)
t2gs <- getBM(attributes = c('ensembl_gene_id', 'transcript_version', 'ensembl_transcript_id', 'ensembl_gene_id', 'external_gene_name', 'description', 'transcript_biotype', 'hgnc_symbol'), mart = ensembl)
t2gs <- dplyr::rename(t2gs, target_id = ensembl_transcript_id,
                      ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
t2gs <- dplyr::select(t2gs, c('target_id', 'ens_gene', 'ext_gene'))
head(t2gs)

# prepare sleuth object metadata
metadata <- read.table("experiment_info_hk.txt", sep = '\t', header = TRUE, stringsAsFactors = FALSE)
head(metadata)
metadata <- dplyr::select(metadata, c('sample', 'cell_line', 'env_o2'))
head(metadata)
metadata <- dplyr::mutate(metadata, 
                          path = file.path('kallisto_quant_v2', sample, "abundance.h5"))
head(metadata)

## initializing sleuth object 
duplicated_targets <- t2gs[duplicated(t2gs$target_id), "target_id"]
t2gs <- t2gs[!duplicated(t2gs$target_id), ]

so <- sleuth_prep(metadata, target_mapping = t2gs,
                  aggregation_column = "ens_gene", extra_bootstrap_summary = TRUE)
head(so$target_mapping)
head(so$obs_raw)
head(so$obs_norm)
head(so$obs_norm_filt)
normalised_counts <- so$obs_norm
est_counts_no_decimal <- normalised_counts$target_id <- sub("\\.\\d+$", "", normalised_counts$target_id)
head(normalised_counts)
normalised_counts <- merge(normalised_counts, t2gs, by.x = "target_id", by.y = "target_id", all.x = TRUE)
head(normalised_counts)

normalised_counts_filtered <- so$obs_norm_filt
filt_counts_no_decimal <- normalised_counts_filtered$target_id <- sub("\\.\\d+$", "", normalised_counts_filtered$target_id)
head(normalised_counts_filtered)
normalised_counts_filtered <- merge(normalised_counts_filtered, t2gs, by.x = "target_id", by.y = "target_id", all.x = TRUE)
head(normalised_counts)

so <- sleuth_prep(metadata, target_mapping = t2gs, aggregation_column = "ext_gene", gene_mode = TRUE, extra_bootstrap_summary = TRUE)

aggregated_data <- so$obs_norm_filt

aggregated_data <- aggregated_data %>%
  mutate(cell_type = case_when(
    sample == 831 ~ "MCF_7",
    sample == 832 ~ "MCF_7",
    sample == 837 ~ "MDA_MB_231",
    sample == 838 ~ "MDA_MB_231",
    sample == 841 ~ "MDA_MB_468",
    sample == 842 ~ "MDA_MB_468",
    sample == 857 ~ "T_47D",
    sample == 858 ~ "T_47D",
    
  ),
  env_o2 = case_when(
    sample == 831 ~ "hypoxia",
    sample == 832 ~ "normoxia",
    sample == 837 ~ "hypoxia",
    sample == 838 ~ "normoxia",
    sample == 841 ~ "hypoxia",
    sample == 842 ~ "normoxia",
    sample == 857 ~ "hypoxia",
    sample == 858 ~ "normoxia",
  ))

## plot sample heatmap to show distances between samples
plot_sample_heatmap(so, use_filtered = TRUE, color_high = "dodgerblue",
                    color_low = "white", x_axis_angle = 50,
                    annotation_cols = setdiff(colnames(so$sample_to_covariates), "sample"),
                    cluster_bool = TRUE)

plot_pca(so, pc_x = 1L, pc_y = 2L, use_filtered = TRUE,
         units = "scaled_reads_per_base", text_labels = TRUE, color_by = "sample",
         point_size = 3, point_alpha = 0.8) + 
  theme(panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))

write.csv(aggregated_data, "G:\\My Drive\\RStudio directories\\RStudio_Housekeeper\\kallisto_sleuth_aggregated_counts.csv", row.names = FALSE)

aggregated_data_unfiltered <- so$obs_norm

aggregated_data_unfiltered <- aggregated_data_unfiltered %>%
  mutate(cell_type = case_when(
    sample == 831 ~ "MCF_7",
    sample == 832 ~ "MCF_7",
    sample == 837 ~ "MDA_MB_231",
    sample == 838 ~ "MDA_MB_231",
    sample == 841 ~ "MDA_MB_468",
    sample == 842 ~ "MDA_MB_468",
    sample == 857 ~ "T_47D",
    sample == 858 ~ "T_47D",
    
  ),
  env_o2 = case_when(
    sample == 831 ~ "hypoxia",
    sample == 832 ~ "normoxia",
    sample == 837 ~ "hypoxia",
    sample == 838 ~ "normoxia",
    sample == 841 ~ "hypoxia",
    sample == 842 ~ "normoxia",
    sample == 857 ~ "hypoxia",
    sample == 858 ~ "normoxia",
  ))


write.csv(aggregated_data_unfiltered, "G:\\My Drive\\RStudio directories\\RStudio_Housekeeper\\full_kallisto_sleuth_aggregated_counts.csv", row.names = FALSE)
