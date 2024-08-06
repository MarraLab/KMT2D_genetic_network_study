#############################################################
# Figure 4
# TCGA COAD/READ cohort - Cibersort immune cell composition analysis
#############################################################
library(pacman)
p_load(tidyverse, survival, lubridate, ggsurvfit, gtsummary, tidycmprsk, condsurv)

# -----------------------------------------------------------------------------
# Path to pan-cancer data
# -----------------------------------------------------------------------------
# Data were downloaded from cBioPortal
tcga_dir <- "cBioPortal/"

# -----------------------------------------------------------------------------
# Load TCGA data
# -----------------------------------------------------------------------------
data_clinical_patient <- read_tsv(paste0(data_key$dir_full[i],"/data_clinical_patient.txt"), skip = 4)
# Created in Figure_5/2_TCGA_COADREAD_KMT2D_MSI_status.R
sample_combined_status <- read_csv(paste0(tcga_dir, "TCGA_COADREAD_sample_status.csv"))
# Downloaded from BioMart
ref <- read_csv("BioMart/All_genes_EnsID_Hugo_UCSC_NCBI_HGNC_ids.csv") %>% janitor::clean_names()
# Cibersort genes:
cibersortx_ref <- read_tsv("Cibersortx_reference_LM22.txt")
# PredictIO 100 gene list received from from Minoru Nakano, BHK lab
load("PredictIO_gene.RData")
# Created in Figure_5/3_TCGA_COADREAD_ICI_response_markers.R
df_tpm_long <- read_tsv(paste0(tcga_dir, "TCGA_COADREAD_KMT2Dgroups_tpm_long.txt"))

# -----------------------------------------------------------------------------
# Extract TCGA sample and generate TPM counts for PredictIO website input
# -----------------------------------------------------------------------------
# Filter for genes in predictIO
ref_predictIO <- ref %>%
  filter(gene_name %in% predictIO_gene$gene) %>%
  dplyr::select(gene_stable_id, gene_name) %>%
  distinct()

df_tpm_long_predictIO <- df_tpm_long %>%
  mutate(EnsID_short = str_split_fixed(Ensembl_ID, "\\.", 2)[, 1]) %>%
  filter(EnsID_short %in% ref_predictIO$gene_stable_id) %>%
  left_join(ref_predictIO %>% dplyr::select(gene_stable_id, gene_name), by = c("EnsID_short" = "gene_stable_id")) %>%
  dplyr::select(gene_name, Tumor_Sample_Barcode, tpm)

# This matrix was submitted to PredictIO.
df_tpm_wide_predictIO <- df_tpm_long_predictIO %>%
  pivot_wider(names_from = Tumor_Sample_Barcode, values_from = tpm) %>%
  write_tsv(paste0(tcga_dir, "TCGA_COADREAD_MSIH_KMT2Dgroups_tpm_predictIO_input.txt"))

# -----------------------------------------------------------------------------
# PredictIO output analysis
# -----------------------------------------------------------------------------
# load results generated from PredictIO's web app & plot
predictIO_res <- read_csv(paste0(tcga_dir, "TCGA_COADREAD_MSIH_KMT2Dgroups_tpm_predictIO_results.csv")) %>%
  mutate(Tumor_Sample_Barcode = str_replace_all(patient_id, "\\.", "-")) %>%
  left_join(sample_combined_status) %>%
  dplyr::select(Tumor_Sample_Barcode, predictio_value, KMT2D_groups) %>%
  distinct() %>%
  mutate(
    KMT2D_groups = fct_relevel(KMT2D_groups, "KMT2D_WT"),
    TCGA_project = project
  )
predictIO_values_all <- bind_rows(predictIO_values_all, predictIO_res)

plot <- ggplot(predictIO_res, aes(x = KMT2D_groups, y = predictio_value)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(height = 0, width = 0.2) +
  theme_bw() +
  theme(text = element_text(size = 12))

pdf(paste0(tcga_dir, "TCGA_COADREAD_MSIH_KMT2Dgroups_tpm_predictIO_results.pdf"), width = 4, height = 6)
print(plot)
dev.off()

# calculate difference between groups
fit <- wilcox.test(predictio_value ~ KMT2D_groups, predictIO_res)