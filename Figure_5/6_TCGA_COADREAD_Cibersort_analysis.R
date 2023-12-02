#############################################################
# Figure 5
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

# -----------------------------------------------------------------------------
# Extract TCGA sample and generate TPM counts for cibersortx website input
# -----------------------------------------------------------------------------
# extract data
# TPM data generated when analysing cibersortx scores
# calculate TPM
df_tpm_long <- read_tsv(paste0(tcga_dir, "TCGA_COADREAD_MSIH_KMT2Dgroups_tpm_long.txt"))

# Filter for genes in cibersort
ref_cibersortx <- ref %>%
  filter(gene_name %in% cibersortx_ref$"Gene symbol") %>%
  dplyr::select(gene_stable_id, gene_name) %>%
  distinct()

df_tpm_long_cibersortx <- df_tpm_long %>%
  mutate(EnsID_short = str_split_fixed(Ensembl_ID, "\\.", 2)[, 1]) %>%
  filter(EnsID_short %in% ref_cibersortx$gene_stable_id) %>%
  left_join(ref_cibersortx %>% dplyr::select(gene_stable_id, gene_name), by = c("EnsID_short" = "gene_stable_id")) %>%
  dplyr::select(gene_name, Tumor_Sample_Barcode, tpm)

# This matrix was submitted to cibersortx.
df_tpm_wide_cibersortx <- df_tpm_long_cibersortx %>%
  pivot_wider(names_from = Tumor_Sample_Barcode, values_from = tpm) %>%
  write_tsv(paste0(tcga_dir, "Cibersort_input/TCGA_COADREAD_MSIH_KMT2Dgroups_tpm_cibersortx_input.txt"))

# -----------------------------------------------------------------------------
# Analyze Cibersortx results
# -----------------------------------------------------------------------------
Pender_colours <- c("#07316c", "#2172b4", "#6caed5", "#00441b", "#006d2c", "#258b42", "#42ae77", "#67c1a5", "#a0dbca", "#cdece6", "#c61a7e", "#f1b5d9", "#4c014c", "#81107c", "#88419d", "#8c6bb0", "#900000", "#c10504", "#2e3336", "#525252", "#737373", "#bdbdbd")
Cibersortish_colours <- c("#fff385", "#fedb05", "#a98d10", "#fba629", "#e56906", "#e9b96d", "#c27c11", "#a2fa4c", "#6fd904", "#3c7902", "#7dafe2", "#3465a4", "#0d54c3", "#b2a2b4", "#c692c0", "#644469", "#ef3537", "#c10504", "#900000", "#cfd5c9", "#7f7f7d", "#2e3336")

# define samples
samples_to_extract <- sample_combined_status %>%
  filter(
    MSI_status %in% c("MSI-H"),
    KMT2D_groups %in% c("KMT2D_WT", "KMT2D_LOF")
  )

cibersortx_relative <- read_csv(paste0(tcga_dir, "Cibersort_output/TCGA_COADREAD_MSIH_KMT2Dgroups_tpm_cibersortx_output_relative.csv")) %>%
  janitor::clean_names() %>%
  dplyr::rename(Tumor_Sample_Barcode = "mixture")
cibersortx_absolute <- read_csv(paste0(tcga_dir, "Cibersort_output/TCGA_COADREAD_MSIH_KMT2Dgroups_tpm_cibersortx_output_absolute.csv")) %>%
  janitor::clean_names() %>%
  dplyr::rename(Tumor_Sample_Barcode = "mixture")

cibersortx_colour_key <- tibble(
  celltype = colnames(cibersortx_relative)[2:23],
  colour = Pender_colours
)

# Generate plot overview (relative)
cibersortx_relative_annot <- left_join(cibersortx_relative, (samples_to_extract %>%
  select(Tumor_Sample_Barcode, PATIENT_ID, KMT2D_groups)))

plot_relative_df <- cibersortx_relative_annot %>%
  select(Tumor_Sample_Barcode:neutrophils, KMT2D_groups) %>%
  arrange(desc(KMT2D_groups), Tumor_Sample_Barcode) %>%
  pivot_longer(-c(Tumor_Sample_Barcode, KMT2D_groups), names_to = "cell_types", values_to = "fraction") %>%
  distinct()

label_order <- plot_relative_df %>%
  select(Tumor_Sample_Barcode, KMT2D_groups) %>%
  distinct() %>%
  arrange(desc(KMT2D_groups), Tumor_Sample_Barcode)
label_colour <- ifelse(label_order$KMT2D_groups == "KMT2D_WT", "black", "purple")

plot_relative <- plot_relative_df %>%
  mutate(Tumor_Sample_Barcode = fct_relevel(Tumor_Sample_Barcode, label_order$Tumor_Sample_Barcode)) %>%
  ggplot(., aes(x = Tumor_Sample_Barcode, y = fraction, fill = cell_types)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(labels = cibersortx_colour_key$celltype, values = cibersortx_colour_key$colour) +
  theme_bw() +
  theme(
    text = element_text(size = 12),
    axis.text.x = element_text(angle = -90, vjust = 0.5, hjust = 0, colour = label_colour)
  )

pdf(paste0(tcga_dir, "Cibersort_output/TCGA_COADREAD_MSIH_KMT2Dgroups_tpm_cibersortx_output_relative.pdf"), width = 12, height = 8)
print(plot_relative)
dev.off()

# Generate plot overview (absolute)
cibersortx_absolute_annot <- left_join(cibersortx_absolute, (samples_to_extract %>% select(Tumor_Sample_Barcode, PATIENT_ID, KMT2D_groups)))

plot_absolute_df <- cibersortx_absolute_annot %>%
  select(Tumor_Sample_Barcode:neutrophils, KMT2D_groups) %>%
  arrange(desc(KMT2D_groups), Tumor_Sample_Barcode) %>%
  pivot_longer(-c(Tumor_Sample_Barcode, KMT2D_groups), names_to = "cell_types", values_to = "score") %>%
  distinct()

# Compare WT vs LOF for each cell type (absolute)
plot_absolute_comparision <- plot_absolute_df %>%
  mutate(
    KMT2D_groups = fct_relevel(KMT2D_groups, "KMT2D_WT"),
    cell_types = fct_relevel(cell_types, unique(cell_types))
  ) %>%
  ggplot(., aes(x = KMT2D_groups, y = score)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(height = 0, width = 0.2) +
  facet_wrap(~cell_types, scales = "free_y") +
  theme_bw() +
  theme(text = element_text(size = 12))

pdf(paste0(tcga_dir, "Cibersort_output/TCGA_COADREAD_MSIH_KMT2Dgroups_tpm_cibersortx_output_absolute_comparision.pdf"), width = 12, height = 10)
print(plot_absolute_comparision)
dev.off()