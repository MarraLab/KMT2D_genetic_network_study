#############################################################
# Figure 4
# TCGA COAD/READ cohort - Annotating KMT2D WT and LOF mutations and MSI status
# TMB, Neoantigen, HRD score, TUBA1B, CD274, CTLA4, M1M2 score,
# Cytolytic activity score, INF-y signature score
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
# TPM were calculated from fpkm counts downloaded from UCSC XenaHub: (because TCGA/GDC people are not organized enough anddata has be be scattered across multiple datasets... I'm not annoyed you are :P )
# https://xenabrowser.net/datapages/?dataset=TCGA-COAD.htseq_fpkm.tsv&host=https%3A%2F%2Fgdc.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
# https://xenabrowser.net/datapages/?cohort=GDC%20TCGA%20Rectal%20Cancer%20(READ)&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
# COAD / READ needs to be combined
# read data
df_coad <- read_tsv(paste0(tcga_dir, "Datasets/TCGA-COAD.htseq_fpkm.tsv"))
df_read <- read_tsv(paste0(tcga_dir, "Datasets/TCGA-READ.htseq_fpkm.tsv"))

# -----------------------------------------------------------------------------
# Prepared RNAseq data
# -----------------------------------------------------------------------------
# define samples
samples_to_extract <- sample_combined_status %>%
  filter(
    TCGA_project == project,
    MSI_status %in% c("MSI-H", "MSS/MSI-L"),
    KMT2D_groups %in% c("KMT2D_WT", "KMT2D_LOF")
  )

# filter and combine two projects
df_coad_long <- df_coad %>%
  pivot_longer(-Ensembl_ID, names_to = "Tumor_Sample_Barcode", values_to = "log2fpkm")
df_read_long <- df_read %>%
  pivot_longer(-Ensembl_ID, names_to = "Tumor_Sample_Barcode", values_to = "log2fpkm")

df_coad_long_subset <- df_coad_long %>%
  mutate(Tumor_Sample_Barcode = str_sub(Tumor_Sample_Barcode, end = -2)) %>%
  filter(Tumor_Sample_Barcode %in% samples_to_extract$Tumor_Sample_Barcode)
df_read_long_subset <- df_read_long %>%
  mutate(Tumor_Sample_Barcode = str_sub(Tumor_Sample_Barcode, end = -2)) %>%
  filter(Tumor_Sample_Barcode %in% samples_to_extract$Tumor_Sample_Barcode)
df_long_subset <- bind_rows(df_coad_long_subset, df_read_long_subset)

# check and remove duplicates - this step will take a while
print("Checking for duplicates.... this can take a while")
test <- df_long_subset %>%
  distinct() %>%
  dplyr::count(Ensembl_ID, Tumor_Sample_Barcode) %>%
  arrange(-n)

if (any(test$n > 1)) {
  multiple <- test %>%
    filter(n > 1) %>%
    select(-n) %>%
    distinct()
  clean <- df_long_subset %>%
    filter(Tumor_Sample_Barcode %in% multiple$Tumor_Sample_Barcode) %>%
    distinct() %>%
    group_by(Ensembl_ID, Tumor_Sample_Barcode) %>%
    filter(row_number() == 1) %>%
    ungroup()

  df_long_subset_clean <- df_long_subset %>%
    filter(!Tumor_Sample_Barcode %in% multiple$Tumor_Sample_Barcode) %>%
    bind_rows(clean)
} else {
  df_long_subset_clean <- df_long_subset
}

# create wide fpkm
df_long_subset_clean_wide <- df_long_subset_clean %>%
  mutate(fpkm = (2^log2fpkm) - 1) %>%
  dplyr::select(-log2fpkm) %>%
  pivot_wider(names_from = Tumor_Sample_Barcode, values_from = fpkm) %>%
  write_tsv(paste0(tcga_dir, "TCGA_COADREAD_KMT2Dgroups_fpkm_wide.txt"))

# calculate TPM
df_tpm_wide <- df_long_subset_clean_wide %>%
  mutate_if(is.numeric, ~ (. / sum(.)) * 10**6) %>%
  mutate_if(is.numeric, ~ round((. / sum(.)) * 10**6, 2)) %>%
  write_tsv(paste0(tcga_dir, "TCGA_COADREAD_KMT2Dgroups_tpm_wide.txt"))

df_tpm_long <- df_tpm_wide %>%
  pivot_longer(-Ensembl_ID, names_to = "Tumor_Sample_Barcode", values_to = "tpm") %>%
  write_tsv(paste0(tcga_dir, "TCGA_COADREAD_KMT2Dgroups_tpm_long.txt"))

# -----------------------------------------------------------------------------
# TMB
# -----------------------------------------------------------------------------
TMB_df <- sample_combined_status %>%
  dplyr::select(-Matched_Norm_Sample_Barcode) %>%
  left_join(data_clinical_sample %>% dplyr::select(Tumor_Sample_Barcode = SAMPLE_ID, TMB_NONSYNONYMOUS)) %>%
  distinct()

# TMB Plot
TMB_MSI_plot <- TMB_df %>%
  filter(
    MSI_status %in% c("MSI-H", "MSS/MSI-L"), 
    KMT2D_groups %in% c("KMT2D_WT", "KMT2D_LOF")) %>%
  mutate(KMT2D_groups = fct_relevel(KMT2D_groups, "KMT2D_WT")) %>%
  ggplot(., aes(y = TMB_NONSYNONYMOUS, x = KMT2D_groups)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(height = 0, width = 0.2) +
  facet_grid(~MSI_status) +
  theme_bw()

pdf(paste0(tcga_dir, "Boxplot_TCGA_COADREAD_TMB_by_groups.pdf"), width = 6, height = 6)
print(TMB_MSI_plot)
dev.off()

# MSI-H ------------------------------------------------------
# T-test
TMB_df %>%
  filter(
    MSI_status == "MSI-H",
    KMT2D_groups %in% c("KMT2D_WT", "KMT2D_LOF")
  ) %>%
  wilcox.test(TMB_NONSYNONYMOUS ~ KMT2D_groups, data = .)

# MSS ------------------------------------------------------
# T-test
TMB_df %>%
  filter(
    MSI_status == "MSS/MSI-L",
    KMT2D_groups %in% c("KMT2D_WT", "KMT2D_LOF")
  ) %>%
  wilcox.test(TMB_NONSYNONYMOUS ~ KMT2D_groups, data = .)

# -----------------------------------------------------------------------------
# Neoantigen score
# -----------------------------------------------------------------------------
# Downloaded from Thorsson et al. 
thorsson_df <- readxl::read_xlsx("/projects/marralab/ytakemon_prj/KMT2D_project/data/Thorsson_etal_Table_S1.xlsx") %>%
  janitor::clean_names() %>%
  dplyr::rename(PATIENT_ID = "tcga_participant_barcode") %>%
  mutate(
    HRD_score = as.numeric(homologous_recombination_defects),
    snv_neoantigens = as.numeric(snv_neoantigens)
  ) %>%
  dplyr::select(
    PATIENT_ID, tcga_study,
    HRD_score,
    snv_neoantigens
  ) %>%
  distinct()
  project <- "TCGA-COADREAD"

# read RNA-seq FPKM data
# TRPKM counts downloaded from UCSC XenaHub: 
# https://xenabrowser.net/datapages/?dataset=TCGA-COAD.htseq_fpkm.tsv&host=https%3A%2F%2Fgdc.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
# https://xenabrowser.net/datapages/?cohort=GDC%20TCGA%20Rectal%20Cancer%20(READ)&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
df_coad <- read_tsv(paste0(tcga_dir, "Datasets/TCGA-COAD.htseq_fpkm.tsv"))
df_read <- read_tsv(paste0(tcga_dir, "Datasets/TCGA-READ.htseq_fpkm.tsv"))

# define samples
samples_to_extract <- sample_combined_status %>%
  filter(
    TCGA_project == project,
    MSI_status %in% c("MSI-H", "MSS/MSI-L"),
    KMT2D_groups %in% c("KMT2D_WT", "KMT2D_LOF")
  )

thorsson_df_subset <- thorsson_df %>%
  filter(PATIENT_ID %in% samples_to_extract$PATIENT_ID) %>%
  left_join(., samples_to_extract) %>%
  mutate(KMT2D_groups = fct_relevel(KMT2D_groups, "KMT2D_WT")) %>%
  dplyr::select(
    PATIENT_ID, MSI_status, KMT2D_groups, everything()
  ) %>%
  distinct() %>%
  group_by(PATIENT_ID) %>%
  arrange(Tumor_Sample_Barcode) %>%
  slice(1) %>%
  ungroup()

# Plot neoantigen scores
plot_neoantigen <- ggplot(thorsson_df_subset, aes(x = KMT2D_groups, y = as.numeric(snv_neoantigens))) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(height = 0, width = 0.2) +
  facet_grid(~MSI_status)+  
  theme_bw() +
  theme(text = element_text(size = 12))

pdf(paste0(tcga_dir, "TCGA_COADREAD_MSIH_KMT2Dgroups_neoantigen_results.pdf"), width = 4, height = 6)
print(plot_neoantigen)
dev.off()

# MSI -------------------------------------------------
thorsson_df_subset %>%
  filter(MSI_status == "MSI-H") %>%
  wilcox.test(snv_neoantigens ~ KMT2D_groups, .)

# MSS -------------------------------------------------
thorsson_df_subset %>%
  filter(MSI_status == "MSS/MSI-L") %>%
  wilcox.test(snv_neoantigens ~ KMT2D_groups, .)

# -----------------------------------------------------------------------------
# HRD score
# -----------------------------------------------------------------------------
# Compare HRD scores
plot_hrd <- ggplot(thorsson_df_subset, aes(x = KMT2D_groups, y = HRD_score)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(height = 0, width = 0.2) +
  facet_grid(~MSI_status) +
  theme_bw() +
  theme(text = element_text(size = 12))

pdf(paste0(tcga_dir, "TCGA_COADREAD_MSIH_KMT2Dgroups_HRD_results.pdf"), width = 6, height = 6)
print(plot_hrd)
dev.off()

# MSI -------------------------------------------------
thorsson_df_subset %>%
  filter(MSI_status == "MSI-H") %>%
  wilcox.test(HRD_score ~ KMT2D_groups, .)

# MSS -------------------------------------------------
thorsson_df_subset %>%
  filter(MSI_status == "MSS/MSI-L") %>%
  wilcox.test(HRD_score ~ KMT2D_groups, .)

# -----------------------------------------------------------------------------
# Calculate additional markers: TUBA1B
# -----------------------------------------------------------------------------
project <- data_key_MSIH$TCGA_project[i]
  print(project)

# define samples
samples_to_extract <- sample_combined_status %>%
  filter(
    MSI_status %in% c("MSI-H", "MSS/MSI-L"),
    KMT2D_groups %in% c("KMT2D_WT", "KMT2D_LOF")) 

# Get fpkm 
fpkm_long <- read_tsv(paste0(tcga_dir, "TCGA_COADREAD_KMT2Dgroups_fpkm_wide.txt")) %>%
  pivot_longer(-Ensembl_ID, names_to = "Tumor_Sample_Barcode", values_to = "fpkm") %>%
  mutate(EnsID_short = str_split_fixed(Ensembl_ID,"\\.",2)[,1]) %>%
  left_join(ref %>% dplyr::select(gene_stable_id, gene_name), by = c("EnsID_short" = "gene_stable_id")) %>%
  left_join(sample_combined_status %>% select(Tumor_Sample_Barcode, KMT2D_groups)) %>% distinct %>%
  filter(Tumor_Sample_Barcode %in% samples_to_extract$Tumor_Sample_Barcode)

# TUBA1B expression comparison plot
TUBA1B_df <- fpkm_long %>%
  filter(gene_name == "TUBA1B", KMT2D_groups %in% c("KMT2D_WT", "KMT2D_LOF"))

TUBA1B_plot <- TUBA1B_df %>%
  mutate(KMT2D_groups = fct_relevel(KMT2D_groups, "KMT2D_WT")) %>%
  ggplot(., aes(x = KMT2D_groups, y = fpkm)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(height = 0, width = 0.2) +
  ylab("TUBA1B FPKM") +
  facet_grid(~MSI_status)+
  theme_bw()+
  theme(text = element_text(size = 12))

pdf(paste0(tcga_dir, "TCGA_COADREAD_KMT2Dgroups_additional_marker_comparision_TUBA1B.pdf"), width = 4, height = 5)
print(TUBA1B_plot)
dev.off()

# MSI p-value
TUBA1B_df %>%
  filter(MSI_status == "MSI-H") %>%
  wilcox.test(fpkm ~ KMT2D_groups, .)

# MSS p-value
TUBA1B_df %>%
  filter(MSI_status == "MSS/MSI-L") %>%
  wilcox.test(fpkm ~ KMT2D_groups, .)

# -----------------------------------------------------------------------------
# Calculate additional markers: CD274
# -----------------------------------------------------------------------------
# CD274 expression comparison (PD-L1) - (FYI - Cibersort matrix doesn't use CD274 expression)
cd274_df <- fpkm_long %>%
  filter(gene_name == "CD274", KMT2D_groups %in% c("KMT2D_WT", "KMT2D_LOF"))

cd274_plot <- cd274_df %>%
  mutate(KMT2D_groups = fct_relevel(KMT2D_groups, "KMT2D_WT")) %>%
  ggplot(., aes(x = KMT2D_groups, y = fpkm)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(height = 0, width = 0.2) +
  facet_grid(~MSI_status)+
  ylab("CD274 FPKM")+
  theme_bw()+
  theme(text = element_text(size = 12))

pdf(paste0(tcga_dir, "TCGA_COADREAD_KMT2Dgroups_additional_marker_comparision_cd274.pdf"), width = 4, height = 5)
print(cd274_plot)
dev.off()

# MSI p-value
cd274_df %>%
  filter(MSI_status == "MSI-H") %>%
  wilcox.test(fpkm ~ KMT2D_groups, .)

# MSS p-value
cd274_df %>%
  filter(MSI_status == "MSS/MSI-L") %>%
  wilcox.test(fpkm ~ KMT2D_groups, .)

# -----------------------------------------------------------------------------
# Calculate additional markers: CTLA4A
# -----------------------------------------------------------------------------
CTLA4_df <- fpkm_long %>%
  filter(gene_name == "CTLA4", KMT2D_groups %in% c("KMT2D_WT", "KMT2D_LOF"))

CTLA4_plot <- CTLA4_df %>%
  mutate(KMT2D_groups = fct_relevel(KMT2D_groups, "KMT2D_WT")) %>%
  ggplot(., aes(x = KMT2D_groups, y = fpkm)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(height = 0, width = 0.2) +
  facet_grid(~MSI_status) +
  ylab("CTLA4 FPKM")+
  theme_bw()+
  theme(text = element_text(size = 12))

pdf(paste0(tcga_dir, "TCGA_COADREAD_KMT2Dgroups_additional_marker_comparision_CTLA4.pdf"), width = 4, height = 5)
print(CTLA4_plot)
dev.off()

# MSI p-value
CTLA4_df %>%
  filter(MSI_status == "MSI-H") %>%
  wilcox.test(fpkm ~ KMT2D_groups, .)

# MSS p-value
CTLA4_df %>%
  filter(MSI_status == "MSS/MSI-L") %>%
  wilcox.test(fpkm ~ KMT2D_groups, .)

fit <- wilcox.test(fpkm ~ KMT2D_groups, CTLA4_df)
n_WT <- CTLA4_df %>% count(KMT2D_groups) %>% 
  filter(str_detect(KMT2D_groups, "KMT2D_WT")) %>% pull(n)
n_LOF <- CTLA4_df %>% count(KMT2D_groups) %>% 
  filter(str_detect(KMT2D_groups, "KMT2D_LOF")) %>% pull(n)
res_CTLA4 <- tibble(project = project, n_WT = n_WT, n_LOF = n_LOF, type = "CTLA4_fpkm", pval = fit$p.value)

# -----------------------------------------------------------------------------
# Calculate additional markers: M1M2 score
# -----------------------------------------------------------------------------
# M1-M2 score (Pender et al paper)
# "The M1-M2 macrophage score for each sample was derived by calculating the mean expression (RPKM) of these 10 genes (CXCL11, IDO1, CCL19, CXCL9, PLA1A, LAMP3, CCR7, APOL6, CXCL10, TNIP3)." 
m1m2 <- c("CXCL11", "IDO1", "CCL19", "CXCL9", "PLA1A", "LAMP3", "CCR7", "APOL6", "CXCL10", "TNIP3")

m1m2_df <- fpkm_long %>%
  filter(gene_name %in% m1m2, KMT2D_groups %in% c("KMT2D_WT", "KMT2D_LOF")) %>%
  distinct() %>%
  group_by(Tumor_Sample_Barcode) %>%
  summarise(m1m2_score = mean(fpkm)) %>%
  ungroup() %>%
  left_join(sample_combined_status %>% select(Tumor_Sample_Barcode, KMT2D_groups)) %>%
  distinct()

m1m2_plot <- m1m2_df %>%
  mutate(KMT2D_groups = fct_relevel(KMT2D_groups, "KMT2D_WT")) %>%
  ggplot(., aes(x = KMT2D_groups, y = m1m2_score)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(height = 0, width = 0.2) +
  facet_grid(~MSI_status) +
  theme_bw() +
  theme(text = element_text(size = 12))

pdf(paste0(tcga_dir, "TCGA_COADREAD_KMT2Dgroups_additional_marker_comparision_m1m2score.pdf"), width = 4, height = 5)
print(m1m2_plot)
dev.off()

# MSI p-value
m1m2_df %>%
  filter(MSI_status == "MSI-H") %>%
  wilcox.test(m1m2_score ~ KMT2D_groups, .)

# MSS p-value
m1m2_df %>%
  filter(MSI_status == "MSS/MSI-L") %>%
  wilcox.test(m1m2_score ~ KMT2D_groups, .)

# -----------------------------------------------------------------------------
# Calculate additional markers: Cytolytic activity score
# -----------------------------------------------------------------------------
# CYT - cytolytic activity score (Pender et al paper via Rooney et al)
# Cytolytic activity (CYT) was calculated as the geometric mean of GZMA and PRF1 (as expressed in TPM, 0.01 offset). 
cyt <- c("GZMA", "PRF1")

cyt_df <- fpkm_long %>%
  filter(gene_name %in% cyt, KMT2D_groups %in% c("KMT2D_WT", "KMT2D_LOF")) %>%
  distinct() %>%
  group_by(Tumor_Sample_Barcode) %>%
  summarise(cyt_score = exp(mean(log(fpkm)))) %>%
  ungroup() %>%
  left_join(sample_combined_status %>% select(Tumor_Sample_Barcode, KMT2D_groups)) %>%
  distinct()

cyt_plot <- cyt_df %>%
  mutate(KMT2D_groups = fct_relevel(KMT2D_groups, "KMT2D_WT")) %>%
  ggplot(., aes(x = KMT2D_groups, y = cyt_score)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(height = 0, width = 0.2) +
  facet_grid(~MSI_status) +
  theme_bw() +
  theme(text = element_text(size = 12))

pdf(paste0(tcga_dir, "TCGA_COADREAD_KMT2Dgroups_additional_marker_comparision_cytscore.pdf"), width = 4, height = 5)
print(cyt_plot)
dev.off()

# MSI p-value
cyt_df %>%
  filter(MSI_status == "MSI-H") %>%
  wilcox.test(cyt_score ~ KMT2D_groups, .)

# MSS p-value
cyt_df %>%
  filter(MSI_status == "MSS/MSI-L") %>%
  wilcox.test(cyt_score ~ KMT2D_groups, .)

# -----------------------------------------------------------------------------
# Calculate additional markers: INF-y signature score
# -----------------------------------------------------------------------------
#https://www.jci.org/articles/view/91190#SEC4 
# After performance of quantile normalization, a log10 transformation was applied, and signature scores were calculated by averaging of the included genes for the IFN-γ 
# A 10-gene “preliminary IFN-γ” signature (IFNG, STAT1, CCR5, CXCL9, CXCL10, CXCL11, IDO1, PRF1, GZMA, and MHCII HLA-DRA)
ifng <- c("IFNG", "STAT1", "CCR5", "CXCL9", "CXCL10", "CXCL11", "IDO1", "PRF1", "GZMA", "MHCII", "HLA-DRA")

ifng_df <- fpkm_long %>%
  filter(gene_name %in% ifng, KMT2D_groups %in% c("KMT2D_WT", "KMT2D_LOF")) %>%
  distinct() %>%
  group_by(Tumor_Sample_Barcode) %>%
  summarise(ifng_score = mean(fpkm)) %>%
  ungroup() %>%
  left_join(sample_combined_status %>% select(Tumor_Sample_Barcode, KMT2D_groups)) %>%
  distinct()

ifng_plot <- ifng_df %>%
  mutate(KMT2D_groups = fct_relevel(KMT2D_groups, "KMT2D_WT")) %>%
  ggplot(., aes(x = KMT2D_groups, y = ifng_score)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(height = 0, width = 0.2) +
  facet_grid(~MSI_status) +
  theme_bw()+
  theme(text = element_text(size = 12))

pdf(paste0(tcga_dir, "TCGA_COADREAD_KMT2Dgroups_additional_marker_comparision_IFNg_score.pdf"), width = 4, height = 5)
  print(ifng_plot)
dev.off()

# MSI p-value
ifng_df %>%
  filter(MSI_status == "MSI-H") %>%
  wilcox.test(cyt_score ~ KMT2D_groups, .)

# MSS p-value
ifng_df %>%
  filter(MSI_status == "MSS/MSI-L") %>%
  wilcox.test(cyt_score ~ KMT2D_groups, .)