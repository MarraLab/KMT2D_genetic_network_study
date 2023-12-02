#############################################################
# Figure 4
# TCGA COAD cohort analysis for MDM2 and TP53 - find samples
#############################################################
library(pacman)
p_load(tidyverse, colorspace, ggrepel)

# -----------------------------------------------------------------------------
# Path to pan-cancer data
# -----------------------------------------------------------------------------
# Data were downloaded from cBioPortal
tcga_dir <- "cBioPortal/"
ref <- read_csv("BioMart/All_genes_EnsID_Hugo_UCSC_NCBI_HGNC_ids.csv") %>% janitor::clean_names()
data_key <- tibble(
  dir_full = list.dirs(paste0(tcga_dir,"Datasets"), full.names = TRUE, recursive = FALSE),
  dir = list.dirs(paste0(tcga_dir,"Datasets"), full.names = FALSE, recursive = FALSE)) %>%
  mutate(TCGA_project = paste0("TCGA-",(str_split_fixed(dir,"_",2)[,1] %>% str_to_upper)))

# -----------------------------------------------------------------------------
# Analyse each TCGA-COAD/READ
# -----------------------------------------------------------------------------
data_key <- data_key %>% filter(TCGA_project == "TCGA-COADREAD")
print(data_key$TCGA_project[i])

# -----------------------------------------------------------------------------
# Pan-cancer data
# -----------------------------------------------------------------------------
data_mut <- read_tsv(paste0(data_key$dir_full[i], "/data_mutations.txt"))
data_cn <- read_tsv(paste0(data_key$dir_full[i], "/data_cna.txt"))
data_clinical_sample <- read_tsv(paste0(data_key$dir_full[i], "/data_clinical_sample.txt"), skip = 4)
data_clinical_patient <- read_tsv(paste0(data_key$dir_full[i], "/data_clinical_patient.txt"), skip = 4)

# -----------------------------------------------------------------------------
# Annotating samples
# -----------------------------------------------------------------------------
# msi scores - determine cut off based on consensus of two scores
# MSIsensor >= 10 is MSI-H http://ascopubs.org/doi/full/10.1200/PO.17.00084
# MANTIS >= 0.4 is MSI-H https://doi.org/10.18632/oncotarget.13918
sample_msi_status <- data_clinical_sample %>%
  mutate(MSI_status = case_when(
    (MSI_SCORE_MANTIS >= 0.4 & MSI_SENSOR_SCORE >= 10) ~ "MSI-H",
    (MSI_SCORE_MANTIS < 0.4 & MSI_SENSOR_SCORE < 10) ~ "MSS/MSI-L",
    !((MSI_SCORE_MANTIS >= 0.4 & MSI_SENSOR_SCORE >= 10) | (MSI_SCORE_MANTIS < 0.4 & MSI_SENSOR_SCORE < 10)) ~ "No_consensus",
    TRUE ~ "Unknown"
  )) %>%
  dplyr::select(Tumor_Sample_Barcode = SAMPLE_ID, MSI_SCORE_MANTIS, MSI_SENSOR_SCORE, MSI_status) %>%
  mutate(TCGA_type = data_key$TCGA_project[i]) %>%
  write_csv(paste0(tcga_dir, "Analysis/KMT2D_project/COAD/", data_key$TCGA_project[i], "_sample_MSI_status.csv"))

# check if MSIH cases are available, if not move to next.
check_msi_avail <- sample_msi_status %>%
  dplyr::count(MSI_status) %>%
  filter(MSI_status == "MSI-H")
if (nrow(check_msi_avail) == 0) {
  next
} else if (check_msi_avail$n < 3) {
  next
}

# mutation annotation -----------------------------------------------
LOF_vc <- c("Nonsense_Mutation", "Frame_Shift_Ins", "Frame_Shift_Del")
# to capture all samples
all_sample_id <- data_mut %>%
  dplyr::select(Tumor_Sample_Barcode, Matched_Norm_Sample_Barcode) %>%
  distinct()

sample_mutation_status_LOF <- data_mut %>%
  dplyr::select(Hugo_Symbol, Tumor_Sample_Barcode, Matched_Norm_Sample_Barcode, Variant_Classification) %>%
  filter(
    Hugo_Symbol %in% c("KMT2D", "MDM2", "TP53"),
    Variant_Classification %in% LOF_vc
  ) %>%
  mutate(LOF_status = Variant_Classification %in% LOF_vc)

sample_mutation_status_other <- data_mut %>%
  dplyr::select(Hugo_Symbol, Tumor_Sample_Barcode, Matched_Norm_Sample_Barcode, Variant_Classification) %>%
  filter(
    Hugo_Symbol %in% c("KMT2D", "MDM2", "TP53"),
    !Variant_Classification %in% c(LOF_vc, "Silent")
  )

sample_mutation_status <- all_sample_id %>%
  mutate(
    KMT2D_mut = case_when(
      Tumor_Sample_Barcode %in% (sample_mutation_status_LOF %>% filter(Hugo_Symbol == "KMT2D") %>% pull(Tumor_Sample_Barcode)) ~ "LOF",
      Tumor_Sample_Barcode %in% (sample_mutation_status_other %>% filter(Hugo_Symbol == "KMT2D") %>% pull(Tumor_Sample_Barcode)) ~ "Other",
      TRUE ~ "WT"
    ),
    MDM2_mut = case_when(
      Tumor_Sample_Barcode %in% (sample_mutation_status_LOF %>% filter(Hugo_Symbol == "MDM2") %>% pull(Tumor_Sample_Barcode)) ~ "LOF",
      Tumor_Sample_Barcode %in% (sample_mutation_status_other %>% filter(Hugo_Symbol == "MDM2") %>% pull(Tumor_Sample_Barcode)) ~ "Other",
      TRUE ~ "WT"
    ),
    TP53_mut = case_when(
      Tumor_Sample_Barcode %in% (sample_mutation_status_LOF %>% filter(Hugo_Symbol == "TP53") %>% pull(Tumor_Sample_Barcode)) ~ "LOF",
      Tumor_Sample_Barcode %in% (sample_mutation_status_other %>% filter(Hugo_Symbol == "TP53") %>% pull(Tumor_Sample_Barcode)) ~ "Other",
      TRUE ~ "WT"
    )
  )

# copy number annotation -----------------------------------------------
cn_key <- tibble(
  cn = c(-2, -1, 0, 1, 2),
  cn_type = c("DeepDel", "ShallowDel", "Neutral", "ShallowAmp", "Amp")
)

sample_cn_status <- data_cn %>%
  filter(Hugo_Symbol %in% c("KMT2D", "MDM2", "TP53")) %>%
  dplyr::select(Hugo_Symbol, starts_with("TCGA-")) %>%
  pivot_longer(-Hugo_Symbol, names_to = "Tumor_Sample_Barcode", values_to = "cn") %>%
  left_join(cn_key) %>%
  dplyr::select(Hugo_Symbol, Tumor_Sample_Barcode, cn_type) %>%
  pivot_wider(names_from = "Hugo_Symbol", values_from = "cn_type") %>%
  dplyr::select(Tumor_Sample_Barcode, KMT2D_cn = KMT2D, WRN_cn = WRN, MDM2_cn = MDM2, TP53_cn = TP53, TUBA1B_cn = TUBA1B, SAA1_cn = SAA1)

# Combine  -----------------------------------------------
sample_combined_status <- sample_mutation_status %>%
  left_join(sample_cn_status) %>%
  left_join(sample_msi_status %>% dplyr::select(Tumor_Sample_Barcode, MSI_status)) %>%
  mutate(
    KMT2D_group = case_when(
      KMT2D_mut == "LOF" | KMT2D_cn == "DeepDel" ~ "KMT2D_LOF",
      KMT2D_mut == "WT" & (KMT2D_cn == "ShallowDel" | KMT2D_cn == "Neutral" | KMT2D_cn == "ShallowAmp" | is.na(KMT2D_cn)) ~ "KMT2D_WT",
      TRUE ~ "Other"
    )
    MDM2_group = case_when(
      MDM2_mut == "LOF" | MDM2_cn == "DeepDel" ~ "MDM2_LOF",
      MDM2_mut == "WT" & (MDM2_cn == "ShallowDel" | MDM2_cn == "Neutral" | MDM2_cn == "ShallowAmp" | is.na(MDM2_cn)) ~ "MDM2_WT",
      TRUE ~ "Other"
    ),
    TP53_group = case_when(
      TP53_mut == "LOF" | TP53_cn == "DeepDel" ~ "TP53_LOF",
      TP53_mut == "WT" & (TP53_cn == "ShallowDel" | TP53_cn == "Neutral" | TP53_cn == "ShallowAmp" | is.na(TP53_cn)) ~ "TP53_WT",
      TRUE ~ "Other"
    )
    Combine_groups = paste0(KMT2D_group, "-", WRN_group),
    PATIENT_ID = str_sub(Tumor_Sample_Barcode, 0, 12)
  ) %>%
  dplyr::select(PATIENT_ID, everything()) %>%
  write_csv(paste0(tcga_dir, "Analysis/KMT2D_project/COAD/", data_key$TCGA_project[i], "_sample_status.csv"))

patient_count <- sample_combined_status %>%
  dplyr::select(PATIENT_ID, MSI_status, Combine_groups) %>%
  distinct() %>%
  dplyr::count(MSI_status, Combine_groups) %>%
  mutate(TCGA_project = data_key$TCGA_project[i])

patient_count_all <- bind_rows(patient_count_all, patient_count)
write_csv(patient_count_all, paste0(tcga_dir, "Analysis/KMT2D_project/COAD/eachTCGA_KMT2D_patient_N.csv"))

#############################################################
# Figure 4
# TCGA COAD cohort analysis for MDM2 and TP53 - RNA
#############################################################
# TPM was processed and download from UCSC XenaHub:
# https://xenabrowser.net/datapages/?dataset=tcga_RSEM_gene_tpm&host=https%3A%2F%2Ftoil.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
df_coad <- read_tsv(paste0(tcga_dir, "Datasets/TCGA-COAD.htseq_fpkm.tsv"))
df_read <- read_tsv(paste0(tcga_dir, "Datasets/TCGA-READ.htseq_fpkm.tsv"))

# -----------------------------------------------------------------------------
# Extract TCGA sample TPM counts for analysis
# -----------------------------------------------------------------------------
samples_to_extract <- sample_combined_status %>%
  filter(
    TCGA_project == "TCGA-COADREAD",
    KMT2D_group %in% c("KMT2D_WT", "KMT2D_LOF")
  )

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
df_coadread_long_subset <- bind_rows(df_coad_long_subset, df_read_long_subset)

df_coadread_long_subset %>%
  dplyr::select(Tumor_Sample_Barcode) %>%
  distinct() %>%
  dim()
# Figure out and clean duplicates
test <- df_coadread_long_subset %>%
  distinct() %>%
  dplyr::count(Ensembl_ID, Tumor_Sample_Barcode) %>%
  arrange(-n)
multiple <- test %>%
  filter(n == 2) %>%
  pull(Tumor_Sample_Barcode) %>%
  unique()
clean <- df_coad_long_subset %>%
  filter(Tumor_Sample_Barcode %in% multiple) %>%
  group_by(Ensembl_ID, Tumor_Sample_Barcode) %>%
  dplyr::slice(1) %>%
  ungroup()

df_coadread_long_subset_clean <- df_coadread_long_subset %>%
  filter(!Tumor_Sample_Barcode %in% multiple) %>%
  bind_rows(clean)
# check
df_coadread_long_subset_clean %>%
  distinct() %>%
  dplyr::count(Ensembl_ID, Tumor_Sample_Barcode) %>%
  arrange(-n) %>%
  filter(n == 2) %>%
  pull(Tumor_Sample_Barcode) %>%
  unique()

df_coadread_long_subset_clean %>%
  dplyr::select(Tumor_Sample_Barcode) %>%
  distinct() %>%
  dim()

df_coadread_long_subset_clean %>%
  mutate(EnsID_short = str_split_fixed(Ensembl_ID, "\\.", 2)[, 1]) %>%
  left_join(ref %>% select(EnsID_short = gene_stable_id, gene_name)) %>%
  write_csv(paste0(tcga_dir, "Analysis/KMT2D_project/COAD/TCGA_COADREAD_log2fpkm.csv"))

df_coadread_long_subset_clean_wide <- df_coadread_long_subset_clean %>%
  mutate(fpkm = (2^log2fpkm) - 1) %>%
  dplyr::select(-log2fpkm) %>%
  pivot_wider(names_from = Tumor_Sample_Barcode, values_from = fpkm)

coadread_tpm_wide <- df_coadread_long_subset_clean_wide %>%
  mutate_if(is.numeric, ~ (. / sum(.)) * 10**6) %>%
  mutate_if(is.numeric, ~ round((. / sum(.)) * 10**6, 2))

coadread_tpm_long <- coadread_tpm_wide %>%
  pivot_longer(-Ensembl_ID, names_to = "Tumor_Sample_Barcode", values_to = "tpm") %>%
  mutate(EnsID_short = str_split_fixed(Ensembl_ID, "\\.", 2)[, 1]) %>%
  left_join(ref %>% select(EnsID_short = gene_stable_id, gene_name)) %>%
  write_csv(paste0(tcga_dir, "Analysis/KMT2D_project/COAD/TCGA_COADREAD_tpm.csv"))

# -----------------------------------------------------------------------------
# Compare MDM2 expression 
# -----------------------------------------------------------------------------
coadread_tpm_long <- read_csv(paste0(tcga_dir, "Analysis/KMT2D_project/COAD/TCGA_COADREAD_tpm.csv"))

MDM2_tpm <- coadread_tpm_long %>%
  filter(gene_name == "MDM2") %>%
  left_join(samples_to_extract %>% select(Tumor_Sample_Barcode, KMT2D_group)) %>%
  select(Tumor_Sample_Barcode, EnsID_short, gene_name, tpm, KMT2D_group) %>%
  mutate(KMT2D_group = fct_relevel(KMT2D_group, "KMT2D_WT")) %>%
  distinct() %>%
  write_csv(paste0(tcga_dir, "Analysis/KMT2D_project/COAD/TCGA_COADREAD_MDM2_tpm.csv"))

pdf(paste0(tcga_dir, "Analysis/KMT2D_project/COAD/TCGA_COADREAD_MDM2_tpm.pdf"), width = 3, height = 5)
ggplot(MDM2_tpm, aes(x = KMT2D_group, y = tpm)) +
  geom_boxplot() +
  geom_jitter(height = 0, width = 0.2, shape = 16) +
  theme_bw() +
  labs(title = "MDM2") +
  theme(text = element_text(size = 12))
dev.off()

# T-test
t.test(tpm ~ KMT2D_group, MDM2_tpm)

# -----------------------------------------------------------------------------
# Compare TP53 expression
# -----------------------------------------------------------------------------
coadread_tpm_long <- read_csv(paste0(tcga_dir, "Analysis/KMT2D_project/COAD/TCGA_COADREAD_tpm.csv"))

TP53_tpm <- coadread_tpm_long %>%
  filter(gene_name == "TP53") %>%
  left_join(samples_to_extract %>% select(Tumor_Sample_Barcode, KMT2D_group)) %>%
  select(Tumor_Sample_Barcode, EnsID_short, gene_name, tpm, KMT2D_group) %>%
  mutate(KMT2D_group = fct_relevel(KMT2D_group, "KMT2D_WT")) %>%
  distinct() %>%
  write_csv(paste0(tcga_dir, "Analysis/KMT2D_project/COAD/TCGA_COADREAD_TP53_tpm.csv"))

pdf(paste0(tcga_dir, "Analysis/KMT2D_project/COAD/TCGA_COADREAD_TP53_tpm.pdf"), width = 3, height = 5)
ggplot(TP53_tpm, aes(x = KMT2D_group, y = tpm)) +
  geom_boxplot() +
  geom_jitter(height = 0, width = 0.2, shape = 16) +
  theme_bw() +
  labs(title = "TP53") +
  theme(text = element_text(size = 12))
dev.off()

# T-test
t.test(tpm ~ KMT2D_group, TP53_tpm)
TP53_tpm %>% dplyr::count(KMT2D_group)

# -----------------------------------------------------------------------------
# Compare downstream targets of TP53 transcriptional repressor - MYC, Cyclin B, VEGF, RAD51, and hTERT
# https://www.nature.com/articles/4402058#Sec7
# https://www.nature.com/articles/cdd2017174#Sec6
# -----------------------------------------------------------------------------
TP53DS_tpm <- coadread_tpm_long %>%
  filter(gene_name %in% c("MYC", "CCNB1", "VEGFA", "RAD51", "TERT", "CDKN1A", "BIRC5", "CDC25B", "CDC25C")) %>%
  left_join(samples_to_extract %>% select(Tumor_Sample_Barcode, KMT2D_group)) %>%
  select(Tumor_Sample_Barcode, EnsID_short, gene_name, tpm, KMT2D_group) %>%
  mutate(KMT2D_group = fct_relevel(KMT2D_group, "KMT2D_WT")) %>%
  distinct() %>%
  write_csv(paste0(tcga_dir, "Analysis/KMT2D_project/COAD/TCGA_COADREAD_TP53DS_tpm.csv"))

pdf(paste0(tcga_dir, "Analysis/KMT2D_project/COAD/TCGA_COADREAD_TP53DS_tpm.pdf"), width = 8, height = 11)
ggplot(TP53DS_tpm, aes(x = KMT2D_group, y = tpm)) +
  geom_boxplot() +
  geom_jitter(height = 0, width = 0.2, shape = 16) +
  facet_wrap(~gene_name, scales = "free_y") +
  theme_bw() +
  labs(title = "TP53 repressed genes") +
  theme(text = element_text(size = 12))
dev.off()

all_pval <- NULL
for (i in c("MYC", "CCNB1", "VEGFA", "RAD51", "TERT", "CDKN1A", "BIRC5", "CDC25B", "CDC25C")) {
  print(i)
  df <- TP53DS_tpm %>%
    filter(gene_name == i)
  pval <- t.test(tpm ~ KMT2D_group, df)$p.value

  res <- tibble(
    gene_name = i,
    pval = pval
  )
  all_pval <- bind_rows(all_pval, res)
}
all_pval <- all_pval %>%
  mutate(padj = p.adjust(pval, method = "BH"))
all_pval %>% arrange(gene_name)
