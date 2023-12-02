#############################################################
# Figure 5
# TCGA COAD/READ cohort - Annotating KMT2D WT and LOF mutations and MSI status
#############################################################
library(pacman)
p_load(tidyverse, colorspace, ggrepel)

# -----------------------------------------------------------------------------
# Path to pan-cancer data
# -----------------------------------------------------------------------------
# Data were downloaded from cBioPortal
tcga_dir <- "cBioPortal/"
data_key <- tibble(
  dir_full = list.dirs(paste0(tcga_dir,"Datasets"), full.names = TRUE, recursive = FALSE),
  dir = list.dirs(paste0(tcga_dir,"Datasets"), full.names = FALSE, recursive = FALSE)) %>%
  mutate(TCGA_project = paste0("TCGA-",(str_split_fixed(dir,"_",2)[,1] %>% str_to_upper)))

# -----------------------------------------------------------------------------
# Analyse each TCGA-tumour group
# -----------------------------------------------------------------------------
data_key <- data_key %>% filter(TCGA_project == "TCGA-COADREAD")
patient_count_all <- sample_msi_status_all <- NULL
print(data_key$TCGA_project[i])

# -----------------------------------------------------------------------------
# Pan-cancer data
# -----------------------------------------------------------------------------
data_mut <- read_tsv(paste0(data_key$dir_full[i],"/data_mutations.txt"))
data_cn <- read_tsv(paste0(data_key$dir_full[i],"/data_cna.txt"))
data_clinical_sample <- read_tsv(paste0(data_key$dir_full[i],"/data_clinical_sample.txt"), skip = 4)
data_clinical_patient <- read_tsv(paste0(data_key$dir_full[i],"/data_clinical_patient.txt"), skip = 4)

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
    TRUE ~ "Unknown")) %>%
  dplyr::select(Tumor_Sample_Barcode = SAMPLE_ID, MSI_SCORE_MANTIS, MSI_SENSOR_SCORE, MSI_status) %>%
  mutate(TCGA_type = data_key$TCGA_project[i]) %>%
  write_csv(paste0(tcga_dir, "Analysis/KMT2D_project/TCGA/", data_key$TCGA_project[i], "_sample_MSI_status.csv"))

msi_plot <- ggplot(sample_msi_status, aes( x = MSI_SCORE_MANTIS, y = MSI_SENSOR_SCORE))+
  geom_point(aes(colour = MSI_status)) +
  geom_hline(yintercept = 10, colour = "grey", linetype = "dashed")+ 
  geom_vline(xintercept = 0.4, colour = "grey", linetype = "dashed")+
  scale_colour_manual(values = c("#f95d6a", "#003f5c", "#008000", "grey"))+
  theme(text = element_text(size = 12))+
  theme_light()

# mutation annotation -----------------------------------------------
LOF_vc <- c("Nonsense_Mutation", "Frame_Shift_Ins", "Frame_Shift_Del")
# to capture all samples 
all_sample_id <- data_mut %>% dplyr::select(Tumor_Sample_Barcode, Matched_Norm_Sample_Barcode) %>% distinct

sample_mutation_status_LOF <- data_mut %>%
  dplyr::select(Hugo_Symbol, Tumor_Sample_Barcode, Matched_Norm_Sample_Barcode, Variant_Classification) %>%
  filter(
    Hugo_Symbol %in% c("KMT2D"),
    Variant_Classification %in% LOF_vc) %>%    
  mutate(LOF_status = Variant_Classification %in% LOF_vc)

sample_mutation_status_other <- data_mut %>%
  dplyr::select(Hugo_Symbol, Tumor_Sample_Barcode, Matched_Norm_Sample_Barcode, Variant_Classification) %>%
  filter(
    Hugo_Symbol %in% c("KMT2D"),
    !Variant_Classification %in% c(LOF_vc, "Silent"))

sample_mutation_status <- all_sample_id %>%
  mutate(
    KMT2D_mut = case_when(
      Tumor_Sample_Barcode %in% (sample_mutation_status_LOF %>% filter(Hugo_Symbol == "KMT2D") %>% pull(Tumor_Sample_Barcode)) ~ "LOF",
      Tumor_Sample_Barcode %in% (sample_mutation_status_other %>% filter(Hugo_Symbol == "KMT2D") %>% pull(Tumor_Sample_Barcode)) ~ "Other",
      TRUE ~ "WT"
    )
  )

# copy number annotation -----------------------------------------------
cn_key <- tibble(
  cn = c(-2, -1, 0, 1, 2),
  cn_type = c("DeepDel","ShallowDel","Neutral","ShallowAmp","Amp"))

sample_cn_status <- data_cn %>%
  filter(Hugo_Symbol %in% c("KMT2D")) %>% dplyr::select(Hugo_Symbol, starts_with("TCGA-")) %>%
  pivot_longer(-Hugo_Symbol, names_to = "Tumor_Sample_Barcode", values_to = "cn") %>%
  left_join(cn_key) %>%
  dplyr::select(Hugo_Symbol, Tumor_Sample_Barcode, cn_type) %>%
  pivot_wider(names_from = "Hugo_Symbol", values_from = "cn_type") %>%
  dplyr::select(Tumor_Sample_Barcode, KMT2D_cn = KMT2D, WRN_cn = WRN)
  
# Combine  -----------------------------------------------
sample_combined_status <- sample_mutation_status %>%
  left_join(sample_cn_status) %>%
  left_join(sample_msi_status %>% dplyr::select(Tumor_Sample_Barcode, MSI_status))  %>%
  mutate(
    KMT2D_group = case_when(
      KMT2D_mut == "LOF" | KMT2D_cn == "DeepDel" ~ "KMT2D_LOF",
      KMT2D_mut == "WT" & (KMT2D_cn == "ShallowDel" | KMT2D_cn == "Neutral"| KMT2D_cn == "ShallowAmp" | is.na(KMT2D_cn)) ~ "KMT2D_WT",
      TRUE ~ "Other")
    PATIENT_ID = str_sub(Tumor_Sample_Barcode, 0, 12)) %>%
    dplyr::select(PATIENT_ID, everything()) %>%
    write_csv(paste0(tcga_dir, "Analysis/KMT2D_project/TCGA/",data_key$TCGA_project[i],"_sample_status.csv"))
