#############################################################
# Figure 5
# TCGA COAD/READ cohort - Mutational signature analysis 
# and chronological age comparison
#############################################################
library(pacman)
p_load(tidyverse)
dir <- "MutationalSignature/"
tcga_dir <- "cBioPortal/"

#############################################################
# Load data
#############################################################
res_act <- read_tsv(paste0(dir, "MutationalSignature/Analysis/TCGA_WES_96/normalised/Assignment_Solution/Activities/Assignment_Solution_Activities.txt")) %>% janitor::clean_names()
tcga_sample_annot <- read_csv(paste0(tcga_dir, "TCGA-COADREAD_sample_status.csv"))

#############################################################
# Clean results table
#############################################################
res_act_clean <- res_act %>%
    filter(str_detect(samples, "ColoRect")) %>%
    mutate(
        cancer_type = str_split_fixed(samples, "::", 2)[, 1],
        patient_id = str_split_fixed(samples, "::", 2)[, 2],
        Tumor_Sample_Barcode = str_sub(patient_id, 1, 15)
    ) %>%
    select(-c(cancer_type, patient_id, samples)) %>%
    select(Tumor_Sample_Barcode, everything()) %>%
    filter(Tumor_Sample_Barcode %in% tcga_sample_annot$Tumor_Sample_Barcode) %>%
    left_join(tcga_sample_annot %>% select(Tumor_Sample_Barcode, MSI_status, KMT2D_groups)) %>%
    distinct() %>%
    pivot_longer(-c(Tumor_Sample_Barcode, MSI_status, KMT2D_groups), names_to = "sbs_type", values_to = "act")

res_act_clean %>%
    select(Tumor_Sample_Barcode, MSI_status, KMT2D_groups) %>%
    distinct() %>%
    filter(KMT2D_groups %in% c("KMT2D_WT", "KMT2D_LOF")) %>%
    count(KMT2D_groups, MSI_status)

res_act %>%
    filter(str_detect(samples, "ColoRect")) %>%
    mutate(
        cancer_type = str_split_fixed(samples, "::", 2)[, 1],
        patient_id = str_split_fixed(samples, "::", 2)[, 2],
        Tumor_Sample_Barcode = str_sub(patient_id, 1, 15)
    ) %>%
    select(-c(cancer_type, patient_id, samples)) %>%
        select(Tumor_Sample_Barcode, everything()) %>%
        filter(Tumor_Sample_Barcode %in% tcga_sample_annot$Tumor_Sample_Barcode) %>%
        left_join(tcga_sample_annot %>% select(Tumor_Sample_Barcode, MSI_status, KMT2D_groups)) %>%
        distinct() %>%
        filter(
            KMT2D_groups %in% c("KMT2D_WT", "KMT2D_LOF"),
            MSI_status %in% c("MSI-H", "MSS/MSI-L")) %>%
        mutate(KMT2D_group = ifelse(KMT2D_groups == "KMT2D_WT", "WT", "LOF")) %>% 
        select(Tumor_Sample_Barcode, MSI_status, KMT2D_group, sbs1:sbs95)

#############################################################
# T-tests
#############################################################
res_act_all <- NULL
for (i in seq(unique(res_act_clean$sbs_type))) {
    # i <- 1
    print(unique(res_act_clean$sbs_type)[i])

    df <- res_act_clean %>%
        filter(
            sbs_type == unique(res_act_clean$sbs_type)[i],
            MSI_status == "MSI-H",
            KMT2D_groups %in%  c("KMT2D_WT", "KMT2D_LOF"),
            !is.na(act),
            !is.nan(act)
        )

    if (all(df$act == 0)) {
        p = 1
    } else {
        fit <- t.test(act ~ KMT2D_groups, df)
        p <- fit$p.value
    }

    res_act <- tibble(
        sbs_type = unique(res_act_clean$sbs_type)[i],
        p_value = p
    )
    res_act_all <- bind_rows(res_act_all, res_act)
}
 res_act_all %>%
     arrange(p_value) %>%
     mutate(padj = p.adjust(p_value, method = "BH", length(p_value)))

#############################################################
# Compare age between groups
#############################################################

# -----------------------------------------------------------------------------
# Path to pan-cancer data
# -----------------------------------------------------------------------------
tcga_dir <- "cBioPortal/"
data_key <- tibble(
  dir_full = list.dirs(paste0(tcga_dir,"Datasets"), full.names = TRUE, recursive = FALSE),
  dir = list.dirs(paste0(tcga_dir,"Datasets"), full.names = FALSE, recursive = FALSE)) %>%
  mutate(TCGA_project = paste0("TCGA-",(str_split_fixed(dir,"_",2)[,1] %>% str_to_upper)))

# -----------------------------------------------------------------------------
# Load data 
# -----------------------------------------------------------------------------
data_clinical_patient <- read_tsv(paste0(data_key$dir_full[i],"/data_clinical_patient.txt"), skip = 4)
sample_combined_status <- read_csv(paste0(tcga_dir, "TCGA_COADREAD_sample_status.csv")) %>%
    left_join(data_clinical_patient %>% select(PATIENT_ID, AGE)) %>%
    filter(KMT2D_groups %in% c("KMT2D_WT", "KMT2D_LOF"))

# -----------------------------------------------------------------------------
# T-tests
# -----------------------------------------------------------------------------
t.test(AGE ~ KMT2D_groups, sample_combined_status %>% filter(MSI_status == "MSI-H"))

t.test(AGE ~ KMT2D_groups, sample_combined_status %>% filter(MSI_status == "MSS/MSI-L"))

pdf(paste0(tcga_dir, "TCGA_COADREAD_compare_age.pdf"), width = 4, height = 4)
sample_combined_status %>%
    filter(MSI_status %in% c("MSS/MSI-L", "MSI-H")) %>%
    mutate(
        KMT2D_groups = fct_relevel(KMT2D_groups, "KMT2D_WT"),
        MSI_status = fct_relevel(MSI_status, "MSS/MSI-L")) %>%
    select(PATIENT_ID, MSI_status, KMT2D_groups, AGE) %>%
    distinct %>%
    ggplot(., aes(x = KMT2D_groups, y = AGE))+
        geom_boxplot(outlier.shape = NA) +
            geom_jitter(width = 0.2, height = 0) +
            facet_wrap(~MSI_status) +
            theme_bw() +
            theme(
                text = element_text(size = 12)
            )
dev.off()