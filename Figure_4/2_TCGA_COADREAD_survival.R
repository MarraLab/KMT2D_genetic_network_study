#############################################################
# Figure 4
# TCGA COAD/READ cohort - Annotating KMT2D WT and LOF mutations and MSI status
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

# -----------------------------------------------------------------------------
# Survival analysis
# -----------------------------------------------------------------------------
# Recoding:
# time: Observed survival time in days
# status: censoring status 1=censored, 2=deceased
# sex: 1=Male, 2=Female
sample_os_survival <- sample_combined_status %>%
  dplyr::select(-Matched_Norm_Sample_Barcode) %>%
  left_join(data_clinical_patient %>% dplyr::select(PATIENT_ID, SEX, OS_STATUS, OS_MONTHS)) %>%
  distinct() %>%
  mutate(
    sex_status = ifelse(SEX == "Female", 2, 1),
    censor_status = ifelse(OS_STATUS == "0:LIVING", 1, 2)
  )

## MSS cases -------------------------------------------------------------
os_surv_plot_wtable <- sample_os_survival %>%
  filter(
    KMT2D_groups %in% c("KMT2D_WT", "KMT2D_LOF"), 
    MSI_status == "MSS/MSI-L") %>%
  survfit2(Surv(OS_MONTHS, censor_status) ~ Combine_groups, data = .) %>%
  ggsurvfit() +
  labs(
    x = "Months",
    y = "Overall survival probability"
  ) +
  add_confidence_interval() +
  add_risktable() +
  theme(text = element_text(size = 12))

pdf(paste0(tcga_dir, "TCGA_COADREAD_MSS_KMT2D_survival_curve_OS_withRisktable.pdf"), width = 6, height = 5)
plot(os_surv_plot_wtable)
dev.off()

# Log-rank test
os_test <- sample_os_survival %>%
  filter(
    KMT2D_groups %in% c("KMT2D_WT", "KMT2D_LOF"),
     MSI_status == "MSS/MSI-L") %>%
  survdiff(Surv(OS_MONTHS, censor_status) ~ Combine_groups, data = .)

## MSI cases -------------------------------------------------------------
os_surv_plot_wtable <- sample_os_survival %>%
  filter(
    KMT2D_groups %in% c("KMT2D_WT", "KMT2D_LOF"), 
    MSI_status == "MSI-H") %>%
  survfit2(Surv(OS_MONTHS, censor_status) ~ Combine_groups, data = .) %>%
  ggsurvfit() +
  labs(
    x = "Months",
    y = "Overall survival probability"
  ) +
  add_confidence_interval() +
  add_risktable() +
  theme(text = element_text(size = 12))

pdf(paste0(tcga_dir, "TCGA_COADREAD_MSIH_KMT2D_survival_curve_OS_withRisktable.pdf"), width = 6, height = 5)
plot(os_surv_plot_wtable)
dev.off()

# Log-rank test
os_test <- sample_os_survival %>%
  filter(
    KMT2D_groups %in% c("KMT2D_WT", "KMT2D_LOF"),
     MSI_status == "MSI-H") %>%
  survdiff(Surv(OS_MONTHS, censor_status) ~ Combine_groups, data = .)
