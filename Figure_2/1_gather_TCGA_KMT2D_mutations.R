#############################################################
# Figure 2A Preparing TCGA+ data to identify KMT2D LOF cases
#############################################################
library(pacman)
p_load(tidyverse, readxl)

# Define paths
dir <- # Path to files 

# TCGA mutation files were downloaded from:
# https://gdc.cancer.gov/about-data/publications/pancanatlas
# Mutations - mc3.v0.2.8.PUBLIC.maf.gz
# TCGA-Clinical Data Resource (CDR) Outcome* - TCGA-CDR-SupplementalTableS1.xlsx

# Read and clean TCGA PanCan open MAF files for downstream use----------
maf_raw <- read_tsv(paste0(dir, "mc3.v0.2.8.PUBLIC.maf"))
clinical_data_raw <- read_xlsx(paste0(dir, "TCGA-CDR-SupplementalTableS1.xlsx"))

maf_annotated <- maf_raw %>%
  select(Hugo_Symbol:INTRON, FILTER) %>%
  mutate(
    Manual_centers = str_split(Tumor_Sample_Barcode, "-", simplify = T)[, 1],
    Patient_ID_pt2 = str_split(Tumor_Sample_Barcode, "-", simplify = T)[, 2],
    Paitnet_ID_pt3 = str_split(Tumor_Sample_Barcode, "-", simplify = T)[, 3],
    bcr_patient_barcode = paste0(Manual_centers, "-", Patient_ID_pt2, "-", Paitnet_ID_pt3)
  ) %>%
  filter(bcr_patient_barcode %in% clinical_data_raw$bcr_patient_barcode) %>%
  select(-c(Patient_ID_pt2, Patient_ID_pt3))

clinical_data_anotated <- clinical_data_raw %>%
  mutate(Manual_centers = str_split(bcr_patient_barcode, "-", simplify = T)[, 1]) %>%
  filter(bcr_patient_barcode %in% maf_annotated$bcr_patient_barcode)

maf_clinical_data <- left_join(maf_annotated, clinical_data_anotated, by = "bcr_patient_barcode") %>%
  rename(sex = gender) %>% 
  select(-c(Manual_centers.y)) %>%
  write_csv(paste0(dir, "processed_TCGA_data/TCGA_maf_clinical_data.csv"))

# Save for next step
saveRDS(maf_clinical_data, paste0(dir, "processed_TCGA_data/TCGA_maf_clinical_data.rds"))

# Identify mutants
muts_of_interest <- c("Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Splice_Site")
# Filter for KMT2D mutations
KMT2D_maf_df <- maf_clinical_data %>%
  filter(
    Hugo_Symbol == "KMT2D",
    Variant_Classification %in% muts_of_interest) %>%
  mutate(Variant_Classification = case_when(
    Variant_Classification %in% c("Frame_Shift_Del", "Frame_Shift_Ins") ~ "Frame_Shift_INDEL",
    Variant_Classification %in% c("In_Frame_Del", "In_Frame_Ins") ~ "In_Frame_INDEL",
    TRUE ~ Variant_Classification
  ))

# hg 18 KMT2D is from 47699025-47735374
# map to hg19
hg.convert <- function(x){
  return(round(((47735374 - x) / (47735374-47699025)) * (49453557-49412758)) + 49412758)
}
my.fixer <- function(x){return(as.numeric(as.character(x)))}

# Morin et al 2011 ; FL; [1] 47  6
# Mutation data from supplemental file: https://doi.org/10.1038/nature10351
morin <- read.delim('morin_2011/downloaded_genes.txt', header = F, sep = '\t',
                    stringsAsFactors = F)
morin <- cbind(morin, var = NA)
morin$var[grep('Frameshift', morin$V3)] <- "Frame_Shift_INDEL"
morin$var[grep('SS ', morin$V3)] <- "Splice_Site"
morin$var[grep('domain', morin$V3)] <- "Missense_Mutation"
morin$var[is.na(morin$var)] <- "Nonsense_Mutation"
morin <- cbind(morin, START_POSITION = do.call(rbind,strsplit(morin$V2, ":"))[,2])
morin$START_POSITION <- my.fixer(morin$START_POSITION)
morin$START_POSITION <- hg.convert(morin$START_POSITION)
morin_df <- morin %>%
  mutate(TCGA_type = "FL") %>%
  select(
    TCGA_type,
    Variant_Classification = var,
    Start = START_POSITION)
morin_all_mut_count <- morin %>% pull(V1) %>% unique %>% length
morin_LOF_mut_count <- morin %>% filter(var %in% c("Nonsense_Mutation","Frame_Shift_INDEL")) %>% pull(V1) %>% unique %>% length

# Parsons et al 2011 ; medulloblastoma (12/53)
# Mutation data from supplemental file: 10.1126/science.1198056
parsons <- read.delim('parsons_2011//downloaded_genes.txt', header = F, sep = '\t',
                     stringsAsFactors = F)
parsons <- parsons[parsons$V2 == "MLL2",]
parsons$V10[parsons$V10 == "INDEL"] <- "Frame_Shift_INDEL"
parsons$V10[parsons$V10 == "Nonsense"] <- "Nonsense_Mutation"
parsons$V10[parsons$V10 == "Missense"] <- "Missense_Mutation"
parsons <- cbind(parsons, START_POSITION = gsub("_.*|[A-Z].*|[a-z].*", "", do.call(rbind,strsplit(parsons$V7, ":"))[,2]))
parsons$START_POSITION <- my.fixer(parsons$START_POSITION)
parsons$START_POSITION <- hg.convert(parsons$START_POSITION)
parsons_df <- parsons %>%
  mutate(TCGA_type = "MED") %>%
  select(
    TCGA_type,
    Variant_Classification = V10,
    Start = START_POSITION)
parsons_all_mut_count <- parsons %>% pull(V6) %>% unique %>% length
parsons_LOF_mut_count <- parsons %>% filter(V10 %in% c("Nonsense_Mutation","Frame_Shift_INDEL")) %>% pull(V6) %>% unique %>% length

# spina et al 2016 ; nodal marginal zone lymphoma (13/51)
# Mutation data from supplemental file: 10.1182/blood-2016-02-696757
spina <- read.delim('spina_2016//downloaded_genes.txt', header = F, sep = '\t',
                    stringsAsFactors = F)
spina <- spina[spina$V9 == "MLL2", ]
spina$V11[grep('frameshift', spina$V11)] <- "Frame_Shift_INDEL"
spina$V11[spina$V11 == "nonsense"] <- "Nonsense_Mutation"
spina$V11[spina$V11 == "missense"] <- "Missense_Mutation"
spina$V11[grep('splice', spina$V11)] <- "Splice_Site"
colnames(spina)[4] <- "START_POSITION"
spina_df <- spina %>%
  mutate(TCGA_type = "NMZL") %>%
  select(
    TCGA_type,
    Variant_Classification = V11,
    Start = START_POSITION)
spina_all_mut_count <- spina %>% pull(V1) %>% unique %>% length
spina_LOF_mut_count <- spina %>% filter(V11 %in% c("Nonsense_Mutation","Frame_Shift_INDEL")) %>% pull(V1) %>% unique %>% length

# bea et al 2013 ; mantle cell lymphoma (4/29)
# Mutation data from supplemental file: 10.1073 / pnas.1314608110
bea <- read.delim('bea_2013//downloaded_genes.txt', header = F, sep = '\t',
                  stringsAsFactors = F)
bea <- bea[bea$V7 == "MLL2", ]
bea$V15[bea$V15 == "non_synonymous"] <- "Missense_Mutation"
bea$V15[bea$V15 == "fs_ins_2"] <- "Frame_Shift_INDEL"
colnames(bea)[3] <- "START_POSITION"
bea$START_POSITION <- my.fixer(bea$START_POSITION)
bea_df <- bea %>%
  mutate(TCGA_type = "MCL") %>%
  select(
    TCGA_type,
    Variant_Classification = V15,
    Start = START_POSITION)
bea_all_mut_count <- bea %>% pull(V1) %>% unique %>% length
bea_LOF_mut_count <- bea %>% filter(V15 %in% c("Nonsense_Mutation","Frame_Shift_INDEL")) %>% pull(V1) %>% unique %>% length

# george et al 2015; small cell lung cancer (/110)
# Mutation data from supplemental file: 10.1038 / nature14664
george_raw <- read_csv("george_2015/downloaded_genes.csv") %>%
  select(1:23) %>%
  mutate(variant_type = case_when(
    Type_1 == "missense" ~ "Missense_Mutation",
    Type_1 == "nonsense" ~ "Nonsense_Mutation",
    Type_1 %in% c("frame_shift_ins", "frame_shift_del") ~ "Frame_Shift_INDEL",
    Type_1 %in% c("splice") ~ "Splice_Site",
    Type_1 %in% c("nonstop") ~ "Nonstop_Mutation",
    TRUE ~ Type_1
  ))
george <- george_raw %>%
  filter(Gene_Hugo %in% c("KMT2D", "MLL2"))
george_df <- george %>%
  filter(variant_type != "silent") %>%
  mutate(TCGA_type = "SCLC") %>%
  select(
    TCGA_type,
    Variant_Classification = variant_type,
    Start)
george_all_patients <- george_raw %>% pull(PAT_ID) %>% unique %>% length
george_all_mut_count <- george %>% pull(PAT_ID) %>% unique %>% length
george_LOF_mut_count <- george %>% filter(variant_type %in% c("Nonsense_Mutation","Frame_Shift_INDEL")) %>% pull(PAT_ID) %>% unique %>% length

# KMT2D mutations by cancer type
Plot_df <- KMT2D_maf_df %>%
  select(
    TCGA_type = type, Start = Start_Position, Variant_Classification) %>%
  bind_rows(., morin_df) %>%
  bind_rows(., parsons_df) %>%
  bind_rows(., tan_df) %>%
  bind_rows(., spina_df) %>%
  bind_rows(., bea_df) %>%
  bind_rows(., george_df) %>%
  mutate(x_map = Start - 49412758) # map x-axis coordinates for each mutation

# map type to proportion of samples mutated
unique_cancer_type <- unique(maf_clinical_data$type)
LOF_muts <- c("Nonsense_Mutation", "Frame_Shift_Del","Frame_Shift_Ins","Nonstop_mutation")

map_prop <- NULL
for(i in 1:length(unique_cancer_type)){
  print(i)
  cancer_type_selected <- unique_cancer_type[i]

  cases_with_all_muts <- maf_clinical_data %>%
    filter(
      type == cancer_type_selected,
      Hugo_Symbol == "KMT2D",
      Variant_Classification %in% muts_of_interest) %>%
    pull(bcr_patient_barcode) %>% unique %>% length

  cases_with_LOF_muts <- maf_clinical_data %>%
    filter(
      type == cancer_type_selected,
      Hugo_Symbol == "KMT2D",
      Variant_Classification %in% LOF_muts) %>%
    pull(bcr_patient_barcode) %>% unique %>% length

  total_cases <- maf_clinical_data %>%
    filter(type == cancer_type_selected) %>%
    pull(bcr_patient_barcode) %>% unique %>% length

  map_prop <- bind_rows(
    map_prop,
    tibble(
      TCGA_type = cancer_type_selected,
      Cases_with_all_muts = cases_with_all_muts,
      Cases_with_LOF_muts = cases_with_LOF_muts,
      Total_cases = total_cases))
}

map_prop <- map_prop %>% bind_rows(. ,
    tibble(
      TCGA_type = "FL",
      Cases_with_all_muts = morin_all_mut_count,
      Cases_with_LOF_muts = morin_LOF_mut_count,
      Total_cases = 35)) %>%
  bind_rows(. ,
    tibble(
      TCGA_type = "MED",
      Cases_with_all_muts = parsons_all_mut_count,
      Cases_with_LOF_muts = parsons_LOF_mut_count,
      Total_cases = 53)) %>%
  bind_rows(. ,
    tibble(
      TCGA_type = "NMZL",
      Cases_with_all_muts = spina_all_mut_count,
      Cases_with_LOF_muts = spina_LOF_mut_count,
      Total_cases = 51)) %>%
  bind_rows(. ,
    tibble(
      TCGA_type = "MCL",
      Cases_with_all_muts = bea_all_mut_count,
      Cases_with_LOF_muts = bea_LOF_mut_count,
      Total_cases = 29)) %>%
  bind_rows(. ,
    tibble(
      TCGA_type = "SCLC",
      Cases_with_all_muts = george_all_mut_count,
      Cases_with_LOF_muts = george_LOF_mut_count,
      Total_cases = george_all_patients)) %>%
  mutate(
    perc = round((Cases_with_all_muts/Total_cases)*100),
    perc_LOF = round((Cases_with_LOF_muts/Cases_with_all_muts)*100),
    perc_LOF = case_when(
      is.nan(perc_LOF) ~ 0,
      TRUE ~ perc_LOF),
    label = paste0(perc, "% (", Cases_with_all_muts ,"/",Total_cases,") | ", perc_LOF,"% (",Cases_with_LOF_muts,"/",Cases_with_all_muts,")")) %>%
  arrange(perc) %>%
  write_csv("PanCan_KMT2D_LOF_mutation_proportions.csv")