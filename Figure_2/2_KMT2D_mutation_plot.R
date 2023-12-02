#############################################################
# Figure 2A Preparing TCGA+ data to identify KMT2D LOF cases
#############################################################
library(pacman)
p_load(tidyverse, GRETTA)

# Define paths
DepMap_dir <- # Define path to DepMap 20Q1 data, see https://github.com/ytakemon/GRETTA for instructions on downloading GRETTA data.
  
## Load files -----------------------------------------------------------------
# see https://github.com/ytakemon/GRETTA for instructions on downloading GRETTA data.
load(paste0(dir,"sample_annot.rda"))
load(paste0(dir, "dep.rda"))

# -----------------------------------------------------------------------------
# Regroup MCL, DLBCL, and FL in to B_NHL
PanCan_KMT2D_muts <- read_csv("PanCan_KMT2D_LOF_mutation_proportions.csv") %>%
  mutate(TCGA_type = case_when(
    TCGA_type %in% c("DLBC", "FL", "MCL", "NMZL") ~ "NHL_B",
    TCGA_type %in% c("COAD", "READ") ~ "COAD/READ",
    TRUE ~ TCGA_type)) %>%
  group_by(TCGA_type) %>%
  summarise(
    Cases_with_all_muts = sum(Cases_with_all_muts),
    Cases_with_LOF_muts = sum(Cases_with_LOF_muts),
    Total_cases = sum(Total_cases)) %>% ungroup %>%
    janitor::adorn_totals("row") %>%
  mutate(
    perc = round((Cases_with_all_muts/Total_cases)*100),
    perc_LOF = round((Cases_with_LOF_muts/Cases_with_all_muts)*100)) %>%
  arrange(-perc)

# Manually annotate cancer type base on Huffogram/TCGA names
depmap_cell_lines <- sample_annot %>%
  filter(DepMap_ID %in% dep$DepMap_ID)

depmap_cell_lines_annoated <- depmap_cell_lines %>% count(disease, disease_subtype) %>%
  #(collapsed bellow)
  mutate(Cancer_type_names = case_when(
    # FL - NA but NHL available
    str_detect(disease, "Lymphoma") & str_detect(disease_subtype, "B-cell") ~ "NHL_B",
    # DLBCL / FL / NHL / NMZL
    # disease_subtype == "Diffuse Large B-cell Lymphoma (DLBCL)" ~ "DLBCL",
    # BLCA
    str_detect(disease, "Bladder") ~ "BLCA",
    # SKCM
    str_detect(disease, "Skin") & str_detect(disease_subtype, "Melanoma") ~ "SKCM",
    # ESCA
    str_detect(disease, "Esophageal") ~ "ESCA",
    # PAAD
    str_detect(disease_subtype, "Ductal Adenocarcinoma") ~ "PAAD",
    # SCLC
    str_detect(disease_subtype, fixed("(SCLC)")) ~ "SCLC",
    # MED
    str_detect(disease_subtype, fixed("Medulloblastoma")) ~ "MED",
    # LUSC
    str_detect(disease_subtype, fixed("(NSCLC), Squamous Cell Carcinoma")) ~ "LUSC",
    # STAD
    str_detect(disease, "Gastric") & str_detect(disease_subtype, "Adenocarcinoma") ~ "STAD",
    # UCEC
    str_detect(disease, "Uterine") & str_detect(disease_subtype, "Endometrial") &
    str_detect(disease_subtype, "arcinoma") ~ "UCEC",
    # HNSC
    str_detect(disease, "Head") & str_detect(disease_subtype, "Squamous") ~ "HNSC",
    # CESC
    str_detect(disease, "Cervical") & str_detect(disease_subtype, "Squamous") ~ "CESC",
    # MCL - NA
    # THYM - NA
    # LIHC
    str_detect(disease, "Liver") & str_detect(disease_subtype, "Hepatocellular") ~ "LIHC",
    # LUAD
    str_detect(disease, "Lung") & str_detect(disease_subtype, "Adenocarcinoma") ~ "LUAD",
    # UCS
    str_detect(disease, "Uterine") & str_detect(disease_subtype, "Carcinoma") &  !str_detect(disease_subtype, "Endometrial") ~ "UCS",
    # CHOL
    str_detect(disease_subtype, "Cholangio") ~ "CHOL",
    # ACC  - NA
    # KIRP
    str_detect(disease, "Kidney") & str_detect(disease_subtype, "Renal Cell Carcinoma") ~ "KIRP",
    # PRAD
    str_detect(disease, "Prostate") ~ "PRAD",
    # KICH - NA
    # KIRC
    str_detect(disease, "Kidney") & str_detect(disease_subtype, "clear cell") ~ "KIRC",
    # TGCT - NA
    # SARC
    str_detect(disease, "Sarcoma") ~ "SARC",
    # LGG
    str_detect(disease_subtype, "Oligodendroglioma") ~ "LGG",
    # UVM
    str_detect(disease_subtype, "Uveal Melanoma") ~ "UVM",
    # GBM
    str_detect(disease_subtype, "Glioblastoma") ~ "GBM",
    # COAD
    str_detect(disease, "Colon") ~ "COAD/READ",
    # BRCA
    str_detect(disease, "Breast") ~ "BRCA",
    # PCPG - NA
    # OV - Ovarian carcinoma
    str_detect(disease, "Ovarian") & str_detect(disease_subtype, "adeno|serous") ~ "OV",
    # LAML
    str_detect(disease_subtype, "(AML)") ~ "LAML",
    # MESO - NA
    # READ
    str_detect(disease, "rectal") ~ "COAD/READ",
    # THCA - NA
    TRUE ~ NA_character_)) %>%
    bind_rows(., tibble(Cancer_type_names = "Total", n = 100)) %>%
    write_csv(paste0(DepMap_dir,"DepMap_cell_lines_TCGA_type_annotated.csv"))

depmap_cell_lines_tally <- depmap_cell_lines_annoated %>%
  group_by(Cancer_type_names) %>%
  summarise(Total_lines = sum(n)) %>%
  arrange(-Total_lines)

Compiled_list <- PanCan_KMT2D_muts %>%
  left_join(depmap_cell_lines_tally, by = c("TCGA_type" = "Cancer_type_names")) %>%
  replace_na(list(Total_lines = 0))

# Count number of DepMap cell lines with KMT2D muts in each TCGA type
DepMap_KMT2D_muts <- select_cell_lines(input_gene = "KMT2D", data_dir = DepMap_dir) %>%
  filter(
    !Group %in% c("Control", "Others", "Amplified"))

DepMap_KMT2D_muts_annotated <- sample_annot %>%
  filter(DepMap_ID %in% DepMap_KMT2D_muts$DepMap_ID) %>%
  left_join(., DepMap_KMT2D_muts %>% select(DepMap_ID, Total_mutations:Group)) %>%
  left_join(.,  depmap_cell_lines_annoated %>% select(-n))

DepMap_KMT2D_muts_tally <- DepMap_KMT2D_muts_annotated %>%
  count(Cancer_type_names, Group) %>%
  pivot_wider(names_from = Group, values_from = n) %>%
  replace_na(list(KMT2D_mut_3 = 0, KMT2D_mut_2 = 0)) %>%
  mutate(Total_KMT2D_muts = KMT2D_mut_2 + KMT2D_mut_3) %>%
  arrange(-Total_KMT2D_muts) %>%
  rename(KMT2D_HetDel = KMT2D_mut_3, KMT2D_THetDel = KMT2D_mut_2)

# Plot -------------------------------------------------------------------------
type_order <- Compiled_list  %>% arrange(-perc) %>% pull(TCGA_type)
plot_df_reg <- Compiled_list %>%
  select(TCGA_type, Percent_mut = perc, Percent_LOF_mut = perc_LOF, DepMap_lines = Total_lines) %>%
  arrange(-Percent_mut) %>%
  mutate(DepMap_lines = case_when(
    TCGA_type == "Total" ~ 100,
    TRUE ~ DepMap_lines)) %>%
  pivot_longer(-TCGA_type, names_to = "data_type", values_to = "values") %>%
  mutate(
    TCGA_type = fct_relevel(TCGA_type, c(type_order)),
    TCGA_type = fct_rev(TCGA_type),
    data_type = fct_relevel(data_type,"Percent_mut","Percent_LOF_mut"))

Selected_groups <- Compiled_list %>%
  select(TCGA_type, Percent_mut = perc, Percent_LOF_mut = perc_LOF, DepMap_lines = Total_lines) %>%
  arrange(-Percent_mut) %>%
  filter(
    Percent_mut >= 10,
    Percent_LOF_mut >= 25,
    DepMap_lines >= 10) %>%
  write_csv(paste0(DepMap_dir, "TCGA_types_passing_threshold.csv"))

label_colour <- plot_df_reg %>%
  mutate(
    TCGA_type_colour = ifelse(TCGA_type %in% Selected_groups$TCGA_type,"dark red","black")) %>%
  select(TCGA_type, TCGA_type_colour) %>% distinct %>%
  pull(TCGA_type_colour) %>% rev

dummy_for_hline <- tibble(
  data_type = c("Percent_mut","Percent_LOF_mut", "DepMap_lines"), hline = c(10,25,10)) %>%
  mutate(data_type = fct_relevel(data_type,"Percent_mut","Percent_LOF_mut"))

pdf(paste0(DepMap_dir, "Compare_TCGA_DepMap_lines.pdf"), height = 6, width = 5)
ggplot(plot_df_reg, aes(x = TCGA_type, y = values, fill = data_type, alpha = data_type)) +
  geom_bar(stat = "identity")+
  # geom_hline(yintercept = 10, colour = "dark grey") +
  facet_grid(~data_type)+
  geom_hline(data = dummy_for_hline, aes(yintercept = hline), colour = "dark grey")+
  coord_flip() +
  scale_fill_manual(values = c("grey", "blue", "blue"))+
  scale_alpha_manual(values = c(1, 0.3, 0.6)) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(colour = label_colour),
    legend.position = "none",
    text = element_text(size = 12)) +
  scale_y_continuous(breaks = seq(0, 100, by = 25)) +
  labs(x = "", y = "")
dev.off()
