#############################################################
# Figure 2 / Supplementary Figure S1
# Comparing KMT2D mRNA/protein and Histone Modification 
# between KMT2D WT & LOF lines
#############################################################
library(pacman)
p_load(tidyverse, GRETTA, ggrepel, janitor)

# Define paths
DepMap_dir <- # Define path to DepMap 20Q1 data, see https://github.com/ytakemon/GRETTA for instructions on downloading GRETTA data.
  
## Load files -----------------------------------------------------------------
# see https://github.com/ytakemon/GRETTA for instructions on downloading GRETTA data.
load(paste0(DepMap_dir, "sample_annot.rda"))

# Supplemental Table S3
Groups <- read_csv(paste0(DepMap_dir, "AllGene_Groups_pan_cancer.csv")) 

# Created in 2_KMT2D_mutation_plot
samples_by_TCGA_type <- read_csv(paste0(DepMap_dir,"DepMap_cell_lines_TCGA_type_annotated.csv")) 
passing_types <- read_csv(paste0(DepMap_dir, "TCGA_types_passing_threshold.csv")) # Created in 2_KMT2D_mutation_plot
# Downloaded from DepMap data portal
histone_marks <- read_csv(paste0(DepMap_dir, "CCLE_GlobalChromatinProfiling_20181130.csv")) %>%
  rename(DepMap_ID = BroadID)

# ------------------------------------------------------------------
# Select cell lines that meet threshold by cancer type
selected_cell_lines <- sample_annot %>%
  left_join(., samples_by_TCGA_type %>% select(-n)) %>%
  filter(
    Cancer_type_names %in% c("Pan_cancer", "BLCA", "COAD", "ESCA", "HNSC", "LUSC", "NHL_B", "SCLC", "STAD", "UCEC", "NHL_B", "COAD/READ"))

# Grab DepMap cell lines with KMT2D mutations
DepMap_KMT2D_muts <- select_cell_lines(input_gene = "KMT2D", data_dir = DepMap_dir) %>%
  filter(
    !Group %in% c("Others", "Amplified")
  )

# Annotate
annotated_lines <- selected_cell_lines %>%
  left_join(DepMap_KMT2D_muts) %>%
  filter(str_detect(Group, "Control|KMT2D_"))

# ----------------------------------------------------------------------------
# Select cell lines that meet threshold by cancer type
selected_cell_lines <- sample_annot %>%
  left_join(.,samples_by_TCGA_type %>% select(-n)) %>%
  filter(
    Cancer_type_names %in% c("Pan_cancer", "BLCA", "COAD", "ESCA", "HNSC", "LUSC", "NHL_B", "SCLC", "STAD", "UCEC", "NHL_B", "COAD/READ"))

# Collect all KMT2D mutations in these cell lines
KMT2D_mutations <- Groups %>%
  filter(DepMap_ID %in% selected_cell_lines$DepMap_ID) %>%
  select(DepMap_ID, Total_mutations:Group)

# Combinselected_cell_linese the two and filter for control and mutants only
annotated_lines <- selected_cell_lines %>% left_join(KMT2D_mutations) %>%
  filter(str_detect(Group, "Control|KMT2D_"))

# -----------------------------------------------------------------------------
# Get expression data for all cell lines
Target_gene <- "KMT2D"
Target_cell_lines <- annotated_lines

# mRNA
KMT2D_rna_expr <- extract_rna(
  input_samples = annotated_lines$DepMap_ID,
  input_genes = "KMT2D",
  data_dir = DepMap_dir
)

# Protein
KMT2D_prot_expr <- extract_prot(
  input_samples = annotated_lines$DepMap_ID,
  input_genes = "KMT2D",
  data_dir = DepMap_dir
)

# Add to annotate data
annotated_lines <- left_join(annotated_lines, KMT2D_rna_expr) %>%
  left_join(., KMT2D_prot_expr)

Target_cell_lines <- Target_cell_lines %>%
  left_join(KMT2D_rna_expr) %>%
  left_join(KMT2D_prot_expr)

# --------------------------------------------------------------------------------
# Compre Mut v Control by disease type
# Combine all mutant groups into one
Target_group <- "KMT2D_LOF"
Target_cell_lines_combined_muts <- Target_cell_lines %>% mutate(
  Cell_line_status = case_when(
    Group == "Control" ~ "Control",
    str_detect(Group, "mut") ~ Target_group)) %>%
  write_csv(paste0(DepMap_dir, "DepMap_cell_lines_TCGA_type_combined_KMT2D_muts_annotated.csv"))

plot_df_subgroups <- Target_cell_lines_combined_muts %>%
  select(DepMap_ID, Cancer_type_names, Cell_line_status, KMT2D_protein, KMT2D_rna)
plot_df_pancan <- Target_cell_lines_combined_muts %>%
  filter(Group %in% c("Control","KMT2D_mut_2")) %>%
  select(DepMap_ID, Cell_line_status, KMT2D_protein, KMT2D_rna) %>%
  mutate(Cancer_type_names = "Pan_cancer")

plot_df <- bind_rows(plot_df_pancan, plot_df_subgroups) %>%
  mutate(Cancer_type_names = fct_relevel(Cancer_type_names, "Pan_cancer", "BLCA", "NHL_B", "COAD/READ", "ESCA", "HNSC", "LUSC", "SCLC", "STAD", "UCEC"))

pdf(paste0(DepMap_dir, "DepMap_cell_lines_TCGA_type_combined_KMT2D_muts_RNA_boxplot.pdf"), width = 7, height = 8)
  ggplot(plot_df, aes(x = Cell_line_status, y = KMT2D_rna))+
    geom_boxplot(
      outlier.shape = NA)+
    geom_jitter(
      shape = 16, colour = "black", alpha = 0.5, 
      width = 0.1, height = 0)+
    facet_wrap(~Cancer_type_names, ncol = 3, nrow = 4)+
    theme_bw()+
    theme(
      text = element_text(size = 12))
dev.off()

pdf(paste0(DepMap_dir, "DepMap_cell_lines_TCGA_type_combined_KMT2D_muts_protein_boxplot.pdf"), width = 7, height = 8)
ggplot(plot_df, aes(x = Cell_line_status, y = KMT2D_protein))+
  geom_boxplot(
    outlier.shape = NA)+
  geom_jitter(
    shape = 16, colour = "black", alpha = 0.5, 
    width = 0.1, height = 0)+
  facet_wrap(~Cancer_type_names, ncol = 3, nrow = 4)+
  theme_bw()+
  theme(
    text = element_text(size = 12))
dev.off()

# T.tests -----------------------------------------------------------------
plot_df_DE_summary <- NULL
for(select_TCGA_type in unique(plot_df$Cancer_type_names)){
  # select_TCGA_type <- "Lung Cancer"
  print(select_TCGA_type)

  # T-test Protein
  if(any(protein_nodup$Gene_Symbol %in% Target_gene) & !is.null(extract_target_protein)){
    plot_df_prot <- plot_df %>%
      filter(
        !is.na(!!sym(paste0(Target_gene,"_protein"))),
        Cancer_type_names %in% select_TCGA_type)

    if(!any(table(plot_df_prot$Cell_line_status) < 2)  & length(table(plot_df_prot$Cell_line_status)) == 2){
      # Dynamic t.tests
      variable <- paste0(Target_gene,"_rna")
      by <- "Cell_line_status"
      exp1 <- expr(!!ensym(variable) ~ !!ensym(by))
      fit_prot <- t.test(formula = eval(exp1), data = plot_df_prot)

      prot_mean <- NULL
      prot_mean$Control <- plot_df_prot %>% filter(Cell_line_status == "Control") %>% pull(!!sym(paste0(Target_gene,"_protein"))) %>% mean
      prot_mean$Mutant <- plot_df_prot %>% filter(Cell_line_status == Target_group) %>% pull(!!sym(paste0(Target_gene,"_protein"))) %>% mean

      prot_N <- NULL
      prot_N$Control <- plot_df_prot %>% filter(Cell_line_status == "Control", !is.na(!!sym(paste0(Target_gene,"_protein")))) %>% nrow()
      prot_N$Mutant <- plot_df_prot %>% filter(Cell_line_status == Target_group, !is.na(!!sym(paste0(Target_gene,"_protein")))) %>% nrow()

    } else {
      fit_prot <- NULL
      fit_prot$p.value <- NA

      prot_mean <- NULL
      prot_mean$Control <- NA
      prot_mean$Mutant <- NA

      prot_N <- NULL
      prot_N$Control <- NA
      prot_N$Mutant <- NA
    }
  } else {
    fit_prot <- NULL
    fit_prot$p.value <- NA

    prot_mean <- NULL
    prot_mean$Control <- NA
    prot_mean$Mutant <- NA

    prot_N <- NULL
    prot_N$Control <- NA
    prot_N$Mutant <- NA
  }

  # RNA -----------------------------------------
  plot_df_rna <- plot_df %>%
    filter(
      !is.na(!!sym(paste0(Target_gene,"_rna"))),
      Cancer_type_names %in% select_TCGA_type)

  if(!any(table(plot_df_rna$Cell_line_status) < 2) & length(table(plot_df_rna$Cell_line_status)) == 2){
    # Dynamic t.tests
    variable <- paste0(Target_gene,"_rna")
    by <- "Cell_line_status"
    exp1 <- expr(!!ensym(variable) ~ !!ensym(by))
    fit_rna <- t.test(formula = eval(exp1), data = plot_df_rna)

    rna_mean <- NULL
    rna_mean$Control <- plot_df_rna %>% filter(Cell_line_status == "Control") %>% pull(!!sym(paste0(Target_gene,"_rna"))) %>% mean
    rna_mean$Mutant <- plot_df_rna %>% filter(Cell_line_status == Target_group) %>% pull(!!sym(paste0(Target_gene,"_rna"))) %>% mean

    rna_N <- NULL
    rna_N$Control <- plot_df_rna %>% filter(Cell_line_status == "Control", !is.na(!!sym(paste0(Target_gene,"_rna")))) %>% nrow()
    rna_N$Mutant <- plot_df_rna %>% filter(Cell_line_status == Target_group, !is.na(!!sym(paste0(Target_gene,"_rna")))) %>% nrow()

  } else {
    fit_rna <- NULL
    fit_rna$p.value <- NA

    rna_mean <- NULL
    rna_mean$Control <- NA
    rna_mean$Mutant <- NA

    rna_N <- NULL
    rna_N$Control <- NA
    rna_N$Mutant <- NA
  }

  DE_summary <- tribble(
    ~RNA_control_mean, ~RNA_mutant_mean,
    ~RNA_control_N, ~RNA_mutant_N, ~RNA_pval,
    ~Protein_control_mean, ~Protein_mutant_mean,
    ~Protein_control_N, ~Protein_mutant_N, ~Protein_pval,
    rna_mean$Control, rna_mean$Mutant,
    rna_N$Control, rna_N$Mutant, fit_rna$p.value,
    prot_mean$Control, prot_mean$Mutant,
    prot_N$Control, prot_N$Mutant, fit_prot$p.value) %>%
    mutate(
      GeneName = Target_gene,
      TCGA_type = select_TCGA_type,
      Mutant_group = !!Target_group) %>%
    select(TCGA_type, GeneName, Mutant_group, everything())

  plot_df_DE_summary <- bind_rows(plot_df_DE_summary, DE_summary)
}
plot_df_DE_summary

# Histone data ------------------------------------------------------
histone_plot_df <- 
  left_join(
    plot_df %>% select(DepMap_ID, Cell_line_status, Cancer_type_names), 
    histone_marks %>% select(DepMap_ID, H3K4me1, H3K4me2, H3K27ac1K36me0, H3K27ac1K36me1, H3K27ac1K36me2, H3K27ac1K36me3)) %>%
  distinct %>%
  pivot_longer(-c(DepMap_ID, Cell_line_status, Cancer_type_names), names_to = "Marks", values_to = "values") %>%
  mutate(
    Marks = fct_relevel(Marks, "H3K4me1", "H3K4me2", "H3K27ac1K36me0", "H3K27ac1K36me1", "H3K27ac1K36me2", "H3K27ac1K36me3"))

pdf(paste0(DepMap_dir, "DepMap_cell_lines_TCGA_type_combined_KMT2D_muts_histones_boxplot.pdf"), width = 8, height = 15)
  ggplot(histone_plot_df, aes(x = Cell_line_status, y = values))+
    geom_boxplot(
      aes(fill = Cell_line_status), outlier.shape = NA)+
    geom_jitter(
      shape = 16, colour = "black", alpha = 0.5, 
      width = 0.1, height = 0)+
    facet_wrap(~Cancer_type_names + Marks, ncol = 6)+
    scale_colour_manual(values = c("#929598", "#7F3F98"))+
    theme_bw()+
    theme(
      text = element_text(size = 12))
dev.off()

# T-tests
histone_plot_DE <- NULL
for(select_TCGA_type in unique(histone_plot_df$Cancer_type_names)){
  for(select_mark in as.character(unique(histone_plot_df$Marks))){
    # select_TCGA_type <- "Pan_cancer"
    # select_mark <- "H3K4me1"

    print(select_TCGA_type)
    print(select_mark)

    # RNA -----------------------------------------
    histone_plot_df_mark <- histone_plot_df %>%
      filter(
        Cancer_type_names == select_TCGA_type,
        Marks == select_mark,
        !is.na(values))

    if(!any(table(histone_plot_df_mark$Cell_line_status) < 2) & length(table(histone_plot_df_mark$Cell_line_status)) == 2){
      fit_mark <- t.test(values ~ Cell_line_status ,histone_plot_df_mark)

      mark_mean <- NULL
      mark_mean$Control <- histone_plot_df_mark %>% filter(Cell_line_status == "Control") %>% pull(values) %>% mean
      mark_mean$Mutant <- histone_plot_df_mark %>% filter(Cell_line_status == Target_group) %>% pull(values) %>% mean

      mark_N <- NULL
      mark_N$Control <- histone_plot_df_mark %>% filter(Cell_line_status == "Control") %>% nrow()
      mark_N$Mutant <- histone_plot_df_mark %>% filter(Cell_line_status == Target_group) %>% nrow()

    } else {
      fit_mark <- NULL
      fit_mark$p.value <- NA

      mark_mean <- NULL
      mark_mean$Control <- NA
      mark_mean$Mutant <- NA

      mark_N <- NULL
      mark_N$Control <- NA
      mark_N$Mutant <- NA
    }
    DE_summary <- tribble(
    ~mark_control_mean, ~mark_mutant_mean,
    ~mark_control_N, ~mark_mutant_N, ~mark_pval,
    mark_mean$Control, mark_mean$Mutant,
    mark_N$Control, mark_N$Mutant, fit_mark$p.value) %>%
    mutate(
      Mark = select_mark,
      TCGA_type = select_TCGA_type) %>%
    select(TCGA_type, Mark, everything())

  histone_plot_DE <- bind_rows(histone_plot_DE, DE_summary)
  }
}
histone_plot_DE