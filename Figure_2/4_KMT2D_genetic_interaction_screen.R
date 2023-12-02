#############################################################
# Figure 2
# In silico screen - Pan Cancer
#############################################################
library(pacman)
p_load(tidyverse, GRETTA)

# Define paths
DepMap_dir <- # Define path to DepMap 20Q1 data, see https://github.com/ytakemon/GRETTA for instructions on downloading GRETTA data.
  
## Load files -----------------------------------------------------------------
# see https://github.com/ytakemon/GRETTA for instructions on downloading GRETTA data.

# Get all pan cancer lines 
KMT2D_groups <- select_cell_lines(input_gene = "KMT2D", data_dir = gretta_data_dir)

KMT2D_mutant_id <- KMT2D_groups %>%
  filter(Group %in% c("KMT2D_THetDel")) %>%
  pull(DepMap_ID)
KMT2D_control_id <- KMT2D_groups %>%
  filter(Group %in% c("Control")) %>%
  pull(DepMap_ID)

# Screen
# This can take several hours depending on number of lines/cores used.
pancancer_screen_results <- GI_screen(
  control_id = KMT2D_control_id,
  mutant_id = KMT2D_mutant_id,
  core_num = 30, # depends on how many cores you have
  output_dir = DepMap_dir, # Will save your results here as well as in the variable
  data_dir = DepMap_dir,
  test = FALSE) %>%
  mutate(
    cancer_type_name = "PanCancer",
    log2FC_by_median = log2(Mutant_median / Control_median),
    log2FC_by_mean = log2(Mutant_mean / Control_mean),
    padj = p.adjust(Pval, "BH", length(Pval)),
    Interaction_score = -log10(Pval) * sin(Mutant_median - Control_median)
  ) %>%
    left_join(., dep_annot %>% select(GeneNameID, GeneNames),
      by = "GeneNameID"
    ) %>%
    select(cancer_type_name, Control_group, Mutant_group, GeneNameID, GeneNames, Control_median:Pval, log2FC_by_median, log2FC_by_mean, everything())

#############################################################
# Figure 2
# In silico screen - Cancer type specific
#############################################################

## Load files ---------------------------------------------------------
# Created in Figure_2/3_DepMap_KMT2D_expression.R
Groups_selected <- read_csv(paste0(DepMap_dir, "DepMap_cell_lines_TCGA_type_combined_KMT2D_muts_annotated.csv"))

# Screen -----------------------------------------------------------
All_cancer_type_screen_output <- NULL
for (cancer_type_name in unique(Groups_selected$cancer_type_name)) {
  print(cancer_type_name)

  # Select cell lines
  Control_group_avail <- Groups_selected %>%
    filter(
      cancer_type_name %in% !!cancer_type_name,
      Cell_line_status == "Control"
    )
  Mutant_groups_avail <- Groups_selected %>%
    filter(
      cancer_type_name %in% !!cancer_type_name,
      Cell_line_status %in% c("KMT2D_HomDel", "KMT2D_THetDel", "KMT2D_HetDel")
    )

  # Screen
  screen_res <- GI_screen(
    control_id = Control_group_avail$DepMap_ID,
    mutant_id = Mutant_groups_avail$DepMap_ID,
    core_num = 30, # depends on how many cores you have
    output_dir = DepMap_dir, # Will save your results here as well as in the variable
    data_dir = DepMap_dir,
    test = FALSE
  )

  # Add cancer type name
  output <- All_res %>%
    mutate(
      cancer_type_name = !!cancer_type_name,
      log2FC_by_median = log2(Mutant_median / Control_median),
      log2FC_by_mean = log2(Mutant_mean / Control_mean),
      padj = p.adjust(Pval, "BH", length(Pval)),
      Interaction_score = -log10(Pval) * sin(Mutant_median - Control_median)
    ) %>%
    left_join(., dep_annot %>% select(GeneNameID, GeneNames),
      by = "GeneNameID"
    ) %>%
    select(cancer_type_name, Control_group, Mutant_group, GeneNameID, GeneNames, Control_median:Pval, log2FC_by_median, log2FC_by_mean, everything())

  # Combine
  All_cancer_type_screen_output <- bind_rows(screen_res, output)
}

bind_rows(pancancer_screen_results, All_cancer_type_screen_output) %>%
  write_csv(.,
    file = paste0(DepMap_dir, "Insilico_genetic_screening_ALL_results.csv")
  )