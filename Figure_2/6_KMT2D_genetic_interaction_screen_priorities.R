#############################################################
# Figure 2
# In silico screen accounting for permutation
#############################################################
library(pacman)
p_load(tidyverse, colorspace, ggrepel)

# Define paths
DepMap_dir <- # Define path to DepMap 20Q1 data, see https://github.com/ytakemon/GRETTA for instructions on downloading GRETTA data.

## Load files -----------------------------------------------------------------
# see https://github.com/ytakemon/GRETTA for instructions on downloading GRETTA data.
# Created in Figure_2/3_DepMap_KMT2D_expression.R
Groups_selected <- read_csv(paste0(DepMap_dir, "DepMap_cell_lines_TCGA_type_combined_KMT2D_muts_annotated.csv"))
KMT2D_groups <- select_cell_lines(input_gene = "KMT2D", data_dir = gretta_data_dir)
Screen_perms <- read_csv(paste0(DepMap_dir, "Screen_10000_perms_results.csv"))
Screen_results <- read_csv(paste0(DepMap_dir, "Insilico_genetic_screening_ALL_results.csv"))

# Append permutated p-values for p-value -------------------------------------------------
# Functions needed 
perm_pvalue <- function(p_value, perm_pvalues) {
  if(is.na(p_value)|p_value == Inf){
    return(NA_integer_)
  }
  x <- table(perm_pvalues <= p_value)
  if(length(x) < 2){
    res <- 1/(x[[1]]+1)
  } else {
    res <- x[["TRUE"]]/(x[["FALSE"]]+ x[["TRUE"]]+1)
  }
  return(res)
}

Combined_screen_res_with_perms <- NULL
for(i in unique(Combined_screen_res$Huffogram_names)){
    # i <- "Pan_cancer"
    screen_res <- Combined_screen_res %>% filter(Huffogram_names == i)

    if(i == "Pan_cancer"){
        perms <- Screen_perms %>% filter(Huffogram_names == "PanCan") %>%
            mutate(
                log2FC_median = log2(Mutant_median/Control_median))
    } else {
        perms <- Screen_perms %>% filter(Huffogram_names == i) %>%
            mutate(
                log2FC_median = log2(Mutant_median/Control_median))
    }

    screen_res_perms <- screen_res %>%
        mutate(
            Pval_perm_pval_correction = map_dbl(Pval, .f = perm_pvalue, perms$Pval))
    
    Combined_screen_res_with_perms <- bind_rows(Combined_screen_res_with_perms, screen_res_perms)
}

Combined_screen_res_with_perms_annot <- Combined_screen_res_with_perms %>% 
    mutate(
        GI_direction = case_when(
            Mutant_median > Control_median ~ "SL",
            Mutant_median < Control_median ~ "AL",
            TRUE ~ "Other"),
        GeneName_GItype = paste0(GeneNames,"_",GI_direction)) %>%
    write_csv(paste0(DepMap_dir, "Screen_results_with_perm_corrected_pvals.csv"))

Combined_screen_res_with_perms_annot <- read_csv(paste0(DepMap_dir, "Screen_results_with_perm_corrected_pvals.csv"))

Combined_screen_res_with_perms_annot_tiers <- Combined_screen_res_with_perms_annot %>%
    mutate(
      Mutant_median = round(Mutant_median, 2),
      Control_median = round(Control_median, 2),
      log2FC_by_median = log2(Mutant_median/Control_median),
      Huffogram_names = fct_relevel(Huffogram_names, "Pan_cancer", "BLCA", "COAD", "ESCA", "HNSC", "LUSC", "NHL_B", "SCLC", "STAD", "UCEC"),
      Signif_tiers = case_when(
          (Pval_perm_pval_correction < 0.01) & (abs(log2FC_by_median) > 2) & (Control_median > 0.5 | Mutant_median > 0.5) ~ "Tier_1",
          (Pval_perm_pval_correction < 0.05) & (abs(log2FC_by_median) > 2) & (Control_median > 0.5 | Mutant_median > 0.5) ~ "Tier_2",
          (Pval_perm_pval_correction < 0.05) & (abs(log2FC_by_median) > 2) ~ "Tier_3",
          TRUE ~ NA_character_),
      to_label = Signif_tiers %in% c("Tier_1", "Tier_2")) %>%
    write_csv(paste0(DepMap_dir, "Screen_results_with_perm_corrected_pvals_tiers_annotated.csv"))

#############################################################
# Figure 2
# Craete upset plot to show overlap of significant GIs
#############################################################
p_load(ComplexHeatmap)
Combined_screen_res_with_perms_annot_tiers <- read_csv(paste0(DepMap_dir, "Screen_results_with_perm_corrected_pvals_tiers_annotated.csv"))

# Selected muts: 
# All muts
plot_list <- list(
  Pan_cancer = Combined_screen_res_with_perms_annot_tiers %>%
    filter(Huffogram_names == "Pan_cancer", !is.na(Signif_tiers)) %>% pull(GeneName_GItype),
  BLCA = Combined_screen_res_with_perms_annot_tiers %>%
    filter(Huffogram_names == "BLCA", !is.na(Signif_tiers)) %>% pull(GeneName_GItype),
  COAD = Combined_screen_res_with_perms_annot_tiers %>%
    filter(Huffogram_names == "COAD", !is.na(Signif_tiers)) %>% pull(GeneName_GItype),
  ESCA = Combined_screen_res_with_perms_annot_tiers %>%
    filter(Huffogram_names == "ESCA", !is.na(Signif_tiers)) %>% pull(GeneName_GItype),
  HNSC = Combined_screen_res_with_perms_annot_tiers %>%
    filter(Huffogram_names == "HNSC", !is.na(Signif_tiers)) %>% pull(GeneName_GItype),
  LUSC = Combined_screen_res_with_perms_annot_tiers %>%
    filter(Huffogram_names == "LUSC", !is.na(Signif_tiers)) %>% pull(GeneName_GItype),
  NHL_B = Combined_screen_res_with_perms_annot_tiers %>%
    filter(Huffogram_names == "NHL_B", !is.na(Signif_tiers)) %>% pull(GeneName_GItype),
  SCLC = Combined_screen_res_with_perms_annot_tiers %>%
    filter(Huffogram_names == "SCLC", !is.na(Signif_tiers)) %>% pull(GeneName_GItype),
  STAD = Combined_screen_res_with_perms_annot_tiers %>%
    filter(Huffogram_names == "STAD", !is.na(Signif_tiers)) %>% pull(GeneName_GItype),
  UCEC = Combined_screen_res_with_perms_annot_tiers %>%
    filter(Huffogram_names == "UCEC", !is.na(Signif_tiers)) %>% pull(GeneName_GItype)
)

plot_matrix <- make_comb_mat(plot_list, mode = "distinct")

ss <- set_size(plot_matrix)
cs <- comb_size(plot_matrix)
cd <- comb_degree(plot_matrix)
ht <- UpSet(plot_matrix,
  set_order =  order(ss, decreasing = T),
  comb_order =  order(cs, cd, decreasing = T),
  pt_size = unit(4, "mm"),
  lwd = 3,
  top_annotation = HeatmapAnnotation(
      "# GIs in group" = anno_barplot(cs,
          ylim = c(0, max(cs)*1.1),
          border = FALSE,
          gp = gpar(fill = "black"),
          height = unit(4, "cm")
      ),
      annotation_name_side = "left",
      annotation_name_rot = 90))

pdf(paste0(DepMap_dir, "UpsetPlot_screen_results_post_perms.pdf"), height = 4, width = 20)
ht <- draw(ht)
od <- column_order(ht)
decorate_annotation("# GIs in group", {
    grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"),
        default.units = "native", just = c("left", "bottom"),
        gp = gpar(fontsize = 10, col = "#404040"), rot = 45)})
dev.off()

#############################################################
# Figure 2
# Craete GI Tiers
#############################################################
DepMap_dir <- "/projects/marralab/ytakemon_prj/DepMap/20Q1/"
Combined_screen_res_with_perms_annot_tiers <- read_csv(paste0(DepMap_dir, "Screen_results_with_perm_corrected_pvals_tiers_annotated.csv"))

# By GI
pdf(paste0(DepMap_dir, "BarPlot_signif_KMT2D_GIs_byGITiers_from_post_perm_list.pdf"), height = 3, width = 8)
Combined_screen_res_with_perms_annot_tiers %>%
  filter(Signif_tiers != "NA") %>%
  count(Huffogram_names, Signif_tiers) %>%
  mutate(
    Signif_tiers = fct_relevel(Signif_tiers, "Tier_3", "Tier_2","Tier_1"),
    Huffogram_names = fct_relevel(Huffogram_names, "Pan_cancer", "BLCA", "COAD", "ESCA", "HNSC", "LUSC", "NHL_B", "SCLC", "STAD", "UCEC")) %>%
  ggplot(aes(x = Huffogram_names, y = log2(n)+1))+
    geom_bar(stat = "identity", aes(fill = Signif_tiers))+
    scale_fill_manual(values = rev(c("#004c6d","#6996b3", "grey")))+
    ylab("# Genetic Interactors") + xlab("")+
    theme_light()+
    theme(
      text = element_text(size = 12))
dev.off()

#############################################################
# Figure 2
# Drug tractability groups
#############################################################
# By new druggability - open target consortium
# Downloaded from OpenTarget consortium
opentarget_df <- read_tsv("OpenTarget/tractability_v23-02.tsv") %>% janitor::clean_names()

compound_tractability <- opentarget_df %>% 
  select(
    symbol, protein_names, top_bucket_sm, category_sm, top_bucket_ab, category_ab, top_bucket_protac, category_protac, top_bucket_othercl, category_othercl) %>%
  mutate(groups = case_when(
    category_sm == "Clinical_Precedence_sm" ~ "Group_1",
    category_ab == "Clinical_Precedence_ab" ~ "Group_1",
    category_protac == "Clinical_Precedence_protac" ~ "Group_1",
    category_othercl == "Clinical_Precedence_othercl" ~ "Group_1",

    category_sm == "Discovery_Precedence_sm" ~ "Group_2",
    category_ab == "Predicted_Tractable_ab_High_confidence_ab" ~ "Group_2",
    category_protac == "Literature_Precedence_protac" ~ "Group_2",
  
    category_sm == "Predicted_Tractable_sm" ~ "Group_3",
    category_ab == "Predicted_Tractable_ab_Medium_to_low_confidence_ab" ~ "Group_3",
    category_protac == "Discovery_Opportunity_protac" ~ "Group_3",
    TRUE ~ NA_character_
  ))

GI_Drug_Tiers <- Combined_screen_res_with_perms_annot_tiers %>%
  filter(Signif_tiers != "NA") %>%
  mutate(Tiers = case_when(
    GeneNames %in% (compound_tractability %>% filter(groups == "Group_1") %>% pull(symbol)) ~ "Tier_1",
    GeneNames %in% (compound_tractability %>% filter(groups == "Group_2") %>% pull(symbol)) ~ "Tier_2",
    GeneNames %in% (compound_tractability %>% filter(groups == "Group_3") %>% pull(symbol)) ~ "Tier_3",
    TRUE ~ "Tier_3"),
    Tiers = fct_relevel(Tiers, "Tier_1","Tier_2")) %>%
  filter(!is.na(Tiers)) %>%
  rename(GI_tiers = Signif_tiers, Drug_tiers = Tiers) %>%
  write_csv(paste0(DepMap_dir, "Signif_KMT2D_GIs_byDrugTractabilityTiers_OpenTarget_from_post_perm_list.csv"))

pdf(paste0(DepMap_dir, "BarPlot_signif_KMT2D_GIs_byDrugTractabilityTiers_OpenTarget_from_post_perm_list.pdf"), height = 3, width = 8)
GI_Drug_Tiers %>%
  count(Huffogram_names, Drug_tiers) %>%
  mutate(
    Drug_tiers = fct_relevel(Drug_tiers, "Tier_3", "Tier_2"),
    Huffogram_names = fct_relevel(Huffogram_names, "Pan_cancer", "BLCA", "COAD", "ESCA", "HNSC", "LUSC", "NHL_B", "SCLC", "STAD", "UCEC")
  ) %>%
  ggplot(aes(x = Huffogram_names, y = log2(n)+1))+
    geom_bar(stat = "identity", aes(fill = Drug_tiers))+
    scale_fill_manual(values = rev(c("#56096d", "#a573b2", "grey")))+
    ylab("# Tractable Targets") + xlab("")+
    theme_light()+
    theme(text = element_text(size = 12)) 
dev.off()

#############################################################
# Figure 2
# Upset plot comparing overlap between GI tiers and Drug tractability groups
#############################################################
# Upset plot 
p_load(ComplexHeatmap)
# Redo plot again to make signif GIs more obvious -----------------------------------------------------------------
Combined_screen_res_with_perms_annot_tiers <- read_csv(paste0(DepMap_dir, "Screen_results_with_perm_corrected_pvals_tiers_annotated.csv"))
GI_Drug_Tiers <- read_csv(paste0(DepMap_dir, "Signif_KMT2D_GIs_byDrugTractabilityTiers_OpenTarget_from_post_perm_list.csv")) %>%
  filter(
    GI_tiers %in% c("Tier_1"),
    Drug_tiers %in% c("Tier_1", "Tier_2")) %>%
  mutate(Best = TRUE) %>%
  select(Huffogram_names, GeneName_GItype, Best)

plot_df <- Combined_screen_res_with_perms_annot_tiers %>%
  left_join(GI_Drug_Tiers) %>%
  mutate(
    log2FC_by_median = ifelse(log2FC_by_median %in% c(Inf, -Inf), 1, log2FC_by_median),
    To_label = case_when(
      Best ~ TRUE,
      TRUE ~ FALSE),
    new_score = -log10(Pval_perm_pval_correction) * (log2FC_by_median),
    new_score = ifelse(is.na(new_score),0, new_score)
    ) %>% 
  group_by(Huffogram_names) %>%
  arrange(new_score) %>% 
  mutate(rank = seq(18333)) %>% ungroup

pdf(paste0(DepMap_dir, "Screen_results_with_perm_corrected_pvals_improved_again.pdf"), width = 7, height = 8)
plot_df %>%
    mutate(
      Huffogram_names = fct_relevel(Huffogram_names, "Pan_cancer", "BLCA", "NHL_B", "COAD", "ESCA", "HNSC", "LUSC",  "SCLC", "STAD", "UCEC"),
      To_colour = ifelse(!is.na(Signif_tiers),"signif","ns")) %>% 
    # filter(Huffogram_names %in% c("Pan_cancer", "COAD")) %>%
    ggplot(., aes(y = new_score, x = rank))+
        geom_point(shape = 16, aes(colour = To_colour))+
        ylab("GI score")+
        xlab("Rank")+
        # facet_wrap(~Huffogram_names)+
        facet_wrap(~Huffogram_names, ncol = 3, nrow = 4)+
        theme(text = element_text(size = 10))+
        theme_light()+
        guides(size = guide_bins(show.limits = TRUE)) +
        geom_label_repel(
            aes(label = ifelse(To_label, GeneNames, "")),
            size = 3,
            direction = "both",
            # alpha = 0.75, 
            fontface = 'bold', 
            color = 'black',
            box.padding   = 0.5,
            point.padding = 0,
            force = 1,
            max.overlaps = Inf,
            label.size = NA, 
            fill = NA,
            min.segment.length = unit(0, 'lines'),
            segment.color = 'grey50')
dev.off()