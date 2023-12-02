#############################################################
# Figure 3
# KMT2D - WRN interaction
#############################################################
library(pacman)
p_load(tidyverse, colorspace, ggrepel)

# Define paths
DepMap_dir <- # Define path to DepMap 20Q1 data, see https://github.com/ytakemon/GRETTA for instructions on downloading GRETTA data.

## Load files -----------------------------------------------------------------
# see https://github.com/ytakemon/GRETTA for instructions on downloading GRETTA data.
load(paste0(DepMap_dir, "dep.rda"))
# Created in Figure_2/3_DepMap_KMT2D_expression.R
Groups_selected <- read_csv(paste0(DepMap_dir, "DepMap_cell_lines_TCGA_type_combined_KMT2D_muts_annotated.csv"))
KMT2D_groups <- select_cell_lines(input_gene = "KMT2D", data_dir = gretta_data_dir)
# Created in Figure_2/4_KMT2D_genetic_interaction_screen.R
Screen_results <- read_csv(paste0(DepMap_dir, "Insilico_genetic_screening_ALL_results.csv"))
# Created in Figure_2/5_KMT2D_genetic_interaction_screen_permutation.R
Screen_perms <- read_csv(paste0(DepMap_dir, "Screen_10000_perms_results.csv"))
# Downloaded from Ghandi et al 2019
MSI_calls <- read_xlsx(paste0(DepMap_dir, "Ghandi_2019_MSI_SuppTable7.xlsx"), sheet = "MSI calls") %>%
  dplyr::rename(
    DepMap_ID = depMapID,
    MS_status = CCLE.MSI.call
  ) %>%
  select(DepMap_ID, MS_status)

#############################################################
# Figure 3
# Compare WRN KO effects on KMT2D WT/ LOF lines
#############################################################
Pan_can_muts_annotated_MS <- All_KMT2D_muts_annotated_MS %>%
  filter(Group == "KMT2D_THetDel" | (Group == "Control")) %>%
  left_join(., MSI_calls, by = "DepMap_ID") %>%
  left_join(., dep %>% select(DepMap_ID, convert_GeneName_to_GeneNameID("WRN")), by = "DepMap_ID") %>%
  dplyr::rename(WRN_dep = convert_GeneName_to_GeneNameID("WRN"),
         Cell_line_status = "Group") %>%
  mutate(MS_status = case_when(
    is.na(MS_status) ~ "NA",
    TRUE ~ MS_status),
    MS_status = fct_relevel(MS_status, "inferred-MSI","inferred-MSS","undetermined"))

# plot
WRN_dep_plot_by_MSI <- ggplot(Pan_can_muts_annotated_MS, aes(x = Cell_line_status, y = WRN_dep, fill = MS_status)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(preserve = "single"))+
  geom_point(colour = "black", alpha = 0.5, 
    position = position_jitterdodge(jitter.width = 0.1))+
  scale_fill_aaas()+
  ylim(c(0,1.2))+
  ylab("WRN lethality probability")+
  theme_light()+
  theme(legend.position="bottom", text = element_text(size = 12))

pdf(paste0(DepMap_dir, "BoxPlot_WRN_lethality_dep_by_KMT2D_THet_MS_status.pdf"), height = 6, width = 7)
print(WRN_dep_plot_by_MSI)
dev.off()

# ANOVA - Tukey
test <- Pan_can_muts_annotated_MS %>%
  mutate(Mut_status_by_MSI = paste0(Cell_line_status, "_", MS_status))

fit <- aov(WRN_dep ~ Mut_status_by_MSI, test)

summary(fit)
TukeyHSD(fit) %>%
  tidy() %>%
  filter(adj.p.value < 0.05)

#############################################################
# Figure 3
# Compare KMT2D KO effects on WRN WT/ LOF lines
#############################################################
WRN_groups <- select_cell_lines(input_gene = "WRN", data_dir = gretta_data_dir)

Pan_can_WRN_muts_annotated_MS <- WRN_groups  %>%
  filter(Group == "WRN_HetDel" | Group == "Control") %>%
  left_join(., MSI_calls, by = "DepMap_ID") %>%
  left_join(., dep %>% select(DepMap_ID, convert_GeneName_to_GeneNameID("KMT2D")), by = "DepMap_ID") %>%
  dplyr::rename(KMT2D_dep = convert_GeneName_to_GeneNameID("KMT2D"),
         Cell_line_status = "Group") %>%
  mutate(MS_status = case_when(
    is.na(MS_status) ~ "NA",
    TRUE ~ MS_status),
    MS_status = fct_relevel(MS_status, "inferred-MSI","inferred-MSS","undetermined"))

Pan_can_WRN_muts_annotated_MS %>%
  count(Cell_line_status, MS_status)
# # A tibble: 7 × 3
#   Cell_line_status MS_status        n
#   <chr>            <fct>        <int>
# 1 Control          inferred-MSI    26
# 2 Control          inferred-MSS   330
# 3 Control          undetermined     3
# 4 Control          NA             135
# 5 WRN_mut_3        inferred-MSI     6
# 6 WRN_mut_3        inferred-MSS   132
# 7 WRN_mut_3        NA              38

# plot
KMT2D_dep_plot_by_MSI <- ggplot(Pan_can_WRN_muts_annotated_MS, aes(x = Cell_line_status, y = KMT2D_dep, fill = MS_status)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(preserve = "single"))+
  geom_point(colour = "black", alpha = 0.5, 
    position = position_jitterdodge(jitter.width = 0.1))+
  scale_fill_aaas()+
  ylim(c(0,1.2))+
  ylab("KMT2D lethality probability")+
  theme_light()+
  theme(legend.position="bottom", text = element_text(size = 12))

pdf(paste0(DepMap_dir, "BoxPlot_KMT2D_lethality_dep_by_WRN_Het_MS_status.pdf"), height = 6, width = 7)
print(KMT2D_dep_plot_by_MSI)
dev.off()

# Plot stats 
test <- Pan_can_WRN_muts_annotated_MS %>%
  mutate(Mut_status_by_MSI = paste0(Cell_line_status,"_",MS_status))

fit <- aov(KMT2D_dep ~ Mut_status_by_MSI, test)
summary(fit)
TukeyHSD(fit) %>%
  tidy() %>%
  filter(adj.p.value < 0.05)