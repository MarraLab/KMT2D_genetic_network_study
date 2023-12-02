#############################################################
# Figure 1D, Overlap between KMT2D co/anti-essential genes and ChIP-MS
#############################################################
# setup
library(pacman)
p_load(tidyverse, clusterProfiler, org.Hs.eg.db, DOSE, enrichplot)
processed_dir <- # Path to fragpipe procesed data

# Data, see Figure_1/KMT2D_ChIP_MS_preprocessing.R and KMT2D_essentiality_network.R
KMT2D_ChIPMS <- read_csv(paste0(processed_dir, "KMT2D_res_positive_in_at_least_2_KO_cleaned.csv")) %>% arrange(gene)
annotated_df # generated in KMT2D_essentiality_network.R

# creat supp table
gene_effs_coess <- All_cor_annotated %>%
  filter(
    Candidate_inflection,
    estimate > 0
  )

gene_effs_antiess <- All_cor_annotated %>%
  filter(
    Candidate_inflection,
    estimate < 0
  )

KMT2D_ChIPMS %>%
  select(ChIP_MS_gene = gene) %>%
  mutate(In_essentiality_network = case_when(
    ChIP_MS_gene %in% gene_effs_coess$GeneName_B ~ "1_co_ess",
    ChIP_MS_gene %in% gene_effs_antiess$GeneName_B ~ "2_anti_ess",
    TRUE ~ NA_character_
  )) %>%
  arrange(In_essentiality_network) %>%
    write_csv(paste0(processed_dir, "/UpsetPlot_KMT2D_ChIPMS_Coess_compare.csv"))

# Upset plot
p_load(ComplexHeatmap)
# Selected muts:
# All muts
plot_list <- list(
  ChIP_MS = KMT2D_ChIPMS$gene,
  Co_ess = All_cor_annotated %>%
    filter(
      Huffogram_names == "Pan_cancer",
      colour_inflection == "Red",
      estimate > 0
    ) %>%
    pull(GeneName_B),
  Anti_ess = All_cor_annotated %>%
    filter(
      Huffogram_names == "Pan_cancer",
      colour_inflection == "Red",
      estimate < 0
    ) %>%
    pull(GeneName_B)
)

plot_matrix <- make_comb_mat(plot_list, mode = "distinct")

ss <- set_size(plot_matrix)
cs <- comb_size(plot_matrix) 
cd <- comb_degree(plot_matrix)
ht <- UpSet(plot_matrix,
  set_order = order(ss, decreasing = T),
  comb_order = order(cs, cd, decreasing = T),
  pt_size = unit(4, "mm"),
  lwd = 3,
  top_annotation = HeatmapAnnotation(
    "# genes in group" = anno_barplot(cs,
      ylim = c(0, max(cs) * 1.1),
      border = FALSE,
      gp = gpar(fill = "black"),
      height = unit(4, "cm")
    ),
    annotation_name_side = "left",
    annotation_name_rot = 90
  )
)

pdf(paste0(processed_dir, "/UpsetPlot_KMT2D_ChIPMS_Coess_compare.pdf"), height = 6, width = 5)
ht <- draw(ht)
od <- column_order(ht)
decorate_annotation("# genes in group", {
  grid.text(cs[od],
    x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"),
    default.units = "native", just = c("left", "bottom"),
    gp = gpar(fontsize = 10, col = "#404040"), rot = 45
  )
})
dev.off()