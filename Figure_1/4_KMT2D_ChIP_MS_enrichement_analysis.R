#############################################################
# Figure 1C, ChIP-MS enrichment analysis
#############################################################
# setup
library(pacman)
p_load(tidyverse, clusterProfiler, org.Hs.eg.db, DOSE, enrichplot)
processed_dir <- # Path to fragpipe procesed data
 
# Data, see Figure_1/KMT2D_ChIP_MS_preprocessing.R
KMT2D_ChIPMS <- read_csv(paste0(processed_dir, "KMT2D_res_positive_in_at_least_2_KO_cleaned.csv")) %>% arrange(gene)

# Get ENTREZ ID for entire set
convert_names_all <- bitr(unique(KMT2D_ChIPMS$gene),
  fromType = "SYMBOL",
  toType = c("ENSEMBL", "ENTREZID"),
  org.Hs.eg.db::org.Hs.eg.db,
  drop = FALSE
)

# GO over-representation analysis based on Boyles et al. 2004, p-value calculated by hyper geometric distribution
# Biological process
ego_BP <- enrichGO(
  gene = unique(convert_names_all$ENTREZID),
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  minGSSize = 5,
  maxGSSize = 100,
  readable = TRUE
)

ego_BP@result %>%
  as_tibble() %>%
  filter(qvalue < 0.05) %>%
  write_csv(paste0(processed_dir, "/KMT2D_ChIPMS_GO_BP_enriched.csv"))

#############################################################
# Figure 1C, ChIP-MS enrichment analysis, calculate Jaccard index
#############################################################
p_unload("all")
library(pacman)
p_load(tidyverse, ggrepel, doMC, clusterProfiler, org.Hs.eg.db, DOSE, enrichplot, ComplexHeatmap, circlize)
registerDoMC(30)

PP_ego <- read_csv(paste0(processed_dir, "/KMT2D_ChIPMS_GO_BP_enriched.csv"))

# Functions
jaccard <- function(a, b) {
  intersection <- length(intersect(a, b))
  union <- length(a) + length(b) - intersection
  return(intersection / union)
}

jaccard_df <- expand_grid(GO1 = PP_ego$ID, GO2 = PP_ego$ID) %>%
  mutate(jaccard_index = NA) %>%
  distinct()

All_res <- NULL
All_res <- foreach(i = 1:nrow(jaccard_df), .combine = bind_rows) %dopar% {
  # All_res <- foreach(i = 1:10, .combine = bind_rows) %dopar% {
  # print count
  if (i == 1) {
    cat(paste0("Processing ", i, " of ", nrow(jaccard_df)), "\n")
  } else if (i == nrow(jaccard_df)) {
    cat(paste0("Processing ", i, " of ", nrow(jaccard_df)), "\n")
  } else if (i %% 1000 == 0) {
    cat(paste0("Processing ", i, " of ", nrow(jaccard_df)), "\n")
  }

  GO1 <- jaccard_df$GO1[i]
  GO2 <- jaccard_df$GO2[i]

  GO1_genes <- PP_ego %>%
    filter(ID == GO1) %>%
    pull(geneID) %>%
    str_split_fixed(., "/", Inf) %>%
    as.vector()
  GO2_genes <- PP_ego %>%
    filter(ID == GO2) %>%
    pull(geneID) %>%
    str_split_fixed(., "/", Inf) %>%
    as.vector()

  # jaccard_df$jaccard_index[i] <- jaccard(GO1_genes, GO2_genes)
  tibble(GO1 = GO1, GO2 = GO2, jaccard_index = jaccard(GO1_genes, GO2_genes))
}

write_csv(All_res, paste0(processed_dir, "/KMT2D_ChIPMS_GO_BP_enriched_full_jaccard.csv"))

jaccard_df <- read_csv(paste0(processed_dir, "/KMT2D_ChIPMS_GO_BP_enriched_full_jaccard.csv"))

jaccard_mat <- jaccard_df %>%
  pivot_wider(names_from = GO2, values_from = jaccard_index) %>%
  dplyr::select(-GO1) %>%
  as.matrix()
rownames(jaccard_mat) <- jaccard_df$GO1 %>% unique()

# Calculate optimal clustering
p_load(cluster, factoextra) # clustering algorithms & visualization
# Gap stat method
gap_stat <- clusGap(df,
  FUN = kmeans, nstart = 25,
  K.max = 15, B = 10000
)
n_clust <- 8

#############################################################
# Figure 1C, ChIP-MS enrichment analysis, plot heatmap
#############################################################
col_fun <- colorRamp2(c(0, 1), c("white", "red"))
hmap <- Heatmap(
  jaccard_mat,
  col = col_fun,
  clustering_distance_rows = "euclidean",
  clustering_method_rows = "complete",
  row_split = n_clust,
  column_split = n_clust,
  top_annotation = HeatmapAnnotation(foo = anno_block(
    gp = gpar(fill = 2),
    labels = c(1:n_clust),
    labels_gp = gpar(col = "white", fontsize = 10)
  )),
  left_annotation = rowAnnotation(foo = anno_block(
    gp = gpar(fill = 2),
    labels = c(1:n_clust),
    labels_gp = gpar(col = "white", fontsize = 10)
  )),
)

pdf(paste0(processed_dir, "/KMT2D_ChIPMS_GO_BP_enriched_full_jaccard_reordered.pdf"), width = 7, height = 6)
draw(hmap)
dev.off()

hmap2 <- draw(hmap)
r_order <- row_order(hmap2) %>% unlist()
c_order <- column_order(hmap2) %>% unlist()
# rownames(jaccard_mat)[row_order(hmap2)]
cluster_annot <- NULL
for (i in 1:length(row_order(hmap2))) {
  # i <- 1
  len <- row_order(hmap2)[[i]] %>% length()
  cluster_annot <- c(cluster_annot, rep(i, len))
}

#############################################################
# Figure 1C, ChIP-MS enrichment analysis, output heatmap annotations
#############################################################
res <- jaccard_mat[r_order, c_order] %>%
  as_tibble(rownames = "GO1") %>%
  mutate(cluster_id = cluster_annot) %>%
  left_join(., PP_ego %>% dplyr::select(GO1 = "ID", Description, geneID)) %>%
  dplyr::select(cluster_id, GO1, Description, geneID) %>%
  write_csv(paste0(processed_dir, "/KMT2D_ChIPMS_GO_BP_enriched_full_jaccard_reordered.csv"))

