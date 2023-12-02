# Figure 1A-B, Mapping KMT2D essentiality network 
library(pacman)
p_load(tidyverse, GRETTA)

#############################################################
# Figure 1A Mapping KMT2D coessential networks
#############################################################

# Define paths, see https://github.com/ytakemon/GRETTA for instructions on downloading GRETTA data.
output_dir <- # Path to output directory
data_dir <- # Path to GRETTA 20Q1 data directory, 

# Calculate coefficients
essentiality_network_df <- coessential_map(
    input_genes = "KMT2D",
    core_num = 20,
    output_dir = output_dir,
    data_dir = data_dir,
    filename = "KMT2D_essentiality_network.csv"
)

# Determine inflection points
inflection_df <- get_inflection_points(essentiality_network_df)
annotated_df <- annotate_coess(essentiality_network_df, inflection_df)

# Plot ranked list
plot_coess(
    result_df = annotated_df,
    inflection_df = inflection_df
)

#############################################################
# Figure 1B Enrichment analysis
#############################################################
# Over-represetation analysis 
# Generage GO enrichment plot
p_unload("all")
library(pacman)
p_load(clusterProfiler, org.Hs.eg.db, DOSE, enrichplot)

# Data:
# Using `annotated_df` from above

# Grab co-essential genes
df_coess <- annotated_df %>%
    filter(
        Candidate == TRUE,
        estimate > 0
    )

# Grab anti-essential genes
df_antiess <- annotated_df %>%
    filter(
        Candidate == T,
        estimate < 0
    )

# Grab ENTREZ ID of the universe
All_IDs <- bitr(unique(annotated_df$GeneName_B),
    fromType = "SYMBOL",
    toType = c("ENSEMBL", "ENTREZID"),
    org.Hs.eg.db::org.Hs.eg.db,
    drop = FALSE
)
# Get ENTREZ IDs of co/anti-essential genes
coess_IDs <- bitr(unique(df_coess$GeneName_B),
    fromType = "SYMBOL",
    toType = c("ENSEMBL", "ENTREZID"),
    org.Hs.eg.db::org.Hs.eg.db,
    drop = FALSE
)
antiess_IDs <- bitr(unique(df_antiess$GeneName_B),
    fromType = "SYMBOL",
    toType = c("ENSEMBL", "ENTREZID"),
    org.Hs.eg.db::org.Hs.eg.db,
    drop = FALSE
)


# GO over-representation analysis based on Boyles et al. 2004, p-value calculated by hyper geometric distribution
ego_coess <- enrichGO(
    gene = unique(coess_IDs$ENTREZID),
    universe = unique(All_IDs$ENTREZID),
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 1,
    minGSSize = 5,
    maxGSSize = nrow(coess_IDs),
    readable = TRUE
)

ego_antiess <- enrichGO(
    gene = unique(antiess_IDs$ENTREZID),
    universe = unique(All_IDs$ENTREZID),
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 1,
    minGSSize = 5,
    maxGSSize = nrow(antiess_IDs),
    readable = TRUE
)

# Collect results
all_egos <-
    bind_rows(
        ego_coess@result %>% tibble() %>%
            mutate(Huffogram_names = i, Essentiality_type = "co_ess")
    ) %>%
    bind_rows(.,
        ego_antiess@result %>% tibble() %>%
            mutate(Huffogram_names = i, Essentiality_type = "anti_ess")
    )

#############################################################
# Figure 1B Enrichment analysis - calculate Jaccard index & create plot
#############################################################
p_unload("all")
library(pacman)
p_load(tidyverse, ggrepel, doMC, clusterProfiler, org.Hs.eg.db, DOSE, enrichplot, ComplexHeatmap, circlize)
registerDoMC(30)

# Data uses `all_egos` from above

# Select significant GO terms
pancan_egos <- all_egos %>%
    filter(
        p.adjust < 0.05, qvalue < 0.05
    )

pancan_egos_coess <- pancan_egos %>% filter(Essentiality_type == "co_ess")
pancan_egos_antiess <- pancan_egos %>% filter(Essentiality_type == "anti_ess")


# Calculate jaccard index --------------------------------------------------------------------
# Define function
jaccard <- function(a, b) {
  intersection <- length(intersect(a, b))
  union <- length(a) + length(b) - intersection
  return(intersection / union)
}

jaccard_df <- expand_grid(GO1 = pancan_egos_coess$geneID, GO2 = pancan_egos_coess$geneID) %>%
  mutate(jaccard_index = NA) %>% distinct

All_res <- NULL
All_res <- foreach(i = 1:nrow(jaccard_df), .combine = bind_rows) %dopar% {
  if (i == 1) {
    cat(paste0("Processing ", i, " of ", nrow(jaccard_df)), "\n")
  } else if (i == nrow(jaccard_df)) {
    cat(paste0("Processing ", i, " of ", nrow(jaccard_df)), "\n")
  } else if (i %% 1000 == 0) {
    cat(paste0("Processing ", i, " of ", nrow(jaccard_df)), "\n")
  }

  GO1 <- jaccard_df$GO1[i]
  GO2 <- jaccard_df$GO2[i]

  GO1_genes <- pancan_egos_coess %>%
    filter(geneID == GO1) %>%
    pull(geneID) %>%
    str_split_fixed(., "/", Inf) %>%
    as.vector()
  GO2_genes <- pancan_egos_coess %>%
    filter(geneID == GO2) %>%
    pull(geneID) %>%
    str_split_fixed(., "/", Inf) %>%
    as.vector()

  # jaccard_df$jaccard_index[i] <- jaccard(GO1_genes, GO2_genes)
  tibble(GO1 = GO1, GO2 = GO2, jaccard_index = jaccard(GO1_genes, GO2_genes))
}

jaccard_mat <- jaccard_df %>%
  pivot_wider(names_from = GO2, values_from = jaccard_index) %>%
  dplyr::select(-GO1) %>%
  as.matrix()
rownames(jaccard_mat) <- jaccard_df$GO1 %>% unique()

# Calculate optimal clustering
p_load(cluster, factoextra) # clustering algorithms & visualization
# Gap stat method
df <- jaccard_mat
gap_stat <- clusGap(df,
  FUN = kmeans, nstart = 25,
  K.max = 15, B = 10000
)
n_clust <- with(gap_stat, maxSE(Tab[, "gap"], Tab[, "SE.sim"]))

# plot heatmap
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

pdf(paste0(output_dir, "Coess_GO_enrichment_clusteredpdf"), width = 7, height = 6)
draw(hmap)
dev.off()

# Annotate df for supplemental table
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

# output heatmap annotations
res <- jaccard_mat[r_order, c_order] %>%
  as_tibble(rownames = "GO1") %>%
  mutate(cluster_id = cluster_annot) %>%
  left_join(., pancan_egos_coess %>% dplyr::select(GO1 = "geneID", Description, ID)) %>%
  dplyr::select(cluster_id, ID, Description, GO1) %>%
  write_csv(paste0(output_dir, "Coess_GO_enrichment_clustered.csv"))
