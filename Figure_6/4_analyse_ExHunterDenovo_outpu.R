#############################################################
# Figure 6
# Analysed Expansion Hunter Denovo results
#############################################################
MSI_KMT2D_LOF <- c("SRR8639191", "SRR8652098", "SRR8652121", "SRR8652136", "SRR8670708", "SRR8652101", "SRR11680468", "SRR11680467", "SRR8670762")
MSS_KMT2D_WT <- c("SRR8670673", "SRR8670700", "SRR8639227")

sample_key <- tibble(
    sample_id = c(MSI_KMT2D_LOF, MSI_KMT2D_WT, MSS_KMT2D_WT),
    cell_line_name = c("CCK81", "LOVO", "MFE319", "NCIH1048", "SW48", "LS180", "RKO", "SKOV3", "SKMEL2", "T84", "NCIH747", "COLO201"),
    cancer_type = c("Colon_Colorectal_Cancer", "Colon_Colorectal_Cancer", "Endometrial_Uterine_Cancer", "Lung_Cancer", "Colon_Colorectal_Cancer", "Colon_Colorectal_Cancer", "Colon_Colorectal_Cancer", "Colon_Colorectal_Cancer", "Colon_Colorectal_Cancer", "Ovarian_Cancer", "Skin_Cancer", "Colon_Colorectal_Cancer", "Colon_Colorectal_Cancer", "Colon_Colorectal_Cancer"),
    group = c(
        rep("MSI_KMT2D_LOF", length(MSI_KMT2D_LOF)),
        rep("MSI_KMT2D_WT", length(MSI_KMT2D_WT)),
        rep("MSS_KMT2D_WT", length(MSS_KMT2D_WT))
    )
)

locus_all <- NULL
for (i in seq(nrow(sample_key))) {
    id <- sample_key$sample_id[i]
    id_group <- sample_key$group[i]
    print(id)

    locus <- read_tsv(paste0(dir, "str_profiles/", id, ".locus.tsv")) %>%
        mutate(
            sample_id = id,
            group = id_group
        )

    locus_all <- bind_rows(locus_all, locus)
}

# Annotate number of AT repeats
# het_str_size:  Finally, the last column provides the estimated repeat expansion size

ATn_locus <- locus_all %>%
    filter(motif == "AT") %>%
    left_join(sample_key) %>%
    mutate(
        sample_id = fct_relevel(sample_id, sample_key$sample_id),
        cell_line_name = fct_relevel(cell_line_name, sample_key$cell_line_name),
    ) %>%
    write_csv(paste0(dir, "/ATn_loci_detected.csv"))

mean_str <- mean(locus_all$het_str_size)
mean_str
# [1] 192.6331

pdf(paste0(dir, "/Est_MS_length.pdf"), height = 5, width = 10)
ATn_locus %>%
    filter(num_anc_irrs > 1) %>%
    ggplot(., aes(x = cell_line_name, y = het_str_size)) +
    geom_boxplot(outlier.shape = NA, aes(fill = group)) +
    geom_jitter(height = 0, width = 0.2, size = 2, shape = 16, aes(colour = group)) +
    geom_hline(yintercept = mean_str) +
    theme_bw() +
    theme(text = element_text(size = 12))
dev.off()