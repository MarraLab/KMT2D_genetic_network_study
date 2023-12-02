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
# Created in Figure_2/3_DepMap_KMT2D_expression.R
Groups_selected <- read_csv(paste0(DepMap_dir, "DepMap_cell_lines_TCGA_type_combined_KMT2D_muts_annotated.csv"))

# Perpare Data -----------------------------------------------------------
# Get all pan cancer lines
KMT2D_groups <- select_cell_lines(input_gene = "KMT2D", data_dir = gretta_data_dir)

Mutant_group <- "KMT2D_LOF"
Control_group <- "Control"
Screen_output <- NULL
for(cancer_type_names in c("Pan_cancer", "BLCA", "NHL_B", "COAD/READ", "ESCA", "HNSC", "LUSC", "SCLC", "STAD", "UCEC")){
  # Context selected
  # cancer_type_names <- "PanCan"
  print(cancer_type_names)

  # Cell lines selection
  if(cancer_type_names == "PanCan"){
    Control_group_avail <- KMT2D_groups %>%
      filter(Group %in% c("Control")) %>%
      pull(DepMap_ID)
    
    Mutant_group_avail <- KMT2D_groups %>%
      filter(Group %in% c("KMT2D_THetDel")) %>%
      pull(DepMap_ID)

  } else {
    Control_group_avail <- Groups_selected %>%
      filter(
        cancer_type_names %in% !!cancer_type_names,
        Cell_line_status == Control_group)

    Mutant_group_avail <- Groups_selected %>%
      filter(
        cancer_type_names %in% !!cancer_type_names,
        Cell_line_status == Mutant_group)
  }

  # In case there are not enough mutant/control cells to perform statistics
  if(nrow(Mutant_group_avail) < 2){
    print(paste0(cancer_type_names, " less than 2 mutant lines"))
    next
  }
  if(nrow(Control_group_avail) < 2){
    print(paste0(cancer_type_names, " less than 2 control lines"))
    next
  }

  # subset 
  select_dep <- dep %>%
    pivot_longer(cols = matches("\\d"), names_to = "GeneNameID", values_to = "DepProb") %>%
    filter(
      !is.na(DepProb)) %>%
    mutate(CellType = case_when(
      DepMap_ID %in% Mutant_group_avail ~ Mutant_group,
      DepMap_ID %in% Control_group_avail ~ Control_group,
      TRUE ~ "Others")) %>%
    filter(CellType != "Others") %>%
    mutate(CellType = fct_relevel(CellType, Control_group, Mutant_group))

  # Begin single gene analysis -----------------------------------------------
  # Begin nested loop
  # 1:length(unique(dep$GeneNameID))
  All_res <- NULL
  All_res <- foreach(each = 1:n_perm, .combine = bind_rows) %dopar% {
  # All_res <- foreach(each = 1:10, .combine = bind_rows) %dopar% {
      if(each == 1){
        cat(paste0("Processing ", each, " of ", n_perm),"\n")
      } else if(each == n_perm){
        cat(paste0("Processing ", each, " of ", n_perm),"\n")
      } else if(each%%100 == 0){
        cat(paste0("Processing ", each, " of ", n_perm),"\n")
      }

      # Create randomly sampled DepProb dataframe
      dummy_geneID <- "A1BG_1"
      df <- select_dep %>% 
        mutate(DepProb_randomize = sample(DepProb, size = length(DepProb), replace = FALSE)) %>%
        filter(GeneNameID == dummy_geneID) %>%
        select(-DepProb) %>%
        rename(DepProb = DepProb_randomize)

      # Begin analysis - analysis was stripped down to only measure the essentials
      if(all(df$DepProb == 0)){
        populate <- rep(0,5)
      # } else if(all(df$DepProb < 0.01)) {
      #   populate <- rep(0,5)

      } else if(all(df$DepProb == 1)){
        populate <- rep(1,5)

      } else {
        # # MWU doesn't handle na or zero's well so
        # # FOR NOW remove zeros.
        # df <- df %>% filter(!is.na(DepProb)) %>%
        #   filter(DepProb != 0)

        stats <- df %>%
          group_by(CellType) %>%
          summarize(Median = median(DepProb, na.rm = TRUE),
                    Mean = mean(DepProb, na.rm = TRUE),
                   .groups = "drop")

        if((any(is.na(stats)) != TRUE) & (nrow(stats) == 2)){

          fit_pval <- wilcox.test(DepProb ~ CellType, df,
                             paired = F,
                             alternative = "two.sided",
                             conf.int = T,
                             na.action = "na.omit")$p.value

        } else if((any(is.na(stats)) == TRUE) & (nrow(stats) == 2)){
          populate <- rep(0,5)
        }

        if((any(is.na(stats)) != TRUE) & (nrow(stats) == 2)){
          populate <- as.numeric(c(unlist(stats)[-c(1,2)], fit_pval))
        } else {
          populate <- rep(0,5)
        }
      }

      tibble(Result = c("Control_median", "Mutant_median", "Control_mean", "Mutant_mean","Pval")) %>%
          mutate(!!sym(geneID) := populate) %>%
          pivot_longer(-Result) %>%
          pivot_wider(names_from = Result, values_from = value) %>%
          rename(GeneNameID = name) %>%
          mutate(
            GeneNameID = paste0("GenePerm_",each),
            Mutant_group = !!Mutant_group,
            Control_group = !!Control_group)
    } # End of for loop

  # Add mutant group name
  output <- All_res %>%
    mutate(
      cancer_type_names = !!cancer_type_names)

  # Combine
  Screen_output <- bind_rows(Screen_output, output)
}

# output
write_csv(Screen_output,
    file = paste0(DepMap_dir, "Screen_",n_perm,"_perms_results.csv"))

end_time <- Sys.time()
cat(paste("Start time", start_time,"\n"))
cat(paste("End time", end_time,"\n"))
cat(end_time - start_time,"\n")
# did you fix 1:length(unique(dep$GeneNameID)) ?
