#############################################################
# Figure 1C, ChIP-MS Fragpipe ouput processing
#############################################################
## Setting up the environment
library('tidyverse')
library('ggplot2')
library('RColorBrewer')

baseRepository <- # FragPipe output data

# Process the data into individual sets per replicate from the raw search output.
samples = tibble('name' = c(rep('control',4),rep('kmt2d',4)),
                 'replicate' = c(rep(seq(1,4,1),2)))

samples$psmPath = file.path(baseRepository,
                            'fragpipeAnalysis',
                            paste(samples$name,'_',samples$replicate,sep=''),
                            'psm.tsv')
samples$quantPath = file.path(baseRepository,
                              'fragpipeAnalysis/quants',
                              paste(samples$name,'_',samples$replicate,'_Matrix.txt',sep=''))
  
all(file.exists(samples$psmPath))
all(file.exists(samples$quantPath))

for (i in 1:nrow(samples)){
  quantTemp = read_tsv(samples$quantPath[i]) %>%
    dplyr::rename(scan = MS2ScanNumber,
                  area = ParentPeakArea) %>%
    dplyr::select(scan, area)

  pepTemp = read_tsv(samples$psmPath[i]) %>%
    dplyr::rename(accession = `Protein ID`,
                  gene = `Gene`,
                  sequence = `Peptide`,
                  mods = `Assigned Modifications`,
                  spectrum = Spectrum) %>%
    dplyr::select(accession, gene, sequence, mods, spectrum) %>%
    dplyr::filter(!grepl('^KRT', gene), !grepl('contam', accession)) %>%
    dplyr::mutate(scan = as.numeric(sub('.*[1-4]\\.(.*)\\.[0-9]+\\.[0-9]+','\\1',spectrum))) %>%
    dplyr::mutate(label = ifelse(grepl('R\\(10', mods), 'heavy',
                                 ifelse(grepl('K\\(8', mods), 'heavy', 'light'))) %>%
    dplyr::mutate(dataset = paste(samples[i,1],'_',samples[i,2],sep=''), numPsm = 1) %>%
    dplyr::select(scan, dataset, accession, gene, sequence, label, numPsm) %>%
    dplyr::mutate(label = factor(label, levels = c('light','heavy'))) %>%
    dplyr::left_join(quantTemp) %>%
    dplyr::group_by(dataset, accession, gene, sequence, label) %>%
    dplyr::summarise(numPsm = sum(numPsm, na.rm = TRUE), area = median(area, na.rm = TRUE)) %>%
    dplyr::mutate(numPeps = 1) %>%
    dplyr::group_by(dataset, accession, gene, label) %>%
    dplyr::summarise(numPsm = sum(numPsm, na.rm = TRUE), numPeps = sum(numPeps, na.rm = TRUE), area = median(area, na.rm = TRUE)) %>%
    tidyr::pivot_wider(id_cols = dataset:gene, names_from = 'label', values_from = c('numPsm','numPeps','area')) %>%
    dplyr::mutate(area_light = ifelse(is.na(area_light), 0, area_light)) %>%
    dplyr::mutate(area_heavy = ifelse(is.na(area_heavy), 0, area_heavy)) %>%
    dplyr::mutate(area_light_log2 = log2(area_light + 1000), area_heavy_log2 = log2(area_heavy + 1000)) %>%
    dplyr::mutate(log2HL = area_heavy_log2 - area_light_log2) %>%
    dplyr::rowwise() %>%
    dplyr::filter(sum(c(numPeps_light,numPeps_heavy), na.rm = TRUE) >= 2) %>%
    dplyr::arrange(-log2HL)

  saveRDS(pepTemp, paste(baseRepository,
                         '/fragpipeAnalysis/processedData/',
                         samples[i,1],
                         '_',
                         samples[i,2],
                         '_processed.rds',
                         sep = ''))

  write.table(pepTemp, paste(baseRepository,
                             '/fragpipeAnalysis/processedData/',
                             samples[i,1],
                             '_',
                             samples[i,2],
                             '_processed.csv',
                             sep = ''),
              quote = FALSE, col.names = TRUE, row.names = FALSE, sep = ',')
}

# Build a combined set with just the fold changes.
samples$rdsFiles = file.path(baseRepository,
                            'fragpipeAnalysis/processedData',
                            paste(samples$name,'_',samples$replicate,'_processed.rds',sep=''))

allRdsFiles = tibble()
for (i in 1:nrow(samples)){
  inputTemp = readRDS(samples$rdsFiles[i]) %>%
    dplyr::select(dataset, accession, gene, log2HL)
  
  allRdsFiles = rbind(allRdsFiles, inputTemp)
}

allRdsFilesWide = allRdsFiles %>%
  tidyr::pivot_wider(id_cols = accession:gene, names_from = 'dataset', values_from = 'log2HL')

saveRDS(allRdsFilesWide,
        paste(baseRepository,
              '/fragpipeAnalysis/processedData/dataset_allSamplesLog2HL.rds', sep =''))
write.table(allRdsFilesWide,
        paste(baseRepository,
              '/fragpipeAnalysis/processedData/dataset_allSamplesLog2HL.csv', sep =''),
        quote = FALSE, col.names = TRUE, row.names = FALSE, sep = ',')

#############################################################
# Figure 1C, ChIP-MS Fragpipe ouput clean up
#############################################################
# setup
library(pacman)
p_load(tidyverse, janitor)
processed_dir <- #Path to processed fragpipe data

# Append peptides
df <- read_csv(paste0(processed_dir, "KMT2D res_positive in at least 2 KO.csv")) %>% # manually selected proteins that were positive in at least two KO samples
    select(accession, gene,
        LogFC_HL_kmt2d_R1KO1 = kmt2d_R1KO1,
        LogFC_HL_kmt2d_R2KO1 = kmt2d_R2KO1,
        LogFC_HL_kmt2d_R1KO3 = kmt2d_R1KO3,
        LogFC_HL_kmt2d_R2KO3 = kmt2d_R2KO3,
        num_reps_with_peptides = "# rep"
    ) %>%
    group_by(accession) %>%
    mutate(
        num_reps_enriched = sum(c(
            (LogFC_HL_kmt2d_R1KO1 > 0),
            (LogFC_HL_kmt2d_R2KO1 > 0),
            (LogFC_HL_kmt2d_R1KO3 > 0),
            (LogFC_HL_kmt2d_R2KO3 > 0)
        ), na.rm = TRUE)
    ) %>%
    ungroup() %>%
    mutate(
        mean_reps_enriched =
            pmap_dbl(
                select(., LogFC_HL_kmt2d_R1KO1, LogFC_HL_kmt2d_R2KO1, LogFC_HL_kmt2d_R1KO3, LogFC_HL_kmt2d_R2KO3),
                function(...) {
                    pos <- c(...)[c(...) > 0]
                    mean(c(pos), na.rm = TRUE)
                }
            )
    ) %>%
    left_join(., read_tsv(paste0(processed_dir, "kmt2d_1/protein.tsv")) %>%
        clean_names() %>%
        select(
            accession = protein_id,
            total_peptides_kmt2d_1 = total_peptides,
            unique_peptides_kmt2d_1 = unique_peptides,
            razor_peptides_kmt2d_1 = razor_peptides
        )) %>%
    left_join(., read_tsv(paste0(processed_dir, "kmt2d_2/protein.tsv")) %>%
        clean_names() %>%
        select(
            accession = protein_id,
            total_peptides_kmt2d_2 = total_peptides,
            unique_peptides_kmt2d_2 = unique_peptides,
            razor_peptides_kmt2d_2 = razor_peptides
        )) %>%
    left_join(., read_tsv(paste0(processed_dir, "kmt2d_3/protein.tsv")) %>%
        clean_names() %>%
        select(
            accession = protein_id,
            total_peptides_kmt2d_3 = total_peptides,
            unique_peptides_kmt2d_3 = unique_peptides,
            razor_peptides_kmt2d_3 = razor_peptides
        )) %>%
    left_join(., read_tsv(paste0(processed_dir, "kmt2d_4/protein.tsv")) %>%
        clean_names() %>%
        select(
            accession = protein_id,
            total_peptides_kmt2d_4 = total_peptides,
            unique_peptides_kmt2d_4 = unique_peptides,
            razor_peptides_kmt2d_4 = razor_peptides
        )) %>%
    write_csv(paste0(processed_dir, "/KMT2D_res_positive_in_at_least_2_KO_cleaned.csv"))

# Following this, manually removed proteins that had zero peptides in all samples