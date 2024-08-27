library(UpSetR)
library(reshape2)
library(dplyr)
library(ggplot2)
library(stringr)
library(RColorBrewer)
library(vegan)

figures_folder <- 'Figures/paper_figures/'

iHMP_upset_plot <- function() {
  taxa_table <- read.csv('iHMP/analysis/data/metaphlan4_taxonomic_profiles.tsv', sep = '\t', skip = 1)
  
  colnames(taxa_table) <- gsub('_taxonomic$', '', colnames(taxa_table))
  taxa_table <- taxa_table[taxa_table$clade_name == 'UNCLASSIFIED' | 
                             (grepl('\\|s__', taxa_table$clade_name) &
                                !grepl('\\|t__', taxa_table$clade_name)),]
  rownames(taxa_table) <- taxa_table$clade_name
  taxa_table$clade_name <- NULL
  
  # Filtered to be common
  keep_taxa <- gsub(".*s__", "", rownames(taxa_table)[rowMeans(taxa_table > 0.1) > 0.01])
  
  all_results <- list.files('iHMP/analysis/results/', full.names = T)
  
  growing_df <- data.frame()
  for (result in all_results[grepl('\\/ibd_', all_results)]) {
    fit_out_joint <- read.csv(result, sep = '\t')
    
    fit_out_joint <- fit_out_joint[!is.na(fit_out_joint$feature),]
    fit_out_joint$feature <- gsub('.*s__', '', fit_out_joint$feature)
    
    if (grepl('Maaslin2', result)) {
      fit_out_joint <- fit_out_joint[fit_out_joint$metadata %in% c('diagnosis', 'dysbiosis_state'),]
      fit_out_joint <- fit_out_joint[!is.na(fit_out_joint$feature),]
      fit_out_joint$metadata_value <- paste0(fit_out_joint$metadata, '_', fit_out_joint$value)
      fit_out_joint$tool <- 'Maaslin2'
    }
    
    if (grepl('Maaslin3', result)) {
      fit_out_joint <- fit_out_joint[fit_out_joint$metadata %in% c('diagnosis', 'dysbiosis_state'),]
      fit_out_joint <- fit_out_joint[!is.na(fit_out_joint$feature),]
      fit_out_joint <- fit_out_joint[is.na(fit_out_joint$error),]
      fit_out_joint$metadata_value <- paste0(fit_out_joint$metadata, '_', fit_out_joint$value)
      fit_out_joint$pval <- fit_out_joint$pval_joint
      fit_out_joint$qval <- fit_out_joint$qval_joint
      if (grepl('Maaslin3CompAdjust', result)) {
        fit_out_joint$tool <- 'MaAsLin 3 Median\nAdjusted'
      } else {
        fit_out_joint$tool <- 'MaAsLin 3 Relative'
      }
    }
    
    if (grepl('ANCOMBC', result)) {
      fit_out_joint <- fit_out_joint[grepl('diagnosis|dysbiosis_state', fit_out_joint$metadata),]
      fit_out_joint <- fit_out_joint[!is.na(fit_out_joint$feature),]
      fit_out_joint <- fit_out_joint[is.na(fit_out_joint$error),]
      fit_out_joint <- fit_out_joint[!is.infinite(fit_out_joint$effect_size),]
      fit_out_joint$metadata_value <- gsub('dysbiosis_state|diagnosis', '', fit_out_joint$metadata)
      fit_out_joint$metadata_value <- ifelse(fit_out_joint$metadata_value %in% c('UC', 'CD'), 
                                             paste0('diagnosis_', fit_out_joint$metadata_value),
                                             paste0('dysbiosis_state_', fit_out_joint$metadata_value))
      fit_out_joint$tool <- 'ANCOMBC'
      fit_out_joint$coef <- fit_out_joint$effect_size
      associations$effect_size <- associations$effect_size / log(2) # Inflate because transformation is originally base e
    }
    
    if (grepl('ALDEx2', result)) {
      fit_out_joint <- fit_out_joint[grepl('diagnosis|dysbiosis_state', fit_out_joint$metadata),]
      fit_out_joint <- fit_out_joint[!is.na(fit_out_joint$feature),]
      fit_out_joint$metadata_value <- gsub('dysbiosis_state|diagnosis', '', fit_out_joint$metadata)
      fit_out_joint$metadata_value <- ifelse(fit_out_joint$metadata_value %in% c('UC', 'CD'), 
                                             paste0('diagnosis_', fit_out_joint$metadata_value),
                                             paste0('dysbiosis_state_', fit_out_joint$metadata_value))
      fit_out_joint$tool <- 'ALDEx2'
      fit_out_joint$coef <- fit_out_joint$effect_size
    }
    
    fit_out_joint <- fit_out_joint[,c('feature', 'metadata_value', 'coef', 'pval', 'qval', 'tool')]
    
    growing_df <- rbind(growing_df, fit_out_joint)
  }
  
  growing_df <- growing_df[growing_df$qval < 0.1 & growing_df$feature %in% keep_taxa,]
  growing_df$association <- paste0(growing_df$feature, '-', growing_df$metadata_value)
  growing_df <- growing_df[abs(growing_df$coef) > 1,]
  growing_df <- growing_df[growing_df$tool != 'MaAsLin 3 Relative',]
  
  # MaAsLin 3 overlap
  print(mean(unique(growing_df$association[growing_df$tool == 'MaAsLin 3 Median\nAdjusted']) %in% 
                 growing_df$association[growing_df$tool == 'Maaslin2']))
  print(mean(unique(growing_df$association[growing_df$tool == 'MaAsLin 3 Median\nAdjusted']) %in% 
                 growing_df$association[growing_df$tool == 'ALDEx2']))
  print(mean(unique(growing_df$association[growing_df$tool == 'MaAsLin 3 Median\nAdjusted']) %in% 
                 growing_df$association[growing_df$tool == 'ANCOMBC']))
  
  # ALDEx2 only
  print(mean(!unique(growing_df$association[growing_df$tool == 'ALDEx2']) %in% 
                 growing_df$association[growing_df$tool != 'ALDEx2']))

  listInput <- list(`MaAsLin 2` = unique(growing_df$association[growing_df$tool == 'Maaslin2']), 
       `MaAsLin 3` = unique(growing_df$association[growing_df$tool == 'MaAsLin 3 Median\nAdjusted']), 
       `ANCOM-BC2` = unique(growing_df$association[growing_df$tool == 'ANCOMBC']),
       ALDEx2 = unique(growing_df$association[growing_df$tool == 'ALDEx2']))
  
  plot_out <- upset(fromList(listInput), order.by = "freq", text.scale = 2.5)
  
  png(file=paste0(figures_folder, 'HMP2_upset.png'), width = 14, height = 6, res = 1000, units = 'in')
  plot_out
  dev.off()
}
iHMP_upset_plot()
  
recreate_summary_plot <- function() {
  # Clean up default HMP2 plot
  # Rename results file with clean titles
  all_results <- read.csv('iHMP/analysis/fit_out_Maaslin3CompAdjust/all_results.tsv', sep='\t')
  all_results <- all_results %>%
    mutate(metadata = case_when(metadata == 'consent_age' ~ 'Age',
                                metadata == 'Antibiotics' ~ 'Abx',
                                metadata == 'diagnosis' ~ 'Diagnosis',
                                metadata == 'dysbiosis_state' ~ 'Dysbiosis',
                                metadata == 'reads_filtered' ~ 'Read depth'),
           value = case_when(value == 'dysbiosis_CD' ~ 'CD',
                             value == 'dysbiosis_UC' ~ 'UC',
                             value == 'Yes' ~ 'Used', # Antibiotics
                             value == 'consent_age' ~ 'Age',
                             value == 'reads_filtered' ~ 'Read depth',
                             TRUE ~ value),
           feature = gsub('.*\\.', '', feature) %>%
             gsub(pattern = 's__', replacement = '') %>% 
             gsub(pattern = '_', replacement = ' ') %>% 
             gsub(pattern = 'sp ', replacement = 'sp. '))
  
  # Write results
  write.table(all_results, 'iHMP/analysis/fit_out_Maaslin3CompAdjust/all_results.tsv', sep='\t', row.names = F)
  
  # Need to create param_list
  taxa_table <- read.csv('iHMP/analysis/data/metaphlan4_taxonomic_profiles.tsv', skip = 1, check.names = F, sep = '\t')
  
  # Reorganize taxa table
  colnames(taxa_table) <- gsub('_taxonomic$', '', colnames(taxa_table))
  taxa_table <- taxa_table[taxa_table$clade_name == 'UNCLASSIFIED' | 
                             (grepl('\\|s__', taxa_table$clade_name) &
                                !grepl('\\|t__', taxa_table$clade_name)),]
  rownames(taxa_table) <- taxa_table$clade_name
  taxa_table$clade_name <- NULL
  
  prepare_metadata <- function(dataset_type) {
    metadata <- read.csv('iHMP/analysis/data/hmp2_metadata_2018-08-20.csv', check.names = F)
    
    # Reorganize metadata table
    for (participant_id in unique(metadata$`Participant ID`)) {
      bmis <- metadata$BMI[metadata$`Participant ID` == participant_id]
      bmis <- bmis[!is.na(bmis)]
      metadata$BMI[metadata$`Participant ID` == participant_id] <- ifelse(length(bmis) > 0, mean(bmis), NA)
      
      smoke_status <- metadata$`smoking status`[metadata$`Participant ID` == participant_id]
      smoke_status <- smoke_status[!is.na(smoke_status)]
      metadata$`smoking status`[metadata$`Participant ID` == participant_id] <- ifelse(length(smoke_status) > 0, smoke_status[1], NA)
    }
    
    if (dataset_type == 'taxa') {
      metadata <- metadata[metadata$data_type == 'metagenomics',]
    } else {
      metadata <- metadata[metadata$data_type == 'metabolomics',]
    }
    rownames(metadata) <- metadata$`External ID`
    metadata <- metadata[,colSums(metadata == '', na.rm = T) != nrow(metadata)]
    keep_cols <- c('External ID', 'Participant ID', 'week_num', 'site_name', 'Age at diagnosis',
                   'Education Level', 'Occupation', 'consent_age', 'diagnosis',
                   colnames(metadata)[c(52:83, 85:111)], 'race', 'sex', 'BMI', 'reads_filtered')
    metadata <- metadata[,keep_cols]
    metadata <- metadata[,colSums(!is.na(metadata)) != 0]
    return(metadata)
  }
  metadata <- prepare_metadata('taxa')
  
  # Calculate dysbiosis score
  veg_dist_out <- vegdist(t(as.matrix(taxa_table)), method="bray")
  veg_dist_out <- as.matrix(veg_dist_out)
  dysbiosis_scores <- vector(length = nrow(veg_dist_out))
  for (i in seq_along(rownames(veg_dist_out))) {
    sample_name <- rownames(veg_dist_out)[i]
    healthy_subset <- colnames(taxa_table)[colnames(taxa_table) %in% 
                                             rownames(metadata[metadata$week_num > 20 & 
                                                                 metadata$diagnosis == 'nonIBD' &
                                                                 metadata$`Participant ID` != metadata[sample_name,]$`Participant ID`,])]
    
    dysbiosis_scores[i] <- median(veg_dist_out[sample_name, healthy_subset])
  }
  names(dysbiosis_scores) <- rownames(veg_dist_out)
  
  # Add the dysbiosis state to the metadata
  dysbiosis_df <- data.frame(sample=names(dysbiosis_scores), dysbiosis_score=dysbiosis_scores)
  metadata <- right_join(dysbiosis_df, metadata, by=c('sample'='External ID'))
  metadata$dysbiosis_state <- metadata$dysbiosis_score > quantile(metadata$dysbiosis_score[metadata$diagnosis == 'nonIBD'], 0.9, na.rm=T)
  metadata$dysbiosis_state <- ifelse(metadata$dysbiosis_state & metadata$diagnosis != 'nonIBD', paste0('dysbiosis_', metadata$diagnosis), 'none')
  metadata$dysbiosis_state <- factor(metadata$dysbiosis_state, levels = c('none', 'dysbiosis_UC', 'dysbiosis_CD'))
  
  dysbiosis_df <- metadata[, c("sample", "dysbiosis_state")]
  dysbiosis_df <- dysbiosis_df[order(dysbiosis_df$sample),]
  dysbiosis_df <- dysbiosis_df[!duplicated(dysbiosis_df$sample),]
  
  metadata <- prepare_metadata('taxa')
  metadata <- right_join(dysbiosis_df, metadata, by=c('sample'='External ID'))
  
  metadata$participant_id <- metadata$`Participant ID`
  metadata$diagnosis <- factor(metadata$diagnosis, levels = c('nonIBD', 'UC', 'CD'))
  
  # Make sure only paired MGX/MBX are used since we need the dysbiosis score
  metadata <- metadata[!is.na(metadata$dysbiosis_state),]
  rownames(metadata) <- metadata$sample
  
  taxa_table <- read.csv('iHMP/analysis/data/metaphlan4_taxonomic_profiles.tsv', skip = 1, check.names = F, sep = '\t')
  
  # Reorganize taxa table
  colnames(taxa_table) <- gsub('_taxonomic$', '', colnames(taxa_table))
  taxa_table <- taxa_table[taxa_table$clade_name == 'UNCLASSIFIED' | 
                             (grepl('\\|s__', taxa_table$clade_name) &
                                !grepl('\\|t__', taxa_table$clade_name)),]
  rownames(taxa_table) <- taxa_table$clade_name
  taxa_table$clade_name <- NULL
  
  metadata <- metadata[metadata$sample %in% colnames(taxa_table),]
  metadata <- metadata[as.numeric(mapvalues(metadata$sample, colnames(taxa_table), 1:ncol(taxa_table))),]
  
  tmp_fit_out <- paste0("iHMP/analysis/fit_out_Maaslin3CompAdjust")
  
  param_list <- list()
  
  keep_taxa <- c("Faecalibacterium prausnitzii", 
                 "Bacteroides uniformis",
                 "Eubacterium rectale",
                 "Phocaeicola vulgatus",
                 "Bacteroides ovatus",
                 "Blautia obeum",
                 "Roseburia inulinivorans",
                 "Roseburia hominis",
                 "Clostridium sp. AT4",
                 "Enterocloster bolteae", 
                 "Ruminococcus gnavus", 
                 "Veillonella parvula",
                 "Clostridium neonatale",
                 "Klebsiella pneumoniae",
                 "Ruminococcus torques",
                 "Escherichia coli",
                 "Akkermansia muciniphila",
                 "Bacteroides fragilis",
                 "Bacteroides thetaiotaomicron",
                 "Dorea formicigenerans")
  
  results_in <- read.csv(file.path(tmp_fit_out, "all_results.tsv"), sep='\t')
  results_in <- results_in[results_in$feature %in% keep_taxa,]
  write.table(results_in, file.path(tmp_fit_out, "all_results.tsv"), sep = '\t', row.names = F)
  
  results_in <- read.csv(file.path(tmp_fit_out, "significant_results.tsv"), sep='\t')
  results_in <- results_in[results_in$feature %in% keep_taxa,]
  write.table(results_in, file.path(tmp_fit_out, "significant_results.tsv"), sep = '\t', row.names = F)
    
  # Set the new heatmap and coefficient plot variables and order them
  maaslin_plot_results_from_output(metadata = metadata, 
                                   output = tmp_fit_out, 
                                   normalization = 'TSS', 
                                   transform = 'LOG', 
                                   plot_summary_plot = T, 
                                   plot_associations = F, 
                                   max_significance = 0.1, 
                                   median_comparison_abundance = T, 
                                   median_comparison_prevalence = F, 
                                   heatmap_vars = c('Diagnosis CD', 'Diagnosis UC', 'Abx Used', 'Age', 'Read depth'),
                                   coef_plot_vars = c('Dysbiosis CD', 'Dysbiosis UC'))
}
recreate_summary_plot()
  
# Diet associations
diet_associations_plot <- function() {
  hmp2_files <- list.files('iHMP/analysis/results/', full.names = T)
  hmp2_files <- hmp2_files[grepl('food_associations', hmp2_files)]
  
  growing_df <- data.frame()
  for (hmp2_file in hmp2_files) {
    in_data <- read.csv(hmp2_file, sep='\t')
    in_data <- in_data[in_data$metadata == 'food_group',]
    in_data$food_group <- gsub('.*food_associations_|\\.tsv', '', hmp2_file)
    in_data$variable_type <- gsub('.*\\/|_food_associations_.*', '', hmp2_file)
    growing_df <- rbind(growing_df, in_data)
  }
  
  signif_level <- 0.1
  
  signif_subset <- growing_df[!is.na(growing_df$qval_individual) & 
               growing_df$qval_individual < signif_level & 
               is.na(growing_df$error),]
  
  signif_subset <- signif_subset[order(signif_subset$qval_individual),]
  signif_subset <- signif_subset[signif_subset$N.not.zero > 0.1 * signif_subset$N,]
  
  all_combinations <- expand.grid(food_group = unique(signif_subset$food_group), 
                                  variable_type = unique(signif_subset$variable_type))
  
  ordered_assoc <- paste0(signif_subset$feature, "_", signif_subset$food_group)[
      signif_subset$variable_type == 'ordered']
  group_assoc <- paste0(signif_subset$feature, "_", signif_subset$food_group)[
      signif_subset$variable_type == 'group']
  ordered_assoc <- unique(ordered_assoc)
  
  # Overlap in group and ordered predictors
  print(mean(ordered_assoc %in% group_assoc))
  print(mean(group_assoc %in% ordered_assoc))
  
  plot_df <- signif_subset %>%
    group_by(food_group, variable_type) %>%
    dplyr::summarize(count = n(), .groups = 'drop') %>%
    right_join(all_combinations, by = c("food_group", "variable_type")) %>%
    replace(is.na(.), 0)
  
  plot_df$food_group <- gsub("_", " ", str_to_title(plot_df$food_group))
  plot_df$food_group <- factor(plot_df$food_group, 
                               levels = unique(plot_df$food_group[order(plot_df$count, decreasing = T)]))
  plot_df <- plot_df %>% 
    mutate(variable_type = case_when(variable_type == 'group' ~ 'Group predictor',
                                     variable_type == 'ordered' ~ 'Ordered predictor'))
  
  plot_out <- ggplot(plot_df, aes(x = food_group, y = count, fill = variable_type)) + 
    geom_bar(position="dodge", stat="identity") + 
    theme_bw() + 
    ylab("Significant associations") + 
    xlab(NULL) + 
    scale_fill_manual(values = brewer.pal(12, "Paired")[9:10]) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          text = element_text(size = 16)) + 
    labs(fill = '')
  ggsave(plot = plot_out, filename = paste0(figures_folder, 'HMP2_diet.png'), width = 12, height = 4)
}
diet_associations_plot()

in_text_vals <- function() {
    associations_in <- read.csv("iHMP/analysis/results//ibd_associations_Maaslin3CompAdjust.tsv", sep = '\t')
    associations_in <- associations_in %>%
        filter(metadata %in% c("diagnosis", "dysbiosis_state"),
               qval_individual < 0.1,
               is.na(error),
               abs(coef) > 1)
    
    print(table(associations_in$association, ifelse(associations_in$coef < 0, "negative", "positive")))
    
    keep_taxa <- c("Faecalibacterium prausnitzii", 
                   "Bacteroides uniformis",
                   "Eubacterium rectale",
                   "Phocaeicola vulgatus",
                   "Bacteroides ovatus",
                   "Blautia obeum",
                   "Roseburia inulinivorans",
                   "Roseburia hominis",
                   "Clostridium sp. AT4",
                   "Enterocloster bolteae", 
                   "Ruminococcus gnavus", 
                   "Veillonella parvula",
                   "Clostridium neonatale",
                   "Klebsiella pneumoniae",
                   "Ruminococcus torques",
                   "Escherichia coli",
                   "Akkermansia muciniphila",
                   "Bacteroides fragilis",
                   "Bacteroides thetaiotaomicron",
                   "Dorea formicigenerans")
    
    associations_in <- read.csv("iHMP/analysis/results//ibd_associations_Maaslin3CompAdjust.tsv", sep = '\t')
    associations_in <- associations_in %>%
        filter(metadata %in% c("diagnosis", "dysbiosis_state"),
               is.na(error))
    
    associations_in$feature <- gsub('.*s__', '', associations_in$feature) %>%
        gsub(pattern = "_", replacement = " ")
    
    associations_in[associations_in$feature %in% keep_taxa & 
                        grepl("Roseburia", associations_in$feature),]
    
}
in_text_vals()


