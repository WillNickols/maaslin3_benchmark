library(UpSetR)
library(reshape2)
library(dplyr)
library(ggplot2)
library(stringr)
library(RColorBrewer)
library(vegan)
library(plyr)
library(maaslin3)

figures_folder <- 'Figures/paper_figures/'

HMP2_upset_plot_v4 <- function() {
  taxa_table <- read.csv('HMP2/analysis/data/metaphlan4_taxonomic_profiles.tsv', sep = '\t', skip = 1)
  
  colnames(taxa_table) <- gsub('_taxonomic$', '', colnames(taxa_table))
  taxa_table <- taxa_table[taxa_table$clade_name == 'UNCLASSIFIED' | 
                             grepl('\\|t__', taxa_table$clade_name),]
  rownames(taxa_table) <- taxa_table$clade_name
  taxa_table$clade_name <- NULL
  
  # Filtered to be common
  keep_taxa <- gsub(".*t__", "", rownames(taxa_table)[rowMeans(taxa_table > 0.1) > 0.01])
  
  all_results <- list.files('HMP2/analysis/results/', full.names = T)
  
  growing_df <- data.frame()
  for (result in all_results[grepl('\\/v4_ibd_', all_results)]) {
    fit_out_joint <- read.csv(result, sep = '\t')
    
    fit_out_joint <- fit_out_joint[!is.na(fit_out_joint$feature),]
    fit_out_joint$feature <- gsub('.*t__', '', fit_out_joint$feature)
    
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
      fit_out_joint$effect_size <- fit_out_joint$effect_size / log(2) # Inflate because transformation is originally base e
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
  
  m3_not_m2 <- unique(growing_df$association[growing_df$tool == 'MaAsLin 3 Median\nAdjusted'])[
      !unique(growing_df$association[growing_df$tool == 'MaAsLin 3 Median\nAdjusted']) %in% 
                 unique(growing_df$association[growing_df$tool == 'Maaslin2'])]
  
  m3_not_m2_df <- growing_df[growing_df$feature %in% gsub('-.*', '', m3_not_m2) & 
                 growing_df$metadata_value %in% gsub('.*-', '', m3_not_m2) &
                 growing_df$tool == 'MaAsLin 3 Median\nAdjusted',]
  
  sum((m3_not_m2_df %>% 
      dplyr::group_by(association) %>%
      dplyr::mutate(agreement = ifelse(n()==2, sign(coef[1]) == sign(coef[2]), NA)) %>%
      dplyr::summarise(mean(agreement)))[,2] == 0, na.rm=T)
  
  # ALDEx2 only
  print(mean(!unique(growing_df$association[growing_df$tool == 'ALDEx2']) %in% 
                 growing_df$association[growing_df$tool != 'ALDEx2']))
  
  table(growing_df$metadata_value[growing_df$tool == 'ALDEx2'])
  tmp_df <- growing_df[growing_df$tool == 'MaAsLin 3 Median\nAdjusted',]
  table(tmp_df$metadata_value[!duplicated(tmp_df$association)])
  
  print(sum(!unique(growing_df$association[growing_df$tool == 'ALDEx2' & 
                                               grepl("dysbiosis", growing_df$metadata_value)]) %in% 
                 growing_df$association[growing_df$tool == 'MaAsLin 3 Median\nAdjusted']))
  print(sum(unique(growing_df$association[growing_df$tool == 'ALDEx2' & 
                                               grepl("dysbiosis", growing_df$metadata_value)]) %in% 
                growing_df$association[growing_df$tool == 'MaAsLin 3 Median\nAdjusted']))
  
  listInput <- list(`MaAsLin 2` = unique(growing_df$association[growing_df$tool == 'Maaslin2']), 
       `MaAsLin 3` = unique(growing_df$association[growing_df$tool == 'MaAsLin 3 Median\nAdjusted']), 
       `ANCOM-BC2` = unique(growing_df$association[growing_df$tool == 'ANCOMBC']),
       ALDEx2 = unique(growing_df$association[growing_df$tool == 'ALDEx2']))
  
  plot_out <- upset(fromList(listInput), order.by = "freq", text.scale = 2.5)
  
  png(file=paste0(figures_folder, 'HMP2_upset_v4.png'), width = 14, height = 6, res = 1000, units = 'in')
  plot_out
  dev.off()
}
HMP2_upset_plot_v4()

HMP2_upset_plot_v3 <- function() {
    taxa_table <- read.csv('HMP2/analysis/data/metaphlan3_taxonomic_profiles.tsv',check.names = F, sep = '\t')
    
    taxa_table <- taxa_table[taxa_table$`#SampleID` == 'UNCLASSIFIED' | 
                                 grepl('\\|s__', taxa_table$`#SampleID`) & 
                                 !grepl('\\|t__', taxa_table$`#SampleID`),]
    rownames(taxa_table) <- taxa_table$`#SampleID`
    taxa_table$`#SampleID` <- NULL
    taxa_table['UNCLASSIFIED',] <- pmax(1 - colSums(taxa_table), 0)
    taxa_table <- taxa_table * 100 # Convert to percents to be consistent with v4
    
    # Filtered to be common
    keep_taxa <- gsub(".*s__", "", rownames(taxa_table)[rowMeans(taxa_table > 0.1) > 0.01])
    
    all_results <- list.files('HMP2/analysis/results/', full.names = T)
    
    growing_df <- data.frame()
    for (result in all_results[grepl('\\/v3_ibd_', all_results)]) {
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
            fit_out_joint$effect_size <- fit_out_joint$effect_size / log(2) # Inflate because transformation is originally base e
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
    
    table(growing_df$metadata_value[growing_df$tool == 'ALDEx2'])
    tmp_df <- growing_df[growing_df$tool == 'MaAsLin 3 Median\nAdjusted',]
    table(tmp_df$metadata_value[!duplicated(tmp_df$association)])
    
    print(sum(!unique(growing_df$association[growing_df$tool == 'ALDEx2' & 
                                                 grepl("dysbiosis", growing_df$metadata_value)]) %in% 
                  growing_df$association[growing_df$tool == 'MaAsLin 3 Median\nAdjusted']))
    
    listInput <- list(`MaAsLin 2` = unique(growing_df$association[growing_df$tool == 'Maaslin2']), 
                      `MaAsLin 3` = unique(growing_df$association[growing_df$tool == 'MaAsLin 3 Median\nAdjusted']), 
                      `ANCOM-BC2` = unique(growing_df$association[growing_df$tool == 'ANCOMBC']),
                      ALDEx2 = unique(growing_df$association[growing_df$tool == 'ALDEx2']))
    
    plot_out <- upset(fromList(listInput), order.by = "freq", text.scale = 2.5)
    
    png(file=paste0(figures_folder, 'HMP2_upset_v3.png'), width = 14, height = 6, res = 1000, units = 'in')
    plot_out
    dev.off()
}
HMP2_upset_plot_v3()

prepare_metadata <- function(dataset_type) {
    metadata <- read.csv('HMP2/analysis/data/hmp2_metadata_2018-08-20.csv', check.names = F)
    
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

process_taxa_table_and_metadata <- function(taxa_table, metadata) {
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
    metadata$dysbiosis_state <- ifelse(metadata$dysbiosis_state, paste0('dysbiosis_', metadata$diagnosis), 'none')
    metadata$dysbiosis_state <- factor(metadata$dysbiosis_state, levels = c('none', 'dysbiosis_nonIBD', 'dysbiosis_UC', 'dysbiosis_CD'))
    
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
    
    return(list(taxa_table = taxa_table, metadata = metadata))
}

# Diet associations
diet_associations_plot_v4 <- function() {
  hmp2_files <- list.files('HMP2/analysis/results/', full.names = T)
  hmp2_files <- hmp2_files[grepl('food_associations', hmp2_files) & grepl('v4', hmp2_files)]
  
  growing_df <- data.frame()
  for (hmp2_file in hmp2_files) {
    in_data <- read.csv(hmp2_file, sep='\t')
    in_data <- in_data[in_data$metadata == 'food_group',]
    in_data$food_group <- gsub('.*food_associations_|\\.tsv', '', hmp2_file)
    in_data$variable_type <- gsub('.*\\/v[0-9]_|_food_associations_.*', '', hmp2_file)
    growing_df <- rbind(growing_df, in_data)
  }
  
  signif_level <- 0.1
  
  signif_subset <- growing_df[!is.na(growing_df$qval_individual) & 
               growing_df$qval_individual < signif_level & 
               is.na(growing_df$error),]
  
  signif_subset <- signif_subset[order(signif_subset$qval_individual),]
  signif_subset <- signif_subset[signif_subset$N_not_zero > 0.1 * signif_subset$N,]
  
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
  print(sum(ordered_assoc %in% group_assoc))
  print(length(ordered_assoc))
  print(sum(group_assoc %in% ordered_assoc))
  print(length(group_assoc))
  
  growing_df[grepl('Roseburia_faecis', growing_df$feature) & growing_df$food_group == 'alcohol',]
  
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
  
  growing_df[grepl("s__Roseburia_faecis", growing_df$feature) & growing_df$food_group == 'alcohol',]
  
  plot_out <- ggplot(plot_df, aes(x = food_group, y = count, fill = variable_type)) + 
    geom_bar(position="dodge", stat="identity") + 
    theme_bw() + 
    ylab("Significant associations") + 
    xlab(NULL) + 
    scale_fill_manual(values = brewer.pal(12, "Paired")[9:10]) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          text = element_text(size = 16)) + 
    labs(fill = '')
  ggsave(plot = plot_out, filename = paste0(figures_folder, 'HMP2_diet_v4.png'), width = 12, height = 4)
}
diet_associations_plot_v4()

diet_associations_plot_v3 <- function() {
    hmp2_files <- list.files('HMP2/analysis/results/', full.names = T)
    hmp2_files <- hmp2_files[grepl('food_associations', hmp2_files) & grepl('v3', hmp2_files)]
    
    growing_df <- data.frame()
    for (hmp2_file in hmp2_files) {
        in_data <- read.csv(hmp2_file, sep='\t')
        in_data <- in_data[in_data$metadata == 'food_group',]
        in_data$food_group <- gsub('.*food_associations_|\\.tsv', '', hmp2_file)
        in_data$variable_type <- gsub('.*\\/v[0-9]_|_food_associations_.*', '', hmp2_file)
        growing_df <- rbind(growing_df, in_data)
    }
    
    signif_level <- 0.1
    
    signif_subset <- growing_df[!is.na(growing_df$qval_individual) & 
                                    growing_df$qval_individual < signif_level & 
                                    is.na(growing_df$error),]
    
    signif_subset <- signif_subset[order(signif_subset$qval_individual),]
    signif_subset <- signif_subset[signif_subset$N_not_zero > 0.1 * signif_subset$N,]
    
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
    print(sum(ordered_assoc %in% group_assoc))
    print(length(ordered_assoc))
    print(sum(group_assoc %in% ordered_assoc))
    print(length(group_assoc))
    
    growing_df[grepl('Roseburia_faecis', growing_df$feature) & growing_df$food_group == 'alcohol',]
    
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
    
    growing_df[grepl("s__Roseburia_faecis", growing_df$feature) & growing_df$food_group == 'alcohol',]
    
    plot_out <- ggplot(plot_df, aes(x = food_group, y = count, fill = variable_type)) + 
        geom_bar(position="dodge", stat="identity") + 
        theme_bw() + 
        ylab("Significant associations") + 
        xlab(NULL) + 
        scale_fill_manual(values = brewer.pal(12, "Paired")[9:10]) + 
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
              text = element_text(size = 16)) + 
        labs(fill = '')
    ggsave(plot = plot_out, filename = paste0(figures_folder, 'HMP2_diet_v3.png'), width = 12, height = 4)
}
diet_associations_plot_v3()

diet_associations_joint <- function() {
    hmp2_files <- list.files('HMP2/analysis/results/', full.names = T)
    hmp2_files <- hmp2_files[grepl('food_associations', hmp2_files) & grepl('v4', hmp2_files)]
    
    growing_df <- data.frame()
    for (hmp2_file in hmp2_files) {
        in_data <- read.csv(hmp2_file, sep='\t')
        in_data <- in_data[in_data$metadata == 'food_group',]
        in_data$food_group <- gsub('.*food_associations_|\\.tsv', '', hmp2_file)
        in_data$variable_type <- gsub('.*\\/v[0-9]_|_food_associations_.*', '', hmp2_file)
        growing_df <- rbind(growing_df, in_data)
    }
    
    v4_results <- growing_df
    v4_results$feature <- gsub(".*\\.s__|\\.t__.*", "", v4_results$feature)
    
    hmp2_files <- list.files('HMP2/analysis/results/', full.names = T)
    hmp2_files <- hmp2_files[grepl('food_associations', hmp2_files) & grepl('v3', hmp2_files)]
    
    growing_df <- data.frame()
    for (hmp2_file in hmp2_files) {
        in_data <- read.csv(hmp2_file, sep='\t')
        in_data <- in_data[in_data$metadata == 'food_group',]
        in_data$food_group <- gsub('.*food_associations_|\\.tsv', '', hmp2_file)
        in_data$variable_type <- gsub('.*\\/v[0-9]_|_food_associations_.*', '', hmp2_file)
        growing_df <- rbind(growing_df, in_data)
    }
    
    growing_df$feature <- gsub(".*\\.s__|\\.t__.*", "", growing_df$feature)
    
    joined_results <- full_join(v4_results, growing_df, c("feature", "metadata", "value", "association", "food_group", "variable_type"))
    
    joined_results_sub <- joined_results[!is.na(joined_results$qval_individual.x) &
                                             !is.na(joined_results$qval_individual.y) & 
                                             joined_results$qval_individual.x < 0.1 &
                                             joined_results$qval_individual.y < 0.1,]
    
    joined_results_sub <- joined_results_sub %>%
        mutate(feature_food_group = interaction(feature, food_group))
    
    joined_results_ordered <- joined_results_sub %>% 
        filter(variable_type == "ordered")
    
    joined_results_group <- joined_results_sub %>% 
        filter(variable_type == "group") %>% 
        pull(feature_food_group)
    
    joined_results_ordered <- joined_results_ordered %>%
        mutate(shape = ifelse(feature_food_group %in% joined_results_group, "Omnibus\nassociation", "No omnibus\nassociation"))
    
    joined_results_ordered$label <- paste0(gsub("_", "\n", joined_results_ordered$feature), ",\n",
                                           gsub("_", " ", joined_results_ordered$food_group))
    
    joined_results_ordered$association <- ifelse(joined_results_ordered$association == 'abundance', 
                                                 "Abundance", "Prevalence")
    
    plot_out <- ggplot(joined_results_ordered, aes(x = coef.x, y = coef.y, shape = shape, color = association)) +
        geom_point(size = 3) +
        ggrepel::geom_text_repel(aes(label = label), seed = 7, max.overlaps = 100, point.padding = 0, force = 100, size = 3, lineheight=.7) +
        scale_shape_manual(values = c(5, 16)) +
        labs(x = "MetaPhlAn 4 coefficient", y = "MetaPhlAn 3 coefficient", color= "Association", shape = "Group testing") +
        scale_color_manual(values = c("Abundance" = "#8B008B", "Prevalence" = "#008B8B")) +
        theme_bw() + 
        geom_abline(slope = 1, intercept = 0, color = 'black', linetype = 'dashed') + 
        ggplot2::theme(
            axis.title = ggplot2::element_text(size = 10),
            axis.text.x = ggplot2::element_text(size = 8),
            axis.text.y = ggplot2::element_text(size = 8),
            legend.title = ggplot2::element_text(size = 10),
            legend.text = ggplot2::element_text(size = 8, 
                                                face = "plain"),
            legend.position = "right",
            legend.direction = 'vertical',
            legend.background = ggplot2::element_rect(
                fill = "transparent"),
            panel.spacing = ggplot2::unit(0, "lines"),
            panel.grid.minor = ggplot2::element_blank(),
            strip.text = ggplot2::element_text(size = 8),
            strip.background = ggplot2::element_rect(
                fill = "transparent")
        ) + 
        scale_x_continuous(breaks = seq(-6, 4, 2), limits = c(-6, 4)) + 
        scale_y_continuous(breaks = seq(-6, 4, 2), limits = c(-6, 4))
    
    ggsave(plot = plot_out, filename = paste0(figures_folder, 'v3_vs_v4_diet.png'), width = 5, height = 4)
    
    # In-text
    v4_results <- v4_results %>%
        filter(is.na(error) & !is.na(qval_individual) & qval_individual < 0.1)
    
    v3_results <- growing_df %>%
        filter(is.na(error) & !is.na(qval_individual) & qval_individual < 0.1)
    
    sum(unique(interaction(v4_results$feature, v4_results$food_group, v4_results$variable_type, v4_results$association)) %in%
             unique(interaction(v3_results$feature, v3_results$food_group, v3_results$variable_type, v3_results$association)))
    
    length(unique(interaction(v3_results$feature, v3_results$food_group, v3_results$variable_type, v3_results$association)))
    length(unique(interaction(v4_results$feature, v4_results$food_group, v4_results$variable_type, v4_results$association)))
}
diet_associations_joint()

dysosmobacter_welbionis_plot <- function() {
    taxa_table <- read.csv('HMP2/analysis/data/metaphlan4_taxonomic_profiles.tsv', skip = 1, check.names = F, sep = '\t')
    
    # Reorganize taxa table
    colnames(taxa_table) <- gsub('_taxonomic$', '', colnames(taxa_table))
    taxa_table <- taxa_table[taxa_table$clade_name == 'UNCLASSIFIED' | 
                                 grepl('\\|t__', taxa_table$clade_name),]
    rownames(taxa_table) <- taxa_table$clade_name
    taxa_table$clade_name <- NULL
    
    metadata <- prepare_metadata('taxa')
    
    taxa_table_and_metadata <- process_taxa_table_and_metadata(taxa_table, metadata)
    taxa_table <- taxa_table_and_metadata[['taxa_table']]
    metadata <- taxa_table_and_metadata[['metadata']]
    
    taxa_sub <- data.frame(sample = colnames(taxa_table),
                           abun = unlist(c(taxa_table['k__Bacteria|p__Firmicutes|c__Clostridia|o__Eubacteriales|f__Oscillospiraceae|g__Dysosmobacter|s__Dysosmobacter_welbionis|t__SGB15078',])) / 100)
    
    metadata_expanded <- left_join(metadata, taxa_sub, by = 'sample')
    metadata_expanded <- metadata_expanded %>%
        mutate(diagnosis_dysbiosis = interaction(diagnosis, dysbiosis_state)) %>%
        mutate(diagnosis_dysbiosis = case_when(diagnosis_dysbiosis == 'nonIBD.none' ~ 'non-IBD',
                                               diagnosis_dysbiosis == 'UC.none' ~ 'UC',
                                               diagnosis_dysbiosis == 'CD.none' ~ 'CD',
                                               diagnosis_dysbiosis == 'nonIBD.dysbiosis_nonIBD' ~ 'non-IBD (Dysbiosis)',
                                               diagnosis_dysbiosis == 'UC.dysbiosis_UC' ~ 'UC (Dysbiosis)',
                                               diagnosis_dysbiosis == 'CD.dysbiosis_CD' ~ 'CD (Dysbiosis)',)) %>%
        mutate(diagnosis_dysbiosis = factor(diagnosis_dysbiosis, levels = 
                                                rev(c('non-IBD', 'UC', 'CD', 'non-IBD (Dysbiosis)', 'UC (Dysbiosis)', 'CD (Dysbiosis)'))))

    associations_in <- read.csv("HMP2/analysis/results/v4_ibd_associations_Maaslin3CompAdjust.tsv", sep = '\t')
    associations_in <- associations_in %>%
        filter(feature == 'k__Bacteria.p__Firmicutes.c__Clostridia.o__Eubacteriales.f__Oscillospiraceae.g__Dysosmobacter.s__Dysosmobacter_welbionis.t__SGB15078',
               metadata %in% c("dysbiosis_state", "diagnosis"))
    signif_df <- data.frame(group1 = case_when(associations_in$value == 'CD' ~ 'non-IBD',
                                               associations_in$value == 'UC' ~ 'non-IBD',
                                               associations_in$value == 'dysbiosis_nonIBD' ~ 'non-IBD',
                                               associations_in$value == 'dysbiosis_UC' ~ 'UC',
                                               associations_in$value == 'dysbiosis_CD' ~ 'CD'),
                            group2 = case_when(associations_in$value == 'CD' ~ 'CD',
                                               associations_in$value == 'UC' ~ 'UC',
                                               associations_in$value == 'dysbiosis_nonIBD' ~ 'non-IBD (Dysbiosis)',
                                               associations_in$value == 'dysbiosis_UC' ~ 'UC (Dysbiosis)',
                                               associations_in$value == 'dysbiosis_CD' ~ 'CD (Dysbiosis)'),
                            qval_individual = ifelse(associations_in$qval_individual < 0.1,
                                                     paste0('q = ', formatC(associations_in$qval_individual, format = "e", digits = 1)),
                                                     'n.s.'),
                            association = associations_in$association)
    
    plot1 <- ggplot(metadata_expanded[metadata_expanded$abun > 0,], 
                    aes(x = diagnosis_dysbiosis, y = abun)) + 
        geom_boxplot(fill = '#8B008B', outlier.shape = NA) + 
        geom_jitter(width = 0.2, alpha = 1, color = "black", size = 0.5) + 
        theme_bw() +
        ggpubr::geom_bracket(
            data = signif_df[signif_df$association == 'abundance',], 
            aes(xmin = group1, xmax = group2, label = qval_individual),
            inherit.aes = FALSE,
            tip.length = 0.01,
            y.position = c(3, 1, 0.001, 2, 4),
            label.size = 3, coord.flip = T
        ) + 
        theme(axis.title.y = element_blank(), 
              axis.text.y = element_blank()) + 
        labs(y = 'Relative abundance') + 
        scale_y_continuous(trans = 'log', breaks = 10^c(-6, -4, -2, 0), limits = c(10^-6, 10^2)) + 
        coord_flip()
    
    metadata_expanded <- metadata_expanded %>%
        mutate(abun_status = ifelse(abun == 0, "Zero", "Non-zero"))
    
    proportion_data <- metadata_expanded %>%
        dplyr::group_by(diagnosis_dysbiosis) %>%
        dplyr::summarise(proportion_nonzero = mean(abun != 0), n = n(), .groups = "drop") %>%
        dplyr::mutate(diagnosis_dysbiosis_name = factor(
            paste0(diagnosis_dysbiosis, '\n(n = ', n, ')'),
            levels = paste0(diagnosis_dysbiosis, '\n(n = ', n, ')')))
    
    signif_df <- signif_df %>%
        mutate(group1 = plyr::mapvalues(group1, 
                                         as.character(proportion_data$diagnosis_dysbiosis), 
                                         as.character(proportion_data$diagnosis_dysbiosis_name)),
               group2 = plyr::mapvalues(group2, 
                                         as.character(proportion_data$diagnosis_dysbiosis), 
                                         as.character(proportion_data$diagnosis_dysbiosis_name)))

    plot2 <- ggplot(proportion_data, aes(x = diagnosis_dysbiosis_name, y = proportion_nonzero)) +
        geom_bar(stat = "identity", fill = "#008B8B") +
        labs(y = "Proportion Present") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        ggpubr::geom_bracket(
            data = signif_df[signif_df$association == 'prevalence',], 
            aes(xmin = group1, xmax = group2, label = qval_individual),
            inherit.aes = FALSE,
            tip.length = 0.01,
            y.position = c(1.4, 1.29, 1.16, 0.9, 1.03), label.size = 3, coord.flip = T
        ) + 
        theme(axis.title.y = element_blank(), panel.grid.minor = element_blank()) + 
        scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1.5)) + 
        coord_flip()
    
    plot_out <- gridExtra::grid.arrange(plot2, plot1, ncol = 2, widths = c(1.5, 1))
    ggsave(plot = plot_out, filename = paste0(figures_folder, 'Dysosmobacter_welbionis.png'), width = 6, height = 3)
}
dysosmobacter_welbionis_plot()

in_text_vals <- function() {
    taxa_table <- read.csv('HMP2/analysis/data/metaphlan4_taxonomic_profiles.tsv', skip = 1, check.names = F, sep = '\t')
    
    # Reorganize taxa table
    colnames(taxa_table) <- gsub('_taxonomic$', '', colnames(taxa_table))
    taxa_table <- taxa_table[taxa_table$clade_name == 'UNCLASSIFIED' | 
                                 grepl('\\|t__', taxa_table$clade_name),]
    rownames(taxa_table) <- taxa_table$clade_name
    taxa_table$clade_name <- NULL
    
    metadata <- prepare_metadata('taxa')
    
    taxa_table_and_metadata <- process_taxa_table_and_metadata(taxa_table, metadata)
    taxa_table <- taxa_table_and_metadata[['taxa_table']]
    metadata <- taxa_table_and_metadata[['metadata']]
    
    print(table(metadata$diagnosis))
    print(metadata %>% 
        dplyr::group_by(diagnosis) %>%
        dplyr::summarize(length(unique(participant_id))))
    
    print(table(metadata$diagnosis, metadata$consent_age <16))
    print(metadata %>% 
              dplyr::group_by(diagnosis, consent_age <16) %>%
              dplyr::summarize(length(unique(participant_id))))
    
    all_results_under16 <- read.csv('HMP2/analysis/results/v4_under16_Maaslin3.tsv', sep='\t')
    all_results_atleast16 <- read.csv('HMP2/analysis/results/v4_atleast16_Maaslin3.tsv', sep='\t')
    
    all_results_under16 <- full_join(all_results_under16[all_results_under16$association == 'abundance',], 
              all_results_under16[all_results_under16$association == 'prevalence',], by = c("feature", "metadata", "value"))
    
    #1. If the association is significant in both abundance/prevalence and matches signs, we take that direction
    #2. If the association is significant in only one, we take that direction
    #3. If the association is significant in both and conflicting, take the smaller q-value
    
    all_results_under16 <- all_results_under16 %>%
        filter(!is.na(metadata) & metadata %in% c('dysbiosis_state', 'diagnosis')) %>%
        filter((qval_individual.x < 0.1 & abs(coef.x) > 1) | (qval_individual.y < 0.1 & abs(coef.y) > 1)) %>%
        mutate(feature = gsub('.*\\.s__', '', feature)) %>%
        filter(!is.na(feature)) %>%
        mutate(association_direction = case_when(qval_individual.x < 0.1 & (qval_individual.y > 0.1 | is.na(qval_individual.y)) & coef.x > 0 ~ "increased",
                                                 qval_individual.x < 0.1 & (qval_individual.y > 0.1 | is.na(qval_individual.y)) & coef.x < 0 ~ "reduced",
                                                 (qval_individual.x > 0.1 | is.na(qval_individual.x)) & qval_individual.y < 0.1 & coef.y > 0 ~ "increased",
                                                 (qval_individual.x > 0.1 | is.na(qval_individual.x)) & qval_individual.y < 0.1 & coef.y < 0 ~ "reduced",
                                                 qval_individual.x < qval_individual.y & coef.x < 0 ~ "reduced",
                                                 qval_individual.x < qval_individual.y & coef.x > 0 ~ "increased",
                                                 qval_individual.x > qval_individual.y & coef.y < 0 ~ "reduced",
                                                 qval_individual.x > qval_individual.y & coef.y > 0 ~ "increased")) %>%
        mutate(qval_individual = pmin(qval_individual.x, qval_individual.y, na.rm = T)) %>%
        select(-qval_individual.x, -qval_individual.y)
    
    all_results_atleast16 <- full_join(all_results_atleast16[all_results_atleast16$association == 'abundance',], 
                                       all_results_atleast16[all_results_atleast16$association == 'prevalence',], by = c("feature", "metadata", "value"))
    
    all_results_atleast16 <- all_results_atleast16 %>%
        filter(!is.na(metadata) & metadata %in% c('dysbiosis_state', 'diagnosis')) %>%
        filter((qval_individual.x < 0.1 & abs(coef.x) > 1) | (qval_individual.y < 0.1 & abs(coef.y) > 1)) %>%
        mutate(feature = gsub('.*\\.s__', '', feature)) %>%
        filter(!is.na(feature)) %>%
        mutate(association_direction = case_when(qval_individual.x < 0.1 & (qval_individual.y > 0.1 | is.na(qval_individual.y)) & coef.x > 0 ~ "increased",
                                                 qval_individual.x < 0.1 & (qval_individual.y > 0.1 | is.na(qval_individual.y)) & coef.x < 0 ~ "reduced",
                                                 (qval_individual.x > 0.1 | is.na(qval_individual.x)) & qval_individual.y < 0.1 & coef.y > 0 ~ "increased",
                                                 (qval_individual.x > 0.1 | is.na(qval_individual.x)) & qval_individual.y < 0.1 & coef.y < 0 ~ "reduced",
                                                 qval_individual.x < qval_individual.y & coef.x < 0 ~ "reduced",
                                                 qval_individual.x < qval_individual.y & coef.x > 0 ~ "increased",
                                                 qval_individual.x > qval_individual.y & coef.y < 0 ~ "reduced",
                                                 qval_individual.x > qval_individual.y & coef.y > 0 ~ "increased")) %>%
        mutate(qval_individual = pmin(qval_individual.x, qval_individual.y, na.rm = T)) %>%
        select(-qval_individual.x, -qval_individual.y)
    
    joined_results <- full_join(all_results_under16, all_results_atleast16, by = c("feature", "metadata", "value"))
    joined_results <- joined_results %>%
        filter(!is.na(metadata) & metadata %in% c('dysbiosis_state', 'diagnosis')) %>%
        filter((qval_individual.x < 0.1 | qval_individual.y)) %>%
        mutate(feature = gsub('.*\\.s__', '', feature)) %>%
        filter(!is.na(feature))
    
    table(!is.na(joined_results$qval_individual.x) & joined_results$qval_individual.x < 0.1 & joined_results$association_direction.x == 'increased', 
          !is.na(joined_results$qval_individual.y) & joined_results$qval_individual.y < 0.1 & joined_results$association_direction.y == 'increased')
    
    table(joined_results$value[!is.na(joined_results$qval_individual.y) & joined_results$qval_individual.y < 0.1 & joined_results$association_direction.y == 'increased'])
    
    table(!is.na(joined_results$qval_individual.x) & joined_results$qval_individual.x < 0.1 & joined_results$association_direction.x == 'reduced', 
          !is.na(joined_results$qval_individual.y) & joined_results$qval_individual.y < 0.1 & joined_results$association_direction.y == 'reduced')
    
    joined_results %>% filter(!is.na(joined_results$qval_individual.x) & joined_results$qval_individual.x < 0.1 & joined_results$association_direction.x == 'increased')
    joined_results %>% filter(!is.na(joined_results$qval_individual.y) & joined_results$qval_individual.y < 0.1 & joined_results$association_direction.y == 'increased')
    joined_results %>% filter(!is.na(joined_results$qval_individual.x) & joined_results$qval_individual.x < 0.1 & joined_results$association_direction.x == 'reduced', 
                              !is.na(joined_results$qval_individual.y) & joined_results$qval_individual.y < 0.1 & joined_results$association_direction.y == 'reduced') %>%
        filter(grepl("vulga", feature))
    
    table(joined_results$qval_individual.x < 0.1 & joined_results$coef.x < 0, joined_results$qval_individual.y < 0.1 & joined_results$coef.y < 0)
    joined_results %>% filter(qval_individual.x < 0.1 & coef.x < 0 & qval_individual.y < 0.1 & coef.y < 0) %>%
        filter(grepl('Phocaeicola_vulgatus', feature))
    
    # MPA4 results
    all_results <- read.csv('HMP2/analysis/results/v4_ibd_associations_Maaslin3CompAdjust.tsv', sep='\t')
    all_results <- all_results %>%
        filter(!is.na(metadata) & metadata %in% c('dysbiosis_state', 'diagnosis')) %>%
        filter((qval_individual < 0.1 & abs(coef) > 1)) %>%
        mutate(feature = gsub('.*\\.s__', '', feature)) %>%
        filter(!is.na(feature))
    
    table(all_results$association, all_results$coef > 0)
    
    all_results <- all_results %>%
        filter(association == 'prevalence') %>%
        filter(value %in% c("dysbiosis_UC", "dysbiosis_CD")) %>%
        filter(coef < 0)
    table(all_results$value)
    
    sort(intersect(all_results$feature[all_results$value == 'dysbiosis_UC'],
              all_results$feature[all_results$value == 'dysbiosis_CD']))
    
    all_results[all_results$feature == 'Dysosmobacter_welbionis.t__SGB15078',]
    
    # MTX
    metadata <- prepare_metadata('taxa')
    
    taxa_table_and_metadata <- process_taxa_table_and_metadata(taxa_table, metadata)
    taxa_table <- taxa_table_and_metadata[['taxa_table']]
    metadata <- taxa_table_and_metadata[['metadata']]
    
    taxa_table <- read.csv('HMP2/analysis/data/metaphlan4_taxonomic_profiles.tsv', skip = 1, check.names = F, sep = '\t')
    
    # Reorganize taxa table
    colnames(taxa_table) <- gsub('_taxonomic$', '', colnames(taxa_table))
    taxa_table <- taxa_table[taxa_table$clade_name == 'UNCLASSIFIED' | 
                                 grepl('\\|t__', taxa_table$clade_name),]
    rownames(taxa_table) <- taxa_table$clade_name
    taxa_table$clade_name <- NULL
    
    metadata <- metadata[metadata$sample %in% colnames(taxa_table),]
    metadata <- metadata[as.numeric(mapvalues(metadata$sample, colnames(taxa_table), 1:ncol(taxa_table))),]
    
    mtx_table <- read.csv('HMP2/analysis/data/pathabundances_3_MTX.tsv', sep = '\t')
    rownames(mtx_table) <- mtx_table$Feature.Sample
    mtx_table$Feature.Sample <- NULL
    
    metadata <- metadata[metadata$sample %in% gsub("_pathabundance.*", "", colnames(mtx_table)),]
    
    table(metadata$diagnosis)
    print(metadata %>% 
              dplyr::group_by(diagnosis) %>%
              dplyr::summarize(length(unique(participant_id))))
    
    results_in <- read.csv("HMP2/analysis/results/mtx_covariate.tsv", sep='\t')
    results_in <- results_in %>%
        filter(is.na(error), qval_individual < 0.1, metadata %in% c("dysbiosis_state", "diagnosis"))
    
    table(results_in$value)
    table(results_in$association)
    
    mean(taxa_table == 0)
    mean(mtx_table == 0)
    
    # F praus oddness
    taxa_table <- read.csv('HMP2/analysis/data/metaphlan4_taxonomic_profiles.tsv', skip = 1, check.names = F, sep = '\t')
    
    # Reorganize taxa table
    colnames(taxa_table) <- gsub('_taxonomic$', '', colnames(taxa_table))
    taxa_table <- taxa_table[taxa_table$clade_name == 'UNCLASSIFIED' | 
                                 grepl('\\|t__', taxa_table$clade_name),]
    rownames(taxa_table) <- taxa_table$clade_name
    taxa_table$clade_name <- NULL
    
    metadata <- prepare_metadata('taxa')
    
    taxa_table_and_metadata <- process_taxa_table_and_metadata(taxa_table, metadata)
    taxa_table <- taxa_table_and_metadata[['taxa_table']]
    metadata <- taxa_table_and_metadata[['metadata']]
    
    taxa_sub <- data.frame(sample = colnames(taxa_table),
                           abun = unlist(c(taxa_table['k__Bacteria|p__Firmicutes|c__Clostridia|o__Eubacteriales|f__Oscillospiraceae|g__Faecalibacterium|s__Faecalibacterium_prausnitzii|t__SGB15326',])) / 100)
    
    metadata_expanded <- left_join(metadata, taxa_sub, by = 'sample')
    metadata_expanded <- metadata_expanded %>%
        mutate(diagnosis_dysbiosis = interaction(diagnosis, dysbiosis_state)) %>%
        mutate(diagnosis_dysbiosis = case_when(diagnosis_dysbiosis == 'nonIBD.none' ~ 'non-IBD',
                                               diagnosis_dysbiosis == 'UC.none' ~ 'UC',
                                               diagnosis_dysbiosis == 'CD.none' ~ 'CD',
                                               diagnosis_dysbiosis == 'nonIBD.dysbiosis_nonIBD' ~ 'non-IBD (Dysbiosis)',
                                               diagnosis_dysbiosis == 'UC.dysbiosis_UC' ~ 'UC (Dysbiosis)',
                                               diagnosis_dysbiosis == 'CD.dysbiosis_CD' ~ 'CD (Dysbiosis)',)) %>%
        mutate(diagnosis_dysbiosis = factor(diagnosis_dysbiosis, levels = 
                                                rev(c('non-IBD', 'UC', 'CD', 'non-IBD (Dysbiosis)', 'UC (Dysbiosis)', 'CD (Dysbiosis)'))))
    table(metadata_expanded$abun > 0, metadata_expanded$diagnosis_dysbiosis)
}
in_text_vals()

recreate_summary_plot_mtx <- function() {
    # Clean up default HMP2 plot
    # Rename results file with clean titles
    all_results <- read.csv('HMP2/analysis/fit_out_Maaslin3_covariate/all_results.tsv', sep='\t')
    all_results <- all_results %>%
        mutate(metadata = case_when(metadata == 'consent_age' ~ 'Age',
                                    metadata == 'Antibiotics' ~ 'Abx',
                                    metadata == 'diagnosis' ~ 'Diagnosis',
                                    metadata == 'dysbiosis_state' ~ 'Dysbiosis',
                                    metadata == 'reads_filtered' ~ 'Read depth'),
               value = case_when(value == 'dysbiosis_CD' ~ 'CD',
                                 value == 'dysbiosis_UC' ~ 'UC',
                                 value == 'dysbiosis_nonIBD' ~ 'non-IBD',
                                 value == 'Yes' ~ 'Used', # Antibiotics
                                 value == 'consent_age' ~ 'Age',
                                 value == 'reads_filtered' ~ 'Read depth',
                                 TRUE ~ value),
               feature = gsub('\\.|\\.', ' ', feature) %>%
                   gsub(pattern = '   ', replacement = ' ', feature) %>%
                   gsub(pattern = '  ', replacement = ' ', feature) %>%
                   gsub(pattern = '_', replacement = ' ') %>% 
                   gsub(pattern = 'sp ', replacement = 'sp. ') %>%
                   gsub(pattern = " g .* s  ", replacement = ', '))
    
    # Write results
    write.table(all_results, 'HMP2/analysis/fit_out_Maaslin3_covariate/all_results.tsv', sep='\t', row.names = F)
    
    all_results <- read.csv('HMP2/analysis/fit_out_Maaslin3_covariate/significant_results.tsv', sep='\t')
    all_results <- all_results %>%
        mutate(metadata = case_when(metadata == 'consent_age' ~ 'Age',
                                    metadata == 'Antibiotics' ~ 'Abx',
                                    metadata == 'diagnosis' ~ 'Diagnosis',
                                    metadata == 'dysbiosis_state' ~ 'Dysbiosis',
                                    metadata == 'reads_filtered' ~ 'Read depth'),
               value = case_when(value == 'dysbiosis_CD' ~ 'CD',
                                 value == 'dysbiosis_UC' ~ 'UC',
                                 value == 'dysbiosis_nonIBD' ~ 'non-IBD',
                                 value == 'Yes' ~ 'Used', # Antibiotics
                                 value == 'consent_age' ~ 'Age',
                                 value == 'reads_filtered' ~ 'Read depth',
                                 TRUE ~ value),
               feature = gsub('\\.|\\.', ' ', feature) %>%
                   gsub(pattern = '   ', replacement = ' ', feature) %>%
                   gsub(pattern = '  ', replacement = ' ', feature) %>%
                   gsub(pattern = '_', replacement = ' ') %>% 
                   gsub(pattern = 'sp ', replacement = 'sp. ') %>%
                   gsub(pattern = " g .* s  ", replacement = ', '))
    
    # Write results
    write.table(all_results, 'HMP2/analysis/fit_out_Maaslin3_covariate/significant_results.tsv', sep='\t', row.names = F)
    
    # Need to create param_list
    taxa_table <- read.csv('HMP2/analysis/data/metaphlan3_taxonomic_profiles.tsv',check.names = F, sep = '\t')
    
    taxa_table <- taxa_table[taxa_table$`#SampleID` == 'UNCLASSIFIED' | 
                                 grepl('\\|s__', taxa_table$`#SampleID`) & 
                                 !grepl('\\|t__', taxa_table$`#SampleID`),]
    rownames(taxa_table) <- taxa_table$`#SampleID`
    taxa_table$`#SampleID` <- NULL
    taxa_table['UNCLASSIFIED',] <- pmax(1 - colSums(taxa_table), 0)
    taxa_table <- taxa_table * 100 # Convert to percents to be consistent with v4
    
    metadata <- prepare_metadata('taxa')
    
    taxa_table_and_metadata <- process_taxa_table_and_metadata(taxa_table, metadata)
    taxa_table <- taxa_table_and_metadata[['taxa_table']]
    metadata <- taxa_table_and_metadata[['metadata']]
    
    taxa_table <- read.csv('HMP2/analysis/data/metaphlan3_taxonomic_profiles.tsv',check.names = F, sep = '\t')
    
    taxa_table <- taxa_table[taxa_table$`#SampleID` == 'UNCLASSIFIED' | 
                                 grepl('\\|s__', taxa_table$`#SampleID`) & 
                                 !grepl('\\|t__', taxa_table$`#SampleID`),]
    rownames(taxa_table) <- taxa_table$`#SampleID`
    taxa_table$`#SampleID` <- NULL
    taxa_table['UNCLASSIFIED',] <- pmax(1 - colSums(taxa_table), 0)
    taxa_table <- taxa_table * 100 # Convert to percents to be consistent with v4
    
    metadata <- metadata[metadata$sample %in% colnames(taxa_table),]
    metadata <- metadata[as.numeric(mapvalues(metadata$sample, colnames(taxa_table), 1:ncol(taxa_table))),]
    
    tmp_fit_out <- paste0("HMP2/analysis/fit_out_Maaslin3_covariate")
    
    results_in <- read.csv(file.path(tmp_fit_out, "all_results.tsv"), sep='\t')
    keep_features <- c(unique(results_in$feature[!(results_in$coef > 5 & results_in$pval_individual < 10^-12 & 
                                                       results_in$model == 'prevalence')][1:20])) # likely prevalence misfits
    results_in <- results_in[results_in$feature %in% keep_features,]
    write.table(results_in, file.path(tmp_fit_out, "all_results.tsv"), sep = '\t', row.names = F)
    
    results_in <- read.csv(file.path(tmp_fit_out, "significant_results.tsv"), sep='\t')
    results_in <- results_in[results_in$feature %in% keep_features,]
    write.table(results_in, file.path(tmp_fit_out, "significant_results.tsv"), sep = '\t', row.names = F)
    
    # Set the new heatmap and coefficient plot variables and order them
    maaslin_plot_results_from_output(metadata = metadata, 
                                     output = tmp_fit_out, 
                                     normalization = 'TSS', 
                                     transform = 'LOG', 
                                     plot_summary_plot = T, 
                                     plot_associations = F, 
                                     max_significance = 0.1, 
                                     median_comparison_abundance = F, 
                                     median_comparison_prevalence = F, 
                                     heatmap_vars = c('Dysbiosis non-IBD', 'Diagnosis CD', 'Diagnosis UC', 'Abx Used', 'Age', 'Read depth'),
                                     coef_plot_vars = c('Dysbiosis CD', 'Dysbiosis UC'))
}
recreate_summary_plot_mtx()

recreate_summary_plot_v4 <- function() {
    # Clean up default HMP2 plot
    # Rename results file with clean titles
    all_results <- read.csv('HMP2/analysis/fit_out_Maaslin3CompAdjust_4/all_results.tsv', sep='\t')
    all_results <- all_results %>%
        mutate(metadata = case_when(metadata == 'consent_age' ~ 'Age',
                                    metadata == 'Antibiotics' ~ 'Abx',
                                    metadata == 'diagnosis' ~ 'Diagnosis',
                                    metadata == 'dysbiosis_state' ~ 'Dysbiosis',
                                    metadata == 'reads_filtered' ~ 'Read depth'),
               value = case_when(value == 'dysbiosis_CD' ~ 'CD',
                                 value == 'dysbiosis_UC' ~ 'UC',
                                 value == 'dysbiosis_nonIBD' ~ 'non-IBD',
                                 value == 'Yes' ~ 'Used', # Antibiotics
                                 value == 'consent_age' ~ 'Age',
                                 value == 'reads_filtered' ~ 'Read depth',
                                 TRUE ~ value),
               feature = gsub('.*\\.s', '', feature) %>%
                   gsub(pattern = 's__', replacement = '') %>% 
                   gsub(pattern = '_', replacement = ' ') %>% 
                   gsub(pattern = 'sp ', replacement = 'sp. ') %>%
                   gsub(pattern = ' group$', replacement = '') %>%
                   gsub(pattern = '\\.t  ', replacement = ' (') %>%
                   paste0(')') %>%
                   trimws())
    
    # Write results
    write.table(all_results, 'HMP2/analysis/fit_out_Maaslin3CompAdjust_4/all_results.tsv', sep='\t', row.names = F)
    
    all_results <- read.csv('HMP2/analysis/fit_out_Maaslin3CompAdjust_4/significant_results.tsv', sep='\t')
    all_results <- all_results %>%
        mutate(metadata = case_when(metadata == 'consent_age' ~ 'Age',
                                    metadata == 'Antibiotics' ~ 'Abx',
                                    metadata == 'diagnosis' ~ 'Diagnosis',
                                    metadata == 'dysbiosis_state' ~ 'Dysbiosis',
                                    metadata == 'reads_filtered' ~ 'Read depth'),
               value = case_when(value == 'dysbiosis_CD' ~ 'CD',
                                 value == 'dysbiosis_UC' ~ 'UC',
                                 value == 'dysbiosis_nonIBD' ~ 'non-IBD',
                                 value == 'Yes' ~ 'Used', # Antibiotics
                                 value == 'consent_age' ~ 'Age',
                                 value == 'reads_filtered' ~ 'Read depth',
                                 TRUE ~ value),
               feature = gsub('.*\\.s', '', feature) %>%
                   gsub(pattern = 's__', replacement = '') %>% 
                   gsub(pattern = '_', replacement = ' ') %>% 
                   gsub(pattern = 'sp ', replacement = 'sp. ') %>%
                   gsub(pattern = ' group$', replacement = '') %>%
                   gsub(pattern = '\\.t  ', replacement = ' (') %>%
                   paste0(')') %>%
                   trimws())
    
    # Write results
    write.table(all_results, 'HMP2/analysis/fit_out_Maaslin3CompAdjust_4/significant_results.tsv', sep='\t', row.names = F)
    
    # Need to create param_list
    taxa_table <- read.csv('HMP2/analysis/data/metaphlan4_taxonomic_profiles.tsv', skip = 1, check.names = F, sep = '\t')
    
    # Reorganize taxa table
    colnames(taxa_table) <- gsub('_taxonomic$', '', colnames(taxa_table))
    taxa_table <- taxa_table[taxa_table$clade_name == 'UNCLASSIFIED' | 
                                 grepl('\\|t__', taxa_table$clade_name),]
    rownames(taxa_table) <- taxa_table$clade_name
    taxa_table$clade_name <- NULL
    
    metadata <- prepare_metadata('taxa')
    
    taxa_table_and_metadata <- process_taxa_table_and_metadata(taxa_table, metadata)
    taxa_table <- taxa_table_and_metadata[['taxa_table']]
    metadata <- taxa_table_and_metadata[['metadata']]
    
    taxa_table <- read.csv('HMP2/analysis/data/metaphlan4_taxonomic_profiles.tsv', skip = 1, check.names = F, sep = '\t')
    
    # Reorganize taxa table
    colnames(taxa_table) <- gsub('_taxonomic$', '', colnames(taxa_table))
    taxa_table <- taxa_table[taxa_table$clade_name == 'UNCLASSIFIED' | 
                                 (grepl('\\|s__', taxa_table$clade_name) &
                                      !grepl('\\|t__', taxa_table$clade_name)),]
    rownames(taxa_table) <- taxa_table$clade_name
    taxa_table$clade_name <- NULL
    
    metadata <- metadata[metadata$sample %in% colnames(taxa_table),]
    metadata <- metadata[as.numeric(mapvalues(metadata$sample, colnames(taxa_table), 1:ncol(taxa_table))),]
    
    tmp_fit_out <- paste0("HMP2/analysis/fit_out_Maaslin3CompAdjust_4/")
    
    keep_taxa <- c("Phocaeicola vulgatus (SGB1814)",
                   "Bacteroides uniformis (SGB1836)",
                   "GGB4802 SGB6641",
                   "Odoribacter splanchnicus (SGB1790)",
                   "Faecalibacterium prausnitzii (SGB15316)", 
                   "Faecalibacterium prausnitzii (SGB15318)", 
                   "Faecalibacterium prausnitzii (SGB15342)",
                   "Faecalibacterium prausnitzii (SGB15332)",
                   "Clostridium sp. AT4 (SGB4753)",
                   "Bacteroides ovatus (SGB1871)",
                   "GGB33469 SGB15236",
                   "Clostridium sp. AF34 10BH (SGB4914)",
                   "Pseudoflavonifractor gallinarum (SGB29328)",
                   "Enterocloster_citroniae (SGB4761)",
                   "Parabacteroides_distasonis (SGB1934)",
                   "Eubacterium rectale (SGB4933)",
                   "Parasutterella excrementihominis (SGB9262)",
                   "Blautia obeum (SGB4811)",
                   "Roseburia inulinivorans (SGB4940)",
                   "Roseburia hominis (SGB4936)",
                   "Enterocloster bolteae (SGB4758)", 
                   "Ruminococcus gnavus (SGB4584)", 
                   "Veillonella parvula (SGB6939)",
                   "Clostridium neonatale (SGB6169)",
                   "Klebsiella pneumoniae (SGB10115)",
                   "Ruminococcus torques (SGB4563)",
                   "Escherichia coli (SGB10068)",
                   "Akkermansia muciniphila (SGB9226)",
                   "Bacteroides fragilis (SGB1855)",
                   "Bacteroides thetaiotaomicron (SGB1861)",
                   "Dorea formicigenerans (SGB4575)")
    
    m2_signif <- read.csv('HMP2/analysis/results/v4_ibd_associations_Maaslin2.tsv', sep='\t')
    
    m2_signif <- m2_signif %>%
        mutate(metadata = case_when(metadata == 'consent_age' ~ 'Age',
                                    metadata == 'Antibiotics' ~ 'Abx',
                                    metadata == 'diagnosis' ~ 'Diagnosis',
                                    metadata == 'dysbiosis_state' ~ 'Dysbiosis',
                                    metadata == 'reads_filtered' ~ 'Read depth'),
               value = case_when(value == 'dysbiosis_CD' ~ 'CD',
                                 value == 'dysbiosis_UC' ~ 'UC',
                                 value == 'dysbiosis_nonIBD' ~ 'non-IBD',
                                 value == 'Yes' ~ 'Used', # Antibiotics
                                 value == 'consent_age' ~ 'Age',
                                 value == 'reads_filtered' ~ 'Read depth',
                                 TRUE ~ value),
               feature = gsub('.*\\.s', '', feature) %>%
                   gsub(pattern = 's__', replacement = '') %>% 
                   gsub(pattern = '_', replacement = ' ') %>% 
                   gsub(pattern = 'sp ', replacement = 'sp. ') %>%
                   gsub(pattern = ' group$', replacement = '') %>%
                   gsub(pattern = '\\.t  ', replacement = ' (') %>%
                   paste0(')') %>%
                   trimws())
    
    m2_signif <- m2_signif %>%
        filter(qval < 0.1 & 
                   metadata == 'Dysbiosis' & value == 'CD')
    
    prev_overlap <- intersect(m2_signif$feature, keep_taxa)
    
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
                                     heatmap_vars = c('Dysbiosis non-IBD', 'Diagnosis CD', 'Diagnosis UC', 'Abx Used', 'Age', 'Read depth'),
                                     coef_plot_vars = c('Dysbiosis CD', 'Dysbiosis UC'))
}
recreate_summary_plot_v4()

recreate_summary_plot_v4_under16 <- function() {
    # Clean up default HMP2 plot
    # Rename results file with clean titles
    all_results <- read.csv('HMP2/analysis/fit_out_Maaslin3_under16_4/all_results.tsv', sep='\t')
    all_results <- all_results %>%
        mutate(metadata = case_when(metadata == 'consent_age' ~ 'Age',
                                    metadata == 'Antibiotics' ~ 'Abx',
                                    metadata == 'diagnosis' ~ 'Diagnosis',
                                    metadata == 'dysbiosis_state' ~ 'Dysbiosis',
                                    metadata == 'reads_filtered' ~ 'Read depth'),
               value = case_when(value == 'dysbiosis_CD' ~ 'CD',
                                 value == 'dysbiosis_UC' ~ 'UC',
                                 value == 'dysbiosis_nonIBD' ~ 'non-IBD',
                                 value == 'Yes' ~ 'Used', # Antibiotics
                                 value == 'consent_age' ~ 'Age',
                                 value == 'reads_filtered' ~ 'Read depth',
                                 TRUE ~ value),
               feature = gsub('.*\\.s', '', feature) %>%
                   gsub(pattern = 's__', replacement = '') %>% 
                   gsub(pattern = '_', replacement = ' ') %>% 
                   gsub(pattern = 'sp ', replacement = 'sp. ') %>%
                   gsub(pattern = ' group$', replacement = '') %>%
                   gsub(pattern = '\\.t  ', replacement = ' (') %>%
                   paste0(')') %>%
                   trimws()) %>%
        mutate(feature = ifelse(grepl('GGB', feature), gsub(' \\(.*', '', feature), feature))
    
    # Write results
    write.table(all_results, 'HMP2/analysis/fit_out_Maaslin3_under16_4/all_results.tsv', sep='\t', row.names = F)
    
    all_results <- read.csv('HMP2/analysis/fit_out_Maaslin3_under16_4/significant_results.tsv', sep='\t')
    all_results <- all_results %>%
        mutate(metadata = case_when(metadata == 'consent_age' ~ 'Age',
                                    metadata == 'Antibiotics' ~ 'Abx',
                                    metadata == 'diagnosis' ~ 'Diagnosis',
                                    metadata == 'dysbiosis_state' ~ 'Dysbiosis',
                                    metadata == 'reads_filtered' ~ 'Read depth'),
               value = case_when(value == 'dysbiosis_CD' ~ 'CD',
                                 value == 'dysbiosis_UC' ~ 'UC',
                                 value == 'dysbiosis_nonIBD' ~ 'non-IBD',
                                 value == 'Yes' ~ 'Used', # Antibiotics
                                 value == 'consent_age' ~ 'Age',
                                 value == 'reads_filtered' ~ 'Read depth',
                                 TRUE ~ value),
               feature = gsub('.*\\.s', '', feature) %>%
                   gsub(pattern = 's__', replacement = '') %>% 
                   gsub(pattern = '_', replacement = ' ') %>% 
                   gsub(pattern = 'sp ', replacement = 'sp. ') %>%
                   gsub(pattern = ' group$', replacement = '') %>%
                   gsub(pattern = '\\.t  ', replacement = ' (') %>%
                   paste0(')') %>%
                   trimws()) %>%
        mutate(feature = ifelse(grepl('GGB', feature), gsub(' \\(.*', '', feature), feature))
    
    # Write results
    write.table(all_results, 'HMP2/analysis/fit_out_Maaslin3_under16_4/significant_results.tsv', sep='\t', row.names = F)
    
    # Need to create param_list
    taxa_table <- read.csv('HMP2/analysis/data/metaphlan4_taxonomic_profiles.tsv', skip = 1, check.names = F, sep = '\t')
    
    # Reorganize taxa table
    colnames(taxa_table) <- gsub('_taxonomic$', '', colnames(taxa_table))
    taxa_table <- taxa_table[taxa_table$clade_name == 'UNCLASSIFIED' | 
                                 grepl('\\|t__', taxa_table$clade_name),]
    rownames(taxa_table) <- taxa_table$clade_name
    taxa_table$clade_name <- NULL
    
    metadata <- prepare_metadata('taxa')
    
    taxa_table_and_metadata <- process_taxa_table_and_metadata(taxa_table, metadata)
    taxa_table <- taxa_table_and_metadata[['taxa_table']]
    metadata <- taxa_table_and_metadata[['metadata']]
    
    taxa_table <- read.csv('HMP2/analysis/data/metaphlan4_taxonomic_profiles.tsv', skip = 1, check.names = F, sep = '\t')
    
    # Reorganize taxa table
    colnames(taxa_table) <- gsub('_taxonomic$', '', colnames(taxa_table))
    taxa_table <- taxa_table[taxa_table$clade_name == 'UNCLASSIFIED' | 
                                 (grepl('\\|s__', taxa_table$clade_name) &
                                      !grepl('\\|t__', taxa_table$clade_name)),]
    rownames(taxa_table) <- taxa_table$clade_name
    taxa_table$clade_name <- NULL
    
    metadata <- metadata[metadata$sample %in% colnames(taxa_table),]
    metadata <- metadata[as.numeric(mapvalues(metadata$sample, colnames(taxa_table), 1:ncol(taxa_table))),]
    
    tmp_fit_out <- paste0("HMP2/analysis/fit_out_Maaslin3_under16_4/")
    
    keep_taxa <- rev(c("Clostridium sp. AT4 (SGB4753)",
                       "Pseudoflavonifractor gallinarum (SGB29328)",
                       "Klebsiella pneumoniae (SGB10115)",
                       "Hungatella hathewayi (SGB4741)",
                       "Roseburia hominis (SGB4936)",
                       "Enterocloster bolteae (SGB4758)", 
                       "Parasutterella excrementihominis (SGB9262)",
                       "Bacteroides ovatus (SGB1871)",
                       "Faecalibacterium prausnitzii (SGB15342)",
                       "Phocaeicola vulgatus (SGB1814)",
                       "Faecalibacterium prausnitzii (SGB15318)"))
    
    output = tmp_fit_out 
    normalization = 'TSS'
    transform = 'LOG'
    plot_summary_plot = T
    plot_associations = F 
    max_significance = 0.1
    median_comparison_abundance = T
    median_comparison_prevalence = F
    heatmap_vars = c('Dysbiosis non-IBD', 'Diagnosis CD', 'Diagnosis UC', 'Abx Used', 'Read depth')
    coef_plot_vars = c('Dysbiosis CD', 'Dysbiosis UC')
    
    merged_results <- read.csv(file.path(tmp_fit_out, "all_results.tsv"), sep='\t')
    merged_results$model[merged_results$model == "abundance"] <- "linear"
    merged_results$model[merged_results$model == "prevalence"] <- "logistic"
    
    merged_results <- maaslin3:::preprocess_merged_results(merged_results)
    
    median_df <- merged_results %>%
        dplyr::group_by(.data$full_metadata_name, .data$model) %>%
        dplyr::summarize(median_val = median(.data$coef), .groups = 'drop')
    
    if (!median_comparison_abundance) {
        median_df$median_val[median_df$model == 'Abundance'] <- 0
    }
    
    if (!median_comparison_prevalence) {
        median_df$median_val[median_df$model == 'Prevalence'] <- 0
    }
    
    results_in <- read.csv(file.path(tmp_fit_out, "all_results.tsv"), sep='\t')
    results_in <- results_in[results_in$feature %in% keep_taxa,]
    write.table(results_in, file.path(tmp_fit_out, "all_results.tsv"), sep = '\t', row.names = F)
    
    results_in <- read.csv(file.path(tmp_fit_out, "significant_results.tsv"), sep='\t')
    results_in <- results_in[results_in$feature %in% keep_taxa,]
    write.table(results_in, file.path(tmp_fit_out, "significant_results.tsv"), sep = '\t', row.names = F)
    
    merged_results <- read.csv(file.path(tmp_fit_out, "all_results.tsv"), sep='\t')
    merged_results$model[merged_results$model == "abundance"] <- "linear"
    merged_results$model[merged_results$model == "prevalence"] <- "logistic"
    
    merged_results <- maaslin3:::preprocess_merged_results(merged_results)
    
    merged_results <-
        merged_results[merged_results$full_metadata_name %in%
                           c(coef_plot_vars, heatmap_vars),]
    
    # Subset associations for plotting
    merged_results_joint_only <-
        unique(merged_results[, c('feature', 'qval_joint')])
    merged_results_joint_only <-
        merged_results_joint_only[
            order(merged_results_joint_only$qval_joint),]
    first_n <- 30
    if (length(unique(merged_results_joint_only$feature)) < first_n) {
        first_n <- length(unique(merged_results_joint_only$feature))
    }
    signif_taxa <-
        unique(merged_results_joint_only$feature)[seq(first_n)]
    
    merged_results_sig <- merged_results %>%
        dplyr::filter(.data$feature %in% signif_taxa)
    
    # Order features
    ord_feature <- keep_taxa
    
    merged_results_sig$feature <-
        factor(merged_results_sig$feature, levels = ord_feature)
    
    # Create coefficient plot
    if (length(coef_plot_vars) > 0 &
        sum(merged_results_sig$full_metadata_name %in% 
            coef_plot_vars) >= 1) {
        coef_plot_data <-
            merged_results_sig[merged_results_sig$full_metadata_name %in% 
                                   coef_plot_vars,]
        
        # Limit plotted coefficients to median +/- 10 times distance to quartiles
        quantile_df <- coef_plot_data %>%
            dplyr::group_by(.data$full_metadata_name) %>%
            dplyr::summarise(
                lower_q = median(.data$coef) - 10 * 
                    (median(.data$coef) - quantile(.data$coef, 0.25)),
                upper_q = median(.data$coef) + 10 * 
                    (quantile(.data$coef, 0.75) - median(.data$coef))
            ) %>%
            data.frame()
        rownames(quantile_df) <- quantile_df$full_metadata_name
        
        # Make sure insignificant coefficients don't distort the plot
        coef_plot_data <-
            coef_plot_data[coef_plot_data$qval_individual < 
                               max_significance |
                               (coef_plot_data$coef > quantile_df[
                                   coef_plot_data$full_metadata_name, 
                                   'lower_q'] &
                                    coef_plot_data$coef < quantile_df[
                                        coef_plot_data$full_metadata_name, 
                                        'upper_q']),]
        
        # Choose breaks for plot
        custom_break_fun <- function(n) {
            return(function(x) {
                extended_breaks <- scales::breaks_extended(n)(x)
                if (max(x) > 0) {
                    extended_breaks <- extended_breaks[
                        extended_breaks <= max(x) * 0.9]
                } else {
                    extended_breaks <- extended_breaks[
                        extended_breaks <= max(x) * 1.1]
                }
                if (min(x) > 0) {
                    extended_breaks <- extended_breaks[
                        extended_breaks >= min(x) * 1.1]
                } else {
                    extended_breaks <- extended_breaks[
                        extended_breaks >= min(x) * 0.9]
                }
                extended_breaks
            })
        }
        
        # Start plot
        p1 <-
            ggplot2::ggplot(coef_plot_data,
                            ggplot2::aes(x = .data$coef, y = .data$feature))
        
        # Show nulls the coefficients were compared against
        if (median_comparison_prevalence |
            median_comparison_abundance) {
            p1 <- p1 +
                ggplot2::guides(linetype = ggplot2::guide_legend(
                    title = 'Null hypothesis', order = 1),
                ) +
                ggplot2::geom_vline(
                    data = median_df[median_df$full_metadata_name %in% 
                                         coef_plot_vars,],
                    ggplot2::aes(
                        xintercept = .data$median_val,
                        linetype = .data$model
                    ),
                    color = "darkgray"
                ) +
                ggplot2::scale_linetype_manual(values = c("Prevalence" = 
                                                              "dashed", 
                                                          "Abundance" = 
                                                              "solid"))
        } else {
            p1 <- p1 +
                ggplot2::geom_vline(
                    ggplot2::aes(xintercept = 0),
                    color = "darkgray",
                    linetype = 'dashed'
                )
        }
        
        # Q-value color scale
        scale_fill_gradient_limits <-
            c(min(max_significance, 10 ^ -20), 1)
        if (min(coef_plot_data$qval_individual) < max_significance) {
            scale_fill_gradient_breaks <-
                c(10 ^ -20, max_significance, 1)
        } else {
            scale_fill_gradient_breaks <- c(max_significance, 1)
        }
        if (min(coef_plot_data$qval_individual) < max_significance) {
            scale_fill_gradient_labels <-
                c(paste0("1e", -20),
                  paste0(max_significance),
                  "1")
        } else {
            scale_fill_gradient_labels <- c(paste0(max_significance),
                                            "1")
        }
        
        # Create the whole plot
        p1 <- p1 +
            ggplot2::geom_errorbar(
                ggplot2::aes(
                    xmin = .data$coef - .data$stderr,
                    xmax = .data$coef + .data$stderr
                ),
                width = 0.2
            ) +
            ggplot2::geom_point(
                data = coef_plot_data[coef_plot_data$model == 
                                          'Prevalence',],
                ggplot2::aes(
                    shape = .data$model,
                    fill = .data$qval_individual
                ),
                size = 4.5,
                color = "black"
            ) +
            ggplot2::scale_fill_gradient(
                low = "#008B8B",
                high = "white",
                limits = scale_fill_gradient_limits,
                breaks = scale_fill_gradient_breaks,
                labels = scale_fill_gradient_labels,
                transform = scales::pseudo_log_trans(sigma = 0.001),
                name = bquote("Prevalence" ~ P["FDR"])
            ) +
            ggnewscale::new_scale_fill() +
            ggplot2::geom_point(
                data = coef_plot_data[coef_plot_data$model == 'Abundance',],
                ggplot2::aes(
                    shape = .data$model,
                    fill = .data$qval_individual
                ),
                size = 4.5,
                color = "black"
            ) +
            ggplot2::scale_fill_gradient(
                low = "#8B008B",
                high = "white",
                limits = scale_fill_gradient_limits,
                breaks = scale_fill_gradient_breaks,
                labels = scale_fill_gradient_labels,
                transform = scales::pseudo_log_trans(sigma = 0.001),
                name = bquote("Abundance" ~ P["FDR"])
            ) +
            ggplot2::scale_x_continuous(
                breaks = custom_break_fun(n = 6),
                limits = c(
                    min(coef_plot_data$coef) - 
                        quantile(coef_plot_data$stderr, 0.8),
                    max(coef_plot_data$coef) + 
                        quantile(coef_plot_data$stderr, 0.8)
                )
            ) +
            ggplot2::scale_shape_manual(name = "Association", values =
                                            c(21, 24)) +
            ggplot2::guides(shape = ggplot2::guide_legend(order = 2), ) +
            ggplot2::labs(x = expression(paste(beta, " coefficient")),  
                          y = "Feature") +
            ggplot2::theme_bw() +
            ggplot2::theme(
                axis.title = ggplot2::element_text(size = 16),
                axis.text.x = ggplot2::element_text(size = 14),
                axis.text.y = ggplot2::element_text(size = 14),
                legend.title = ggplot2::element_text(size = 16),
                legend.text = ggplot2::element_text(size = 14, 
                                                    face = "plain"),
                legend.position = "right",
                legend.background = ggplot2::element_rect(
                    fill = "transparent"),
                panel.spacing = ggplot2::unit(0, "lines"),
                panel.grid.minor = ggplot2::element_blank(),
                strip.text = ggplot2::element_text(size = 14),
                strip.background = ggplot2::element_rect(
                    fill = "transparent")
            ) +
            ggplot2::facet_wrap(
                ~ factor(full_metadata_name, 
                         levels = unique(coef_plot_vars)),
                scales = 'free_x',
                ncol = length(coef_plot_vars)
            )
    } else {
        p1 <- NULL
    }
    
    # Create heatmap plot
    if (length(heatmap_vars) > 0 &
        sum(merged_results_sig$full_metadata_name %in% heatmap_vars) >= 1) {
        
        merged_results_sig$sig_star <- cut(merged_results_sig$qval_individual, 
                                           breaks = c(-Inf, max_significance/10, max_significance, 
                                                      Inf), label = c("**", "*", ""))
        coefficient_thresh <- 2
        coef_breaks <- c(-coefficient_thresh, -coefficient_thresh/2, 
                         0, coefficient_thresh/2, coefficient_thresh, Inf)
        threshold_set <- c(paste0("(-Inf,", -1 * coefficient_thresh, "]"), 
                           paste0("(", -1 * coefficient_thresh, ",", -1/2 * coefficient_thresh, "]"), 
                           paste0("(", -1/2 * coefficient_thresh, ",0]"), 
                           paste0("(0,", 1/2 * coefficient_thresh, "]"), 
                           paste0("(", 1/2 * coefficient_thresh, ",", 1 * coefficient_thresh, "]"), 
                           paste0("(", 1 * coefficient_thresh, ",Inf)"))
        threshold_indices <- vapply(merged_results_sig$coef, function(value) {
            which(value < coef_breaks)[1]
        }, FUN.VALUE = 0)
        merged_results_sig <- merged_results_sig %>% dplyr::mutate(coef_cat = threshold_set[threshold_indices])
        merged_results_sig$coef_cat <- factor(merged_results_sig$coef_cat, 
                                              levels = threshold_set)
        scale_fill_values <- rev((RColorBrewer::brewer.pal(n = 6, 
                                                           name = "RdBu")))
        names(scale_fill_values) <- threshold_set
        heatmap_data <- merged_results_sig[merged_results_sig$full_metadata_name %in% 
                                               heatmap_vars, ]
        grid <- expand.grid(feature = unique(heatmap_data$feature), 
                            full_metadata_name = unique(heatmap_data$full_metadata_name), 
                            model = unique(heatmap_data$model))
        grid$model <- factor(grid$model, levels = c('Prevalence', 'Abundance'))
        heatmap_data <- merge(grid, heatmap_data, by = c("feature", 
                                                         "full_metadata_name", "model"), all.x = TRUE)
        heatmap_data$coef[is.na(heatmap_data$coef)] <- NA
        p2 <- ggplot2::ggplot(heatmap_data, ggplot2::aes(x = factor(.data$full_metadata_name, 
                                                                    unique(heatmap_vars)), y = .data$feature)) + 
            ggplot2::geom_tile(data = heatmap_data, 
                               ggplot2::aes(fill = .data$coef_cat), colour = "white", 
                               linewidth = 0.2) + 
            ggplot2::scale_fill_manual(name = "Beta coefficient", 
                                       na.value = "#EEEEEE", values = scale_fill_values) + 
            ggplot2::geom_text(ggplot2::aes(label = .data$sig_star, 
                                            color = .data$sig_star), size = 6, vjust = 0.75, hjust = 0.5, 
                               key_glyph = ggplot2::draw_key_blank) + 
            ggplot2::scale_color_manual(name = bquote("Covariates" ~ 
                                                          P["FDR"]), breaks = c("*", "**", ""), values = c("black", 
                                                                                                           "black", "black"), labels = c(paste0("* < ", round(max_significance, 
                                                                                                                                                              3)), paste0("** < ", round(max_significance/10, 5)), 
                                                                                                                                         ""), drop = FALSE) + 
            ggplot2::labs(x = "", y = "Feature", 
                          caption = "") + 
            ggplot2::theme_bw() + 
            ggplot2::theme(axis.title = ggplot2::element_text(size = 16), 
                           axis.text.x = ggplot2::element_text(size = 14, angle = 90, 
                                                               vjust = 0.5, hjust = 1), legend.title = ggplot2::element_text(size = 16), 
                           legend.text = ggplot2::element_text(size = 14, face = "plain"), 
                           legend.position = "bottom", legend.background = ggplot2::element_rect(fill = "transparent"), 
                           panel.spacing = ggplot2::unit(0, "lines"), panel.grid.minor = ggplot2::element_blank(), 
                           strip.text = ggplot2::element_text(size = 14), strip.background = ggplot2::element_rect(fill = "transparent")) + 
            ggplot2::guides(fill = ggplot2::guide_legend(order = 1), 
                            color = ggplot2::guide_legend(order = 2), ) + ggplot2::facet_grid(~model, 
                                                                                              labeller = ggplot2::labeller(model = c(abundance = "Abundance", 
                                                                                                                                     prevalence = "Prevalence")))
        
        if (!is.null(p1)) {
            p2 <- p2 + ggplot2::theme(
                axis.text.y = ggplot2::element_blank(),
                axis.title.y = ggplot2::element_blank(),
                axis.ticks.y = ggplot2::element_blank(),
            )
        }
        
    } else {
        p2 <- NULL
    }
    
    # Combine plots
    if (!is.null(p1) & !is.null(p2)) {
        final_plot <- patchwork::wrap_plots(
            p1,
            p2,
            ncol = 3,
            widths = c(max(
                0, length(coef_plot_vars) * (max(15, max(
                    nchar(as.character(coef_plot_vars))
                ))) / 15 - 2
            ) + 2,
            max(0, length(
                heatmap_vars
            ) / 4 - 2) + 2,
            0.5),
            guides = 'collect'
        ) + patchwork::plot_layout(guides = "collect") & theme(legend.position = 'right', legend.direction = 'vertical')
        
    } else if (is.null(p1) & !is.null(p2)) {
        final_plot <- p2
    } else if (!is.null(p1) & is.null(p2)) {
        final_plot <- p1
    } else {
        final_plot <- NULL
    }
    
    summary_plot_file <- 'HMP2/analysis/fit_out_Maaslin3_under16_4/figures/summary_plot.pdf'
    figures_folder <- 'HMP2/analysis/fit_out_Maaslin3_under16_4/figures'
    
    # Save plot
    if (!is.null(final_plot)) {
        height_out <-
            5 + max(first_n / 5 - 5, 0) + max(nchar(c(
                as.character(coef_plot_vars),
                as.character(heatmap_vars)
            ))) / 10
        width_out <-
            6.5 + max(nchar(merged_results$feature)) / 12 +
            (length(coef_plot_vars) * (max(20, max(
                nchar(as.character(coef_plot_vars))
            ))) / 20) * 2.5 +
            length(heatmap_vars) * 0.25
        
        ggplot2::ggsave(
            summary_plot_file,
            plot = final_plot,
            height = height_out,
            width = width_out
        )
        png_file <-
            file.path(figures_folder, "summary_plot.png")
        ggplot2::ggsave(png_file,
                        plot = final_plot,
                        height = height_out,
                        width = width_out)
    }
}
recreate_summary_plot_v4_under16()

recreate_summary_plot_v4_atleast16 <- function() {
    # Clean up default HMP2 plot
    # Rename results file with clean titles
    all_results <- read.csv('HMP2/analysis/fit_out_Maaslin3_atleast16_4/all_results.tsv', sep='\t')
    all_results <- all_results %>%
        mutate(metadata = case_when(metadata == 'consent_age' ~ 'Age',
                                    metadata == 'Antibiotics' ~ 'Abx',
                                    metadata == 'diagnosis' ~ 'Diagnosis',
                                    metadata == 'dysbiosis_state' ~ 'Dysbiosis',
                                    metadata == 'reads_filtered' ~ 'Read depth'),
               value = case_when(value == 'dysbiosis_CD' ~ 'CD',
                                 value == 'dysbiosis_UC' ~ 'UC',
                                 value == 'dysbiosis_nonIBD' ~ 'non-IBD',
                                 value == 'Yes' ~ 'Used', # Antibiotics
                                 value == 'consent_age' ~ 'Age',
                                 value == 'reads_filtered' ~ 'Read depth',
                                 TRUE ~ value),
               feature = gsub('.*\\.s', '', feature) %>%
                   gsub(pattern = 's__', replacement = '') %>% 
                   gsub(pattern = '_', replacement = ' ') %>% 
                   gsub(pattern = 'sp ', replacement = 'sp. ') %>%
                   gsub(pattern = ' group$', replacement = '') %>%
                   gsub(pattern = '\\.t  ', replacement = ' (') %>%
                   paste0(')') %>%
                   trimws()) %>%
        mutate(feature = ifelse(grepl('GGB', feature), gsub(' \\(.*', '', feature), feature))
    
    # Write results
    write.table(all_results, 'HMP2/analysis/fit_out_Maaslin3_atleast16_4/all_results.tsv', sep='\t', row.names = F)
    
    all_results <- read.csv('HMP2/analysis/fit_out_Maaslin3_atleast16_4/significant_results.tsv', sep='\t')
    all_results <- all_results %>%
        mutate(metadata = case_when(metadata == 'consent_age' ~ 'Age',
                                    metadata == 'Antibiotics' ~ 'Abx',
                                    metadata == 'diagnosis' ~ 'Diagnosis',
                                    metadata == 'dysbiosis_state' ~ 'Dysbiosis',
                                    metadata == 'reads_filtered' ~ 'Read depth'),
               value = case_when(value == 'dysbiosis_CD' ~ 'CD',
                                 value == 'dysbiosis_UC' ~ 'UC',
                                 value == 'dysbiosis_nonIBD' ~ 'non-IBD',
                                 value == 'Yes' ~ 'Used', # Antibiotics
                                 value == 'consent_age' ~ 'Age',
                                 value == 'reads_filtered' ~ 'Read depth',
                                 TRUE ~ value),
               feature = gsub('.*\\.s', '', feature) %>%
                   gsub(pattern = 's__', replacement = '') %>% 
                   gsub(pattern = '_', replacement = ' ') %>% 
                   gsub(pattern = 'sp ', replacement = 'sp. ') %>%
                   gsub(pattern = ' group$', replacement = '') %>%
                   gsub(pattern = '\\.t  ', replacement = ' (') %>%
                   paste0(')') %>%
                   trimws()) %>%
        mutate(feature = ifelse(grepl('GGB', feature), gsub(' \\(.*', '', feature), feature))
    
    # Write results
    write.table(all_results, 'HMP2/analysis/fit_out_Maaslin3_atleast16_4/significant_results.tsv', sep='\t', row.names = F)
    
    # Need to create param_list
    taxa_table <- read.csv('HMP2/analysis/data/metaphlan4_taxonomic_profiles.tsv', skip = 1, check.names = F, sep = '\t')
    
    # Reorganize taxa table
    colnames(taxa_table) <- gsub('_taxonomic$', '', colnames(taxa_table))
    taxa_table <- taxa_table[taxa_table$clade_name == 'UNCLASSIFIED' | 
                                 grepl('\\|t__', taxa_table$clade_name),]
    rownames(taxa_table) <- taxa_table$clade_name
    taxa_table$clade_name <- NULL
    
    metadata <- prepare_metadata('taxa')
    
    taxa_table_and_metadata <- process_taxa_table_and_metadata(taxa_table, metadata)
    taxa_table <- taxa_table_and_metadata[['taxa_table']]
    metadata <- taxa_table_and_metadata[['metadata']]
    
    taxa_table <- read.csv('HMP2/analysis/data/metaphlan4_taxonomic_profiles.tsv', skip = 1, check.names = F, sep = '\t')
    
    # Reorganize taxa table
    colnames(taxa_table) <- gsub('_taxonomic$', '', colnames(taxa_table))
    taxa_table <- taxa_table[taxa_table$clade_name == 'UNCLASSIFIED' | 
                                 (grepl('\\|s__', taxa_table$clade_name) &
                                      !grepl('\\|t__', taxa_table$clade_name)),]
    rownames(taxa_table) <- taxa_table$clade_name
    taxa_table$clade_name <- NULL
    
    metadata <- metadata[metadata$sample %in% colnames(taxa_table),]
    metadata <- metadata[as.numeric(mapvalues(metadata$sample, colnames(taxa_table), 1:ncol(taxa_table))),]
    
    tmp_fit_out <- paste0("HMP2/analysis/fit_out_Maaslin3_atleast16_4/")
    
    keep_taxa <- rev(c("Clostridium sp. AT4 (SGB4753)",
                       "Pseudoflavonifractor gallinarum (SGB29328)",
                       "Klebsiella pneumoniae (SGB10115)",
                       "Hungatella hathewayi (SGB4741)",
                       "Roseburia hominis (SGB4936)",
                       "Enterocloster bolteae (SGB4758)", 
                       "Parasutterella excrementihominis (SGB9262)",
                       "Bacteroides ovatus (SGB1871)",
                       "Faecalibacterium prausnitzii (SGB15342)",
                       "Phocaeicola vulgatus (SGB1814)",
                       "Faecalibacterium prausnitzii (SGB15318)"))
    
    output = tmp_fit_out 
    normalization = 'TSS'
    transform = 'LOG'
    plot_summary_plot = T
    plot_associations = F 
    max_significance = 0.1
    median_comparison_abundance = T
    median_comparison_prevalence = F
    heatmap_vars = c('Dysbiosis non-IBD', 'Diagnosis CD', 'Diagnosis UC', 'Abx Used', 'Read depth')
    coef_plot_vars = c('Dysbiosis CD', 'Dysbiosis UC')
    
    merged_results <- read.csv(file.path(tmp_fit_out, "all_results.tsv"), sep='\t')
    merged_results$model[merged_results$model == "abundance"] <- "linear"
    merged_results$model[merged_results$model == "prevalence"] <- "logistic"
    
    merged_results <- maaslin3:::preprocess_merged_results(merged_results)
    
    median_df <- merged_results %>%
        dplyr::group_by(.data$full_metadata_name, .data$model) %>%
        dplyr::summarize(median_val = median(.data$coef), .groups = 'drop')
    
    if (!median_comparison_abundance) {
        median_df$median_val[median_df$model == 'Abundance'] <- 0
    }
    
    if (!median_comparison_prevalence) {
        median_df$median_val[median_df$model == 'Prevalence'] <- 0
    }
    
    results_in <- read.csv(file.path(tmp_fit_out, "all_results.tsv"), sep='\t')
    results_in <- results_in[results_in$feature %in% keep_taxa,]
    write.table(results_in, file.path(tmp_fit_out, "all_results.tsv"), sep = '\t', row.names = F)
    
    results_in <- read.csv(file.path(tmp_fit_out, "significant_results.tsv"), sep='\t')
    results_in <- results_in[results_in$feature %in% keep_taxa,]
    write.table(results_in, file.path(tmp_fit_out, "significant_results.tsv"), sep = '\t', row.names = F)
    
    merged_results <- read.csv(file.path(tmp_fit_out, "all_results.tsv"), sep='\t')
    merged_results$model[merged_results$model == "abundance"] <- "linear"
    merged_results$model[merged_results$model == "prevalence"] <- "logistic"
    
    merged_results <- maaslin3:::preprocess_merged_results(merged_results)
    
    merged_results <-
        merged_results[merged_results$full_metadata_name %in%
                           c(coef_plot_vars, heatmap_vars),]
    
    # Subset associations for plotting
    merged_results_joint_only <-
        unique(merged_results[, c('feature', 'qval_joint')])
    merged_results_joint_only <-
        merged_results_joint_only[
            order(merged_results_joint_only$qval_joint),]
    first_n <- 30
    if (length(unique(merged_results_joint_only$feature)) < first_n) {
        first_n <- length(unique(merged_results_joint_only$feature))
    }
    signif_taxa <-
        unique(merged_results_joint_only$feature)[seq(first_n)]
    
    merged_results_sig <- merged_results %>%
        dplyr::filter(.data$feature %in% signif_taxa)
    
    # Order features
    ord_feature <- keep_taxa
    
    merged_results_sig$feature <-
        factor(merged_results_sig$feature, levels = ord_feature)
    
    # Create coefficient plot
    if (length(coef_plot_vars) > 0 &
        sum(merged_results_sig$full_metadata_name %in% 
            coef_plot_vars) >= 1) {
        coef_plot_data <-
            merged_results_sig[merged_results_sig$full_metadata_name %in% 
                                   coef_plot_vars,]
        
        # Limit plotted coefficients to median +/- 10 times distance to quartiles
        quantile_df <- coef_plot_data %>%
            dplyr::group_by(.data$full_metadata_name) %>%
            dplyr::summarise(
                lower_q = median(.data$coef) - 10 * 
                    (median(.data$coef) - quantile(.data$coef, 0.25)),
                upper_q = median(.data$coef) + 10 * 
                    (quantile(.data$coef, 0.75) - median(.data$coef))
            ) %>%
            data.frame()
        rownames(quantile_df) <- quantile_df$full_metadata_name
        
        # Make sure insignificant coefficients don't distort the plot
        coef_plot_data <-
            coef_plot_data[coef_plot_data$qval_individual < 
                               max_significance |
                               (coef_plot_data$coef > quantile_df[
                                   coef_plot_data$full_metadata_name, 
                                   'lower_q'] &
                                    coef_plot_data$coef < quantile_df[
                                        coef_plot_data$full_metadata_name, 
                                        'upper_q']),]
        
        # Choose breaks for plot
        custom_break_fun <- function(n) {
            return(function(x) {
                extended_breaks <- scales::breaks_extended(n)(x)
                if (max(x) > 0) {
                    extended_breaks <- extended_breaks[
                        extended_breaks <= max(x) * 0.9]
                } else {
                    extended_breaks <- extended_breaks[
                        extended_breaks <= max(x) * 1.1]
                }
                if (min(x) > 0) {
                    extended_breaks <- extended_breaks[
                        extended_breaks >= min(x) * 1.1]
                } else {
                    extended_breaks <- extended_breaks[
                        extended_breaks >= min(x) * 0.9]
                }
                extended_breaks
            })
        }
        
        # Start plot
        p1 <-
            ggplot2::ggplot(coef_plot_data,
                            ggplot2::aes(x = .data$coef, y = .data$feature))
        
        # Show nulls the coefficients were compared against
        if (median_comparison_prevalence |
            median_comparison_abundance) {
            p1 <- p1 +
                ggplot2::guides(linetype = ggplot2::guide_legend(
                    title = 'Null hypothesis', order = 1),
                ) +
                ggplot2::geom_vline(
                    data = median_df[median_df$full_metadata_name %in% 
                                         coef_plot_vars,],
                    ggplot2::aes(
                        xintercept = .data$median_val,
                        linetype = .data$model
                    ),
                    color = "darkgray"
                ) +
                ggplot2::scale_linetype_manual(values = c("Prevalence" = 
                                                              "dashed", 
                                                          "Abundance" = 
                                                              "solid"))
        } else {
            p1 <- p1 +
                ggplot2::geom_vline(
                    ggplot2::aes(xintercept = 0),
                    color = "darkgray",
                    linetype = 'dashed'
                )
        }
        
        # Q-value color scale
        scale_fill_gradient_limits <-
            c(min(max_significance, 10 ^ -20), 1)
        if (min(coef_plot_data$qval_individual) < max_significance) {
            scale_fill_gradient_breaks <-
                c(10 ^ -20, max_significance, 1)
        } else {
            scale_fill_gradient_breaks <- c(max_significance, 1)
        }
        if (min(coef_plot_data$qval_individual) < max_significance) {
            scale_fill_gradient_labels <-
                c(paste0("1e", -20),
                  paste0(max_significance),
                  "1")
        } else {
            scale_fill_gradient_labels <- c(paste0(max_significance),
                                            "1")
        }
        
        # Create the whole plot
        p1 <- p1 +
            ggplot2::geom_errorbar(
                ggplot2::aes(
                    xmin = .data$coef - .data$stderr,
                    xmax = .data$coef + .data$stderr
                ),
                width = 0.2
            ) +
            ggplot2::geom_point(
                data = coef_plot_data[coef_plot_data$model == 
                                          'Prevalence',],
                ggplot2::aes(
                    shape = .data$model,
                    fill = .data$qval_individual
                ),
                size = 4.5,
                color = "black"
            ) +
            ggplot2::scale_fill_gradient(
                low = "#008B8B",
                high = "white",
                limits = scale_fill_gradient_limits,
                breaks = scale_fill_gradient_breaks,
                labels = scale_fill_gradient_labels,
                transform = scales::pseudo_log_trans(sigma = 0.001),
                name = bquote("Prevalence" ~ P["FDR"])
            ) +
            ggnewscale::new_scale_fill() +
            ggplot2::geom_point(
                data = coef_plot_data[coef_plot_data$model == 'Abundance',],
                ggplot2::aes(
                    shape = .data$model,
                    fill = .data$qval_individual
                ),
                size = 4.5,
                color = "black"
            ) +
            ggplot2::scale_fill_gradient(
                low = "#8B008B",
                high = "white",
                limits = scale_fill_gradient_limits,
                breaks = scale_fill_gradient_breaks,
                labels = scale_fill_gradient_labels,
                transform = scales::pseudo_log_trans(sigma = 0.001),
                name = bquote("Abundance" ~ P["FDR"])
            ) +
            ggplot2::scale_x_continuous(
                breaks = custom_break_fun(n = 6),
                limits = c(
                    min(coef_plot_data$coef) - 
                        quantile(coef_plot_data$stderr, 0.8),
                    max(coef_plot_data$coef) + 
                        quantile(coef_plot_data$stderr, 0.8)
                )
            ) +
            ggplot2::scale_shape_manual(name = "Association", values =
                                            c(21, 24)) +
            ggplot2::guides(shape = ggplot2::guide_legend(order = 2), ) +
            ggplot2::labs(x = expression(paste(beta, " coefficient")),  
                          y = "Feature") +
            ggplot2::theme_bw() +
            ggplot2::theme(
                axis.title = ggplot2::element_text(size = 16),
                axis.text.x = ggplot2::element_text(size = 14),
                axis.text.y = ggplot2::element_text(size = 14),
                legend.title = ggplot2::element_text(size = 16),
                legend.text = ggplot2::element_text(size = 14, 
                                                    face = "plain"),
                legend.position = "bottom",
                legend.background = ggplot2::element_rect(
                    fill = "transparent"),
                panel.spacing = ggplot2::unit(0, "lines"),
                panel.grid.minor = ggplot2::element_blank(),
                strip.text = ggplot2::element_text(size = 14),
                strip.background = ggplot2::element_rect(
                    fill = "transparent")
            ) +
            ggplot2::facet_wrap(
                ~ factor(full_metadata_name, 
                         levels = unique(coef_plot_vars)),
                scales = 'free_x',
                ncol = length(coef_plot_vars)
            )
    } else {
        p1 <- NULL
    }
    
    # Create heatmap plot
    if (length(heatmap_vars) > 0 &
        sum(merged_results_sig$full_metadata_name %in% heatmap_vars) >= 1) {
        
        merged_results_sig$sig_star <- cut(merged_results_sig$qval_individual, 
                                           breaks = c(-Inf, max_significance/10, max_significance, 
                                                      Inf), label = c("**", "*", ""))
        coefficient_thresh <- 2
        coef_breaks <- c(-coefficient_thresh, -coefficient_thresh/2, 
                         0, coefficient_thresh/2, coefficient_thresh, Inf)
        threshold_set <- c(paste0("(-Inf,", -1 * coefficient_thresh, "]"), 
                           paste0("(", -1 * coefficient_thresh, ",", -1/2 * coefficient_thresh, "]"), 
                           paste0("(", -1/2 * coefficient_thresh, ",0]"), 
                           paste0("(0,", 1/2 * coefficient_thresh, "]"), 
                           paste0("(", 1/2 * coefficient_thresh, ",", 1 * coefficient_thresh, "]"), 
                           paste0("(", 1 * coefficient_thresh, ",Inf)"))
        threshold_indices <- vapply(merged_results_sig$coef, function(value) {
            which(value < coef_breaks)[1]
        }, FUN.VALUE = 0)
        merged_results_sig <- merged_results_sig %>% dplyr::mutate(coef_cat = threshold_set[threshold_indices])
        merged_results_sig$coef_cat <- factor(merged_results_sig$coef_cat, 
                                              levels = threshold_set)
        scale_fill_values <- rev((RColorBrewer::brewer.pal(n = 6, 
                                                           name = "RdBu")))
        names(scale_fill_values) <- threshold_set
        heatmap_data <- merged_results_sig[merged_results_sig$full_metadata_name %in% 
                                               heatmap_vars, ]
        grid <- expand.grid(feature = unique(heatmap_data$feature), 
                            full_metadata_name = unique(heatmap_data$full_metadata_name), 
                            model = unique(heatmap_data$model))
        grid$model <- factor(grid$model, levels = c('Prevalence', 'Abundance'))
        heatmap_data <- merge(grid, heatmap_data, by = c("feature", 
                                                         "full_metadata_name", "model"), all.x = TRUE)
        heatmap_data$coef[is.na(heatmap_data$coef)] <- NA
        p2 <- ggplot2::ggplot(heatmap_data, ggplot2::aes(x = factor(.data$full_metadata_name, 
                                                                    unique(heatmap_vars)), y = .data$feature)) + 
            ggplot2::geom_tile(data = heatmap_data, 
                               ggplot2::aes(fill = .data$coef_cat), colour = "white", 
                               linewidth = 0.2) + ggplot2::scale_fill_manual(name = "Beta coefficient",
                                                                             na.value = "#EEEEEE", values = scale_fill_values) + 
            ggplot2::geom_text(ggplot2::aes(label = .data$sig_star, color = .data$sig_star), size = 6, vjust = 0.75, hjust = 0.5,
                               key_glyph = ggplot2::draw_key_blank) + 
            ggplot2::scale_color_manual(name = bquote("Covariates" ~ P["FDR"]), 
                                        breaks = c("*", "**", ""), values = c("black", "black", "black"), 
                                        labels = c(paste0("* < ", round(max_significance, 3)), 
                                                   paste0("** < ", round(max_significance/10, 5)), 
                                                   ""), drop = FALSE) + ggplot2::labs(x = "", y = "Feature", 
                                                                                      caption = "") + 
            ggplot2::theme_bw() + 
            ggplot2::theme(axis.title = ggplot2::element_text(size = 16), 
                           axis.text.x = ggplot2::element_text(size = 14, angle = 90, vjust = 0.5, hjust = 1), 
                           legend.title = ggplot2::element_text(size = 16), 
                           legend.text = ggplot2::element_text(size = 14, face = "plain"), 
                           legend.position = "bottom", 
                           legend.background = ggplot2::element_rect(fill = "transparent"), 
                           panel.spacing = ggplot2::unit(0, "lines"), panel.grid.minor = ggplot2::element_blank(), 
                           strip.text = ggplot2::element_text(size = 14), strip.background = ggplot2::element_rect(fill = "transparent")) + 
            ggplot2::guides(fill = ggplot2::guide_legend(order = 1), 
                            color = ggplot2::guide_legend(order = 2), ) + 
            ggplot2::facet_grid(~model, 
                                labeller = ggplot2::labeller(model = c(abundance = "Abundance", 
                                                                       prevalence = "Prevalence")))
        
        if (!is.null(p1)) {
            p2 <- p2 + ggplot2::theme(
                axis.text.y = ggplot2::element_blank(),
                axis.title.y = ggplot2::element_blank(),
                axis.ticks.y = ggplot2::element_blank(),
            )
        }
        
    } else {
        p2 <- NULL
    }
    
    # Combine plots
    if (!is.null(p1) & !is.null(p2)) {
        final_plot <- patchwork::wrap_plots(
            p1,
            p2,
            ncol = 3,
            widths = c(max(
                0, length(coef_plot_vars) * (max(15, max(
                    nchar(as.character(coef_plot_vars))
                ))) / 15 - 2
            ) + 2,
            max(0, length(
                heatmap_vars
            ) / 4 - 2) + 2,
            0.5),
            guides = 'collect'
        ) + patchwork::plot_layout(guides = "collect") & theme(legend.position = 'right', legend.direction = 'vertical')
    } else if (is.null(p1) & !is.null(p2)) {
        final_plot <- p2
    } else if (!is.null(p1) & is.null(p2)) {
        final_plot <- p1
    } else {
        final_plot <- NULL
    }
    
    summary_plot_file <- 'HMP2/analysis/fit_out_Maaslin3_atleast16_4/figures/summary_plot.pdf'
    figures_folder <- 'HMP2/analysis/fit_out_Maaslin3_atleast16_4/figures'
    
    # Save plot
    if (!is.null(final_plot)) {
        height_out <-
            5 + max(first_n / 5 - 5, 0) + max(nchar(c(
                as.character(coef_plot_vars),
                as.character(heatmap_vars)
            ))) / 10
        width_out <-
            6.5 + max(nchar(merged_results$feature)) / 12 +
            (length(coef_plot_vars) * (max(20, max(
                nchar(as.character(coef_plot_vars))
            ))) / 20) * 2.5 +
            length(heatmap_vars) * 0.25
        
        ggplot2::ggsave(
            summary_plot_file,
            plot = final_plot,
            height = height_out,
            width = width_out
        )
        png_file <-
            file.path(figures_folder, "summary_plot.png")
        ggplot2::ggsave(png_file,
                        plot = final_plot,
                        height = height_out,
                        width = width_out)
    }
}
recreate_summary_plot_v4_atleast16()






