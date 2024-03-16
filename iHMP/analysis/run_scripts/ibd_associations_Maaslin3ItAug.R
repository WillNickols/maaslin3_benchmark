#!/usr/bin/env Rscript
package_vec = c("reshape2", "ggplot2", "optparse", 
                "parallel", "stringi", "doParallel", "plyr", "tidyr", 'dplyr', 'vegan')
invisible(suppressPackageStartupMessages(lapply(package_vec, require, character.only = TRUE)))

# Command Line Usage
option_list = list(
  make_option(
    c("-c", "--nCores"), default=4, # c stands for how many cores to be used in the analysis
    type = "integer"),
  make_option(
    c("-w", "--workingDirectory"), # w stands for working directory
    type = "character"),
  make_option(
    c("-a", "--analysisDirectory"), # a stands for analysis directory
    type = "character"),
  make_option(
    c("-d", "--dataset"), # d stands for dataset
    type = "character"))
opt <- parse_args(OptionParser(option_list=option_list), positional_arguments = TRUE)
workingDirectory <- opt$options$workingDirectory
analysisDirectory <- opt$options$analysisDirectory
nCores<- opt$options$nCores

# Read in data
if (opt$options$dataset == 'taxa') {
  taxa_table <- read.csv('data/metaphlan4_taxonomic_profiles.tsv', skip = 1, check.names = F, sep = '\t')
  
  # Reorganize taxa table
  colnames(taxa_table) <- gsub('_taxonomic$', '', colnames(taxa_table))
  taxa_table <- taxa_table[taxa_table$clade_name == 'UNCLASSIFIED' | 
                             (grepl('\\|s__', taxa_table$clade_name) &
                                !grepl('\\|t__', taxa_table$clade_name)),]
  rownames(taxa_table) <- taxa_table$clade_name
  taxa_table$clade_name <- NULL
  
} else {
  Maaslin3_path <- paste0(gsub("/$", "", workingDirectory), "/Maaslin3/R/")
  for (R_file in dir(Maaslin3_path, pattern = "*.R$")) {
    source(file.path(Maaslin3_path, R_file))
  }
  
  # Read in data
  mbx_table <- read.csv("data/intensities_hmp2.csv")
  annotations <- read.csv("data/annotations_hmp2.csv")
  
  if (all(abs(annotations$prev - rowMeans(!is.na(mbx_table))) < 0.001)) {
    mbx_table <- mbx_table[annotations$prim_feature == 'primary' & annotations$Metabolite != '',]
    rownames(mbx_table) <- paste0(annotations$Metabolite[annotations$prim_feature == 'primary' & annotations$Metabolite != ''], '_',
                                  annotations$HMDB.ID[annotations$prim_feature == 'primary' & annotations$Metabolite != ''], '_',
                                  annotations$Method[annotations$prim_feature == 'primary' & annotations$Metabolite != ''])
    mbx_table[is.na(mbx_table)] <- 0
  }
  taxa_table <- mbx_table
}
metadata <- read.csv('data/hmp2_metadata_2018-08-20.csv', check.names = F)

# Reorganize metadata table
for (participant_id in unique(metadata$`Participant ID`)) {
  bmis <- metadata$BMI[metadata$`Participant ID` == participant_id]
  bmis <- bmis[!is.na(bmis)]
  metadata$BMI[metadata$`Participant ID` == participant_id] <- ifelse(length(bmis) > 0, mean(bmis), NA)
  
  smoke_status <- metadata$`smoking status`[metadata$`Participant ID` == participant_id]
  smoke_status <- smoke_status[!is.na(smoke_status)]
  metadata$`smoking status`[metadata$`Participant ID` == participant_id] <- ifelse(length(smoke_status) > 0, smoke_status[1], NA)
}

if (opt$options$dataset == 'taxa') {
  metadata <- metadata[metadata$data_type == 'metagenomics',]
} else {
  metadata <- metadata[metadata$data_type == 'metabolomics',]
}
rownames(metadata) <- metadata$`External ID`
metadata <- metadata[,colSums(metadata == '', na.rm = T) != nrow(metadata)]
keep_cols <- c('External ID', 'Participant ID', 'week_num', 'site_name', 'Age at diagnosis',
               'Education Level', 'Occupation', 'consent_age', 'diagnosis',
               colnames(metadata)[c(52:83, 85:111)], 'race', 'sex', 'BMI')
metadata <- metadata[,keep_cols]
metadata <- metadata[,colSums(!is.na(metadata)) != 0]

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
metadata <- full_join(dysbiosis_df, metadata, by=c('sample' = 'External ID'))
metadata$dysbiosis_state <- metadata$dysbiosis_score > quantile(metadata$dysbiosis_score[metadata$diagnosis == 'nonIBD'], 0.9, na.rm=T)
metadata$dysbiosis_state <- ifelse(metadata$dysbiosis_state & metadata$diagnosis != 'nonIBD', paste0('dysbiosis_', metadata$diagnosis), 'none')
metadata$dysbiosis_state <- factor(metadata$dysbiosis_state, levels = c('none', 'dysbiosis_UC', 'dysbiosis_CD'))
rownames(metadata) <- metadata$sample

# Run Maaslin 3
Maaslin3_path <- paste0(gsub("/$", "", workingDirectory), "/Maaslin3/R/")
for (R_file in dir(Maaslin3_path, pattern = "*.R$")) {
  source(file.path(Maaslin3_path, R_file))
}

metadata$participant_id <- metadata$`Participant ID`
metadata$diagnosis <- factor(metadata$diagnosis, levels = c('nonIBD', 'UC', 'CD'))

if (opt$options$dataset == 'taxa') {
  tmp_fit_out <- paste0(gsub("/$", "", analysisDirectory), "/tmp_fit_out_Maaslin3ItAug")
} else {
  tmp_fit_out <- paste0(gsub("/$", "", analysisDirectory), "/tmp_fit_out_Maaslin3ItAug_mbx")
}

param_list <- list(input_data = taxa_table, input_metadata = metadata, min_abundance = 0, min_prevalence = 0, output = tmp_fit_out, 
                   min_variance = 0, normalization = 'TSS', transform = 'LOG', analysis_method = 'LM', 
                   formula = '~ diagnosis + dysbiosis_state + Antibiotics + consent_age + (1 | participant_id)', 
                   save_scatter = FALSE, 
                   save_models = F, plot_heatmap = F, plot_scatter = F, max_significance = 0.1, augment = T, iterative_mode = T, cores=nCores)
fit_out <- Maaslin3(param_list)

unlink(tmp_fit_out, recursive = T)

fit_out_lm <- fit_out$fit_data_non_zero$results
fit_out_lm <- fit_out_lm[c("feature", "metadata", "value", "coef", "pval_single", "error", "qval_single", "pval_joint", "qval_joint")]
fit_out_lm$association <- "abundance"

fit_out_binary <- fit_out$fit_data_binary$results
fit_out_binary <- fit_out_binary[c("feature", "metadata", "value", "coef", "pval_single", "error", "qval_single", "pval_joint", "qval_joint")]
fit_out_binary$association <- "prevalence"

fit_out_joint <- full_join(fit_out_lm, fit_out_binary, by = colnames(fit_out_lm))
fit_out_joint <- fit_out_joint[!is.na(fit_out_joint$feature),]

dir.create(paste0(gsub("/$", "", analysisDirectory), "/results"))
write.table(fit_out_joint, paste0(gsub("/$", "", analysisDirectory), "/results/", ifelse(opt$options$dataset == 'taxa', '', 'mbx_'), "ibd_associations_Maaslin3ItAug.tsv"), row.names = F, sep='\t')


