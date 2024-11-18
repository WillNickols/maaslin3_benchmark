#!/usr/bin/env Rscript
package_vec = c("reshape2", "ggplot2", "optparse", 
                "parallel", "stringi", "doParallel", "plyr", "tidyr", 'dplyr', 'vegan', 'maaslin3')
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

mtx_in <- read.csv("data/pathabundances_3_MTX.tsv", sep = '\t')
rownames(mtx_in) <- mtx_in$Feature.Sample
mtx_in$Feature.Sample <- NULL
mgx_in <- read.csv("data/pathabundances_3_MGX.tsv", sep = '\t')
rownames(mgx_in) <- mgx_in$Feature.Sample
mgx_in$Feature.Sample <- NULL
mgx_in <- mgx_in[,colnames(mtx_in)]
mgx_in <- mgx_in[rowSums(mgx_in != 0) != 0,]

colnames(mtx_in) <- gsub("_pathabundance_cpm", "", colnames(mtx_in))
colnames(mgx_in) <- gsub("_pathabundance_cpm", "", colnames(mgx_in))

mgx_in <- t(mgx_in)
mtx_in <- t(mtx_in)

taxa_table <- read.csv('data/metaphlan3_taxonomic_profiles.tsv',check.names = F, sep = '\t')

taxa_table <- taxa_table[taxa_table$`#SampleID` == 'UNCLASSIFIED' | 
                             grepl('\\|s__', taxa_table$`#SampleID`) & 
                             !grepl('\\|t__', taxa_table$`#SampleID`),]
rownames(taxa_table) <- taxa_table$`#SampleID`
taxa_table$`#SampleID` <- NULL
taxa_table['UNCLASSIFIED',] <- pmax(1 - colSums(taxa_table), 0)
taxa_table <- taxa_table * 100 # Convert to percents to be consistent with v4

prepare_metadata <- function(dataset_type) {
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

taxa_table <- read.csv('data/metaphlan3_taxonomic_profiles.tsv',check.names = F, sep = '\t')

taxa_table <- taxa_table[taxa_table$`#SampleID` == 'UNCLASSIFIED' | 
                             grepl('\\|s__', taxa_table$`#SampleID`) & 
                             !grepl('\\|t__', taxa_table$`#SampleID`),]
rownames(taxa_table) <- taxa_table$`#SampleID`
taxa_table$`#SampleID` <- NULL
taxa_table['UNCLASSIFIED',] <- pmax(1 - colSums(taxa_table), 0)
taxa_table <- taxa_table * 100 # Convert to percents to be consistent with v4

metadata <- metadata[metadata$sample %in% colnames(taxa_table),]
metadata <- metadata[as.numeric(mapvalues(metadata$sample, colnames(taxa_table), 1:ncol(taxa_table))),]

tmp_fit_out <- paste0(gsub("/$", "", analysisDirectory), "/fit_out_Maaslin3_", opt$options$dataset)

preprocess_out <- preprocess_dna_mtx(mgx_in, mtx_in)

set.seed(1)

if (opt$options$dataset == 'ratio') {
    mtx_mgx_ratios <- preprocess_out$rna_table/2^preprocess_out$dna_table
    
    # Set to NA since no normalization to take care of this
    mtx_mgx_ratios[mtx_mgx_ratios == 0] <- NA

    fit_out <- maaslin3::maaslin3(
        input_data = mtx_mgx_ratios,
        input_metadata = metadata,
        output = tmp_fit_out,
        normalization = 'NONE', 
        median_comparison_abundance = F,
        formula = 'diagnosis + dysbiosis_state + Antibiotics + consent_age + reads_filtered + (1 | participant_id)', 
        cores = nCores, 
        warn_prevalence = F)
} else {
    fit_out <- maaslin3::maaslin3(
        input_data = preprocess_out$rna_table,
        input_metadata = metadata,
        output = tmp_fit_out,
        normalization = 'NONE',
        median_comparison_abundance = F,
        formula = 'diagnosis + dysbiosis_state + Antibiotics + consent_age + reads_filtered + (1 | participant_id)',
        feature_specific_covariate = preprocess_out$dna_table,
        feature_specific_covariate_name = 'DNA',
        feature_specific_covariate_record = FALSE,
        cores = nCores)
}

fit_out_lm <- fit_out$fit_data_abundance$results
fit_out_lm <- fit_out_lm[c("feature", "metadata", "value", "coef", "pval_individual", "error", "qval_individual", "pval_joint", "qval_joint", "N", "N_not_zero")]
fit_out_lm$association <- "abundance"

fit_out_binary <- fit_out$fit_data_prevalence$results
fit_out_binary <- fit_out_binary[c("feature", "metadata", "value", "coef", "pval_individual", "error", "qval_individual", "pval_joint", "qval_joint", "N", "N_not_zero")]
fit_out_binary$association <- "prevalence"

fit_out_joint <- full_join(fit_out_lm, fit_out_binary, by = colnames(fit_out_lm))
fit_out_joint <- fit_out_joint[!is.na(fit_out_joint$feature),]

dir.create(paste0(gsub("/$", "", analysisDirectory), "/results"))
write.table(fit_out_joint, paste0(gsub("/$", "", analysisDirectory), "/results/mtx_", 
                                  ifelse(opt$options$dataset == 'ratio', 'ratio', 'covariate'), ".tsv"), 
            row.names = F, sep='\t')


