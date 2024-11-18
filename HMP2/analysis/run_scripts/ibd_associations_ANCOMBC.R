#!/usr/bin/env Rscript
package_vec = c("reshape2", "ggplot2", "optparse", 
                "parallel", "stringi", "doParallel", "plyr", "tidyr", 'dplyr', 'vegan',
                'ANCOMBC', 'TreeSummarizedExperiment')
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
      type = "character"),
  make_option(
      c("-v", "--version"), # d stands for dataset
      type = "integer"))
opt <- parse_args(OptionParser(option_list=option_list), positional_arguments = TRUE)
workingDirectory <- opt$options$workingDirectory
analysisDirectory <- opt$options$analysisDirectory
nCores <- opt$options$nCores

if (opt$options$version == 3) {
    taxa_table <- read.csv('data/metaphlan3_taxonomic_profiles.tsv',check.names = F, sep = '\t')
    
    taxa_table <- taxa_table[taxa_table$`#SampleID` == 'UNCLASSIFIED' | 
                                 grepl('\\|s__', taxa_table$`#SampleID`) & 
                                 !grepl('\\|t__', taxa_table$`#SampleID`),]
    rownames(taxa_table) <- taxa_table$`#SampleID`
    taxa_table$`#SampleID` <- NULL
    taxa_table['UNCLASSIFIED',] <- pmax(1 - colSums(taxa_table), 0)
    taxa_table <- taxa_table * 100 # Convert to percents to be consistent with v4
} else if (opt$options$version == 4) {
    taxa_table <- read.csv('data/metaphlan4_taxonomic_profiles.tsv', skip = 1, check.names = F, sep = '\t')
    
    # Reorganize taxa table
    colnames(taxa_table) <- gsub('_taxonomic$', '', colnames(taxa_table))
    taxa_table <- taxa_table[taxa_table$clade_name == 'UNCLASSIFIED' | 
                                 grepl('\\|t__', taxa_table$clade_name),]
    rownames(taxa_table) <- taxa_table$clade_name
    taxa_table$clade_name <- NULL
} else {
    stop("--version not valid")
}

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

metadata <- prepare_metadata(opt$options$dataset)
metadata <- right_join(dysbiosis_df, metadata, by=c('sample'='External ID'))

metadata$participant_id <- metadata$`Participant ID`
metadata$diagnosis <- factor(metadata$diagnosis, levels = c('nonIBD', 'UC', 'CD'))

# Make sure only paired MGX/MBX are used since we need the dysbiosis score
metadata <- metadata[!is.na(metadata$dysbiosis_state),]
rownames(metadata) <- metadata$sample

if (opt$options$dataset == 'taxa') {
    if (opt$options$version == 3) {
        taxa_table <- read.csv('data/metaphlan3_taxonomic_profiles.tsv',check.names = F, sep = '\t')
        
        taxa_table <- taxa_table[taxa_table$`#SampleID` == 'UNCLASSIFIED' | 
                                     grepl('\\|s__', taxa_table$`#SampleID`) & 
                                     !grepl('\\|t__', taxa_table$`#SampleID`),]
        rownames(taxa_table) <- taxa_table$`#SampleID`
        taxa_table$`#SampleID` <- NULL
        taxa_table['UNCLASSIFIED',] <- pmax(1 - colSums(taxa_table), 0)
        taxa_table <- taxa_table * 100 # Convert to percents to be consistent with v4
    } else if (opt$options$version == 4) {
        taxa_table <- read.csv('data/metaphlan4_taxonomic_profiles.tsv', skip = 1, check.names = F, sep = '\t')
        
        # Reorganize taxa table
        colnames(taxa_table) <- gsub('_taxonomic$', '', colnames(taxa_table))
        taxa_table <- taxa_table[taxa_table$clade_name == 'UNCLASSIFIED' | 
                                     grepl('\\|t__', taxa_table$clade_name),]
        rownames(taxa_table) <- taxa_table$clade_name
        taxa_table$clade_name <- NULL
    } else {
        stop("--version not valid")
    }
    
} else {
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

metadata <- metadata[metadata$sample %in% colnames(taxa_table),]
metadata <- metadata[as.numeric(mapvalues(metadata$sample, colnames(taxa_table), 1:ncol(taxa_table))),]

# Convert to pseudo-read counts
taxa_table_in <- round(taxa_table / min(taxa_table[taxa_table > 0]))

taxa_table_in <- t(taxa_table_in)
taxa_table_in <- taxa_table_in[rownames(taxa_table_in) %in% rownames(metadata),]
taxa_table_in <- taxa_table_in[match(rownames(metadata), rownames(taxa_table_in)),]
taxa_table <- t(taxa_table_in)

assays = S4Vectors::SimpleList(counts = as.matrix(taxa_table))
smd = S4Vectors::DataFrame(metadata)
both_tse <- TreeSummarizedExperiment(
  assays = assays,
  colData = smd)
if (opt$options$dataset == 'taxa') {
  ancombc_out <- ancombc2(both_tse, 
                          fix_formula = 'diagnosis + dysbiosis_state + Antibiotics + consent_age + reads_filtered',
                          rand_formula = '(1|participant_id)', 
                          p_adj_method = "BH", 
                          alpha = 0.1, 
                          n_cl = nCores, prv_cut = 0)
} else {
  ancombc_out <- ancombc2(both_tse, 
                          fix_formula = 'diagnosis + dysbiosis_state + Antibiotics + consent_age',
                          rand_formula = '(1|participant_id)', 
                          p_adj_method = "BH", 
                          alpha = 0.1, 
                          n_cl = nCores, prv_cut = 0)
}

glm.test <- ancombc_out$res
glm.test <- glm.test[,grepl("^taxon|^lfc_|^p_|^passed_ss_", colnames(glm.test))]
glm.test <- glm.test[,!grepl("Intercept", colnames(glm.test))]

glm.test <- reshape2::melt(glm.test, id.vars = c("taxon"))
glm.test$metric <- gsub("_.*", "", glm.test$variable)
glm.test$variable <- gsub("^[^_]*_|^passed_ss_", "", glm.test$variable)
glm.test <- reshape2::dcast(formula = taxon + variable ~ metric, glm.test)

# struc_zeros <- ancombc_out$zero_ind[ancombc_out$zero_ind[,2] != ancombc_out$zero_ind[,3],]
# glm.test <- rbind(glm.test, data.frame(taxon = struc_zeros$taxon,
#                                        variable = 'dysbiosis_state',
#                                        lfc = ifelse(struc_zeros[,3] == F, Inf, -Inf),
#                                        p = 0,
#                                        passed = 1))

outputs <- data.frame(feature = glm.test$taxon,
                      metadata = glm.test$variable,
                      effect_size = glm.test$lfc,
                      pval = glm.test$p,
                      qval = p.adjust(glm.test$p, method = "BH"),
                      error = ifelse(glm.test$passed == 1, NA, "sensitivity failed"),
                      associations = ifelse(is.infinite(glm.test$lfc), "prevalence", "abundance"))

dir.create(paste0(gsub("/$", "", analysisDirectory), "/results"))
write.table(outputs, paste0(gsub("/$", "", analysisDirectory), "/results/", ifelse(opt$options$dataset == 'taxa', '', 'mbx_'), "v", opt$options$version, "_ibd_associations_ANCOMBC.tsv"), row.names = F, sep='\t')


