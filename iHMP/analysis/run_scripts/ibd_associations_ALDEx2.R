#!/usr/bin/env Rscript
package_vec = c("reshape2", "ggplot2", "optparse", 
                "parallel", "stringi", "doParallel", "plyr", "tidyr", 'dplyr', 'vegan',
                'ALDEx2')
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
    mbx_table <- t(TSSnorm(t(mbx_table)))
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

metadata$participant_id <- metadata$`Participant ID`
metadata$diagnosis <- factor(metadata$diagnosis, levels = c('nonIBD', 'UC', 'CD'))

metadata <- metadata[metadata$sample %in% colnames(taxa_table),]
mm <- model.matrix(formula(paste0("~ diagnosis + dysbiosis_state + Antibiotics + consent_age")), metadata)

taxa_table_in <- round(taxa_table / min(taxa_table[taxa_table > 0]))

taxa_table_in <- t(taxa_table_in)
taxa_table_in <- taxa_table_in[rownames(taxa_table_in) %in% rownames(mm),]
taxa_table_in <- t(taxa_table_in)

aldex_clr_out <- aldex.clr(taxa_table_in, mm, denom="all", useMC = F)
glm.test <- aldex.glm(aldex_clr_out)

glm.test <- glm.test[,grepl("Est$|pval$", colnames(glm.test))]

if ('ID' %in% colnames(metadata)) {
  glm.test <- glm.test[,!grepl(paste0(c("Intercept", unique(metadata$ID)), collapse = "|"), colnames(glm.test))]
} else {
  glm.test <- glm.test[,!grepl("Intercept", colnames(glm.test))]
}
glm.test$Feature <- rownames(glm.test)
glm.test <- reshape2::melt(glm.test, id.vars = c("Feature"))
glm.test$metric <- gsub(".*\\:", "", glm.test$variable)
glm.test$variable <- gsub("\\:.*", "", glm.test$variable)
glm.test <- reshape2::dcast(formula = Feature + variable ~ metric, glm.test)
sink()

outputs <- data.frame(feature = glm.test$Feature,
                      metadata = glm.test$variable,
                      effect_size = glm.test$Est,
                      pval = glm.test$pval,
                      qval = p.adjust(glm.test$pval, method = "BH"),
                      associations = "abundance")

dir.create(paste0(gsub("/$", "", analysisDirectory), "/results"))
write.table(outputs, paste0(gsub("/$", "", analysisDirectory), "/results/", ifelse(opt$options$dataset == 'taxa', '', 'mbx_'), "ibd_associations_ALDEx2.tsv"), row.names = F, sep='\t')


