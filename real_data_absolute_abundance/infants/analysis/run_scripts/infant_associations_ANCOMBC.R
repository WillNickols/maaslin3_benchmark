#!/usr/bin/env Rscript
package_vec = c("reshape2", "ggplot2", "optparse", 
                "parallel", "stringi", "doParallel", "plyr", "tidyr", 'dplyr', 'vegan',
                'ANCOMBC', 'TreeSummarizedExperiment')
invisible(suppressPackageStartupMessages(lapply(package_vec, require, character.only = TRUE)))

analysisDirectory <- "real_data_absolute_abundance/infants/analysis"

# Read in data
taxa_table <- read.csv(paste0(gsub('/$', '', analysisDirectory), '/data/41586_2021_3241_MOESM4_ESM.csv'), skip = 1, check.names = F, sep = ',')
# Drop OTU_ID and confidence
taxa_table <- taxa_table[,grepl("NICU|OTU_ID", colnames(taxa_table))]
taxa_table <- taxa_table[rowSums(is.na(taxa_table)) == 0,]

rownames(taxa_table) <- taxa_table$OTU_ID
taxa_table$OTU_ID <- NULL
metadata <- data.frame(sample = colnames(taxa_table), 
                       infant = as.character(gsub("^[^_]*_(.*?)_.+", "\\1", colnames(taxa_table))),
                       day = as.numeric(gsub(".*_([0-9]+)(?:re.*)?$", "\\1", colnames(taxa_table))))

read_depth_file <- read.csv(paste0(gsub('/$', '', analysisDirectory), '/data/filereport_read_run_PRJEB36435.tsv'), sep='\t')
read_depth_file$sample_name <- gsub('_R2.*', '', gsub('.*/', '', read_depth_file$submitted_ftp))
read_depth_file <- read_depth_file[grepl('NICU_', read_depth_file$sample_name) & 
                                     grepl('_bac16S$', read_depth_file$sample_name),]
read_depth_file$sample_name <- gsub('_bac16.*', '', read_depth_file$sample_name)
read_depth_file <- read_depth_file[,c('sample_name', 'read_count')]

metadata <- left_join(metadata, read_depth_file, by=c('sample' = 'sample_name'))
rownames(metadata) <- metadata$sample
metadata$sample <- NULL
colnames(metadata) <- c('infant', 'day', 'read_depth')

taxa_table <- apply(t(taxa_table), 1, function(x) {x / sum(x)}) # Convert to relative abundance

for (col in colnames(metadata)) {
  if (is.numeric(metadata[,col])) {
    metadata[,col] <- scale(metadata[,col])
  }
}

assays = S4Vectors::SimpleList(counts = as.matrix(taxa_table))
smd = S4Vectors::DataFrame(metadata)
both_tse <- TreeSummarizedExperiment(
  assays = assays,
  colData = smd)
ancombc_out <- ancombc2(both_tse, fix_formula = 'day + read_depth',
                        rand_formula = '(1|infant)', p_adj_method = "BH", 
                        alpha = 0.1, n_cl = 1,
                        prv_cut = 0)

glm.test <- ancombc_out$res
glm.test <- glm.test[,grepl("^taxon|^lfc_|^p_|^passed_ss_", colnames(glm.test))]
glm.test <- glm.test[,!grepl("Intercept", colnames(glm.test))]

glm.test <- reshape2::melt(glm.test, id.vars = c("taxon"))
glm.test$metric <- gsub("_.*", "", glm.test$variable)
glm.test$variable <- gsub("^[^_]*_|^passed_ss_", "", glm.test$variable)
glm.test <- reshape2::dcast(formula = taxon + variable ~ metric, glm.test)

outputs <- data.frame(feature = glm.test$taxon,
                      metadata = glm.test$variable,
                      effect_size = glm.test$lfc,
                      pval = glm.test$p,
                      qval = p.adjust(glm.test$p, method = "BH"),
                      error = ifelse(glm.test$passed == 1, NA, "sensitivity failed"),
                      associations = ifelse(is.infinite(glm.test$lfc), "prevalence", "abundance"))

dir.create(paste0(gsub("/$", "", analysisDirectory), "/results"))
write.table(outputs, paste0(gsub("/$", "", analysisDirectory), "/results/", "infant_associations_ANCOMBC.tsv"), row.names = F, sep='\t')


