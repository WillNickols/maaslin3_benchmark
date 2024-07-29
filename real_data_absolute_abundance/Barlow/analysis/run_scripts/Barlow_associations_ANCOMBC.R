#!/usr/bin/env Rscript
package_vec = c("reshape2", "ggplot2", "optparse", 
                "parallel", "stringi", "doParallel", "plyr", "tidyr", 'dplyr', 'vegan',
                'ANCOMBC', 'TreeSummarizedExperiment')
invisible(suppressPackageStartupMessages(lapply(package_vec, require, character.only = TRUE)))

analysisDirectory <- "real_data_absolute_abundance/Barlow/analysis"

# Read in data
taxa_table <- read.csv(paste0(gsub('/$', '', analysisDirectory), '/data/Absolute_Abundance_Table.csv'), check.names = F, sep = ',')
rownames(taxa_table) <- taxa_table[,1]
taxa_table[,1] <- NULL

metadata <- taxa_table[,c('Diet', 'Site', 'Day', 'mouse', 'Cage')]
metadata <- metadata[metadata$Site == 'Stool',]
metadata$Site <- NULL
metadata$Cage <- NULL

read_depths <- read.csv(paste0(gsub('/$', '', analysisDirectory), '/data/filereport_read_run_PRJNA575097.tsv'), sep='\t')
metadata$sample_alias <- gsub(' ', '_', rownames(metadata))
metadata <- left_join(metadata, read_depths, by = c('sample_alias'))
rownames(metadata) <- gsub('_', ' ', metadata$sample_alias)
metadata$sample_alias <- NULL

taxa_table <- taxa_table[,!colnames(taxa_table) %in% c('Diet', 'Site', 'Day', 'mouse', 'Cage')]
taxa_table <- t(apply(taxa_table, 1, function(x) {x / sum(x)})) # Convert to relative abundance
colnames(taxa_table) <- make.names(colnames(taxa_table))
taxa_table <- t(taxa_table[rownames(taxa_table) %in% rownames(metadata),])

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
ancombc_out <- ancombc2(both_tse, fix_formula = 'Diet + Day + read_count', #rand_formula = '(1|mouse)', # Doesn't work for some reason
                        p_adj_method = "BH", 
                        alpha = 0.1, n_cl = 1, prv_cut = 0)

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
write.table(outputs, paste0(gsub("/$", "", analysisDirectory), "/results/", "Barlow_associations_ANCOMBC.tsv"), row.names = F, sep='\t')


