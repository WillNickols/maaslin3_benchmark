#!/usr/bin/env Rscript
package_vec = c("reshape2", "ggplot2", "optparse", 
                "parallel", "stringi", "doParallel", "plyr", "tidyr", 'dplyr', 'vegan',
                'maaslin3')
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

for (col in colnames(metadata)) {
  if (is.numeric(metadata[,col])) {
    metadata[,col] <- scale(metadata[,col])
  }
}

tmp_fit_out <- paste0(gsub("/$", "", analysisDirectory), "/tmp_fit_out_Maaslin3")

param_list <- list(input_data = taxa_table, 
                   input_metadata = metadata, 
                   output = tmp_fit_out, 
                   normalization = 'TSS', 
                   transform = 'LOG',
                   standardize = F,
                   formula = '~ Diet + Day + (1|mouse)', # No read count because samples fixed at 19500 read depth
                   plot_summary_plot = T, plot_associations = T, max_significance = 0.1, 
                   augment = T, median_comparison_abundance=T, median_comparison_prevalence=F)
fit_out <- maaslin3::maaslin3(param_list)

unlink(tmp_fit_out, recursive = T)

fit_out_lm <- fit_out$fit_data_abundance$results
fit_out_lm <- fit_out_lm[c("feature", "metadata", "value", "coef", "pval_individual", "error", "qval_individual", "pval_joint", "qval_joint")]
fit_out_lm$association <- "abundance"

fit_out_binary <- fit_out$fit_data_prevalence$results
fit_out_binary <- fit_out_binary[c("feature", "metadata", "value", "coef", "pval_individual", "error", "qval_individual", "pval_joint", "qval_joint")]
fit_out_binary$association <- "prevalence"

fit_out_joint <- full_join(fit_out_lm, fit_out_binary, by = colnames(fit_out_lm))
fit_out_joint <- fit_out_joint[!is.na(fit_out_joint$feature),]

dir.create(paste0(gsub("/$", "", analysisDirectory), "/results"))
write.table(fit_out_joint, paste0(gsub("/$", "", analysisDirectory), "/results/", "Barlow_associations_Maaslin3CompAdjust.tsv"), row.names = F, sep='\t')
