#!/usr/bin/env Rscript
package_vec = c("reshape2", "ggplot2", "optparse", 
                "parallel", "stringi", "doParallel", "plyr", "tidyr", 'dplyr', 'vegan',
                'maaslin3')
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

for (col in colnames(metadata)) {
  if (is.numeric(metadata[,col])) {
    metadata[,col] <- scale(metadata[,col])
  }
}

tmp_fit_out <- paste0(gsub("/$", "", analysisDirectory), "/tmp_fit_out_Maaslin3")

param_list <- list(input_data = taxa_table, 
                   input_metadata = metadata, 
                   output = tmp_fit_out, 
                   normalization = 'NONE', 
                   transform = 'LOG',
                   standardize = F,
                   formula = '~ day + read_depth + (1|infant)',
                   plot_summary_plot = T, plot_associations = T, max_significance = 0.1, 
                   augment = T, median_comparison_abundance=F, median_comparison_prevalence=F)
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
write.table(fit_out_joint, paste0(gsub("/$", "", analysisDirectory), "/results/", "infant_associations_Maaslin3Inferred.tsv"), row.names = F, sep='\t')


