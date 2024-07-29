#!/usr/bin/env Rscript
package_vec = c("reshape2", "ggplot2", "optparse", 
                "parallel", "stringi", "doParallel", "plyr", "tidyr", 'dplyr', 'vegan',
                'Maaslin2')
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

tmp_fit_out <- paste0(gsub("/$", "", analysisDirectory), "/tmp_fit_out_Maaslin2")
dir.create(tmp_fit_out, showWarnings = F)

fit_out <- Maaslin2::Maaslin2(taxa_table, metadata, min_abundance = 0, min_prevalence = 0, output = tmp_fit_out, 
                              min_variance = 0, normalization = 'TSS', transform = 'log', analysis_method = 'LM', 
                              random_effects = c('mouse'), fixed_effects = c("Diet", "Day"), # No read count because samples fixed at 19500 read depth
                              save_scatter = FALSE, save_models = F, plot_heatmap = F, plot_scatter = F, standardize = F,
                              max_significance = 0.1)$results

unlink(tmp_fit_out, recursive = T)

fit_out <- data.frame(feature = fit_out$feature,
           metadata = fit_out$metadata,
           value = fit_out$value,
           coef = fit_out$coef,
           pval = fit_out$pval,
           qval = fit_out$qval,
           association = "abundance")

dir.create(paste0(gsub("/$", "", analysisDirectory), "/results"))
write.table(fit_out, paste0(gsub("/$", "", analysisDirectory), "/results/", "Barlow_associations_Maaslin2.tsv"), row.names = F, sep='\t')


