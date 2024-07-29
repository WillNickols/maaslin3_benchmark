#!/usr/bin/env Rscript
package_vec = c("reshape2", "ggplot2", "optparse", 
                "parallel", "stringi", "doParallel", "plyr", "tidyr", 'dplyr', 'vegan',
                'Maaslin2')
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

tmp_fit_out <- paste0(gsub("/$", "", analysisDirectory), "/tmp_fit_out_Maaslin2")

fit_out <- Maaslin2::Maaslin2(taxa_table, metadata, min_abundance = 0, min_prevalence = 0, output = tmp_fit_out, 
                              min_variance = 0, normalization = 'TSS', transform = 'log', analysis_method = 'LM', 
                              random_effects = c("infant"), fixed_effects = c('day', 'read_depth'), 
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
write.table(fit_out, paste0(gsub("/$", "", analysisDirectory), "/results/", "infant_associations_Maaslin2.tsv"), row.names = F, sep='\t')


