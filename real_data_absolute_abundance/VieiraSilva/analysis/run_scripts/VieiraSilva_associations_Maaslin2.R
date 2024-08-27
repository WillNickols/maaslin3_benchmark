#!/usr/bin/env Rscript
package_vec = c("reshape2", "ggplot2", "optparse", 
                "parallel", "stringi", "doParallel", "plyr", "tidyr", 'dplyr', 'vegan',
                'Maaslin2')
analysisDirectory <- "real_data_absolute_abundance/VieiraSilva/analysis"

taxa_table <- read.csv(paste0(gsub('/$', '', analysisDirectory), '/data/QMP.matrix.tsv'), check.names = F, sep = '\t')
rownames(taxa_table) <- taxa_table[,1]
taxa_table[,1] <- NULL
taxa_table <- t(apply(taxa_table, 1, function(x) {x / sum(x)})) # Convert to relative abundance

metadata <- read.csv(paste0(gsub('/$', '', analysisDirectory), '/data/41564_2019_483_MOESM3_ESM.csv'), check.names = F, sep = ',')
metadata <- metadata[,c(1:5, 7, 8)]
colnames(metadata) <- c("sample", "diagnosis", "age", "gender", "bmi", "fc", "crp")
metadata <- metadata[rowSums(is.na(metadata)) == 0,] # Drops 20
rownames(metadata) <- metadata$sample
metadata$sample <- NULL
metadata <- metadata[metadata$diagnosis != "UC",] # Only 4, too few
metadata$diagnosis[metadata$diagnosis == "PSC"] <- 'PSC-only'
metadata$diagnosis <- factor(metadata$diagnosis, levels = c('mHC', 'PSC', 'PSC-UC', 'CD', 'PSC-CD'))
metadata$gender <- factor(metadata$gender, levels = c("F", "M"))

read_counts <- apply(taxa_table, 1, function(x) {1 / min(x[x>0])})

# read counts times relative abundances give all integers, so these were the original read counts
if (sum(abs(taxa_table * read_counts - round(taxa_table * read_counts)) > 0.01) > 0) {
    stop("Read counts went wrong")
}

read_count_df <- data.frame(sample = names(read_counts), read_count = read_counts)
metadata$sample <- rownames(metadata)
metadata <- full_join(metadata, read_count_df, by = c("sample"))
rownames(metadata) <- metadata$sample
metadata$sample <- NULL

for (col in colnames(metadata)) {
    if (is.numeric(metadata[,col])) {
        metadata[,col] <- scale(metadata[,col])
    }
}

tmp_fit_out <- paste0(gsub("/$", "", analysisDirectory), "/tmp_fit_out_Maaslin2")

fit_out <- Maaslin2::Maaslin2(taxa_table, metadata, min_abundance = 0, min_prevalence = 0, output = tmp_fit_out, 
                              min_variance = 0, normalization = 'TSS', transform = 'log', analysis_method = 'LM', 
                              fixed_effects = c('diagnosis', 'age', 'gender', 'bmi', 'read_count'), 
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
write.table(fit_out, paste0(gsub("/$", "", analysisDirectory), "/results/", "VieiraSilva_associations_Maaslin2.tsv"), row.names = F, sep='\t')


