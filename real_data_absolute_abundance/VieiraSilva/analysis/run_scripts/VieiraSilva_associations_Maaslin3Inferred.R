#!/usr/bin/env Rscript
package_vec = c("reshape2", "ggplot2", "optparse", 
                "parallel", "stringi", "doParallel", "plyr", "tidyr", 'dplyr', 'vegan',
                'maaslin3')
invisible(suppressPackageStartupMessages(lapply(package_vec, require, character.only = TRUE)))

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

# Read taxa table again as absolute abundances
taxa_table <- read.csv(paste0(gsub('/$', '', analysisDirectory), '/data/QMP.matrix.tsv'), check.names = F, sep = '\t')
rownames(taxa_table) <- taxa_table[,1]
taxa_table[,1] <- NULL

tmp_fit_out <- paste0(gsub("/$", "", analysisDirectory), "/tmp_fit_out_Maaslin3")

fit_out <- maaslin3::maaslin3(input_data = taxa_table, 
                   input_metadata = metadata, 
                   output = tmp_fit_out, 
                   normalization = 'NONE', 
                   transform = 'LOG',
                   standardize = F,
                   formula = '~ diagnosis + age + gender + bmi + read_count',
                   plot_summary_plot = T, plot_associations = T, max_significance = 0.1, 
                   augment = T, median_comparison_abundance=F, median_comparison_prevalence=F)

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
write.table(fit_out_joint, paste0(gsub("/$", "", analysisDirectory), "/results/", "VieiraSilva_associations_Maaslin3Inferred.tsv"), row.names = F, sep='\t')


