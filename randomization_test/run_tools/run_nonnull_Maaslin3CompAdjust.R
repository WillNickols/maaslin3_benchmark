remove(list = ls())
package_vec = c("devtools", "pkgmaker", "optparse", 
                "parallel", "stringi", "doParallel", "plyr", "tidyr", "maaslin3", "dplyr")
invisible(suppressPackageStartupMessages(lapply(package_vec, require, character.only = TRUE)))

option_list = list(
  make_option(
    c("-w", "--workingDirectory"), # w stands for working directory
    type = "character"),
  make_option(
    c("-d", "--dataset"),
    type = "character"))
opt <- parse_args(OptionParser(option_list=option_list), positional_arguments = TRUE)

setwd(opt$options$workingDirectory)

dataset <- opt$options$dataset

# Don't randomize
if (!grepl("OTU|asv|otu", readLines(paste0('randomization_test/raw_data/', dataset, '_ASVs_table.tsv'), n = 1))) {
  abun_in <- read.csv(paste0('randomization_test/raw_data/', dataset, '_ASVs_table.tsv'), sep='\t', skip = 1, row.names = 1)
} else {
  abun_in <- read.csv(paste0('randomization_test/raw_data/', dataset, '_ASVs_table.tsv'), sep='\t', row.names = 1)
}

metadata_in <- read.csv(paste0('randomization_test/raw_data/', dataset, '_metadata.tsv'), sep='\t')
if (ncol(metadata_in) == 2) {
  colnames(metadata_in) <- c("sample", "metadatum")
  rownames(metadata_in) <- make.names(metadata_in$sample)
  metadata_in$sample <- NULL
} else {
  rownames(metadata_in) <- make.names(rownames(metadata_in))
  colnames(metadata_in) <- "metadatum"
}

intersected_samples <- intersect(colnames(abun_in), rownames(metadata_in))

abun_in <- abun_in[, intersected_samples, drop = F]
metadata_in <- metadata_in[intersected_samples, , drop = F]

if (dim(abun_in)[2] != dim(metadata_in)[1]) {
  stop("Columns and rows not equal")
}

if (dim(abun_in)[2] == 0) {
  stop("All samples gone")
}

if (dim(metadata_in)[2] != 1) {
  stop("More than 1 metadatum")
}

dir.create(paste0('randomization_test/data_non_null/'), showWarnings = F, recursive = T)
write.table(abun_in, paste0('randomization_test/data_non_null/', dataset, '_ASVs_table.tsv'), sep='\t', row.names = T)
write.table(metadata_in, paste0('randomization_test/data_non_null/', dataset, '_metadata.tsv'), sep='\t', row.names = T)

# Run MaAsLin 3 on all datasets
growing_outputs <- data.frame()
abundance <- read.csv(paste0('randomization_test/data_non_null/', dataset, '_ASVs_table.tsv'), sep='\t')
metadata <- read.csv(paste0('randomization_test/data_non_null/', dataset, '_metadata.tsv'), sep='\t')

# Get read depths
abundance <- abundance[apply(abundance, 1, var) != 0,]
read_depths <- data.frame(sample = names(colSums(abundance)),
                          read_depth = colSums(abundance))
metadata$sample <- rownames(metadata)
metadata <- left_join(metadata, read_depths, by = 'sample')
rownames(metadata) <- metadata$sample
metadata$sample <- NULL
metadata$metadatum <- factor(metadata$metadatum)

tmp_fit_out <- paste0('randomization_test', "/tmp_out_non_null/", dataset)
dir.create(tmp_fit_out, recursive = T, showWarnings = F)

set.seed(1)
fit_out <- maaslin3::maaslin3(input_data = abundance, 
                   input_metadata = metadata, 
                   output = tmp_fit_out, 
                   normalization = 'TSS', 
                   transform = 'LOG',
                   formula = 'metadatum', 
                   median_comparison_abundance = T, 
                   median_comparison_prevalence = F,
                   plot_summary_plot = F, 
                   plot_associations = F, 
                   max_significance = 0.1)

unlink(tmp_fit_out, recursive = T)

fit_out_lm <- fit_out$fit_data_abundance$results
fit_out_lm <- fit_out_lm[c("feature", "metadata", "coef", "pval_individual", "error", "qval_individual", "pval_joint", "qval_joint")]
fit_out_lm$association <- "abundance"

fit_out_binary <- fit_out$fit_data_prevalence$results
fit_out_binary <- fit_out_binary[c("feature", "metadata", "coef", "pval_individual", "error", "qval_individual", "pval_joint", "qval_joint")]
fit_out_binary$association <- "prevalence"

fit_out <- full_join(fit_out_lm, fit_out_binary, by = colnames(fit_out_lm))

outputs <- data.frame(pval = fit_out$pval_individual,
                      qval = fit_out$qval_individual,
                      pval_joint = fit_out$pval_joint,
                      qval_joint = fit_out$qval_joint,
                      error = fit_out$error,
                      associations = fit_out$association)
outputs <- outputs[is.na(outputs$error),]
outputs$error <- NULL
outputs$dataset <- dataset

growing_outputs <- rbind(growing_outputs, outputs)

dir.create('randomization_test/associations/', showWarnings = F, recursive = T)
write.table(growing_outputs,
            file = paste0('randomization_test/associations/', dataset, "_non_null.tsv"),
            sep = '\t',
            row.names = F)

