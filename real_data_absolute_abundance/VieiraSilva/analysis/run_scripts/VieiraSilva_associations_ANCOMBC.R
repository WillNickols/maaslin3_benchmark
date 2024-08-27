#!/usr/bin/env Rscript
package_vec = c("reshape2", "ggplot2", "optparse", 
                "parallel", "stringi", "doParallel", "plyr", "tidyr", 'dplyr', 'vegan',
                'ANCOMBC', 'TreeSummarizedExperiment')
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

taxa_table <- taxa_table[rownames(metadata)[rownames(metadata) 
                                            %in% rownames(taxa_table)],]
metadata <- metadata[rownames(taxa_table),]

assays = S4Vectors::SimpleList(counts = t(as.matrix(taxa_table)))
smd = S4Vectors::DataFrame(metadata)
both_tse <- TreeSummarizedExperiment(
  assays = assays,
  colData = smd)
ancombc_out <- ancombc2(both_tse, fix_formula = 'diagnosis + age + gender + bmi + read_count',
                        p_adj_method = "BH", 
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
write.table(outputs, paste0(gsub("/$", "", analysisDirectory), "/results/", "VieiraSilva_associations_ANCOMBC.tsv"), row.names = F, sep='\t')


