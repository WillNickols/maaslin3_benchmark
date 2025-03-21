#!/usr/bin/env Rscript
package_vec = c("reshape2", "ggplot2", "optparse", 
                "parallel", "stringi", "doParallel", "plyr", "tidyr", 'dplyr', 'vegan',
                'ALDEx2')
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
unscaled <- rowSums(taxa_table)
taxa_table <- t(apply(taxa_table, 1, function(x) {x / sum(x)})) # Convert to relative abundance
colnames(taxa_table) <- make.names(colnames(taxa_table))
taxa_table <- t(taxa_table[rownames(taxa_table) %in% rownames(metadata),])
unscaled <- unscaled[names(unscaled) %in% rownames(metadata)]

for (col in colnames(metadata)) {
  if (is.numeric(metadata[,col])) {
    metadata[,col] <- scale(metadata[,col])
  }
}

# No read count because samples fixed at 19500 read depth
mm <- model.matrix(formula(paste0("~ Day + Diet + mouse")), metadata)

# This gets back the counts table since all samples are at 19500 read depth
taxa_table_in <- round(taxa_table / min(taxa_table[taxa_table > 0]))

taxa_table_in <- t(taxa_table_in)
taxa_table_in <- taxa_table_in[rownames(taxa_table_in) %in% rownames(mm),]
taxa_table_in <- t(taxa_table_in)

aldex_clr_out <- aldex.clr(taxa_table_in, mm, denom="all", useMC = F, gamma = matrix(unscaled, nrow = length(unscaled), ncol = 128))
glm.test <- aldex.glm(aldex_clr_out, verbose = T)

glm.test <- glm.test[,grepl("Est$|pval$", colnames(glm.test))]

glm.test <- glm.test[,!grepl(paste0(c("Intercept", unique(metadata$mouse)), collapse = "|"), colnames(glm.test))]
glm.test$Feature <- rownames(glm.test)
glm.test <- reshape2::melt(glm.test, id.vars = c("Feature"))
glm.test$metric <- gsub(".*\\:", "", glm.test$variable)
glm.test$variable <- gsub("\\:.*", "", glm.test$variable)
glm.test <- reshape2::dcast(formula = Feature + variable ~ metric, glm.test)
sink()

outputs <- data.frame(feature = glm.test$Feature,
                      metadata = glm.test$variable,
                      effect_size = glm.test$Est,
                      pval = glm.test$pval,
                      qval = p.adjust(glm.test$pval, method = "BH"),
                      associations = "abundance")

dir.create(paste0(gsub("/$", "", analysisDirectory), "/results"))
write.table(outputs, paste0(gsub("/$", "", analysisDirectory), "/results/", "Barlow_associations_ALDEx2_scale.tsv"), row.names = F, sep='\t')


