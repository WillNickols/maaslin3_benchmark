#!/usr/bin/env Rscript
package_vec = c("reshape2", "ggplot2", "optparse", 
                "parallel", "stringi", "doParallel", "plyr", "tidyr", 'dplyr', 'vegan',
                'ALDEx2')
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

mm <- model.matrix(formula(paste0("~ day + infant + read_depth")), metadata)

# Convert to relative abundances then counts since ALDEx2 needs that
taxa_table_in <- apply(t(taxa_table), 1, function(x) {x / sum(x)})
scaling_factors <- read_depth_file$read_count[match(colnames(taxa_table_in), read_depth_file$sample_name)]
taxa_table_in <- round(sweep(taxa_table_in, 2, scaling_factors, `*`) + 0.5)

taxa_table_in <- t(taxa_table_in)
taxa_table_in <- taxa_table_in[rownames(taxa_table_in) %in% rownames(mm),]
taxa_table_in <- t(taxa_table_in)

aldex_clr_out <- aldex.clr(taxa_table_in, mm, denom="all", useMC = F, gamma = 0.5)
glm.test <- aldex.glm(aldex_clr_out, verbose = T)

glm.test <- glm.test[,grepl("Est$|pval$", colnames(glm.test))]

glm.test <- glm.test[,!grepl(paste0(c("Intercept", unique(metadata$infant)), collapse = "|"), colnames(glm.test))]
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
write.table(outputs, paste0(gsub("/$", "", analysisDirectory), "/results/", "infant_associations_ALDEx2.tsv"), row.names = F, sep='\t')


