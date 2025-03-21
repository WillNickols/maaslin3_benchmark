#!/usr/bin/env Rscript
package_vec = c("reshape2", "ggplot2", "optparse", 
                "parallel", "stringi", "doParallel", "plyr", "tidyr", 'dplyr', 'vegan',
                'ALDEx2')
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

mm <- model.matrix(formula(paste0("~ diagnosis + age + gender + bmi + read_count")), metadata)

# Convert to relative abundances then counts since ALDEx2 needs that
taxa_table_in <- round(taxa_table * read_counts)
taxa_table_in <- taxa_table_in[rownames(taxa_table_in) %in% rownames(mm),]
taxa_table_in <- taxa_table_in[rownames(mm),]
mm <- mm[rownames(taxa_table_in),]
taxa_table_in <- t(taxa_table_in)

aldex_clr_out <- aldex.clr(taxa_table_in, mm, denom="all", useMC = F, gamma = 0.5)
glm.test <- aldex.glm(aldex_clr_out, verbose = T)

glm.test <- glm.test[,grepl("Est$|pval$", colnames(glm.test))]

glm.test <- glm.test[,!grepl(paste0(c("Intercept"), collapse = "|"), colnames(glm.test))]
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
write.table(outputs, paste0(gsub("/$", "", analysisDirectory), "/results/", "VieiraSilva_associations_ALDEx2.tsv"), row.names = F, sep='\t')


