Maaslin3_path <- "/Users/williamnickols/Documents/GitHub/Maaslin3/R/"

for (R_file in dir(Maaslin3_path, pattern = "*.R$")) {
  source(file.path(Maaslin3_path, R_file))
}

taxa_table <- read.csv('hmp_aerobiosis_small.txt', sep = '\t', header = F)
metadata <- t(taxa_table[1:3,])
colnames(metadata) <- metadata[1,]
metadata <- metadata[-1,]
rownames(metadata) <- paste0("Sample", 1:nrow(metadata))
metadata <- data.frame(metadata)
metadata$oxygen_availability <- factor(metadata$oxygen_availability, levels = c('Low_O2', 'Mid_O2', 'High_O2'))

taxa_table <- taxa_table[-c(1:3),]
rownames_taxa_table <- taxa_table[,1]
taxa_table <- apply(taxa_table[,-1], 2, function(x){x <- as.numeric(x); x <- ifelse(x == min(x), 0, x)})
rownames(taxa_table) <- rownames_taxa_table
colnames(taxa_table) <- paste0("Sample", 1:nrow(metadata))
taxa_table <- data.frame(taxa_table)

param_list <- list(input_data = taxa_table, input_metadata = metadata, min_abundance = 0, min_prevalence = 0, output = 'output/', 
                   min_variance = 0, normalization = 'NONE', transform = 'LOG', analysis_method = 'LM', 
                   formula = '~ oxygen_availability + (1 | subject_id)', 
                   save_scatter = FALSE, save_models = FALSE, plot_heatmap = F, plot_scatter = F, 
                   max_significance = 0.1, augment = TRUE, iterative_mode = F, cores=6)
fit_out <- Maaslin3(param_list)
maaslin_write_results_lefse_format(fit_out)







