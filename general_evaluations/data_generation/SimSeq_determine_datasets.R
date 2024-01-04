# Run from the folder containing the datasets of Nearing et al. 2023

folders <- list.files(".")
output_df <- data.frame(matrix(nrow = 0, ncol = 5))
for (folder in folders) {
  print(folder)
  files <- list.files(folder)
  genus_table_name <- files[grepl("_genus_table.tsv$", files)]
  taxa_table <- read.csv(paste0(folder, "/", genus_table_name), skip = 1, sep = "\t")
  print(dim(taxa_table))
  metadata_name <- files[grepl("_meta.tsv$|_metadata.tsv$", files)]
  if (length(metadata_name) != 0) {
    metadata <- read.csv(paste0(folder, "/", metadata_name), sep = "\t")
    print(dim(metadata))
    output_df <- rbind(output_df, c(folder, dim(taxa_table)[1], sort(table(metadata[,ncol(metadata)])), mean(taxa_table == 0)))
  }
}

colnames(output_df) <- c("Dataset", "Taxa", "Group 1", "Group 2", "Zeros")
output_df[,2:5] <- apply(output_df[,2:5], 2, as.numeric)

# Final datasets
output_df[output_df$`Group 1` > 50 & output_df$`Group 2` > 100 & output_df$Zeros < 0.9,]$Dataset
