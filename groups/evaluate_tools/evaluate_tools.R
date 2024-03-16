#!/usr/bin/env Rscript
# Group

remove(list = ls())
package_vec = c("reshape2", "ggplot2", "optparse", 
                "parallel", "stringi", "doParallel", "plyr", "tidyr")
invisible(suppressPackageStartupMessages(lapply(package_vec, require, character.only = TRUE)))

# Command Line Usage
option_list = list(
  make_option(
    c("-l", "--parameters"),
    type = "character"))
opt <- parse_args(OptionParser(option_list=option_list), positional_arguments = TRUE)

# Extract Parameters
parameters<-normalizePath(opt$options$parameters)

filename_pieces <- unlist(strsplit(parameters, '/'))
workingDirectory <- paste0(filename_pieces[1:(length(filename_pieces) - 3)], collapse = '/')
generator <- gsub('.txt', '', filename_pieces[length(filename_pieces)])
source(paste0(workingDirectory, '/library/run_evaluation_helpers.R'))

nIterations = 100

parameter_text <- readLines(parameters)

param_list = list()
for (parameter_line in parameter_text) {
  parts = trimws(unlist(strsplit(parameter_line, ':')))
  vals = unlist(strsplit(parts[2], ' '))
  param_list[[parts[1]]] = vals
}
  
metadata_types = param_list[['metadataType']]
param_list[['metadataType']] <- NULL

params_for_files <- list()
for (param in names(param_list)) {
  first_params <- c()
  for (param_tmp in names(param_list)) {
    first_params <- c(first_params, param_list[[param_tmp]][1])
  }
  for (j in seq_along(param_list[[param]])) {
    params_for_files <- c(params_for_files, list(first_params))
    names(params_for_files[[length(params_for_files)]]) <- names(param_list)
    params_for_files[[length(params_for_files)]][param] <- param_list[[param]][j]
  }
}

params_for_files <- unique(params_for_files)

param_list_final = list()
for (param_file in params_for_files) {
  for (metadata_type in metadata_types) {
    tmp_param_file <- param_file
    tmp_param_file['metadataType'] <- metadata_type
    if (metadata_type == 'binary') {
      tmp_param_file['nMetadata'] <- '1'
      tmp_param_file['nPerSubject'] <- '1'
    }
    param_list_final <- c(param_list_final, list(tmp_param_file))
  }
}

general_results_df <- data.frame(matrix(nrow = 0, ncol = 0))
maaslin3_results_df <- data.frame(matrix(nrow = 0, ncol = 0))
metrics <- c("Precision", "Recall", "Precision\n(common taxa)", "Recall\n(common taxa)", 
             "AUC", "AUC (Common taxa)", "Issue proportion")
for (param_list in param_list_final) {
  print(param_list)
  if (param_list[['nPerSubject']]!='1'){
    inputSubString<-'RandomEffect'
  } else {
    inputSubString<-'noRandomEffect'
  }
  
  options("scipen"=10)
  inputString<-paste(inputSubString, 
                     param_list[['metadataType']], 
                     param_list[['nSubjects']], 
                     param_list[['nPerSubject']], 
                     param_list[['nMicrobes']], 
                     param_list[['spikeMicrobes']], 
                     param_list[['nMetadata']], 
                     param_list[['effectSize']], 
                     param_list[['effectPos']], 
                     param_list[['readDepth']], 
                     param_list[['nGroups']],
                     param_list[['nLevels']],
                     sep='_')
  options("scipen"=5)
  
  inputDirectory <- file.path(workingDirectory, 'Input', 'group_evaluations', generator)
  this_params_folder <- file.path(inputDirectory, inputString)
  outputDirectory <- file.path(workingDirectory, 'Output', 'group_evaluations', generator)
  tmp_files <- list.files(outputDirectory)
  tmp_files <- tmp_files[grepl(paste0(inputString, "_"), tmp_files)]
  tools <- gsub(paste0(inputString, "_"), '', tmp_files, perl = T)
  this_output_folder <- file.path(outputDirectory, inputString)
  
  for (i in 1:nIterations) {
    possible_error <- tryCatch({
      metadata <- read.csv(paste0(this_params_folder, "/metadata_", i, ".tsv"), sep = "\t")
      metadata <- metadata[,colnames(metadata) != 'ID']
      abundance <- read.csv(paste0(this_params_folder, "/abundance_", i, ".tsv"), sep = "\t")
      truth <- read.csv(paste0(this_params_folder, "/truth_", i, ".tsv"), sep = "\t")
    }, error = function(err) {
      err
    })
    if(inherits(possible_error, "error")) next
    
    truth <- prepare_truth_groups(truth, generator)
    for (tool in tools) {
      possible_error <- tryCatch({
        associations <- read.csv(paste0(this_output_folder, '_', tool, "/associations_", i, ".tsv"), sep = "\t")
      }, error = function(err) {
        err
      })
      
      if(inherits(possible_error, "error")) next
      
      associations_out <- prepare_associations_maaslin3_group(associations, tool)
      associations_out <- associations_out[!paste0(associations_out$taxon, '_', associations_out$metadata) %in% 
                                             paste0(truth$taxon, '_', truth$org_metadata)[truth$effect_size == 0],]
      associations_out$metadata <- gsub('[0-9]', '', associations_out$metadata)
      associations_out <- associations_out[associations_out$metadata %in% truth$metadata & !is.na(associations_out$signif),]
      
      new_row <- c(unweighted_precision_recall_maaslin3(truth, associations_out, abundance, metadata),
                   weighted_precision_recall_maaslin3(truth, associations_out, abundance, metadata),
                   pval_auc_maaslin3(truth, associations_out, abundance, metadata),
                   weighted_pval_auc_maaslin3(truth, associations_out, abundance, metadata),
                   issue_prop(associations, tool))
      new_row <- lapply(new_row, function(x) {x})
      names(new_row) <- metrics
      
      new_row[['tool']] <- tool
      new_row[['iter']] <- i
      new_row <- c(new_row, param_list)
      maaslin3_results_df <- plyr::rbind.fill(maaslin3_results_df, data.frame(new_row, check.names = F))
    }
  }
}

names_to_iter <- names(param_list)[names(param_list) != 'metadataType']
for (metadata_type in metadata_types) {
  figures_folder <- paste0(workingDirectory, '/Figures/group_evaluations/', generator, '/', metadata_type, '/')
  dir.create(figures_folder, recursive = T, showWarnings = F)
  for (param_name in names_to_iter) {
    default_vals <- first_params
    names(default_vals) <- names_to_iter
    default_vals <- default_vals[-which(names(default_vals) == param_name)]
    
    results_subset <- maaslin3_results_df[apply(maaslin3_results_df[,names(default_vals)], 1, function(x) {all(x == default_vals)}),]
    results_subset <- results_subset[results_subset$metadataType == metadata_type,]
    results_subset <- results_subset[,c(metrics, param_name, 'tool', 'iter')]
    
    melted_df <- melt(results_subset, id.vars = c('tool', 'iter', param_name))
    melted_df[,param_name] <- factor(melted_df[,param_name], levels = as.character(sort(unique(as.numeric(melted_df[,param_name])))))
    plot_out <- ggplot(melted_df, aes(x = get(param_name), y = value, fill = tool)) + 
      geom_boxplot(position = position_dodge(preserve = "single")) + 
      facet_wrap(~variable, scales = 'free') + 
      theme_bw() + 
      ylim(c(0,1)) +
      xlab(param_name) + 
      ylab('') + 
      theme(text=element_text(size=21))
    ggsave(paste0(figures_folder, 'Maaslin3_', param_name, '.png'),
           plot = plot_out, width = 18, height = 9)
  }
}






