remove(list = ls())
package_vec = c("reshape2", "ggplot2", "optparse", 
                "parallel", "stringi", "doParallel", "plyr", "tidyr")
invisible(suppressPackageStartupMessages(lapply(package_vec, require, character.only = TRUE)))

parameters <- "~/Documents/GitHub/maaslin3_benchmark/unscaled/data_generation/SD2.txt"

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
for (param_list in param_list_final) { #c(2, 14, 17, 20, 23, 26)
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
                     sep='_')
  options("scipen"=5)
  
  inputDirectory <- file.path(workingDirectory, 'Input', 'unscaled', generator)
  this_params_folder <- file.path(inputDirectory, inputString)
  outputDirectory <- file.path(workingDirectory, 'Output', 'unscaled', generator)
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
    
    truth <- prepare_truth(truth, generator)
    for (tool in tools) {
      possible_error <- tryCatch({
        associations <- read.csv(paste0(this_output_folder, '_', tool, "/associations_", i, ".tsv"), sep = "\t")
      }, error = function(err) {
        err
      })
      
      if(inherits(possible_error, "error")) next
      
      new_row <- c(unweighted_precision_recall(truth, prepare_associations_general(associations, tool, generator), abundance, metadata),
                   weighted_precision_recall(truth, prepare_associations_general(associations, tool, generator), abundance, metadata),
                   effect_size_error(truth, prepare_associations_abundance(associations, tool, generator)),
                   effect_size_correlation(truth, prepare_associations_abundance(associations, tool, generator)),
                   pval_auc(truth, prepare_associations_general(associations, tool, generator), abundance, metadata),
                   weighted_pval_auc(truth, prepare_associations_general(associations, tool, generator), abundance, metadata),
                   issue_prop(associations, tool))
      new_row <- lapply(new_row, function(x) {x})
      names(new_row) <- metrics
      
      new_row[['tool']] <- tool
      new_row[['iter']] <- i
      new_row <- c(new_row, param_list)
      general_results_df <- plyr::rbind.fill(general_results_df, data.frame(new_row, check.names = F))
      
      if (grepl('Maaslin3', tool) & generator == 'SD2') {
        tmp_truth <- truth[truth$associations == 'abundance',]
        tmp_associations <- prepare_associations_maaslin3(associations, tool)
        tmp_associations <- tmp_associations[tmp_associations$association == 'abundance',]
        new_row <- c(unweighted_precision_recall_maaslin3(tmp_truth, tmp_associations, abundance, metadata),
                     weighted_precision_recall_maaslin3(tmp_truth, tmp_associations, abundance, metadata),
                     effect_size_error_maaslin3(tmp_truth, tmp_associations),
                     effect_size_correlation_maaslin3(tmp_truth, tmp_associations),
                     pval_auc_maaslin3(tmp_truth, tmp_associations, abundance, metadata),
                     weighted_pval_auc_maaslin3(tmp_truth, tmp_associations, abundance, metadata),
                     issue_prop(associations[associations$associations == 'abundance',], tool))
        new_row <- lapply(new_row, function(x) {x})
        names(new_row) <- metrics
        
        new_row[['tool']] <- tool
        new_row[['iter']] <- i
        new_row[['association_type']] <- 'abundance'
        new_row <- c(new_row, param_list)
        maaslin3_results_df <- plyr::rbind.fill(maaslin3_results_df, data.frame(new_row, check.names = F))
        
        tmp_truth <- truth[truth$associations == 'prevalence',]
        tmp_associations <- prepare_associations_maaslin3(associations, tool)
        tmp_associations <- tmp_associations[tmp_associations$association == 'prevalence',]
        new_row <- c(unweighted_precision_recall_maaslin3(tmp_truth, tmp_associations, abundance, metadata),
                     weighted_precision_recall_maaslin3(tmp_truth, tmp_associations, abundance, metadata),
                     effect_size_error_maaslin3(tmp_truth, tmp_associations),
                     effect_size_correlation_maaslin3(tmp_truth, tmp_associations),
                     pval_auc_maaslin3(tmp_truth, tmp_associations, abundance, metadata),
                     weighted_pval_auc_maaslin3(tmp_truth, tmp_associations, abundance, metadata),
                     issue_prop(associations[associations$associations == 'prevalence',], tool))
        new_row <- lapply(new_row, function(x) {x})
        names(new_row) <- metrics
        
        new_row[['tool']] <- tool
        new_row[['iter']] <- i
        new_row[['association_type']] <- 'prevalence'
        new_row <- c(new_row, param_list)
        maaslin3_results_df <- plyr::rbind.fill(maaslin3_results_df, data.frame(new_row, check.names = F))
      }
    }
  }
}

names_to_iter <- names(param_list)[names(param_list) != 'metadataType']
for (metadata_type in metadata_types) {
  figures_folder <- paste0(workingDirectory, '/Figures/unscaled/', generator, '/', metadata_type, '/')
  dir.create(figures_folder, recursive = T, showWarnings = F)
  for (param_name in names_to_iter) {
    default_vals <- first_params
    names(default_vals) <- names_to_iter
    to_drop <- which(names(default_vals) == param_name)
    if (length(to_drop) > 0) {
      default_vals <- default_vals[-to_drop]
    }
    
    results_subset <- general_results_df[apply(general_results_df[,names(default_vals)], 1, function(x) {all(x == default_vals)}),]
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
    ggsave(paste0(figures_folder, 'general_', param_name, '.png'),
           plot = plot_out, width = 18, height = 9)
  }
}

#####################
# In-text variation #
#####################

out_diffs <- c()
for (param_name in c("effectSize", "nSubjects", "nPerSubject", "effectPos", "nMetadata", "nMicrobes", "readDepth", "spikeMicrobes")) {
  default_vals <- first_params
  names(default_vals) <- names_to_iter
  default_vals <- default_vals[-which(names(default_vals) == param_name)]
  
  results_subset <- general_results_df[apply(general_results_df[,names(default_vals)], 1, function(x) {all(x == default_vals)}),]
  results_subset <- results_subset[results_subset$metadataType == metadata_type,]
  results_subset <- results_subset[,c(metrics, param_name, 'tool', 'iter')]
  
  melted_df <- melt(results_subset, id.vars = c('tool', 'iter', param_name))
  
  melted_df[,param_name] <- factor(melted_df[,param_name], levels = as.character(sort(unique(as.numeric(melted_df[,param_name])))))
  melted_df <- melted_df[melted_df$variable != "Issue proportion",]
  melted_df$variable <- case_when(melted_df$variable == 'Precision' ~ 'Precision',
                                  melted_df$variable == 'Recall' ~ 'Recall',
                                  melted_df$variable == 'Precision\n(common taxa)' ~ 'Precision (common taxa)',
                                  melted_df$variable == 'Recall\n(common taxa)' ~ 'Recall (common taxa)',
                                  melted_df$variable == 'AUC' ~ 'AUC',
                                  melted_df$variable == 'AUC (Common taxa)' ~ 'AUC (common taxa)',
                                  melted_df$variable == 'Relative effect error' ~ 'Relative effect error',
                                  melted_df$variable == 'Relative shrinkage error' ~ 'Relative shrinkage error',
                                  melted_df$variable == 'Effect size\nSpearman cor.' ~ 'Effect size Spearman cor.')
  melted_df$variable <- factor(melted_df$variable, levels = c('Precision', 
                                                              'Recall',
                                                              'AUC',
                                                              'Precision (common taxa)',
                                                              'Recall (common taxa)',
                                                              'AUC (common taxa)',
                                                              'Relative effect error',
                                                              'Relative shrinkage error',
                                                              'Effect size Spearman cor.'))
  melted_df$tool <- case_when(melted_df$tool == 'ALDEx2' ~ 'ALDEx2',
                              melted_df$tool == 'ANCOMBC' ~ 'ANCOM-BC2',
                              melted_df$tool == 'Maaslin2' ~ 'MaAsLin 2',
                              melted_df$tool == 'Maaslin3' ~ 'MaAsLin 3 Base',
                              melted_df$tool == 'Maaslin3itaug' ~ 'MaAsLin 3',
                              melted_df$tool == 'Maaslin3unscaled' ~ 'MaAsLin 3 Inferred\nAbundance',
                              melted_df$tool == 'Maaslin3correct' ~ 'MaAsLin 3 New\nCorrection')
  
  first_agg <- aggregate(value ~ variable + get(param_name) + tool, data = melted_df, FUN = mean)
  diffs <- abs(aggregate(value ~ variable + tool, data = first_agg, FUN = max)$value - 
    aggregate(value ~ variable + tool, data = first_agg, FUN = min)$value)
  print(mean(diffs))
  out_diffs <- c(out_diffs, diffs)
}

figure_1 <- function() {
  metadata_type <- "MVAref"
  param_name <- "nSubjects"
  
  general_results_df <- data.frame(matrix(nrow = 0, ncol = 0))
  maaslin3_results_df <- data.frame(matrix(nrow = 0, ncol = 0))
  metrics <- c("Precision", "Recall", "Precision\n(common taxa)", "Recall\n(common taxa)", 
               "Relative effect error", "Relative shrinkage error", "Effect size\nSpearman cor.", "AUC", "AUC (Common taxa)", 
               "Issue proportion")
  for (param_list in param_list_final[c(2, 10, 12, 14, 16, 18)]) { #c(2, 14, 17, 20, 23, 26)
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
                       sep='_')
    options("scipen"=5)
    
    inputDirectory <- file.path(workingDirectory, 'Input', 'unscaled', generator)
    this_params_folder <- file.path(inputDirectory, inputString)
    outputDirectory <- file.path(workingDirectory, 'Output', 'unscaled', generator)
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
      
      truth <- prepare_truth(truth, generator)
      for (tool in tools) {
        possible_error <- tryCatch({
          associations <- read.csv(paste0(this_output_folder, '_', tool, "/associations_", i, ".tsv"), sep = "\t")
        }, error = function(err) {
          err
        })
        
        if(inherits(possible_error, "error")) next
        
        new_row <- c(unweighted_precision_recall(truth, prepare_associations_general(associations, tool, generator), abundance, metadata),
                     weighted_precision_recall(truth, prepare_associations_general(associations, tool, generator), abundance, metadata),
                     effect_size_error(truth, prepare_associations_abundance(associations, tool, generator)),
                     effect_size_correlation(truth, prepare_associations_abundance(associations, tool, generator)),
                     pval_auc(truth, prepare_associations_general(associations, tool, generator), abundance, metadata),
                     weighted_pval_auc(truth, prepare_associations_general(associations, tool, generator), abundance, metadata),
                     issue_prop(associations, tool))
        new_row <- lapply(new_row, function(x) {x})
        names(new_row) <- metrics
        
        new_row[['tool']] <- tool
        new_row[['iter']] <- i
        new_row <- c(new_row, param_list)
        general_results_df <- plyr::rbind.fill(general_results_df, data.frame(new_row, check.names = F))
        
        if (grepl('Maaslin3', tool) & generator == 'SD2') {
          tmp_truth <- truth[truth$associations == 'abundance',]
          tmp_associations <- prepare_associations_maaslin3(associations, tool)
          tmp_associations <- tmp_associations[tmp_associations$association == 'abundance',]
          new_row <- c(unweighted_precision_recall_maaslin3(tmp_truth, tmp_associations, abundance, metadata),
                       weighted_precision_recall_maaslin3(tmp_truth, tmp_associations, abundance, metadata),
                       effect_size_error_maaslin3(tmp_truth, tmp_associations),
                       effect_size_correlation_maaslin3(tmp_truth, tmp_associations),
                       pval_auc_maaslin3(tmp_truth, tmp_associations, abundance, metadata),
                       weighted_pval_auc_maaslin3(tmp_truth, tmp_associations, abundance, metadata),
                       issue_prop(associations[associations$associations == 'abundance',], tool))
          new_row <- lapply(new_row, function(x) {x})
          names(new_row) <- metrics
          
          new_row[['tool']] <- tool
          new_row[['iter']] <- i
          new_row[['association_type']] <- 'abundance'
          new_row <- c(new_row, param_list)
          maaslin3_results_df <- plyr::rbind.fill(maaslin3_results_df, data.frame(new_row, check.names = F))
          
          tmp_truth <- truth[truth$associations == 'prevalence',]
          tmp_associations <- prepare_associations_maaslin3(associations, tool)
          tmp_associations <- tmp_associations[tmp_associations$association == 'prevalence',]
          new_row <- c(unweighted_precision_recall_maaslin3(tmp_truth, tmp_associations, abundance, metadata),
                       weighted_precision_recall_maaslin3(tmp_truth, tmp_associations, abundance, metadata),
                       effect_size_error_maaslin3(tmp_truth, tmp_associations),
                       effect_size_correlation_maaslin3(tmp_truth, tmp_associations),
                       pval_auc_maaslin3(tmp_truth, tmp_associations, abundance, metadata),
                       weighted_pval_auc_maaslin3(tmp_truth, tmp_associations, abundance, metadata),
                       issue_prop(associations[associations$associations == 'prevalence',], tool))
          new_row <- lapply(new_row, function(x) {x})
          names(new_row) <- metrics
          
          new_row[['tool']] <- tool
          new_row[['iter']] <- i
          new_row[['association_type']] <- 'prevalence'
          new_row <- c(new_row, param_list)
          maaslin3_results_df <- plyr::rbind.fill(maaslin3_results_df, data.frame(new_row, check.names = F))
        }
      }
    }
  }
  
  names_to_iter <- names(param_list)[names(param_list) != 'metadataType']
  figures_folder <- paste0(workingDirectory, '/Figures/thesis_figures/')
  dir.create(figures_folder, recursive = T, showWarnings = F)
  
  default_vals <- first_params
  names(default_vals) <- names_to_iter
  default_vals <- default_vals[-which(names(default_vals) == param_name)]
  
  results_subset <- general_results_df[apply(general_results_df[,names(default_vals)], 1, function(x) {all(x == default_vals)}),]
  results_subset <- results_subset[results_subset$metadataType == metadata_type,]
  results_subset <- results_subset[,c(metrics, param_name, 'tool', 'iter')]
  
  melted_df <- melt(results_subset, id.vars = c('tool', 'iter', param_name))
  melted_df[,param_name] <- factor(melted_df[,param_name], levels = as.character(sort(unique(as.numeric(melted_df[,param_name])))))
  melted_df <- melted_df[melted_df$variable != "Issue proportion",]
  melted_df$variable <- case_when(melted_df$variable == 'Precision' ~ 'Precision',
                                  melted_df$variable == 'Recall' ~ 'Recall',
                                  melted_df$variable == 'Precision\n(common taxa)' ~ 'Precision (common taxa)',
                                  melted_df$variable == 'Recall\n(common taxa)' ~ 'Recall (common taxa)',
                                  melted_df$variable == 'AUC' ~ 'AUC',
                                  melted_df$variable == 'AUC (Common taxa)' ~ 'AUC (common taxa)',
                                  melted_df$variable == 'Relative effect error' ~ 'Relative effect error',
                                  melted_df$variable == 'Relative shrinkage error' ~ 'Relative shrinkage error',
                                  melted_df$variable == 'Effect size\nSpearman cor.' ~ 'Effect size Spearman cor.')
  melted_df$variable <- factor(melted_df$variable, levels = c('Precision', 
                                                              'Recall',
                                                              'AUC',
                                                              'Precision (common taxa)',
                                                              'Recall (common taxa)',
                                                              'AUC (common taxa)',
                                                              'Relative effect error',
                                                              'Relative shrinkage error',
                                                              'Effect size Spearman cor.'))
  melted_df$tool <- case_when(melted_df$tool == 'ALDEx2' ~ 'ALDEx2',
                              melted_df$tool == 'ANCOMBC' ~ 'ANCOM-BC2',
                              melted_df$tool == 'Maaslin2' ~ 'MaAsLin 2',
                              melted_df$tool == 'Maaslin3' ~ 'MaAsLin 3 Base',
                              melted_df$tool == 'Maaslin3itaug' ~ 'MaAsLin 3',
                              melted_df$tool == 'Maaslin3unscaled' ~ 'MaAsLin 3 Inferred\nAbundance',
                              melted_df$tool == 'Maaslin3correct' ~ 'MaAsLin 3 Sparsity\nCorrected')
  melted_df$tool <- factor(melted_df$tool, levels = c("ALDEx2", "ANCOM-BC2", "MaAsLin 2", "MaAsLin 3 Base", "MaAsLin 3", "MaAsLin 3 Inferred\nAbundance", "MaAsLin 3 Sparsity\nCorrected"))
  
  # In-text
  aggregate(value ~ tool, data = melted_df[melted_df$variable == 'Precision' & melted_df$nSubjects %in% c(20,50),], FUN = mean)
  aggregate(value ~ tool, data = melted_df[melted_df$variable == 'Precision' & melted_df$nSubjects %in% 1000,], FUN = mean)
  aggregate(value ~ tool, data = melted_df[melted_df$variable == 'Recall' & melted_df$nSubjects %in% 50,], FUN = mean)
  aggregate(value ~ tool, data = melted_df[melted_df$variable == 'Recall' & melted_df$nSubjects %in% 500,], FUN = mean)
  aggregate(value ~ tool, data = melted_df[melted_df$variable == 'AUC' & melted_df$nSubjects %in% 20,], FUN = mean)
  mean((aggregate(value ~ tool + nSubjects, data = melted_df[melted_df$variable == 'Recall (common taxa)',], FUN = mean) - 
    aggregate(value ~ tool + nSubjects, data = melted_df[melted_df$variable == 'Recall',], FUN = mean))$value)
  aggregate(value ~ tool, data = melted_df[melted_df$variable == 'Relative shrinkage error' & melted_df$nSubjects %in% 1000,], FUN = mean)
  aggregate(value ~ tool, data = melted_df[melted_df$variable == 'Effect size Spearman cor.' & melted_df$nSubjects %in% 100,], FUN = mean)
  
  plot_out <- ggplot(melted_df, aes(x = get(param_name), y = value, fill = tool)) + 
    geom_boxplot(position = position_dodge(preserve = "single")) + 
    facet_wrap(~variable, scales = 'free') + 
    theme_bw() + 
    xlab("Number of Subjects") + 
    ylab('') + 
    theme(text=element_text(size=21),
          legend.position = 'bottom') + 
    labs(fill = 'Model') + 
    scale_fill_manual(values=c("#104E8B", "#458B00", "#EEAD0E", "#EE7600", "#8B1A1A", "#68228B", "#9852CB"),
                      breaks=c("ALDEx2", "ANCOM-BC2", "MaAsLin 2", "MaAsLin 3 Base", "MaAsLin 3", "MaAsLin 3 Inferred\nAbundance", "MaAsLin 3 Sparsity\nCorrected"))
  ggsave(paste0(figures_folder, 'fig_1.png'),
         plot = plot_out, width = 12, height = 12)
}

figure_2 <- function() {
  metadata_type <- "MVAref"
  param_name <- "nPerSubject"
  
  general_results_df <- data.frame(matrix(nrow = 0, ncol = 0))
  maaslin3_results_df <- data.frame(matrix(nrow = 0, ncol = 0))
  metrics <- c("Precision", "Recall", "Precision\n(common taxa)", "Recall\n(common taxa)", 
               "Relative effect error", "Relative shrinkage error", "Effect size\nSpearman cor.", "AUC", "AUC (Common taxa)", 
               "Issue proportion")
  for (param_list in param_list_final[c(2, 42, 44, 46)]) { #c(2, 62, 65, 68)
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
                       sep='_')
    options("scipen"=5)
    
    inputDirectory <- file.path(workingDirectory, 'Input', 'unscaled', generator)
    this_params_folder <- file.path(inputDirectory, inputString)
    outputDirectory <- file.path(workingDirectory, 'Output', 'unscaled', generator)
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
      
      truth <- prepare_truth(truth, generator)
      for (tool in tools) {
        possible_error <- tryCatch({
          associations <- read.csv(paste0(this_output_folder, '_', tool, "/associations_", i, ".tsv"), sep = "\t")
        }, error = function(err) {
          err
        })
        
        if(inherits(possible_error, "error")) next
        
        new_row <- c(unweighted_precision_recall(truth, prepare_associations_general(associations, tool, generator), abundance, metadata),
                     weighted_precision_recall(truth, prepare_associations_general(associations, tool, generator), abundance, metadata),
                     effect_size_error(truth, prepare_associations_abundance(associations, tool, generator)),
                     effect_size_correlation(truth, prepare_associations_abundance(associations, tool, generator)),
                     pval_auc(truth, prepare_associations_general(associations, tool, generator), abundance, metadata),
                     weighted_pval_auc(truth, prepare_associations_general(associations, tool, generator), abundance, metadata),
                     issue_prop(associations, tool))
        new_row <- lapply(new_row, function(x) {x})
        names(new_row) <- metrics
        
        new_row[['tool']] <- tool
        new_row[['iter']] <- i
        new_row <- c(new_row, param_list)
        general_results_df <- plyr::rbind.fill(general_results_df, data.frame(new_row, check.names = F))
        
        if (grepl('Maaslin3', tool) & generator == 'SD2') {
          tmp_truth <- truth[truth$associations == 'abundance',]
          tmp_associations <- prepare_associations_maaslin3(associations, tool)
          tmp_associations <- tmp_associations[tmp_associations$association == 'abundance',]
          new_row <- c(unweighted_precision_recall_maaslin3(tmp_truth, tmp_associations, abundance, metadata),
                       weighted_precision_recall_maaslin3(tmp_truth, tmp_associations, abundance, metadata),
                       effect_size_error_maaslin3(tmp_truth, tmp_associations),
                       effect_size_correlation_maaslin3(tmp_truth, tmp_associations),
                       pval_auc_maaslin3(tmp_truth, tmp_associations, abundance, metadata),
                       weighted_pval_auc_maaslin3(tmp_truth, tmp_associations, abundance, metadata),
                       issue_prop(associations[associations$associations == 'abundance',], tool))
          new_row <- lapply(new_row, function(x) {x})
          names(new_row) <- metrics
          
          new_row[['tool']] <- tool
          new_row[['iter']] <- i
          new_row[['association_type']] <- 'abundance'
          new_row <- c(new_row, param_list)
          maaslin3_results_df <- plyr::rbind.fill(maaslin3_results_df, data.frame(new_row, check.names = F))
          
          tmp_truth <- truth[truth$associations == 'prevalence',]
          tmp_associations <- prepare_associations_maaslin3(associations, tool)
          tmp_associations <- tmp_associations[tmp_associations$association == 'prevalence',]
          new_row <- c(unweighted_precision_recall_maaslin3(tmp_truth, tmp_associations, abundance, metadata),
                       weighted_precision_recall_maaslin3(tmp_truth, tmp_associations, abundance, metadata),
                       effect_size_error_maaslin3(tmp_truth, tmp_associations),
                       effect_size_correlation_maaslin3(tmp_truth, tmp_associations),
                       pval_auc_maaslin3(tmp_truth, tmp_associations, abundance, metadata),
                       weighted_pval_auc_maaslin3(tmp_truth, tmp_associations, abundance, metadata),
                       issue_prop(associations[associations$associations == 'prevalence',], tool))
          new_row <- lapply(new_row, function(x) {x})
          names(new_row) <- metrics
          
          new_row[['tool']] <- tool
          new_row[['iter']] <- i
          new_row[['association_type']] <- 'prevalence'
          new_row <- c(new_row, param_list)
          maaslin3_results_df <- plyr::rbind.fill(maaslin3_results_df, data.frame(new_row, check.names = F))
        }
      }
    }
  }
  
  names_to_iter <- names(param_list)[names(param_list) != 'metadataType']
  figures_folder <- paste0(workingDirectory, '/Figures/thesis_figures/')
  dir.create(figures_folder, recursive = T, showWarnings = F)
  
  default_vals <- first_params
  names(default_vals) <- names_to_iter
  default_vals <- default_vals[-which(names(default_vals) == param_name)]
  
  results_subset <- general_results_df[apply(general_results_df[,names(default_vals)], 1, function(x) {all(x == default_vals)}),]
  results_subset <- results_subset[results_subset$metadataType == metadata_type,]
  results_subset <- results_subset[,c(metrics, param_name, 'tool', 'iter')]
  
  melted_df <- melt(results_subset, id.vars = c('tool', 'iter', param_name))
  melted_df[,param_name] <- factor(melted_df[,param_name], levels = as.character(sort(unique(as.numeric(melted_df[,param_name])))))
  melted_df <- melted_df[melted_df$variable != "Issue proportion",]
  melted_df$variable <- case_when(melted_df$variable == 'Precision' ~ 'Precision',
                                  melted_df$variable == 'Recall' ~ 'Recall',
                                  melted_df$variable == 'Precision\n(common taxa)' ~ 'Precision (common taxa)',
                                  melted_df$variable == 'Recall\n(common taxa)' ~ 'Recall (common taxa)',
                                  melted_df$variable == 'AUC' ~ 'AUC',
                                  melted_df$variable == 'AUC (Common taxa)' ~ 'AUC (common taxa)',
                                  melted_df$variable == 'Relative effect error' ~ 'Relative effect error',
                                  melted_df$variable == 'Relative shrinkage error' ~ 'Relative shrinkage error',
                                  melted_df$variable == 'Effect size\nSpearman cor.' ~ 'Effect size Spearman cor.')
  melted_df$variable <- factor(melted_df$variable, levels = c('Precision', 
                                                              'Recall',
                                                              'AUC',
                                                              'Precision (common taxa)',
                                                              'Recall (common taxa)',
                                                              'AUC (common taxa)',
                                                              'Relative effect error',
                                                              'Relative shrinkage error',
                                                              'Effect size Spearman cor.'))
  melted_df$tool <- case_when(melted_df$tool == 'ALDEx2' ~ 'ALDEx2',
                              melted_df$tool == 'ANCOMBC' ~ 'ANCOM-BC2',
                              melted_df$tool == 'Maaslin2' ~ 'MaAsLin 2',
                              melted_df$tool == 'Maaslin3' ~ 'MaAsLin 3 Base',
                              melted_df$tool == 'Maaslin3itaug' ~ 'MaAsLin 3',
                              melted_df$tool == 'Maaslin3unscaled' ~ 'MaAsLin 3 Inferred\nAbundance',
                              melted_df$tool == 'Maaslin3correct' ~ 'MaAsLin 3 Sparsity\nCorrected')
  
  melted_df$tool <- factor(melted_df$tool, levels = c("ALDEx2", "ANCOM-BC2", "MaAsLin 2", "MaAsLin 3 Base", "MaAsLin 3", "MaAsLin 3 Inferred\nAbundance", "MaAsLin 3 Sparsity\nCorrected"))
  
  # In-text
  aggregate(value ~ tool + nPerSubject, data = melted_df[melted_df$variable == 'Precision',], FUN = mean)
  aggregate(value ~ tool + nPerSubject, data = melted_df[melted_df$variable == 'Recall',], FUN = mean)
  aggregate(value ~ tool + nPerSubject, data = melted_df[melted_df$variable == 'Recall (common taxa)',], FUN = mean)
  
  plot_out <- ggplot(melted_df, aes(x = get(param_name), y = value, fill = tool)) + 
    geom_boxplot(position = position_dodge(preserve = "single")) + 
    facet_wrap(~variable, scales = 'free') + 
    theme_bw() + 
    xlab("Samples per Subject") + 
    ylab('') + 
    theme(text=element_text(size=21),
          legend.position = 'bottom') + 
    labs(fill = 'Model') + 
    scale_fill_manual(values=c("#104E8B", "#458B00", "#EEAD0E", "#EE7600", "#8B1A1A", "#68228B", "#9852CB"),
                      breaks=c("ALDEx2", "ANCOM-BC2", "MaAsLin 2", "MaAsLin 3 Base", "MaAsLin 3", "MaAsLin 3 Inferred\nAbundance", "MaAsLin 3 Sparsity\nCorrected"))
  ggsave(paste0(figures_folder, 'fig_2.png'),
         plot = plot_out, width = 12, height = 12)
}

figure_3 <- function() {
  metadata_type <- "MVAref"
  param_name <- "effectSize"
  
  general_results_df <- data.frame(matrix(nrow = 0, ncol = 0))
  maaslin3_results_df <- data.frame(matrix(nrow = 0, ncol = 0))
  metrics <- c("Precision", "Recall", "Precision\n(common taxa)", "Recall\n(common taxa)", 
               "Relative effect error", "Relative shrinkage error", "Effect size\nSpearman cor.", "AUC", "AUC (Common taxa)", 
               "Issue proportion")
  for (param_list in param_list_final[c(2, 24, 26, 28, 30)]) { #c(2, 35, 38, 41, 44)
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
                       sep='_')
    options("scipen"=5)
    
    inputDirectory <- file.path(workingDirectory, 'Input', 'unscaled', generator)
    this_params_folder <- file.path(inputDirectory, inputString)
    outputDirectory <- file.path(workingDirectory, 'Output', 'unscaled', generator)
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
      
      truth <- prepare_truth(truth, generator)
      for (tool in tools) {
        possible_error <- tryCatch({
          associations <- read.csv(paste0(this_output_folder, '_', tool, "/associations_", i, ".tsv"), sep = "\t")
        }, error = function(err) {
          err
        })
        
        if(inherits(possible_error, "error")) next
        
        new_row <- c(unweighted_precision_recall(truth, prepare_associations_general(associations, tool, generator), abundance, metadata),
                     weighted_precision_recall(truth, prepare_associations_general(associations, tool, generator), abundance, metadata),
                     effect_size_error(truth, prepare_associations_abundance(associations, tool, generator)),
                     effect_size_correlation(truth, prepare_associations_abundance(associations, tool, generator)),
                     pval_auc(truth, prepare_associations_general(associations, tool, generator), abundance, metadata),
                     weighted_pval_auc(truth, prepare_associations_general(associations, tool, generator), abundance, metadata),
                     issue_prop(associations, tool))
        new_row <- lapply(new_row, function(x) {x})
        names(new_row) <- metrics
        
        new_row[['tool']] <- tool
        new_row[['iter']] <- i
        new_row <- c(new_row, param_list)
        general_results_df <- plyr::rbind.fill(general_results_df, data.frame(new_row, check.names = F))
        
        if (grepl('Maaslin3', tool) & generator == 'SD2') {
          tmp_truth <- truth[truth$associations == 'abundance',]
          tmp_associations <- prepare_associations_maaslin3(associations, tool)
          tmp_associations <- tmp_associations[tmp_associations$association == 'abundance',]
          new_row <- c(unweighted_precision_recall_maaslin3(tmp_truth, tmp_associations, abundance, metadata),
                       weighted_precision_recall_maaslin3(tmp_truth, tmp_associations, abundance, metadata),
                       effect_size_error_maaslin3(tmp_truth, tmp_associations),
                       effect_size_correlation_maaslin3(tmp_truth, tmp_associations),
                       pval_auc_maaslin3(tmp_truth, tmp_associations, abundance, metadata),
                       weighted_pval_auc_maaslin3(tmp_truth, tmp_associations, abundance, metadata),
                       issue_prop(associations[associations$associations == 'abundance',], tool))
          new_row <- lapply(new_row, function(x) {x})
          names(new_row) <- metrics
          
          new_row[['tool']] <- tool
          new_row[['iter']] <- i
          new_row[['association_type']] <- 'abundance'
          new_row <- c(new_row, param_list)
          maaslin3_results_df <- plyr::rbind.fill(maaslin3_results_df, data.frame(new_row, check.names = F))
          
          tmp_truth <- truth[truth$associations == 'prevalence',]
          tmp_associations <- prepare_associations_maaslin3(associations, tool)
          tmp_associations <- tmp_associations[tmp_associations$association == 'prevalence',]
          new_row <- c(unweighted_precision_recall_maaslin3(tmp_truth, tmp_associations, abundance, metadata),
                       weighted_precision_recall_maaslin3(tmp_truth, tmp_associations, abundance, metadata),
                       effect_size_error_maaslin3(tmp_truth, tmp_associations),
                       effect_size_correlation_maaslin3(tmp_truth, tmp_associations),
                       pval_auc_maaslin3(tmp_truth, tmp_associations, abundance, metadata),
                       weighted_pval_auc_maaslin3(tmp_truth, tmp_associations, abundance, metadata),
                       issue_prop(associations[associations$associations == 'prevalence',], tool))
          new_row <- lapply(new_row, function(x) {x})
          names(new_row) <- metrics
          
          new_row[['tool']] <- tool
          new_row[['iter']] <- i
          new_row[['association_type']] <- 'prevalence'
          new_row <- c(new_row, param_list)
          maaslin3_results_df <- plyr::rbind.fill(maaslin3_results_df, data.frame(new_row, check.names = F))
        }
      }
    }
  }
  
  names_to_iter <- names(param_list)[names(param_list) != 'metadataType']
  figures_folder <- paste0(workingDirectory, '/Figures/thesis_figures/')
  dir.create(figures_folder, recursive = T, showWarnings = F)
  
  default_vals <- first_params
  names(default_vals) <- names_to_iter
  default_vals <- default_vals[-which(names(default_vals) == param_name)]
  
  results_subset <- general_results_df[apply(general_results_df[,names(default_vals)], 1, function(x) {all(x == default_vals)}),]
  results_subset <- results_subset[results_subset$metadataType == metadata_type,]
  results_subset <- results_subset[,c(metrics, param_name, 'tool', 'iter')]
  
  melted_df <- melt(results_subset, id.vars = c('tool', 'iter', param_name))
  melted_df[,param_name] <- factor(melted_df[,param_name], levels = as.character(sort(unique(as.numeric(melted_df[,param_name])))))
  melted_df <- melted_df[melted_df$variable != "Issue proportion",]
  melted_df$variable <- case_when(melted_df$variable == 'Precision' ~ 'Precision',
                                  melted_df$variable == 'Recall' ~ 'Recall',
                                  melted_df$variable == 'Precision\n(common taxa)' ~ 'Precision (common taxa)',
                                  melted_df$variable == 'Recall\n(common taxa)' ~ 'Recall (common taxa)',
                                  melted_df$variable == 'AUC' ~ 'AUC',
                                  melted_df$variable == 'AUC (Common taxa)' ~ 'AUC (common taxa)',
                                  melted_df$variable == 'Relative effect error' ~ 'Relative effect error',
                                  melted_df$variable == 'Relative shrinkage error' ~ 'Relative shrinkage error',
                                  melted_df$variable == 'Effect size\nSpearman cor.' ~ 'Effect size Spearman cor.')
  melted_df$variable <- factor(melted_df$variable, levels = c('Precision', 
                                                              'Recall',
                                                              'AUC',
                                                              'Precision (common taxa)',
                                                              'Recall (common taxa)',
                                                              'AUC (common taxa)',
                                                              'Relative effect error',
                                                              'Relative shrinkage error',
                                                              'Effect size Spearman cor.'))
  melted_df$tool <- case_when(melted_df$tool == 'ALDEx2' ~ 'ALDEx2',
                              melted_df$tool == 'ANCOMBC' ~ 'ANCOM-BC2',
                              melted_df$tool == 'Maaslin2' ~ 'MaAsLin 2',
                              melted_df$tool == 'Maaslin3' ~ 'MaAsLin 3 Base',
                              melted_df$tool == 'Maaslin3itaug' ~ 'MaAsLin 3',
                              melted_df$tool == 'Maaslin3unscaled' ~ 'MaAsLin 3 Inferred\nAbundance',
                              melted_df$tool == 'Maaslin3correct' ~ 'MaAsLin 3 Sparsity\nCorrected')
  
  melted_df$tool <- factor(melted_df$tool, levels = c("ALDEx2", "ANCOM-BC2", "MaAsLin 2", "MaAsLin 3 Base", "MaAsLin 3", "MaAsLin 3 Inferred\nAbundance", "MaAsLin 3 Sparsity\nCorrected"))
  
  # In-text
  aggregate(value ~ tool, data = melted_df[melted_df$variable == 'Relative shrinkage error' & melted_df$effectSize == 0.5,], FUN = mean)
  aggregate(value ~ tool + effectSize, data = melted_df[melted_df$variable == 'Precision',], FUN = mean)
  aggregate(value ~ tool + effectSize, data = melted_df[melted_df$variable == 'Recall',], FUN = mean)
  aggregate(value ~ tool + effectSize, data = melted_df[melted_df$variable == 'Relative effect error',], FUN = mean)
  aggregate(value ~ tool + effectSize, data = melted_df[melted_df$variable == 'Effect size Spearman cor.',], FUN = mean)
  
  plot_out <- ggplot(melted_df, aes(x = get(param_name), y = value, fill = tool)) + 
    geom_boxplot(position = position_dodge(preserve = "single")) + 
    facet_wrap(~variable, scales = 'free') + 
    theme_bw() + 
    xlab("Effect Size") + 
    ylab('') + 
    theme(text=element_text(size=21),
          legend.position = 'bottom') + 
    labs(fill = 'Model') + 
    scale_fill_manual(values=c("#104E8B", "#458B00", "#EEAD0E", "#EE7600", "#8B1A1A", "#68228B", "#9852CB"),
                      breaks=c("ALDEx2", "ANCOM-BC2", "MaAsLin 2", "MaAsLin 3 Base", "MaAsLin 3", "MaAsLin 3 Inferred\nAbundance", "MaAsLin 3 Sparsity\nCorrected"))
  ggsave(paste0(figures_folder, 'fig_3.png'),
         plot = plot_out, width = 12, height = 12)
}

figure_4 <- function() {
  metadata_type <- "MVAref"
  param_name <- "nSubjects"
  
  general_results_df <- data.frame(matrix(nrow = 0, ncol = 0))
  maaslin3_results_df <- data.frame(matrix(nrow = 0, ncol = 0))
  metrics <- c("Precision", "Recall", "Precision\n(common taxa)", "Recall\n(common taxa)", 
               "Relative effect error", "Relative shrinkage error", "Effect size\nSpearman cor.", "AUC", "AUC (Common taxa)", 
               "Issue proportion")
  for (param_list in param_list_final[c(2, 10, 12, 14, 16, 18)]) { #c(2, 14, 17, 20, 23, 26)
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
                       sep='_')
    options("scipen"=5)
    
    inputDirectory <- file.path(workingDirectory, 'Input', 'unscaled', generator)
    this_params_folder <- file.path(inputDirectory, inputString)
    outputDirectory <- file.path(workingDirectory, 'Output', 'unscaled', generator)
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
      
      truth <- prepare_truth(truth, generator)
      for (tool in tools) {
        possible_error <- tryCatch({
          associations <- read.csv(paste0(this_output_folder, '_', tool, "/associations_", i, ".tsv"), sep = "\t")
        }, error = function(err) {
          err
        })
        
        if(inherits(possible_error, "error")) next
        
        new_row <- c(unweighted_precision_recall(truth, prepare_associations_general(associations, tool, generator), abundance, metadata),
                     weighted_precision_recall(truth, prepare_associations_general(associations, tool, generator), abundance, metadata),
                     effect_size_error(truth, prepare_associations_abundance(associations, tool, generator)),
                     effect_size_correlation(truth, prepare_associations_abundance(associations, tool, generator)),
                     pval_auc(truth, prepare_associations_general(associations, tool, generator), abundance, metadata),
                     weighted_pval_auc(truth, prepare_associations_general(associations, tool, generator), abundance, metadata),
                     issue_prop(associations, tool))
        new_row <- lapply(new_row, function(x) {x})
        names(new_row) <- metrics
        
        new_row[['tool']] <- tool
        new_row[['iter']] <- i
        new_row <- c(new_row, param_list)
        general_results_df <- plyr::rbind.fill(general_results_df, data.frame(new_row, check.names = F))
        
        if (grepl('Maaslin3', tool) & generator == 'SD2') {
          tmp_truth <- truth[truth$associations == 'abundance',]
          tmp_associations <- prepare_associations_maaslin3(associations, tool)
          tmp_associations <- tmp_associations[tmp_associations$association == 'abundance',]
          new_row <- c(unweighted_precision_recall_maaslin3(tmp_truth, tmp_associations, abundance, metadata),
                       weighted_precision_recall_maaslin3(tmp_truth, tmp_associations, abundance, metadata),
                       effect_size_error_maaslin3(tmp_truth, tmp_associations),
                       effect_size_correlation_maaslin3(tmp_truth, tmp_associations),
                       pval_auc_maaslin3(tmp_truth, tmp_associations, abundance, metadata),
                       weighted_pval_auc_maaslin3(tmp_truth, tmp_associations, abundance, metadata),
                       issue_prop(associations[associations$associations == 'abundance',], tool))
          new_row <- lapply(new_row, function(x) {x})
          names(new_row) <- metrics
          
          new_row[['tool']] <- tool
          new_row[['iter']] <- i
          new_row[['association_type']] <- 'abundance'
          new_row <- c(new_row, param_list)
          maaslin3_results_df <- plyr::rbind.fill(maaslin3_results_df, data.frame(new_row, check.names = F))
          
          tmp_truth <- truth[truth$associations == 'prevalence',]
          tmp_associations <- prepare_associations_maaslin3(associations, tool)
          tmp_associations <- tmp_associations[tmp_associations$association == 'prevalence',]
          new_row <- c(unweighted_precision_recall_maaslin3(tmp_truth, tmp_associations, abundance, metadata),
                       weighted_precision_recall_maaslin3(tmp_truth, tmp_associations, abundance, metadata),
                       effect_size_error_maaslin3(tmp_truth, tmp_associations),
                       effect_size_correlation_maaslin3(tmp_truth, tmp_associations),
                       pval_auc_maaslin3(tmp_truth, tmp_associations, abundance, metadata),
                       weighted_pval_auc_maaslin3(tmp_truth, tmp_associations, abundance, metadata),
                       issue_prop(associations[associations$associations == 'prevalence',], tool))
          new_row <- lapply(new_row, function(x) {x})
          names(new_row) <- metrics
          
          new_row[['tool']] <- tool
          new_row[['iter']] <- i
          new_row[['association_type']] <- 'prevalence'
          new_row <- c(new_row, param_list)
          maaslin3_results_df <- plyr::rbind.fill(maaslin3_results_df, data.frame(new_row, check.names = F))
        }
      }
    }
  }
  
  names_to_iter <- names(param_list)[names(param_list) != 'metadataType']
  figures_folder <- paste0(workingDirectory, '/Figures/thesis_figures/')
  dir.create(figures_folder, recursive = T, showWarnings = F)
  
  default_vals <- first_params
  names(default_vals) <- names_to_iter
  default_vals <- default_vals[-which(names(default_vals) == param_name)]
  
  results_subset <- maaslin3_results_df[apply(maaslin3_results_df[,names(default_vals)], 1, function(x) {all(x == default_vals)}),]
  results_subset <- results_subset[results_subset$metadataType == metadata_type,]
  results_subset <- results_subset[,c(metrics, param_name, 'tool', 'iter', 'association_type')]
  
  melted_df <- melt(results_subset, id.vars = c('tool', 'iter', 'association_type', param_name))
  melted_df[,param_name] <- factor(melted_df[,param_name], levels = as.character(sort(unique(as.numeric(melted_df[,param_name])))))
  melted_df <- melted_df[melted_df$variable != "Issue proportion",]
  melted_df$variable <- case_when(melted_df$variable == 'Precision' ~ 'Precision',
                                  melted_df$variable == 'Recall' ~ 'Recall',
                                  melted_df$variable == 'Precision\n(common taxa)' ~ 'Precision (common taxa)',
                                  melted_df$variable == 'Recall\n(common taxa)' ~ 'Recall (common taxa)',
                                  melted_df$variable == 'AUC' ~ 'AUC',
                                  melted_df$variable == 'AUC (Common taxa)' ~ 'AUC (common taxa)',
                                  melted_df$variable == 'Relative effect error' ~ 'Relative effect error',
                                  melted_df$variable == 'Relative shrinkage error' ~ 'Relative shrinkage error',
                                  melted_df$variable == 'Effect size\nSpearman cor.' ~ 'Effect size Spearman cor.')
  melted_df$variable <- factor(melted_df$variable, levels = c('Precision', 
                                                              'Recall',
                                                              'AUC',
                                                              'Precision (common taxa)',
                                                              'Recall (common taxa)',
                                                              'AUC (common taxa)',
                                                              'Relative effect error',
                                                              'Relative shrinkage error',
                                                              'Effect size Spearman cor.'))
  melted_df$tool <- case_when(melted_df$tool == 'ALDEx2' ~ 'ALDEx2',
                              melted_df$tool == 'ANCOMBC' ~ 'ANCOM-BC2',
                              melted_df$tool == 'Maaslin2' ~ 'MaAsLin 2',
                              melted_df$tool == 'Maaslin3' ~ 'MaAsLin 3 Base',
                              melted_df$tool == 'Maaslin3itaug' ~ 'MaAsLin 3',
                              melted_df$tool == 'Maaslin3unscaled' ~ 'MaAsLin 3 Inferred\nAbundance',
                              melted_df$tool == 'Maaslin3correct' ~ 'MaAsLin 3 Sparsity\nCorrected')
  
  melted_df$tool <- paste0(melted_df$tool, " ", ifelse(melted_df$association_type == 'abundance', "Abundance", "Prevalence"))
  
  tool_vec <- c("ALDEx2", "ANCOM-BC2", "MaAsLin 2", "MaAsLin 3 Base", 
                "MaAsLin 3", "MaAsLin 3 Inferred\nAbundance", "MaAsLin 3 Sparsity\nCorrected")
  
  melted_df$tool <- factor(melted_df$tool, levels = c(paste0(tool_vec, " ", "Abundance"),
                                                      paste0(tool_vec, " ", "Prevalence")))
  
  # In-text
  aggregate(value ~ tool, data = melted_df[melted_df$variable == 'Precision' & melted_df$nSubjects == 1000,], FUN = mean)
  aggregate(value ~ tool, data = melted_df[melted_df$variable == 'AUC' & melted_df$nSubjects == 100,], FUN = mean)
  aggregate(value ~ tool, data = melted_df[melted_df$variable == 'Relative shrinkage error' & melted_df$nSubjects == 1000,], FUN = mean)
  aggregate(value ~ tool + nSubjects, data = melted_df[melted_df$variable == 'Effect size Spearman cor.',], FUN = mean)
  
  plot_out <- ggplot(melted_df, aes(x = get(param_name), y = value, fill = tool)) + 
    geom_boxplot(position = position_dodge(preserve = "single")) + 
    facet_wrap(~variable, scales = 'free') + 
    theme_bw() + 
    xlab("Number of Subjects") + 
    ylab('') + 
    theme(text=element_text(size=21),
          legend.position = 'bottom') + 
    labs(fill = 'Model') + 
    scale_fill_manual(values=c("#104E8B", "#458B00", "#EEAD0E", "#EE7600", "#8B1A1A", "#68228B",  "#9852CB",
                               "#88A7C5", "#A2C580", "#F7D687", "#F7BB80", "#C58D8D", "#B491C5",  "#B852CB"),
                      breaks=c(paste0(tool_vec, " ", "Abundance"),
                               paste0(tool_vec, " ", "Prevalence")))
  ggsave(paste0(figures_folder, 'fig_4.png'),
         plot = plot_out, width = 12, height = 12)
}

figure_5 <- function() {
  metadata_type <- "MVAref"
  param_name <- "nPerSubject"
  
  general_results_df <- data.frame(matrix(nrow = 0, ncol = 0))
  maaslin3_results_df <- data.frame(matrix(nrow = 0, ncol = 0))
  metrics <- c("Precision", "Recall", "Precision\n(common taxa)", "Recall\n(common taxa)", 
               "Relative effect error", "Relative shrinkage error", "Effect size\nSpearman cor.", "AUC", "AUC (Common taxa)", 
               "Issue proportion")
  for (param_list in param_list_final[c(2, 42, 44, 46)]) { #c(2, 62, 65, 68)
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
                       sep='_')
    options("scipen"=5)
    
    inputDirectory <- file.path(workingDirectory, 'Input', 'unscaled', generator)
    this_params_folder <- file.path(inputDirectory, inputString)
    outputDirectory <- file.path(workingDirectory, 'Output', 'unscaled', generator)
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
      
      truth <- prepare_truth(truth, generator)
      for (tool in tools) {
        possible_error <- tryCatch({
          associations <- read.csv(paste0(this_output_folder, '_', tool, "/associations_", i, ".tsv"), sep = "\t")
        }, error = function(err) {
          err
        })
        
        if(inherits(possible_error, "error")) next
        
        new_row <- c(unweighted_precision_recall(truth, prepare_associations_general(associations, tool, generator), abundance, metadata),
                     weighted_precision_recall(truth, prepare_associations_general(associations, tool, generator), abundance, metadata),
                     effect_size_error(truth, prepare_associations_abundance(associations, tool, generator)),
                     effect_size_correlation(truth, prepare_associations_abundance(associations, tool, generator)),
                     pval_auc(truth, prepare_associations_general(associations, tool, generator), abundance, metadata),
                     weighted_pval_auc(truth, prepare_associations_general(associations, tool, generator), abundance, metadata),
                     issue_prop(associations, tool))
        new_row <- lapply(new_row, function(x) {x})
        names(new_row) <- metrics
        
        new_row[['tool']] <- tool
        new_row[['iter']] <- i
        new_row <- c(new_row, param_list)
        general_results_df <- plyr::rbind.fill(general_results_df, data.frame(new_row, check.names = F))
        
        if (grepl('Maaslin3', tool) & generator == 'SD2') {
          tmp_truth <- truth[truth$associations == 'abundance',]
          tmp_associations <- prepare_associations_maaslin3(associations, tool)
          tmp_associations <- tmp_associations[tmp_associations$association == 'abundance',]
          new_row <- c(unweighted_precision_recall_maaslin3(tmp_truth, tmp_associations, abundance, metadata),
                       weighted_precision_recall_maaslin3(tmp_truth, tmp_associations, abundance, metadata),
                       effect_size_error_maaslin3(tmp_truth, tmp_associations),
                       effect_size_correlation_maaslin3(tmp_truth, tmp_associations),
                       pval_auc_maaslin3(tmp_truth, tmp_associations, abundance, metadata),
                       weighted_pval_auc_maaslin3(tmp_truth, tmp_associations, abundance, metadata),
                       issue_prop(associations[associations$associations == 'abundance',], tool))
          new_row <- lapply(new_row, function(x) {x})
          names(new_row) <- metrics
          
          new_row[['tool']] <- tool
          new_row[['iter']] <- i
          new_row[['association_type']] <- 'abundance'
          new_row <- c(new_row, param_list)
          maaslin3_results_df <- plyr::rbind.fill(maaslin3_results_df, data.frame(new_row, check.names = F))
          
          tmp_truth <- truth[truth$associations == 'prevalence',]
          tmp_associations <- prepare_associations_maaslin3(associations, tool)
          tmp_associations <- tmp_associations[tmp_associations$association == 'prevalence',]
          new_row <- c(unweighted_precision_recall_maaslin3(tmp_truth, tmp_associations, abundance, metadata),
                       weighted_precision_recall_maaslin3(tmp_truth, tmp_associations, abundance, metadata),
                       effect_size_error_maaslin3(tmp_truth, tmp_associations),
                       effect_size_correlation_maaslin3(tmp_truth, tmp_associations),
                       pval_auc_maaslin3(tmp_truth, tmp_associations, abundance, metadata),
                       weighted_pval_auc_maaslin3(tmp_truth, tmp_associations, abundance, metadata),
                       issue_prop(associations[associations$associations == 'prevalence',], tool))
          new_row <- lapply(new_row, function(x) {x})
          names(new_row) <- metrics
          
          new_row[['tool']] <- tool
          new_row[['iter']] <- i
          new_row[['association_type']] <- 'prevalence'
          new_row <- c(new_row, param_list)
          maaslin3_results_df <- plyr::rbind.fill(maaslin3_results_df, data.frame(new_row, check.names = F))
        }
      }
    }
  }
  
  names_to_iter <- names(param_list)[names(param_list) != 'metadataType']
  figures_folder <- paste0(workingDirectory, '/Figures/thesis_figures/')
  dir.create(figures_folder, recursive = T, showWarnings = F)
  
  default_vals <- first_params
  names(default_vals) <- names_to_iter
  default_vals <- default_vals[-which(names(default_vals) == param_name)]
  
  results_subset <- maaslin3_results_df[apply(maaslin3_results_df[,names(default_vals)], 1, function(x) {all(x == default_vals)}),]
  results_subset <- results_subset[results_subset$metadataType == metadata_type,]
  results_subset <- results_subset[,c(metrics, param_name, 'tool', 'iter', 'association_type')]
  
  melted_df <- melt(results_subset, id.vars = c('tool', 'iter', 'association_type', param_name))
  melted_df[,param_name] <- factor(melted_df[,param_name], levels = as.character(sort(unique(as.numeric(melted_df[,param_name])))))
  melted_df <- melted_df[melted_df$variable != "Issue proportion",]
  melted_df$variable <- case_when(melted_df$variable == 'Precision' ~ 'Precision',
                                  melted_df$variable == 'Recall' ~ 'Recall',
                                  melted_df$variable == 'Precision\n(common taxa)' ~ 'Precision (common taxa)',
                                  melted_df$variable == 'Recall\n(common taxa)' ~ 'Recall (common taxa)',
                                  melted_df$variable == 'AUC' ~ 'AUC',
                                  melted_df$variable == 'AUC (Common taxa)' ~ 'AUC (common taxa)',
                                  melted_df$variable == 'Relative effect error' ~ 'Relative effect error',
                                  melted_df$variable == 'Relative shrinkage error' ~ 'Relative shrinkage error',
                                  melted_df$variable == 'Effect size\nSpearman cor.' ~ 'Effect size Spearman cor.')
  melted_df$variable <- factor(melted_df$variable, levels = c('Precision', 
                                                              'Recall',
                                                              'AUC',
                                                              'Precision (common taxa)',
                                                              'Recall (common taxa)',
                                                              'AUC (common taxa)',
                                                              'Relative effect error',
                                                              'Relative shrinkage error',
                                                              'Effect size Spearman cor.'))
  melted_df$tool <- case_when(melted_df$tool == 'ALDEx2' ~ 'ALDEx2',
                              melted_df$tool == 'ANCOMBC' ~ 'ANCOM-BC2',
                              melted_df$tool == 'Maaslin2' ~ 'MaAsLin 2',
                              melted_df$tool == 'Maaslin3' ~ 'MaAsLin 3 Base',
                              melted_df$tool == 'Maaslin3itaug' ~ 'MaAsLin 3',
                              melted_df$tool == 'Maaslin3unscaled' ~ 'MaAsLin 3 Inferred\nAbundance',
                              melted_df$tool == 'Maaslin3correct' ~ 'MaAsLin 3 Sparsity\nCorrected')
  
  melted_df$tool <- paste0(melted_df$tool, " ", ifelse(melted_df$association_type == 'abundance', "Abundance", "Prevalence"))
  
  tool_vec <- c("ALDEx2", "ANCOM-BC2", "MaAsLin 2", "MaAsLin 3 Base", 
                "MaAsLin 3", "MaAsLin 3 Inferred\nAbundance", "MaAsLin 3 Sparsity\nCorrected")
  
  melted_df$tool <- factor(melted_df$tool, levels = c(paste0(tool_vec, " ", "Abundance"),
                                                      paste0(tool_vec, " ", "Prevalence")))
  
  melted_df <- melted_df[(!melted_df$variable %in% c("Relative effect error", "Relative shrinkage error")) |
    (melted_df$value < 2 & melted_df$variable %in% c("Relative effect error", "Relative shrinkage error")),] # Remove obviously wrong coefficients
  
  melted_df <- melted_df[!is.na(melted_df$variable),]
  
  # In-text
  aggregate(value ~ tool, data = melted_df[melted_df$variable == 'Precision' & melted_df$nPerSubject == 20,], FUN = mean)
  aggregate(value ~ tool, data = melted_df[melted_df$variable == 'Recall' & melted_df$nPerSubject == 20,], FUN = mean)
  aggregate(value ~ tool + nPerSubject, data = melted_df[melted_df$variable == 'AUC',], FUN = mean)
  aggregate(value ~ tool + nPerSubject, data = melted_df[melted_df$variable == 'Effect size Spearman cor.',], FUN = mean)
  aggregate(value ~ tool, data = melted_df[melted_df$variable == 'Relative shrinkage error' & melted_df$nPerSubject == 20,], FUN = mean)
  
  plot_out <- ggplot(melted_df, aes(x = get(param_name), y = value, fill = tool)) + 
    geom_boxplot(position = position_dodge(preserve = "single")) + 
    facet_wrap(~variable, scales = 'free') + 
    theme_bw() + 
    xlab("Samples per Subject") + 
    ylab('') + 
    theme(text=element_text(size=21),
          legend.position = 'bottom') + 
    labs(fill = 'Model') + 
    scale_fill_manual(values=c("#104E8B", "#458B00", "#EEAD0E", "#EE7600", "#8B1A1A", "#68228B",  "#9852CB",
                               "#88A7C5", "#A2C580", "#F7D687", "#F7BB80", "#C58D8D", "#B491C5",  "#B852CB"),
                      breaks=c(paste0(tool_vec, " ", "Abundance"),
                               paste0(tool_vec, " ", "Prevalence")))
  ggsave(paste0(figures_folder, 'fig_5.png'),
         plot = plot_out, width = 12, height = 12)
}

figure_6 <- function() {
  metadata_type <- "MVAref"
  param_name <- "effectSize"
  
  general_results_df <- data.frame(matrix(nrow = 0, ncol = 0))
  maaslin3_results_df <- data.frame(matrix(nrow = 0, ncol = 0))
  metrics <- c("Precision", "Recall", "Precision\n(common taxa)", "Recall\n(common taxa)", 
               "Relative effect error", "Relative shrinkage error", "Effect size\nSpearman cor.", "AUC", "AUC (Common taxa)", 
               "Issue proportion")
  for (param_list in param_list_final[c(2, 24, 26, 28, 30)]) { #c(2, 35, 38, 41, 44)
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
                       sep='_')
    options("scipen"=5)
    
    inputDirectory <- file.path(workingDirectory, 'Input', 'unscaled', generator)
    this_params_folder <- file.path(inputDirectory, inputString)
    outputDirectory <- file.path(workingDirectory, 'Output', 'unscaled', generator)
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
      
      truth <- prepare_truth(truth, generator)
      for (tool in tools) {
        possible_error <- tryCatch({
          associations <- read.csv(paste0(this_output_folder, '_', tool, "/associations_", i, ".tsv"), sep = "\t")
        }, error = function(err) {
          err
        })
        
        if(inherits(possible_error, "error")) next
        
        new_row <- c(unweighted_precision_recall(truth, prepare_associations_general(associations, tool, generator), abundance, metadata),
                     weighted_precision_recall(truth, prepare_associations_general(associations, tool, generator), abundance, metadata),
                     effect_size_error(truth, prepare_associations_abundance(associations, tool, generator)),
                     effect_size_correlation(truth, prepare_associations_abundance(associations, tool, generator)),
                     pval_auc(truth, prepare_associations_general(associations, tool, generator), abundance, metadata),
                     weighted_pval_auc(truth, prepare_associations_general(associations, tool, generator), abundance, metadata),
                     issue_prop(associations, tool))
        new_row <- lapply(new_row, function(x) {x})
        names(new_row) <- metrics
        
        new_row[['tool']] <- tool
        new_row[['iter']] <- i
        new_row <- c(new_row, param_list)
        general_results_df <- plyr::rbind.fill(general_results_df, data.frame(new_row, check.names = F))
        
        if (grepl('Maaslin3', tool) & generator == 'SD2') {
          tmp_truth <- truth[truth$associations == 'abundance',]
          tmp_associations <- prepare_associations_maaslin3(associations, tool)
          tmp_associations <- tmp_associations[tmp_associations$association == 'abundance',]
          new_row <- c(unweighted_precision_recall_maaslin3(tmp_truth, tmp_associations, abundance, metadata),
                       weighted_precision_recall_maaslin3(tmp_truth, tmp_associations, abundance, metadata),
                       effect_size_error_maaslin3(tmp_truth, tmp_associations),
                       effect_size_correlation_maaslin3(tmp_truth, tmp_associations),
                       pval_auc_maaslin3(tmp_truth, tmp_associations, abundance, metadata),
                       weighted_pval_auc_maaslin3(tmp_truth, tmp_associations, abundance, metadata),
                       issue_prop(associations[associations$associations == 'abundance',], tool))
          new_row <- lapply(new_row, function(x) {x})
          names(new_row) <- metrics
          
          new_row[['tool']] <- tool
          new_row[['iter']] <- i
          new_row[['association_type']] <- 'abundance'
          new_row <- c(new_row, param_list)
          maaslin3_results_df <- plyr::rbind.fill(maaslin3_results_df, data.frame(new_row, check.names = F))
          
          tmp_truth <- truth[truth$associations == 'prevalence',]
          tmp_associations <- prepare_associations_maaslin3(associations, tool)
          tmp_associations <- tmp_associations[tmp_associations$association == 'prevalence',]
          new_row <- c(unweighted_precision_recall_maaslin3(tmp_truth, tmp_associations, abundance, metadata),
                       weighted_precision_recall_maaslin3(tmp_truth, tmp_associations, abundance, metadata),
                       effect_size_error_maaslin3(tmp_truth, tmp_associations),
                       effect_size_correlation_maaslin3(tmp_truth, tmp_associations),
                       pval_auc_maaslin3(tmp_truth, tmp_associations, abundance, metadata),
                       weighted_pval_auc_maaslin3(tmp_truth, tmp_associations, abundance, metadata),
                       issue_prop(associations[associations$associations == 'prevalence',], tool))
          new_row <- lapply(new_row, function(x) {x})
          names(new_row) <- metrics
          
          new_row[['tool']] <- tool
          new_row[['iter']] <- i
          new_row[['association_type']] <- 'prevalence'
          new_row <- c(new_row, param_list)
          maaslin3_results_df <- plyr::rbind.fill(maaslin3_results_df, data.frame(new_row, check.names = F))
        }
      }
    }
  }
  
  names_to_iter <- names(param_list)[names(param_list) != 'metadataType']
  figures_folder <- paste0(workingDirectory, '/Figures/thesis_figures/')
  dir.create(figures_folder, recursive = T, showWarnings = F)
  
  default_vals <- first_params
  names(default_vals) <- names_to_iter
  default_vals <- default_vals[-which(names(default_vals) == param_name)]
  
  results_subset <- maaslin3_results_df[apply(maaslin3_results_df[,names(default_vals)], 1, function(x) {all(x == default_vals)}),]
  results_subset <- results_subset[results_subset$metadataType == metadata_type,]
  results_subset <- results_subset[,c(metrics, param_name, 'tool', 'iter', 'association_type')]
  
  melted_df <- melt(results_subset, id.vars = c('tool', 'iter', 'association_type', param_name))
  melted_df[,param_name] <- factor(melted_df[,param_name], levels = as.character(sort(unique(as.numeric(melted_df[,param_name])))))
  melted_df <- melted_df[melted_df$variable != "Issue proportion",]
  melted_df$variable <- case_when(melted_df$variable == 'Precision' ~ 'Precision',
                                  melted_df$variable == 'Recall' ~ 'Recall',
                                  melted_df$variable == 'Precision\n(common taxa)' ~ 'Precision (common taxa)',
                                  melted_df$variable == 'Recall\n(common taxa)' ~ 'Recall (common taxa)',
                                  melted_df$variable == 'AUC' ~ 'AUC',
                                  melted_df$variable == 'AUC (Common taxa)' ~ 'AUC (common taxa)',
                                  melted_df$variable == 'Relative effect error' ~ 'Relative effect error',
                                  melted_df$variable == 'Relative shrinkage error' ~ 'Relative shrinkage error',
                                  melted_df$variable == 'Effect size\nSpearman cor.' ~ 'Effect size Spearman cor.')
  melted_df$variable <- factor(melted_df$variable, levels = c('Precision', 
                                                              'Recall',
                                                              'AUC',
                                                              'Precision (common taxa)',
                                                              'Recall (common taxa)',
                                                              'AUC (common taxa)',
                                                              'Relative effect error',
                                                              'Relative shrinkage error',
                                                              'Effect size Spearman cor.'))
  melted_df$tool <- case_when(melted_df$tool == 'ALDEx2' ~ 'ALDEx2',
                              melted_df$tool == 'ANCOMBC' ~ 'ANCOM-BC2',
                              melted_df$tool == 'Maaslin2' ~ 'MaAsLin 2',
                              melted_df$tool == 'Maaslin3' ~ 'MaAsLin 3 Base',
                              melted_df$tool == 'Maaslin3itaug' ~ 'MaAsLin 3',
                              melted_df$tool == 'Maaslin3unscaled' ~ 'MaAsLin 3 Inferred\nAbundance',
                              melted_df$tool == 'Maaslin3correct' ~ 'MaAsLin 3 Sparsity\nCorrected')
  
  melted_df$tool <- paste0(melted_df$tool, " ", ifelse(melted_df$association_type == 'abundance', "Abundance", "Prevalence"))
  
  tool_vec <- c("ALDEx2", "ANCOM-BC2", "MaAsLin 2", "MaAsLin 3 Base", 
                "MaAsLin 3", "MaAsLin 3 Inferred\nAbundance", "MaAsLin 3 Sparsity\nCorrected")
  
  melted_df$tool <- factor(melted_df$tool, levels = c(paste0(tool_vec, " ", "Abundance"),
                                                      paste0(tool_vec, " ", "Prevalence")))
  
  # In-text
  aggregate(value ~ tool, data = melted_df[melted_df$variable == 'Precision' & melted_df$effectSize == 2,], FUN = mean)
  aggregate(value ~ tool, data = melted_df[melted_df$variable == 'Recall' & melted_df$effectSize == 2,], FUN = mean)
  aggregate(value ~ tool, data = melted_df[melted_df$variable == 'Relative shrinkage error' & melted_df$effectSize == 2,], FUN = mean)
  aggregate(value ~ tool, data = melted_df[melted_df$variable == 'Precision' & melted_df$effectSize == 5,], FUN = mean)
  aggregate(value ~ tool, data = melted_df[melted_df$variable == 'Recall' & melted_df$effectSize == 5,], FUN = mean)
  aggregate(value ~ tool, data = melted_df[melted_df$variable == 'Relative effect error' & melted_df$effectSize == 5,], FUN = mean)
  aggregate(value ~ tool, data = melted_df[melted_df$variable == 'Precision' & melted_df$effectSize == 10,], FUN = mean)
  aggregate(value ~ tool, data = melted_df[melted_df$variable == 'Relative shrinkage error' & melted_df$effectSize == 10,], FUN = mean)
  
  plot_out <- ggplot(melted_df, aes(x = get(param_name), y = value, fill = tool)) + 
    geom_boxplot(position = position_dodge(preserve = "single")) + 
    facet_wrap(~variable, scales = 'free') + 
    theme_bw() + 
    xlab("Effect Size") + 
    ylab('') + 
    theme(text=element_text(size=21),
          legend.position = 'bottom') + 
    labs(fill = 'Model') + 
    scale_fill_manual(values=c("#104E8B", "#458B00", "#EEAD0E", "#EE7600", "#8B1A1A", "#68228B",  "#9852CB",
                               "#88A7C5", "#A2C580", "#F7D687", "#F7BB80", "#C58D8D", "#B491C5",  "#B852CB"),
                      breaks=c(paste0(tool_vec, " ", "Abundance"),
                               paste0(tool_vec, " ", "Prevalence")))
  ggsave(paste0(figures_folder, 'fig_6.png'),
         plot = plot_out, width = 12, height = 12)
}

figure_1()
figure_2()
figure_3()
figure_4()
figure_5()
figure_6()

parameters <- "~/Documents/GitHub/maaslin3_benchmark/general_evaluations/data_generation/ANCOM_BC_generator.txt"

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

figure_7 <- function() {
  metadata_type <- "soil"
  param_name <- "nSubjects"
  
  general_results_df <- data.frame(matrix(nrow = 0, ncol = 0))
  maaslin3_results_df <- data.frame(matrix(nrow = 0, ncol = 0))
  metrics <- c("Precision", "Recall", "Precision\n(common taxa)", "Recall\n(common taxa)", 
               "Relative effect error", "Relative shrinkage error", "Effect size\nSpearman cor.", "AUC", "AUC (Common taxa)", 
               "Issue proportion")
  for (param_list in param_list_final[c(2, 8, 10, 12, 14, 16)]) {
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
                       sep='_')
    options("scipen"=5)
    
    inputDirectory <- file.path(workingDirectory, 'Input', 'general_evaluations', generator)
    this_params_folder <- file.path(inputDirectory, inputString)
    outputDirectory <- file.path(workingDirectory, 'Output', 'general_evaluations', generator)
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
      
      truth <- prepare_truth(truth, generator)
      for (tool in tools) {
        possible_error <- tryCatch({
          associations <- read.csv(paste0(this_output_folder, '_', tool, "/associations_", i, ".tsv"), sep = "\t")
        }, error = function(err) {
          err
        })
        
        if(inherits(possible_error, "error")) next
        
        new_row <- c(unweighted_precision_recall(truth, prepare_associations_general(associations, tool, generator), abundance, metadata),
                     weighted_precision_recall(truth, prepare_associations_general(associations, tool, generator), abundance, metadata),
                     effect_size_error(truth, prepare_associations_abundance(associations, tool, generator)),
                     effect_size_correlation(truth[truth$associations == 'abundance',], prepare_associations_abundance(associations, tool, generator)),
                     pval_auc(truth, prepare_associations_general(associations, tool, generator), abundance, metadata),
                     weighted_pval_auc(truth, prepare_associations_general(associations, tool, generator), abundance, metadata),
                     issue_prop(associations, tool))
        new_row <- lapply(new_row, function(x) {x})
        names(new_row) <- metrics
        
        new_row[['tool']] <- tool
        new_row[['iter']] <- i
        new_row <- c(new_row, param_list)
        general_results_df <- plyr::rbind.fill(general_results_df, data.frame(new_row, check.names = F))
      }
    }
  }
  
  names_to_iter <- names(param_list)[names(param_list) != 'metadataType']
  figures_folder <- paste0(workingDirectory, '/Figures/thesis_figures/')
  dir.create(figures_folder, recursive = T, showWarnings = F)
  
  default_vals <- first_params
  names(default_vals) <- names_to_iter
  default_vals <- default_vals[-which(names(default_vals) == param_name)]
  
  results_subset <- general_results_df[apply(general_results_df[,names(default_vals)], 1, function(x) {all(x == default_vals)}),]
  results_subset <- results_subset[results_subset$metadataType == metadata_type,]
  results_subset <- results_subset[,c(metrics, param_name, 'tool', 'iter')]
  
  melted_df <- melt(results_subset, id.vars = c('tool', 'iter', param_name))
  melted_df[,param_name] <- factor(melted_df[,param_name], levels = as.character(sort(unique(as.numeric(melted_df[,param_name])))))
  melted_df <- melted_df[melted_df$variable != "Issue proportion",]
  melted_df$variable <- case_when(melted_df$variable == 'Precision' ~ 'Precision',
                                  melted_df$variable == 'Recall' ~ 'Recall',
                                  melted_df$variable == 'Precision\n(common taxa)' ~ 'Precision (common taxa)',
                                  melted_df$variable == 'Recall\n(common taxa)' ~ 'Recall (common taxa)',
                                  melted_df$variable == 'AUC' ~ 'AUC',
                                  melted_df$variable == 'AUC (Common taxa)' ~ 'AUC (common taxa)',
                                  melted_df$variable == 'Relative effect error' ~ 'Relative effect error',
                                  melted_df$variable == 'Relative shrinkage error' ~ 'Relative shrinkage error',
                                  melted_df$variable == 'Effect size\nSpearman cor.' ~ 'Effect size Spearman cor.')
  melted_df$variable <- factor(melted_df$variable, levels = c('Precision', 
                                                              'Recall',
                                                              'AUC',
                                                              'Precision (common taxa)',
                                                              'Recall (common taxa)',
                                                              'AUC (common taxa)',
                                                              'Relative effect error',
                                                              'Relative shrinkage error',
                                                              'Effect size Spearman cor.'))
  melted_df$tool <- case_when(melted_df$tool == 'ALDEx2' ~ 'ALDEx2',
                              melted_df$tool == 'ANCOMBC' ~ 'ANCOM-BC2',
                              melted_df$tool == 'Maaslin2' ~ 'MaAsLin 2',
                              melted_df$tool == 'Maaslin3' ~ 'MaAsLin 3 Base',
                              melted_df$tool == 'Maaslin3itaug' ~ 'MaAsLin 3',
                              melted_df$tool == 'Maaslin3unscaled' ~ 'MaAsLin 3 Inferred\nAbundance')
  melted_df$tool <- factor(melted_df$tool, levels = c("ALDEx2", "ANCOM-BC2", "MaAsLin 2", "MaAsLin 3 Base", "MaAsLin 3", "MaAsLin 3 Inferred\nAbundance"))
  
  # In-text
  aggregate(value ~ tool + nSubjects, data = melted_df[melted_df$variable == 'Precision',], FUN = mean)
  aggregate(value ~ tool + nSubjects, data = melted_df[melted_df$variable == 'Precision (common taxa)',], FUN = mean)
  aggregate(value ~ tool + nSubjects, data = melted_df[melted_df$variable == 'Recall',], FUN = mean)
  aggregate(value ~ tool + nSubjects, data = melted_df[melted_df$variable == 'Recall (common taxa)',], FUN = mean)
  aggregate(value ~ tool + nSubjects, data = melted_df[melted_df$variable == 'AUC',], FUN = mean)
  aggregate(value ~ tool + nSubjects, data = melted_df[melted_df$variable == 'Relative shrinkage error',], FUN = mean)
  aggregate(value ~ tool + nSubjects, data = melted_df[melted_df$variable == 'Effect size Spearman cor.',], FUN = mean)
  
  plot_out <- ggplot(melted_df, aes(x = get(param_name), y = value, fill = tool)) + 
    geom_boxplot(position = position_dodge(preserve = "single")) + 
    facet_wrap(~variable, scales = 'free') + 
    theme_bw() + 
    xlab("Number of Subjects") + 
    ylab('') + 
    theme(text=element_text(size=21),
          legend.position = 'bottom') + 
    labs(fill = 'Model') + 
    scale_fill_manual(values=c("#104E8B", "#458B00", "#EEAD0E", "#EE7600", "#8B1A1A", "#68228B"),
                      breaks=c("ALDEx2", "ANCOM-BC2", "MaAsLin 2", "MaAsLin 3 Base", "MaAsLin 3", "MaAsLin 3 Inferred\nAbundance"))
  ggsave(paste0(figures_folder, 'fig_7.png'),
         plot = plot_out, width = 12, height = 12)
}
figure_7()

parameters <- "~/Documents/GitHub/maaslin3_benchmark/general_evaluations/data_generation/SimSeq.txt"

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

figure_8 <- function() {
  metadata_type <- "all"
  param_name <- "nSubjects"
  
  general_results_df <- data.frame(matrix(nrow = 0, ncol = 0))
  maaslin3_results_df <- data.frame(matrix(nrow = 0, ncol = 0))
  metrics <- c("Precision", "Recall", "Precision\n(common taxa)", "Recall\n(common taxa)", 
               "Relative effect error", "Relative shrinkage error", "Effect size\nSpearman cor.", "AUC", "AUC (Common taxa)", 
               "Issue proportion")
  for (param_list in param_list_final[c(1:5, 11:20)]) {
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
                       sep='_')
    options("scipen"=5)
    
    inputDirectory <- file.path(workingDirectory, 'Input', 'general_evaluations', generator)
    this_params_folder <- file.path(inputDirectory, inputString)
    outputDirectory <- file.path(workingDirectory, 'Output', 'general_evaluations', generator)
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
      
      truth <- prepare_truth(truth, generator)
      for (tool in tools) {
        possible_error <- tryCatch({
          associations <- read.csv(paste0(this_output_folder, '_', tool, "/associations_", i, ".tsv"), sep = "\t")
        }, error = function(err) {
          err
        })
        
        if(inherits(possible_error, "error")) next
        
        new_row <- c(unweighted_precision_recall(truth, prepare_associations_general(associations, tool, generator), abundance, metadata),
                     weighted_precision_recall(truth, prepare_associations_general(associations, tool, generator), abundance, metadata),
                     effect_size_error(truth, prepare_associations_abundance(associations, tool, generator)),
                     effect_size_correlation(truth, prepare_associations_abundance(associations, tool, generator)),
                     pval_auc(truth, prepare_associations_general(associations, tool, generator), abundance, metadata),
                     weighted_pval_auc(truth, prepare_associations_general(associations, tool, generator), abundance, metadata),
                     issue_prop(associations, tool))
        new_row <- lapply(new_row, function(x) {x})
        names(new_row) <- metrics
        
        new_row[['tool']] <- tool
        new_row[['iter']] <- i
        new_row <- c(new_row, param_list)
        general_results_df <- plyr::rbind.fill(general_results_df, data.frame(new_row, check.names = F))
      }
    }
  }
  
  names_to_iter <- names(param_list)[names(param_list) != 'metadataType']
  figures_folder <- paste0(workingDirectory, '/Figures/thesis_figures/')
  dir.create(figures_folder, recursive = T, showWarnings = F)
  
  default_vals <- first_params
  names(default_vals) <- names_to_iter
  default_vals <- default_vals[-which(names(default_vals) == param_name)]
  
  results_subset <- general_results_df[apply(general_results_df[,names(default_vals)], 1, function(x) {all(x == default_vals)}),]
  #results_subset <- results_subset[results_subset$metadataType == metadata_type,]
  results_subset <- results_subset[,c(metrics, param_name, 'tool', 'iter')]
  
  melted_df <- melt(results_subset, id.vars = c('tool', 'iter', param_name))
  melted_df[,param_name] <- factor(melted_df[,param_name], levels = as.character(sort(unique(as.numeric(melted_df[,param_name])))))
  melted_df <- melted_df[melted_df$variable != "Issue proportion",]
  melted_df$variable <- case_when(melted_df$variable == 'Precision' ~ 'Precision',
                                  melted_df$variable == 'Recall' ~ 'Recall',
                                  melted_df$variable == 'Precision\n(common taxa)' ~ 'Precision (common taxa)',
                                  melted_df$variable == 'Recall\n(common taxa)' ~ 'Recall (common taxa)',
                                  melted_df$variable == 'AUC' ~ 'AUC',
                                  melted_df$variable == 'AUC (Common taxa)' ~ 'AUC (common taxa)',
                                  melted_df$variable == 'Relative effect error' ~ 'Relative effect error',
                                  melted_df$variable == 'Relative shrinkage error' ~ 'Relative shrinkage error',
                                  melted_df$variable == 'Effect size\nSpearman cor.' ~ 'Effect size Spearman cor.')
  melted_df$variable <- factor(melted_df$variable, levels = c('Precision', 
                                                              'Recall',
                                                              'AUC',
                                                              'Precision (common taxa)',
                                                              'Recall (common taxa)',
                                                              'AUC (common taxa)',
                                                              'Relative effect error',
                                                              'Relative shrinkage error',
                                                              'Effect size Spearman cor.'))
  melted_df$tool <- case_when(melted_df$tool == 'ALDEx2' ~ 'ALDEx2',
                              melted_df$tool == 'ANCOMBC' ~ 'ANCOM-BC2',
                              melted_df$tool == 'Maaslin2' ~ 'MaAsLin 2',
                              melted_df$tool == 'Maaslin3' ~ 'MaAsLin 3 Base',
                              melted_df$tool == 'Maaslin3itaug' ~ 'MaAsLin 3',
                              melted_df$tool == 'Maaslin3unscaled' ~ 'MaAsLin 3 Inferred\nAbundance')
  melted_df$tool <- factor(melted_df$tool, levels = c("ALDEx2", "ANCOM-BC2", "MaAsLin 2", "MaAsLin 3 Base", "MaAsLin 3", "MaAsLin 3 Inferred\nAbundance"))
  
  melted_df <- melted_df[melted_df$variable %in% c("Precision", "Recall", "AUC", "Precision (common taxa)", "Recall (common taxa)", "AUC (common taxa)"),]
  
  # In-text
  aggregate(value ~ tool + nSubjects, data = melted_df[melted_df$variable == 'Precision',], FUN = mean)
  aggregate(value ~ tool + nSubjects, data = melted_df[melted_df$variable == 'Recall',], FUN = mean)
  aggregate(value ~ tool + nSubjects, data = melted_df[melted_df$variable == 'AUC',], FUN = mean)
  
  plot_out <- ggplot(melted_df, aes(x = get(param_name), y = value, fill = tool)) + 
    geom_boxplot(position = position_dodge(preserve = "single")) + 
    facet_wrap(~variable, scales = 'free') + 
    theme_bw() + 
    xlab("Number of Subjects") + 
    ylab('') + 
    theme(text=element_text(size=21),
          legend.position = 'bottom') + 
    labs(fill = 'Model') + 
    scale_fill_manual(values=c("#104E8B", "#458B00", "#EEAD0E", "#EE7600", "#8B1A1A", "#68228B"),
                      breaks=c("ALDEx2", "ANCOM-BC2", "MaAsLin 2", "MaAsLin 3 Base", "MaAsLin 3", "MaAsLin 3 Inferred\nAbundance"))
  ggsave(paste0(figures_folder, 'fig_8.png'),
         plot = plot_out, width = 12, height = 8.5)
}
figure_8()

parameters <- "~/Documents/GitHub/maaslin3_benchmark/omps/data_generation/SD2.txt"

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

figure_9 <- function() {
  metadata_type <- 'MVBomp'
  param_name <- 'nSubjects'
  
  general_results_df <- data.frame(matrix(nrow = 0, ncol = 0))
  maaslin3_results_df <- data.frame(matrix(nrow = 0, ncol = 0))
  metrics <- c("Precision", "Recall", "Precision\n(common taxa)", "Recall\n(common taxa)", 
               "AUC", "AUC (Common taxa)", "Relative effect error", "Relative shrinkage error", "Effect size\nSpearman cor.",
               "Issue proportion")
  for (param_list in param_list_final[c(2, 4, 6, 8, 10)]) {
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
                       param_list[['nOmps']],
                       param_list[['nLevels']],
                       sep='_')
    options("scipen"=5)
    
    inputDirectory <- file.path(workingDirectory, 'Input', 'omp_evaluations', generator)
    this_params_folder <- file.path(inputDirectory, inputString)
    outputDirectory <- file.path(workingDirectory, 'Output', 'omp_evaluations', generator)
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
      truth <- truth[truth$effect_size != 0,]
      for (tool in tools) {
        possible_error <- tryCatch({
          associations <- read.csv(paste0(this_output_folder, '_', tool, "/associations_", i, ".tsv"), sep = "\t")
        }, error = function(err) {
          err
        })
        
        if(inherits(possible_error, "error")) next
        
        if (tool %in% c("Maaslin3ItAugLinear", "Maaslin3ItAugGomp")) {
          associations_out <- prepare_associations_maaslin3(associations, tool)
          associations_out <- associations_out[associations_out$metadata %in% truth$metadata & !is.na(associations_out$signif),]
          
          new_row <- c(unweighted_precision_recall_maaslin3(truth, associations_out, abundance, metadata),
                       weighted_precision_recall_maaslin3(truth, associations_out, abundance, metadata),
                       pval_auc_maaslin3(truth, associations_out, abundance, metadata),
                       weighted_pval_auc_maaslin3(truth, associations_out, abundance, metadata),
                       c(NA, NA),
                       NA,
                       issue_prop(associations[associations$metadata %in% truth$metadata,], tool))
        }
        
        if (tool %in% c("Maaslin3ItAugOneHot", "Maaslin3ItAugOmp")) {
          truth_out <- truth
          truth_out$metadata <- truth_out$org_metadata
          truth_out$org_metadata <- NULL
          associations_out <- prepare_associations_maaslin3(associations, tool, remove_possible_error = T)
          associations_out <- associations_out[associations_out$metadata %in% truth_out$metadata & !is.na(associations_out$signif),]
          truth_out <- truth_out[truth_out$effect_size != 0,]
          
          new_row <- c(unweighted_precision_recall_maaslin3(truth_out, associations_out, abundance, metadata),
                       weighted_precision_recall_maaslin3(truth_out, associations_out, abundance, metadata),
                       pval_auc_maaslin3(truth_out, associations_out, abundance, metadata),
                       weighted_pval_auc_maaslin3(truth_out, associations_out, abundance, metadata),
                       effect_size_error_maaslin3(truth_out[abs(truth_out$effect_size) > 0.5,], associations_out),
                       effect_size_correlation_maaslin3(truth_out, associations_out),
                       issue_prop(associations[associations$metadata %in% truth$org_metadata,], tool))
        }
        
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
  figures_folder <- paste0(workingDirectory, '/Figures/thesis_figures/')
  dir.create(figures_folder, recursive = T, showWarnings = F)
  
  default_vals <- first_params
  names(default_vals) <- names_to_iter
  default_vals <- default_vals[-which(names(default_vals) == param_name)]
  
  results_subset <- maaslin3_results_df[apply(maaslin3_results_df[,names(default_vals)], 1, function(x) {all(x == default_vals)}),]
  results_subset <- results_subset[results_subset$metadataType == metadata_type,]
  results_subset <- results_subset[,c(metrics, param_name, 'tool', 'iter')]
  
  melted_df <- melt(results_subset, id.vars = c('tool', 'iter', param_name))
  melted_df[,param_name] <- factor(melted_df[,param_name], levels = as.character(sort(unique(as.numeric(melted_df[,param_name])))))
  melted_df <- melted_df[melted_df$variable != "Issue proportion",]
  melted_df$variable <- case_when(melted_df$variable == 'Precision' ~ 'Precision',
                                  melted_df$variable == 'Recall' ~ 'Recall',
                                  melted_df$variable == 'Precision\n(common taxa)' ~ 'Precision (common taxa)',
                                  melted_df$variable == 'Recall\n(common taxa)' ~ 'Recall (common taxa)',
                                  melted_df$variable == 'AUC' ~ 'AUC',
                                  melted_df$variable == 'AUC (Common taxa)' ~ 'AUC (common taxa)',
                                  melted_df$variable == 'Relative effect error' ~ 'Relative effect error',
                                  melted_df$variable == 'Relative shrinkage error' ~ 'Relative shrinkage error',
                                  melted_df$variable == 'Effect size\nSpearman cor.' ~ 'Effect size Spearman cor.')
  melted_df$variable <- factor(melted_df$variable, levels = c('Precision', 
                                                              'Recall',
                                                              'AUC',
                                                              'Precision (common taxa)',
                                                              'Recall (common taxa)',
                                                              'AUC (common taxa)',
                                                              'Relative effect error',
                                                              'Relative shrinkage error',
                                                              'Effect size Spearman cor.'))
  melted_df$tool <- case_when(melted_df$tool == 'Maaslin3ItAugGomp' ~ 'MaAsLin 3 GOMP',
                              melted_df$tool == 'Maaslin3ItAugLinear' ~ 'MaAsLin 3 Linear',
                              melted_df$tool == 'Maaslin3ItAugOmp' ~ 'MaAsLin 3 OMP',
                              melted_df$tool == 'Maaslin3ItAugOneHot' ~ 'MaAsLin 3 One-hot')
  
  # In-text
  aggregate(value ~ tool + nSubjects, data = melted_df[melted_df$variable == 'Precision',], FUN = mean)
  aggregate(value ~ tool + nSubjects, data = melted_df[melted_df$variable == 'Relative shrinkage error',], FUN = mean)
  aggregate(value ~ tool + nSubjects, data = melted_df[melted_df$variable == 'Effect size Spearman cor.',], FUN = mean)
  
  plot_out <- ggplot(melted_df, aes(x = get(param_name), y = value, fill = tool)) + 
    geom_boxplot(position = position_dodge(preserve = "single")) + 
    facet_wrap(~variable, scales = 'free') + 
    theme_bw() + 
    xlab('Number of Subjects') + 
    ylab('') + 
    theme(text=element_text(size=21),
          legend.position = 'bottom') + 
    labs(fill = 'Model') + 
    scale_fill_brewer(palette="Paired")
  ggsave(paste0(figures_folder, 'fig_9.png'),
         plot = plot_out, width = 12, height = 10)
}

figure_10 <- function() {
  metadata_type <- 'MVBomp'
  param_name <- 'nPerSubject'
  
  general_results_df <- data.frame(matrix(nrow = 0, ncol = 0))
  maaslin3_results_df <- data.frame(matrix(nrow = 0, ncol = 0))
  metrics <- c("Precision", "Recall", "Precision\n(common taxa)", "Recall\n(common taxa)", 
               "AUC", "AUC (Common taxa)", "Relative effect error", "Relative shrinkage error", "Effect size\nSpearman cor.",
               "Issue proportion")
  for (param_list in param_list_final[c(2, 32, 34)]) {
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
                       param_list[['nOmps']],
                       param_list[['nLevels']],
                       sep='_')
    options("scipen"=5)
    
    inputDirectory <- file.path(workingDirectory, 'Input', 'omp_evaluations', generator)
    this_params_folder <- file.path(inputDirectory, inputString)
    outputDirectory <- file.path(workingDirectory, 'Output', 'omp_evaluations', generator)
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
      truth <- truth[truth$effect_size != 0,]
      for (tool in tools) {
        possible_error <- tryCatch({
          associations <- read.csv(paste0(this_output_folder, '_', tool, "/associations_", i, ".tsv"), sep = "\t")
        }, error = function(err) {
          err
        })
        
        if(inherits(possible_error, "error")) next
        
        if (tool %in% c("Maaslin3ItAugLinear", "Maaslin3ItAugGomp")) {
          associations_out <- prepare_associations_maaslin3(associations, tool)
          associations_out <- associations_out[associations_out$metadata %in% truth$metadata & !is.na(associations_out$signif),]
          
          new_row <- c(unweighted_precision_recall_maaslin3(truth, associations_out, abundance, metadata),
                       weighted_precision_recall_maaslin3(truth, associations_out, abundance, metadata),
                       pval_auc_maaslin3(truth, associations_out, abundance, metadata),
                       weighted_pval_auc_maaslin3(truth, associations_out, abundance, metadata),
                       c(NA, NA),
                       NA,
                       issue_prop(associations[associations$metadata %in% truth$metadata,], tool))
        }
        
        if (tool %in% c("Maaslin3ItAugOneHot", "Maaslin3ItAugOmp")) {
          truth_out <- truth
          truth_out$metadata <- truth_out$org_metadata
          truth_out$org_metadata <- NULL
          associations_out <- prepare_associations_maaslin3(associations, tool, remove_possible_error = T)
          associations_out <- associations_out[associations_out$metadata %in% truth_out$metadata & !is.na(associations_out$signif),]
          truth_out <- truth_out[truth_out$effect_size != 0,]
          
          new_row <- c(unweighted_precision_recall_maaslin3(truth_out, associations_out, abundance, metadata),
                       weighted_precision_recall_maaslin3(truth_out, associations_out, abundance, metadata),
                       pval_auc_maaslin3(truth_out, associations_out, abundance, metadata),
                       weighted_pval_auc_maaslin3(truth_out, associations_out, abundance, metadata),
                       effect_size_error_maaslin3(truth_out[abs(truth_out$effect_size) > 0.5,], associations_out),
                       effect_size_correlation_maaslin3(truth_out, associations_out),
                       issue_prop(associations[associations$metadata %in% truth$org_metadata,], tool))
        }
        
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
  figures_folder <- paste0(workingDirectory, '/Figures/thesis_figures/')
  dir.create(figures_folder, recursive = T, showWarnings = F)
  
  default_vals <- first_params
  names(default_vals) <- names_to_iter
  default_vals <- default_vals[-which(names(default_vals) == param_name)]
  
  results_subset <- maaslin3_results_df[apply(maaslin3_results_df[,names(default_vals)], 1, function(x) {all(x == default_vals)}),]
  results_subset <- results_subset[results_subset$metadataType == metadata_type,]
  results_subset <- results_subset[,c(metrics, param_name, 'tool', 'iter')]
  
  melted_df <- melt(results_subset, id.vars = c('tool', 'iter', param_name))
  melted_df[,param_name] <- factor(melted_df[,param_name], levels = as.character(sort(unique(as.numeric(melted_df[,param_name])))))
  melted_df <- melted_df[melted_df$variable != "Issue proportion",]
  melted_df$variable <- case_when(melted_df$variable == 'Precision' ~ 'Precision',
                                  melted_df$variable == 'Recall' ~ 'Recall',
                                  melted_df$variable == 'Precision\n(common taxa)' ~ 'Precision (common taxa)',
                                  melted_df$variable == 'Recall\n(common taxa)' ~ 'Recall (common taxa)',
                                  melted_df$variable == 'AUC' ~ 'AUC',
                                  melted_df$variable == 'AUC (Common taxa)' ~ 'AUC (common taxa)',
                                  melted_df$variable == 'Relative effect error' ~ 'Relative effect error',
                                  melted_df$variable == 'Relative shrinkage error' ~ 'Relative shrinkage error',
                                  melted_df$variable == 'Effect size\nSpearman cor.' ~ 'Effect size Spearman cor.')
  melted_df$variable <- factor(melted_df$variable, levels = c('Precision', 
                                                              'Recall',
                                                              'AUC',
                                                              'Precision (common taxa)',
                                                              'Recall (common taxa)',
                                                              'AUC (common taxa)',
                                                              'Relative effect error',
                                                              'Relative shrinkage error',
                                                              'Effect size Spearman cor.'))
  melted_df$tool <- case_when(melted_df$tool == 'Maaslin3ItAugGomp' ~ 'MaAsLin 3 GOMP',
                              melted_df$tool == 'Maaslin3ItAugLinear' ~ 'MaAsLin 3 Linear',
                              melted_df$tool == 'Maaslin3ItAugOmp' ~ 'MaAsLin 3 OMP',
                              melted_df$tool == 'Maaslin3ItAugOneHot' ~ 'MaAsLin 3 One-hot')
  
  # In-text
  aggregate(value ~ tool + nPerSubject, data = melted_df[melted_df$variable == 'Precision',], FUN = mean)
  aggregate(value ~ tool + nPerSubject, data = melted_df[melted_df$variable == 'Precision (common taxa)',], FUN = mean)
  
  plot_out <- ggplot(melted_df, aes(x = get(param_name), y = value, fill = tool)) + 
    geom_boxplot(position = position_dodge(preserve = "single")) + 
    facet_wrap(~variable, scales = 'free') + 
    theme_bw() + 
    xlab('Samples per Subject') + 
    ylab('') + 
    theme(text=element_text(size=21),
          legend.position = 'bottom') + 
    labs(fill = 'Model') + 
    scale_fill_brewer(palette="Paired")
  ggsave(paste0(figures_folder, 'fig_10.png'),
         plot = plot_out, width = 12, height = 10)
}

figure_11 <- function() {
  metadata_type <- 'MVBomp'
  param_name <- 'nLevels'
  
  general_results_df <- data.frame(matrix(nrow = 0, ncol = 0))
  maaslin3_results_df <- data.frame(matrix(nrow = 0, ncol = 0))
  metrics <- c("Precision", "Recall", "Precision\n(common taxa)", "Recall\n(common taxa)", 
               "AUC", "AUC (Common taxa)", "Relative effect error", "Relative shrinkage error", "Effect size\nSpearman cor.",
               "Issue proportion")
  for (param_list in param_list_final[c(2, 28, 30)]) {
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
                       param_list[['nOmps']],
                       param_list[['nLevels']],
                       sep='_')
    options("scipen"=5)
    
    inputDirectory <- file.path(workingDirectory, 'Input', 'omp_evaluations', generator)
    this_params_folder <- file.path(inputDirectory, inputString)
    outputDirectory <- file.path(workingDirectory, 'Output', 'omp_evaluations', generator)
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
      truth <- truth[truth$effect_size != 0,]
      for (tool in tools) {
        possible_error <- tryCatch({
          associations <- read.csv(paste0(this_output_folder, '_', tool, "/associations_", i, ".tsv"), sep = "\t")
        }, error = function(err) {
          err
        })
        
        if(inherits(possible_error, "error")) next
        
        if (tool %in% c("Maaslin3ItAugLinear", "Maaslin3ItAugGomp")) {
          associations_out <- prepare_associations_maaslin3(associations, tool)
          associations_out <- associations_out[associations_out$metadata %in% truth$metadata & !is.na(associations_out$signif),]
          
          new_row <- c(unweighted_precision_recall_maaslin3(truth, associations_out, abundance, metadata),
                       weighted_precision_recall_maaslin3(truth, associations_out, abundance, metadata),
                       pval_auc_maaslin3(truth, associations_out, abundance, metadata),
                       weighted_pval_auc_maaslin3(truth, associations_out, abundance, metadata),
                       c(NA, NA),
                       NA,
                       issue_prop(associations[associations$metadata %in% truth$metadata,], tool))
        }
        
        if (tool %in% c("Maaslin3ItAugOneHot", "Maaslin3ItAugOmp")) {
          truth_out <- truth
          truth_out$metadata <- truth_out$org_metadata
          truth_out$org_metadata <- NULL
          associations_out <- prepare_associations_maaslin3(associations, tool, remove_possible_error = T)
          associations_out <- associations_out[associations_out$metadata %in% truth_out$metadata & !is.na(associations_out$signif),]
          truth_out <- truth_out[truth_out$effect_size != 0,]
          
          new_row <- c(unweighted_precision_recall_maaslin3(truth_out, associations_out, abundance, metadata),
                       weighted_precision_recall_maaslin3(truth_out, associations_out, abundance, metadata),
                       pval_auc_maaslin3(truth_out, associations_out, abundance, metadata),
                       weighted_pval_auc_maaslin3(truth_out, associations_out, abundance, metadata),
                       effect_size_error_maaslin3(truth_out[abs(truth_out$effect_size) > 0.5,], associations_out),
                       effect_size_correlation_maaslin3(truth_out, associations_out),
                       issue_prop(associations[associations$metadata %in% truth$org_metadata,], tool))
        }
        
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
  figures_folder <- paste0(workingDirectory, '/Figures/thesis_figures/')
  dir.create(figures_folder, recursive = T, showWarnings = F)
  
  default_vals <- first_params
  names(default_vals) <- names_to_iter
  default_vals <- default_vals[-which(names(default_vals) == param_name)]
  
  results_subset <- maaslin3_results_df[apply(maaslin3_results_df[,names(default_vals)], 1, function(x) {all(x == default_vals)}),]
  results_subset <- results_subset[results_subset$metadataType == metadata_type,]
  results_subset <- results_subset[,c(metrics, param_name, 'tool', 'iter')]
  
  melted_df <- melt(results_subset, id.vars = c('tool', 'iter', param_name))
  melted_df[,param_name] <- factor(melted_df[,param_name], levels = as.character(sort(unique(as.numeric(melted_df[,param_name])))))
  melted_df <- melted_df[melted_df$variable != "Issue proportion",]
  melted_df$variable <- case_when(melted_df$variable == 'Precision' ~ 'Precision',
                                  melted_df$variable == 'Recall' ~ 'Recall',
                                  melted_df$variable == 'Precision\n(common taxa)' ~ 'Precision (common taxa)',
                                  melted_df$variable == 'Recall\n(common taxa)' ~ 'Recall (common taxa)',
                                  melted_df$variable == 'AUC' ~ 'AUC',
                                  melted_df$variable == 'AUC (Common taxa)' ~ 'AUC (common taxa)',
                                  melted_df$variable == 'Relative effect error' ~ 'Relative effect error',
                                  melted_df$variable == 'Relative shrinkage error' ~ 'Relative shrinkage error',
                                  melted_df$variable == 'Effect size\nSpearman cor.' ~ 'Effect size Spearman cor.')
  melted_df$variable <- factor(melted_df$variable, levels = c('Precision', 
                                                              'Recall',
                                                              'AUC',
                                                              'Precision (common taxa)',
                                                              'Recall (common taxa)',
                                                              'AUC (common taxa)',
                                                              'Relative effect error',
                                                              'Relative shrinkage error',
                                                              'Effect size Spearman cor.'))
  melted_df$tool <- case_when(melted_df$tool == 'Maaslin3ItAugGomp' ~ 'MaAsLin 3 GOMP',
                              melted_df$tool == 'Maaslin3ItAugLinear' ~ 'MaAsLin 3 Linear',
                              melted_df$tool == 'Maaslin3ItAugOmp' ~ 'MaAsLin 3 OMP',
                              melted_df$tool == 'Maaslin3ItAugOneHot' ~ 'MaAsLin 3 One-hot')
  
  plot_out <- ggplot(melted_df, aes(x = get(param_name), y = value, fill = tool)) + 
    geom_boxplot(position = position_dodge(preserve = "single")) + 
    facet_wrap(~variable, scales = 'free') + 
    theme_bw() + 
    xlab('Number of Levels') + 
    ylab('') + 
    theme(text=element_text(size=21),
          legend.position = 'bottom') + 
    labs(fill = 'Model') + 
    scale_fill_brewer(palette="Paired")
  ggsave(paste0(figures_folder, 'fig_11.png'),
         plot = plot_out, width = 12, height = 10)
}

parameters <- "~/Documents/GitHub/maaslin3_benchmark/groups/data_generation/SD2.txt"

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

figure_12 <- function() {
  metadata_type <- 'MVBgroup'
  param_name <- 'nSubjects'
  
  general_results_df <- data.frame(matrix(nrow = 0, ncol = 0))
  maaslin3_results_df <- data.frame(matrix(nrow = 0, ncol = 0))
  metrics <- c("Precision", "Recall", "Precision\n(common taxa)", "Recall\n(common taxa)", 
               "AUC", "AUC (Common taxa)", "Issue proportion")
  for (param_list in param_list_final[c(2, 4, 6, 8, 10)]) {
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
  figures_folder <- paste0(workingDirectory, '/Figures/thesis_figures/')
  dir.create(figures_folder, recursive = T, showWarnings = F)
  
  default_vals <- first_params
  names(default_vals) <- names_to_iter
  default_vals <- default_vals[-which(names(default_vals) == param_name)]
  
  results_subset <- maaslin3_results_df[apply(maaslin3_results_df[,names(default_vals)], 1, function(x) {all(x == default_vals)}),]
  results_subset <- results_subset[results_subset$metadataType == metadata_type,]
  results_subset <- results_subset[,c(metrics, param_name, 'tool', 'iter')]
  
  melted_df <- melt(results_subset, id.vars = c('tool', 'iter', param_name))
  melted_df[,param_name] <- factor(melted_df[,param_name], levels = as.character(sort(unique(as.numeric(melted_df[,param_name])))))
  melted_df <- melted_df[melted_df$variable != "Issue proportion",]
  melted_df$variable <- case_when(melted_df$variable == 'Precision' ~ 'Precision',
                                  melted_df$variable == 'Recall' ~ 'Recall',
                                  melted_df$variable == 'Precision\n(common taxa)' ~ 'Precision (common taxa)',
                                  melted_df$variable == 'Recall\n(common taxa)' ~ 'Recall (common taxa)',
                                  melted_df$variable == 'AUC' ~ 'AUC',
                                  melted_df$variable == 'AUC (Common taxa)' ~ 'AUC (common taxa)',
                                  melted_df$variable == 'Relative effect error' ~ 'Relative effect error',
                                  melted_df$variable == 'Relative shrinkage error' ~ 'Relative shrinkage error',
                                  melted_df$variable == 'Effect size\nSpearman cor.' ~ 'Effect size Spearman cor.')
  melted_df$variable <- factor(melted_df$variable, levels = c('Precision', 
                                                              'Recall',
                                                              'AUC',
                                                              'Precision (common taxa)',
                                                              'Recall (common taxa)',
                                                              'AUC (common taxa)',
                                                              'Relative effect error',
                                                              'Relative shrinkage error',
                                                              'Effect size Spearman cor.'))
  melted_df$tool <- case_when(melted_df$tool == 'Maaslin3ItAugGroup' ~ 'MaAsLin 3 Group',
                              melted_df$tool == 'Maaslin3ItAugOneHot' ~ 'MaAsLin 3 One-hot')
  
  # In-text
  aggregate(value ~ tool + nSubjects, data = melted_df[melted_df$variable == 'Precision',], FUN = mean)
  aggregate(value ~ tool + nSubjects, data = melted_df[melted_df$variable == 'Recall',], FUN = mean)
  aggregate(value ~ tool + nSubjects, data = melted_df[melted_df$variable == 'AUC',], FUN = mean)
  
  plot_out <- ggplot(melted_df, aes(x = get(param_name), y = value, fill = tool)) + 
    geom_boxplot(position = position_dodge(preserve = "single")) + 
    facet_wrap(~variable, scales = 'free') + 
    theme_bw() + 
    xlab('Number of Subjects') + 
    ylab('') + 
    theme(text=element_text(size=21),
          legend.position = 'bottom') + 
    labs(fill = 'Model') +
    scale_fill_brewer(palette="Pastel2")

  ggsave(paste0(figures_folder, 'fig_12.png'),
         plot = plot_out, width = 12, height = 7)
}

#########################
# Read depth evaluation #
#########################

figure_13 <- function() {
  metadata_type <- "MVAref"
  param_name <- "readDepth"
  
  general_results_df <- data.frame(matrix(nrow = 0, ncol = 0))
  maaslin3_results_df <- data.frame(matrix(nrow = 0, ncol = 0))
  metrics <- c("Precision", "Recall", "Precision\n(common taxa)", "Recall\n(common taxa)", 
               "Relative effect error", "Relative shrinkage error", "Effect size\nSpearman cor.", "AUC", "AUC (Common taxa)", 
               "Issue proportion")
  for (param_list in param_list_final[c(2, 4, 6, 8)]) { #c(2, 14, 17, 20, 23, 26)
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
                       sep='_')
    options("scipen"=5)
    
    inputDirectory <- file.path(workingDirectory, 'Input', 'unscaled', generator)
    this_params_folder <- file.path(inputDirectory, inputString)
    outputDirectory <- file.path(workingDirectory, 'Output', 'unscaled', generator)
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
      
      truth <- prepare_truth(truth, generator)
      for (tool in tools) {
        possible_error <- tryCatch({
          associations <- read.csv(paste0(this_output_folder, '_', tool, "/associations_", i, ".tsv"), sep = "\t")
        }, error = function(err) {
          err
        })
        
        if(inherits(possible_error, "error")) next
        
        new_row <- c(unweighted_precision_recall(truth, prepare_associations_general(associations, tool, generator), abundance, metadata),
                     weighted_precision_recall(truth, prepare_associations_general(associations, tool, generator), abundance, metadata),
                     effect_size_error(truth, prepare_associations_abundance(associations, tool, generator)),
                     effect_size_correlation(truth, prepare_associations_abundance(associations, tool, generator)),
                     pval_auc(truth, prepare_associations_general(associations, tool, generator), abundance, metadata),
                     weighted_pval_auc(truth, prepare_associations_general(associations, tool, generator), abundance, metadata),
                     issue_prop(associations, tool))
        new_row <- lapply(new_row, function(x) {x})
        names(new_row) <- metrics
        
        new_row[['tool']] <- tool
        new_row[['iter']] <- i
        new_row <- c(new_row, param_list)
        general_results_df <- plyr::rbind.fill(general_results_df, data.frame(new_row, check.names = F))
        
        if (grepl('Maaslin3', tool) & generator == 'SD2') {
          tmp_truth <- truth[truth$associations == 'abundance',]
          tmp_associations <- prepare_associations_maaslin3(associations, tool)
          tmp_associations <- tmp_associations[tmp_associations$association == 'abundance',]
          new_row <- c(unweighted_precision_recall_maaslin3(tmp_truth, tmp_associations, abundance, metadata),
                       weighted_precision_recall_maaslin3(tmp_truth, tmp_associations, abundance, metadata),
                       effect_size_error_maaslin3(tmp_truth, tmp_associations),
                       effect_size_correlation_maaslin3(tmp_truth, tmp_associations),
                       pval_auc_maaslin3(tmp_truth, tmp_associations, abundance, metadata),
                       weighted_pval_auc_maaslin3(tmp_truth, tmp_associations, abundance, metadata),
                       issue_prop(associations[associations$associations == 'abundance',], tool))
          new_row <- lapply(new_row, function(x) {x})
          names(new_row) <- metrics
          
          new_row[['tool']] <- tool
          new_row[['iter']] <- i
          new_row[['association_type']] <- 'abundance'
          new_row <- c(new_row, param_list)
          maaslin3_results_df <- plyr::rbind.fill(maaslin3_results_df, data.frame(new_row, check.names = F))
          
          tmp_truth <- truth[truth$associations == 'prevalence',]
          tmp_associations <- prepare_associations_maaslin3(associations, tool)
          tmp_associations <- tmp_associations[tmp_associations$association == 'prevalence',]
          new_row <- c(unweighted_precision_recall_maaslin3(tmp_truth, tmp_associations, abundance, metadata),
                       weighted_precision_recall_maaslin3(tmp_truth, tmp_associations, abundance, metadata),
                       effect_size_error_maaslin3(tmp_truth, tmp_associations),
                       effect_size_correlation_maaslin3(tmp_truth, tmp_associations),
                       pval_auc_maaslin3(tmp_truth, tmp_associations, abundance, metadata),
                       weighted_pval_auc_maaslin3(tmp_truth, tmp_associations, abundance, metadata),
                       issue_prop(associations[associations$associations == 'prevalence',], tool))
          new_row <- lapply(new_row, function(x) {x})
          names(new_row) <- metrics
          
          new_row[['tool']] <- tool
          new_row[['iter']] <- i
          new_row[['association_type']] <- 'prevalence'
          new_row <- c(new_row, param_list)
          maaslin3_results_df <- plyr::rbind.fill(maaslin3_results_df, data.frame(new_row, check.names = F))
        }
      }
    }
  }
  
  names_to_iter <- names(param_list)[names(param_list) != 'metadataType']
  figures_folder <- paste0(workingDirectory, '/Figures/thesis_figures/')
  dir.create(figures_folder, recursive = T, showWarnings = F)
  
  default_vals <- first_params
  names(default_vals) <- names_to_iter
  default_vals <- default_vals[-which(names(default_vals) == param_name)]
  
  results_subset <- maaslin3_results_df[apply(maaslin3_results_df[,names(default_vals)], 1, function(x) {all(x == default_vals)}),]
  results_subset <- results_subset[results_subset$metadataType == metadata_type,]
  results_subset <- results_subset[,c(metrics, param_name, 'tool', 'iter', 'association_type')]
  
  melted_df <- melt(results_subset, id.vars = c('tool', 'iter', 'association_type', param_name))
  melted_df[,param_name] <- factor(melted_df[,param_name], levels = as.character(sort(unique(as.numeric(melted_df[,param_name])))))
  melted_df <- melted_df[melted_df$variable != "Issue proportion",]
  melted_df$variable <- case_when(melted_df$variable == 'Precision' ~ 'Precision',
                                  melted_df$variable == 'Recall' ~ 'Recall',
                                  melted_df$variable == 'Precision\n(common taxa)' ~ 'Precision (common taxa)',
                                  melted_df$variable == 'Recall\n(common taxa)' ~ 'Recall (common taxa)',
                                  melted_df$variable == 'AUC' ~ 'AUC',
                                  melted_df$variable == 'AUC (Common taxa)' ~ 'AUC (common taxa)',
                                  melted_df$variable == 'Relative effect error' ~ 'Relative effect error',
                                  melted_df$variable == 'Relative shrinkage error' ~ 'Relative shrinkage error',
                                  melted_df$variable == 'Effect size\nSpearman cor.' ~ 'Effect size Spearman cor.')
  melted_df$variable <- factor(melted_df$variable, levels = c('Precision', 
                                                              'Recall',
                                                              'AUC',
                                                              'Precision (common taxa)',
                                                              'Recall (common taxa)',
                                                              'AUC (common taxa)',
                                                              'Relative effect error',
                                                              'Relative shrinkage error',
                                                              'Effect size Spearman cor.'))
  melted_df$tool <- case_when(melted_df$tool == 'ALDEx2' ~ 'ALDEx2',
                              melted_df$tool == 'ANCOMBC' ~ 'ANCOM-BC2',
                              melted_df$tool == 'Maaslin2' ~ 'MaAsLin 2',
                              melted_df$tool == 'Maaslin3' ~ 'MaAsLin 3 Base',
                              melted_df$tool == 'Maaslin3itaug' ~ 'MaAsLin 3',
                              melted_df$tool == 'Maaslin3unscaled' ~ 'MaAsLin 3 Inferred\nAbundance')
  
  melted_df$tool <- paste0(melted_df$tool, " ", ifelse(melted_df$association_type == 'abundance', "Abundance", "Prevalence"))
  
  tool_vec <- c("ALDEx2", "ANCOM-BC2", "MaAsLin 2", "MaAsLin 3 Base", 
                "MaAsLin 3", "MaAsLin 3 Inferred\nAbundance")
  
  melted_df$tool <- factor(melted_df$tool, levels = c(paste0(tool_vec, " ", "Abundance"),
                                                      paste0(tool_vec, " ", "Prevalence")))
  
  plot_out <- ggplot(melted_df, aes(x = get(param_name), y = value, fill = tool)) + 
    geom_boxplot(position = position_dodge(preserve = "single")) + 
    facet_wrap(~variable, scales = 'free') + 
    theme_bw() + 
    xlab("Number of Subjects") + 
    ylab('') + 
    theme(text=element_text(size=21),
          legend.position = 'bottom') + 
    labs(fill = 'Model') + 
    scale_fill_manual(values=c("#104E8B", "#458B00", "#EEAD0E", "#EE7600", "#8B1A1A", "#68228B",
                               "#88A7C5", "#A2C580", "#F7D687", "#F7BB80", "#C58D8D", "#B491C5"),
                      breaks=c(paste0(tool_vec, " ", "Abundance"),
                               paste0(tool_vec, " ", "Prevalence")))
}














