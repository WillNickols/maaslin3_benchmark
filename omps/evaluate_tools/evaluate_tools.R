remove(list = ls())
package_vec = c("reshape2", "ggplot2", "optparse", 
                "parallel", "stringi", "doParallel", "plyr", "tidyr")
invisible(suppressPackageStartupMessages(lapply(package_vec, require, character.only = TRUE)))
source('library/run_evaluation_helpers.R')

ordered_sample_size_figure <- function() {
  parameters <- "omps/data_generation/SD2.txt"
  
  filename_pieces <- unlist(strsplit(parameters, '/'))
  generator <- gsub('.txt', '', filename_pieces[length(filename_pieces)])
  
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
  
  metadata_type <- "MVAomp"
  param_name <- "nSubjects"
  
  general_results_df <- data.frame(matrix(nrow = 0, ncol = 0))
  metrics <- c("Precision", "Recall", "Precision\n(common taxa)", "Recall\n(common taxa)", 
               "AUC", "AUC (Common taxa)", "Relative effect error", "Relative shrinkage error", "Effect size\nSpearman cor.",
               "Issue proportion")
  for (param_list in param_list_final[c(1, 3, 5, 7, 9)]) {
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
    
    inputDirectory <- file.path('Input', 'omp_evaluations', generator)
    this_params_folder <- file.path(inputDirectory, inputString)
    outputDirectory <- file.path('Output', 'omp_evaluations', generator)
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
        
        truth_out <- truth
        truth_out$metadata <- truth_out$org_metadata
        truth_out$org_metadata <- NULL
        associations_out <- prepare_associations_maaslin3(associations, tool, remove_possible_error = T)
        associations_out <- associations_out[associations_out$metadata %in% truth_out$metadata & !is.na(associations_out$signif),]

        new_row <- c(unweighted_precision_recall_maaslin3(truth_out, associations_out, abundance, metadata),
                     weighted_precision_recall_maaslin3(truth_out, associations_out, abundance, metadata),
                     pval_auc_maaslin3(truth_out, associations_out, abundance, metadata),
                     weighted_pval_auc_maaslin3(truth_out, associations_out, abundance, metadata),
                     effect_size_error_maaslin3(truth_out[abs(truth_out$effect_size) > 0.5,], associations_out),
                     effect_size_correlation_maaslin3(truth_out, associations_out),
                     issue_prop(associations[associations$metadata %in% truth$org_metadata,], tool))
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
  figures_folder <- paste0('Figures/paper_figures/')
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
                                                              'Relative shrinkage error',
                                                              'Effect size Spearman cor.'))
  melted_df <- melted_df[!is.na(melted_df$variable),]
  melted_df$tool <- case_when(melted_df$tool == 'Maaslin3' ~ 'MaAsLin 3 (No adjustment)',
                              melted_df$tool == 'Maaslin3CompAdjust' ~ 'MaAsLin 3')
  melted_df$tool <- factor(melted_df$tool, levels = c('MaAsLin 3 (No adjustment)', 'MaAsLin 3'))
  melted_df <- melted_df[!is.na(melted_df$tool),]
  
  plot_out <- ggplot(melted_df, aes(x = get(param_name), y = value, fill = tool)) + 
    geom_boxplot(position = position_dodge(preserve = "single")) + 
    facet_wrap(~variable, scales = 'free') + 
    theme_bw() + 
    xlab("Number of Subjects") + 
    ylab('') + 
    theme(text=element_text(size=21),
          legend.position = 'bottom',
          strip.background = element_rect(fill = "gray95")) + 
    labs(fill = 'Model') + 
    scale_fill_manual(values=c("#EE7600", "#68228B"),
                      breaks=c('MaAsLin 3 (No adjustment)', 'MaAsLin 3'))
  ggsave(paste0(figures_folder, 'ordered_sample_size_figure.png'),
         plot = plot_out, width = 8, height = 8)
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
  melted_df$tool <- case_when(melted_df$tool == 'Maaslin3' ~ 'MaAsLin 3 Relative',
                              melted_df$tool == 'Maaslin3CompAdjust' ~ 'MaAsLin 3 Median\nAdjusted')
  melted_df$tool <- factor(melted_df$tool, levels = c('MaAsLin 3 Relative', 'MaAsLin 3 Median\nAdjusted'))
  
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
