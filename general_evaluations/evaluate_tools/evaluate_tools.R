remove(list = ls())
package_vec = c("reshape2", "ggplot2", "optparse", 
                "parallel", "stringi", "doParallel", "plyr", "tidyr")
invisible(suppressPackageStartupMessages(lapply(package_vec, require, character.only = TRUE)))
source('library/run_evaluation_helpers.R')

ancom_gen_figure <- function() {
  parameters <- "general_evaluations/data_generation/ANCOM_BC_generator.txt"
  
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
  
  metadata_type <- "soil"
  param_name <- "nSubjects"
  
  general_results_df <- data.frame(matrix(nrow = 0, ncol = 0))
  metrics <- c("Precision", "Recall", "Precision\n(common taxa)", "Recall\n(common taxa)", 
               "Relative effect error", "Relative shrinkage error", "Effect size\nSpearman cor.", "AUC", "AUC (Common taxa)", 
               "Issue proportion")
  
  # Just look at variations in sample number
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
    
    inputDirectory <- file.path('Input', 'general_evaluations', generator)
    this_params_folder <- file.path(inputDirectory, inputString)
    outputDirectory <- file.path('Output', 'general_evaluations', generator)
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
  figures_folder <- 'Figures/paper_figures/'
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
  melted_df$tool <- case_when(melted_df$tool == 'ALDEx2' ~ 'ALDEx2',
                              melted_df$tool == 'ANCOMBC' ~ 'ANCOM-BC2',
                              melted_df$tool == 'Maaslin2' ~ 'MaAsLin 2',
                              melted_df$tool == 'Maaslin3CompAdjust' ~ 'MaAsLin 3')
  melted_df$tool <- factor(melted_df$tool, levels = c("ALDEx2", "ANCOM-BC2", "MaAsLin 2", "MaAsLin 3"))
  melted_df <- melted_df[!is.na(melted_df$tool),]
  
  # In-text numbers
  print(melted_df %>%
      dplyr::filter(variable == 'Precision') %>%
      dplyr::group_by(tool, nSubjects) %>%
      dplyr::summarize(mean(value, na.rm=T)), n = 100)
  
  print(melted_df %>%
            dplyr::filter(variable == 'Effect size Spearman cor.') %>%
            dplyr::group_by(tool, nSubjects) %>%
            dplyr::summarize(mean(value, na.rm=T)), n = 100)
  
  melted_df %>%
            dplyr::filter(variable == 'Relative shrinkage error', nSubjects == "1000") %>%
            dplyr::group_by(tool, nSubjects) %>%
            dplyr::summarize(mean(value, na.rm=T))
  
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
    scale_fill_manual(values=c("#104E8B", "#458B00", "#EEAD0E", "#68228B"),
                      breaks=c("ALDEx2", "ANCOM-BC2", "MaAsLin 2", "MaAsLin 3")) + 
    geom_point(aes(x = get(param_name), y = 0), alpha = 0)
  ggsave(paste0(figures_folder, 'ancom_gen_figure.png'),
         plot = plot_out, width = 8, height = 8)
}
ancom_gen_figure()

# simseq_figure <- function() {
#   parameters <- "general_evaluations/data_generation/SimSeq.txt"
#   
#   filename_pieces <- unlist(strsplit(parameters, '/'))
#   generator <- gsub('.txt', '', filename_pieces[length(filename_pieces)])
#   
#   nIterations = 100
#   
#   parameter_text <- readLines(parameters)
#   
#   param_list = list()
#   for (parameter_line in parameter_text) {
#     parts = trimws(unlist(strsplit(parameter_line, ':')))
#     vals = unlist(strsplit(parts[2], ' '))
#     param_list[[parts[1]]] = vals
#   }
#   
#   metadata_types = param_list[['metadataType']]
#   param_list[['metadataType']] <- NULL
#   
#   params_for_files <- list()
#   for (param in names(param_list)) {
#     first_params <- c()
#     for (param_tmp in names(param_list)) {
#       first_params <- c(first_params, param_list[[param_tmp]][1])
#     }
#     for (j in seq_along(param_list[[param]])) {
#       params_for_files <- c(params_for_files, list(first_params))
#       names(params_for_files[[length(params_for_files)]]) <- names(param_list)
#       params_for_files[[length(params_for_files)]][param] <- param_list[[param]][j]
#     }
#   }
#   
#   params_for_files <- unique(params_for_files)
#   
#   param_list_final = list()
#   for (param_file in params_for_files) {
#     for (metadata_type in metadata_types) {
#       tmp_param_file <- param_file
#       tmp_param_file['metadataType'] <- metadata_type
#       if (metadata_type == 'binary') {
#         tmp_param_file['nMetadata'] <- '1'
#         tmp_param_file['nPerSubject'] <- '1'
#       }
#       param_list_final <- c(param_list_final, list(tmp_param_file))
#     }
#   }
#   
#   metadata_type <- "all"
#   param_name <- "nSubjects"
#   
#   general_results_df <- data.frame(matrix(nrow = 0, ncol = 0))
#   metrics <- c("Precision", "Recall", "Precision\n(common taxa)", "Recall\n(common taxa)", 
#                "Relative effect error", "Relative shrinkage error", "Effect size\nSpearman cor.", "AUC", "AUC (Common taxa)", 
#                "Issue proportion")
#   for (param_list in param_list_final[c(1:5, 11:20)]) { # Select only relevant settings
#     print(param_list)
#     if (param_list[['nPerSubject']]!='1'){
#       inputSubString<-'RandomEffect'
#     } else {
#       inputSubString<-'noRandomEffect'
#     }
#     
#     options("scipen"=10)
#     inputString<-paste(inputSubString, 
#                        param_list[['metadataType']], 
#                        param_list[['nSubjects']], 
#                        param_list[['nPerSubject']], 
#                        param_list[['nMicrobes']], 
#                        param_list[['spikeMicrobes']], 
#                        param_list[['nMetadata']], 
#                        param_list[['effectSize']], 
#                        param_list[['effectPos']], 
#                        param_list[['readDepth']], 
#                        sep='_')
#     options("scipen"=5)
#     
#     inputDirectory <- file.path('Input', 'general_evaluations', generator)
#     this_params_folder <- file.path(inputDirectory, inputString)
#     outputDirectory <- file.path('Output', 'general_evaluations', generator)
#     tmp_files <- list.files(outputDirectory)
#     tmp_files <- tmp_files[grepl(paste0(inputString, "_"), tmp_files)]
#     tools <- gsub(paste0(inputString, "_"), '', tmp_files, perl = T)
#     this_output_folder <- file.path(outputDirectory, inputString)
#     
#     for (i in 1:nIterations) {
#       possible_error <- tryCatch({
#         metadata <- read.csv(paste0(this_params_folder, "/metadata_", i, ".tsv"), sep = "\t")
#         metadata <- metadata[,colnames(metadata) != 'ID']
#         abundance <- read.csv(paste0(this_params_folder, "/abundance_", i, ".tsv"), sep = "\t")
#         truth <- read.csv(paste0(this_params_folder, "/truth_", i, ".tsv"), sep = "\t")
#       }, error = function(err) {
#         err
#       })
#       if(inherits(possible_error, "error")) next
#       
#       truth <- prepare_truth(truth, generator)
#       for (tool in tools) {
#         possible_error <- tryCatch({
#           associations <- read.csv(paste0(this_output_folder, '_', tool, "/associations_", i, ".tsv"), sep = "\t")
#         }, error = function(err) {
#           err
#         })
#         
#         if(inherits(possible_error, "error")) next
#         
#         new_row <- c(unweighted_precision_recall(truth, prepare_associations_general(associations, tool, generator), abundance, metadata),
#                      weighted_precision_recall(truth, prepare_associations_general(associations, tool, generator), abundance, metadata),
#                      effect_size_error(truth, prepare_associations_abundance(associations, tool, generator)),
#                      effect_size_correlation(truth, prepare_associations_abundance(associations, tool, generator)),
#                      pval_auc(truth, prepare_associations_general(associations, tool, generator), abundance, metadata),
#                      weighted_pval_auc(truth, prepare_associations_general(associations, tool, generator), abundance, metadata),
#                      issue_prop(associations, tool))
#         new_row <- lapply(new_row, function(x) {x})
#         names(new_row) <- metrics
#         
#         new_row[['tool']] <- tool
#         new_row[['iter']] <- i
#         new_row <- c(new_row, param_list)
#         general_results_df <- plyr::rbind.fill(general_results_df, data.frame(new_row, check.names = F))
#       }
#     }
#   }
#   
#   names_to_iter <- names(param_list)[names(param_list) != 'metadataType']
#   figures_folder <- paste0('Figures/paper_figures/')
#   dir.create(figures_folder, recursive = T, showWarnings = F)
#   
#   default_vals <- first_params
#   names(default_vals) <- names_to_iter
#   default_vals <- default_vals[-which(names(default_vals) == param_name)]
#   
#   results_subset <- general_results_df[apply(general_results_df[,names(default_vals)], 1, function(x) {all(x == default_vals)}),]
#   results_subset <- results_subset[,c(metrics, param_name, 'tool', 'iter')]
#   
#   melted_df <- melt(results_subset, id.vars = c('tool', 'iter', param_name))
#   melted_df[,param_name] <- factor(melted_df[,param_name], levels = as.character(sort(unique(as.numeric(melted_df[,param_name])))))
#   melted_df <- melted_df[melted_df$variable != "Issue proportion",]
#   melted_df$variable <- case_when(melted_df$variable == 'Precision' ~ 'Precision',
#                                   melted_df$variable == 'Recall' ~ 'Recall',
#                                   melted_df$variable == 'Precision\n(common taxa)' ~ 'Precision (common taxa)',
#                                   melted_df$variable == 'Recall\n(common taxa)' ~ 'Recall (common taxa)',
#                                   melted_df$variable == 'AUC' ~ 'AUC',
#                                   melted_df$variable == 'AUC (Common taxa)' ~ 'AUC (common taxa)',
#                                   melted_df$variable == 'Relative effect error' ~ 'Relative effect error',
#                                   melted_df$variable == 'Relative shrinkage error' ~ 'Relative shrinkage error',
#                                   melted_df$variable == 'Effect size\nSpearman cor.' ~ 'Effect size Spearman cor.')
#   melted_df$variable <- factor(melted_df$variable, levels = c('Precision', 
#                                                               'Recall'))
#   melted_df <- melted_df[!is.na(melted_df$variable),]
#   melted_df$tool <- case_when(melted_df$tool == 'ALDEx2' ~ 'ALDEx2',
#                               melted_df$tool == 'ANCOMBC' ~ 'ANCOM-BC2',
#                               melted_df$tool == 'Maaslin2' ~ 'MaAsLin 2',
#                               melted_df$tool == 'Maaslin3CompAdjust' ~ 'MaAsLin 3')
#   melted_df$tool <- factor(melted_df$tool, levels = c("ALDEx2", "ANCOM-BC2", "MaAsLin 2", "MaAsLin 3"))
#   melted_df <- melted_df[!is.na(melted_df$tool),]
#   
#   plot_out <- ggplot(melted_df, aes(x = get(param_name), y = value, fill = tool)) + 
#     geom_boxplot(position = position_dodge(preserve = "single")) + 
#     facet_wrap(~variable, scales = 'free') + 
#     theme_bw() + 
#     xlab("Number of Subjects") + 
#     ylab('') + 
#     theme(text=element_text(size=21),
#           legend.position = 'bottom',
#           strip.background = element_rect(fill = "gray95")) + 
#     labs(fill = 'Model') + 
#     scale_fill_manual(values=c("#104E8B", "#458B00", "#EEAD0E", "#68228B"),
#                       breaks=c("ALDEx2", "ANCOM-BC2", "MaAsLin 2", "MaAsLin 3"))
#   ggsave(paste0(figures_folder, 'simseq_figure.png'),
#          plot = plot_out, width = 8, height = 4)
# }
# simseq_figure()









