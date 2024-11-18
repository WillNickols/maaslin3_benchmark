remove(list = ls())
package_vec = c("reshape2", "ggplot2", "optparse", 
                "parallel", "stringi", "doParallel", "plyr", "tidyr")
invisible(suppressPackageStartupMessages(lapply(package_vec, require, character.only = TRUE)))
source('library/run_evaluation_helpers.R')

SD2_read_depth_figure_maaslin3 <- function() {
  parameters <- "deep_sequencing/data_generation/SD2.txt"
  
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
  
  metadata_type <- "MVAref"
  param_name <- "readDepth"
  
  general_results_df <- data.frame(matrix(nrow = 0, ncol = 0))
  maaslin3_results_df <- data.frame(matrix(nrow = 0, ncol = 0))
  metrics <- c("Precision", "Recall", "Precision\n(common taxa)", "Recall\n(common taxa)", 
               "Relative effect error", "Relative shrinkage error", "Effect size\nSpearman cor.", "AUC", "AUC (Common taxa)", 
               "Issue proportion")
  for (param_list in param_list_final[c(2, 4, 6, 8)]) {
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
    
    inputDirectory <- file.path('Input', 'deep_sequencing', generator)
    this_params_folder <- file.path(inputDirectory, inputString)
    outputDirectory <- file.path('Output', 'deep_sequencing', generator)
    tmp_files <- list.files(outputDirectory)
    tmp_files <- tmp_files[grepl(paste0(inputString, "_"), tmp_files)]
    tools <- gsub(paste0(inputString, "_"), '', tmp_files, perl = T)
    this_output_folder <- file.path(outputDirectory, inputString)
    
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
    
    inputDirectory <- file.path('Input', 'deep_sequencing', generator)
    this_params_folder <- file.path(inputDirectory, inputString)
    outputDirectory <- file.path('Output', 'deep_sequencing', generator)
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
  figures_folder <- paste0('Figures/paper_figures/')
  dir.create(figures_folder, recursive = T, showWarnings = F)
  
  default_vals <- first_params
  names(default_vals) <- names_to_iter
  default_vals <- default_vals[-which(names(default_vals) == param_name)]
  
  results_subset <- maaslin3_results_df[apply(maaslin3_results_df[,names(default_vals)], 1, function(x) {all(x == default_vals)}),]
  results_subset <- results_subset[results_subset$metadataType == metadata_type,]
  results_subset <- results_subset[,c(metrics, param_name, 'tool', 'iter', 'association_type')]
  
  melted_df <- melt(results_subset, id.vars = c('tool', 'iter', 'association_type', param_name))
  melted_df[,param_name] <- format(as.numeric(melted_df[,param_name]), scientific = TRUE)
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
                              melted_df$tool == 'Maaslin3CompAdjust' ~ 'MaAsLin 3',
                              melted_df$tool == 'Maaslin3' ~ 'MaAsLin 3\n(No adjustment)',
                              melted_df$tool == 'Maaslin3unscaled' ~ 'MaAsLin 3\n(Spike-in)')
  melted_df$tool <- paste0(melted_df$tool, " ", ifelse(melted_df$association_type == 'abundance', "Abundance", "Prevalence"))
  
  tool_vec <- c("MaAsLin 3\n(No adjustment)", "MaAsLin 3", "MaAsLin 3\n(Spike-in)")
  
  melted_df$tool <- factor(melted_df$tool, levels = c(paste0(tool_vec, " ", "Abundance"),
                                                      paste0(tool_vec[2], " ", "Prevalence")))
  melted_df <- melted_df[!is.na(melted_df$tool),]
  
  # In-text numbers
  melted_df %>%
      dplyr::filter(variable == 'Precision') %>%
      dplyr::group_by(tool, readDepth) %>%
      dplyr::summarize(mean(value, na.rm=T))
  
  melted_df %>%
      dplyr::filter(variable == 'Recall') %>%
      dplyr::group_by(tool, readDepth) %>%
      dplyr::summarize(mean(value, na.rm=T))
  
  melted_df %>%
      dplyr::filter(variable == 'Relative shrinkage error') %>%
      dplyr::group_by(tool, readDepth) %>%
      dplyr::summarize(mean(value, na.rm=T))
  
  melted_df %>%
      dplyr::filter(variable == 'Effect size Spearman cor.') %>%
      dplyr::group_by(tool, readDepth) %>%
      dplyr::summarize(mean(value, na.rm=T))
  
  plot_out <- ggplot(melted_df, aes(x = get(param_name), y = value, fill = tool)) + 
    geom_boxplot(position = position_dodge(preserve = "single")) + 
    facet_wrap(~variable, scales = 'free', ncol = 4) + 
    theme_bw() + 
    xlab("Read depth") + 
    ylab('') + 
    theme(text=element_text(size=21),
          legend.position = 'bottom',
          strip.background = element_rect(fill = "gray95")) + 
    labs(fill = 'Model') + 
    scale_fill_manual(values=c("#EE7600", "#68228B", "#8B1A1A",
                               "#F7BB80", "#996515", "#C58D8D"),
                      breaks=c(paste0(tool_vec, " ", "Abundance"),
                               paste0(tool_vec, " ", "Prevalence"))) + 
      geom_point(aes(x = get(param_name), y = 0), alpha = 0)
  ggsave(paste0(figures_folder, 'SD2_read_depth_figure_maaslin3.png'),
         plot = plot_out, width = 16, height = 4)
  
  melted_df$readDepth <- case_when(melted_df$readDepth == '5e+04' ~ '50K',
                                   melted_df$readDepth == '5e+07' ~ '50M',
                                   TRUE ~ NA)
  melted_df$readDepth <- factor(melted_df$readDepth, levels = c('50K', '50M'))
  melted_df <- melted_df[!grepl("No adjustment", melted_df$tool),]
  
  melted_df$variable <- case_when(melted_df$variable == 'Precision' ~ 'Precision',
                                  melted_df$variable == 'Recall' ~ 'Recall',
                                  melted_df$variable == 'Relative shrinkage error' ~ 'Relative shrinkage\nerror',
                                  melted_df$variable == 'Effect size Spearman cor.' ~ 'Effect size\nSpearman cor.')
  melted_df$variable <- factor(melted_df$variable, levels = c('Precision', 
                                                              'Recall',
                                                              'Relative shrinkage\nerror',
                                                              'Effect size\nSpearman cor.'))
  
  tool_vec <- c("MaAsLin 3 (No adjustment)", "MaAsLin 3", "MaAsLin 3 (Spike-in)")
  melted_df$tool <- ifelse(melted_df$tool == 'MaAsLin 3\n(Spike-in) Abundance', 
                           'MaAsLin 3 (Spike-in) Abundance',
                           as.character(melted_df$tool))
  
  melted_df$tool <- factor(melted_df$tool, levels = c(paste0(tool_vec, " ", "Abundance"),
                                                      paste0(tool_vec[2], " ", "Prevalence")))
  
  plot_out <- ggplot(melted_df[melted_df$readDepth %in% c("50K", "50M"),], aes(x = get(param_name), y = value, fill = tool)) + 
      geom_boxplot(position = position_dodge(preserve = "single")) + 
      facet_wrap(~variable, scales = 'free', ncol = 4) + 
      theme_bw() + 
      xlab("Read depth") + 
      ylab('') + 
      theme(text=element_text(size=21),
            legend.position = 'bottom',
            strip.background = element_rect(fill = "gray95")) + 
      labs(fill = 'Model') + 
      scale_fill_manual(values=c("#EE7600", "#68228B", "#8B1A1A",
                                 "#F7BB80", "#996515", "#C58D8D"),
                        breaks=c(paste0(tool_vec, " ", "Abundance"),
                                 paste0(tool_vec, " ", "Prevalence"))) + 
      geom_point(aes(x = get(param_name), y = 0), alpha = 0)
  ggsave(paste0(figures_folder, 'SD2_read_depth_subfigure_maaslin3.png'),
         plot = plot_out, width = 12, height = 4)
}
SD2_read_depth_figure_maaslin3()

SD2_read_depth_figure <- function() {
  parameters <- "deep_sequencing/data_generation/SD2.txt"
  
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
  
  metadata_type <- "MVAref"
  param_name <- "readDepth"
  
  general_results_df <- data.frame(matrix(nrow = 0, ncol = 0))
  maaslin3_results_df <- data.frame(matrix(nrow = 0, ncol = 0))
  metrics <- c("Precision", "Recall", "Precision\n(common taxa)", "Recall\n(common taxa)", 
               "Relative effect error", "Relative shrinkage error", "Effect size\nSpearman cor.", "AUC", "AUC (Common taxa)", 
               "Issue proportion")
  for (param_list in param_list_final[c(2, 4, 6, 8)]) {
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
    
    inputDirectory <- file.path('Input', 'deep_sequencing', generator)
    this_params_folder <- file.path(inputDirectory, inputString)
    outputDirectory <- file.path('Output', 'deep_sequencing', generator)
    tmp_files <- list.files(outputDirectory)
    tmp_files <- tmp_files[grepl(paste0(inputString, "_"), tmp_files)]
    tools <- gsub(paste0(inputString, "_"), '', tmp_files, perl = T)
    this_output_folder <- file.path(outputDirectory, inputString)
    
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
    
    inputDirectory <- file.path('Input', 'deep_sequencing', generator)
    this_params_folder <- file.path(inputDirectory, inputString)
    outputDirectory <- file.path('Output', 'deep_sequencing', generator)
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
  figures_folder <- paste0('Figures/paper_figures/')
  dir.create(figures_folder, recursive = T, showWarnings = F)
  
  default_vals <- first_params
  names(default_vals) <- names_to_iter
  default_vals <- default_vals[-which(names(default_vals) == param_name)]
  
  results_subset <- general_results_df[apply(general_results_df[,names(default_vals)], 1, function(x) {all(x == default_vals)}),]
  results_subset <- results_subset[results_subset$metadataType == metadata_type,]
  results_subset <- results_subset[,c(metrics, param_name, 'tool', 'iter')]
  
  melted_df <- melt(results_subset, id.vars = c('tool', 'iter', param_name))
  melted_df[,param_name] <- format(as.numeric(melted_df[,param_name]), scientific = TRUE)
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
  
  # Remove few extreme ANCOM values that make plot worse
  melted_df <- melted_df[!is.na(melted_df$value) & melted_df$value <= 1,]
  
  plot_out <- ggplot(melted_df, aes(x = get(param_name), y = value, fill = tool)) + 
    geom_boxplot(position = position_dodge(preserve = "single")) + 
    facet_wrap(~variable, scales = 'free', ncol = 4) + 
    theme_bw() + 
    xlab("Read depth") + 
    ylab('') + 
    theme(text=element_text(size=21),
          legend.position = 'bottom',
          strip.background = element_rect(fill = "gray95")) + 
    labs(fill = 'Model') + 
    scale_fill_manual(values=c("#104E8B", "#458B00", "#EEAD0E", "#68228B"),
                      breaks=c("ALDEx2", "ANCOM-BC2", "MaAsLin 2", "MaAsLin 3")) + 
      geom_point(aes(x = get(param_name), y = 0), alpha = 0)
  ggsave(paste0(figures_folder, 'SD2_read_depth_figure.png'),
         plot = plot_out, width = 16, height = 4)
  
  melted_df$readDepth <- case_when(melted_df$readDepth == '5e+04' ~ '50K',
                                   melted_df$readDepth == '5e+07' ~ '50M',
                                   TRUE ~ NA)
  melted_df$readDepth <- factor(melted_df$readDepth, levels = c('50K', '50M'))
  melted_df <- melted_df[!grepl("No adjustment", melted_df$tool),]
  
  melted_df$variable <- case_when(melted_df$variable == 'Precision' ~ 'Precision',
                                  melted_df$variable == 'Recall' ~ 'Recall',
                                  melted_df$variable == 'Relative shrinkage error' ~ 'Relative shrinkage\nerror',
                                  melted_df$variable == 'Effect size Spearman cor.' ~ 'Effect size\nSpearman cor.')
  melted_df$variable <- factor(melted_df$variable, levels = c('Precision', 
                                                              'Recall',
                                                              'Relative shrinkage\nerror',
                                                              'Effect size\nSpearman cor.'))
  
  # In-text numbers
  melted_df %>%
      dplyr::filter(variable == 'Precision', readDepth %in% c("50K", "50M")) %>%
      dplyr::group_by(tool, readDepth) %>%
      dplyr::summarize(mean(value, na.rm=T))
  
  melted_df %>%
      dplyr::filter(variable == 'Recall', readDepth %in% c("50K", "50M")) %>%
      dplyr::group_by(tool, readDepth) %>%
      dplyr::summarize(mean(value, na.rm=T))
  
  melted_df %>%
      dplyr::filter(variable == 'Relative shrinkage\nerror', readDepth %in% c("50K", "50M")) %>%
      dplyr::group_by(tool, readDepth) %>%
      dplyr::summarize(mean(value, na.rm=T))
  
  plot_out <- ggplot(melted_df[melted_df$readDepth %in% c("50K", "50M"),], aes(x = get(param_name), y = value, fill = tool)) + 
      geom_boxplot(position = position_dodge(preserve = "single")) + 
      facet_wrap(~variable, scales = 'free', ncol = 4) + 
      theme_bw() + 
      xlab("Read depth") + 
      ylab('') + 
      theme(text=element_text(size=21),
            legend.position = 'bottom',
            strip.background = element_rect(fill = "gray95")) + 
      labs(fill = 'Model') + 
      scale_fill_manual(values=c("#104E8B", "#458B00", "#EEAD0E", "#68228B"),
                        breaks=c("ALDEx2", "ANCOM-BC2", "MaAsLin 2", "MaAsLin 3")) + 
      geom_point(aes(x = get(param_name), y = 0), alpha = 0)
  ggsave(paste0(figures_folder, 'SD2_read_depth_subfigure.png'),
         plot = plot_out, width = 12, height = 4)
}
SD2_read_depth_figure()

SD2_read_depth_roc_figure_maaslin3 <- function() {
    parameters <- "deep_sequencing/data_generation/SD2.txt"
    
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
    
    metadata_type <- "MVAref"
    param_name <- "readDepth"
    
    general_results_df <- data.frame(matrix(nrow = 0, ncol = 0))
    maaslin3_results_df <- data.frame(matrix(nrow = 0, ncol = 0))

    metrics <- c("Precision", "Recall")
    for (param_list in param_list_final[c(8)]) {
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
        
        inputDirectory <- file.path('Input', 'deep_sequencing', generator)
        this_params_folder <- file.path(inputDirectory, inputString)
        outputDirectory <- file.path('Output', 'deep_sequencing', generator)
        tmp_files <- list.files(outputDirectory)
        tmp_files <- tmp_files[grepl(paste0(inputString, "_"), tmp_files)]
        tools <- gsub(paste0(inputString, "_"), '', tmp_files, perl = T)
        this_output_folder <- file.path(outputDirectory, inputString)
        
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
        
        inputDirectory <- file.path('Input', 'deep_sequencing', generator)
        this_params_folder <- file.path(inputDirectory, inputString)
        outputDirectory <- file.path('Output', 'deep_sequencing', generator)
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
                
                if (grepl('Maaslin3', tool) & generator == 'SD2') {
                    tmp_truth <- truth[truth$associations == 'abundance',]
                    tmp_associations <- prepare_associations_maaslin3(associations, tool)
                    tmp_associations <- tmp_associations[tmp_associations$association == 'abundance',]
                    
                    thresholds = c(-1, 2^(c(seq(-300, -20, 5), seq(-20, 1, 0.5))), 0.1)
                    prec_recall <- t(sapply(thresholds, unweighted_precision_recall_maaslin3, truth = tmp_truth,
                                            associations = tmp_associations, abundance = abundance,
                                            metadata = metadata))
                    prec_recall <- data.frame(prec_recall)
                    colnames(prec_recall) <- c("Precision", "Recall")
                    prec_recall$thresholds <- thresholds
                    prec_recall$tool <- tool
                    prec_recall$iter <- i
                    prec_recall$association_type <- 'abundance'
                    prec_recall$readDepth <- param_list['readDepth']
                    
                    maaslin3_results_df <- plyr::rbind.fill(maaslin3_results_df, prec_recall)
                    
                    tmp_truth <- truth[truth$associations == 'prevalence',]
                    tmp_associations <- prepare_associations_maaslin3(associations, tool)
                    tmp_associations <- tmp_associations[tmp_associations$association == 'prevalence',]
                    
                    thresholds = c(-1, 2^(c(seq(-300, -20, 5), seq(-20, 1, 0.5))), 0.1)
                    prec_recall <- t(sapply(thresholds, unweighted_precision_recall_maaslin3, truth = tmp_truth,
                                            associations = tmp_associations, abundance = abundance,
                                            metadata = metadata))
                    prec_recall <- data.frame(prec_recall)
                    colnames(prec_recall) <- c("Precision", "Recall")
                    prec_recall$thresholds <- thresholds
                    prec_recall$tool <- tool
                    prec_recall$iter <- i
                    prec_recall$association_type <- 'prevalence'
                    prec_recall$readDepth <- param_list['readDepth']
                    
                    maaslin3_results_df <- plyr::rbind.fill(maaslin3_results_df, prec_recall)
                }
            }
        }
    }
    
    parameters <- "unscaled/data_generation/SD2.txt"
    
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
    
    metadata_type <- "MVAref"
    param_name <- "readDepth"
    
    metrics <- c("Precision", "Recall")
    for (param_list in param_list_final[c(18)]) {
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
        
        inputDirectory <- file.path('Input', 'unscaled', generator)
        this_params_folder <- file.path(inputDirectory, inputString)
        outputDirectory <- file.path('Output', 'unscaled', generator)
        tmp_files <- list.files(outputDirectory)
        tmp_files <- tmp_files[grepl(paste0(inputString, "_"), tmp_files)]
        tools <- gsub(paste0(inputString, "_"), '', tmp_files, perl = T)
        this_output_folder <- file.path(outputDirectory, inputString)
        
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
        
        inputDirectory <- file.path('Input', 'deep_sequencing', generator)
        this_params_folder <- file.path(inputDirectory, inputString)
        outputDirectory <- file.path('Output', 'deep_sequencing', generator)
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
                
                if (grepl('Maaslin3', tool) & generator == 'SD2') {
                    tmp_truth <- truth[truth$associations == 'abundance',]
                    tmp_associations <- prepare_associations_maaslin3(associations, tool)
                    tmp_associations <- tmp_associations[tmp_associations$association == 'abundance',]
                    
                    thresholds = c(-1, 2^(c(seq(-300, -20, 5), seq(-20, 1, 0.5))), 0.1)
                    prec_recall <- t(sapply(thresholds, unweighted_precision_recall_maaslin3, truth = tmp_truth,
                                            associations = tmp_associations, abundance = abundance,
                                            metadata = metadata))
                    prec_recall <- data.frame(prec_recall)
                    colnames(prec_recall) <- c("Precision", "Recall")
                    prec_recall$thresholds <- thresholds
                    prec_recall$tool <- tool
                    prec_recall$iter <- i
                    prec_recall$association_type <- 'abundance'
                    prec_recall$readDepth <- param_list['readDepth']
                    
                    maaslin3_results_df <- plyr::rbind.fill(maaslin3_results_df, prec_recall)
                    
                    tmp_truth <- truth[truth$associations == 'prevalence',]
                    tmp_associations <- prepare_associations_maaslin3(associations, tool)
                    tmp_associations <- tmp_associations[tmp_associations$association == 'prevalence',]
                    
                    thresholds = c(-1, 2^(c(seq(-300, -20, 5), seq(-20, 1, 0.5))), 0.1)
                    prec_recall <- t(sapply(thresholds, unweighted_precision_recall_maaslin3, truth = tmp_truth,
                                            associations = tmp_associations, abundance = abundance,
                                            metadata = metadata))
                    prec_recall <- data.frame(prec_recall)
                    colnames(prec_recall) <- c("Precision", "Recall")
                    prec_recall$thresholds <- thresholds
                    prec_recall$tool <- tool
                    prec_recall$iter <- i
                    prec_recall$association_type <- 'prevalence'
                    prec_recall$readDepth <- param_list['readDepth']
                    
                    maaslin3_results_df <- plyr::rbind.fill(maaslin3_results_df, prec_recall)
                }
            }
        }
    }
    
    names_to_iter <- names(param_list)[names(param_list) != 'metadataType']
    figures_folder <- paste0('Figures/paper_figures/')
    dir.create(figures_folder, recursive = T, showWarnings = F)
    
    melted_df <- maaslin3_results_df %>%
        dplyr::group_by(tool, association_type, thresholds, readDepth) %>%
        dplyr::summarise(Recall = mean(Recall, na.rm=T), Precision = mean(Precision, na.rm=T))
    
    melted_df$tool <- case_when(melted_df$tool == 'ALDEx2' ~ 'ALDEx2',
                                melted_df$tool == 'ANCOMBC' ~ 'ANCOM-BC2',
                                melted_df$tool == 'Maaslin2' ~ 'MaAsLin 2',
                                melted_df$tool == 'Maaslin3CompAdjust' ~ 'MaAsLin 3',
                                melted_df$tool == 'Maaslin3' ~ 'MaAsLin 3\n(No adjustment)',
                                melted_df$tool == 'Maaslin3unscaled' ~ 'MaAsLin 3\n(Spike-in)')
    melted_df$tool <- paste0(melted_df$tool, " ", ifelse(melted_df$association_type == 'abundance', "Abundance", "Prevalence"))
    
    tool_vec <- c("MaAsLin 3\n(No adjustment)", "MaAsLin 3", "MaAsLin 3\n(Spike-in)")
    
    melted_df$tool <- factor(melted_df$tool, levels = c(paste0(tool_vec, " ", "Abundance"),
                                                        paste0(tool_vec[2], " ", "Prevalence")))
    melted_df <- melted_df[!is.na(melted_df$tool),]
    melted_df$alpha <- melted_df$readDepth == max(melted_df$readDepth)
    melted_df$alpha <- ifelse(melted_df$alpha, 1, 0.5)
    
    plot_out <- ggplot(melted_df, 
           aes(x = Recall, y = Precision, color = tool,
               group = interaction(tool, readDepth))) +
        scale_alpha_identity() + 
        geom_line(linewidth = 1, aes(alpha = alpha)) +
        labs(
            title = "Increase read depth",
            x = "Recall",
            y = "Precision",
            color = "Model"
        ) +
        theme_bw() + 
        scale_x_continuous(limits = c(0, 1)) + 
        scale_y_continuous(limits = c(0,1)) + 
        theme(text=element_text(size=21), 
              title = element_text(size = 18),
              legend.position = 'none',
              strip.background = element_rect(fill = "gray95")) + 
        scale_color_manual(values=c("#EE7600", "#68228B", "#8B1A1A",
                                   "#F7BB80", "#996515", "#C58D8D"),
                          breaks=c(paste0(tool_vec, " ", "Abundance"),
                                   paste0(tool_vec, " ", "Prevalence"))) + 
        geom_point(data = melted_df[melted_df$thresholds == 0.1,], 
                   aes(x = Recall, y = Precision), 
                   size = 2) + 
        geom_segment(
            data = melted_df[melted_df$thresholds == 0.1,] %>%
                dplyr::group_by(tool) %>%
                dplyr::arrange(readDepth) %>%
                dplyr::mutate(xend = lead(Recall),
                       yend = lead(Precision),
                       xend = Recall + (lead(Recall) - Recall) * 0.96,  # Scale xend
                       yend = Precision + (lead(Precision) - Precision) * 0.96  # Scale yend
                ) %>%
                dplyr::filter(!is.na(xend)),
            aes(xend = xend, yend = yend),
            arrow = arrow(length = unit(0.2, "cm")),
        )
    
    ggsave(paste0(figures_folder, 'SD2_read_depth_roc_figure_maaslin3.png'),
           plot = plot_out, width = 5, height = 5)
}
# SD2_read_depth_roc_figure_maaslin3()
