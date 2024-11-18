remove(list = ls())
package_vec = c("reshape2", "ggplot2", "optparse", 
                "parallel", "stringi", "doParallel", "plyr", "tidyr", "dplyr")
invisible(suppressPackageStartupMessages(lapply(package_vec, require, character.only = TRUE)))
source('library/run_evaluation_helpers.R')

SD2_sample_size_figure <- function() {
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
  param_name <- "nSubjects"
  
  general_results_df <- data.frame(matrix(nrow = 0, ncol = 0))
  metrics <- c("Precision", "Recall", "Precision\n(common taxa)", "Recall\n(common taxa)", 
               "Relative effect error", "Relative shrinkage error", "Effect size\nSpearman cor.", "AUC", "AUC (Common taxa)", 
               "Issue proportion")
  for (param_list in param_list_final[c(2, 10, 12, 14, 16, 18)]) {
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
  figures_folder <- paste0('Figures/paper_figures/')
  dir.create(figures_folder, recursive = T, showWarnings = F)
  
  default_vals <- first_params
  names(default_vals) <- names_to_iter
  default_vals <- default_vals[-which(names(default_vals) == param_name)]
  
  results_subset <- general_results_df[apply(general_results_df[,names(default_vals)], 1, function(x) {all(x == default_vals)}),]
  results_subset <- results_subset[results_subset$metadataType == metadata_type,]
  results_subset <- results_subset[,c(metrics, param_name, 'tool', 'iter')]
  results_subset$F1 <- 2 / (1 / results_subset$Precision + 1 / results_subset$Recall)
  
  melted_df <- melt(results_subset, id.vars = c('tool', 'iter', param_name))
  melted_df[,param_name] <- factor(melted_df[,param_name], levels = as.character(sort(unique(as.numeric(melted_df[,param_name])))))
  melted_df <- melted_df[melted_df$variable != "Issue proportion",]
  melted_df$variable <- case_when(melted_df$variable == 'Precision' ~ 'Precision',
                                  melted_df$variable == 'Recall' ~ 'Recall',
                                  melted_df$variable == 'F1' ~ 'F1',
                                  melted_df$variable == 'Precision\n(common taxa)' ~ 'Precision (common taxa)',
                                  melted_df$variable == 'Recall\n(common taxa)' ~ 'Recall (common taxa)',
                                  melted_df$variable == 'AUC' ~ 'AUC',
                                  melted_df$variable == 'AUC (Common taxa)' ~ 'AUC (common taxa)',
                                  melted_df$variable == 'Relative effect error' ~ 'Relative effect error',
                                  melted_df$variable == 'Relative shrinkage error' ~ 'Relative shrinkage error',
                                  melted_df$variable == 'Effect size\nSpearman cor.' ~ 'Effect size correlation')
  melted_df$variable <- factor(melted_df$variable, levels = c('F1',
                                                              'Relative shrinkage error',
                                                              'Effect size correlation'))
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
  ggsave(paste0(figures_folder, 'SD2_sample_size_figure.png'),
         plot = plot_out, width = 12, height = 4.5)
}
SD2_sample_size_figure()

SD2_sample_size_sup_figure <- function() {
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
    param_name <- "nSubjects"
    
    general_results_df <- data.frame(matrix(nrow = 0, ncol = 0))
    metrics <- c("Precision", "Recall",
                 "Relative effect error", "Relative shrinkage error", "Effect size\nSpearman cor.",
                 "Issue proportion")
    for (param_list in param_list_final[c(2, 10, 12, 14, 16, 18)]) {
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
                             effect_size_error(truth, prepare_associations_abundance(associations, tool, generator)),
                             effect_size_correlation(truth, prepare_associations_abundance(associations, tool, generator)),
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
    melted_df$tool <- case_when(melted_df$tool == 'ALDEx2' ~ 'ALDEx2',
                                melted_df$tool == 'ANCOMBC' ~ 'ANCOM-BC2',
                                melted_df$tool == 'Maaslin2' ~ 'MaAsLin 2',
                                melted_df$tool == 'Maaslin3CompAdjust' ~ 'MaAsLin 3')
    melted_df$tool <- factor(melted_df$tool, levels = c("ALDEx2", "ANCOM-BC2", "MaAsLin 2", "MaAsLin 3"))
    melted_df <- melted_df[!is.na(melted_df$tool),]
    
    # In-text numbers
    melted_df %>%
        dplyr::filter(variable == 'Precision', tool %in% c("MaAsLin 2", "MaAsLin 3")) %>%
        dplyr::group_by(tool, nSubjects) %>%
        dplyr::summarize(mean(value, na.rm=T))
    
    melted_df %>%
        dplyr::filter(variable == 'Precision', tool %in% c("ALDEx2", "MaAsLin 3")) %>%
        dplyr::group_by(tool, nSubjects) %>%
        dplyr::summarize(mean(value, na.rm=T))
    
    melted_df %>%
        dplyr::filter(variable == 'Recall', nSubjects %in% c(20, 50)) %>%
        dplyr::group_by(tool, nSubjects) %>%
        dplyr::summarize(mean(value, na.rm=T))
    
    melted_df %>%
        dplyr::filter(variable == 'Precision', nSubjects %in% c("20", "50")) %>%
        dplyr::group_by(tool, nSubjects) %>%
        dplyr::summarize(mean(value, na.rm=T))
    
    cor_vals <- melted_df %>%
        dplyr::filter(variable == 'Effect size Spearman cor.', nSubjects != "20") %>%
        dplyr::group_by(tool, nSubjects) %>%
        dplyr::summarize(mean_out = mean(value, na.rm=T))
    
    range(cor_vals$mean_out, na.rm=T)
    
    melted_df %>%
        dplyr::filter(variable == 'Effect size Spearman cor.', nSubjects %in% c("1000")) %>%
        dplyr::group_by(tool, nSubjects) %>%
        dplyr::summarize(mean_out = mean(value, na.rm=T))
    
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
    ggsave(paste0(figures_folder, 'SD2_sample_size_sup_figure.png'),
           plot = plot_out, width = 8, height = 8)
}
SD2_sample_size_sup_figure()

SD2_repeated_samples_figure <- function() {
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
  param_name <- "nPerSubject"
  
  general_results_df <- data.frame(matrix(nrow = 0, ncol = 0))
  metrics <- c("Precision", "Recall", "Precision\n(common taxa)", "Recall\n(common taxa)", 
               "Relative effect error", "Relative shrinkage error", "Effect size\nSpearman cor.", "AUC", "AUC (Common taxa)", 
               "Issue proportion")
  for (param_list in param_list_final[c(2, 40, 42, 44)]) {
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
  melted_df$tool <- case_when(melted_df$tool == 'ALDEx2' ~ 'ALDEx2',
                              melted_df$tool == 'ANCOMBC' ~ 'ANCOM-BC2',
                              melted_df$tool == 'Maaslin2' ~ 'MaAsLin 2',
                              melted_df$tool == 'Maaslin3CompAdjust' ~ 'MaAsLin 3')
  melted_df$tool <- factor(melted_df$tool, levels = c("ALDEx2", "ANCOM-BC2", "MaAsLin 2", "MaAsLin 3"))
  melted_df <- melted_df[!is.na(melted_df$tool),]
  
  # In-text numbers
  melted_df %>%
      dplyr::filter(variable == 'Precision', tool %in% c("ANCOM-BC2")) %>%
      dplyr::group_by(tool, nPerSubject) %>%
      dplyr::summarize(mean(value, na.rm=T))
  
  melted_df %>%
      dplyr::filter(variable == 'Recall', tool %in% c("ANCOM-BC2")) %>%
      dplyr::group_by(tool, nPerSubject) %>%
      dplyr::summarize(mean(value, na.rm=T))
  
  melted_df %>%
      dplyr::filter(variable == 'Precision') %>%
      dplyr::group_by(tool, nPerSubject) %>%
      dplyr::summarize(mean(value, na.rm=T))
  
  melted_df %>%
      dplyr::filter(variable == 'Recall') %>%
      dplyr::group_by(tool, nPerSubject) %>%
      dplyr::summarize(mean(value, na.rm=T))
  
  melted_df %>%
      dplyr::filter(variable == 'Relative shrinkage error') %>%
      dplyr::group_by(tool, nPerSubject) %>%
      dplyr::summarize(mean(value, na.rm=T))
  
  melted_df %>%
      dplyr::filter(variable == 'Effect size Spearman cor.') %>%
      dplyr::group_by(tool, nPerSubject) %>%
      dplyr::summarize(mean(value, na.rm=T))
  
  # Remove few extreme ANCOM values that make plot worse
  melted_df <- melted_df[!is.na(melted_df$value) & melted_df$value <= 1,]
  
  plot_out <- ggplot(melted_df, aes(x = get(param_name), y = value, fill = tool)) + 
    geom_boxplot(position = position_dodge(preserve = "single")) + 
    facet_wrap(~variable, scales = 'free', ncol = 4) + 
    theme_bw() + 
    xlab("Samples per Subject") + 
    ylab('') + 
    theme(text=element_text(size=21),
          legend.position = 'bottom',
          strip.background = element_rect(fill = "gray95")) + 
    labs(fill = 'Model') + 
    scale_fill_manual(values=c("#104E8B", "#458B00", "#EEAD0E", "#68228B"),
                      breaks=c("ALDEx2", "ANCOM-BC2", "MaAsLin 2", "MaAsLin 3")) + 
    geom_point(aes(x = get(param_name), y = 0), alpha = 0)
  ggsave(paste0(figures_folder, 'SD2_repeated_samples_figure.png'),
         plot = plot_out, width = 16, height = 4)
}
SD2_repeated_samples_figure()

SD2_repeated_samples_investigation <- function() {
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
    param_name <- "nPerSubject"
    
    general_results_df <- data.frame(matrix(nrow = 0, ncol = 0))
    metrics <- c("prop_small_precision")
    for (param_list in param_list_final[c(2, 40, 42, 44)]) {
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
            for (tool in tools[grepl('Maaslin3', tools)]) {
                possible_error <- tryCatch({
                    associations <- read.csv(paste0(this_output_folder, '_', tool, "/associations_", i, ".tsv"), sep = "\t")
                }, error = function(err) {
                    err
                })
                
                if(inherits(possible_error, "error")) next
                
                tmp_metric <- function(truth, associations, abundance, metadata, threshold = 0.1) {
                    if (nrow(associations) == 0) {
                        return(c(NA, NA))
                    }
                    truth_match_vec <- unique(paste0(truth$taxon, '_', truth$metadata))
                    associations_match_vec <- paste0(associations$taxon, '_', associations$metadata)

                    tmp_associations <- associations[
                        !associations_match_vec %in% truth_match_vec & associations$signif < threshold,]
                    
                    if (nrow(tmp_associations) > 0) {
                        return(mean(abs(tmp_associations$effect_size) < 1 & tmp_associations$association == 'prevalence'))
                    } else {
                        return(NA)
                    }
                }

                new_row <- tmp_metric(truth, prepare_associations_maaslin3(associations, tool), abundance, metadata)
                new_row <- lapply(new_row, function(x) {x})
                names(new_row) <- 'prop_small_precision'
                
                new_row[['tool']] <- tool
                new_row[['iter']] <- i
                new_row <- c(new_row, param_list)
                general_results_df <- plyr::rbind.fill(general_results_df, data.frame(new_row, check.names = F))
            }
        }
    }
    
    names_to_iter <- names(param_list)[names(param_list) != 'metadataType']

    default_vals <- first_params
    names(default_vals) <- names_to_iter
    default_vals <- default_vals[-which(names(default_vals) == param_name)]
    
    results_subset <- general_results_df[apply(general_results_df[,names(default_vals)], 1, function(x) {all(x == default_vals)}),]
    results_subset <- results_subset[results_subset$metadataType == metadata_type,]
    results_subset <- results_subset[,c(metrics, param_name, 'tool', 'iter')]

    melted_df <- melt(results_subset, id.vars = c('tool', 'iter', param_name))
    melted_df[,param_name] <- factor(melted_df[,param_name], levels = as.character(sort(unique(as.numeric(melted_df[,param_name])))))
    melted_df$tool <- case_when(melted_df$tool == 'ALDEx2' ~ 'ALDEx2',
                                melted_df$tool == 'ANCOMBC' ~ 'ANCOM-BC2',
                                melted_df$tool == 'Maaslin2' ~ 'MaAsLin 2',
                                melted_df$tool == 'Maaslin3CompAdjust' ~ 'MaAsLin 3')
    melted_df$tool <- factor(melted_df$tool, levels = c("ALDEx2", "ANCOM-BC2", "MaAsLin 2", "MaAsLin 3"))
    melted_df <- melted_df[!is.na(melted_df$tool),]
    
    # In-text numbers
    melted_df %>%
        dplyr::filter(variable == 'prop_small_precision') %>%
        dplyr::group_by(tool, nPerSubject) %>%
        dplyr::summarize(mean(value, na.rm=T))
}
# SD2_repeated_samples_investigation()

SD2_repeated_samples_thresholded_figure <- function() {
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
    param_name <- "nPerSubject"
    
    general_results_df <- data.frame(matrix(nrow = 0, ncol = 0))
    metrics <- c("Precision", "Recall", "Precision\n(common taxa)", "Recall\n(common taxa)", 
                 "Relative effect error", "Relative shrinkage error", "Effect size\nSpearman cor.", "AUC", "AUC (Common taxa)", 
                 "Issue proportion")
    for (param_list in param_list_final[c(2, 40, 42, 44)]) {
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
                
                new_row <- c(unweighted_precision_recall(truth, prepare_associations_general(associations, tool, generator, 1), abundance, metadata),
                             weighted_precision_recall(truth, prepare_associations_general(associations, tool, generator, 1), abundance, metadata),
                             effect_size_error(truth, prepare_associations_abundance(associations, tool, generator, 1)),
                             effect_size_correlation(truth, prepare_associations_abundance(associations, tool, generator, 1)),
                             pval_auc(truth, prepare_associations_general(associations, tool, generator, 1), abundance, metadata),
                             weighted_pval_auc(truth, prepare_associations_general(associations, tool, generator, 1), abundance, metadata),
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
    melted_df$tool <- case_when(melted_df$tool == 'ALDEx2' ~ 'ALDEx2',
                                melted_df$tool == 'ANCOMBC' ~ 'ANCOM-BC2',
                                melted_df$tool == 'Maaslin2' ~ 'MaAsLin 2',
                                melted_df$tool == 'Maaslin3CompAdjust' ~ 'MaAsLin 3')
    melted_df$tool <- factor(melted_df$tool, levels = c("ALDEx2", "ANCOM-BC2", "MaAsLin 2", "MaAsLin 3"))
    melted_df <- melted_df[!is.na(melted_df$tool),]
    
    melted_df %>%
        dplyr::filter(variable == 'Precision') %>%
        dplyr::group_by(tool, nPerSubject) %>%
        dplyr::summarize(mean(value, na.rm=T))
    
    melted_df %>%
        dplyr::filter(variable == 'Recall') %>%
        dplyr::group_by(tool, nPerSubject) %>%
        dplyr::summarize(mean(value, na.rm=T))
    
    melted_df %>%
        dplyr::filter(variable == 'Relative shrinkage error') %>%
        dplyr::group_by(tool, nPerSubject) %>%
        dplyr::summarize(mean(value, na.rm=T))
    
    melted_df %>%
        dplyr::filter(variable == 'Effect size Spearman cor.') %>%
        dplyr::group_by(tool, nPerSubject) %>%
        dplyr::summarize(mean(value, na.rm=T))
    
    # Remove few extreme ANCOM values that make plot worse
    melted_df <- melted_df[!is.na(melted_df$value) & melted_df$value <= 1,]
    
    plot_out <- ggplot(melted_df, aes(x = get(param_name), y = value, fill = tool)) + 
        geom_boxplot(position = position_dodge(preserve = "single")) + 
        facet_wrap(~variable, scales = 'free', ncol=4) + 
        theme_bw() + 
        xlab("Samples per Subject") + 
        ylab('') + 
        theme(text=element_text(size=21),
              legend.position = 'bottom',
              strip.background = element_rect(fill = "gray95")) + 
        labs(fill = 'Model') + 
        scale_fill_manual(values=c("#104E8B", "#458B00", "#EEAD0E", "#68228B"),
                          breaks=c("ALDEx2", "ANCOM-BC2", "MaAsLin 2", "MaAsLin 3")) + 
        geom_point(aes(x = get(param_name), y = 0), alpha = 0)
    ggsave(paste0(figures_folder, 'SD2_repeated_samples_thresholded_figure.png'),
           plot = plot_out, width = 16, height = 4)
}
SD2_repeated_samples_thresholded_figure()

SD2_effect_size_figure <- function() {
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
  param_name <- "effectSize"
  
  general_results_df <- data.frame(matrix(nrow = 0, ncol = 0))
  metrics <- c("Precision", "Recall", "Precision\n(common taxa)", "Recall\n(common taxa)", 
               "Relative effect error", "Relative shrinkage error", "Effect size\nSpearman cor.", "AUC", "AUC (Common taxa)", 
               "Issue proportion")
  for (param_list in param_list_final[c(2, 24, 26, 28, 30)]) {
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
  melted_df$tool <- case_when(melted_df$tool == 'ALDEx2' ~ 'ALDEx2',
                              melted_df$tool == 'ANCOMBC' ~ 'ANCOM-BC2',
                              melted_df$tool == 'Maaslin2' ~ 'MaAsLin 2',
                              melted_df$tool == 'Maaslin3CompAdjust' ~ 'MaAsLin 3')
  melted_df$tool <- factor(melted_df$tool, levels = c("ALDEx2", "ANCOM-BC2", "MaAsLin 2", "MaAsLin 3"))
  melted_df <- melted_df[!is.na(melted_df$tool),]
  
  melted_df <- melted_df[melted_df$effectSize != "0.5",]
  melted_df$effectSize <- case_when(melted_df$effectSize == '1' ~ "0.5 - 1",
                                    melted_df$effectSize == '2' ~ "1 - 2",
                                    melted_df$effectSize == '5' ~ "2.5 - 5",
                                    melted_df$effectSize == '10' ~ "5 - 10")
  
  # In-text numbers
  melted_df %>%
      dplyr::filter(variable == 'Precision') %>%
      dplyr::group_by(tool, effectSize) %>%
      dplyr::summarize(mean(value, na.rm=T))
  
  melted_df %>%
      dplyr::filter(variable == 'Recall') %>%
      dplyr::group_by(tool, effectSize) %>%
      dplyr::summarize(mean(value, na.rm=T))
  
  plot_out <- ggplot(melted_df, aes(x = get(param_name), y = value, fill = tool)) + 
    geom_boxplot(position = position_dodge(preserve = "single")) + 
    facet_wrap(~variable, scales = 'free') + 
    theme_bw() + 
    xlab("Effect Size") + 
    ylab('') + 
    theme(text=element_text(size=21),
          legend.position = 'bottom',
          strip.background = element_rect(fill = "gray95")) + 
    labs(fill = 'Model') + 
    scale_fill_manual(values=c("#104E8B", "#458B00", "#EEAD0E", "#68228B"),
                      breaks=c("ALDEx2", "ANCOM-BC2", "MaAsLin 2", "MaAsLin 3")) + 
      geom_point(aes(x = get(param_name), y = 0), alpha = 0)
  ggsave(paste0(figures_folder, 'SD2_effect_size_figure.png'),
         plot = plot_out, width = 9, height = 8)
}
SD2_effect_size_figure()

SD2_sample_size_figure_maaslin3 <- function() {
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
    
    inputDirectory <- file.path('Input', 'unscaled', generator)
    this_params_folder <- file.path(inputDirectory, inputString)
    outputDirectory <- file.path('Output', 'unscaled', generator)
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
      dplyr::filter(association_type == 'abundance', variable == 'Precision', nSubjects == "1000") %>%
      dplyr::group_by(tool) %>%
      dplyr::summarise(mean(value, na.rm=T))
  
  melted_df %>%
      dplyr::filter(association_type == 'abundance', variable == 'Recall', nSubjects == "1000") %>%
      dplyr::group_by(tool) %>%
      dplyr::summarise(mean(value, na.rm=T))
  
  plot_out <- ggplot(melted_df, aes(x = get(param_name), y = value, fill = tool)) + 
    geom_boxplot(position = position_dodge(preserve = "single")) + 
    facet_wrap(~variable, scales = 'free', ncol = 4) + 
    theme_bw() + 
    xlab("Number of Subjects") + 
    ylab('') + 
    theme(text=element_text(size=21),
          legend.position = 'bottom',
          strip.background = element_rect(fill = "gray95")) + 
    labs(fill = 'Model') + 
    scale_fill_manual(values=c("#EE7600", "#68228B", "#8B1A1A",
                               "#F7BB80", "#996515", "#C58D8D"),
                      breaks=c(paste0(tool_vec, " ", "Abundance"),
                               paste0(tool_vec, " ", "Prevalence")))
  ggsave(paste0(figures_folder, 'SD2_sample_size_figure_maaslin3.png'),
         plot = plot_out, width = 16, height = 4)
  
  plot_out <- ggplot(melted_df[melted_df$variable %in% c("Precision", "Recall"),], 
                     aes(x = get(param_name), y = value, fill = tool)) + 
      geom_boxplot(position = position_dodge(preserve = "single")) + 
      facet_wrap(~variable, scales = 'free', ncol = 4) + 
      theme_bw() + 
      xlab("Number of Subjects") + 
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
  ggsave(paste0(figures_folder, 'SD2_sample_size_subfigure_maaslin3.png'),
         plot = plot_out, width = 10.5, height = 5)
}
SD2_sample_size_figure_maaslin3()

SD2_threshold_figure_maaslin3 <- function() {
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
    
    inputDirectory <- file.path('Input', 'unscaled', generator)
    this_params_folder <- file.path(inputDirectory, inputString)
    outputDirectory <- file.path('Output', 'unscaled', generator)
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
        
        new_row <- c(unweighted_precision_recall(truth, prepare_associations_general(associations, tool, generator, threshold = 1), abundance, metadata),
                     weighted_precision_recall(truth, prepare_associations_general(associations, tool, generator, threshold = 1), abundance, metadata),
                     effect_size_error(truth, prepare_associations_abundance(associations, tool, generator, threshold = 1)),
                     effect_size_correlation(truth, prepare_associations_abundance(associations, tool, generator, threshold = 1)),
                     pval_auc(truth, prepare_associations_general(associations, tool, generator, threshold = 1), abundance, metadata),
                     weighted_pval_auc(truth, prepare_associations_general(associations, tool, generator, threshold = 1), abundance, metadata),
                     issue_prop(associations, tool))
        new_row <- lapply(new_row, function(x) {x})
        names(new_row) <- metrics
        
        new_row[['tool']] <- tool
        new_row[['iter']] <- i
        new_row <- c(new_row, param_list)
        general_results_df <- plyr::rbind.fill(general_results_df, data.frame(new_row, check.names = F))
        
        if (grepl('Maaslin3', tool) & generator == 'SD2') {
          tmp_truth <- truth[truth$associations == 'abundance',]
          tmp_associations <- prepare_associations_maaslin3(associations, tool, threshold = 1)
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
          tmp_associations <- prepare_associations_maaslin3(associations, tool, threshold = 1)
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
                              melted_df$tool == 'Maaslin3CompAdjust' ~ 'MaAsLin 3',
                              melted_df$tool == 'Maaslin3' ~ 'MaAsLin 3\n(No adjustment)',
                              melted_df$tool == 'Maaslin3unscaled' ~ 'MaAsLin 3\n(Spike-in)')
  melted_df$tool <- paste0(melted_df$tool, " ", ifelse(melted_df$association_type == 'abundance', "Abundance", "Prevalence"))
  
  tool_vec <- c("MaAsLin 3\n(No adjustment)", "MaAsLin 3", "MaAsLin 3\n(Spike-in)")
  
  melted_df$tool <- factor(melted_df$tool, levels = c(paste0(tool_vec, " ", "Abundance"),
                                                      paste0(tool_vec[2], " ", "Prevalence")))
  melted_df <- melted_df[!is.na(melted_df$tool),]
  
  plot_out <- ggplot(melted_df, aes(x = get(param_name), y = value, fill = tool)) + 
    geom_boxplot(position = position_dodge(preserve = "single")) + 
    facet_wrap(~variable, scales = 'free', ncol = 4) + 
    theme_bw() + 
    xlab("Number of Subjects") + 
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
  ggsave(paste0(figures_folder, 'SD2_threshold_figure_maaslin3.png'),
         plot = plot_out, width = 16, height = 4)
}
SD2_threshold_figure_maaslin3()

SD2_threshold_figure <- function() {
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
  param_name <- "nSubjects"
  
  general_results_df <- data.frame(matrix(nrow = 0, ncol = 0))
  metrics <- c("Precision", "Recall", "Precision\n(common taxa)", "Recall\n(common taxa)", 
               "Relative effect error", "Relative shrinkage error", "Effect size\nSpearman cor.", "AUC", "AUC (Common taxa)", 
               "Issue proportion")
  for (param_list in param_list_final[c(2, 10, 12, 14, 16, 18)]) {
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
        
        new_row <- c(unweighted_precision_recall(truth, prepare_associations_general(associations, tool, generator, threshold = 1), abundance, metadata),
                     weighted_precision_recall(truth, prepare_associations_general(associations, tool, generator, threshold = 1), abundance, metadata),
                     effect_size_error(truth, prepare_associations_abundance(associations, tool, generator, threshold = 1)),
                     effect_size_correlation(truth, prepare_associations_abundance(associations, tool, generator, threshold = 1)),
                     pval_auc(truth, prepare_associations_general(associations, tool, generator, threshold = 1), abundance, metadata),
                     weighted_pval_auc(truth, prepare_associations_general(associations, tool, generator, threshold = 1), abundance, metadata),
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
    xlab("Number of Subjects") + 
    ylab('') + 
    theme(text=element_text(size=21),
          legend.position = 'bottom',
          strip.background = element_rect(fill = "gray95")) + 
    labs(fill = 'Model') + 
    scale_fill_manual(values=c("#104E8B", "#458B00", "#EEAD0E", "#68228B"),
                      breaks=c("ALDEx2", "ANCOM-BC2", "MaAsLin 2", "MaAsLin 3")) + 
      geom_point(aes(x = get(param_name), y = 0), alpha = 0)
  ggsave(paste0(figures_folder, 'SD2_threshold_figure.png'),
         plot = plot_out, width = 16, height = 4)
}
SD2_threshold_figure()

SD2_threshold_pr_figure_maaslin3 <- function() {
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
    
    general_results_df <- data.frame(matrix(nrow = 0, ncol = 0))
    maaslin3_results_df <- data.frame(matrix(nrow = 0, ncol = 0))
    
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
                    prec_recall$threshold = 0
                    
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
                    prec_recall$threshold = 0
                    
                    maaslin3_results_df <- plyr::rbind.fill(maaslin3_results_df, prec_recall)
                }
                
                if (grepl('Maaslin3', tool) & generator == 'SD2') {
                    tmp_truth <- truth[truth$associations == 'abundance',]
                    tmp_associations <- prepare_associations_maaslin3(associations, tool, threshold = 1)
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
                    prec_recall$threshold = 1
                    
                    maaslin3_results_df <- plyr::rbind.fill(maaslin3_results_df, prec_recall)
                    
                    tmp_truth <- truth[truth$associations == 'prevalence',]
                    tmp_associations <- prepare_associations_maaslin3(associations, tool, threshold = 1)
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
                    prec_recall$threshold = 1
                    
                    maaslin3_results_df <- plyr::rbind.fill(maaslin3_results_df, prec_recall)
                }
            }
        }
    }
    
    names_to_iter <- names(param_list)[names(param_list) != 'metadataType']
    figures_folder <- paste0('Figures/paper_figures/')
    dir.create(figures_folder, recursive = T, showWarnings = F)
    
    melted_df <- maaslin3_results_df %>%
        group_by(tool, association_type, thresholds, threshold) %>%
        summarise(Recall = mean(Recall, na.rm=T), Precision = mean(Precision, na.rm=T))
    
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
    melted_df$alpha <- ifelse(melted_df$threshold == 0, 0.5, 1)
    
    plot_out <- ggplot(melted_df, 
                       aes(x = Recall, y = Precision, color = tool,
                           group = interaction(tool, threshold))) +
        scale_alpha_identity() + 
        geom_line(linewidth = 1, aes(alpha = alpha)) +
        labs(
            title = "Threshold coefficients at 1",
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
                group_by(tool) %>%
                arrange(threshold) %>%
                mutate(xend = lead(Recall),
                       yend = lead(Precision),
                       xend = Recall + (lead(Recall) - Recall) * 0.95,  # Scale xend
                       yend = Precision + (lead(Precision) - Precision) * 0.95  # Scale yend
                ) %>%
                filter(!is.na(xend)),
            aes(xend = xend, yend = yend),
            arrow = arrow(length = unit(0.2, "cm")),
        )
    
    ggsave(paste0(figures_folder, 'SD2_threshold_pr_figure_maaslin3.png'),
           plot = plot_out, width = 5, height = 5)
}
# SD2_threshold_pr_figure_maaslin3()

SD2_improvement_pr_figure_maaslin3 <- function() {
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
    
    general_results_df <- data.frame(matrix(nrow = 0, ncol = 0))
    maaslin3_results_df <- data.frame(matrix(nrow = 0, ncol = 0))
    
    metrics <- c("Precision", "Recall")
    for (param_list in param_list_final[c(2, 18)]) {
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
                    prec_recall$nSubjects = param_list[['nSubjects']]
                    prec_recall$allow_abun_to_prev = F
                    
                    maaslin3_results_df <- plyr::rbind.fill(maaslin3_results_df, prec_recall)
                    
                    tmp_truth <- truth[truth$associations == 'prevalence',]
                    tmp_associations <- prepare_associations_maaslin3(associations, tool, allow_abun_to_prev = F)
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
                    prec_recall$nSubjects = param_list[['nSubjects']]
                    prec_recall$allow_abun_to_prev = F
                    
                    maaslin3_results_df <- plyr::rbind.fill(maaslin3_results_df, prec_recall)
                    
                    tmp_truth <- truth[truth$associations == 'prevalence',]
                    tmp_associations <- prepare_associations_maaslin3(associations, tool, allow_abun_to_prev = T)
                    tmp_associations <- tmp_associations[tmp_associations$association == 'prevalence',]
                    
                    thresholds = c(-1, 2^(c(seq(-300, -20, 2), seq(-20, 1, 0.2))), 0.1, 1, 2)
                    prec_recall <- t(sapply(thresholds, unweighted_precision_recall_maaslin3, truth = tmp_truth,
                                            associations = tmp_associations, abundance = abundance,
                                            metadata = metadata))
                    prec_recall <- data.frame(prec_recall)
                    colnames(prec_recall) <- c("Precision", "Recall")
                    prec_recall$thresholds <- thresholds
                    prec_recall$tool <- tool
                    prec_recall$iter <- i
                    prec_recall$association_type <- 'prevalence'
                    prec_recall$nSubjects = param_list[['nSubjects']]
                    prec_recall$allow_abun_to_prev = T
                    
                    maaslin3_results_df <- plyr::rbind.fill(maaslin3_results_df, prec_recall)
                }
            }
        }
    }
    
    names_to_iter <- names(param_list)[names(param_list) != 'metadataType']
    figures_folder <- paste0('Figures/paper_figures/')
    dir.create(figures_folder, recursive = T, showWarnings = F)
    
    melted_df <- maaslin3_results_df %>%
        dplyr::group_by(tool, association_type, thresholds, nSubjects, allow_abun_to_prev) %>%
        dplyr::summarise(Recall = mean(Recall, na.rm=T), Precision = mean(Precision, na.rm=T))
    
    melted_df$tool <- case_when(melted_df$tool == 'ALDEx2' ~ 'ALDEx2',
                                melted_df$tool == 'ANCOMBC' ~ 'ANCOM-BC2',
                                melted_df$tool == 'Maaslin2' ~ 'MaAsLin 2',
                                melted_df$tool == 'Maaslin3CompAdjust' ~ 'MaAsLin 3',
                                melted_df$tool == 'Maaslin3' ~ 'MaAsLin 3 (No median adjustment)',
                                melted_df$tool == 'Maaslin3unscaled' ~ 'MaAsLin 3 (Spike-in)',
                                melted_df$tool == 'Maaslin3NoAugment' ~ 'MaAsLin 3 (No augmentation)')

    tool_vec <- c("MaAsLin 3 (No median adjustment)", "MaAsLin 3", "MaAsLin 3 (Spike-in)")
    
    melted_df$tool <- factor(melted_df$tool, levels = tool_vec)
    melted_df <- melted_df[!is.na(melted_df$tool),]
    melted_df$linetype <- factor(ifelse(melted_df$nSubjects == 100, "100", "1000"))
    melted_df <- melted_df[melted_df$association_type == 'abundance',]
    
    plot_out <- ggplot(melted_df, 
                       aes(x = Recall, y = Precision, color = tool, linetype = linetype,
                           group = interaction(tool, nSubjects))) +
        scale_linetype_manual(values = c("100" = "twodash", "1000" = "solid"),
                              guide = "legend") + 
        geom_line(linewidth = 1) +
        labs(
            title = "Abundance association detection",
            x = "Recall",
            y = "Precision",
            color = "Model",
            linetype = "Number of Subjects"
        ) +
        theme_bw() + 
        scale_x_continuous(limits = c(0, 1)) + 
        scale_y_continuous(limits = c(0,1)) + 
        theme(text=element_text(size=21),
              title = element_text(size = 18),
              legend.position = 'right',
              legend.direction = 'vertical',
              strip.background = element_rect(fill = "gray95")) + 
        scale_color_manual(values=c("#EE7600", "#68228B", "#8B1A1A"),
                           breaks=tool_vec)
    
    ggsave(paste0(figures_folder, 'SD2_improvement_abundance_pr_figure_maaslin3.png'),
           plot = plot_out, width = 9, height = 5)
    
    melted_df <- maaslin3_results_df %>%
        dplyr::group_by(tool, association_type, thresholds, nSubjects, allow_abun_to_prev) %>%
        dplyr::summarise(Recall = mean(Recall, na.rm=T), Precision = mean(Precision, na.rm=T))
    
    melted_df$tool <- case_when(melted_df$tool == 'ALDEx2' ~ 'ALDEx2',
                                melted_df$tool == 'ANCOMBC' ~ 'ANCOM-BC2',
                                melted_df$tool == 'Maaslin2' ~ 'MaAsLin 2',
                                melted_df$tool == 'Maaslin3CompAdjust' ~ 'MaAsLin 3',
                                melted_df$tool == 'Maaslin3' ~ 'MaAsLin 3 (No adjustment)',
                                melted_df$tool == 'Maaslin3unscaled' ~ 'MaAsLin 3 (Spike-in)',
                                melted_df$tool == 'Maaslin3NoAugment' ~ 'MaAsLin 3 (No augmentation)')
    melted_df$tool <- ifelse(melted_df$allow_abun_to_prev, 
                             paste0(melted_df$tool, " (All associations)"),
                             paste0(melted_df$tool, " (Passed screen)"))
    
    melted_df$tool <- factor(melted_df$tool, levels = c(paste0(tool_vec[2], " ", "(All associations)"),
                                                        paste0(tool_vec[2], " ", "(Passed screen)"),
                                                        paste0('MaAsLin 3 (No augmentation)', " ", "(All associations)")))
    
    melted_df <- melted_df[!is.na(melted_df$tool),]
    melted_df <- melted_df %>%
        mutate(tool = case_when(tool == 'MaAsLin 3 (All associations)' ~ 'MaAsLin 3 (No prevalence screen)',
                                tool == 'MaAsLin 3 (Passed screen)' ~ 'MaAsLin 3',
                                tool == 'MaAsLin 3 (No augmentation) (All associations)' ~ 'MaAsLin 3 (No augmentation)') %>%
                   factor(levels = c('MaAsLin 3 (No augmentation)', 'MaAsLin 3 (No prevalence screen)', 'MaAsLin 3')),
               )
    melted_df$linetype <- factor(ifelse(melted_df$nSubjects == 100, "100", "1000"))
    melted_df <- melted_df[melted_df$association_type == 'prevalence',]
    
    plot_out <- ggplot(melted_df, 
                       aes(x = Recall, y = Precision, color = tool, linetype = linetype,
                           group = interaction(tool, nSubjects))) +
        scale_linetype_manual(values = c("100" = "twodash", "1000" = "solid"),
                              guide = "legend") + 
        geom_line(linewidth = 1) +
        labs(
            title = "Prevalence association detection",
            x = "Recall",
            y = "Precision",
            color = "Model",
            linetype = "Number of Subjects"
        ) +
        theme_bw() + 
        scale_x_continuous(limits = c(0, 1)) + 
        scale_y_continuous(limits = c(0,1)) + 
        theme(text=element_text(size=21),
              title = element_text(size = 18),
              legend.position = 'right',
              legend.direction = 'vertical',
              strip.background = element_rect(fill = "gray95")) + 
        scale_color_manual(values=c("#D53E4F", "#26867d", "#68228B"),
                           breaks=c('MaAsLin 3 (No augmentation)', 'MaAsLin 3 (No prevalence screen)', 'MaAsLin 3'))
    
    ggsave(paste0(figures_folder, 'SD2_improvement_precision_pr_figure_maaslin3.png'),
           plot = plot_out, width = 9, height = 5)
}
SD2_improvement_pr_figure_maaslin3()

SD2_sample_size_prevalence_only_recall <- function() {
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
    param_name <- "nSubjects"

    general_results_df <- data.frame(matrix(nrow = 0, ncol = 0))
    metrics <- c("Precision", "Recall")
    for (param_list in param_list_final[c(2, 10, 12, 14, 16, 18)]) {
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
            truth <- truth[truth$associations == 'prevalence',]
            for (tool in tools) {
                possible_error <- tryCatch({
                    associations <- read.csv(paste0(this_output_folder, '_', tool, "/associations_", i, ".tsv"), sep = "\t")
                }, error = function(err) {
                    err
                })

                if(inherits(possible_error, "error")) next

                new_row <- c(unweighted_precision_recall(truth, prepare_associations_general(associations, tool, generator), abundance, metadata))
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
    melted_df$tool <- case_when(melted_df$tool == 'ALDEx2' ~ 'ALDEx2',
                                melted_df$tool == 'ANCOMBC' ~ 'ANCOM-BC2',
                                melted_df$tool == 'Maaslin2' ~ 'MaAsLin 2',
                                melted_df$tool == 'Maaslin3CompAdjust' ~ 'MaAsLin 3')
    melted_df$tool <- factor(melted_df$tool, levels = c("ALDEx2", "ANCOM-BC2", "MaAsLin 2", "MaAsLin 3"))
    melted_df <- melted_df[!is.na(melted_df$tool),]

    # In-text numbers
    melted_df %>%
        dplyr::filter(variable == 'Recall' & nSubjects == "1000") %>%
        dplyr::group_by(tool, nSubjects) %>%
        dplyr::summarize(mean(value, na.rm=T))
}
SD2_sample_size_prevalence_only_recall()


