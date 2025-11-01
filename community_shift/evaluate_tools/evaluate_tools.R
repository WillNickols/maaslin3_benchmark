remove(list = ls())
package_vec = c("reshape2", "ggplot2", "optparse", 
                "parallel", "stringi", "doParallel", "plyr", "tidyr")
invisible(suppressPackageStartupMessages(lapply(package_vec, require, character.only = TRUE)))
source('library/run_evaluation_helpers.R')

SD2_spike_microbe_figure <- function() {
  parameters <- "community_shift/data_generation/SD2.txt"
  
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
      if (tmp_param_file['depthConfound'] == 'TRUE'){
          tmp_param_file['nMetadata'] = '1'
      }
      param_list_final <- c(param_list_final, list(tmp_param_file))
    }
  }
  
  metadata_type <- "MVAref"
  param_name <- "spikeMicrobes"
  
  general_results_df <- data.frame(matrix(nrow = 0, ncol = 0))
  metrics <- c("Precision", "Recall",
               "Effect size mean diff", "Effect size correlation")
  for (param_list in param_list_final[c(2, 4, 6, 8, 10, 12, 14, 16, 18)]) {
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
                       param_list[['depthConfound']],
                       param_list[['propAbun']],
                       param_list[['zeroInflate']],
                       sep='_')
    options("scipen"=5)
    
    inputDirectory <- file.path('Input', 'community_shift', generator)
    this_params_folder <- file.path(inputDirectory, inputString)
    outputDirectory <- file.path('Output', 'community_shift', generator)
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
                     effect_size_mean_diff(truth, prepare_associations_abundance(associations, tool, generator)),
                     effect_size_correlation(truth, prepare_associations_abundance(associations, tool, generator)))

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
                                  melted_df$variable == 'Effect size mean diff' ~ 'Effect size bias',
                                  melted_df$variable == 'Effect size correlation' ~ 'Effect size correlation')
  melted_df$variable <- factor(melted_df$variable, levels = c('Precision', 
                                                              'Recall',
                                                              'Effect size bias', 
                                                              'Effect size correlation'))
  melted_df <- melted_df[!is.na(melted_df$variable),]
  melted_df$tool <- case_when(melted_df$tool == 'ALDEx2' ~ 'ALDEx2',
                              melted_df$tool == 'ALDEx2_scale' ~ 'ALDEx2\n(scale informed)',
                              melted_df$tool == 'ANCOMBC' ~ 'ANCOM-BC2',
                              melted_df$tool == 'DESeq2' ~ 'DESeq2',
                              melted_df$tool == 'edgeR' ~ 'edgeR',
                              melted_df$tool == 'Maaslin2' ~ 'MaAsLin 2',
                              melted_df$tool == 'Maaslin3CompAdjust' ~ 'MaAsLin 3',
                              melted_df$tool == 'Maaslin3unscaled' ~ 'MaAsLin 3\n(Spike-in)')
  melted_df$tool <- factor(melted_df$tool, 
    levels = c("ALDEx2", 'ALDEx2\n(scale informed)', "ANCOM-BC2", 'DESeq2', 
               'edgeR', "MaAsLin 2", "MaAsLin 3", "MaAsLin 3\n(Spike-in)"))
  melted_df <- melted_df[!is.na(melted_df$tool),]
  
  # Make plotting dimensions reasonable
  melted_df <- melted_df[!is.na(melted_df$value) & melted_df$value < 5 & melted_df$value > -6,]
  d <- melted_df[melted_df$spikeMicrobes %in% c("0.1", "0.3", "0.5", "0.7", "0.9"),]
  
  # In-text numbers
  melted_df %>%
      dplyr::filter(variable == 'Effect size bias' & spikeMicrobes == "0.9") %>%
      dplyr::group_by(tool) %>%
      dplyr::summarize(mean(value, na.rm=T))
  
  plot_out <- ggplot(melted_df, aes(x = get(param_name), y = value, fill = tool)) + 
    geom_boxplot(position = position_dodge(preserve = "single")) + 
    facet_wrap(~variable, scales = 'free') + 
    theme_bw() + 
    xlab("Proportion of associations non-null") + 
    ylab('') + 
    theme(text=element_text(size=21),
          legend.position = 'right',
          strip.background = element_rect(fill = "gray95")) + 
    labs(fill = 'Model') + 
    scale_fill_manual(values=c("#104E8B", "#5B9BD5", "#458B00", "darkgray", "lightgray", "#EEAD0E", "#68228B", "#8B1A1A"),
                      breaks=c("ALDEx2", 'ALDEx2\n(scale informed)', "ANCOM-BC2", 'DESeq2', 
                               'edgeR', "MaAsLin 2", "MaAsLin 3", "MaAsLin 3\n(Spike-in)"))
  ggsave(paste0(figures_folder, 'SD2_spike_microbe_figure_sup.png'),
         plot = plot_out, width = 14, height = 8)
  
  
  plot_out <- ggplot(melted_df %>% filter(grepl('MaAsLin|ALDEx|ANCOM', tool) & 
                    spikeMicrobes %in% seq(0.1, 0.9, 0.2)), 
                     aes(x = get(param_name), y = value, fill = tool)) + 
      geom_boxplot(position = position_dodge(preserve = "single")) + 
      facet_wrap(~variable, scales = 'free', ncol = 4) + 
      theme_bw() + 
      xlab("Proportion of associations non-null") + 
      ylab('') + 
      theme(text=element_text(size=21),
            legend.position = 'right',
            strip.background = element_rect(fill = "gray95")) + 
      labs(fill = 'Model') + 
      scale_fill_manual(values=c("#104E8B", "#5B9BD5", "#458B00", "darkgray", "lightgray", "#EEAD0E", "#68228B", "#8B1A1A"),
                        breaks=c("ALDEx2", 'ALDEx2\n(scale informed)', "ANCOM-BC2", 'DESeq2', 
                                 'edgeR', "MaAsLin 2", "MaAsLin 3", "MaAsLin 3\n(Spike-in)"))
  ggsave(paste0(figures_folder, 'SD2_spike_microbe_figure.svg'),
         plot = plot_out, width = 16, height = 4)
}
SD2_spike_microbe_figure()