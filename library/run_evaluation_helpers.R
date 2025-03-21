library(dplyr)
library(pROC)

prepare_truth <- function(truth, generator) {
  if (generator == 'SD2') {
    outputs <- data.frame(taxon = truth$feature_spiked,
                          metadata = truth$metadata_datum,
                          effect_size = truth$effect_size / ifelse(truth$associated_property == 'abundance', log(2), 1), # SD2 is on exp rather than base 2
                          associations = truth$associated_property)
    return(outputs)
  }
  if (generator == 'SimSeq') {
    outputs <- data.frame(taxon = truth$feature_spiked,
                          metadata = truth$metadata_datum,
                          effect_size = truth$effect_size,
                          associations = truth$associated_property)
    return(outputs)
  }
  if (generator == 'ANCOM_BC_generator') {
    outputs <- data.frame(taxon = truth$feature_spiked,
                          metadata = truth$metadata_datum,
                          effect_size = pmax(log2(truth$effect_size), -10),
                          associations = truth$associated_property)
    outputs <- outputs[outputs$effect_size != 0,]
    outputs$effect_size <- -1 * outputs$effect_size
    outputs$associations = ifelse(outputs$effect_size == 10, "prevalence", "abundance")
    return(outputs)
  }
}

prepare_truth_groups <- function(truth, generator) {
  if (generator == 'SD2') {
    truth <- truth[grepl('[0-9]', substr(truth$metadata_datum, nchar(truth$metadata_datum), nchar(truth$metadata_datum))),]
    outputs <- data.frame(taxon = truth$feature_spiked,
                          metadata = truth$metadata_datum,
                          effect_size = truth$effect_size / ifelse(truth$associated_property == 'abundance', log(2), 1),
                          associations = truth$associated_property,
                          org_metadata = truth$metadata_datum)
    outputs$metadata <- gsub('[0-9]', '', outputs$metadata)
    outputs <- outputs[outputs$effect_size != 0,]
    outputs <- unique(outputs)
    return(outputs)
  }
}

allowed_errors <- c(NA)

prepare_associations_general <- function(associations, tool, generator = 'SD2', threshold = -1) {
    if (tool %in% c('ALDEx2', 'ALDEx2_scale', 'ANCOMBC', 'edgeR')) {
        if (tool == 'ANCOMBC') {
            associations <- associations[is.na(associations$error),]
            if (generator != 'ANCOM_BC_generator') {
                associations <- associations[associations$associations != 'prevalence',]
            }
            associations$effect_size <- associations$effect_size / log(2) # Inflate because transformation is originally base e
        }
        if (generator == 'SimSeq' | generator == 'ANCOM_BC_generator') {
            outputs <- data.frame(taxon = associations$taxon,
                                  metadata = gsub('[0-9]$', '', associations$metadata),
                                  effect_size = associations$effect_size,
                                  signif = associations$qval)
        } else {
            outputs <- data.frame(taxon = associations$taxon,
                                  metadata = gsub("([0-9])[0-9]$", "\\1", associations$metadata),
                                  effect_size = associations$effect_size,
                                  signif = associations$qval)
        }
        if (nrow(outputs) == 0) {
            outputs <- data.frame(matrix(nrow = 0, ncol = 4))
            colnames(outputs) <- c("taxon", "metadata", "effect_size", "signif")
        }
        return(outputs[rowSums(is.na(outputs)) == 0 & outputs$metadata != 'read_depth' & 
                           abs(outputs$effect_size) > threshold,])
    }
    if (tool %in% c('Maaslin2')) {
        outputs <- data.frame(taxon = associations$taxon,
                              metadata = associations$metadata,
                              effect_size = associations$effect_size,
                              signif = associations$qval)
        return(outputs[rowSums(is.na(outputs)) == 0 & outputs$metadata != 'read_depth' & 
                           abs(outputs$effect_size) > threshold,])
    }
    if (grepl('Maaslin3', tool)) {
        associations <- associations[associations$error %in% allowed_errors,]
        outputs <- data.frame(taxon = associations$taxon,
                              metadata = associations$metadata,
                              effect_size = associations$effect_size,
                              signif = associations$qval_joint)
        return(outputs[rowSums(is.na(outputs)) == 0 & outputs$metadata != 'read_depth' & 
                           abs(outputs$effect_size) > threshold,])
    }
    if (tool == 'DESeq2') {
        if (generator == 'SimSeq' | generator == 'ANCOM_BC_generator') {
            outputs <- data.frame(taxon = associations$taxon,
                                  metadata = sub("_.*", "", associations$metadata),
                                  effect_size = associations$effect_size,
                                  signif = associations$qval)
        } else {
            outputs <- data.frame(taxon = associations$taxon,
                                  metadata = sub("^(.*?_[0-9]+).*", "\\1", associations$metadata),
                                  effect_size = associations$effect_size,
                                  signif = associations$qval)
        }
        if (nrow(outputs) == 0) {
            outputs <- data.frame(matrix(nrow = 0, ncol = 4))
            colnames(outputs) <- c("taxon", "metadata", "effect_size", "signif")
        }
        return(outputs[rowSums(is.na(outputs)) == 0 & outputs$metadata != 'read_depth' & 
                           abs(outputs$effect_size) > threshold,])
    }
}

prepare_associations_abundance <- function(associations, tool, generator = 'SD2', threshold = -1) {
  if (grepl('Maaslin3', tool)) {
    associations <- associations[associations$error %in% allowed_errors,]
    associations <- associations[associations$associations == 'abundance',]
    outputs <- data.frame(taxon = associations$taxon,
                          metadata = associations$metadata,
                          effect_size = associations$effect_size,
                          signif = associations$qval)
    return(outputs[rowSums(is.na(outputs)) == 0 & outputs$metadata != 'read_depth' & 
                     abs(outputs$effect_size) > threshold,])
  } else if (tool == 'ANCOMBC') {
    associations <- associations[associations$associations == 'abundance',]
    return(prepare_associations_general(associations, tool, generator))
  }
  else {
    return(prepare_associations_general(associations, tool, generator))
  }
}

prepare_associations_maaslin3 <- function(associations, tool, allow_abun_to_prev = F, threshold = -1) {
  if (!grepl('Maaslin3', tool)) {
    stop("Only works with Maaslin3")
  }
    
  if (allow_abun_to_prev) {
      associations$error[grepl("Prevalence association possibly induced", 
                               associations$error)] <- NA
  }
    
  associations <- associations[associations$error %in% allowed_errors,]
  # associations <- associations[abs(associations$effect_size) > 1,]
  outputs <- data.frame(taxon = associations$taxon,
                        metadata = associations$metadata,
                        effect_size = associations$effect_size,
                        signif = associations$qval,
                        association = associations$associations)
  return(outputs[rowSums(is.na(outputs)) == 0 & outputs$metadata != 'read_depth' & 
                   abs(outputs$effect_size) > threshold,])
}

prepare_associations_maaslin3_group <- function(associations, tool) {
  if (!grepl('Maaslin3', tool)) {
    stop("Only works with Maaslin3")
  }
  associations <- associations[associations$error %in% allowed_errors,]
  outputs <- data.frame(taxon = associations$taxon,
                        metadata = associations$metadata,
                        effect_size = associations$effect_size,
                        signif = associations$qval,
                        association = associations$associations)
  return(outputs)
}

unweighted_precision_recall <- function(truth, associations, abundance, metadata, threshold = 0.1) {
  if (nrow(associations) == 0) {
    return(c(NA, NA))
  }
  # all_possible_pairs <- apply(expand.grid(rownames(abundance), colnames(metadata)),
  #                             1, paste0, collapse = '_')
  truth_match_vec <- unique(paste0(truth$taxon, '_', truth$metadata))
  associations_match_vec <- paste0(associations$taxon, '_', associations$metadata)
  associations_match_vec <- unique(associations_match_vec[associations$signif < threshold])
  
  precision = sum(associations_match_vec %in% truth_match_vec) / length(associations_match_vec)
  recall = sum(associations_match_vec %in% truth_match_vec) / length(truth_match_vec)
  
  return(c(precision, recall))
}

effect_size_fp <- function(truth, associations, abundance, metadata, threshold = 0.1) {
    if (nrow(associations) == 0) {
        return(c(NA, NA))
    }
    # all_possible_pairs <- apply(expand.grid(rownames(abundance), colnames(metadata)),
    #                             1, paste0, collapse = '_')
    truth_match_vec <- unique(paste0(truth$taxon, '_', truth$metadata))
    associations_match_vec <- paste0(associations$taxon, '_', associations$metadata)
    associations_match_vec <- unique(associations_match_vec[associations$signif < threshold])
    
    effect_size_fp <- mean(abs(associations$effect_size[(!paste0(associations$taxon, '_', associations$metadata) %in% truth_match_vec) & 
                                                       associations$signif < threshold]))
    
    if (length(effect_size_fp) == 0) {
        effect_size_fp <- NA
    }
    
    return(effect_size_fp)
}

effect_size_tp <- function(truth, associations, abundance, metadata, threshold = 0.1) {
    if (nrow(associations) == 0) {
        return(c(NA, NA))
    }
    # all_possible_pairs <- apply(expand.grid(rownames(abundance), colnames(metadata)),
    #                             1, paste0, collapse = '_')
    truth_match_vec <- unique(paste0(truth$taxon, '_', truth$metadata))
    associations_match_vec <- paste0(associations$taxon, '_', associations$metadata)
    associations_match_vec <- unique(associations_match_vec[associations$signif < threshold])
    
    effect_size_tp <- mean(abs(associations$effect_size[(paste0(associations$taxon, '_', associations$metadata) %in% truth_match_vec) & 
                                                            associations$signif < threshold]))
    
    if (length(effect_size_tp) == 0) {
        effect_size_tp <- NA
    }
    
    return(effect_size_tp)
}


TSSnorm = function(features) {
  # Convert to Matrix from Data Frame
  features_norm = as.matrix(features)
  X_mask <- ifelse(features > 0, 1, 0)
  dd <- colnames(features_norm)
  
  ##############
  # From vegan #
  ##############
  
  x <- as.matrix(features_norm)
  if (any(x < 0, na.rm = TRUE)) {
    k <- min(x, na.rm = TRUE)
    warning("input data contains negative entries: result may be non-sense")
  } else {
    k <- .Machine$double.eps
  }
  
  MARGIN <- 2
  
  tmp <- pmax(k, apply(x, MARGIN, sum, na.rm = TRUE))
  x <- sweep(x, MARGIN, tmp, "/")
  attr <- list(total = tmp, margin = MARGIN)
  if (any(is.nan(x))) 
    warning("result contains NaN, perhaps due to impossible mathematical\n
            operation\n")
  
  x <- ifelse(X_mask, x, 0)
  
  # Convert back to data frame
  features_TSS <- as.data.frame(x)
  
  # Rename the True Positive Features - Same Format as Before
  colnames(features_TSS) <- dd
  
  # Return
  return(features_TSS)
}

weighted_precision_recall <- function(truth, associations, abundance, metadata, threshold = 0.1, abundance_threshold = 0.001, prevalence_threshold = 0.1) {
  if (nrow(associations) == 0) {
    return(c(NA, NA))
  }
  
  abundance <- TSSnorm(abundance)
  allowed_taxa <- rownames(abundance)[rowMeans(abundance > abundance_threshold, na.rm = T) > prevalence_threshold]
  
  truth <- truth[truth$taxon %in% allowed_taxa,]
  associations <- associations[associations$taxon %in% allowed_taxa,]
  
  truth_match_vec <- unique(paste0(truth$taxon, '_', truth$metadata))
  associations_match_vec <- paste0(associations$taxon, '_', associations$metadata)
  associations_match_vec <- unique(associations_match_vec[associations$signif < threshold])
  
  precision = sum(associations_match_vec %in% truth_match_vec) / length(associations_match_vec)
  recall = sum(associations_match_vec %in% truth_match_vec) / length(truth_match_vec)
  
  return(c(precision, recall))
}

pval_auc <- function(truth, associations, abundance, metadata) {
  if (nrow(associations) == 0) {
    return(c(NA))
  }
  
  if (nrow(associations) == 0) {
    return(NA)
  }
  truth_match_vec <- unique(paste0(truth$taxon, '_', truth$metadata))
  associations_match_vec <- paste0(associations$taxon, '_', associations$metadata)
  associations$correct <- associations_match_vec %in% truth_match_vec
  if (length(unique(associations$correct)) < 2) {
    return(NA)
  }
  associations <- unique(associations[,c('taxon', 'metadata', 'signif', 'correct')])
  roc_curve <- roc(associations$correct, associations$signif, direction = "<", levels = c(TRUE, FALSE))
  return(c(auc(roc_curve)))
}

weighted_pval_auc <- function(truth, associations, abundance, metadata, abundance_threshold = 0.001, prevalence_threshold = 0.1) {
  if (nrow(associations) == 0) {
    return(NA)
  }
  
  abundance <- TSSnorm(abundance)
  allowed_taxa <- rownames(abundance)[rowMeans(abundance > abundance_threshold, na.rm = T) > prevalence_threshold]
  
  truth <- truth[truth$taxon %in% allowed_taxa,]
  associations <- associations[associations$taxon %in% allowed_taxa,]
  if (nrow(associations) == 0) {
    return(NA)
  }
  
  truth_match_vec <- unique(paste0(truth$taxon, '_', truth$metadata))
  associations_match_vec <- paste0(associations$taxon, '_', associations$metadata)
  associations$correct <- associations_match_vec %in% truth_match_vec
  if (length(unique(associations$correct)) < 2) {
    return(NA)
  }
  associations <- unique(associations[,c('taxon', 'metadata', 'signif', 'correct')])
  roc_curve <- roc(associations$correct, associations$signif, direction = "<", levels = c(TRUE, FALSE))
  return(c(auc(roc_curve)))
}

effect_size_error <- function(truth, associations, threshold = 0.1) {
  if (nrow(associations) == 0) {
    return(c(NA, NA))
  }
  
  truth <- truth[truth$associations == 'abundance',]
  
  truth_match_vec <- unique(paste0(truth$taxon, '_', truth$metadata))
  associations_match_vec <- paste0(associations$taxon, '_', associations$metadata)
  associations_match_vec <- associations_match_vec[associations$signif < threshold]
  associations <- associations[associations$signif < threshold,]
  
  associations <- associations[associations_match_vec %in% truth_match_vec, -4]
  truth <- truth[truth_match_vec %in% associations_match_vec, -4]

  colnames(truth) <- c(colnames(truth)[-3], "truth")
  colnames(associations) <- c(colnames(associations)[-3], "association")
  
  if (nrow(truth) > 0 & nrow(associations) > 0) {
    joined_df <- full_join(truth, associations, by=c('taxon', 'metadata'))
  } else {
    return(c(NA, NA))
  }

  first_return <- abs(joined_df$truth - joined_df$association) / abs(joined_df$truth)
  second_return <- (abs(joined_df$association) - abs(joined_df$truth)) / abs(joined_df$truth)
  for_return <- c(mean(first_return[!is.infinite(first_return)]),
                  mean(second_return[!is.infinite(second_return)]))
  
  return(for_return)
}

effect_size_correlation <- function(truth, associations, threshold = 0.1) {
  truth <- truth[truth$associations == 'abundance',]
  
  if (nrow(associations) == 0) {
    return(NA)
  }
  
  truth_match_vec <- paste0(truth$taxon, '_', truth$metadata)
  associations_match_vec <- paste0(associations$taxon, '_', associations$metadata)
  associations_match_vec <- associations_match_vec[associations$signif < threshold]
  associations <- associations[associations$signif < threshold,]
  
  associations <- associations[associations_match_vec %in% truth_match_vec, -4]
  truth <- truth[truth_match_vec %in% associations_match_vec, -4]
  
  colnames(truth) <- c(colnames(truth)[-3], "truth")
  colnames(associations) <- c(colnames(associations)[-3], "association")
  
  if (nrow(truth) > 0 & nrow(associations) > 0) {
    joined_df <- full_join(truth, associations, by=c('taxon', 'metadata'))
  } else {
    return(NA)
  }
  
  if (nrow(joined_df) < 3) {
    return(NA)
  }
  
  joined_df = joined_df %>%
    dplyr::group_by(metadata) %>%
    dplyr::summarise(correlation = cor(truth, association, use = 'pairwise.complete.obs', method = 'spearman'))
  
  return(mean(joined_df$correlation, na.rm=T))
}

unweighted_precision_recall_maaslin3 <- function(truth, associations, abundance, metadata, threshold = 0.1) {
  if (nrow(associations) == 0) {
    return(c(NA, NA))
  }
  
  truth_match_vec <- unique(paste0(truth$taxon, '_', truth$metadata, '_', truth$associations))
  associations_match_vec <- paste0(associations$taxon, '_', associations$metadata, '_', associations$association)
  associations_match_vec <- unique(associations_match_vec[associations$signif < threshold])
  
  precision = sum(associations_match_vec %in% truth_match_vec) / length(associations_match_vec)
  recall = sum(associations_match_vec %in% truth_match_vec) / length(truth_match_vec)
  
  return(c(precision, recall))
}

effect_size_fp_maaslin3 <- function(truth, associations, abundance, metadata, threshold = 0.1) {
    if (nrow(associations) == 0) {
        return(c(NA, NA))
    }
    # all_possible_pairs <- apply(expand.grid(rownames(abundance), colnames(metadata)),
    #                             1, paste0, collapse = '_')
    truth_match_vec <- unique(paste0(truth$taxon, '_', truth$metadata, '_', truth$associations))
    associations_match_vec <- paste0(associations$taxon, '_', associations$metadata, '_', associations$association)

    effect_size_fp <- mean(abs(associations$effect_size[(!associations_match_vec %in% truth_match_vec) & 
                                                       associations$signif < threshold]))
    
    if (length(effect_size_fp) == 0) {
        effect_size_fp <- NA
    }
    
    
    return(effect_size_fp)
}

effect_size_tp_maaslin3 <- function(truth, associations, abundance, metadata, threshold = 0.1) {
    if (nrow(associations) == 0) {
        return(c(NA, NA))
    }
    # all_possible_pairs <- apply(expand.grid(rownames(abundance), colnames(metadata)),
    #                             1, paste0, collapse = '_')
    truth_match_vec <- unique(paste0(truth$taxon, '_', truth$metadata, '_', truth$associations))
    associations_match_vec <- paste0(associations$taxon, '_', associations$metadata, '_', associations$association)
    
    effect_size_tp <- mean(abs(associations$effect_size[(associations_match_vec %in% truth_match_vec) & 
                                                            associations$signif < threshold]))
    
    if (length(effect_size_tp) == 0) {
        effect_size_tp <- NA
    }
    
    
    return(effect_size_tp)
}

weighted_precision_recall_maaslin3 <- function(truth, associations, abundance, metadata, threshold = 0.1, abundance_threshold = 0.001, prevalence_threshold = 0.1) {
  if (nrow(associations) == 0) {
    return(c(NA, NA))
  }
  
  
  abundance <- TSSnorm(abundance)
  allowed_taxa <- rownames(abundance)[rowMeans(abundance > abundance_threshold, na.rm = T) > prevalence_threshold]

  truth <- truth[truth$taxon %in% allowed_taxa,]
  associations <- associations[associations$taxon %in% allowed_taxa,]
  
  truth_match_vec <- unique(paste0(truth$taxon, '_', truth$metadata, '_', truth$associations))
  associations_match_vec <- paste0(associations$taxon, '_', associations$metadata, '_', associations$association)
  associations_match_vec <- unique(associations_match_vec[associations$signif < threshold])
  
  precision = sum(associations_match_vec %in% truth_match_vec) / length(associations_match_vec)
  recall = sum(associations_match_vec %in% truth_match_vec) / length(truth_match_vec)
  
  return(c(precision, recall))
}

pval_auc_maaslin3 <- function(truth, associations, abundance, metadata) {
  if (nrow(associations) == 0) {
    return(NA)
  }
  
  if (nrow(associations) == 0) {
    return(NA)
  }
  truth_match_vec <- unique(paste0(truth$taxon, '_', truth$metadata, '_', truth$associations))
  associations_match_vec <- paste0(associations$taxon, '_', associations$metadata, '_', associations$association)
  associations$correct <- associations_match_vec %in% truth_match_vec
  if (length(unique(associations$correct)) < 2) {
    return(NA)
  }
  roc_curve <- roc(associations$correct, associations$signif, direction = "<", levels = c(TRUE, FALSE))
  return(c(auc(roc_curve)))
}

weighted_pval_auc_maaslin3 <- function(truth, associations, abundance, metadata, abundance_threshold = 0.001, prevalence_threshold = 0.1) {
  if (nrow(associations) == 0) {
    return(NA)
  }
  
  abundance <- TSSnorm(abundance)
  allowed_taxa <- rownames(abundance)[rowMeans(abundance > abundance_threshold, na.rm = T) > prevalence_threshold]
  
  truth <- truth[truth$taxon %in% allowed_taxa,]
  associations <- associations[associations$taxon %in% allowed_taxa,]
  if (nrow(associations) == 0) {
    return(NA)
  }
  
  truth_match_vec <- unique(paste0(truth$taxon, '_', truth$metadata, '_', truth$associations))
  associations_match_vec <- paste0(associations$taxon, '_', associations$metadata, '_', associations$association)
  associations$correct <- associations_match_vec %in% truth_match_vec
  if (length(unique(associations$correct)) < 2) {
    return(NA)
  }
  roc_curve <- roc(associations$correct, associations$signif, direction = "<", levels = c(TRUE, FALSE))
  return(c(auc(roc_curve)))
}

effect_size_error_maaslin3 <- function(truth, associations, threshold = 0.1) {
  if (nrow(associations) == 0) {
    return(c(NA, NA))
  }
  
  truth_match_vec <- paste0(truth$taxon, '_', truth$metadata, '_', truth$associations)
  associations_match_vec <- paste0(associations$taxon, '_', associations$metadata, '_', associations$association)
  associations_match_vec <- associations_match_vec[associations$signif < threshold]
  associations <- associations[associations$signif < threshold,]
  
  associations <- associations[associations_match_vec %in% truth_match_vec, -c(4)]
  truth <- truth[truth_match_vec %in% associations_match_vec,]
  
  colnames(truth) <- c(colnames(truth)[-c(3,4)], "truth", "associations")
  colnames(associations) <- c(colnames(associations)[-c(3,4)], "association", "associations")
  
  if (nrow(truth) > 0 & nrow(associations) > 0) {
    joined_df <- full_join(truth, associations, by=c('taxon', 'metadata', 'associations'))
  } else {
    return(c(NA, NA))
  }
  
  return_val <- c(mean(abs(joined_df$truth - joined_df$association) / abs(joined_df$truth)),
                  mean((abs(joined_df$association) - abs(joined_df$truth)) / abs(joined_df$truth)))
  return(return_val)
}

effect_size_correlation_maaslin3 <- function(truth, associations, threshold = 0.1) {
  if (nrow(associations) == 0) {
    return(NA)
  }
  
  truth_match_vec <- paste0(truth$taxon, '_', truth$metadata, '_', truth$associations)
  associations_match_vec <- paste0(associations$taxon, '_', associations$metadata, '_', associations$association)
  associations_match_vec <- associations_match_vec[associations$signif < threshold]
  associations <- associations[associations$signif < threshold,]
  
  associations <- associations[associations_match_vec %in% truth_match_vec, -c(4)]
  truth <- truth[truth_match_vec %in% associations_match_vec,]
  
  colnames(truth) <- c(colnames(truth)[-c(3,4)], "truth", "associations")
  colnames(associations) <- c(colnames(associations)[-c(3,4)], "association", "associations")
  
  if (nrow(truth) > 0 & nrow(associations) > 0) {
    joined_df <- full_join(truth, associations, by=c('taxon', 'metadata', 'associations'))
  } else {
    return(NA)
  }
  
  if (nrow(joined_df) < 3) {
    return(NA)
  }
  
  joined_df = joined_df %>%
    dplyr::group_by(metadata) %>%
    dplyr::summarise(correlation = cor(truth, association, use = 'pairwise.complete.obs', method = 'spearman'))
  
  return(mean(joined_df$correlation, na.rm=T))
}

issue_prop <- function(associations, tool) {
  if (tool %in% c('ALDEx2', 'ALDEx2_scale', 'ANCOMBC', 'DESeq2', 'edgeR')) {
    outputs <- data.frame(taxon = associations$taxon,
                          metadata = associations$metadata,
                          effect_size = associations$effect_size,
                          signif = associations$qval)
    return(mean(is.na(outputs$signif)))
  }
  if (tool %in% c('Maaslin2')) {
    outputs <- data.frame(taxon = associations$taxon,
                          metadata = associations$metadata,
                          effect_size = associations$effect_size,
                          signif = associations$qval)
    return(mean(is.na(outputs$signif)))
  }
  if (grepl('Maaslin3', tool)) {
    associations <- associations[associations$error %in% allowed_errors,]
    outputs <- data.frame(taxon = associations$taxon,
                          metadata = associations$metadata,
                          effect_size = associations$effect_size,
                          signif = associations$qval_joint)
    return(mean(!is.na(associations$error) | is.na(associations$pval)))
  }
}

effect_size_mean_diff <- function(truth, associations) {
  if (nrow(associations) == 0) {
    return(c(NA))
  }
  
  truth <- truth[truth$associations == 'abundance',]
  
  truth_match_vec <- unique(paste0(truth$taxon, '_', truth$metadata))
  associations_match_vec <- paste0(associations$taxon, '_', associations$metadata)

  associations <- associations[associations_match_vec %in% truth_match_vec, -4]
  truth <- truth[truth_match_vec %in% associations_match_vec, -4]
  
  colnames(truth) <- c(colnames(truth)[-3], "truth")
  colnames(associations) <- c(colnames(associations)[-3], "association")
  
  if (nrow(truth) > 0 & nrow(associations) > 0) {
    joined_df <- full_join(truth, associations, by=c('taxon', 'metadata'))
  } else {
    return(c(NA, NA))
  }
  
  first_return <- mean(joined_df$association - joined_df$truth, na.rm=T)

  return(first_return)
}
