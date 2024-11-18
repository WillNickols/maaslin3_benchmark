#!/usr/bin/env Rscript

remove(list = ls())
package_vec = c("devtools", "pkgmaker", "optparse", 
                "parallel", "stringi", "doParallel", "plyr", "tidyr", "maaslin3")
invisible(suppressPackageStartupMessages(lapply(package_vec, require, character.only = TRUE)))

# Command Line Usage
option_list = list(
  make_option(
    c("-l", "--RandomEffect"), default=FALSE, # l stands for longitudinal design
    action = "store_true"),
  make_option(
    c("-d", "--metadataType"), # d stands for design matrix
    type = "character"),
  make_option(
    c("-n", "--nSubjects"), # n stands for number of subjects
    type = "integer"),
  make_option(
    c("-b", "--nPerSubject"), 
    type = "integer",
    default=1), # b stands for block size
  make_option(
    c("-f", "--nMicrobes"), # f stands for number of features
    type = "integer"),
  make_option(
    c("-s", "--spikeMicrobes"), # s stands for spike
    type = "numeric"),
  make_option(
    c("-m", "--nMetadata"), # m stands for number of metadata
    type = "integer"),
  make_option(
    c("-e", "--effectSize"), # e stands for effect size
    type = "numeric"),
  make_option(
    c("--effectPos"), default=0.5,
    type = "numeric"),
  make_option(
    c("-q", "--readDepth"), default=500000, # q stands for quantity and quality of sequencing depth or library size
    type = "integer"),
  make_option(
    c("-g", "--noParallel"), default=FALSE, # g stands for grid
    action = "store_true"),
  make_option(
    c("-t", "--nIterations"), default=100, # t stands for how many times an experiment is repeated
    type = "integer"),
  make_option(
    c("-r", "--rSeed"), default=1234, # r stands for reproducibility index (random seed)
    type = "integer"),
  make_option(
    c("-c", "--nCores"), default=4, # c stands for how many cores to be used in the analysis
    type = "integer"),
  make_option(
    c("-w", "--workingDirectory"), # w stands for working directory
    type = "character"),
  make_option(
    c("--generator"),
    type = "character"))
opt <- parse_args(OptionParser(option_list=option_list), positional_arguments = TRUE)

# Extract Parameters
RandomEffect<-opt$options$RandomEffect # High-level parameter
metadataType <- opt$options$metadataType # High-level parameter
nSubjects <- opt$options$nSubjects # Low-level parameter
nPerSubject<-opt$options$nPerSubject # Low-level parameter
nMicrobes <- opt$options$nMicrobes # Low-level parameter
spikeMicrobes <- opt$options$spikeMicrobes # Low-level parameter
nMetadata<- opt$options$nMetadata # Low-level parameter
effectSize<- opt$options$effectSize # Low-level parameter
effectPos<- opt$options$effectPos # Low-level parameter
readDepth<- opt$options$readDepth # Default parameter
noParallel<-opt$options$noParallel # Default parameter
nIterations<- opt$options$nIterations # Default parameter
rSeed<- opt$options$rSeed # Default parameter
nCores<- opt$options$nCores # Default parameter
workingDirectory <- opt$options$workingDirectory # Default parameter
generator <- opt$options$generator

if (RandomEffect==TRUE){
  inputSubString<-'RandomEffect'
} else {
  inputSubString<-'noRandomEffect'
}

options("scipen"=10)
inputString<-paste(inputSubString, metadataType, nSubjects, nPerSubject, nMicrobes, spikeMicrobes, nMetadata, effectSize, effectPos, readDepth, sep='_')
options("scipen"=5)

# Create Input Directory
inputDirectory <- file.path(workingDirectory, 'Input', 'deep_sequencing', generator)
this_params_folder <- file.path(inputDirectory, inputString)

outputs_already_exist <- TRUE
if (!dir.exists(this_params_folder)) {
  outputs_already_exist <- FALSE
} else {
  for (i in 1:length(nIterations)) {
    if (!file.exists(paste0(this_params_folder, "/metadata_", i, ".tsv")) |
        !file.exists(paste0(this_params_folder, "/abundance_", i, ".tsv")) |
        !file.exists(paste0(this_params_folder, "/truth_", i, ".tsv"))) {
      outputs_already_exist <- FALSE
    }
  }
}

# Load Dataset
if (!outputs_already_exist) stop('The input file does not exist. Generate the dataset first.') 

# How Many Datasets to Run In A List 
reps<-1:nIterations

# Create Output Directory
outputDirectory <- file.path(workingDirectory, 'Output', 'deep_sequencing', generator)
outputString <- paste(inputString, "Maaslin3CompAdjust", sep="_")
this_output_folder <- file.path(outputDirectory, outputString)

outputs_already_exist <- TRUE
if (!dir.exists(this_output_folder)) {
  outputs_already_exist <- FALSE
  dir.create(this_output_folder, recursive = T, showWarnings = F)
} else {
  for (i in 1:length(nIterations)) {
    if (!file.exists(paste0(this_output_folder, "/associations_", i, ".tsv"))) {
      outputs_already_exist <- FALSE
    }
  }
}

# Run the Model Only if the Output Does Not Exist
if (!outputs_already_exist){
  no_cores <- nCores
  cl <- makeCluster(no_cores)
  registerDoParallel(cl)
  
  f <- foreach(i = reps, .packages = c("MASS", "stringi", 'dplyr', 'plyr', 'lme4')) %dopar% {
    metadata <- read.csv(paste0(this_params_folder, "/metadata_", i, ".tsv"), sep = "\t")
    abundance <- read.csv(paste0(this_params_folder, "/abundance_", i, ".tsv"), sep = "\t")
    truth <- read.csv(paste0(this_params_folder, "/truth_", i, ".tsv"), sep = "\t")
    
    # Remove spike-in
    abundance <- abundance[1:(nrow(abundance) - 1),]
    
    # Get read depths
    abundance <- abundance[apply(abundance, 1, var) != 0,]
    read_depths <- data.frame(sample = names(colSums(abundance)),
                              read_depth = colSums(abundance))
    metadata$sample <- rownames(metadata)
    metadata <- left_join(metadata, read_depths, by = 'sample')
    rownames(metadata) <- metadata$sample
    metadata$sample <- NULL
    
    # Convert to relative abundance
    abundance <- t(t(abundance) / colSums(abundance))
    
    for (column in colnames(metadata)) {
      if (length(unique(metadata[,column])) == 2) {
        metadata[,column] <- as.factor(metadata[,column])
      }
    }
    
    if ("ID" %in% colnames(metadata)) {
      ID <- metadata$ID
    } else {
      ID <- rownames(metadata)
    }
    
    tmp_fit_out <- paste0(this_output_folder, "/tmp_out_", i)
    dir.create(tmp_fit_out, recursive = T)
    
    sink('/dev/null')
    set.seed(1)
    if(length(ID)==length(unique(ID))){
      fit_out <- maaslin3::maaslin3(input_data = abundance, 
                         input_metadata = metadata, 
                         output = tmp_fit_out, 
                         normalization = 'TSS', 
                         transform = 'LOG',
                         fixed_effects = colnames(metadata)[!colnames(metadata) %in% c("ID", "read_depth")], 
                         median_comparison_abundance = T, 
                         median_comparison_prevalence = F,
                         subtract_median = T,
                         plot_summary_plot = F, 
                         plot_associations = F, 
                         max_significance = 0.1)
    } else{
      fit_out <- maaslin3::maaslin3(input_data = abundance, 
                         input_metadata = metadata, 
                         output = tmp_fit_out, 
                         normalization = 'TSS', 
                         transform = 'LOG',
                         fixed_effects = colnames(metadata)[!colnames(metadata) %in% c("ID", "read_depth")],
                         random_effects = "ID", 
                         median_comparison_abundance = T, 
                         median_comparison_prevalence = F,
                         subtract_median = T,
                         plot_summary_plot = F, 
                         plot_associations = F, 
                         max_significance = 0.1)
    }
    sink()
    
    unlink(tmp_fit_out, recursive = T)
    
    fit_out_lm <- fit_out$fit_data_abundance$results
    fit_out_lm <- fit_out_lm[c("feature", "metadata", "coef", "pval_individual", "error", "qval_individual", "pval_joint", "qval_joint")]
    fit_out_lm$association <- "abundance"
    
    fit_out_binary <- fit_out$fit_data_prevalence$results
    fit_out_binary <- fit_out_binary[c("feature", "metadata", "coef", "pval_individual", "error", "qval_individual", "pval_joint", "qval_joint")]
    fit_out_binary$association <- "prevalence"
    
    fit_out <- full_join(fit_out_lm, fit_out_binary, by = colnames(fit_out_lm))
    
    outputs <- data.frame(taxon = fit_out$feature,
                          metadata = fit_out$metadata,
                          effect_size = fit_out$coef,
                          pval = fit_out$pval_individual,
                          qval = fit_out$qval_individual,
                          pval_joint = fit_out$pval_joint,
                          qval_joint = fit_out$qval_joint,
                          error = fit_out$error,
                          associations = fit_out$association)
    return(outputs)
  }
  
  for (i in 1:nIterations) {
    write.table(f[[i]],
                file = paste0(this_output_folder, "/associations_", i, ".tsv"),
                sep = '\t',
                row.names = F)
  }
  
  # Stop the Cluster 
  stopCluster(cl)
}

