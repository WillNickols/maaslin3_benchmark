#!/usr/bin/env Rscript

remove(list = ls())
package_vec = c("devtools", "pkgmaker", "optparse", "ALDEx2",
                "parallel", "stringi", "doParallel", "plyr", "tidyr")
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
inputDirectory <- file.path(workingDirectory, 'Input', 'general_evaluations', generator)
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
outputDirectory <- file.path(workingDirectory, 'Output', 'general_evaluations', generator)
outputString <- paste(inputString, "ALDEx2", sep="_")
this_output_folder <- file.path(outputDirectory, outputString)

outputs_already_exist <- TRUE
if (!dir.exists(this_output_folder)) {
  outputs_already_exist <- FALSE
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
  
  f <- foreach(i = reps, .packages = c("MASS", "stringi", 'dplyr', 'plyr', 'ALDEx2')) %dopar% {
    metadata <- read.csv(paste0(this_params_folder, "/metadata_", i, ".tsv"), sep = "\t")
    abundance <- read.csv(paste0(this_params_folder, "/abundance_", i, ".tsv"), sep = "\t")
    truth <- read.csv(paste0(this_params_folder, "/truth_", i, ".tsv"), sep = "\t")
    
    for (col in names(metadata)[names(metadata) != "ID"]) {
      if (all(metadata[[col]] %in% c(0, 1))) {
        metadata[[col]] <- as.character(metadata[[col]])
      } else {
        metadata[[col]] <- as.numeric(metadata[[col]])
      }
    }
    
    if ("ID" %in% colnames(metadata)) {
      ID <- metadata$ID
    } else {
      ID <- rownames(metadata)
    }
    
    sink('/dev/null')
    if ('ID' %in% colnames(metadata)) {
      if (length(unique(metadata$ID)) == length(metadata$ID)) { # No random effect
        mm <- model.matrix(formula(paste0("~ ", 
                                          paste0(colnames(metadata)[colnames(metadata) != "ID"], 
                                                 collapse = " + "), 
                                          collapse = "")),metadata)
      } else { # Random effect
        mm <- model.matrix(formula(paste0("~ ", 
                                          paste0(colnames(metadata), 
                                                 collapse = " + "), 
                                          collapse = "")),metadata)
      }
    } else { # No random effect
      mm <- model.matrix(formula(paste0("~ ", 
                                        paste0(colnames(metadata)[colnames(metadata) != "ID"], 
                                               collapse = " + "), 
                                        collapse = "")),metadata)
    }
    
    aldex_clr_out <- aldex.clr(abundance, mm, denom="all")
    glm.test <- aldex.glm(aldex_clr_out)
    
    glm.test <- glm.test[,grepl("Est$|pval$", colnames(glm.test))]
    
    if ('ID' %in% colnames(metadata)) {
      glm.test <- glm.test[,!grepl(paste0(c("Intercept", unique(metadata$ID)), collapse = "|"), colnames(glm.test))]
    } else {
      glm.test <- glm.test[,!grepl("Intercept", colnames(glm.test))]
    }
    glm.test$Feature <- rownames(glm.test)
    glm.test <- reshape2::melt(glm.test, id.vars = c("Feature"))
    glm.test$metric <- gsub(".*\\:", "", glm.test$variable)
    glm.test$variable <- gsub("\\:.*", "", glm.test$variable)
    glm.test <- reshape2::dcast(formula = Feature + variable ~ metric, glm.test)
    sink()
    
    outputs <- data.frame(taxon = glm.test$Feature,
                          metadata = glm.test$variable,
                          effect_size = glm.test$Est,
                          pval = glm.test$pval,
                          qval = p.adjust(glm.test$pval, method = "BH"),
                          associations = "abundance")
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
