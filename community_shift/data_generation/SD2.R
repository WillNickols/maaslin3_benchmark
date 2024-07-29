#!/usr/bin/env Rscript

remove(list = ls())
package_vec = c("devtools", "SparseDOSSA2", "pkgmaker", "optparse", 
                "parallel", "stringi", "doParallel", "dplyr")
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

# Load Utility Functions
pkgmaker::source_files(paste(workingDirectory, '/library/SD2_helper_functions.R', sep=''))

if (RandomEffect==TRUE){
  inputSubString<-'RandomEffect'
} else {
  inputSubString<-'noRandomEffect'
}

options("scipen"=10)
inputString<-paste(inputSubString, metadataType, nSubjects, nPerSubject, nMicrobes, spikeMicrobes, nMetadata, effectSize, effectPos, readDepth, sep='_')
options("scipen"=5)

# Create Input Directory
inputDirectory <- file.path(workingDirectory, 'Input', 'community_shift', 'SD2')
this_params_folder <- file.path(inputDirectory, inputString)

outputs_already_exist <- TRUE
if (!dir.exists(this_params_folder)) {
  outputs_already_exist <- FALSE
} else {
  for (i in 1:length(nIterations)) {
    if (!file.exists(paste0(this_params_folder, "/metadata_", i, ".tsv")) |
        !file.exists(paste0(this_params_folder, "/abundance_", i, ".tsv")) |
        !file.exists(paste0(this_params_folder, "/truth_", i, ".tsv")) |
        !file.exists(paste0(this_params_folder, "/unscaled_", i, ".tsv"))) {
      outputs_already_exist <- FALSE
    }
  }
}

if (!outputs_already_exist){
  dir.create(this_params_folder, recursive = T)
  
  simlist<-trigger_sparseDOSSA2_Simulator_unscaled(RandomEffect=RandomEffect,
                                          metadataType=metadataType,
                                          nSubjects=nSubjects,
                                          nPerSubject=nPerSubject,
                                          nMicrobes=nMicrobes,
                                          spikeMicrobes = spikeMicrobes, 
                                          nMetadata = nMetadata, 
                                          effectSize = effectSize,
                                          effectPos = effectPos,
                                          readDepth = readDepth,
                                          noParallel = noParallel,
                                          nIterations=nIterations,
                                          rSeed=rSeed,
                                          nCores=nCores)
  
  for (i in 1:length(simlist)) {
    current_sim <- simlist[[i]]
    write.table(current_sim$metadata, 
                file = paste0(this_params_folder, "/metadata_", i, ".tsv"), 
                sep = '\t')
    write.table(current_sim$features, 
                file = paste0(this_params_folder, "/abundance_", i, ".tsv"), 
                sep = '\t')
    write.table(current_sim$truth, 
                file = paste0(this_params_folder, "/truth_", i, ".tsv"), 
                sep = '\t',
                row.names = F)
    write.table(current_sim$abundance_unscaled, 
                file = paste0(this_params_folder, "/unscaled_", i, ".tsv"), 
                sep = '\t')
  }
} else {
  print("Parameter combination already exists")
}

#######################################
# Delete Temporary sparseDOSSA Files  #
#######################################

if (file.exists("SyntheticMicrobiome-Counts.pcl")) file.remove("SyntheticMicrobiome-Counts.pcl")
if (file.exists("SyntheticMicrobiome.pcl")) file.remove("SyntheticMicrobiome.pcl")
if (file.exists("SyntheticMicrobiomeParameterFile.txt")) file.remove("SyntheticMicrobiomeParameterFile.txt")
