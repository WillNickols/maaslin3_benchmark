#!/usr/bin/env Rscript

remove(list = ls())
package_vec = c("devtools", "pkgmaker", "optparse", "ANCOMBC", "TreeSummarizedExperiment",
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
inputDirectory <- file.path(workingDirectory, 'Input', 'unscaled', generator)
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
outputDirectory <- file.path(workingDirectory, 'Output', 'unscaled', generator)
outputString <- paste(inputString, "ANCOMBC", sep="_")
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
  
  f <- foreach(i = reps, .packages = c("MASS", "stringi", 'dplyr', 'plyr', 'ANCOMBC', 'TreeSummarizedExperiment')) %dopar% {
    metadata <- read.csv(paste0(this_params_folder, "/metadata_", i, ".tsv"), sep = "\t")
    abundance <- read.csv(paste0(this_params_folder, "/abundance_", i, ".tsv"), sep = "\t")
    truth <- read.csv(paste0(this_params_folder, "/truth_", i, ".tsv"), sep = "\t")
    
    abundance <- abundance[apply(abundance, 1, var) != 0,]
    
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
        fix_formula <- paste0(colnames(metadata)[colnames(metadata) != "ID"], collapse = " + ")
        random_formula <- NULL
      } else { # Random effect
        fix_formula <- paste0(colnames(metadata)[colnames(metadata) != "ID"], collapse = " + ")
        random_formula <- "(1|ID)"
      }
    } else { # No random effect
      fix_formula <- paste0(colnames(metadata)[colnames(metadata) != "ID"], collapse = " + ")
      random_formula <- NULL
    }
    
    binary_cols <- colnames(metadata)[apply(metadata, 2, function(x){length(unique(x)) == 2})]
    
    ancombc_out <- tryCatch({
      assays = S4Vectors::SimpleList(counts = as.matrix(abundance))
      smd = S4Vectors::DataFrame(metadata)
      both_tse <- TreeSummarizedExperiment(
        assays = assays,
        colData = smd)
      ancombc_out <- ancombc2(both_tse, fix_formula = fix_formula,
                              rand_formula = random_formula, p_adj_method = "BH", alpha = 0.1, struc_zero = T,
                              group = binary_cols[1])
      ancombc_out
    }, error = function(err) {
      if (grepl('nPlease remove', err)) {
        tryCatch({
          features_to_drop <- unlist(strsplit(gsub("\nPlease remove.*", "", gsub(".*taxa\\:\n", "", err)), ", "))
          abundance <- abundance[!(rownames(abundance) %in% features_to_drop),]
          assays = S4Vectors::SimpleList(counts = as.matrix(abundance))
          smd = S4Vectors::DataFrame(metadata)
          both_tse <- TreeSummarizedExperiment(
            assays = assays,
            colData = smd)
          ancombc_out <- ancombc2(both_tse, fix_formula = fix_formula,
                                  rand_formula = random_formula, p_adj_method = "BH", alpha = 0.1, struc_zero = T,
                                  group = binary_cols[1])
          return(ancombc_out)
        }, error = function(err) {
          return(NULL)
        })
      } else {
        return(NULL)
      }
    })
    
    if (is.null(ancombc_out)) {
      outputs <- data.frame(taxon = character(0),
                            metadata = character(0),
                            effect_size = character(0),
                            pval = character(0),
                            qval = character(0),
                            associations = character(0))
    }
    
    outputs <- tryCatch({
      glm.test <- ancombc_out$res
      glm.test <- glm.test[,grepl("^taxon|^lfc_|^p_|^passed_ss_", colnames(glm.test))]
      glm.test <- glm.test[,!grepl("Intercept", colnames(glm.test))]
      
      glm.test <- reshape2::melt(glm.test, id.vars = c("taxon"))
      glm.test$metric <- gsub("_.*", "", glm.test$variable)
      glm.test$variable <- gsub("^[^_]*_|^passed_ss_", "", glm.test$variable)
      glm.test <- reshape2::dcast(formula = taxon + variable ~ metric, glm.test)
      
      struc_zeros <- ancombc_out$zero_ind[ancombc_out$zero_ind[,2] != ancombc_out$zero_ind[,3],]
      glm.test <- rbind(glm.test, data.frame(taxon = struc_zeros$taxon,
                                             variable = binary_cols[1],
                                             lfc = ifelse(struc_zeros[,3] == F, Inf, -Inf),
                                             p = 0,
                                             passed = 1))
      
      outputs <- data.frame(taxon = glm.test$taxon,
                            metadata = glm.test$variable,
                            effect_size = glm.test$lfc,
                            pval = glm.test$p,
                            qval = p.adjust(glm.test$p, method = "BH"),
                            error = ifelse(glm.test$passed == 1, NA, "sensitivity failed"),
                            associations = ifelse(is.infinite(glm.test$lfc), "prevalence", "abundance"))
      
      outputs
    }, error = function(err) {
      return(NULL)
    })
    
    sink()
    
    if (is.null(ancombc_out) | is.null(outputs)) {
      return(data.frame(taxon = character(0),
                        metadata = character(0),
                        effect_size = character(0),
                        pval = character(0),
                        qval = character(0),
                        error = character(0),
                        associations = character(0)))
    }
    
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
