#!/usr/bin/env Rscript

remove(list = ls())
package_vec = c("devtools", "pkgmaker", "optparse", "ALDEx2",
                "parallel", "stringi", "doParallel", "plyr", "tidyr", "dplyr")
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
        c("--depthConfound"), default=FALSE,
        type = "logical"),
    make_option(
        c("--propAbun"), default=0.5,
        type = "numeric"),
    make_option(
        c("--zeroInflate"), default=TRUE,
        type = "logical"),
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
depthConfound<- opt$options$depthConfound
propAbun<- opt$options$propAbun
zeroInflate<- opt$options$zeroInflate
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
inputString<-paste(inputSubString, metadataType, nSubjects, nPerSubject, 
                   nMicrobes, spikeMicrobes, nMetadata, effectSize, effectPos, 
                   readDepth, depthConfound, propAbun, zeroInflate, sep='_')
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
outputString <- paste(inputString, "ALDEx2_scale", sep="_")
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
    
    f <- list()
    for (i in reps) {
        metadata <- read.csv(paste0(this_params_folder, "/metadata_", i, ".tsv"), sep = "\t")
        abundance <- read.csv(paste0(this_params_folder, "/abundance_", i, ".tsv"), sep = "\t")
        truth <- read.csv(paste0(this_params_folder, "/truth_", i, ".tsv"), sep = "\t")
        unscaled <- read.csv(paste0(this_params_folder, "/unscaled_", i, ".tsv"), sep = "\t")
        
        if (colnames(unscaled) != 'total') {
            rows_to_drop <- colnames(abundance)[abundance[colnames(unscaled),] == 0]
            metadata <- metadata[!rownames(metadata) %in% rows_to_drop, , drop = F]
            abundance <- abundance[,!colnames(abundance) %in% rows_to_drop, drop = F]
            unscaled <- unscaled[!rownames(unscaled) %in% rows_to_drop, , drop = F]
            
            # Convert to total
            unscaled <- unscaled[,1] / abundance[nrow(abundance),]
        }
        unscaled <- c(unlist(unscaled))
        
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
        
        # Don't convert to relative abundance - ALDEx needs raw reads
        # abundance <- t(t(abundance) / colSums(abundance))
        
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
        
        if(depthConfound) {
            fixed_effects <- colnames(metadata)[!colnames(metadata) %in% c("ID")]
        } else {
            fixed_effects <- colnames(metadata)[!colnames(metadata) %in% c("ID", "read_depth")]
        }
        
        sink('/dev/null')
        if ('ID' %in% colnames(metadata) & length(unique(metadata$ID)) != length(metadata$ID)) {
            mm <- model.matrix(formula(paste0("~ ", 
                                              paste0(c(fixed_effects, "ID"), 
                                                     collapse = " + "), 
                                              collapse = "")),metadata)
        } else {
            mm <- model.matrix(formula(paste0("~ ", 
                                              paste0(fixed_effects, 
                                                     collapse = " + "), 
                                              collapse = "")),metadata)
        }
        
        aldex_clr_out <- aldex.clr(abundance, mm, denom="all", gamma = matrix(unscaled, nrow = length(unscaled), ncol = 128))
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
        f[[i]] <- outputs
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
