##############################
## Synthetic Data Generation #
##############################

trigger_sparseDOSSA2_Simulator<-function(noZeroInflate=FALSE,
                                        RandomEffect=FALSE,
                                        metadataType,
                                        nSubjects,
                                        nPerSubject,
                                        nMicrobes,
                                        spikeMicrobes,
                                        nMetadata,
                                        effectSize,
                                        effectPos,
                                        readDepth = 50000,
                                        nIterations = 100,
                                        noParallel = FALSE,
                                        rSeed = 1234,
                                        nCores = 4){
  
  # Create Replicates 
  reps = 1:nIterations
  
  ########################
  # Catch Obvious Errors #
  ########################
  
  # Check Character Values
  if (!metadataType %in% c('UVA', 'UVB', 'MVA', 'MVB', 'binary'))
    stop('Must be one of the following: UVA, UVB, MVA, MVB, or binary.')
  
  # Check Positive Integer Values
  if (round(nSubjects) != nSubjects || 
      nSubjects<0 ||
      round(nPerSubject) != nPerSubject || 
      nPerSubject<0 ||
      round(nMicrobes) != nMicrobes || 
      nMicrobes<0 ||
      round(nMetadata) != nMetadata ||
      nMetadata<0 ||
      round(readDepth) != readDepth || 
      readDepth<0)
    stop('nSubjects/nPerSubject/nMicrobes/nMetadata/readDepth must be positive integers.')
  
  # Check Proportion Values
  if (spikeMicrobes>1 || spikeMicrobes<=0)
    stop('spikeMicrobes must be in (0, 1].')
  
  # Check Illegal Combinations 
  if(RandomEffect==TRUE && nPerSubject==1)
    stop('nPerSubject must be greater 1 when RandomEffect is TRUE.')
  
  if(RandomEffect==FALSE && nPerSubject>1)
    stop('nPerSubject must be equal to  1 when RandomEffect is FALSE.')
  
  if(metadataType %in% c('UVA', 'UVB') && (nMetadata!=1))
    stop('nMetadata must be equal to 1 when metadataType is UVA or UVB.')
  
  if(metadataType %in% c('MVA', 'MVB') && nMetadata==1)
    stop('nMetadata must be greater than 1 when metadataType is MVA or MVB')
  
  if(metadataType == 'binary') {
    if(RandomEffect) { stop('RandomEffect cannot be true when metadataType is binary.') }
    if(nMetadata != 1) {stop('nMetadata must be equal to 1 when metadataType is binary.')}
  }
  
  # Define the Simulation Parameters Combinations
  simparams = apply(expand.grid(metadataType,
                                nSubjects,
                                nPerSubject,
                                nMicrobes, 
                                spikeMicrobes, 
                                nMetadata,
                                effectSize,
                                effectPos,
                                readDepth,
                                reps), 1, paste, collapse = '_')
  
  # Define the Labels to Go with Each Element of the Simulation Parameter
  simparamslabels = c("metadataType","nSubjects", "nPerSubject", "nMicrobes", "spikeMicrobes", "nMetadata", "effectSize", "effectPos", "readDepth", "rep")
  
  # Track Start Time
  cat(c("Job started at:",date()), "\n")
  start.time <- Sys.time()
  
  # Set Reproducibility Seed
  set.seed(rSeed) 
  
  # Call Grid Computing Only When Specified
  
  if (noParallel){
    
    # Call SparseDOSSA Wrapper (noParallel)
    simlist <- sparseDOSSA2_Wrapper_noParallel(simparams, simparamslabels, noZeroInflate=noZeroInflate)
  }
  else{
    
    # Set Up Clustering Environment
    no_cores <- nCores 
    cl <- makeCluster(no_cores)
    registerDoParallel(cl)
    
    ####################
    # Data Generation #
    ###################
    
    # Call SparseDOSSA Wrapper 
    simlist <- sparseDOSSA2_Wrapper(simparams, simparamslabels, noZeroInflate=noZeroInflate)
    
    # Stop the Cluster 
    stopCluster(cl)
  }
  
  # Set Names
  if (RandomEffect==TRUE) {
    simnames<- paste('RandomEffect', simparams, sep='_')
  } else {
    simnames<- paste('noRandomEffect', simparams, sep='_')
  }
  names(simlist) <- simnames
  
  # Track End Time
  cat(c("Job ended at:",date()), "\n")
  stop.time <- Sys.time()
  minutes <- round(difftime(stop.time, start.time, units="min"), 3)
  cat("Computational time:", minutes, "minutes \n")
  
  # Return
  return(simlist)
}

trigger_sparseDOSSA2_Simulator_omps<-function(noZeroInflate=FALSE,
                                         RandomEffect=FALSE,
                                         metadataType,
                                         nSubjects,
                                         nPerSubject,
                                         nMicrobes,
                                         spikeMicrobes,
                                         nMetadata,
                                         effectSize,
                                         effectPos,
                                         nOmps = nOmps,
                                         nLevels = nLevels,
                                         readDepth = 50000,
                                         nIterations = 100,
                                         noParallel = FALSE,
                                         rSeed = 1234,
                                         nCores = 4){
  
  # Create Replicates 
  reps = 1:nIterations
  
  ########################
  # Catch Obvious Errors #
  ########################
  
  # Check Character Values
  if (!metadataType %in% c('MVAomp', 'MVBomp'))
    stop('Must be one of the following: MVAomp MVBomp')
  
  # Check Positive Integer Values
  if (round(nSubjects) != nSubjects || 
      nSubjects<0 ||
      round(nPerSubject) != nPerSubject || 
      nPerSubject<0 ||
      round(nMicrobes) != nMicrobes || 
      nMicrobes<0 ||
      round(nMetadata) != nMetadata ||
      nMetadata<0 ||
      round(readDepth) != readDepth || 
      readDepth<0)
    stop('nSubjects/nPerSubject/nMicrobes/nMetadata/readDepth must be positive integers.')
  
  # Check Proportion Values
  if (spikeMicrobes>1 || spikeMicrobes<=0)
    stop('spikeMicrobes must be in (0, 1].')
  
  # Check Illegal Combinations 
  if(RandomEffect==TRUE && nPerSubject==1)
    stop('nPerSubject must be greater 1 when RandomEffect is TRUE.')
  
  if(RandomEffect==FALSE && nPerSubject>1)
    stop('nPerSubject must be equal to  1 when RandomEffect is FALSE.')
  
  # Define the Simulation Parameters Combinations
  simparams = apply(expand.grid(metadataType,
                                nSubjects,
                                nPerSubject,
                                nMicrobes, 
                                spikeMicrobes, 
                                nMetadata,
                                effectSize,
                                effectPos,
                                readDepth,
                                nOmps,
                                nLevels,
                                reps), 1, paste, collapse = '_')
  
  # Define the Labels to Go with Each Element of the Simulation Parameter
  simparamslabels = c("metadataType","nSubjects", "nPerSubject", "nMicrobes", "spikeMicrobes", "nMetadata", "effectSize", "effectPos", "readDepth", "nOmps", "nLevels", "rep")
  
  # Track Start Time
  cat(c("Job started at:",date()), "\n")
  start.time <- Sys.time()
  
  # Set Reproducibility Seed
  set.seed(rSeed) 
  
  # Set Up Clustering Environment
  no_cores <- nCores 
  cl <- makeCluster(no_cores)
  registerDoParallel(cl)
  
  ####################
  # Data Generation #
  ###################
  
  # Call SparseDOSSA Wrapper 
  simlist <- sparseDOSSA2_Wrapper_omps(simparams, simparamslabels, noZeroInflate=noZeroInflate)
  
  # Stop the Cluster 
  stopCluster(cl)
  
  # Set Names
  if (RandomEffect==TRUE) {
    simnames<- paste('RandomEffect', simparams, sep='_')
  } else {
    simnames<- paste('noRandomEffect', simparams, sep='_')
  }
  names(simlist) <- simnames
  
  # Track End Time
  cat(c("Job ended at:",date()), "\n")
  stop.time <- Sys.time()
  minutes <- round(difftime(stop.time, start.time, units="min"), 3)
  cat("Computational time:", minutes, "minutes \n")
  
  # Return
  return(simlist)
}

trigger_sparseDOSSA2_Simulator_groups<-function(noZeroInflate=FALSE,
                                              RandomEffect=FALSE,
                                              metadataType,
                                              nSubjects,
                                              nPerSubject,
                                              nMicrobes,
                                              spikeMicrobes,
                                              nMetadata,
                                              effectSize,
                                              effectPos,
                                              nGroups = nGroups,
                                              nLevels = nLevels,
                                              readDepth = 50000,
                                              nIterations = 100,
                                              noParallel = FALSE,
                                              rSeed = 1234,
                                              nCores = 4){
  
  # Create Replicates 
  reps = 1:nIterations
  
  ########################
  # Catch Obvious Errors #
  ########################
  
  # Check Character Values
  if (!metadataType %in% c('MVAgroup', 'MVBgroup'))
    stop('Must be one of the following: MVAgroup MVBgroup')
  
  # Check Positive Integer Values
  if (round(nSubjects) != nSubjects || 
      nSubjects<0 ||
      round(nPerSubject) != nPerSubject || 
      nPerSubject<0 ||
      round(nMicrobes) != nMicrobes || 
      nMicrobes<0 ||
      round(nMetadata) != nMetadata ||
      nMetadata<0 ||
      round(readDepth) != readDepth || 
      readDepth<0)
    stop('nSubjects/nPerSubject/nMicrobes/nMetadata/readDepth must be positive integers.')
  
  # Check Proportion Values
  if (spikeMicrobes>1 || spikeMicrobes<=0)
    stop('spikeMicrobes must be in (0, 1].')
  
  # Check Illegal Combinations 
  if(RandomEffect==TRUE && nPerSubject==1)
    stop('nPerSubject must be greater 1 when RandomEffect is TRUE.')
  
  if(RandomEffect==FALSE && nPerSubject>1)
    stop('nPerSubject must be equal to  1 when RandomEffect is FALSE.')
  
  # Define the Simulation Parameters Combinations
  simparams = apply(expand.grid(metadataType,
                                nSubjects,
                                nPerSubject,
                                nMicrobes, 
                                spikeMicrobes, 
                                nMetadata,
                                effectSize,
                                effectPos,
                                readDepth,
                                nGroups,
                                nLevels,
                                reps), 1, paste, collapse = '_')
  
  # Define the Labels to Go with Each Element of the Simulation Parameter
  simparamslabels = c("metadataType","nSubjects", "nPerSubject", "nMicrobes", "spikeMicrobes", "nMetadata", "effectSize", "effectPos", "readDepth", "nGroups", "nLevels", "rep")
  
  # Track Start Time
  cat(c("Job started at:",date()), "\n")
  start.time <- Sys.time()
  
  # Set Reproducibility Seed
  set.seed(rSeed) 
  
  # Set Up Clustering Environment
  no_cores <- nCores 
  cl <- makeCluster(no_cores)
  registerDoParallel(cl)
  
  ####################
  # Data Generation #
  ###################
  
  # Call SparseDOSSA Wrapper 
  simlist <- sparseDOSSA2_Wrapper_groups(simparams, simparamslabels, noZeroInflate=noZeroInflate)
  
  # Stop the Cluster 
  stopCluster(cl)
  
  # Set Names
  if (RandomEffect==TRUE) {
    simnames<- paste('RandomEffect', simparams, sep='_')
  } else {
    simnames<- paste('noRandomEffect', simparams, sep='_')
  }
  names(simlist) <- simnames
  
  # Track End Time
  cat(c("Job ended at:",date()), "\n")
  stop.time <- Sys.time()
  minutes <- round(difftime(stop.time, start.time, units="min"), 3)
  cat("Computational time:", minutes, "minutes \n")
  
  # Return
  return(simlist)
}

trigger_sparseDOSSA2_Simulator_unscaled<-function(noZeroInflate=FALSE,
                                                  RandomEffect=FALSE,
                                                  metadataType,
                                                  nSubjects,
                                                  nPerSubject,
                                                  nMicrobes,
                                                  spikeMicrobes,
                                                  nMetadata,
                                                  effectSize,
                                                  effectPos,
                                                  readDepth = 50000,
                                                  nIterations = 100,
                                                  noParallel = FALSE,
                                                  rSeed = 1234,
                                                  nCores = 4){
  
  # Create Replicates 
  reps = 1:nIterations
  
  ########################
  # Catch Obvious Errors #
  ########################
  
  # Check Character Values
  if (!metadataType %in% c('MVAtotal', 'MVAref'))
    stop('Must be one of the following: MVAtotal MVAref')
  
  # Check Positive Integer Values
  if (round(nSubjects) != nSubjects || 
      nSubjects<0 ||
      round(nPerSubject) != nPerSubject || 
      nPerSubject<0 ||
      round(nMicrobes) != nMicrobes || 
      nMicrobes<0 ||
      round(nMetadata) != nMetadata ||
      nMetadata<0 ||
      round(readDepth) != readDepth || 
      readDepth<0)
    stop('nSubjects/nPerSubject/nMicrobes/nMetadata/readDepth must be positive integers.')
  
  # Check Proportion Values
  if (spikeMicrobes>1 || spikeMicrobes<=0)
    stop('spikeMicrobes must be in (0, 1].')
  
  # Check Illegal Combinations 
  if(RandomEffect==TRUE && nPerSubject==1)
    stop('nPerSubject must be greater 1 when RandomEffect is TRUE.')
  
  if(RandomEffect==FALSE && nPerSubject>1)
    stop('nPerSubject must be equal to  1 when RandomEffect is FALSE.')
  
  # Define the Simulation Parameters Combinations
  simparams = apply(expand.grid(metadataType,
                                nSubjects,
                                nPerSubject,
                                nMicrobes, 
                                spikeMicrobes, 
                                nMetadata,
                                effectSize,
                                effectPos,
                                readDepth,
                                reps), 1, paste, collapse = '_')
  
  # Define the Labels to Go with Each Element of the Simulation Parameter
  simparamslabels = c("metadataType","nSubjects", "nPerSubject", "nMicrobes", "spikeMicrobes", "nMetadata", "effectSize", "effectPos", "readDepth", "rep")
  
  # Track Start Time
  cat(c("Job started at:",date()), "\n")
  start.time <- Sys.time()
  
  # Set Reproducibility Seed
  set.seed(rSeed) 
  
  # Set Up Clustering Environment
  no_cores <- nCores 
  cl <- makeCluster(no_cores)
  registerDoParallel(cl)
  
  ####################
  # Data Generation #
  ###################
  
  # Call SparseDOSSA Wrapper 
  simlist <- sparseDOSSA2_Wrapper_unscaled(simparams, simparamslabels, noZeroInflate=noZeroInflate)
  
  # Stop the Cluster 
  stopCluster(cl)
  
  # Set Names
  if (RandomEffect==TRUE) {
    simnames<- paste('RandomEffect', simparams, sep='_')
  } else {
    simnames<- paste('noRandomEffect', simparams, sep='_')
  }
  names(simlist) <- simnames
  
  # Track End Time
  cat(c("Job ended at:",date()), "\n")
  stop.time <- Sys.time()
  minutes <- round(difftime(stop.time, start.time, units="min"), 3)
  cat("Computational time:", minutes, "minutes \n")
  
  # Return
  return(simlist)
}

sparseDOSSA2_Wrapper<-function(simparams, simparamslabels, noZeroInflate){
  f<-foreach(i = simparams, .packages = c("SparseDOSSA2", "MASS", "stringi", 'dplyr'),
             .export = c("generateMetadata")) %dopar% {
               
               # Extract Parameter Strings
               params = strsplit(i, '_')[[1]]
               names(params) <- simparamslabels
               
               # Extract Relevant Parameters
               metadataType = as.character(params["metadataType"]) # Type of Metadata
               nSubjects <- as.numeric(params["nSubjects"])  # Number of Subjects
               nPerSubject <- as.numeric(params["nPerSubject"])  # Number of Samples Per Subject
               nSamples<-round(nSubjects*nPerSubject) # Number of Samples
               nMicrobes <- as.numeric(params["nMicrobes"])  # Number of Microbes
               spikeMicrobes <- as.numeric(params["spikeMicrobes"]) # Proportion of Spiked-in Microbes
               nMetadata<-as.numeric(params["nMetadata"])  # Number of Metadata
               effectSize<-as.character(params["effectSize"]) # Effect Size
               effectPos<-as.numeric(params["effectPos"]) # Effect Size
               readDepth<-as.numeric(params["readDepth"]) # Library Size
               
               # Initialize
               DD = NULL
               
               # sparseDOSSA Error Control 
               tryAgain = TRUE
               infiniteloopcounter = 1
               while (tryAgain & infiniteloopcounter < 5) {
                 
                 # Generate Metadata
                 FF<-generateMetadata(metadataType=metadataType, 
                                      nSubjects=nSubjects, 
                                      nPerSubject=nPerSubject, 
                                      nMetadata=nMetadata)
                 
                 # Extract Relevant Information
                 UserMetadata<-FF$UserMetadata; 

                 spike_metadata_df <- data.frame(matrix(ncol = 4, nrow = 0))
                 while (nrow(spike_metadata_df) < spikeMicrobes * nMicrobes * nMetadata) {
                   # Multiply by 10 to allow subsetting to avoid duplication
                   feature_spiked <- sample(paste0("Feature", 1:nMicrobes), spikeMicrobes * nMicrobes * nMetadata * 10, replace = T)
                   metadata_datum <- sample(1:nMetadata, spikeMicrobes * nMicrobes * nMetadata * 10, replace = T)
                   associated_property <- sample(c("abundance", "prevalence"), spikeMicrobes * nMicrobes * nMetadata * 10, replace = T)
                   spike_metadata_df <- data.frame(metadata_datum, feature_spiked, associated_property)
                   spike_metadata_df <- distinct(spike_metadata_df)
                 }
                 spike_metadata_df <- spike_metadata_df[1:floor(spikeMicrobes * nMicrobes * nMetadata),]
                 
                 spike_metadata_df$effect_size <- runif(nrow(spike_metadata_df), as.numeric(effectSize) * 0.5, as.numeric(effectSize)) *
                   sample(c(1, -1), nrow(spike_metadata_df), replace = T, prob = c(effectPos, 1-effectPos))
                 
                 # Generate sparseDOSSA Synthetic Abundances
                 DD<-SparseDOSSA2::SparseDOSSA2(template="Stool",
                                                n_sample = nSamples,
                                                n_feature = nMicrobes,
                                                spike_metadata = spike_metadata_df,
                                                metadata_matrix = t(UserMetadata),
                                                median_read_depth = readDepth,
                                                verbose=F)
                 
                 if (is.null(DD) | inherits(DD, "try-error")) {
                   tryAgain = TRUE
                   infiniteloopcounter = infiniteloopcounter + 1
                 } else {
                   tryAgain = FALSE
                 }
               }
               if (infiniteloopcounter >= 5) {
                 stop("Consistent error found during simulation. Need to investigate cause.")
               }
               
               sparsedossa_results <- DD$simulated_data
               significant_features <- DD$spike_metadata$feature_metadata_spike_df
               significant_features$metadata_datum <- paste0("Metadata_", significant_features$metadata_datum)
               sparsedossa_metadata <- data.frame(DD$spike_metadata$metadata_matrix)
               colnames(sparsedossa_metadata) <- paste0("Metadata_", 1:ncol(sparsedossa_metadata))
               ID <- rep(paste('Subject', 1:nSubjects, sep=''), each = nPerSubject)
               sparsedossa_metadata$ID <- ID
               rownames(sparsedossa_metadata) <- paste0("Sample", 1:nrow(sparsedossa_metadata))
               
               # Return
               return(list(metadata=sparsedossa_metadata, 
                           features=sparsedossa_results, 
                           truth=significant_features, 
                           true_tax=DD$params$feature_param, 
                           ID=ID, 
                           libSize=colSums(sparsedossa_results)))
             }
  return(f)
}

sparseDOSSA2_Wrapper_noParallel<-function(simparams, simparamslabels, noZeroInflate){
  
  # Intitialize
  pclList<-list()
  
  # Repeated Loop 
  for(i in simparams){
    
    # Extract Parameter Strings
    params = strsplit(i, '_')[[1]]
    names(params) <- simparamslabels
    
    # Extract Relevant Parameters
    metadataType = as.character(params["metadataType"]) # Type of Metadata
    nSubjects <- as.numeric(params["nSubjects"])  # Number of Subjects
    nPerSubject <- as.numeric(params["nPerSubject"])  # Number of Samples Per Subject
    nSamples<-round(nSubjects*nPerSubject) # Number of Samples
    nMicrobes <- as.numeric(params["nMicrobes"])  # Number of Microbes
    spikeMicrobes <- as.numeric(params["spikeMicrobes"]) # Proportion of Spiked-in Microbes
    nMetadata<-as.numeric(params["nMetadata"])  # Number of Metadata
    effectSize<-as.character(params["effectSize"]) # Effect Size
    effectPos<-as.numeric(params["effectPos"]) # Effect Size
    readDepth<-as.numeric(params["readDepth"]) # Library Size
    
    # Initialize
    DD = NULL
    
    # sparseDOSSA Error Control 
    tryAgain = TRUE
    infiniteloopcounter = 1
    while (tryAgain & infiniteloopcounter < 5) {
      
      # Generate Metadata
      FF<-generateMetadata(metadataType=metadataType, 
                           nSubjects=nSubjects, 
                           nPerSubject=nPerSubject, 
                           nMetadata=nMetadata)
      
      # Extract Relevant Information
      UserMetadata<-FF$UserMetadata; 
      
      spike_metadata_df <- data.frame(matrix(ncol = 4, nrow = 0))
      while (nrow(spike_metadata_df) < spikeMicrobes * nMicrobes * nMetadata) {
        # Multiply by 10 to allow subsetting to avoid duplication
        feature_spiked <- sample(paste0("Feature", 1:nMicrobes), spikeMicrobes * nMicrobes * nMetadata * 10, replace = T)
        metadata_datum <- sample(1:nMetadata, spikeMicrobes * nMicrobes * nMetadata * 10, replace = T)
        associated_property <- sample(c("abundance", "prevalence"), spikeMicrobes * nMicrobes * nMetadata * 10, replace = T)
        spike_metadata_df <- data.frame(metadata_datum, feature_spiked, associated_property)
        spike_metadata_df <- distinct(spike_metadata_df)
      }
      spike_metadata_df <- spike_metadata_df[1:floor(spikeMicrobes * nMicrobes * nMetadata),]
      
      spike_metadata_df$effect_size <- runif(nrow(spike_metadata_df), as.numeric(effectSize) * 0.5, as.numeric(effectSize)) *
        sample(c(1, -1), nrow(spike_metadata_df), replace = T, prob = c(effectPos, 1-effectPos))
      
      # Generate sparseDOSSA Synthetic Abundances
      DD<-SparseDOSSA2::SparseDOSSA2(template="Stool",
                                     n_sample = nSamples,
                                     n_feature = nMicrobes,
                                     spike_metadata = spike_metadata_df,
                                     metadata_matrix = t(UserMetadata),
                                     median_read_depth = readDepth,
                                     verbose=F)
      
      if (is.null(DD) | inherits(DD, "try-error")) {
        tryAgain = TRUE
        infiniteloopcounter = infiniteloopcounter + 1
      } else {
        tryAgain = FALSE
      }
    }
    if (infiniteloopcounter >= 5) {
      stop("Consistent error found during simulation. Need to investigate cause.")
    }
    
    sparsedossa_results <- DD$simulated_data
    significant_features <- DD$spike_metadata$feature_metadata_spike_df
    sparsedossa_metadata <- data.frame(DD$spike_metadata$metadata_matrix)
    colnames(sparsedossa_metadata) <- paste0("Metadata_", 1:ncol(sparsedossa_metadata))
    ID <- rep(paste('Subject', 1:nSubjects, sep=''), each = nPerSubject)
    sparsedossa_metadata$ID <- ID
    
    # Save
    pclList[[i]] <- list(metadata=sparsedossa_metadata, 
                         features=sparsedossa_results, 
                         truth=significant_features, 
                         true_tax=DD$params$feature_param, 
                         ID=ID, 
                         libSize=colSums(sparsedossa_results))
  }
  
  # Return
  return(pclList)
}

sparseDOSSA2_Wrapper_omps <- function(simparams, simparamslabels, noZeroInflate){
  f<-foreach(i = simparams, .packages = c("SparseDOSSA2", "MASS", "stringi", 'plyr', 'dplyr', 'MCMCprecision'),
             .export = c("generateMetadata_omps")) %dopar% {
               
               # Extract Parameter Strings
               params = strsplit(i, '_')[[1]]
               names(params) <- simparamslabels
               
               # Extract Relevant Parameters
               metadataType = as.character(params["metadataType"]) # Type of Metadata
               nSubjects <- as.numeric(params["nSubjects"])  # Number of Subjects
               nPerSubject <- as.numeric(params["nPerSubject"])  # Number of Samples Per Subject
               nSamples<-round(nSubjects*nPerSubject) # Number of Samples
               nMicrobes <- as.numeric(params["nMicrobes"])  # Number of Microbes
               spikeMicrobes <- as.numeric(params["spikeMicrobes"]) # Proportion of Spiked-in Microbes
               nMetadata<-as.numeric(params["nMetadata"])  # Number of Metadata
               effectSize<-as.character(params["effectSize"]) # Effect Size
               effectPos<-as.numeric(params["effectPos"]) # Proportion positive
               nOmps<-as.numeric(params["nOmps"]) # Effect Size
               nLevels<-as.numeric(params["nLevels"]) # Effect Size
               readDepth<-as.numeric(params["readDepth"]) # Library Size
               
               # Initialize
               DD = NULL
               
               # sparseDOSSA Error Control 
               tryAgain = TRUE
               infiniteloopcounter = 1
               while (tryAgain & infiniteloopcounter < 5) {
                 
                 # Generate Metadata
                 FF<-generateMetadata_omps(metadataType=metadataType, 
                                      nSubjects=nSubjects, 
                                      nPerSubject=nPerSubject, 
                                      nMetadata=nMetadata,
                                      nOmps=nOmps, 
                                      nLevels=nLevels)
                 
                 # Extract Relevant Information
                 UserMetadata <- FF$UserMetadata; columns_not_to_level <- FF$columns_not_to_level; unaltered_metadata <- FF$original_metadata
                 
                 spike_metadata_df <- data.frame(matrix(ncol = 4, nrow = 0))
                 totalMetadata <- nMetadata + nOmps
                 while (nrow(spike_metadata_df) < spikeMicrobes * nMicrobes * totalMetadata) {
                   # Multiply by 10 to allow subsetting to avoid duplication
                   feature_spiked <- sample(paste0("Feature", 1:nMicrobes), spikeMicrobes * nMicrobes * totalMetadata * 10, replace = T)
                   metadata_datum <- sample(1:totalMetadata, spikeMicrobes * nMicrobes * totalMetadata * 10, replace = T)
                   associated_property <- sample(c("abundance", "prevalence"), spikeMicrobes * nMicrobes * totalMetadata * 10, replace = T)
                   spike_metadata_df <- data.frame(metadata_datum, feature_spiked, associated_property)
                   spike_metadata_df <- distinct(spike_metadata_df)
                 }
                 spike_metadata_df <- spike_metadata_df[1:floor(spikeMicrobes * nMicrobes * totalMetadata),]
                 spike_metadata_df$metadata_datum <- LETTERS[spike_metadata_df$metadata_datum]
                 
                 spike_metadata_df$effect_size <- runif(nrow(spike_metadata_df), as.numeric(effectSize) * 0.5, as.numeric(effectSize)) *
                   sample(c(1, -1), nrow(spike_metadata_df), replace = T, prob = c(effectPos, 1-effectPos))
                 
                 new_spike_metadata_df <- data.frame(matrix(ncol = ncol(spike_metadata_df), nrow = 0))
                 for (i in 1:nrow(spike_metadata_df)) {
                   if (!spike_metadata_df[i, 'metadata_datum'] %in% columns_not_to_level) {
                     new_addition <- do.call(rbind, replicate(nLevels - 1, spike_metadata_df[i,], simplify = FALSE))
                     new_addition$metadata_datum <- paste0(new_addition$metadata_datum, 2:nLevels)
                     effect_divisions <- MCMCprecision::rdirichlet(1, rep(1, nLevels - 1))
                     for (i in seq_along(effect_divisions[-1])) {
                       if (sample(c(0,1), 1) == 1) {
                         effect_divisions[i + 1] = effect_divisions[i + 1] + effect_divisions[i]
                         effect_divisions[i] <- 0
                       }
                     }
                     new_addition$effect_size <- c(new_addition$effect_size[1] * effect_divisions)
                     new_spike_metadata_df <- rbind(new_spike_metadata_df, new_addition)
                   } else {
                     new_spike_metadata_df <- rbind(new_spike_metadata_df, spike_metadata_df[i,])
                   }
                 }
                 
                 new_spike_metadata_df$metadata_datum <- as.numeric(mapvalues(new_spike_metadata_df$metadata_datum, 
                                                                              rownames(UserMetadata), 
                                                                              1:nrow(UserMetadata)))
                 
                 # Generate sparseDOSSA Synthetic Abundances
                 DD<-SparseDOSSA2::SparseDOSSA2(template="Stool",
                                                n_sample = nSamples,
                                                n_feature = nMicrobes,
                                                spike_metadata = new_spike_metadata_df,
                                                metadata_matrix = t(UserMetadata),
                                                median_read_depth = readDepth,
                                                verbose=F)
                 
                 if (is.null(DD) | inherits(DD, "try-error")) {
                   tryAgain = TRUE
                   infiniteloopcounter = infiniteloopcounter + 1
                 } else {
                   tryAgain = FALSE
                 }
               }
               if (infiniteloopcounter >= 5) {
                 stop("Consistent error found during simulation. Need to investigate cause.")
               }
               
               sparsedossa_results <- DD$simulated_data
               significant_features <- DD$spike_metadata$feature_metadata_spike_df
               significant_features$metadata_datum <- paste0("Metadata_", rownames(UserMetadata)[significant_features$metadata_datum])
               sparsedossa_metadata <- data.frame(unaltered_metadata)
               colnames(sparsedossa_metadata) <- paste0("Metadata_", colnames(sparsedossa_metadata))
               ID <- rep(paste('Subject', 1:nSubjects, sep=''), each = nPerSubject)
               sparsedossa_metadata$ID <- ID
               rownames(sparsedossa_metadata) <- paste0("Sample", 1:nrow(sparsedossa_metadata))
               colnames(sparsedossa_results) <- paste0("Sample", 1:ncol(sparsedossa_results))
               
               # Return
               return(list(metadata=sparsedossa_metadata, 
                           features=sparsedossa_results, 
                           truth=significant_features, 
                           true_tax=DD$params$feature_param, 
                           ID=ID, 
                           libSize=colSums(sparsedossa_results)))
             }
  return(f)
}

sparseDOSSA2_Wrapper_groups <- function(simparams, simparamslabels, noZeroInflate){
  f<-foreach(i = simparams, .packages = c("SparseDOSSA2", "MASS", "stringi", 'plyr', 'dplyr', 'MCMCprecision'),
             .export = c("generateMetadata_groups")) %dopar% {
               
               # Extract Parameter Strings
               params = strsplit(i, '_')[[1]]
               names(params) <- simparamslabels
               
               # Extract Relevant Parameters
               metadataType = as.character(params["metadataType"]) # Type of Metadata
               nSubjects <- as.numeric(params["nSubjects"])  # Number of Subjects
               nPerSubject <- as.numeric(params["nPerSubject"])  # Number of Samples Per Subject
               nSamples<-round(nSubjects*nPerSubject) # Number of Samples
               nMicrobes <- as.numeric(params["nMicrobes"])  # Number of Microbes
               spikeMicrobes <- as.numeric(params["spikeMicrobes"]) # Proportion of Spiked-in Microbes
               nMetadata<-as.numeric(params["nMetadata"])  # Number of Metadata
               effectSize<-as.character(params["effectSize"]) # Effect Size
               effectPos<-as.numeric(params["effectPos"]) # Effect Size
               nGroups<-as.numeric(params["nGroups"]) # Effect Size
               nLevels<-as.numeric(params["nLevels"]) # Effect Size
               readDepth<-as.numeric(params["readDepth"]) # Library Size
               
               # Initialize
               DD = NULL
               
               # sparseDOSSA Error Control 
               tryAgain = TRUE
               infiniteloopcounter = 1
               while (tryAgain & infiniteloopcounter < 5) {
                 
                 # Generate Metadata
                 FF<-generateMetadata_groups(metadataType=metadataType, 
                                           nSubjects=nSubjects, 
                                           nPerSubject=nPerSubject, 
                                           nMetadata=nMetadata,
                                           nGroups=nGroups, 
                                           nLevels=nLevels)
                 
                 # Extract Relevant Information
                 UserMetadata <- FF$UserMetadata; columns_not_to_level <- FF$columns_not_to_level; unaltered_metadata <- FF$original_metadata
                 
                 spike_metadata_df <- data.frame(matrix(ncol = 4, nrow = 0))
                 totalMetadata <- nMetadata + nGroups
                 while (nrow(spike_metadata_df) < spikeMicrobes * nMicrobes * totalMetadata) {
                   # Multiply by 10 to allow subsetting to avoid duplication
                   feature_spiked <- sample(paste0("Feature", 1:nMicrobes), spikeMicrobes * nMicrobes * totalMetadata * 10, replace = T)
                   metadata_datum <- sample(1:totalMetadata, spikeMicrobes * nMicrobes * totalMetadata * 10, replace = T)
                   associated_property <- sample(c("abundance", "prevalence"), spikeMicrobes * nMicrobes * totalMetadata * 10, replace = T)
                   spike_metadata_df <- data.frame(metadata_datum, feature_spiked, associated_property)
                   spike_metadata_df <- distinct(spike_metadata_df)
                 }
                 spike_metadata_df <- spike_metadata_df[1:floor(spikeMicrobes * nMicrobes * totalMetadata),]
                 spike_metadata_df$metadata_datum <- LETTERS[spike_metadata_df$metadata_datum]
                 
                 spike_metadata_df$effect_size <- runif(nrow(spike_metadata_df), as.numeric(effectSize) * 0.5, as.numeric(effectSize)) *
                   sample(c(1, -1), nrow(spike_metadata_df), replace = T, prob = c(effectPos, 1-effectPos))
                 
                 new_spike_metadata_df <- data.frame(matrix(ncol = ncol(spike_metadata_df), nrow = 0))
                 for (i in 1:nrow(spike_metadata_df)) {
                   if (!spike_metadata_df[i, 'metadata_datum'] %in% columns_not_to_level) {
                     new_addition <- do.call(rbind, replicate(nLevels - 1, spike_metadata_df[i,], simplify = FALSE))
                     new_addition$metadata_datum <- paste0(new_addition$metadata_datum, 2:nLevels)
                     effect_divisions <- MCMCprecision::rdirichlet(1, rep(1, nLevels - 1))
                     for (i in seq_along(effect_divisions[-1])) {
                       if (sample(c(0,1), 1) == 1) {
                         effect_divisions[i + 1] = effect_divisions[i + 1] + effect_divisions[i]
                         effect_divisions[i] <- 0
                       }
                     }
                     new_addition$effect_size <- c(new_addition$effect_size[1] * effect_divisions)
                     new_spike_metadata_df <- rbind(new_spike_metadata_df, new_addition)
                   } else {
                     new_spike_metadata_df <- rbind(new_spike_metadata_df, spike_metadata_df[i,])
                   }
                 }
                 
                 new_spike_metadata_df$metadata_datum <- as.numeric(mapvalues(new_spike_metadata_df$metadata_datum, 
                                                                              rownames(UserMetadata), 
                                                                              1:nrow(UserMetadata)))
                 
                 # Generate sparseDOSSA Synthetic Abundances
                 DD<-SparseDOSSA2::SparseDOSSA2(template="Stool",
                                                n_sample = nSamples,
                                                n_feature = nMicrobes,
                                                spike_metadata = new_spike_metadata_df,
                                                metadata_matrix = t(UserMetadata),
                                                median_read_depth = readDepth,
                                                verbose=F)
                 
                 if (is.null(DD) | inherits(DD, "try-error")) {
                   tryAgain = TRUE
                   infiniteloopcounter = infiniteloopcounter + 1
                 } else {
                   tryAgain = FALSE
                 }
               }
               if (infiniteloopcounter >= 5) {
                 stop("Consistent error found during simulation. Need to investigate cause.")
               }
               
               sparsedossa_results <- DD$simulated_data
               significant_features <- DD$spike_metadata$feature_metadata_spike_df
               significant_features$metadata_datum <- paste0("Metadata_", rownames(UserMetadata)[significant_features$metadata_datum])
               sparsedossa_metadata <- data.frame(unaltered_metadata)
               colnames(sparsedossa_metadata) <- paste0("Metadata_", colnames(sparsedossa_metadata))
               ID <- rep(paste('Subject', 1:nSubjects, sep=''), each = nPerSubject)
               sparsedossa_metadata$ID <- ID
               rownames(sparsedossa_metadata) <- paste0("Sample", 1:nrow(sparsedossa_metadata))
               colnames(sparsedossa_results) <- paste0("Sample", 1:ncol(sparsedossa_results))
               
               # Return
               return(list(metadata=sparsedossa_metadata, 
                           features=sparsedossa_results, 
                           truth=significant_features, 
                           true_tax=DD$params$feature_param, 
                           ID=ID, 
                           libSize=colSums(sparsedossa_results)))
             }
  return(f)
}

sparseDOSSA2_Wrapper_unscaled <- function(simparams, simparamslabels, noZeroInflate){
  args <- commandArgs(trailingOnly = FALSE)
  
  script_index <- grep("--file=", args)
  if (length(script_index) > 0) {
    script_filepath <- sub("--file=", "", args[script_index])
  } else {
    script_filepath <- args[1]
  }
  this_file <- normalizePath(script_filepath)
  
  source(file.path(paste0(unlist(strsplit(this_file, "/", perl = TRUE))[1:(length(unlist(strsplit(this_file, "/", perl = TRUE))) - 3)], collapse = '/'), 'library/SD_spike.R'))
  f<-foreach(i = simparams, .packages = c("SparseDOSSA2", "MASS", "stringi", 'plyr', 'dplyr'),
             .export = c("generateMetadata_unscaled", "SparseDOSSA2_spike")) %dopar% {
               
               # Extract Parameter Strings
               params = strsplit(i, '_')[[1]]
               names(params) <- simparamslabels
               
               # Extract Relevant Parameters
               metadataType = as.character(params["metadataType"]) # Type of Metadata
               nSubjects <- as.numeric(params["nSubjects"])  # Number of Subjects
               nPerSubject <- as.numeric(params["nPerSubject"])  # Number of Samples Per Subject
               nSamples<-round(nSubjects*nPerSubject) # Number of Samples
               nMicrobes <- as.numeric(params["nMicrobes"])  # Number of Microbes
               spikeMicrobes <- as.numeric(params["spikeMicrobes"]) # Proportion of Spiked-in Microbes
               nMetadata<-as.numeric(params["nMetadata"])  # Number of Metadata
               effectSize<-as.character(params["effectSize"]) # Effect Size
               effectPos<-as.numeric(params["effectPos"]) # Effect Size
               readDepth<-as.numeric(params["readDepth"]) # Library Size
               
               # Initialize
               DD = NULL
               
               # sparseDOSSA Error Control 
               tryAgain = TRUE
               infiniteloopcounter = 1
               while (tryAgain & infiniteloopcounter < 5) {
                 
                 # Generate Metadata
                 FF<-generateMetadata_unscaled(metadataType=metadataType, 
                                             nSubjects=nSubjects, 
                                             nPerSubject=nPerSubject, 
                                             nMetadata=nMetadata)
                 
                 # Extract Relevant Information
                 UserMetadata<-FF$UserMetadata; 
                 
                 spike_metadata_df <- data.frame(matrix(ncol = 4, nrow = 0))
                 while (nrow(spike_metadata_df) < spikeMicrobes * nMicrobes * nMetadata) {
                   # Multiply by 10 to allow subsetting to avoid duplication
                   feature_spiked <- sample(paste0("Feature", 1:nMicrobes), spikeMicrobes * nMicrobes * nMetadata * 10, replace = T)
                   metadata_datum <- sample(1:nMetadata, spikeMicrobes * nMicrobes * nMetadata * 10, replace = T)
                   associated_property <- sample(c("abundance", "prevalence"), spikeMicrobes * nMicrobes * nMetadata * 10, replace = T)
                   spike_metadata_df <- data.frame(metadata_datum, feature_spiked, associated_property)
                   spike_metadata_df <- distinct(spike_metadata_df)
                 }
                 spike_metadata_df <- spike_metadata_df[1:floor(spikeMicrobes * nMicrobes * nMetadata),]
                 
                 spike_metadata_df$effect_size <- runif(nrow(spike_metadata_df), as.numeric(effectSize) * 0.5, as.numeric(effectSize)) *
                   sample(c(1, -1), nrow(spike_metadata_df), replace = T, prob = c(effectPos, 1-effectPos))
                 
                 # Generate sparseDOSSA Synthetic Abundances
                 DD<-SparseDOSSA2_spike(template="Stool",
                                                n_sample = nSamples,
                                                n_feature = nMicrobes,
                                                spike_metadata = spike_metadata_df,
                                                metadata_matrix = t(UserMetadata),
                                                median_read_depth = readDepth,
                                                verbose=F)
                 
                 if (is.null(DD) | inherits(DD, "try-error")) {
                   tryAgain = TRUE
                   infiniteloopcounter = infiniteloopcounter + 1
                 } else {
                   tryAgain = FALSE
                 }
               }
               if (infiniteloopcounter >= 5) {
                 stop("Consistent error found during simulation. Need to investigate cause.")
               }
               
               sparsedossa_results <- DD$simulated_data
               significant_features <- DD$spike_metadata$feature_metadata_spike_df
               significant_features$metadata_datum <- paste0("Metadata_", significant_features$metadata_datum)
               sparsedossa_metadata <- data.frame(DD$spike_metadata$metadata_matrix)
               colnames(sparsedossa_metadata) <- paste0("Metadata_", 1:ncol(sparsedossa_metadata))
               ID <- rep(paste('Subject', 1:nSubjects, sep=''), each = nPerSubject)
               sparsedossa_metadata$ID <- ID
               rownames(sparsedossa_metadata) <- paste0("Sample", 1:nrow(sparsedossa_metadata))
               
               if (metadataType == 'MVAtotal') {
                 abundance_unscaled <- data.frame('total' = colSums(DD$simulated_matrices$a_spiked))
                 rownames(abundance_unscaled) <- colnames(DD$simulated_matrices$a_spiked)
               } else if (metadataType == 'MVAref') {
                 abundance_unscaled <- data.frame('reference' = DD$spike_in_abs_abun)
                 colnames(abundance_unscaled) <- rownames(sparsedossa_results)[nrow(sparsedossa_results)]
                 rownames(abundance_unscaled) <- colnames(DD$simulated_matrices$a_spiked)
               }
               
               # Return
               return(list(metadata=sparsedossa_metadata, 
                           features=sparsedossa_results, 
                           truth=significant_features, 
                           true_tax=DD$params$feature_param, 
                           ID=ID, 
                           libSize=colSums(sparsedossa_results),
                           abundance_unscaled = abundance_unscaled))
             }
  return(f)
}

#####################
# Generate Metadata #
#####################

generateMetadata<-function(metadataType, 
                           nSubjects, 
                           nPerSubject, 
                           nMetadata){
  
  # Calculate Number of Samples
  nSamples = round(nSubjects*nPerSubject)
  
  # Create Blocking Variable
  if (nPerSubject==1){  # NO RANDOM EFFECTS 
    subjectRandomEffects <- as.matrix(rnorm(nSubjects,mean=0,sd=0))
  }
  if (nPerSubject>1){  # SUBJECT-SPECIFIC RANDOM EFFECTS
    subjectRandomEffects <- as.matrix(rnorm(nSubjects,mean=0,sd=1))
  }
  BLOCK <- as.vector(matrix(subjectRandomEffects, nrow=nPerSubject, ncol=length(subjectRandomEffects), byrow=TRUE))
  
  # Specify Mean and Covariance Structure
  mu<-rep(0,nMetadata)
  cov<-diag(1,nMetadata, nMetadata)
  
  if (metadataType == 'MVB'){
    for (i in 1:nMetadata){
      for (j in 1:nMetadata){
        if(i!=j) cov[i,j]=0.5**(abs(i-j)) # AR(1)
      }
    }
  }
  
  # Generate from MVN
  fakeMetadata<-as.matrix(MASS::mvrnorm(n=nSamples, mu, cov))
  
  # Transpose and Add Blocking Structure
  finalMetadata<-apply(fakeMetadata, 2, function(x) x+BLOCK)
  
  #############################
  # Modularize Specific Cases #
  #############################
  
  # Multivariable Scenario - Dichotomize Half of the Features
  if (metadataType %in% c('MVA', 'MVB')){
    t_UserMetadata<-apply(finalMetadata, 2, function(x) ifelse(x>median(x), 1, 0))
    columns_not_to_binarize<-sample(1:nMetadata, nMetadata/2)
    t_UserMetadata[,columns_not_to_binarize]<-finalMetadata[, columns_not_to_binarize]
    UserMetadata<-t(t_UserMetadata)
  }
  
  # Univariate Binary
  else if (metadataType %in% c('UVB', 'binary')){
    UserMetadata<-t(apply(finalMetadata, 2, function(x) ifelse(x>median(x), 1, 0)))
  }
  # Univariate Continuous
  else if (metadataType == 'UVA') {
    UserMetadata<-t(finalMetadata)
  } else {
    stop("Unrecognized metadataType")
  }

  # Return 
  return(list(UserMetadata=UserMetadata))
}

generateMetadata_omps<-function(metadataType, 
                           nSubjects, 
                           nPerSubject, 
                           nMetadata,
                           nOmps, 
                           nLevels) {
  
  # Calculate Number of Samples
  nSamples = round(nSubjects*nPerSubject)
  
  # Create Blocking Variable
  if (nPerSubject==1){  # NO RANDOM EFFECTS 
    subjectRandomEffects <- as.matrix(rnorm(nSubjects,mean=0,sd=0))
  }
  if (nPerSubject>1){  # SUBJECT-SPECIFIC RANDOM EFFECTS
    subjectRandomEffects <- as.matrix(rnorm(nSubjects,mean=0,sd=1))
  }
  BLOCK <- as.vector(matrix(subjectRandomEffects, nrow=nPerSubject, ncol=length(subjectRandomEffects), byrow=TRUE))
  
  totalMetadata <- nMetadata + nOmps
  
  # Specify Mean and Covariance Structure
  mu<-rep(0,totalMetadata)
  cov<-diag(1,totalMetadata, totalMetadata)
  
  if (metadataType == 'MVBomp'){
    for (i in 1:totalMetadata){
      for (j in 1:totalMetadata){
        if(i!=j) cov[i,j]=0.5**(abs(i-j)) # AR(1)
      }
    }
  }
  
  # Generate from MVN
  fakeMetadata<-as.matrix(MASS::mvrnorm(n=nSamples, mu, cov))
  
  # Transpose and Add Blocking Structure
  finalMetadata<-apply(fakeMetadata, 2, function(x) x+BLOCK)
  
  #############################
  # Modularize Specific Cases #
  #############################
  
  # Multivariable Scenario - Dichotomize Half of the Features
  t_UserMetadata <- finalMetadata
  for (i in 1:nrow(finalMetadata)) {
    for (j in 1:ncol(finalMetadata)) {
      t_UserMetadata[i,j] <- max(c(which(t_UserMetadata[i,j] > quantile(finalMetadata[,j], seq(0, 1, 1/nLevels))), 1))
    }
  }
  columns_not_to_level<-sample(1:totalMetadata, nMetadata)
  t_UserMetadata[,columns_not_to_level]<-finalMetadata[, columns_not_to_level]
  
  t_UserMetadata <- data.frame(t_UserMetadata)
  if (length(columns_not_to_level) == 0) {
    t_UserMetadata <- data.frame(apply(t_UserMetadata, 2, as.factor))
  } else {
    t_UserMetadata[,-columns_not_to_level] <- data.frame(apply(t_UserMetadata[,-columns_not_to_level, drop=F], 2, as.factor))
  }
  for (col_num in setdiff(1:ncol(t_UserMetadata), columns_not_to_level)) {
      t_UserMetadata[,col_num] <- as.factor(t_UserMetadata[,col_num])
  }
  colnames(t_UserMetadata) <- LETTERS[1:ncol(t_UserMetadata)]
  
  thermometer_encode <- function(factor_col) {
      levels <- levels(factor_col)
      encoded_matrix <- sapply(2:length(levels), function(i) as.numeric(as.integer(factor_col) >= i))
      colnames(encoded_matrix) <- levels[-1]
      return(encoded_matrix)
  }
  
  t_UserMetadata_final <- t_UserMetadata %>%
      dplyr::mutate(across(where(is.factor), thermometer_encode))
  
  t_UserMetadata_final <- model.matrix(formula(' ~ .'), t_UserMetadata_final)[,-1, drop=F]
  
  UserMetadata<-t(t_UserMetadata_final)
  
  # Return 
  return(list(UserMetadata=UserMetadata, columns_not_to_level=LETTERS[columns_not_to_level], original_metadata = t_UserMetadata))
}

generateMetadata_groups<-function(metadataType, 
                                nSubjects, 
                                nPerSubject, 
                                nMetadata,
                                nGroups, 
                                nLevels){
  
  # Calculate Number of Samples
  nSamples = round(nSubjects*nPerSubject)
  
  # Create Blocking Variable
  if (nPerSubject==1){  # NO RANDOM EFFECTS 
    subjectRandomEffects <- as.matrix(rnorm(nSubjects,mean=0,sd=0))
  }
  if (nPerSubject>1){  # SUBJECT-SPECIFIC RANDOM EFFECTS
    subjectRandomEffects <- as.matrix(rnorm(nSubjects,mean=0,sd=1))
  }
  BLOCK <- as.vector(matrix(subjectRandomEffects, nrow=nPerSubject, ncol=length(subjectRandomEffects), byrow=TRUE))
  
  totalMetadata <- nMetadata + nGroups
  
  # Specify Mean and Covariance Structure
  mu<-rep(0,totalMetadata)
  cov<-diag(1,totalMetadata, totalMetadata)
  
  if (metadataType == 'MVBgroup'){
    for (i in 1:totalMetadata){
      for (j in 1:totalMetadata){
        if(i!=j) cov[i,j]=0.5**(abs(i-j)) # AR(1)
      }
    }
  }
  
  # Generate from MVN
  fakeMetadata<-as.matrix(MASS::mvrnorm(n=nSamples, mu, cov))
  
  # Transpose and Add Blocking Structure
  finalMetadata<-apply(fakeMetadata, 2, function(x) x+BLOCK)
  
  #############################
  # Modularize Specific Cases #
  #############################
  
  # Multivariable Scenario - Dichotomize Half of the Features
  t_UserMetadata <- finalMetadata
  for (i in 1:nrow(finalMetadata)) {
      for (j in 1:ncol(finalMetadata)) {
          t_UserMetadata[i,j] <- max(c(which(t_UserMetadata[i,j] > quantile(finalMetadata[,j], seq(0, 1, 1/nLevels))), 1))
      }
  }
  columns_not_to_level<-sample(1:totalMetadata, nMetadata)
  t_UserMetadata[,columns_not_to_level]<-finalMetadata[, columns_not_to_level]
  
  t_UserMetadata <- data.frame(t_UserMetadata)
  if (length(columns_not_to_level) == 0) {
      t_UserMetadata <- data.frame(apply(t_UserMetadata, 2, as.factor))
  } else {
      t_UserMetadata[,-columns_not_to_level] <- data.frame(apply(t_UserMetadata[,-columns_not_to_level, drop=F], 2, as.factor))
  }
  for (col_num in setdiff(1:ncol(t_UserMetadata), columns_not_to_level)) {
      t_UserMetadata[,col_num] <- as.factor(t_UserMetadata[,col_num])
  }
  colnames(t_UserMetadata) <- LETTERS[1:ncol(t_UserMetadata)]
  
  thermometer_encode <- function(factor_col) {
      levels <- levels(factor_col)
      encoded_matrix <- sapply(2:length(levels), function(i) as.numeric(as.integer(factor_col) >= i))
      colnames(encoded_matrix) <- levels[-1]
      return(encoded_matrix)
  }
  
  t_UserMetadata_final <- t_UserMetadata %>%
      dplyr::mutate(across(where(is.factor), thermometer_encode))
  
  t_UserMetadata_final <- model.matrix(formula(' ~ .'), t_UserMetadata_final)[,-1, drop=F]
  
  UserMetadata<-t(t_UserMetadata_final)
  
  # Return 
  return(list(UserMetadata=UserMetadata, columns_not_to_level=LETTERS[columns_not_to_level], original_metadata = t_UserMetadata))
}

generateMetadata_unscaled<-function(metadataType, 
                           nSubjects, 
                           nPerSubject, 
                           nMetadata){
  
  # Calculate Number of Samples
  nSamples = round(nSubjects*nPerSubject)
  
  # Create Blocking Variable
  if (nPerSubject==1){  # NO RANDOM EFFECTS 
    subjectRandomEffects <- as.matrix(rnorm(nSubjects,mean=0,sd=0))
  }
  if (nPerSubject>1){  # SUBJECT-SPECIFIC RANDOM EFFECTS
    subjectRandomEffects <- as.matrix(rnorm(nSubjects,mean=0,sd=1))
  }
  BLOCK <- as.vector(matrix(subjectRandomEffects, nrow=nPerSubject, ncol=length(subjectRandomEffects), byrow=TRUE))
  
  # Specify Mean and Covariance Structure
  mu<-rep(0,nMetadata)
  cov<-diag(1,nMetadata, nMetadata)
  
  if (metadataType == 'MVB'){
    for (i in 1:nMetadata){
      for (j in 1:nMetadata){
        if(i!=j) cov[i,j]=0.5**(abs(i-j)) # AR(1)
      }
    }
  }
  
  # Generate from MVN
  fakeMetadata<-as.matrix(MASS::mvrnorm(n=nSamples, mu, cov))
  
  # Transpose and Add Blocking Structure
  finalMetadata<-apply(fakeMetadata, 2, function(x) x+BLOCK)
  
  #############################
  # Modularize Specific Cases #
  #############################
  
  # Multivariable Scenario - Dichotomize Half of the Features
  if (metadataType %in% c('MVAtotal', 'MVAref')){
    t_UserMetadata<-apply(finalMetadata, 2, function(x) ifelse(x>median(x), 1, 0))
    columns_not_to_binarize<-sample(1:nMetadata, nMetadata/2)
    t_UserMetadata[,columns_not_to_binarize]<-finalMetadata[, columns_not_to_binarize]
    UserMetadata<-t(t_UserMetadata)
  } else {
    stop("Unrecognized metadataType")
  }
  
  # Return 
  return(list(UserMetadata=UserMetadata))
}
