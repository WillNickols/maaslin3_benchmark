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
                                        readDepth = 500000,
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


