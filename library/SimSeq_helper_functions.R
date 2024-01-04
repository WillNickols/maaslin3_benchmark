trigger_SimSeq_Simulator<-function(metadataType,
                                   nSubjects,
                                   nMicrobes,
                                   spikeMicrobes, 
                                   nIterations,
                                   rSeed,
                                   nCores,
                                   workingDirectory){
  
  # Create Replicates 
  reps = 1:nIterations
  
  ########################
  # Catch Obvious Errors #
  ########################
  
  # Check Character Values
  if (!metadataType %in% c("ArcticTransects", "cdi_schubert", "crc_baxter", "edd_singh", "ob_goodrich"))
    stop('Must be one of the following: ArcticTransects, cdi_schubert, crc_baxter, edd_singh, or ob_goodrich')
  
  # Check Positive Integer Values
  if (round(nSubjects) != nSubjects || 
      nSubjects<0 ||
      round(nMicrobes) != nMicrobes || 
      nMicrobes<0)
    stop('nSubjects/nMicrobes must be positive integers.')
  
  # Check Proportion Values
  if (spikeMicrobes>1 || spikeMicrobes<=0)
    stop('spikeMicrobes must be in (0, 1].')
  
  # Define the Simulation Parameters Combinations
  simparams = apply(expand.grid(nSubjects,
                                nMicrobes, 
                                spikeMicrobes, 
                                reps), 1, paste, collapse = '_')
  
  # Define the Labels to Go with Each Element of the Simulation Parameter
  simparamslabels = c("nSubjects", "nMicrobes", "spikeMicrobes", "rep")
  
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
  simlist <- SimSeq_Wrapper(simparams, simparamslabels, metadataType, workingDirectory)
  
  # Stop the Cluster 
  stopCluster(cl)
  
  # Set Names
  simnames<- paste('noRandomEffect', simparams, sep='_')
  
  names(simlist) <- simnames
  
  # Return
  return(simlist)
}

SimSeq_Wrapper <- function(simparams, simparamslabels, metadataType, workingDirectory){
  f<-foreach(i = simparams, .packages = c("SimSeq", "MASS", "stringi", 'dplyr', 'plyr'),
             .export = c("sim_data", "parse_input_data")) %dopar% {
               
               # Extract Parameter Strings
               params = strsplit(i, '_')[[1]]
               names(params) <- simparamslabels
               
               # Extract Relevant Parameters
               nSubjects <- as.numeric(params["nSubjects"])  # Number of Subjects
               nMicrobes <- as.numeric(params["nMicrobes"])  # Number of Microbes
               spikeMicrobes <- as.numeric(params["spikeMicrobes"]) # Proportion of Spiked-in Microbes
               rep <- as.numeric(params["rep"])
               
               set.seed(rep)

               return(sim_data(metadataType, nSubjects, nMicrobes, spikeMicrobes, workingDirectory))
             }
  return(f)
}

parse_input_data <- function(metadataType, workingDirectory) {
  template_folder = paste0(workingDirectory, "/general_evaluations/data_generation/simseq_templates/", metadataType, "/")
  true_tax <- read.csv(paste0(template_folder, metadataType, '_genus_table.tsv'), sep = '\t')
  true_metadata <- read.csv(paste0(template_folder, metadataType, '_metadata.tsv'), sep='\t')
  
  rownames(true_tax) <- paste0("Feature", 1:nrow(true_tax))
  true_tax <- as.matrix(true_tax[,-1])
  colnames(true_tax) <- gsub("^X", "", colnames(true_tax))
  true_tax <- true_tax[,colnames(true_tax) %in% true_metadata$sample]
  true_metadata <- true_metadata[true_metadata$sample %in% colnames(true_tax),]
  extract_vals <- mapvalues(true_metadata$sample, colnames(true_tax), 1:ncol(true_tax))
  extract_vals <- as.numeric(extract_vals)
  
  true_metadata <- true_metadata[order(true_metadata$sample),]
  true_tax <- true_tax[,order(colnames(true_tax))]

  treatment <- true_metadata$variable == names(sort(table(true_metadata$variable))[1])
  
  return(list(true_tax = true_tax, treatment = treatment))
}

sim_data <- function(metadataType, nSubjects, nMicrobes, spikeMicrobes, workingDirectory) {
  input_list <- parse_input_data(metadataType, workingDirectory)
  
  true_tax <- input_list[[1]]
  treatment <- input_list[[2]]
  
  sim_data <- SimData(true_tax,
                      treatment,
                      replic=NULL,
                      sort.method = "unpaired",
                      k.ind = as.integer(nSubjects / 2),
                      n.genes=as.integer(nMicrobes),
                      n.diff=as.integer(nMicrobes * spikeMicrobes),
                      norm.factors=apply(true_tax, 2, function(x) quantile(x[x > 0], 0.75)), #CSS
                      samp.independent=FALSE,
                      genes.select=NULL,
                      genes.diff=NULL,
                      switch.trt=FALSE,
                      probs=NULL,
                      weights=NULL,
                      exact=FALSE,
                      power=1)
  
  taxa_table <- sim_data$counts
  treatment_new <- sim_data$treatment
  
  colnames(taxa_table) <- paste0("Sample", 1:ncol(taxa_table))
  rownames(taxa_table) <- paste0("Feature", 1:nrow(taxa_table))
  
  true_diff_abun <- rownames(taxa_table)[sim_data$genes.subset %in% sim_data$DE.genes]
  
  metadata <- data.frame(group = treatment_new)
  rownames(metadata) <- paste0("Sample", 1:ncol(taxa_table))
  
  truth <- data.frame(metadata_datum = "group",
                      feature_spiked = true_diff_abun,
                      associated_property = "abundance_prevalence",
                      effect_size = NA)
  
  return(list(metadata=metadata, 
              features=taxa_table, 
              truth=truth))
}

