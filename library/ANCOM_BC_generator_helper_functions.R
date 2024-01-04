trigger_ANCOM_BC_Simulator <- function(metadataType=metadataType,
                                       nSubjects=nSubjects,
                                       nMicrobes=nMicrobes,
                                       spikeMicrobes = spikeMicrobes, 
                                       effectSize = effectSize,
                                       effectPos = effectPos,
                                       readDepth = readDepth,
                                       nIterations=nIterations,
                                       nCores=nCores) {
  # Create Replicates 
  reps = 1:nIterations
  
  ########################
  # Catch Obvious Errors #
  ########################
  
  # Check Character Values
  if (!metadataType %in% c('pois_gam', 'soil'))
    stop('Must be one of the following: pois_gam, soil')
  
  # Check Positive Integer Values
  if (round(nSubjects) != nSubjects || 
      nSubjects<0 ||
      round(nMicrobes) != nMicrobes || 
      nMicrobes<0 ||
      round(readDepth) != readDepth || 
      readDepth<0)
    stop('nSubjects/nMicrobes/readDepth must be positive integers.')
  
  # Check Proportion Values
  if (spikeMicrobes>1 || spikeMicrobes<=0)
    stop('spikeMicrobes must be in (0, 1].')
  
  # Define the Simulation Parameters Combinations
  simparams = apply(expand.grid(nSubjects,
                                nMicrobes,
                                spikeMicrobes, 
                                effectSize,
                                effectPos,
                                readDepth,
                                reps), 1, paste, collapse = '_')
  
  # Define the Labels to Go with Each Element of the Simulation Parameter
  simparamslabels = c("nSubjects", "nMicrobes", "spikeMicrobes", "effectSize", "effectPos", "readDepth", "rep")
  
  # Set Up Clustering Environment
  no_cores <- nCores 
  cl <- makeCluster(no_cores)
  registerDoParallel(cl)
  
  ####################
  # Data Generation #
  ###################
  
  # Call ANCOM-BC Simulator 
  simlist <- ANCOM_BC_Simulator_wrapper(simparams, simparamslabels, metadataType)

  # Stop the Cluster 
  stopCluster(cl)
  
  # Set Names
  simnames<- paste('noRandomEffect', simparams, sep='_')
  names(simlist) <- simnames
  
  # Return
  return(simlist)
}

ANCOM_BC_Simulator_wrapper <- function(simparams, simparamslabels, metadataType) {
  f<-foreach(i = simparams, .packages = c("phyloseq", "microbiome", "MASS", "stringi", 'dplyr'), 
             .export = c("abn.tab.gen.soil", "abn.tab.gen.pois.gam")) %dopar% {
              # Extract Parameter Strings
              params = strsplit(i, '_')[[1]]
              names(params) <- simparamslabels
              
              # Extract Relevant Parameters
              n.samp.grp1 <- ceiling(as.numeric(params["nSubjects"]) / 2)  # Number of Subjects
              n.samp.grp2 <- floor(as.numeric(params["nSubjects"]) / 2)  # Number of Subjects
              n.taxa <- as.numeric(params["nMicrobes"])  # Number of Microbes
              prop.diff <- as.numeric(params["spikeMicrobes"]) # Proportion of Spiked-in Microbes
              effectSize<-as.character(params["effectSize"]) # Effect Size
              effectPos<-as.numeric(params["effectPos"]) # Effect Size
              readDepth<-as.numeric(params["readDepth"]) # Library Size
              abn.seed = as.numeric(params["rep"])
              obs.seed = as.numeric(params["rep"]) + 1
              struc.zero.prop=0.2
              
              if (metadataType == "pois_gam") {
                single_sim <- abn.tab.gen.pois.gam(n.taxa, 
                                                n.samp.grp1, 
                                                n.samp.grp2, 
                                                prop.diff, 
                                                abn.seed, 
                                                obs.seed,
                                                effectPos,
                                                effectSize,
                                                readDepth)
              } else {
                single_sim <- abn.tab.gen.soil(n.taxa, 
                                            n.samp.grp1, 
                                            n.samp.grp2, 
                                            prop.diff, 
                                            abn.seed, 
                                            obs.seed,
                                            effectPos,
                                            effectSize,
                                            readDepth)
              }
              metadata = data.frame(single_sim$grp)
              colnames(metadata) <- "group"
              features <- single_sim$obs.abn
              rownames(features) <- paste0("Feature", rownames(features))
              truth = data.frame("metadata_datum" = c("group"), 
                                 "feature_spiked" = paste0("Feature", names(single_sim$effect.size)),
                                 "associated_property" = c("abundance"),
                                 "effect_size" = unname(single_sim$effect.size))
              return(list(metadata=metadata, 
                          features=features, 
                          truth=truth))
  }
  return(f)
}

# Soil
abn.tab.gen.soil = function(n.taxa, 
                            n.samp.grp1, 
                            n.samp.grp2, 
                            prop.diff, 
                            abn.seed, 
                            obs.seed,
                            effectPos,
                            effectSize,
                            readDepth,
                            struc.zero.prop=0.2){
  data("GlobalPatterns")
  pseq = GlobalPatterns
  # Simulations were evaluated for soil environments
  meta.data = meta(pseq)
  pseq.subset = subset_samples(pseq, SampleType == "Soil")
  # Prune taxa
  pseq.prune = prune_taxa(taxa_sums(pseq.subset) > 50, pseq.subset)
  template = taxa_sums(pseq.prune) * ceiling(readDepth / 10000) # Ensure sampling still works
  
  # Template for absolute abundance in the ecosystem
  set.seed(abn.seed)
  abn.temp = sample(template, n.taxa)
  taxa.id = names(abn.temp); n.samp = n.samp.grp1 + n.samp.grp2
  
  # Which taxa are differentially abundant
  diff.ind=rep(0, n.taxa)
  # Group1 is higher than group2
  diff1.ind=sample(c(1:n.taxa), floor(n.taxa*prop.diff), replace=FALSE)
  diff.ind[diff1.ind]=1
  # Group2 is higher than group1
  wt=1-effectPos
  diff2.ind=sample(diff1.ind, wt*length(diff1.ind), replace=FALSE)
  diff.ind[diff2.ind]=2
  # Structural zeros
  diff3.ind=sample(which(diff.ind==0), struc.zero.prop*length(which(diff.ind==0)), replace = FALSE)
  diff.ind[diff3.ind]=-1
  
  # Effect size
  effect.size=rep(1, n.taxa)
  effect.size[diff1.ind]=2^(runif(length(diff1.ind), as.numeric(effectSize) * 0.5, as.numeric(effectSize)))
  effect.size[diff2.ind]=2^(-runif(length(diff2.ind), as.numeric(effectSize) * 0.5, as.numeric(effectSize)))
  effect.size[diff3.ind]=0
  names(effect.size)=taxa.id
  
  # Absolute abundance in the ecosystem of two groups
  abn.grp1 = round(abn.temp * effect.size)
  abn.grp2 = round(abn.temp)
  abn.mat = cbind(abn.grp1, abn.grp2)
  
  abn.total = colSums(abn.mat)
  names(abn.total) = c("grp1", "grp2")
  
  mean_depth <- readDepth / max(abn.total)
  mult_fact <- 1 / (mean_depth * 2 / log(5) * 4 / 11 * 10)
  # Match means to original set-up
  #depth=1/sample(c(runif(n.samp, 10, 50), runif(n.samp, 100, 500)), n.samp, replace = T)
  depth=1/sample(c(runif(n.samp, mult_fact, 5 * mult_fact), runif(n.samp, 10 * mult_fact, 50 * mult_fact)), 
                 n.samp, replace = T)
  obs.total=round(max(abn.total)*depth)
  names(obs.total) = paste0("sub", seq(n.samp))
  
  # Specimen abundance
  set.seed(obs.seed)
  obs.mat = matrix(NA, nrow = n.taxa, ncol  =  n.samp)
  for (i in 1:n.samp.grp1) {
    obs.mat[, i] = phyloseq:::rarefaction_subsample(x = abn.mat[, 1], sample.size = obs.total[i])
  }
  for (i in (n.samp.grp1 + 1):n.samp) {
    obs.mat[, i] = phyloseq:::rarefaction_subsample(x = abn.mat[, 2], sample.size = obs.total[i])
  }
  
  # Prepare output data sets
  abn.dat = data.frame(abn.mat, row.names  =  NULL)
  rownames(abn.dat) = taxa.id
  colnames(abn.dat) = c("grp1", "grp2")
  
  obs.dat = data.frame(obs.mat, row.names  =  NULL)
  rownames(obs.dat) = taxa.id
  colnames(obs.dat) = paste0("sub", seq(n.samp))
  
  grp.ind = c(rep(1, n.samp.grp1), rep(2, n.samp.grp2))
  names(grp.ind) = paste0("sub", seq(n.samp))
  
  names(diff.ind) = taxa.id
  
  c.mult = c(obs.total[1:n.samp.grp1]/abn.total[1], 
             obs.total[(n.samp.grp1 + 1):n.samp]/abn.total[2])
  names(c.mult) = paste0("sub", seq(n.samp))
  
  test.data = list(abn.dat, obs.dat, effect.size, grp.ind, 
                   diff.ind, c.mult, abn.total, obs.total)
  names(test.data) = c("pop.abn", "obs.abn", "effect.size", "grp", 
                       "diff.taxa", "mult", "abn.total", "obs.total")
  return(test.data)
}

# Poisson-gamma two group
abn.tab.gen.pois.gam=function(n.taxa, 
                              n.samp.grp1, 
                              n.samp.grp2, 
                              prop.diff, 
                              abn.seed, 
                              obs.seed,
                              effectPos,
                              effectSize,
                              readDepth,
                              struc.zero.prop=0.2) {
  # From sim_high_var
  low.abn=50 * ceiling(readDepth / 50000)
  med.abn=200 * ceiling(readDepth / 50000)
  high.abn=10000 * ceiling(readDepth / 50000)
  out.zero.prop=0.05
  
  # Total number of samples
  n.samp=n.samp.grp1+n.samp.grp2
  
  set.seed(abn.seed)
  low.prop=0.6 # Proportion of low abundance 
  med.prop=0.3 # Proportion of medium abundance
  hi.prop=0.1  # Proportion of high abundance
  # Indices for taxa abundance 
  index=sample(c(1, 2, 3), n.taxa, replace = T, prob = c(low.prop, med.prop, hi.prop)) 
  
  # Poisson parameters
  lambda=rep(NA, n.taxa)
  lambda[which(index==1)]=rgamma(length(which(index==1)), shape=low.abn, rate=1)
  lambda[which(index==2)]=rgamma(length(which(index==2)), shape=med.abn, rate=1)
  lambda[which(index==3)]=rgamma(length(which(index==3)), shape=high.abn, rate=1)
  
  # Construct unbalanced microbial in the ecosystem
  
  # Which taxa are differentially abundant
  diff.ind=rep(0, n.taxa)
  # Group1 is higher than group2
  diff1.ind=sample(c(1:n.taxa), floor(n.taxa*prop.diff), replace=FALSE)
  diff.ind[diff1.ind]=1
  # Group2 is higher than group1
  wt=1-effectPos
  diff2.ind=sample(diff1.ind, wt*length(diff1.ind), replace=FALSE)
  diff.ind[diff2.ind]=2
  # Structural zeros
  diff3.ind=sample(which(diff.ind==0), struc.zero.prop*length(which(diff.ind==0)), replace = FALSE)
  diff.ind[diff3.ind]=-1
  
  # Effect size
  effect.size=rep(1, n.taxa)
  effect.size[diff1.ind]=2^(runif(length(diff1.ind), as.numeric(effectSize) * 0.5, as.numeric(effectSize)))
  effect.size[diff2.ind]=2^(-runif(length(diff2.ind), as.numeric(effectSize) * 0.5, as.numeric(effectSize)))
  effect.size[diff3.ind]=0
  names(effect.size)=paste0("taxon", seq(n.taxa))
  
  # Mean absolute abundance in the ecosystem
  temp.grp1=round(lambda*effect.size)
  temp.grp2=round(lambda)
  for (i in which(effect.size!=1)) {
    if(temp.grp1[i]==temp.grp2[i]) temp.grp1[i]=temp.grp1[i]+1
  }
  temp.dat=data.frame(temp.grp1, temp.grp2, effect.size)
  rownames(temp.dat)=paste0("taxon", seq(n.taxa))
  
  # Absolute abundance in ecosystems
  abn.mat=matrix(0, ncol=n.samp, nrow=n.taxa)
  for(i in 1:n.taxa){
    abn.mat[i, ]=c(rpois(n.samp.grp1, temp.grp1[i]), rpois(n.samp.grp2, temp.grp2[i]))
  }
  # Outlier zeros
  out.ind=rep(0, n.taxa); out.ind[sample(seq(n.taxa), out.zero.prop*n.taxa, replace = F)]=1
  names(out.ind)=paste0("taxon", seq(n.taxa))
  abn.mat[which(out.ind==1), sample(seq(n.samp), out.zero.prop*n.samp, replace = F)]=0
  
  # Microbial load
  abn.total=colSums(abn.mat)
  names(abn.total)=paste0("sub", seq(n.samp))
  abn.prob.mat=t(t(abn.mat)/abn.total)
  
  # library size
  mean_depth <- readDepth / max(abn.total)
  mult_fact <- 1 / (mean_depth * 2 / log(5) * 4 / 11 * 10)
  # Match means to original set-up
  #depth=1/sample(c(runif(n.samp, 10, 50), runif(n.samp, 100, 500)), n.samp, replace = T)
  depth=1/sample(c(runif(n.samp, mult_fact, 5 * mult_fact), runif(n.samp, 10 * mult_fact, 50 * mult_fact)), 
                 n.samp, replace = T)
  obs.total=round(max(abn.total)*depth)
  names(obs.total)=paste0("sub", seq(n.samp))
  
  # Absolute abundance in samples
  set.seed(obs.seed)
  obs.list=lapply(1:n.samp, function(i)
    phyloseq:::rarefaction_subsample(x=abn.mat[, i], sample.size=obs.total[i]))
  obs.mat=Reduce('cbind', obs.list)
  
  # Prepare outputs
  abn.dat=data.frame(abn.mat, row.names = NULL)
  rownames(abn.dat)=paste0("taxon", seq(n.taxa))
  colnames(abn.dat)=paste0("sub", seq(n.samp))
  
  obs.dat=data.frame(obs.mat, row.names = NULL)
  rownames(obs.dat)=paste0("taxon", seq(n.taxa))
  colnames(obs.dat)=paste0("sub", seq(n.samp))
  
  abn.prob.dat=data.frame(abn.prob.mat, row.names = NULL)
  rownames(abn.prob.dat)=paste0("taxon", seq(n.taxa))
  colnames(abn.prob.dat)=paste0("sub", seq(n.samp))
  
  grp.ind=c(rep(1, n.samp.grp1), rep(2, n.samp.grp2))
  names(grp.ind)=paste0("sub", seq(n.samp))
  
  names(diff.ind)=paste0("taxon", seq(n.taxa))
  
  # Sampling fractions
  c.mult=obs.total/abn.total
  names(c.mult)=paste0("sub", seq(n.samp))
  
  test.data=list(temp.dat, abn.dat, obs.dat, effect.size, grp.ind, 
                 diff.ind, out.ind, c.mult, abn.total, obs.total)
  names(test.data)=c("mean.eco.abn", "eco.abn", "obs.abn", "effect.size", "grp", 
                     "diff.taxa", "outlier", "samp.frac", "micro.load", "lib.size")
  return(test.data)
}
