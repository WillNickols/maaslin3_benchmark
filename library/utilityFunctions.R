

#########################
# Summarization Scripts #
#########################

# Calculated Evaluation Metris for A Single Data Frame

eval_res_list = function(resi, alpha=0.05) {
  
  # Replace Missing Q-values to the Highest Possible Value 1.0
  resi[is.na(resi[, "qval"]), "qval"] <- 1
  
  # Evaluate Detection Performance
  
  time= mean(resi[,"time"], na.rm=TRUE)
  wh.pred = (resi[, "qval"] < alpha)
  wh.pos = which(wh.pred)
  wh.neg = which(!wh.pred)
  wh.TP = grep("[[:print:]]+\\_TP$", resi[, "pairwiseAssociation"])
  FPs = sum(!wh.pos %in% wh.TP)
  TPs = sum(wh.pos %in% wh.TP)
  TNs = sum(!wh.neg %in% wh.TP)
  FNs = sum(wh.neg %in% wh.TP)
  
  
  # Sensitivity: True Positives Divided by All Positives (Sum of True
  # Positives and False Negatives)
  Sensitivity = TPs/(TPs + FNs)
  
  # Specificity: True Negatives Divided by All Negatives (Sum of True
  # Negatives and False Positives)
  Specificity = TNs/(TNs + FPs)
  
  # False Discovery Rate: False Positives Divided by All Detected Positives
  FDR = if ((TPs + FPs) == 0) 
    0 else FPs/(TPs + FPs)
  # If no true positives, return NA's for irrelevant measures
  
  # FScore
  FScore<-2*TPs/(2*TPs + FNs + FPs) 
  
  # Matthew's Correlation Coefficient
  numerator <- (TPs * TNs - FPs * FNs)
  denominator <- sqrt((TPs + FPs)*(TPs + FNs)*(TNs + FPs)*(TNs + FNs))
  if(denominator == 0) denominator <- 1
  MCC <- numerator/denominator
  
  # AUC, pAUC (FPR < 0.20)
  wh.truth = (1:nrow(resi) %in% wh.TP)
  if (all(wh.truth=='TRUE')) {
    AUC=1
    pAUC=1
    fAUC=1
  } else if (all(wh.truth=='FALSE')) {
    AUC=0
    pAUC=0
  } else {
    pred <- prediction(as.numeric(wh.pred), factor(wh.truth, levels=c("TRUE", "FALSE")))
    AUC = performance(pred, "auc")@y.values[[1]]
    pAUC = performance(pred, "auc", fpr.stop=0.20)@y.values[[1]][[1]]
  }
  
  # Departure from Uniformity Under the Null
  AucAocVals <- AucAocFun(resi[!wh.truth, "pval"], plotIt = FALSE, pch = 20, type = "l")
  totalArea = AucAocVals["conservArea"] + AucAocVals["liberalArea"]
  names(totalArea) = "totalArea"  
  
  # Return
  return(c(Sensitivity = Sensitivity,
           Specificity=Specificity,
           FDR=FDR,
           FScore=FScore,
           MCC=MCC,
           AUC=AUC,
           pAUC=pAUC,
           AucAocVals["conservArea"],
           AucAocVals["liberalArea"],
           totalArea,
           time=time))
}

# Extract Parameters from A List of Results

make_power_df = function(reslist, simparamslabels) {
  powerdf = ldply(reslist, .parallel=TRUE)
  colnames(powerdf)[1] <- "Combinations"
  paramdf = ldply(strsplit(powerdf[, "Combinations"], "_"), .parallel=TRUE)
  colnames(paramdf) <- simparamslabels
  powerdf = cbind(powerdf, paramdf)
  return(powerdf)
}


# Organize All Lists into A Coherent Data Frame

list.perfdf<-function(resultslist, simparams, simparamslabels){
  foreach(resultslist = resultslist, .export=c("make_power_df", "eval_res_list", "AucAocFun"),
          .packages = c("ROCR", "plyr"), .errorhandling = 'remove') %dopar% {
  perflist <- lapply(resultslist, eval_res_list)
  if (is.null(names(resultslist))) {
    names(perflist) <- simparams
  } else {
    names(perflist) <- names(resultslist)
  }
  perfdf = make_power_df(perflist, simparamslabels)
  return(perfdf)
}}



#######################################
# Arc-Sine Square Root Transformation #
#######################################

# Arc Sine Square Root Transformation
AST<-function(x){
  return(sign(x)*asin(sqrt(abs(x))))
}

########################
# Logit Transformation #
########################

# LOGIT Transformation 
LOGIT<-function(x){
  y<-car::logit(x, adjust=0)
  y[!is.finite(y)]<-0
  return(y)
}

# Shifted Logit Transformation (Lukens et al, 2014, Nature)
# LOGIT_S<-function(x){
#   y<-0.5*log(x/(1-x)) + 10
#   y[!is.finite(y)]<-0
#   return(y)
# }

######################
# Log Transformation #
######################

# LOG Transformation
LOG<-function(x){
  return(log(x+1))
}

# Apply SimpleFilter (Prevalence and Abundance Filtering) To A Dataset

SimpleFilter = function(physeq, Threshold_Abundance, Threshold_Prevalence) {
  features<-physeq$features
  filtered_features<-features[,colSums(features > Threshold_Abundance) > nrow(features)*Threshold_Prevalence] 
  physeq$features<-filtered_features
  return(physeq)
}

# Apply SimpleFilter (Prevalence and Abundance Filtering) To A List of Datasets

list.SimpleFilter<- function(simlist, Threshold_Abundance, Threshold_Prevalence){
  return(lapply(simlist, SimpleFilter, Threshold_Abundance, Threshold_Prevalence))
}

# Calculate No. of TP Features in A Dataset

TruePositives<-function(physeq, breakitdown) {
  features<-physeq$features
  d<-colnames(features)[grep("[[:print:]]+\\_TP$", colnames(features))]
  if(breakitdown==TRUE) {
    return(d)
  } else {
    return(length(d))
  }
}

# Calculate No. of TP Features in A List of Datasets

list.TruePositives<- function(simlist, breakitdown){
  return(lapply(simlist, TruePositives, breakitdown))
}

# Calculate Prevalence of All Features in A Dataset

calculatePrevalence<-function(physeq) {
  features<-physeq$features
  d<-colSums(features > Threshold_Abundance)/nrow(features)
  return(d)
}

# Calculate Prevalence of All Features in A List of Datasets

list.calculatePrevalence<- function(simlist){
  return(lapply(simlist, calculatePrevalence))
}

# Check if a metadata is binary 

is.binary <- function(v) {
  x <- unique(v)
  dim(x)[1] - sum(is.na(x)) == 2L
}


# Function from Hawinkel et al. (2017)
## A function to look at the distribution of the P-values under the null
## distribution 

### Compute AUC and AOC of P-value distribution to be applied on the Monte
### Carlo distribution of each OTU separately
AucAocFun <- function(pVals, maxVal = 0.25, plotIt = FALSE, ...) {
  ## sum of height differences between the curve and the y=x line, half maximum
  ## = 0.5 * max(pVals) * length(pVals) this is the case of 100% rejection with
  ## any alpha value _onlyTotArea_ controls if only the total area is returned,
  ## and not also the Area Over and Under the y=x line.  halfMaxArea <- 0.5 *
  ## max(maxVal) * length(pVals)
  halfMaxArea <- 0.5 * length(pVals)
  
  pVals <- pVals[!is.na(pVals)]
  estimated <- sort(pVals)
  theoretic <- seq_along(pVals)/length(pVals)
  
  if (plotIt) {
    plot(theoretic, estimated, ...)
    abline(0, 1)
  } else {
  }
  
  diffPVals <- theoretic - estimated
  indConserv <- theoretic <= estimated
  conservArea <- sum(-diffPVals[indConserv])/halfMaxArea
  liberalArea <- sum(diffPVals[!indConserv])/halfMaxArea
  
  c(conservArea = conservArea, liberalArea = liberalArea, totalArea = liberalArea + 
      conservArea)
  
}  # END - function: aucAocFun, AUC/AOC calculation for p-values distribution
