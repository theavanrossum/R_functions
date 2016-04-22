#'Calculates the number of features observed (presence/absence) if a 
#'subset of the total samples is considered. Repeats this calculation
#'for each subset (niter) and for a series of increasingly large subsets
#'@param inputMatrix numeric matrix with features as rows and samples as columns
#'@param subsampStep step size between subsampling levels
#'@param start number of samples to start with (minimum subsample size)
#'@param niter number of iterations to perform for each subsampling level
#'@return two column data frame with subsample size and number of features observed (column names: subsampleSize, pcCount)
featureCountRarefaction <- function(
  inputMatrix, subsampStart = 1, subsampStep = 1 , niter = 10
){
  stopifnot(is.matrix(inputMatrix) & is.numeric(inputMatrix))
  stopifnot(ncol(inputMatrix) >=1 & nrow(inputMatrix) >=1)

  mpa <- inputMatrix
  mpa[mpa != 0] <- 1 # code all abundances as present/absent
  nsamples <- ncol(mpa)

  subsampSizeVec <- c()
  pcCountVec <- c()
  
  # randomly subsample the samples (without replacement)
  subsampSize <- subsampStart
  while(subsampSize <= nsamples){
    print(paste("Sampling size:",subsampSize))
    # repeat each subsampling niter times
    for( i in 1:niter){
      
      mpaSub <- mpa[,sample(1:nsamples, subsampSize, replace=F)]
      # get the number of PCs with non-zero abundance
      pcCount <- sum(rowSums(mpaSub != 0)>0)
      
      pcCountVec <- c(pcCountVec, pcCount)
      subsampSizeVec <- c(subsampSizeVec,subsampSize)
    }
    subsampSize <- subsampSize+subsampStep
  }
  
  
  pcCountRare <- cbind.data.frame(
    subsampleSize = subsampSizeVec, 
    pcCount = pcCountVec)
  return(pcCountRare)
}