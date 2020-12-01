makeXvX_gam2 <- function(covData = NA , dIData , listfefdKnots , incInts = TRUE , intMthd = 0 , colnamesXcns = NULL , nDiscPts = 100 , lnTfmdData = FALSE){
  if(identical(covData , NA)){ 
    nCovs <- 0  
    p <- 1
    return(list('X' = matrix(1 , nrow(dIData) , 1) , 'vXU' = NA , 'iU' = NA , 'namesX' = 'const'))     
  }else{}
  
  nCovs <- ncol(covData)
  n <- nrow(covData)
  
  for(nm in names(listfefdKnots)){
    if(!is.element(nm , names(covData))){ stop(paste0('Error - ' , nm , ' not found in covariate data!')) }else{}
  }
  if(!is.element('dIMidPts' , names(covData))){ print(paste0('Attention - dIMidPts not found in covariate data or listfefdKnots, so trend function will not vary depth! Just checking that this is what you want')) }else{}

  if(is.null(colnamesXcns)){ stop('Error - enter colnamesXcns for makeXvX_gam2 function!') }else{}
  
  p <- length(colnamesXcns)
  colsd <- which(grepl('dIMidPts_KNOT' , colnamesXcns))
  colsnod <- setdiff(seq(p) , colsd)

  if(lnTfmdData & (length(colsd) > 0)){
    
    maxTmp <- 1000000
    nPerBatch <- floor(maxTmp / (nDiscPts * p * p))
    if(nPerBatch < 1){ nPerBatch <- 1 }else if(nPerBatch > n){ nPerBatch <- n }else{}
    nBatches <- ceiling(n / nPerBatch)
    
    fnVarApply <- function(idx , X4Apply){ var(X4Apply[idx,,drop=FALSE]) }
    fnMeanApply <- function(idx , X4Apply){ colMeans(X4Apply[idx,,drop=FALSE]) }
    
### initialise mats...    
    Xcns <- matrix(NA , n , p)
    vXUcns <- matrix(NA , n * p , p)

    ########################################################
    ### loop over each batch, and calc vXU using disc pts...
    ########################################################
    for (iBatch in 1:nBatches){
      
      if(iBatch < nBatches){ 
        i <- seq((iBatch-1)*nPerBatch+1 , iBatch*nPerBatch) 
      }else{ 
        i <- seq((iBatch-1)*nPerBatch+1 , n) 
      }
      
      if(iBatch == 1){
        i4Apply <- matrix(seq(nDiscPts * length(i)) , nDiscPts , length(i))
      }else if(iBatch == nBatches){
        i4Apply <- matrix(seq(nDiscPts * length(i)) , nDiscPts , length(i))
      }else{} # in between 1 and nBatches, no need to update it. 
      
      linTmp <- seq(0 , 1 , 1 / (nDiscPts-1))
      iRep <- rep(i , each = nDiscPts)
      
      if(nDiscPts > 1){
        dDiscPts <- dIData[iRep,1] + rep(linTmp , length(i)) * (dIData[iRep,2] - dIData[iRep,1])
      }else{
        dDiscPts <- 0.5 * (dIData[i,1] + dIData[i,2])
      }          
      
      covDataThis <- covData[iRep,,drop=FALSE]
      covDataThis$dIMidPts <- dDiscPts

      ### if we don't pass in dIData, then makeXcns will use the given dIMidPts on pt support... 
      XTmp <- makeXcns(dfCovs = covData , listfefdKnots = listfefdKnots , incInts = incInts , intMthd = intMthd , colnamesX = colnamesXcns)
      
      mXTmp <- t(apply(i4Apply , 2 , fnMeanApply , X4Apply = XTmp))

      Xcns[i,] <- as.matrix(mXTmp)
      
      ### remove means from XTmp...
      XTmp <- XTmp - Xcns[iRep,,drop=FALSE] 
      
      iThis <- rep((i-1)*p , each = p) + rep(seq(1,p) , length(i))
      
      vXTmp <- apply(i4Apply , 2 , fnVarApply , X4Apply = XTmp)
      vXTmp <- matrix(vXTmp , length(i) * p , p , byrow = TRUE)
      vXUcns[iThis,] <- vXTmp
    }
  
    iU <- colsd  
    iK <- colsnod
  ###############################################
  ### remove any variances from vXU for cols without d...
  ###############################################
    if(length(iK) > 0){
      vXUcns <- matrix(vXUcns[,-iK] , nrow = n * p)
      iRemove <- kronecker(seq(0 , n - 1) * p , matrix(1 , length(iK) , 1)) + kronecker(matrix(1 , n , 1) , iK)        
      vXUcns <- vXUcns[-iRemove,,drop=FALSE]
    }else{}    

    ### just to make sure things are matrices...    
    Xcns <- matrix(Xcns , nrow = n , ncol = p)
    vXUcns <- matrix(vXUcns , nrow = n * length(iU) , ncol = length(iU))

  }else{
    Xcns <- makeXcns(dfCovs = covData , dIData = dIData , listfefdKnots = listfefdKnots , incInts = incInts , intMthd = intMthd , colnamesX = colnamesXcns)
    vXUcns <- NA
    iU <- NA
  }
  
  return(list('X' = Xcns , 'vXU' = vXUcns , 'iU' = iU , 'namesX' = colnames(Xcns)))     
}
    