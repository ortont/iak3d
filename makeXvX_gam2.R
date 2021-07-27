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
    
#################################################################################################
### fn to wrap up all of the gam2 stuff needed before the fitIAK3D fn is called. 
#################################################################################################
gam2IAK3DInit <- function(dIFit , covsFit , scaleCovs = TRUE , nIntKnotsd = 1 , nIntKnotss = 1 , incInts = NULL , intMthd = 0 , 
                          q4BdryKnotsd = c(0.05 , 0.95) , q4BdryKnotss = c(0.01 , 0.99)){
  ###################################################################################
  ### set up for fitting a spline model.
  ### nIntKnotsd : number of internal knots for the spline function (nat spline, clamped to have grad=0 at upper bdry) of depth; if this is complex enough, probably no need for the depth component in prod-sum covariance model
  ### nIntKnotss : number of internal knots for the spline functions (nat spline, clamped to have grad=0 at upper and lower bdries) of covariates
  ### incIntsWithd : include interactions between depth and spatial covariates (but here not between different spatial covariates)
  #####################
  ### if incInts = TRUE, then include all interactions of basis fns.
  ### if incInts = list of two vectors, then all variables in first vector interact with all variables in the second vector
  #####################
  ### if variable v has basis fns v1 - v5 and w has w1 - w5
  ### intMthd = 0 : interactions between v and w are v1*w0, v2*w0, ..., v5*w0, v0*w1, v0*w2,..., v0*w5 (ie 10 basis fns)    
  ### where v0 is spline type basis of v with 1 df, same for w0
  ### intMthd = 1 : interactions between v and w are all vi * wj (ie 25 basis fns)
  #####################
  ### q4BdryKnotsd are the quantiles of the dIMidPts data that are used to define the boundary knots
  ### q4BdryKnotsd are the quantiles of the spatial covariates data that are used to define the boundary knots
  ### if these are outside 0-1, then they define an extension of the range (-0.25 = min - 0.25 * range; 1.25 = max + (1.25 - 1) * range)
  ### at present, only one q4BdryKnotss applies to all spatial covariates. should generalise at some stage.
  ###################################################################################
  
  if((!is.numeric(q4BdryKnotsd)) | (length(q4BdryKnotsd) != 2)){ stop('Error - q4BdryKnotsd must be numeric vector of length 2!') }else{}
  if((!is.numeric(q4BdryKnotss)) | (length(q4BdryKnotss) != 2)){ stop('Error - q4BdryKnotss must be numeric vector of length 2!') }else{}

  modelX <- list('type' = 'gam2')
  
  ### don't include depth here.   
  spatialCovs <- setdiff(names(covsFit) , 'dIMidPts')
  
  iCat <- iCts <- c()
  if(length(spatialCovs) > 0){
    for(inm in 1:length(spatialCovs)){
      if(is.factor(covsFit[[spatialCovs[inm]]])){
        iCat <- c(iCat , inm)
      }else if(is.numeric(covsFit[[spatialCovs[inm]]])){
        iCts <- c(iCts , inm)
      }else{
        stop('Unknown data type detected in covsFit!')
      }
    }
  }else{}
  
  if(scaleCovs & (length(iCts) > 0)){

    ### to work with scaled covariates
    if(is.list(incInts) && length(incInts) == 2){
      for(i in 1:2){
        for(j in iCts){
          incInts[[i]] <- gsub(spatialCovs[j] , paste0(spatialCovs[j] , '_SCALED') , incInts[[i]])
        }
      }
    }else{}

    spatialCovs[iCts] <- paste0(spatialCovs[iCts] , '_SCALED')

  }else{
    ### to work with unscaled (raw) covariates
    spatialCovs <- spatialCovs # no change here
  } 

  ### add any scaled variables (_SCALED) to covs dfs...
  if(length(iCts) > 0){
    # tmp <- addScaledCovs2df(dfFit = covsFit , dfPred = covsVal , covNames = spatialCovs)
    # covsFit <- tmp$dfFit
    # covsVal <- tmp$dfPred
    
    tmp <- addScaledCovs2df(dfFit = covsFit , covNames = spatialCovs)
    covsFit <- tmp$dfFit
    modelX$scalePars_m <- tmp$scalePars_m
    modelX$scalePars_sd <- tmp$scalePars_sd
  }else{
    modelX['scalePars_m'] <- list(NULL)
    modelX['scalePars_sd'] <- list(NULL)
  }
  
  # for bdry knot positions for depth fn, use 5th ptile (not clamped) to stop it being too wiggly at the surface.   
  # and 95th ptile for other boundary - this is clamped, which will force a plateau
  # for other (spatial) covariates, boundary knots at 1st and 99th ptiles, which together with clamping at both ends will stop extrapolation.
  # q4BdryKnots <- c(0.05 , rep(0.01 , length(spatialCovs))) 
  # q4BdryKnots <- cbind(q4BdryKnots , 1 - q4BdryKnots)
  
  q4BdryKnots <- cbind(c(q4BdryKnotsd[1] , rep(q4BdryKnotss[1] , length(spatialCovs))) , 
                       c(q4BdryKnotsd[2] , rep(q4BdryKnotss[2] , length(spatialCovs))))

  nIntKnots <- c(nIntKnotsd , rep(nIntKnotss , length(spatialCovs)))
  
  sTypeSpatialCovs <- rep('nsclug' , length(spatialCovs))
  if(length(iCat) > 0){ sTypeSpatialCovs[iCat] <- 'cat' }else{}
  sType <- c('nscug' , sTypeSpatialCovs)
  
  modelX$listfefdKnots <- makelistfefdKnots(dfFit = covsFit , covNames = c('dIMidPts' , spatialCovs) , nIntKnots = nIntKnots , q4BdryKnots = q4BdryKnots , sType = sType)
  if(is.null(incInts)){
    ### default is to include interactions between depth and the other spatial covariates...  
    modelX$incInts <- list('dIMidPts' , spatialCovs)
  }else{
    ### incInts should be passed in as list with two sets of covariates that interact.    
    modelX$incInts <- incInts 
  }

  modelX$intMthd <- intMthd
  
  ### just to get names of columns...  
  XcnsTmp <- makeXcns(dfCovs = covsFit , dIData = dIFit , listfefdKnots = modelX$listfefdKnots , incInts = modelX$incInts , colnamesX = NULL , intMthd = modelX$intMthd) # intMthd = 1 for now. 0 = simpler
  
  modelX$colnamesX <- colnames(XcnsTmp)
  
  # return(list('modelX' = modelX , 'covsFit' = covsFit , 'covsVal' = covsVal))
  return(list('modelX' = modelX , 'covsFit' = covsFit))
}
  