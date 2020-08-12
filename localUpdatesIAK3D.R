###################################################################
### for local updating, need to redefine some matrices that are in the fitted object...
### if removeLocal, data coinciding with local data are removed from lmmFit
### else, are local data are added to lmmFit
###################################################################
lmmUpdateLocal <- function(lmmFit , xLocal = NULL , dILocal = NULL , zLocal = NULL , covsLocal = NULL , removeLocalData = FALSE , attachBigMats = TRUE , updateLmmPars = TRUE , mindFromLegData = 0.001){

#  mindFromLegData is for checking local data aren't already in legacy data. (0.001 = 1 m in km)
  
### to update stored variables from older versions of code
  lmmFit <- updateSplineNames(lmmFit)
  lmmFit <- updateCubistModel(lmmFit)

  if(!is.null(xLocal)){
    iKeep <- which((!is.na(zLocal)) & (rowSums(is.na(covsLocal)) == 0) & (apply(xyDist(lmmFit$xData , xLocal) , 2 , min) > mindFromLegData)) # last one is a check that not already in leg data   
    
    if(length(iKeep) == 0){ return(lmmFit) }else{}

    xLocal <- xLocal[iKeep,,drop=FALSE]
    dILocal <- dILocal[iKeep,,drop=FALSE]
    zLocal <- zLocal[iKeep]
    covsLocal <- covsLocal[iKeep,,drop=FALSE]
  
    if(removeLocalData){
### check if any local data are already in legacy data, and if so, remove from legacy...    
      iRm <- which(apply(xyDist(lmmFit$xData , xLocal) , 1 , min) < 0.001) 
      if(length(iRm) == 0){ return(lmmFit) }else{}
      lmmFit$xData <- lmmFit$xData[-iRm,,drop=FALSE]
      lmmFit$dIData <- lmmFit$dIData[-iRm,,drop=FALSE]
      lmmFit$zData <- lmmFit$zData[-iRm]
      lmmFit$covsData <- lmmFit$covsData[-iRm,,drop=FALSE]
### only do for new data...
      if (!is.element('iU' , names(lmmFit))){
        lmmFit$iU <- NA
      }else{}
      
      XLims <- matrix(NA , 2 , ncol(lmmFit$XData))
      XLims[1,] <- Inf
      XLims[2,] <- -Inf
 
      lmmFit$XData <- lmmFit$XData[-iRm,,drop=FALSE]
      if (is.element('vXU' , names(lmmFit))){
        iRmvXU <- rows2vXURows(iRm , iUUpdate)
        lmmFit$vXU <- lmmFit$vXU[-iRmvXU,,drop=FALSE]
      }else{
        lmmFit$vXU <- list(NULL)
      }
      
    }else{
### add the local data to the legacy data...    
      lmmFit$xData <- rbind(lmmFit$xData , xLocal)
      lmmFit$dIData <- rbind(lmmFit$dIData , dILocal)
      lmmFit$zData <- c(lmmFit$zData , zLocal)
      lmmFit$covsData <- rbind(lmmFit$covsData , covsLocal)
### only do for new data...
      if (!is.element('iU' , names(lmmFit))){
        lmmFit$iU <- NA
      }else{}
      
      XLims <- matrix(NA , 2 , ncol(lmmFit$XData))
      XLims[1,] <- Inf
      XLims[2,] <- -Inf
 
      tmp <- makeXvX(covData = covsLocal , dIData = dILocal , modelX = lmmFit$modelX , allKnotsd = lmmFit$allKnotsd , iU = lmmFit$iU , nDiscPts = 1000 , lnTfmdData = lmmFit$lnTfmdData , XLims = XLims)
      lmmFit$XData <- rbind(lmmFit$XData , tmp$X)
      if (is.element('vXU' , names(lmmFit))){
        lmmFit$vXU <- rbind(lmmFit$vXU , tmp$vXU)
      }else{
        lmmFit$vXU <- list(NULL)
      }
    }

###################################
### update the compLikMats...  
###################################
    if(lmmFit$compLikMats$compLikOptn > 0){
      vcentres <- matrix(NA , length(lmmFit$compLikMats$listBlocks) , ncol(lmmFit$compLikMats$listBlocks[[1]]$centroid))
      for (i in 1:length(lmmFit$compLikMats$listBlocks)){ vcentres[i,] <- lmmFit$compLikMats$listBlocks[[i]]$centroid }
      optnTmp <- lmmFit$compLikMats$compLikOptn
      lmmFit$compLikMats <- NULL # delete old vs.
      lmmFit$compLikMats <- setVoronoiBlocks(x = lmmFit$xData , nPerBlock = NA , plotVor = F , vcentres = vcentres)
      lmmFit$compLikMats$compLikOptn <- optnTmp
    }else{}

###################################
### update the setupMats...  
### could be speeded up in future to just update the matrices rather than remake them.
###################################
    lmmFit$setupMats <- NULL # delete any old vs.
    
    if (lmmFit$compLikMats$compLikOptn == 0){
      if(attachBigMats){
        lmmFit$setupMats <- setupIAK3D(lmmFit$xData , lmmFit$dIData , nDscPts = 0)
      }else{
        lmmFit$setupMats <- list('xData' = lmmFit$xData , 'dIData' = lmmFit$dIData)
      }
    }else{
      
      lmmFit$setupMats <- list()
      ### order is all adj subset pairs, then all individual subsets (which can be used to get non-adj subset pairs)
      for (i in 1:nrow(lmmFit$compLikMats$subsetPairsAdj)){
        iThis <- c(lmmFit$compLikMats$listBlocks[[lmmFit$compLikMats$subsetPairsAdj[i,1]]]$i , lmmFit$compLikMats$listBlocks[[lmmFit$compLikMats$subsetPairsAdj[i,2]]]$i)
        if(attachBigMats){
          lmmFit$setupMats[[i]] <- setupIAK3D(lmmFit$xData[iThis,,drop = FALSE] , lmmFit$dIData[iThis,,drop = FALSE] , nDscPts = 0)
        }else{ # attach locations and depths as setupMats...
          lmmFit$setupMats[[i]] <- list('xData' = lmmFit$xData[iThis,,drop = FALSE] , 'dIData' = lmmFit$dIData[iThis,,drop = FALSE])
        }
      }
      ### now all individual subsets...
      for (i in 1:length(lmmFit$compLikMats$listBlocks)){
        iThis <- lmmFit$compLikMats$listBlocks[[i]]$i
        if(attachBigMats){
          lmmFit$setupMats[[nrow(lmmFit$compLikMats$subsetPairsAdj)+i]] <- setupIAK3D(lmmFit$xData[iThis,,drop = FALSE] , lmmFit$dIData[iThis,,drop = FALSE] , nDscPts = 0)
        }else{
          lmmFit$setupMats[[nrow(lmmFit$compLikMats$subsetPairsAdj)+i]] <- list('xData' = lmmFit$xData[iThis,,drop = FALSE] , 'dIData' = lmmFit$dIData[iThis,,drop = FALSE])
        }
      }
    }

    if(updateLmmPars){
      
      printnllTime <<- TRUE
      verboseOptim <<- TRUE  

### this will probably update beta - which I don't think I really want to do (though could do)
### probably just want to update list iCX things?
      if (lmmFit$compLikMats$compLikOptn == 0){
        lmmFitUpdate <- nllIAK3D(pars = lmmFit$pars , zData = lmmFit$zData , XData = lmmFit$XData , vXU = lmmFit$vXU , iU = lmmFit$iU , modelx = lmmFit$modelx , nud = lmmFit$nud ,
                    sdfdType_cd1 = lmmFit$sdfdType_cd1 , sdfdType_cxd0 = lmmFit$sdfdType_cxd0 , sdfdType_cxd1 = lmmFit$sdfdType_cxd1 ,
                    cmeOpt = lmmFit$cmeOpt , prodSum = lmmFit$prodSum , setupMats = lmmFit$setupMats , parBnds = lmmFit$parBnds , useReml = lmmFit$useReml , lnTfmdData = lmmFit$lnTfmdData , rtnAll = T , attachBigMats = attachBigMats)
      }else{
        lmmFitUpdate <- nllIAK3D_CL(pars = lmmFit$pars , zData = lmmFit$zData , XData = lmmFit$XData , modelx = lmmFit$modelx , nud = lmmFit$nud ,
                    sdfdType_cd1 = lmmFit$sdfdType_cd1 , sdfdType_cxd0 = lmmFit$sdfdType_cxd0 , sdfdType_cxd1 = lmmFit$sdfdType_cxd1 ,
                    cmeOpt = lmmFit$cmeOpt , prodSum = lmmFit$prodSum , setupMats = lmmFit$setupMats , parBnds = lmmFit$parBnds , useReml = lmmFit$useReml , compLikMats = lmmFit$compLikMats , rtnAll = T , attachBigMats = attachBigMats)
      }
      
### write over lmmFit elements...
      for(i in 1:length(lmmFitUpdate)){
        if(!is.null(lmmFitUpdate[[i]])){
          lmmFit[[names(lmmFitUpdate)[i]]] <- lmmFitUpdate[[i]]
        }else{
          lmmFit[names(lmmFitUpdate)[i]] <- list(NULL)
        }
      }
    }else{}
  }else{} # no need to update anything else as no local data given

  return(lmmFit)

}

###################################################################
### for local updating, need to redefine some matrices that are in the fitted object...
### if removeLocal, data coinciding with local data are removed from lmmFit
### else, are local data are added to lmmFit
###################################################################
lmmUpdateLocal.DEFUNCT <- function(lmmFit , xLocal = NULL , dILocal = NULL , zLocal = NULL , covsLocal = NULL , removeLocalData = FALSE , attachBigMats = TRUE , mindFromLegData = 0.001){
  
  #  mindFromLegData is for checking local data aren't already in legacy data. (0.001 = 1 m in km)
  
  ### to update stored variables from older versions of code
  lmmFit <- updateSplineNames(lmmFit)
  lmmFit <- updateCubistModel(lmmFit)
  
  if(!is.null(xLocal)){
    iKeep <- which((!is.na(zLocal)) & (rowSums(is.na(covsLocal)) == 0) & (apply(xyDist(lmmFit$xData , xLocal) , 2 , min) > mindFromLegData)) # last one is a check that not already in leg data   
    
    if(length(iKeep) == 0){ return(lmmFit) }else{}
    
    xLocal <- xLocal[iKeep,,drop=FALSE]
    dILocal <- dILocal[iKeep,,drop=FALSE]
    zLocal <- zLocal[iKeep]
    covsLocal <- covsLocal[iKeep,,drop=FALSE]
    
    if(removeLocalData){
      ### check if any local data are already in legacy data, and if so, remove from legacy...    
      iRm <- which(apply(xyDist(lmmFit$xData , xLocal) , 1 , min) < 0.001) 
      if(length(iRm) == 0){ return(lmmFit) }else{}
      lmmFit$xData <- lmmFit$xData[-iRm,,drop=FALSE]
      lmmFit$dIData <- lmmFit$dIData[-iRm,,drop=FALSE]
      lmmFit$zData <- lmmFit$zData[-iRm]
      lmmFit$covsData <- lmmFit$covsData[-iRm,,drop=FALSE]
      ### only do for new data...
      if (is.element('iU' , names(lmmFit))){
        iUUpdate <- lmmFit$iU
      }else{
        iUUpdate <- NA
      }
      
      XLims <- matrix(NA , 2 , ncol(lmmFit$XData))
      XLims[1,] <- Inf
      XLims[2,] <- -Inf
      
      XUpdate <- lmmFit$XData[-iRm,,drop=FALSE]
      if (is.element('vXU' , names(lmmFit))){
        iRmvXU <- rows2vXURows(iRm , iUUpdate)
        vXUUpdate <- lmmFit$vXU[-iRmvXU,,drop=FALSE]
      }else{
        vXUUpdate <- c()
      }
      
    }else{
      ### add the local data to the legacy data...    
      lmmFit$xData <- rbind(lmmFit$xData , xLocal)
      lmmFit$dIData <- rbind(lmmFit$dIData , dILocal)
      lmmFit$zData <- c(lmmFit$zData , zLocal)
      lmmFit$covsData <- rbind(lmmFit$covsData , covsLocal)
      ### only do for new data...
      if (is.element('iU' , names(lmmFit))){
        iUUpdate <- lmmFit$iU
      }else{
        iUUpdate <- NA
      }
      
      XLims <- matrix(NA , 2 , ncol(lmmFit$XData))
      XLims[1,] <- Inf
      XLims[2,] <- -Inf
      
      tmp <- makeXvX(covData = covsLocal , dIData = dILocal , modelX = lmmFit$modelX , allKnotsd = lmmFit$allKnotsd , iU = iUUpdate , nDiscPts = 1000 , lnTfmdData = lmmFit$lnTfmdData , XLims = XLims)
      XUpdate <- rbind(lmmFit$XData , tmp$X)
      if (is.element('vXU' , names(lmmFit))){
        vXUUpdate <- rbind(lmmFit$vXU , tmp$vXU)
      }else{
        vXUUpdate <- c()
      }
    }
    
    ###################################
    ### update the compLikMats...  
    ###################################
    if(lmmFit$compLikMats$compLikOptn > 0){
      vcentres <- matrix(NA , length(lmmFit$compLikMats$listBlocks) , ncol(lmmFit$compLikMats$listBlocks[[1]]$centroid))
      for (i in 1:length(lmmFit$compLikMats$listBlocks)){ vcentres[i,] <- lmmFit$compLikMats$listBlocks[[i]]$centroid }
      optnTmp <- lmmFit$compLikMats$compLikOptn
      lmmFit$compLikMats <- NULL # delete old vs.
      lmmFit$compLikMats <- setVoronoiBlocks(x = lmmFit$xData , nPerBlock = NA , plotVor = F , vcentres = vcentres)
      lmmFit$compLikMats$compLikOptn <- optnTmp
    }else{}
    
    ###################################
    ### update the setupMats...  
    ### could be speeded up in future to just update the matrices rather than remake them.
    ###################################
    lmmFit$setupMats <- NULL # delete any old vs.
    
    if (lmmFit$compLikMats$compLikOptn == 0){
      setupMats <- setupIAK3D(lmmFit$xData , lmmFit$dIData , nDscPts = 0)
    }else{
      setupMats <- list()
      ### order is all adj subset pairs, then all individual subsets (which can be used to get non-adj subset pairs)
      for (i in 1:nrow(lmmFit$compLikMats$subsetPairsAdj)){
        iThis <- c(lmmFit$compLikMats$listBlocks[[lmmFit$compLikMats$subsetPairsAdj[i,1]]]$i , lmmFit$compLikMats$listBlocks[[lmmFit$compLikMats$subsetPairsAdj[i,2]]]$i)
        setupMats[[i]] <- setupIAK3D(lmmFit$xData[iThis,,drop = FALSE] , lmmFit$dIData[iThis,,drop = FALSE] , nDscPts = 0)
      }
      ### now all individual subsets...      
      for (i in 1:length(lmmFit$compLikMats$listBlocks)){
        iThis <- lmmFit$compLikMats$listBlocks[[i]]$i
        setupMats[[nrow(lmmFit$compLikMats$subsetPairsAdj)+i]] <- setupIAK3D(lmmFit$xData[iThis,,drop = FALSE] , lmmFit$dIData[iThis,,drop = FALSE] , nDscPts = 0)
      }
    }
    
    printnllTime <<- TRUE
    verboseOptim <<- TRUE  
    
    ### this will probably update beta - which I don't think I really want to do (though could do)
    ### probably just want to update list iCX things?
    if (lmmFit$compLikMats$compLikOptn == 0){
      lmmFitUpdate <- nllIAK3D(pars = lmmFit$pars , zData = lmmFit$zData , XData = XUpdate , vXU = vXUUpdate , iU = iUUpdate , modelx = lmmFit$modelx , nud = lmmFit$nud ,
                               sdfdType_cd1 = lmmFit$sdfdType_cd1 , sdfdType_cxd0 = lmmFit$sdfdType_cxd0 , sdfdType_cxd1 = lmmFit$sdfdType_cxd1 ,
                               cmeOpt = lmmFit$cmeOpt , prodSum = lmmFit$prodSum , setupMats = setupMats , parBnds = lmmFit$parBnds , useReml = lmmFit$useReml , lnTfmdData = lmmFit$lnTfmdData , rtnAll = T , attachBigMats = attachBigMats)
    }else{
      lmmFitUpdate <- nllIAK3D_CL(pars = lmmFit$pars , zData = lmmFit$zData , XData = XUpdate , modelx = lmmFit$modelx , nud = lmmFit$nud ,
                                  sdfdType_cd1 = lmmFit$sdfdType_cd1 , sdfdType_cxd0 = lmmFit$sdfdType_cxd0 , sdfdType_cxd1 = lmmFit$sdfdType_cxd1 ,
                                  cmeOpt = lmmFit$cmeOpt , prodSum = lmmFit$prodSum , setupMats = setupMats , parBnds = lmmFit$parBnds , useReml = lmmFit$useReml , compLikMats = lmmFit$compLikMats , rtnAll = T , attachBigMats = attachBigMats)
    }
    
    if(attachBigMats){
      lmmFit$setupMats <- setupMats
    }else{}
    remove(setupMats)
    
    ### write over lmmFit elements...
    for(i in 1:length(lmmFitUpdate)){
      if(!is.null(lmmFitUpdate[[i]])){
        lmmFit[[names(lmmFitUpdate)[i]]] <- lmmFitUpdate[[i]]
      }else{
        lmmFit[names(lmmFitUpdate)[i]] <- list(NULL)
      }
    }
    
  }else{} # no need to update anything else as no local data given
  
  return(lmmFit)
  
}

################################################################
### function to change names for spline columns to new version...
################################################################
updateSplineNames <- function(lmmFit){
  iSpline <- which(lmmFit$namesX == 'dSpline')
  if(length(iSpline) > 0){
    lmmFit$namesX[iSpline] <- paste0('dSpline.' , seq(length(iSpline)))
  }else{}

  return(lmmFit)  

}

################################################################
### function to change stored variables in the cubist model to those used in new version...
################################################################
updateCubistModel <- function(lmmFit){

  if(is(lmmFit$modelX , 'cubist')){
    
### add a couple of things (if missing) that weren't kept in older version of code...    
    if(!is.element('allKnotsd' , names(lmmFit$modelX))){
      if(length(lmmFit$allKnotsd) > 0){
        lmmFit$modelX$allKnotsd <- lmmFit$allKnotsd
      }else{ # assigning as c() doesn't work, so: 
        lmmFit$modelX['allKnotsd'] <- list(NULL)
      }
    }else{}
    if(!is.element('namesDataFit' , names(lmmFit$modelX))){
      lmmFit$modelX$namesDataFit <- names(lmmFit$covsData)
    }else{}
    
    cmTmp <- cubist2X(cubistModel = lmmFit$modelX , dataFit = lmmFit$covsData , zFit = lmmFit$zData , allKnotsd = lmmFit$allKnotsd)
    lmmFit$modelX <- cmTmp$cubistModel

    if(!is.element('namesX' , names(lmmFit$modelX))){
      lmmFit$modelX$namesX <- lmmFit$namesX
    }else{}  
    
    if(!is.element('comms4XCubist' , names(lmmFit$modelX))){
      lmmFit$modelX$comms4XCubist <- rep(1 , length(lmmFit$modelX$names4XCubist))
    }else{}
  }else{}

  return(lmmFit)  

}
