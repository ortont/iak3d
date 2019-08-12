### for local updating, need to redefine some matrices that are in the fitted object...
lmmUpdateLocal <- function(lmmFit , xLocal = NULL , dILocal = NULL , zLocal = NULL , covsLocal = NULL){

### to update stored variables from older versions of code
  lmmFit <- updateSplineNames(lmmFit)
  lmmFit <- updateCubistModel(lmmFit)

  if(!is.null(xLocal)){
    iKeep <- which((!is.na(zLocal)) & (rowSums(is.na(covsLocal)) == 0))
    xLocal <- xLocal[iKeep,,drop=FALSE]
    dILocal <- dILocal[iKeep,,drop=FALSE]
    zLocal <- zLocal[iKeep]
    covsLocal <- covsLocal[iKeep,,drop=FALSE]
  
    lmmFit$x <- rbind(lmmFit$x , xLocal)
    lmmFit$dI <- rbind(lmmFit$dI , dILocal)
    lmmFit$z <- c(lmmFit$z , zLocal)
    lmmFit$covs <- rbind(lmmFit$covs , covsLocal)
  
### only do for new data...
    if (is.element('iU' , names(lmmFit))){
      iUUpdate <- lmmFit$iU
    }else{
      iUUpdate <- NA
    }
        
    XLims <- matrix(NA , 2 , ncol(lmmFit$X))
    XLims[1,] <- Inf
    XLims[2,] <- -Inf
 
    tmp <- makeXvX(covData = covsLocal , dI = dILocal , modelX = lmmFit$modelX , allKnotsd = lmmFit$allKnotsd , iU = iUUpdate , nDiscPts = 1000 , lnTfmdData = lmmFit$lnTfmdData , XLims = XLims)
    XUpdate <- rbind(lmmFit$X , tmp$X)
    if (is.element('vXU' , names(lmmFit))){
      vXUUpdate <- rbind(lmmFit$vXU , tmp$vXU)
    }else{
      vXUUpdate <- c()
    }
  
###################################
### update the compLikMats...  
###################################
    if(lmmFit$compLikMats$compLikOptn > 0){
      vcentres <- matrix(NA , length(lmmFit$compLikMats$listBlocks) , ncol(lmmFit$compLikMats$listBlocks[[1]]$centroid))
      for (i in 1:length(lmmFit$compLikMats$listBlocks)){ vcentres[i,] <- lmmFit$compLikMats$listBlocks[[i]]$centroid }
      optnTmp <- lmmFit$compLikMats$compLikOptn
      lmmFit$compLikMats <- setVoronoiBlocks(x = lmmFit$x , nPerBlock = NA , plotVor = F , vcentres = vcentres)
      lmmFit$compLikMats$compLikOptn <- optnTmp
    }else{}

###################################
### update the setupMats...  
### could be speeded up in future to just update the matrices rather than remake them.
###################################
    if (lmmFit$compLikMats$compLikOptn == 0){
        lmmFit$setupMats <- setupIAK3D(lmmFit$x , lmmFit$dI , nDscPts = 0)
    }else{
        lmmFit$setupMats <- list()
### order is all adj subset pairs, then all individual subsets (which can be used to get non-adj subset pairs)
        for (i in 1:nrow(lmmFit$compLikMats$subsetPairsAdj)){
          iThis <- c(lmmFit$compLikMats$listBlocks[[lmmFit$compLikMats$subsetPairsAdj[i,1]]]$i , lmmFit$compLikMats$listBlocks[[lmmFit$compLikMats$subsetPairsAdj[i,2]]]$i)
          lmmFit$setupMats[[i]] <- setupIAK3D(lmmFit$x[iThis,,drop = FALSE] , lmmFit$dI[iThis,,drop = FALSE] , nDscPts = 0)
        }
### now all individual subsets...      
        for (i in 1:length(lmmFit$compLikMats$listBlocks)){
          iThis <- lmmFit$compLikMats$listBlocks[[i]]$i
          lmmFit$setupMats[[nrow(lmmFit$compLikMats$subsetPairsAdj)+i]] <- setupIAK3D(lmmFit$x[iThis,,drop = FALSE] , lmmFit$dI[iThis,,drop = FALSE] , nDscPts = 0)
        }
    }
  
    printnllTime <<- TRUE
    verboseOptim <<- TRUE  

### this will probably update beta - which I don't think I really want to do (though could do)
### probably just want to update list iCX things?
    if (lmmFit$compLikMats$compLikOptn == 0){
      lmmFitUpdate <- nllIAK3D(pars = lmmFit$pars , z = lmmFit$z , X = XUpdate , vXU = vXUUpdate , iU = iUUpdate , modelx = lmmFit$modelx , nud = lmmFit$nud , 
                  sdfdType_cd1 = lmmFit$sdfdType_cd1 , sdfdType_cxd0 = lmmFit$sdfdType_cxd0 , sdfdType_cxd1 = lmmFit$sdfdType_cxd1 , 
                  cmeOpt = lmmFit$cmeOpt , prodSum = lmmFit$prodSum , setupMats = lmmFit$setupMats , parBnds = lmmFit$parBnds , useReml = lmmFit$useReml , lnTfmdData = lmmFit$lnTfmdData , rtnAll = T)	
    }else{
      lmmFitUpdate <- nllIAK3D_CL(pars = lmmFit$pars , z = lmmFit$z , X = XUpdate , modelx = lmmFit$modelx , nud = lmmFit$nud , 
                  sdfdType_cd1 = lmmFit$sdfdType_cd1 , sdfdType_cxd0 = lmmFit$sdfdType_cxd0 , sdfdType_cxd1 = lmmFit$sdfdType_cxd1 , 
                  cmeOpt = lmmFit$cmeOpt , prodSum = lmmFit$prodSum , setupMats = lmmFit$setupMats , parBnds = lmmFit$parBnds , useReml = lmmFit$useReml , compLikMats = lmmFit$compLikMats , rtnAll = T)	
    }

### write over lmmFit elements...
    for(i in 1:length(lmmFitUpdate)){
      lmmFit[[names(lmmFitUpdate)[i]]] <- lmmFitUpdate[[i]]
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

  if(class(lmmFit$modelX) == 'cubist'){
    
### add a couple of things (if missing) that weren't kept in older version of code...    
    if(!is.element('allKnotsd' , names(lmmFit$modelX))){
      if(length(lmmFit$allKnotsd) > 0){
        lmmFit$modelX$allKnotsd <- lmmFit$allKnotsd
      }else{ # assigning as c() doesn't work, so: 
        lmmFit$modelX['allKnotsd'] <- list(NULL)
      }
    }else{}
    if(!is.element('namesDataFit' , names(lmmFit$modelX))){
      lmmFit$modelX$namesDataFit <- names(lmmFit$covs)
    }else{}
    
    cmTmp <- cubist2X(cubistModel = lmmFit$modelX , dataFit = lmmFit$covs , zFit = lmmFit$z , allKnotsd = lmmFit$allKnotsd)
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
