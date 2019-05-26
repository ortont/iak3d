makeXvX <- function(covData = NA , dI , modelX , allKnotsd = c() , iU = NA , nDiscPts = 100 , lnTfmdData = FALSE , XLims = NULL){
###########################################################################
### name convention...
###########################################################################
### for constant, use 'const'
### for depth on its own, use 'd'
### for depth ^ 2, use 'd2'
### for interaction between d and covariate 'dem' (for eg), use 'd.dem' 
### for interaction between d ^ 2 and covariate 'dem' (for eg), use 'd2.dem' 
### for a spatial categorical covariate with nc classes, 
### include nc - 1 columns all with the same name
###########################################################################
### if modelX is integer/numeric = 0/1/2... 
### 	then set up X with all covariates and their i-order interactions with depth
### 	ie, if 0, then no interactions, if 1 then interactions with depth, if 1 then interactions with depth ^ 2 
### if modelX is a character vector...
###	then use the names conventions above to exctract appropriate design matrices
### if modelX is a cubist model...
### 	then use the cubist splits and rules to set design matrices
###########################################################################
    n <- dim(dI)[[1]]
    if((length(covData) == 1) && is.na(covData)){
        nCovs <- 0
    }else{
        nCovs <- dim(covData)[[2]]
    }

    if(length(allKnotsd) > 0){ incdSpline <- TRUE }else{ incdSpline <- FALSE }
    
    modelXIn <- modelX
    
    if(is.integer(modelX) | is.numeric(modelX)){
        modelType <- 'mlr'
        namesX <- c('const')
#############################################################
### if modelX = 0, just covData
### if modelX = 0.5, covData + d
### if modelX = 1, covData + d + d.covData
### if modelX = 1.5, covData + d + d.covData + d2
### if modelX = 2, covData + d + d.covData + d2 + d2.covData
### if modelX = 2.5, covData + d + d.covData + d2 + d2.covData + d3
### if modelX = 3, covData + d + d.covData + d2 + d2.covData + d3 + d3.covData
#############################################################
	  for (o in 0:floor(modelX)){
        if (o == 0){ 
		  dNameInts <- '' 
  		  dNameUni <- 'd' 
	    }else if(o == 1){
		  dNameInts <- 'd.' 
		  dNameUni <- 'd2' 
	    }else if(o == 2){
		  dNameInts <- 'd2.' 
		  dNameUni <- 'd3' 
	    }else if(o == 3){
		  dNameInts <- 'd3.' 
		  dNameUni <- NA #ie 3.5 not programmed 
	    }else{
          stop('Error - the makeXvX function has only been coded for up to cubic functions of depth - generalise the code to get higher-order polynomials, or try a spline / Cubist model!')
        }

        if(nCovs > 0){
          for (j in 1:nCovs){
            if(is.factor(covData[[j]])){
              for(l in 1:(length(levels(covData[[j]])) - 1)){
                namesX <- c(namesX , paste0(dNameInts , names(covData)[j]))
              }
            }else{
                namesX <- c(namesX , paste0(dNameInts , names(covData)[j]))
            }
          }
        }else{}

        if(modelX > o){
		  namesX <- c(namesX , dNameUni)
        }else{}
      }

	  p <- length(namesX)

    }else if(is.character(modelX)){
      modelType <- 'mlr'
  	  namesX <- modelX
	  p <- length(namesX)

    }else if(class(modelX) == "cubist"){
########################################
### assume cubist model was fitted with covariates covData
### and that 'dIMidPts' is included in covData
### note inconsistency with approach for linear models, 
### where covData was just spatial covariates
### done it this way, so that covData is what was presented to cubist
### check depth is in covData...
### and check the name of the depth interval mid points is dIMidPts
########################################
      modelType <- 'cubist'

      idInCovData <- which(names(covData) == 'dIMidPts')
      if(length(idInCovData) == 0){
        stop('Refit the cubist model, but make the name of the depth variable dIMidPts')
      }else{}
 
### just to get the names and number of parameters...
      namesX <- modelX$names4XCubist
      p <- sum(modelX$pVec)

    }else{}

###############################################
### if allKnotsd is not empty, then it gives knots (internal and boundary) 
### for a cubic spline for the general depth function 
### (ie an approximation to the depth-correlated random effect)
### in this case, remove all polynomial terms (d, d2,...), from fixed effects
### but keep constant as bs specified with intercept = F.
###############################################
    if(length(allKnotsd) > 0){

#        if(length(allKnotsd) < 3){ stop('Error - must input two boundary knots and at least one internal knot in allKnotsd!') }else{}
        if(length(allKnotsd) < 2){ stop('Error - must input at least two boundary knots and at least one internal knot in allKnotsd!') }else{}
        
        intKnots <- allKnotsd[-1]
        intKnots <- intKnots[-length(intKnots)]
        bdryKnots <- c(allKnotsd[1] , allKnotsd[length(allKnotsd)])

        if(is.integer(modelXIn) | is.numeric(modelXIn) | class(modelXIn) == "cubist"){
          XSplineCheck <- bs(x = c(0.1 , 0.2) , knots = intKnots , degree = 3 , intercept = F , Boundary.knots = bdryKnots)
          namesX <- c(namesX , rep('dSpline' , ncol(XSplineCheck))) 
          p <- length(namesX)
        }else{} # if came in as character vector, this has already been done.
        
        if (modelType == 'mlr'){
          iRemove <- which((namesX == 'd') | (namesX == 'd2') | (namesX == 'd3'))    
          if(length(iRemove) > 0){
            namesX <- namesX[-iRemove]
            p <- p - length(iRemove)
          }else{}        
        }else{

        }
    }else{}

###############################################
### initialize... 
###############################################
    X <- matrix(NA , n , p)
    if(lnTfmdData){ 
      vXU <- matrix(NA , n * p , p)
    }else{
      vXU <- NA
    }    

    if(length(allKnotsd) > 0){
      ipSpline <- which(namesX == 'dSpline')
    }else{}
    
    if(is.null(XLims)){
      setXLims <- TRUE
      XLims <- matrix(NA , 2 , p)
      XLims[1,] <- Inf
      XLims[2,] <- -Inf
### don't apply to the spline functions though...   
      if(length(allKnotsd) > 0){
        XLims[1,ipSpline] <- -Inf
        XLims[2,ipSpline] <- Inf
      }else{}
    }else{
      setXLims <- FALSE
    }

# number of parameters associated with each name
# 1 for continuous variable, nClasses - 1 for categorical, put in first column with that name (rest are NA)
    pX <- NA * integer(p) 

###############################################
### if we have a multiple linear regression model...
###############################################
    if (modelType == 'mlr'){
### first put in with d = 1
        iconst <- which(namesX  == 'const')    
        id <- which(namesX == 'd')    
        id2 <- which(namesX == 'd2')    
        id3 <- which(namesX == 'd3')    
        idInt <- which(substr(namesX , 1 , 2) == 'd.')    
        id2Int <- which(substr(namesX , 1 , 3) == 'd2.')    
        id3Int <- which(substr(namesX , 1 , 3) == 'd3.')    
        idSpline <- which(namesX == 'dSpline')    

        namesX_d <- namesX
        namesX_d[c(id,id2,id3,idSpline)] <- 'const'
        namesX_d[idInt] <- substr(namesX_d[idInt] , 3 , nchar(namesX_d[idInt]))
        namesX_d[id2Int] <- substr(namesX_d[id2Int] , 4 , nchar(namesX_d[id2Int]))
        namesX_d[id3Int] <- substr(namesX_d[id3Int] , 4 , nchar(namesX_d[id3Int]))

        X_d <- matrix(NA , n , p)
        for (j in 1:p){
          if (namesX_d[j] == 'const'){
            X_d[,j] <- 1
            pX[j] <- 1
          }else{
            covDataThis <- covData[[namesX_d[j]]]
            if(is.factor(covDataThis)){
### check name against previous name...
### treatment coding: first additional column for factor variable is 2nd level...
              if(namesX[j] != namesX[j-1]){ 
                l <- 2 
                pX[j] <- length(levels(covDataThis)) - 1
              }else{}
              X_d[,j] <- 0
              X_d[which(covDataThis == levels(covDataThis)[l]),j] <- 1
              l <- l + 1
            }else{
              X_d[,j] <- covData[[namesX_d[j]]]
              pX[j] <- 1
            }
          }
        }

        for(i in 1:n){
          XTmp <- kronecker(matrix(X_d[i,] , nrow = 1) , matrix(1 , nDiscPts , 1))
          if(nDiscPts > 1){
              dDiscPts <- seq(dI[i,1] , dI[i,2] , (dI[i,2] - dI[i,1]) / (nDiscPts - 1))
          }else{
              dDiscPts <- 0.5 * (dI[i,1] + dI[i,2])
          }
          if(length(id) > 0){ XTmp[,id] <- XTmp[,id] * dDiscPts }else{}
          if(length(id2) > 0){ XTmp[,id2] <- XTmp[,id2] * (dDiscPts ^ 2) }else{}
          if(length(id3) > 0){ XTmp[,id3] <- XTmp[,id3] * (dDiscPts ^ 3) }else{}
          if(length(idInt) > 0){ XTmp[,idInt] <- XTmp[,idInt] * kronecker(dDiscPts , matrix(1 , 1 , length(idInt))) }else{}
          if(length(id2Int) > 0){ XTmp[,id2Int] <- XTmp[,id2Int] * kronecker(dDiscPts ^ 2 , matrix(1 , 1 , length(id2Int))) }else{}
          if(length(id3Int) > 0){ XTmp[,id3Int] <- XTmp[,id3Int] * kronecker(dDiscPts ^ 3 , matrix(1 , 1 , length(id3Int))) }else{}
          if(length(idSpline) > 0){ 
            XTmp[,idSpline] <- XTmp[,idSpline] * bs(x = dDiscPts , knots = intKnots , degree = 3 , intercept = F , Boundary.knots = bdryKnots) 
          }else{}

          if(setXLims){
            XLims[1,] <- apply(rbind(XTmp , XLims[1,]) , 2 , minNonZero)
            XLims[2,] <- apply(rbind(XTmp , XLims[2,]) , 2 , maxNonZero)
          }else{
### apply the given constraints to the point-support design matrix...
            XTmp <- matrix(mapply(replaceLT , x = t(XTmp) , llim = XLims[1,] , replaceZeros = FALSE , SIMPLIFY = TRUE) , nrow = nrow(XTmp) , ncol = ncol(XTmp) , byrow = TRUE)
            XTmp <- matrix(mapply(replaceGT , x = t(XTmp) , ulim = XLims[2,] , replaceZeros = FALSE , SIMPLIFY = TRUE) , nrow = nrow(XTmp) , ncol = ncol(XTmp) , byrow = TRUE)
          }

          X[i,] <- colMeans(XTmp)

          if(lnTfmdData){ 
            iThis <- seq((i - 1) * p + 1 , i * p)
            vXU[iThis,] <- var(XTmp)
          }else{}
          
        }

###############################################
### or if we have a cubist model...
###############################################
    }else if(modelType == 'cubist'){

#n x p x nDiscPts <= maxTmp

        maxTmp <- 10000000

        if(lnTfmdData){ 
          nPerBatch <- floor(maxTmp / (nDiscPts * p * p))
        }else{
          nPerBatch <- floor(maxTmp / (nDiscPts * p))
        }

        if(nPerBatch < 1){ 
          nPerBatch <- 1 
        }else if(nPerBatch > n){
          nPerBatch <- n
        }else{}
        nBatches <- ceiling(n / nPerBatch)

        fnVarApply <- function(idx , X4Apply){ var(X4Apply[idx,,drop=FALSE]) }
        fnMeanApply <- function(idx , X4Apply){ colMeans(X4Apply[idx,,drop=FALSE]) }

        for (iBatch in 1:nBatches){

 #         print(paste0('Making X for Batch ' , iBatch , ' of ' , nBatches))

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
            dDiscPts <- dI[iRep,1] + rep(linTmp , length(i)) * (dI[iRep,2] - dI[iRep,1])
          }else{
            dDiscPts <- 0.5 * (dI[i,1] + dI[i,2])
          }          

          covDataThis <- covData[iRep,,drop=FALSE]
          covDataThis$dIMidPts <- dDiscPts
          
          XTmp <- cubist2X(cubistModel = modelX , dataFit = covDataThis , incdSpline = incdSpline)
          XTmp <- XTmp$X

          if(length(allKnotsd) > 0){
              XTmp <- cbind(XTmp , bs(x = dDiscPts , knots = intKnots , degree = 3 , intercept = F , Boundary.knots = c(allKnotsd[1] , allKnotsd[length(allKnotsd)])) )
          }else{}

          if(setXLims){
            XLims[1,] <- apply(rbind(XTmp , XLims[1,]) , 2 , minNonZero)
            XLims[2,] <- apply(rbind(XTmp , XLims[2,]) , 2 , maxNonZero)
          }else{
### apply the given constraints to the point-support design matrix...
            XTmp <- matrix(mapply(replaceLT , x = t(XTmp) , llim = XLims[1,] , replaceZeros = FALSE , SIMPLIFY = TRUE) , nrow = nrow(XTmp) , ncol = ncol(XTmp) , byrow = TRUE)
            XTmp <- matrix(mapply(replaceGT , x = t(XTmp) , ulim = XLims[2,] , replaceZeros = FALSE , SIMPLIFY = TRUE) , nrow = nrow(XTmp) , ncol = ncol(XTmp) , byrow = TRUE)
          }

          mXTmp <- t(apply(i4Apply , 2 , fnMeanApply , X4Apply = XTmp))
          X[i,] <- as.matrix(mXTmp)

          if(lnTfmdData){ 
### remove means from XTmp...
            XTmp <- XTmp - X[iRep,,drop=FALSE] 
            
            iThis <- rep((i-1)*p , each = p) + rep(seq(1,p) , length(i))

            vXTmp <- apply(i4Apply , 2 , fnVarApply , X4Apply = XTmp)
            vXTmp <- matrix(vXTmp , length(i) * p , p , byrow = TRUE)
            vXU[iThis,] <- vXTmp

          }else{}
        }
        
### all variables within Cubist rules are continuous, so...
        pX <- 1 + integer(p)


    }else{}


###############################################
### finally, remove any variances close to zero...
###############################################
    if(lnTfmdData){ 
      if(is.na(iU[1])){
### use names to determine which will vary within sample support
        if(modelType == 'mlr'){
          iU <- which(namesX == 'd' | namesX == 'd2' | substr(namesX , 1 , 2) == 'd.' | substr(namesX , 1 , 3) == 'd2.' | substr(namesX , 1 , 3) == 'd3.' | namesX == 'dSpline') 
          iK <- setdiff(seq(p) , iU)
          
        }else if(modelType == 'cubist'){
### iU for cubist model is all columns. add spline cols if necessary...
### even if d not part of rule conditions or predictor, 
### d can affect the prob of being in a rule (eg if 3 rules are apt for d = 0.1 but only 2 for d = 0.2)
### so safest to assume all cols affects by d
          iU <- seq(p)
          iK <- c()
        
        }else{
          stop('Invalid model type!')
        }


      }else{
        iK <- setdiff(seq(p) , iU)
      }

      if(length(iK) > 0){
        vXU <- matrix(vXU[,-iK] , nrow = n * p)
        iRemove <- kronecker(seq(0 , n - 1) * p , matrix(1 , length(iK) , 1)) + kronecker(matrix(1 , n , 1) , iK)        
        vXU <- vXU[-iRemove,]
      }else{}    

    }else{
      iU <- NA
    }
      
### just to make sure things are matrices...    
    X <- matrix(X , nrow = n , ncol = p)
    if(lnTfmdData){ 
      vXU <- matrix(vXU , nrow = n * length(iU) , ncol = length(iU))
    }else{}

    return(list('X' = X , 'vXU' = vXU , 'iU' = iU , 'namesX' = namesX , 'pX' = pX , 'XLims' = XLims))     
}

#########################################################
### scale covariates to std normals and remember scaling parameters...
#########################################################
scaleCovs <- function(covsData , scalePars = matrix(NA)){
    if(is.na(scalePars[1,1])){
      calcScaleParsNow <- T
      scalePars <- matrix(NA , 2 , dim(covsData)[[2]])
      scalePars <- as.data.frame(scalePars)
      names(scalePars) <- names(covsData)
    }else{
      calcScaleParsNow <- F
    }
    
    for (j in 1:dim(covsData)[[2]]){
      if(class(covsData[[j]]) != 'factor'){
        if(calcScaleParsNow){
          jThis <- j
          scalePars[1,jThis] <- mean(covsData[[j]])
          scalePars[2,jThis] <- sqrt(var((covsData[[j]])))
        }else{
          jThis <- which(names(scalePars) == names(covsData)[j])
          if(length(jThis) != 1){         
            print(paste0('Error with scaling for ' , names(covsMap)[j]))
            stop('Check the covariate names in fitting, scaling and mapping data frames!')
          }else{}
        }
        covsData[[j]] <- (covsData[[j]] - scalePars[1,jThis]) / scalePars[2,jThis]
      }else{}
    }
    
    return(list('covsData' = covsData , 'scalePars' = scalePars))
}


getIndexLower <- function(n , incDiag = FALSE){
  z <- sequence(n)
  if(incDiag){
    idxL <- cbind(row = unlist(lapply(1:n, function(x) x:n), use.names = FALSE),
                  col = rep(z, times = rev(z)))
  }else{
    idxL <- cbind(row = unlist(lapply(2:n, function(x) x:n), use.names = FALSE),
                  col = rep(z[-length(z)], times = rev(tail(z, -1))-1))
  }
  return(idxL)
}

getIndexUpper <- function(n , incDiag = FALSE){
  z <- sequence(n)
  if(incDiag){
    idxL <- cbind(row = rep(z, times = rev(z)),
                  col = unlist(lapply(1:n, function(x) x:n), use.names = FALSE))
  }else{
    idxL <- cbind(row = rep(z[-length(z)], times = rev(tail(z, -1))-1),
                  col = unlist(lapply(2:n, function(x) x:n), use.names = FALSE))
  }
  return(idxL)
}

replaceLT <- function(x , llim , replaceZeros = TRUE){ 
  if(replaceZeros){
    x[which(x < llim)] <- llim 
  }else{
    x[which(x < llim & x != 0)] <- llim 
  }
  return(x)
}
replaceGT <- function(x , ulim , replaceZeros = TRUE){ 
  if(replaceZeros){
    x[which(x > ulim)] <- ulim 
  }else{
    x[which(x > ulim & x != 0)] <- ulim 
  }
  return(x)
}


minNonZero <- function(x , all0Val = Inf){
  xNZ <- x[x!=0]
  if(length(xNZ) == 0){
    return(all0Val)
  }else{
    return(min(xNZ))  
  }
}

maxNonZero <- function(x , all0Val = -Inf){
  xNZ <- x[x!=0]
  if(length(xNZ) == 0){
    return(all0Val)
  }else{
    return(max(xNZ))  
  }
}

