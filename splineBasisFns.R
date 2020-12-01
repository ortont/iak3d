####################################################################################
####################################################################################
### AS SET UP HERE, TO FIT, ALL BDRYKNOTS AND INTKNOTS MUST BE WITHIN RANGE OF DATA 
### (IE AT LEAST ONE DATA PT BEYOND BDRYKNOTS)
####################################################################################
####################################################################################


####################################################################################
### fn to get positions for knots based on equal quantiles of data within bdryKnots
####################################################################################
getQuantileKnots <- function(x , bdryKnots , nIntKnots){
  if (nIntKnots == 0){
    return(c())
  }else{
    xIn <- x[x >= bdryKnots[1] & x <= bdryKnots[2]]
    qStep <- 1 / (nIntKnots + 1)
    intKnots <- quantile(xIn , seq(qStep , 1 - qStep , qStep))
    return(intKnots)
  }
}

####################################################################################
### fn to make basis fn for cubic spline
####################################################################################
makeXbs <- function(x , knots){
  if(is.null(dim(x)) || (ncol(x) == 1)){
    ### x just vector of known vals of x    
    n <- length(x)
    iaX <- FALSE
  }else{
    ### x should have 2 cols, first is lb, 2nd is ub for x
    n <- nrow(x)
    iaX <- TRUE
  }
  
  X <- matrix(1 , n , 4 + length(knots))
  
  if(iaX){
    X[,2] <- 0.5 * (x[,1] + x[,2])
    X[,3] <- (x[,2] ^ 3 - x[,1] ^ 3) / (3 * (x[,2] - x[,1]))
    X[,4] <- (x[,2] ^ 4 - x[,1] ^ 4) / (4 * (x[,2] - x[,1]))
  }else{
    X[,2] <- x
    X[,3] <- x ^ 2
    X[,4] <- x ^ 3
  }
  
  if(length(knots) > 0){
    if(iaX){
      ### inc-avd vs... 
      for (ik in 1:length(knots)){
        tmpL <- (x[,1] - knots[ik])
        tmpL[tmpL < 0] <- 0
        tmpU <- (x[,2] - knots[ik])
        tmpU[tmpU < 0] <- 0
        X[,4+ik] <- (tmpU ^ 4 - tmpL ^ 4) / (4 * (x[,2] - x[,1]))
      }
    }else{
      ### pt-supp vs... 
      for (ik in 1:length(knots)){
        tmp <- (x - knots[ik])
        tmp[tmp < 0] <- 0
        X[,4+ik] <- tmp ^ 3
      }
    }
  }else{}
  
  return(X)  
}

####################################################################################
### fn to make basis fn for natural cubic spline
### if intKnots is empty, will give 2 cols, 1 and x (ie a linear fn would be fitted)
####################################################################################
makeXns <- function(x , bdryKnots , intKnots){
  if(is.null(dim(x)) || (ncol(x) == 1)){
    ### x just vector of known vals of x    
    n <- length(x)
    iaX <- FALSE
  }else{
    ### x should have 2 cols, first is lb, 2nd is ub for x
    n <- nrow(x)
    iaX <- TRUE
  }
  
  X <- matrix(1 , n , 2 + length(intKnots))
  if(iaX){
    X[,2] <- 0.5 * (x[,1] + x[,2])
  }else{
    X[,2] <- x
  }
  
  if(length(intKnots) > 0){
    if(iaX){
      ### inc-avd vs... 
      tmpL <- (x[,1] - bdryKnots[1])
      tmpL[tmpL < 0] <- 0
      tmpU <- (x[,2] - bdryKnots[1])
      tmpU[tmpU < 0] <- 0
      Ex_xL3 <- 0.25 * (tmpU ^ 4 - tmpL ^ 4) / (x[,2] - x[,1])

      tmpL <- (x[,1] - bdryKnots[2])
      tmpL[tmpL < 0] <- 0
      tmpU <- (x[,2] - bdryKnots[2])
      tmpU[tmpU < 0] <- 0
      Ex_xU3 <- 0.25 * (tmpU ^ 4 - tmpL ^ 4) / (x[,2] - x[,1])

      for (ik in 1:length(intKnots)){
        tmpL <- (x[,1] - intKnots[ik])
        tmpL[tmpL < 0] <- 0
        tmpU <- (x[,2] - intKnots[ik])
        tmpU[tmpU < 0] <- 0
        
        X[,2+ik] <- 0.25 * (tmpU ^ 4 - tmpL ^ 4) / (x[,2] - x[,1]) - ((bdryKnots[2] - intKnots[ik]) / (bdryKnots[2] - bdryKnots[1])) * Ex_xL3 + ((bdryKnots[1] - intKnots[ik]) / (bdryKnots[2] - bdryKnots[1])) * Ex_xU3
      }
      
    }else{
      ### pt-supp vs...      
      tmp <- (x - bdryKnots[1])
      tmp[tmp < 0] <- 0
      x_xL3 <- (tmp ^ 3)
      
      tmp <- (x - bdryKnots[2])
      tmp[tmp < 0] <- 0
      x_xU3 <- (tmp ^ 3) 
      for (ik in 1:length(intKnots)){
        tmp <- (x - intKnots[ik])
        tmp[tmp < 0] <- 0
        X[,2+ik] <- tmp ^ 3 - ((bdryKnots[2] - intKnots[ik]) / (bdryKnots[2] - bdryKnots[1])) * x_xL3 + ((bdryKnots[1] - intKnots[ik]) / (bdryKnots[2] - bdryKnots[1])) * x_xU3
      }
    }
  }else{}
  
  return(X)  
}

####################################################################################
### fn to make basis fn for natural cubic spline with gradient clamped to 0 beyond upper bdryKnot
### if intKnots is empty, will give 1 col of 1s (ie a constant would be fitted)
####################################################################################
makeXnsClampedUG0 <- function(x , bdryKnots , intKnots , paramVs = 0){
  
  if(is.null(dim(x)) || (ncol(x) == 1)){
### x just vector of known vals of x    
    n <- length(x)
    iaX <- FALSE
  }else{
### x should have 2 cols, first is lb, 2nd is ub for x
    n <- nrow(x)
    iaX <- TRUE
  }
  
  if(paramVs == 0){ # a, g1, ... , gK
    X <- matrix(1 , n , 1 + length(intKnots))
    
    if(length(intKnots) > 0){
      if(iaX){
        ### inc-avd vs...      
        tmpL <- (x[,1] - bdryKnots[1])
        tmpL[tmpL < 0] <- 0
        tmpU <- (x[,2] - bdryKnots[1])
        tmpU[tmpU < 0] <- 0
        Ex_xL3 <- 0.25 * (tmpU ^ 4 - tmpL ^ 4) / (x[,2] - x[,1])
        
        tmpL <- (x[,1] - bdryKnots[2])
        tmpL[tmpL < 0] <- 0
        tmpU <- (x[,2] - bdryKnots[2])
        tmpU[tmpU < 0] <- 0
        Ex_xU3 <- 0.25 * (tmpU ^ 4 - tmpL ^ 4) / (x[,2] - x[,1])
        
        Ex <- 0.5 * (x[,1] + x[,2])
        
        for (ik in 1:length(intKnots)){
          tmpL <- (x[,1] - intKnots[ik])
          tmpL[tmpL < 0] <- 0
          tmpU <- (x[,2] - intKnots[ik])
          tmpU[tmpU < 0] <- 0
          
          X[,1+ik] <- 0.25 * (tmpU ^ 4 - tmpL ^ 4) / (x[,2] - x[,1]) - ((bdryKnots[2] - intKnots[ik]) / (bdryKnots[2] - bdryKnots[1])) * Ex_xL3 + 
            ((bdryKnots[1] - intKnots[ik]) / (bdryKnots[2] - bdryKnots[1])) * Ex_xU3 + 
            3 * (bdryKnots[2] - intKnots[ik]) * (intKnots[ik] - bdryKnots[1]) * Ex
          
        }
      }else{
        ### pt-supp vs...      
        tmp <- (x - bdryKnots[1])
        tmp[tmp < 0] <- 0
        x_xL3 <- (tmp ^ 3) 
        
        tmp <- (x - bdryKnots[2])
        tmp[tmp < 0] <- 0
        x_xU3 <- (tmp ^ 3) 
        
        for (ik in 1:length(intKnots)){
          tmp <- (x - intKnots[ik])
          tmp[tmp < 0] <- 0
          X[,1+ik] <- tmp ^ 3 - ((bdryKnots[2] - intKnots[ik]) / (bdryKnots[2] - bdryKnots[1])) * x_xL3 + 
            ((bdryKnots[1] - intKnots[ik]) / (bdryKnots[2] - bdryKnots[1])) * x_xU3 + 
            3 * (bdryKnots[2] - intKnots[ik]) * (intKnots[ik] - bdryKnots[1]) * x
        }
      }
    }else{}
    
  }else if(paramVs == 1){ # a, g1, ... , gK-1, w; where w is value at x = xU
    nIntKnots <- length(intKnots)
    
    X <- matrix(1 , n , 1 + nIntKnots)

    if(nIntKnots > 0){
      if(iaX){
        ### inc-avd vs...      
        tmpL <- (x[,1] - bdryKnots[1])
        tmpL[tmpL < 0] <- 0
        tmpU <- (x[,2] - bdryKnots[1])
        tmpU[tmpU < 0] <- 0
        Ex_xL3 <- 0.25 * (tmpU ^ 4 - tmpL ^ 4) / (x[,2] - x[,1])
        
        tmpL <- (x[,1] - bdryKnots[2])
        tmpL[tmpL < 0] <- 0
        tmpU <- (x[,2] - bdryKnots[2])
        tmpU[tmpU < 0] <- 0
        Ex_xU3 <- 0.25 * (tmpU ^ 4 - tmpL ^ 4) / (x[,2] - x[,1])
        
        tmpL <- (x[,1] - intKnots[nIntKnots])
        tmpL[tmpL < 0] <- 0
        tmpU <- (x[,2] - intKnots[nIntKnots])
        tmpU[tmpU < 0] <- 0
        Ex_xK3 <- 0.25 * (tmpU ^ 4 - tmpL ^ 4) / (x[,2] - x[,1])

        Ex <- 0.5 * (x[,1] + x[,2])
        
        Ehxk <- Ex_xK3 - ((bdryKnots[2] - intKnots[nIntKnots]) / (bdryKnots[2] - bdryKnots[1])) * Ex_xL3 + 
          ((bdryKnots[1] - intKnots[nIntKnots]) / (bdryKnots[2] - bdryKnots[1])) * Ex_xU3 + 
          3 * (bdryKnots[2] - intKnots[nIntKnots]) * (intKnots[nIntKnots] - bdryKnots[1]) * Ex

        denomTmp <- ((bdryKnots[2] - intKnots[nIntKnots]) * (intKnots[nIntKnots] - bdryKnots[1]) * (intKnots[nIntKnots] + bdryKnots[1] + bdryKnots[2]))

        X[,1] <- 1 - Ehxk / denomTmp
        
        if(nIntKnots > 1){
          for (ik in 1:(nIntKnots-1)){
            tmpL <- (x[,1] - intKnots[ik])
            tmpL[tmpL < 0] <- 0
            tmpU <- (x[,2] - intKnots[ik])
            tmpU[tmpU < 0] <- 0
            
            qtmp <- ((bdryKnots[2] - intKnots[ik]) * (intKnots[ik] - bdryKnots[1]) * (intKnots[ik] + bdryKnots[1] + bdryKnots[2])) / denomTmp

            X[,1+ik] <- 0.25 * (tmpU ^ 4 - tmpL ^ 4) / (x[,2] - x[,1]) - ((bdryKnots[2] - intKnots[ik]) / (bdryKnots[2] - bdryKnots[1])) * Ex_xL3 + 
              ((bdryKnots[1] - intKnots[ik]) / (bdryKnots[2] - bdryKnots[1])) * Ex_xU3 + 
              3 * (bdryKnots[2] - intKnots[ik]) * (intKnots[ik] - bdryKnots[1]) * Ex - qtmp * Ehxk
          }
        }else{}
        
        X[,1+nIntKnots] <- Ehxk / denomTmp
        
      }else{
        ### pt-supp vs...      
        tmp <- (x - bdryKnots[1])
        tmp[tmp < 0] <- 0
        x_xL3 <- (tmp ^ 3) 
        
        tmp <- (x - bdryKnots[2])
        tmp[tmp < 0] <- 0
        x_xU3 <- (tmp ^ 3) 
        
        tmp <- (x - intKnots[nIntKnots])
        tmp[tmp < 0] <- 0
        x_xK3 <- (tmp ^ 3) 
        
        hxk <- x_xK3 - ((bdryKnots[2] - intKnots[nIntKnots]) / (bdryKnots[2] - bdryKnots[1])) * x_xL3 + 
          ((bdryKnots[1] - intKnots[nIntKnots]) / (bdryKnots[2] - bdryKnots[1])) * x_xU3 + 
          3 * (bdryKnots[2] - intKnots[nIntKnots]) * (intKnots[nIntKnots] - bdryKnots[1]) * x
        
        denomTmp <- ((bdryKnots[2] - intKnots[nIntKnots]) * (intKnots[nIntKnots] - bdryKnots[1]) * (intKnots[nIntKnots] + bdryKnots[1] + bdryKnots[2]))
        
        X[,1] <- 1 - hxk / denomTmp
        
        if(nIntKnots > 1){
          for (ik in 1:(nIntKnots-1)){
            tmp <- (x - intKnots[ik])
            tmp[tmp < 0] <- 0
            qtmp <- ((bdryKnots[2] - intKnots[ik]) * (intKnots[ik] - bdryKnots[1]) * (intKnots[ik] + bdryKnots[1] + bdryKnots[2])) / denomTmp

            X[,1+ik] <- tmp ^ 3 - ((bdryKnots[2] - intKnots[ik]) / (bdryKnots[2] - bdryKnots[1])) * x_xL3 + 
              ((bdryKnots[1] - intKnots[ik]) / (bdryKnots[2] - bdryKnots[1])) * x_xU3 + 
              3 * (bdryKnots[2] - intKnots[ik]) * (intKnots[ik] - bdryKnots[1]) * x - qtmp * hxk
          }
        }else{}

        X[,1+nIntKnots] <- hxk / denomTmp

      }
    }else{}
    
        
  }else{
    stop('Unknonwn parameterisation option for spline function!')
  }

  return(X)  
}

####################################################################################
### fn to make basis fn for natural cubic spline with gradient clamped to 0 beyond upper bdryKnot and below lower bdryKnot
### if intKnots is empty, then error.
### if intKnots length 1, will give 1 col of 1s (ie a constant would be fitted)
####################################################################################
makeXnsClampedLUG0 <- function(x , bdryKnots , intKnots){

  nIntKnots <- length(intKnots)
  
  if(nIntKnots == 0){ stop('Attention - for double-clamped natural cubic spline, <=1 internal knot will give constant model)') }else{}
  if(nIntKnots == 0){ intKnots <- mean(bdryKnots) }else{} # just in case above line were to get turned off.

  if(is.null(dim(x)) || (ncol(x) == 1)){
    ### x just vector of known vals of x    
    n <- length(x)
    iaX <- FALSE
  }else{
    ### x should have 2 cols, first is lb, 2nd is ub for x
    n <- nrow(x)
    iaX <- TRUE
  }
  
  X <- matrix(1 , n , nIntKnots)
  
  if(nIntKnots > 1){
    if(iaX){
      ### inc-avd vs...      
      tmpL <- (x[,1] - bdryKnots[1])
      tmpL[tmpL < 0] <- 0
      tmpU <- (x[,2] - bdryKnots[1])
      tmpU[tmpU < 0] <- 0
      Ex_xL3 <- 0.25 * (tmpU ^ 4 - tmpL ^ 4) / (x[,2] - x[,1])
      
      tmpL <- (x[,1] - bdryKnots[2])
      tmpL[tmpL < 0] <- 0
      tmpU <- (x[,2] - bdryKnots[2])
      tmpU[tmpU < 0] <- 0
      Ex_xU3 <- 0.25 * (tmpU ^ 4 - tmpL ^ 4) / (x[,2] - x[,1])
      
      tmpL <- (x[,1] - intKnots[nIntKnots])
      tmpL[tmpL < 0] <- 0
      tmpU <- (x[,2] - intKnots[nIntKnots])
      tmpU[tmpU < 0] <- 0
      Ex_xK3 <- 0.25 * (tmpU ^ 4 - tmpL ^ 4) / (x[,2] - x[,1])
      
      for (ik in 1:(nIntKnots-1)){
        tmpL <- (x[,1] - intKnots[ik])
        tmpL[tmpL < 0] <- 0
        tmpU <- (x[,2] - intKnots[ik])
        tmpU[tmpU < 0] <- 0
        
        X[,1+ik] <- 0.25 * (tmpU ^ 4 - tmpL ^ 4) / (x[,2] - x[,1]) - (((bdryKnots[2] - intKnots[ik]) * (bdryKnots[1] - intKnots[ik])) / ((bdryKnots[2] - intKnots[nIntKnots]) * (bdryKnots[1] - intKnots[nIntKnots]))) * Ex_xK3 + 
          (((bdryKnots[2] - intKnots[ik]) * (intKnots[nIntKnots] - intKnots[ik])) / ((bdryKnots[2] - bdryKnots[1]) * (bdryKnots[1] - intKnots[nIntKnots]))) * Ex_xL3 + 
          (((bdryKnots[1] - intKnots[ik]) * (intKnots[ik] - intKnots[nIntKnots])) / ((bdryKnots[2] - bdryKnots[1]) * (bdryKnots[2] - intKnots[nIntKnots]))) * Ex_xU3
      }
    }else{
      ### pt-supp vs...      
      tmp <- (x - bdryKnots[1])
      tmp[tmp < 0] <- 0
      x_xL3 <- (tmp ^ 3) 
      
      tmp <- (x - bdryKnots[2])
      tmp[tmp < 0] <- 0
      x_xU3 <- (tmp ^ 3) 

      tmp <- (x - intKnots[nIntKnots])
      tmp[tmp < 0] <- 0
      x_xK3 <- (tmp ^ 3) 
      
      for (ik in 1:(nIntKnots-1)){
        tmp <- (x - intKnots[ik])
        tmp[tmp < 0] <- 0
        X[,1+ik] <- tmp ^ 3 - (((bdryKnots[2] - intKnots[ik]) * (bdryKnots[1] - intKnots[ik])) / ((bdryKnots[2] - intKnots[nIntKnots]) * (bdryKnots[1] - intKnots[nIntKnots]))) * x_xK3 + 
          (((bdryKnots[2] - intKnots[ik]) * (intKnots[nIntKnots] - intKnots[ik])) / ((bdryKnots[2] - bdryKnots[1]) * (bdryKnots[1] - intKnots[nIntKnots]))) * x_xL3 + 
          (((bdryKnots[1] - intKnots[ik]) * (intKnots[ik] - intKnots[nIntKnots])) / ((bdryKnots[2] - bdryKnots[1]) * (bdryKnots[2] - intKnots[nIntKnots]))) * x_xU3
      }
    }
  }else{}
  
  return(X)  
}


#####################################################
### fn to make des mtx...
### fn nanme not great, as this is for multiple variables, whereas the other makeX fns above were for univariate. 
### update naming sometime.
#####################################################
makeXcns <- function(dfCovs , dIData = NULL , listfefdKnots , incInts = TRUE , colnamesX = NULL , intMthd = 0){
  
  if(length(listfefdKnots) == 0){
    X <- matrix(1 , nrow(dfCovs) , 1)
    colnames(X) <- 'const'
    return(X)
  }
  rcondLim <- 1E-12 # if rcond < this lim, will be determined as too close to singular, so cols will not be included.
  
  #####################
  ### if colnamesX is not given, then we are making X with the fitting dataset - so colin cols have to be removed
  ### if colnamesX is given, then we have already made with the fitting dataset, so just include those named cols 
  #####################
  ### if incInts = TRUE, then include all interactions of basis fns.
  ### if incInts = list of two vectors, then all variables in first vector interact with all variables in the second vector
  #####################
  ### if variable v has basis fns v1 - v5 and w has w1 - w5
  ### intMthd = 0 : interactions between v and w are v1*w0, v2*w0, ..., v5*w0, v0*w1, v0*w2,..., v0*w5 (ie 10 basis fns)    
  ### where v0 is spline type basis of v with 1 df, same for w0
  ### intMthd = 1 : interactions between v and w are all vi * wj (ie 25 basis fns)
  #####################
  covNames <- names(listfefdKnots)
  
  listXMain <- list()
  if(intMthd == 0){ listXMain0 <- list() }else{}
  for (nm in covNames){
    
    if((nm == 'dIMidPts') & (!is.null(dIData))){
      xThis <- dIData
    }else{
      xThis <- dfCovs[[nm]]
    }
    
    if(listfefdKnots[[nm]]$sType == 'nscug'){
      listXMain[[nm]] <- makeXnsClampedUG0(x = xThis , bdryKnots = listfefdKnots[[nm]]$bdryKnots , intKnots = listfefdKnots[[nm]]$intKnots , paramVs = 0)[,-1,drop=FALSE]
      colnames(listXMain[[nm]]) <- paste0(nm , '_KNOT' , seq(ncol(listXMain[[nm]])))
      if(intMthd == 0){ 
### spline with 1 basis fn requires 1 int knot...        
        listXMain0[[nm]] <- makeXnsClampedUG0(x = xThis , bdryKnots = listfefdKnots[[nm]]$bdryKnots , intKnots = median(listfefdKnots[[nm]]$intKnots) , paramVs = 0)[,-1,drop=FALSE]
        colnames(listXMain0[[nm]]) <- paste0(nm , '_KNOT0')
      }else{}
    }else if(listfefdKnots[[nm]]$sType == 'nsclug'){
      listXMain[[nm]] <- makeXnsClampedLUG0(x = xThis , bdryKnots = listfefdKnots[[nm]]$bdryKnots , intKnots = listfefdKnots[[nm]]$intKnots)[,-1,drop=FALSE]
      colnames(listXMain[[nm]]) <- paste0(nm , '_KNOT' , seq(ncol(listXMain[[nm]])))
      if(intMthd == 0){ 
### spline with 1 basis fn requires 2 int knots...        
        listXMain0[[nm]] <- makeXnsClampedLUG0(x = xThis , bdryKnots = listfefdKnots[[nm]]$bdryKnots , intKnots = quantile(listfefdKnots[[nm]]$intKnots , c(1/3,2/3)))[,-1,drop=FALSE]
        colnames(listXMain0[[nm]]) <- paste0(nm , '_KNOT0')
      }else{}
    }else if(listfefdKnots[[nm]]$sType == 'ns'){
      listXMain[[nm]] <- makeXns(x = xThis , bdryKnots = listfefdKnots[[nm]]$bdryKnots , intKnots = listfefdKnots[[nm]]$intKnots)[,-1,drop=FALSE]
      colnames(listXMain[[nm]]) <- paste0(nm , '_KNOT' , seq(ncol(listXMain[[nm]])))
      if(intMthd == 0){ 
### spline with 1 basis fn requires 0 int knots...        
        listXMain0[[nm]] <- makeXns(x = xThis , bdryKnots = listfefdKnots[[nm]]$bdryKnots , intKnots = c())[,-1,drop=FALSE]
        colnames(listXMain0[[nm]]) <- paste0(nm , '_KNOT0')
      }else{}
    }else if(listfefdKnots[[nm]]$sType == 'cat'){
      # qwe2 = model.matrix(as.formula(paste0('~' , nm)) , data = dfCovs , xlev = 'C1')[,-1]
### will use the same naming for categorical variable, using _KNOT1, etc.
### when testing for drops later, will check if cts/cat, then decide what is appropriate to test for dropping from model...
      fnTmp <- function(i){ vecRtn <- numeric(nrow(dfCovs)) ; vecRtn[dfCovs[[nm]] == listfefdKnots[[nm]]$levels[i]] <- 1 ; return(vecRtn) }
      listXMain[[nm]] <- do.call(cbind , lapply(seq(2,length(listfefdKnots[[nm]]$levels)) , fnTmp))
      colnames(listXMain[[nm]]) <- paste0(nm , '_KNOTC' , seq(2 , ncol(listXMain[[nm]]) + 1)) # for cats, use _KNOTSC[2 - nlevels]
    }else{
      stop('Error - unrecognised sType entered into makeXcns!')
    }
  }
    
  ### mtcs of interactive effects...
  if(is.list(incInts) || incInts){
    listXInt <- list()
    for (inm in 1:length(covNames)){
      for(jnm in inm:length(covNames)){
        
        if(is.list(incInts)){
          incThisInt <- (inm != jnm) & ((is.element(covNames[inm] , incInts[[1]]) & is.element(covNames[jnm] , incInts[[2]])) | (is.element(covNames[inm] , incInts[[2]]) & is.element(covNames[jnm] , incInts[[1]])))
        }else{
          incThisInt <- (inm != jnm)
        }
        
        if(incThisInt){
          if(intMthd == 0){ 
### a mtx with cols = products of all cols of X1 with simplest basis fn of X2, and all cols of X2 with simplest basis fn of X1  
            XIntTest <- cbind(listXMain[[inm]] * matrix(listXMain0[[jnm]] , nrow(listXMain[[inm]]) , ncol(listXMain[[inm]])) , 
                              listXMain[[jnm]] * matrix(listXMain0[[inm]] , nrow(listXMain[[jnm]]) , ncol(listXMain[[jnm]])))
            colnames(XIntTest) <- c(paste0(colnames(listXMain[[inm]]) , '___' , colnames(listXMain0[[jnm]])) ,
                                    paste0(colnames(listXMain[[jnm]]) , '___' , colnames(listXMain0[[inm]])))
          }else{
### a mtx with cols = all products of cols in X1 and X2  
            XIntTest <- do.call(cbind , lapply(lapply(seq_len(ncol(listXMain[[inm]])), function(i) listXMain[[inm]][,i]) , function(mtx){ mtx * listXMain[[jnm]] }))
            colnames(XIntTest) <- paste0(rep(colnames(listXMain[[inm]]) , each = ncol(listXMain[[jnm]])) , '___' , rep(colnames(listXMain[[jnm]]) , ncol(listXMain[[inm]])))
          }
          
          if(is.null(colnamesX)){
            ### test whether data allow interactive basis fns to be included...
            XTest <- cbind(matrix(1 , nrow(dfCovs) , 1) , listXMain[[inm]] , listXMain[[jnm]])
            betaTest <- try(solve(t(XTest) %*% XTest , matrix(1 , ncol(XTest) , 1)), silent = TRUE)
            if(is.character(betaTest)){ stop(paste0('Error - system singular with basis fns ' , covNames[inm] , ' and ' , covNames[jnm])) }else{}

            rcondTest <- rcond(t(XTest) %*% XTest)
            if(rcondTest < rcondLim){ stop(paste0('Error - system rcond too small with basis fns ' , covNames[inm] , ' and ' , covNames[jnm])) }else{}

            incCol <- logical(ncol(XIntTest))
            for(icolAdd in 1:ncol(XIntTest)){
              XTestThis <- cbind(XTest , XIntTest[,icolAdd,drop=FALSE])
              betaTest <- try(solve(t(XTestThis) %*% XTestThis , matrix(1 , ncol(XTestThis) , 1)), silent = TRUE)
              rcondTest <- rcond(t(XTestThis) %*% XTestThis)

              if(is.character(betaTest) | (rcondTest < rcondLim)){ 
                print(paste0('Attention - system singular when interaction added for ' , colnames(XIntTest)[icolAdd] , '! So not including this column.')) 
                incCol[icolAdd] <- FALSE
              }else{
                incCol[icolAdd] <- TRUE
                XTest <- XTestThis
              }
            }
            
          }else{
            incCol <- is.element(colnames(XIntTest) , colnamesX)
          }
          
          listXInt[[paste0(covNames[inm] , '_' , covNames[jnm])]] <- XIntTest[,incCol,drop=FALSE]
        }else{}
      }
    }
  }else{}
  
  ### X for main effects...
  Xcns <- cbind(matrix(1 , nrow(dfCovs) , 1) , do.call(cbind , listXMain))
  
  ### X for interactive effects...
  if(is.list(incInts) || incInts){
    Xcns <- cbind(Xcns , do.call(cbind , listXInt))
  }else{}
  
  colnames(Xcns)[1] <- 'const'
  
  if(is.null(colnamesX)){
    # print('Not doing any checks of colinearity in the full mtx.')
    ### another check for colin with the full mtx...
    colsInc <- rep(TRUE , ncol(Xcns))
    for(i in 1:ncol(Xcns)){
      XXTest <- t(Xcns[,seq(1,i)[colsInc[1:i]],drop=FALSE]) %*% Xcns[,seq(1,i)[colsInc[1:i]],drop=FALSE]
      betaTest <- try(solve(XXTest , matrix(1 , sum(colsInc[1:i]) , 1)), silent = TRUE)
      rcondTest <- rcond(XXTest)
      if(is.character(betaTest) | (rcondTest < rcondLim)){
        print(paste0('Colinearity found when adding column ' , colnames(Xcns)[i] , ', so this will not be included.'))
        colsInc[i] <- FALSE
      }else{}
    }
    Xcns <- Xcns[,colsInc,drop=FALSE]
    
  }else{
    ### only keep the rqd cols...  
    Xcns <- Xcns[,is.element(colnames(Xcns) , colnamesX),drop=FALSE]
  }
  
  return(Xcns)  
}

#####################################################
### fn to make list, each element is a set of cols that could be tested for dropping...
#####################################################
getPossibleReductions <- function(XFull , allowKnotRemoval = TRUE){
  varsWithKnots <- unique(unlist(strsplit(colnames(XFull) , '___')))
  varsWithKnots <- varsWithKnots[which(grepl('_KNOT' , varsWithKnots))]
  varsWithoutKnots <- unique(unlist(lapply(strsplit(varsWithKnots , '_KNOT') , '[[' , 1)))

  catvarsWithKnots <- varsWithKnots[which(grepl('_KNOTC' , varsWithKnots))]
  ctsvarsWithKnots <- setdiff(varsWithKnots , catvarsWithKnots)
  catvarsWithoutKnots <- unique(unlist(lapply(strsplit(catvarsWithKnots , '_KNOT') , '[[' , 1)))
  ctsvarsWithoutKnots <- setdiff(varsWithoutKnots , catvarsWithoutKnots)
  
  ### add to list cols with poss knots to rm
  if(allowKnotRemoval){
    fnTmp <- function(vk){ which(grepl(vk , colnames(XFull))) }
    listRm <- lapply(ctsvarsWithKnots , fnTmp)
  }else{
    listRm <- list()
  }

  ### add to list cols to remove a whole covariate from X...
  fnTmp <- function(nm){ which(grepl(paste0(nm , '_KNOT') , colnames(XFull))) }
  lTmp <- lapply(varsWithoutKnots , fnTmp)
  lTmp <- lTmp[unlist(lapply(lTmp , length)) > 0]
  if(length(lTmp) > 0){
    listRm <- c(listRm , lTmp)
  }else{}
  
  ### add to list cols to remove a single basis fn from an interaction between 2 basis fns of cts vars...
  if(length(catvarsWithKnots) > 0){
    lTmp <- as.list(which((unlist(lapply(strsplit(colnames(XFull) , '___') , length)) == 2) & (!grepl(paste(catvarsWithKnots , collapse = '|') , colnames(XFull)))))
  }else{
    lTmp <- as.list(which(unlist(lapply(strsplit(colnames(XFull) , '___') , length)) == 2))
  }
  lTmp <- lTmp[unlist(lapply(lTmp , length)) > 0]
  if(length(lTmp) > 0){
    listRm <- c(listRm , lTmp)
  }else{}

  ### add to list cols with any interactions to fully rm...
  if(length(varsWithoutKnots) > 1){
    listVarInts <- combn(varsWithoutKnots , 2 , simplify = FALSE)
    fnTmp <- function(vVec){ which(grepl(paste0(vVec[1] , '_KNOT') , colnames(XFull)) & grepl(paste0(vVec[2] , '_KNOT') , colnames(XFull))) }
    
    lTmp <- lapply(listVarInts , fnTmp)
    lTmp <- lTmp[unlist(lapply(lTmp , length)) > 0]
    if(length(lTmp) > 0){
      listRm <- c(listRm , lTmp)
    }else{}
  }else{}
  
  return(listRm)  
}

#####################################################
### fn to perform backward stepwise selection, using Wald tests...
### permissible drops include 1. any knot; 2. any basis column from an interaction; 3. any complete interaction.
#####################################################
stepBackXcns <- function(Xcns , zData , alpha = 0.15 , iAX = NULL , iAz = NULL){
  
  ### in phase 1, try dropping basis fns of interactions, or complete interactions...
  ### in phase 2, also allow complete removal of knots...
  # listInfo <- list()
  
  if(!is.null(iAX)){ 
    ziAz <- t(zData) %*% iAz 
  }else{
    zz <- sum(zData ^ 2)
  }
  
  continueDropping <- TRUE
  allowKnotRemoval <- FALSE
  while(continueDropping){
    
    if(is.null(iAX)){
      # tmp <- nllLm(zData = zData , XData = Xcns , REML = TRUE)
      # betahat <- tmp$betahat
      # vbetahat <- tmp$vbetahat

###      
      XX <- t(Xcns) %*% Xcns
      Xz <- t(Xcns) %*% zData
      betahat <- try(solve(XX , Xz) , silent = TRUE)

### add small increment to diag of XX (all cols of X should be stdzd), until system solvable...            
      ridgeParAdd <- 1E-12
      diagXX0 <- diag(XX)
      while(is.character(betahat)){
        diag(XX) <- diagXX0 + ridgeParAdd
        betahat <- try(solve(XX , Xz) , silent = TRUE)
        if(is.character(betahat)){ ridgeParAdd <- ridgeParAdd * 10 }else{}
      }

      sigma2hat <- as.numeric(zz - 2 * t(betahat) %*% Xz + t(betahat) %*% XX %*% betahat) / (nrow(Xcns) - ncol(Xcns))
      vbetahat <- sigma2hat * solve(XX)
      
    }else{
      ### these things could be done more efficiently, but for now...      
      XiAX <- t(Xcns) %*% iAX
      XiAz <- t(Xcns) %*% iAz
      betahat <- try(solve(XiAX , XiAz) , silent = TRUE)
      
      ### keep adding small increments to diag of XX (all cols of X should be stdzd), until system solvable...            
      ridgeParAdd <- 1E-12
      while(is.character(betahat)){
        diag(XiAX) <- diag(XiAX) + ridgeParAdd
        betahat <- try(solve(XiAX , XiAz) , silent = TRUE)
      }
      
      sigma2hat <- as.numeric(ziAz - 2 * t(betahat) %*% XiAz + t(betahat) %*% XiAX %*% betahat) / (nrow(Xcns) - ncol(Xcns))
      vbetahat <- sigma2hat * solve(XiAX)
    }
    
    listPossRedns <- getPossibleReductions(Xcns , allowKnotRemoval = allowKnotRemoval)
    
    pVec <- unlist(lapply(listPossRedns , waldTest , betahat = betahat , vbetahat = vbetahat))
    # infoThis <- list(head(Xcns) , pVec)
    # listInfo <- c(listInfo , infoThis)
    if(max(pVec) > alpha){
      imax <- which.max(pVec)
      Xcns <- Xcns[,-listPossRedns[[imax]],drop=FALSE]
      if(!is.null(iAX)){ 
        iAX <- iAX[,-listPossRedns[[imax]],drop=FALSE]
      }else{}
    }else{
      if(allowKnotRemoval){
        continueDropping <- FALSE
      }else{
        ### switch to phase 2, but don't stop yet...
        allowKnotRemoval <- TRUE
      }
    }
  }

### p values from final tests for each individual basis fn...
  pVec <- unlist(lapply(as.list(seq(nrow(vbetahat))) , waldTest , betahat = betahat , vbetahat = vbetahat))
  names(pVec) <- colnames(Xcns)
  
  # return(list('Xcns' = Xcns , 'betahat' = betahat , 'vbetahat' = vbetahat , 'listInfo' = listInfo))
  return(list('Xcns' = Xcns , 'betahat' = betahat , 'vbetahat' = vbetahat , 'pVec' = pVec))
}

##########################################################################
### add scaled versions of covariates to fitting and calibration data...
##########################################################################
addScaledCovs2df <- function(dfFit , dfPred = NULL , covNames){
  for(nm in covNames){
    if(substr(nm , nchar(nm) - 6 , nchar(nm)) == '_SCALED'){
      
### get and check the original name, which must be in the df...      
      nmunscaled <- substr(nm , 1 , nchar(nm) - 7)
      if(!is.element(nmunscaled , names(dfFit))){ stop(paste0('Error - ' , nmunscaled , ' not found in dfFit!')) }else{}
      if(!is.null(dfPred)){
        if(!is.element(nmunscaled , names(dfPred))){ stop(paste0('Error - ' , nmunscaled , ' not found in dfPred!')) }else{}
      }else{}
      
      mTmp <- mean(dfFit[[nmunscaled]])
      sdTmp <- sd(dfFit[[nmunscaled]])
      dfFit[[nm]] <- (dfFit[[nmunscaled]] - mTmp) / sdTmp
      if(!is.null(dfPred)){
        dfPred[[nm]] <- (dfPred[[nmunscaled]] - mTmp) / sdTmp
      }else{}
      
    }else{
      if(!is.element(nm , names(dfFit))){ stop(paste0('Error - ' , nm , ' not found in dfFit!')) }else{}
      if(!is.null(dfPred)){
        if(!is.element(nm , names(dfPred))){ stop(paste0('Error - ' , nm , ' not found in dfPred!')) }else{}
      }else{}
    }
  }
  
  return(list('dfFit' = dfFit , 'dfPred' = dfPred))
}

##############################################################
### make a list with initial knots for each covariate...
##############################################################
makelistfefdKnots <- function(dfFit , covNames , q4BdryKnots = c(0.05 , 0.95) , nIntKnots = 5 , sType = 'nscug'){

  ############################################
  ### sType is the type of spline - can be:
  ###             nscug (natural spline, clamped to have gradient 0 at upper boundary knot)  
  ###             nsclug (natural spline, clamped to have gradient 0 at upper boundary knot and at lower boundary knot)  
  ###             ns (natural spline, no clamping)  
  ############################################
  listfefdKnots <- list()
  if(length(covNames) > 0){
    for(inm in 1:length(covNames)){
      nm <- covNames[inm]
      fefdKnotsTmp <- list()
      
      if(is.null(dfFit[[nm]])){ stop(paste0('Error - ' , nm , ' not found in data.frame!')) }else{}
      
      if(is(dfFit[[nm]] , 'numeric')){
        ### knots for a cts variable...      
        if(length(sType) > 1){
          if(length(sType) != length(covNames)){ stop('Error - sType entered as vector with wrong length - should be same length as covNames!') }else{}
          fefdKnotsTmp$sType <- sType[inm]
        }else{
          fefdKnotsTmp$sType <- sType
        }
        if(is.matrix(q4BdryKnots)){
          if((ncol(q4BdryKnots) != 2) | (nrow(q4BdryKnots) != length(covNames))){ stop('Error - q4BdryKnots entered as matrix with wrong dimensions - should have 2 cols and number of rows = length(covNames)!') }else{}
          fefdKnotsTmp$bdryKnots <- quantile(dfFit[[nm]] , q4BdryKnots[inm,])
        }else{
          fefdKnotsTmp$bdryKnots <- quantile(dfFit[[nm]] , q4BdryKnots)
        }
        if(length(nIntKnots) > 1){
          if(length(nIntKnots) != length(covNames)){ stop('Error - nIntKnots entered as vector with wrong length - should be same length as covNames!') }else{}
          fefdKnotsTmp$intKnots <- getQuantileKnots(dfFit[[nm]] , bdryKnots = fefdKnotsTmp$bdryKnots , nIntKnots = nIntKnots[inm])
        }else{
          fefdKnotsTmp$intKnots <- getQuantileKnots(dfFit[[nm]] , bdryKnots = fefdKnotsTmp$bdryKnots , nIntKnots = nIntKnots)
        }
        
        ### check for too few knots...      
        if((fefdKnotsTmp$sType == 'nsclug') & (length(fefdKnotsTmp$intKnots) < 2)){ stop('Error - nsclug (clamped grad at upper and lower boundary knots) should have at least 2 internal knots, otherwise it is constant.') }else{}
        if((fefdKnotsTmp$sType == 'nscug') & (length(fefdKnotsTmp$intKnots) < 1)){ stop('Error - nscug (clamped grad at upper boundary knot) should have at least 1 internal knot, otherwise it is constant.') }else{}
        
      }else{
        ### cats to make basis fns for a cat variable...      
        if((length(sType) > 1) & (sType[inm] != 'cat')){ stop('Error - sType was not entered as cat for a factor variable - check this!') }else{}
        fefdKnotsTmp$sType <- 'cat'
        fefdKnotsTmp$levels <- levels(dfFit[[nm]])
      }      
      
      listfefdKnots[[nm]] <- fefdKnotsTmp
    }
  }else{}
  
  return(listfefdKnots)  
}

####################################################################################
### fn to make basis fn for semi-natural cubic spline with gradient clamped to 0 beyond upper bdryKnot
### semi-natural because only upper bdry included
####################################################################################
# makeXsnsClampedUG0 <- function(x , bdryKnotU , intKnots){
#   
#   if(is.null(dim(x)) || (ncol(x) == 1)){
#     ### x just vector of known vals of x    
#     n <- length(x)
#     iaX <- FALSE
#   }else{
#     ### x should have 2 cols, first is lb, 2nd is ub for x
#     n <- nrow(x)
#     iaX <- TRUE
#   }
#   
#   X <- matrix(1 , n , 2 + length(intKnots))
#   
#   if(length(intKnots) > 0){
#     if(iaX){
#       
#       Ex <- rowMeans(x)
#       Ex2 <- (x[,2] ^ 3 - x[,1] ^ 3) / (3 * (x[,2] - x[,1]))
#       Ex3 <- (x[,2] ^ 4 - x[,1] ^ 4) / (4 * (x[,2] - x[,1]))
#       
#       ### inc-avd vs...      
#       tmpL <- (x[,1] - bdryKnotU)
#       tmpL[tmpL < 0] <- 0
#       tmpU <- (x[,2] - bdryKnotU)
#       tmpU[tmpU < 0] <- 0
# 
#       X[,2] <- 0.25 * (tmpU ^ 4 - tmpL ^ 4) / (x[,2] - x[,1]) + 3 * bdryKnotU * Ex2 - Ex3 - 3 * (bdryKnotU ^ 2) * Ex
# 
#       for (ik in 1:length(intKnots)){
#         tmpL <- (x[,1] - intKnots[ik])
#         tmpL[tmpL < 0] <- 0
#         tmpU <- (x[,2] - intKnots[ik])
#         tmpU[tmpU < 0] <- 0
#         
#         X[,2+ik] <- 0.25 * (tmpU ^ 4 - tmpL ^ 4) / (x[,2] - x[,1]) + 3 * intKnots[ik] * Ex2 - Ex3 - 3 * (intKnots[ik] ^ 2) * Ex
#       }
#     }else{
#       ### pt-supp vs...      
#       tmp <- (x - bdryKnots[2])
#       tmp[tmp < 0] <- 0
#       alphaU <- (tmp ^ 3) / (bdryKnots[2] - bdryKnots[1])
#       X[,2] <- tmp ^ 3 + 3 * bdryKnotU * (x ^ 2) - (x ^ 3) - 3 * (bdryKnotU ^ 2) * x
#       
#       for (ik in 1:length(intKnots)){
#         tmp <- (x - intKnots[ik])
#         tmp[tmp < 0] <- 0
#         X[,2+ik] <- tmp ^ 3 + 3 * intKnots[ik] * (x ^ 2) - (x ^ 3) - 3 * (intKnots[ik] ^ 2) * x
#       }
#     }
#   }else{}
#   
#   return(X)  
# }


