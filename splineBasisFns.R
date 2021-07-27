####################################################################################
####################################################################################
### AS SET UP HERE, TO FIT, ALL BDRYKNOTS AND INTKNOTS MUST BE WITHIN RANGE OF DATA 
### (IE AT LEAST ONE DATA PT BEYOND BDRYKNOTS)
### FOR FUTURE - ADD DEGREE AS ARGUMENT IN ALL BASIS FUNCTION FUNCTIONS.
### DEGREE = 1 (PIECEWISE LINEAR) MIGHT BE QUITE USEFUL
####################################################################################
####################################################################################


####################################################################################
### fn to get positions for knots based on equal quantiles of data within bdryKnots
### (or if bdry knots outside range of data, then within the range of data)
####################################################################################
getQuantileKnots <- function(x , bdryKnots , nIntKnots){
  if (nIntKnots == 0){
    return(c())
  }else{
    minx <- min(x)
    maxx <- max(x)
    if(bdryKnots[1] < minx){ 
      xlb <- minx
    }else{
      xlb <- bdryKnots[1]
    }
    if(bdryKnots[2] > maxx){ 
      xub <- maxx
    }else{
      xub <- bdryKnots[2]
    }
    
    xIn <- x[x >= xlb & x <= xub]
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
  
  if(nIntKnots == 2 & bdryKnots[1] == -Inf & bdryKnots[2] == Inf){
### way to code a linear function. will ignore internal knot positions.
    if(iaX){
      X[,2] <- 0.5 * (x[,1] + x[,2])
    }else{
      X[,2] <- x 
    }
    return(X)  
  }else{}
  
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
### fn name not great, as this is for multiple variables, whereas the other makeX fns above were for univariate. 
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
      if(intMthd == 0){
        listXMain0[[nm]] <- listXMain[[nm]]
        colnames(listXMain0[[nm]]) <- paste0(nm , '_KNOTC0L' , seq(ncol(listXMain[[nm]]))) # should not be used, inserted as place holder as idx in list used.
      }else{}
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
          if((intMthd == 0) & (listfefdKnots[[covNames[inm]]]$sType != 'cat') & (listfefdKnots[[covNames[jnm]]]$sType != 'cat') & (ncol(listXMain[[inm]]) > 1) & (ncol(listXMain[[jnm]]) > 1)){ 
            ### cat variables or variables defined with 1 basis fn are not included in this way, but done using all prods in the else bit.
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
          
          # if(is.null(colnamesX)){
          #   incCol <- logical(ncol(XIntTest))
          # }else{
          #   incCol <- is.element(colnames(XIntTest) , colnamesX)
          # }
          # 
          # listXInt[[paste0(covNames[inm] , '_' , covNames[jnm])]] <- XIntTest[,incCol,drop=FALSE]
          
          listXInt[[paste0(covNames[inm] , '_' , covNames[jnm])]] <- XIntTest
        }else{}
      }
    }
  }else{}
  
  ### X for main effects...
  Xcns <- cbind(matrix(1 , nrow(dfCovs) , 1) , do.call(cbind , listXMain))
  colnames(Xcns)[1] <- 'const'
  
  ### check the main effects, and remove cols as required...
  if(is.null(colnamesX)){
    
    tmp <- reduceXUsingrcond(XIn = Xcns , rqr_p_LT_n = TRUE , rcondLim = rcondLim)
    Xcns <- tmp$XOut  
    colsDropFromMainEffs <- tmp$colsDropped

  }else{
    Xcns <- Xcns[,is.element(colnames(Xcns) , colnamesX),drop=FALSE]
  }

  ### X for interactive effects...
  if(is.list(incInts) || incInts){
    pMain <- ncol(Xcns)

    # Xcns <- cbind(Xcns , do.call(cbind , listXInt))
    XcnsInt <- do.call(cbind , listXInt)
    # colnames(Xcns)[1] <- 'const'
    
    if(is.null(colnamesX)){
      
      ### drop from int mtx any with same name as dropped from main effs
      ### (ie knot dropped from main -> drop from int)
      if(length(colsDropFromMainEffs) > 0){
        iDrop <- which(grepl(colsDropFromMainEffs , colnames(XcnsInt)))
        if(length(iDrop) > 0){
          XcnsInt <- XcnsInt[,-iDrop,drop=FALSE]
        }else{}
      }else{}
      
      Xcns <- cbind(Xcns , XcnsInt)
      
      tmp <- reduceXUsingrcond(XIn = Xcns , rqr_p_LT_n = TRUE , rcondLim = rcondLim)
      Xcns <- tmp$XOut  
      colsDropFromIntEffs <- tmp$colsDropped
      
      # XcnsTest <- Xcns
      # 
      # ### stdise...
      # mX <- colMeans(XcnsTest[,-1,drop=FALSE])
      # sdX <- apply(XcnsTest[,-1,drop=FALSE] , 2 , sd)
      # XcnsTest[,-1] <- (XcnsTest[,-1,drop=FALSE] - matrix(mX , nrow(XcnsTest) , ncol(XcnsTest)-1 , byrow = TRUE)) / matrix(sdX , nrow(XcnsTest) , ncol(XcnsTest)-1 , byrow = TRUE)
      # 
      # XXTest <- t(XcnsTest) %*% XcnsTest
      # rcondTest <- rcond(XXTest)
      # colsDropFromMainEffs <- c()
      # while((rcondTest < rcondLim) & (ncol(XXTest) > pMain)){
      #   ### drop the column from X that gives largest rcond...      
      #   rTmp <- unlist(lapply(seq((pMain+1):ncol(XXTest)) , rcondSubMtx , XXData = XXTest))
      #   
      #   iDrop <- which.max(rTmp)
      #   colsDropFromMainEffs <- c(colsDropFromMainEffs , colnames(Xcs)[iDrop])
      #   Xcns <- Xcns[,-iDrop,drop=FALSE]
      #   XcnsTest <- XcnsTest[,-iDrop,drop=FALSE]
      #   XXTest <- XXTest[-iDrop,-iDrop,drop=FALSE]
      #   
      #   # print(paste0('Colinearity found when adding column ' , colnames(Xcns)[i] , ', so this will not be included.'))
      #   # colsInc[i] <- FALSE
      # }
      
    }else{
      
      Xcns <- cbind(Xcns , XcnsInt)
      
      ### only keep the rqd cols...  
      Xcns <- Xcns[,is.element(colnames(Xcns) , colnamesX),drop=FALSE]
    }
    
  }else{}
  
  return(Xcns)  
}

#############################################################################
### drop rows/cols from XX and return rcond...
#############################################################################
rcondSubMtx <- function(colsDrop = c() , XXData){
  if(length(colsDrop) > 0){
    XXData <- XXData[-colsDrop,-colsDrop,drop=FALSE]
  }else{}
  return(rcond(XXData))  
}

###########################################################################
### stepwise greedy algorithm to remove columns of X so that t(X) %*% X is invertible
### based on increasing condition number. 
### Should reduce colinearity issues while retaining input variable meaning
### But could remove useful predictors 
###########################################################################
reduceXUsingrcond <- function(XIn , rqr_p_LT_n = TRUE , rcondLim = 1E-12){
  # if rqr_p_LT_n, then we need p < n so that we have a variance
  # else, we just need p <= n; would give prediction (interpolation) but not a variance.
  
  if(rqr_p_LT_n){
    pmax <- nrow(XIn) - 1
  }else{
    pmax <- nrow(XIn)
  }
  
  ### stdise...
  mX <- colMeans(XIn[,-1,drop=FALSE])
  sdX <- apply(XIn[,-1,drop=FALSE] , 2 , sd)
  XTest <- XIn
  XTest[,-1] <- (XIn[,-1,drop=FALSE] - matrix(mX , nrow(XIn) , ncol(XIn)-1 , byrow = TRUE)) / matrix(sdX , nrow(XIn) , ncol(XIn)-1 , byrow = TRUE)
  
  XXTest <- t(XTest) %*% XTest
  rcondTest <- rcond(XXTest)
  colsDropped <- c()
  while(((rcondTest < rcondLim) & (ncol(XXTest) > 1)) | (ncol(XXTest) > pmax)){
    ### drop the column from X that gives largest rcond...      
    rTmp <- unlist(lapply(seq(2:ncol(XXTest)) , rcondSubMtx , XXData = XXTest))
    
    icolDrop <- 1 + which.max(rTmp)
    colsDropped <- c(colsDropped , colnames(XTest)[icolDrop])
    XIn <- XIn[,-icolDrop,drop=FALSE]
    XTest <- XTest[,-icolDrop,drop=FALSE]
    XXTest <- XXTest[-icolDrop,-icolDrop,drop=FALSE]
    rcondTest <- rcond(XXTest)
  }
  
  return(list('XOut' = XIn, 'colsDropped' = colsDropped))  
}

#####################################################
### fn to make list, each element is a set of cols that could be tested for dropping...
#####################################################
getPossibleReductions <- function(XFull , allowKnotRemoval = TRUE , allowVariableRemoval = TRUE){
  varsWithKnots <- unique(unlist(strsplit(colnames(XFull) , '___')))
  varsWithKnots <- varsWithKnots[which(grepl('_KNOT' , varsWithKnots))]
  varsWithoutKnots <- unique(unlist(lapply(strsplit(varsWithKnots , '_KNOT') , '[[' , 1)))
  
  catvarsWithKnots <- varsWithKnots[which(grepl('_KNOTC' , varsWithKnots))]
  ctsvarsWithKnots <- setdiff(varsWithKnots , catvarsWithKnots)
### just the variable names, with the '_KNOTi' bit removed...
  catvarsWithoutKnots <- unique(unlist(lapply(strsplit(catvarsWithKnots , '_KNOT') , '[[' , 1)))
  ctsvarsWithoutKnots <- setdiff(varsWithoutKnots , catvarsWithoutKnots)
  
  ### add to list cols with poss knots to rm
  ### but make sure removing does not leave interaction term but no main effect.
  if(allowKnotRemoval){
    fnTmp <- function(vk){ which(grepl(vk , colnames(XFull))) }
    # listRm <- lapply(ctsvarsWithKnots , fnTmp)
    
    if(length(ctsvarsWithKnots) > 0){
      varNamesTmp <- unlist(lapply(strsplit(ctsvarsWithKnots , '_KNOT') , '[[' , 1))
      iKnotsTmp <- as.numeric(unlist(lapply(strsplit(ctsvarsWithKnots , '_KNOT') , '[[' , 2)))
      listRm <- list() ; iCount <- 1
      for(i in 1:length(ctsvarsWithKnots)){
        iTmp <- which(varNamesTmp == varNamesTmp[i])
### add to list unless KNOTi (i > 0) and exactly 2 instances of variable found, one with KNOTi and the other with KNOT0
        if(!((length(iTmp) == 2) && (iKnotsTmp[i] > 0) && ((iKnotsTmp[iTmp[1]] == 0) | (iKnotsTmp[iTmp[2]] == 0)))){
          listRm[[iCount]] <- fnTmp(ctsvarsWithKnots[i])
          iCount <- iCount + 1
        }else{}
      }      
    }else{
      listRm <- list()
    }
  }else{
    listRm <- list()
  }
  
  ### add to list cols to remove a whole covariate from X...
  if(allowVariableRemoval){
    fnTmp <- function(nm){ which(grepl(paste0(nm , '_KNOT') , colnames(XFull))) }
    lTmp <- lapply(varsWithoutKnots , fnTmp)
    lTmp <- lTmp[unlist(lapply(lTmp , length)) > 0]
    if(length(lTmp) > 0){
      listRm <- c(listRm , lTmp)
    }else{}
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
### fn to get the names of all spatial variables in the model...
#####################################################
getSpatVars <- function(colnamesX){
  varsWithKnots <- unique(unlist(strsplit(colnamesX , '___')))
  varsWithKnots <- varsWithKnots[which(grepl('_KNOT' , varsWithKnots))]
  varsWithoutKnots <- unique(unlist(lapply(strsplit(varsWithKnots , '_KNOT') , '[[' , 1)))
  spatVars <- setdiff(varsWithoutKnots , 'dIMidPts')
  
  return(spatVars)
}

#####################################################
### fn to get the number of knots for the fn of dIMidPts...
### (this is the number of columns with unique _KNOT numbers)
### (don't include _KNOT0 in the count)
#####################################################
getnKnotsd <- function(colnamesX){
  varsWithKnots <- unique(unlist(strsplit(colnamesX , '___')))
  
  varsWithKnots <- gsub('_KNOT0' , '_JUNKFORNOW' , varsWithKnots)

  varsWithKnots <- varsWithKnots[which(grepl('_KNOT' , varsWithKnots))]
  varsWithKnotsd <- varsWithKnots[which(substr(varsWithKnots , 1 , 13) == 'dIMidPts_KNOT')]
  nKnotsd <- length(unique(varsWithKnotsd))
  
  return(nKnotsd)
}

#####################################################
### fn to get the number of knots for the fn of each spatial covariate...
### (this is the number of columns with unique _KNOT numbers)
### (don't include _KNOT0 in the count)
#####################################################
getnKnotss <- function(colnamesX , allSpatVars){
  if(length(allSpatVars) == 0){ return(c()) }else{}
  
  varsWithKnots <- unique(unlist(strsplit(colnamesX , '___')))
  
  varsWithKnots <- gsub('_KNOT0' , '_JUNKFORNOW' , varsWithKnots)
  
  varsWithKnots <- varsWithKnots[which(grepl('_KNOT' , varsWithKnots))]

  nKnotssVec <- NA * integer(length(allSpatVars)) 
  names(nKnotssVec) <- allSpatVars
  for(i in 1:length(allSpatVars)){
    varsWithKnotssThis <- varsWithKnots[which(substr(varsWithKnots , 1 , nchar(allSpatVars[i]) + 5) == paste0(allSpatVars[i] , '_KNOT'))]
    nKnotssVec[i] <- length(unique(varsWithKnotssThis))
  }

  return(nKnotssVec)
}

#####################################################
### fn to get list of the vars (with knots included in names) in each column...
#####################################################
getListVarsWithKnots <- function(colnamesX){
  listVarsWithKnots <- strsplit(colnamesX , '___')
  
  return(listVarsWithKnots)
}

#####################################################
### fn to get list of the vars (without knots included in names) in each column...
#####################################################
getListVarsWithoutKnots <- function(colnamesX){
  listVarsWithKnots <- strsplit(colnamesX , '___')
  
  tmp <- strsplit(unlist(listVarsWithKnots) , '_KNOT')
  tmp <- tmp[unlist(lapply(tmp , length)) == 2]
  if(length(tmp) == 0){
    return(listVarsWithKnots)
  }else{}
  
  allKnotNos <- as.numeric(unique(unlist(lapply(tmp , '[[' , 2))))
  listVarsWithoutKnots <- listVarsWithKnots
  for(ik in allKnotNos){
    listVarsWithoutKnots <- lapply(listVarsWithoutKnots , gsub , pattern = paste0('_KNOT' , ik) , replacement = '')
  }

  return(listVarsWithoutKnots)
}

#####################################################
### fn to perform backward stepwise selection, using Wald tests...
### permissible drops include 1. any knot; 2. any basis column from an interaction; 3. any complete interaction.
### pVec/aicVec in returned list are for removal of each individual basis fn from the final model.
### pHistory/aicHistory track through the process, first val is from full model for aic, for p second value is p of the first dropped term
### function wraps up stepBackXcnsWorker which does the work
### if no condits applied on complexity (maxnSpatVars = Inf , maxnKnotsd = Inf , maxnKnotss = Inf)
### then only one call to the worker
### if some conditions applied, then first call results in a subset of variables/knots that satisty conditions
### worker is then re-called with this subset (perhaps useful knots were removed early on in the first call)
### (the restart is from reduced pool which satisfies just the constraint on number of variables included)
#####################################################
stepBackXcns <- function(Xcns , zData , alpha = 0.15 , iAX = NULL , iAz = NULL , XziAXz = NULL , sumn_CLEidsvik = NULL , optStat4Drop = 1 , allowVariableRemoval = TRUE , 
                         maxnSpatVars = Inf , maxnKnotsd = Inf , maxnKnotss = Inf){

  # optStat4Drop = 1 : aic of models fitted by reml (nested one evaluated in framework of larger one, see Marchant et al., 2009/welham and thompson, 1997)
  # optStat4Drop = 2 : aic of models fitted by ml (nested one evaluated in framework of larger one)
  # optStat4Drop = 3 : p values based on Wald test (with alpha the sig level)
  
  if(!is.null(sumn_CLEidsvik)){
    if(optStat4Drop != 2){ stop('Error in input to stepBackXcns - if you want to use Eidsvik approx, use optStat4Drop=2 to indicate AIC with ML') }else{}
  }else{}

  resBS <- stepBackXcnsWorker(Xcns = Xcns , zData = zData , alpha = alpha , iAX = iAX , iAz = iAz , XziAXz = XziAXz , sumn_CLEidsvik = sumn_CLEidsvik , optStat4Drop = optStat4Drop , allowVariableRemoval = allowVariableRemoval , 
                           maxnSpatVars = maxnSpatVars , maxnKnotsd = maxnKnotsd , maxnKnotss = maxnKnotss)

  print('Done the initial step back.')
  
  ### if some conditions applied, then run again with reduced initial pool...
  ### (just reduced to satisfy constraint on number of variables included)
  if(min(c(maxnSpatVars , maxnKnotsd , maxnKnotss)) < Inf){
    spatVarsInReducedModel <- getSpatVars(colnames(resBS$Xcns))
    listVarsInReducedModel <- getListVarsWithoutKnots(colnames(resBS$Xcns))
    
    listVarsInFullModel <- getListVarsWithoutKnots(colnames(Xcns))
      
    listCheck <- unique(unlist(listVarsInReducedModel))
    listCheck <- setdiff(listCheck , 'const')

    if(length(listCheck) > 1){
      iCombos <- combn(length(listCheck) , 2)
      listCheck1 <- matrix(NA , ncol(iCombos) , 2)
      listCheck1[,1] <- listCheck[iCombos[1,]]
      listCheck1[,2] <- listCheck[iCombos[2,]]
      listCheck1 <- rbind(listCheck1 , listCheck1[,c(2,1),drop=FALSE]) # gets either order.
      
      listCheck1 <- lapply(seq_len(nrow(listCheck1)), function(i) listCheck1[i,]) # convert to list.
      listCheck <- c(as.list(listCheck) , listCheck1)
      
    }else{
      listCheck <- as.list(listCheck)
    }
    listCheck <- c('const' , listCheck)
    
    colsInc4Restart <- which(is.element(listVarsInFullModel , listCheck))

    if(!is.null(XziAXz)){
      XziAXz4Restart <- XziAXz[c(colsInc4Restart,nrow(XziAXz)),c(colsInc4Restart,nrow(XziAXz)),drop=FALSE]
    }else{
      XziAXz4Restart <- NULL
    }
    
    print('Now to try a second step back...')
    
    ### run again with just these cols to start off...
    resBS <- stepBackXcnsWorker(Xcns = Xcns[,colsInc4Restart,drop=FALSE] , zData = zData , alpha = alpha , 
                          iAX = iAX[,colsInc4Restart,drop=FALSE] , iAz = iAz , XziAXz = XziAXz4Restart , sumn_CLEidsvik = sumn_CLEidsvik , optStat4Drop = optStat4Drop , allowVariableRemoval = allowVariableRemoval , 
                          maxnSpatVars = maxnSpatVars , maxnKnotsd = maxnKnotsd , maxnKnotss = maxnKnotss)

  }else{}   
  
  return(resBS)
}

###############################################################################
### this is the function associated with stepBackXcns that does the actual work...
###############################################################################
stepBackXcnsWorker <- function(Xcns , zData , alpha = 0.15 , iAX = NULL , iAz = NULL , XziAXz = NULL , sumn_CLEidsvik = NULL , optStat4Drop = 1 , allowVariableRemoval = TRUE , 
                         maxnSpatVars = Inf , maxnKnotsd = Inf , maxnKnotss = Inf){
  
  if(!is.null(sumn_CLEidsvik)){
    if(optStat4Drop != 2){ stop('Error in input to stepBackXcns - if you want to use Eidsvik approx, use optStat4Drop=2 to indicate AIC with ML') }else{}
  }else{}
  
  # optStat4Drop = 1 : aic of models fitted by reml (nested one evaluated in framework of larger one, see Marchant et al., 2009/welham and thompson, 1997)
  # optStat4Drop = 2 : aic of models fitted by ml (nested one evaluated in framework of larger one)
  # optStat4Drop = 3 : p values based on Wald test (with alpha the sig level)
  
  colsDropped <- c()
  
  p <- ncol(Xcns)
  
  ### all of the spatial variables in the model at the start of the algorithm...
  allSpatVars <- getSpatVars(colnames(Xcns)) 
  
  if(is.null(iAX) & is.null(XziAXz)){ 
    XiAX <- t(Xcns) %*% Xcns
    XiAz <- t(Xcns) %*% zData
    ziAz <- sum(zData ^ 2)
  }else{
    if(is.null(XziAXz)){
      XiAX <- t(Xcns) %*% iAX
      XiAz <- t(Xcns) %*% iAz
      ziAz <- t(zData) %*% iAz 
    }else{
      XiAX <- XziAXz[1:p,1:p,drop=FALSE]
      XiAz <- XziAXz[1:p,p+1,drop=FALSE]
      ziAz <- XziAXz[p+1,p+1]
    }
  }    
  
  if(ncol(Xcns) == 1){ 
    betahat <- solve(XiAX , XiAz)
    if(is.null(sumn_CLEidsvik)){
      sigma2hat <- as.numeric(ziAz - 2 * t(betahat) %*% XiAz + t(betahat) %*% XiAX %*% betahat) / (nrow(Xcns) - ncol(Xcns))
      vbetahat <- sigma2hat * solve(XiAX)
    }else{
      # the real X iA X / n is approx sum(Xj iAj Xj) / sum(nj)      
      sigma2hat <- as.numeric(ziAz - 2 * t(betahat) %*% XiAz + t(betahat) %*% XiAX %*% betahat) / sumn_CLEidsvik
      vbetahat <- sigma2hat * solve(XiAX * nrow(Xcns) / sumn_CLEidsvik)
    }
    return(list('Xcns' = Xcns , 'betahat' = betahat , 'vbetahat' = vbetahat , 'pVec' = NA , 'aicVec' = NA , 'pHistory' = NA , 'aicHistory' = NA , 'colsDropped' = colsDropped)) 
  }else{}
  
  #######################################################
  ### check initial model is solvable...   
  ### if not, then go through cols in order and only include if solvable.
  #######################################################
  betahat <- try(solve(XiAX , XiAz) , silent = TRUE)
  sigma2hat <- try(as.numeric(ziAz - 2 * t(betahat) %*% XiAz + t(betahat) %*% XiAX %*% betahat) / (nrow(Xcns) - ncol(Xcns)) , silent = TRUE)
  
  if((ncol(Xcns) >= nrow(Xcns)) | is.character(betahat) | is.character(sigma2hat) | (is.numeric(sigma2hat) && (sigma2hat < (1E-12 * var(zData))))){
    print('Error - initial matrix gives very small condition number...')
    stop('Try removing columns before passing to the stepBackXcns function (perhaps using reduceXUsingrcond function) to remove colinearity.')
    
    # colsInc <- 1
    # print('Initial model not solvable, so dropping columns based on order of entry so that system is solvable.')
    # print('NOTE - THIS SHOULD NOT HAPPEN ANYMORE BUT IT HAS!')
    # print('Dropped columns are:')
    # 
    # colsIn <- colnames(Xcns)
    # tmp <- reduceXUsingrcond(XIn = Xcns , rqr_p_LT_n = TRUE , rcondLim = rcondLim)
    # Xcns <- tmp$XOut  
    # colsDropped <- tmp$colsDropped # overwriting the blank initialised colsDropped
    # colsInc <- setdiff(colsIn , colsDropped)
    # rm(colsIn)
    # 
    # print(colsDropped)
    # print('')
    # Xcns <- Xcns[,colsInc,drop=FALSE]
    # if(!is.null(iAX)){
    #   iAX <- iAX[colsInc,,drop=FALSE]
    # }else{}
    # XiAX <- XiAX[colsInc,colsInc,drop=FALSE]
    # XiAz <- XiAz[colsInc,,drop=FALSE]
  }else{}    
  
  if(optStat4Drop == 1 | optStat4Drop == 2){
    aicFull <- aic4Drop(colsDrop = c() , XFull = Xcns , XiAXFull = XiAX , XiAzFull = XiAz , ziAz = ziAz , optStat4Drop = optStat4Drop , sumn_CLEidsvik = sumn_CLEidsvik)
    
    pVec <- NA
    pHistory <- NA
    aicHistory <- aicFull
    
  }else if(optStat4Drop == 3){
    ### based on init (legalized) model...  
    betahat <- try(solve(XiAX , XiAz) , silent = TRUE)
    sigma2hat <- as.numeric(ziAz - 2 * t(betahat) %*% XiAz + t(betahat) %*% XiAX %*% betahat) / (nrow(Xcns) - ncol(Xcns))
    vbetahat <- sigma2hat * solve(XiAX)
    
    aicVec <- NA
    pHistory <- NA # starts with NA for consistent vector with aic, in which first element is aic of full model.
    aicHistory <- NA
  }else{}
  
  ### in phase 1, try dropping basis fns of interactions, or complete interactions...
  ### in phase 2, also allow complete removal of knots...
  continueDropping <- TRUE
  allowKnotRemoval <- FALSE
  while(continueDropping){
    
    listPossRedns <- getPossibleReductions(Xcns , allowKnotRemoval = allowKnotRemoval , allowVariableRemoval = allowVariableRemoval)
    
    if(length(listPossRedns) > 0){
      if(optStat4Drop == 1 | optStat4Drop == 2){
        aicFull <- aic4Drop(c() , XFull = Xcns , XiAXFull = XiAX , XiAzFull = XiAz , ziAz = ziAz , optStat4Drop = optStat4Drop , sumn_CLEidsvik = sumn_CLEidsvik)
        aicVec <- unlist(lapply(listPossRedns , aic4Drop , XFull = Xcns , XiAXFull = XiAX , XiAzFull = XiAz , ziAz = ziAz , optStat4Drop = optStat4Drop , sumn_CLEidsvik = sumn_CLEidsvik))
        
        if(min(aicVec) < aicFull){
          colsDrop <- listPossRedns[[which.min(aicVec)]]
          colsDrop_to_satisfy_constraints <- c()
          aicBest <- min(aicVec)
          aicHistory <- c(aicHistory , aicBest)
        }else{
          colsDrop <- c()
          colsDrop_to_satisfy_constraints <- listPossRedns[[which.min(aicVec)]] # this only needs defining if colsDrop empty (ie dropping doesn't improve fit, but needed to get to model satisfying complexity constraints)
        }
      }else if(optStat4Drop == 3){
        betahat <- try(solve(XiAX , XiAz) , silent = TRUE)
        sigma2hat <- as.numeric(ziAz - 2 * t(betahat) %*% XiAz + t(betahat) %*% XiAX %*% betahat) / (nrow(Xcns) - ncol(Xcns))
        vbetahat <- sigma2hat * solve(XiAX)
        
        pVec <- unlist(lapply(listPossRedns , waldTest , betahat = betahat , vbetahat = vbetahat))
        if(max(pVec) > alpha){
          colsDrop <- listPossRedns[[which.max(pVec)]]
          colsDrop_to_satisfy_constraints <- c()
          pHistory <- c(pHistory , max(pVec))
        }else{
          colsDrop <- c()
          colsDrop_to_satisfy_constraints <- listPossRedns[[which.max(pVec)]]
        }
      }else{
        stop('Error - unknown option for optStat4Drop!')
      }
    }else{
      # print('No more possible reductions.')
      colsDrop <- c()
      colsDrop_to_satisfy_constraints <- c()
    }
    
    if(length(colsDrop) > 0){
      colsDropped <- c(colsDropped , colnames(Xcns)[colsDrop])
      
      Xcns <- Xcns[,-colsDrop,drop=FALSE]
      
      XiAX <- XiAX[-colsDrop,-colsDrop,drop=FALSE]
      XiAz <- XiAz[-colsDrop,,drop=FALSE]
      
      if(!is.null(iAX)){ 
        iAX <- iAX[,-colsDrop,drop=FALSE]
      }else{}
    }else{
      
      if(allowKnotRemoval){
        # continueDropping <- FALSE # was just this line in the if(allowKnotRemoval) block.
        
        nSpatVars <- length(getSpatVars(colnames(Xcns)))
        nKnotsd <- getnKnotsd(colnames(Xcns))
        nKnotssVec <- getnKnotss(colnames(Xcns) , allSpatVars = allSpatVars)
        
        if((length(colsDrop_to_satisfy_constraints) > 0) & 
           ((nSpatVars > maxnSpatVars) | (nKnotsd > maxnKnotsd) | (max(nKnotssVec) > maxnKnotss))){
          colsDropped <- c(colsDropped , colnames(Xcns)[colsDrop_to_satisfy_constraints])
          
          Xcns <- Xcns[,-colsDrop_to_satisfy_constraints,drop=FALSE]
          
          XiAX <- XiAX[-colsDrop_to_satisfy_constraints,-colsDrop_to_satisfy_constraints,drop=FALSE]
          XiAz <- XiAz[-colsDrop_to_satisfy_constraints,,drop=FALSE]
          
          if(!is.null(iAX)){ 
            iAX <- iAX[,-colsDrop_to_satisfy_constraints,drop=FALSE]
          }else{}
          
        }else{
          continueDropping <- FALSE
        }
        
      }else{
        ### switch to phase 2, but don't stop yet...
        allowKnotRemoval <- TRUE
      }
    }
  }
  
  if(optStat4Drop == 1 | optStat4Drop == 2){
    aicVec <- unlist(lapply(as.list(seq(nrow(XiAX))) , aic4Drop , XFull = Xcns , XiAXFull = XiAX , XiAzFull = XiAz , ziAz = ziAz , optStat4Drop = optStat4Drop , sumn_CLEidsvik = sumn_CLEidsvik))
    ### also include in aicVec    
    names(aicVec) <- colnames(Xcns)
    aicFinalModel <- aic4Drop(c() , XFull = Xcns , XiAXFull = XiAX , XiAzFull = XiAz , ziAz = ziAz , optStat4Drop = optStat4Drop , sumn_CLEidsvik = sumn_CLEidsvik)
  }else if(optStat4Drop == 3){
    ### p values from final tests for each individual basis fn...
    pVec <- unlist(lapply(as.list(seq(nrow(vbetahat))) , waldTest , betahat = betahat , vbetahat = vbetahat))
    names(pVec) <- colnames(Xcns)
    ### get ml aic...    
    aicFinalModel <- aic4Drop(c() , XFull = Xcns , XiAXFull = XiAX , XiAzFull = XiAz , ziAz = ziAz , optStat4Drop = 2 , sumn_CLEidsvik = sumn_CLEidsvik)
  }else{
    stop('Error - unknown option for optStat4Drop!')
  }  
  
  betahat <- try(solve(XiAX , XiAz) , silent = TRUE)
  sigma2hat <- as.numeric(ziAz - 2 * t(betahat) %*% XiAz + t(betahat) %*% XiAX %*% betahat) / (nrow(Xcns) - ncol(Xcns))
  vbetahat <- sigma2hat * solve(XiAX)
  
  return(list('Xcns' = Xcns , 'betahat' = betahat , 'vbetahat' = vbetahat , 'pVec' = pVec , 'aicVec' = aicVec , 'pHistory' = pHistory , 'aicHistory' = aicHistory , 'aicFinalModel' = aicFinalModel , 'colsDropped' = colsDropped))
}

##########################################################################
### fn to calc term based on aic of models for testing drop...
### iA assumed fixed, reml, within a 'full' model, X (welham/thompson method).
### lndetC = lndetA + n * logsigma2hat
### because A not being changed, lndetA can be ignored here.
### and in w/t framework, lndetXiAX (if A fixed) is same for both models, so not required
### constant also not included.
##########################################################################
aic4Drop <- function(colsDrop , XFull , XiAXFull , XiAzFull , ziAz , optStat4Drop = 1 , sumn_CLEidsvik = NULL){
  if(length(colsDrop) == 0){ 
    X <- XFull
    XiAX <- XiAXFull
    XiAz <- XiAzFull
  }else{
    X <- XFull[,-colsDrop,drop=FALSE]
    XiAX <- XiAXFull[-colsDrop,-colsDrop,drop=FALSE]
    XiAz <- XiAzFull[-colsDrop,,drop=FALSE]
  }
  betahat <- try(solve(XiAX , XiAz) , silent = TRUE)
  
  if(is.character(betahat)){
    return(9E99)
  }else{}
  
  if(optStat4Drop == 1){
    sigma2hat <- as.numeric(ziAz - 2 * t(betahat) %*% XiAz + t(betahat) %*% XiAX %*% betahat) / (nrow(X) - ncol(X))
    nll <- 0.5 * ((nrow(X) - ncol(XFull)) * log(sigma2hat) + as.numeric(ziAz - t(XiAz) %*% betahat) / sigma2hat)
    aic <- 2 * nll + 2 * ncol(X)
  }else if(optStat4Drop == 2){
    if(is.null(sumn_CLEidsvik)){
      sigma2hat <- as.numeric(ziAz - 2 * t(betahat) %*% XiAz + t(betahat) %*% XiAX %*% betahat) / nrow(X)
      nll <- 0.5 * (nrow(X) * log(sigma2hat) + as.numeric(ziAz - t(XiAz) %*% betahat) / sigma2hat)
      aic <- 2 * nll + 2 * ncol(X)
    }else{
# inconsistent with the above is.null(sumn_CLEidsvik) version, but won't ever be compared.  
# in this case, XiAXFull , XiAzFull , ziAz entered as sums over all pairs of blocks
      sigma2hat <- as.numeric(ziAz - 2 * t(betahat) %*% XiAz + t(betahat) %*% XiAX %*% betahat) / sumn_CLEidsvik
      nll <- 0.5 * sumn_CLEidsvik * (log(sigma2hat) + 1)
      aic <- 2 * nll + 2 * ncol(X)
    }
    
  }else{
    stop('Unknown optStat4Drop in aic function (should be 1 or 2)!')
  }
  
  return(aic)
}

##########################################################################
### add scaled versions of covariates to fitting and calibration data...
##########################################################################
addScaledCovs2df_OLD <- function(dfFit = NULL , dfPred = NULL , covNames , scalePars_m = NULL , scalePars_sd = NULL){
  
  for(nm in covNames){
    if(substr(nm , nchar(nm) - 6 , nchar(nm)) == '_SCALED'){
      
      ### get and check the original name, which must be in the df...      
      nmunscaled <- substr(nm , 1 , nchar(nm) - 7)
      if(!is.element(nmunscaled , names(dfFit))){ stop(paste0('Error - ' , nmunscaled , ' not found in dfFit!')) }else{}
      if(!is.null(dfPred)){
        if(nrow(dfPred) > 0){
          if(!is.element(nmunscaled , names(dfPred))){ stop(paste0('Error - ' , nmunscaled , ' not found in dfPred!')) }else{}
        }else{}
      }else{}
      
      mTmp <- mean(dfFit[[nmunscaled]])
      sdTmp <- sd(dfFit[[nmunscaled]])
      dfFit[[nm]] <- (dfFit[[nmunscaled]] - mTmp) / sdTmp
      if(!is.null(dfPred)){
        if(nrow(dfPred) > 0){
          dfPred[[nm]] <- (dfPred[[nmunscaled]] - mTmp) / sdTmp
        }else{}
      }else{}
      
    }else{
      if(!is.element(nm , names(dfFit))){ stop(paste0('Error - ' , nm , ' not found in dfFit!')) }else{}
      if(!is.null(dfPred)){
        if(nrow(dfPred) > 0){
          if(!is.element(nm , names(dfPred))){ stop(paste0('Error - ' , nm , ' not found in dfPred!')) }else{}
        }else{}
      }else{}
    }
  }
  
  ### if dfPred was passed in as df with no rows, update the col names by copying from dfFit...
  if(!is.null(dfPred)){
    if(nrow(dfPred) == 0){
      dfPred <- dfFit[c(),,drop=FALSE]
    }else{}
  }else{}
  
  return(list('dfFit' = dfFit , 'dfPred' = dfPred))
}

##########################################################################
### add scaled versions of covariates to fitting and calibration data...
##########################################################################
addScaledCovs2df <- function(dfFit = NULL , covNames = NULL , scalePars_m = NULL , scalePars_sd = NULL){

  if(is.null(scalePars_m)){
    getScalingParsNow <- TRUE
    scalePars_m <- scalePars_sd <- list()
  }else{
    getScalingParsNow <- FALSE
    if(!is.null(covNames)){ print('No need to pass in covNames to addScaledCovs2df if passing in scalePars_m and scalePars_sd.') }else{}
    covNames <- paste0(names(scalePars_m) , '_SCALED') 
  }
  
  # if(length(covNames) > 0)

  for(nm in covNames){
    if(substr(nm , nchar(nm) - 6 , nchar(nm)) == '_SCALED'){
      
      ### get and check the original name, which must be in the df...      
      nmunscaled <- substr(nm , 1 , nchar(nm) - 7)
      if(!is.element(nmunscaled , names(dfFit))){ stop(paste0('Error - ' , nmunscaled , ' not found in dfFit!')) }else{}

      if(getScalingParsNow){
        mTmp <- mean(dfFit[[nmunscaled]])
        sdTmp <- sd(dfFit[[nmunscaled]])
        scalePars_m[[nmunscaled]] <- mTmp
        scalePars_sd[[nmunscaled]] <- sdTmp
      }else{
        mTmp <- scalePars_m[[nmunscaled]]
        sdTmp <- scalePars_sd[[nmunscaled]]
        if(is.null(mTmp)){ stop(paste0('Error - scaling parameters for ' , nmunscaled , ' not found in scalePars_m!')) }else{}
        if(is.null(sdTmp)){ stop(paste0('Error - scaling parameters for ' , nmunscaled , ' not found in scalePars_sd!')) }else{}
      }      
      dfFit[[nm]] <- (dfFit[[nmunscaled]] - mTmp) / sdTmp

    }else{
      if(!is.element(nm , names(dfFit))){ stop(paste0('Error - ' , nm , ' not found in dfFit!')) }else{}
    }
  }
  
  return(list('dfFit' = dfFit , 'scalePars_m' = scalePars_m , 'scalePars_sd' = scalePars_sd))
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
  ### if q4BdryKnots are in 0-1, they are quantiles of the input covariate values
  ### if q4BdryKnots are <0 or >1, they are used to define bdry knots by extending range of the input covariate values
  ###             eg -0.25 means min(x) - 0.25 * (max(x) - min(x))
  ###             eg 1.25 means max(x) + 0.25 * (max(x) - min(x))
  ### this option / quantile wrapped up in quantile4BdryKnotsFn function
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
          fefdKnotsTmp$bdryKnots <- quantile4BdryKnotsFn(dfFit[[nm]] , q4BdryKnots[inm,])
        }else{
          fefdKnotsTmp$bdryKnots <- quantile4BdryKnotsFn(dfFit[[nm]] , q4BdryKnots)
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

############################################
### function to return boundary knots as either quantiles or extended ranges. 
### if q4BdryKnots are in 0-1, they are quantiles of the input covariate values
### if q4BdryKnots are <0 or >1, they are used to define bdry knots by extending range of the input covariate values
###             eg -0.25 means min(x) - 0.25 * (max(x) - min(x))
###             eg 1.25 means max(x) + 0.25 * (max(x) - min(x)) [ie (1.25 - 1) * the range]
############################################
quantile4BdryKnotsFn <- function(xIn , q4BdryKnots){

  if(q4BdryKnots[1] >= q4BdryKnots[2]){ stop('Error - q4BdryKnots[1] must be smaller than q4BdryKnots[2]!') }else{}
  if(q4BdryKnots[1] >= 1){ stop('Error - q4BdryKnots[1] must be smaller than 1 - that would be lower boundary at the maximum of the data!') }else{}
  if(q4BdryKnots[2] <= 0){ stop('Error - q4BdryKnots[2] must be larger than 0 - that would be upper boundary at the minimum of the data!') }else{}
  
  bdryKnots <- NA * numeric(2)
  if(all(q4BdryKnots >= 0 & q4BdryKnots <= 1)){
    bdryKnots <- quantile(xIn , q4BdryKnots)
  }else if(q4BdryKnots[1] < 0 & q4BdryKnots[2] <= 1){
    maxxIn <- max(xIn) 
    minxIn <- min(xIn)
    bdryKnots[1] <- minxIn + q4BdryKnots[1] * (maxxIn - minxIn)
    bdryKnots[2] <- quantile(xIn , q4BdryKnots[2])
  }else if(q4BdryKnots[1] >= 0 & q4BdryKnots[2] > 1){
    maxxIn <- max(xIn) 
    minxIn <- min(xIn)
    bdryKnots[1] <- quantile(xIn , q4BdryKnots[1])
    bdryKnots[2] <- maxxIn + (q4BdryKnots[2] - 1) * (maxxIn - minxIn)
  }else if(q4BdryKnots[1] < 0 & q4BdryKnots[2] > 1){
    maxxIn <- max(xIn) 
    minxIn <- min(xIn)
    rangexIn <- maxxIn - minxIn
    bdryKnots[1] <- minxIn + q4BdryKnots[1] * rangexIn
    bdryKnots[2] <- maxxIn + (q4BdryKnots[2] - 1) * rangexIn
  }else{
    stop('Error - this should never happen - check what happened in quantile4BdryKnotsFn!')
  }
  
  return(bdryKnots)  
}

##########################################################
### an algorithm to suggest variables that might be useful for explaining within-field variation
### for (i) topsoil and (ii) subsoil?
##########################################################
selectSpatVars4IAK3D <- function(dfFit , responseVar , covariateVars , siteVar = NULL , 
                                 nSpatVarsSelect , dITopSub = data.frame('dU' = c(0 , 1) , 'dL' = c(0.1 , 1.1) , stringsAsFactors = FALSE)){
  
  if(is.null(dfFit$profID)){ stop('In selectSpatVars4IAK3D, dfFit must have profID as a column!') }else{}
  
  if(nSpatVarsSelect >= length(covariateVars)){
    print(paste0('Requested to select ' , nSpatVarsSelect ,  ' variables in function selectSpatVars4IAK3D, but only ' , length(covariateVars) , ' variables passed into function! So they are all selected!'))
    return(covariateVars)
  }else{}
  
  dfFit <- dfFit[which(!is.na(dfFit[[soilPropModelTfmd]])),,drop=FALSE]
  if(nrow(dfFit) == 0){ stop('Error - no valid response data found in selectSpatVars4IAK3D!') }else{}
  
  if(any(is.na(dfFit[,covariateVars,drop=FALSE]))){ stop('In selectSpatVars4IAK3D, dfFit must have no NA values for covariates!') }else{}
  
  ###############################################
  ### harmonise soil data to the top and subsoil depths...
  ###############################################
  dfFit_Harm <- harmonizeMPS(profIDData = as.character(dfFit$profID) , dIData = dfFit[,c('dU' , 'dL'),drop=FALSE] , zData = dfFit[[responseVar]] , dfSpatialData4Copy = dfFit[,c(siteVar , covariateVars),drop=FALSE] , 
                             dIStd = dITopSub , vlow = -Inf , vhigh = Inf , nmSoilProp = responseVar , singlesByRegression = TRUE)$dfHrmnzdData
  
  dfadjR2 <- data.frame('covariate' = covariateVars , 'adjR2_Top' = NA , 'adjR2_Sub' = NA , stringsAsFactors = FALSE)
  for(idepth in 1:2){
    rowsTmp <- which((!is.na(dfFit_Harm[[responseVar]])) & 
                       (round(dfFit_Harm$dU , digits = 2) == round(dITopSub$dU[idepth] , digits = 2)) & 
                       (round(dfFit_Harm$dL , digits = 2) == round(dITopSub$dL[idepth] , digits = 2)))
    for(i in 1:nrow(dfadjR2)){
      covTmp <- dfFit_Harm[[dfadjR2$covariate[i]]][rowsTmp]
      covTmp <- (covTmp - mean(covTmp)) / sd(covTmp)
      Xcns <- makeXnsClampedLUG0(covTmp , bdryKnots = quantile(covTmp , c(0.01 , 0.99)) , intKnots = quantile(covTmp , c(0.25 , 0.5 , 0.75)))
      basisFnsTmp <- Xcns[,-1,drop=FALSE]
      namesbasisFnsTmp <- paste0('bf_' , seq(ncol(basisFnsTmp)))
      dfTmp <- cbind(dfFit_Harm[rowsTmp,responseVar,drop=FALSE] , basisFnsTmp)
      names(dfTmp)[2:ncol(dfTmp)] <- namesbasisFnsTmp
      if(!is.null(siteVar)){
        dfTmp$siteID <- dfFit_Harm[[siteVar]][rowsTmp]
        formulaTmp <- as.formula(paste0(responseVar , ' ~ siteID + ' , paste(namesbasisFnsTmp , collapse = ' + ')))
      }else{
        formulaTmp <- as.formula(paste0(responseVar , ' ~ ' , paste(namesbasisFnsTmp , collapse = ' + ')))
      }
      
      lmTmp <- lm(formulaTmp , data = dfTmp)
      sumlmTmp <- summary(lmTmp)
      if(idepth == 1){
        dfadjR2$adjR2_Top[i] <- sumlmTmp$adj.r.squared
      }else{
        dfadjR2$adjR2_Sub[i] <- sumlmTmp$adj.r.squared
      }
    }
  }
  
  ################################
  # select variables:
  # one with highest overall adjR2, the other with the highest adjR2 for the other depth...
  # iterate until we have nSpatVarsSelect (if odd, will have one too many, so take first nSpatVarsSelect)
  ################################
  varsInc <- c()
  for(i in 1:ceiling(nSpatVarsSelect/2)){
    maxadjR2 <- max(dfadjR2[,c('adjR2_Top' , 'adjR2_Sub'),drop=FALSE])
    if(any(dfadjR2$adjR2_Top == maxadjR2)){
      ivarsInc <- which.max(dfadjR2$adjR2_Top)
      adjR2Tmp <- dfadjR2$adjR2_Sub
      adjR2Tmp[ivarsInc] <- -Inf
      dfadjR2$adjR2_Top[ivarsInc] <- -Inf
      dfadjR2$adjR2_Sub[ivarsInc] <- -Inf
      ivarsInc <- c(ivarsInc , which.max(adjR2Tmp))
      varsInc <- c(varsInc , dfadjR2$covariate[ivarsInc])
    }else{
      ivarsInc <- which.max(dfadjR2$adjR2_Sub)
      adjR2Tmp <- dfadjR2$adjR2_Top
      adjR2Tmp[ivarsInc] <- -Inf
      dfadjR2$adjR2_Top[ivarsInc] <- -Inf
      dfadjR2$adjR2_Sub[ivarsInc] <- -Inf
      ivarsInc <- c(ivarsInc , which.max(adjR2Tmp))
      varsInc <- c(varsInc , dfadjR2$covariate[ivarsInc])
    }
  }
  varsInc <- varsInc[1:nSpatVarsSelect]
  
  return(varsInc)
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



