gam2X <- function(gamModel , dataFit = NULL , colsRmv = NULL , allKnotsd = NULL){

######################################################################################
### gamModel is a fitted gam (mgcv) model.
###                to include depth in the gamModel, it should be called dIMidPts in the covariates and the model
###                when iak uses it, basis functions that include depth will be averaged over depth intervals
###                (which is a bit different to evaluating basis functions at dIMidPts)
### dataFit is either the covariate data that were used to fit the model 
###                or the covariate data for the prediction locations
### colsRmv is a vector with the indices of the gamModel design matrix to be removed
###                eg, if gamModel has te(x,y,dIMidPts) terms, then remove these 
###                because iak code models this variation with prod-sum covariance model
### allKnotsd are the knots (both boundaries and any interior knots) for a spline function of depth only
###                don't include this if gamModel already has a smooth of depth only, ie s(dIMidPts)
###                or if cov model for iak includes the d sum term
######################################################################################
  
### if gam already has colsRmv attached, check agrees with the argument (if given)
  if(is.element('colsRmv' , names(gamModel))){
    if(!is.null(colsRmv)){
      if(!identical(colsRmv , gamModel$colsRmv)){
        stop('Error - colsRmv entered as argument to gam2X is not the same as the version attached to the gam!')
      }else{}
    }else{}
  }else{
### attach the given colsRmv now...
    if(is.null(colsRmv)){
      gamModel['colsRmv'] <- list(NULL)
    }else{
      gamModel$colsRmv <- colsRmv
    }
  }

### if gam already has allKnotsd attached, check agrees with the argument (if given)
  if(is.element('allKnotsd' , names(gamModel))){
    if(!is.null(allKnotsd)){
      if(!identical(allKnotsd , gamModel$allKnotsd)){
        stop('Error - allKnotsd entered as argument to gam2X is not the same as the version attached to the gam!')
      }else{}
    }else{}
  }else{
### attach the given allKnotsd now...
    if(is.null(allKnotsd)){
      gamModel['allKnotsd'] <- list(NULL)
    }else{
      gamModel$allKnotsd <- allKnotsd
    }
  }

### make X for the given gam model...
  if(is.null(dataFit)){
    X <- model.matrix(gamModel)
  }else{
    X <- model.matrix(gamModel , newdata = dataFit)
  }

### remove colsRmv from the X matrix...
  if(!is.null(gamModel$colsRmv)){

#########################################################################
### will have to define colsRmv yourself and attach it to the fitted gam: 
### if there is a te(x,y,dIMidPts) term in the gam, remove this 
### as the IAK method models this as spatially/vertically correlated with prod-sum covariance model
### if the iak code is going to be called with allKnotsd given, then also remove any smooths of dIMidPts
#########################################################################
  
    X <- X[,-gamModel$colsRmv,drop=FALSE]

  }else{}
  
###########################################################    
### if spline of depth included, then cbind this...    
### done in the same way as for the cubist2X code.
###########################################################    
  if(length(gamModel$allKnotsd) > 0){
      
      intKnots <- gamModel$allKnotsd[-1]
      intKnots <- intKnots[-length(intKnots)]
      bdryKnots <- c(gamModel$allKnotsd[1] , gamModel$allKnotsd[length(gamModel$allKnotsd)])

      XSpline <- bs(dataFit$dIMidPts , knots = intKnots , degree = 3 , intercept = F , Boundary.knots = bdryKnots)
      colnames(XSpline) <- paste0('dSpline.' , seq(ncol(XSpline)))
      X <- cbind(X , XSpline)

  }else{}
    
  gamModel$namesX <- colnames(X)
  
  return(list('X' = X , 'gamModel' = gamModel))

}


