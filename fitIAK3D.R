fitIAK3D <- function(x , dI , z , covs , modelX , modelx = 'matern' , nud = 0.5 , allKnotsd = c() , 
                sdfdType_cd1 = 0 , sdfdType_cxd0 = 0 , sdfdType_cxd1 = 0 , cmeOpt = 0 , prodSum = TRUE , 
                lnTfmdData = FALSE , useReml = TRUE , compLikMats = list('compLikOptn' = 0) , namePlot = NA , lmmFit = list() , rqrBTfmdPreds = TRUE , parsInit = NULL){

########################################################
### if x or dI were dataframes, convert to matrices here.
### and make sure all are numeric...
########################################################
    if(!is.matrix(x)){
        x <- as.matrix(x)
    }else{}
    if(!is.matrix(dI)){
        dI <- as.matrix(dI)
    }else{}

    xCopy <- matrix(NA , nrow(x) , ncol(x))
    for (i in 1:ncol(x)){ xCopy[,i] <- as.numeric(x[,i]) }
    x <- xCopy
    remove(xCopy)
    
    dICopy <- matrix(NA , nrow(dI) , ncol(dI))
    for (i in 1:ncol(dI)){ dICopy[,i] <- as.numeric(dI[,i]) }
    dI <- dICopy
    remove(dICopy)

#########################################################
### error check for duplicated x,dI data, if no measurement error being included...
#########################################################
    if(cmeOpt == 0){
      iDuplicated <- which(duplicated(cbind(x , dI)))
      if(length(iDuplicated) > 0){ stop(paste0('Error - the data at positions ' , iDuplicated , ' are duplicated!')) }else{} 
    }else{}
    
#########################################################
### error check for any NAs; removing them if found...
#########################################################
    iNA <- which(is.na(z) | rowSums(is.na(covs)) > 0)
    if(length(iNA) > 0){ 
      print('Attention! Some NAs found in data ; removing them.')
      x <- x[-iNA,drop = FALSE]
      dI <- dI[-iNA,drop = FALSE]
      z <- z[iNA]
    }else{}

#########################################################
### round depths to nearest cm (assuming depths were inputted in m)
#########################################################
    dI <- round(dI , digits = 2)
    
################################################
### if compLikMats is not null, fit assuming additive on input scale (ie put lnTfmdData = FALSE)
### ie approx involves approx by arithmetic averaging on log-transformed data + comp lik.
### but allow rqrBTfmdPreds to remain as inputted
################################################
    if(compLikMats$compLikOptn > 0){
      if (lnTfmdData){
        lnTfmdData <- FALSE
        print('Attention - averaging assumed to occur on the transformed scale.')
      }else{}
    }else{}

#######################################
### for Normal data, fit beta and cxd[d=0] automatically
#######################################
### for lnNormal data, fit beta by newton-raphson
#######################################
    if(lnTfmdData){
### in this case, only allowed reml if depth not in fixed effects, 
### so assume has to be ml...
        if(useReml){
            print('Attention! When underlying point-support variable assumed lognormal, cannot use REML because beta parameters contribute non-linearly to likelihood function!')
            print('Therefore switching option to use ML.')
            useReml <- F
        }else{}
    }else{}

###################################################
### if using Eidsvik, make sure ML, not REML...
###################################################
    if(compLikMats$compLikOptn == 2 || compLikMats$compLikOptn == 3){
      if(useReml){
        print('Attention - do not useReml with Eidsvik CL - appropriate for ML only! Changing option here.')
        useReml <- FALSE
      }else{}
    }else{}
    
################################################
### if not prodSum, will be product covariance function
### in which case, put sdfdType_cd1 <- -9 
################################################
    if(!prodSum){
        if (modelx != 'matern'){ stop('Error - for product covariance, specify matern for modelx!') }else{}
        sdfdType_cd1 <- -9 
    }else{}

######################################################
### set up fixed-effect design matrices...
######################################################
    tmp <- makeXvX(covData = covs , dI = dI , modelX = modelX , allKnotsd = allKnotsd , nDiscPts = 1000 , lnTfmdData = lnTfmdData)
    X <- tmp$X
    vXU <- tmp$vXU
    iU <- tmp$iU
    namesX <- tmp$namesX
    XLims <- tmp$XLims

######################################################
### setup random-effects design matrices...
######################################################
    if (compLikMats$compLikOptn == 0){
      setupMats <- setupIAK3D(x , dI , nDscPts = 0)
    }else{
      setupMats <- list()
### order is all adj subset pairs, then all individual subsets (which can be used to get non-adj subset pairs)
      for (i in 1:nrow(compLikMats$subsetPairsAdj)){
        iThis <- c(compLikMats$listBlocks[[compLikMats$subsetPairsAdj[i,1]]]$i , compLikMats$listBlocks[[compLikMats$subsetPairsAdj[i,2]]]$i)
        setupMats[[i]] <- setupIAK3D(x[iThis,,drop = FALSE] , dI[iThis,,drop = FALSE] , nDscPts = 0)
      }
### now all individual subsets...      
      for (i in 1:length(compLikMats$listBlocks)){
        iThis <- compLikMats$listBlocks[[i]]$i
        setupMats[[nrow(compLikMats$subsetPairsAdj)+i]] <- setupIAK3D(x[iThis,,drop = FALSE] , dI[iThis,,drop = FALSE] , nDscPts = 0)
      }
    }
      
######################################################
### set sensible initial parameters and evaluate nll...
######################################################
    verboseOptim <<- T

    if(is.null(lmmFit$parBnds)){
        parBnds <- setParBndsIAK3D(x = x , dI = dI , setupMats = setupMats , compLikMats = compLikMats)
    }else{
        parBnds <- lmmFit$parBnds
    }

    if(is.null(lmmFit$pars)){
    
      if(is.null(parsInit)){
        tmpInit <- setInitsIAK3D(x = x , dI = dI , z = z , X = X , vXU = vXU , iU = iU , modelx = modelx , nud = nud , 
            sdfdType_cd1 = sdfdType_cd1 , sdfdType_cxd0 = sdfdType_cxd0 , sdfdType_cxd1 = sdfdType_cxd1 , 
            cmeOpt = cmeOpt , prodSum = prodSum , setupMats = setupMats , parBnds = parBnds , lnTfmdData = lnTfmdData , lmmFit = lmmFit , compLikMats = compLikMats)
        parsInit <- tmpInit$pars

### for exact lik, run through setInits function again, with the C from above used to give better initial beta (then better pars)...
        updateInits <- FALSE
        if(updateInits & compLikMats$compLikOptn == 0){
          parsBTfmd <- readPars(pars = parsInit , modelx = modelx , nud = nud , 
                sdfdType_cd1 = sdfdType_cd1 , sdfdType_cxd0 = sdfdType_cxd0 , sdfdType_cxd1 = sdfdType_cxd1 , 
                cmeOpt = cmeOpt , prodSum = prodSum , parBnds = parBnds , lnTfmdData = lnTfmdData)

          tmp <- setCIAK3D(parsBTfmd = parsBTfmd , modelx = modelx , 
                sdfdType_cd1 = sdfdType_cd1 , sdfdType_cxd0 = sdfdType_cxd0 , sdfdType_cxd1 = sdfdType_cxd1 , 
                cmeOpt = cmeOpt , setupMats = setupMats)
          initC <- tmp$C

          tmpInit <- setInitsIAK3D(x = x , dI = dI , z = z , X = X , vXU = vXU , iU = iU , modelx = modelx , nud = nud , 
            sdfdType_cd1 = sdfdType_cd1 , sdfdType_cxd0 = sdfdType_cxd0 , sdfdType_cxd1 = sdfdType_cxd1 , 
            cmeOpt = cmeOpt , prodSum = prodSum , setupMats = setupMats , parBnds = parBnds , lnTfmdData = lnTfmdData , lmmFit = lmmFit , compLikMats = compLikMats , initC = initC)
          parsInit <- tmpInit$pars
        }else{}
        
      }else{}
        
      verboseOptim <<- F
    }else{
      parsInit <- lmmFit$pars
    }

    if(is.null(lmmFit$fitRange)){  
      fitRange <- setFitRangeIAK3D(pars = parsInit , modelx = modelx , 
            sdfdType_cd1 = sdfdType_cd1 , sdfdType_cxd0 = sdfdType_cxd0 , sdfdType_cxd1 = sdfdType_cxd1 , 
            cmeOpt = cmeOpt , prodSum = prodSum)            

      fitRangeStat <- setFitRangeIAK3D(pars = parsInit , modelx = modelx , 
            sdfdType_cd1 = 0 , sdfdType_cxd0 = 0 , sdfdType_cxd1 = 0 , 
            cmeOpt = cmeOpt , prodSum = prodSum)            
    }else{
      fitRange <- lmmFit$fitRange
    }

######################################################
### optimize parameters...
######################################################
    if(is.null(lmmFit$pars)){
### print out init pars and nll...
      print('At initial parameters...')
      verboseOptim <<- T
      if (compLikMats$compLikOptn == 0){
        nllInit <- nllIAK3D(pars = parsInit , z = z , X = X , vXU = vXU , iU = iU , modelx = modelx , nud = nud , 
            sdfdType_cd1 = sdfdType_cd1 , sdfdType_cxd0 = sdfdType_cxd0 , sdfdType_cxd1 = sdfdType_cxd1 , 
            cmeOpt = cmeOpt , prodSum = prodSum , setupMats = setupMats , parBnds = parBnds , useReml = useReml , lnTfmdData = lnTfmdData , rtnAll = T)	
      }else{
        nllInit <- nllIAK3D_CL(pars = parsInit , z = z , X = X , modelx = modelx , nud = nud , 
            sdfdType_cd1 = sdfdType_cd1 , sdfdType_cxd0 = sdfdType_cxd0 , sdfdType_cxd1 = sdfdType_cxd1 , 
            cmeOpt = cmeOpt , prodSum = prodSum , setupMats = setupMats , parBnds = parBnds , useReml = useReml , compLikMats = compLikMats , rtnAll = T)	
      }
      verboseOptim <<- F

### and fit...
      print('Now fitting...')
      
      if (compLikMats$compLikOptn == 0){
        parsFit <- optimIt(par = parsInit , fn = nllIAK3D , fitRange = fitRange , 
            z = z , X = X , vXU = vXU , iU = iU , modelx = modelx , nud = nud , 
            sdfdType_cd1 = sdfdType_cd1 , sdfdType_cxd0 = sdfdType_cxd0 , sdfdType_cxd1 = sdfdType_cxd1 , 
            cmeOpt = cmeOpt , prodSum = prodSum , setupMats = setupMats , parBnds = parBnds , useReml = useReml , lnTfmdData = lnTfmdData)
      }else{
        parsFit <- optimIt(par = parsInit , fn = nllIAK3D_CL , fitRange = fitRange , 
            z = z , X = X , modelx = modelx , nud = nud , 
            sdfdType_cd1 = sdfdType_cd1 , sdfdType_cxd0 = sdfdType_cxd0 , sdfdType_cxd1 = sdfdType_cxd1 , 
            cmeOpt = cmeOpt , prodSum = prodSum , setupMats = setupMats , parBnds = parBnds , useReml = useReml , compLikMats = compLikMats)
      }
      
### note that parameters have now been fitted...
      print('Parameters fitted...')

    }else{
### or to just take the values in inits, which were the transformed values taken from the inputted lmmFit...
      parsFit <- list('par' = parsInit)
    }

    verboseOptim <<- T
    if (compLikMats$compLikOptn == 0){
      lmmFit <- nllIAK3D(pars = parsFit$par , z = z , X = X , vXU = vXU , iU = iU , modelx = modelx , nud = nud , 
                sdfdType_cd1 = sdfdType_cd1 , sdfdType_cxd0 = sdfdType_cxd0 , sdfdType_cxd1 = sdfdType_cxd1 , 
                cmeOpt = cmeOpt , prodSum = prodSum , setupMats = setupMats , parBnds = parBnds , useReml = useReml , lnTfmdData = lnTfmdData , rtnAll = T)	
    }else{
      lmmFit <- nllIAK3D_CL(pars = parsFit$par , z = z , X = X , modelx = modelx , nud = nud , 
                sdfdType_cd1 = sdfdType_cd1 , sdfdType_cxd0 = sdfdType_cxd0 , sdfdType_cxd1 = sdfdType_cxd1 , 
                cmeOpt = cmeOpt , prodSum = prodSum , setupMats = setupMats , parBnds = parBnds , useReml = useReml , compLikMats = compLikMats , rtnAll = T)	
    }
    verboseOptim <<- F

############################################
### also attach everything that was used to call the function to lmmFit...
############################################
    lmmFit <- attachSettings2lmmFit(lmmFit , x , dI , z , covs , modelX , namesX , modelx , nud , allKnotsd , XLims , 
               sdfdType_cd1 , sdfdType_cxd0 , sdfdType_cxd1 , cmeOpt , prodSum , setupMats , lnTfmdData , useReml , parBnds , fitRange , compLikMats)

    if (!is.na(namePlot)){
##############################################
### a quick prediction routine...
### predict through the profile for all data locations and for one distant location (with average of covs)
##############################################
######################################################
### set up fixed-effect design matrices...
######################################################
#      dIPred <- cbind(seq(0 , 1.99 , 0.01) , seq(0.01 , 2 , 0.01))
#      dIPred <- cbind(seq(0 , 1.98 , 0.02) , seq(0.02 , 2 , 0.02))
#      dIPred <- cbind(seq(0 , 1.9 , 0.1) , seq(0.1 , 2 , 0.1))
#      dIPred <- cbind(seq(0 , 14.98 , 0.02) , seq(0.02 , 15 , 0.02))
      dIPred <- cbind(seq(0 , parBnds$maxd-0.02 , 0.02) , seq(0.02 , parBnds$maxd , 0.02))

      dIPred <- round(dIPred , digits = 2)
      ndIPred <- nrow(dIPred)

      iTmp <- which(!duplicated(x))
      
      maxnProfPlot <- 200
      if(length(iTmp) > maxnProfPlot){
### take a subsample of 200 profiles to plot...      
### so that repeatable, don't make it random...
        iTmp2 <- round(seq(1 , length(iTmp) , length(iTmp) / maxnProfPlot))
        iTmp <- iTmp[iTmp2]
        profNamesPlot <- paste0('Profile ' , iTmp2)
      }else{
        profNamesPlot <- paste0('Profile ' , seq(length(iTmp)))
      }
      
      covsPred <- covs[iTmp,,drop=FALSE]
      xPred <- x[iTmp,,drop=FALSE]

### add a distant location to the prediction locations...
      xPredDistant <- matrix(c(9E99 , 9E99) , 1 , 2)
#      covsPredDistant <- colMeansAndModes(covsPred)
#      covsPredDistant <- colMeansAndModes(covs)


### initialise with first row...      
      covsPredDistant <- covs[1,,drop=FALSE]
### don't include the depths in finding medoid...
      idCovs <- which(names(covs) == 'dIMidPts')
      if(length(idCovs) > 0){ covsPredDistant[1,idCovs] <- NA }else{}
      if(length(idCovs) < ncol(covs)){ covsPredDistant[1,-idCovs] <- medoid(covs[,-idCovs,drop=FALSE]) }else{}
            
      nxPred <- dim(xPred)[[1]]
        
### call the predict function...
      tmp <- predictIAK3D(xMap = rbind(xPred , xPredDistant) , dIMap = dIPred , covsMap = rbind(covsPred , covsPredDistant) , lmmFit = lmmFit , rqrBTfmdPreds = rqrBTfmdPreds , constrainX4Pred = constrainX4Pred)

      zPred <- tmp$zMap[,1:nxPred]
      vPred <- tmp$vMap[,1:nxPred]
      pi90UPred <- tmp$pi90UMap[,1:nxPred]
      pi90LPred <- tmp$pi90LMap[,1:nxPred]

      zPredDistant <- tmp$zMap[,nxPred+1]
      vPredDistant <- tmp$vMap[,nxPred+1]

### to back transform data for mapping...
      if(lnTfmdData & rqrBTfmdPreds){
###(were inputted and used all the way through on the log scale, 
### but for plot put back to real scale)        
        z <- exp(z)    
      }else{}

#################################################    
### make a pdf with the distant profile prediction (page 1) and all data profiles (6 per page thereafter)...
### note these are not validation predictions, they are predicted at the data profiles given the data for the same profiles
### also note, they are a bi-product of the method (ie you can use the method to predict at profiles where we have data), 
### not to be confused with the spline-then-krige type approach where similar plots may be produced in the first step of analysis
#################################################    
      tmp <- plotProfilesIAK3D(namePlot = namePlot , x = x , dI = dI , z = z , 
                xPred = xPred , dIPred = dIPred , zPred = zPred , pi90LPred = pi90LPred , pi90UPred = pi90UPred , zPredDistant = zPredDistant , profNames = profNamesPlot)

    }else{
      zPredDistant <- vPredDistant <- zPred <- vPred <- pi90LPred <- pi90UPred <- xPred <- dIPred <- NA
    }

#####################################################
### and return...
#####################################################
    return(list('parsFit' = parsFit , 'lmmFit' = lmmFit , 
                'zPredDistant' = zPredDistant , 'vPredDistant' = vPredDistant , 
                'zPred' = zPred , 'vPred' = vPred , 'pi90LPred' = pi90LPred , 'pi90UPred' = pi90UPred , 
                'xPred' = xPred , 'dIPred' = dIPred))
}

###################################################################
### the nll function that is to be minimized...
###################################################################
nllIAK3D <- function(pars , z , X , vXU , iU , modelx , nud ,  
                     sdfdType_cd1 , sdfdType_cxd0 , sdfdType_cxd1 , cmeOpt , prodSum , setupMats , parBnds , useReml , lnTfmdData , rtnAll = F , forCompLik = FALSE){

### default values for return...
    badnll <- 9E99
    cxdhat <- NA
    WiAW <- NA
    lndetA <- NA

    if(printnllTime){
      ptm <- proc.time()
    }else{}

    if(exists("parsTrace4Optim") & (!forCompLik)){
    
      if(is.null(ncol(parsTrace4Optim)) || (ncol(parsTrace4Optim) == length(pars))){

        if(length(parsTrace4Optim) > 0){
          parsTmp <- rbind(parsTrace4Optim , pars)
          nllTmp <- c(nllTrace4Optim , NA) # for this trace, use NA for failed evaluation.
        }else{
          parsTmp <- matrix(pars , nrow = 1)
          nllTmp <- NA # for this trace, use NA for failed evaluation.
        }
        assign("parsTrace4Optim" , parsTmp , envir = .GlobalEnv)
        assign("nllTrace4Optim" , nllTmp , envir = .GlobalEnv)
      }else{}
      
    }else{}
    
    parsBTfmd <- readPars(pars = pars , modelx = modelx , nud = nud , 
                sdfdType_cd1 = sdfdType_cd1 , sdfdType_cxd0 = sdfdType_cxd0 , sdfdType_cxd1 = sdfdType_cxd1 , 
                cmeOpt = cmeOpt , prodSum = prodSum , parBnds = parBnds , lnTfmdData = lnTfmdData)

    tmp <- setCIAK3D(parsBTfmd = parsBTfmd , modelx = modelx , 
                sdfdType_cd1 = sdfdType_cd1 , sdfdType_cxd0 = sdfdType_cxd0 , sdfdType_cxd1 = sdfdType_cxd1 , 
                cmeOpt = cmeOpt , setupMats = setupMats)

    A <- tmp$C
    sigma2Vec <- tmp$sigma2Vec

    parsOK <- T
    if(max(is.na(A)) == 1){ parsOK <- F }else{}
    if(is.infinite(parsBTfmd$cx0) | is.infinite(parsBTfmd$cx1) | is.infinite(parsBTfmd$cd1) | is.infinite(parsBTfmd$cxd0) | is.infinite(parsBTfmd$cxd1)){ parsOK <- F }else{}
   
    if(parsOK){
        n <- length(z)
        X <- as.matrix(X , nrow = n)
        p <- dim(X)[[2]]
        W <- cbind(X,z)

### just to make sure numerical errors have not made it non-symmetric...  		
        A <- 0.5 * (A + t(A))

        if(lnTfmdData){
            tmp <- nrUpdatesIAK3DlnN(z = z , X = X , vXU = vXU , iU = iU , C = A , sigma2Vec = sigma2Vec)
            nll <- tmp$nll
            betahat <- tmp$betahat
            if(rtnAll){ 
#                vbetahat <- solve(tmp$fim) 
### or, alternatively, 
                vXUTmp <- matrix(0 , n * p , p)
                vXUTmp[kronecker(seq(0 , (n - 1) * p , p) , matrix(1 , length(iU) , 1)) + rep(iU , n) , iU] <- vXU
                tmp <- gradHessIAK3DlnN(z = z , X = X , vXU = vXUTmp , iU = seq(p) , betaU = betahat , C = A , lndetC = 0 , TK  = tmp$iC , iXKiCXKXKiC = c() , sigma2Vec = sigma2Vec)

                vbetahat <- solve(tmp$fim) 

            }else{ 
                vbetahat <- NA 
            } 
            cxdhat <- 1
            C <- A
        }else{
          tmp <- lndetANDinvCb(A , W)
          
          if(is.character(tmp$cholC)){
            print('Inversion of A failed!')
            print(paste0('nll = ' , badnll , '; nmPars = ' , paste(round(pars , digits = 3) , collapse = ', ')))
            nll <- badnll
            cxdhat <- NA
            if(forCompLik){
              return(list('pars' = pars , 'parsBTfmd' = parsBTfmd , 'sigma2Vec' = sigma2Vec , 'WiAW' = matrix(NA , p+1 , p+1) , 'lndetA' = NA))
            }else{}
            
            if(rtnAll){ 
### -999 * A returned because it is not the covariance matrix yet, as cxhat was not computed. Still allows A to be retrieved if required though. 
        	    return(list('nll' = nll , 'pars' = pars , 'parsBTfmd' = parsBTfmd , 'betahat' = NA , 'vbetahat' = NA , 'cxdhat' = NA , 
                            'C' = -999 * A , 'sigma2Vec' = sigma2Vec , 'X' = X , 'vXU' = vXU , 'iU' = iU))
            }else{
        	    return(nll)
            }
          }else{}
          
          lndetA <- tmp$lndetC
          iAW <- tmp$invCb

          WiAW <- t(W) %*% iAW
          WiAW  <- 0.5 * (WiAW + t(WiAW))

          if(forCompLik){
            return(list('pars' = pars , 'parsBTfmd' = parsBTfmd , 'sigma2Vec' = sigma2Vec , 'WiAW' = WiAW , 'lndetA' = lndetA , 'A' = A , 'iAW' = iAW))
          }else{}

          XiAX <- WiAW[1:p , 1:p , drop = FALSE]
          XiAz <- WiAW[1:p , p+1 , drop = FALSE]
          ziAz <- as.numeric(WiAW[p+1 , p+1])

          tmp <- lndetANDinvCb(XiAX , XiAz)
          if(is.character(tmp$cholC)){
              print('Inversion of t(X) iA X failed!')
              print(paste0('nll = ' , badnll , '; nmPars = ' , paste(round(pars , digits = 3) , collapse = ', ')))
              nll <- badnll
              if(rtnAll){ 
### -999 * A returned because it is not the covariance matrix yet, as cxhat was not computed. Still allows A to be retrieved if required though. 
        	    return(list('nll' = nll , 'pars' = pars , 'parsBTfmd' = parsBTfmd , 'betahat' = NA , 'vbetahat' = NA , 'cxdhat' = NA , 
                            'C' = -999 * A , 'sigma2Vec' = sigma2Vec , 'X' = X , 'vXU' = vXU , 'iU' = iU))
              }else{
        	    return(nll)
              }
          }else{}
            
          lndetXiAX <- tmp$lndetC
          betahat <- tmp$invCb
          resiAres <- as.numeric(ziAz - 2 * t(betahat) %*% XiAz + t(betahat) %*% XiAX %*% betahat)

          if(useReml){
              cxdhat <- resiAres / (n - p)
          }else{
              cxdhat <- resiAres / n
          }
### update the values in parsBTfmd so that it makes sense when returned to user.          
          parsBTfmd$cx0 <- cxdhat * parsBTfmd$cx0
          parsBTfmd$cx1 <- cxdhat * parsBTfmd$cx1
          parsBTfmd$cd1 <- cxdhat * parsBTfmd$cd1
          parsBTfmd$cxd0 <- cxdhat * parsBTfmd$cxd0
          parsBTfmd$cxd1 <- cxdhat * parsBTfmd$cxd1
          if(cmeOpt == 1){ parsBTfmd$cme <- cxdhat * parsBTfmd$cme }else{}
                
          lndetC <- lndetA + n * log(cxdhat)
          lndetXiCX <- lndetXiAX - p * log(cxdhat)
          resiCres <- resiAres / cxdhat

          if(rtnAll){ 
              vbetahat <- cxdhat * solve(XiAX) 
              C <- cxdhat * A
              sigma2Vec <- cxdhat * sigma2Vec 
          }else{}
          
          if(useReml){
              nll <- 0.5 * ((n - p) * log(2 * pi) + lndetC + lndetXiCX + resiCres) 
          }else{
              nll <- 0.5 * (n * log(2 * pi) + lndetC + resiCres) 
          }
          nll <- as.numeric(nll)

      }
      
    }else{
        nll <- badnll
        if(forCompLik){
          return(list('pars' = pars , 'parsBTfmd' = parsBTfmd , 'sigma2Vec' = sigma2Vec , 'WiAW' = matrix(NA , p+1 , p+1) , 'lndetA' = NA))
        }else{}
    }

    if(is.na(nll) | is.nan(nll) | is.infinite(nll)){ nll <- badnll }else{}

    if(!lnTfmdData){
        txtOut <- paste0('nll = ' , round(nll , digits = 3) , 
                     '; nmPars = ' , paste(round(pars , digits = 3) , collapse = ', ') , 
                     ', cxdhat = ' , round(cxdhat , digits = 4))
    }else{
        txtOut <- paste0('nll = ' , round(nll , digits = 3) , 
                     '; nmPars = ' , paste(round(pars , digits = 3) , collapse = ', '))
    }

    if(verboseOptim){    
      print(txtOut)    
    }else{}

    if(printnllTime){
      print('time for nll evaluation was:')
      print(proc.time() - ptm)
    }else{}
    
    if(exists("parsTrace4Optim") & (length(nllTrace4Optim) > 0) & (!forCompLik)){
      nllTmp <- nllTrace4Optim # update the final value, now fn has been successful.
      if(nll != badnll){
        nllTmp[length(nllTmp)] <- nll
        assign("nllTrace4Optim" , nllTmp , envir = .GlobalEnv)
      }else{}
    }else{}

    if(rtnAll){ 
### added extra info for quicker prediction, 28/2/19...
        if(!lnTfmdData){
          muhat <- X %*% betahat 
        }else{
          muhat <- setmuIAK3D(X = X , vXU = vXU , iU = iU , beta = betahat , diagC = diag(C) , sigma2Vec = sigma2Vec) 
        }
        tmp <- lndetANDinvCb(C , cbind(X , z - muhat))
        
    	return(list('nll' = nll , 'pars' = pars , 'parsBTfmd' = parsBTfmd , 'betahat' = betahat , 'vbetahat' = vbetahat , 'cxdhat' = cxdhat , 
                  'C' = C , 'sigma2Vec' = sigma2Vec , 'X' = X , 'vXU' = vXU , 'iU' = iU , 'iCX' = tmp$invCb[,1:p,drop=FALSE] , 'iCz_muhat' = tmp$invCb[,p+1,drop=FALSE] , 'iC' = chol2inv(tmp$cholC)))
                  
    }else{
    	return(nll)
    }
}	

setCIAK3D <- function(parsBTfmd , modelx , 
                sdfdType_cd1 , sdfdType_cxd0 , sdfdType_cxd1 , cmeOpt , setupMats){

### to be included in args: , spliti = NA
### if spliti is given as an integer, rqr cov mtx between first spliti points and the remaining points of z.
### else full cov mtx for all of z rqd.

### defined for the unique locations...
    phix0 <- 1 * (setupMats$Dx == 0)
    if (modelx == 'matern'){
        phix1 <- maternCov(setupMats$Dx , c(1 , parsBTfmd$ax , parsBTfmd$nux))
    }else if(modelx == 'nugget'){
        phix1 <- matrix(0 , dim(setupMats$Dx)[[1]] , dim(setupMats$Dx)[[1]])
    }else{
        stop('Not ready!')
    }

### if any components are stationary, function only neds to be run once...
    iStat <- which(c(sdfdType_cd1 , sdfdType_cxd0 , sdfdType_cxd1) == 0)
    if(length(iStat) > 0){
        tmp <- iaCovMatern(setupMats$dIU , parsBTfmd$ad , parsBTfmd$nud , c(1 , 0 , 1) , 0 , setupMats$dIUabcd , setupMats$dIUiUElements)
        phidStat <- tmp$avCovMtx
        avVardStat <- tmp$avVarVec
    }else{}

    if(sdfdType_cd1 == 0){
        phid_cd1 <- phidStat
        avVard_cd1 <- avVardStat
    }else if(sdfdType_cd1 == -9){
        phid_cd1 <- 0
        avVard_cd1 <- 0
    }else{
        tmp <- iaCovMatern(setupMats$dIU , parsBTfmd$ad , parsBTfmd$nud , parsBTfmd$sdfdPars_cd1 , sdfdType_cd1 , setupMats$dIUabcd , setupMats$dIUiUElements)
        phid_cd1 <- tmp$avCovMtx
        avVard_cd1 <- tmp$avVarVec
    }

    if(sdfdType_cxd0 == 0){
        phid_cxd0 <- phidStat
        avVard_cxd0 <- avVardStat
    }else{
        tmp <- iaCovMatern(setupMats$dIU , parsBTfmd$ad , parsBTfmd$nud , parsBTfmd$sdfdPars_cxd0 , sdfdType_cxd0 , setupMats$dIUabcd , setupMats$dIUiUElements)
        phid_cxd0 <- tmp$avCovMtx
        avVard_cxd0 <- tmp$avVarVec
    }

    if(sdfdType_cxd1 == 0){
        phid_cxd1 <- phidStat
        avVard_cxd1 <- avVardStat
    }else{
        tmp <- iaCovMatern(setupMats$dIU , parsBTfmd$ad , parsBTfmd$nud , parsBTfmd$sdfdPars_cxd1 , sdfdType_cxd1 , setupMats$dIUabcd , setupMats$dIUiUElements)
        phid_cxd1 <- tmp$avCovMtx
        avVard_cxd1 <- tmp$avVarVec
    }
    
    parsOK <- T
    if((max(is.na(phid_cd1)) == 1) | (max(is.na(phid_cxd0)) == 1) | (max(is.na(phid_cxd1)) == 1)){ parsOK <- F } else{}

    if(parsOK){
        Kphix0K <- setupMats$Kx %*% phix0 %*% t(setupMats$Kx)
        Cx0 <- parsBTfmd$cx0 * Kphix0K 
        Kphix1K <- setupMats$Kx %*% phix1 %*% t(setupMats$Kx)
        Cx1 <- parsBTfmd$cx1 * Kphix1K 

        if(length(iStat) > 0){
            KphidKStat <- setupMats$Kd %*% phidStat %*% t(setupMats$Kd)	
        }else{}
		
        if(sdfdType_cd1 == 0){
            Cd <- parsBTfmd$cd1 * KphidKStat
        }else if(sdfdType_cd1 == -9){
            Cd <- 0
        }else{
            Cd <- parsBTfmd$cd1 * (setupMats$Kd %*% phid_cd1 %*% t(setupMats$Kd))
        }
        if(sdfdType_cxd0 == 0){
            Cxd0 <- parsBTfmd$cxd0 * Kphix0K * KphidKStat  
        }else{
            Cxd0 <- parsBTfmd$cxd0 * Kphix0K * (setupMats$Kd %*% phid_cxd0 %*% t(setupMats$Kd))
        }
        if(sdfdType_cxd1 == 0){
            Cxd1 <- parsBTfmd$cxd1 * Kphix1K * KphidKStat 
        }else{
            Cxd1 <- parsBTfmd$cxd1 * Kphix1K * (setupMats$Kd %*% phid_cxd1 %*% t(setupMats$Kd))
        }

        n <- nrow(setupMats$Kx)
        C <- Cx0 + Cx1 + Cd + Cxd0 + Cxd1 + parsBTfmd$cme * diag(n)

        sigma2Vec <- parsBTfmd$cx0 + parsBTfmd$cx1 + parsBTfmd$cd1 * avVard_cd1 + parsBTfmd$cxd0 * avVard_cxd0 + parsBTfmd$cxd1 * avVard_cxd1 + parsBTfmd$cme
        sigma2Vec <- setupMats$Kd %*%sigma2Vec 

    }else{        
        C <- NA
        sigma2Vec <- NA
        print('Bad parameters...not sure how this can happen')
        print(parsBTfmd)
    }

    return(list('C' = C , 'sigma2Vec' = sigma2Vec))
#	return(list('C' = C , 'sigma2Vec' = sigma2Vec , 'dC' = list(Cx0 , Cx1 , Cd , Cxd0 , Cxd1) , 'phi' = list(phix0 , phix1 , phid_cd1 , phid_cxd0 , phid_cxd1)))
}

setmuIAK3D <- function(X , vXU , iU , beta , diagC , sigma2Vec){
###########################################################
### sigma2Vec is a length-n vector with the average variances (not average covariances) in each sample support
### with stationary variances, these would all be equal to the sum of all variance parameters. 
### with non-stationary variances, these will vary.
### this function is only used in the lognormal setting
###########################################################
    betavXbeta <- calcbetavXbeta(vXU , beta[iU])
    
    mu <- X %*% beta + 0.5 * (sigma2Vec - diagC + betavXbeta)

    return(mu)
}
	
calcbetavXbeta <- function(vX , beta){
#################################################
### vX is n * p x p    
### vX could contain covariances for all covariates, but this would not be required,
### since only those that vary in sample supports (indexed by iU) contribute to this term
### so usually vXU and betaU will be passed in to function
### this function is only used in the lognormal setting
#################################################
    p <- dim(vX)[[2]]
    n <- dim(vX)[[1]] / p
    vXbeta <- vX %*% beta

    i1 <- kronecker(seq(n) , matrix(1 , p , 1))
    j1 <- seq(n * p)

    tKTEMP <- sparseMatrix(i = i1 , j = j1 , x = 1)
    betavXbeta <- tKTEMP %*% (kronecker(matrix(1, n , 1) , beta) * vXbeta)

    return(betavXbeta) 
}

    
readPars <- function(pars , modelx , nud , sdfdType_cd1 , sdfdType_cxd0 , sdfdType_cxd1 , cmeOpt , prodSum , parBnds , lnTfmdData){

    inext <- 1
    parsOK <- T
    if (modelx == 'matern'){
      ax <- logitab(pars[inext] , parBnds$axBnds[1] , parBnds$axBnds[2] , invt = T); inext <- inext + 1  
      nux <- logitab(pars[inext] , parBnds$nuxBnds[1] , parBnds$nuxBnds[2] , invt = T) ; inext <- inext + 1 
    }else if(modelx == 'nugget'){
### nothing to read in.    
        ax <- nux <- NA
    }else{
        stop('Build in options for other spatial correlation models')
    }

    ad <- logitab(pars[inext] , parBnds$adBnds[1] , parBnds$adBnds[2] , invt = T); inext <- inext + 1  

    if(prodSum){    
      cx0 <- 1E-8 + exp(pars[inext]) ; inext <- inext + 1 ; 
      if(modelx != 'nugget'){
        cx1 <- 1E-8 + exp(pars[inext]) ; inext <- inext + 1 ; 
      }else{
        cx1 <- 0        
      }
    }else{
      cx0 <- 0
      cx1 <- 0
    }
    if(sdfdType_cd1 != -9){        
      cd1 <- 1E-8 + exp(pars[inext]) ; inext <- inext + 1 ; 
    }else{
      cd1 <- 0
    }

    if(modelx != 'nugget'){
      if(!lnTfmdData){
        cxd1 <- logitab(pars[inext] , parBnds$sxd1Bnds[1] , parBnds$sxd1Bnds[2] , invt = T) ; inext <- inext + 1 
        cxd0 <- 1 - cxd1
      }else{ # for lnN methods, read in all variance parameters (as lnci)...
        cxd0 <- 1E-8 + exp(pars[inext]) ; inext <- inext + 1 
        cxd1 <- 1E-8 + exp(pars[inext]) ; inext <- inext + 1 
      }
    }else{
      if(!lnTfmdData){
        cxd0 <- 1
        cxd1 <- 0
      }else{
        cxd0 <- 1E-8 + exp(pars[inext]) ; inext <- inext + 1 
        cxd1 <- 0
      }        
    }
   
    tmp <- sdfdRead(pars , sdfdType_cd1 , inext , parBnds = parBnds)
    sdfdPars_cd1 <- tmp$sdfdPars
    inext <- tmp$inext

    tmp <- sdfdRead(pars , sdfdType_cxd0 , inext , parBnds = parBnds)
    sdfdPars_cxd0 <- tmp$sdfdPars
    inext <- tmp$inext

    tmp <- sdfdRead(pars , sdfdType_cxd1 , inext , parBnds = parBnds)
    sdfdPars_cxd1 <- tmp$sdfdPars
    inext <- tmp$inext

    if (cmeOpt == 1){
      cme <- exp(pars[inext])
    }else{
      cme <- 0
    }

    return(list('ax' = ax , 'nux' = nux , 'ad' = ad , 'nud' = nud , 
                'cx0' = cx0 , 'cx1' = cx1 , 'cd1' = cd1 , 'cxd0' = cxd0 , 'cxd1' = cxd1 ,
                'sdfdPars_cd1' = sdfdPars_cd1 , 'sdfdPars_cxd0' = sdfdPars_cxd0 , 'sdfdPars_cxd1' = sdfdPars_cxd1 ,
                'cme' = cme))
}

##########################################
### the non-stationary standard deviation function, evaluated at depths d
### used if approximating average covariances numerically with dicretization points
### not being used in the current code, but left in as might be used for checking in future.
##########################################
sdfd <- function(d , sdfdPars , sdfdType){

    if (sdfdType == 0){
        sdfdVal <- matrix(1 , length(d) , 1)
    }else if(sdfdType == -1){
        sdfdVal <- (1 - sdfdPars[1]) + sdfdPars[1] * exp(-sdfdPars[2] * d)
    }else if (sdfdType == -9){
        sdfdVal <- matrix(1 , length(d) , 1)
    }else{
        stop('sdfd option not yet programmed')
    }

### check sdfdDsc, return just NA if not ok...
    if (min(sdfdVal) < 0){ sdfdVal <- NA }else{}
	
    return(sdfdVal)
}

sdfdRead <- function(pars , sdfdType , inext , parBnds){
    if (sdfdType == 0){
        sdfdPars <- c()
    }else if(sdfdType == -1){
#############################################
### sdfd <- (1 - tau1) + tau1 * exp(-(tau2 * d))
### to ensure (1 - tau1) + tau1 * exp(-(tau2 * d)) > 0, 
### make tau2 > 0 and furthermore in bounds stated in parBnds$tau2Bnds...param by par2 = logitab(tau2 , tau2min, tau2max)
### and to make sure it is positive for all values in range d = 0 to d = dmax...
### make (1 - tau1) + tau1 * exp(-(tau2 * dmax)) > 0...param by par1 = log[(1 - tau1) + tau1 * exp(-(tau2 * dmax))]
### inverse of latter is tau1 = (1 - exp(par1)) / (1 - exp(-tau2 * dmax))
#############################################
        sdfdPars <- tauTfm(parsIn = pars[inext:(inext+1)] , parBnds = parBnds , invt = T)

        inext <- inext + 2
    }else if(sdfdType == -9){
        sdfdPars <- c()
    }else{
        stop('sdfd option not yet programmed')
    }
	
    return(list('sdfdPars' = sdfdPars , 'inext' = inext))
}


tauTfm <- function(parsIn , parBnds , invt = F){
  nparsIn <- length(parsIn)
  if(nparsIn == 2){
    parsOut <- c(NA , NA)
    if(!invt){
      parsOut[1] <- log(1 - parsIn[1] * (1 - exp(-parsIn[2] * parBnds$maxd)))
      parsOut[2] <- logitab(parsIn[2] , parBnds$tau2Bnds[1]  , parBnds$tau2Bnds[2])
    }else{
      parsOut[2] <- logitab(parsIn[2] , parBnds$tau2Bnds[1]  , parBnds$tau2Bnds[2] , invt = T)
      parsOut[1] <- (1 - exp(parsIn[1])) / (1 - exp(-parsOut[2] * parBnds$maxd))
    }
  }else if(nparsIn == 0){
    parsOut <- c()  
  }else{
    stop('Error - not ready for transforming a tau vector with this number of parameters!')
  }  
  return(parsOut)
}


attachSettings2lmmFit <- function(lmmFit , x , dI , z , covs , modelX , namesX , modelx , nud , allKnotsd , XLims , 
               sdfdType_cd1 , sdfdType_cxd0 , sdfdType_cxd1 , cmeOpt , prodSum , setupMats , lnTfmdData , useReml , parBnds , fitRange , compLikMats){
    lmmFit$x <- x
    lmmFit$dI <- dI
    lmmFit$z <- z
    lmmFit$covs <- covs
    lmmFit$modelX <- modelX
    lmmFit$namesX <- namesX
    lmmFit$modelx <- modelx
    lmmFit$nud <- nud
    if(length(allKnotsd) > 0){
      lmmFit$allKnotsd <- allKnotsd
    }else{ # assigning as c() doesn't work, so: 
      lmmFit['allKnotsd'] <- list(NULL)
    }
    lmmFit$XLims <- XLims
    lmmFit$sdfdType_cd1 <- sdfdType_cd1
    lmmFit$sdfdType_cxd0 <- sdfdType_cxd0
    lmmFit$sdfdType_cxd1 <- sdfdType_cxd1
    lmmFit$cmeOpt <- cmeOpt
    lmmFit$prodSum <- prodSum
    lmmFit$setupMats <- setupMats
    lmmFit$lnTfmdData <- lnTfmdData
    lmmFit$useReml <- useReml
    lmmFit$parBnds <- parBnds
    lmmFit$fitRange <- fitRange
    lmmFit$compLikMats <- compLikMats
    
    return(lmmFit)
}

#########################################################################
### function to iterate optim, each time starting from previous finishing point, until no more change...
########################################################################
### to put this into the global environment, so that it can be seen in function and in optimIt without passing...
verboseOptim <<- F

optimIt <- function(par , fn , methodOptim = c("Nelder-Mead" , "L-BFGS-B") , vecFixedIn = logical(length(par)) , fitRange = matrix(NA , length(par) , 2) , ...){
# fitRange is npar x 2 with lower and upper values
# iterate between first NM and then L-BFGS-B (or other specified algorithms) until no more improvement.
    verboseOptim <<- T
    eval(fn(par,...))
    verboseOptim <<- F
    
    prevSSE = 9E9
    tolIt <- 1E-4 # stop when we improve by tol or less.
    stillImproving = TRUE
    listRes <- list()
    it <- 1
    parStore <- c()
    vecFixed <- vecFixedIn

    while (stillImproving){
        iMethodThis <- it %% length(methodOptim)
        if(iMethodThis == 0){ iMethodThis <- length(methodOptim) }else{}
        if(is.element(methodOptim[iMethodThis] , c("L-BFGS-B" , "Brent"))){
            lower = fitRange[,1] 
            upper = fitRange[,2] 
        }else{
            lower <- -Inf
            upper <- Inf
        }    
    
        res <- optifix(par = par , fn = fn , fixed = vecFixed , method = methodOptim[iMethodThis] , lower = lower , upper = upper , ... = ...)
        listRes[[it]] <- res

        verboseOptim <<- T
        eval(fn(res$fullpars,...))
        verboseOptim <<- F

        parInitsNext <- res$fullpars
### also test if coefficient of variation of last 3 vals is low...       
        parStore <- cbind(parStore , res$fullpars)
        it <- it + 1
        
        if (res$value < (prevSSE - tolIt)){
            stillImproving = TRUE
            prevSSE <- res$value
            par <- parInitsNext
        }else{
            stillImproving = FALSE
        }
    }

    res$allRes <- listRes
    res$parStore <- parStore
            
    return(res)
}

##########################################################
### function to calculate lndetC and invC b (or just invC if b not given) using cholesky...
##########################################################
lndetANDinvCb <- function(C , b = NA){

    if(is.numeric(C)){
      cholC <- try(chol(C) , silent = TRUE)
    }else{
      n <- nrow(C)
### some issue with dgeMatrix matrices, chol(C) was different to cholC[1:n1:n], so...    
      cholC <- try(chol(C[1:n,1:n,drop = FALSE]) , silent = TRUE)
    }
    if (is.character(cholC)){
        lndetC <- invCb <- NA
    }else{
        lndetC <- 2 * sum(log(diag(cholC)))
        if (is.na(as.numeric(b)[1])){        
            invCb <- chol2inv(cholC)
        }else{
### test if C is sparse...
          if(is.matrix(C)){    
### seems to be most efficient way if dense...
            invCb<- backsolve(cholC , forwardsolve(t(cholC) , b))   
          }else{
### but above doesn't work for sparse matrices, so...          
            invCb <- chol2inv(cholC) %*% b 
          }
        }
   }

   return(list('lndetC' = lndetC , 'invCb' = invCb, 'cholC' = cholC))
}

##########################################################
### logit transform function for range [a - b]...
##########################################################
logitab <- function(z , a = 0 , b = 1 , invt = F){
	if(!invt){
		p <- (z - a) / (b - a)
		z <- log(p / (1 - p))
	}else{
		p <- 1 / (1 + exp(-z))
		z <- p * (b - a) + a
	}
	return(z)
}   

##########################################################
### a transform function for the tau parameters...
### sdfd(d=maxd) / sdfd(d=0) = [1 - tau1 ( 1 - exp{-tau2 maxd})] / 1
### bounds can be put on this ratio (sdfdmaxdBnds)
##########################################################
tauTfm <- function(parsIn , parBnds , invt = F){
  nparsIn <- length(parsIn)
  if(nparsIn == 2){
    parsOut <- c(NA , NA)
    if(!invt){
      parsOut[1] <- log(1 - parsIn[1] * (1 - exp(-parsIn[2] * parBnds$maxd)))
      parsOut[2] <- logitab(parsIn[2] , parBnds$tau2Bnds[1]  , parBnds$tau2Bnds[2])
    }else{
      parsOut[2] <- logitab(parsIn[2] , parBnds$tau2Bnds[1]  , parBnds$tau2Bnds[2] , invt = T)
      parsOut[1] <- (1 - exp(parsIn[1])) / (1 - exp(-parsOut[2] * parBnds$maxd))
    }
  }else if(nparsIn == 0){
    parsOut <- c()  
  }else{
    stop('Error - not ready for transforming a tau vector with this number of parameters!')
  }  
  return(parsOut)
}

########################################################################
### sdfd(d=maxd) / sdfd(d=0) = [1 - tau1 ( 1 - exp{-tau2 maxd})] / 1
### bounds can be put on this ratio (sdfdmaxdBnds)
### only used in defining initial parameters... 
### in DEFUNCT version,
### parsIn here are on the transformed scale...
### with tfmed pars, sd(d) = 1 - [1-exp(tau1Tfmd)] * [1 - exp(-tau2 d)] / [1 - exp(-tau2 maxd)]
### and              sd(maxd) = exp(tau1Tfmd)    
########################################################################
contraintau1_DEFUNCT <- function(tau1Tfmd , maxd , sdfmaxdBnds = c(0 , Inf)){
    sdTmpmaxd <- exp(tau1Tfmd)
    if(sdTmpmaxd < sdfmaxdBnds[1]){
      tau1Tfmd <- log(sdfmaxdBnds[1])
    }else if(sdTmpmaxd > sdfmaxdBnds[2]){
      tau1Tfmd <- log(sdfmaxdBnds[2])
    }else{}
    return(tau1Tfmd)
}

#############################################
### in this version, tau is on real scale...
#############################################
contraintau <- function(tau , parBnds , sdfmaxdBnds = c(0 , Inf)){
### check ok for tau2...
    if(tau[2] <= parBnds$tau2Bnds[1]){ 
        tau[2] <- parBnds$tau2Bnds[1] + 0.01 * (parBnds$tau2Bnds[2] - parBnds$tau2Bnds[1]) 
    }else if(tau[2] >= parBnds$tau2Bnds[1]){ 
        tau[2] <- parBnds$tau2Bnds[1] + 0.99 * (parBnds$tau2Bnds[2] - parBnds$tau2Bnds[1]) 
    }else{}

### check ok for tau1...
    sdTmpmaxd <- 1 - tau[1] * (1 - exp(-tau[2] * parBnds$maxd))

    if(sdTmpmaxd < sdfmaxdBnds[1]){ 
        tau[1] <- (1 - sdfmaxdBnds[1]) / (1 - exp(-tau[2] * parBnds$maxd)) 
    }else if(sdTmpmaxd > sdfmaxdBnds[2]){ 
        tau[1] <- (1 - sdfmaxdBnds[2]) / (1 - exp(-tau[2] * parBnds$maxd)) 
    }else{}

    return(tau)
}



############################################################
### uses modes of any factors, means of any non-factors...
############################################################
colMeansAndModes <- function(dfIn , na.rm = FALSE){
  if(ncol(dfIn) > 0){
    jFactor <- c()
    for (j in 1:ncol(dfIn)){
      if(class(dfIn[,j]) == 'factor'){
        jFactor <- c(jFactor , j)
      }else{}
    }
    jNonFactor <- setdiff(1:ncol(dfIn) , jFactor)

    cm <- dfIn[1,]
    if(length(jNonFactor) > 0){
      cm[jNonFactor] <- colMeans(dfIn[,jNonFactor] , na.rm = na.rm)
    }else{}

    if(length(jFactor) > 0){  
      for(j in 1:length(jFactor)){  
        lThis <- levels(dfIn[,jFactor[j]])
        cm[jFactor[j]] <- lThis[which.max(tabulate(match(dfIn[,jFactor[j]] , lThis)))]
      }
    }else{}
  }else{
    cm <- NULL
  }
  return(cm)
}

###########################################################################
### a function to get the sample that is closest to the 'centre'...
###########################################################################
medoid <- function(dfIn){

  if(ncol(dfIn) > 0){
    dfInOrig <- dfIn
 
########################################################################  
### standardise numeric covariates and take abs distance from centre...
### for factors, use 1 - count / countmax, stdardized to have mean = avD
########################################################################  
    avD <- NA * numeric(ncol(dfIn))
    for (j in 1:ncol(dfIn)){
      if(is.numeric(dfIn[,j])){
        dfIn[,j] <- dfIn[,j,drop=FALSE] - mean(dfIn[,j])
        dfIn[,j] <- abs(dfIn[,j,drop=FALSE]) / sqrt(var(dfIn[,j]))

        avD[j] <- mean(dfIn[,j])
      
      }else if(is.factor(dfIn[,j])){
### nothing for now.    
        tableTmp <- table(dfInOrig[,j])
        dfIn[,j] <- NA
        dfIn[,j] <- as.numeric(dfIn[,j])
        for(i in 1:nrow(tableTmp)){
          iThis <- which(dfInOrig[,j] == rownames(tableTmp)[i])
          if(length(iThis) > 0){
            dfIn[iThis,j] <- 1 - tableTmp[i]/max(tableTmp)
          }else{}
        }
      }else{
        stop('Error in medoid function - currently only written for numeric or factor variables!')
      }
    }
  
    if(all(is.na(avD))){ 
      mavD <- 1 
    }else{
      mavD <- mean(avD , na.rm = TRUE)
    }
  
############################################################
### define each sample with a distance from centre
### distance = sum of all univariate distances (city block metric)
### distance for factor is 0 if most common class
############################################################
    jFactors <- which(is.na(avD))
    if(length(jFactors) > 0){
      for (j in jFactors){
        dfIn[,j] <- dfIn[,j] * mavD / mean(dfIn[,j])
      }
    }

### get row with min sum of dists and return that row of the original df...
    imed <- which.min(rowSums(dfIn))

    return(dfInOrig[imed,,drop=FALSE])
  }else{
    return(dfIn[1,,drop=FALSE])
  }
}
 
