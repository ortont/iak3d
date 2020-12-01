selectCovIAK3D <- function(xData , dIData , zData , covsData , modelX , modelx , nud = NULL , sdfdTypeANDcmeInit = c() , allKnotsd = c() , prodSum = TRUE , lnTfmdData , useReml , compLikMats = NULL , rqrBTfmdPreds = TRUE , dirPlot = getwd()){

#########################################################
### error check for duplicated xData,dIData data, if no measurement error being included...
#########################################################
    iDuplicated <- which(duplicated(cbind(xData , dIData)))
    if(length(iDuplicated) > 0){ stop(paste0('Error - the data at positions ' , iDuplicated , ' are duplicated! Crashing for now - but could generalise model selection to force measurement error to be included.')) }else{} 

#########################################################
### error check for any NAs; removing them if found...
#########################################################
    iNA <- which(is.na(zData) | rowSums(is.na(covsData)) > 0)
    if(length(iNA) > 0){ 
      print('Attention! Some NAs found in data ; removing them.')
      xData <- xData[-iNA,drop = FALSE]
      dIData <- dIData[-iNA,drop = FALSE]
      covsData <- covsData[-iNA,drop = FALSE]
      zData <- zData[iNA]
    }else{}
    
    #################################################
    ### the initial model, from input if sdfdTypeANDcmeInit is length 4, else default...
    ### the components represent: sdfdType_cd1 , sdfdType_cxd0 , sdfdType_cxd1 , cmeOpt 
    ###     0 = stationary (ie constant with depth)
    ###    -1 = non-stationary, variance modelled as exponential function of depth
    ###     1 = included (for measurement error) or non-stationary and modelled as linear function of depth for other components (not inlcuded as option yet) 
    ###    -9 = variance componenet not included
    #################################################
    if(length(sdfdTypeANDcmeInit) == 4){
      sdfdTypeANDcmeOpt <- sdfdTypeANDcmeInit
      if(!prodSum){
        sdfdTypeANDcmeOpt[1] <- -9
      }else{}
      
    }else{
      ### start with these values by default.
      if (modelx == 'nugget'){
        #        sdfdTypeANDcmeOpt <- c(-1 , -1 , 0 , 1)  
        sdfdTypeANDcmeOpt <- c(0 , -1 , 0 , 1)  
      }else{
        #        sdfdTypeANDcmeOpt <- c(-1 , -1 , -1 , 1)  
        sdfdTypeANDcmeOpt <- c(0 , -1 , -1 , 1)  
      }
      if(!prodSum){
        sdfdTypeANDcmeOpt[1] <- -9
      }else{}
    }
    
    #################################################
    ### first select nud, if not specified...
    #################################################
    if(is.null(nud)){
      nudVec <- c(0.5 , 1.5 , 2.5)
      lmm.fit.Selectnud <- list()
      for (i in 1:3){
        print(paste0('nud = ' , nudVec[i] , ', sdfdType_cd1 = ' , sdfdTypeANDcmeOpt[1] , ', sdfdType_cxd0 = ' , sdfdTypeANDcmeOpt[2] , 
                     ', sdfdType_cxd1 = ' , sdfdTypeANDcmeOpt[3] , ', cmeOpt = ' , sdfdTypeANDcmeOpt[4] , '...'))
        
        namePlot <- paste0(dirPlot , '/lmm.fit.Selectnud' , floor(nudVec[i]) , '.pdf')
        lmm.fit.Selectnud[[i]] <- fitIAK3D(xData = xData , dIData = dIData , zData = zData , covsData = covsData , modelX = modelX , modelx = modelx , nud = nudVec[i] , allKnotsd = allKnotsd ,
                                           sdfdType_cd1 = sdfdTypeANDcmeOpt[1] , sdfdType_cxd0 = sdfdTypeANDcmeOpt[2] , sdfdType_cxd1 = sdfdTypeANDcmeOpt[3] , cmeOpt = sdfdTypeANDcmeOpt[4] , prodSum = prodSum , 
                                           lnTfmdData = lnTfmdData , useReml = useReml , compLikMats = compLikMats , namePlot = namePlot , rqrBTfmdPreds = rqrBTfmdPreds)
      }
      
      iBest <- which.min(c(lmm.fit.Selectnud[[1]]$lmmFit$nll , lmm.fit.Selectnud[[2]]$lmmFit$nll , lmm.fit.Selectnud[[3]]$lmmFit$nll))
      nudSelected <- nudVec[iBest]
    }else{
      
      namePlot <- paste0(dirPlot , '/lmm.fit.Selectnud' , floor(nud) , '.pdf')
      iBest <- ceiling(nud)
      lmm.fit.Selectnud <- list()
      lmm.fit.Selectnud[[1]] <- NA
      lmm.fit.Selectnud[[2]] <- NA
      lmm.fit.Selectnud[[3]] <- NA
      lmm.fit.Selectnud[[iBest]] <- fitIAK3D(xData = xData , dIData = dIData , zData = zData , covsData = covsData , modelX = modelX , modelx = modelx , nud = nud , allKnotsd = allKnotsd ,
                                             sdfdType_cd1 = sdfdTypeANDcmeOpt[1] , sdfdType_cxd0 = sdfdTypeANDcmeOpt[2] , sdfdType_cxd1 = sdfdTypeANDcmeOpt[3] , cmeOpt = sdfdTypeANDcmeOpt[4] , prodSum = prodSum , 
                                             lnTfmdData = lnTfmdData , useReml = useReml , compLikMats = compLikMats , namePlot = namePlot , rqrBTfmdPreds = rqrBTfmdPreds)
      nudSelected <- nud
      print(paste0('nud = ' , nudSelected , ' selected.'))
    }
    
    ### to be updated as (if) submodels get accepted...
    lmm.fit.Selected <- lmm.fit.Selectnud[[iBest]]$lmmFit
    lmm.fit.Selected.AllInfo <- lmm.fit.Selectnud[[iBest]]
    
    sdfdTypeANDcmeOptSelected <- sdfdTypeANDcmeOpt
    ### values that will be updated for comparison...
    nllCurrent <- lmm.fit.Selectnud[[iBest]]$lmmFit$nll
    
    ### counting nud as a parameter, full model (if matern) should have 16 covariance parameters
    ### ie 1 (nud) + 3 (other corr pars) + 5 (var pars) + 6 (non-stat pars) + 1 (me par) covariance parameters
    if(lnTfmdData){
      pCurrent <- length(lmm.fit.Selectnud[[iBest]]$lmmFit$betahat + length(lmm.fit.Selectnud[[iBest]]$parsFit$par) + 1)
    }else{
      pCurrent <- length(lmm.fit.Selectnud[[iBest]]$lmmFit$betahat + length(lmm.fit.Selectnud[[iBest]]$parsFit$par) + 2) # as 1 var param was done automatically
    }
    aicCurrent <- 2 * nllCurrent + 2 * pCurrent
    
    ##################################################
    ### with nud selected, see if we can drop terms from the fixed-effect model...
    ### though if modelX was a cubist model, will have to make some changes to the function, 
    ### so for the moment, only do this if modelX was a linear model
    ##################################################
    
    ##################################################
    ### now see if we can drop terms (non-stationary variances or measurement error component) from the covariance model...
    ### fit up to 4 submodels and drop terms until no further reductions in AIC...
    ##################################################
    continueDropping <- T
    
    ### list of all models fitted on the way. lmm.fit.Drop[[1]] will itself be a list of the <= 4 models fitted in the first step
    lmm.fit.Drop <- list()
    fitWarnings <- nllDrop <- pDrop <- aicDrop <- matrix(NA , 4 , 4)
    iStep <- 1
    while(continueDropping){
      
      lmm.fit.Drop[[iStep]] <- list()
      for (i in 1:4){
        if(!is.element(sdfdTypeANDcmeOpt[i] , c(0 , -9))){
          sdfdTypeANDcmeOptThis <- sdfdTypeANDcmeOpt
          sdfdTypeANDcmeOptThis[i] <- 0
          
          print(paste0('nud = ' , nudSelected , ', sdfdType_cd1 = ' , sdfdTypeANDcmeOptThis[1] , ', sdfdType_cxd0 = ' , sdfdTypeANDcmeOptThis[2] , 
                       ', sdfdType_cxd1 = ' , sdfdTypeANDcmeOptThis[3] , ', cmeOpt = ' , sdfdTypeANDcmeOptThis[4] , '...'))
          
          
          namePlot <- paste0(dirPlot , '/lmm.fit.DropStep' , iStep , '_Drop' , i , '.pdf')
          lmm.fit.Drop[[iStep]][[i]] <- fitIAK3D(xData = xData , dIData = dIData , zData = zData , covsData = covsData , modelX = modelX , modelx = modelx , nud = nudSelected , allKnotsd = allKnotsd ,
                                                 sdfdType_cd1 = sdfdTypeANDcmeOptThis[1] , sdfdType_cxd0 = sdfdTypeANDcmeOptThis[2] , sdfdType_cxd1 = sdfdTypeANDcmeOptThis[3] , cmeOpt = sdfdTypeANDcmeOptThis[4] , prodSum = prodSum , 
                                                 lnTfmdData = lnTfmdData , useReml = useReml  , compLikMats = compLikMats , namePlot = namePlot , rqrBTfmdPreds = rqrBTfmdPreds)
          nllDrop[i,iStep] <- lmm.fit.Drop[[iStep]][[i]]$lmmFit$nll
          
          if(nllDrop[i,iStep] < (nllCurrent - 1E-4)){
            print('Failure of the nesting model fit!!!! This simpler model has smaller nll than the more complex one. Check this!')
            
            ### perhaps in this case we should return and refit the nesting model, and check whether that was the best
            fitWarnings[i,iStep] <- 1
          }else{
            fitWarnings[i,iStep] <- 0
          }
          
          if(i == 4){
            pDrop[i,iStep] <- pCurrent - 1
          }else{
            pDrop[i,iStep] <- pCurrent - 2
          }    
          aicDrop[i,iStep] <- 2 * nllDrop[i,iStep] + 2 * pDrop[i,iStep]
        }else{}
      }
      
      ### compare all to the current best via AIC...
      if (min(aicDrop[,iStep] , na.rm = T) < aicCurrent){
        ### iBest is the parameter which when dropped leaves the best aic.        
        iBest <- which.min(aicDrop[,iStep])
        sdfdTypeANDcmeOpt[iBest] <- 0
        nllCurrent <- nllDrop[iBest,iStep]
        pCurrent <- pDrop[iBest,iStep]
        aicCurrent <- aicDrop[iBest,iStep]
        lmm.fit.Selected <- lmm.fit.Drop[[iStep]][[iBest]]$lmmFit
        lmm.fit.Selected.AllInfo <- lmm.fit.Drop[[iStep]][[iBest]]
        sdfdTypeANDcmeOptSelected <- sdfdTypeANDcmeOpt
        namesTmp <- c('sdfdType_cd1' , 'sdfdType_cxd0' , 'sdfdType_cxd1' , 'cmeOpt')
        print(paste0('Dropping ' , namesTmp[iBest] , ' from model.'))
      }else{
        continueDropping <- F
      }
      
      if(all(is.element(sdfdTypeANDcmeOpt , c(0 , -9)))){
        continueDropping <- F
      }else{}
      
      iStep <- iStep + 1
    }
    
    ### return...
    return(list('lmmSelected' = lmm.fit.Selected , 'lmm.fit.Selected.AllInfo' = lmm.fit.Selected.AllInfo , 'nud' = nudSelected , 'sdfdTypeANDcmeOptSelected' = sdfdTypeANDcmeOptSelected , 
                'nllDrop' = nllDrop , 'pDrop' = pDrop , 'aicDrop' = aicDrop , 'fitWarnings' = fitWarnings ,
                'lmm.fit.Selectnud' = lmm.fit.Selectnud , 'lmm.fit.Drop' = lmm.fit.Drop))
}

################################################################################
### here the models with sdfdTypeANDcmeOpt[1] == -9 are also included 
### (-9 = no depth rand eff trend ; useful if spline included for approx with big datasets), 
### which weren't searched for in the above model selection algorithm.
### function useful for parallelizing model selection.
################################################################################
listAllSubmodels <- function(modelx , nud = NULL , sdfdTypeANDcmeInit = c() , allKnotsd = c() , prodSum = TRUE){

#################################################
### the initial model, from input if sdfdTypeANDcmeInit is length 4, else default...
### the components represent: sdfdType_cd1 , sdfdType_cxd0 , sdfdType_cxd1 , cmeOpt 
###     0 = stationary (ie constant with depth)
###    -1 = non-stationary, variance modelled as exponential function of depth
###     1 = included (for measurement error) or non-stationary and modelled as linear function of depth for other components (not inlcuded as option yet) 
###    -9 = variance componenet not included
#################################################
    if(length(sdfdTypeANDcmeInit) == 4){
      sdfdTypeANDcmeOpt <- sdfdTypeANDcmeInit
      if(!prodSum){
        sdfdTypeANDcmeOpt[1] <- -9
      }else{}

    }else{
### start with these values by default.
      if (modelx == 'nugget'){
       sdfdTypeANDcmeOpt <- c(0 , -1 , 0 , 1)  
      }else{
        sdfdTypeANDcmeOpt <- c(0 , -1 , -1 , 1)  
      }
      if(!prodSum){
        sdfdTypeANDcmeOpt[1] <- -9
      }else{}
    }

    sdfdTypeANDcmeOptAll <- matrix(sdfdTypeANDcmeOpt , nrow = 1)
    
    if(sdfdTypeANDcmeOpt[1] == -1){
      sdfdTypeANDcmeOptTmp <- sdfdTypeANDcmeOptAll

      sdfdTypeANDcmeOptTmp[,1] <- 0
      sdfdTypeANDcmeOptAll <- rbind(sdfdTypeANDcmeOptAll , sdfdTypeANDcmeOptTmp)
      
      sdfdTypeANDcmeOptTmp[,1] <- -9
      sdfdTypeANDcmeOptAll <- rbind(sdfdTypeANDcmeOptAll , sdfdTypeANDcmeOptTmp)
    }else if(sdfdTypeANDcmeOpt[1] == 0){
      sdfdTypeANDcmeOptTmp <- sdfdTypeANDcmeOptAll

      sdfdTypeANDcmeOptTmp[,1] <- -9
      sdfdTypeANDcmeOptAll <- rbind(sdfdTypeANDcmeOptAll , sdfdTypeANDcmeOptTmp)
    }else{}

    if(sdfdTypeANDcmeOpt[2] == -1){
      sdfdTypeANDcmeOptTmp <- sdfdTypeANDcmeOptAll

      sdfdTypeANDcmeOptTmp[,2] <- 0
      sdfdTypeANDcmeOptAll <- rbind(sdfdTypeANDcmeOptAll , sdfdTypeANDcmeOptTmp)
    }else{}

    if(sdfdTypeANDcmeOpt[3] == -1){
      sdfdTypeANDcmeOptTmp <- sdfdTypeANDcmeOptAll

      sdfdTypeANDcmeOptTmp[,3] <- 0
      sdfdTypeANDcmeOptAll <- rbind(sdfdTypeANDcmeOptAll , sdfdTypeANDcmeOptTmp)
    }else{}

    if(sdfdTypeANDcmeOpt[4] == 1){
      sdfdTypeANDcmeOptTmp <- sdfdTypeANDcmeOptAll

      sdfdTypeANDcmeOptTmp[,4] <- 0
      sdfdTypeANDcmeOptAll <- rbind(sdfdTypeANDcmeOptAll , sdfdTypeANDcmeOptTmp)
    }else{}

    if(is.null(nud)){
      
      nudAll <- c(rep(0.5 , nrow(sdfdTypeANDcmeOptAll)) , rep(1.5 , nrow(sdfdTypeANDcmeOptAll)) , rep(2.5 , nrow(sdfdTypeANDcmeOptAll)))
      sdfdTypeANDcmeOptAll <- rbind(sdfdTypeANDcmeOptAll , sdfdTypeANDcmeOptAll , sdfdTypeANDcmeOptAll)

    }else{
      nudAll <- rep(nud , nrow(sdfdTypeANDcmeOptAll))
    }

    return(list('sdfdTypeANDcmeOptAll' = sdfdTypeANDcmeOptAll , 'nudAll' = nudAll))
} 

#############################################################
### a backwards elimination procedure for fixed-effect model based on AIC and independent residuals...
### to do, allow lmmFit object to be passed in, keep all cov params fixed throughout, but refit beta
### (a conditional AIC? Wald test with alpha = 0.15 prob quicker and v similar)
### assumes additive normal (ie arithmetic averaging) effects on given scale 
#############################################################
selectXAicIAK3D <- function(xData , dIData , zData , covsData , modelXInit , allKnotsd = c()){

    tmp <- makeXvX(covData = covsData , dIData = dIData , modelX = modelXInit, allKnotsd = allKnotsd , nDiscPts = 10 , lnTfmdData = FALSE)

    namesXCurrent <- tmp$namesX
    XCurrent <- tmp$X
    pXCurrent <- tmp$pX

    tmp <- nllLm(zData = zData , XData = XCurrent , REML = FALSE)
    nllCurrent <- tmp$nll
    pCurrent <- length(tmp$betahat)

    continueRemoving <- T

    nllStore <- c()
    namesXStore <- list() ; iNext <- 1
    while (continueRemoving){
        tmp <- stepBackAic(namesX = namesXCurrent , pX = pXCurrent , allKnotsd = allKnotsd , zData = zData , dIData = dIData , covsData = covsData , nllCurrent = nllCurrent , pCurrent = pCurrent)
        if(length(tmp$iRemove) > 0){
            namesXCurrent <- tmp$namesX 
            pXCurrent <- tmp$pX
            nllCurrent <- tmp$nllNew 
            pCurrent <- tmp$pNew 
            XCurrent <- tmp$XNew
            nllStore <- c(nllStore , nllCurrent)
            namesXStore[[iNext]] <- namesXCurrent 
            iNext <- iNext + 1
        }else{
            continueRemoving <- F
        }
    }
    return(namesXCurrent)
}

#############################################################
### one step of the backwards elimination procedure for fixed-effect model based on AIC and independent residuals...
### assumes additive normal (ie arithmetic averaging) effects on given scale 
#############################################################
stepBackAic <- function(namesX , pX , allKnotsd , zData , dIData , covsData , nllCurrent , pCurrent){
### the nll and number of parameters in the current model 
### fit all immediate submodels; aic for backward elimination
    aicCurrent <- 2 * nllCurrent + 2 * pCurrent

    nllAll <- pAll <- NA * numeric(length(namesX))
    crAll <- canRemove(namesX) 
    XList <- list()
    for (i in 1:length(namesX)){
        if((crAll[i]) & (!is.na(pX[i]))){
            iRmvTmp <- which(namesX == namesX[i])
            namesXThis <- namesX[-iRmvTmp]
            pXThis <- pX[-iRmvTmp]

            tmp <- makeXvX(covData = covsData , dIData = dIData , modelX = namesX , allKnotsd = allKnotsd , nDiscPts = 10 , lnTfmdData = FALSE)
            XList[[i]] <- tmp$X
            tmp <- nllLm(zData = zData , XData = XList[[i]] , REML = F)
            nllAll[i] <- tmp$nll
            pAll[i] <- pCurrent - pX[i]
        }else{}
    }
    aicAll <- 2 * nllAll + 2 * pAll
    if(all(is.na(aicAll))){
        iRemove <- c()
        nllNew <- nllCurrent
        pNew <- pCurrent
        XNew <- NA # whatever it was before.
    }else{
        if(min(aicAll , na.rm = T) < aicCurrent){
          iRemove <- which.min(aicAll)
### get all of this name, in case its a categorical variable being removed...
          nameXRmv <- namesX[iRemove]
          iRemove <- which(namesX == nameXRmv)

          namesX <- namesX[-iRemove]
          pX <- pX[-iRemove]
          nllNew <- nllAll[iRemove[1]]
          pNew <- pAll[iRemove[1]]
          XNew <- XList[[iRemove[1]]]
        }else{
          iRemove <- c()
          nllNew <- nllCurrent
          pNew <- pCurrent
          XNew <- NA # whatever it was before.
        }
    }
    return(list('namesX' = namesX , 'pX' = pX , 'iRemove' = iRemove , 'nllNew' = nllNew , 'pNew' = pNew , 'XNew' = XNew))
}

#############################################################
### a backwards elimination procedure for fixed-effect model based on Wald tests and independent residuals...
### can also be used with covariance parameters in lmmFit used to model residuals...
### using alpha = 0.15 should give something similar to AIC
#############################################################
selectXWaldIAK3D <- function(xData , dIData , zData , covsData , modelXInit = c() , allKnotsd = c() , lmmFit = list() , alpha = 0.05){

    if(!is.null(lmmFit$parsBTfmd)){
### a lmm has been fitted - get modelXInit from lmmFit, stop if modelXInit is given as well
        if (length(modelXInit) > 0){
            stop('When giving a lmmFit to selectXAicIAK3D, do not also give a modelXInit - this will be taken from the lmmFit object!')
        }else{}
        if (length(allKnotsd) > 0){
            stop('When giving a lmmFit to selectXAicIAK3D, do not also give an allKnotsd - this will be taken from the lmmFit object!')
        }else{}
        modelXInit <- lmmFit$modelX
        allKnotsd <- lmmFit$allKnotsd
        lnTfmdData <- lmmFit$lnTfmdData    
    }else{
        lnTfmdData <- FALSE
    }

    tmp <- makeXvX(covData = covsData , dIData = dIData , modelX = modelXInit, allKnotsd = allKnotsd , nDiscPts = 10 , lnTfmdData = lnTfmdData)

    namesXCurrent <- tmp$namesX
    pXCurrent <- tmp$pX

    continueRemoving <- T

    namesXStore <- list() ; iNext <- 1
    pValStore <- list()
    while (continueRemoving){
        tmp <- stepBackWald(namesX = namesXCurrent , pX = pXCurrent, allKnotsd = allKnotsd , zData = zData , dIData = dIData , covsData = covsData , lmmFit = lmmFit , alpha = alpha)

        if(length(tmp$iRemove) > 0){
            namesXCurrent <- tmp$namesX 
            pXCurrent <- tmp$pX
            namesXStore[[iNext]] <- namesXCurrent 
            pValStore[[iNext]]<- tmp$pValAll
            iNext <- iNext + 1
        }else{
            continueRemoving <- F
        }
    }
    return(namesXCurrent)
}

#############################################################
### one step of the backwards elimination procedure for fixed-effect model based on Wald tests and independent residuals...
#############################################################
stepBackWald <- function(namesX , pX , allKnotsd , zData , dIData , covsData , lmmFit , alpha){
### fit the current model; wald tests for backward elimination
    if(is.null(lmmFit$parsBTfmd)){
### in this case, fit a lm to give betahat and vbetahat
        tmp <- makeXvX(covData = covsData , dIData = dIData , modelX = namesX , allKnotsd = allKnotsd , nDiscPts = 10 , lnTfmdData = FALSE)
        XThis <- tmp$X

        tmp <- nllLm(zData = zData , XData = XThis , REML = T)
        betahat <- tmp$betahat
        vbetahat <- tmp$vbetahat
    }else{
### in this case, use fitted model in lmmFit, and calculate give betahat and vbetahat
        lmmFit$modelX <- namesX
### to make sure these will be re calculated in the fit function (which will just refit betahat, conditional on cov parameters.
        lmmFit$X <- NULL
        lmmFit$vXU <- NULL
        lmmFit$iU <- NULL

        tmp <- fitIAK3D(xData = lmmFit$xData , dIData = lmmFit$dIData , zData = lmmFit$zData , covsData = lmmFit$covsData , modelX = namesX , modelx = lmmFit$modelx , nud = lmmFit$nud , allKnotsd = lmmFit$allKnotsd , 
			sdfdType_cd1 = lmmFit$sdfdType_cd1 , sdfdType_cxd0 = lmmFit$sdfdType_cxd0 , sdfdType_cxd1 = lmmFit$sdfdType_cxd1 , cmeOpt = lmmFit$cmeOpt , prodSum = lmmFit$prodSum , 
                  lnTfmdData = lmmFit$lnTfmdData  , useReml = lmmFit$useReml , compLikMats = lmmFit$compLikMats , lmmFit = lmmFit)

        lmmFit <- tmp$lmmFit

        betahat <- lmmFit$betahat
        vbetahat <- lmmFit$vbetahat
    }

    crAll <- canRemove(namesX) 
    pValAll <- NA * numeric(length(namesX))
    for (i in 1:length(namesX)){
        if((crAll[i]) & (!is.na(pX[i]))){
          iRmvTmp <- which(namesX == namesX[i])

          pValAll[i] <- waldTest(betahat = betahat , vbetahat = vbetahat , ip0 = iRmvTmp)
        }else{}
    }

    if(all(is.na(pValAll))){
        iRemove <- c()
    }else{
        if(max(pValAll , na.rm = T) > alpha){
          iRemove <- which.max(pValAll)

### get all of this name, in case its a categorical variable being removed...
          nameXRmv <- namesX[iRemove]
          iRemove <- which(namesX == nameXRmv)

          namesX <- namesX[-iRemove]
          pX <- pX[-iRemove]
        }else{
          iRemove <- c()
        }
    }
#####################################
### and return...
#####################################
    return(list('namesX' = namesX , 'pX' = pX , 'iRemove' = iRemove , 'pValAll' = pValAll))
}

##########################################################################
### function to give the Wald test p value, given betahat, vbetahat...
### and the which elements of betahat we want to test are equal to 0....
### must specifiy exactly one of ip0 or L
##########################################################################
waldTest <- function(betahat , vbetahat , ip0 = NULL , L = NULL){
  if(is.null(ip0) & is.null(L)){ stop('Error, specify exactly one of ip0 or L for Wald test!') }else{} 
  if((!is.null(ip0)) & (!is.null(L))){ stop('Error, specify exactly one of ip0 or L for Wald test!') }else{} 

  if(is.null(L)){
    L <- matrix(0 , length(ip0) , length(betahat))
    for (j in 1:length(ip0)){
      L[j , ip0[j]] <- 1
    }
    diffp <- length(ip0)
  }else{
    diffp <- dim(L)[[1]]
  }

  v <- L %*% vbetahat %*% t(L)
  Lbetahat <- L %*% betahat
  WaldStat <- as.numeric(t(Lbetahat) %*% solve(v , Lbetahat))
  
  pValWT <- 1 - pchisq(WaldStat , diffp)
  
  return(pValWT)
}

#############################################################
### the nll of the linear model, and its associated parameters...
#############################################################
nllLm <- function(zData , XData , REML = FALSE){
    n <- length(zData)
    p <- dim(XData)[[2]]
    XX <- t(XData) %*% XData 
    cholXX <- try(chol(XX) , silent = TRUE)
    if(is.character(cholXX) | min(eigen(XX)$value) <= 0){
      return(list('nll' = NA , 'betahat' = matrix(NA , p , 1) , 'vbetahat' = matrix(NA , p , p) , 'sigma2hat' = NA))
    }else{}
    vbetahat <- chol2inv(cholXX)
    betahat <- matrix(vbetahat %*% (t(XData) %*% zData) , ncol = 1)

    lndetXX <- determinant(XX , logarithm = TRUE)
    if(lndetXX$sign < 0){
      return(list('nll' = NA , 'betahat' = matrix(NA , p , 1) , 'vbetahat' = matrix(NA , p , p) , 'sigma2hat' = NA))
    }else{}
    lndetXX <- as.numeric(lndetXX$modulus)
    
    res <- zData - XData %*% betahat
    if(REML){
        sigma2hat <- as.numeric(t(res) %*% res) / (n - p)
        nll <- 0.5 * ((n - p) * log(2 * pi) + (n - p) * log(sigma2hat) + lndetXX + n - p)
    }else{
        sigma2hat <- as.numeric(t(res) %*% res) / n
        nll <- 0.5 * n * (log(2 * pi) + log(sigma2hat) + 1)
    }
    if(sigma2hat < 0){
      return(list('nll' = NA , 'betahat' = matrix(NA , p , 1) , 'vbetahat' = matrix(NA , p , p) , 'sigma2hat' = NA)) 
    }else{}
    vbetahat <- vbetahat * sigma2hat

    return(list('nll' = nll , 'betahat' = betahat , 'vbetahat' = vbetahat , 'sigma2hat' = sigma2hat))
}

#############################################################
### calculate a logical vector saying which of the variables in 'namesX' can be removed from the model
### assumes terms up to d2 (i.e. d ^ 2 or spatialVariable * (d ^ 2)) included. 
#############################################################
canRemove <- function(namesX){
  canRemove <- logical(length(namesX))
# start at 2, as can't remove the constant.
  for(j in 2:length(namesX)){
    if(namesX[j] == 'd2'){
### search for any terms that start with d2.; if found, can't remove d2
        if(length(which(substr(namesX , 1 , 3) == 'd2.')) == 0){
            canRemove[j] <- T
        }else{}
    }else if(substr(namesX[j] , 1 , 3) == 'd2.'){
### can remove any d2 interaction terms.
        canRemove[j] <- T
    }else if(namesX[j] == 'd'){
### search for any terms that start with d.; if found, can't remove d
        if(length(which(substr(namesX , 1 , 2) == 'd.')) == 0){
            canRemove[j] <- T
        }else{}
    }else if(substr(namesX[j] , 1 , 2) == 'd.'){
        sptlNameThis <- substr(namesX[j] , 3 , nchar(namesX[j]))
### search for d2.name  
        if(length(which(namesX == paste0('d2.', sptlNameThis))) == 0){
            canRemove[j] <- T
        }else{}
    }else if(namesX[j] == 'dSpline'){
### cannot remove spline terms...    
    }else{
### just the name of a variable, so search for 'd.variable' or 'd2.variable' (shouldn't need the latter search)         
        if(length(which((namesX == paste0('d.', namesX[j])) | (namesX == paste0('d2.', namesX[j])))) == 0){
            canRemove[j] <- T
        }else{}
    }
  }
  return(canRemove)
}


####################################################
### function to select knots for general depth trend...
####################################################
selectKnots <- function(dIData , zData , covsData , modelX , plotSplines = FALSE , degreeSpline = 3 , usedMdPts = FALSE){

### calc initial residuals based on independence + no spline...
    if(usedMdPts){
        tmp <- makeXvX(covData = covsData , dIData = dIData , modelX = modelX , allKnotsd = c() , nDiscPts = 1 , lnTfmdData = FALSE)
    }else{
        tmp <- makeXvX(covData = covsData , dIData = dIData , modelX = modelX , allKnotsd = c() , nDiscPts = 10 , lnTfmdData = FALSE)
    }
    lmNoSpline <- lmGivenX(zData = zData , XData = tmp$X , method = 'REML')
    zRes <- zData - tmp$X %*% lmNoSpline$betahat

### unique values of dIMdPt, rounded to nearest 5 cm...
    dIMdPt <- rowMeans(dIData)
    dIMdPt <- round(dIMdPt * 20) / 20
    
    dIMdPtU <- unique(dIMdPt)
    dIMdPtU <- dIMdPtU[order(dIMdPtU)]

    knotsInit <- 0.5 * (dIMdPtU[-1] + dIMdPtU[-length(dIMdPtU)])
    
    zResdIMdPtU <- NA * numeric(length(dIMdPtU))
    for(i in 1:length(dIMdPtU)){
        iThis <- which(dIMdPt == dIMdPtU[i])
        zResdIMdPtU[i] <- mean(zRes[iThis])
    }
    
    bdryKnots <- c(0 , max(dIData[,2]))
    dPred <- seq(bdryKnots[1] , bdryKnots[2] , (bdryKnots[2] - bdryKnots[1]) / 100)

    knotsInit <- knotsInit[seq(2 , length(knotsInit) - 1 , 2)]

    if(usedMdPts){
        XInit <- bs(x = dIMdPt , knots = knotsInit , degree = degreeSpline , intercept = T , Boundary.knots = bdryKnots)
    }else{
        allKnotsd <- c(bdryKnots[1] , knotsInit , bdryKnots[2])
        tmp <- makeXvX(covData = NULL , dIData = dIData , modelX = 'const' , allKnotsd = allKnotsd , nDiscPts = 10 , lnTfmdData = FALSE)
        XInit <- tmp$X
    }
    lmInit <- lmGivenX(zData = zRes , XData = XInit , method = 'REML')
    aicInit <- lmInit$AIC

    if(plotSplines){
        XPred <- bs(x = dPred , knots = knotsInit , degree = degreeSpline , intercept = T , Boundary.knots = bdryKnots)
        zPred <- XPred %*% lmInit$betahat
        dev.new()
        plot(zRes , -dIMdPt , ylim = c(-bdryKnots[2] , 0))
        points(zResdIMdPtU , -dIMdPtU , col = 'magenta' , pch = 4)
        lines(zPred , -dPred , col = 'black' , lty = 1 , lwd = 2)
    }else{}

    knotsCurrent <- knotsInit
    aicCurrent <- aicInit
    
    stillRemoving <- TRUE
    while(stillRemoving){
        aicTmp <- NA * numeric(length(knotsCurrent))
        for(j in 1:length(knotsCurrent)){
          if(usedMdPts){
            XRed <- bs(x = dIMdPt , knots = knotsCurrent[-j] , degree = degreeSpline , intercept = T , Boundary.knots = bdryKnots)
          }else{
            allKnotsd <- c(bdryKnots[1] , knotsCurrent[-j] , bdryKnots[2])
            tmp <- makeXvX(covData = NULL , dIData = dIData , modelX = 'const' , allKnotsd = allKnotsd , nDiscPts = 10 , lnTfmdData = FALSE)
            XRed <- tmp$X
          }
          lmRed <- lmGivenX(zData = zRes , XData = XRed , method = 'REML' , XFULL = XInit)
          aicTmp[j] <- lmRed$AIC
        }
        
        if(min(aicTmp , na.rm = TRUE) < aicCurrent){
            iRemove <- which.min(aicTmp) 
            knotsCurrent <- knotsCurrent[-iRemove]
            aicCurrent <- aicTmp[iRemove]
            if(length(knotsCurrent) < 2){
                stillRemoving <- FALSE
                print('Removed all but one knot - stopping process now!')
            }else{}        
        }else{
            stillRemoving <- FALSE
        }    
    }

    if(plotSplines){
        if(usedMdPts){
          XRed <- bs(x = dIMdPt , knots = knotsCurrent , degree = degreeSpline , intercept = T , Boundary.knots = bdryKnots)
        }else{      
          allKnotsd <- c(bdryKnots[1] , knotsCurrent , bdryKnots[2])
          tmp <- makeXvX(covData = NULL , dIData = dIData , modelX = 'const' , allKnotsd = allKnotsd , nDiscPts = 10 , lnTfmdData = FALSE)
          XRed <- tmp$X
        }  
        lmRed <- lmGivenX(zData = zRes , XData = XRed , method = 'REML' , XFULL = XInit)
        
        XPred <- bs(x = dPred , knots = knotsCurrent , degree = degreeSpline , intercept = T , Boundary.knots = bdryKnots)
        zPred <- XPred %*% lmRed$betahat
        lines(zPred , -dPred , col = 'red' , lty = 2 , lwd = 2)
    }else{}
    
    allKnotsd <- c(bdryKnots[1] , knotsCurrent , bdryKnots[2])

    return(allKnotsd)
}

####################################################
### function to fit a lm given the design matrix XData
### also allows comparison of REML values betwen different XData's if both nested in XFULL...
####################################################
lmGivenX <- function(zData , XData , method = 'ML' , XFULL = NULL){

    zData <- matrix(zData , ncol = 1)
    n <- length(zData)
    p <- ncol(XData)
    XX <- t(XData) %*% XData
    Xz <- matrix(t(XData) %*% zData , ncol = 1)
    
# #    tmp0 <- try(solve(XX , Xz) , silent = TRUE)
# #    tmp <- lndetANDinvCb(XX , Xz)
# #    if(is.character(tmp$cholC) | is.character(tmp0)){
#     betahat <- try(solve(XX , Xz) , silent = TRUE)
#     if(is.character(betahat)){
#       return(list('nll' = NA , 'AIC' = NA , 'betahat' = matrix(NA , p , 1)))
#     }else{}
# #    lndetXX <- tmp$lndetC
    # lndetXX <- as.numeric(determinant(XX , logarithm = TRUE)$modulus)
    
    betahat <- lndetANDinvCb(XX , Xz)
    lndetXX <- betahat$lndetC
    betahat <- betahat$invCb
    
    if(is.na(lndetXX) || (is.infinite(lndetXX))){
      return(list('nll' = NA , 'AIC' = NA , 'betahat' = matrix(NA , p , 1)))
    }else{}
    betahat <- matrix(betahat , ncol = 1)
    
    zRes <- zData - XData %*% betahat
    if(method == 'ML'){
        if(!is.null(XFULL)){ stop('Error - XFULL should only be included for comparison based on REML fits!') }else{}
        sigma2hat <- as.numeric(t(zRes) %*% zRes / n)
        nll <- 0.5 * n * (log(2 * pi) + log(sigma2hat) + 1)
        
        vbetahat <- sigma2hat * solve(XX)
        
    }else if(method == 'REML'){
        sigma2hat <- as.numeric(t(zRes) %*% zRes / (n - p))
        vbetahat <- sigma2hat * solve(XX)
        
        if(is.null(XFULL)){
            nll <- 0.5 * (n - p) * (log(2 * pi) + log(sigma2hat) + 1) + 0.5 * lndetXX
        }else{
### could add check here that XData is nested in XFULL. But for now assuming this is so. 
            pFULL <- ncol(XFULL)
            XXFULL <- t(XFULL) %*% XFULL
            XzFULL <- t(XFULL) %*% zData
#            tmp <- lndetANDinvCb(XXFULL , XzFULL)
#            lndetXXFULL <- tmp$lndetC
             lndetXXFULL <- as.numeric(determinant(XXFULL , logarithm = TRUE)$modulus)

            nll <- 0.5 * ((n - pFULL) * log(2 * pi) + (n - pFULL) * log(sigma2hat) + (n - p)) + 0.5 * lndetXXFULL
        }
    }else{
        stop('Method should be ML or REML!')
    }

    return(list('nll' = nll , 'AIC' = 2 * nll + 2 * p , 'AICc' = 2 * nll + 2 * p + 2 * p * (p + 1) / (n - p - 1) , 'betahat' = betahat , 'vbetahat' = vbetahat))
}
