###################################################################
### the nll function via composite likelihood - nll called for all pairs of adjacent subsets...
### non-adjacent subsets included but with cov[subset i , subset j] = 0
### should have no depth sum component in cov model (include via spline instead)
###################################################################
nllIAK3D_CL <- function(pars , zData , XData , modelx , nud ,  
                     sdfdType_cd1 , sdfdType_cxd0 , sdfdType_cxd1 , cmeOpt , prodSum , setupMats = NULL , parBnds , useReml , compLikMats , rtnAll = F , attachBigMats = TRUE){

    if(exists("printnllTime") && printnllTime){
      ptm <- proc.time()
    }else{}

    if(exists("lnTfmdData") && lnTfmdData){
      stop('Error - lnTfmdData entered as TRUE for composite likelihood - this option not available with composite likelihood!')
    }else{
      lnTfmdData <- FALSE
    }

    if(exists("parsTrace4Optim")){
      parsTmp <- rbind(parsTrace4Optim , pars)
      assign("parsTrace4Optim" , parsTmp , envir = .GlobalEnv)

      nllTmp <- c(nllTrace4Optim , NA) # for this trace, use NA for failed evaluation.
      assign("nllTrace4Optim" , nllTmp , envir = .GlobalEnv)
    }else{}
    
    badnll <- 9E99
    cxdhat <- NA
       
    n <- length(zData)
    p <- ncol(XData)

    XziAXz_SUM <- matrix(0 , ncol(XData)+1 , ncol(XData)+1)
    lndetA_SUM <- n_SUM <- 0
    if(rtnAll){ listiCkl_kl <- listCkl_kl <- listiAXkl_kl <- listiAzkl_kl <- list() }else{}

    if(compLikMats$compLikOptn < 4){
      for(i in 1:nrow(compLikMats$subsetPairsAdj)){
        iThis <- c(compLikMats$listBlocks[[compLikMats$subsetPairsAdj[i,1]]]$i , compLikMats$listBlocks[[compLikMats$subsetPairsAdj[i,2]]]$i)
        tmp <- nllIAK3D(pars = pars , zData = zData[iThis] , XData = XData[iThis,,drop = FALSE] , vXU = NA , iU = NA , modelx = modelx , nud = nud , 
                        sdfdType_cd1 = sdfdType_cd1 , sdfdType_cxd0 = sdfdType_cxd0 , sdfdType_cxd1 = sdfdType_cxd1 , 
                        cmeOpt = cmeOpt , prodSum = prodSum , setupMats = setupMats[[i]] , parBnds = parBnds , useReml = useReml , lnTfmdData = lnTfmdData , rtnAll = TRUE , forCompLik = TRUE)	      
        XziAXz_SUM <- XziAXz_SUM + tmp$WiAW
        lndetA_SUM <- lndetA_SUM + tmp$lndetA
        n_SUM <- n_SUM + length(iThis)
        parsBTfmd <- tmp$parsBTfmd
        sigma2Vec <- tmp$sigma2Vec
        
        if(rtnAll){ 
          if(attachBigMats){
            listCkl_kl[[i]] <- tmp$A
          #        listiCkl_kl[[i]] <- lndetANDinvCb(tmp$A)$invCb
            listiCkl_kl[[i]] <- chol2inv(chol(tmp$A)) # should make code better as already solved a system in nll function.
          }else{}
          # listiAXkl_kl[[i]] <- tmp$iAW[,1:p,drop=FALSE]
          # listiAzkl_kl[[i]] <- tmp$iAW[,p+1,drop=FALSE]
        }else{}
      }
    }else{}

    if(compLikMats$compLikOptn == 1 || compLikMats$compLikOptn == 3 || compLikMats$compLikOptn == 4){
      nNonadjSubsets <- nrow(compLikMats$subsetPairsNonadj)
  
      if(nNonadjSubsets > 0 || compLikMats$compLikOptn == 4){  
        nadjSubsets <- nrow(compLikMats$subsetPairsAdj)
        XziAXz_BLOCKS <- list() # mayeb try as array in future?
        lndetA_BLOCKS <- NA * numeric(length(compLikMats$listBlocks))
        for(i in 1:length(compLikMats$listBlocks)){
          iThis <- compLikMats$listBlocks[[i]]$i
          tmp <- nllIAK3D(pars = pars , zData = zData[iThis] , XData = XData[iThis,,drop = FALSE] , vXU = NA , iU = NA , modelx = modelx , nud = nud , 
                sdfdType_cd1 = sdfdType_cd1 , sdfdType_cxd0 = sdfdType_cxd0 , sdfdType_cxd1 = sdfdType_cxd1 , 
                cmeOpt = cmeOpt , prodSum = prodSum , setupMats = setupMats[[i+nadjSubsets]] , parBnds = parBnds , useReml = useReml , lnTfmdData = lnTfmdData , rtnAll = TRUE , forCompLik = TRUE)	      
          XziAXz_BLOCKS[[i]] <- tmp$WiAW
          lndetA_BLOCKS[i] <- tmp$lndetA
          if(compLikMats$compLikOptn == 4){
            parsBTfmd <- tmp$parsBTfmd
            sigma2Vec <- tmp$sigma2Vec
          }else{}
          
          if(rtnAll){ # these ones are for a single block, not for a pair 
            if(attachBigMats){
              listCkl_kl[[i+nadjSubsets]] <- tmp$A
              #            listiCkl_kl[[i+nadjSubsets]] <- lndetANDinvCb(tmp$A)$invCb
              listiCkl_kl[[i+nadjSubsets]] <- chol2inv(chol(tmp$A))
            }else{}
            # listiAXkl_kl[[i+nadjSubsets]] <- tmp$iAW[,1:p,drop=FALSE]
            # listiAzkl_kl[[i+nadjSubsets]] <- tmp$iAW[,p+1,drop=FALSE]
          }else{}
        }

        if(compLikMats$compLikOptn == 1 || compLikMats$compLikOptn == 3){
          for(i in 1:nrow(compLikMats$subsetPairsNonadj)){
            iBlockThis <- compLikMats$subsetPairsNonadj[i,1]
            jBlockThis <- compLikMats$subsetPairsNonadj[i,2]
    
            XziAXz_SUM <- XziAXz_SUM + XziAXz_BLOCKS[[iBlockThis]] + XziAXz_BLOCKS[[jBlockThis]]
            lndetA_SUM <- lndetA_SUM + lndetA_BLOCKS[iBlockThis] + lndetA_BLOCKS[jBlockThis]
          }
        }else if(compLikMats$compLikOptn == 4){
          for(i in 1:length(compLikMats$listBlocks)){
            XziAXz_SUM <- XziAXz_SUM + XziAXz_BLOCKS[[i]]
            lndetA_SUM <- lndetA_SUM + lndetA_BLOCKS[i]
          }
          
        }else{}
      }else{}

### with all subset pairs included...
      if(compLikMats$compLikOptn == 3){
        n_SUM <- n * (length(compLikMats$listBlocks) - 1)
      }else{}

    }else{}

    if(compLikMats$compLikOptn == 1){
### effectively likelihood for nSubsets-1 copies of the data, so:
      XziAXz <- XziAXz_SUM / (length(compLikMats$listBlocks)-1)
      lndetA <- lndetA_SUM / (length(compLikMats$listBlocks)-1)
    }else{
### eidsvik doesn't average, just sums, so...
      XziAXz <- XziAXz_SUM 
      lndetA <- lndetA_SUM 
    }
        
    if(p > 0){
      betahat <- lndetANDinvCb(XziAXz[1:p,1:p,drop = FALSE] , XziAXz[1:p,p+1,drop = FALSE])
      lndetXiAX <- betahat$lndetC
      betahat <- betahat$invCb
      if(is.na(lndetXiAX) || is.infinite(lndetXiAX)){
        return(badnll)
      }else{}

    }else{
      betahat <- c()
      lndetXiAX <- c()
    }
      
    resiAres <- as.numeric(XziAXz[p+1,p+1] - 2 * t(betahat) %*% XziAXz[1:p,p+1,drop=FALSE] + t(betahat) %*% XziAXz[1:p,1:p,drop=FALSE] %*% betahat)

    if(compLikMats$compLikOptn == 2 || compLikMats$compLikOptn == 3){
      if(useReml){ stop('Error - Eidsvik was only defined for ML, not REML estimation! Work on this!') }else{}
      cxdhat <- resiAres / n_SUM 
    }else{
      if(useReml){
        cxdhat <- resiAres / (n - p)
      }else{
        cxdhat <- resiAres / n 
      }
    }
### update the values in parsBTfmd so that it makes sense when returned to user.          
    parsBTfmd$cx0 <- cxdhat * parsBTfmd$cx0
    parsBTfmd$cx1 <- cxdhat * parsBTfmd$cx1
    parsBTfmd$cd1 <- cxdhat * parsBTfmd$cd1
    parsBTfmd$cxd0 <- cxdhat * parsBTfmd$cxd0
    parsBTfmd$cxd1 <- cxdhat * parsBTfmd$cxd1
    if(cmeOpt == 1){ parsBTfmd$cme <- cxdhat * parsBTfmd$cme }else{}
                
    if(compLikMats$compLikOptn == 2 || compLikMats$compLikOptn == 3){
      lndetC <- lndetA + n_SUM * log(cxdhat)
    }else{
      lndetC <- lndetA + n * log(cxdhat)
    }
    lndetXiCX <- lndetXiAX - p * log(cxdhat)
    resiCres <- resiAres / cxdhat

    if(rtnAll){ 
        vbetahat <- cxdhat * solve(XziAXz[1:p,1:p,drop = FALSE]) 
        sigma2Vec <- cxdhat * sigma2Vec 
        if(attachBigMats){
          listiCReskl_kl <- list()
          for (i in 1:length(listCkl_kl)){
            listCkl_kl[[i]] <- cxdhat * listCkl_kl[[i]]
            listiCkl_kl[[i]] <- listiCkl_kl[[i]] / cxdhat
            # listiCReskl_kl[[i]] <- (listiAzkl_kl[[i]] - listiAXkl_kl[[i]] %*% betahat) / cxdhat
          }
        }else{}
    }else{}

    if(useReml){
        nll <- 0.5 * ((n - p) * log(2 * pi) + lndetC + lndetXiCX + resiCres) 
    }else{
        nll <- 0.5 * (n * log(2 * pi) + lndetC + resiCres) 
    }
    nll <- as.numeric(nll)

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

    if(exists("printnllTime") && printnllTime){
      print('time for nll complik evaluation was:')
      print(proc.time() - ptm)
    }else{}
    
    if(exists("parsTrace4Optim")){
      nllTmp <- nllTrace4Optim # update the final value, now fn has been successful.
      if(nll != badnll){
        nllTmp[length(nllTmp)] <- nll
        assign("nllTrace4Optim" , nllTmp , envir = .GlobalEnv)
      }else{}
    }else{}

    if(rtnAll){ 
      if(attachBigMats){
        return(list('nll' = nll , 'pars' = pars , 'parsBTfmd' = parsBTfmd , 'betahat' = betahat , 'vbetahat' = vbetahat , 'cxdhat' = cxdhat , 
                    'sigma2Vec' = sigma2Vec , 'XData' = XData , 'listCkl_kl' = listCkl_kl , 'listiCkl_kl' = listiCkl_kl))
      }else{
        return(list('nll' = nll , 'pars' = pars , 'parsBTfmd' = parsBTfmd , 'betahat' = betahat , 'vbetahat' = vbetahat , 'cxdhat' = cxdhat , 
                    'sigma2Vec' = sigma2Vec , 'XData' = XData))
      }
    }else{
    	return(nll)
    }
}


### this function not used currently. replaced by ByBlock functions.
predMatsIAK3D_CL <- function(z_muhat , XData , xMap , dIMap , iData = seq(length(z_muhat)) , lmmFit){

  print('Running predMatsIAK3D_CL function')

  oldSetC <- TRUE # new way to be checked.
  
    n <- length(z_muhat)
    p <- ncol(XData)
    nxMap <- nrow(xMap)

    diagCkhiCChk_SUM <- matrix(0 , nxMap , 1)
    CkhiCX_SUM <- matrix(0 , nxMap , p)
    CkhiCz_muhat_SUM <- matrix(0 , nxMap , 1)

    for(i in 1:nrow(lmmFit$compLikMats$subsetPairsAdj)){
      iThis <- c(lmmFit$compLikMats$listBlocks[[lmmFit$compLikMats$subsetPairsAdj[i,1]]]$i , lmmFit$compLikMats$listBlocks[[lmmFit$compLikMats$subsetPairsAdj[i,2]]]$i)

      if(length(iData) < n){ ### allows a subset of complete dataset to be used...
        iThis <- intersect(iThis , iData)
      }else{}

      if(oldSetC){
        # setupMatsMap <- setupIAK3D(xData = rbind(lmmFit$xData[iThis,,drop=FALSE] , xMap) , dIData = rbind(as.matrix(lmmFit$dIData[iThis,,drop=FALSE]) , dIMap) , nDscPts = 0)
        setupMatsMap <- setupIAK3D(xData = rbind(lmmFit$xData[iThis,,drop=FALSE] , xMap) , dIData = rbind(as.matrix(lmmFit$dIData[iThis,,drop=FALSE]) , dIMap) , nDscPts = 0 , 
                                   sdfdType_cd1 = lmmFit$sdfdType_cd1 , sdfdType_cxd0 = lmmFit$sdfdType_cxd0 , sdfdType_cxd1 = lmmFit$sdfdType_cxd1 , sdfdKnots = lmmFit$sdfdKnots)
        
        tmp <- setCIAK3D(parsBTfmd = lmmFit$parsBTfmd , modelx = lmmFit$modelx ,
                sdfdType_cd1 = lmmFit$sdfdType_cd1 , sdfdType_cxd0 = lmmFit$sdfdType_cxd0 , sdfdType_cxd1 = lmmFit$sdfdType_cxd1 ,
                cmeOpt = lmmFit$cmeOpt , setupMats = setupMatsMap)

        if(i == 1){ diagCkk <- diag(tmp$C[(length(iThis)+1):(length(iThis)+nxMap),(length(iThis)+1):(length(iThis)+nxMap),drop = FALSE]) }else{}

        ChkThis <- tmp$C[1:length(iThis),(length(iThis)+1):(length(iThis)+nxMap),drop = FALSE]
        ChThis <- tmp$C[1:length(iThis),1:length(iThis),drop = FALSE]
      }else{
        if(i == 1){ 
          # setupMatsMap <- setupIAK3D(xData = xMap , dIData = dIMap , nDscPts = 0)
          setupMatsMap <- setupIAK3D(xData = xMap , dIData = dIMap , nDscPts = 0 , 
                                     sdfdType_cd1 = lmmFit$sdfdType_cd1 , sdfdType_cxd0 = lmmFit$sdfdType_cxd0 , sdfdType_cxd1 = lmmFit$sdfdType_cxd1 , sdfdKnots = lmmFit$sdfdKnots)

          tmp <- setCIAK3D(parsBTfmd = lmmFit$parsBTfmd , modelx = lmmFit$modelx , 
                           sdfdType_cd1 = lmmFit$sdfdType_cd1 , sdfdType_cxd0 = lmmFit$sdfdType_cxd0 , sdfdType_cxd1 = lmmFit$sdfdType_cxd1 , 
                           cmeOpt = lmmFit$cmeOpt , setupMats = setupMatsMap)
          
          diagCkk <- diag(tmp$C)
          rm(setupMatsMap , tmp)
        }else{}
        
        # setupMatsMap <- setupIAK3D2(xData = lmmFit$xData[iThis,,drop=FALSE] , dIData = as.matrix(lmmFit$dIData[iThis,,drop=FALSE]) ,
        #                             xData2 = xMap , dIData2 = dIMap)
        setupMatsMap <- setupIAK3D2(xData = lmmFit$xData[iThis,,drop=FALSE] , dIData = as.matrix(lmmFit$dIData[iThis,,drop=FALSE]) ,
                                    xData2 = xMap , dIData2 = dIMap , sdfdType_cd1 = lmmFit$sdfdType_cd1 , sdfdType_cxd0 = lmmFit$sdfdType_cxd0 , sdfdType_cxd1 = lmmFit$sdfdType_cxd1 , sdfdKnots = lmmFit$sdfdKnots)
        
        ChkThis <- setCIAK3D2(parsBTfmd = lmmFit$parsBTfmd , modelx = lmmFit$modelx ,
                              sdfdType_cd1 = lmmFit$sdfdType_cd1 , sdfdType_cxd0 = lmmFit$sdfdType_cxd0 , sdfdType_cxd1 = lmmFit$sdfdType_cxd1 ,
                              cmeOpt = lmmFit$cmeOpt , setupMats = setupMatsMap)
        
        # ChThis <- tmp$C[1:length(iThis),1:length(iThis),drop = FALSE]
        ### Ch should be in the stored lmmFit object?
        print('Testing ChThis is correct')
        # setupMatsMap <- setupIAK3D(xData = rbind(lmmFit$xData[iThis,,drop=FALSE] , xMap) , dIData = rbind(as.matrix(lmmFit$dIData[iThis,,drop=FALSE]) , dIMap) , nDscPts = 0)
        setupMatsMap <- setupIAK3D(xData = rbind(lmmFit$xData[iThis,,drop=FALSE] , xMap) , dIData = rbind(as.matrix(lmmFit$dIData[iThis,,drop=FALSE]) , dIMap) , nDscPts = 0 , 
                                   sdfdType_cd1 = lmmFit$sdfdType_cd1 , sdfdType_cxd0 = lmmFit$sdfdType_cxd0 , sdfdType_cxd1 = lmmFit$sdfdType_cxd1 , sdfdKnots = lmmFit$sdfdKnots)
        
        tmp <- setCIAK3D(parsBTfmd = lmmFit$parsBTfmd , modelx = lmmFit$modelx ,
                         sdfdType_cd1 = lmmFit$sdfdType_cd1 , sdfdType_cxd0 = lmmFit$sdfdType_cxd0 , sdfdType_cxd1 = lmmFit$sdfdType_cxd1 ,
                         cmeOpt = lmmFit$cmeOpt , setupMats = setupMatsMap)
        ChThisLONG <- tmp$C[1:length(iThis),1:length(iThis),drop = FALSE]
        
        ChThis <- lmmFit$listCkl_kl[[i]]
        
        if(max(abs(ChThis - ChThisLONG)) > 1E-8){
          stop('Looks like an error in loading ChThis - ordering?')
        }
      }
      
#      tmp <- lndetANDinvCb(ChThis , ChkThis)
#      iCChkThis <- tmp$invCb
      iCChkThis <- solve(ChThis , ChkThis)
      
      diagCkhiCChk_SUM <- diagCkhiCChk_SUM + matrix(colSums(ChkThis * iCChkThis) , ncol = 1)

      CkhiCX_SUM <- CkhiCX_SUM + t(iCChkThis) %*% XData[iThis,,drop=FALSE]
      CkhiCz_muhat_SUM <- CkhiCz_muhat_SUM + t(iCChkThis) %*% z_muhat[iThis]
    }

    nadjSubsets <- nrow(lmmFit$compLikMats$subsetPairsAdj)
    diagCkhiCChk_BLOCKS <- matrix(NA , nxMap , length(lmmFit$compLikMats$listBlocks))
    CkhiCX_BLOCKS <- list()
    CkhiCz_muhat_BLOCKS <- matrix(NA , nxMap , length(lmmFit$compLikMats$listBlocks))
    for(i in 1:length(lmmFit$compLikMats$listBlocks)){
      iThis <- lmmFit$compLikMats$listBlocks[[i]]$i

      if(length(iData) < n){ ### allows a subset of complete dataset to be used...
        iThis <- intersect(iThis , iData)
      }else{}

      if(oldSetC){
        # setupMatsMap <- setupIAK3D(xData = rbind(lmmFit$xData[iThis,,drop=FALSE] , xMap) , dIData = rbind(as.matrix(lmmFit$dIData[iThis,,drop=FALSE]) , dIMap) , nDscPts = 0)
        setupMatsMap <- setupIAK3D(xData = rbind(lmmFit$xData[iThis,,drop=FALSE] , xMap) , dIData = rbind(as.matrix(lmmFit$dIData[iThis,,drop=FALSE]) , dIMap) , nDscPts = 0 , 
                                   sdfdType_cd1 = lmmFit$sdfdType_cd1 , sdfdType_cxd0 = lmmFit$sdfdType_cxd0 , sdfdType_cxd1 = lmmFit$sdfdType_cxd1 , sdfdKnots = lmmFit$sdfdKnots)
        
        tmp <- setCIAK3D(parsBTfmd = lmmFit$parsBTfmd , modelx = lmmFit$modelx ,
                sdfdType_cd1 = lmmFit$sdfdType_cd1 , sdfdType_cxd0 = lmmFit$sdfdType_cxd0 , sdfdType_cxd1 = lmmFit$sdfdType_cxd1 ,
                cmeOpt = lmmFit$cmeOpt , setupMats = setupMatsMap)

        ChkThis <- tmp$C[1:length(iThis),(length(iThis)+1):(length(iThis)+nxMap),drop = FALSE]
        ChThis <- tmp$C[1:length(iThis),1:length(iThis),drop = FALSE]
      }else{
        # setupMatsMap <- setupIAK3D2(xData = lmmFit$xData[iThis,,drop=FALSE] , dIData = as.matrix(lmmFit$dIData[iThis,,drop=FALSE]) , 
        #                             xData2 = xMap , dIData2 = dIMap)
        setupMatsMap <- setupIAK3D2(xData = lmmFit$xData[iThis,,drop=FALSE] , dIData = as.matrix(lmmFit$dIData[iThis,,drop=FALSE]) , 
                                    xData2 = xMap , dIData2 = dIMap , sdfdType_cd1 = lmmFit$sdfdType_cd1 , sdfdType_cxd0 = lmmFit$sdfdType_cxd0 , sdfdType_cxd1 = lmmFit$sdfdType_cxd1 , sdfdKnots = lmmFit$sdfdKnots)
        ChkThis <- setCIAK3D2(parsBTfmd = lmmFit$parsBTfmd , modelx = lmmFit$modelx , 
                              sdfdType_cd1 = lmmFit$sdfdType_cd1 , sdfdType_cxd0 = lmmFit$sdfdType_cxd0 , sdfdType_cxd1 = lmmFit$sdfdType_cxd1 , 
                              cmeOpt = lmmFit$cmeOpt , setupMats = setupMatsMap)
        
        print('Testing old and long way...')    
        # setupMatsMap <- setupIAK3D(xData = rbind(lmmFit$xData[iThis,,drop=FALSE] , xMap) , dIData = rbind(as.matrix(lmmFit$dIData[iThis,,drop=FALSE]) , dIMap) , nDscPts = 0)
        setupMatsMap <- setupIAK3D(xData = rbind(lmmFit$xData[iThis,,drop=FALSE] , xMap) , dIData = rbind(as.matrix(lmmFit$dIData[iThis,,drop=FALSE]) , dIMap) , nDscPts = 0 , 
                                   sdfdType_cd1 = lmmFit$sdfdType_cd1 , sdfdType_cxd0 = lmmFit$sdfdType_cxd0 , sdfdType_cxd1 = lmmFit$sdfdType_cxd1 , sdfdKnots = lmmFit$sdfdKnots)
        
        tmp <- setCIAK3D(parsBTfmd = lmmFit$parsBTfmd , modelx = lmmFit$modelx ,
                         sdfdType_cd1 = lmmFit$sdfdType_cd1 , sdfdType_cxd0 = lmmFit$sdfdType_cxd0 , sdfdType_cxd1 = lmmFit$sdfdType_cxd1 ,
                         cmeOpt = lmmFit$cmeOpt , setupMats = setupMatsMap)
        
        ChThisLONG <- tmp$C[1:length(iThis),1:length(iThis),drop = FALSE]
        
        ChThis <- lmmFit$listCkl_kl[[i]]
        
        if(max(abs(ChThis - ChThisLONG)) > 1E-8){ stop('Looking error like in ChThis.') }else{}
      }
      
#      tmp <- lndetANDinvCb(ChThis , ChkThis)
#      iCChkThis <- tmp$invCb
      iCChkThis <- solve(ChThis , ChkThis)
      
      diagCkhiCChk_BLOCKS[,i] <- matrix(colSums(ChkThis * iCChkThis), ncol = 1)

      CkhiCX_BLOCKS[[i]] <- t(iCChkThis) %*% XData[iThis,,drop=FALSE]
      CkhiCz_muhat_BLOCKS[,i] <- matrix(t(iCChkThis) %*% z_muhat[iThis] , ncol = 1)
    }

    for(i in 1:nrow(lmmFit$compLikMats$subsetPairsNonadj)){
      iBlockThis <- lmmFit$compLikMats$subsetPairsNonadj[i,1]
      jBlockThis <- lmmFit$compLikMats$subsetPairsNonadj[i,2]
    
      diagCkhiCChk_SUM <- diagCkhiCChk_SUM + diagCkhiCChk_BLOCKS[,iBlockThis] + diagCkhiCChk_BLOCKS[,jBlockThis]
      CkhiCX_SUM <- CkhiCX_SUM + CkhiCX_BLOCKS[[iBlockThis]] + CkhiCX_BLOCKS[[jBlockThis]]
      CkhiCz_muhat_SUM <- CkhiCz_muhat_SUM + CkhiCz_muhat_BLOCKS[,iBlockThis] + CkhiCz_muhat_BLOCKS[,jBlockThis]
    }

### effectively likelihood for nSubsets-1 copies of the data, so:
    diagCkhiCChk <- diagCkhiCChk_SUM / (length(lmmFit$compLikMats$listBlocks)-1)
    CkhiCX <- CkhiCX_SUM / (length(lmmFit$compLikMats$listBlocks)-1)
    CkhiCz_muhat <- CkhiCz_muhat_SUM / (length(lmmFit$compLikMats$listBlocks)-1)

    return(list('diagCkk' = diagCkk , 'diagCkhiCChk' = diagCkhiCChk , 'CkhiCX' = CkhiCX , 'CkhiCz_muhat' = CkhiCz_muhat))
}

predMatsIAK3D_CL_ByBlock <- function(z_muhat , XData , xMap , dIMap , iData = seq(length(z_muhat)) , lmmFit){ # , listCkl_kl , listiCkl_kl

  print('Running predMatsIAK3D_CL_ByBlock function')
  
### get the blocks for the prediction locations...
  nBlocks <- length(lmmFit$compLikMats$listBlocks)
  allCentroids <- matrix(NA , nBlocks , 2)
  for(i in 1:nBlocks){ allCentroids[i,] <- lmmFit$compLikMats$listBlocks[[i]]$centroid }
  DTmp <- xyDist(xMap , allCentroids)
  blockMap <- apply(DTmp , 1 , which.min)

### also pass in what block each prediction location is in...
  nBlocksMap <- max(blockMap)
  nxMap <- length(blockMap)
  
  n <- length(z_muhat)
  p <- ncol(XData)
  
  A0 <- b0 <- matrix(NA , nxMap , 1)
#  diagCkk <- diagCkhiCChk <- CkhiCX <- CkhiCz_muhat <- list()

  diagCkk <- matrix(NA , nxMap , 1)
  diagCkhiCChk <- matrix(NA , nxMap , 1)
  CkhiCX <- matrix(NA , nxMap , p)
  CkhiCz_muhat <- matrix(NA , nxMap , 1)

  for (iBlock in 1:nBlocksMap){

    iMapThis <- which(blockMap == iBlock)
    xMapThis <- xMap[iMapThis,,drop=FALSE]
    dIMapThis <- dIMap[iMapThis,,drop=FALSE]

    nxMapThis <- length(iMapThis)

    if(nxMapThis > 0){

### what are the neighbouring blocks?
      iBlocksAdj <- c()
      iTmp <- which(lmmFit$compLikMats$subsetPairsAdj[,1] == iBlock)
      if(length(iTmp) > 0){ iBlocksAdj <- c(iBlocksAdj , lmmFit$compLikMats$subsetPairsAdj[iTmp,2]) }else{}
      iTmp <- which(lmmFit$compLikMats$subsetPairsAdj[,2] == iBlock)
      if(length(iTmp) > 0){ iBlocksAdj <- c(iBlocksAdj , lmmFit$compLikMats$subsetPairsAdj[iTmp,1]) }else{}
      
      if(length(iBlocksAdj) == 0){ stop('Error - every block must have a neighbour - check this!') }else{}

### what are the data for the prediciton block?
      iDatak <- lmmFit$compLikMats$listBlocks[[iBlock]]$i
      if(length(iData) < n){ ### allows a subset of complete dataset to be used...
          iDatak <- intersect(iDatak , iData)
      }else{}
      
      
      diagCkhiCChk_SUM <- matrix(0 , nxMapThis , 1)
      CkhiCX_SUM <- matrix(0 , nxMapThis , p)
      CkhiCz_muhat_SUM <- matrix(0 , nxMapThis , 1)

      for(i in 1:length(iBlocksAdj)){
### what are the data for this neighbouring block?
        iDatal <- lmmFit$compLikMats$listBlocks[[iBlocksAdj[i]]]$i
        if(length(iData) < n){ ### allows a subset of complete dataset to be used...
          iDatal <- intersect(iDatal , iData)
        }else{}

        iThis <- c(iDatak , iDatal)

### note different order to other stuff (preds first here) for consistency with Eidsvik...
        oldSetC <- TRUE # new way to be checked.
        if(oldSetC){
          # setupMatsMap <- setupIAK3D(xData = rbind(xMapThis , lmmFit$xData[iThis,,drop=FALSE]) , dIData = rbind(dIMapThis , as.matrix(lmmFit$dIData[iThis,,drop=FALSE])) , nDscPts = 0)
          setupMatsMap <- setupIAK3D(xData = rbind(xMapThis , lmmFit$xData[iThis,,drop=FALSE]) , dIData = rbind(dIMapThis , as.matrix(lmmFit$dIData[iThis,,drop=FALSE])) , nDscPts = 0 , 
                                     sdfdType_cd1 = lmmFit$sdfdType_cd1 , sdfdType_cxd0 = lmmFit$sdfdType_cxd0 , sdfdType_cxd1 = lmmFit$sdfdType_cxd1 , sdfdKnots = lmmFit$sdfdKnots)
          
          tmp <- setCIAK3D(parsBTfmd = lmmFit$parsBTfmd , modelx = lmmFit$modelx ,
                  sdfdType_cd1 = lmmFit$sdfdType_cd1 , sdfdType_cxd0 = lmmFit$sdfdType_cxd0 , sdfdType_cxd1 = lmmFit$sdfdType_cxd1 ,
                  cmeOpt = lmmFit$cmeOpt , setupMats = setupMatsMap)

          if(i == 1){
            diagCkk[iMapThis,1] <- diag(tmp$C[1:nxMapThis,1:nxMapThis,drop = FALSE])
          }else{}
          ChkThis <- tmp$C[(nxMapThis+1):(nxMapThis+length(iThis)),1:nxMapThis,drop = FALSE]
          ChThis <- tmp$C[(nxMapThis+1):(nxMapThis+length(iThis)),(nxMapThis+1):(nxMapThis+length(iThis)),drop = FALSE]
        }else{
          if(i == 1){ 
            # setupMatsMap <- setupIAK3D(xData = xMapThis , dIData = dIMapThis , nDscPts = 0)
            setupMatsMap <- setupIAK3D(xData = xMapThis , dIData = dIMapThis , nDscPts = 0 , 
                                       sdfdType_cd1 = lmmFit$sdfdType_cd1 , sdfdType_cxd0 = lmmFit$sdfdType_cxd0 , sdfdType_cxd1 = lmmFit$sdfdType_cxd1 , sdfdKnots = lmmFit$sdfdKnots)
            
            tmp <- setCIAK3D(parsBTfmd = lmmFit$parsBTfmd , modelx = lmmFit$modelx , 
                             sdfdType_cd1 = lmmFit$sdfdType_cd1 , sdfdType_cxd0 = lmmFit$sdfdType_cxd0 , sdfdType_cxd1 = lmmFit$sdfdType_cxd1 , 
                             cmeOpt = lmmFit$cmeOpt , setupMats = setupMatsMap)
            
            diagCkk[iMapThis,1] <- diag(tmp$C) 
          }else{}
          # setupMatsMap <- setupIAK3D2(xData = lmmFit$xData[iThis,,drop=FALSE] , dIData = as.matrix(lmmFit$dIData[iThis,,drop=FALSE]) , 
          #                             xData2 = xMapThis , dIData2 = dIMapThis)
          setupMatsMap <- setupIAK3D2(xData = lmmFit$xData[iThis,,drop=FALSE] , dIData = as.matrix(lmmFit$dIData[iThis,,drop=FALSE]) , 
                                      xData2 = xMapThis , dIData2 = dIMapThis , sdfdType_cd1 = lmmFit$sdfdType_cd1 , sdfdType_cxd0 = lmmFit$sdfdType_cxd0 , sdfdType_cxd1 = lmmFit$sdfdType_cxd1 , sdfdKnots = lmmFit$sdfdKnots)
          
          ChkThis <- setCIAK3D2(parsBTfmd = lmmFit$parsBTfmd , modelx = lmmFit$modelx ,
                                sdfdType_cd1 = lmmFit$sdfdType_cd1 , sdfdType_cxd0 = lmmFit$sdfdType_cxd0 , sdfdType_cxd1 = lmmFit$sdfdType_cxd1 ,
                                cmeOpt = lmmFit$cmeOpt , setupMats = setupMatsMap)
          # ChThis <- tmp$C[(nxMapThis+1):(nxMapThis+length(iThis)),(nxMapThis+1):(nxMapThis+length(iThis)),drop = FALSE]
          
          print('testing old long way')
          # setupMatsMap <- setupIAK3D(xData = rbind(xMapThis , lmmFit$xData[iThis,,drop=FALSE]) , dIData = rbind(dIMapThis , as.matrix(lmmFit$dIData[iThis,,drop=FALSE])) , nDscPts = 0)
          setupMatsMap <- setupIAK3D(xData = rbind(xMapThis , lmmFit$xData[iThis,,drop=FALSE]) , dIData = rbind(dIMapThis , as.matrix(lmmFit$dIData[iThis,,drop=FALSE])) , nDscPts = 0 , 
                                     sdfdType_cd1 = lmmFit$sdfdType_cd1 , sdfdType_cxd0 = lmmFit$sdfdType_cxd0 , sdfdType_cxd1 = lmmFit$sdfdType_cxd1 , sdfdKnots = lmmFit$sdfdKnots)
          
          tmp <- setCIAK3D(parsBTfmd = lmmFit$parsBTfmd , modelx = lmmFit$modelx ,
                           sdfdType_cd1 = lmmFit$sdfdType_cd1 , sdfdType_cxd0 = lmmFit$sdfdType_cxd0 , sdfdType_cxd1 = lmmFit$sdfdType_cxd1 ,
                           cmeOpt = lmmFit$cmeOpt , setupMats = setupMatsMap)
          
          ChThisOLDLONG <- tmp$C[(nxMapThis+1):(nxMapThis+length(iThis)),(nxMapThis+1):(nxMapThis+length(iThis)),drop = FALSE]
          ChThis <- lmmFit$listCkl_kl[[i]]
          if(max(abs(ChThis - ChThisOLDLONG)) > 1E-8){ stop('error in chThis') }else{}
        }
        
#        tmp <- lndetANDinvCb(ChThis , ChkThis)
#        iCChkThis <- tmp$invCb
        iCChkThis <- solve(ChThis , ChkThis)
        
        diagCkhiCChkThis <- matrix(colSums(ChkThis * iCChkThis) , ncol = 1)
        
        diagCkhiCChk_SUM <- diagCkhiCChk_SUM + diagCkhiCChkThis
        CkhiCX_SUM <- CkhiCX_SUM + t(iCChkThis) %*% XData[iThis,,drop=FALSE]
        CkhiCz_muhat_SUM <- CkhiCz_muhat_SUM + t(iCChkThis) %*% z_muhat[iThis]
      }

### effectively likelihood for nadj copies of the data, so:
      diagCkhiCChk[iMapThis,1] <- matrix(diagCkhiCChk_SUM / length(iBlocksAdj) , nrow = length(iMapThis))
      CkhiCX[iMapThis,] <- matrix(CkhiCX_SUM / length(iBlocksAdj) , nrow = length(iMapThis))
      CkhiCz_muhat[iMapThis,1] <- matrix(CkhiCz_muhat_SUM / length(iBlocksAdj) , nrow = length(iMapThis))
    }else{}
  }

  return(list('diagCkk' = diagCkk , 'diagCkhiCChk' = diagCkhiCChk , 'CkhiCX' = CkhiCX , 'CkhiCz_muhat' = CkhiCz_muhat))
}

predMatsIAK3D_CLEV_ByBlock <- function(z_muhat , XData , xMap , dIMap , iData = seq(length(z_muhat)) , lmmFit){

  print('Running predMatsIAK3D_CLEV_ByBlock function')
  
### get the blocks for the prediction locations...
  nBlocks <- length(lmmFit$compLikMats$listBlocks)
  allCentroids <- matrix(NA , nBlocks , 2)
  for(i in 1:nBlocks){ allCentroids[i,] <- lmmFit$compLikMats$listBlocks[[i]]$centroid }
  DTmp <- xyDist(xMap , allCentroids)
  blockMap <- apply(DTmp , 1 , which.min)
  rm(DTmp)

### also pass in what block each prediction location is in...
  nBlocksMap <- max(blockMap)
  nxMap <- length(blockMap)
  
  n <- length(z_muhat)
  p <- ncol(XData)
  
  zk <- vk <- matrix(NA , nxMap , 1)

  for (iBlock in 1:nBlocksMap){

    iMapThis <- which(blockMap == iBlock)
    xMapThis <- xMap[iMapThis,,drop=FALSE]
    dIMapThis <- dIMap[iMapThis,,drop=FALSE]

    nxMapThis <- length(iMapThis)

    if(nxMapThis > 0){

### what are the neighbouring blocks?
      iBlocksAdj <- c()
      iTmp <- which(lmmFit$compLikMats$subsetPairsAdj[,1] == iBlock)
      if(length(iTmp) > 0){ iBlocksAdj <- c(iBlocksAdj , lmmFit$compLikMats$subsetPairsAdj[iTmp,2]) }else{}
      iTmp <- which(lmmFit$compLikMats$subsetPairsAdj[,2] == iBlock)
      if(length(iTmp) > 0){ iBlocksAdj <- c(iBlocksAdj , lmmFit$compLikMats$subsetPairsAdj[iTmp,1]) }else{}
      
      if(length(iBlocksAdj) == 0){ stop('Error - every block must have a neighbour - check this!') }else{}

### what are the data for the prediciton block?
      iDatak <- lmmFit$compLikMats$listBlocks[[iBlock]]$i
      if(length(iData) < n){ ### allows a subset of complete dataset to be used...
          iDatak <- intersect(iDatak , iData)
      }else{}

      i0 <- 1:nxMapThis
      i1 <- (nxMapThis+1):(nxMapThis+length(iDatak))

      i2List <- list()
      iDatalList <- list()
      iDatalAll <- c()
      counter <- 1
      for(i in 1:length(iBlocksAdj)){
### what are the data for this neighbouring block?
          iDatal <- lmmFit$compLikMats$listBlocks[[iBlocksAdj[i]]]$i
          if(length(iData) < n){ ### allows a subset of complete dataset to be used...
            iDatal <- intersect(iDatal , iData)
          }else{}
          iDatalList[[i]] <- iDatal
          iDatalAll <- c(iDatalAll , iDatal)
          i2List[[i]] <- (nxMapThis+length(iDatak)+counter):(nxMapThis+length(iDatak)+counter+length(iDatal)-1)
          counter <- counter+length(iDatal)
      }

      iThis <- c(iDatak , iDatalAll)
### note different order to other stuff (preds first here) for consistency with Eidsvik...
      # setupMatsMap <- setupIAK3D(xData = rbind(xMapThis , lmmFit$xData[iThis,,drop=FALSE]) , dIData = rbind(dIMapThis , as.matrix(lmmFit$dIData[iThis,,drop=FALSE])) , nDscPts = 0)
      setupMatsMap <- setupIAK3D(xData = rbind(xMapThis , lmmFit$xData[iThis,,drop=FALSE]) , dIData = rbind(dIMapThis , as.matrix(lmmFit$dIData[iThis,,drop=FALSE])) , nDscPts = 0 , 
                                 sdfdType_cd1 = lmmFit$sdfdType_cd1 , sdfdType_cxd0 = lmmFit$sdfdType_cxd0 , sdfdType_cxd1 = lmmFit$sdfdType_cxd1 , sdfdKnots = lmmFit$sdfdKnots)

      C0klAll <- setCIAK3D(parsBTfmd = lmmFit$parsBTfmd , modelx = lmmFit$modelx , 
                sdfdType_cd1 = lmmFit$sdfdType_cd1 , sdfdType_cxd0 = lmmFit$sdfdType_cxd0 , sdfdType_cxd1 = lmmFit$sdfdType_cxd1 , 
                cmeOpt = lmmFit$cmeOpt , setupMats = setupMatsMap)$C
      rm(setupMatsMap)

      C0k_0k <- C0klAll[c(i0,i1) , c(i0,i1) , drop = FALSE]

      if((!is.null(lmmFit$listCkl_kl)) & (!is.null(lmmFit$listiCkl_kl))){
        loadFromFit <- TRUE
      }else{
        loadFromFit <- FALSE
      }
      
      C0k_lList <- Q0kl_i0_i2_List <- list()
      A0_SUM <- matrix(0 , nxMapThis , nxMapThis)
      b0_SUM <- matrix(0 , nxMapThis , 1)
      Bk0_SUM <- matrix(0 , nxMapThis , nxMapThis + length(iDatak))
      for(i in 1:length(iBlocksAdj)){
          iDatal <- iDatalList[[i]]
          i2 <- i2List[[i]]

          if(loadFromFit){
            ijTmp <- c(iBlock , iBlocksAdj[i])
            iPair <- which(lmmFit$compLikMats$subsetPairsAdj[,1] == min(ijTmp) & lmmFit$compLikMats$subsetPairsAdj[,2] == max(ijTmp))
            if(length(iPair) != 1){ stop('Error - not found a unique subset pair!') }else{}
            if(iBlock > iBlocksAdj[i]){
              iOrderPair <- c(seq(length(iDatal) + 1 , length(iDatal) + length(iDatak)) , seq(1 , length(iDatal)))  
              Ckl_klThis <- lmmFit$listCkl_kl[[iPair]][iOrderPair , iOrderPair , drop=FALSE]
              iCkl_klThis <- lmmFit$listiCkl_kl[[iPair]][iOrderPair , iOrderPair , drop=FALSE]
            }else{
              Ckl_klThis <- lmmFit$listCkl_kl[[iPair]]
              iCkl_klThis <- lmmFit$listiCkl_kl[[iPair]] 
            }
          }else{
            Ckl_klThis <- C0klAll[c(i1,i2) , c(i1,i2) , drop = FALSE]
          }
          
          Ckl_0This <- C0klAll[c(i1,i2) , i0 , drop = FALSE]

          if(loadFromFit){
            iCkl_kl_Ckl_0This <- iCkl_klThis %*% Ckl_0This
          }else{
#            iCkl_kl_Ckl_0This <- lndetANDinvCb(Ckl_klThis , Ckl_0This)$invCb 
            iCkl_kl_Ckl_0This <- solve(Ckl_klThis , Ckl_0This)
          }

          C0IFklThis <- C0k_0k[i0,i0,drop=FALSE] - C0klAll[i0 , c(i1,i2) , drop = FALSE] %*% iCkl_kl_Ckl_0This
### ensure it is symmetric...          
          C0IFklThis <- 0.5 * (C0IFklThis + t(C0IFklThis))

#          iC0IFklThis <- lndetANDinvCb(C0IFklThis)$invCb
          iC0IFklThis <- chol2inv(chol(C0IFklThis))
          Q0klThis_i0_i12 <- -iC0IFklThis %*% t(iCkl_kl_Ckl_0This)

          Q0klThis_i0_i0 <- iC0IFklThis

          Q0kl_i0_i2_List[[i]] <- Q0klThis_i0_i12[,(length(i1)+1):(length(i1)+length(i2)),drop=FALSE]
          C0k_lList[[i]] <- C0klAll[c(i0,i1),i2,drop=FALSE]
          
          A0This <- iC0IFklThis
          b0This <- -Q0klThis_i0_i12 %*% z_muhat[c(iDatak , iDatal),drop=FALSE]
          Bk0This <- cbind(A0This , Q0klThis_i0_i12[,1:length(i1),drop=FALSE])

          A0_SUM <- A0_SUM + A0This
          b0_SUM <- b0_SUM + b0This
          Bk0_SUM <- Bk0_SUM + Bk0This
      }

### for Eidsvik implemented with all pairs in likelihood function, non-adj pairs assumed independent...
      if(lmmFit$compLikMats$compLikOptn == 3){
        nnonadjBlocks <- length(lmmFit$compLikMats$listBlocks) - length(iBlocksAdj) - 1
        nadjSubsets <- nrow(lmmFit$compLikMats$subsetPairsAdj) # not to be confused with above which was refering to adjacency to the prediction block
        if(loadFromFit){
          iCk_k <- lmmFit$listiCkl_kl[[nadjSubsets+iBlock]]
        }else{
#          iCk_k <- lndetANDinvCb(C0k_0k[i1,i1,drop=FALSE])$invCb # could load this from saved objects.
          iCk_k <- chol2inv(chol(C0k_0k[i1,i1,drop=FALSE])) # could load this from saved objects.
        }
        C0_kiCk_k <- C0k_0k[i0,i1,drop=FALSE] %*% iCk_k 

        C0_0_IF_k <- C0k_0k[i0,i0,drop=FALSE] - C0_kiCk_k %*% C0k_0k[i1,i0,drop=FALSE]
### to ensure C0_0_IF_k is symmetric (possibly not due to numerical errors)...
        C0_0_IF_k <- 0.5 * (C0_0_IF_k + t(C0_0_IF_k))

#        iC0_0_IF_k <- lndetANDinvCb(C0_0_IF_k)$invCb
        iC0_0_IF_k <- chol2inv(chol(C0_0_IF_k))
        iC0_0_IF_k_C0_kiCk_k <- iC0_0_IF_k %*% C0_kiCk_k
        
        A0_SUM <- A0_SUM + nnonadjBlocks * iC0_0_IF_k
        b0_SUM <- b0_SUM + nnonadjBlocks * iC0_0_IF_k_C0_kiCk_k %*% z_muhat[iDatak,drop=FALSE]

        Bk0_SUM[,1:nxMapThis] <- Bk0_SUM[,1:nxMapThis,drop=FALSE] + nnonadjBlocks * iC0_0_IF_k
        Bk0_SUM[,(nxMapThis+1):(nxMapThis+length(iDatak))] <- Bk0_SUM[,(nxMapThis+1):(nxMapThis+length(iDatak)),drop=FALSE] - nnonadjBlocks * iC0_0_IF_k_C0_kiCk_k
      }else{}

      J0_SUM <- Bk0_SUM %*% C0k_0k %*% t(Bk0_SUM)
      for(i in 1:length(iBlocksAdj)){
        i2i <- i2List[[i]]
        J0_SUM <- J0_SUM + 2 * Bk0_SUM %*% C0k_lList[[i]] %*% t(Q0kl_i0_i2_List[[i]])

        for(j in 1:length(iBlocksAdj)){
          i2j <- i2List[[j]]
          J0_SUM <- J0_SUM + Q0kl_i0_i2_List[[i]] %*% C0klAll[i2i , i2j , drop = FALSE] %*% t(Q0kl_i0_i2_List[[j]])
        }
      }

### to ensure A0_SUM and J0_SUM are symmetric (possibly not due to numerical errors)...
      A0_SUM <- 0.5 * (A0_SUM + t(A0_SUM))
      J0_SUM <- 0.5 * (J0_SUM + t(J0_SUM))

#      invA0 <- lndetANDinvCb(A0_SUM)$invCb
      invA0 <- chol2inv(chol(A0_SUM))

      zk[iMapThis,1] <- matrix(invA0 %*% b0_SUM , ncol = 1)

      vkTmp <- invA0 %*% J0_SUM %*% invA0
      
      vk[iMapThis,1] <- matrix(diag(vkTmp) , nrow = length(iMapThis))
    }else{}
  }

### remember, Xk betahat will have to be added to zk after returning.
  return(list('zk' = zk , 'vk' = vk))
}

setVoronoiBlocks <- function(x , nPerBlock = 50 , plotVor = F , vcentres = NULL){
###########################################################
### the approximate number of locations in each block
### though all blocks will probably have slightly different numbers
###########################################################
  xU <- x[!duplicated(x),]
  nU <- nrow(xU)
  
  if(is.null(vcentres)){
    nBlocks <- floor(nU / nPerBlock)
    nPerSmallBlock <- nPerBlock
    nPerBigBlock <- nPerBlock + 1
  
    nSmallBlocks <- nBlocks * (nPerSmallBlock+1) - nU
    nBigBlocks <- nBlocks - nSmallBlocks 
  }else{
    nBlocks <- nrow(vcentres)
  }
  
  xUTmp <- xU
  listBlocks <- list()
  
  if(is.null(vcentres)){
    vcentres <- matrix(NA , nBlocks , 2)
    for(i in 1:nBlocks){
      listThis <- list()      
### get the longest axis: NS or EW...      
      DNS <- max(xUTmp[,2]) - min(xUTmp[,2])
      DEW <- max(xUTmp[,1]) - min(xUTmp[,1])
      if(DNS > DEW){
          iTmp <- which.max(xUTmp[,2])
      }else{
          iTmp <- which.max(xUTmp[,1])
      }      
      
      pointNThis <- xUTmp[iTmp,,drop=FALSE]

### get the closest points...
      DThis <- xyDist(xUTmp , pointNThis)
      oThis <- order(DThis)
      if(i <= nBigBlocks){
        iThis <- oThis[1:min(nPerBigBlock , nrow(xUTmp))]
      }else if(i < nBlocks){
        iThis <- oThis[1:min(nPerSmallBlock , nrow(xUTmp))]
      }else{ # the rest...
        iThis <- oThis[1:nrow(xUTmp)]
      }    

      listThis$centroid <- matrix(colMeans(xUTmp[iThis,]) , nrow = 1)
      listBlocks[[i]] <- listThis
      vcentres[i,] <- listThis$centroid
    
      xUTmp <- xUTmp[-iThis,,drop=FALSE]
    }
  }else{
    for(i in 1:nBlocks){
      listThis <- list()
      listThis$centroid <- vcentres[i,,drop=FALSE] 
      listBlocks[[i]] <- listThis
    }  
  }

### apply above centroids in a voronoi... 
  D2Centroids <- xyDist(xU , vcentres)
  iBlock <- apply(D2Centroids , 1 , which.min)
  
  for(i in 1:nBlocks){
    iBThis <- which(iBlock == i)
    listBlocks[[i]]$xU <- xU[iBThis,,drop=FALSE]
    iThis <- which(is.element(x[,1] , listBlocks[[i]]$xU[,1]) & is.element(x[,2] , listBlocks[[i]]$xU[,2]))
    listBlocks[[i]]$i <- iThis
  }

######################################################
### use deldir to get neighbours...
######################################################
  if(plotVor){
    dev.new()
    tmp <- deldir(x = vcentres[,1] , y = vcentres[,2] , plotit = T , wl = 'te' , asp = 1) # rw = c(0 , 100 , 0 , 100) , 
    points(x , col = 'red' , pch = '.')
#    points(vcentres[,1] , vcentres[,2] , col = 'blue')
    for(i in 1:nBlocks){ text(vcentres[i,1] , vcentres[i,2] , nrow(listBlocks[[i]]$xU) , col = 'blue') }
  }else{
    tmp <- deldir(x = vcentres[,1] , y = vcentres[,2])
  }
  adjMtx <- sparseMatrix(i = c(tmp$dirsgs$ind1 , tmp$dirsgs$ind2) , j = c(tmp$dirsgs$ind2 , tmp$dirsgs$ind1) , x = 1)

#######################################################
### split the pairs into adjacent subset pairs and non-adjacent subset pairs...
#######################################################
  subsetPairsAdj <- c()
  subsetPairsNonadj <- c()
  for (i in 1:(length(listBlocks)-1)){
    for(j in (i+1):length(listBlocks)){
      if(adjMtx[i,j] == 1){
        subsetPairsAdj <- rbind(subsetPairsAdj , c(i,j))
      }else{
        subsetPairsNonadj <- rbind(subsetPairsNonadj , c(i,j))
      }
    }
  }

  return(list('listBlocks' = listBlocks , 'subsetPairsAdj' = subsetPairsAdj , 'subsetPairsNonadj' = subsetPairsNonadj))
}


setVoronoiBlocksWrap <- function(x , nPerBlock = 50 , plotVor = F , optnBalance = 1){

  compLikMats <- setVoronoiBlocks(x = x , nPerBlock = nPerBlock , plotVor = plotVor)

  if(optnBalance == 1){
    compLikMats <- balanceNeighbours(x = x , nPerBlock = nPerBlock , plotVor = plotVor , 
          listBlocks = compLikMats$listBlocks , subsetPairsAdj = compLikMats$subsetPairsAdj , subsetPairsNonadj = compLikMats$subsetPairsNonadj)
  }else if(optnBalance == 2){
    compLikMats <- balanceNeighbours2(x = x , nPerBlock = nPerBlock , plotVor = plotVor , listBlocks = compLikMats$listBlocks)
  }else{
    stop('Error - enter valid option for balancing neighbours!')
  }  
  
  return(compLikMats)
}


######################################################
### to try and improve, balance between uneven neighbours...
######################################################
balanceNeighbours <- function(x , nPerBlock = 50 , plotVor = F , listBlocks , subsetPairsAdj , subsetPairsNonadj){
  
  maxPermissDiff <- nPerBlock
  maxnUpdates <- 10
  
  carryOnEvening <- TRUE  
  adjDiffs <- NA * numeric(nrow(subsetPairsAdj))
  for(i in 1:nrow(subsetPairsAdj)){
    iThis <- subsetPairsAdj[i,1]
    jThis <- subsetPairsAdj[i,2]
    
    adjDiffs[i] <- nrow(listBlocks[[iThis]]$xU) - nrow(listBlocks[[jThis]]$xU)
  }

  vcentres <- matrix(NA , length(listBlocks) , 2)
  for(i in 1:length(listBlocks)){ 
    vcentres[i,] <- listBlocks[[i]]$centroid 
  }

  nUpdates <- 1
  while(carryOnEvening && max(abs(adjDiffs)) > maxPermissDiff){
    iTmp <- which.max(abs(adjDiffs))

    iThis <- subsetPairsAdj[iTmp,1]
    jThis <- subsetPairsAdj[iTmp,2]

    xUi <- listBlocks[[iThis]]$xU
    xUj <- listBlocks[[jThis]]$xU
       
### drop the boundary, resplit based on direction of max sep dist
    xUTmp <- rbind(xUi , xUj)
    DTmp <- xyDist(xUTmp , xUTmp)
    iDistant <- which(DTmp == max(DTmp) , arr.ind = TRUE)[1,1]

    iTmp <- order(DTmp[iDistant,])
    iTmp <- iTmp[1:ceiling(length(iTmp) / 2)]

### get which is more similar to the previous arrangement (just for some neatness and consistency, not really needed)    
    D2i <- xyDist(xUTmp[iDistant,,drop=FALSE] , listBlocks[[iThis]]$centroid)
    D2j <- xyDist(xUTmp[iDistant,,drop=FALSE] , listBlocks[[jThis]]$centroid)

    if(D2i < D2j){
      listBlocks[[iThis]]$centroid <- matrix(colMeans(xUTmp[iTmp,,drop=FALSE]) , nrow = 1)
      vcentres[iThis,] <- listBlocks[[iThis]]$centroid
      
      listBlocks[[jThis]]$centroid <- matrix(colMeans(xUTmp[-iTmp,,drop=FALSE]) , nrow = 1)
      vcentres[jThis,] <- listBlocks[[iThis]]$centroid
      
    }else{ # it more closely resembles block jThis 
      listBlocks[[jThis]]$centroid <- matrix(colMeans(xUTmp[iTmp,,drop=FALSE]) , nrow = 1)
      vcentres[jThis,] <- listBlocks[[jThis]]$centroid
    
      listBlocks[[iThis]]$centroid <- matrix(colMeans(xUTmp[-iTmp,,drop=FALSE]) , nrow = 1)
      vcentres[iThis,] <- listBlocks[[iThis]]$centroid
    }

    tmp <- setVoronoiBlocks(x = x , nPerBlock = nPerBlock , plotVor = FALSE , vcentres = vcentres)
    listBlocks <- tmp$listBlocks 
    subsetPairsAdj <- tmp$subsetPairsAdj 
    subsetPairsNonadj <- tmp$subsetPairsNonadj

    if(nUpdates > maxnUpdates){
      carryOnEvening <- FALSE
    }else{}

    nUpdates <- nUpdates + 1
    adjDiffs <- NA * numeric(nrow(subsetPairsAdj))
    for(i in 1:nrow(subsetPairsAdj)){
      iThis <- subsetPairsAdj[i,1]
      jThis <- subsetPairsAdj[i,2]
    
      adjDiffs[i] <- nrow(listBlocks[[iThis]]$xU) - nrow(listBlocks[[jThis]]$xU)
    }

  }

### to plot now, if rqd...
  if(plotVor){
    tmp <- setVoronoiBlocks(x = x , nPerBlock = nPerBlock , plotVor = TRUE , vcentres = vcentres)
  }else{}
  
  return(list('listBlocks' = listBlocks , 'subsetPairsAdj' = subsetPairsAdj , 'subsetPairsNonadj' = subsetPairsNonadj))
}

######################################################
### to try and improve, balance between uneven neighbours...
######################################################
calcnPerCluster <- function(vcentres , xU){
  ncentres <- nrow(vcentres)
  
  DTmp <- xyDist(xU,vcentres)
  
  iclstr <- apply(DTmp , 1 , which.min)

  tallyFn <- function(i , vecData){ length(which(vecData == i)) }

  nPerCluster <- sapply(seq(ncentres) , tallyFn , vecData = iclstr)

  return(nPerCluster) 
}

balanceNeighbours2 <- function(x , nPerBlock = 50 , plotVor = F , listBlocks){

### must improve within the past nItsOfNoChange4Opt
### start nReps times...
  nReps <- 5
  nItsOfNoChange4Opt <- 200
  nEpochs <- 15

  ncentres <- length(listBlocks)
  xU <- x[!duplicated(x),,drop=FALSE]

  ofVec <- NA * numeric(nReps)
  vcentresList <- list()

  for(irep in 1:nReps){
    vcentres <- matrix(NA , ncentres , 2)
    for(i in 1:ncentres){ 
      vcentres[i,] <- listBlocks[[i]]$centroid 
    }

    stddevProp <- max(c(max(xU[,1]) - min(xU[,1]) , max(xU[,2]) - min(x[,2]))) / 10 # max extent / 10
  
    nPerClusterCurrent <- calcnPerCluster(vcentres , xU = xU)
#  ofCurrent <- max(nPerClusterCurrent) - min(nPerClusterCurrent) 
    ofCurrent <- - min(nPerClusterCurrent) 
  
    print(paste0('Initial of = ' , ofCurrent))
  
    for (iEpoch in 1:nEpochs){
      print(paste0('Epoch ' , iEpoch , '...'))
      iCount <- 0
  
      while(iCount <= nItsOfNoChange4Opt){

### move the centres with the largest and smallest numbers...  
#      iMove <- c(which.min(nPerClusterCurrent) , which.max(nPerClusterCurrent))
### move all...   
#        iMove <- seq(ncentres)
### move smallest 3 and largest 3...
        iMove <- unique(c(order(nPerClusterCurrent)[1:3] , order(-nPerClusterCurrent)[1:3])) # unique in case of 6 or fewer subsets or repeatet numbers.
   
        vcentresProp <- vcentres
        vcentresProp[iMove,] <- vcentresProp[iMove,] + stddevProp * rnorm(2*length(iMove))  

        nPerClusterProp <- calcnPerCluster(vcentresProp , xU = xU)
#      ofProp <- max(nPerClusterProp) - min(nPerClusterProp)
        ofProp <-  - min(nPerClusterProp)
    
        if(ofProp < ofCurrent){
          vcentres <- vcentresProp
          nPerClusterCurrent <- nPerClusterProp
          ofCurrent <- ofProp
          iCount <- 0
          print(paste0('of reduced to ' , ofCurrent , ' in Epoch ' , iEpoch))
        }else{
          iCount <- iCount + 1
        }
      }
    
      stddevProp <- stddevProp / 2
    }
    
    vcentresList[[irep]] <- vcentres
    ofVec[irep] <- ofCurrent
  }
  
  ibest <- which.min(ofVec)
  vcentres <- vcentresList[[ibest]]
  print(paste0('Best was from rep ' , ibest , ' which gave of value of ' , ofVec[ibest]))
  
  tmp <- setVoronoiBlocks(x = x , nPerBlock = nPerBlock , plotVor = plotVor , vcentres = vcentres)
  listBlocks <- tmp$listBlocks 
  subsetPairsAdj <- tmp$subsetPairsAdj 
  subsetPairsNonadj <- tmp$subsetPairsNonadj
  
  return(list('listBlocks' = listBlocks , 'subsetPairsAdj' = subsetPairsAdj , 'subsetPairsNonadj' = subsetPairsNonadj))

}
