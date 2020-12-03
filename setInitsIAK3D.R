setParBndsIAK3D <- function(xData , dIData , setupMats , compLikMats = list('compLikOptn' = 0)){
  
  if (compLikMats$compLikOptn == 0){
    DxTmp <- setupMats$Dx[lower.tri(setupMats$Dx)]
    minDx <- min(DxTmp)
    maxDx <- max(DxTmp)
    tmp <- round(setupMats$dIU , digits = 2)
    DTmp <- xyDist(rowMeans(tmp) , rowMeans(tmp))
    DTmp <- round(DTmp , digits = 2)
    allDdIU <- DTmp[DTmp > 0]
    
    allmaxd <- setupMats$maxd
  }else{    
    setupMatsFULL <- setupIAK3D(xData , dIData , nDscPts = 0 , partSetup = TRUE)
    allmaxd <- setupMatsFULL$maxd
    
    tmp <- round(setupMatsFULL$dIU , digits = 2)
    DTmp <- xyDist(rowMeans(tmp) , rowMeans(tmp))
    DTmp <- round(DTmp , digits = 2)
    allDdIU <- DTmp[DTmp > 0]
    
    DxTmp <- xyDist(setupMatsFULL$xU , setupMatsFULL$xU)
    DxTmp <- DxTmp[lower.tri(DxTmp)]
    minDx <- min(DxTmp)
    maxDx <- max(DxTmp)
  }
  
  #################################################
  ### parBnds are used to transform parameters
  #################################################
  ### make range between 1.5 * minD and 0.5 * maxD
  parBnds <- list()
  parBnds$axBnds <- c(1.5 * minDx / 3 , 0.5 * maxDx / 3)
  parBnds$nuxBnds <- c(0.05 , 20)
  
  #    parBnds$adBnds <- c(0.05 , 1) # so that range is between 15 cm and 3 m
  parBnds$adBnds <- c(max(0.05 , min(allDdIU)/3) , 0.5 * max(allDdIU) / 3) # so that range is between at least 15 cm and 1/2 max depth diff m
  
  parBnds$tau2Bnds <- c(0.01 , 20) # tau2 = 0.01 gives range of 300 m! tau2 = 20 gives range of 15 cm. The former allows a linear function to be represented.
  parBnds$sxd1Bnds <- c(1E-8 , 1 - 1E-8) # so that at d = dmax, the sd is between 0.1 and 10 x the sd at d = 0.
  
  ### added 3/5/19, 
  #    parBnds$tau1Bnds
  
  parBnds$maxd <- allmaxd
  
  return(parBnds)
}

setInitsIAK3D <- function(xData , dIData , zData , XData , vXU , iU , modelx , nud ,  
                          sdfdType_cd1 , sdfdType_cxd0 , sdfdType_cxd1 , prodSum , 
                          cmeOpt , setupMats , parBnds , lnTfmdData , lmmFit , compLikMats = list('compLikOptn' = 0) , initC = NULL){
  
  #     if(max(c(sdfdType_cd1 , sdfdType_cxd0 , sdfdType_cxd1)) > 0){
  # ### just a quick go, check no of params and use const vec of 0.2...      
  #       nPars0 <- getnParsIAK3D(modelx = modelx , sdfdType_cd1 = sdfdType_cd1 , sdfdType_cxd0 = sdfdType_cxd0 , sdfdType_cxd1 = sdfdType_cxd1 , 
  #                               prodSum = prodSum , cmeOpt = cmeOpt , lnTfmdData = lnTfmdData)
  #       return(list('pars' = rep(0.2 , nPars0)))
  #     }else{}
  
  if (compLikMats$compLikOptn == 0){
    allKx <- setupMats$Kx
    allKd <- setupMats$Kd
    alldIU <- setupMats$dIU
  }else{
    setupMatsFULL <- setupIAK3D(xData , dIData , nDscPts = 0 , partSetup = TRUE)
    allKx <- setupMatsFULL$Kx
    allKd <- setupMatsFULL$Kd
    alldIU <- setupMatsFULL$dIU
  }
  
  #################################################
  ### fitRange (at end of this function) is used by the optimization routine and is for the transformed parameters.
  #################################################
  if(is.null(lmmFit$parsBTfmd)){ # define inits (else model already fitted)
    
    ### for inits for a newton-raphson routine to fit the variance parameters...
    if(is.null(initC)){
      betaInit <- solve(t(XData) %*% XData , t(XData) %*% zData)
    }else{
      tmp <- lndetANDinvCb(initC , cbind(zData,XData))
      betaInit <- solve(t(XData) %*% tmp$invCb[,-1,drop=FALSE] , t(XData) %*% tmp$invCb[,1,drop=FALSE])
    }
    resInit <- zData - XData %*% betaInit
    
    varResInit <- var(as.numeric(resInit))
    
    ###############################################
    ### first, fit model to profile-average data, crudely averaged at 5cm-rounded midpoints... (I think by this I meant averages for each depth in profile)
    ### this should give inits for depth-wise sum component, cd1
    ###############################################
    zTmp <- lndetANDinvCb(t(allKd) %*% allKd , t(allKd) %*% resInit)$invCb
    dIUTmp <- rowMeans(alldIU)
    dIUTmp <- round(dIUTmp * 20) / 20
    dIUTmpU <- unique(dIUTmp)
    zTmpU <- NA * numeric(length(dIUTmpU))
    for(i in 1:length(dIUTmpU)){ 
      iThis <- which(dIUTmp == dIUTmpU[i])
      zTmpU[i] <- mean(zTmp[iThis])
    }
    ### initially, make bounds a bit more estrictive than usual to make sure not starting at extremes. 
    minaTmp <- parBnds$adBnds[1] + 0.2 * (parBnds$adBnds[2] - parBnds$adBnds[1])
    maxaTmp <- parBnds$adBnds[2] - 0.4 * (parBnds$adBnds[2] - parBnds$adBnds[1])
    
    DTmp <- xyDist(dIUTmpU , dIUTmpU)
    
    tmp <- optimIt(par = c(0.5 , 0.5 * (minaTmp + maxaTmp)), fn = fLMM2 , methodOptim = c('Nelder-Mead') , c = dIUTmpU , z = zTmpU , X = matrix(1 , length(zTmpU) , 1) , D = DTmp , 
                   covModel = paste0('matern' , nud) , nSpatStructs = 1 , returnAll = F , optionML = F , verbose = F , forCompLik = FALSE , mina = minaTmp , maxa = maxaTmp)
    
    ad <- tmp$par[2]
    
    ###############################################
    ### now, fit spatial models to three depths:, A, B , C...
    ### use midpoints to assign, max of one datum for each depth per location.
    ### this should give inits for depth-wise sum component, cd1
    ###############################################
    dMidpnts <- rowMeans(dIData)
    
    ### only one from each location...
    dAB <- quantile(dMidpnts , 1/3) # the cutoff depth between A and B
    dBC <- quantile(dMidpnts , 2/3) # the cutoff depth between B and C
    dA <- quantile(dMidpnts , 1/6) # initial representative depth for A 
    dB <- quantile(dMidpnts , 0.5) # initial representative depth for B
    dC <- quantile(dMidpnts , 5/6) # initial representative depth for C
    
    iA <- iB <- iC <- c()
    #      for(i in 1:dim(setupMats$xU)[[1]]){
    for(i in 1:ncol(allKx)){
      iAThis <- which((allKx[,i] == 1) & (dMidpnts <= dAB))
      if(length(iAThis) == 0){
      }else if(length(iAThis) == 1){
        iA <- c(iA , iAThis)
      }else{
        iAThis <- iAThis[which.min(abs(dMidpnts[iAThis] - dA))]
        iA <- c(iA , iAThis)
      }
      
      iBThis <- which((allKx[,i] == 1) & (dMidpnts > dAB) & (dMidpnts <= dBC))
      if(length(iBThis) == 0){
      }else if(length(iBThis) == 1){
        iB <- c(iB , iBThis)
      }else{
        iBThis <- iBThis[which.min(abs(dMidpnts[iBThis] - dB))]
        iB <- c(iB , iBThis)
      }
      
      iCThis <- which((allKx[,i] == 1) & (dMidpnts > dBC))
      if(length(iCThis) == 0){
      }else if(length(iCThis) == 1){
        iC <- c(iC , iCThis)
      }else{
        iCThis <- iCThis[which.min(abs(dMidpnts[iCThis] - dC))]
        iC <- c(iC , iCThis)
      }
    }
    
    ### update representative depths with means...
    dMidpntsA <- dMidpnts[iA]
    xA <- xData[iA,]
    zTmpA <- resInit[iA]
    dA <- mean(dMidpntsA)
    
    dMidpntsB <- dMidpnts[iB]
    xB <- xData[iB,]
    zTmpB <- resInit[iB]
    dB <- mean(dMidpntsB)
    
    dMidpntsC <- dMidpnts[iC]
    xC <- xData[iC,]
    zTmpC <- resInit[iC]
    dC <- mean(dMidpntsC)
    
    if(modelx == 'nugget'){
      covInitxdA <- list('c0' = var(zTmpA) , 'c1' = 0 , 'a' = NA , 'nu' = NA)
      covInitxdB <- list('c0' = var(zTmpB) , 'c1' = 0 , 'a' = NA , 'nu' = NA)
      covInitxdC <- list('c0' = var(zTmpC) , 'c1' = 0 , 'a' = NA , 'nu' = NA)
    }else{
      minaTmp <- parBnds$ax[1] + 0.25 * (parBnds$ax[2] - parBnds$ax[1])
      maxaTmp <- parBnds$ax[1] + 0.75 * (parBnds$ax[2] - parBnds$ax[1])
      
      DTmp <- xyDist(xA , xA)
      tmp <- optimIt(par = c(0.5 , 0.5 * (minaTmp + maxaTmp)) , fn = fLMM2 , methodOptim = c('Nelder-Mead') , c = xA , z = zTmpA , X = matrix(1 , length(zTmpA) , 1) , D = DTmp , 
                     covModel = 'matern0.5' , nSpatStructs = 1 , returnAll = F , optionML = F , verbose = F , forCompLik = FALSE , mina = minaTmp , maxa = maxaTmp)
      tmp <- fLMM2(tmp$par , c = xA , z = zTmpA , X = matrix(1 , length(zTmpA) , 1) , D = DTmp , 
                   covModel = 'matern0.5' , nSpatStructs = 1 , returnAll = T , optionML = F , verbose = F , forCompLik = FALSE , mina = minaTmp , maxa = maxaTmp)
      covInitxdA <- list('c0' = tmp$sigma2hat * (1 - tmp$pars[1]) , 'c1' = tmp$sigma2hat * tmp$pars[1] , 'a' = tmp$pars[2] , 'nu' = 0.5)
      
      ax <- covInitxdA$a
      nux <- covInitxdA$nu
      vTmp <- (covInitxdA$c0 + covInitxdA$c1)
      if((covInitxdA$c1 / vTmp) < 0.05){ covInitxdA$c1 <- 0.05 * vTmp ; covInitxdA$c0 <- 0.95 * vTmp }else{}
      if((covInitxdA$c0 / vTmp) < 0.05){ covInitxdA$c0 <- 0.05 * vTmp ; covInitxdA$c1 <- 0.95 * vTmp }else{}
      
      DTmp <- xyDist(xB , xB)
      fitRangeTmp <- matrix(NA , 2 , 2)
      fitRangeTmp[1,] <- c(0 , 1)
      mina4B <- ax - 1
      maxa4B <- ax + 1
      tmp <- optimIt(par = c(0.5 , ax) , fn = fLMM2 , methodOptim = c('Brent') , fitRange = fitRangeTmp , vecFixedIn = c(F , T) , c = xB , z = zTmpB , X = matrix(1 , length(zTmpB) , 1) , D = DTmp , 
                     covModel = 'matern0.5' , nSpatStructs = 1 , returnAll = F , optionML = F , verbose = F , forCompLik = FALSE , mina = mina4B , maxa = maxa4B)
      
      tmp <- fLMM2(c(tmp$par , ax) , c = xB , z = zTmpB , X = matrix(1 , length(zTmpB) , 1) , D = DTmp , 
                   covModel = 'matern0.5' , nSpatStructs = 1 , returnAll = T , optionML = F , verbose = F , forCompLik = FALSE , mina = mina4B , maxa = maxa4B)
      covInitxdB <- list('c0' = tmp$sigma2hat * (1 - tmp$pars[1]) , 'c1' = tmp$sigma2hat * tmp$pars[1] , 'a' = tmp$pars[2] , 'nu' = 0.5)
      
      vTmp <- (covInitxdB$c0 + covInitxdB$c1)
      if((covInitxdB$c1 / vTmp) < 0.05){ covInitxdB$c1 <- 0.05 * vTmp ; covInitxdB$c0 <- 0.95 * vTmp }else{}
      if((covInitxdB$c0 / vTmp) < 0.05){ covInitxdB$c0 <- 0.05 * vTmp ; covInitxdB$c1 <- 0.95 * vTmp }else{}
      
      DTmp <- xyDist(xC , xC)
      fitRangeTmp <- matrix(NA , 2 , 2)
      fitRangeTmp[1,] <- c(0 , 1)
      mina4C <- ax - 1
      maxa4C <- ax + 1
      tmp <- optimIt(par = c(0.5 , ax) , fn = fLMM2 , methodOptim = c('Brent') , fitRange = fitRangeTmp , vecFixedIn = c(F , T) , c = xC , z = zTmpC , X = matrix(1 , length(zTmpC) , 1) , D = DTmp , 
                     covModel = 'matern0.5' , nSpatStructs = 1 , returnAll = F , optionML = F , verbose = F , forCompLik = FALSE , mina = mina4C , maxa = maxa4C)
      tmp <- fLMM2(c(tmp$par , ax) , c = xC , z = zTmpC , X = matrix(1 , length(zTmpC) , 1) , D = DTmp , 
                   covModel = 'matern0.5' , nSpatStructs = 1 , returnAll = T , optionML = F , verbose = F , forCompLik = FALSE , mina = mina4C , maxa = maxa4C)
      covInitxdC <- list('c0' = tmp$sigma2hat * (1 - tmp$pars[1]) , 'c1' = tmp$sigma2hat * tmp$pars[1] , 'a' = tmp$pars[2] , 'nu' = 0.5)
      
      vTmp <- (covInitxdC$c0 + covInitxdC$c1)
      if((covInitxdC$c1 / vTmp) < 0.05){ covInitxdC$c1 <- 0.05 * vTmp ; covInitxdC$c0 <- 0.95 * vTmp }else{}
      if((covInitxdC$c0 / vTmp) < 0.05){ covInitxdC$c0 <- 0.05 * vTmp ; covInitxdC$c1 <- 0.95 * vTmp }else{}
    }
    
    ### for now, if any splines being used, just use inits of ax,nux,ad...            
    if(max(c(sdfdType_cd1 , sdfdType_cxd0 , sdfdType_cxd1)) > 0){
      ### just a quick go, check no of params and use const vec of 0.2...      
      nPars0 <- getnParsIAK3D(modelx = modelx , sdfdType_cd1 = sdfdType_cd1 , sdfdType_cxd0 = sdfdType_cxd0 , sdfdType_cxd1 = sdfdType_cxd1 , 
                              prodSum = prodSum , cmeOpt = cmeOpt , lnTfmdData = lnTfmdData)
      if(modelx == 'nugget'){
        return(list('pars' = c(logitab(ad , parBnds$adBnds[1] , parBnds$adBnds[2]) , rep(0.2 , nPars0 - 1))))
      }else{
        return(list('pars' = c(logitab(ax , parBnds$axBnds[1] , parBnds$axBnds[2]) , 
                               logitab(nux , parBnds$nuxBnds[1] , parBnds$nuxBnds[2]) , 
                               logitab(ad , parBnds$adBnds[1] , parBnds$adBnds[2]) , 
                               rep(0.2 , nPars0 - 3))))
      }
    }else{}
    
    ######################################################
    ### use topsoil c0A to give cx0 + cxd0 + cd1 
    ### to allow fitting initial tau1s below, have different values of cxd0 and cd1...
    ### and topsoil c1A to give cx1 + cxd1 (assume initial equal variances for 2 components)
    ######################################################
    if(prodSum){
      if (sdfdType_cd1 != -9){
        cx0 <- covInitxdA$c0 / 6
        cxd0 <- covInitxdA$c0 / 2
        cd1 <- covInitxdA$c0 / 3
      }else{
        cx0 <- covInitxdA$c0 / 3
        cxd0 <- covInitxdA$c0 * 2 / 3
        cd1 <- 0
      }
    }else{
      cx0 <- 0
      cxd0 <- covInitxdA$c0
      cd1 <- 0
    }
    
    if(prodSum){
      cx1 <- cxd1 <- covInitxdA$c1 / 2
    }else{
      cx1 <- 0
      cxd1 <- covInitxdA$c1 
    }
    cSum <- cx0 + cx1 + cd1 + cxd0 + cxd1
    cx0 <- cx0 * varResInit / cSum
    cx1 <- cx1 * varResInit / cSum
    cd1 <- cd1 * varResInit / cSum
    cxd0 <- cxd0 * varResInit / cSum
    cxd1 <- cxd1 * varResInit / cSum
    
    ### must have at least a bit of each term...
    lbTmp <- 0.05 * cSum
    
    if(prodSum){ cx0 <- max(cx0 , lbTmp) }else{}
    if(prodSum){ cx1 <- max(cx1 , lbTmp) }else{}
    if (sdfdType_cd1 != -9){ cd1 <- max(cd1 , lbTmp) }else{}
    cxd0 <- max(cxd0 , lbTmp)
    cxd1 <- max(cxd1 , lbTmp)
    
    ### put all variance into nugget parameters if fitting nugget model...
    if (modelx == 'nugget'){
      cx0 <- cx0 + cx1
      cx1 <- 0
      cxd0 <- cxd0 + cxd1
      cxd1 <- 0
    }else{}
    
    if((!is.element(sdfdType_cd1 , c(0 , -9))) | (sdfdType_cxd0 != 0) | (sdfdType_cxd1 != 0)){
      #######################################################
      ### assume initially tau2d1 = tau2xd0 = tau2xd1 = 2
      #######################################################
      tau2d1 <- tau2xd0 <- tau2xd1 <- 2
      
      #######################################################
      ### use c0A, c1A, c0B, c1B, c0C and c1C to fit initial tau1 parameters...
      #######################################################
      ### c0A = cx0 + cd1 * (1 - tau1d1 * (1 - exp(-tau2d1 * dA))) + cxd0 * (1 - tau1xd0 * (1 - exp(-tau2xd0 * dA)))
      ### c1A = cx1 + cxd1 * (1 - tau1xd1 * (1 - exp(-tau2xd1 * dA)))
      ### c0B = cx0 + cd1 * (1 - tau1d1 * (1 - exp(-tau2d1 * dB))) + cxd0 * (1 - tau1xd0 * (1 - exp(-tau2xd0 * dB))) 
      ### c1B = cx1 + cxd1 * (1 - tau1xd1 * (1 - exp(-tau2xd1 * dB)))
      ### c0C = cx0 + cd1 * (1 - tau1d1 * (1 - exp(-tau2d1 * dC))) + cxd0 * (1 - tau1xd0 * (1 - exp(-tau2xd0 * dC))) 
      ### c1C = cx1 + cxd1 * (1 - tau1xd1 * (1 - exp(-tau2xd1 * dC)))
      ### ie ...
      ### c0A = cx0 + cxd0  + cd1 + tau1d1 * cd1 * (exp(-tau2d1 * dA) - 1) + tau1xd0 * cxd0 * (exp(-tau2xd0 * dA) - 1)
      ### c1A = cx1 + cxd1 + tau1xd1 * cxd1 * (exp(-tau2xd1 * dA) - 1)
      ### c0B = cx0 + cxd0  + cd1 + tau1d1 * cd1 * (exp(-tau2d1 * dB) - 1) + tau1xd0 * cxd0 * (exp(-tau2xd0 * dB) - 1)
      ### c1B = cx1 + cxd1 + tau1xd1 * cxd1 * (exp(-tau2xd1 * dB) - 1)
      ### c0C = cx0 + cxd0  + cd1 + tau1d1 * cd1 * (exp(-tau2d1 * dC) - 1) + tau1xd0 * cxd0 * (exp(-tau2xd0 * dC) - 1)
      ### c1C = cx1 + cxd1 + tau1xd1 * cxd1 * (exp(-tau2xd1 * dC) - 1)
      ### ie ...
      ### c0A - cx0 - cxd0 - cd1 = tau1d1 * cd1 * (exp(-tau2d1 * dA) - 1) + tau1xd0 * cxd0 * (exp(-tau2xd0 * dA) - 1)
      ### c1A - cx1 - cxd1 = tau1xd1 * cxd1 * (exp(-tau2xd1 * dA) - 1)
      ### c0B - cx0 - cxd0 - cd1 = tau1d1 * cd1 * (exp(-tau2d1 * dB) - 1) + tau1xd0 * cxd0 * (exp(-tau2xd0 * dB) - 1)
      ### c1B - cx1  cxd1 = tau1xd1 * cxd1 * (exp(-tau2xd1 * dB) - 1)
      ### c0C - cx0 - cxd0 - cd1 = tau1d1 * cd1 * (exp(-tau2d1 * dC) - 1) + tau1xd0 * cxd0 * (exp(-tau2xd0 * dC) - 1) 
      ### c1C - cx1 - cxd1 = tau1xd1 * cxd1 * (exp(-tau2xd1 * dC) - 1)
      
      zTmp <- matrix(c(sqrt(covInitxdA$c0) - sqrt(cx0) - sqrt(cxd0) - sqrt(cd1) , sqrt(covInitxdA$c1) - sqrt(cx1) - sqrt(cxd1) , 
                       sqrt(covInitxdB$c0) - sqrt(cx0) - sqrt(cxd0) - sqrt(cd1) , sqrt(covInitxdB$c1) - sqrt(cx1) - sqrt(cxd1) , 
                       sqrt(covInitxdC$c0) - sqrt(cx0) - sqrt(cxd0) - sqrt(cd1) , sqrt(covInitxdC$c1) - sqrt(cx1) - sqrt(cxd1)) , 6 , 1)
      
      ### fitted params will be tau1d1 , tau1xd0 , tau1xd1
      XTmp <- matrix(0 , 6 , 3)
      XTmp[1,1] <- sqrt(cd1) * (exp(-tau2d1 * dA) - 1)
      XTmp[1,2] <- sqrt(cxd0) * (exp(-tau2xd0 * dA) - 1)
      XTmp[2,3] <- sqrt(cxd1) * (exp(-tau2xd1 * dA) - 1)
      XTmp[3,1] <- sqrt(cd1) * (exp(-tau2d1 * dB) - 1)
      XTmp[3,2] <- sqrt(cxd0) * (exp(-tau2xd0 * dB) - 1)
      XTmp[4,3] <- cxd1 * (exp(-tau2xd1 * dB) - 1)
      XTmp[5,1] <- cd1 * (exp(-tau2d1 * dC) - 1)
      XTmp[5,2] <- cxd0 * (exp(-tau2xd0 * dC) - 1)
      XTmp[6,3] <- cxd1 * (exp(-tau2xd1 * dC) - 1)
      
      iTmpd1 <- 1
      iTmpxd0 <- 2
      iTmpxd1 <- 3
      
      iDlt <- c()
      jDlt <- c()
      
      if(is.element(sdfdType_cd1 , c(0 , -9))){
        iDlt <- c(iDlt , c())
        jDlt <- c(jDlt , 1)
        
        iTmpd1 <- c()
        iTmpxd0 <- iTmpxd0 - 1
        iTmpxd1 <- iTmpxd1 - 1
      }else{}
      
      if(sdfdType_cxd0 == 0){
        iDlt <- c(iDlt , c())
        jDlt <- c(jDlt , 2)
        
        iTmpxd0 <- c()
        iTmpxd1 <- iTmpxd1 - 1
      }else{}
      
      if((sdfdType_cd1 == 0) & (sdfdType_cxd0 == 0)){ 
        # also delete rows of XTmp and zTmp if both stat...
        iDlt <- c(iDlt , c(1,3,5))
      }else{}
      
      if(sdfdType_cxd1 == 0){
        iDlt <- c(iDlt , c(2,4,6))
        jDlt <- c(jDlt , 3)
        
        iTmpxd1 <- c()
      }else{}
      
      if(length(iDlt) > 0){
        XTmp <- XTmp[-iDlt,]
        zTmp <- zTmp[-iDlt]
      }else{}
      
      if(length(jDlt) > 0){
        XTmp <- XTmp[,-jDlt]
      }else{}
      
      parTmp <- solve(t(XTmp) %*% XTmp , t(XTmp) %*% zTmp)
      tau1d1 <- parTmp[iTmpd1]
      tau1xd0 <- parTmp[iTmpxd0]
      tau1xd1 <- parTmp[iTmpxd1]
    }else{}      
    
    if(is.element(sdfdType_cd1 , c(0 , -9))){
      sdfdPars_cd1 <- c()
      tau12d12 <- c() 
    }else if(sdfdType_cd1 == -1){
      
      tau12d1 <- contraintau(c(tau1d1,tau2d1) , parBnds = parBnds , sdfmaxdBnds = c(0.5 , 2))
      sdfdPars_cd1 <- tauTfm(parsIn = tau12d1 , parBnds = parBnds, invt = F)
      
    }else{
      stop('Not ready for this standxard deviation function yet!')
    } 
    
    if(sdfdType_cxd0 == 0){
      sdfdPars_cxd0 <- c() 
      tau12xd0 <- c()
    }else if(sdfdType_cxd0 == -1){
      
      tau12xd0 <- contraintau(c(tau1xd0,tau2xd0) , parBnds = parBnds , sdfmaxdBnds = c(0.5 , 2))
      sdfdPars_cxd0 <- tauTfm(parsIn = tau12xd0 , parBnds = parBnds, invt = F)
      
    }else{
      stop('Not ready for this standxard deviation function yet!')
    } 
    
    if(sdfdType_cxd1 == 0){
      sdfdPars_cxd1 <- c() 
      tau12xd1 <- c()
    }else if(sdfdType_cxd1 == -1){
      
      tau12xd1 <- contraintau(c(tau1xd1,tau2xd1) , parBnds = parBnds , sdfmaxdBnds = c(0.5 , 2))
      sdfdPars_cxd1 <- tauTfm(parsIn = tau12xd1 , parBnds = parBnds, invt = F)
      
    }else{
      stop('Not ready for this standxard deviation function yet!')
    } 
    
    if(cmeOpt == 1){
      lncme <- log(varResInit * 0.2) 
    }else{
      lncme <- c()
    }
    
    ###############################################
    ### one run of the nll function with lnTfmdData = F to get better initial sigma2, 
    ### only useful if lnTfmdData = T
    ###############################################
    if (lnTfmdData){
      if(sdfdType_cd1 == -9){
        lncd1ParTmp <- c()        
      }else{
        lncd1ParTmp <- log(cd1/(cxd0+cxd1))
      }
      
      if(modelx == 'matern'){
        ### beta,cxd[d=0] auto    
        parsTmp <- c(logitab(ax , parBnds$axBnds[1] , parBnds$axBnds[2]) , logitab(nux , parBnds$nuxBnds[1] , parBnds$nuxBnds[2]) , 
                     logitab(ad , parBnds$adBnds[1] , parBnds$adBnds[2]) , 
                     log(cx0/(cxd0+cxd1)) , log(cx1/(cxd0+cxd1)) , lncd1ParTmp , logitab(cxd1/(cxd0+cxd1)) ,
                     sdfdPars_cd1 , sdfdPars_cxd0 , sdfdPars_cxd1 , lncme - log(cxd0+cxd1))
      }else if(modelx == 'nugget'){
        ### beta,cxd[d=0] auto    
        parsTmp <- c(logitab(ad , parBnds$adBnds[1] , parBnds$adBnds[2]) , 
                     log((cx0 + cx1) / (cxd0 + cxd1)) , lncd1ParTmp ,  
                     sdfdPars_cd1 , sdfdPars_cxd0 , log(0.2))
      }else{}
      
      print('Running with lnTfmdData = F to get better initial variance...')
      verboseOptim <<- F
      if (compLikMats$compLikOptn == 0){
        tmp <- nllIAK3D(pars = parsTmp , zData = zData , XData = XData , vXU = vXU , iU = iU , modelx = modelx , nud = nud , 
                        sdfdType_cd1 = sdfdType_cd1 , sdfdType_cxd0 = sdfdType_cxd0 , sdfdType_cxd1 = sdfdType_cxd1 , 
                        cmeOpt = cmeOpt , prodSum = prodSum , setupMats = setupMats , parBnds = parBnds , useReml = F , lnTfmdData = F , rtnAll = T)	
      }else{
        tmp <- nllIAK3D_CL(pars = parsTmp , zData = zData , XData = XData , modelx = modelx , nud = nud , 
                           sdfdType_cd1 = sdfdType_cd1 , sdfdType_cxd0 = sdfdType_cxd0 , sdfdType_cxd1 = sdfdType_cxd1 , 
                           cmeOpt = cmeOpt , prodSum = prodSum , setupMats = setupMats , parBnds = parBnds , useReml = F , compLikMats = compLikMats , rtnAll = T)	
      }
      verboseOptim <<- T
      
      ### better estimate of cxd is cxdhat...
      rTmp <- tmp$cxdhat / (cxd0 + cxd1)
      cx0 <- cx0 * rTmp
      cx1 <- cx1 * rTmp
      cd1 <- cd1 * rTmp
      cxd0 <- cxd0 * rTmp
      cxd1 <- cxd1 * rTmp
    }else{}
    
  }else{
    ### parameters already fitted in lmmFit...
    ax <- lmmFit$parsBTfmd$ax
    nux <- lmmFit$parsBTfmd$nux
    ad <- lmmFit$parsBTfmd$ad
    nud <- lmmFit$parsBTfmd$nud
    cx0 <- lmmFit$parsBTfmd$cx0
    cx1 <- lmmFit$parsBTfmd$cx1
    cd1 <- lmmFit$parsBTfmd$cd1
    cxd0 <- lmmFit$parsBTfmd$cxd0
    cxd1 <- lmmFit$parsBTfmd$cxd1
    
    sdfdPars_cd1 <- tauTfm(parsIn = lmmFit$parsBTfmd$sdfdPars_cd1 , parBnds = parBnds , invt = F)
    sdfdPars_cxd0 <- tauTfm(parsIn = lmmFit$parsBTfmd$sdfdPars_cxd0 , parBnds = parBnds , invt = F)
    sdfdPars_cxd1 <- tauTfm(parsIn = lmmFit$parsBTfmd$sdfdPars_cxd1 , parBnds = parBnds , invt = F)
    
    if(cmeOpt == 1){
      lncme <- log(lmmFit$parsBTfmd$cme)
    }else{
      lncme <- c()
    }
  }
  
  ############################################################
  ### finally set the parameters in the pars vector...
  ############################################################
  if(sdfdType_cd1 == -9){
    lncd1ParTmp <- c()        
  }else{
    if(!lnTfmdData){
      lncd1ParTmp <- log(cd1/(cxd0+cxd1))
    }else{
      lncd1ParTmp <- log(cd1)
    }
  }
  
  if(!prodSum){
    lncx0ParTmp <- c()
    lncx1ParTmp <- c()
  }else{
    if(!lnTfmdData){
      lncx0ParTmp <- log(cx0/(cxd0+cxd1))
      lncx1ParTmp <- log(cx1/(cxd0+cxd1))
    }else{
      lncx0ParTmp <- log(cx0)
      lncx1ParTmp <- log(cx1)
    }
  }
  
  if(modelx == 'matern'){
    if(!lnTfmdData){
      ### beta,cxd[d=0] auto    
      parsStat1 <- c(logitab(ax , parBnds$axBnds[1] , parBnds$axBnds[2]) , logitab(nux , parBnds$nuxBnds[1] , parBnds$nuxBnds[2]) , 
                     logitab(ad , parBnds$adBnds[1] , parBnds$adBnds[2]) , 
                     lncx0ParTmp , lncx1ParTmp , lncd1ParTmp , logitab(cxd1/(cxd0+cxd1)))
      parsNonStat <- c(sdfdPars_cd1 , sdfdPars_cxd0 , sdfdPars_cxd1)
      parsStat2 <- c(lncme - log(cxd0+cxd1))
      parNames <- c('ax.lt' , 'nux.lt' , 'ad.lt')
      if(length(lncx0ParTmp) > 0){ parNames <- c(parNames , 'cx0OVERcxd.l') }else{}
      if(length(lncx1ParTmp) > 0){ parNames <- c(parNames , 'cx1OVERcxd.l') }else{}
      if(length(lncd1ParTmp) > 0){ parNames <- c(parNames , 'cd1OVERcxd.l') }else{}
      parNames <- c(parNames , 'cxd1OVERcxd.lt') 
      if(length(sdfdPars_cd1) > 0){ parNames <- c(parNames , c('taud1.1.tfm' , 'taud1.2.tfm')) }else{}
      if(length(sdfdPars_cxd0) > 0){ parNames <- c(parNames , c('tauxd0.1.tfm' , 'tauxd0.2.tfm')) }else{}
      if(length(sdfdPars_cxd1) > 0){ parNames <- c(parNames , c('tauxd1.1.tfm' , 'tauxd1.2.tfm')) }else{}
      if(length(lncme) > 0){ parNames <- c(parNames , 'cmeOVERcxd.l') }else{}
    }else{
      ### beta by nr
      parsStat1 <- c(logitab(ax , parBnds$axBnds[1] , parBnds$axBnds[2]) , logitab(nux , parBnds$nuxBnds[1] , parBnds$nuxBnds[2]) , 
                     logitab(ad , parBnds$adBnds[1] , parBnds$adBnds[2]) , 
                     lncx0ParTmp , lncx1ParTmp , lncd1ParTmp , log(cxd0) , log(cxd1))
      parsNonStat <- c(sdfdPars_cd1 , sdfdPars_cxd0 , sdfdPars_cxd1)
      parsStat2 <- lncme
      parNames <- c('ax.lt' , 'nux.lt' , 'ad.lt')
      if(length(lncx0ParTmp) > 0){ parNames <- c(parNames , 'cx0.l') }else{}
      if(length(lncx1ParTmp) > 0){ parNames <- c(parNames , 'cx1.l') }else{}
      if(length(lncd1ParTmp) > 0){ parNames <- c(parNames , 'cd1.l') }else{}
      if(length(cxd0) > 0){ parNames <- c(parNames , 'cxd0.l') }else{}
      if(length(cxd1) > 0){ parNames <- c(parNames , 'cxd1.l') }else{}
      if(length(sdfdPars_cd1) > 0){ parNames <- c(parNames , c('taud1.1.tfm' , 'taud1.2.tfm')) }else{}
      if(length(sdfdPars_cxd0) > 0){ parNames <- c(parNames , c('tauxd0.1.tfm' , 'tauxd0.2.tfm')) }else{}
      if(length(sdfdPars_cxd1) > 0){ parNames <- c(parNames , c('tauxd1.1.tfm' , 'tauxd1.2.tfm')) }else{}
      if(length(lncme) > 0){ parNames <- c(parNames , 'cme.l') }else{}
    }
    
  }else if(modelx == 'nugget'){
    if(!lnTfmdData){
      ### beta,cxd[d=0] auto    
      parsStat1 <- c(logitab(ad , parBnds$adBnds[1] , parBnds$adBnds[2]) , 
                     log((cx0 + cx1) / (cxd0 + cxd1)) , lncd1ParTmp)
      parsNonStat <- c(sdfdPars_cd1 , sdfdPars_cxd0)
      parsStat2 <- lncme - log(cxd0+cxd1)
      parNames <- c('ad.lt')
      if((length(cx0)+length(cx1)) > 0){ parNames <- c(parNames , 'cx0OVERcxd.l') }else{}
      if(length(lncd1ParTmp) > 0){ parNames <- c(parNames , 'cd1OVERcxd.l') }else{}
      if(length(sdfdPars_cd1) > 0){ parNames <- c(parNames , c('taud1.1.tfm' , 'taud1.2.tfm')) }else{}
      if(length(sdfdPars_cxd0) > 0){ parNames <- c(parNames , c('tauxd0.1.tfm' , 'tauxd0.2.tfm')) }else{}
      if(length(lncme) > 0){ parNames <- c(parNames , 'cmeOVERcxd.l') }else{}
    }else{
      ### beta by nr
      parsStat1 <- c(logitab(ad , parBnds$adBnds[1] , parBnds$adBnds[2]) , 
                     log(cx0 + cx1) , lncd1ParTmp , log(cxd0 + cxd1))
      parsNonStat <- c(sdfdPars_cd1 , sdfdPars_cxd0)
      parsStat2 <- lncme
      parNames <- c('ad.lt')
      if(length(cx0) > 0){ parNames <- c(parNames , 'cx0.l') }else{}
      if(length(lncd1ParTmp) > 0){ parNames <- c(parNames , 'cd1.l') }else{}
      if(length(cxd0) > 0){ parNames <- c(parNames , 'cxd0.l') }else{}
      
      if(length(sdfdPars_cd1) > 0){ parNames <- c(parNames , c('taud1.1.tfm' , 'taud1.2.tfm')) }else{}
      if(length(sdfdPars_cxd0) > 0){ parNames <- c(parNames , c('tauxd0.1.tfm' , 'tauxd0.2.tfm')) }else{}
      if(length(lncme) > 0){ parNames <- c(parNames , 'cme.l') }else{}
    }
    
  }
  
  ### put the above-created subvectors together. 
  ### iStat, iTau1 and iTau2 are the indices of these components in the full pars vector.    
  pars <- c(parsStat1 , parsNonStat , parsStat2)
  parsStat <- c(parsStat1 , parsStat2)
  
  iStat <- seq(length(parsStat1))
  if (length(parsStat2) > 0){ iStat <- c(iStat , length(pars)) }else{} # add ME to the stat pars...
  
  if(length(parsNonStat) > 0){
    iTau1 <- length(parsStat1) + seq(1 , length(parsNonStat) , 2)
    iTau2 <- iTau1 + 1
  }else{
    iTau1 <- iTau2 <- c()        
  }
  
  return(list('pars' = pars , 'parsStat' = parsStat , 'iStat' = iStat , 'iTau1' = iTau1 , 'iTau2' = iTau2))
}

getParNamesIAK3D <- function(modelx , sdfdType_cd1 , sdfdType_cxd0 , sdfdType_cxd1 , prodSum , cmeOpt , lnTfmdData){
  
#  if(max(c(sdfdType_cd1 , sdfdType_cxd0 , sdfdType_cxd1)) > 0){ stop('Update getParNamesIAK3D function for sdfdType > 0!') }else{}
  
  if(modelx == 'matern'){
    if(!lnTfmdData){
      ### beta,cxd[d=0] auto    
      parNames <- c('ax.lt' , 'nux.lt' , 'ad.lt')
      if(prodSum){ parNames <- c(parNames , 'cx0OVERcxd.l') }else{}
      if(prodSum){ parNames <- c(parNames , 'cx1OVERcxd.l') }else{}
      if(prodSum & (sdfdType_cd1 != -9)){ parNames <- c(parNames , 'cd1OVERcxd.l') }else{}
      parNames <- c(parNames , 'cxd1OVERcxd.lt')
      if(sdfdType_cd1 == -1){ parNames <- c(parNames , c('taud1.1.tfm' , 'taud1.2.tfm')) }else{}
      if(sdfdType_cd1 > 0){ parNames <- c(parNames , paste0('taud1.' , seq(sdfdType_cd1))) }else{}
      if(sdfdType_cxd0 == -1){ parNames <- c(parNames , c('tauxd0.1.tfm' , 'tauxd0.2.tfm')) }else{}
      if(sdfdType_cxd0 > 0){ parNames <- c(parNames , paste0('tauxd0.' , seq(sdfdType_cxd0))) }else{}
      if(sdfdType_cxd1 == -1){ parNames <- c(parNames , c('tauxd1.1.tfm' , 'tauxd1.2.tfm')) }else{}
      if(sdfdType_cxd1 > 0){ parNames <- c(parNames , paste0('tauxd1.' , seq(sdfdType_cxd1))) }else{}
      if(cmeOpt == 1){ parNames <- c(parNames , 'cmeOVERcxd.l') }else{}
    }else{
      ### beta by nr
      parNames <- c('ax.lt' , 'nux.lt' , 'ad.lt')
      if(prodSum){ parNames <- c(parNames , 'cx0.l') }else{}
      if(prodSum){ parNames <- c(parNames , 'cx1.l') }else{}
      if(prodSum & (sdfdType_cd1 != -9)){ parNames <- c(parNames , 'cd1.l') }else{}
      if(sdfdType_cxd0 != -9){ parNames <- c(parNames , 'cxd0.l') }else{}
      if(sdfdType_cxd1 != -9){ parNames <- c(parNames , 'cxd1.l') }else{}
      if(sdfdType_cd1 == -1){ parNames <- c(parNames , c('taud1.1.tfm' , 'taud1.2.tfm')) }else{}
      if(sdfdType_cd1 > 0){ parNames <- c(parNames , paste0('taud1.' , seq(sdfdType_cd1))) }else{}
      if(sdfdType_cxd0 == -1){ parNames <- c(parNames , c('tauxd0.1.tfm' , 'tauxd0.2.tfm')) }else{}
      if(sdfdType_cxd0 > 0){ parNames <- c(parNames , paste0('tauxd0.' , seq(sdfdType_cxd0))) }else{}
      if(sdfdType_cxd1 == -1){ parNames <- c(parNames , c('tauxd1.1.tfm' , 'tauxd1.2.tfm')) }else{}
      if(sdfdType_cxd1 > 0){ parNames <- c(parNames , paste0('tauxd1.' , seq(sdfdType_cxd1))) }else{}
      if(cmeOpt == 1){ parNames <- c(parNames , 'cme.l') }else{}
    }
    
  }else if(modelx == 'nugget'){
    if(!lnTfmdData){
      ### beta,cxd[d=0] auto    
      parNames <- c('ad.lt')
      if(prodSum){ parNames <- c(parNames , 'cx0OVERcxd.l') }else{}
      if(prodSum & (sdfdType_cd1 != -9)){ parNames <- c(parNames , 'cd1OVERcxd.l') }else{}
      if(sdfdType_cd1 == -1){ parNames <- c(parNames , c('taud1.1.tfm' , 'taud1.2.tfm')) }else{}
      if(sdfdType_cd1 > 0){ parNames <- c(parNames , paste0('taud1.' , seq(sdfdType_cd1))) }else{}
      if(sdfdType_cxd0 == -1){ parNames <- c(parNames , c('tauxd0.1.tfm' , 'tauxd0.2.tfm')) }else{}
      if(sdfdType_cxd0 > 0){ parNames <- c(parNames , paste0('tauxd0.' , seq(sdfdType_cxd0))) }else{}
      if(cmeOpt == 1){ parNames <- c(parNames , 'cmeOVERcxd.l') }else{}
    }else{
      ### beta by nr
      parNames <- c('ad.lt' , 'cx0.l')
      if(prodSum){ parNames <- c(parNames , 'cx0.l') }else{}
      if(prodSum & (sdfdType_cd1 != -9)){ parNames <- c(parNames , 'cd1.l') }else{}
      parNames <- c(parNames , 'cxd0.l')
      if(sdfdType_cd1 == -1){ parNames <- c(parNames , c('taud1.1.tfm' , 'taud1.2.tfm')) }else{}
      if(sdfdType_cd1 > 0){ parNames <- c(parNames , paste0('taud1.' , seq(sdfdType_cd1))) }else{}
      if(sdfdType_cxd0 == -1){ parNames <- c(parNames , c('tauxd0.1.tfm' , 'tauxd0.2.tfm')) }else{}
      if(sdfdType_cxd0 > 0){ parNames <- c(parNames , paste0('tauxd0.' , seq(sdfdType_cxd0))) }else{}
      if(cmeOpt == 1){ parNames <- c(parNames , 'cme.l') }else{}
    }
    
  }
  
  return(parNames)
}

setFitRangeIAK3D <- function(pars , modelx , sdfdType_cd1 , sdfdType_cxd0 , sdfdType_cxd1 , prodSum , cmeOpt , lnTfmdData){
  
  ### logitab tfms are given the range -20 , 20
  ### log tfms are given the range -20 , Inf
  
  if(modelx == 'matern'){
    if(!lnTfmdData){
      
      ### beta,cxd[d=0] auto    
      ### parsStat1 <- logitab[c(ax , nux , ad)] , lncx0ParTmp , lncx1ParTmp , lncd1ParTmp , logitab(cxd1/(cxd0+cxd1))
      ### parsNonStat <- c(sdfdPars_cd1 , sdfdPars_cxd0 , sdfdPars_cxd1)
      ### parsStat2 <- c(lncme - log(cxd0+cxd1))
      
      fitRange <- matrix(NA , length(pars) , 2) ; fitRange[,1] <- -Inf ; fitRange[,2] <- Inf ; 
      inext <- 1
      fitRange[inext,] <- c(-20 , 20) ; inext <- inext + 1
      fitRange[inext,] <- c(-20 , 20) ; inext <- inext + 1
      fitRange[inext,] <- c(-20 , 20) ; inext <- inext + 1
      
      if(prodSum){
        fitRange[inext:(inext+1),1] <- -20 ; inext <- inext + 2 # the variance ratios for cx0 and cx1
      }else{}
      if(sdfdType_cd1 != -9){
        fitRange[inext,1] <- -20 ; inext <- inext + 1 # the variance ratio for cd1
      }else{}
      fitRange[inext,] <- c(-20 , 20) ; inext <- inext + 1 # the 1 variance logit, (cxd1 / (cxd0 + cxd1)) 
      
      if(sdfdType_cd1 == -1){ 
        fitRange[inext,1] <- -20 ; inext <- inext + 1 # log of sd fn at dmax
        fitRange[inext,] <- c(-20 , 20) ; inext <- inext + 1 # logit of tau2
      }else if(sdfdType_cd1 > 0){
        # inext <- inext + sdfdType_cd1 # for spline, leave unbounded (checks in fn)
        fitRange[inext + sdfdType_cd1 - 1,] <- c(0 , Inf)
        inext <- inext + sdfdType_cd1 # for spline, bound final param below by 0, others leave unbounded (checks in fn)
      }else{}        
      if(sdfdType_cxd0 == -1){ 
        fitRange[inext,1] <- -20 ; inext <- inext + 1 # log of sd fn at dmax
        fitRange[inext,] <- c(-20 , 20) ; inext <- inext + 1 # logit of tau2
      }else if(sdfdType_cxd0 > 0){
        # inext <- inext + sdfdType_cxd0 # for spline, leave unbounded (checks in fn)
        fitRange[inext + sdfdType_cxd0 - 1,] <- c(0 , Inf)
        inext <- inext + sdfdType_cxd0 # for spline, bound final param below by 0, others leave unbounded (checks in fn)
      }else{}        
      if(sdfdType_cxd1 == -1){ 
        fitRange[inext,1] <- -20 ; inext <- inext + 1 # log of sd fn at dmax
        fitRange[inext,] <- c(-20 , 20) ; inext <- inext + 1 # logit of tau2
      }else if(sdfdType_cxd1 > 0){
        # inext <- inext + sdfdType_cxd1 # for spline, leave unbounded (checks in fn)
        fitRange[inext + sdfdType_cxd1 - 1,] <- c(0 , Inf)
        inext <- inext + sdfdType_cxd1 # for spline, bound final param below by 0, others leave unbounded (checks in fn)
      }else{}        
      if(cmeOpt == 1){
        fitRange[inext,1] <- -20 ; inext <- inext + 1
      }else{}
      
    }else{
      
      ### beta by nr
      ### beta by nr
      ### parsStat1 <- logitab[ax , nux , ad] , lncx0ParTmp , lncx1ParTmp , lncd1ParTmp , log(cxd0) , log(cxd1)
      ### parsNonStat <- c(sdfdPars_cd1 , sdfdPars_cxd0 , sdfdPars_cxd1)
      ### parsStat2 <- lncme
      
      fitRange <- matrix(NA , length(pars) , 2) ; fitRange[,1] <- -Inf ; fitRange[,2] <- Inf ; 
      inext <- 1
      fitRange[inext,] <- c(-20 , 20) ; inext <- inext + 1
      fitRange[inext,] <- c(-20 , 20) ; inext <- inext + 1
      fitRange[inext,] <- c(-20 , 20) ; inext <- inext + 1
      if(prodSum){
        fitRange[inext:(inext+1),1] <- -20 ; inext <- inext + 2 # for cx0 and cx1 variances.
      }else{}
      
      if(prodSum & (sdfdType_cd1 != -9)){
        fitRange[inext,1] <- -20 ; inext <- inext + 1 # the cd1 variance.
      }else{}
      
      fitRange[inext:(inext+1),1] <- -20 ; inext <- inext + 2 # for the cxd0 and cxd1 variances.
      
      if(sdfdType_cd1 == -1){ 
        fitRange[inext,1] <- -20 ; inext <- inext + 1 # log of sd fn at dmax
        fitRange[inext,] <- c(-20 , 20) ; inext <- inext + 1 # logit of tau2
      }else if(sdfdType_cd1 > 0){
        # inext <- inext + sdfdType_cd1 # for spline, leave unbounded (checks in fn)
        fitRange[inext + sdfdType_cd1 - 1,] <- c(0 , Inf)
        inext <- inext + sdfdType_cd1 # for spline, bound final param below by 0, others leave unbounded (checks in fn)
      }else{}        
      if(sdfdType_cxd0 == -1){ 
        fitRange[inext,1] <- -20 ; inext <- inext + 1 # log of sd fn at dmax
        fitRange[inext,] <- c(-20 , 20) ; inext <- inext + 1 # logit of tau2
      }else if(sdfdType_cxd0 > 0){
        # inext <- inext + sdfdType_cxd0 # for spline, leave unbounded (checks in fn)
        fitRange[inext + sdfdType_cxd0 - 1,] <- c(0 , Inf)
        inext <- inext + sdfdType_cxd0 # for spline, bound final param below by 0, others leave unbounded (checks in fn)
      }else{}        
      if(sdfdType_cxd1 == -1){ 
        fitRange[inext,1] <- -20 ; inext <- inext + 1 # log of sd fn at dmax
        fitRange[inext,] <- c(-20 , 20) ; inext <- inext + 1 # logit of tau2
      }else if(sdfdType_cxd1 > 0){
        # inext <- inext + sdfdType_cxd1 # for spline, leave unbounded (checks in fn)
        fitRange[inext + sdfdType_cxd1 - 1,] <- c(0 , Inf)
        inext <- inext + sdfdType_cxd1 # for spline, bound final param below by 0, others leave unbounded (checks in fn)
      }else{}        
      if(cmeOpt == 1){
        fitRange[inext,1] <- -20 ; inext <- inext + 1
      }else{}
    }
    
  }else if(modelx == 'nugget'){
    
    if(!lnTfmdData){
      ### beta,cxd[d=0] auto    
      ### parsStat1 <- logitab[ad] , log((cx0 + cx1) / (cxd0 + cxd1)) , lncd1ParTmp
      ### parsNonStat <- c(sdfdPars_cd1 , sdfdPars_cxd0)
      ### parsStat2 <- lncme - log(cxd0+cxd1)
      fitRange <- matrix(NA , length(pars) , 2) ; fitRange[,1] <- -Inf ; fitRange[,2] <- Inf ; 
      inext <- 1
      fitRange[inext,] <- c(-20 , 20) ; inext <- inext + 1
      
      if(prodSum){
        fitRange[inext,1] <- -20 ; inext <- inext + 1 # the variance ratio for cx0 
      }else{}
      if(prodSum & (sdfdType_cd1 != -9)){
        fitRange[inext,1] <- -20 ; inext <- inext + 1 # the variance ratio for cd1
      }else{}
      
      if(sdfdType_cd1 == -1){ 
        fitRange[inext,1] <- -20 ; inext <- inext + 1 # log of sd fn at dmax
        fitRange[inext,] <- c(-20 , 20) ; inext <- inext + 1 # logit of tau2
      }else if(sdfdType_cd1 > 0){
        # inext <- inext + sdfdType_cd1 # for spline, leave unbounded (checks in fn)
        fitRange[inext + sdfdType_cd1 - 1,] <- c(0 , Inf)
        inext <- inext + sdfdType_cd1 # for spline, bound final param below by 0, others leave unbounded (checks in fn)
      }else{}        
      if(sdfdType_cxd0 == -1){ 
        fitRange[inext,1] <- -20 ; inext <- inext + 1 # log of sd fn at dmax
        fitRange[inext,] <- c(-20 , 20) ; inext <- inext + 1 # logit of tau2
      }else if(sdfdType_cxd0 > 0){
        # inext <- inext + sdfdType_cxd0 # for spline, leave unbounded (checks in fn)
        fitRange[inext + sdfdType_cxd0 - 1,] <- c(0 , Inf)
        inext <- inext + sdfdType_cxd0 # for spline, bound final param below by 0, others leave unbounded (checks in fn)
      }else{}        
      if(cmeOpt == 1){
        fitRange[inext,1] <- -20 ; inext <- inext + 1
      }else{}
      
    }else{
      ### beta by nr
      ### parsStat1 <- logitab[ad] , log(cx0 + cx1) , lncd1ParTmp , log(cxd0 + cxd1)
      ### parsNonStat <- c(sdfdPars_cd1 , sdfdPars_cxd0)
      ### parsStat2 <- lncme
      
      fitRange <- matrix(NA , length(pars) , 2) ; fitRange[,1] <- -Inf ; fitRange[,2] <- Inf ; 
      inext <- 1
      fitRange[inext,] <- c(-20 , 20) ; inext <- inext + 1
      
      if(prodSum){
        fitRange[inext,1] <- -20 ; inext <- inext + 1 # the variance for cx0 
      }else{}
      if(prodSum & (sdfdType_cd1 != -9)){
        fitRange[inext,1] <- -20 ; inext <- inext + 1 # the variance for cd1
      }else{}
      fitRange[inext,1] <- -20 ; inext <- inext + 1 # the variance for cxd0
      
      if(sdfdType_cd1 == -1){ 
        fitRange[inext,1] <- -20 ; inext <- inext + 1 # log of sd fn at dmax
        fitRange[inext,] <- c(-20 , 20) ; inext <- inext + 1 # logit of tau2
      }else if(sdfdType_cd1 > 0){
        # inext <- inext + sdfdType_cd1 # for spline, leave unbounded (checks in fn)
        fitRange[inext + sdfdType_cd1 - 1,] <- c(0 , Inf)
        inext <- inext + sdfdType_cd1 # for spline, bound final param below by 0, others leave unbounded (checks in fn)
      }else{}   
      if(sdfdType_cxd0 == -1){ 
        fitRange[inext,1] <- -20 ; inext <- inext + 1 # log of sd fn at dmax
        fitRange[inext,] <- c(-20 , 20) ; inext <- inext + 1 # logit of tau2
      }else if(sdfdType_cxd0 > 0){
        # inext <- inext + sdfdType_cxd0 # for spline, leave unbounded (checks in fn)
        fitRange[inext + sdfdType_cxd1 - 1,] <- c(0 , Inf)
        inext <- inext + sdfdType_cxd1 # for spline, bound final param below by 0, others leave unbounded (checks in fn)
      }else{}        
      if(cmeOpt == 1){
        fitRange[inext,1] <- -20 ; inext <- inext + 1
      }else{}
    }
    
  }
  
  ### tidy up, in case just using function to find out how many parameters to be fitted numerically...
  fitRange <- fitRange[1:(inext-1),,drop=FALSE]
  
  return(fitRange)
}

getnParsIAK3D <- function(modelx , sdfdType_cd1 , sdfdType_cxd0 , sdfdType_cxd1 , 
                          prodSum , cmeOpt , lnTfmdData){
  
  # ax,nux,ad,
  # cx0.,cx1.,cd1.,[cxd0.,]sxd1.,
  # taud.1,taud.2,
  # tauxd0.1,tauxd0.2,
  # tauxd1.1,tauxd1.2,
  # cme.
  
  tmp <- setFitRangeIAK3D(pars = NA * numeric(25) , modelx , sdfdType_cd1 , sdfdType_cxd0 , sdfdType_cxd1 , prodSum , cmeOpt , lnTfmdData)
  
  return(nrow(tmp))
}
