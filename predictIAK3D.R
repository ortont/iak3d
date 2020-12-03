predictIAK3D <- function(xMap , dIMap , covsMap , lmmFit , rqrBTfmdPreds = TRUE , constrainX4Pred = FALSE){

########################################################
### if xMap or dIMap were dataframes, convert to matrices here.
### and make sure all are numeric...
########################################################
    if(!is.matrix(xMap)){
        xMap <- as.matrix(xMap)
    }else{}
    if(!is.matrix(dIMap)){
        dIMap <- as.matrix(dIMap)
    }else{}

    xMapCopy <- matrix(NA , nrow(xMap) , ncol(xMap))
    for (i in 1:ncol(xMap)){ xMapCopy[,i] <- as.numeric(xMap[,i]) }
    xMap <- xMapCopy
    remove(xMapCopy)
    
    dIMapCopy <- matrix(NA , nrow(dIMap) , ncol(dIMap))
    for (i in 1:ncol(dIMap)){ dIMapCopy[,i] <- as.numeric(dIMap[,i]) }
    dIMap <- dIMapCopy
    remove(dIMapCopy)
    
    ##############################################
    ### make maps for depth intervals dIMap...
    ### predict for all depths for all data locations 
    ##############################################
    dIMap <- round(dIMap , digits = 2)
    dIMap <- matrix(dIMap , ncol = 2)
    ndIMap <- nrow(dIMap)
    
    ############################################################
    ### remove missing covariates after defining XMap, 
    ### as perhaps XMap will not need all covariates...
    ############################################################
    ndimTmp <- ncol(lmmFit$xData)
    if(is.null(ndimTmp)){ ndimTmp <- 1 }else{}
    
    xMap <- matrix(xMap , ncol = ndimTmp)
    nxMap <- nrow(xMap)
    
    ### up to 5000 locations at a time?...
    #    nxPerBatch <- 5000
    nxPerBatch <- 2000
    nBatches <- ceiling(nxMap / nxPerBatch)
    
    betahat <- lmmFit$betahat
    vbetahat <- lmmFit$vbetahat
    nData <- length(lmmFit$zData)
    p <- ncol(lmmFit$XData)
    pU <- length(lmmFit$iU)
    
    if(constrainX4Pred){
      XLims4Pred <- lmmFit$XLims
    }else{
      XLims4Pred <- matrix(NA , 2 , p)
      XLims4Pred[1,] <- -Inf
      XLims4Pred[2,] <- Inf
    }

    if(lmmFit$compLikMats$compLikOptn == 0){
### make or load C and delete from lmmFit...    
      if(is.element('C' , names(lmmFit))){
        C <- lmmFit$C
        lmmFit$C <- NULL
      }else{
        print('Before making C for data in predictIAK3D:')
        gc() ; # try to make sure garbege collected before making C (as runs out of mem making C in shiny app)
        print(mem_used())
        if(is.element('setupMats' , names(lmmFit))){
          C <- setCIAK3D(parsBTfmd = lmmFit$parsBTfmd , modelx = lmmFit$modelx , 
                         sdfdType_cd1 = lmmFit$sdfdType_cd1 , sdfdType_cxd0 = lmmFit$sdfdType_cxd0 , sdfdType_cxd1 = lmmFit$sdfdType_cxd1 , 
                         cmeOpt = lmmFit$cmeOpt , setupMats = lmmFit$setupMats)$C
        }else{
          C <- setCIAK3D(parsBTfmd = lmmFit$parsBTfmd , modelx = lmmFit$modelx , 
                         sdfdType_cd1 = lmmFit$sdfdType_cd1 , sdfdType_cxd0 = lmmFit$sdfdType_cxd0 , sdfdType_cxd1 = lmmFit$sdfdType_cxd1 , 
                         cmeOpt = lmmFit$cmeOpt , setupMats = list('xData' = lmmFit$xData , 'dIData' = lmmFit$dIData))$C
        }
        print('Made C for data in predictIAK3D. Now:')
        print(mem_used())
      }
    }else{}
    
    if(!lmmFit$lnTfmdData){
      muhat <- lmmFit$XData %*% betahat 
    }else{
      muhat <- setmuIAK3D(XData = lmmFit$XData , vXU = lmmFit$vXU , iU = lmmFit$iU , beta = betahat , diagC = diag(C) , sigma2Vec = lmmFit$sigma2Vec) 
    }
    
    if(lmmFit$compLikMats$compLikOptn == 0){
      if(is.element('iCX'  , names(lmmFit)) & is.element('iCz_muhat'  , names(lmmFit)) & is.element('iC'  , names(lmmFit))){
        iCX <- lmmFit$iCX
        iCz_muhat <- lmmFit$iCz_muhat
        iC <- lmmFit$iC
      }else{
        C <- chol2inv(chol(C))
        iC <- C # maybe better for memory to make as C then copy to iC then remove? 
        remove(C)
        print('Inverted C for data in predictIAK3D. Now:')
        print(mem_used())

        iCXResTmp <- matrix(iC %*% cbind(lmmFit$XData , lmmFit$zData - muhat) , ncol = p+1)
        iCX <- iCXResTmp[,1:p,drop=FALSE]
        iCz_muhat <- iCXResTmp[,p+1,drop=FALSE]
        remove(iCXResTmp)
      }
    }else{}
    
    print(paste0('Predicting for ' , ndIMap , ' depths and ' , nBatches , ' batches...'))
    ptm <- proc.time()
    
    zMap <- vMap <- matrix(NA , ndIMap , nxMap)
    for (i in 1:ndIMap){
      for(iBatch in 1:nBatches){
        iThis <- ((iBatch - 1) * nxPerBatch + 1):(min(iBatch * nxPerBatch , nxMap))
        nxMapThis <- length(iThis)
        
        dIMapThis <- kronecker(matrix(dIMap[i,,drop=FALSE] , 1 , 2) , matrix(1 , nxMapThis , 1))
        
        if(identical(lmmFit$modelX$type , 'gam2')){
          tmp <- makeXvX_gam2(covData = covsMap[iThis,,drop = FALSE] , dIData = dIMapThis , listfefdKnots = lmmFit$modelX$listfefdKnots , incInts = lmmFit$modelX$incInts , intMthd = lmmFit$modelX$intMthd , colnamesXcns = lmmFit$modelX$colnamesX , nDiscPts = 1000 , lnTfmdData = lmmFit$lnTfmdData)
        }else{
          tmp <- makeXvX(covData = covsMap[iThis,,drop = FALSE] , dIData = dIMapThis , modelX = lmmFit$modelX , allKnotsd = lmmFit$allKnotsd , iU = lmmFit$iU , nDiscPts = 10 , lnTfmdData = lmmFit$lnTfmdData , XLims = XLims4Pred)
        }
        XMapThis <- as.matrix(tmp$X)
        vXUMapThis <- as.matrix(tmp$vXU)
        
        iOKInThis <- which(!is.na(rowMeans(XMapThis)))
        nxMapThis <- length(iOKInThis)
        
        zMapThis <- vMapThis <- NA * numeric(length(iThis))
        if(nxMapThis > 0){
          XMapThis <- XMapThis[iOKInThis,  , drop = FALSE]
          
          if(lmmFit$lnTfmdData){
            iTmp <- kronecker(((iOKInThis - 1) * pU) , matrix(1 , pU , 1)) + kronecker(matrix(1 , nxMapThis , 1) , seq(pU))
            vXUMapThis <- vXUMapThis[iTmp, , drop = FALSE]
          }else{}
          
          if(lmmFit$compLikMats$compLikOptn == 0){
            
            # setupMatsMap <- setupIAK3D(xData = xMap[iThis[iOKInThis],,drop=FALSE] , dIData = dIMapThis[iOKInThis,,drop=FALSE] , nDscPts = 0)
            
            setupMatsMap <- setupIAK3D(xData = xMap[iThis[iOKInThis],,drop=FALSE] , dIData = dIMapThis[iOKInThis,,drop=FALSE] , nDscPts = 0 , 
                                       sdfdType_cd1 = lmmFit$sdfdType_cd1 , sdfdType_cxd0 = lmmFit$sdfdType_cxd0 , sdfdType_cxd1 = lmmFit$sdfdType_cxd1 , sdfdKnots = lmmFit$sdfdKnots)

            tmp <- setCIAK3D(parsBTfmd = lmmFit$parsBTfmd , modelx = lmmFit$modelx , 
                             sdfdType_cd1 = lmmFit$sdfdType_cd1 , sdfdType_cxd0 = lmmFit$sdfdType_cxd0 , sdfdType_cxd1 = lmmFit$sdfdType_cxd1 , 
                             cmeOpt = lmmFit$cmeOpt , setupMats = setupMatsMap)
            rm(setupMatsMap) 
            
            sigma2Veck <- tmp$sigma2Vec
            Ckk <- diag(tmp$C)
            rm(tmp)
            
            # setupMatsMap <- setupIAK3D2(xData = xMap[iThis[iOKInThis],,drop=FALSE] , dIData = dIMapThis[iOKInThis,,drop=FALSE] , 
            #                             xData2 = lmmFit$xData , dIData2 = as.matrix(lmmFit$dIData))
            setupMatsMap <- setupIAK3D2(xData = xMap[iThis[iOKInThis],,drop=FALSE] , dIData = dIMapThis[iOKInThis,,drop=FALSE] , 
                                        xData2 = lmmFit$xData , dIData2 = as.matrix(lmmFit$dIData) , 
                                        sdfdType_cd1 = lmmFit$sdfdType_cd1 , sdfdType_cxd0 = lmmFit$sdfdType_cxd0 , sdfdType_cxd1 = lmmFit$sdfdType_cxd1 , sdfdKnots = lmmFit$sdfdKnots)

            Ckh <- setCIAK3D2(parsBTfmd = lmmFit$parsBTfmd , modelx = lmmFit$modelx , 
                              sdfdType_cd1 = lmmFit$sdfdType_cd1 , sdfdType_cxd0 = lmmFit$sdfdType_cxd0 , sdfdType_cxd1 = lmmFit$sdfdType_cxd1 , 
                              cmeOpt = lmmFit$cmeOpt , setupMats = setupMatsMap)
            rm(setupMatsMap) 

            Ckhtmp <- Ckh %*% cbind(iCX , iCz_muhat , iC)
            CkhiCX <- Ckhtmp[,1:p , drop = FALSE]
            CkhiCz_muhat <- Ckhtmp[,p+1 , drop = FALSE]
            CkhiC <- Ckhtmp[,(p+2):(dim(Ckhtmp)[[2]]) , drop = FALSE]
            rm(Ckhtmp) 
            
            CkhiCChk <- rowSums(CkhiC * Ckh)
          }else{
            ### define stuff via CL method...          
            if(lmmFit$compLikMats$compLikOptn == 2 || lmmFit$compLikMats$compLikOptn == 3){
              tmp <- predMatsIAK3D_CLEV_ByBlock(z_muhat = lmmFit$zData - muhat , XData = lmmFit$XData , 
                                                xMap = xMap[iThis[iOKInThis],,drop=FALSE] , dIMap = dIMapThis[iOKInThis,,drop=FALSE], lmmFit = lmmFit)
              zkTmp <- tmp$zk
              vkTmp <- tmp$vk
              rm(tmp)
            }else{
              tmp <- predMatsIAK3D_CL_ByBlock(z_muhat = lmmFit$zData - muhat , XData = lmmFit$XData , 
                                              xMap = xMap[iThis[iOKInThis],,drop=FALSE] , dIMap = dIMapThis[iOKInThis,,drop=FALSE], lmmFit = lmmFit)
              
              Ckk <- tmp$diagCkk 
              CkhiCChk <- tmp$diagCkhiCChk 
              CkhiCX <- tmp$CkhiCX 
              CkhiCz_muhat <- tmp$CkhiCz_muhat
              rm(tmp)
            }
          }
          
          if(!lmmFit$lnTfmdData){
            muhatMapThis <- XMapThis %*% betahat
          }else{
            muhatMapThis <- setmuIAK3D(XData = XMapThis , vXU = vXUMapThis , iU = lmmFit$iU , beta = betahat , diagC = Ckk , sigma2Vec = sigma2Veck)
          }
          
          if(lmmFit$compLikMats$compLikOptn == 2 || lmmFit$compLikMats$compLikOptn == 3){
            zMapThis[iOKInThis] <- as.numeric(muhatMapThis) + as.numeric(zkTmp)     
            vMapThis[iOKInThis] <- as.numeric(vkTmp)
          }else{
            
            zMapThis[iOKInThis] <- as.numeric(muhatMapThis) + as.numeric(CkhiCz_muhat)     
            if(!lmmFit$lnTfmdData){
              tmp <- XMapThis - CkhiCX
              if(lmmFit$useReml){
                vMapThis[iOKInThis] <- Ckk - CkhiCChk + rowSums(tmp * t(vbetahat %*% t(tmp)))
              }else{
                vMapThis[iOKInThis] <- Ckk - CkhiCChk #TEMP WITHOUT beta unc...
              }

            }else{
              ### don't include uncerainty due to fixed effects as not correct formula in lognomral case. Could do as TS-FIM approx in future. 
              vMapThis[iOKInThis] <- Ckk - CkhiCChk
            }
          }
          
        }else{}
        
        zMap[i,iThis] <- as.numeric(zMapThis)
        vMap[i,iThis] <- as.numeric(vMapThis)
        
        if((i == 1) & (iBatch == 1)){
          tTmp <- proc.time() - ptm
            tTmp <- ceiling(tTmp[3] * nBatches * ndIMap / 60)
            print(paste0('...should take about ' , tTmp , ' minutes...'))
        }else{}

      } # end of iBatch loop

    } # end of i loop
    
    pi90LMap <- zMap - 1.645 * sqrt(vMap)
    pi90UMap <- zMap + 1.645 * sqrt(vMap)
    
### to back transform...
    if(lmmFit$lnTfmdData & rqrBTfmdPreds){
      zMap <- exp(zMap + 0.5 * vMap)
        
      pi90LMap <- exp(pi90LMap)
      pi90UMap <- exp(pi90UMap)

### because zPred and zPredDistant are now back-transformed,...
      vMap <- (exp(vMap) - 1) * (zMap ^ 2)
    }else{}

    return(list('zMap' = zMap , 'vMap' = vMap , 'pi90LMap' = pi90LMap , 'pi90UMap' = pi90UMap))
}


profilePredictIAK3D <- function(xMap , dIMap , covsMap , iData = seq(length(lmmFit$zData)) , lmmFit , rqrBTfmdPreds = TRUE , constrainX4Pred = FALSE){
########################################################
### if xMap or dIMap were dataframes, convert to matrices here.
### and make sure all are numeric...
########################################################
    if(!is.matrix(xMap)){
        xMap <- as.matrix(xMap)
    }else{}
    if(!is.matrix(dIMap)){
        dIMap <- as.matrix(dIMap)
    }else{}

    xMapCopy <- matrix(NA , nrow(xMap) , ncol(xMap))
    for (i in 1:ncol(xMap)){ xMapCopy[,i] <- as.numeric(xMap[,i]) }
    xMap <- xMapCopy
    remove(xMapCopy)
    
    dIMapCopy <- matrix(NA , nrow(dIMap) , ncol(dIMap))
    for (i in 1:ncol(dIMap)){ dIMapCopy[,i] <- as.numeric(dIMap[,i]) }
    dIMap <- dIMapCopy
    remove(dIMapCopy)

##############################################
### a version to predict through the profile for one location
### iData allows a subset of the full dataset to be used for prediction
##############################################
    dIMap <- round(dIMap , digits = 2)
    ndIMap <- dim(dIMap)[[1]]

    nData <- length(iData)

    betahat <- lmmFit$betahat
    vbetahat <- lmmFit$vbetahat

    p <- dim(lmmFit$XData)[[2]]
    pU <- length(lmmFit$iU)

    if(constrainX4Pred){
      XLims4Pred <- lmmFit$XLims
    }else{
      XLims4Pred <- matrix(NA , 2 , p)
      XLims4Pred[1,] <- -Inf
      XLims4Pred[2,] <- Inf
    }
    
    if(!lmmFit$lnTfmdData){
      muhat <- lmmFit$XData %*% betahat 
    }else{
      muhat <- setmuIAK3D(XData = lmmFit$XData , vXU = lmmFit$vXU , iU = lmmFit$iU , beta = betahat , diagC = diag(lmmFit$C) , sigma2Vec = lmmFit$sigma2Vec) 
    }
 
    if(lmmFit$compLikMats$compLikOptn == 0){
      if((length(iData) == length(lmmFit$zData)) & is.element('iCX'  , names(lmmFit)) & is.element('iCz_muhat'  , names(lmmFit)) & is.element('iC'  , names(lmmFit))){
        iCX <- lmmFit$iCX
        iCz_muhat <- lmmFit$iCz_muhat
        iC <- lmmFit$iC
      }else{
#        tmp <- lndetANDinvCb(lmmFit$C[iData,iData,drop=FALSE] , cbind(lmmFit$XData[iData,,drop=FALSE] , lmmFit$zData[iData] - muhat[iData]))
#        iCX <- tmp$invCb[,1:p,drop=FALSE]
#        iCz_muhat <- tmp$invCb[,p+1,drop=FALSE]
#        iC <- chol2inv(tmp$cholC)

        iC <- chol2inv(chol(lmmFit$C[iData,iData,drop=FALSE]))
        iCXResTmp <- matrix(iC %*% cbind(lmmFit$XData[iData,,drop=FALSE] , lmmFit$zData[iData] - muhat[iData]) , ncol = p+1)
        iCX <- iCXResTmp[,1:p,drop=FALSE]
        iCz_muhat <- iCXResTmp[,p+1,drop=FALSE]
        remove(iCXResTmp)
      }
    }else{}
    
### repeat the prediction location and covariates,,,
    xMap <- matrix(xMap , nrow = 1)
    xMap <- xMap[integer(ndIMap) + 1,,drop=FALSE]
    covsMap <- covsMap[integer(ndIMap) + 1,,drop = FALSE]

    if(identical(lmmFit$modelX$type , 'gam2')){
      tmp <- makeXvX_gam2(covData = covsMap , dIData = dIMap , listfefdKnots = lmmFit$modelX$listfefdKnots , incInts = lmmFit$modelX$incInts , intMthd = lmmFit$modelX$intMthd , colnamesXcns = lmmFit$modelX$colnamesX , nDiscPts = 10 , lnTfmdData = lmmFit$lnTfmdData)
    }else{
      tmp <- makeXvX(covData = covsMap , dIData = dIMap , modelX = lmmFit$modelX , allKnotsd = lmmFit$allKnotsd , iU = lmmFit$iU , nDiscPts = 10 , lnTfmdData = lmmFit$lnTfmdData , XLims = XLims4Pred)
    }
    XMap <- as.matrix(tmp$X)
    if(lmmFit$lnTfmdData){
      vXUMap <- as.matrix(tmp$vXU)
    }else{}

    if(lmmFit$compLikMats$compLikOptn == 0){
      # setupMatsMap <- setupIAK3D(xData = xMap , dIData = dIMap , nDscPts = 0)
      setupMatsMap <- setupIAK3D(xData = xMap , dIData = dIMap , nDscPts = 0 , 
                                 sdfdType_cd1 = lmmFit$sdfdType_cd1 , sdfdType_cxd0 = lmmFit$sdfdType_cxd0 , sdfdType_cxd1 = lmmFit$sdfdType_cxd1 , sdfdKnots = lmmFit$sdfdKnots)
      
      tmp <- setCIAK3D(parsBTfmd = lmmFit$parsBTfmd , modelx = lmmFit$modelx , 
                       sdfdType_cd1 = lmmFit$sdfdType_cd1 , sdfdType_cxd0 = lmmFit$sdfdType_cxd0 , sdfdType_cxd1 = lmmFit$sdfdType_cxd1 , 
                       cmeOpt = lmmFit$cmeOpt , setupMats = setupMatsMap)
      rm(setupMatsMap) 
      
      sigma2Veck <- tmp$sigma2Vec
      Ckk <- diag(tmp$C)
      
      # setupMatsMap <- setupIAK3D2(xData = xMap , dIData = dIMap , 
      #                             xData2 = lmmFit$xData[iData,,drop=FALSE] , dIData2 = as.matrix(lmmFit$dIData[iData,,drop=FALSE]))
      setupMatsMap <- setupIAK3D2(xData = xMap , dIData = dIMap , 
                                  xData2 = lmmFit$xData[iData,,drop=FALSE] , dIData2 = as.matrix(lmmFit$dIData[iData,,drop=FALSE]) , 
                                  sdfdType_cd1 = lmmFit$sdfdType_cd1 , sdfdType_cxd0 = lmmFit$sdfdType_cxd0 , sdfdType_cxd1 = lmmFit$sdfdType_cxd1 , sdfdKnots = lmmFit$sdfdKnots)

      Ckh <- setCIAK3D2(parsBTfmd = lmmFit$parsBTfmd , modelx = lmmFit$modelx , 
                        sdfdType_cd1 = lmmFit$sdfdType_cd1 , sdfdType_cxd0 = lmmFit$sdfdType_cxd0 , sdfdType_cxd1 = lmmFit$sdfdType_cxd1 , 
                        cmeOpt = lmmFit$cmeOpt , setupMats = setupMatsMap)
      
      rm(tmp,setupMatsMap) 

      Ckhtmp <- Ckh %*% cbind(iCX , iCz_muhat , iC)
      CkhiCX <- Ckhtmp[,1:p]
      CkhiCz_muhat <- Ckhtmp[,p+1]
      CkhiC <- Ckhtmp[,(p+2):(dim(Ckhtmp)[[2]])]
      rm(Ckhtmp) ; gc()
      CkhiCChk <- rowSums(CkhiC * Ckh)

    }else{
### define stuff via CL method...          
      if(lmmFit$compLikMats$compLikOptn == 2 || lmmFit$compLikMats$compLikOptn == 3){
      
        tmp <- predMatsIAK3D_CLEV_ByBlock(z_muhat = lmmFit$zData - muhat , XData = lmmFit$XData , 
              xMap = xMap , dIMap = dIMap , iData = iData , lmmFit = lmmFit)
        zkTmp <- tmp$zk
        vkTmp <- tmp$vk
      }else{
### define stuff via CL method...          
        tmp <- predMatsIAK3D_CL_ByBlock(z_muhat = lmmFit$zData - muhat , XData = lmmFit$XData , 
                xMap = xMap , dIMap = dIMap , iData = iData , lmmFit = lmmFit)
        Ckk <- tmp$diagCkk 
        CkhiCChk <- tmp$diagCkhiCChk 
        CkhiCX <- tmp$CkhiCX 
        CkhiCz_muhat <- tmp$CkhiCz_muhat
      }
    }

    if(!lmmFit$lnTfmdData){
        muhatMap <- XMap %*% betahat
    }else{
        muhatMap <- setmuIAK3D(XData = XMap , vXU = vXUMap , iU = lmmFit$iU , beta = betahat , diagC = Ckk , sigma2Vec = sigma2Veck)
    }
    
    if(lmmFit$compLikMats$compLikOptn == 2 || lmmFit$compLikMats$compLikOptn == 3){
      zMap <- as.numeric(muhatMap) + as.numeric(zkTmp)     
      vMap <- as.numeric(vkTmp)
    }else{
      zMap <- muhatMap + CkhiCz_muhat     
      if(!lmmFit$lnTfmdData){
        tmp <- XMap - CkhiCX
        if(lmmFit$useReml){
          vMap <- Ckk - CkhiCChk + rowSums(tmp * t(vbetahat %*% t(tmp)))
        }else{
          vMap <- Ckk - CkhiCChk # TEMP WO BETA UNC.
        }
      }else{
### don't include uncerainty due to fixed effects as not correct formula in lognomral case. Could do as TS-FIM approx in future. 
### though maybe with vbetahat based on FIM it works ok, see one of Gerards paper's /Ben's uncertainty paper
        vMap <- Ckk - CkhiCChk
      }
    }       
    pi90LMap <- zMap - 1.645 * sqrt(vMap)
    pi90UMap <- zMap + 1.645 * sqrt(vMap)
    
### to back transform...
    if(lmmFit$lnTfmdData & rqrBTfmdPreds){
      zMap <- exp(zMap + 0.5 * vMap)
        
      pi90LMap <- exp(pi90LMap)
      pi90UMap <- exp(pi90UMap)

### because zPred and zPredDistant are now back-transformed,...
      vMap <- (exp(vMap) - 1) * (zMap ^ 2)
    }else{}

    return(list('zMap' = zMap , 'vMap' = vMap , 'pi90LMap' = pi90LMap , 'pi90UMap' = pi90UMap , 'XMap' = XMap , 'muhatMap' = muhatMap))
}

xValIAK3D <- function(lmmFit , removeAllWithin = 0 , namePlot = 'xvPlots.pdf' , rqrBTfmdPreds = TRUE){
####################################################
### make full profile predictions for plotting 
### and predictions at data supports for validation
### remove full profiles at a time, and also any other profiles within 
### a distance of removeAllWithin of the validation location
###
### rqrBTfmdPreds is just for the plots...
####################################################
    iTmp <- which(!duplicated(lmmFit$xData))
    covsPred <- lmmFit$covsData[iTmp,,drop = FALSE]
    xPred <- lmmFit$xData[iTmp,,drop = FALSE]

    nxPred <- dim(xPred)[[1]]
    nData <- length(lmmFit$zData)
    p <- dim(lmmFit$XData)[[2]]

    dIPredPlot <- cbind(seq(0 , 1.98 , 0.02) , seq(0.02 , 2 , 0.02))
    ndIPredPlot <- dim(dIPredPlot)[[1]]

    print(paste0('Cross-validating for ' , nxPred , 'locations...'))
    ptm <- proc.time()

    zhatPlot <- vhatPlot <- pi90LPlot <- pi90UPlot <- matrix(NA , ndIPredPlot , nxPred)
    zhatxv <- vhatxv <- pi90Lxv <- pi90Uxv <- NA * numeric(length(nData))
    for (i in 1:nxPred){
        xPredThis <- matrix(xPred[i,] , nrow = 1)
        covsPredThis <- covsPred[i,,drop = FALSE]

### get all the dIData to be predicted for this location...
        iPredThis <- which((lmmFit$xData[,1] == xPredThis[1]) & (lmmFit$xData[,2] == xPredThis[2]))
        dIPredThis <- as.matrix(lmmFit$dI[iPredThis,])

### get the prediction data, removing any locations closer than removeAllWithin 
        DTmp <- xyDist(lmmFit$xData , xPredThis)
        DTmp[iPredThis] <- -999 # just to make sure they will be removed.
        iDataThis <- which(DTmp >= removeAllWithin)

### use the profilePredict function (call without back-transform, to get stdzd sqd errs on transformed scale)...
### note that for xval, columns of XData are not constrained. 
        tmp <- profilePredictIAK3D(xMap = xPredThis , dIMap = rbind(matrix(dIPredPlot , ncol = 2) , matrix(dIPredThis , ncol = 2)) ,
                    covsMap = covsPredThis , iData = iDataThis , lmmFit = lmmFit , rqrBTfmdPreds = FALSE)

### extract the results...
### i won't back-transform xv values for calculating sspe stats.
        zhatxv[iPredThis] <- tmp$zMap[(ndIPredPlot+1):(ndIPredPlot+length(iPredThis))]
        vhatxv[iPredThis] <- tmp$vMap[(ndIPredPlot+1):(ndIPredPlot+length(iPredThis))]
        pi90Lxv[iPredThis] <- tmp$pi90LMap[(ndIPredPlot+1):(ndIPredPlot+length(iPredThis))]
        pi90Uxv[iPredThis] <- tmp$pi90UMap[(ndIPredPlot+1):(ndIPredPlot+length(iPredThis))]

        zhatPlot[,i] <- tmp$zMap[1:ndIPredPlot]
        vhatPlot[,i] <- tmp$vMap[1:ndIPredPlot]
        pi90LPlot[,i] <- tmp$pi90LMap[1:ndIPredPlot]
        pi90UPlot[,i] <- tmp$pi90UMap[1:ndIPredPlot]
        if(i == 1){
            tTmp <- proc.time() - ptm
            tTmp <- ceiling(tTmp[3] * nxPred / 60)
            print(paste0('...should take about ' , tTmp , ' minutes...'))
        }else{}
    }

##################################################
### copmute stdzd sqd errs on transformed scale for validating prediction variances...
##################################################
    stdzdSqdErrs <- ((zhatxv - lmmFit$zData) ^2) / vhatxv

    zData <- lmmFit$zData
##################################################
### back-transform predictions and 90% pis here, if rqd... 
##################################################
    if(lmmFit$lnTfmdData & rqrBTfmdPreds){
      zhatxv <- exp(zhatxv + 0.5 * vhatxv)        
      pi90Lxv<- exp(pi90Lxv)
      pi90Uxv<- exp(pi90Uxv)
      vhatxv <- (exp(vhatxv) - 1) * (zhatxv ^ 2)

      zhatPlot <- exp(zhatPlot + 0.5 * vhatPlot)        
      pi90LPlot <- exp(pi90LPlot)
      pi90UPlot <- exp(pi90UPlot)
      vhatPlot <- (exp(vhatPlot) - 1) * (zhatPlot ^ 2)

      zData <- exp(zData)
    }else{}

#####################################################
### calculate the validation statistics to assess prediction quality...
#####################################################
    errs <- zhatxv - zData
    sqdErrs <- errs ^ 2

    errs0 <- mean(zData) - zData
    sqdErrs0 <- errs0 ^ 2

############################
### make a matrix of xv results: 
###    overall results, results split by depth 
###    (based on midpoints, for three depths; < 0.2 , -0.5 , >= 0.5)
### Columns are:
###     bias, rmse, r2 , mean stdzd sqd errs, median stdzd sqd errs, n, 
###     approx bounds for mean stdzd sqd errs (l90 and u90) if fair
###     approx bounds for median stdzd sqd errs (l90 and u90) if fair
############################
    dLims <- cbind(c(-Inf , 0.2 , 0.5) , c(0.2 , 0.5 , Inf))
    xvStats <- matrix(NA , 4 , 9)
    rownames(xvStats) <- c('Overall' , 'dMidpnt < 0.2' , '0.2 <= dMidpnt < 0.5' , 'dMidpnt >= 0.5')
    colnames(xvStats) <- c('Bias' , 'RMSE' , 'Mean SSE' , 'Median SSE' , 'n' , 
            'Approx theretical L90 MSSE' , 'Approx theretical U90 MSSE' ,
            'Approx theretical L90 MedSSE' , 'Approx theretical U90 MedSSE')
### overall stats
    xvStats[1,1] <- mean(errs)
    xvStats[1,2] <- sqrt(mean(sqdErrs))
    xvStats[1,3] <- mean(stdzdSqdErrs)
    xvStats[1,4] <- median(stdzdSqdErrs)
    xvStats[1,5] <- length(errs)

    xvStats[1,6] <- 1 - 1.645 * sqrt(2/length(errs))
    xvStats[1,7] <- 1 + 1.645 * sqrt(2/length(errs))

    sdmed <- sqrt(1/(8*((length(errs)-1)/2)*0.2219)) # 0.2219 is chisqpdf(0.455, 1df) ^ 2
    xvStats[1,8] <- 0.455 - 1.645 * sdmed
    xvStats[1,9] <- 0.455 + 1.645 * sdmed

### by depth, using midpoints to classify...
    for (id in 1:3){    
        iThis <- which((rowMeans(lmmFit$dI) >= dLims[id,1]) & (rowMeans(lmmFit$dI) < dLims[id,2]))

        xvStats[id+1,1] <- mean(errs[iThis])
        xvStats[id+1,2] <- sqrt(mean(sqdErrs[iThis]))
        xvStats[id+1,3] <- mean(stdzdSqdErrs[iThis])
        xvStats[id+1,4] <- median(stdzdSqdErrs[iThis])
        xvStats[id+1,5] <- length(iThis)

        xvStats[id+1,6] <- 1 - 1.645 * sqrt(2/length(iThis))
        xvStats[id+1,7] <- 1 + 1.645 * sqrt(2/length(iThis))

        sdmed <- sqrt(1/(8*((length(iThis)-1)/2)*0.2219)) # 0.2219 is chisqpdf(0.455, 1df) ^ 2
        xvStats[id+1,8] <- 0.455 - 1.645 * sdmed
        xvStats[id+1,9] <- 0.455 + 1.645 * sdmed
    }

########################################
### plot...
########################################
    tmp <- plotProfilesIAK3D(namePlot = namePlot , xData = lmmFit$xData , dIData = lmmFit$dIData, zData = lmmFit$zData , 
                xPred = xPred , dIPred = dIPredPlot , zPred = zhatPlot , pi90LPred = pi90LPlot , pi90UPred = pi90UPlot , 
                zhatxv = zhatxv , pi90Lxv = pi90Lxv , pi90Uxv = pi90Uxv)

########################################
### and return...
########################################
    return(list('xvStats' = xvStats , 'zhatxv' = zhatxv , 'vhatxv' = vhatxv , 'pi90Lxv' = pi90Lxv , 'pi90Uxv' = pi90Uxv ,
                'zhatPlot' = zhatPlot , 'vhatPlot' = vhatPlot , 'pi90LPlot' = pi90LPlot  , 'pi90UPlot' = pi90UPlot)) 
}

plotProfilesIAK3D <- function(namePlot = 'profilePlots.pdf' , xData , dIData , zData , xPred = NULL , dIPred = NULL , zPred = NULL  , pi90LPred = NULL , pi90UPred = NULL , 
                              dIStd = NULL , zStd = NULL , pi90LStd = NULL , pi90UStd = NULL , 
                              zPredDistant = NULL , zhatxv = NULL , pi90Lxv = NULL , pi90Uxv = NULL , profNames = NULL , xlim = NULL , xlab = NULL){
#################################################    
### make a pdf with the distant profile prediction (page 1) and all data profiles (6 per page thereafter)...
### note these are not validation predictions, they are predicted at the data profiles given the data for the same profiles
### also note, they are a bi-product of the method (ie you can use the method to predict at profiles where we have data), 
### not to be confused with the spline-then-krige type approach where similar plots may be produced in the first step of analysis
###
### xData, dIData and zData are the raw data (all horizons)
### xPred, etc are unique locations and assciated profile predictions.
#################################################    
########################################################
### if xData or dIData were dataframes, convert to matrices here.
### and make sure dI numeric, and x not a factor...
########################################################
    if(!is.matrix(xData)){
        xData <- as.matrix(xData)
    }else{}
    if(!is.matrix(dIData)){
        dIData <- as.matrix(dIData)
    }else{}

    for (i in 1:ncol(xData)){ if(is.factor(xData[,i])){ stop('Enter xData as either coordinates or as character vec of ids') }else{} }

    dICopy <- matrix(NA , nrow(dIData) , ncol(dIData))
    for (i in 1:ncol(dIData)){ dICopy[,i] <- as.numeric(dIData[,i]) }
    dIData <- dICopy
    remove(dICopy)

########################################################
### if xData or dIData were dataframes, convert to matrices here.
### and make sure all are numeric...
########################################################
    if(!is.null(xPred)){
      if(!is.matrix(xPred)){
        xPred <- as.matrix(xPred)
      }else{}
      if(!is.matrix(dIPred)){
        dIPred <- as.matrix(dIPred)
      }else{}

      for (i in 1:ncol(xPred)){ if(is.factor(xPred[,i])){ stop('Enter xPred as either coordinates or as character vec of ids') }else{} }

      dIPredCopy <- matrix(NA , nrow(dIPred) , ncol(dIPred))
      for (i in 1:ncol(dIPred)){ dIPredCopy[,i] <- as.numeric(dIPred[,i]) }
      dIPred <- dIPredCopy
      remove(dIPredCopy)
    }else{}

    ndx <- ncol(xData)
    
    nPerPage <- 6
    
    if (is.null(xPred)){
      xPred <- xData[!duplicated(xData) , ,drop = FALSE]
      plotPreds <- FALSE
      zPredDistant <- NULL
    }else{
      plotPreds <- TRUE
    }
    nPages <- ceiling(dim(xPred)[[1]] / nPerPage)

    if(plotPreds){
      maxd <- max(dIPred[,2])
    }else{
      maxd <- max(dIData[,2])
    }
    ylim  <- c(-maxd , 0)
    if(maxd <= 0.5){
      yaxTck <- seq(-0.5 , 0 , 0.1)
      yaxTckLbls <- c('0.5' , '0.4' , '0.3' , '0.2' , '0.1' , '0.0')
    }else if(maxd <= 1.0){
      yaxTck <- seq(-1 , 0 , 0.2)
      yaxTckLbls <- c('1.0' , '0.8' , '0.6' , '0.4' , '0.2' , '0.0')
    }else if(maxd <= 2.0){
      yaxTck <- seq(-2 , 0 , 0.5)
      yaxTckLbls <- c('2.0' , '1.5' , '1.0' , '0.5' , '0.0')
    }else if(maxd <= 5.0){
      yaxTck <- seq(-5 , 0 , 1)
      yaxTckLbls <- c('5.0' , '4.0' , '3.0' , '2.0' , '1.0' , '0.0')
    }else if(maxd <= 10.0){
      yaxTck <- seq(-10 , 0 , 2)
      yaxTckLbls <- c('10.0' , '8.0' , '6.0' , '4.0' , '2.0' , '0.0')
    }else if(maxd <= 20.0){
      yaxTck <- seq(-20 , 0 , 5)
      yaxTckLbls <- c('20' , '15' , '10' , '5' , '0')
    }else{}
    
    if(is.null(xlim)){
      rangez <- max(zData) - min(zData)
      xlim <- c(min(zData) - 0.01 * rangez , max(zData) + 0.01 * rangez)
      if(xlim[2] - xlim[1] < 0.1){ xlim <- c(xlim[1] - 0.1 , xlim[2] + 0.1) }else{} 
    }else{}
    if(is.null(xlab)){
      xlab <- 'z'
    }else{}
    
    pdf(file = namePlot)
    if(!is.null(zPredDistant)){
      if(is.character(xlim) && xlim == 'flex'){
        xlimThis <- c(min(c(zData,zPredDistant)) , max(c(zData,zPredDistant)))
        if(xlimThis[2] - xlimThis[1] < 0.1){ xlimThis <- c(xlimThis[1] - 0.1 , xlimThis[2] + 0.1) }else{} 
      }else{
        xlimThis <- xlim
      }
      plot(zPredDistant , -rowMeans(dIPred) , xlim = xlimThis , ylim = ylim , type = 'l' , lwd = 2 , col = 'red' ,  
            xlab = xlab , ylab = 'depth, m' , yaxt = 'n' , main = 'Distant profile')
### add all data to this plot...                  
      for (j in 1:length(zData)){
        lines(zData[j] * c(1,1) , -dIData[j,] , lty = 1 , lwd = 1 , col = 'black')
      }
### re-add the prediction line to make it on top of the data in this case...
      lines(zPredDistant , -rowMeans(dIPred) , lwd = 2 , col = 'red')

      axis(side = 2 , at = yaxTck , labels = yaxTckLbls , las = 2)
    }else{}

    for (ip in 1:nPages){
      par(mfrow = c(2,3))
      for (i in 1:nPerPage){
        iProfThis <- (ip - 1) * nPerPage + i
        if(iProfThis <= dim(xPred)[[1]]){
          if(is.null(profNames)){
            nameThis <- paste0('Profile ' , iProfThis)
          }else{
            nameThis <- profNames[iProfThis]
          }      
### get data for this profile...
          if (ndx == 1){
            iThis <- which(xData[,1] == xPred[iProfThis,1])
          }else if(ndx == 2){
            iThis <- which((xData[,1] == xPred[iProfThis,1]) & (xData[,2] == xPred[iProfThis,2]))       
          }else if(ndx == 3){
            iThis <- which((xData[,1] == xPred[iProfThis,1]) & (xData[,2] == xPred[iProfThis,2]) & (xData[,3] == xPred[iProfThis,3]))       
          }else{
            stop('x is 4d or more? Check this.')            
          }
          
          if(is.character(xlim) && xlim == 'flex'){
            ### justuse data and preds to define xlim...              
            if(plotPreds){
              if(!is.null(dIStd)){
                xlimThis <- c(min(c(zPred[,iProfThis],zData[iThis],zStd[iProfThis,]),na.rm=T) , max(c(zPred[,iProfThis],zData[iThis],zStd[iProfThis,]),na.rm=T))
              }else{
                xlimThis <- c(min(c(zPred[,iProfThis],zData[iThis]),na.rm=T) , max(c(zPred[,iProfThis],zData[iThis]),na.rm=T))
              }
            }else{
              xlimThis <- c(min(zData[iThis]) , max(zData[iThis]))
            }
          }else{
            xlimThis <- xlim
          }
          if(xlimThis[2] - xlimThis[1] < 0.1){ xlimThis[1] <- xlimThis[1] - 0.1 ; xlimThis[2] <- xlimThis[2] + 0.1 }else{}

          if(plotPreds){
            plot(zPred[,iProfThis] , -rowMeans(dIPred) , xlim = xlimThis , ylim = ylim , type = 'l' , lwd = 2 , col = 'red' , 
              xlab = xlab , ylab = 'depth, m' , yaxt = 'n' , main = nameThis)
            if((!is.null(pi90LPred)) & (!is.null(pi90UPred))){
              lines(pi90LPred[,iProfThis] , -rowMeans(dIPred) , lty = 2 , lwd = 1 , col = 'red')
              lines(pi90UPred[,iProfThis] , -rowMeans(dIPred) , lty = 2 , lwd = 1 , col = 'red')
            }else{}
            
            if(!is.null(dIStd)){
              for(j in 1:nrow(dIStd)){
                lines(zStd[iProfThis,j] * c(1,1) , -dIStd[j,] , lty = 1 , lwd = 2 , col = 'green')
                if((!is.null(pi90LStd)) & (!is.null(pi90UStd))){
                  lines(pi90LStd[iProfThis,j] * c(1,1) , -dIStd[j,] , lty = 2 , lwd = 1 , col = 'green')
                  lines(pi90UStd[iProfThis,j] * c(1,1) , -dIStd[j,] , lty = 2 , lwd = 1 , col = 'green')
                }else{}
              }
            }else{}
          }else{
            plot(c() , c() , xlim = xlimThis , ylim = ylim , xlab = 'z' , ylab = 'depth, m' , yaxt = 'n' , main = nameThis)
          } 
          
          for (j in 1:length(iThis)){
            lines(zData[iThis[j]] * c(1,1) , -dIData[iThis[j],] , lty = 1 , lwd = 3 , col = 'black')
            if((!is.null(zhatxv)) && (length(zhatxv) >= iThis[j]) && (!is.na(zhatxv[iThis[j]]))){
              lines(zhatxv[iThis[j]] * c(1,1) , -dIData[iThis[j],] , lty = 1 , lwd = 2 , col = 'cyan')
              
              if((!is.null(pi90Lxv)) & (!is.null(pi90Uxv))){
                lines(pi90Lxv[iThis[j]] * c(1,1) , -dIData[iThis[j],] , lty = 2 , lwd = 1 , col = 'cyan')
                lines(pi90Uxv[iThis[j]] * c(1,1) , -dIData[iThis[j],] , lty = 2 , lwd = 1 , col = 'cyan')
              }else{}
            }else{}
          }
          axis(side = 2 , at = yaxTck , labels = yaxTckLbls , las = 2)
        }else{}
      }
    }
    dev.off()

    return()
}

#############################################################
### calculate lins CCC...
#############################################################
linsCCC <- function(o , p , na.rm = FALSE){
  
  if(na.rm){
    iNA <- which(is.na(o) | is.na(p))
    
    if(length(iNA) > 0){
      o <- o[-iNA]
      p <- p[-iNA]
    }else{}
  }else{}
  
  mp <- mean(p)
  mo <- mean(o)
  sp <- mean((p - mp) ^ 2)
  so <- mean((o - mo) ^ 2)
  spo <- mean((p-mp) * (o - mo))
  
  lCCC <- 2 * spo / (sp + so + (mp - mo) ^ 2)
  
  return(lCCC)
}

##########################################
### calc val stats...
##########################################
calcValStats <- function(zVal , dIVal , zkVal , vkVal , layerMidPts = c(0.025 , 0.1 , 0.225 , 0.45 , 0.8 , 1.5) , printValStats = TRUE){
  
  
  dTmp <- xyDist(layerMidPts , rowMeans(dIVal))
  ilayer <- apply(dTmp , 2 , which.min)

  inPI <- numeric(nrow(dIVal))
  inPI[which((zVal >= (zkVal - 1.64 * sqrt(vkVal))) & (zVal <= (zkVal + 1.64 * sqrt(vkVal))))] <- 1
  
  valStatsAllLayers <- data.frame('layerMidPts' = layerMidPts , 'bias' = NA , 'rmse' = NA ,  'R2' = NA , 'ccc' = NA , 'pInPI90' = NA)
  valStatsTot <- data.frame('bias' = NA , 'rmse' = NA ,  'ccc' = NA , 'pInP90' = NA)
  for(i in 1:length(layerMidPts)){
    iThis <- which(ilayer == i & rowMeans(dIVal) <= 2) # constrain so that v deep val data not incd
  
    zkThis <- zkVal[iThis]
    zThis <- zVal[iThis]
  
    valStatsAllLayers$bias[i] <- mean(zkThis - zThis)
    valStatsAllLayers$rmse[i] <- sqrt(mean((zkThis - zThis) ^ 2))
    valStatsAllLayers$R2[i] <- 1 - sum((zkThis - zThis) ^ 2) / sum((mean(zThis) - zThis) ^ 2)
    valStatsAllLayers$ccc[i] <- linsCCC(zkThis , zThis)
    valStatsAllLayers$pInPI90[i] <- mean(inPI[iThis]) 
  }

  valStatsTot$bias <- mean(zkVal - zVal)
  valStatsTot$rmse <- sqrt(mean((zkVal - zVal) ^ 2))
  valStatsTot$R2 <- 1 - sum((zkVal - zVal) ^ 2) / sum((mean(zVal) - zVal) ^ 2)
  valStatsTot$ccc <- linsCCC(zkVal , zVal)
  valStatsTot$pInPI90 <- mean(inPI) 

##############################
### print out val stats    
##############################
  if(printValStats){
    print('Overall prob in PI90:')
    print(valStatsTot$pInPI90)
    print('Prob in PI90 for each of the 6 layers:')
    print(valStatsAllLayers$pInPI90)
  
    print('Overall RMSE:')
    print(valStatsAllLayers$rmse)
    print('RMSE for each of the 6 layers:')
    print(valStatsTot$rmse)
  
    print('Overall CCC:')
    print(valStatsAllLayers$ccc)
    print('CCC for each of the 6 layers:')
    print(valStatsTot$ccc)
  }

  return(list('valStatsAllLayers' = valStatsAllLayers , 'valStatsTot' = valStatsTot))
}
