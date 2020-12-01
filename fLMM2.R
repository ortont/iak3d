fitLMM2 <- function(c , z , X , covModel = 'exponential' , nSpatStructs = 1 , incNugget = T , optionML = F , verbose = FALSE , mina = NULL , maxa = NULL , parsInit = NULL , attachBigMats = TRUE){

  ina <- which(is.na(z))
  if(length(ina) > 0){
    c <- c[-ina,,drop=FALSE]
    z <- z[-ina]
    X <- X[-ina,,drop=FALSE]
  }else{}
  
  D <- xyDist(c , c)

### note parameterisation in fLMM2 fn for exp model is with a = range (not 3a = range)
  if(is.null(mina)){ mina <- min(D[D > 0]) }else{}
  if(is.null(maxa)){ maxa <- max(D) / 2 }else{}

  if(is.null(parsInit)){
    if(incNugget){
      if(nSpatStructs == 1){
        parsInit <- c(0.6 , 0.5 * mina + 0.5 * maxa)
      }else if(nSpatStructs == 2){
        parsInit <- c(0.3 , 0.9 * mina + 0.1 * maxa , 0.3 , 0.3 * mina + 0.7 * maxa)
      }else{
        stop('Error - generalise fitLMM2 for other values of nSpatStructs!')
      }
    }else{
      if(nSpatStructs == 1){
        parsInit <- c(0.5 * mina + 0.5 * maxa)
      }else if(nSpatStructs == 2){
        parsInit <- c(0.5 , 0.9 * mina + 0.1 * maxa , 0.3 * mina + 0.7 * maxa)
      }else{
        stop('Error - generalise fitLMM2 for other values of nSpatStructs!')
      }
    }
  }else{}
  
#####################################################    
### run to get initial lik val...
#####################################################    
  ftLMM2Init <- fLMM2(pars = parsInit , c = c , z = z , X = X , D = D , 
                    covModel = covModel , nSpatStructs = nSpatStructs , incNugget = incNugget , returnAll = FALSE , optionML = optionML , 
                    verbose = TRUE , forCompLik = FALSE , mina = mina , maxa = maxa)

#####################################################    
### fit model with temporal correlation accounted for...
#####################################################    
  badNll <- 9E99
  if(length(parsInit) > 1){
    tmp <- optimIt(par = parsInit , fn = fLMM2 , methodOptim = c('Nelder-Mead') , c = c , z = z , X = X , D = D ,
               covModel = covModel , nSpatStructs = nSpatStructs , incNugget = incNugget , returnAll = FALSE , optionML = optionML , 
               verbose = verbose , forCompLik = FALSE , mina = mina , maxa = maxa , badNll = badNll)
    parsFit <- tmp$par
  }else{    
    tmp <- optimize(f = fLMM2 , lower = mina , upper = maxa , c = c , z = z , X = X , D = D ,
               covModel = covModel , nSpatStructs = nSpatStructs , incNugget = incNugget , returnAll = FALSE , optionML = optionML , 
               verbose = verbose , forCompLik = FALSE , mina = mina , maxa = maxa , badNll = badNll)
    parsFit <- tmp$minimum
  }
  
#####################################################    
### run again to get final lik val + other info...
#####################################################    
  ftLMM2Fit <- fLMM2(pars = parsFit , c = c , z = z , X = X , D = D , 
                      covModel = covModel , nSpatStructs = nSpatStructs , incNugget = incNugget , returnAll = TRUE , optionML = optionML , 
                      verbose = TRUE , forCompLik = FALSE , mina = mina , maxa = maxa , attachBigMats = attachBigMats)
  
  return(ftLMM2Fit)
}

fLMM2 <- function(pars , c , z , X , D , covModel , nSpatStructs , incNugget = T , returnAll , optionML , verbose , forCompLik = FALSE , mina = NULL , maxa = NULL , badNll = 9E99 , attachBigMats = TRUE){
    N <- length(z)
    np <- dim(X)[2]

    if(missing(returnAll)){ returnAll <- F }else{} # default is to return just nll
    if(missing(optionML)){ optionML <- F }else{} # default is REML
    if(missing(verbose)){ verbose <- T }else{} # default is to print the current nll

### a1 must be at least the minimum sep dist...
    if(is.null(mina)){ mina <- min(D + 9E9 * diag(N)) }else{}
### a2 (or a1 if 1 spat struct) must be less than the following value...
    if(is.null(maxa)){ maxa <- max(D) / 2 }else{}
    
### assume nugget. no matern yet.
### pars will be..
### if nSpatStructs = 1...s,a
### if nSpatStructs = 2, will be ...s1,a1,s2,a2
    if(incNugget){
      if (nSpatStructs == 1){
        s1 <- pars[1] ; a1 <- pars[2] ; s2 <- 0 ; a2 <- a1 ; s0 <- 1 - s1
        if((s1 < 0) | (s1 > 1) | (a1 < mina) | (a1 > maxa)){ paramsOK <- F }else{ paramsOK = T } 
      }else{
        s1 <- pars[1] ; a1 <- pars[2]; s2 <- pars[3] ; a2 <- pars[4] ; s0 <- 1 - (s1 + s2)
        if((s1 < 0) | (s1 > 1) | (s2 < 0) | (s2 > 1) | ((s1 + s2) > 1) | (a2 < a1)  | (a1 < mina) | (a2 > maxa)){ paramsOK <- F }else{ paramsOK = T } 
      }
    }else{
      if (nSpatStructs == 1){
        s1 <- 1 ; a1 <- pars[1] ; s2 <- 0 ; a2 <- a1 ; s0 <- 0
        if((s1 < 0) | (s1 > 1) | (a1 < mina) | (a1 > maxa)){ paramsOK <- F }else{ paramsOK = T } 
      }else{
        s1 <- pars[1] ; a1 <- pars[2]; s2 <- 1 - s1 ; a2 <- pars[3] ; s0 <- 0
        if((s1 < 0) | (s1 > 1) | (s2 < 0) | (s2 > 1) | ((s1 + s2) > 1) | (a2 < a1)  | (a1 < mina) | (a2 > maxa)){ paramsOK <- F }else{ paramsOK = T } 
      }
    }
    
    parsOut <- paste0('s1=' , round(s1, digits = 3) , ', a1=' , round(a1, digits = 2))
    if(nSpatStructs == 2){ parsOut <- paste0(parsOut , ', s2=' , round(s2, digits = 3) , ', a2=' , round(a2, digits = 2)) }else{}

### default lists to be returned when some error arises...
    if(forCompLik){
      if(returnAll){
        listOut <- list('zinvAz' = NA , 'zinvAX' = NA , 'XinvAX' = NA , 'lndetA' = NA , 'parsOut' = parsOut , 
                        'invAX' = NA , 'invAz' = NA , 'iA' = NA , 'onesinvAones' = NA , 'onesinvAz' = NA , 'onesinvAX' = NA)
      }else{
        listOut <- list('zinvAz' = NA , 'zinvAX' = NA , 'XinvAX' = NA , 'lndetA' = NA , 'parsOut' = parsOut , 
                            'onesinvAones' = NA , 'onesinvAz' = NA , 'onesinvAX' = NA)
      }
    }else{
      if(returnAll){
        listOut <- list('nll' = badNll , 'sigma2hat' = NA , 'betahat' = NA , 'vbetahat' = NA , 'pars' = pars , 'covParams4Kriging' = NA , 'C' = NA , 'iC' = NA)
      }else{             
#        listOut <- NA
#        listOut <- 9E99
        listOut <- badNll
      }
    }
    
    if(paramsOK){
        C <- defineCLMM2(c0 = s0 , c1 = s1 , a1 = a1 , c2 = s2 , a2 = a2 , D = D , covModel = covModel , nSpatStructs = nSpatStructs)
        
        oneszX <- cbind(matrix(1 , N , 1) , z , X)

        # return(list(C,oneszX))
        
        tmp <- lndetANDinvCb(C , oneszX)
         
        if(returnAll){ cholA <- chol(C) }else{}
        lndetC <- tmp$lndetC
        if(is.na(lndetC)){             
            printNll(nll = NA , parsOut = parsOut , verbose = verbose)
            return(listOut) 
        }else{}
        iConeszX <- tmp$invCb

        oneszXinvConeszX <- t(oneszX) %*% iConeszX
        
        onesinvCones <- oneszXinvConeszX[1,1 , drop = FALSE]
        onesinvCz <- oneszXinvConeszX[1,2 , drop = FALSE]
        onesinvCX <- oneszXinvConeszX[1,-c(1,2) , drop = FALSE]
        
        zinvCz <- oneszXinvConeszX[2,2 , drop = FALSE]
        zinvCX <- oneszXinvConeszX[2,-c(1,2) , drop = FALSE]
        XinvCX <- oneszXinvConeszX[-c(1,2),-c(1,2) , drop = FALSE]

        if(forCompLik){
            if(returnAll){
                invCX <- iConeszX[,-c(1,2) , drop = FALSE]
                invCz <- iConeszX[,2 , drop = FALSE]
                return(list('zinvAz' = zinvCz , 'zinvAX' = zinvCX , 'XinvAX' = XinvCX , 'lndetA' = lndetC , 'parsOut' = parsOut , 
                        'invAX' = invCX , 'invAz' = invCz , 'iA' = chol2inv(cholA) , 
                        'onesinvAones' = onesinvCones , 'onesinvAz' = onesinvCz , 'onesinvAX' = onesinvCX))        
            }else{
                return(list('zinvAz' = zinvCz , 'zinvAX' = zinvCX , 'XinvAX' = XinvCX , 'lndetA' = lndetC , 'parsOut' = parsOut , 
                            'onesinvAones' = onesinvCones , 'onesinvAz' = onesinvCz , 'onesinvAX' = onesinvCX))        
            }
        }else{}
        
        tmp <- lndetANDinvCb(XinvCX , t(zinvCX))
        lndetXinvCX <- tmp$lndetC
        if(is.na(lndetXinvCX)){     
          printNll(nll = NA , parsOut = parsOut , verbose = verbose)
          return(listOut) 
        }else{}
        betahat <- tmp$invCb
        
        resiCres <- as.numeric(zinvCz - 2 * zinvCX %*% betahat + t(betahat) %*% XinvCX %*% betahat)
        
        if (optionML){
            sigma2hat <- resiCres / N
            lndetXinvCX <- 0
            np <- 0 # not really, but just because it appears in const of nll
        }else{
            sigma2hat <- resiCres / (N - np)
            lndetXinvCX <- lndetXinvCX - np * log(sigma2hat)
        }        

        if(returnAll){ 
          XiAX <- XinvCX
          XiAz <- t(zinvCX)
          ziAz <- zinvCz
          vbetahat <- sigma2hat * chol2inv(chol(XinvCX)) 
        }else{}
        lndetC <- lndetC + N * log(sigma2hat)
        resiCres <- resiCres / sigma2hat
        
        nll <- -0.5 * (N-np) * log(2 * pi) - 0.5 * lndetC - 0.5 * lndetXinvCX - 0.5 * resiCres
        nll <- -nll
    }else{
        return(listOut)
    }   

    printNll(nll = nll , parsOut = parsOut , verbose = verbose)
    
    if(is.na(nll)){ nll <- badNll }else{}

    if(returnAll){
      covParams4Kriging <- c((1 - s1) * sigma2hat , s1 * sigma2hat , a1 , s2 * sigma2hat , a2)
      if(attachBigMats){
        C <- sigma2hat * C    
        iC <- chol2inv(cholA) / sigma2hat
        return(list('nll' = nll , 'sigma2hat' = sigma2hat , 'betahat' = betahat , 'vbetahat' = vbetahat , 'pars' = pars , 
                    'covModel' = covModel , 'covParams4Kriging' = covParams4Kriging , 'nSpatStructs' = nSpatStructs , 
                    'C' = C , 'iC' = iC , 'XiAX' = XiAX , 'XiAz' = XiAz , 'ziAz' = ziAz))
      }else{
        return(list('nll' = nll , 'sigma2hat' = sigma2hat , 'betahat' = betahat , 'vbetahat' = vbetahat , 'pars' = pars , 
                    'covModel' = covModel , 'covParams4Kriging' = covParams4Kriging , 'nSpatStructs' = nSpatStructs , 
                    'XiAX' = XiAX , 'XiAz' = XiAz , 'ziAz' = ziAz))
      }
    }else{             
        return(nll)
    }
}

#############################################################################
### a composite likelihood approximation, based on blocks, 
### correlation within blocks, but each block assumed independent of each other 
### blocks defined in a list, the jth element of the list having 
### the indices of the data to be used in block j
### possibly two levels of blocks, bricked, so that both sets include all data 
### (ie each data point appears once in level 1 blocks and once in level 2 blocks)
#############################################################################
compLikLMM2 <- function(pars , c , z , X , DBlocks , blocks , nBlocks1 = length(blocks) , covModel , nSpatStructs , incNugget = TRUE , returnAll = FALSE , optionML = FALSE , verbose = TRUE , mina = NULL , maxa = NULL , partOfBiggerModel = F , badNll = 9E99){

    n <- dim(X)[[1]]
    p <- dim(X)[[2]]
       
    c <- as.matrix(c , nrow = n) 

    nLevels <- checkBlocks(blocks , nBlocks1 , n)

    sumlndetA <- sumziAz <- 0
    sumziAX <- matrix(0 , 1 , p)
    sumXiAX <- matrix(0 , p , p)

### for the approach with a station-specific random effect for Robin's data...
### only relevant if partOfBiggerModel is TRUE...
    sum1iA1 <- 0
    sum1iAz <- 0
    sum1iAX <- matrix(0 , 1 , p)
    sumln1PLUironesiCones <- 0 
       
    if (returnAll){
        iCz <- matrix(NA , n , nLevels)
        iCX <- list()
        iCX[[1]] <- matrix(NA , n , p)
        iCBlocks <- list()
        if(nLevels == 2){
            iCX[[2]] <- iCX[[1]]
        }else{}
    }else{}

    if(is.null(maxa)){
        maxa <- max(unlist(lapply(DBlocks , max)))
    }else{}
    
    if(is.null(mina)){
        fnTmp <- function(D){ min(D + 9E99 * diag(dim(D)[[1]])) }
        mina <- unlist(lapply(DBlocks , fnTmp))
        mina <- min(mina[mina > 0])
    }else{}
    
    paramsOK <- T

    for (i in 1:length(blocks)){    
      if(paramsOK){  
          tmp <- fLMM2(pars = pars , c = c[blocks[[i]],] , z = z[blocks[[i]]] , X = X[blocks[[i]],] , D = DBlocks[[i]] , 
                    covModel = covModel , nSpatStructs = nSpatStructs , incNugget = incNugget , returnAll = returnAll , optionML = optionML , verbose = T , forCompLik = TRUE , mina = mina , maxa = maxa  , badNll = badNll)

          if(i == 1){ 
            parsOut <- tmp$parsOut 
            
### will always come into this bit first time around, so now set default return list...                
            if(partOfBiggerModel){
              listOut <- list('sumlndetA' = NA , 'sumziAz' = NA , 'sumziAX' = NA , 'sumXiAX' = NA , 
                    'sum1iA1' = NA , 'sum1iAz' = NA , 'sum1iAX' = NA)
            }else{
              if(returnAll){
                listOut <- list('nll' = badNll , 'sigma2hat' = NA , 'betahat' = NA , 'vbetahat' = NA , 
                    'iCX' = NA , iCz = NA , 'iCres' = NA , 'iCBlocks' = NA , 'pars' = pars , 'covParams4Kriging' = NA)
              }else{
#                listOut <- NA
#                listOut <- 9E99
                listOut <- badNll
              }
            }
          }else{}

### what level is this block?...
          if(i <= nBlocks1){
            iLevelThis <- 1
          }else{
            iLevelThis <- 2
          }           

          if(!is.na(tmp$lndetA)){
            sumlndetA <- sumlndetA + tmp$lndetA
            sumziAz <- sumziAz + tmp$zinvAz
            sumziAX <- sumziAX + tmp$zinvAX
            sumXiAX <- sumXiAX + tmp$XinvAX
            
### for the approach with a station-specific random effect for Robin's data...
            if(partOfBiggerModel){
                sum1iA1 <- sum1iA1 + tmp$onesinvAones
                sum1iAz <- sum1iAz + tmp$onesinvAz
                sum1iAX <- sum1iAX + tmp$onesinvAX
            }else{}
            
            if(returnAll){
              iCBlocks[[i]] <- tmp$iA
              iCX[[iLevelThis]][blocks[[i]],] <- tmp$invAX
              iCz[blocks[[i]] , iLevelThis] <- tmp$invAz
              if(i == 1){ covParams4Kriging <- tmp$covParams4Kriging }else{}
            }else{}
          }else{
            paramsOK <- F
          }
      }else{}
    }
            
### dividing by nLevels here so that lieklihood of data (not nLevels reps of it) will be approximated...      
    sumlndetA <- sumlndetA / nLevels
    sumziAz <- sumziAz / nLevels
    sumziAX <- sumziAX / nLevels
    sumXiAX <- sumXiAX / nLevels
     
    sum1iA1 <- sum1iA1 / nLevels
    sum1iAz <- sum1iAz / nLevels
    sum1iAX <- sum1iAX / nLevels
     
    if(partOfBiggerModel){
      if(paramsOK){
        return(list('sumlndetA' = sumlndetA , 'sumziAz' = sumziAz , 'sumziAX' = sumziAX , 'sumXiAX' = sumXiAX , 
                    'sum1iA1' = sum1iA1 , 'sum1iAz' = sum1iAz , 'sum1iAX' = sum1iAX))
      }else{
        return(listOut)
      }
    }else{}     
     
      
    if(paramsOK){
        tmp <- lndetANDinvCb(sumXiAX , t(sumziAX))
        betahat <- as.matrix(tmp$invCb)
        lndetXiAX <- tmp$lndetC

#        save(sumXiAX , file = "/scratch/rsc6/ortont/Data/robinData/output/sumXiAX.RData")
#        save(sumziAX , file = "/scratch/rsc6/ortont/Data/robinData/output/sumziAX.RData")
#        print(lndetXiAX)
#        print(dim(sumziAX))
#        print(dim(betahat))

        if(is.na(lndetXiAX)){
          printNll(nll = NA , parsOut = parsOut , verbose = verbose)
          return(listOut)
        }else{}        

        resiAres <- as.numeric(sumziAz - 2 * sumziAX %*% betahat + t(betahat) %*% sumXiAX %*% betahat)
        
        if(optionML){
            sigma2hat <- resiAres / n
        }else{
            sigma2hat <- resiAres / (n - p)
        }
            
        if(returnAll){
            vbetahat <- sigma2hat * solve(sumXiAX) 
            iCz <- iCz / sigma2hat

            iCres <- NA * iCz
            for(iLevel in 1:nLevels){ 
                iCX[[iLevel]] <- iCX[[iLevel]] / sigma2hat 
                iCres[,iLevel] <- iCz[,iLevel] - iCX[[iLevel]] %*% betahat
            }
            
            for (i in 1:length(blocks)){ iCBlocks[[i]] <- iCBlocks[[i]] / sigma2hat }

        }else{}

        if (optionML){
            nll <- 0.5 * (n * log(2 * pi) + n * log(sigma2hat) + sumlndetA + n) # [last term because (1 / sigma2hat) * resiAres  = n]
        }else{
            nll <- 0.5 * ((n - p) * log(2 * pi) + (n - p) * log(sigma2hat) + sumlndetA + lndetXiAX + n - p) # [last term because (1 / sigma2hat) * resiAres  = n - p]
        }
    }else{
#        nll <- NA
#        nll <- 9E99
        nll <- badNll
    }

    printNll(nll = nll , parsOut = parsOut , verbose = verbose)
    
    if(is.na(nll)){ nll <- badNll }else{}

    if(returnAll){
      if(paramsOK){
        return(list('nll' = nll , 'sigma2hat' = sigma2hat , 'betahat' = betahat , 'vbetahat' = vbetahat , 
                    'iCX' = iCX , iCz = iCz , 'iCres' = iCres , 'iCBlocks' = iCBlocks , 
                    'pars' = pars , 'covParams4Kriging' = covParams4Kriging))
      }else{
        return(list('nll' = badNll , 'sigma2hat' = NA , 'betahat' = NA , 'vbetahat' = NA , 
                    'iCX' = NA , iCz = NA , 'iCres' = NA , 'iCBlocks' = NA , 
                    'pars' = pars , 'covParams4Kriging' = NA))
      }
    }else{             
        return(nll)
    }
    
}    


krigingLMM2 <- function(ck , Xk , c , z , X , D , covModel , nSpatStructs , covParams , betahat , predVarsRqd , fullCovMtx = FALSE , blockLength = 0){

    ndimc <- ncol(c)
    nk <- nrow(ck)
    if(is.null(ndimc)){
      ndimc <- 1
      nk <- length(ck)
    }else{}
    
    ### prediction block centred at 0,0
    if(blockLength > 0){
      if(ncol(Xk) > 1){ stop('Generalise this bit of krigingLMM2 for block kriging with trend!') }else{}
      if(ndimc>2){ stop('Generalise this bit of krigingLMM2 for 3d blocks!') }else{}
      xBlock0 <- seq(0,blockLength,length=5)
      if(ndimc == 1){
        xBlock0 <- matrix(xBlock0 , ncol = 1)
      }else if(ndimc == 2){
        xBlock0 <- cbind(rep(xBlock0 , each = length(xBlock0)) , rep(xBlock0 , length(xBlock0)))
      }else{}
      xBlock0 <- xBlock0 - blockLength / 2
    }else{
      xBlock0 <- matrix(0 , 1 , ndimc)
    }
    nkPerBlock <- nrow(xBlock0)
    
    N <- length(z)
    np <- ncol(X)
    if(is.null(np)){ np <- 1 }else{}
    
    c0 <- covParams[1] ; c1 <- covParams[2] ; a1 <- covParams[3] ; 
    if (nSpatStructs == 1){
        c2 <- 0 ; a2 <- a1
    }else{
        c2 <- covParams[4] ; a2 <- covParams[5]
    }

    print(paste0('c0 = ' , c0))
    print(paste0('c1 = ' , c1))
    print(paste0('a1 = ' , a1))
    if (nSpatStructs == 2){
      print(paste0('c2 = ' , c2))
      print(paste0('a2 = ' , a2))
    }else{}
    
    C <- defineCLMM2(c0 = c0 , c1 = c1 , a1 = a1 , c2 = c2 , a2 = a2 , D = D , covModel = covModel , nSpatStructs = nSpatStructs)

    zX <- cbind(z , X)
    tmp <- lndetANDinvCb(C , zX)

    lndetC <- tmp$lndetC
    invCz <- tmp$invCb[,1 , drop = FALSE]
    invCX <- tmp$invCb[,2:(1+np) , drop = FALSE]
      
    # invC <- chol2inv(tmp$cholC)
    invC <- solve(C)
    invCRes <- invCz - invCX %*% betahat    
        
    XinvCX <- t(X) %*% invCX
    invXinvCX <- solve(XinvCX)
        
    if(blockLength == 0){
      nPerChunk <- 5000
      nChunks <- ceiling(nk / nPerChunk)
    }else{
      nPerChunk <- 1
      nChunks <- nk
    }

    pred <- matrix(NA,nk,1)
    predVars <- matrix(NA,nk,1)

    if(fullCovMtx){ 
        if(nChunks > 1){
            stop('Returning full prediction covariance matrix only possible if less than 5000 prediction locations given!')
        }else{}
        predVars <- matrix(NA,nk,nk)
    }else{}

###############################################
### for each chunk...
###############################################
	for (iChunk in 1:nChunks){
    
	  if(is.element(iChunk , round(quantile(seq(nChunks) , seq(0.05 , 0.95 , 0.05))))){ print(paste0(iChunk , ' of ' , nChunks)) }else{}

		if (iChunk < nChunks){
			ikThis <- seq( (iChunk - 1) * nPerChunk + 1 , iChunk * nPerChunk)	
		}else{
			ikThis <- seq( (iChunk - 1) * nPerChunk + 1 , nk)
		}
    nkThis <- length(ikThis)
    
    if(ndimc == 1){
      ckThis <- ck[ikThis]
    }else{
      ckThis <- ck[ikThis,  , drop = FALSE]
    }

    if(blockLength > 0){
      if(ndimc == 1){
        ckThis <- ckThis + xBlock0
      }else{
        ckTmp <- xBlock0
        for(j in 1:ndimc){
          ckTmp[,j] <- ckTmp[,j] + ckThis[1,j]
        }
        ckThis <- ckTmp 
        rm(ckTmp)
      }
    }else{}

 		Dkh <- rdist(ckThis,c)
    Ckh <- defineCLMM2(c0 = c0 , c1 = c1 , a1 = a1 , c2 = c2 , a2 = a2 , D = Dkh , covModel = covModel , nSpatStructs = nSpatStructs)

    if(blockLength == 0){
      pred[ikThis] <- Xk[ikThis, , drop = FALSE] %*% betahat + Ckh %*% invCRes
    }else{
### needs to be generalised in trend is not a const mean.      
      pred[ikThis] <- Xk[ikThis, , drop = FALSE] %*% betahat + mean(Ckh %*% invCRes)
    }
    
### also calculate the prediction variances; but only if rqd...
    if(predVarsRqd){
      if(fullCovMtx | (blockLength > 0)){
        Dk <- rdist(ckThis,ckThis)
        Ck <- defineCLMM2(c0 = c0 , c1 = c1 , a1 = a1 , c2 = c2 , a2 = a2 , D = Dk , covModel = covModel , nSpatStructs = nSpatStructs)
        Xk_CkhiChXh <- Xk[ikThis,] - Ckh %*% invCX
        varTrend <- (Xk_CkhiChXh %*% invXinvCX) %*% t(Xk_CkhiChXh)
        if(blockLength == 0){
          predVars[ikThis,ikThis] <- Ck - (Ckh %*% invC) %*% t(Ckh) + varTrend
        }else{
          predVars[ikThis] <- mean(Ck - (Ckh %*% invC) %*% t(Ckh) + varTrend)
        }
      }else{
        Ck <- c0 + c1 + c2
        Xk_CkhiChXh <- Xk[ikThis,] - Ckh %*% invCX
        varTrend <- rowSums((Xk_CkhiChXh %*% invXinvCX) * Xk_CkhiChXh)
        predVars[ikThis] <- Ck - rowSums((Ckh %*% invC) * Ckh) + varTrend
      }
    }
  } 
    
  return(list('pred' = pred , 'predVars' = predVars))
}

###################################################################
### for the compLik approach...
### only for point-support predoiction.
###################################################################
compLikKrigingLMM2 <- function(ck , Xk , blocksk , c , z , X , blocks , DBlocks , nBlocks1 = length(blocks) , iCres , iCX , vbetahat , covModel , nSpatStructs , covParams , betahat , nsdsPlot){
    muk <- Xk %*% betahat

    X <- as.matrix(X)
    Xk <- as.matrix(Xk)

    c <- as.matrix(c) 
    ck <- as.matrix(ck) 

    n <- dim(X)[[1]]
    p <- dim(X)[[2]]

    nLevels <- checkBlocks(blocks , nBlocks1 , n)
      
    nk <- dim(ck)[1]

    pred <- predVars <- predVarsTrend <- predVarsSK <- matrix(NA , nk , nLevels)

    if(length(covParams) == 3) { covParams <- c(covParams , 0 , NA) }else{}
    if(length(covParams) != 5) { stop('Error! covParams not correctly inputted for kriging routine!') }else{}
    
    Ck <- covParams[1] + covParams[2] + covParams[4]
    for (i in 1:length(blocks)){  
    
### what level is this block?...
        if(i <= nBlocks1){
            iLevelThis <- 1
        }else{
            iLevelThis <- 2
        }           
    
        DkhThis <- rdist(ck[blocksk[[i]],] , c[blocks[[i]],])
        CkhThis <- defineCLMM2(c0 = covParams[1] , c1 = covParams[2] , a1 = covParams[3] , c2 = covParams[4] , a2 = covParams[5] , 
              D = DkhThis , covModel = covModel , nSpatStructs = nSpatStructs)

        ChThis <- defineCLMM2(c0 = covParams[1] , c1 = covParams[2] , a1 = covParams[3] , c2 = covParams[4] , a2 = covParams[5] , 
              D = DBlocks[[i]] , covModel = covModel , nSpatStructs = nSpatStructs)

        tmp <- lndetANDinvCb(ChThis , t(CkhThis))
        iChChkThis <- tmp$invCb

        pred[blocksk[[i]] , iLevelThis] <- muk[blocksk[[i]]] + CkhThis %*% iCres[blocks[[i]],iLevelThis]
        vkSK <- Ck - rowSums(CkhThis * t(iChChkThis))

        Xk_CkhiCXThis <- Xk[blocksk[[i]],] - CkhThis %*% iCX[[iLevelThis]][blocks[[i]],]
        tmp <- vbetahat %*% t(Xk_CkhiCXThis)
        predVars[blocksk[[i]] , iLevelThis] <- vkSK + rowSums(Xk_CkhiCXThis * t(tmp)) 

        predVarsSK[blocksk[[i]] , iLevelThis] <- vkSK

        tmp <- vbetahat %*% t(Xk[blocksk[[i]],])
        predVarsTrend[blocksk[[i]] , iLevelThis] <- rowSums(Xk[blocksk[[i]],] * t(tmp)) 
    }

##################################################################
### if we have two levels, then our estimate is N(zk;m1,sig1) * N(zk;m2,sig2)
### this product gives a N(zk;m3,sig3) distribution (up to constant terms) with 
###   m3 = inv(inv[sig1] + inv[sig2]) (inv[sig1] m1 + inv[sig2] m2)   
###   sig3 = inv(inv[sig1] + inv[sig2]) 
### see matrix cookbook.
##################################################################
    if(nLevels == 2){
        tmp1 <- rowSums(pred / predVars)
        
        predVars <- 1 / rowSums(1 / predVars)
        pred <- predVars * tmp1

        predVars <- 2 * predVars
### not sure about these 2 for 2 levels, but...        
        predVarsSK <- 2 / rowSums(1 / predVarsSK)
        predVarsTrend <- 2 / rowSums(1 / predVarsTrend)
    }else{}

### for plotting, make polygons that can be shaded to show prediction intervals...
    polyPredVarsa <- cbind(ck , pred - nsdsPlot * sqrt(predVars)) 
    polyPredVarsb <- cbind(ck , pred + nsdsPlot * sqrt(predVars)) 
    polyPredVarsb <- polyPredVarsb[nrow(polyPredVarsb):1,]
    polyPredVars <- rbind(polyPredVarsa , polyPredVarsb , polyPredVarsa[1,])

    polyTrendVarsa <- cbind(ck , muk - nsdsPlot * sqrt(predVarsTrend)) 
    polyTrendVarsb <- cbind(ck , muk + nsdsPlot * sqrt(predVarsTrend)) 
    polyTrendVarsb <- polyTrendVarsb[nrow(polyTrendVarsb):1,]
    polyTrendVars <- rbind(polyTrendVarsa , polyTrendVarsb , polyTrendVarsa[1,])

    return(list('pred' = pred , 'predVars' = predVars , 'predTrend' = muk , 'predVarsTrend' = predVarsTrend , 'predVarsSK' = predVarsSK , 
                'polyPredVars' = polyPredVars , 'polyTrendVars' = polyTrendVars))
}

####################################################
### cross-validation, given estimated (trend and covariance) parameters 
### if iSubsets is given, it is a list of subsets 
### for each of which a prediction of the subset's average is required
### for the final subset only, the full prediction vector and covariance matrix is also returned
####################################################
XVLMM2 <- function(c , z , X , iC , betahat , vbetahat , nsdsPlot , iSubsets = NULL){ 

######################################################
### Write C as block matrix: C = [Chh , Chk ; Ckh , Ckk] (k is single prediction location, h is all other data)
### Write its inverse in block form as iC = [iChh , iChk ; iCkh , iCkk]
### Then from equations for inverse of partitioned matrix:
###             iChh = inv[Chh - Chk inv(Ckk) Ckh]
###             iChk = -inv[Chh] Chk inv[Ckk - Ckh inv(Chh) Chk]
###             iCkh = -inv[Ckk] Ckh inv[Chh - Chk inv(Ckk) Ckh]
###             iCkk = inv[Ckk - Ckh inv(Chh) Chk]
### For cross validation, given C and betahat with uncertainty vbetahat, we need...
###    zk = muk + Ckh inv[Chh] (zh - muh)
###    vk = Ckk - Ckh inv[Chh] Chk + (Xk - Ckh inv[Chh] Xh) vbetahat t(Xk - Ckh inv[Chh] Xh)
### From our four equations, we have:
###    iChk  inv[iCkk] = -inv[Chh] Chk inv[Ckk - Ckh inv(Chh) Chk] [Ckk - Ckh inv(Chh) Chk]
###                = -inv[Chh] Chk
### So:
###    inv[Chh] Chk = -iChk inv[iCkk]
### Also, the SK variance (the first two terms of vk) are given by:
###    Ckk - Ckh inv[Chh] Chk = inv[iCkk]
### Other terms are calculated from these two terms.
#############################################################
    c <- as.matrix(c) 
    X <- as.matrix(X)
    mu <- X %*% betahat
    res <- z - mu

    n <- dim(X)[[1]]
    if(is.null(iSubsets)){
        nSubsets <- n
        nPerSubset <- integer(n) + 1
    }else{
### a list whose jth element contains the indices of data to be removed in the jth iteration
### and for which an average will be required...
        nSubsets <- length(iSubsets)
        nPerSubset <- unlist(lapply(iSubsets , length))
    }
    
    zk <- vk <- zbar <- matrix(NA,nSubsets,1)
    for (i in 1:nSubsets){
        if(is.null(iSubsets)){
            iThis <- i
            nThis <- 1
            zbar[i] <- z[iThis]
        }else{
            iThis <- iSubsets[[i]]
            nThis <- length(iThis)
            zbar[i] <- mean(z[iThis])
        }

        if(length(iThis) > 0){
          CkhiChh <- -solve(matrix(iC[iThis,iThis] , nThis , nThis) , matrix(iC[iThis,-iThis] , nThis , n - nThis)) 

          tmp <- CkhiChh %*% cbind(res[-iThis] , X[-iThis,])
          CkhiChhzh <- tmp[,1]
          CkhiChhXh <- tmp[,-1]
        
          Xk_CkhiChhXh <- matrix(X[iThis,] - CkhiChhXh , nrow = nThis)
          vSK <- solve(iC[iThis,iThis])

          zkFull <- mu[iThis] + CkhiChhzh
          vkFull <- vSK + Xk_CkhiChhXh %*% vbetahat %*% t(Xk_CkhiChhXh)
        
          zk[i] <- mean(zkFull)
          vk[i] <- mean(vkFull)
        }else{}
    }
      
    err <- (zbar - zk)
    serr <- err ^ 2
    sse <- serr / vk

### for plotting, make polygons that can be shaded to show prediction intervals...
### only if point-support predictions
    if((is.null(iSubsets)) & (dim(c)[[2]] == 1)){
        polyvka <- cbind(c , zk - nsdsPlot * sqrt(vk)) 
        polyvkb <- cbind(c , zk + nsdsPlot * sqrt(vk)) 
        polyvkb <- polyvkb[nrow(polyvkb):1,]
        polyvk <- rbind(polyvka , polyvkb , polyvka[1,])
    }else{
        polyvk <- NA
    }
    
    return(list('zk' = zk , 'vk' = vk , 'err' = err , 'serr' = serr , 'sse' = sse , 'polyvk' = polyvk , 
                'zbar' = zbar , 'nPerSubset' = nPerSubset , 'zkFullFinal' = zkFull , 'vkFullFinal' = vkFull))
}

####################################################
### cross-validation based on the composite likelihood, given estimated (trend and covariance) parameters 
####################################################
compLikXVLMM2 <- function(c , z , X , blocks , iCBlocks , nBlocks1 = length(blocks) , betahat , vbetahat , nsdsPlot){ 
    c <- as.matrix(c) 
    X <- as.matrix(X)

    n <- dim(X)[[1]]

    nLevels <- checkBlocks(blocks , nBlocks1 , n)

    zk <- vk <- matrix(NA , n , nLevels)
    for (i in 1:length(blocks)){  
### what level is this block?...
        if(i <= nBlocks1){
            iLevelThis <- 1
        }else{
            iLevelThis <- 2
        }           
    
### for 'point-support' prediction
        tmp <- XVLMM2(c = c[blocks[[i]],], z = z[blocks[[i]]] , X = X[blocks[[i]],] , iC = iCBlocks[[i]] , betahat = betahat , vbetahat = vbetahat , nsdsPlot = nsdsPlot) 

        zk[blocks[[i]],iLevelThis] <- tmp$zk 
        vk[blocks[[i]],iLevelThis] <- tmp$vk
    }

    zkAllLevels <- zk
    vkAllLevels <- vk

##################################################################
### if we have two levels, then our estimate is [N(zk;m1,sig1) * N(zk;m2,sig2)] ^ (1 / nLevels)
### this product gives a N(zk;m3,sig3) distribution (up to constant terms) with 
###   m3 = inv(inv[sig1] + inv[sig2]) (inv[sig1] m1 + inv[sig2] m2)   
###   sig3 = inv(inv[sig1] + inv[sig2]) 
### to include the power, var = nLevels * sig3 
### see matrix cookbook.
##################################################################
    if(nLevels == 2){
        iMin <- cbind(seq(n) , apply(vk , 1 , which.min))
        vk.Min <- vk[iMin]
        zk.Min <- zk[iMin]
        
        tmp1 <- rowSums(zk / vk)
        
        vk <- 1 / rowSums(1 / vk)
        zk <- vk * tmp1
        
        vk <- nLevels * vk # because based on likelihood for two copies of data.
    }else{
        zk.Min <- vk
        vk.Min <- vk
    }

    err <- (z - zk)
    serr <- err ^ 2
    sse <- serr / vk

### for plotting, make polygons that can be shaded to show prediction intervals...
    if(dim(c)[[2]] == 1){
        polyvka <- cbind(c , zk - nsdsPlot * sqrt(vk)) 
        polyvkb <- cbind(c , zk + nsdsPlot * sqrt(vk)) 
        polyvkb <- polyvkb[nrow(polyvkb):1,]
        polyvk <- rbind(polyvka , polyvkb , polyvka[1,])
    }else{
        polyvk <- NA
    }
    
    return(list('zk' = zk , 'vk' = vk , 'err' = err , 'serr' = serr , 'sse' = sse , 'polyvk' = polyvk , 'zkAllLevels' = zkAllLevels , 'vkAllLevels' = vkAllLevels , 'zk.Min' = zk.Min , 'vk.Min' = vk.Min , 'zbar' = z))

}

####################################################
### cross-validation based on the composite likelihood, given estimated (trend and covariance) parameters 
### but for a complete subset, for which a prediction of the mean of the entire subset is required. 
####################################################
compLikXVLMM2.Subsets <- function(c , z , X , blocks , iCBlocks , nBlocks1 = length(blocks) , betahat , vbetahat , nsdsPlot , iSubsets){ 
    c <- as.matrix(c) 
    X <- as.matrix(X)
    
    n <- dim(X)[[1]]
    nLevels <- checkBlocks(blocks , nBlocks1 , n)
    nSubsets <- length(iSubsets)
    
    zbar <- zk <- vk <- zk.Min <- vk.Min <- NA * numeric(nSubsets)
    for (isub in 1:nSubsets){
      if((isub %% 100) == 0){
          print(paste0('Cross validating for subset ' , isub , ' of ' , nSubsets , '...'))
      }else{}
      
      iSubset <- iSubsets[[isub]]
      nInSubset <- length(iSubset)

      if(nInSubset > 0){
### for the mulivariate prediction at all locations in subset, store for each level...
        zkFull <- matrix(0 , nInSubset , nLevels)
        vkFull <- list()
        for(i in 1:nLevels){ vkFull[[i]] <- matrix(0 , nInSubset , nInSubset) }

### which blocks overlap with this subset?...
        iBlocksThis <- lapply(blocks , intersect , iSubset)
        iBlocksThis <- which(unlist(lapply(iBlocksThis , length)) > 0)

        for (i in iBlocksThis){
   
### what level is this block?...
          if(i <= nBlocks1){
            iLevelThis <- 1
          }else{
            iLevelThis <- 2
          }           

### to which points will this block contribute?...
          iPtsThisInFull <- which(is.element(iSubset , blocks[[i]]))
        
          if(length(iPtsThisInFull) > 0){
            iPtsThis <- iSubset[iPtsThisInFull]
        
### what are the indices of these points within this block?...
### renumber non-empty iSubsetsThis to go from 1 - nDataThisBlock
            iminThisBlock <- min(blocks[[i]])
          
            iSubsetThis <- list()
            iSubsetThis[[1]] <- iPtsThis - iminThisBlock + 1
    
            tmp <- XVLMM2(c = c[blocks[[i]],], z = z[blocks[[i]]] , X = X[blocks[[i]],] , iC = iCBlocks[[i]] , betahat = betahat , vbetahat = vbetahat , nsdsPlot = nsdsPlot , iSubset = iSubsetThis) 

            zkFull[iPtsThisInFull,iLevelThis] <- tmp$zkFullFinal
            vkFull[[iLevelThis]][iPtsThisInFull,iPtsThisInFull] <- tmp$vkFullFinal
          }else{}
        }

##################################################################
### if we have two levels, then our estimate is [N(zk;m1,sig1) * N(zk;m2,sig2)] ^ (1/nLevels)
### this product gives a N(zk;m3,sig3) distribution (up to constant terms) with 
###   m3 = inv(inv[sig1] + inv[sig2]) (inv[sig1] m1 + inv[sig2] m2)   
###   sig3 = inv(inv[sig1] + inv[sig2])
### to include the power, var = nLevels * sig3 
### see matrix cookbook.
##################################################################
        zkFullFinal.AllLevels <- zkFull
        vkFullFinal.AllLevels <- vkFull
        if(nLevels == 2){
          tmp1 <- lndetANDinvCb(vkFull[[1]] , zkFull[,1])
          tmp2 <- lndetANDinvCb(vkFull[[2]] , zkFull[,2])

#          tmp3 <- lndetANDinvCb(chol2inv(tmp1$cholC) + chol2inv(tmp2$cholC) , tmp1$invCb +tmp2$invCb)
### no longer returning also cholC, so:
          m3Tmp <- chol2inv(chol(vkFull[[1]])) + chol2inv(chol(vkFull[[2]]))
          tmp3 <- lndetANDinvCb(m3Tmp , tmp1$invCb +tmp2$invCb)

          zkFull <- tmp3$invCb
#          vkFull <- nLevels * chol2inv(tmp3$cholC) # because based on likelihood for two copies of data
### no longer returning also cholC, so:
          vkFull <- nLevels * chol2inv(chol(m3Tmp)) # because based on likelihood for two copies of data
        
### alt - which level gave block with the best coverage of prediction block?
          n0.L1 <- length(which(vkFullFinal.AllLevels[[1]] == 0))
          n0.L2 <- length(which(vkFullFinal.AllLevels[[2]] == 0))
          if((n0.L1 == 0) & (n0.L2 == 0)){
### both complete, so choose the one with the minimum average variance(ie best prediction)...        
            vTmp <- unlist(lapply(vkFullFinal.AllLevels , mean))
            iMin <- which.min(vTmp)
          }else if(n0.L1 < n0.L2){
            iMin <- 1
          }else{
            iMin <- 2
          }
          zkFull.Min <- zkFullFinal.AllLevels[,iMin]
          vkFull.Min <- vkFullFinal.AllLevels[[iMin]]
        }else{
          vkFull <- vkFull[[1]]
          vkFull.Min <- vkFull[[1]]
          zkFull.Min <- zkFull
        }

##############################################################
### predictions of subset average, and relevant stats...
##############################################################
        zbar[isub] <- mean(z[iSubset])
    
        zk[isub] <- mean(zkFull)
        vk[isub] <- mean(vkFull)
      
        zk.Min[isub] <- mean(zkFull.Min)
        vk.Min[isub] <- mean(vkFull.Min)

      }else{} # close of "if (nInSubset > 0)"

    } # close of loop over subsets.
    
    err <- zbar - zk 
    serr <- err ^ 2
    sse <- serr / vk

##############################################################
### and return, including the full prediction and covariance matrix for just the final subset...    
##############################################################
    return(list('zk' = zk , 'vk' = vk , 'err' = err , 'serr' = serr , 'sse' = sse , 'zbar' = zbar , 
                'zkFullFinal' = zkFull , 'vkFullFinal' = vkFull , 'zkFullFinal.AllLevels' = zkFullFinal.AllLevels , 'vkFullFinal.AllLevels' = vkFullFinal.AllLevels , 
                'zkFull.Min' = zkFull.Min , 'vkFull.Min' = vkFull.Min , 'zk.Min' = zk.Min , 'vk.Min' = vk.Min))

}

####################################################
### cross-validation based on the composite likelihood, given estimated (trend and covariance) parameters 
### but for complete subsets, for which a prediction of the mean of the entire subset is required. 
####################################################
compLikXVLMM2.Subsets.1Level <- function(c , z , X , blocks , iCBlocks , betahat , vbetahat , nsdsPlot , iSubsets){ 
    c <- as.matrix(c) 
    X <- as.matrix(X)

    n <- dim(X)[[1]]
    
    nBlocks <- length(blocks)
    nSubsets <- length(iSubsets)
    nPerSubset <- unlist(lapply(iSubsets , length))
    nLevels <- checkBlocks(blocks , nBlocks , n)
    if(nLevels != 1){ stop('Must have 1 level of blocks for this routine!') }else{}

    zk <- vk <- numeric(nSubsets)
    for (i in 1:length(blocks)){  
### for 'block-support' prediction (bad name, don't confuse with the blocks of composite likelihood)
### XV routine returns 'block' averages, we will use to make block sums then divide at the end...
### only those in this block...    
        
### to which subsets will this block contribute?...
        iSubsetsThis <- lapply(iSubsets , intersect , blocks[[i]])
        iWhichSubsetsThis <- which(unlist(lapply(iSubsetsThis , length)) > 0)
        
### renumber non-empty iSubsetsThis to go from 1 - nDataThisBlock
        iminThisBlock <- min(blocks[[i]])
        iSubsetsThis <- iSubsetsThis[iWhichSubsetsThis]
        iSubsetsThis <- lapply(iSubsetsThis , function(x , c) x - c , iminThisBlock - 1)
    
        tmp <- XVLMM2(c = c[blocks[[i]],], z = z[blocks[[i]]] , X = X[blocks[[i]],] , iC = iCBlocks[[i]] , betahat = betahat , vbetahat = vbetahat , nsdsPlot = nsdsPlot , iSubsets = iSubsetsThis) 

        zk[iWhichSubsetsThis] <- zk[iWhichSubsetsThis] + tmp$zk * tmp$nPerSubset 
        vk[iWhichSubsetsThis] <- vk[iWhichSubsetsThis] + tmp$vk * (tmp$nPerSubset ^ 2)
    }

    zbar <- unlist(lapply(iSubsets , function(iSubsets , z){ if(length(iSubsets) == 0){ return(NA) }else{ return(mean(z[iSubsets])) } } , z = z))
    
##################################################
### zk[i] is now sum[zk predictions in subset i]
### vk[i] is now sum[covariances in subset i]
### so to get predictions and variance of subset average...
##################################################
    iGT0 <- which(nPerSubset > 0)
    zk[iGT0] <- zk[iGT0] / nPerSubset[iGT0]
    vk[iGT0] <- vk[iGT0] / (nPerSubset[iGT0] ^ 2)
    
    err <- zbar - zk 
    serr <- err ^ 2
    sse <- serr / vk

##############################################################
### and return, including the full prediction and covariance matrix for just the final subset...    
##############################################################
    return(list('zk' = zk , 'vk' = vk , 'err' = err , 'serr' = serr , 'sse' = sse , 'zbar' = zbar))

}


##########################################################################
### function to define C for these functions... 
##########################################################################
defineCLMM2 <- function(c0 , c1 , a1 , c2 , a2 , D , covModel , nSpatStructs){

    iDPossErrs <- which((D > 0) & (D <=1E-10))
    if(length(iDPossErrs) > 0){ print('WARNING - SOME SMALL VALUES IN D COULD BE COLOCATED - CHECK DEFINITION OF D!') }else{}

    if(covModel == 'exponential'){ 
        C <-  c0 * (D == 0) + c1 * exp(-3 * D / a1) 
        if(nSpatStructs == 2){ C <- C + c2 * exp(-3 * D / a2) }else{}
    }else if(covModel == 'spherical'){ 
        DOVERa1 <- D / a1
        DOVERa2 <- D / a2
        C1 <- 1 - (1.5 * DOVERa1 - 0.5 * (DOVERa1 ^ 3) )
        C1[which(DOVERa1 > 1)] <- 0
        C <- c0 * (D == 0) + c1 * C1

        if(nSpatStructs == 2){ 
            C2 <- 1 - (1.5 * DOVERa2 - 0.5 * (DOVERa2 ^ 3) )
            C2[which(DOVERa2 > 1)] <- 0
            C <- C + c2 * C2
        }else{}

    }else if (covModel == 'gaussian'){ 
        C <-  c0 * (D == 0) + c1 * exp(-((sqrt(3) * D / a1) ^ 2)) 
        if (nSpatStructs == 2){ C <- C + c2 * exp(-((sqrt(3) * D / a2) ^ 2)) }else{}

    }else if (covModel == 'exponentialgaussian'){
### exp for short range, gaussian for long range...
        if (nSpatStructs != 2){ stop('Error, for exponentialgaussian cov model, must put nSpatStructs = 2') }else{}        
        C <-  c0 * (D == 0) + c1 * exp(-3 * D / a1) 
        C <- C + c2 * exp(-((sqrt(3) * D / a2) ^ 2)) 

    }else if(substr(covModel , 1 , 6) == 'matern'){ 
### note - with nu = 0.5 will be different to exponential, cos of different parameterization. 
### sensible in future to align all parameterizations (st a is effective range)
### although if using in iak work, attention that iaCov functions are based on the distance (not eff range) parameterization
        nu <- as.numeric(substr(covModel , 7 , nchar(covModel)))
        C <-  c0 * (D == 0) + c1 * maternCov4fLMM2(D = D , pars = c(1 , a1 , nu)) 
        if(nSpatStructs == 2){ C <- C + c2 * maternCov4fLMM2(D = D , pars = c(1 , a2 , nu)) }else{}

    }else{}

    return(C)
}

######################################################################
### to wrap up the printout...
######################################################################
printNll <- function(nll , parsOut , verbose = T){
    if(verbose){
        print(paste0('nll=' , round(nll , digits = 3) , '; ' , parsOut)) 
    }else{}
}

######################################################################
### some blocking functions...
######################################################################
getnBlocks1 <- function(blocks , n){
    nPerBlock <- unlist(lapply(blocks , length))
    nBlocks1 <- which(cumsum(nPerBlock) == n)
    return(nBlocks1)
}

checkBlocks <- function(blocks , nBlocks1 , n){

    nBlocks <- length(blocks)
    nBlocks2 <- nBlocks - nBlocks1

    if(nBlocks2 == 0){ nLevels <- 1 }else{ nLevels <- 2 }

    tmp <- unlist(blocks)
    if((length(tmp) == n) & (nLevels == 1)){
### all ok for a single level of blocks.
    }else if((length(tmp) == (2 * n)) & (nLevels == 2)){
### assume first nBlocks1 blocks containing all n data, as do the remaingin nBlocks2 blocks
        nBlocks1Check <- getnBlocks1(blocks , n)
        if((length(nBlocks1Check) != 1) || (nBlocks1Check != nBlocks1)){ 
            stop('Some error in 2-level blocking - first nBlocks1 blocks containing all n data, as do the remaingin nBlocks2 blocks!') 
        }else{
### all ok for a two levels of blocks.
        }
        
#        stop('Not ready for this yet! To be finished!')
        
    }else{
        stop('Some error in blocking - union of blocks must contain each datapoint once or each datapoint twice!')    
    }    
    
    return(nLevels)
}

##################################################################
### the matern covariance function...
##################################################################
maternCov4fLMM2 <- function(D , pars){
	c1 <- pars[1]
	a <- pars[2]
	nu <- pars[3]

	if((c1 > 0) & (a > 0) & (nu >= 0.05)  & (nu <= 20)){

		iD0 <- which(D == 0)
		iDGT0 <- which(D > 0)

### range is approx a * 3...this is from wiki, 
### and is i think what stein's parameterization was supposed to be.
		sqrt2nuOVERa <- sqrt(2 * nu) / a 
		Dsqrt2nuOVERa <- D * sqrt2nuOVERa  

		bes <- 0 * D
		bes[iDGT0] <- besselK(Dsqrt2nuOVERa[iDGT0] , nu) 
      
		lnconstmatern <- NA * Dsqrt2nuOVERa
	    lnconstmatern[iDGT0] <- nu * log(Dsqrt2nuOVERa[iDGT0]) - (nu - 1) * log(2) - lgamma(nu)
      
		realmin <- 3.448490e-304 
		ibesGT0 <- which(bes > realmin)

        C <- 0 * D # initiate.
        C[ibesGT0] <- c1 * exp(lnconstmatern[ibesGT0]+log(bes[ibesGT0]))
        C[iD0] <- c1
        C[which(is.infinite(bes))] <- c1

	}else{
		C <- NA
	}

	return(C)
}

