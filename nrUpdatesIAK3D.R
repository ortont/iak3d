nrUpdatesIAK3D <- function(dA , lnvParInits , W , smallV = NA){

    printNrProgress <- F

    if(is.na(smallV)){
        smallV <- log(var(W[,dim(W)[[2]]]) / 10000)
    }else{}

    lnvPars <- lnvParInits
    tmp <- gradHessIAK3D(dA , lnvPars , W)

    nllNew <- tmp$nll
    gradNew <- tmp$grad
    hessNew <- tmp$hess
    hessNew <- 0.5 * (hessNew + t(hessNew))
    fimNew <- tmp$fim
    fimNew <- 0.5 * (fimNew + t(fimNew))
    betahatNew <- tmp$betahat
    vbetahatNew <- tmp$vbetahat
    CNew <- tmp$C
    
    if(printNrProgress){        
        print('Initial values...')
        print(lnvPars)
        print(nllNew)
    }else{}
    
    tolNr <- 1E-6
    noits <- 0
    nllPrev <- Inf
    while((noits < 20) & ((nllNew < (nllPrev - tolNr)) | (noits < 4))){ # do at least 4 its
### copy 'New' to 'Prev' values...
        nllPrev <- nllNew

### check current hess (or fim)...
        if (noits < 2){
### at least 2 fisher scoring updates...        
#            checkHess <- lndetANDinvCb(fimNew , gradNew)
            checkHess <- try(solve(fimNew , gradNew) , silent = TRUE)
            if(printNrProgress){ print('Fisher-scoring update...') }else{}
        }else{
#            checkHess <- lndetANDinvCb(hessNew , gradNew)
            checkHess <- try(solve(hessNew , gradNew) , silent = TRUE)
            if(is.character(checkHess$cholC)){
### try another fisher-scoring update...            
                if(printNrProgress){ print('Fisher-scoring update...') }else{}
#                checkHess <- lndetANDinvCb(fimNew , gradNew)
                checkHess <- try(solve(fimNew , gradNew) , silent = TRUE)
            }else{
                if(printNrProgress){ print('Newton-Raphson update...') }else{}
            }
        }
        
#        if(!is.character(checkHess$cholC)){
        if(!is.character(checkHess)){
#################################################           
### updates of very small variances can be slow because of flat likelihood.
### speed this up by moving to -20 to compare... 
#################################################           
#            iSmallV <- which((lnvPars < smallV) & (checkHess$invCb > 0))
            iSmallV <- which((lnvPars < smallV) & (checkHess > 0))
            
#            lnvPars <- lnvPars - checkHess$invCb
            lnvPars <- lnvPars - checkHess
            lnvPars[iSmallV] <- -20
            
            tmp <- gradHessIAK3D(dA , lnvPars , W)

            nllNew <- tmp$nll
            gradNew <- tmp$grad
            hessNew <- tmp$hess
            hessNew <- 0.5 * (hessNew + t(hessNew))
            fimNew <- tmp$fim
            fimNew <- 0.5 * (fimNew + t(fimNew))
            betahatNew <- tmp$betahat
            vbetahatNew <- tmp$vbetahat
            CNew <- tmp$C

            if(printNrProgress){        
                print(paste0('After ' , noits + 1 , ' updates...'))
                print(lnvPars)
                print(nllNew)
            }else{}
        }else{
            print('Exiting Newton-Raphson because neither Hessian or FIM are positive definite!')
            return(list('lnvPars' = lnvPars , 'nll' = nllNew , 'grad' = gradNew , 'hess' = hessNew , 'fim' = fimNew , 
                        'betahat' = betahatNew , 'vbetahat' = vbetahatNew , 'C' = CNew , 'noits' = noits , 'errFlag' = 1))
        }
        noits <- noits + 1
    }

    if(noits == 20){     
        print('Exiting Newton-Raphson because not reached optimum in 20 iterations!')
        errFlag <- 2 
    }else{ 
        errFlag <- 0 
    }
    
    if(printNrProgress){        
        print('Final values...')
        print(lnvPars)
        print(nllNew)
    }else{}
    
##############################################################
### since this is being used in nelder mead, round final nll to tolNr * 10 to avoid confusing the nm routine...    
##############################################################
    rndnll <- tolNr * 10
    nllNew <- round(nllNew / rndnll) * rndnll
    
    return(list('lnvPars' = lnvPars , 'nll' = nllNew , 'grad' = gradNew , 'hess' = hessNew  , 'fim' = fimNew, 
                'betahat' = betahatNew , 'vbetahat' = vbetahatNew , 'C' = CNew , 'noits' = noits , 'errFlag' = errFlag))
}

##########################################################
### calculate grad hess and fim for the updates...
##########################################################
gradHessIAK3D <- function(dA , lnvPars , W){
	n <- dim(W)[[1]]
	p <- dim(W)[[2]] - 1

###########################################        
### dT <- -T dC T	
### diC <- -iC dC iC
### dlndetC <- tr(iC dC)
### dlndetXiCX <- tr(iXiCX dXiCX)
### dXiCX <- -X iC dC iC X 
###########################################        
### dA was defined with all variance parameters = 1.    
    dC <- dA
    for (i in 1:length(dC)){ dC[[i]] <- exp(lnvPars[i]) * dC[[i]] }
    C <- dC[[1]] + dC[[2]] + dC[[3]] + dC[[4]] + dC[[5]]
    dCAll <- cbind(dC[[1]] , dC[[2]] , dC[[3]] , dC[[4]] , dC[[5]])   
    
    if(length(dC) == 6){ 
        C <- C + dC[[6]] 
        dCAll <- cbind(dCAll , dC[[6]])
    }else{}
   
	C <- 0.5 * (C + t(C))

#    return(list('dA' = dA , 'dC' = dC , 'lnvPars' = lnvPars , 'C' = C))

#	tmp <- lndetANDinvCb(C , W)
#    if(is.character(tmp$cholC)){

	iC <- try(chol2inv(chol(C)) , silent = TRUE)
    if(is.character(iC)){
        return(list('nll' = NA , 'grad' = NA , 'hess' = NA , 'fim' = NA , 'C' = C))
    }else{}
    
#	lndetC <- tmp$lndetC
#	iCW <- tmp$invCb
#
#	iCX <- as.matrix(iCW[, 1:p])
#	iCz <- as.matrix(iCW[, p+1])
#    iC <- chol2inv(tmp$cholC)

	iCW <- matrix(iC %*% W , ncol = p+1)
	iCX <- as.matrix(iCW[, 1:p])
	iCz <- as.matrix(iCW[, p+1])
	lndetC <- as.numeric(determinant(C , logarithm = TRUE)$modulus)

	WiCW <- t(W) %*% iCW
	WiCW  <- 0.5 * (WiCW + t(WiCW))

	XiCX <- WiCW[1:p , 1:p]
	XiCz <- WiCW[1:p , p+1]
	ziCz <- as.numeric(WiCW[p+1 , p+1])

#	tmp <- lndetANDinvCb(XiCX , XiCz)
#	lndetXiCX <- tmp$lndetC
#	betahat <- tmp$invCb
#    iXiCX <- chol2inv(tmp$cholC)

	iXiCX <- try(chol2inv(chol(XiCX)) , silent = TRUE)
    if(is.character(iXiCX)){
        return(list('nll' = NA , 'grad' = NA , 'hess' = NA , 'fim' = NA , 'C' = C))
    }else{}
	betahat <- matrix(iXiCX %*% XiCz , ncol = 1)
	lndetXiCX <- as.numeric(determinant(XiCX , logarithm = TRUE)$modulus)

    T <- iC - iCX %*% iXiCX %*% t(iCX)
	Tz <- as.matrix(iCz - iCX %*% betahat)
	zTz <- as.numeric(ziCz - 2 * t(betahat) %*% XiCz + t(betahat) %*% XiCX %*% betahat)

    TdCAll <- T %*% dCAll
    dCAllTz <- t(dCAll) %*% Tz
    zTdCTz <- t(Tz) %*% matrix(dCAllTz , nrow = n)

    grad <- trTdC <- NA * numeric(length(dC))
    hess <- fim <- matrix(NA , length(dC) , length(dC))
	for (i in 1:length(dC)){
        iThis <- ((i - 1) * n + 1):(i * n)

        trTdC[i] <- sum(diag(TdCAll[,iThis]))
 		grad[i] <- 0.5 * (trTdC[i] - zTdCTz[i])
        for (j in i:length(dC)){
            jThis <- ((j - 1) * n + 1):(j * n)
            TdCiTdCjThis <- TdCAll[,iThis] %*% TdCAll[,jThis]
            if (i == j){
                hess[i,j] <- as.numeric(0.5 * (trTdC[i] - sum(diag(TdCiTdCjThis)) - zTdCTz[i] + 2 * t(W[,p+1]) %*% TdCiTdCjThis %*% Tz))
                fim[i,j] <- as.numeric(0.5 * sum(diag(TdCiTdCjThis)))
            }else{
                hess[i,j] <- as.numeric(0.5 * (-sum(diag(TdCiTdCjThis)) + 2 * t(W[,p+1]) %*% TdCiTdCjThis %*% Tz))
                hess[j,i] <- hess[i,j]
                fim[i,j] <- as.numeric(0.5 * sum(diag(TdCiTdCjThis)))
                fim[j,i] <- fim[i,j]
            }
        }
	}
    
	nll <- 0.5 * ((n - p) * log(2 * pi) + lndetC + lndetXiCX + zTz) 
    
    return(list('nll' = nll , 'grad' = grad , 'hess' = hess , 'fim' = fim , 'betahat' = betahat , 'vbetahat' = iXiCX , 'C' = C))
}
