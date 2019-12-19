nrUpdatesIAK3DlnN <- function(zData , XData , vXU , iU , C , sigma2Vec){
 
    printNrProgress <- F

    pU <- length(iU)
    p <- dim(XData)[[2]]
    pK <- p - pU
    n <- length(zData)
    iK <- setdiff(seq(p) , iU)

    if(pU > 0){
        W <- cbind(XData , zData)          
    }else{
        W <- cbind(XData , zData - 0.5 * (sigma2Vec - diag(C)))
    }

#  	tmp <- lndetANDinvCb(C , W)
#	lndetC <- tmp$lndetC
#	iCW <- tmp$invCb
#    iC <- chol2inv(tmp$cholC)

	iC <- try(chol2inv(chol(C)) , silent = TRUE)
    if(is.character(iC)){
        return(list('betahatNew' = NA , 'nll' = NA , 'grad' = NA , 'hess' = NA  , 'fim' = NA, 
                'noits' = 0 , 'errFlag' = -123))
    }else{}
	iCW <- matrix(iC %*% W , ncol = p+1)
	lndetC <- as.numeric(determinant(C , logarithm = TRUE)$modulus)

    WiCW <- t(W) %*% iCW
	WiCW  <- 0.5 * (WiCW + t(WiCW))

    XiCX <- WiCW[1:p , 1:p,drop = FALSE]
    XiCz <- WiCW[1:p , p+1,drop = FALSE]
	ziCz <- as.numeric(WiCW[p+1 , p+1,drop = FALSE])

    if(pU == 0){
#     	tmp <- lndetANDinvCb(XiCX , XiCz) # the zData here already is zData - 0.5 * (sigma2 - diagC), and note that betavXUbeta = 0
#	    betahatNew <- tmp$invCb
     	betahatNew <- solve(XiCX , XiCz) # the zData here already is zData - 0.5 * (sigma2 - diagC), and note that betavXUbeta = 0
   		fimNew <- -XiCX
        hessNew <- gradNew <- NA # not needed, but hess = fim and grad = t(XData) iC (zData - XData beta - 0.5 *(sigma2 - diagC  + betavXbeta)), with betavXbeta = 0
        resiCres <- ziCz - 2 * t(betahatNew) %*% XiCz + t(betahatNew) %*% XiCX %*% betahatNew 
        
        nllNew <- 0.5 * (n * log(2 * pi) + lndetC + resiCres) 
        
        return(list('betahatNew' = betahatNew , 'nll' = nllNew , 'grad' = gradNew , 'hess' = hessNew  , 'fim' = fimNew, 
                'noits' = 0 , 'errFlag' = 0))
    }else{}


###################################
### define initial parameters...
###################################
#  	tmp <- lndetANDinvCb(XiCX , XiCz)
#	betaUNew <- tmp$invCb[iU]

  	betaUNew <- matrix(solve(XiCX , XiCz)[iU] , ncol = 1)

    if(pK > 0){
        XKiCXK <- WiCW[iK,iK,drop = FALSE]
        iXKiCXKXKiC <- solve(XKiCXK , t(iCW[,iK,drop = FALSE]))    
        TK <- iC - iCW[,iK,drop = FALSE] %*% iXKiCXKXKiC
    }else{
        XKiCXK <- c()
        iXKiCXKXKiC <- c()    
        TK <- iC
    }
    
###################################
### calc grad, hess and nll at inits...
###################################
    tmp <- gradHessIAK3DlnN(zData = zData , XData = XData , vXU = vXU , iU = iU , betaU = betaUNew , 
                            C = C , lndetC = lndetC , TK = TK , iXKiCXKXKiC = iXKiCXKXKiC , sigma2Vec = sigma2Vec)

    nllNew <- tmp$nll
    gradNew <- tmp$grad
    hessNew <- tmp$hess
    hessNew <- 0.5 * (hessNew + t(hessNew))
    fimNew <- tmp$fim
    fimNew <- 0.5 * (fimNew + t(fimNew))
    betahatNew <- tmp$betahat
    
    if(printNrProgress){        
        print('Initial values...')
        print(betahatNew)
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
#            betaUNew <- betaUNew - checkHess$invCb
            betaUNew <- betaUNew - checkHess
            
            tmp <- gradHessIAK3DlnN(zData = zData , XData = XData , vXU = vXU , iU = iU , betaU = betaUNew , 
                            C = C , lndetC = lndetC , TK = TK , iXKiCXKXKiC = iXKiCXKXKiC , sigma2Vec = sigma2Vec)

            nllNew <- tmp$nll
            gradNew <- tmp$grad
            hessNew <- tmp$hess
            hessNew <- 0.5 * (hessNew + t(hessNew))
            fimNew <- tmp$fim
            fimNew <- 0.5 * (fimNew + t(fimNew))
            betahatNew <- tmp$betahat
    
            if(printNrProgress){        
                print(paste0('After ' , noits + 1 , ' updates...'))
                print(betahatNew)
                print(nllNew)
            }else{}
        }else{
            print('Exiting Newton-Raphson because neither Hessian or FIM are positive definite!')
            return(list('betahat' = betahatNew , 'nll' = nllNew , 'grad' = gradNew , 'hess' = hessNew , 'fim' = fimNew , 
                        'noits' = noits , 'errFlag' = 1))
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
        print(betahatNew)
        print(nllNew)
    }else{}
    
##############################################################
### since this is being used in nelder mead, round final nll to tolNr * 10 to avoid confusing the nm routine...    
##############################################################
    rndnll <- tolNr * 10
    nllNew <- round(nllNew / rndnll) * rndnll
    
    return(list('betahatNew' = betahatNew , 'nll' = nllNew , 'grad' = gradNew , 'hess' = hessNew  , 'fim' = fimNew, 'iC' = iC ,
                'noits' = noits , 'errFlag' = errFlag))
}


gradHessIAK3DlnN <- function(zData , XData , vXU , iU , betaU , C , lndetC , TK , iXKiCXKXKiC , sigma2Vec){
###############################################################
### U represents the columns of XData that vary in the sample supports
### betaU the current values of these parameters
### vXU has the covariances of these covariates in their supports
### K represents the colums of XData that are always constant in the sample supports
### TK = iC - iC XK inv(XK iC XK) XK iC is the usual T matrix, 
### but defined with XK rather than the full XData
### iXKiCXKXKiC is inv(XK iC XK) XK iC
###############################################################
    pU <- length(iU)
    p <- dim(XData)[[2]]
    pK <- p - pU
    n <- length(zData)
    iK <- setdiff(seq(p) , iU)
    
    betavXbeta <- calcbetavXbeta(vXU , betaU)

    if(length(iU) > 1){
        XbetaU <- XData[,iU] %*% betaU 
    }else{
        XbetaU <- XData[,iU] * betaU 
    } 
    resU <- zData - XbetaU - 0.5 * (sigma2Vec - diag(C) + betavXbeta)

    betahat <- matrix(NA , p , 1)
    if(pU > 0){ betahat[iU] <- as.numeric(betaU) }else{}

    if (pK > 0){
        betaKhat <- iXKiCXKXKiC %*% resU
        betahat[iK] <- as.numeric(betaKhat)
    }else{
        betaKhat <- c()
    }

    vXUbetaU <- rowSums(vXU * kronecker(matrix(1 , n * pU , 1) , matrix(betaU , nrow = 1)))
    dbetavXbeta = matrix(vXUbetaU , n , pU , byrow = T)

    XUPLUSdbetavXbeta <- XData[,iU,drop = FALSE] + dbetavXbeta

    tmp <- TK %*% cbind(resU , XUPLUSdbetavXbeta)
    TKresU <- tmp[,1,drop = FALSE]
    TKXUPLUSdbetavXbeta <- tmp[,2:dim(tmp)[[2]],drop = FALSE]

    nll <- as.numeric(0.5 * (n * log(2 * pi) + lndetC + t(resU) %*% TKresU))

    grad <- -t(XUPLUSdbetavXbeta) %*% TKresU

    tmp <- matrix(as.numeric(t(vXU)), pU * pU , n)
    vXTKresU <- matrix(tmp %*% TKresU , pU , pU)
    XUPLUSdbetavXbetaTKXUPLUSdbetavXbeta <- t(XUPLUSdbetavXbeta) %*% TKXUPLUSdbetavXbeta
    hess <- -vXTKresU + XUPLUSdbetavXbetaTKXUPLUSdbetavXbeta

    fim <- XUPLUSdbetavXbetaTKXUPLUSdbetavXbeta

    return(list('nll' = nll , 'grad' = grad , 'hess' = hess , 'fim' = fim , 'betahat' = betahat))
}
