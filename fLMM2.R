################################################################
### a wrapper for the two different versions of fitLMM2
### (NoPrior or Bayes)
################################################################
fitLMM2 <- function(c , z , X , nAliquots = 1 , rdmEffFactors = NULL , beta0 = NULL , pbeta0 = NULL , 
                    covModel = 'exponential' , nSpatStructs = 1 , incNugget = T , 
                    optionML = F , verbose = FALSE , mina = NULL , maxa = NULL , parsInit = NULL , 
                    attachBigMats = TRUE , fixa2 = FALSE , parsFit = NULL){

  if(((is.null(beta0)) & (!is.null(pbeta0))) | ((!is.null(beta0)) & (is.null(pbeta0)))){
    stop('Error - you have defined only one of beta0 and pbeta0 - define neither or both!')
  }else{}
  if(!is.null(rdmEffFactors)){
    if(is.matrix(rdmEffFactors)){ rdmEffFactors <- as.data.frame(rdmEffFactors) }else{}
    
    if(is.data.frame(rdmEffFactors)){
      if(nrow(rdmEffFactors) != length(z)){ stop('Error - rdmEffFactors must be a matrix with the same number of rows as you have data!') }else{} 
    }else{
      if(length(rdmEffFactors) == length(z)){ 
        rdmEffFactors <- data.frame('rdmEff1' = rdmEffFactors) 
      }else{
        stop('Error - rdmEffFactors must be a matrix with the same number of rows as you have data!')        
      }
    }
    
    ### make sure all are factors...
    for(i in 1:ncol(rdmEffFactors)){
      if(!is.factor(rdmEffFactors[,i])){ rdmEffFactors[,i] <- factor(rdmEffFactors[,i]) }else{}
    }
  }else{}
  
  if(is.null(beta0) & is.null(pbeta0) & is.null(rdmEffFactors)){
    listRtn <- fitLMM2_NoPrior(c = c , z = z , X = X , nAliquots = nAliquots , 
                               covModel = covModel , nSpatStructs = nSpatStructs , incNugget = incNugget , 
                               optionML = optionML , verbose = verbose , mina = mina , maxa = maxa , parsInit = parsInit , 
                               attachBigMats = attachBigMats , fixa2 = fixa2 , parsFit = parsFit)
  }else if(((!is.null(beta0)) & (!is.null(pbeta0))) | (!is.null(rdmEffFactors))){
    listRtn <- fitLMM2_Bayes(c = c , z = z , X = X , nAliquots = nAliquots , rdmEffFactors = rdmEffFactors , beta0 = beta0 , pbeta0 = pbeta0 , 
                             covModel = covModel , nSpatStructs = nSpatStructs , incNugget = incNugget , 
                             optionML = optionML , verbose = verbose , mina = mina , maxa = maxa , parsInit = parsInit , 
                             attachBigMats = attachBigMats , fixa2 = fixa2 , parsFit = parsFit)
  }else{
    stop('Error - you have defined only one of beta0 and pbeta0 - define neither or both!')
  }
  
  return(listRtn)
}

################################################################
### previously called fitLMM2, now fitLMM2_NoPrior, 
### fitLMM2 function written over with a more general wrapper, above...
################################################################
fitLMM2_NoPrior <- function(c , z , X , nAliquots = 1 , 
                    covModel = 'exponential' , nSpatStructs = 1 , incNugget = T , 
                    optionML = F , verbose = FALSE , mina = NULL , maxa = NULL , parsInit = NULL , 
                    attachBigMats = TRUE , fixa2 = FALSE , parsFit = NULL){

  ### nAliquots used as if these cores were sampled a negligible distance from each other
  ### so diagonal of cov mtx for entries where nAliquots > 1 = c1 + c2 + c0/nAliquots
  if(length(nAliquots) > 1){ if(length(nAliquots) != length(z)){ stop('Error - nAliquots must be same length as z for fitLMM2!') }else{} }else{}
  if(fixa2 & nSpatStructs == 1){ stop('Error - asked to fixa2, but only put nSpatStructs = 1!') }else{}
  
  ina <- which(is.na(z) | (rowSums(is.na(X)) > 0))
  if(length(ina) > 0){
    print('Attention - some NA values in the data, getting removed for fitting LMM...')
    c <- c[-ina,,drop=FALSE]
    z <- z[-ina]
    X <- X[-ina,,drop=FALSE]
    if(length(nAliquots) > 1){ nAliquots <- nAliquots[-ina] }else{}
  }else{}
  
  if(nSpatStructs == 0){
    if(!incNugget){ stop('Error - put incNugget = T for nSpatStructs = 0! (covModel will be ignored)') }else{}

    XX <- t(X) %*% X
    Xz <- t(X) %*% z
    zz <- sum(z * z)
    
    n <- length(z)
    p <- ncol(X)

    if((length(nAliquots) > 1) || (nAliquots > 1)){
      A <- sparseMatrix(i = seq(n) , j = seq(n) , x = 1 / nAliquots , symmetric = T)

      badNll <- 9E99
      pars <- c()
      parsOut <- ''
      listOut <- list('nll' = badNll , 'sigma2hat' = NA , 'betahat' = NA , 'vbetahat' = NA , 'pars' = pars , 'covParams4Kriging' = NA , 'C' = NA , 'iC' = NA)

      listRtn <- fLMM2_Given_A(A = A , z = z , X = X , nAliquots = nAliquots , 
                               returnAll = T , optionML = optionML , verbose = verbose , forCompLik = FALSE , 
                               badNll = badNll , attachBigMats = attachBigMats , 
                               pars = pars , parsOut = parsOut , listOut = listOut)

      return(listRtn)
    }else{
      if(p > 1){
        cholXX <- try(chol(XX) , silent = TRUE)
        if(is.character(cholXX) | min(eigen(XX)$value) <= 0){
          return(list('nll' = NA  , 'sigma2hat' = NA, 'betahat' = matrix(NA , p , 1) , 'vbetahat' = matrix(NA , p , p) , 'pars' = c() , 
                      'covModel' = covModel , 'covParams4Kriging' = rep(NA , 5) , 'nSpatStructs' = nSpatStructs , 
                      'XiAX' = XX , 'XiAz' = Xz , 'ziAz' = zz , 'lndetA' = 0 , 
                      'iA' = sparseMatrix(i = seq(n) , j = seq(n) , x = 1 , symmetric = T) , 
                      'iAX' = X , 'iAz' = z))
        }else{}
        vbetahat <- chol2inv(cholXX)
        betahat <- matrix(vbetahat %*% Xz , ncol = 1)
        
        lndetXX <- determinant(XX , logarithm = TRUE)
        if(lndetXX$sign < 0){
          return(list('nll' = NA  , 'sigma2hat' = NA, 'betahat' = matrix(NA , p , 1) , 'vbetahat' = matrix(NA , p , p) , 'pars' = c() , 
                      'covModel' = covModel , 'covParams4Kriging' = rep(NA , 5) , 'nSpatStructs' = nSpatStructs , 
                      'XiAX' = XX , 'XiAz' = Xz , 'ziAz' = zz , 'lndetA' = 0 , 
                      'iA' = sparseMatrix(i = seq(n) , j = seq(n) , x = 1 , symmetric = T) , 
                      'iAX' = X , 'iAz' = z))
        }else{}
        lndetXX <- as.numeric(lndetXX$modulus)
      }else{
        vbetahat <- matrix(1 / XX , 1 , 1)
        betahat <- matrix(vbetahat %*% Xz , ncol = 1)
        lndetXX <- as.numeric(log(XX))
      }
      
      res <- z - X %*% betahat
      if(!optionML){
        sigma2hat <- as.numeric(t(res) %*% res) / (n - p)
        if((!is.na(sigma2hat)) && (sigma2hat > 0)){
          nll <- 0.5 * ((n - p) * log(2 * pi) + (n - p) * log(sigma2hat) + lndetXX + n - p)
        }else{
          nll <- NA
          sigma2hat <- NA
        }
      }else{
        sigma2hat <- as.numeric(t(res) %*% res) / n
        if((!is.na(sigma2hat)) && (sigma2hat > 0)){
          nll <- 0.5 * n * (log(2 * pi) + log(sigma2hat) + 1)
        }else{
          nll <- NA
          sigma2hat <- NA
        }
      }
      
      vbetahat <- vbetahat * sigma2hat
      
      return(list('nll' = nll , 'sigma2hat' = sigma2hat , 'betahat' = betahat , 'vbetahat' = vbetahat , 'pars' = c() , 
                  'covModel' = covModel , 'covParams4Kriging' = c(sigma2hat , 0 , 1 , 0 , 2) , 'nSpatStructs' = nSpatStructs , 
                  'XiAX' = XX , 'XiAz' = Xz , 'ziAz' = zz , 'lndetA' = 0 , 
                  'iA' = sparseMatrix(i = seq(n) , j = seq(n) , x = 1 , symmetric = T) , 
                  'iAX' = X , 'iAz' = z))
    }
    
  }else{}
  
  D <- xyDist(c , c)

### note parameterisation in fLMM2 fn for exp model is with a = range (not 3a = range)
  if(is.null(mina)){ mina <- min(D[D > 0]) }else{}
  if(is.null(maxa)){ maxa <- max(D) / 2 }else{}

  if(covModel == 'spherical' | substr(covModel , 1 , 8) == 'wendland' | substr(covModel , 1 , 9) == 'gwendland'){
    ### use sparse matrices (if v big, will have to rethink this bit)
    D[D > maxa] <- 0 # not really D=0, but cov fn will copy the 0s to use as covariances
    D <- as(D , "dsCMatrix")
  }else{}    

  ### check D for any really small dists...
  if(length(which((D > 0) & (D <= 1E-10))) > 0){ print('WARNING - SOME SMALL VALUES IN D COULD BE COLOCATED - CHECK DEFINITION OF D!') }else{}

  if(is.null(parsFit)){
    ############################################
    ### model to be fitted now..
    ############################################

    ############################################
    ### set initial pars...
    ############################################
    if(is.null(parsInit)){
      parsInit <- getParsInit_fLMM2(covModel = covModel , nSpatStructs = nSpatStructs , incNugget = incNugget , 
                                    mina = mina , maxa = maxa , fixa2 = fixa2)
    }else{}
    
    #####################################################    
    ### run to get initial lik val...
    #####################################################    
    ftLMM2Init <- fLMM2(pars = parsInit , c = c , z = z , X = X , nAliquots = nAliquots , D = D , 
                        covModel = covModel , nSpatStructs = nSpatStructs , incNugget = incNugget , returnAll = FALSE , optionML = optionML , 
                        verbose = verbose , forCompLik = FALSE , mina = mina , maxa = maxa)
    
    #####################################################    
    ### fit model with temporal correlation accounted for...
    #####################################################    
    badNll <- 9E99
    if(length(parsInit) > 1){
      
      if(fixa2){
        # at most 3 params to optimize with a2 fixed. 
        # maybe no need to iterate with optimIt...
        vecFixed <- rep(F , length(parsInit))
        vecFixed[length(parsInit)] <- T
        tmp <- optifix(par = parsInit , fn = fLMM2 , fixed = vecFixed , method = 'Nelder-Mead' , 
                       c = c , z = z , X = X , nAliquots = nAliquots , D = D ,
                       covModel = covModel , nSpatStructs = nSpatStructs , incNugget = incNugget , returnAll = FALSE , optionML = optionML ,
                       verbose = verbose , forCompLik = FALSE , mina = mina , maxa = maxa , badNll = badNll)
        
        parsFit <- tmp$fullpars
        
      }else{
        tmp <- optimIt(par = parsInit , fn = fLMM2 , methodOptim = c('Nelder-Mead') , 
                       c = c , z = z , X = X , nAliquots = nAliquots , D = D ,
                       covModel = covModel , nSpatStructs = nSpatStructs , incNugget = incNugget , returnAll = FALSE , optionML = optionML ,
                       verbose = verbose , forCompLik = FALSE , mina = mina , maxa = maxa , badNll = badNll)
        parsFit <- tmp$par
      }
    }else{    
      tmp <- optimize(f = fLMM2 , lower = mina , upper = maxa , c = c , z = z , X = X , nAliquots = nAliquots , D = D ,
                      covModel = covModel , nSpatStructs = nSpatStructs , incNugget = incNugget , returnAll = FALSE , optionML = optionML , 
                      verbose = verbose , forCompLik = FALSE , mina = mina , maxa = maxa , badNll = badNll)
      parsFit <- tmp$minimum
    }
  }else{} # parsFit given as input to fn

#####################################################    
### run again to get final lik val + other info...
#####################################################    
  ftLMM2Fit <- fLMM2(pars = parsFit , c = c , z = z , X = X , nAliquots = nAliquots , D = D , 
                      covModel = covModel , nSpatStructs = nSpatStructs , incNugget = incNugget , returnAll = TRUE , optionML = optionML , 
                      verbose = verbose , forCompLik = FALSE , mina = mina , maxa = maxa , attachBigMats = attachBigMats)
  
  return(ftLMM2Fit)
}

###########################################################################
### the fit function specifically for Bayes
### or when rdmEffFactors is not null (rdm Effects are added)
###########################################################################
fitLMM2_Bayes <- function(c , z , X , nAliquots = 1 , rdmEffFactors = NULL , beta0 = NULL , pbeta0 = NULL , 
                          covModel = 'exponential' , nSpatStructs = 1 , incNugget = T , 
                          optionML = F , verbose = FALSE , mina = NULL , maxa = NULL , parsInit = NULL , 
                          attachBigMats = TRUE , fixa2 = FALSE , parsFit = NULL){
  
  ### nAliquots used as if these cores were sampled a negligible distance from each other
  ### so diagonal of cov mtx for entries where nAliquots > 1 = c1 + c2 + c0/nAliquots
  if(length(nAliquots) > 1){ if(length(nAliquots) != length(z)){ stop('Error - nAliquots must be same length as z for fitLMM2!') }else{} }else{}
  if(fixa2 & nSpatStructs == 1){ stop('Error - asked to fixa2, but only put nSpatStructs = 1!') }else{}
  
  ina <- which(is.na(z) | (rowSums(is.na(X)) > 0))
  if(length(ina) > 0){
    print('Attention - some NA values in the data, getting removed for fitting LMM...')
    c <- c[-ina,,drop=FALSE]
    z <- z[-ina]
    X <- X[-ina,,drop=FALSE]
    if(length(nAliquots) > 1){ nAliquots <- nAliquots[-ina] }else{}
  }else{}
  
  n <- length(z)
  p <- ncol(X)
  
  ### put prior info matrices together in list...
  if(!is.null(beta0)){ 
    beta0 <- matrix(beta0 , ncol = 1)
    if(is.null(pbeta0)){ stop('Error - you have entered a prior beta as beta0, but not an associated prior precision, pbeta0!') }else{}
    
    priorInfo <- list('beta0' = beta0 , 'pbeta0' = pbeta0 , 
                      'pbeta0beta0' = pbeta0 %*% beta0)
    priorInfo$beta0pbeta0beta0 <- t(priorInfo$beta0) %*% priorInfo$pbeta0beta0
  }else{
    if(!is.null(pbeta0)){ stop('Error - you have entered a prior precision pbeta as pbeta0, but not an associated prior estimate, beta0!') }else{}
    priorInfo <- NULL
  }
  
  if(nSpatStructs == 0){
    stop('Error - fitLMM2 in fLMM2_Bayes.R only set up for including a spatial structure.')
  }else{}
  
  D <- xyDist(c , c)
  
  ### note parameterisation in fLMM2 fn for exp model is with a = range (not 3a = range)
  if(is.null(mina)){ mina <- min(D[D > 0]) }else{}
  if(is.null(maxa)){ maxa <- max(D) / 2 }else{}
  
  if(covModel == 'spherical' | substr(covModel , 1 , 8) == 'wendland' | substr(covModel , 1 , 9) == 'gwendland'){
    ### use sparse matrices (if v big, will have to rethink this bit)
    D[D > maxa] <- 0 # not really D=0, but cov fn will copy the 0s to use as covariances
    D <- as(D , "dsCMatrix")
  }else{}    
  
  ### check D for any really small dists...
  if(length(which((D > 0) & (D <= 1E-10))) > 0){ print('WARNING - SOME SMALL VALUES IN D COULD BE COLOCATED - CHECK DEFINITION OF D!') }else{}
  
  if(is.null(parsFit)){
    ############################################
    ### model to be fitted now..
    ############################################
    
    ############################################
    ### set initial pars...
    ############################################
    if(is.null(parsInit)){
      parsInit <- getParsInit_fLMM2_Bayes(covModel = covModel , nSpatStructs = nSpatStructs , incNugget = incNugget , 
                                          mina = mina , maxa = maxa , fixa2 = fixa2 , z = z , X = X , rdmEffFactors = rdmEffFactors)
    }else{}
    
    #####################################################    
    ### run to get initial lik val...
    #####################################################    
    ftLMM2Init <- fLMM2_Bayes(pars = parsInit , c = c , z = z , X = X , nAliquots = nAliquots ,rdmEffFactors = rdmEffFactors , D = D , priorInfo = priorInfo ,
                              covModel = covModel , nSpatStructs = nSpatStructs , incNugget = incNugget , returnAll = FALSE , optionML = optionML , 
                              verbose = verbose , forCompLik = FALSE , mina = mina , maxa = maxa)
    
    #####################################################    
    ### fit model with temporal correlation accounted for...
    #####################################################    
    badNll <- 9E99
    if(length(parsInit) > 1){
      
      if(fixa2){
        # at most 3 params to optimize with a2 fixed. 
        # maybe no need to iterate with optimIt...
        vecFixed <- rep(F , length(parsInit))
        vecFixed[length(parsInit)] <- T
        tmp <- optifix(par = parsInit , fn = fLMM2_Bayes , fixed = vecFixed , method = 'Nelder-Mead' , 
                       c = c , z = z , X = X , nAliquots = nAliquots , rdmEffFactors = rdmEffFactors , D = D , priorInfo = priorInfo ,
                       covModel = covModel , nSpatStructs = nSpatStructs , incNugget = incNugget , returnAll = FALSE , optionML = optionML ,
                       verbose = verbose , forCompLik = FALSE , mina = mina , maxa = maxa , badNll = badNll)
        
        parsFit <- tmp$fullpars
        
      }else{
        tmp <- optimIt(par = parsInit , fn = fLMM2_Bayes , methodOptim = c('Nelder-Mead') , 
                       c = c , z = z , X = X , nAliquots = nAliquots , rdmEffFactors = rdmEffFactors , D = D , priorInfo = priorInfo ,
                       covModel = covModel , nSpatStructs = nSpatStructs , incNugget = incNugget , returnAll = FALSE , optionML = optionML ,
                       verbose = verbose , forCompLik = FALSE , mina = mina , maxa = maxa , badNll = badNll)
        parsFit <- tmp$par
      }
    }else{    
      tmp <- optimize(f = fLMM2_Bayes , lower = mina , upper = maxa , c = c , z = z , X = X , nAliquots = nAliquots , rdmEffFactors = rdmEffFactors , D = D , priorInfo = priorInfo ,
                      covModel = covModel , nSpatStructs = nSpatStructs , incNugget = incNugget , returnAll = FALSE , optionML = optionML , 
                      verbose = verbose , forCompLik = FALSE , mina = mina , maxa = maxa , badNll = badNll)
      parsFit <- tmp$minimum
    }
  }else{} # parsFit given as input to fn
  
  #####################################################    
  ### run again to get final lik val + other info...
  #####################################################    
  ftLMM2Fit <- fLMM2_Bayes(pars = parsFit , c = c , z = z , X = X , nAliquots = nAliquots , rdmEffFactors = rdmEffFactors , D = D , priorInfo = priorInfo ,
                           covModel = covModel , nSpatStructs = nSpatStructs , incNugget = incNugget , returnAll = TRUE , optionML = optionML , 
                           verbose = verbose , forCompLik = FALSE , mina = mina , maxa = maxa , attachBigMats = attachBigMats)
  
  return(ftLMM2Fit)
}

##########################################################
### function to get initial params for optim for fLMM2...
##########################################################
getParsInit_fLMM2 <- function(covModel = NULL , nSpatStructs = NULL , incNugget = NULL , mina = NULL , maxa = NULL , fixa2 = NULL){
  if(is.null(covModel) | is.null(nSpatStructs) | is.null(incNugget)){
    stop('Error - getParsInit_fLMM2 only makes sense if all input arguments defined!')
  }else{}

  if(nSpatStructs == 0){
    parsInit <- c()
  }else{
    if(is.null(mina) | is.null(maxa) | is.null(fixa2)){
      stop('Error - getParsInit_fLMM2 only makes sense if all input arguments defined!')
    }else{}
    
    if(incNugget){
      if(nSpatStructs == 1){
        parsInit <- c(0.6 , 0.5 * mina + 0.5 * maxa)
      }else if(nSpatStructs == 2){
        parsInit <- c(0.3 , 0.9 * mina + 0.1 * maxa , 0.3 , 0.3 * mina + 0.7 * maxa)
        if(fixa2){ parsInit[4] <- maxa }else{}
      }else{
        stop('Error - generalise fitLMM2 for other values of nSpatStructs!')
      }
    }else{
      if(nSpatStructs == 1){
        parsInit <- c(0.5 * mina + 0.5 * maxa)
      }else if(nSpatStructs == 2){
        parsInit <- c(0.5 , 0.9 * mina + 0.1 * maxa , 0.3 * mina + 0.7 * maxa)
        if(fixa2){ parsInit[3] <- maxa }else{}
      }else{
        stop('Error - generalise fitLMM2 for other values of nSpatStructs!')
      }
    }
  }  

  return(parsInit)
}

##########################################################
### function to calc the nll for fLMM2...
##########################################################
fLMM2 <- function(pars , c , z , X , nAliquots = 1 , D , covModel , nSpatStructs , incNugget = T , returnAll , optionML , verbose , forCompLik = FALSE , mina = NULL , maxa = NULL , badNll = 9E99 , attachBigMats = TRUE , timeCodeNow = FALSE){
    N <- length(z)
    np <- dim(X)[2]

    if(length(nAliquots) > 1){ if(length(nAliquots) != length(z)){ stop('Error - nAliquots must be same length as z for fitLMM2!') }else{} }else{}

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
      if(timeCodeNow){ start_time <- Sys.time() }else{}
      # C <- defineCLMM2(c0 = s0 , c1 = s1 , a1 = a1 , c2 = s2 , a2 = a2 , D = D , nAliquots = nAliquots , covModel = covModel , nSpatStructs = nSpatStructs)
      A <- defineCLMM2(c0 = s0 , c1 = s1 , a1 = a1 , c2 = s2 , a2 = a2 , D = D , nAliquots = nAliquots , covModel = covModel , nSpatStructs = nSpatStructs)
      if(timeCodeNow){ print('Time to define cor mtx A:') ; print(Sys.time() - start_time) }else{}
    }else{
      return(listOut)
    }   

    if(timeCodeNow){ start_time <- Sys.time() }else{}
    listRtn <- fLMM2_Given_A(A = A , z = z , X = X , nAliquots = nAliquots , 
                             returnAll = returnAll , optionML = optionML , verbose = verbose , forCompLik = FALSE , 
                             badNll = badNll , attachBigMats = attachBigMats , 
                             pars = pars , parsOut = parsOut , listOut = listOut)
    if(timeCodeNow){ print('Time to calculate nll:') ; print(Sys.time() - start_time) }else{}
    
    if(returnAll){
      ### update the covMode,covParams4Kriging,nSpatStructs...
      listRtn$covModel <- covModel
      listRtn$covParams4Kriging <- c((1 - s1 - s2) * listRtn$sigma2hat , s1 * listRtn$sigma2hat , a1 , s2 * listRtn$sigma2hat , a2)
      listRtn$nSpatStructs <- nSpatStructs
    }else{}
    
    return(listRtn)
}

###################################################################################
### with the correlation matrix A defined, return the fLMM2 output...
###################################################################################
fLMM2_Given_A <- function(A , z , X , nAliquots = 1 , returnAll , optionML , verbose , forCompLik = FALSE , badNll = 9E99 , attachBigMats = TRUE , 
                          pars , parsOut , listOut){
  N <- length(z)
  np <- dim(X)[2]
  
  if(length(nAliquots) > 1){ if(length(nAliquots) != length(z)){ stop('Error - nAliquots must be same length as z for fitLMM2!') }else{} }else{}
  if(missing(returnAll)){ returnAll <- F }else{} # default is to return just nll
  if(missing(optionML)){ optionML <- F }else{} # default is REML
  if(missing(verbose)){ verbose <- T }else{} # default is to print the current nll
  
  oneszX <- cbind(matrix(1 , N , 1) , z , X)
  
  tmp <- lndetANDinvCb(A , oneszX)
  
  if(returnAll){ cholA <- chol(A) }else{}

  lndetA <- tmp$lndetC
  if(is.na(lndetA)){             
    printNll(nll = NA , parsOut = parsOut , verbose = verbose)
    return(listOut) 
  }else{}
  
  ### NOTE CONFUSING NOTATION - iConeszX etc are actually iAoneszX at this stage!
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
      return(list('zinvAz' = zinvCz , 'zinvAX' = zinvCX , 'XinvAX' = XinvCX , 'lndetA' = lndetA , 'parsOut' = parsOut , 
                  'invAX' = invCX , 'invAz' = invCz , 'iA' = chol2inv(cholA) , 
                  'onesinvAones' = onesinvCones , 'onesinvAz' = onesinvCz , 'onesinvAX' = onesinvCX))        
    }else{
      return(list('zinvAz' = zinvCz , 'zinvAX' = zinvCX , 'XinvAX' = XinvCX , 'lndetA' = lndetA , 'parsOut' = parsOut , 
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
    if((!is.na(sigma2hat)) && (sigma2hat > 0)){
      lndetXinvCX <- lndetXinvCX - np * log(sigma2hat)
    }else{
      lndetXinvCX <- NA
    }
  }        
  
  if(returnAll){ 
    XiAX <- XinvCX
    XiAz <- t(zinvCX)
    ziAz <- zinvCz
    vbetahat <- sigma2hat * chol2inv(chol(XinvCX)) 
  }else{}
  if((!is.na(sigma2hat)) && (sigma2hat > 0)){
    lndetC <- lndetA + N * log(sigma2hat)
    resiCres <- resiCres / sigma2hat
  }else{
    lndetC <- NA
    resiCres <- NA
  }
  nll <- -0.5 * (N-np) * log(2 * pi) - 0.5 * lndetC - 0.5 * lndetXinvCX - 0.5 * resiCres
  nll <- -nll

  printNll(nll = nll , parsOut = parsOut , verbose = verbose)
  
  if(is.na(nll)){ nll <- badNll }else{}
  
  if(returnAll){
    # covParams4Kriging <- c((1 - s1) * sigma2hat , s1 * sigma2hat , a1 , s2 * sigma2hat , a2)
    # covParams4Kriging <- c((1 - s1 - s2) * sigma2hat , s1 * sigma2hat , a1 , s2 * sigma2hat , a2)
    # covModel,covParams4Kriging,nSpatStructs just being reserved a space in the list for now...
    if(attachBigMats){
      C <- sigma2hat * A    
      iA <- chol2inv(cholA)
      iC <- iA / sigma2hat
      return(list('nll' = nll , 'sigma2hat' = sigma2hat , 'betahat' = betahat , 'vbetahat' = vbetahat , 'pars' = pars , 
                  'covModel' = NA , 'covParams4Kriging' = NA , 'nSpatStructs' = NA , 
                  'C' = C , 'iC' = iC , 'XiAX' = XiAX , 'XiAz' = XiAz , 'ziAz' = ziAz , 'lndetA' = lndetA , 
                  'iA' = iA , 'iAX' = iConeszX[,-c(1,2),drop=FALSE] , 'iAz' = iConeszX[,2]))
    }else{
      return(list('nll' = nll , 'sigma2hat' = sigma2hat , 'betahat' = betahat , 'vbetahat' = vbetahat , 'pars' = pars , 
                  'covModel' = NA , 'covParams4Kriging' = NA , 'nSpatStructs' = NA , 
                  'XiAX' = XiAX , 'XiAz' = XiAz , 'ziAz' = ziAz))
    }
  }else{             
    return(nll)
  }
}

###################################################################################
### a further wrap up, with WiAW and lndetA defined, return the fLMM2 output...
### returnOpt = 0 for nll; 1 for nll & sigma2hat ; 2 for also betahat & vbetahat
###################################################################################
fLMM2_Given_WiAW <- function(XiAX , XiAz , ziAz , lndetA , n , optionML , badNll = 9E99 , returnOpt = 0){
  
  np <- ncol(XiAX)
  if(is.null(np)){ np <- 1 }else{}

  tmp <- lndetANDinvCb(XiAX , XiAz)
  lndetXinvCX <- tmp$lndetC
  if(is.na(lndetXinvCX)){     
    if(returnOpt == 0){
      return(badNll)
    }else if(returnOpt == 1){
      return(list('nll' = badNll , 'sigma2hat' = NA)) 
    }else if(returnOpt == 2){
      return(list('nll' = badNll , 'betahat' = NA , 'vbetahat' = NA , 'sigma2hat' = NA)) 
    }else{
      stop('Unknown returnOpt!')
    }
  }else{}
  betahat <- tmp$invCb

  resiCres <- as.numeric(ziAz - 2 * t(XiAz) %*% betahat + t(betahat) %*% XiAX %*% betahat)
  
  if (optionML){
    sigma2hat <- resiCres / n
    lndetXinvCX <- 0
    np <- 0 # not really, but just because it appears in const of nll
  }else{
    sigma2hat <- resiCres / (n - np)
    if((!is.na(sigma2hat)) && (sigma2hat > 0)){
      lndetXinvCX <- lndetXinvCX - np * log(sigma2hat)
    }else{
      lndetXinvCX <- NA
    }
  }        
  
  if(returnOpt > 1){ 
    vbetahat <- sigma2hat * chol2inv(chol(XiAX)) 
  }else{}

  if((!is.na(sigma2hat)) && (sigma2hat > 0)){
    lndetC <- lndetA + n * log(sigma2hat)
    resiCres <- resiCres / sigma2hat
  }else{
    lndetC <- NA
    resiCres <- NA
  }
  nll <- -0.5 * (n-np) * log(2 * pi) - 0.5 * lndetC - 0.5 * lndetXinvCX - 0.5 * resiCres
  nll <- -nll

  if(is.na(nll)){
    nll <- badNll
  }else{}
  
  if(returnOpt == 0){
    return(nll)
  }else if(returnOpt == 1){
    return(list('nll' = nll , 'sigma2hat' = sigma2hat)) 
  }else if(returnOpt == 2){
    return(list('nll' = nll , 'betahat' = betahat , 'vbetahat' = vbetahat , 'sigma2hat' = sigma2hat)) 
  }else{
    stop('Unknown returnOpt!')
  }
}

##########################################################
### function to get initial params for optim for fLMM2_Bayes 
### (ie prior info given, so we will fit all covariance parameters in optim, 
###  rather than fitting total variance analytically when no prior given)...
##########################################################
getParsInit_fLMM2_Bayes <- function(covModel = NULL , nSpatStructs = NULL , incNugget = NULL , mina = NULL , maxa = NULL , fixa2 = NULL , z , X , rdmEffFactors = NULL){
  if(is.null(covModel) | is.null(nSpatStructs) | is.null(incNugget)){
    stop('Error - getParsInit_fLMM2 only makes sense if all input arguments defined!')
  }else{}
  
  betahatInit <- solve(t(X) %*% X , t(X) %*% z)
  resInit <- z - X %*% betahatInit
  varInit <- sum(resInit ^ 2) / (length(z) - 1)
  
  if(nSpatStructs == 0){
    parsInit <- c(varInit)
  }else{
    if(is.null(mina) | is.null(maxa) | is.null(fixa2)){
      stop('Error - getParsInit_fLMM2 only makes sense if all input arguments defined!')
    }else{}
    
    if(incNugget){
      if(nSpatStructs == 1){
        parsInit <- c(varInit * 0.4 , varInit * 0.6 , 0.5 * mina + 0.5 * maxa)
      }else if(nSpatStructs == 2){
        parsInit <- c(varInit * 0.4 , varInit * 0.3 , 0.9 * mina + 0.1 * maxa , varInit * 0.3 , 0.3 * mina + 0.7 * maxa)
        if(fixa2){ parsInit[5] <- maxa }else{}
      }else{
        stop('Error - generalise fitLMM2 for other values of nSpatStructs!')
      }
    }else{
      if(nSpatStructs == 1){
        parsInit <- c(varInit , 0.5 * mina + 0.5 * maxa)
      }else if(nSpatStructs == 2){
        parsInit <- c(varInit * 0.5 , 0.9 * mina + 0.1 * maxa , varInit * 0.5 , 0.3 * mina + 0.7 * maxa)
        if(fixa2){ parsInit[4] <- maxa }else{}
      }else{
        stop('Error - generalise fitLMM2 for other values of nSpatStructs!')
      }
    }
  }  

  if(!is.null(rdmEffFactors)){
    parsInit <- c(parsInit , rep(varInit * 0.1 , ncol(rdmEffFactors)))
  }else{}

  return(parsInit)
}

##########################################################
### function to calc the nll for fLMM2 with priorInfo given...
### or with rdmEffFactors.
### This version optimises all cov parameters numerically
### rather than no prior version, which optimises one variance parameter analytically
##########################################################
fLMM2_Bayes <- function(pars , c , z , X , nAliquots = 1 , rdmEffFactors = NULL , D , priorInfo = NULL , covModel , nSpatStructs , incNugget = T , returnAll , optionML , verbose , forCompLik = FALSE , mina = NULL , maxa = NULL , badNll = 9E99 , attachBigMats = TRUE , timeCodeNow = FALSE){

  if(forCompLik){ stop('Error - fLMM2_Bayes not been written for comp lik!') }else{}
  
  N <- length(z)
  np <- dim(X)[2]
  
  if(length(nAliquots) > 1){ if(length(nAliquots) != length(z)){ stop('Error - nAliquots must be same length as z for fitLMM2!') }else{} }else{}
  
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
      c0 <- pars[1] ; c1 <- pars[2] ; a1 <- pars[3] ; c2 <- 0 ; a2 <- a1 ; iParNext <- 4
      if((c0 < 0) | (c1 < 0) | (a1 < mina) | (a1 > maxa)){ paramsOK <- F }else{ paramsOK = T } 
    }else{
      c0 <- pars[1] ; c1 <- pars[2] ; a1 <- pars[3]; c2 <- pars[4] ; a2 <- pars[5] ; iParNext <- 6
      if((c0 < 0) | (c1 < 0) | (c2 < 0) | (a2 < a1)  | (a1 < mina) | (a2 > maxa)){ paramsOK <- F }else{ paramsOK = T } 
    }
  }else{
    if (nSpatStructs == 1){
      c0 <- 0 ; c1 <- pars[1] ; a1 <- pars[2] ; c2 <- 0 ; a2 <- a1 ; iParNext <- 3
      if((c1 < 0) | (c2 < 0) | (a1 < mina) | (a1 > maxa)){ paramsOK <- F }else{ paramsOK = T } 
    }else{
      c0 <- 0 ; c1 <- pars[1] ; a1 <- pars[2]; c2 <- pars[3] ; a2 <- pars[4] ; iParNext <- 5
      if((c1 < 0) | (c2 < 0) | (a2 < a1)  | (a1 < mina) | (a2 > maxa)){ paramsOK <- F }else{ paramsOK = T } 
    }
  }

  if(!is.null(rdmEffFactors)){
    cRdmEff <- pars[iParNext:(iParNext + ncol(rdmEffFactors) - 1)]
    if(min(cRdmEff) < 0){ paramsOK <- F }else{}
  }else{
    cRdmEff <- c()
  }
  
  parsOut <- paste0('c0=' , round(c0, digits = 3) , ', c1=' , round(c1, digits = 3) , ', a1=' , round(a1, digits = 2))
  if(nSpatStructs == 2){ parsOut <- paste0(parsOut , ', c2=' , round(c2, digits = 3) , ', a2=' , round(a2, digits = 2)) }else{}
  if(!is.null(rdmEffFactors)){
    parsOut <- paste0(parsOut , ', cRdmEff=[' , paste(round(cRdmEff, digits = 3) , collapse = '; ') , ']')
  }else{}
  
  ### default lists to be returned when some error arises...
  if(returnAll){
    listOut <- list('nll' = badNll , 'betahat' = NA , 'vbetahat' = NA , 'pars' = pars , 'covParams4Kriging' = NA , 'C' = NA , 'iC' = NA)
  }else{             
    listOut <- badNll
  }

  if(paramsOK){
    if(timeCodeNow){ start_time <- Sys.time() }else{}
    C <- defineCLMM2(c0 = c0 , c1 = c1 , a1 = a1 , c2 = c2 , a2 = a2 , cRdmEff = cRdmEff , D = D , nAliquots = nAliquots , rdmEffFactors = rdmEffFactors , covModel = covModel , nSpatStructs = nSpatStructs)
    if(timeCodeNow){ print('Time to define cov mtx C:') ; print(Sys.time() - start_time) }else{}
  }else{
    return(listOut)
  }   
  
  if(timeCodeNow){ start_time <- Sys.time() }else{}
  listRtn <- fLMM2_Given_C_Bayes(C = C , z = z , X = X , nAliquots = nAliquots , priorInfo = priorInfo , 
                                 returnAll = returnAll , optionML = optionML , verbose = verbose , forCompLik = FALSE , 
                                 badNll = badNll , attachBigMats = attachBigMats , 
                                 pars = pars , parsOut = parsOut , listOut = listOut)
  if(timeCodeNow){ print('Time to calculate nll:') ; print(Sys.time() - start_time) }else{}
  
  if(returnAll){
    ### update the covMode,covParams4Kriging,nSpatStructs...
    listRtn$covModel <- covModel
    listRtn$covParams4Kriging <- c(c0 , c1 , a1 , c2 , a2)
    if(length(cRdmEff) > 0){ listRtn$covParams4Kriging <- c(listRtn$covParams4Kriging , cRdmEff) }else{}
    listRtn$nSpatStructs <- nSpatStructs
  }else{}
  
  return(listRtn)
}

###################################################################################
### with the covariance matrix C defined, return the fLMM2 output (with priorInfo given)...
###################################################################################
fLMM2_Given_C_Bayes <- function(C , z , X , nAliquots = 1 , priorInfo = NULL , returnAll , optionML , verbose , forCompLik = FALSE , badNll = 9E99 , attachBigMats = TRUE , 
                                pars , parsOut , listOut){
  N <- length(z)
  np <- dim(X)[2]
  
  if(length(nAliquots) > 1){ if(length(nAliquots) != length(z)){ stop('Error - nAliquots must be same length as z for fitLMM2!') }else{} }else{}
  if(missing(returnAll)){ returnAll <- F }else{} # default is to return just nll
  if(missing(optionML)){ optionML <- F }else{} # default is REML
  if(missing(verbose)){ verbose <- T }else{} # default is to print the current nll
  
  oneszX <- cbind(matrix(1 , N , 1) , z , X)
  
  tmp <- lndetANDinvCb(C , oneszX)
  
  if(returnAll){ cholC <- chol(C) }else{}
  
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
  
  if(!is.null(priorInfo)){
    zinvCz <- zinvCz + priorInfo$beta0pbeta0beta0
    zinvCX <- zinvCX + t(priorInfo$pbeta0beta0)
    XinvCX <- XinvCX + priorInfo$pbeta0
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
    nll <- -0.5 * N * log(2 * pi) - 0.5 * lndetC - 0.5 * resiCres
  }else{
    nll <- -0.5 * (N-np) * log(2 * pi) - 0.5 * lndetC - 0.5 * lndetXinvCX - 0.5 * resiCres
  }
  nll <- -nll
  
  printNll(nll = nll , parsOut = parsOut , verbose = verbose)
  
  if(is.na(nll)){ nll <- badNll }else{}
  
  if(returnAll){
    vbetahat <- chol2inv(chol(XinvCX)) 
    if(attachBigMats){
      iC <- chol2inv(cholC)
      if((!is.null(priorInfo)) && (!is.null(priorInfo$beta0)) && (max(abs(priorInfo$beta0)) > 0)){ print('Note - the returned matrices XiCX etc are actually XiCX + pbeta0 etc (ie Bayes updated versions)') }else{}
      
      listRtn <- list('nll' = nll , 'betahat' = betahat , 'vbetahat' = vbetahat , 'pars' = pars , 
                      'covModel' = NA , 'covParams4Kriging' = NA , 'nSpatStructs' = NA , 
                      'C' = C , 'iC' = iC , 'XiCX' = XinvCX , 'XiCz' = t(zinvCX) , 'ziCz' = zinvCz , 'lndetC' = lndetC , 
                      'iC' = iC , 'iCX' = iConeszX[,-c(1,2),drop=FALSE] , 'iCz' = iConeszX[,2])
      
      ### add other info, for consistency with the noPrior version...
      if((length(nAliquots) == 1) & (nAliquots > 1)){
        stop('ATTENTION - nAliquots > 1, SO CANNOT TELL AT THIS STAGE WHAT THE VARIANCE IS!')
      }else if((length(nAliquots) == 1) & (nAliquots == 1)){
        i1 <- 1
      }else if((length(nAliquots) > 1) & (min(nAliquots) > 1)){
        stop('ATTENTION - nAliquots > 1, SO CANNOT TELL AT THIS STAGE WHAT THE VARIANCE IS!')
      }else if((length(nAliquots) > 1) & (min(nAliquots) == 1)){
        i1 <- which(nAliquots == 1)[1]
      }else{
        stop('ALL POSSIBILITIES EXHAUSTED!')
      }
      listRtn$sigma2hat <- listRtn$C[i1,i1]
      
      listRtn$iA <- listRtn$iC / listRtn$sigma2hat
      listRtn$iAX <- listRtn$iCX / listRtn$sigma2hat
      listRtn$iAz <- listRtn$iCz / listRtn$sigma2hat

      listRtn$XiAX <- listRtn$XiCX / listRtn$sigma2hat
      listRtn$XiAz <- listRtn$XiCz / listRtn$sigma2hat
      listRtn$ziAz <- listRtn$ziCz / listRtn$sigma2hat
      
      listRtn$lndetA <- listRtn$lndetC - N * log(listRtn$sigma2hat)

      return(listRtn)
    }else{
      if((!is.null(priorInfo)) && (!is.null(priorInfo$beta0)) && (max(abs(priorInfo$beta0)) > 0)){ print('Note - the returned matrices XiCX etc are actually XiCX + pbeta0 etc (ie Bayes updated versions)') }else{}
      return(list('nll' = nll , 'betahat' = betahat , 'vbetahat' = vbetahat , 'pars' = pars , 
                  'covModel' = NA , 'covParams4Kriging' = NA , 'nSpatStructs' = NA , 
                  'XiCX' = XinvCX , 'XiCz' = t(zinvCX) , 'ziCz' = zinvCz))
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

        if((!is.na(sigma2hat)) && (sigma2hat > 0)){
          
          if (optionML){
              nll <- 0.5 * (n * log(2 * pi) + n * log(sigma2hat) + sumlndetA + n) # [last term because (1 / sigma2hat) * resiAres  = n]
          }else{
              nll <- 0.5 * ((n - p) * log(2 * pi) + (n - p) * log(sigma2hat) + sumlndetA + lndetXiAX + n - p) # [last term because (1 / sigma2hat) * resiAres  = n - p]
          }
          
        }else{

          nll <- badNll
          
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

###################################################################
### kringing for a fitted LMM2 model...
### same version of this function for Bayes/no prior given
### must pass in either 
###    [betahat/vbetahat & not beta0/pbeta0 = no prior used or prior info already built into betahat/vbetahat - betahat/vbetahat fixed]
###    [not betahat/vbetahat & not beta0/pbeta0 = no prior used - betahat/vbetahat re-estimated]
###    [not betahat/vbetahat & beta0/pbeta0 = prior used - betahat/vbetahat re-estimated]
### if all betahat/vbetahat & beta0/pbeta0 passed in = error!
###################################################################
krigingLMM2 <- function(ck , Xk , nAliquotsk = 1 , rdmEffFactorsk = NULL , c , z = NULL , X = NULL , nAliquots = 1 , rdmEffFactors = NULL , D = NULL , beta0 = NULL , pbeta0 = NULL , covModel , nSpatStructs , covParams , 
                              betahat = NULL , vbetahat = NULL , invC = NULL , invCz = NULL , invCX = NULL , 
                              predVarsRqd , fullCovMtx = FALSE , blockLength = 0 , verbose = TRUE , optReturn = 0 , nPerChunk = 5000){
  
  ### if betahat passed in, that is used
  ### if not (ie NULL), then betahat estimated from the data passed into fn.
  ### could be different if for eg the data passed in to this function are not identical to the data used to fit the covariance model
  ### optReturn = 0 ; preds and predvars (or predVar mtx)
  ### optReturn = 1 ; also CkIFh_SK = Ck - Ckh iCh Chk and CkhiChzh = Ckh iCh zh and CkhiChXh = Ckh iCh Xh
  
  if(is.null(z) | is.null(X) | is.null(D)){
    ### in this case, must pass in all of betahat = NULL , vbetahat = NULL , invC = NULL , invCz = NULL , invCX = NULL
    if(is.null(betahat) | is.null(vbetahat) | is.null(invC) | is.null(invCz) | is.null(invCX)){
      stop('Error in krigingLMM2 - if not passing in z/X/D, then must pass in all of betahat/vbetahat/invC/invCz/invCX!')
    }else{}
  }else{}
  
  if((!is.null(betahat)) & (!is.null(beta0))){
    stop('Error - you have passed in an estimated betahat as well as a prior for beta (beta0) - unclear whether you have already incorporated prior info in the given betahat, so stopping now!')
  }else{}
  
  if((!is.null(beta0)) & (is.null(pbeta0))){ stop('Error - beta0 given but not pbeta0!') }else{}
  if((is.null(beta0)) & (!is.null(pbeta0))){ stop('Error - pbeta0 given but not beta0!') }else{}
  
  if(is.data.frame(c)){
    c <- as.matrix(c)
  }else{}
  
  if(is.data.frame(ck)){
    ck <- as.matrix(ck)
  }else{}
  
  ndimc <- ncol(c)
  nk <- nrow(ck)
  if(is.null(ndimc)){
    ndimc <- 1
    nk <- length(ck)
  }else{}

  ### check and format rdmEffFactors and rdmEffFactorsk...
  if((is.null(rdmEffFactors)) & !is.null(rdmEffFactorsk)){ stop('Error - rdmEffFactorsk given but not rdmEffFactors!') }else{}
  if((!is.null(rdmEffFactors)) & is.null(rdmEffFactorsk)){ stop('Error - rdmEffFactors given but not rdmEffFactorsk!') }else{}
  if(!is.null(rdmEffFactors)){
    if(is.matrix(rdmEffFactors)){ rdmEffFactors <- as.data.frame(rdmEffFactors) }else{}
    if(is.data.frame(rdmEffFactors)){
      if(nrow(rdmEffFactors) != length(z)){ stop('Error - rdmEffFactors must be a matrix with the same number of rows as you have data!') }else{} 
    }else{
      if(length(rdmEffFactors) == length(z)){ 
        rdmEffFactors <- data.frame('rdmEff1' = rdmEffFactors) 
      }else{
        stop('Error - rdmEffFactors must be a matrix with the same number of rows as you have data!')        
      }
    }
    
    if(is.matrix(rdmEffFactorsk)){ rdmEffFactorsk <- as.data.frame(rdmEffFactorsk) }else{}
    
    if(is.data.frame(rdmEffFactorsk)){
      if(nrow(rdmEffFactorsk) != nrow(ck)){ stop('Error - rdmEffFactors must be a matrix with the same number of rows as you have data!') }else{} 
    }else{
      if(length(rdmEffFactorsk) == nrow(ck)){ 
        rdmEffFactorsk <- data.frame('rdmEff1' = rdmEffFactorsk) 
      }else{
        stop('Error - rdmEffFactors must be a matrix with the same number of rows as you have data!')        
      }
    }
    ### make sure all are factors, and that rdmEffFactorsk has the same levels as rdmEffFactors...
    for(i in 1:ncol(rdmEffFactors)){
      if(!is.factor(rdmEffFactors[,i])){ rdmEffFactors[,i] <- factor(rdmEffFactors[,i]) }else{}
      rdmEffFactorsk[,i] <- factor(rdmEffFactorsk[,i] , levels = levels(rdmEffFactors[,i])) 
    }
  }else{}
  
  ### prediction block centred at 0,0
  if(blockLength > 0){
    if(ncol(Xk) > 1){ print('Attention - krigingLMM2 is doing block kriging with trend - trend assumed constant within blocks!') }else{}
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

  if(!is.null(rdmEffFactors)){ 
    if(length(covParams) != (5 + ncol(rdmEffFactors))){ stop('Error - rdmEffFactors given, but parameters in covParams do not line up properly! (Should be 5 spatial covariance model parameters + the rest as other random effect parameters.)') }else{}
    cRdmEff <- covParams[6:length(covParams)] 
  }else{
    cRdmEff <- c()
  }
  
  if(verbose){
    print(paste0('c0 = ' , c0))
    print(paste0('c1 = ' , c1))
    print(paste0('a1 = ' , a1))
    if (nSpatStructs == 2){
      print(paste0('c2 = ' , c2))
      print(paste0('a2 = ' , a2))
    }else{}
    if(length(cRdmEff) > 0){
      print(paste0('cRdmEff = [' , paste(cRdmEff , collapse = ' , ') , ']'))
    }else{}
  }else{}
  
  if(is.null(invC) | is.null(invCz) | is.null(invCX) | is.null(betahat) | is.null(vbetahat)){
    C <- defineCLMM2(c0 = c0 , c1 = c1 , a1 = a1 , c2 = c2 , a2 = a2 , cRdmEff = cRdmEff , D = D , nAliquots = nAliquots , rdmEffFactors = rdmEffFactors , covModel = covModel , nSpatStructs = nSpatStructs)

    zX <- cbind(z , X)
    tmp <- lndetANDinvCb(C , zX)
    
    lndetC <- tmp$lndetC
    invCz <- tmp$invCb[,1 , drop = FALSE]
    invCX <- tmp$invCb[,2:(1+np) , drop = FALSE]
    
    # invC <- chol2inv(tmp$cholC)
    invC <- solve(C)
    
    XinvCX <- t(X) %*% invCX
    XinvCz <- t(X) %*% invCz
    
    if(!is.null(beta0)){
      XinvCX <- XinvCX + pbeta0
      XinvCz <- XinvCz + pbeta0 %*% beta0
    }else{}
    
    vbetahat <- solve(XinvCX)
    betahat <- vbetahat %*% XinvCz

  }else{}
  
  invCRes <- invCz - invCX %*% betahat    
  
  if(blockLength == 0){
    nChunks <- ceiling(nk / nPerChunk)
  }else{
    nPerChunk <- 1
    nChunks <- nk
  }
  
  pred <- matrix(NA,nk,1)
  if(predVarsRqd){ predVars <- matrix(NA,nk,1) }else{}
  if(optReturn == 1){
    if(predVarsRqd){ CkIFh_SK <- matrix(NA,nk,1) }else{}
    CkhiChzh <- matrix(NA,nk,1)
    CkhiChXh <- matrix(NA,nk,np)
  }else{}
  
  if(fullCovMtx){ 
    if(nChunks > 1){
      stop('Returning full prediction covariance matrix only possible if less than 5000 prediction locations given!')
    }else{}
    if(predVarsRqd){
      predVars <- matrix(NA,nk,nk)
      if(optReturn == 1){
        CkIFh_SK <- matrix(NA,nk,nk)
      }else{}
    }else{}
  }else{}
  
  ###############################################
  ### for each chunk...
  ###############################################
  for (iChunk in 1:nChunks){
    
    if(verbose && is.element(iChunk , round(quantile(seq(nChunks) , seq(0.05 , 0.95 , 0.05))))){ print(paste0(iChunk , ' of ' , nChunks)) }else{}
    
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
    
    Dkh <- xyDist(ckThis,c)
    
    ### pass in nAliquots = 1 here, as this is between prediction locations and data locations, so nAliquots not used.
    Ckh <- defineCLMM2(c0 = c0 , c1 = c1 , a1 = a1 , c2 = c2 , a2 = a2 , cRdmEff = cRdmEff , D = Dkh , nAliquots = 1 , 
                       rdmEffFactors = rdmEffFactorsk[ikThis,,drop=FALSE] , rdmEffFactors2 = rdmEffFactors , covModel = covModel , nSpatStructs = nSpatStructs)

    if(blockLength == 0){
      if(optReturn > 0){
        CkhiChzh[ikThis] <- Ckh %*% invCz
        CkhiChXh[ikThis,] <- Ckh %*% invCX
      }else{}
      pred[ikThis] <- Xk[ikThis, , drop = FALSE] %*% betahat + Ckh %*% invCRes
    }else{
      if(optReturn > 0){
        CkhiChzh[ikThis] <- mean(Ckh %*% invCz)
        CkhiChXh[ikThis,] <- colMeans(Ckh %*% invCX)
      }else{}
      ### orig note : needs to be generalised in trend is not a const mean.      
      ### updated : if we have trend, this will assume trend constant within the block supports.      
      pred[ikThis] <- Xk[ikThis, , drop = FALSE] %*% betahat + mean(Ckh %*% invCRes)
    }
    
    ### also calculate the prediction variances; but only if rqd...
    if(predVarsRqd){
      if(fullCovMtx | (blockLength > 0)){
        Dk <- xyDist(ckThis,ckThis)
        
        ### pass in nAliquots = nAliquotsk here, for prediction of composite of nAliquotsk aliquots.
        Ck <- defineCLMM2(c0 = c0 , c1 = c1 , a1 = a1 , c2 = c2 , a2 = a2 , cRdmEff = cRdmEff , D = Dk , nAliquots = nAliquotsk , rdmEffFactors = rdmEffFactorsk[ikThis,,drop=FALSE] , covModel = covModel , nSpatStructs = nSpatStructs)

        # Xk_CkhiChXh <- Xk[ikThis,] - Ckh %*% invCX
        if(length(ikThis) < nrow(ckThis)){
          Xk_CkhiChXh <- Xk[rep(ikThis , nrow(ckThis)),,drop=FALSE] - Ckh %*% invCX
        }else{
          Xk_CkhiChXh <- Xk[ikThis,,drop=FALSE] - Ckh %*% invCX
        }
        
        varTrend <- (Xk_CkhiChXh %*% vbetahat) %*% t(Xk_CkhiChXh)
        
        if(blockLength == 0){
          if(optReturn == 0){
            ### using as.matrix to avoid errors if matrices are sparse...
            predVars[ikThis,ikThis] <- as.matrix(Ck - (Ckh %*% invC) %*% t(Ckh) + varTrend)
          }else{
            ### using as.matrix to avoid errors if matrices are sparse...
            CkIFh_SK[ikThis,ikThis] <- as.matrix(Ck - (Ckh %*% invC) %*% t(Ckh))
            predVars[ikThis,ikThis] <- as.matrix(CkIFh_SK[ikThis,ikThis] + varTrend)
          }
        }else{
          if(optReturn == 0){
            predVars[ikThis] <- mean(Ck - (Ckh %*% invC) %*% t(Ckh)) + mean(varTrend)
          }else{
            CkIFh_SK[ikThis] <- mean(Ck - (Ckh %*% invC) %*% t(Ckh))
            predVars[ikThis] <- CkIFh_SK[ikThis] + mean(varTrend)
          }
          # predVars[ikThis] <- mean(Ck - (Ckh %*% invC) %*% t(Ckh) + varTrend)
        }
      }else{
        Ck <- c0 + c1 + c2 + sum(cRdmEff)
        Xk_CkhiChXh <- Xk[ikThis,,drop=FALSE] - Ckh %*% invCX
        varTrend <- rowSums((Xk_CkhiChXh %*% vbetahat) * Xk_CkhiChXh)

        if(optReturn == 0){
          predVars[ikThis] <- Ck - rowSums((Ckh %*% invC) * Ckh) + varTrend
        }else{
          CkIFh_SK[ikThis] <- Ck - rowSums((Ckh %*% invC) * Ckh)
          predVars[ikThis] <- CkIFh_SK[ikThis] + varTrend
        }
      }
    }
  } 
  
  if(optReturn == 0){
    return(list('pred' = pred , 'predVars' = predVars))
  }else if(optReturn == 1){
    return(list('pred' = pred , 'predVars' = predVars , 
                'CkIFh_SK' = CkIFh_SK , 'CkhiChzh' = CkhiChzh , 'CkhiChXh' = CkhiChXh))
  }else{
    stop('Unknown optReturn for krigingLMM2!')
  }
}

###################################################################
### for the compLik approach...
### only for point-support prediction.
###################################################################
compLikKrigingLMM2 <- function(ck , Xk , blocksk , c , z , X , nAliquots = 1 , blocks , DBlocks , nBlocks1 = length(blocks) , iCres , iCX , vbetahat , covModel , nSpatStructs , covParams , betahat , nsdsPlot){
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
    
      # DkhThis <- rdist(ck[blocksk[[i]],] , c[blocks[[i]],])
      DkhThis <- xyDist(ck[blocksk[[i]],] , c[blocks[[i]],])
      ### pass in nAliquots = 1 here, as this is between prediction locations and data locations, so nAliquots not used.
        CkhThis <- defineCLMM2(c0 = covParams[1] , c1 = covParams[2] , a1 = covParams[3] , c2 = covParams[4] , a2 = covParams[5] , 
              D = DkhThis , nAliquots = 1 , covModel = covModel , nSpatStructs = nSpatStructs)

        if(length(nAliquots) > 1){ nAliquotsThis <- nAliquots[blocks[[i]]] }else{ nAliquotsThis <- nAliquots }
        ChThis <- defineCLMM2(c0 = covParams[1] , c1 = covParams[2] , a1 = covParams[3] , c2 = covParams[4] , a2 = covParams[5] , 
              D = DBlocks[[i]] , nAliquots = nAliquotsThis , covModel = covModel , nSpatStructs = nSpatStructs)

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
### if betahat and vbetahat are not given (NULL), then these are re-estimated each time a subset is removed
### note - nAliquots not needed as argument here, as already incorporated in iC
### same version of this function for Bayes/no prior given
### must pass in either 
###    [betahat/vbetahat & not beta0/pbeta0 = no prior used or prior info already built into betahat/vbetahat - betahat/vbetahat fixed]
###    [not betahat/vbetahat & not beta0/pbeta0 = no prior used - betahat/vbetahat re-estimated each time]
###    [not betahat/vbetahat & beta0/pbeta0 = prior used - betahat/vbetahat re-estimated each time]
### if all betahat/vbetahat & beta0/pbeta0 passed in = error!
####################################################
XVLMM2 <- function(c , z , X , iC , betahat = NULL , vbetahat = NULL , beta0 = NULL , pbeta0 = NULL , nsdsPlot = 1.96 , iSubsets = NULL , rtnFullCovMtxAsList = FALSE){ 
  
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
  if(!is.null(betahat)){
    if(!is.null(beta0)){ stop('Error - you have defined an estimate of betahat, and also a prior for beta, beta0! Define either neither (to re-estimate each time), just beta0 (to re-estimate each time with prior beta0, precision pbeta0), or just betahat (fixed throughout)!') }else{}
    mu <- X %*% betahat
    res <- z - mu
  }else{}
  
  ### put prior info matrices together in list...
  if(!is.null(beta0)){ 
    beta0 <- matrix(beta0 , ncol = 1)
    if(is.null(pbeta0)){ stop('Error - you have entered a prior beta as beta0, but not an associated prior precision, pbeta0!') }else{}
    
    priorInfo <- list('beta0' = beta0 , 'pbeta0' = pbeta0 , 
                      'pbeta0beta0' = pbeta0 %*% beta0)
    priorInfo$beta0pbeta0beta0 <- t(priorInfo$beta0) %*% priorInfo$pbeta0beta0
  }else{
    if(!is.null(pbeta0)){ stop('Error - you have entered a prior precision pbeta as pbeta0, but not an associated prior estimate, beta0!') }else{}
    priorInfo <- NULL
  }
  
  n <- dim(X)[[1]]
  p <- dim(X)[[2]]
  
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
  if(rtnFullCovMtxAsList){
    zkList <- vkList <- list()
  }else{}
  
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
      if(is.null(betahat)){
        iCThis <- iC[-iThis,-iThis,drop=FALSE] - iC[-iThis,iThis,drop=FALSE] %*% solve(iC[iThis,iThis,drop=FALSE]) %*% iC[iThis,-iThis,drop=FALSE]
        iCXzhThis <- iCThis %*% cbind(X[-iThis,,drop=FALSE] , z[-iThis])
        
        XzhiCXzhThis <- t(cbind(X[-iThis,,drop=FALSE] , z[-iThis])) %*% iCXzhThis
        if(is.null(priorInfo)){
          vbetahatThis <- solve(XzhiCXzhThis[1:p,1:p,drop=FALSE])
          betahatThis <- vbetahatThis %*% XzhiCXzhThis[1:p,p+1,drop=FALSE]
        }else{
          XiCXTmp <- XzhiCXzhThis[1:p,1:p,drop=FALSE] + priorInfo$pbeta0
          XiCzTmp <- XzhiCXzhThis[1:p,p+1,drop=FALSE] + priorInfo$pbeta0beta0
          
          vbetahatThis <- solve(XiCXTmp)
          betahatThis <- solve(XiCXTmp , XiCzTmp)
        }
        
        mukThis <- X[iThis,,drop=FALSE] %*% betahatThis
        reshThis <- z[-iThis] - X[-iThis,,drop=FALSE] %*% betahatThis
      }else{
        mukThis <- mu[iThis]
        reshThis <- res[-iThis]
        vbetahatThis <- vbetahat
      }
      CkhiChh <- -solve(iC[iThis,iThis,drop=FALSE] , iC[iThis,-iThis,drop=FALSE]) 
      
      tmp <- CkhiChh %*% cbind(reshThis , X[-iThis,,drop=FALSE])
      CkhiChhresh <- tmp[,1,drop=FALSE]
      CkhiChhXh <- tmp[,-1,drop=FALSE]
      
      Xk_CkhiChhXh <- X[iThis,,drop=FALSE] - CkhiChhXh
      vSK <- solve(iC[iThis,iThis,drop=FALSE])
      
      zkFull <- mukThis + CkhiChhresh
      vkFull <- vSK + Xk_CkhiChhXh %*% vbetahatThis %*% t(Xk_CkhiChhXh)
      
      if(rtnFullCovMtxAsList){
        zkList[[i]] <- zkFull
        vkList[[i]] <- as.matrix(vkFull)
      }else{}
      
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
  
  listRtn <- list('zk' = zk , 'vk' = vk , 'err' = err , 'serr' = serr , 'sse' = sse , 'polyvk' = polyvk , 
                  'zbar' = zbar , 'nPerSubset' = nPerSubset , 'zkFullFinal' = zkFull , 'vkFullFinal' = vkFull)
  if(rtnFullCovMtxAsList){
    listRtn$zkList <- zkList
    listRtn$vkList <- vkList
  }else{}
  
  return(listRtn)
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
### if rdmEffFactors2 is given in addition to rdmEffFactors, then 
### the fn gives cov between distinct data (e.g. between data and prediction locations); (D would be non-symmetric)
### else, assumed symmetric, and rdmEffFactors2 = rdmEffFactors
##########################################################################
defineCLMM2 <- function(c0 , c1 , a1 , c2 , a2 , cRdmEff = c() , D , nAliquots = 1 , rdmEffFactors = NULL , rdmEffFactors2 = NULL , covModel , nSpatStructs){

    # iDPossErrs <- which((D > 0) & (D <=1E-10))
    # if(length(iDPossErrs) > 0){ print('WARNING - SOME SMALL VALUES IN D COULD BE COLOCATED - CHECK DEFINITION OF D!') }else{}

    if(ncol(D) != nrow(D)){
      if(length(nAliquots) > 1){ 
        stop('Error - put nAliquots = 1 when defining non-symmetric covariance matrix - nAliquots only used when defining symmetric covariance matrices')
      }else if(nAliquots > 1){
        stop('Error - put nAliquots = 1 when defining non-symmetric covariance matrix - nAliquots only used when defining symmetric covariance matrices')
      }else{}
    }else{}

    if(covModel == 'exponential'){ 
        C <-  c1 * exp(-3 * D / a1) 
        if(nSpatStructs == 2){ C <- C + c2 * exp(-3 * D / a2) }else{}
    }else if(covModel == 'spherical'){ 
      if(is(D , 'dsCMatrix')){
        C <- D
        C@x <- c1 * (1 - (1.5 * (D@x / a1) - 0.5 * ((D@x / a1) ^ 3) ))
        C@x[D@x > a1] <- 0

        if(nSpatStructs == 2){ 
          C@x <- C@x + c2 * (1 - (1.5 * (D@x / a2) - 0.5 * ((D@x / a2) ^ 3) ))
          C@x[D@x > a2] <- 0
        }else{}

      }else{
        DOVERa1 <- D / a1
        DOVERa2 <- D / a2
        C1 <- 1 - (1.5 * DOVERa1 - 0.5 * (DOVERa1 ^ 3) )
        C1[which(DOVERa1 > 1)] <- 0
        C <- c1 * C1
        
        if(nSpatStructs == 2){ 
          C2 <- 1 - (1.5 * DOVERa2 - 0.5 * (DOVERa2 ^ 3) )
          C2[which(DOVERa2 > 1)] <- 0
          C <- C + c2 * C2
        }else{}
      }

    }else if (covModel == 'gaussian'){ 
        C <-  c1 * exp(-((sqrt(3) * D / a1) ^ 2)) 
        if (nSpatStructs == 2){ C <- C + c2 * exp(-((sqrt(3) * D / a2) ^ 2)) }else{}

    }else if (covModel == 'exponentialgaussian'){
### exp for short range, gaussian for long range...
        if (nSpatStructs != 2){ stop('Error, for exponentialgaussian cov model, must put nSpatStructs = 2') }else{}        
        C <-  c1 * exp(-3 * D / a1) 
        C <- C + c2 * exp(-((sqrt(3) * D / a2) ^ 2)) 

    }else if(substr(covModel , 1 , 6) == 'matern'){ 
### note - with nu = 0.5 will be different to exponential, cos of different parameterization. 
### sensible in future to align all parameterizations (st a is effective range)
### although if using in iak work, attention that iaCov functions are based on the distance (not eff range) parameterization
        nu <- as.numeric(substr(covModel , 7 , nchar(covModel)))
        C <-  c1 * maternCov4fLMM2(D = D , pars = c(1 , a1 , nu)) 
        if(nSpatStructs == 2){ C <- C + c2 * maternCov4fLMM2(D = D , pars = c(1 , a2 , nu)) }else{}

    }else if(substr(covModel , 1 , 8) == 'wendland'){ 
      nu <- as.numeric(substr(covModel , 9 , nchar(covModel)))
      if(is(D , 'dsCMatrix')){
        C <- D
        C@x <- wendlandCov4fLMM2(D = D@x , pars = c(c1 , a1 , nu))
        if(nSpatStructs == 2){ 
          C@x <- C@x + wendlandCov4fLMM2(D = D@x , pars = c(c2 , a2 , nu))
        }else{}
        
      }else{
        C <- wendlandCov4fLMM2(D = D , pars = c(c1 , a1 , nu))
        if(nSpatStructs == 2){ 
          C <- C + wendlandCov4fLMM2(D = D , pars = c(c2 , a2 , nu))
        }else{}
      }
    }else if(substr(covModel , 1 , 9) == 'gwendland'){ 
      nu <- as.numeric(substr(covModel , 10 , nchar(covModel)))
      if(is(D , 'dsCMatrix')){
        C <- D
        C@x <- gwendlandCov4fLMM2(D = D@x , pars = c(c1 , a1 , nu))
        if(nSpatStructs == 2){ 
          C@x <- C@x + gwendlandCov4fLMM2(D = D@x , pars = c(c2 , a2 , nu))
        }else{}
        
      }else{
        C <- gwendlandCov4fLMM2(D = D , pars = c(c1 , a1 , nu))
        if(nSpatStructs == 2){ 
          C <- C + gwendlandCov4fLMM2(D = D , pars = c(c2 , a2 , nu))
        }else{}
      }
    }else{}

    ### add nugget, divided by nAliquots 
    if(max(nAliquots) > 1){
      ### assuming symmetric D with the only D = 0 values on diagonal
      diag(C) <- c1 + c2 + c0/nAliquots
    }else if(is(D , 'dsCMatrix')){
      ### sparse symmetric, must be for rows-cols = same locations
      diag(C) <- c1 + c2 + c0/nAliquots
    }else{
      ### not necessarily symmetric - could be between data and prediction locations, with potential co-location
      C <- C + c0 * (D == 0)
    }

    if(!is.null(rdmEffFactors)){
      for(i in 1:ncol(rdmEffFactors)){
        for(j in 1:nlevels(rdmEffFactors[,i])){
          iThis <- which(rdmEffFactors[,i] == levels(rdmEffFactors[,i])[j])
          if(is.null(rdmEffFactors2)){
            iThis2 <- iThis
          }else{
            iThis2 <- which(rdmEffFactors2[,i] == levels(rdmEffFactors[,i])[j])
          }
          if((length(iThis) > 0) & (length(iThis2) > 0)){
            C[iThis,iThis2] <- C[iThis,iThis2] + cRdmEff[i]
          }else{}
        }
      }
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

##################################################################
### the wendland covariance function...
##################################################################
wendlandCov4fLMM2 <- function(D , pars){
  if(!exists('Wendland')){
    stop('Error - source iaCovMatern or load the fields library to load the Wendland function')
  }else{}
  
  c1 <- pars[1]
  a <- pars[2] # the range, or theta
  nu <- pars[3] # 
  if((c1 > 0) & (a > 0) & (nu >= 1)  & (nu <= 200)){
    
    if(nu == round(nu)){
      C <- Wendland(D, theta = a , k = nu, dimension = 2)
    }else{
      C <- (ceiling(nu) - nu) * Wendland(D, theta = a , k = floor(nu), dimension = 2)
      C <- C + (nu - floor(nu)) * Wendland(D, theta = a , k = ceiling(nu), dimension = 2)
    }

    C <- c1 * C
    
  }else{
    C <- NA * D
  }  
  
  return(C)
}

##################################################################
### a version of the generalized wendland covariance function...
### with parameter mu restricted to its lower boundary value (mu = ((nd+1) / 2 + nu))
##################################################################
gwendlandCov4fLMM2 <- function(D , pars , nd = 2){

  require(gsl) # for the hyperg_2F1 function (hypergeometric gaussian function)
  # from belivaqua 2022
  # nd = number of spatial dims
  # here, putting mu = (nd+1) / 2 + nu, which is the equal to its boundary constraint
  # seems to allow mu to still give variety of smoothness and doesn't get numerical errors in initial tests (which i got with mu greather than this boundary)
  
  c1 <- pars[1]
  a <- pars[2] # the range for compact support - beta in Bevilacqua
  nu <- pars[3] # smoothness,
  
  # mu <- pars[3] # ? , pos def if mu >= (nd+1) / 2 + nu, nd = number of spatial dims
  # & (mu >= ((nd+1) / 2 + nu))
  mu <- ((nd+1) / 2 + nu)

  if((c1 > 0) & (a > 0) & (nu >= 0)  & (nu <= 200)){
    
    if(nu == 0){
      C <- 0 * D
      C[D > 0 & D < a] <- (1 - D[D > 0 & D < a] / a) ^ mu
      C[D == 0] <- 1
    }else if (nu > 0){
      C <- 0 * D
      K <- exp(lngamma(nu) + lngamma(2 * nu + mu + 1) - lngamma(2 * nu) - lngamma(nu + mu + 1) - (mu + 1) * log(2))
      C[D > 0 & D < a] <- K * ((1 - (D[D > 0 & D < a] / a) ^ 2) ^ (nu + mu)) * hyperg_2F1(mu/2,(mu+1)/2,nu+mu+1,1 - (D[D > 0 & D < a] / a) ^ 2)
      C[D == 0] <- 1
    }else{
      C <- NA * D
    }
    C <- c1 * C
    
  }else{
    C <- NA * D
  }  

  return(C)
}

#############################################
### natural spline forward selection (NSFS)
### can add cts variable
###     add knot to cts variable (as nat spline)
###     multiply two columns of current X
###     add cat variable
###     mutiply cat variable (all cols) by existing cts col
### perhaps do with ML
### 1. fit model X0 with const mean + 2 structs, a2<= 0.5 * maxD
### 2. select 3 candidate additions using LRT with fixed s0,s1,s2,a1,a2; free sigma2hat
### 3. properly refit these models and use LRT 
### 4. repeat 2 and 3 - would expect a2 to decrease as more variance explained, so always put maxa = a2 of previous fit
### [or fix a2 at eg 100 km throughout?]
#############################################
modelSelectLMM2_NSFS <- function(c , z , dfCovs = NULL , namesCovsStart = NULL , nAliquots = 1 , rdmEffFactors = NULL , covModel = 'exponential' , nSpatStructs = 1 , incNugget = T , optionML = F , verbose_NSFS = FALSE , mina = NULL , maxa = NULL , parsInit = NULL , maxnIntKnots = 3 , maxnInts = 2 , prefer_fewer_covs = FALSE){

  fixa2 <- F

  if(is.null(dfCovs)){ stop('Error - must include dfCovs in call to modelSelectLMM2_NSFS!') }else{}
  if(!is.data.frame(dfCovs)){ stop('Error - dfCovs must be a data.frame in call to modelSelectLMM2_NSFS!') }else{}
  if((!is.null(namesCovsStart)) && (length(namesCovsStart) > 0) && (!is.element(namesCovsStart , names(dfCovs)))){ stop('Error - namesCovsStart must be in dfCovs!') }else{}

  #########################################
  ### nAliquots used as if these cores were sampled a negligible distance from each other
  ### so diagonal of cov mtx for entries where nAliquots > 1 = c1 + c2 + c0/nAliquots
  #########################################
  if(length(nAliquots) > 1){ if(length(nAliquots) != length(z)){ stop('Error - nAliquots must be same length as z for fitLMM2!') }else{} }else{}
  
  #########################################
  ### check rdm eff factors...
  #########################################
  if(!is.null(rdmEffFactors)){
    if(is.matrix(rdmEffFactors)){ rdmEffFactors <- as.data.frame(rdmEffFactors) }else{}
    
    if(is.data.frame(rdmEffFactors)){
      if(nrow(rdmEffFactors) != length(z)){ stop('Error - rdmEffFactors must be a matrix with the same number of rows as you have data!') }else{} 
    }else{
      if(length(rdmEffFactors) == length(z)){ 
        rdmEffFactors <- data.frame('rdmEff1' = rdmEffFactors) 
      }else{
        stop('Error - rdmEffFactors must be a matrix with the same number of rows as you have data!')        
      }
    }
    
    ### make sure all are factors...
    for(i in 1:ncol(rdmEffFactors)){
      if(!is.factor(rdmEffFactors[,i])){ rdmEffFactors[,i] <- factor(rdmEffFactors[,i]) }else{}
    }
    
    if(any(is.na(rdmEffFactors))){ stop('Error - if specifying rdmEffFactors, make sure there are no NAs!') }else{}
  }else{}

  #########################################
  ### get rid of any incomplete rows...
  #########################################
  ina <- which(is.na(z))
  ina <- unique(c(ina , which((rowSums(is.na(dfCovs)) > 0))))

  if(length(ina) > 0){
    print('Attention - some NA values in the data, getting removed for fitting LMM...')
    c <- c[-ina,,drop=FALSE]
    z <- z[-ina]
    dfCovs <- dfCovs[-ina,,drop=FALSE] 
    if(length(nAliquots) > 1){ nAliquots <- nAliquots[-ina] }else{}
    if(!is.null(rdmEffFactors)){ rdmEffFactors <- rdmEffFactors[-ina,,drop=FALSE] }else{}
  }else{}
  
  n <- length(z)

  #########################################
  ### check/preprocess data for NSFS...
  #########################################
  tmp <- preprocessdfCovs_NSFS(dfCovs = dfCovs)
  dfCovsCts <- tmp$dfCovsCts 
  dfCovsCat <- tmp$dfCovsCat 
  dfCovsCatBin <- tmp$dfCovsCatBin
  mdfCovsCtsIn <- tmp$mdfCovsCtsIn
  sddfCovsCtsIn <- tmp$sddfCovsCtsIn
  covsCat <- tmp$covsCat

  #########################################
  ### set up for fLMM2 with spatial covariance...
  #########################################
  if(nSpatStructs > 0){
    D <- xyDist(c , c)
    
    ### note parameterisation in fLMM2 fn for exp model is with a = range (not 3a = range)
    if(is.null(mina)){ mina <- min(D[D > 0]) }else{}
    if(is.null(maxa)){ maxa <- max(D) / 2 }else{}
    
    if(covModel == 'spherical' | substr(covModel , 1 , 8) == 'wendland' | substr(covModel , 1 , 9) == 'gwendland'){
      ### use sparse matrices (if v big, will have to rethink this bit)
      D[D > maxa] <- 0 # not really D=0, but cov fn will copy the 0s to use as covariances
      D <- as(D , "dsCMatrix")
    }else{}    
    
    ### check D for any really small dists...
    if(length(which((D > 0) & (D <= 1E-10))) > 0){ print('WARNING - SOME SMALL VALUES IN D COULD BE COLOCATED - CHECK DEFINITION OF D!') }else{}
  }else{}
  
  if(is.null(rdmEffFactors)){ 
    parsInit <- getParsInit_fLMM2(covModel = covModel , nSpatStructs = nSpatStructs , incNugget = incNugget , 
                                  mina = mina , maxa = maxa , fixa2 = fixa2)
  }else{
    parsInit <- getParsInit_fLMM2_Bayes(covModel = covModel , nSpatStructs = nSpatStructs , incNugget = incNugget , 
                                        mina = mina , maxa = maxa , fixa2 = fixa2 , z = z , X = matrix(1 , n , 1) , rdmEffFactors = rdmEffFactors)
  }

  ##########################################################  
  ### fit initial model with const mean...
  ##########################################################  
  if(is.null(namesCovsStart)){
    X0 <- matrix(1 , n , 1)
    colnames(X0) <- 'constant'
  }else{

    if(!is.null(dfCovsCts)){
      namesCovsCtsStart <- intersect(namesCovsStart , names(dfCovsCts))
    }else{
      namesCovsCtsStart <- c()
    }
    if(!is.null(dfCovsCat)){
      namesCovsCatStart <- intersect(namesCovsStart , names(dfCovsCat))
    }else{
      namesCovsCatStart <- c()
    }
    
    ### if we have any categorical variables requested to be (fully) in model, then use 'indicator' coding (ie no column of 1s, just columns with binary indicators)
    ### else, first column will be 1s
    if(!is.null(is.null(namesCovsCatStart))){
      X0 <- c()
      for(i in 1:length(namesCovsCatStart)){
        colsTmp <- which(is.element(names(dfCovsCatBin) , paste0(namesCovsCatStart[i] , levels(dfCovsCat[,i]))))
        X0This <- as.matrix(dfCovsCatBin[,colsTmp])
        colnames(X0This) <- names(dfCovsCatBin)[colsTmp]
        X0 <- cbind(X0 , X0This)
        min_colSums_X0This <- min(colSums(X0This))
        if(min_colSums_X0This == 0){
          stop('Error - the X0 initialised with the given namesCovsStart gives some categories with 0 cases! Remove these classes.')
        }else if(min_colSums_X0This == 1){
          stop('Error - the X0 initialised with the given namesCovsStart gives some categories with only 1 case! Merge these classes with another similar class.')
        }else{}
      }
    }else{
      X0 <- matrix(1 , n , 1)
      colnames(X0) <- 'constant'
    }    
    
    if(!is.null(is.null(namesCovsCtsStart))){
      X0This <- as.matrix(dfCovsCts[,namesCovsCtsStart])
      colnames(X0This) <- namesCovsCtsStart
      X0 <- cbind(X0 , X0This)
    }else{}
  }

  ftlmm20 <- fitLMM2(c = c , z = z , X = X0 , nAliquots = nAliquots , rdmEffFactors = rdmEffFactors ,
                 covModel = covModel , nSpatStructs = nSpatStructs , incNugget = incNugget , 
                 optionML = optionML , verbose = F , mina = mina , maxa = maxa , parsInit = parsInit , 
                 attachBigMats = T , fixa2 = fixa2)

  XCurrent <- X0 
  ftlmm2Current <- ftlmm20
  
  ################################################  
  ### lists to keep track of covariates and knots in each column of X...
  ### -9 indicates not relevant for knots
  ################################################  
  if((ncol(X0) == 1) && (colnames(X0)[1] == 'constant')){
    listCovsXCurrent <- list('constant' = 'constant')
    listqIntKnotsCurrent <- list('constant' = -9)
    
    isAddedCts <- rep(F , length(names(dfCovsCts)))
    isAddedCatBin <- rep(F , length(names(dfCovsCatBin)))
  }else{
    listCovsXCurrent <- list()
    listqIntKnotsCurrent <- list()
    for (i in 1:ncol(X0)){
      listCovsXCurrent[[colnames(X0)[i]]] <- colnames(X0)[i]
      listqIntKnotsCurrent[[colnames(X0)[i]]] <- -9
    }
    isAddedCts <- rep(F , length(names(dfCovsCts)))
    isAddedCatBin <- rep(F , length(names(dfCovsCatBin)))
    isAddedCts[is.element(names(dfCovsCts)  , colnames(X0))] <- T
    isAddedCatBin[is.element(names(isAddedCatBin)  , colnames(X0))] <- T
  }

  continueAdding <- TRUE
  while(continueAdding){
    tmp <- getAllPossColsToAdd_NSFS(XCurrent = XCurrent , listCovsXCurrent = listCovsXCurrent , listqIntKnotsCurrent = listqIntKnotsCurrent , 
                                    dfCovsCts = dfCovsCts , dfCovsCatBin = dfCovsCatBin , isAddedCts = isAddedCts , isAddedCatBin = isAddedCatBin , 
                                    maxnIntKnots = maxnIntKnots , maxnInts = maxnInts , covsCat = covsCat)
    dfXADDPoss <- tmp$dfXADDPoss 
    listCovsXADD <- tmp$listCovsXADD
    vecqKnotsXADD <- tmp$vecqKnotsXADD
    vecPoolXADD <- tmp$vecPool

    if(ncol(dfXADDPoss) > 0){
      nllVec <- delaicVec <- pvalVec <- rep(NA , ncol(dfXADDPoss))

      XADDiAXADD <- as.numeric(colSums(as.matrix(dfXADDPoss , nrow = n) * (ftlmm2Current$iA %*% as.matrix(dfXADDPoss , nrow = n))))
      XADDiAXCurrent <- t(as.matrix(dfXADDPoss)) %*% ftlmm2Current$iAX
      XADDiAz <- t(as.matrix(dfXADDPoss)) %*% ftlmm2Current$iAz

      for(j in 1:ncol(dfXADDPoss)){
        XiAXThis <- rbind(cbind(ftlmm2Current$XiAX , t(XADDiAXCurrent[j,,drop=FALSE])) , 
                          cbind(XADDiAXCurrent[j,,drop=FALSE] , XADDiAXADD[j]))
        XiAzThis <- rbind(ftlmm2Current$XiAz , XADDiAz[j,,drop=FALSE])

        ### make this function also return sigma2hat and nll...
        tmp <- fLMM2_Given_WiAW(XiAX = XiAXThis , XiAz = XiAzThis , ziAz = ftlmm2Current$ziAz , 
                                      lndetA = ftlmm2Current$lndetA , n = n , optionML = optionML , badNll = 9E99 , returnOpt = 1)
        nllVec[j] <- tmp$nll
        sigma2hatTEST <- tmp$sigma2hat
        
        delnll <- delnllFixedA_NSFS(sigma2hatF = sigma2hatTEST , sigma2hatR = ftlmm2Current$sigma2hat , 
                                    npF = ncol(XCurrent) + 1 , npR = ncol(XCurrent) , n = n , optionML = optionML)
        delaicVec[j] <- 2 * delnll + 2 * 1
      }

      if(!all(is.na(delaicVec))){
        if(prefer_fewer_covs){
          if(length(vecPoolXADD) != length(delaicVec)){ stop('Ouch - vecPool not right!') }else{}
          iP1 <- which(vecPoolXADD == 1)
          if(length(iP1) > 0){
            delaicVecTMP <- delaicVec
            delaicVecTMP[-iP1] <- 9E98
            delaicVecTMP[is.na(delaicVecTMP)|is.nan(delaicVecTMP)] <- 9E99
            if(min(delaicVecTMP) < 0){
              ### min from pool 1, if < 0
              iminaic <- which.min(delaicVecTMP)
            }else{
              # min overall...
              iminaic <- which.min(delaicVec)
            }
          }else{
            iminaic <- which.min(delaicVec)
          }
        }else{
          ### any the best according to smallest aic
          iminaic <- which.min(delaicVec)
        }

        if(verbose_NSFS){
          print(paste0('Model with ' , 
                       names(dfXADDPoss)[iminaic] , 
                       ' suggested as the best to add (change in AIC of ' , 
                       round(delaicVec[iminaic]) , 
                       '). Testing this properly...'))
        }else{}       
        ### refit that one properly?...
        if(nSpatStructs > 0){
          ### updating the maxa to make use of sparse savings...
          if(is.null(rdmEffFactors)){
            ipara2 <- length(ftlmm2Current$pars)
          }else{
            ipara2 <- length(ftlmm2Current$pars) - ncol(rdmEffFactors)
          } 
          maxa <- ftlmm2Current$pars[ipara2]
        }else{}
        
        ftlmm2New <- fitLMM2(c = c , z = z , X = cbind(XCurrent , dfXADDPoss[,iminaic]) , nAliquots = nAliquots , rdmEffFactors = rdmEffFactors , 
                             covModel = covModel , nSpatStructs = nSpatStructs , incNugget = incNugget , 
                             optionML = optionML , verbose = F , mina = mina , maxa = maxa , parsInit = ftlmm2Current$pars , 
                             attachBigMats = T , fixa2 = fixa2)

        if(!is.null(ftlmm2New$XiAX)){
          lndetXiAXNew <- determinant(ftlmm2New$XiAX , logarithm = TRUE)
          if(lndetXiAXNew$sign < 0){ 
            lndetXiAXNew <- NA
          }else{
            lndetXiAXNew <- as.numeric(lndetXiAXNew$modulus)
          }
          
          XFiACurrentXF <- rbind(cbind(ftlmm2Current$XiAX , t(XADDiAXCurrent[iminaic,,drop=FALSE])) , 
                                 cbind(XADDiAXCurrent[iminaic,,drop=FALSE] , XADDiAXADD[iminaic]))
          lndetXFiACurrentXF <- determinant(XFiACurrentXF , logarithm = TRUE)
          if(lndetXFiACurrentXF$sign < 0){ 
            lndetXFiACurrentXF <- NA
          }else{
            lndetXFiACurrentXF <- as.numeric(lndetXFiACurrentXF$modulus)
          }

          ### change in nll with different params in reduced and full models.  
          delnllTest <- delnll_NSFS(sigma2hatF = ftlmm2New$sigma2hat , sigma2hatR = ftlmm2Current$sigma2hat , 
                                    lndetAF  = ftlmm2New$lndetA , lndetAR = ftlmm2Current$lndetA , 
                                    lndetXFiAFXF = lndetXiAXNew , lndetXFiARXF = lndetXFiACurrentXF , 
                                    npF = ncol(XCurrent) + 1 , npR = ncol(XCurrent) , n = n , optionML = optionML)
          delaicTest <- 2 * delnllTest + 2 * 1
          
          if(is.na(delaicTest) | is.nan(delaicTest)){
            lndetXiAXNew <- NA
            delnllTest <- Inf
            delaicTest <- Inf
          }else{}

        }else{
          lndetXiAXNew <- NA
          delnllTest <- Inf
          delaicTest <- Inf
        }
        
        if(verbose_NSFS){
          txtTmp <- paste0('When correlation structure for the full model was refitted, change in AIC was ' , 
                           round(delaicTest))
          if(delaicTest >= 0){
            txtTmp <- paste0(txtTmp , '. So stopping now.')
          }else{
            txtTmp <- paste0(txtTmp , '. So adding this to model...')
          }
          print(txtTmp)
        }else{}
        
      }else{
        # so that will stop now.
        delaicTest <- Inf
        if(verbose_NSFS){
          print('All of the additions give some colinearity, so stopping.')
        }else{}
      }
      
      ### any decreases in aic?
      if(delaicTest < 0){
        newNames <- c(colnames(XCurrent) , names(dfXADDPoss)[iminaic])
        XCurrent <- cbind(XCurrent , dfXADDPoss[,iminaic]) 
        colnames(XCurrent) <- newNames
        
        ### add to lists of covs/knots
        listCovsXCurrent[[ncol(XCurrent)]] <- listCovsXADD[[iminaic]]
        if(vecqKnotsXADD[iminaic] == -9){
          listqIntKnotsCurrent[[ncol(XCurrent)]] <- vecqKnotsXADD[iminaic]
        }else{
          qTmp <- colnames(XCurrent)[ncol(XCurrent)]
          if(!grepl('_NS_' , qTmp)){ stop('Error - should have found _NS_ in new colname!') }else{}
          qTmp <- as.numeric(strsplit(qTmp , '_NS_')[[1]][-1])
          listqIntKnotsCurrent[[ncol(XCurrent)]] <- qTmp
        }

        ### update the 'Current' model...
        ftlmm2Current <- ftlmm2New
        
        if(length(isAddedCts) > 0){
          isAddedCts[is.element(names(dfCovsCts) , listCovsXADD[[iminaic]])] <- TRUE
        }else{}
        if(length(isAddedCatBin) > 0){
          isAddedCatBin[is.element(names(dfCovsCatBin) , listCovsXADD[[iminaic]])] <- TRUE
        }else{}

        if(verbose_NSFS){
          print('Added to give model with:')
          print(colnames(XCurrent))
        }else{}
      }else{
        continueAdding <- FALSE
      }
      
    }else{
      continueAdding <- FALSE
    }
  }

  if(verbose_NSFS){
    print('All done! Final model has:')
    print(colnames(XCurrent))
  }else{}
  
  return(list('ftlmm2Final' = ftlmm2Current , 'XFinal' = XCurrent , 
              'listCovsXFinal' = listCovsXCurrent , 
              'listqIntKnotsFinal' = listqIntKnotsCurrent))
}

# ##################################################################
# ### test version where cat variables are added fully or interacted fully
# ### instead of adding a single binary column or its interation.
# ### NOT DONE!
# ##################################################################
# modelSelectLMM2_NSFS_TEST <- function(c , z , dfCovs = NULL , namesCovsStart = NULL , nAliquots = 1 , covModel = 'exponential' , nSpatStructs = 1 , incNugget = T , optionML = F , verbose_NSFS = FALSE , mina = NULL , maxa = NULL , parsInit = NULL , maxnIntKnots = 3 , maxnInts = 2){
#   
#   fixa2 <- F
#   
#   if(is.null(dfCovs)){ stop('Error - must include dfCovs in call to modelSelectLMM2_NSFS!') }else{}
#   if(!is.data.frame(dfCovs)){ stop('Error - dfCovs must be a data.frame in call to modelSelectLMM2_NSFS!') }else{}
#   if((!is.null(namesCovsStart)) && (length(namesCovsStart) > 0) && (!is.element(namesCovsStart , names(dfCovs)))){ stop('Error - namesCovsStart must be in dfCovs!') }else{}
#   
#   #########################################
#   ### nAliquots used as if these cores were sampled a negligible distance from each other
#   ### so diagonal of cov mtx for entries where nAliquots > 1 = c1 + c2 + c0/nAliquots
#   #########################################
#   if(length(nAliquots) > 1){ if(length(nAliquots) != length(z)){ stop('Error - nAliquots must be same length as z for fitLMM2!') }else{} }else{}
# 
#   #########################################
#   ### get rid of any incomplete rows...
#   #########################################
#   ina <- which(is.na(z))
#   ina <- unique(c(ina , which((rowSums(is.na(dfCovs)) > 0))))
#   
#   if(length(ina) > 0){
#     print('Attention - some NA values in the data, getting removed for fitting LMM...')
#     c <- c[-ina,,drop=FALSE]
#     z <- z[-ina]
#     dfCovs <- dfCovs[-ina,,drop=FALSE] 
#     if(length(nAliquots) > 1){ nAliquots <- nAliquots[-ina] }else{}
#   }else{}
#   
#   n <- length(z)
#   
#   #########################################
#   ### check/preprocess data for NSFS...
#   #########################################
#   tmp <- preprocessdfCovs_NSFS(dfCovs = dfCovs)
#   dfCovsCts <- tmp$dfCovsCts 
#   dfCovsCat <- tmp$dfCovsCat 
#   dfCovsCatBin <- tmp$dfCovsCatBin
#   mdfCovsCtsIn <- tmp$mdfCovsCtsIn
#   sddfCovsCtsIn <- tmp$sddfCovsCtsIn
#   covsCat <- tmp$covsCat
#   
#   #########################################
#   ### set up for fLMM2 with spatial covariance...
#   #########################################
#   if(nSpatStructs > 0){
#     D <- xyDist(c , c)
#     
#     ### note parameterisation in fLMM2 fn for exp model is with a = range (not 3a = range)
#     if(is.null(mina)){ mina <- min(D[D > 0]) }else{}
#     if(is.null(maxa)){ maxa <- max(D) / 2 }else{}
#     
#     if(covModel == 'spherical' | substr(covModel , 1 , 8) == 'wendland' | substr(covModel , 1 , 9) == 'gwendland'){
#       ### use sparse matrices (if v big, will have to rethink this bit)
#       D[D > maxa] <- 0 # not really D=0, but cov fn will copy the 0s to use as covariances
#       D <- as(D , "dsCMatrix")
#     }else{}    
#     
#     ### check D for any really small dists...
#     if(length(which((D > 0) & (D <= 1E-10))) > 0){ print('WARNING - SOME SMALL VALUES IN D COULD BE COLOCATED - CHECK DEFINITION OF D!') }else{}
#   }else{}
#   
#   parsInit <- getParsInit_fLMM2(covModel = covModel , nSpatStructs = nSpatStructs , incNugget = incNugget , 
#                                 mina = mina , maxa = maxa , fixa2 = fixa2)
# 
#   ##########################################################  
#   ### fit initial model with const mean...
#   ##########################################################  
#   if(is.null(namesCovsStart)){
#     X0 <- matrix(1 , n , 1)
#     colnames(X0) <- 'constant'
#   }else{
#     
#     if(!is.null(dfCovsCts)){
#       namesCovsCtsStart <- intersect(namesCovsStart , names(dfCovsCts))
#     }else{
#       namesCovsCtsStart <- c()
#     }
#     if(!is.null(dfCovsCat)){
#       namesCovsCatStart <- intersect(namesCovsStart , names(dfCovsCat))
#     }else{
#       namesCovsCatStart <- c()
#     }
#     
#     ### if we have any categorical variables requested to be (fully) in model, then use 'indicator' coding (ie no column of 1s, just columns with binary indicators)
#     ### else, first column will be 1s
#     if(!is.null(is.null(namesCovsCatStart))){
#       X0 <- c()
#       for(i in 1:length(namesCovsCatStart)){
#         colsTmp <- which(is.element(names(dfCovsCatBin) , paste0(namesCovsCatStart[i] , levels(dfCovsCat[,i]))))
#         X0This <- as.matrix(dfCovsCatBin[,colsTmp])
#         colnames(X0This) <- names(dfCovsCatBin)[colsTmp]
#         X0 <- cbind(X0 , X0This)
#         min_colSums_X0This <- min(colSums(X0This))
#         if(min_colSums_X0This == 0){
#           stop('Error - the X0 initialised with the given namesCovsStart gives some categories with 0 cases! Remove these classes.')
#         }else if(min_colSums_X0This == 1){
#           stop('Error - the X0 initialised with the given namesCovsStart gives some categories with only 1 case! Merge these classes with another similar class.')
#         }else{}
#       }
#     }else{
#       X0 <- matrix(1 , n , 1)
#       colnames(X0) <- 'constant'
#     }    
#     
#     if(!is.null(is.null(namesCovsCtsStart))){
#       X0This <- as.matrix(dfCovsCts[,namesCovsCtsStart])
#       colnames(X0This) <- namesCovsCtsStart
#       X0 <- cbind(X0 , X0This)
#     }else{}
#   }
#   
#   ftlmm20 <- fitLMM2(c = c , z = z , X = X0 , nAliquots = nAliquots , 
#                      covModel = covModel , nSpatStructs = nSpatStructs , incNugget = incNugget , 
#                      optionML = optionML , verbose = F , mina = mina , maxa = maxa , parsInit = parsInit , 
#                      attachBigMats = T , fixa2 = fixa2)
#   
#   XCurrent <- X0 
#   ftlmm2Current <- ftlmm20
#   
#   ################################################  
#   ### lists to keep track of covariates and knots in each column of X...
#   ### -9 indicates not relevant for knots
#   ################################################  
#   if((ncol(X0) == 1) && (colnames(X0)[1] == 'constant')){
#     listCovsXCurrent <- list('constant' = 'constant')
#     listqIntKnotsCurrent <- list('constant' = -9)
#     
#     isAddedCts <- rep(F , length(names(dfCovsCts)))
#     isAddedCatBin <- rep(F , length(names(dfCovsCatBin)))
#   }else{
#     listCovsXCurrent <- list()
#     listqIntKnotsCurrent <- list()
#     for (i in 1:ncol(X0)){
#       listCovsXCurrent[[colnames(X0)[i]]] <- colnames(X0)[i]
#       listqIntKnotsCurrent[[colnames(X0)[i]]] <- -9
#     }
#     isAddedCts <- rep(F , length(names(dfCovsCts)))
#     isAddedCatBin <- rep(F , length(names(dfCovsCatBin)))
#     isAddedCts[is.element(names(dfCovsCts)  , colnames(X0))] <- T
#     isAddedCatBin[is.element(names(isAddedCatBin)  , colnames(X0))] <- T
#   }
#   
#   
#   continueAdding <- TRUE
#   while(continueAdding){
#     tmp <- getAllPossColsToAdd_NSFS(XCurrent = XCurrent , listCovsXCurrent = listCovsXCurrent , listqIntKnotsCurrent = listqIntKnotsCurrent , 
#                                     dfCovsCts = dfCovsCts , dfCovsCatBin = dfCovsCatBin , isAddedCts = isAddedCts , isAddedCatBin = isAddedCatBin , 
#                                     maxnIntKnots = maxnIntKnots , maxnInts = maxnInts , covsCat = covsCat)
#     dfXADDPoss <- tmp$dfXADDPoss 
#     listCovsXADD <- tmp$listCovsXADD
#     vecqKnotsXADD <- tmp$vecqKnotsXADD
#     
#     if(verbose_NSFS){
#       print('We will test the following:')
#       print(names(dfXADDPoss))
#     }else{}
#     
#     if(ncol(dfXADDPoss) > 0){
#       nllVec <- delaicVec <- rep(NA , ncol(dfXADDPoss))
#       
#       XADDiAXADD <- as.numeric(colSums(as.matrix(dfXADDPoss , nrow = n) * (ftlmm2Current$iA %*% as.matrix(dfXADDPoss , nrow = n))))
#       XADDiAXCurrent <- t(as.matrix(dfXADDPoss)) %*% ftlmm2Current$iAX
#       XADDiAz <- t(as.matrix(dfXADDPoss)) %*% ftlmm2Current$iAz
#       
#       for(j in 1:ncol(dfXADDPoss)){
#         XiAXThis <- rbind(cbind(ftlmm2Current$XiAX , t(XADDiAXCurrent[j,,drop=FALSE])) , 
#                           cbind(XADDiAXCurrent[j,,drop=FALSE] , XADDiAXADD[j]))
#         XiAzThis <- rbind(ftlmm2Current$XiAz , XADDiAz[j,,drop=FALSE])
#         
#         ### make this function also return sigma2hat and nll...
#         tmp <- fLMM2_Given_WiAW(XiAX = XiAXThis , XiAz = XiAzThis , ziAz = ftlmm2Current$ziAz , 
#                                 lndetA = ftlmm2Current$lndetA , n = n , optionML = optionML , badNll = 9E99 , returnOpt = 1)
#         nllVec[j] <- tmp$nll
#         sigma2hatTEST <- tmp$sigma2hat
#         
#         delnll <- delnllFixedA_NSFS(sigma2hatF = sigma2hatTEST , sigma2hatR = ftlmm2Current$sigma2hat , 
#                                     npF = ncol(XCurrent) + 1 , npR = ncol(XCurrent) , n = n , optionML = optionML)
#         delaicVec[j] <- 2 * delnll + 2 * 1
#       }
#       
#       if(!all(is.na(delaicVec))){
#         ### any the best according to smallest aic
#         iminaic <- which.min(delaicVec)
#         
#         if(verbose_NSFS){
#           print(paste0('Model with ' , 
#                        names(dfXADDPoss)[iminaic] , 
#                        ' suggested as the best to add (change in AIC of ' , 
#                        round(delaicVec[iminaic]) , 
#                        '). Testing this properly...'))
#         }else{}       
#         ### refit that one properly?...
#         if(nSpatStructs > 0){
#           ### updating the maxa to make use of sparse savings...
#           maxa <- ftlmm2Current$pars[length(ftlmm2Current$pars)]
#         }else{}
#         ftlmm2New <- fitLMM2(c = c , z = z , X = cbind(XCurrent , dfXADDPoss[,iminaic]) , nAliquots = nAliquots , 
#                              covModel = covModel , nSpatStructs = nSpatStructs , incNugget = incNugget , 
#                              optionML = optionML , verbose = F , mina = mina , maxa = maxa , parsInit = ftlmm2Current$pars , 
#                              attachBigMats = T , fixa2 = fixa2)
#         
#         if(!is.null(ftlmm2New$XiAX)){
#           lndetXiAXNew <- determinant(ftlmm2New$XiAX , logarithm = TRUE)
#           if(lndetXiAXNew$sign < 0){ 
#             lndetXiAXNew <- NA
#           }else{
#             lndetXiAXNew <- as.numeric(lndetXiAXNew$modulus)
#           }
#           
#           XFiACurrentXF <- rbind(cbind(ftlmm2Current$XiAX , t(XADDiAXCurrent[iminaic,,drop=FALSE])) , 
#                                  cbind(XADDiAXCurrent[iminaic,,drop=FALSE] , XADDiAXADD[iminaic]))
#           lndetXFiACurrentXF <- determinant(XFiACurrentXF , logarithm = TRUE)
#           if(lndetXFiACurrentXF$sign < 0){ 
#             lndetXFiACurrentXF <- NA
#           }else{
#             lndetXFiACurrentXF <- as.numeric(lndetXFiACurrentXF$modulus)
#           }
#           
#           ### change in nll with different params in reduced and full models.  
#           delnllTest <- delnll_NSFS(sigma2hatF = ftlmm2New$sigma2hat , sigma2hatR = ftlmm2Current$sigma2hat , 
#                                     lndetAF  = ftlmm2New$lndetA , lndetAR = ftlmm2Current$lndetA , 
#                                     lndetXFiAFXF = lndetXiAXNew , lndetXFiARXF = lndetXFiACurrentXF , 
#                                     npF = ncol(XCurrent) + 1 , npR = ncol(XCurrent) , n = n , optionML = optionML)
#           delaicTest <- 2 * delnllTest + 2 * 1
#           
#           if(is.na(delaicTest) | is.nan(delaicTest)){
#             lndetXiAXNew <- NA
#             delnllTest <- Inf
#             delaicTest <- Inf
#           }else{}
#           
#         }else{
#           lndetXiAXNew <- NA
#           delnllTest <- Inf
#           delaicTest <- Inf
#         }
#         
#         if(verbose_NSFS){
#           txtTmp <- paste0('When correlation structure for the full model was refitted, change in AIC was ' , 
#                            round(delaicTest))
#           if(delaicTest >= 0){
#             txtTmp <- paste0(txtTmp , '. So stopping now.')
#           }else{
#             txtTmp <- paste0(txtTmp , '. So adding this to model...')
#           }
#           print(txtTmp)
#         }else{}
#         
#       }else{
#         # so that will stop now.
#         delaicTest <- Inf
#         if(verbose_NSFS){
#           print('All of the additions give some colinearity, so stopping.')
#         }else{}
#       }
#       
#       ### any decreases in aic?
#       if(delaicTest < 0){
#         newNames <- c(colnames(XCurrent) , names(dfXADDPoss)[iminaic])
#         XCurrent <- cbind(XCurrent , dfXADDPoss[,iminaic]) 
#         colnames(XCurrent) <- newNames
#         
#         ### add to lists of covs/knots
#         listCovsXCurrent[[ncol(XCurrent)]] <- listCovsXADD[[iminaic]]
#         if(vecqKnotsXADD[iminaic] == -9){
#           listqIntKnotsCurrent[[ncol(XCurrent)]] <- vecqKnotsXADD[iminaic]
#         }else{
#           qTmp <- colnames(XCurrent)[ncol(XCurrent)]
#           if(!grepl('_NS_' , qTmp)){ stop('Error - should have found _NS_ in new colname!') }else{}
#           qTmp <- as.numeric(strsplit(qTmp , '_NS_')[[1]][-1])
#           listqIntKnotsCurrent[[ncol(XCurrent)]] <- qTmp
#         }
#         
#         ### update the 'Current' model...
#         ftlmm2Current <- ftlmm2New
#         
#         if(length(isAddedCts) > 0){
#           isAddedCts[is.element(names(dfCovsCts) , listCovsXADD[[iminaic]])] <- TRUE
#         }else{}
#         if(length(isAddedCatBin) > 0){
#           isAddedCatBin[is.element(names(dfCovsCatBin) , listCovsXADD[[iminaic]])] <- TRUE
#         }else{}
#         
#         if(verbose_NSFS){
#           print('Added to give model with:')
#           print(colnames(XCurrent))
#         }else{}
#       }else{
#         continueAdding <- FALSE
#       }
#       
#     }else{
#       continueAdding <- FALSE
#     }
#   }
#   
#   if(verbose_NSFS){
#     print('All done! Final model has:')
#     print(colnames(XCurrent))
#   }else{}
#   
#   return(list('ftlmm2Final' = ftlmm2Current , 'XFinal' = XCurrent , 
#               'listCovsXFinal' = listCovsXCurrent , 
#               'listqIntKnotsFinal' = listqIntKnotsCurrent))
# }
# 
##################################################################################################
### calculate change in nll (nllFull - nllRed) from Full and nested models
### assume corr mtx unchanged, A.
### if model fitted by reml : Welham and Thomson method for comparing reml vals.
### ie for computing reml of reduced model, use full X in lndet(X iCR X)
### so fit represents a constrained fit of full model, with betai =  0 constrained. (or any lin constraints?)
##################################################################################################
delnllFixedA_NSFS <- function(sigma2hatF , sigma2hatR , npF , npR , n , optionML = FALSE){
  
  if((!is.na(sigma2hatF)) && (!is.na(sigma2hatR)) && (sigma2hatF > 0) && (sigma2hatR > 0)){
    if(optionML){
      delnll <- 0.5 * n * (log(sigma2hatF) - log(sigma2hatR))
    }else{
      delnll <- 0.5 * ((n - npF) * (log(sigma2hatF) - log(sigma2hatR)) + (npR - npF))
    }
  }else{
    delnll <- NA
  }
  
  return(delnll)
}

##################################################################################################
### calculate change in nll (nllFull - nllRed) from Full and nested models
### corr mtx, A, refitted.
### if model fitted by reml : Welham and Thomson method for comparing reml vals.
### ie for computing reml of reduced model, use full X in lndet(X iCR X)
### so fit represents a constrained fit of full model, with betai =  0 constrained. (or any lin constraints?)
##################################################################################################
delnll_NSFS <- function(sigma2hatF , sigma2hatR , lndetAF , lndetAR , lndetXFiAFXF , lndetXFiARXF , npF , npR , n , optionML = FALSE){

  if((!is.na(sigma2hatF)) && (!is.na(sigma2hatR)) && (sigma2hatF > 0) && (sigma2hatR > 0)){
    if(optionML){
      delnll <- 0.5 * (n * (log(sigma2hatF) - log(sigma2hatR)) + lndetAF - lndetAR)
    }else{
      delnll <- 0.5 * ((n - npF) * (log(sigma2hatF) - log(sigma2hatR)) + 
                          lndetAF - lndetAR + lndetXFiAFXF - lndetXFiARXF +
                          (npR - npF))
    }
  }else{
    delnll <- NA
  }  
  return(delnll)
}

##########################################################
### get a df with all possible cols to add to XCurrent
### either lin fn of cts or cat-bin
### or extra knot for cts var
### or multiply two cols
###
### columns for binary of factor are labelled covlevel (ie simple append)
### columns for linear of cts variable are labelled cov
### columns for ns of cts variable are labelled cov_NS_0.5,cov_NS_0.5_0.75,cov_NS_0.5_0.75_0.625, etc...
### columns for interaction of two columns are made by joining with '___'
###
### also making to return a vec with 'pool' = 1 if in the first pool, 2 if in the second
### pool 1 for just using variables already in model (ie interact/spline)
### pool 2 for adding in a new variable
### other cats from a factor variable also put in pool 1.
### covsCat is the vector with the names of the original categorical covariates,
### so can put other non-added cats in pool 1.
##########################################################
getAllPossColsToAdd_NSFS <- function(XCurrent , listCovsXCurrent , listqIntKnotsCurrent , dfCovsCts , dfCovsCatBin , isAddedCts , isAddedCatBin , maxnIntKnots , maxnInts , covsCat = NULL){
  
  dfXADDPoss <- data.frame('junk' = rep(NA , nrow(XCurrent)) , stringsAsFactors = F)[,c(),drop=FALSE]
  vecPool <- c()
  
  ### new covariates as lin fns...
  if((!is.null(dfCovsCts)) && (ncol(dfCovsCts) > 0) && (sum(isAddedCts) < ncol(dfCovsCts))){
    dfXADDPoss <- cbind(dfXADDPoss , dfCovsCts[,!isAddedCts,drop=FALSE])
    vecPool <- c(vecPool , rep(2 , sum(!isAddedCts)))
  }else{}

  ### new cat binary covariates as lin fns...
  if((!is.null(dfCovsCatBin)) && (ncol(dfCovsCatBin) > 0) && (sum(isAddedCatBin) < ncol(dfCovsCatBin))){
    dfXADDPoss <- cbind(dfXADDPoss , dfCovsCatBin[,!isAddedCatBin,drop=FALSE])
    catVarsAdded <- unique(covsCat[which(isAddedCatBin)])
    
    jTmp <- which(!isAddedCatBin)
    for(j in jTmp){
      if(is.element(covsCat[j] , catVarsAdded)){
        vecPool <- c(vecPool , 1)
      }else{
        vecPool <- c(vecPool , 2)
      }
    }
  }else{}

  ### dummy values for knots for the linear fns...  
  vecqKnotsXADD <- rep(-9 , ncol(dfXADDPoss))
  listCovsXADD <- as.list(names(dfXADDPoss))

  ### columns to complexify ns fns... 
  ### remember to also check not too many knots...
  cols4NS <- which((colnames(XCurrent) != 'constant') & 
                     (!is.element(colnames(XCurrent) , names(dfCovsCatBin))) & 
                     (!grepl('___' , colnames(XCurrent))))  
  if((!is.null(dfCovsCts)) && (ncol(dfCovsCts) > 0) && (length(cols4NS) > 0)){
    for(i in cols4NS){
      covNameThis <- listCovsXCurrent[[i]]
      covValsThis <- dfCovsCts[,covNameThis]
      ### number of columns with _NS_ and no ___ (interaction) is number of internal knots.
      countColsTmp <- length(which(listCovsXCurrent == covNameThis)) - 1 # -1 as one is the lin term
      if(countColsTmp < maxnIntKnots){
        dfTmp <- getPossKnotSegmentsToAdd_NSFS(qKnotsIn = listqIntKnotsCurrent[[i]] , maxnIntKnots = maxnIntKnots)
        if(nrow(dfTmp) > 0){
          vecqKnotsXADD <- c(vecqKnotsXADD , dfTmp$qintKnots)
          for(j in 1:nrow(dfTmp)){
            vecPool <- c(vecPool , 1)
            nameTmp <- paste0(colnames(XCurrent)[i] , '_NS_' , dfTmp$qintKnots[j])
            dfXADDPoss[,nameTmp] <- makeXns(covValsThis , 
                                            bdryKnots = quantile(covValsThis , c(dfTmp$qbdryKnotsL[j] , dfTmp$qbdryKnotsU[j])) , 
                                            intKnots = quantile(covValsThis , dfTmp$qintKnots[j]))[,3]
            listCovsXADD[[length(listCovsXADD) + 1]] <- covNameThis
          }
        }else{}
      }else{}
    }
  }else{}

  ### or multiply together two existing columns...
  ### as long as they are not based on the same coavriate
  if(ncol(XCurrent) > 1){
    for(i in 1:ncol(XCurrent)){
      for(j in 1:ncol(XCurrent)){
        if((length(intersect(listCovsXCurrent[[i]] , listCovsXCurrent[[j]])) == 0) & 
           (length(c(listCovsXCurrent[[i]] , listCovsXCurrent[[j]])) <= maxnInts) &
           (colnames(XCurrent)[i] != 'constant') & (colnames(XCurrent)[j] != 'constant')){
          vecPool <- c(vecPool , 1)
          dfXADDPoss[,paste0(colnames(XCurrent)[i] , '___' , colnames(XCurrent)[j])] <- XCurrent[,i] * XCurrent[,j]
          listCovsXADD[[length(listCovsXADD) + 1]] <- c(listCovsXCurrent[[i]] , listCovsXCurrent[[j]])
          vecqKnotsXADD <- c(vecqKnotsXADD , -9)
        }else{}
      }
    }
  }else{}

  ### get rid of anything already in model...
  ### and anything with too many interactions (not sure why these seem to have got through?)...
  ### counting interactions as the number of '___' + 1 ; so eg, a___b has nInts = 2, a has nInts = 1
  iDel <- which(is.element(names(dfXADDPoss) , colnames(XCurrent)))
  intCount <- unlist(lapply(lapply(gregexpr('___' , names(dfXADDPoss)) , 'as.numeric') , length)) + 1
  i1 <- which(unlist(lapply(gregexpr('___' , names(dfXADDPoss)) , '[' , 1) ) == -1)
  intCount[i1] <- 1
  iDel <- c(iDel , which(intCount > maxnInts))
  if(length(iDel) > 0){
    vecPool <- vecPool[-iDel]
    dfXADDPoss <- dfXADDPoss[,-iDel,drop=FALSE]
    listCovsXADD <- listCovsXADD[-iDel]
    vecqKnotsXADD <- vecqKnotsXADD[-iDel]
  }else{}
  
  return(list('dfXADDPoss' = dfXADDPoss , 'listCovsXADD' = listCovsXADD , 'vecqKnotsXADD' = vecqKnotsXADD , 'vecPool' = vecPool))  
}

# ############################################################
# ### test version using complete categorical variables, not their binaries
# ### NOT DONE!!!
# ############################################################
# getAllPossColsToAdd_NSFS_TEST <- function(XCurrent , listCovsXCurrent , listqIntKnotsCurrent , dfCovsCts , list_dfCovsCat , isAddedCts , isAddedCat , maxnIntKnots , maxnInts){
#   
#   ### all names of the purely categorical columns...
#   allNamesCovsCatBin <- c()
#   if(length(list_dfCovsCat) > 0){
#     for(i in 1:length(list_dfCovsCat)){
#       allNamesCovsCatBin <- c(allNamesCovsCatBin , names(list_dfCovsCat[[i]]))
#     }
#   }else{}
#   
#   # dfXADDPoss <- data.frame('junk' = rep(NA , nrow(XCurrent)) , stringsAsFactors = F)[,c(),drop=FALSE]
#   listXADDPoss <- list()
#   # namesXADDPoss <- c()
#   npXADDPoss <- c()
#   
#   ### new covariates as lin fns...
#   if((!is.null(dfCovsCts)) && (ncol(dfCovsCts) > 0) && (sum(isAddedCts) < ncol(dfCovsCts))){
#     for(i in 1:length(isAddedCts)){
#       if(!isAddedCts[i]){
#         listXADDPoss[[names(dfCovsCts)[i]]] <- dfCovsCts[,i,drop=FALSE]
#         # namesXADDPoss <- c(namesXADDPoss , names(dfCovsCts)[i])
#         npXADDPoss <- c(npXADDPoss , 1)
#       }else{}
#     }
#   }else{}
#   
#   ### new cat covariates as lin fns...
#   # if((!is.null(dfCovsCatBin)) && (ncol(dfCovsCatBin) > 0) && (sum(isAddedCatBin) < ncol(dfCovsCatBin))){
#   if((!is.null(list_dfCovsCat)) && (length(list_dfCovsCat) > 0) && (sum(isAddedCat) < length(list_dfCovsCat))){
#     for(i in 1:length(isAddedCat)){
#       if(!isAddedCat[i]){
#         listXADDPoss[[names(list_dfCovsCat)[i]]] <- list_dfCovsCat[[i]]
#         # namesXADDPoss <- c(namesXADDPoss , names(list_dfCovsCat)[i])
#         npXADDPoss <- c(npXADDPoss , ncol(list_dfCovsCat[[i]]))
#       }else{}
#     }
#   }else{}
#   
#   ### dummy values for knots for the linear fns...  
#   vecqKnotsXADD <- rep(-9 , length(listXADDPoss))
#   listCovsXADD <- as.list(names(listXADDPoss))
#   
#   ### columns to complexify ns fns... 
#   ### remember to also check not too many knots...
#   cols4NS <- which((colnames(XCurrent) != 'constant') & 
#                      (!is.element(colnames(XCurrent) , allNamesCovsCatBin)) & 
#                      (!grepl('___' , colnames(XCurrent))))  
#   if((!is.null(dfCovsCts)) && (ncol(dfCovsCts) > 0) && (length(cols4NS) > 0)){
#     for(i in cols4NS){
#       covNameThis <- listCovsXCurrent[[i]]
#       covValsThis <- dfCovsCts[,covNameThis]
#       ### number of columns with _NS_ and no ___ (interaction) is number of internal knots.
#       countColsTmp <- length(which(listCovsXCurrent == covNameThis)) - 1 # -1 as one is the lin term
#       if(countColsTmp < maxnIntKnots){
#         dfTmp <- getPossKnotSegmentsToAdd_NSFS(qKnotsIn = listqIntKnotsCurrent[[i]] , maxnIntKnots = maxnIntKnots)
#         if(nrow(dfTmp) > 0){
#           vecqKnotsXADD <- c(vecqKnotsXADD , dfTmp$qintKnots)
#           for(j in 1:nrow(dfTmp)){
#             nameTmp <- paste0(colnames(XCurrent)[i] , '_NS_' , dfTmp$qintKnots[j])
#             listXADDPoss[[nameTmp]] <- makeXns(covValsThis , 
#                                             bdryKnots = quantile(covValsThis , c(dfTmp$qbdryKnotsL[j] , dfTmp$qbdryKnotsU[j])) , 
#                                             intKnots = quantile(covValsThis , dfTmp$qintKnots[j]))[,3]
#             npXADDPoss <- c(npXADDPoss , 1)
#             listCovsXADD[[length(listCovsXADD) + 1]] <- covNameThis
#           }
#         }else{}
#       }else{}
#     }
#   }else{}
#   
#   ### or multiply together two existing columns...
#   ### as long as they are not based on the same coavriate
#   if(ncol(XCurrent) > 1){
#     for(i in 1:ncol(XCurrent)){
#       for(j in 1:ncol(XCurrent)){
#         if((length(intersect(listCovsXCurrent[[i]] , listCovsXCurrent[[j]])) == 0) & 
#            (length(c(listCovsXCurrent[[i]] , listCovsXCurrent[[j]])) <= maxnInts) &
#            (colnames(XCurrent)[i] != 'constant') & (colnames(XCurrent)[j] != 'constant')){
#           dfXADDPoss[,paste0(colnames(XCurrent)[i] , '___' , colnames(XCurrent)[j])] <- XCurrent[,i] * XCurrent[,j]
#           listCovsXADD[[length(listCovsXADD) + 1]] <- c(listCovsXCurrent[[i]] , listCovsXCurrent[[j]])
#           vecqKnotsXADD <- c(vecqKnotsXADD , -9)
#         }else{}
#       }
#     }
#   }else{}
#   
#   ### get rid of anything already in model...
#   ### and anything with too many interactions (not sure why these seem to have got through?)...
#   ### counting interactions as the number of '___' + 1 ; so eg, a___b has nInts = 2, a has nInts = 1
#   iDel <- which(is.element(names(dfXADDPoss) , colnames(XCurrent)))
#   intCount <- unlist(lapply(lapply(gregexpr('___' , names(dfXADDPoss)) , 'as.numeric') , length)) + 1
#   i1 <- which(unlist(lapply(gregexpr('___' , names(dfXADDPoss)) , '[' , 1) ) == -1)
#   intCount[i1] <- 1
#   iDel <- c(iDel , which(intCount > maxnInts))
#   if(length(iDel) > 0){
#     dfXADDPoss <- dfXADDPoss[,-iDel,drop=FALSE]
#     listCovsXADD <- listCovsXADD[-iDel]
#     vecqKnotsXADD <- vecqKnotsXADD[-iDel]
#   }else{}
#   
#   return(list('dfXADDPoss' = dfXADDPoss , 'listCovsXADD' = listCovsXADD , 'vecqKnotsXADD' = vecqKnotsXADD))  
# }
# 
##################################################
### currently model has knots at qKnotsIn
### suggest possible additional segments by splitting these with natural spline 
##################################################
getPossKnotSegmentsToAdd_NSFS <- function(qKnotsIn , maxnIntKnots = Inf){
  if(length(qKnotsIn) >= maxnIntKnots){
    return(data.frame('qbdryKnotsL' =  NA , 'qbdryKnotsU' = NA , 'qintKnots' = NA)[c(),,drop=FALSE]) 
  }else{}
  
  if(length(qKnotsIn) > 1){
    qKnotsIn <- qKnotsIn[order(qKnotsIn)]
  }else{}
  
  if((length(qKnotsIn) == 1) && (qKnotsIn == -9)){
    dfqKnotsOut <- data.frame('qbdryKnotsL' =  0 , 
                              'qbdryKnotsU' = 1)
    dfqKnotsOut$qintKnots <- 0.5 
  }else{
    dfqKnotsOut <- data.frame('qbdryKnotsL' =  c(0 , qKnotsIn) , 
                              'qbdryKnotsU' = c(qKnotsIn , 1))
    dfqKnotsOut$qintKnots <- 0.5 * (dfqKnotsOut$qbdryKnotsL + dfqKnotsOut$qbdryKnotsU)
  }
  
  return(dfqKnotsOut)
}

##################################################
### function to make X based on the names colnamesX_NSFS
### dfCovsCts and dfCovsCat MUST be the data passed to modelSelectLMM2_NSFS
### if dfCovsCts_Pred and dfCovsCat_Pred are null, then X for the calibration data (dfCovsCts , dfCovsCat) is returned
### if dfCovsCts_Pred and dfCovsCat_Pred are given, then X for the prediction data (dfCovsCts_Pred , dfCovsCat_Pred) is returned
##################################################
makeX_NSFS <- function(colnamesX_NSFS , z = NULL , dfCovs = NULL , dfCovs_Pred = NULL){
  
  if(is.null(z)){
    stop('Error - must enter z (the original calibration response data) to makeX_NSFS, just in case some rows were removed because of z = NA values when model was fitted!')
  }else{}
  
  #########################################
  ### get rid of any incomplete rows...
  #########################################
  ina <- which(is.na(z))
  if(!is.null(dfCovs)){
    ina <- unique(c(ina , which((rowSums(is.na(dfCovs)) > 0))))
  }else{}

  if(length(ina) > 0){
    print('Attention - some NA values in the data, getting removed for fitting LMM...')
    z <- z[-ina]
    if(!is.null(dfCovs)){ dfCovs <- dfCovs[-ina,,drop=FALSE] }else{}
  }else{}
  
  n <- length(z)
  
  #########################################
  ### check/preprocess data for NSFS...
  ### including scaling of cts vars and converting cat to binary vars 
  #########################################
  tmp <- preprocessdfCovs_NSFS(dfCovs = dfCovs , dfCovs_Pred = dfCovs_Pred)
  dfCovsCts <- tmp$dfCovsCts 
  dfCovsCat <- tmp$dfCovsCat 
  dfCovsCatBin <- tmp$dfCovsCatBin
  mdfCovsCtsIn <- tmp$mdfCovsCtsIn
  sddfCovsCtsIn <- tmp$sddfCovsCtsIn
  dfCovsCts_Pred <- tmp$dfCovsCts_Pred 
  dfCovsCat_Pred <- tmp$dfCovsCat_Pred
  dfCovsCatBin_Pred <- tmp$dfCovsCatBin_Pred
  covsCat <- tmp$covsCat
  
  ##############################
  ### initialise XOut...
  ##############################
  if(is.null(dfCovsCts_Pred) & is.null(dfCovsCat_Pred)){
    XOut <- matrix(NA , n , length(colnamesX_NSFS))
    n_Pred <- 0
    
    dfCovsCts_Out <- dfCovsCts
    dfCovsCatBin_Out <- dfCovsCatBin
  }else{
    if(!is.null(dfCovsCts_Pred)){
      n_Pred <- nrow(dfCovsCts_Pred)
    }else if(!is.null(dfCovsCat_Pred)){
      n_Pred <- nrow(dfCovsCat_Pred)
    }else{}
    XOut <- matrix(NA , n_Pred , length(colnamesX_NSFS))

    dfCovsCts_Out <- dfCovsCts_Pred
    dfCovsCatBin_Out <- dfCovsCatBin_Pred
  } 
  colnames(XOut) <- colnamesX_NSFS
  
  for(j in 1:length(colnamesX_NSFS)){
    if(colnamesX_NSFS[j] == 'constant'){
      XOut[,j] <- 1
    }else if(is.element(colnamesX_NSFS[j] , names(dfCovsCts))){
      ### a linear fn of cts var...
      XOut[,j] <- dfCovsCts_Out[,colnamesX_NSFS[j]]
    }else if(is.element(colnamesX_NSFS[j] , names(dfCovsCatBin))){
      ### a linear fn of cat bin var...
      XOut[,j] <- dfCovsCatBin_Out[,colnamesX_NSFS[j]]
    }else if(grepl('___' , colnamesX_NSFS[j])){
      ### an interaction of two preceding columns...
      colnames4Int <- strsplit(colnamesX_NSFS[j] , '___')[[1]]
      XOut[,j] <- apply(XOut[,colnames4Int,drop=FALSE] , 1 , prod)
    }else{
      ### must be a NS of a single variable
      covNameThis <- strsplit(colnamesX_NSFS[j] , '_NS_')[[1]]
      qKnotsThis <- as.numeric(covNameThis[-1])
      covNameThis <- covNameThis[1]
      covValsThis <- dfCovsCts_Out[,covNameThis]
      covVals4qThis <- dfCovsCts[,covNameThis]
      
      if(length(qKnotsThis) == 1){
        qIntKnotThis <- qKnotsThis
        qbdryKnotsThis <- c(0,1)
      }else{
        qIntKnotThis <- qKnotsThis[length(qKnotsThis)]

        ### bdry for this are the ones immediately above/below qIntKnotThis
        qbdryKnotsThis <- qKnotsThis[1:(length(qKnotsThis) - 1)]
        qbdryKnotsThis <- qbdryKnotsThis[order(qbdryKnotsThis)]
        qbdryKnotsThis <- c(0 , qbdryKnotsThis , 1)
        qbdryKnotsThis <- c(max(qbdryKnotsThis[qbdryKnotsThis < qIntKnotThis]) , 
                            min(qbdryKnotsThis[qbdryKnotsThis > qIntKnotThis]))
      }
      
      XOut[,j] <- makeXns(covValsThis , 
                          bdryKnots = quantile(covVals4qThis , qbdryKnotsThis) , 
                          intKnots = quantile(covVals4qThis , qIntKnotThis))[,3]
    }
  }  
  
  return(XOut)
}

#########################################################################
### function to check/preprocess data for NSFS...
#########################################################################
preprocessdfCovs_NSFS <- function(dfCovs = NULL , dfCovs_Pred = NULL){
  
  if(is.null(dfCovs)){ stop('Error - must enter dfCovs for preprocessdfCovs_NSFS!') }else{}
    
  #########################################
  ### check input variable names...  
  #########################################
  if(any(grepl('___' , names(dfCovs)))){
    stop('Error - do not include any variable names with ___ (triple underscore, for indicating interaction) in!')
  }else{}
  if(any(grepl('_NS_' , names(dfCovs)))){
    stop('Error - do not include any variable names with _NS_ (for indicating natural spline of preceding variable, with knots at following quantiles) in!')
  }else{}

  if(!is.data.frame(dfCovs)){ stop('Error - dfCovs must be entered as data.frame!') }else{}
  if((!is.null(dfCovs_Pred)) && (is.null(ncol(dfCovs_Pred)))){ stop('Error - dfCovs_Pred must be entered as data.frame!') }else{}

  if(!is.null(dfCovs_Pred)){
    if(!identical(names(dfCovs) , names(dfCovs_Pred))){ stop('Error - dfCovs_Pred must be entered with identical column names to dfCovs!') }else{}
  }else{}

  ##############################################################
  ### split into cts and cat variables...
  ### and make sure factor levels are defined the same in dfCovs and dfCovs_Pred
  ### this converts to character then to factor.
  ##############################################################
  namesCts <- c() ; namesCat <- c()
  for(covname in names(dfCovs)){
    if(is.numeric(dfCovs[[covname]])){
      namesCts <- c(namesCts , covname)
    }else{
      namesCat <- c(namesCat , covname)
      dfCovs[[covname]] <- factor(as.character(dfCovs[[covname]]))
      if(!is.null(dfCovs_Pred)){
        dfCovs_Pred[[covname]] <- factor(as.character(dfCovs_Pred[[covname]]) , levels = levels(dfCovs[[covname]]))
      }else{}
    }
  }

  if(length(namesCts) > 0){
    dfCovsCts <- dfCovs[,namesCts,drop=FALSE]
    if(!is.null(dfCovs_Pred)){
      dfCovsCts_Pred <- dfCovs_Pred[,namesCts,drop=FALSE]
    }else{
      dfCovsCts_Pred <- NULL
    }
  }else{
    dfCovsCts <- NULL
    dfCovsCts_Pred <- NULL
  }
  
  if(length(namesCat) > 0){
    dfCovsCat <- dfCovs[,namesCat,drop=FALSE]
    if(!is.null(dfCovs_Pred)){
      dfCovsCat_Pred <- dfCovs_Pred[,namesCat,drop=FALSE]
    }else{
      dfCovsCat_Pred <- NULL
    }
  }else{
    dfCovsCat <- NULL
    dfCovsCat_Pred <- NULL
  }

  ##########################################################  
  ### stdz cts vars - should keep track...
  ##########################################################  
  if(!is.null(dfCovsCts)){
    mdfCovsCtsIn <- colMeans(dfCovsCts)
    sddfCovsCtsIn <- apply(dfCovsCts , 2 , sd)
    for (j in 1:ncol(dfCovsCts)){
      dfCovsCts[,j] <- (dfCovsCts[,j] - mdfCovsCtsIn[j]) / sddfCovsCtsIn[j] 
      if(!is.null(dfCovsCts_Pred)){
        dfCovsCts_Pred[,j] <- (dfCovsCts_Pred[,j] - mdfCovsCtsIn[j]) / sddfCovsCtsIn[j] 
      }else{}
    }
  }else{}
  
  ################################################  
  ### convert categorical variables to binary...
  ################################################  
  covsCat <- c()
  if(!is.null(dfCovsCat)){
    namesTmp <- c()
    dfCovsCatBin <- dfCovsCat[,c(),drop=F]
    if(!is.null(dfCovsCat_Pred)){ dfCovsCatBin_Pred <- dfCovsCat_Pred[,c(),drop=F] }else{}
    
    for(i in 1:ncol(dfCovsCat)){
      namesTmp <- c(namesTmp , paste0(names(dfCovsCat)[i] , levels(dfCovsCat[,i])))
      for(j in 1:nlevels(dfCovsCat[,i])){
        covsCat <- c(covsCat , names(dfCovsCat)[i])
        dfCovsCatBin <- cbind(dfCovsCatBin , as.numeric(dfCovsCat[,i] == levels(dfCovsCat[,i])[j]))
        if(!is.null(dfCovsCat_Pred)){
          dfCovsCatBin_Pred <- cbind(dfCovsCatBin_Pred , as.numeric(dfCovsCat_Pred[,i] == levels(dfCovsCat[,i])[j]))
        }else{}
      }
    }
    names(dfCovsCatBin) <- namesTmp
    if(!is.null(dfCovsCat_Pred)){ names(dfCovsCatBin_Pred) <- namesTmp }else{}
    rm(namesTmp)
    
    if(is.null(dfCovsCat_Pred)){
      dfCovsCatBin_Pred <- NULL
    }else{}
    
  }else{
    dfCovsCatBin <- NULL
    dfCovsCatBin_Pred <- NULL
  }
 
  return(list('dfCovsCts' = dfCovsCts , 'dfCovsCat' = dfCovsCat , 'dfCovsCatBin' = dfCovsCatBin , 'mdfCovsCtsIn' = mdfCovsCtsIn , 'sddfCovsCtsIn' = sddfCovsCtsIn , 
              'dfCovsCts_Pred' = dfCovsCts_Pred , 'dfCovsCat_Pred' = dfCovsCat_Pred , 'dfCovsCatBin_Pred' = dfCovsCatBin_Pred , 'covsCat' = covsCat)) 
}

