##' 
##' optifix. Optimise with fixed parameters 
##' 
##' its like optim, but with fixed parameters. 
##' 
##' specify a second argument 'fixed', a vector of TRUE/FALSE values. 
##' If TRUE, the corresponding parameter in fn() is fixed. Otherwise its 
##' variable and optimised over. 
##' 
##' The return thing is the return thing from optim() but with a couple of extra 
##' bits - a vector of all the parameters and a vector copy of the 'fixed' argument. 
##' 
##' Written by Barry Rowlingson  October 2011 
##' 
##' This file released under a CC By-SA license: 
##' http://creativecommons.org/licenses/by-sa/3.0/ 
##' 
##' and must retain the text: "Originally written by Barry Rowlingson" in comments. 
##' 
##' eg 
##' tmp <- optifix(c(1.00026,0.3),c(TRUE, FALSE), fr,method="Nelder-Mead")
##'
optifix <- function(par, fixed, fn, gr = NULL, ..., 
		method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN"), 
		lower = -Inf, upper = Inf, control = list(), hessian = FALSE){ 
	
	force(fn) 
	force(fixed) 
	.npar=length(par) 
	.fixValues = par[fixed] 
	.parStart = par[!fixed] 

	.fn <- function(par,...){ 
		.par = rep(NA,sum(!fixed)) 
		.par[!fixed] = par 
		.par[fixed] = .fixValues 

		fn(.par,...) 
	} 

	if(!is.null(gr)){ 
		.gr <- function(par,...){ 
			.gpar = rep(NA,sum(!fixed)) 
			.gpar[!fixed] = par 
			.gpar[fixed] = .fixValues 
			gr(.gpar,...)[!fixed] 
		} 
	}else{ 
		.gr <- NULL 
	} 
	
#	.opt = optim(.parStart,.fn,.gr,...,
#			method=method,lower=lower,control=control,hessian=hessian) 
### TO update (22/12/17) to subset lower/upper and to pass on upper...
      if(length(lower) > 1){
      	.opt = optim(.parStart,.fn,.gr,...,
			method=method,lower=lower[!fixed],upper = upper[!fixed] , control=control,hessian=hessian) 
      }else{
      	.opt = optim(.parStart,.fn,.gr,...,
			method=method,lower=lower,upper = upper , control=control,hessian=hessian) 
      }
	.opt$fullpars = rep(NA,sum(!fixed)) 
	.opt$fullpars[fixed]=.fixValues 
	.opt$fullpars[!fixed]=.opt$par 
	.opt$fixed = fixed 
	return(.opt) 
} 



Nelder_Mead2fix <- function(par, fn, fixed, ..., 
                    lower = -Inf, upper = Inf, control = list()){ 
  
  force(fn) 
  force(fixed) 
  .npar=length(par) 
  .fixValues = par[fixed] 
  .parStart = par[!fixed] 
  
  .fn <- function(par,...){ 
    .par = rep(NA,sum(!fixed)) 
    .par[!fixed] = par 
    .par[fixed] = .fixValues 
    
    fn(.par,...) 
  } 
  
  if(length(lower) > 1){
    .opt = Nelder_Mead2(par = .parStart, fn = .fn,...=...,
                 lower=lower[!fixed],upper = upper[!fixed] , control=control) 
  }else{
    .opt = Nelder_Mead2(par=.parStart,fn=.fn,...=...,
                 lower=lower,upper = upper , control=control) 
  }
  .opt$fullpars = rep(NA,sum(!fixed)) 
  .opt$fullpars[fixed]=.fixValues 
  .opt$fullpars[!fixed]=.opt$par 
  .opt$fixed = fixed 
  return(.opt) 
} 

Nelder_Mead2 <- function (par, fn, ... , lower = rep.int(-Inf, n), upper = rep.int(Inf, 
                                                             n), control = list()) 
{
  n <- length(par)
  if (is.null(xst <- control[["xst"]])) 
    xst <- rep.int(0.02, n)
  if (is.null(xt <- control[["xt"]])) 
    xt <- xst * 5e-04
  control[["xst"]] <- control[["xt"]] <- NULL
  if (is.null(verbose <- control[["verbose"]])) 
    verbose <- 0
  control[["verbose"]] <- NULL
  if (is.null(control[["iprint"]])) {
    control[["iprint"]] <- switch(as.character(min(as.numeric(verbose), 
                                                   3L)), `0` = 0, `1` = 20, `2` = 10, 
                                  `3` = 1)
  }
  # stopifnot(is.function(fn), length(formals(fn)) == 1L, (n <- length(par <- as.numeric(par))) == 
  #             length(lower <- as.numeric(lower)), length(upper <- as.numeric(upper)) == 
  #             n, length(xst <- as.numeric(xst)) == n, all(xst != 0), 
  #           length(xt <- as.numeric(xt)) == n)
  stopifnot(is.function(fn), (n <- length(par <- as.numeric(par))) == 
              length(lower <- as.numeric(lower)), length(upper <- as.numeric(upper)) == 
              n, length(xst <- as.numeric(xst)) == n, all(xst != 0), 
            length(xt <- as.numeric(xt)) == n)
  nM <- NelderMead$new(lower = lower, upper = upper, x0 = par, 
                       xst = xst, xt = xt)
  cc <- do.call(function(iprint = 0L, maxfun = 10000L, FtolAbs = 1e-05, 
                         FtolRel = 1e-15, XtolRel = 1e-07, MinfMax = -.Machine$double.xmax, 
                         warnOnly = FALSE, ...) {
    if (length(list(...)) > 0) 
      warning("unused control arguments ignored")
    list(iprint = iprint, maxfun = maxfun, FtolAbs = FtolAbs, 
         FtolRel = FtolRel, XtolRel = XtolRel, MinfMax = MinfMax, 
         warnOnly = warnOnly)
  }, control)
  nM$setFtolAbs(cc$FtolAbs)
  nM$setFtolRel(cc$FtolRel)
  nM$setIprint(cc$iprint)
  nM$setMaxeval(cc$maxfun)
  nM$setMinfMax(cc$MinfMax)
  it <- 0
  repeat {
    it <- it + 1
    nMres <- nM$newf(fn(nM$xeval(),...=...))
    if (nMres != 0L) 
      break
  }
  cmsg <- "reached max evaluations"
  if (nMres == -4) {
    cmsg <- warning(sprintf("failure to converge in %d evaluations", 
                            cc$maxfun))
    nMres <- 4
  }
  msgvec <- c("nm_forced", "cannot generate a feasible simplex", 
              "initial x is not feasible", "active", "objective function went below allowed minimum", 
              "objective function values converged to within tolerance", 
              "parameter values converged to within tolerance", 
              cmsg)
  if (nMres < 0) {
    (if (cc$warnOnly) 
      warning
     else stop)(msgvec[nMres + 4])
  }
  list(fval = nM$value(), par = nM$xpos(), convergence = pmin(0, 
                                                              nMres), NM.result = nMres, message = msgvec[nMres + 4], 
       control = c(cc, xst = xst, xt = xt), feval = it)
}

#########################################################################
### function to iterate optim, each time starting from previous finishing point, until no more change...
########################################################################
### to put this into the global environment, so that it can be seen in function and in optimIt without passing...
verboseOptim <<- F

#optimIt <- function(par , fn , gr = NULL , methodOptim = c("Nelder-Mead" , "L-BFGS-B") , vecFixedIn = logical(length(par)) , fitRange = matrix(NA , length(par) , 2) , ...){
optimIt <- function(par , fn , gr = NULL , methodOptim = c("L-BFGS-B" , "Nelder-Mead") , vecFixedIn = logical(length(par)) , fitRange = matrix(NA , length(par) , 2) , tolIt = 1E-4 , ...){
  # fitRange is npar x 2 with lower and upper values
  # iterate between first NM and then L-BFGS-B (or other specified algorithms) until no more improvement.
  verboseOptim <<- T
  eval(fn(par,...))
  verboseOptim <<- F
  
  fitRange[is.na(fitRange[,1]),1] <- -Inf
  fitRange[is.na(fitRange[,2]),2] <- Inf
  
  warnIn <- options()$warn
  
  # prevOF = 9E9
  prevOF = eval(fn(par,...))
  # stop when we improve by tol or less.
  stillImproving = TRUE
  listRes <- list()
  it <- 1
  parStore <- c()
  vecFixed <- vecFixedIn
  
  while (stillImproving){
    iMethodThis <- it %% length(methodOptim)
    if(iMethodThis == 0){ iMethodThis <- length(methodOptim) }else{}
    
    # set max n fn evals a bit less than the default, because we will iterate
    if(all(!vecFixedIn)){ # to avoid an extra level of wrapping, straight to optim if nothing fixed
      if(methodOptim[iMethodThis] == 'Nelder-Mead'){
        options(warn = -1) # to ignore warnings from Nelder_Mead2 about lack of convergence.
        res <- Nelder_Mead2(par = par , fn = fn , lower = fitRange[,1] , upper = fitRange[,2] , control = list(maxfun = 350 , warnOnly = TRUE) , ... = ...)
        options(warn = warnIn)
        newOF <- res$fval
        res$value <- res$fval
      }else{
        res <- optim(par = par , fn = fn , gr = gr , method = methodOptim[iMethodThis] , lower = fitRange[,1] , upper = fitRange[,2] , control = list(maxit = 70) , ... = ...)
        newOF <- res$value
      }
      verboseOptim <<- T
      eval(fn(res$par,...))
      verboseOptim <<- F
      parInitsNext <- res$par
      parStore <- cbind(parStore , res$par)
    }else{
      if(methodOptim[iMethodThis] == 'Nelder-Mead'){
        options(warn = -1) # to ignore warnings from Nelder_Mead2 about lack of convergence.
        res <- Nelder_Mead2fix(par = par , fn = fn , fixed = vecFixed , lower = fitRange[,1] , upper = fitRange[,2] , control = list(maxfun = 350 , warnOnly = TRUE) , ... = ...)
        options(warn = warnIn)
        newOF <- res$fval
        res$value <- res$fval
      }else{
        res <- optifix(par = par , fn = fn , gr = gr , fixed = vecFixed , method = methodOptim[iMethodThis] , lower = fitRange[,1] , upper = fitRange[,2] , control = list(maxit = 70) , ... = ...)
        newOF <- res$value
      }
      verboseOptim <<- T
      eval(fn(res$fullpars,...))
      verboseOptim <<- F
      parInitsNext <- res$fullpars
      parStore <- cbind(parStore , res$fullpars)
      res$par <- res$fullpars
    }
    listRes[[it]] <- res
    
    it <- it + 1
    
    if (newOF < (prevOF - tolIt)){
      stillImproving = TRUE
      prevOF <- newOF
      par <- parInitsNext
    }else{
      stillImproving = FALSE
    }
  }
  
  res$allRes <- listRes
  res$parStore <- parStore
  
  ###23/04/2020, changing to return newres...
  newres <- list('par' = res$par , 'value' = newOF , 'allRes' = listRes , 'parStore' = parStore)
  
  return(newres)
}


##########################################################
### this version quicker version of calculating lndet and invCb (compared with my previous function); 
### solve makes factorisation of C as side effect. 
### determinant uses that if chol factorisation (for dgeMatrix or matrix)
### for dspMatrix, factoristaion is BunchKaufman, not used by determinant fn, so manual calc...
### not returning chol to save mem
##########################################################
lndetANDinvCb <- function(C , b = NULL){
  if(is.null(dim(C))){ stop('Error - enter matrix for lndetANDinvCb_NEW') }else{}
  if(!exists('methodSolveTEST')){ methodSolveTEST <- 1 }else{}
  
  if(methodSolveTEST == 0){
    if(is.null(b)){
      invCb <- try(solve(C) , silent = TRUE)
    }else{
      invCb <- try(solve(C , b) , silent = TRUE)
    }
    if (is.character(invCb)){
      lndetC <- invCb <- NA
    }else{
      # if(class(C) == 'dspMatrix'){
      if(is(C , 'dspMatrix')){
        ### side effect of solve should have added BunchKaufman factorization (not used by determinant fn, so doing it here)...     
        lndetC <- sum(log(diag(C@factors$pBunchKaufman)))
      }else{
        ### for matrix, side effect of solve should have added chol factorization, used by determinant fn...     
        lndetC <- as.numeric(determinant(C , logarithm = TRUE)$modulus)
      }
    }
  }else if(methodSolveTEST == 1){
    
    invCb <- try(chol(C) , silent = TRUE) # note - not yet invCb, cholC at mo - but will be a bit later...
    if (is.character(invCb)){
      lndetC <- invCb <- NA
    }else{
      ### calc lndetC from the chol...
      lndetC <- 2 * sum(log(diag(invCb))) #
      
      if(is.null(b)){
        invCb <- chol2inv(invCb)
      }else{
        invCb<- backsolve(invCb , forwardsolve(t(invCb) , b))   
      }
    }
    
  }else{
    stop('Error - enter valid methodSolveTEST!')
  }
  
  
  return(list('lndetC' = lndetC , 'invCb' = invCb))
}

# ##########################################################
# ### function to calculate lndetC and invC b (or just invC if b not given) using cholesky...
# ##########################################################
# lndetANDinvCb_OLD <- function(C , b = NA){
#   
#   if(is.numeric(C)){
#     cholC <- try(chol(C) , silent = TRUE)
#   }else{
#     n <- nrow(C)
#     ### some issue with dgeMatrix matrices, chol(C) was different to cholC[1:n1:n], so...    
#     cholC <- try(chol(C[1:n,1:n,drop = FALSE]) , silent = TRUE)
#   }
#   if (is.character(cholC)){
#     lndetC <- invCb <- NA
#   }else{
#     lndetC <- 2 * sum(log(diag(cholC)))
#     if (is.na(as.numeric(b)[1])){        
#       invCb <- chol2inv(cholC)
#     }else{
#       ### test if C is sparse...
#       if(is.matrix(C)){    
#         ### seems to be most efficient way if dense...
#         invCb<- backsolve(cholC , forwardsolve(t(cholC) , b))   
#       }else{
#         ### but above doesn't work for sparse matrices, so...          
#         invCb <- chol2inv(cholC) %*% b 
#       }
#     }
#   }
#   
#   return(list('lndetC' = lndetC , 'invCb' = invCb, 'cholC' = cholC))
# }

# ##########################################################
# ### this version quicker; solve makes factorisation of C as side effect. 
# ### determinant uses that if chol factorisation (for dgeMatrix or matrix)
# ### for dspMatrix, factoristaion is BunchKaufman, not used by determinant fn, so manual calc...
# ### not returning chol to save mem
# ##########################################################
# lndetANDinvCb_TMP <- function(C , b = NULL){
#   if(is.null(dim(C))){ stop('Error - enter matrix for lndetANDinvCb_NEW') }else{}
#   if(is.null(b)){
#     invCb <- try(solve(C) , silent = TRUE)
#   }else{
#     invCb <- try(solve(C , b) , silent = TRUE)
#   }
#   if (is.character(invCb)){
#     lndetC <- invCb <- NA
#   }else{
#     # if(class(C) == 'dspMatrix'){
#     if(is(C , 'dspMatrix')){
#       ### side effect of solve should have added BunchKaufman factorization (not used by determinant fn, so doing it here)...     
#       lndetC <- sum(log(diag(C@factors$pBunchKaufman)))
#     }else{
#       ### for matrix, side effect of solve should have added chol factorization, used by determinant fn...     
#       lndetC <- as.numeric(determinant(C , logarithm = TRUE)$modulus)
#     }
#   }
#   
#   return(list('lndetC' = lndetC , 'invCb' = invCb))
# }
# 
