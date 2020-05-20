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

