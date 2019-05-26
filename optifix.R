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



