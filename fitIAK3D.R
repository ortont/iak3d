fitIAK3D <- function(xData , dIData , zData , covsData , modelX , modelx = 'matern' , nud = 0.5 , allKnotsd = c() , 
                sdfdType_cd1 = 0 , sdfdType_cxd0 = 0 , sdfdType_cxd1 = 0 , cmeOpt = 0 , sdfdKnots = NULL , prodSum = TRUE , 
                lnTfmdData = FALSE , useReml = TRUE , compLikMats = list('compLikOptn' = 0) , namePlot = NA , 
                lmmFit = list() , rqrBTfmdPreds = TRUE , parsInit = NULL , attachBigMats = TRUE , testMC = TRUE){

  
  # if(!identical(class(xData) , 'matrix')){ stop('Stopping - for fitIAKD3D function, enter xData as a matrix') }else{}
  # if(!identical(class(dIData) , 'matrix')){ stop('Stopping - for fitIAKD3D function, enter dIData as a matrix with 2 columns') }else{}
  
########################################################
### if xData or dIData were dataframes, convert to matrices here.
### and make sure all are numeric...
########################################################
    if(!is.matrix(xData)){
        xData <- as.matrix(xData)
    }else{}
    if(!is.matrix(dIData)){
        dIData <- as.matrix(dIData)
    }else{}

    xCopy <- matrix(NA , nrow(xData) , ncol(xData))
    for (i in 1:ncol(xData)){ xCopy[,i] <- as.numeric(xData[,i]) }
    xData <- xCopy
    remove(xCopy)
    
    dICopy <- matrix(NA , nrow(dIData) , ncol(dIData))
    for (i in 1:ncol(dIData)){ dICopy[,i] <- as.numeric(dIData[,i]) }
    dIData <- dICopy
    remove(dICopy)

#########################################################
### error check for duplicated xData,dIData data, if no measurement error being included...
#########################################################
    if(cmeOpt == 0){
      iDuplicated <- which(duplicated(cbind(xData , dIData)))
      if(length(iDuplicated) > 0){ stop(paste0('Error - the data at positions ' , iDuplicated , ' are duplicated!')) }else{} 
    }else{}
    
#########################################################
### error check for any NAs; removing them if found...
#########################################################
    iNA <- which(is.na(zData) | rowSums(is.na(covsData)) > 0)
    if(length(iNA) > 0){ 
      print('Attention! Some NAs found in data ; removing them.')
      xData <- xData[-iNA,drop = FALSE]
      dIData <- dIData[-iNA,drop = FALSE]
      covsData <- covsData[-iNA,drop = FALSE]
      zData <- zData[iNA]
    }else{}

#########################################################
### round depths to nearest cm (assuming depths were inputted in m)
#########################################################
    dIData <- round(dIData , digits = 2)
    if(min(dIData[,2] - dIData[,1]) <= 0){ stop('Error - some data entered with lower depth less than or equal to upper depth. Check profiles!') }else{}
    
################################################
### if compLikMats is not null, fit assuming additive on input scale (ie put lnTfmdData = FALSE)
### ie approx involves approx by arithmetic averaging on log-transformed data + comp lik.
### but allow rqrBTfmdPreds to remain as inputted
################################################
    if(compLikMats$compLikOptn > 0){
      if (lnTfmdData){
        lnTfmdData <- FALSE
        print('Attention - averaging assumed to occur on the transformed scale.')
      }else{}
    }else{}

#######################################
### for Normal data, fit beta and cxd[d=0] automatically
#######################################
### for lnNormal data, fit beta by newton-raphson
#######################################
    if(lnTfmdData){
### in this case, only allowed reml if depth not in fixed effects, 
### so assume has to be ml...
        if(useReml){
            print('Attention! When underlying point-support variable assumed lognormal, cannot use REML because beta parameters contribute non-linearly to likelihood function!')
            print('Therefore switching option to use ML.')
            useReml <- F
        }else{}
    }else{}

###################################################
### if using Eidsvik, make sure ML, not REML...
###################################################
    if(compLikMats$compLikOptn == 2 || compLikMats$compLikOptn == 3){
      if(useReml){
        print('Attention - do not useReml with Eidsvik CL - appropriate for ML only! Changing option here.')
        useReml <- FALSE
      }else{}
    }else{}
    
################################################
### if not prodSum, will be product covariance function
### in which case, put sdfdType_cd1 <- -9 
################################################
    if(!prodSum){
        if (modelx != 'matern'){ stop('Error - for product covariance, specify matern for modelx!') }else{}
        sdfdType_cd1 <- -9 
    }else{}

################################################
### if using gam2, all knot info will be in modelX$listfefdKnots...
### other info will be in modelX$incInts and modelX$colnamesXcns
################################################
    if(identical(modelX$type , 'gam2')){
      if(length(allKnotsd) > 0){ stop('Error - if using gam2 for trend, enter knot info in modelX$listfefdKnots (you can use function makelistfefdKnots to do this); do not use allKnotsd.\n
                                      Note, modelX should also set modelX$incInts and modelX$comnamesXcns')}else{}
    }else{}

################################################
### for sdfd as spline fn (if sdfdType > 0) set knots if not given...
################################################
    if(is.null(sdfdKnots)){
      sdfdKnots <- setKnots4sdfd(dIData , sdfdType_cd1 = sdfdType_cd1 , sdfdType_cxd0 = sdfdType_cxd0 , sdfdType_cxd1 = sdfdType_cxd1)
    }else{
      if(sdfdType_cd1 > 0 | sdfdType_cxd0 > 0 | sdfdType_cxd1 > 0){
        if(is.null(sdfdKnots$bdryKnots) | ((sdfdType_cd1 > 0) & (length(sdfdKnots$intKnots_cd1) != sdfdType_cd1)) |
         ((sdfdType_cxd0 > 0) & (length(sdfdKnots$intKnots_cxd0) != sdfdType_cxd0)) | 
         ((sdfdType_cxd1 > 0) & (length(sdfdKnots$intKnots_cxd1) != sdfdType_cxd1))){
          stop('Error - sdfdKnots incorrectly defined for sdfdType - check this!')
        }else{}
      }else{}
    }
    
######################################################
### set up fixed-effect design matrices...
######################################################
    if(identical(modelX$type , 'gam2')){
      tmp <- makeXvX_gam2(covData = covsData , dIData = dIData , listfefdKnots = modelX$listfefdKnots , incInts = modelX$incInts , intMthd = modelX$intMthd , colnamesXcns = modelX$colnamesX , nDiscPts = 1000 , lnTfmdData = lnTfmdData)
      XData <- tmp$X
      vXU <- tmp$vXU
      iU <- tmp$iU
      namesX <- tmp$namesX
      XLims <- NA
    }else{
      tmp <- makeXvX(covData = covsData , dIData = dIData , modelX = modelX , allKnotsd = allKnotsd , nDiscPts = 1000 , lnTfmdData = lnTfmdData)
      XData <- tmp$X
      vXU <- tmp$vXU
      iU <- tmp$iU
      namesX <- tmp$namesX
      XLims <- tmp$XLims
    }

    remove(tmp) ; gc(verbose = FALSE)

############################################################
### check if XData has a lot of zeros, and if so, make it sparse...
############################################################
#    propnX0 <- sum(XData == 0) / (ncol(XData) * nrow(XData))
#    if(propnX0 > 0.8){
#      XData <- Matrix(XData)
#    }else{}
#    if(length(vXU) > 1){
#      propnvXU0 <- sum(vXU == 0) / (ncol(vXU) * nrow(vXU))
#      if(propnvXU0 > 0.8){
#        vXU <- Matrix(vXU)
#      }else{}
#    }else{} 
#
######################################################
### setup random-effects design matrices...
######################################################
    if (compLikMats$compLikOptn == 0){
      # setupMats <- setupIAK3D(xData , dIData , nDscPts = 0) 
      
      setupMats <- setupIAK3D(xData , dIData , nDscPts = 0 ,
                              sdfdType_cd1 = sdfdType_cd1 , sdfdType_cxd0 = sdfdType_cxd0 , sdfdType_cxd1 = sdfdType_cxd1 , sdfdKnots = sdfdKnots)
    }else{
      # setupMats <- setupIAK3D_CL(xData , dIData , nDscPts = 0 , compLikMats = compLikMats)
      setupMats <- setupIAK3D_CL(xData , dIData , nDscPts = 0 , compLikMats = compLikMats ,
                                 sdfdType_cd1 = sdfdType_cd1 , sdfdType_cxd0 = sdfdType_cxd0 , sdfdType_cxd1 = sdfdType_cxd1 , sdfdKnots = sdfdKnots)
    }
      
######################################################
### set sensible initial parameters and evaluate nll...
######################################################
    verboseOptim <<- T

    if(is.null(lmmFit$parBnds)){
        parBnds <- setParBndsIAK3D(xData = xData , dIData = dIData , setupMats = setupMats , compLikMats = compLikMats)
    }else{
        parBnds <- lmmFit$parBnds
    }

    if(is.null(lmmFit$pars)){
    
      if(is.null(parsInit)){
        parsInit <- setInitsIAK3D(xData = xData , dIData = dIData , zData = zData , XData = XData , vXU = vXU , iU = iU , modelx = modelx , nud = nud , 
            sdfdType_cd1 = sdfdType_cd1 , sdfdType_cxd0 = sdfdType_cxd0 , sdfdType_cxd1 = sdfdType_cxd1 , 
            cmeOpt = cmeOpt , prodSum = prodSum , setupMats = setupMats , parBnds = parBnds , lnTfmdData = lnTfmdData , lmmFit = lmmFit , compLikMats = compLikMats)$pars

### for exact lik, run through setInits function again, with the C from above used to give better initial beta (then better pars)...
        updateInits <- FALSE
        if(updateInits & compLikMats$compLikOptn == 0){
          parsBTfmd <- readPars(pars = parsInit , modelx = modelx , nud = nud , 
                sdfdType_cd1 = sdfdType_cd1 , sdfdType_cxd0 = sdfdType_cxd0 , sdfdType_cxd1 = sdfdType_cxd1 , 
                cmeOpt = cmeOpt , prodSum = prodSum , parBnds = parBnds , lnTfmdData = lnTfmdData)

          initC <- setCIAK3D(parsBTfmd = parsBTfmd , modelx = modelx , 
              sdfdType_cd1 = sdfdType_cd1 , sdfdType_cxd0 = sdfdType_cxd0 , sdfdType_cxd1 = sdfdType_cxd1 , 
              cmeOpt = cmeOpt , setupMats = setupMats)$C

          parsInit <- setInitsIAK3D(xData = xData , dIData = dIData , zData = zData , XData = XData , vXU = vXU , iU = iU , modelx = modelx , nud = nud , 
            sdfdType_cd1 = sdfdType_cd1 , sdfdType_cxd0 = sdfdType_cxd0 , sdfdType_cxd1 = sdfdType_cxd1 , 
            cmeOpt = cmeOpt , prodSum = prodSum , setupMats = setupMats , parBnds = parBnds , lnTfmdData = lnTfmdData , lmmFit = lmmFit , compLikMats = compLikMats , initC = initC)$pars
        }else{}
        
      }else{}
        
      verboseOptim <<- F
    }else{
      parsInit <- lmmFit$pars
    }

    if(is.null(lmmFit$fitRange)){  
      fitRange <- setFitRangeIAK3D(pars = parsInit , modelx = modelx , 
            sdfdType_cd1 = sdfdType_cd1 , sdfdType_cxd0 = sdfdType_cxd0 , sdfdType_cxd1 = sdfdType_cxd1 , 
            cmeOpt = cmeOpt , prodSum = prodSum , lnTfmdData = lnTfmdData)            
    }else{
      fitRange <- lmmFit$fitRange
    }

######################################################
### optimize parameters...
######################################################
    if(is.null(lmmFit$pars)){
### print out init pars and nll...
      print('At initial parameters...')
      start_time <- Sys.time()
      verboseOptim <<- T
      if (compLikMats$compLikOptn == 0){
        nllInit <- nllIAK3D(pars = parsInit , zData = zData , XData = XData , vXU = vXU , iU = iU , modelx = modelx , nud = nud , 
            sdfdType_cd1 = sdfdType_cd1 , sdfdType_cxd0 = sdfdType_cxd0 , sdfdType_cxd1 = sdfdType_cxd1 , 
            cmeOpt = cmeOpt , prodSum = prodSum , setupMats = setupMats , parBnds = parBnds , useReml = useReml , lnTfmdData = lnTfmdData , rtnAll = F)	
      }else{
        nllInit <- nllIAK3D_CL(pars = parsInit , zData = zData , XData = XData , modelx = modelx , nud = nud , 
            sdfdType_cd1 = sdfdType_cd1 , sdfdType_cxd0 = sdfdType_cxd0 , sdfdType_cxd1 = sdfdType_cxd1 , 
            cmeOpt = cmeOpt , prodSum = prodSum , setupMats = setupMats , parBnds = parBnds , useReml = useReml , compLikMats = compLikMats , rtnAll = F)	
      }
      verboseOptim <<- F
      end_time <- Sys.time()
      nllTime <- end_time - start_time
      
### and fit...
      if (compLikMats$compLikOptn == 0){

        if(testMC){
          print('Testing whether parallel or series version of gradient function is quicker...')
          start_time <- Sys.time()
          gr <- gradnllIAK3D(pars = parsInit , zData = zData , XData = XData , vXU = vXU , iU = iU , modelx = modelx , nud = nud , 
                             sdfdType_cd1 = sdfdType_cd1 , sdfdType_cxd0 = sdfdType_cxd0 , sdfdType_cxd1 = sdfdType_cxd1 , 
                             cmeOpt = cmeOpt , prodSum = prodSum , setupMats = setupMats , parBnds = parBnds , useReml = useReml , lnTfmdData = lnTfmdData)	
          end_time <- Sys.time()
          grTime <- end_time - start_time
          
          start_time <- Sys.time()
          gr.mc <- gradnllIAK3D.mc(pars = parsInit , zData = zData , XData = XData , vXU = vXU , iU = iU , modelx = modelx , nud = nud , 
                                   sdfdType_cd1 = sdfdType_cd1 , sdfdType_cxd0 = sdfdType_cxd0 , sdfdType_cxd1 = sdfdType_cxd1 , 
                                   cmeOpt = cmeOpt , prodSum = prodSum , setupMats = setupMats , parBnds = parBnds , useReml = useReml , lnTfmdData = lnTfmdData)	
          end_time <- Sys.time()
          grTime.mc <- end_time - start_time
          
          if(grTime.mc < grTime){
            print('Multi-core version of grad was quicker, so will use this to compute gradients of nll using all cores for parameter fitting (see function gradnllIAK3D.mc to change this).')
            gr4Optim <- gradnllIAK3D.mc
          }else{
            print('Single-core version of grad was quicker, so will use this to compute gradients of nll for parameter fitting.')
            gr4Optim <- gradnllIAK3D
          }
        }else{
          gr4Optim <- gradnllIAK3D
        }
        print('Now fitting...')
        parsFit <- optimIt(par = parsInit , fn = nllIAK3D , gr = gr4Optim , fitRange = fitRange ,
            zData = zData , XData = XData , vXU = vXU , iU = iU , modelx = modelx , nud = nud ,
            sdfdType_cd1 = sdfdType_cd1 , sdfdType_cxd0 = sdfdType_cxd0 , sdfdType_cxd1 = sdfdType_cxd1 ,
            cmeOpt = cmeOpt , prodSum = prodSum , setupMats = setupMats , parBnds = parBnds , useReml = useReml , lnTfmdData = lnTfmdData)
        ### all args taken from this env now.
        # parsFit <- optimIt(par = parsInit , fn = nllIAK3D , gr = gr4Optim , fitRange = fitRange)
      }else{
        
        if(testMC){
          print('Testing whether parallel or series version of gradient function is quicker...')
          start_time <- Sys.time()
          gr_CL <- gradnllIAK3D_CL(pars = parsInit , zData = zData , XData = XData , modelx = modelx , nud = nud , 
                                 sdfdType_cd1 = sdfdType_cd1 , sdfdType_cxd0 = sdfdType_cxd0 , sdfdType_cxd1 = sdfdType_cxd1 , 
                                 cmeOpt = cmeOpt , prodSum = prodSum , setupMats = setupMats , parBnds = parBnds , useReml = useReml , compLikMats = compLikMats)	
          end_time <- Sys.time()
          grTime <- end_time - start_time
          
          start_time <- Sys.time()
          gr_CL <- gradnllIAK3D_CL.mc(pars = parsInit , zData = zData , XData = XData , modelx = modelx , nud = nud , 
                                   sdfdType_cd1 = sdfdType_cd1 , sdfdType_cxd0 = sdfdType_cxd0 , sdfdType_cxd1 = sdfdType_cxd1 , 
                                   cmeOpt = cmeOpt , prodSum = prodSum , setupMats = setupMats , parBnds = parBnds , useReml = useReml , compLikMats = compLikMats)	
          end_time <- Sys.time()
          grTime.mc <- end_time - start_time
          
          if(grTime.mc < grTime){
            print('Multi-core version of grad was quicker, so will use this to compute gradients of nll using all cores for parameter fitting (see function gradnllIAK3D_CL.mc to change this).')
            gr4Optim <- gradnllIAK3D_CL.mc
          }else{
            print('Single-core version of grad was quicker, so will use this to compute gradients of nll for parameter fitting.')
            gr4Optim <- gradnllIAK3D_CL
          }
        }else{
          gr4Optim <- gradnllIAK3D_CL
        }
        print('Now fitting...')
        parsFit <- optimIt(par = parsInit , fn = nllIAK3D_CL , gr = gr4Optim , fitRange = fitRange ,
            zData = zData , XData = XData , modelx = modelx , nud = nud ,
            sdfdType_cd1 = sdfdType_cd1 , sdfdType_cxd0 = sdfdType_cxd0 , sdfdType_cxd1 = sdfdType_cxd1 ,
            cmeOpt = cmeOpt , prodSum = prodSum , setupMats = setupMats , parBnds = parBnds , useReml = useReml , compLikMats = compLikMats)
        ### args now taken from this env.
        # parsFit <- optimIt(par = parsInit , fn = nllIAK3D_CL , gr = gr4Optim , fitRange = fitRange)
      }
      
### note that parameters have now been fitted...
      print('Parameters fitted...')

    }else{
### or to just take the values in inits, which were the transformed values taken from the inputted lmmFit...
      parsFit <- list('par' = parsInit)
    }

    verboseOptim <<- T
    if (compLikMats$compLikOptn == 0){
      lmmFit <- nllIAK3D(pars = parsFit$par , zData = zData , XData = XData , vXU = vXU , iU = iU , modelx = modelx , nud = nud , 
                sdfdType_cd1 = sdfdType_cd1 , sdfdType_cxd0 = sdfdType_cxd0 , sdfdType_cxd1 = sdfdType_cxd1 , 
                cmeOpt = cmeOpt , prodSum = prodSum , setupMats = setupMats , parBnds = parBnds , useReml = useReml , lnTfmdData = lnTfmdData , rtnAll = T , attachBigMats = attachBigMats)	
    }else{
      lmmFit <- nllIAK3D_CL(pars = parsFit$par , zData = zData , XData = XData , modelx = modelx , nud = nud , 
                sdfdType_cd1 = sdfdType_cd1 , sdfdType_cxd0 = sdfdType_cxd0 , sdfdType_cxd1 = sdfdType_cxd1 , 
                cmeOpt = cmeOpt , prodSum = prodSum , setupMats = setupMats , parBnds = parBnds , useReml = useReml , compLikMats = compLikMats , rtnAll = T , attachBigMats = attachBigMats)	
    }
    verboseOptim <<- F

############################################
### also attach everything that was used to call the function to lmmFit...
############################################
    lmmFit <- attachSettings2lmmFit(lmmFit = lmmFit , xData = xData , dIData = dIData , zData = zData , covsData = covsData , 
                                    modelX = modelX , namesX = namesX , modelx = modelx , nud = nud , allKnotsd = allKnotsd , XLims = XLims , 
                                    sdfdType_cd1 = sdfdType_cd1 , sdfdType_cxd0 = sdfdType_cxd0 , sdfdType_cxd1 = sdfdType_cxd1 , cmeOpt = cmeOpt , 
                                    prodSum = prodSum , sdfdKnots = sdfdKnots , setupMats = setupMats , lnTfmdData = lnTfmdData , 
                                    useReml = useReml , parBnds = parBnds , fitRange = fitRange , compLikMats = compLikMats)

print('After attaching...')
print(lmmFit$sdfdKnots)
print(head(lmmFit$setupMats$XsdfdSplineU_cxd0))

    if (!is.na(namePlot)){
##############################################
### a quick prediction routine...
### predict through the profile for all data locations and for one distant location (with average of covsData)
##############################################
######################################################
### set up fixed-effect design matrices...
######################################################
#      dIPred <- cbind(seq(0 , 1.99 , 0.01) , seq(0.01 , 2 , 0.01))
#      dIPred <- cbind(seq(0 , 1.98 , 0.02) , seq(0.02 , 2 , 0.02))
#      dIPred <- cbind(seq(0 , 1.9 , 0.1) , seq(0.1 , 2 , 0.1))
#      dIPred <- cbind(seq(0 , 14.98 , 0.02) , seq(0.02 , 15 , 0.02))
      dIPred <- cbind(seq(0 , parBnds$maxd-0.02 , 0.02) , seq(0.02 , parBnds$maxd , 0.02))

### below 1m, only do every 5th one of these depth intervals...
      iKeep1 <- which(dIPred[,1] < 1)
      iKeep2 <- which(dIPred[,1] >= 1)
      if(length(iKeep2) > 0){ iKeep2 <- iKeep2[seq(1,length(iKeep2),5)] }else{}
      dIPred <- dIPred[c(iKeep1,iKeep2),,drop=FALSE]
      rm(iKeep1,iKeep2)

      dIPred <- round(dIPred , digits = 2)
      ndIPred <- nrow(dIPred)

      iTmp <- which(!duplicated(xData))
      
      maxnProfPlot <- 100
      if(length(iTmp) > maxnProfPlot){
### take a subsample of 200 profiles to plot...      
### so that repeatable, don't make it random...
        iTmp2 <- round(seq(1 , length(iTmp) , length(iTmp) / maxnProfPlot))
        iTmp <- iTmp[iTmp2]
        profNamesPlot <- paste0('Profile ' , iTmp2)
      }else{
        profNamesPlot <- paste0('Profile ' , seq(length(iTmp)))
      }
      
      covsPred <- covsData[iTmp,,drop=FALSE]
      xPred <- xData[iTmp,,drop=FALSE]

### add a distant location to the prediction locations...
      xPredDistant <- matrix(c(9E99 , 9E99) , 1 , 2)
#      covsPredDistant <- colMeansAndModes(covsPred)
#      covsPredDistant <- colMeansAndModes(covsData)


### initialise with first row...      
      covsPredDistant <- covsData[1,,drop=FALSE]
### don't include the depths in finding medoid...
      idCovs <- which(names(covsData) == 'dIMidPts')
      if(length(idCovs) > 0){ covsPredDistant[1,idCovs] <- NA }else{}
      if(length(idCovs) < ncol(covsData)){ covsPredDistant[1,-idCovs] <- medoid(covsData[,-idCovs,drop=FALSE]) }else{}
            
      nxPred <- dim(xPred)[[1]]
        
      # return(list('xTmp' = rbind(xPred , xPredDistant) , 'dIMap' = dIPred , 'covsMap' = rbind(covsPred , covsPredDistant) , 'lmmFit' = lmmFit))
      
### call the predict function...
      tmp <- predictIAK3D(xMap = rbind(xPred , xPredDistant) , dIMap = dIPred , covsMap = rbind(covsPred , covsPredDistant) , lmmFit = lmmFit , rqrBTfmdPreds = rqrBTfmdPreds , constrainX4Pred = FALSE)

      zPred <- tmp$zMap[,1:nxPred]
      vPred <- tmp$vMap[,1:nxPred]
      pi90UPred <- tmp$pi90UMap[,1:nxPred]
      pi90LPred <- tmp$pi90LMap[,1:nxPred]

      zPredDistant <- tmp$zMap[,nxPred+1]
      vPredDistant <- tmp$vMap[,nxPred+1]

### to back transform data for mapping...
      if(lnTfmdData & rqrBTfmdPreds){
###(were inputted and used all the way through on the log scale, 
### but for plot put back to real scale)        
        zData <- exp(zData)    
      }else{}

#################################################    
### make a pdf with the distant profile prediction (page 1) and all data profiles (6 per page thereafter)...
### note these are not validation predictions, they are predicted at the data profiles given the data for the same profiles
### also note, they are a bi-product of the method (ie you can use the method to predict at profiles where we have data), 
### not to be confused with the spline-then-krige type approach where similar plots may be produced in the first step of analysis
#################################################    
      tmp <- plotProfilesIAK3D(namePlot = namePlot , xData = xData , dIData = dIData , zData = zData , 
                xPred = xPred , dIPred = dIPred , zPred = zPred , pi90LPred = pi90LPred , pi90UPred = pi90UPred , zPredDistant = zPredDistant , profNames = profNamesPlot)

    }else{
      zPredDistant <- vPredDistant <- zPred <- vPred <- pi90LPred <- pi90UPred <- xPred <- dIPred <- NA
    }

#####################################################
### and return...
#####################################################
    return(list('parsFit' = parsFit , 'lmmFit' = lmmFit , 
                'zPredDistant' = zPredDistant , 'vPredDistant' = vPredDistant , 
                'zPred' = zPred , 'vPred' = vPred , 'pi90LPred' = pi90LPred , 'pi90UPred' = pi90UPred , 
                'xPred' = xPred , 'dIPred' = dIPred))
}

###################################################################
### the nll function that is to be minimized...
###################################################################
nllIAK3D <- function(pars , zData , XData , vXU , iU , modelx , nud ,  
                     sdfdType_cd1 , sdfdType_cxd0 , sdfdType_cxd1 , cmeOpt , prodSum , setupMats = NULL , parBnds , useReml , lnTfmdData , rtnAll = F , forCompLik = FALSE , attachBigMats = FALSE){

### default values for return...
    badnll <- 9E99
    cxdhat <- NA
    WiAW <- NA
    lndetA <- NA

    p <- ncol(XData)
    n <- length(zData)
    
    if(exists("printnllTime") && printnllTime){
      ptm <- proc.time()
    }else{}
    
    if(exists("parsTrace4Optim") && (!forCompLik)){
    
      if(is.null(ncol(parsTrace4Optim)) || (ncol(parsTrace4Optim) == length(pars))){

        if(length(parsTrace4Optim) > 0){
          parsTmp <- rbind(parsTrace4Optim , pars)
          nllTmp <- c(nllTrace4Optim , NA) # for this trace, use NA for failed evaluation.
        }else{
          parsTmp <- matrix(pars , nrow = 1)
          nllTmp <- NA # for this trace, use NA for failed evaluation.
        }
        assign("parsTrace4Optim" , parsTmp , envir = .GlobalEnv)
        assign("nllTrace4Optim" , nllTmp , envir = .GlobalEnv)
      }else{}
      
    }else{}

    parsBTfmd <- readPars(pars = pars , modelx = modelx , nud = nud , 
                sdfdType_cd1 = sdfdType_cd1 , sdfdType_cxd0 = sdfdType_cxd0 , sdfdType_cxd1 = sdfdType_cxd1 , 
                cmeOpt = cmeOpt , prodSum = prodSum , parBnds = parBnds , lnTfmdData = lnTfmdData)

    if(any(is.na(c(parsBTfmd$sdfdPars_cd1,parsBTfmd$sdfdPars_cxd0,parsBTfmd$sdfdPars_cxd1))) || 
       any(is.nan(c(parsBTfmd$sdfdPars_cd1,parsBTfmd$sdfdPars_cxd0,parsBTfmd$sdfdPars_cxd1))) ||
       any(is.infinite(c(parsBTfmd$sdfdPars_cd1,parsBTfmd$sdfdPars_cxd0,parsBTfmd$sdfdPars_cxd1)))){
      
      parsOK <- FALSE
      A <- NA
      sigma2Vec <- NA
      p <- ncol(XData)
      if(is.null(p)){ p <- 1 }else{}
      
    }else{
      
      tmp <- setCIAK3D(parsBTfmd = parsBTfmd , modelx = modelx , 
                         sdfdType_cd1 = sdfdType_cd1 , sdfdType_cxd0 = sdfdType_cxd0 , sdfdType_cxd1 = sdfdType_cxd1 , 
                         cmeOpt = cmeOpt , setupMats = setupMats)

      A <- tmp$C
      sigma2Vec <- tmp$sigma2Vec
      remove(tmp) ; 

      parsOK <- T
      if(max(is.na(A)) == 1){ parsOK <- F }else{}
      if(is.infinite(parsBTfmd$cx0) | is.infinite(parsBTfmd$cx1) | is.infinite(parsBTfmd$cd1) | is.infinite(parsBTfmd$cxd0) | is.infinite(parsBTfmd$cxd1)){ parsOK <- F }else{}
    }
    
    if(parsOK){
        n <- length(zData)
        XData <- as.matrix(XData , nrow = n)
        p <- ncol(XData)
        W <- cbind(XData,zData)

### just to make sure numerical errors have not made it non-symmetric...  		
        if(!is(A , 'dspMatrix')){ A <- forceSymmetric(A) }else{}

        if(lnTfmdData){
            tmp <- nrUpdatesIAK3DlnN(zData = zData , XData = XData , vXU = vXU , iU = iU , C = A , sigma2Vec = sigma2Vec)
            nll <- tmp$nll
            betahat <- tmp$betahat
            if(rtnAll){ 
#                vbetahat <- solve(tmp$fim) 
### or, alternatively, 
                vXUTmp <- matrix(0 , n * p , p)
                vXUTmp[kronecker(seq(0 , (n - 1) * p , p) , matrix(1 , length(iU) , 1)) + rep(iU , n) , iU] <- vXU
                tmp <- gradHessIAK3DlnN(zData = zData , XData = XData , vXU = vXUTmp , iU = seq(p) , betaU = betahat , C = A , lndetC = 0 , TK  = tmp$iC , iXKiCXKXKiC = c() , sigma2Vec = sigma2Vec)

                vbetahat <- solve(tmp$fim) 

            }else{ 
                vbetahat <- NA 
            } 
            cxdhat <- 1
            C <- A
        }else{

         iAW <- lndetANDinvCb(A , W)
         lndetA <- iAW$lndetC
         iAW <- iAW$invCb
         if(is.na(lndetA) | is.infinite(lndetA)){
            print('Inversion of A failed!')
            print(paste0('nll = ' , badnll , '; nmPars = ' , paste(round(pars , digits = 3) , collapse = ', ')))
            nll <- badnll
            cxdhat <- NA
            if(forCompLik){
              return(list('pars' = pars , 'parsBTfmd' = parsBTfmd , 'sigma2Vec' = sigma2Vec , 'WiAW' = matrix(NA , p+1 , p+1) , 'lndetA' = NA))
            }else{}
            
            if(rtnAll){ 
### -999 * A returned because it is not the covariance matrix yet, as cxhat was not computed. Still allows A to be retrieved if required though. 
              if(attachBigMats){
                return(list('nll' = nll , 'pars' = pars , 'parsBTfmd' = parsBTfmd , 'betahat' = NA , 'vbetahat' = NA , 'cxdhat' = NA , 
                            'C' = -999 * A , 'sigma2Vec' = sigma2Vec , 'XData' = XData , 'vXU' = vXU , 'iU' = iU))
              }else{
                return(list('nll' = nll , 'pars' = pars , 'parsBTfmd' = parsBTfmd , 'betahat' = NA , 'vbetahat' = NA , 'cxdhat' = NA , 
                            'sigma2Vec' = sigma2Vec , 'XData' = XData , 'vXU' = vXU , 'iU' = iU))
              }
            }else{
        	    return(nll)
            }
          }else{}
          
          iAW <- matrix(iAW , ncol = ncol(W))

          WiAW <- t(W) %*% iAW
          if(!is(WiAW , 'dspMatrix')){ WiAW <- forceSymmetric(WiAW) }else{}
          
          if(forCompLik){
            return(list('pars' = pars , 'parsBTfmd' = parsBTfmd , 'sigma2Vec' = sigma2Vec , 'WiAW' = WiAW , 'lndetA' = lndetA , 'A' = A , 'iAW' = iAW))
          }else{}

          XiAX <- WiAW[1:p , 1:p , drop = FALSE]
          XiAz <- WiAW[1:p , p+1 , drop = FALSE]
          ziAz <- as.numeric(WiAW[p+1 , p+1])

          betahat <- lndetANDinvCb(XiAX , XiAz)
          lndetXiAX <- betahat$lndetC
          betahat <- betahat$invCb
          if(is.na(lndetXiAX) | is.infinite(lndetXiAX)){
            
              print('Inversion of t(XData) iA XData failed!')
              print(paste0('nll = ' , badnll , '; nmPars = ' , paste(round(pars , digits = 3) , collapse = ', ')))
              nll <- badnll
              if(rtnAll){ 
### -999 * A returned because it is not the covariance matrix yet, as cxhat was not computed. Still allows A to be retrieved if required though. 
                if(attachBigMats){
                  return(list('nll' = nll , 'pars' = pars , 'parsBTfmd' = parsBTfmd , 'betahat' = NA , 'vbetahat' = NA , 'cxdhat' = NA , 
                            'C' = -999 * A , 'sigma2Vec' = sigma2Vec , 'XData' = XData , 'vXU' = vXU , 'iU' = iU))
                }else{
                  return(list('nll' = nll , 'pars' = pars , 'parsBTfmd' = parsBTfmd , 'betahat' = NA , 'vbetahat' = NA , 'cxdhat' = NA , 
                              'sigma2Vec' = sigma2Vec , 'XData' = XData , 'vXU' = vXU , 'iU' = iU))
                }
              }else{
        	    return(nll)
              }
          }else{}
            
          betahat <- matrix(betahat , ncol = 1)
          
          resiAres <- as.numeric(ziAz - 2 * t(betahat) %*% XiAz + t(betahat) %*% XiAX %*% betahat)

          if(useReml){
              cxdhat <- resiAres / (n - p)
          }else{
              cxdhat <- resiAres / n
          }
### update the values in parsBTfmd so that it makes sense when returned to user.          
          parsBTfmd$cx0 <- cxdhat * parsBTfmd$cx0
          parsBTfmd$cx1 <- cxdhat * parsBTfmd$cx1
          parsBTfmd$cd1 <- cxdhat * parsBTfmd$cd1
          parsBTfmd$cxd0 <- cxdhat * parsBTfmd$cxd0
          parsBTfmd$cxd1 <- cxdhat * parsBTfmd$cxd1
          if(cmeOpt == 1){ parsBTfmd$cme <- cxdhat * parsBTfmd$cme }else{}
                
          lndetC <- lndetA + n * log(cxdhat)
          lndetXiCX <- lndetXiAX - p * log(cxdhat)
          resiCres <- resiAres / cxdhat

          if(rtnAll){ 
              vbetahat <- cxdhat * solve(XiAX) 
              C <- cxdhat * A
              sigma2Vec <- cxdhat * sigma2Vec 
          }else{}
          
          if(useReml){
              nll <- 0.5 * ((n - p) * log(2 * pi) + lndetC + lndetXiCX + resiCres) 
          }else{
              nll <- 0.5 * (n * log(2 * pi) + lndetC + resiCres) 
          }
          nll <- as.numeric(nll)

      }
      
    }else{
        nll <- badnll
        if(forCompLik){
          return(list('pars' = pars , 'parsBTfmd' = parsBTfmd , 'sigma2Vec' = sigma2Vec , 'WiAW' = matrix(NA , p+1 , p+1) , 'lndetA' = NA))
        }else{}
    }

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
      print('time for nll evaluation was:')
      print(proc.time() - ptm)
    }else{}
    
    if(exists("parsTrace4Optim") && (length(nllTrace4Optim) > 0) && (!forCompLik)){
      nllTmp <- nllTrace4Optim # update the final value, now fn has been successful.
      if(nll != badnll){
        nllTmp[length(nllTmp)] <- nll
        assign("nllTrace4Optim" , nllTmp , envir = .GlobalEnv)
      }else{}
    }else{}

    if(rtnAll){ 
### added extra info for quicker prediction, 28/2/19...
        if(!lnTfmdData){
          muhat <- XData %*% betahat 
        }else{
          muhat <- setmuIAK3D(XData = XData , vXU = vXU , iU = iU , beta = betahat , diagC = diag(C) , sigma2Vec = sigma2Vec) 
        }

#        tmp <- lndetANDinvCb(C , cbind(XData , zData - muhat))
#    	return(list('nll' = nll , 'pars' = pars , 'parsBTfmd' = parsBTfmd , 'betahat' = betahat , 'vbetahat' = vbetahat , 'cxdhat' = cxdhat , 
#                  'C' = C , 'sigma2Vec' = sigma2Vec , 'XData' = XData , 'vXU' = vXU , 'iU' = iU , 'iCX' = tmp$invCb[,1:p,drop=FALSE] , 'iCz_muhat' = tmp$invCb[,p+1,drop=FALSE] , 'iC' = chol2inv(tmp$cholC)))

      if(attachBigMats){
        
        iC <- chol2inv(chol(C))
        iCXRes <- matrix(iC %*% cbind(XData , zData - muhat) , ncol = p+1)
        
      	return(list('nll' = nll , 'pars' = pars , 'parsBTfmd' = parsBTfmd , 'betahat' = betahat , 'vbetahat' = vbetahat , 'cxdhat' = cxdhat , 
                  'C' = C , 'sigma2Vec' = sigma2Vec , 'XData' = XData , 'vXU' = vXU , 'iU' = iU , 'iCX' = iCXRes[,1:p,drop=FALSE] , 'iCz_muhat' = iCXRes[,p+1,drop=FALSE] , 'iC' = iC))
      }else{
        return(list('nll' = nll , 'pars' = pars , 'parsBTfmd' = parsBTfmd , 'betahat' = betahat , 'vbetahat' = vbetahat , 'cxdhat' = cxdhat , 
                    'sigma2Vec' = sigma2Vec , 'XData' = XData , 'vXU' = vXU , 'iU' = iU))
      }                
    }else{
    	return(nll)
    }
}	

### updated from the above function, 28/2/20-5/3/20, to use dspMatrix class for more efficient memory.
setCIAK3D <- function(parsBTfmd , modelx , 
                        sdfdType_cd1 , sdfdType_cxd0 , sdfdType_cxd1 , cmeOpt , setupMats){

### check sdfdTypes - if any > 0, use the approx version...  
  if(max(c(sdfdType_cd1 , sdfdType_cxd0 , sdfdType_cxd1)) > 0){
    return(setCApproxIAK3D(parsBTfmd = parsBTfmd , modelx = modelx , 
                     sdfdType_cd1 = sdfdType_cd1 , sdfdType_cxd0 = sdfdType_cxd0 , sdfdType_cxd1 = sdfdType_cxd1 , cmeOpt = cmeOpt , setupMats = setupMats))
  }else{}
  
### if setupMats doesn't include 'Kx', it should include 'xData' and 'dIData', so that proper setupMats can be made now...
  if(is.null(setupMats$Kx)){
    if(max(c(sdfdType_cd1 , sdfdType_cxd0 , sdfdType_cxd1)) > 0){ stop('Not ready to run setCIAK3D without passing in setupMats yet, because this requires passing in sdfdKnots - update the function!') }else{}
    setupMats <- setupIAK3D(xData = setupMats$xData , dIData = as.data.frame(setupMats$dIData) , nDscPts = 0 , 
                            sdfdType_cd1 = sdfdType_cd1 , sdfdType_cxd0 = sdfdType_cxd0 , sdfdType_cxd1 = sdfdType_cxd1) 
  }else{}
  
  n <- nrow(setupMats$Kx)
  C <- sparseMatrix(i=seq(n),j=seq(n),x=parsBTfmd$cme , symmetric = TRUE) # initiate with just meas err on diag. , dims = as.integer(c(n,n))

  ### adding Cx0...
  C <- C + parsBTfmd$cx0 * forceSymmetric(setupMats$Kx %*% t(setupMats$Kx))

  if (modelx == 'matern'){
    phix1 <- maternCov(setupMats$Dx , c(1 , parsBTfmd$ax , parsBTfmd$nux))
  }else if(modelx == 'nugget'){
    phix1 <- 0
  }else{
    stop('Not ready!')
  }
  
  ### if any components are stationary, function only neds to be run once...
  iStat <- which(c(sdfdType_cd1 , sdfdType_cxd0 , sdfdType_cxd1) == 0)
  if(length(iStat) > 0){
    tmp <- iaCovMatern(dIData = setupMats$dIU , ad = parsBTfmd$ad , nud = parsBTfmd$nud , sdfdPars = c(1 , 0 , 1) , sdfdType = 0 , abcd = setupMats$dIUabcd , iUElements = setupMats$dIUiUElements)
    phidStat <- tmp$avCovMtx
    avVardStat <- tmp$avVarVec
    remove(tmp) ;
  }else{}
  
  if(sdfdType_cd1 == 0){
    phid_cd1 <- phidStat
    avVard_cd1 <- avVardStat
  }else if(sdfdType_cd1 == -9){
    phid_cd1 <- 0
    avVard_cd1 <- 0
  }else{
    tmp <- iaCovMatern(dIData = setupMats$dIU , ad = parsBTfmd$ad , nud = parsBTfmd$nud , sdfdPars = parsBTfmd$sdfdPars_cd1 , sdfdType = sdfdType_cd1 , abcd = setupMats$dIUabcd , iUElements = setupMats$dIUiUElements)
    phid_cd1 <- tmp$avCovMtx
    avVard_cd1 <- tmp$avVarVec
    remove(tmp) ; 
  }
  
  if(sdfdType_cxd0 == 0){
    phid_cxd0 <- phidStat
    avVard_cxd0 <- avVardStat
  }else{
    tmp <- iaCovMatern(dIData = setupMats$dIU , ad = parsBTfmd$ad , nud = parsBTfmd$nud , sdfdPars = parsBTfmd$sdfdPars_cxd0 , sdfdType = sdfdType_cxd0 , abcd = setupMats$dIUabcd , iUElements = setupMats$dIUiUElements)
    phid_cxd0 <- tmp$avCovMtx
    avVard_cxd0 <- tmp$avVarVec
    remove(tmp) ; 
  }

  if(modelx == 'nugget'){
    phid_cxd1 <- 0
    avVard_cxd1 <- 0
  }else{
    if(sdfdType_cxd1 == 0){
      phid_cxd1 <- phidStat
      avVard_cxd1 <- avVardStat
    }else{
      tmp <- iaCovMatern(dIData = setupMats$dIU , ad = parsBTfmd$ad , nud = parsBTfmd$nud , sdfdPars = parsBTfmd$sdfdPars_cxd1 , sdfdType = sdfdType_cxd1 , abcd = setupMats$dIUabcd , iUElements = setupMats$dIUiUElements)
      phid_cxd1 <- tmp$avCovMtx
      avVard_cxd1 <- tmp$avVarVec
      remove(tmp) ; 
    }
  }
  
  parsOK <- T
  if((max(is.na(phid_cd1)) == 1) | (max(is.na(phid_cxd0)) == 1) | (max(is.na(phid_cxd1)) == 1)){ parsOK <- F } else{}
  
  if(parsOK){
    ### adding Cxd0...
    C <- C + parsBTfmd$cxd0 * sparseMatrix(i=setupMats$summKxKx$i , j=setupMats$summKxKx$j,x = phid_cxd0@x[setupMats$utriKdIdxdKd][setupMats$summKxKx$idxUtri],symmetric = TRUE)
    remove(phid_cxd0)

    ### without using diagBlocks fn...    
    if(length(iStat) > 0){
      KphidKStat <- new("dspMatrix" , Dim = as.integer(c(nrow(setupMats$Kd),nrow(setupMats$Kd))), x =  phidStat@x[setupMats$utriKdIdxdKd])
    }else{}

### adding Cd...
    if(sdfdType_cd1 == 0){
      C <- C + parsBTfmd$cd1 * KphidKStat
    }else if(sdfdType_cd1 == -9){
      Cd <- 0
    }else{
      C <- C + new("dspMatrix" , Dim = as.integer(c(nrow(setupMats$Kd),nrow(setupMats$Kd))), x =  phid_cd1@x[setupMats$utriKdIdxdKd])
    }
    remove(phid_cd1)

    if(modelx == 'nugget'){
      if(exists('KphidKStat')){ remove(KphidKStat) ;  }else{}
      Cx1 <- Cxd1 <- 0
    }else{
      if(sdfdType_cxd1 != 0){ # if not rqd, free up mem...
        if(exists('KphidKStat')){ remove(KphidKStat) ;  }else{}
      }else{}
      
      Kphix1K <- new("dspMatrix" , Dim = as.integer(c(nrow(setupMats$Kx),nrow(setupMats$Kx))), x =  phix1@x[setupMats$utriKxIdxxKx])
      rm(phix1)
### adding Cx1...
      C <- C + parsBTfmd$cx1 * Kphix1K

### adding Cxd1...
      if(sdfdType_cxd1 == 0){
        C <- C + parsBTfmd$cxd1 * Kphix1K * KphidKStat
      }else{
        if(exists('KphidKStat')){ remove(KphidKStat) ;  }else{}
        C <- C + parsBTfmd$cxd1 * Kphix1K * new("dspMatrix" , Dim = as.integer(c(nrow(setupMats$Kd),nrow(setupMats$Kd))), x =  phid_cxd1@x[setupMats$utriKdIdxdKd])
      }

      remove(Kphix1K) ; 
    }
    if(exists('KphidKStat')){ remove(KphidKStat) ;  }else{}
    if(exists('phid_cxd1')){ remove(phid_cxd1) ;  }else{}

    ### add the variances...
    sigma2Vec <- parsBTfmd$cxd0 * avVard_cxd0 + parsBTfmd$cxd1 * avVard_cxd1 + parsBTfmd$cme + parsBTfmd$cx0 + parsBTfmd$cx1 + parsBTfmd$cd1 * avVard_cd1

    ### expand for all sampled depths (from the unique ones)...    
    sigma2Vec <- setupMats$Kd %*% sigma2Vec

  }else{        
    C <- NA
    sigma2Vec <- NA
    #    print('Bad parameters...not sure how this can happen')
#    print(parsBTfmd)
  }
  
  return(list('C' = C , 'sigma2Vec' = sigma2Vec))
  #	return(list('C' = C , 'sigma2Vec' = sigma2Vec , 'dC' = list(Cx0 , Cx1 , Cd , Cxd0 , Cxd1) , 'phi' = list(phix0 , phix1 , phid_cd1 , phid_cxd0 , phid_cxd1)))
}



################################################################
### use setCIAK3D2 for covs between different sets of data (setupMats will also contain xU2, dIU2, etc)
################################################################
setCIAK3D2 <- function(parsBTfmd , modelx , 
                       sdfdType_cd1 , sdfdType_cxd0 , sdfdType_cxd1 , cmeOpt , setupMats){

  ### check sdfdTypes - if any > 0, use the approx version...  
  if(max(c(sdfdType_cd1 , sdfdType_cxd0 , sdfdType_cxd1)) > 0){
    return(setCApproxIAK3D2(parsBTfmd = parsBTfmd , modelx = modelx , 
                     sdfdType_cd1 = sdfdType_cd1 , sdfdType_cxd0 = sdfdType_cxd0 , sdfdType_cxd1 = sdfdType_cxd1 , cmeOpt = cmeOpt , setupMats = setupMats))
  }else{}
  
  ### if setupMats doesn't include 'Kx', it should include 'xData', 'dIData', 'xData2' and 'dIData2', so that proper setupMats can be made now...
  if(is.null(setupMats$Kx)){
    if(max(c(sdfdType_cd1 , sdfdType_cxd0 , sdfdType_cxd1)) > 0){ stop('Not ready to run setCIAK3D2 without passing in setupMats yet, because this requires passing in sdfdKnots - update the function!') }else{}
    setupMats <- setupIAK3D2(xData = setupMats$xData , dIData = setupMats$dIData , xData2 = setupMats$xData2 , dIData2 = setupMats$dIData2 , 
                            sdfdType_cd1 = sdfdType_cd1 , sdfdType_cxd0 = sdfdType_cxd0 , sdfdType_cxd1 = sdfdType_cxd1) 
  }else{}
  
  ### starting with Cx0...
  ### defined for the unique locations...
  ijTmp <- which(setupMats$Dx == 0 , arr.ind = TRUE)
  if(nrow(ijTmp) > 0){
    phix0 <- sparseMatrix(i = ijTmp[,1] , j = ijTmp[,2] , x = 1 , dims = c(nrow(setupMats$Dx) , ncol(setupMats$Dx)))
    C <- parsBTfmd$cx0 * setupMats$Kx %*% phix0 %*% t(setupMats$Kx2)
  }else{
    phix0 <- sparseMatrix(i = 1 , ,j=1 , x=0 , dims = c(nrow(setupMats$Dx),ncol(setupMats$Dx)))
    C <- sparseMatrix(i=1,j=1,x=0,dims = c(nrow(setupMats$Kx),nrow(setupMats$Kx2)))
  }
  
  if (modelx == 'matern'){
    phix1 <- maternCov(setupMats$Dx , c(1 , parsBTfmd$ax , parsBTfmd$nux))
  }else if(modelx == 'nugget'){
    phix1 <- 0
  }else{
    stop('Not ready!')
  }
  
  ### if any components are stationary, function only neds to be run once...
  iStat <- which(c(sdfdType_cd1 , sdfdType_cxd0 , sdfdType_cxd1) == 0)
  if(length(iStat) > 0){
    phidStat <- iaCovMatern2(dIData = setupMats$dIU , dIData2 = setupMats$dIU2 , ad = parsBTfmd$ad , nud = parsBTfmd$nud , sdfdPars = c(1 , 0 , 1) , sdfdType = 0)
  }else{}
  
  if(sdfdType_cd1 == 0){
    phid_cd1 <- phidStat
  }else if(sdfdType_cd1 == -9){
    phid_cd1 <- 0
  }else{
    phid_cd1 <- iaCovMatern2(dIData = setupMats$dIU , dIData2 = setupMats$dIU2 , ad = parsBTfmd$ad , nud = parsBTfmd$nud , sdfdPars = parsBTfmd$sdfdPars_cd1 , sdfdType = sdfdType_cd1)
  }
  
  if(sdfdType_cxd0 == 0){
    phid_cxd0 <- phidStat
  }else{
    phid_cxd0 <- iaCovMatern2(dIData = setupMats$dIU , dIData2 = setupMats$dIU2 , ad = parsBTfmd$ad , nud = parsBTfmd$nud , sdfdPars = parsBTfmd$sdfdPars_cxd0 , sdfdType = sdfdType_cxd0)
  }
  
  if(modelx == 'nugget'){
    phid_cxd1 <- 0
  }else{
    if(sdfdType_cxd1 == 0){
      phid_cxd1 <- phidStat
    }else{
      phid_cxd1 <- iaCovMatern2(dIData = setupMats$dIU , dIData2 = setupMats$dIU2 , ad = parsBTfmd$ad , nud = parsBTfmd$nud , sdfdPars = parsBTfmd$sdfdPars_cxd1 , sdfdType = sdfdType_cxd1)
    }
  }
  
  parsOK <- T
  if((max(is.na(phid_cd1)) == 1) | (max(is.na(phid_cxd0)) == 1) | (max(is.na(phid_cxd1)) == 1)){ parsOK <- F } else{}
  
  if(parsOK){
    
### adding Cxd0...    
    if(length(iStat) > 0){
      
      if(modelx == 'nugget' & sdfdType_cd1 == -9){
        KphidKStat <- (setupMats$Kx %*% phix0 %*% t(setupMats$Kx2)) * (setupMats$Kd %*% phidStat %*% t(setupMats$Kd2))
        C <- C + parsBTfmd$cxd0 * (setupMats$Kx %*% phix0 %*% t(setupMats$Kx2)) * (setupMats$Kd %*% phid_cxd0 %*% t(setupMats$Kd2))
      }else{        
        KphidKStat <- setupMats$Kd %*% phidStat %*% t(setupMats$Kd2)	
        C <- C + parsBTfmd$cxd0 * (setupMats$Kx %*% phix0 %*% t(setupMats$Kx2)) * (setupMats$Kd %*% phid_cxd0 %*% t(setupMats$Kd2))
      }
    }else{
      C <- C + parsBTfmd$cxd0 * (setupMats$Kx %*% phix0 %*% t(setupMats$Kx2)) * (setupMats$Kd %*% phid_cxd0 %*% t(setupMats$Kd2))
    }
    
### adding Cd...    
    if(sdfdType_cd1 == 0){
      C <- C + parsBTfmd$cd1 * KphidKStat
    }else if(sdfdType_cd1 == -9){
      Cd <- 0
    }else{
      C <- C + parsBTfmd$cd1 * (setupMats$Kd %*% phid_cd1 %*% t(setupMats$Kd2))
    }

    if(modelx == 'nugget'){
      Cx1 <- Cxd1 <- 0
    }else{
      Kphix1K <- setupMats$Kx %*% phix1 %*% t(setupMats$Kx2)
      ### adding Cx1...    
      C <- C + parsBTfmd$cx1 * Kphix1K 
      ### adding Cxd1...    
      if(sdfdType_cxd1 == 0){
        C <- C + parsBTfmd$cxd1 * Kphix1K * KphidKStat 
      }else{
        C <- C + parsBTfmd$cxd1 * Kphix1K * (setupMats$Kd %*% phid_cxd1 %*% t(setupMats$Kd2))
      }
      remove(Kphix1K) ; 
    }
    if(exists('KphidKStat')){ remove(KphidKStat) ;  }else{}
    
    ### note, no meas err being added to cov between distinct data...
    
  }else{        
    C <- NA
#    print('Bad parameters...not sure how this can happen')
#    print(parsBTfmd)
  }
  
  return(C)
}


setmuIAK3D <- function(XData , vXU , iU , beta , diagC , sigma2Vec){
###########################################################
### sigma2Vec is a length-n vector with the average variances (not average covariances) in each sample support
### with stationary variances, these would all be equal to the sum of all variance parameters. 
### with non-stationary variances, these will vary.
### this function is only used in the lognormal setting
###########################################################
    betavXbeta <- calcbetavXbeta(vXU , beta[iU])
    
    mu <- XData %*% beta + 0.5 * (sigma2Vec - diagC + betavXbeta)

    return(mu)
}

######################################################################
### function for calculating (Kx %*% t(Kx)) * (Kd %*% M %*% t(Kd)) 
### where one entry of 1 in each row of Kx and in each row of Kd
### ie only the diagonal blocks of (Kd %*% M %*% t(Kd)) rqd
### useful for nugget.
######################################################################
diagBlocksKdMKd <- function(setupMats , M){  
  newcoladds <- as.numeric(setupMats$Kx %*% seq(0,ncol(setupMats$Kx)-1) * ncol(setupMats$Kd))
  ijx <- summary(setupMats$Kd)
  ijx <- ijx[order(ijx$i),]
  ijx$j <- ijx$j + newcoladds
  KdB <- sparseMatrix(i = ijx$i , j = ijx$j , x = ijx$x , dims = c(nrow(setupMats$Kd) , ncol(setupMats$Kd) * ncol(setupMats$Kx)))

### M is dense matrix, as is Kd %*% M...
### this bit could be improved so that sorting isn't needed, and so that M could be a list of similarly sized matrices.
  newcoladds <- rep(setupMats$Kx %*% seq(0,ncol(setupMats$Kx)-1) * ncol(setupMats$Kd) , each = ncol(setupMats$Kd))
  
  if(is.list(M)){

    iTmp <- rep(seq(nrow(setupMats$Kd)) , each = ncol(setupMats$Kd)) 
    jTmp <- rep(seq(ncol(setupMats$Kd)), nrow(setupMats$Kd)) + newcoladds 
                 
    KdMKd <- list()
    for(i in 1:length(M)){
      KdMKd[[i]] <- sparseMatrix(i = iTmp , j = jTmp , x = as.numeric(t(setupMats$Kd %*% M[[i]])) , dims = c(nrow(setupMats$Kd) , ncol(setupMats$Kd) * ncol(setupMats$Kx)))
      KdMKd[[i]] <- KdMKd[[i]] %*% t(KdB)
    }  
  }else{
    KdMKd <- sparseMatrix(i = rep(seq(nrow(setupMats$Kd)) , each = ncol(setupMats$Kd)) , 
                          j = rep(seq(ncol(setupMats$Kd)), nrow(setupMats$Kd)) + newcoladds , 
                          x = as.numeric(t(setupMats$Kd %*% M)) , dims = c(nrow(setupMats$Kd) , ncol(setupMats$Kd) * ncol(setupMats$Kx)))
    KdMKd <- KdMKd %*% t(KdB)
  }
  return(KdMKd)
}
	

calcbetavXbeta <- function(vX , beta){
#################################################
### vX is n * p x p    
### vX could contain covariances for all covariates, but this would not be required,
### since only those that vary in sample supports (indexed by iU) contribute to this term
### so usually vXU and betaU will be passed in to function
### this function is only used in the lognormal setting
#################################################
    p <- dim(vX)[[2]]
    n <- dim(vX)[[1]] / p
    vXbeta <- vX %*% beta

    i1 <- kronecker(seq(n) , matrix(1 , p , 1))
    j1 <- seq(n * p)

    tKTEMP <- sparseMatrix(i = i1 , j = j1 , x = 1)
    betavXbeta <- tKTEMP %*% (kronecker(matrix(1, n , 1) , beta) * vXbeta)

    return(betavXbeta) 
}

    
readPars <- function(pars , modelx , nud , sdfdType_cd1 , sdfdType_cxd0 , sdfdType_cxd1 , cmeOpt , prodSum , parBnds , lnTfmdData){

    inext <- 1
    parsOK <- T
    if (modelx == 'matern'){
      ax <- logitab(pars[inext] , parBnds$axBnds[1] , parBnds$axBnds[2] , invt = T); inext <- inext + 1  
      nux <- logitab(pars[inext] , parBnds$nuxBnds[1] , parBnds$nuxBnds[2] , invt = T) ; inext <- inext + 1 
    }else if(modelx == 'nugget'){
### nothing to read in.    
        ax <- nux <- NA
    }else{
        stop('Build in options for other spatial correlation models')
    }

    ad <- logitab(pars[inext] , parBnds$adBnds[1] , parBnds$adBnds[2] , invt = T); inext <- inext + 1  

    if(prodSum){    
      cx0 <- 1E-8 + exp(pars[inext]) ; inext <- inext + 1 ; 
      if(modelx != 'nugget'){
        cx1 <- 1E-8 + exp(pars[inext]) ; inext <- inext + 1 ; 
      }else{
        cx1 <- 0        
      }
    }else{
      cx0 <- 0
      cx1 <- 0
    }
    if(sdfdType_cd1 != -9){        
      cd1 <- 1E-8 + exp(pars[inext]) ; inext <- inext + 1 ; 
    }else{
      cd1 <- 0
    }

    if(modelx != 'nugget'){
      if(!lnTfmdData){
        cxd1 <- logitab(pars[inext] , parBnds$sxd1Bnds[1] , parBnds$sxd1Bnds[2] , invt = T) ; inext <- inext + 1 
        cxd0 <- 1 - cxd1
      }else{ # for lnN methods, read in all variance parameters (as lnci)...
        cxd0 <- 1E-8 + exp(pars[inext]) ; inext <- inext + 1 
        cxd1 <- 1E-8 + exp(pars[inext]) ; inext <- inext + 1 
      }
    }else{
      if(!lnTfmdData){
        cxd0 <- 1
        cxd1 <- 0
      }else{
        cxd0 <- 1E-8 + exp(pars[inext]) ; inext <- inext + 1 
        cxd1 <- 0
      }        
    }
   
    tmp <- sdfdRead(pars , sdfdType_cd1 , inext , parBnds = parBnds)
    sdfdPars_cd1 <- tmp$sdfdPars
    inext <- tmp$inext

    tmp <- sdfdRead(pars , sdfdType_cxd0 , inext , parBnds = parBnds)
    sdfdPars_cxd0 <- tmp$sdfdPars
    inext <- tmp$inext

    tmp <- sdfdRead(pars , sdfdType_cxd1 , inext , parBnds = parBnds)
    sdfdPars_cxd1 <- tmp$sdfdPars
    inext <- tmp$inext

    if (cmeOpt == 1){
      cme <- exp(pars[inext])
    }else{
      cme <- 0
    }

    return(list('ax' = ax , 'nux' = nux , 'ad' = ad , 'nud' = nud , 
                'cx0' = cx0 , 'cx1' = cx1 , 'cd1' = cd1 , 'cxd0' = cxd0 , 'cxd1' = cxd1 ,
                'sdfdPars_cd1' = sdfdPars_cd1 , 'sdfdPars_cxd0' = sdfdPars_cxd0 , 'sdfdPars_cxd1' = sdfdPars_cxd1 ,
                'cme' = cme))
}

##########################################
### the non-stationary standard deviation function, evaluated at depths d
### used if approximating average covariances numerically with dicretization points
### not being used in the current code, but left in as might be used for checking in future.
##########################################
sdfd <- function(d , sdfdPars , sdfdType , XsdfdSpline){

    if (sdfdType == 0){
        sdfdVal <- matrix(1 , length(d) , 1)
    }else if(sdfdType == -1){
        sdfdVal <- (1 - sdfdPars[1]) + sdfdPars[1] * exp(-sdfdPars[2] * d)
    }else if (sdfdType == -9){
        sdfdVal <- matrix(1 , length(d) , 1)
    }else if(sdfdType > 0){
      # sdfdVal <- XsdfdSpline %*% matrix(sdfdPars , ncol = 1) # assuming XsdfdSpline was created on point support
      sdfdVal <- XsdfdSpline %*% matrix(c(1 , sdfdPars) , ncol = 1) # assuming XsdfdSpline was created on point support
    }else{
        stop('sdfd option not yet programmed')
    }

### check sdfdDsc, return just NA if not ok...
    if (min(sdfdVal) < 0){ sdfdVal <- NA }else{}
	
    return(sdfdVal)
}

##########################################
### increment-averaged version of the sdfd
##########################################
iasdfd <- function(dI , sdfdPars , sdfdType , XsdfdSpline){
  
  if (sdfdType == 0){
    sdfdVal <- matrix(1 , nrow(dI) , 1)
  }else if(sdfdType == -1){
    # sdfdVal <- (1 - sdfdPars[1]) + sdfdPars[1] * exp(-sdfdPars[2] * d)
    sdfdVal <- (1 - sdfdPars[1]) - (sdfdPars[1] / sdfdPars[2]) * (exp(-sdfdPars[2] * dI[,2]) - exp(-sdfdPars[2] * dI[,1])) / (dI[,2] - dI[,1])
  }else if (sdfdType == -9){
    sdfdVal <- matrix(1 , nrow(d) , 1)
  }else if(sdfdType > 0){
    # sdfdVal <- XsdfdSpline %*% matrix(sdfdPars , ncol = 1)  # assuming XsdfdSpline was created on ia support
    sdfdVal <- XsdfdSpline %*% matrix(c(1 , sdfdPars) , ncol = 1)  # assuming XsdfdSpline was created on ia support
  }else{
    stop('sdfd option not yet programmed')
  }
  
  ### check sdfdDsc, return just NA if not ok...
  # if (min(sdfdVal) < 0){ sdfdVal <- NA }else{}
  
  return(sdfdVal)
}


sdfdRead <- function(pars , sdfdType , inext , parBnds){
    if (sdfdType == 0){
        sdfdPars <- c()
    }else if(sdfdType == -1){
#############################################
### sdfd <- (1 - tau1) + tau1 * exp(-(tau2 * d))
### to ensure (1 - tau1) + tau1 * exp(-(tau2 * d)) > 0, 
### make tau2 > 0 and furthermore in bounds stated in parBnds$tau2Bnds...param by par2 = logitab(tau2 , tau2min, tau2max)
### and to make sure it is positive for all values in range d = 0 to d = dmax...
### make (1 - tau1) + tau1 * exp(-(tau2 * dmax)) > 0...param by par1 = log[(1 - tau1) + tau1 * exp(-(tau2 * dmax))]
### inverse of latter is tau1 = (1 - exp(par1)) / (1 - exp(-tau2 * dmax))
#############################################
        sdfdPars <- tauTfm(parsIn = pars[inext:(inext+1)] , parBnds = parBnds , invt = T)

        inext <- inext + 2
    }else if(sdfdType > 0){
### natural spline, clamped to grad = 0 at upper bdry knot; sdfdType internal knots (sdfdType free pars if forced to 1 at d=0)       

      sdfdPars <- pars[inext:(inext+sdfdType-1)] 
      inext <- inext + sdfdType
      
      # stop('Not ready yet.')      
      
    }else if(sdfdType == -9){
        sdfdPars <- c()
    }else{
        stop('sdfd option not yet programmed')
    }
	
    return(list('sdfdPars' = sdfdPars , 'inext' = inext))
}


tauTfm <- function(parsIn , parBnds , invt = F){
  nparsIn <- length(parsIn)
  if(nparsIn == 2){
    parsOut <- c(NA , NA)
    if(!invt){
      parsOut[1] <- log(1 - parsIn[1] * (1 - exp(-parsIn[2] * parBnds$maxd)))
      parsOut[2] <- logitab(parsIn[2] , parBnds$tau2Bnds[1]  , parBnds$tau2Bnds[2])
    }else{
      parsOut[2] <- logitab(parsIn[2] , parBnds$tau2Bnds[1]  , parBnds$tau2Bnds[2] , invt = T)
      parsOut[1] <- (1 - exp(parsIn[1])) / (1 - exp(-parsOut[2] * parBnds$maxd))
    }
  }else if(nparsIn == 0){
    parsOut <- c()  
  }else{
    stop('Error - not ready for transforming a tau vector with this number of parameters!')
  }  
  return(parsOut)
}

updateLmmFitNames <- function(lmmFit){
  ix <- which(names(lmmFit) == 'x')
  if (length(ix) == 1){ names(lmmFit)[ix] <- 'xData' }else{}
  idI <- which(names(lmmFit) == 'dI')
  if (length(idI) == 1){ names(lmmFit)[idI] <- 'dIData' }else{}
  iz <- which(names(lmmFit) == 'z')
  if (length(iz) == 1){ names(lmmFit)[iz] <- 'zData' }else{}
  icovs <- which(names(lmmFit) == 'covs')
  if (length(icovs) == 1){ names(lmmFit)[icovs] <- 'covsData' }else{}
  iX <- which(names(lmmFit) == 'X')
  if (length(iX) == 1){ names(lmmFit)[iX] <- 'XData' }else{}
  
  # if(class(lmmFit$modelX) == "cubist"){
  if(is(lmmFit$modelX , "cubist")){
    lmmFit$modelX$allKnotsd <- lmmFit$allKnotsd
    if(length(lmmFit$modelX$allKnotsd) > 0){
      nparSpline <- ncol(lmmFit$XData) - length(lmmFit$modelX$names4XCubist)
      lmmFit$modelX$namesX <- c(lmmFit$modelX$names4XCubist , paste0('dSpline.' , seq(nparSpline)))
    }else{
      lmmFit$modelX$namesX <- lmmFit$modelX$names4XCubist
    }
    lmmFit$modelX$namesDataFit <- names(lmmFit$covsData)
    lmmFit$modelX$comms4XCubist <- rep(1 , length(lmmFit$modelX$names4XCubist))
  }else{}

  return(lmmFit)
}

attachSettings2lmmFit <- function(lmmFit , xData , dIData , zData , covsData , modelX , namesX , modelx , nud , allKnotsd , XLims , 
               sdfdType_cd1 , sdfdType_cxd0 , sdfdType_cxd1 , cmeOpt , prodSum , sdfdKnots , setupMats , lnTfmdData , useReml , parBnds , fitRange , compLikMats){
    lmmFit$xData <- xData
    lmmFit$dIData <- dIData
    lmmFit$zData <- zData
    lmmFit$covsData <- covsData
    lmmFit$modelX <- modelX
    lmmFit$namesX <- namesX
    lmmFit$modelx <- modelx
    lmmFit$nud <- nud
    if(length(allKnotsd) > 0){
      lmmFit$allKnotsd <- allKnotsd
    }else{ # assigning as c() doesn't work, so: 
      lmmFit['allKnotsd'] <- list(NULL)
    }
    lmmFit$XLims <- XLims
    lmmFit$sdfdType_cd1 <- sdfdType_cd1
    lmmFit$sdfdType_cxd0 <- sdfdType_cxd0
    lmmFit$sdfdType_cxd1 <- sdfdType_cxd1
    lmmFit$cmeOpt <- cmeOpt
    lmmFit$prodSum <- prodSum
    if(!is.null(sdfdKnots)){
      lmmFit$sdfdKnots <- sdfdKnots
    }else{  
      lmmFit['sdfdKnots'] <- list(NULL)
    }
    lmmFit$setupMats <- setupMats
    lmmFit$lnTfmdData <- lnTfmdData
    lmmFit$useReml <- useReml
    lmmFit$parBnds <- parBnds
    lmmFit$fitRange <- fitRange
    lmmFit$compLikMats <- compLikMats
    
    return(lmmFit)
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
### function to calculate lndetC and invC b (or just invC if b not given) using cholesky...
##########################################################
lndetANDinvCb_OLD <- function(C , b = NA){

    if(is.numeric(C)){
      cholC <- try(chol(C) , silent = TRUE)
    }else{
      n <- nrow(C)
### some issue with dgeMatrix matrices, chol(C) was different to cholC[1:n1:n], so...    
      cholC <- try(chol(C[1:n,1:n,drop = FALSE]) , silent = TRUE)
    }
    if (is.character(cholC)){
        lndetC <- invCb <- NA
    }else{
        lndetC <- 2 * sum(log(diag(cholC)))
        if (is.na(as.numeric(b)[1])){        
            invCb <- chol2inv(cholC)
        }else{
### test if C is sparse...
          if(is.matrix(C)){    
### seems to be most efficient way if dense...
            invCb<- backsolve(cholC , forwardsolve(t(cholC) , b))   
          }else{
### but above doesn't work for sparse matrices, so...          
            invCb <- chol2inv(cholC) %*% b 
          }
        }
   }

   return(list('lndetC' = lndetC , 'invCb' = invCb, 'cholC' = cholC))
}

##########################################################
### this version quicker; solve makes factorisation of C as side effect. 
### determinant uses that if chol factorisation (for dgeMatrix or matrix)
### for dspMatrix, factoristaion is BunchKaufman, not used by determinant fn, so manual calc...
### not returning chol to save mem
##########################################################
lndetANDinvCb <- function(C , b = NULL){
  if(is.null(dim(C))){ stop('Error - enter matrix for lndetANDinvCb_NEW') }else{}
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
  
  return(list('lndetC' = lndetC , 'invCb' = invCb))
}

##########################################################
### logit transform function for range [a - b]...
##########################################################
logitab <- function(z , a = 0 , b = 1 , invt = F){
	if(!invt){
		p <- (z - a) / (b - a)
		z <- log(p / (1 - p))
	}else{
		p <- 1 / (1 + exp(-z))
		z <- p * (b - a) + a
	}
	return(z)
}   

##########################################################
### a transform function for the tau parameters...
### sdfd(d=maxd) / sdfd(d=0) = [1 - tau1 ( 1 - exp{-tau2 maxd})] / 1
### bounds can be put on this ratio (sdfdmaxdBnds)
##########################################################
tauTfm <- function(parsIn , parBnds , invt = F){
  nparsIn <- length(parsIn)
  if(nparsIn == 2){
    parsOut <- c(NA , NA)
    if(!invt){
      parsOut[1] <- log(1 - parsIn[1] * (1 - exp(-parsIn[2] * parBnds$maxd)))
      parsOut[2] <- logitab(parsIn[2] , parBnds$tau2Bnds[1]  , parBnds$tau2Bnds[2])
    }else{
      parsOut[2] <- logitab(parsIn[2] , parBnds$tau2Bnds[1]  , parBnds$tau2Bnds[2] , invt = T)
      parsOut[1] <- (1 - exp(parsIn[1])) / (1 - exp(-parsOut[2] * parBnds$maxd))
    }
  }else if(nparsIn == 0){
    parsOut <- c()  
  }else{
    stop('Error - not ready for transforming a tau vector with this number of parameters!')
  }  
  return(parsOut)
}

########################################################################
### sdfd(d=maxd) / sdfd(d=0) = [1 - tau1 ( 1 - exp{-tau2 maxd})] / 1
### bounds can be put on this ratio (sdfdmaxdBnds)
### only used in defining initial parameters... 
### in DEFUNCT version,
### parsIn here are on the transformed scale...
### with tfmed pars, sd(d) = 1 - [1-exp(tau1Tfmd)] * [1 - exp(-tau2 d)] / [1 - exp(-tau2 maxd)]
### and              sd(maxd) = exp(tau1Tfmd)    
########################################################################
contraintau1_DEFUNCT <- function(tau1Tfmd , maxd , sdfmaxdBnds = c(0 , Inf)){
    sdTmpmaxd <- exp(tau1Tfmd)
    if(sdTmpmaxd < sdfmaxdBnds[1]){
      tau1Tfmd <- log(sdfmaxdBnds[1])
    }else if(sdTmpmaxd > sdfmaxdBnds[2]){
      tau1Tfmd <- log(sdfmaxdBnds[2])
    }else{}
    return(tau1Tfmd)
}

#############################################
### in this version, tau is on real scale...
#############################################
contraintau <- function(tau , parBnds , sdfmaxdBnds = c(0 , Inf)){
### check ok for tau2...
    if(tau[2] <= parBnds$tau2Bnds[1]){ 
        tau[2] <- parBnds$tau2Bnds[1] + 0.01 * (parBnds$tau2Bnds[2] - parBnds$tau2Bnds[1]) 
    }else if(tau[2] >= parBnds$tau2Bnds[1]){ 
        tau[2] <- parBnds$tau2Bnds[1] + 0.99 * (parBnds$tau2Bnds[2] - parBnds$tau2Bnds[1]) 
    }else{}

### check ok for tau1...
    sdTmpmaxd <- 1 - tau[1] * (1 - exp(-tau[2] * parBnds$maxd))

    if(sdTmpmaxd < sdfmaxdBnds[1]){ 
        tau[1] <- (1 - sdfmaxdBnds[1]) / (1 - exp(-tau[2] * parBnds$maxd)) 
    }else if(sdTmpmaxd > sdfmaxdBnds[2]){ 
        tau[1] <- (1 - sdfmaxdBnds[2]) / (1 - exp(-tau[2] * parBnds$maxd)) 
    }else{}

    return(tau)
}



############################################################
### uses modes of any factors, means of any non-factors...
############################################################
colMeansAndModes <- function(dfIn , na.rm = FALSE){
  if(ncol(dfIn) > 0){
    jFactor <- c()
    for (j in 1:ncol(dfIn)){
      # if(class(dfIn[,j]) == 'factor'){
      if(is(dfIn[,j] , 'factor')){
        jFactor <- c(jFactor , j)
      }else{}
    }
    jNonFactor <- setdiff(1:ncol(dfIn) , jFactor)

    cm <- dfIn[1,]
    if(length(jNonFactor) > 0){
      cm[jNonFactor] <- colMeans(dfIn[,jNonFactor] , na.rm = na.rm)
    }else{}

    if(length(jFactor) > 0){  
      for(j in 1:length(jFactor)){  
        lThis <- levels(dfIn[,jFactor[j]])
        cm[jFactor[j]] <- lThis[which.max(tabulate(match(dfIn[,jFactor[j]] , lThis)))]
      }
    }else{}
  }else{
    cm <- NULL
  }
  return(cm)
}

###########################################################################
### a function to get the sample that is closest to the 'centre'...
###########################################################################
medoid <- function(dfIn){

  if(ncol(dfIn) > 0){
    dfInOrig <- dfIn
 
########################################################################  
### standardise numeric covariates and take abs distance from centre...
### for factors, use 1 - count / countmax, stdardized to have mean = avD
########################################################################  
    avD <- NA * numeric(ncol(dfIn))
    for (j in 1:ncol(dfIn)){
      if(is.numeric(dfIn[,j])){
        dfIn[,j] <- dfIn[,j,drop=FALSE] - mean(dfIn[,j])
        dfIn[,j] <- abs(dfIn[,j,drop=FALSE]) / sqrt(var(dfIn[,j]))

        avD[j] <- mean(dfIn[,j])
      
      }else if(is.factor(dfIn[,j])){
### nothing for now.    
        tableTmp <- table(dfInOrig[,j])
        dfIn[,j] <- NA
        dfIn[,j] <- as.numeric(dfIn[,j])
        for(i in 1:nrow(tableTmp)){
          iThis <- which(dfInOrig[,j] == rownames(tableTmp)[i])
          if(length(iThis) > 0){
            dfIn[iThis,j] <- 1 - tableTmp[i]/max(tableTmp)
          }else{}
        }
      }else{
        stop('Error in medoid function - currently only written for numeric or factor variables!')
      }
    }
  
    if(all(is.na(avD))){ 
      mavD <- 1 
    }else{
      mavD <- mean(avD , na.rm = TRUE)
    }
  
############################################################
### define each sample with a distance from centre
### distance = sum of all univariate distances (city block metric)
### distance for factor is 0 if most common class
############################################################
    jFactors <- which(is.na(avD))
    if(length(jFactors) > 0){
      for (j in jFactors){
        dfIn[,j] <- dfIn[,j] * mavD / mean(dfIn[,j])
      }
    }

### get row with min sum of dists and return that row of the original df...
    imed <- which.min(rowSums(dfIn))

    return(dfInOrig[imed,,drop=FALSE])
  }else{
    return(dfIn[1,,drop=FALSE])
  }
}
 

###################################################################
### compute gradients of nll by finite diff with mc...
###################################################################
gradnllIAK3D.mc <- function(pars , zData , XData , vXU , iU , modelx , nud ,  
                     sdfdType_cd1 , sdfdType_cxd0 , sdfdType_cxd1 , cmeOpt , prodSum , setupMats , parBnds , useReml , lnTfmdData){

  if(length(pars) == 1){ stop('Error - do not use gradnllIAK3D.mc with 1 parameter!') }else{}

  delVec <- rep(1E-6 , length(pars))
  numCores <- min(detectCores() , 8)
  
  parsMtx <- matrix(pars , length(pars) , length(pars) + 1)
  parsMtx[1:length(pars),1:length(pars)] <- parsMtx[1:length(pars),1:length(pars)] +  diag(delVec)
  parsList <- split(parsMtx, rep(1:ncol(parsMtx), each = nrow(parsMtx)))

  if (Sys.info()[1] == "Windows"){
    listRqd4nll <- c("setupIAK3D" , "nllIAK3D" , "readPars" , "setCIAK3D" , "setCApproxIAK3D" ,  
                          "printnllTime" , "verboseOptim" , "logitab" , "sdfdRead" , "tauTfm" , "iasdfd" ,
                          "maternCov" , "iaCovMatern" , "intRy" , "intRy12" , "intRy34" , "intTy" , "gammafnInt" , "poly2Exp" , "diagBlocksKdMKd" ,
                          "sparseMatrix" , "t" , "summary" , "diag" , "determinant" , "lndetANDinvCb" , "xyDist")
    cl <- makeCluster(numCores)
    clusterExport(cl, listRqd4nll)
    clusterEvalQ(cl, library("Matrix"))
    
    outtmp <- clusterApply(cl=cl, x=parsList, fun=nllIAK3D , 
                           zData = zData , XData = XData , vXU = vXU , iU = iU , modelx = modelx , nud = nud ,  
                           sdfdType_cd1 = sdfdType_cd1 , sdfdType_cxd0 = sdfdType_cxd0 , sdfdType_cxd1 = sdfdType_cxd1 , cmeOpt = cmeOpt , 
                           prodSum = prodSum , setupMats = setupMats , parBnds = parBnds , useReml = useReml , lnTfmdData = lnTfmdData , rtnAll = F , forCompLik = FALSE)

    stopCluster(cl) 
  }else{
    outtmp <- mclapply(X = parsList , FUN = nllIAK3D , mc.cores = numCores , 
                       zData = zData , XData = XData , vXU = vXU , iU = iU , modelx = modelx , nud = nud ,  
                       sdfdType_cd1 = sdfdType_cd1 , sdfdType_cxd0 = sdfdType_cxd0 , sdfdType_cxd1 = sdfdType_cxd1 , cmeOpt = cmeOpt , 
                       prodSum = prodSum , setupMats = setupMats , parBnds = parBnds , useReml = useReml , lnTfmdData = lnTfmdData , rtnAll = F , forCompLik = FALSE)
  }

  nllVec <- unlist(outtmp)

  grad <- (nllVec[1:length(pars)] - nllVec[length(pars)+1]) / delVec
  grad[is.na(grad)] <- 0 # set any NAs to grad 0

  return(grad)
}

gradnllIAK3D_CL.mc <- function(pars , zData , XData , modelx , nud ,  
                     sdfdType_cd1 , sdfdType_cxd0 , sdfdType_cxd1 , cmeOpt , prodSum , setupMats , parBnds , useReml , compLikMats){

  if(length(pars) == 1){ stop('Error - do not use gradnllIAK3D.mc with 1 parameter!') }else{}

  delVec <- rep(1E-6 , length(pars))
  numCores <- min(detectCores() , 8)
  
  parsMtx <- matrix(pars , length(pars) , length(pars) + 1)
  parsMtx[1:length(pars),1:length(pars)] <- parsMtx[1:length(pars),1:length(pars)] +  diag(delVec)
  parsList <- split(parsMtx, rep(1:ncol(parsMtx), each = nrow(parsMtx)))

  if (Sys.info()[1] == "Windows"){
    listRqd4nll <- c("setupIAK3D" , "nllIAK3D" , "nllIAK3D_CL" , "readPars" , "setCIAK3D" , "setCApproxIAK3D" , 
                           "printnllTime" , "verboseOptim" , "logitab" , "sdfdRead" , "tauTfm" , "iasdfd" ,
                           "maternCov" , "iaCovMatern" , "intRy" , "intRy12" , "intRy34" , "intTy" , "gammafnInt" , "poly2Exp" , "diagBlocksKdMKd" ,
                           "sparseMatrix" , "t" , "summary" , "diag" , "determinant" , "lndetANDinvCb" , "xyDist")
    cl <- makeCluster(numCores)
    clusterExport(cl, listRqd4nll)
    clusterEvalQ(cl, library("Matrix"))
    
    outtmp <- clusterApply(cl=cl, x=parsList, fun=nllIAK3D_CL , 
                       zData = zData , XData = XData , modelx = modelx , nud = nud ,  
                       sdfdType_cd1 = sdfdType_cd1 , sdfdType_cxd0 = sdfdType_cxd0 , sdfdType_cxd1 = sdfdType_cxd1 , cmeOpt = cmeOpt , 
                       prodSum = prodSum , setupMats = setupMats , parBnds = parBnds , useReml = useReml , compLikMats = compLikMats , rtnAll = F)

    stopCluster(cl) 
  }else{
    outtmp <- mclapply(X = parsList , FUN = nllIAK3D_CL , mc.cores = numCores , 
                       zData = zData , XData = XData , modelx = modelx , nud = nud ,  
                       sdfdType_cd1 = sdfdType_cd1 , sdfdType_cxd0 = sdfdType_cxd0 , sdfdType_cxd1 = sdfdType_cxd1 , cmeOpt = cmeOpt , 
                       prodSum = prodSum , setupMats = setupMats , parBnds = parBnds , useReml = useReml , compLikMats = compLikMats , rtnAll = F)
  }

  nllVec <- unlist(outtmp)

  grad <- (nllVec[1:length(pars)] - nllVec[length(pars)+1]) / delVec
  grad[is.na(grad)] <- 0 # set any NAs to grad 0

  return(grad)
}

###################################################################
### compute gradients of nll by finite diff without mc!...
###################################################################
gradnllIAK3D <- function(pars , zData , XData , vXU , iU , modelx , nud ,  
                            sdfdType_cd1 , sdfdType_cxd0 , sdfdType_cxd1 , cmeOpt , prodSum , setupMats , parBnds , useReml , lnTfmdData){
  
  if(length(pars) == 1){ stop('Error - do not use gradnllIAK3D with 1 parameter!') }else{}
  
  delVec <- rep(1E-6 , length(pars))

  parsMtx <- matrix(pars , length(pars) , length(pars) + 1)
  parsMtx[1:length(pars),1:length(pars)] <- parsMtx[1:length(pars),1:length(pars)] +  diag(delVec)
  parsList <- split(parsMtx, rep(1:ncol(parsMtx), each = nrow(parsMtx)))
  
  outtmp <- lapply(parsList , nllIAK3D , 
                       zData = zData , XData = XData , vXU = vXU , iU = iU , modelx = modelx , nud = nud ,  
                       sdfdType_cd1 = sdfdType_cd1 , sdfdType_cxd0 = sdfdType_cxd0 , sdfdType_cxd1 = sdfdType_cxd1 , cmeOpt = cmeOpt , 
                       prodSum = prodSum , setupMats = setupMats , parBnds = parBnds , useReml = useReml , lnTfmdData = lnTfmdData , rtnAll = F , forCompLik = FALSE)

  nllVec <- unlist(outtmp)
  
  grad <- (nllVec[1:length(pars)] - nllVec[length(pars)+1]) / delVec
  grad[is.na(grad)] <- 0 # set any NAs to grad 0
  
  return(grad)
}

gradnllIAK3D_CL <- function(pars , zData , XData , modelx , nud ,  
                               sdfdType_cd1 , sdfdType_cxd0 , sdfdType_cxd1 , cmeOpt , prodSum , setupMats , parBnds , useReml , compLikMats){
  
  if(length(pars) == 1){ stop('Error - do not use gradnllIAK3D.mc with 1 parameter!') }else{}
  
  delVec <- rep(1E-6 , length(pars))

  parsMtx <- matrix(pars , length(pars) , length(pars) + 1)
  parsMtx[1:length(pars),1:length(pars)] <- parsMtx[1:length(pars),1:length(pars)] +  diag(delVec)
  parsList <- split(parsMtx, rep(1:ncol(parsMtx), each = nrow(parsMtx)))
  
  outtmp <- lapply(parsList , nllIAK3D_CL , 
                       zData = zData , XData = XData , modelx = modelx , nud = nud ,  
                       sdfdType_cd1 = sdfdType_cd1 , sdfdType_cxd0 = sdfdType_cxd0 , sdfdType_cxd1 = sdfdType_cxd1 , cmeOpt = cmeOpt , 
                       prodSum = prodSum , setupMats = setupMats , parBnds = parBnds , useReml = useReml , compLikMats = compLikMats , rtnAll = F)

  nllVec <- unlist(outtmp)
  
  grad <- (nllVec[1:length(pars)] - nllVec[length(pars)+1]) / delVec
  grad[is.na(grad)] <- 0 # set any NAs to grad 0
  
  return(grad)
}

###############################################################################
### uses C = avCor * avsdfd1 * avsdfd2, 
### rather than the proper vs which uses C = av[Cor * sdfd1 * sdfd2]  
### allows other sdfd fns to be used that haven't had proper int avd cov fn worked out
### while should be quicker than discretisation approx
###############################################################################
setCApproxIAK3D <- function(parsBTfmd , modelx , 
                            sdfdType_cd1 , sdfdType_cxd0 , sdfdType_cxd1 , cmeOpt , setupMats){
  
  ### if setupMats doesn't include 'Kx', it should include 'xData' and 'dIData', so that proper setupMats can be made now...
  if(is.null(setupMats$Kx)){
    if(max(c(sdfdType_cd1 , sdfdType_cxd0 , sdfdType_cxd1)) > 0){ stop('Not ready to run setCIAK3D without passing in setupMats yet, because this requires passing in sdfdKnots - update the function!') }else{}
    setupMats <- setupIAK3D(xData = setupMats$xData , dIData = as.data.frame(setupMats$dIData) , nDscPts = 0 , 
                            sdfdType_cd1 = sdfdType_cd1 , sdfdType_cxd0 = sdfdType_cxd0 , sdfdType_cxd1 = sdfdType_cxd1) 
  }else{}
  
  n <- nrow(setupMats$Kx)
  C <- sparseMatrix(i=seq(n),j=seq(n),x=parsBTfmd$cme , symmetric = TRUE) # initiate with just meas err on diag. , dims = as.integer(c(n,n))
  
  ### adding Cx0...
  C <- C + parsBTfmd$cx0 * forceSymmetric(setupMats$Kx %*% t(setupMats$Kx))
  
  if (modelx == 'matern'){
    phix1 <- maternCov(setupMats$Dx , c(1 , parsBTfmd$ax , parsBTfmd$nux))
  }else if(modelx == 'nugget'){
    phix1 <- 0
  }else{
    stop('Not ready!')
  }
  
  ### get stationary cov mtx...
  tmp <- iaCovMatern(dIData = setupMats$dIU , ad = parsBTfmd$ad , nud = parsBTfmd$nud , sdfdPars = c(1 , 0 , 1) , sdfdType = 0 , abcd = setupMats$dIUabcd , iUElements = setupMats$dIUiUElements)
  phidStat <- tmp$avCovMtx
  avVardStat <- tmp$avVarVec
  remove(tmp) ;
  
  if(sdfdType_cd1 == 0){
    phid_cd1 <- phidStat
    avVard_cd1 <- avVardStat
    avsdfd_cd1 <- sqrt(avVard_cd1)
  }else if(sdfdType_cd1 == -9){
    phid_cd1 <- 0
    avVard_cd1 <- 0
    avsdfd_cd1 <- 0
  }else{
    
    avsdfd_cd1 <- iasdfd(dI = setupMats$dIU , sdfdPars = parsBTfmd$sdfdPars_cd1 , sdfdType = sdfdType_cd1 , XsdfdSpline = setupMats$XsdfdSplineU_cd1)
    avVard_cd1 <- avsdfd_cd1 ^ 2
    phid_cd1 <- phidStat
    phid_cd1@x <- phid_cd1@x * avsdfd_cd1[rep(seq(nrow(setupMats$dIU)) , times = seq(1,nrow(setupMats$dIU)))] * avsdfd_cd1[unlist(lapply(seq(nrow(setupMats$dIU)) , function(i){ seq(1,i) }))]
    
    # tmp <- iaCovMatern(dIData = setupMats$dIU , ad = parsBTfmd$ad , nud = parsBTfmd$nud , sdfdPars = parsBTfmd$sdfdPars_cd1 , sdfdType = sdfdType_cd1 , abcd = setupMats$dIUabcd , iUElements = setupMats$dIUiUElements)
    # phid_cd1 <- tmp$avCovMtx
    # avVard_cd1 <- tmp$avVarVec
    # remove(tmp) ; 
  }
  
  if(sdfdType_cxd0 == 0){
    phid_cxd0 <- phidStat
    avVard_cxd0 <- avVardStat
    avsdfd_cxd0 <- sqrt(avVard_cxd0)
  }else{
    
    avsdfd_cxd0 <- iasdfd(dI = setupMats$dIU , sdfdPars = parsBTfmd$sdfdPars_cxd0 , sdfdType = sdfdType_cxd0 , XsdfdSpline = setupMats$XsdfdSplineU_cxd0)
    avVard_cxd0 <- avsdfd_cxd0 ^ 2
    phid_cxd0 <- phidStat
    phid_cxd0@x <- phid_cxd0@x * avsdfd_cxd0[rep(seq(nrow(setupMats$dIU)) , times = seq(1,nrow(setupMats$dIU)))] * avsdfd_cxd0[unlist(lapply(seq(nrow(setupMats$dIU)) , function(i){ seq(1,i) }))]
    
    # tmp <- iaCovMatern(dIData = setupMats$dIU , ad = parsBTfmd$ad , nud = parsBTfmd$nud , sdfdPars = parsBTfmd$sdfdPars_cxd0 , sdfdType = sdfdType_cxd0 , abcd = setupMats$dIUabcd , iUElements = setupMats$dIUiUElements)
    # phid_cxd0 <- tmp$avCovMtx
    # avVard_cxd0 <- tmp$avVarVec
    # remove(tmp) ; 
    
  }
  
  if(modelx == 'nugget'){
    phid_cxd1 <- 0
    avVard_cxd1 <- 0
    avsdfd_cxd1 <- 0
  }else{
    if(sdfdType_cxd1 == 0){
      phid_cxd1 <- phidStat
      avVard_cxd1 <- avVardStat
      avsdfd_cxd1 <- sqrt(avVard_cxd1)
    }else{
      avsdfd_cxd1 <- iasdfd(dI = setupMats$dIU , sdfdPars = parsBTfmd$sdfdPars_cxd1 , sdfdType = sdfdType_cxd1 , XsdfdSpline = setupMats$XsdfdSplineU_cxd1)
      avVard_cxd1 <- avsdfd_cxd1 ^ 2
      phid_cxd1 <- phidStat
      phid_cxd1@x <- phid_cxd1@x * avsdfd_cxd1[rep(seq(nrow(setupMats$dIU)) , times = seq(1,nrow(setupMats$dIU)))] * avsdfd_cxd1[unlist(lapply(seq(nrow(setupMats$dIU)) , function(i){ seq(1,i) }))]
      
      # tmp <- iaCovMatern(dIData = setupMats$dIU , ad = parsBTfmd$ad , nud = parsBTfmd$nud , sdfdPars = parsBTfmd$sdfdPars_cxd1 , sdfdType = sdfdType_cxd1 , abcd = setupMats$dIUabcd , iUElements = setupMats$dIUiUElements)
      # phid_cxd1 <- tmp$avCovMtx
      # avVard_cxd1 <- tmp$avVarVec
      # remove(tmp) ; 
    }
  }
  
  parsOK <- T
  if((max(is.na(phid_cd1)) == 1) | (max(is.na(phid_cxd0)) == 1) | (max(is.na(phid_cxd1)) == 1)){ parsOK <- F } else{}
  if((min(avsdfd_cd1) < 0) | (min(avsdfd_cxd0) < 0) | (min(avsdfd_cxd1) < 0)){ parsOK <- FALSE }else{}
  
  if(parsOK){
    ### adding Cxd0...
    C <- C + parsBTfmd$cxd0 * sparseMatrix(i=setupMats$summKxKx$i , j=setupMats$summKxKx$j,x = phid_cxd0@x[setupMats$utriKdIdxdKd][setupMats$summKxKx$idxUtri],symmetric = TRUE)
    remove(phid_cxd0)
    
    ### without using diagBlocks fn...    
    KphidKStat <- new("dspMatrix" , Dim = as.integer(c(nrow(setupMats$Kd),nrow(setupMats$Kd))), x =  phidStat@x[setupMats$utriKdIdxdKd])
    
    ### adding Cd...
    if(sdfdType_cd1 == 0){
      C <- C + parsBTfmd$cd1 * KphidKStat
    }else if(sdfdType_cd1 == -9){
      Cd <- 0
    }else{
      C <- C + new("dspMatrix" , Dim = as.integer(c(nrow(setupMats$Kd),nrow(setupMats$Kd))), x =  phid_cd1@x[setupMats$utriKdIdxdKd])
    }
    remove(phid_cd1)
    
    if(modelx == 'nugget'){
      if(exists('KphidKStat')){ remove(KphidKStat) ;  }else{}
      Cx1 <- Cxd1 <- 0
    }else{
      if(sdfdType_cxd1 != 0){ # if not rqd, free up mem...
        if(exists('KphidKStat')){ remove(KphidKStat) ;  }else{}
      }else{}
      
      Kphix1K <- new("dspMatrix" , Dim = as.integer(c(nrow(setupMats$Kx),nrow(setupMats$Kx))), x =  phix1@x[setupMats$utriKxIdxxKx])
      rm(phix1)
      ### adding Cx1...
      C <- C + parsBTfmd$cx1 * Kphix1K
      
      ### adding Cxd1...
      if(sdfdType_cxd1 == 0){
        C <- C + parsBTfmd$cxd1 * Kphix1K * KphidKStat
      }else{
        if(exists('KphidKStat')){ remove(KphidKStat) ;  }else{}
        C <- C + parsBTfmd$cxd1 * Kphix1K * new("dspMatrix" , Dim = as.integer(c(nrow(setupMats$Kd),nrow(setupMats$Kd))), x =  phid_cxd1@x[setupMats$utriKdIdxdKd])
      }
      
      remove(Kphix1K) ; 
    }
    if(exists('KphidKStat')){ remove(KphidKStat) ;  }else{}
    if(exists('phid_cxd1')){ remove(phid_cxd1) ;  }else{}
    
    sigma2Vec <- parsBTfmd$cx0 + parsBTfmd$cx1 + parsBTfmd$cd1 * avVard_cd1 + parsBTfmd$cxd0 * avVard_cxd0 + parsBTfmd$cxd1 * avVard_cxd1 + parsBTfmd$cme
    sigma2Vec <- setupMats$Kd %*%sigma2Vec 
    
  }else{        
    C <- NA
    sigma2Vec <- NA
#    print('Bad parameters...not sure how this can happen')
#    print(parsBTfmd)
  }
  
  return(list('C' = C , 'sigma2Vec' = sigma2Vec))
  #	return(list('C' = C , 'sigma2Vec' = sigma2Vec , 'dC' = list(Cx0 , Cx1 , Cd , Cxd0 , Cxd1) , 'phi' = list(phix0 , phix1 , phid_cd1 , phid_cxd0 , phid_cxd1)))
}

################################################################
### use setCApproxIAK3D2 for covs between different sets of data (setupMats will also contain xU2, dIU2, etc)
### uses C = avCor * avsdfd1 * avsdfd2, 
### rather than the proper vs which uses C = av[Cor * sdfd1 * sdfd2]  
### allows other sdfd fns to be used that haven't had proper int avd cov fn worked out
### while should be quicker than discretisation approx
################################################################
setCApproxIAK3D2 <- function(parsBTfmd , modelx , 
                              sdfdType_cd1 , sdfdType_cxd0 , sdfdType_cxd1 , cmeOpt , setupMats){
  
  ### if setupMats doesn't include 'Kx', it should include 'xData', 'dIData', 'xData2' and 'dIData2', so that proper setupMats can be made now...
  if(is.null(setupMats$Kx)){
    if(max(c(sdfdType_cd1 , sdfdType_cxd0 , sdfdType_cxd1)) > 0){ stop('Not ready to run setCIAK3D without passing in setupMats yet, because this requires passing in sdfdKnots - update the function!') }else{}
    setupMats <- setupIAK3D2(xData = setupMats$xData , dIData = setupMats$dIData , xData2 = setupMats$xData2 , dIData2 = setupMats$dIData2 , 
                            sdfdType_cd1 = sdfdType_cd1 , sdfdType_cxd0 = sdfdType_cxd0 , sdfdType_cxd1 = sdfdType_cxd1) 
  }else{}
  
  ### starting with Cx0...
  ### defined for the unique locations...
  ijTmp <- which(setupMats$Dx == 0 , arr.ind = TRUE)
  if(nrow(ijTmp) > 0){
    phix0 <- sparseMatrix(i = ijTmp[,1] , j = ijTmp[,2] , x = 1 , dims = c(nrow(setupMats$Dx) , ncol(setupMats$Dx)))
    C <- parsBTfmd$cx0 * setupMats$Kx %*% phix0 %*% t(setupMats$Kx2)
  }else{
    phix0 <- sparseMatrix(i = 1 , ,j=1 , x=0 , dims = c(nrow(setupMats$Dx),ncol(setupMats$Dx)))
    C <- sparseMatrix(i=1,j=1,x=0,dims = c(nrow(setupMats$Kx),nrow(setupMats$Kx2)))
  }
  
  if (modelx == 'matern'){
    phix1 <- maternCov(setupMats$Dx , c(1 , parsBTfmd$ax , parsBTfmd$nux))
  }else if(modelx == 'nugget'){
    phix1 <- 0
  }else{
    stop('Not ready!')
  }
  
  ### get stationary cov mtx...
  phidStat <- iaCovMatern2(dIData = setupMats$dIU , dIData2 = setupMats$dIU2 , ad = parsBTfmd$ad , nud = parsBTfmd$nud , sdfdPars = c(1 , 0 , 1) , sdfdType = 0)
  
  if(sdfdType_cd1 == 0){
    phid_cd1 <- phidStat
    avsdfd_cd11 <- avsdfd_cd12 <- 1
  }else if(sdfdType_cd1 == -9){
    phid_cd1 <- 0
    avsdfd_cd11 <- avsdfd_cd12 <- 0
  }else{
    
    avsdfd_cd11 <- iasdfd(dI = setupMats$dIU , sdfdPars = parsBTfmd$sdfdPars_cd1 , sdfdType = sdfdType_cd1 , XsdfdSpline = setupMats$XsdfdSplineU_cd1)
    avsdfd_cd12 <- iasdfd(dI = setupMats$dIU2 , sdfdPars = parsBTfmd$sdfdPars_cd1 , sdfdType = sdfdType_cd1 , XsdfdSpline = setupMats$XsdfdSplineU_cd12)
    phid_cd1 <- phidStat * matrix(avsdfd_cd11 , length(avsdfd_cd11) , length(avsdfd_cd12)) * matrix(avsdfd_cd12 , length(avsdfd_cd11) , length(avsdfd_cd12) , byrow = T) 
    
    # phid_cd1 <- iaCovMatern2(dIData = setupMats$dIU , dIData2 = setupMats$dIU2 , ad = parsBTfmd$ad , nud = parsBTfmd$nud , sdfdPars = parsBTfmd$sdfdPars_cd1 , sdfdType = sdfdType_cd1)
  }
  
  if(sdfdType_cxd0 == 0){
    phid_cxd0 <- phidStat
    avsdfd_cxd01 <- avsdfd_cxd02 <- 1
  }else{
    avsdfd_cxd01 <- iasdfd(dI = setupMats$dIU , sdfdPars = parsBTfmd$sdfdPars_cxd0 , sdfdType = sdfdType_cxd0 , XsdfdSpline = setupMats$XsdfdSplineU_cxd0)
    avsdfd_cxd02 <- iasdfd(dI = setupMats$dIU2 , sdfdPars = parsBTfmd$sdfdPars_cxd0 , sdfdType = sdfdType_cxd0 , XsdfdSpline = setupMats$XsdfdSplineU_cxd02)
    phid_cxd0 <- phidStat * matrix(avsdfd_cxd01 , length(avsdfd_cxd01) , length(avsdfd_cxd02)) * matrix(avsdfd_cxd02 , length(avsdfd_cxd01) , length(avsdfd_cxd02) , byrow = T) 
    
    # phid_cxd0 <- iaCovMatern2(dIData = setupMats$dIU , dIData2 = setupMats$dIU2 , ad = parsBTfmd$ad , nud = parsBTfmd$nud , sdfdPars = parsBTfmd$sdfdPars_cxd0 , sdfdType = sdfdType_cxd0)
  }
  
  if(modelx == 'nugget'){
    phid_cxd1 <- 0
    avsdfd_cxd11 <- avsdfd_cxd12 <- 0
  }else{
    if(sdfdType_cxd1 == 0){
      phid_cxd1 <- phidStat
      avsdfd_cxd11 <- avsdfd_cxd12 <- 1
    }else{
      avsdfd_cxd11 <- iasdfd(dI = setupMats$dIU , sdfdPars = parsBTfmd$sdfdPars_cxd1 , sdfdType = sdfdType_cxd1 , XsdfdSpline = setupMats$XsdfdSplineU_cxd1)
      avsdfd_cxd12 <- iasdfd(dI = setupMats$dIU2 , sdfdPars = parsBTfmd$sdfdPars_cxd1 , sdfdType = sdfdType_cxd1 , XsdfdSpline = setupMats$XsdfdSplineU_cxd12)
      phid_cxd1 <- phidStat * matrix(avsdfd_cxd11 , length(avsdfd_cxd11) , length(avsdfd_cxd12)) * matrix(avsdfd_cxd12 , length(avsdfd_cxd11) , length(avsdfd_cxd12) , byrow = T) 
      
      # phid_cxd1 <- iaCovMatern2(dIData = setupMats$dIU , dIData2 = setupMats$dIU2 , ad = parsBTfmd$ad , nud = parsBTfmd$nud , sdfdPars = parsBTfmd$sdfdPars_cxd1 , sdfdType = sdfdType_cxd1)
    }
  }
  
  parsOK <- T
  if((max(is.na(phid_cd1)) == 1) | (max(is.na(phid_cxd0)) == 1) | (max(is.na(phid_cxd1)) == 1)){ parsOK <- F } else{}
  if(min(c(avsdfd_cd11 , avsdfd_cd12 , avsdfd_cxd01 , avsdfd_cxd02 , avsdfd_cxd11 , avsdfd_cxd12)) < 0){ parsOK <- F } else{}
  
  if(parsOK){
    
    ### adding Cxd0...    
    if(modelx == 'nugget' & sdfdType_cd1 == -9){
      KphidKStat <- (setupMats$Kx %*% phix0 %*% t(setupMats$Kx2)) * (setupMats$Kd %*% phidStat %*% t(setupMats$Kd2))
      C <- C + parsBTfmd$cxd0 * (setupMats$Kx %*% phix0 %*% t(setupMats$Kx2)) * (setupMats$Kd %*% phid_cxd0 %*% t(setupMats$Kd2))
    }else{        
      KphidKStat <- setupMats$Kd %*% phidStat %*% t(setupMats$Kd2)	
      C <- C + parsBTfmd$cxd0 * (setupMats$Kx %*% phix0 %*% t(setupMats$Kx2)) * (setupMats$Kd %*% phid_cxd0 %*% t(setupMats$Kd2))
    }
    
    ### adding Cd...    
    if(sdfdType_cd1 == 0){
      C <- C + parsBTfmd$cd1 * KphidKStat
    }else if(sdfdType_cd1 == -9){
      Cd <- 0
    }else{
      C <- C + parsBTfmd$cd1 * (setupMats$Kd %*% phid_cd1 %*% t(setupMats$Kd2))
    }
    
    if(modelx == 'nugget'){
      Cx1 <- Cxd1 <- 0
    }else{
      Kphix1K <- setupMats$Kx %*% phix1 %*% t(setupMats$Kx2)
      ### adding Cx1...    
      C <- C + parsBTfmd$cx1 * Kphix1K 
      ### adding Cxd1...    
      if(sdfdType_cxd1 == 0){
        C <- C + parsBTfmd$cxd1 * Kphix1K * KphidKStat 
      }else{
        C <- C + parsBTfmd$cxd1 * Kphix1K * (setupMats$Kd %*% phid_cxd1 %*% t(setupMats$Kd2))
      }
      remove(Kphix1K) ; 
    }
    if(exists('KphidKStat')){ remove(KphidKStat) ;  }else{}
    
    ### note, no meas err being added to cov between distinct data...
    
  }else{        
    C <- NA
#    print('Bad parameters...not sure how this can happen')
#    print(parsBTfmd)
  }
  
  return(C)
}

