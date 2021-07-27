##############################################################
### source('/export/home/ortont/scripts/iakTests/egEdgeroi.R')
### source('U://scripts/iakTests/egEdgeroi.R')
##############################################################
rm(list=ls())
assign("last.warning", NULL, envir = baseenv())

##############################################################
### full model selection took ~ 19 hours.
##############################################################
fitCubistModelNow <- TRUE
fitModelNow <- TRUE
plotVargiogramFit <- TRUE
valNow <- TRUE
val4PlotNow <- TRUE
mapNow <- FALSE 

##############################################################
### load libraries and source files...
##############################################################
library(sp)
library(raster)
if(Sys.info()[1] != 'Windows'){ rasterOptions(tmpdir="/scratch/tmp/") }else{}
library(rgdal)
library(Cubist)
library(mgcv)

library(Matrix)
library(MASS)
library(splines)
library(deldir)
library(lme4)

library(aqp)
library(GSIF)
library(ithir)
library(parallel)


if(Sys.info()[1] == 'Windows'){
  wDir <- 'U://scripts/'
  iakDir <- 'U://scripts/iakTests'
  lmm2Dir <- 'U://scripts/fLMM2'
  dataDir <- 'Z://ortont/Data/edgeroiEg'
}else{
  wDir <- '/export/home/ortont/scripts/'
  iakDir <- '/export/home/ortont/scripts/iakTests'
  lmm2Dir <- '/export/home/ortont/scripts/fLMM2'
  dataDir <- '/scratch/rsc2/ortont/Data/soilConstraints2/soilData/ecMap2000s'
}

setwd(wDir)
source(paste0(iakDir , '/fitIAK3D.R'))
source(paste0(iakDir , '/setInitsIAK3D.R'))
source(paste0(iakDir , '/iaCovMatern.R'))
source(paste0(iakDir , '/nrUpdatesIAK3D.R'))
source(paste0(iakDir , '/nrUpdatesIAK3DlnN.R'))
source(paste0(iakDir , '/makeXvX.R'))
source(paste0(iakDir , '/optifix.R'))
source(paste0(iakDir , '/cubist2XIAK3D.R')) 
source(paste0(iakDir , '/gam2XIAK3D.R')) 
source(paste0(iakDir , '/predictIAK3D.R')) 
source(paste0(iakDir , '/modelSelectIAK3D.R')) 
source(paste0(iakDir , '/compositeLikIAK3D.R')) 
source(paste0(iakDir , '/mpspline.source.R')) 
source(paste0(iakDir , '/getEdgeroiData.R')) 
source(paste0(lmm2Dir , '/fLMM2.R')) 

source(paste0(iakDir , '/splineBasisFns.R'))
source(paste0(iakDir , '/makeXvX_gam2.R'))

printnllTime <<- FALSE

crsAusAlbers <- CRS("+proj=aea +lat_1=-18 +lat_2=-36 +lat_0=0 +lon_0=132 +x_0=0 +y_0=0 +ellps=GRS80 +datum=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
crsLongLat <- CRS('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')

##########################################################
### set some options for iak...
##########################################################
modelX <- 'cubist' # 'gam2' or 'cubist'. For linear model and how to code, see top of makeXvX.R - e.g. 1 for a model with first order d/spat cov interactions
constrainX4Pred <- FALSE # use the limits of each column of X (non-zeros only) from the fitting data to constrain the columns of X for prediction?
prodSum <- TRUE # product sum model? If FALSE, product model is used.  
modelx <- 'matern' # matern (range = a / 3) or wendland or spherical or nugget. wendland or spherical (range = a) have C = 0 beyond range so can be computable with bigger datasets. 
minRange <- NA # NA will give default minRange (1.5 * min sep dist)
maxRange <- NA # NA will give default maxRange (0.5 * max sep dist)

# nud <- NULL # use this if running the modelSelectIAK3D function.
nud <- 0.5
lnTfmdData <- FALSE 
rqrBTfmdPreds <- FALSE # only relevant if data were log-transformed - back transforms via exp(z + 0.5 * v), ie to minimize expected sqd err
useReml <- TRUE

allKnotsd <- c() # can use this to give spline of depth in combination with using cubist for trend. if gam2 for trend keep empty
if(length(allKnotsd) > 0){ incdSpline <- TRUE }else{ incdSpline <- FALSE }

# sdfdTypeANDcmeInit <- c(0 , -1 , -1 , 1) # stationary d component in prodSum model (0), non-stat cxd0 with exp fn (-1) and cxd1 with exp fn (-1), meas err included (1)
sdfdTypeANDcmeInit <- c(0 , 2 , 2 , 1) # stationary d component in prodSum model (0), non-stat cxd0 with spline fn (2 = 2 extra pars fitted per non-stat fn of depth = same complexity as -1) and cxd1 with spline fn (2), meas err included (1)
# sdfdTypeANDcmeInit <- c(0 , 0 , 0 , 1) # stationary d component in prodSum model (0), stat cxd0 (0) and cxd1 (0), meas err included (1)
if((modelX == 'gam2') | (length(allKnotsd) > 0)){
  sdfdTypeANDcmeInit[1] <- -9 # no d component in prodSum model because spline used instead (-9)
}else{}

if(identical(modelX , 'cubist')){
  ### nRules : number of rules for Cubist model. If NA, then a x val routine to select nRules will be used. 
  ### refineCubistModel : use a stepwise algorithm (refineXIAK3D in cubist2XIAK3D.R) to remove predictors from the fitted Cubist model? # If NA, a x val routine will be used to select T/F
  optionsModelX <- list('refineCubistModel' = TRUE , 'nRules' = 5 , 
                        'reduceXAfterInitFit' = FALSE)
}else if(identical(modelX , 'gam2')){
  optionsModelX <- list('nIntKnotsd' = 1 , 'nIntKnotss' = 2 , 'incInts' = FALSE , 'intMthd' = 0 ,
                        'reduceXAfterInitFit' = FALSE , 'optStat4Drop' = 1 , 'alpha' = NA ,
                        'maxnSpatVars' = 3 , 'maxnKnotsd' = 2 , 'maxnKnotss' = 2 ,
                        'q4BdryKnotsd' = c(0.05 , 0.95) , 'q4BdryKnotss' = c(-0.2 , 1.2))
  
  # ### Inf for bdry knots + 2 internal knots = linear function. See makeXnsClampedLUG0
  # optionsModelX <- list('nIntKnotsd' = 2 , 'nIntKnotss' = 2 , 'incInts' = FALSE , 'intMthd' = 0 ,
  #                       'reduceXAfterInitFit' = FALSE , 'optStat4Drop' = 1 , 'alpha' = NA ,
  #                       'maxnSpatVars' = 3 , 'maxnKnotsd' = 2 , 'maxnKnotss' = 2 ,
  #                       'q4BdryKnotsd' = c(0.05 , 0.95) , 'q4BdryKnotss' = c(-Inf , Inf)) # -Inf , Inf with maxnKnotss = 2 gives linear fns of spatial covariates
  
  ### or with interactions between spatial covs and depth...
  # optionsModelX <- list('nIntKnotsd' = 2 , 'nIntKnotss' = 2 , 'incInts' = list('dIMidPts' , setdiff(covNames4IAK3D , 'dIMidPts')) , 'intMthd' = 0 ,
  #                       'reduceXAfterInitFit' = FALSE , 'optStat4Drop' = 1 , 'alpha' = NA ,
  #                       'maxnSpatVars' = 3 , 'maxnKnotsd' = 2 , 'maxnKnotss' = 2 , 
  #                       'q4BdryKnotsd' = c(0.05 , 0.95) , 'q4BdryKnotss' = c(-0.2 , 1.2))
  
}else if(is.numeric(modelX) | (length(modelX) > 1 & is.character(modelX))){ # assume this will be linear model. see top of makeXvX for how to specify
  # this needs checking as hasn't been used for a while
  covsFit$dIMidPts <- NULL
  covsVal$dIMidPts <- NULL
  optionsModelX <- list()
}else{
  stop('Check code is still working for this type of modelX!')
}

useCL <- FALSE # use the composite likelihood approx?

##############################################
### load the edgeroi dataset (from GSIF package) and put into format for iak3d...
##############################################
tmp <- getEdgeroiData()
cFit <- tmp$cFit
dIFit <- tmp$dIFit
covsFit <- tmp$covsFit
zFit <- tmp$zFit
profIDFit <- tmp$profIDFit

cVal <- tmp$cVal
dIVal <- tmp$dIVal
covsVal <- tmp$covsVal
zVal <- tmp$zVal
profIDVal <- tmp$profIDVal

rList <- tmp$rList

### to do a temp test with a cat var...
# covsFit$radK_TMP <- NA
# covsFit$radK_TMP[covsFit$radK < 0.5] <- 'K_1'
# covsFit$radK_TMP[covsFit$radK >= 0.5 & covsFit$radK < 1] <- 'K_2'
# covsFit$radK_TMP[covsFit$radK >= 1 & covsFit$radK < 1.5] <- 'K_3'
# covsFit$radK_TMP[covsFit$radK >= 1.5] <- 'K_4'
# covsFit$radK <- factor(covsFit$radK_TMP , levels = paste0('K_' , seq(4)))
# covsFit$radK_TMP <- NULL
# 
# covsVal$radK_TMP <- NA
# covsVal$radK_TMP[covsVal$radK < 0.5] <- 'K_1'
# covsVal$radK_TMP[covsVal$radK >= 0.5 & covsVal$radK < 1] <- 'K_2'
# covsVal$radK_TMP[covsVal$radK >= 1 & covsVal$radK < 1.5] <- 'K_3'
# covsVal$radK_TMP[covsVal$radK >= 1.5] <- 'K_4'
# covsVal$radK <- factor(covsVal$radK_TMP , levels = paste0('K_' , seq(4)))
# covsVal$radK_TMP <- NULL

##################################################
### define (or load) the matrices for composite likelihood
### for looking at different models with same dataset, make sure these matrices are always the same (ie done once then reloaded)
##################################################
if(useCL){
  ### if fitting different models, save out the compLikMats so that it is the same each time.
  #  compLikMats <- setVoronoiBlocksWrap(x = cFit , nPerBlock = 50 , plotVor = T , optnBalance = 2)
  ### also attach compLikOptn
  #  compLikMats$compLikOptn <- 2
  ## save(compLikMats , file = paste0(dataDir , '/compLikMats.RData'))
  # saveRDS(compLikMats , file = paste0(dataDir , '/compLikMats.RData'))
  
  ## load(file = paste0(dataDir , '/compLikMats.RData'))
  
  compLikMats <- readRDS(file = paste0(dataDir , '/compLikMats.RData'))
  
}else{
  compLikMats <- list('compLikOptn' = 0) # use exact likelihood
}

#################################################################################################
### set knots for sdfd spline fn (if used)
#################################################################################################
sdfdKnots <- setKnots4sdfd(dIFit , sdfdType_cd1 = sdfdTypeANDcmeInit[1] , sdfdType_cxd0 = sdfdTypeANDcmeInit[2] , sdfdType_cxd1 = sdfdTypeANDcmeInit[3])


lmmFitFile <- paste0(dataDir , '/lmm.fit.selected.RData') # for the fitted model
nmplt <- paste0(dataDir , '/plot.selected.gam2.pdf') # for a plot with the internal 'predictions' = predictions through profiles of sampled profiles (not validation, can be a check of what's going on)

if(fitModelNow){

### this bit can be run to do the model selection - took around 19 hours. Could be set up to be done in parallel though.
  # start_time <- Sys.time()
  # modSelectOutput <- selectCovIAK3D(xData = cFit , dIData = dIFit , zData = zFit , covsData = covsFit , modelX = cmFit , modelx = 'matern' , nud = nud , 
  #                                   sdfdTypeANDcmeInit = sdfdTypeANDcmeInit , allKnotsd = allKnotsd , prodSum = prodSum , lnTfmdData = lnTfmdData , 
  #                                   useReml = useReml , compLikMats = compLikMats , rqrBTfmdPreds = rqrBTfmdPreds , dirPlot = dataDir)
  # sdfdTypeANDcmeInit <- modSelectOutput$sdfdTypeANDcmeOptSelected
  # nud <- modSelectOutput$nud
  #   
  # end_time <- Sys.time()
  # print('Time to fit was:')
  # print(end_time - start_time)
  # 
  # save(modSelectOutput , file = paste0(dataDir , '/modSelectOutput.RData'))
  # 
  # lmm.fit.selected <- modSelectOutput$lmmSelected
  # 
  # save(lmm.fit.selected , file = lmmFitFile)
  # 
  # stop('What model was selected?')
  #     
  
### refit cubist model as lmm...if selectCovIAK3D was run don't need to do this bit
  start_time <- Sys.time()
  
  tmp <- fitIAK3D(xData = cFit , dIData = dIFit , zData = zFit , covsData = covsFit , modelX = modelX , modelx = modelx , nud = nud , allKnotsd = allKnotsd , 
                  sdfdType_cd1 = sdfdTypeANDcmeInit[1] , sdfdType_cxd0 = sdfdTypeANDcmeInit[2] , sdfdType_cxd1 = sdfdTypeANDcmeInit[3] , 
                  cmeOpt = sdfdTypeANDcmeInit[4] , sdfdKnots = sdfdKnots , minRange = minRange , maxRange = maxRange , 
                  prodSum = prodSum , lnTfmdData = lnTfmdData , useReml = useReml , 
                  optionsModelX = optionsModelX , compLikMats = compLikMats , namePlot = nmplt , rqrBTfmdPreds = rqrBTfmdPreds) 

  end_time <- Sys.time()
  print('Time to fit was:')
  print(end_time - start_time)

  lmm.fit.selected <- tmp$lmmFit

  save(lmm.fit.selected , file = lmmFitFile)

}else{
  load(file = lmmFitFile)
}

stop('done')

###########################################################################
### some plots of the fitted covariance model...
###########################################################################
if(plotVargiogramFit){
  dIPlot <- data.frame('dU' = c(0 , 20 , 50 , 90 , 150 , 190)/100 , 'dL' = c(10 , 30 , 60 , 100 , 160 , 200)/100)
  hx <- seq(0 , 20 , 1)
  pdf(file = paste0(dataDir , '/varioFitgam22.pdf'))
  tmp <- plotCovx(lmm.fit = lmm.fit.selected , hx = hx , dIPlot = dIPlot , addExpmntlV = TRUE , hzntlUnits = 'km')
  dev.off()
  
  hdPlot <- seq(0 , 2 , 0.01)
  pdf(file = paste0(dataDir , '/cordFit.pdf'))
  qwe <- plotCord(lmm.fit = lmm.fit.selected , hdPlot = hdPlot, vrtclUnits = 'm')
  dev.off()
  
  dTmp <- seq(0 , 2 , 0.1)
  dIPlot <- data.frame('dU' = dTmp[-length(dTmp)] , 'dL' = dTmp[-1])
  pdf(file = paste0(dataDir , '/covardFit.pdf'))
  qwe <- plotCovd(lmm.fit = lmm.fit.selected , dIPlot = dIPlot , vrtclUnits = 'm')
  dev.off()
  
  ### plot of the variances...
  pdf(file = paste0(dataDir , '/varComps.pdf'))
  dPlot <- seq(0 , 2 , 0.01)
  plotVarComps(lmm.fit = lmm.fit.selected , dPlot = dPlot)
  dev.off()

}else{}

#########################################################
### validation bit...
#########################################################
fnamezkVal <- paste0(dataDir , '/zkVal.RData')
fnamevkVal <- paste0(dataDir , '/vkVal.RData')
namePlot = paste0(dataDir , '/plotVal.pdf')

if(valNow){

  nVal <- nrow(cVal)
  iU <- which(!duplicated(cVal))
  cValU <- cVal[iU,,drop=FALSE]
  covsValU <- covsVal[iU,,drop=FALSE]
  zkVal <- vkVal <- NA * numeric(nVal)
  for(i in 1:nrow(cValU)){
      iTmp <- which(cVal[,1] == cValU[i,1] & cVal[,2] == cValU[i,2])
      tmp <- profilePredictIAK3D(xMap = cValU[i,,drop=FALSE] , covsMap = covsVal[iTmp,,drop=FALSE] , dIMap = dIVal[iTmp,,drop=FALSE] , lmmFit = lmm.fit.selected , rqrBTfmdPreds = rqrBTfmdPreds , constrainX4Pred = constrainX4Pred)
      zkVal[iTmp] <- tmp$zMap 
      vkVal[iTmp] <- tmp$vMap 
  }

  dIPred <- cbind(seq(0 , 1.98 , 0.02) , seq(0.02 , 2 , 0.02))
  tmp <- predictIAK3D(xMap = cValU , dIMap = dIPred , covsMap = covsValU , lmmFit = lmm.fit.selected , rqrBTfmdPreds = rqrBTfmdPreds , constrainX4Pred = constrainX4Pred)

  zkProfPred <- tmp$zMap
  vkProfPred <- tmp$vMap
  pi90LkProfPred <- tmp$pi90LMap
  pi90UkProfPred <- tmp$pi90UMap

  ### calc and print val stats...
  tmp <- calcValStats(zVal = zVal , dIVal = dIVal , zkVal = zkVal , vkVal = vkVal , layerMidPts = c(0.025 , 0.1 , 0.225 , 0.45 , 0.8 , 1.5) , printValStats = TRUE)
  valStatsAllLayers <- tmp$valStatsAllLayers 
  valStatsTot <- tmp$valStatsTot
  
  tmp <- plotProfilesIAK3D(namePlot = namePlot , xData = cVal , dIData = dIVal , zData = zVal , 
                  xPred = cValU , dIPred = dIPred , zPred = zkProfPred , pi90LPred = pi90LkProfPred , pi90UPred = pi90UkProfPred , 
                  zhatxv = zkVal , pi90Lxv = zkVal - 1.64 * sqrt(vkVal) , pi90Uxv = zkVal + 1.64 * sqrt(vkVal)) 

  save(zkVal , file = fnamezkVal)
  save(vkVal , file = fnamevkVal)

}else{

  load(file = fnamezkVal)
  load(file = fnamevkVal)

}
 
### keeping this bit separate. 
if(val4PlotNow){

#> sample(60 , 6)
# [1]  6 19 49 41  3 24

  rand6ForPlot <- c(6 , 19 , 49 , 41 , 3 , 24)

  i4PlotU <- which(!duplicated(cVal))[rand6ForPlot]
  cVal4PlotU <- cVal[i4PlotU,,drop=FALSE]
  covsVal4PlotU <- covsVal[i4PlotU,,drop=FALSE]

  iVal4Plot <- c()
  for(i in 1:nrow(cVal4PlotU)){
      iTmp <- which(cVal[,1] == cVal4PlotU[i,1] & cVal[,2] == cVal4PlotU[i,2])
      iVal4Plot <- c(iVal4Plot , iTmp)
  }

  cVal4Plot <- cVal[iVal4Plot,,drop=FALSE]
  dIVal4Plot <- dIVal[iVal4Plot,,drop=FALSE]
  covsVal4Plot <- covsVal[iVal4Plot,,drop=FALSE]
  zVal4Plot <- zVal[iVal4Plot]

  zkVal4Plot <- vkVal4Plot <- NA * numeric(length(zVal4Plot))
  for(i in 1:nrow(cVal4PlotU)){
      iTmp <- which(cVal4Plot[,1] == cVal4PlotU[i,1] & cVal4Plot[,2] == cVal4PlotU[i,2])
      tmp <- profilePredictIAK3D(xMap = cVal4PlotU[i,,drop=FALSE] , covsMap = covsVal4Plot[iTmp,,drop=FALSE] , dIMap = dIVal4Plot[iTmp,,drop=FALSE] , lmmFit = lmm.fit.selected , rqrBTfmdPreds = rqrBTfmdPreds , constrainX4Pred = constrainX4Pred)
      zkVal4Plot[iTmp] <- tmp$zMap 
      vkVal4Plot[iTmp] <- tmp$vMap 
  }

  dIPred <- cbind(seq(0 , 1.98 , 0.02) , seq(0.02 , 2 , 0.02))
  tmp <- predictIAK3D(xMap = cVal4PlotU , dIMap = dIPred , covsMap = covsVal4PlotU , lmmFit = lmm.fit.selected , rqrBTfmdPreds = rqrBTfmdPreds , constrainX4Pred = constrainX4Pred)

  zkProfPred <- tmp$zMap
  vkProfPred <- tmp$vMap
  pi90LkProfPred <- tmp$pi90LMap
  pi90UkProfPred <- tmp$pi90UMap

  zVal4Plot_PLOT <- zVal4Plot

  zkVal4Plot_PLOT <- zkVal4Plot
  pi90LkVal4Plot_PLOT <- zkVal4Plot - 1.64 * sqrt(vkVal4Plot)
  pi90UkVal4Plot_PLOT <- zkVal4Plot + 1.64 * sqrt(vkVal4Plot)

  zkProfPred_PLOT <- zkProfPred
  pi90LkProfPred_PLOT <- pi90LkProfPred
  pi90UkProfPred_PLOT <- pi90UkProfPred

  xlab <- "Clay percent"
  vecTmp <- c(zVal4Plot_PLOT , zkVal4Plot_PLOT , as.numeric(zkProfPred_PLOT))
  xlim <- c(min(vecTmp) , max(vecTmp))

  namePlot = paste0(dataDir , '/plotVal4Plot.pdf')

  tmp <- plotProfilesIAK3D(namePlot = namePlot , xData = cVal4Plot , dIData = dIVal4Plot , zData = zVal4Plot_PLOT , 
                  xPred = cVal4PlotU , dIPred = dIPred , zPred = zkProfPred_PLOT , pi90LPred = pi90LkProfPred_PLOT , pi90UPred = pi90UkProfPred_PLOT , 
                  zhatxv = zkVal4Plot_PLOT , pi90Lxv = pi90LkVal4Plot_PLOT , pi90Uxv = pi90UkVal4Plot_PLOT , 
                  profNames = paste0('Profile ' , rand6ForPlot) , xlim = xlim , xlab = xlab) 

}else{}
### plot pred vs obs for 6 gsm layers...

stop('Fitting and validation done!')

########################################################################
### set up covaraites for mapping, based on the grid that all covariates are on in rList...
########################################################################
if(mapNow){

#  dIMap <- data.frame('dU' = c(0 , 0.05 , 0.15 , 0.3 , 0.6 , 1.0) , 'dL' = c(0.05 , 0.15 , 0.3 , 0.6 , 1.0 , 2.0))
  dIMap <- data.frame('dU' = c(0 , 0.15 , 0.6) , 'dL' = c(0.05 , 0.3 , 1.0))

  firstRow <- 1
  lastRow <- 10
  
### get which rows are needed in this batch...
  rowsToDo <- firstRow:lastRow

  for(irow in rowsToDo){

### cell centres for this row...
    xVecMap <- seq(xFromCol(rList[[1]] , 1) , xFromCol(rList[[1]] , ncol(rList[[1]])) , res(rList[[1]])[1])
    yVecMap <- yFromRow(rList[[1]] , irow)
    
### get coordinates and covariates for this row...
    cMap <- data.frame('Eastings' = xVecMap , 'Northings' = yVecMap)
    
### define cMap and extract covsMap for this row...
    covsMap <- data.frame(matrix(NA , ncol(rList[[1]]) , length(rList)))
    for (icov in 1:length(rList)){ 
      covsMap[,icov] <- extract(rList[[icov]] , cMap) 
    }
    names(covsMap) <- names(rList)
    iIn <- which(!is.na(rowSums((covsMap))))

    covsMap[['dIMidPts']] <- NA
    
################################################
### calculate the predictions for the mapping depth...
################################################
    lmm.map.tmp <- predictIAK3D(xMap = cMap[iIn,,drop=FALSE] , covsMap = covsMap[iIn,,drop=FALSE] , dIMap = dIMap , lmmFit = lmm.fit.selected , rqrBTfmdPreds = rqrBTfmdPreds , constrainX4Pred = constrainX4Pred)
    lmm.map <- list()
    for (j in 1:length(lmm.map.tmp)){
      lmm.map[[names(lmm.map.tmp)[j]]] <- matrix(NA , nrow(dIMap) , ncol(rList[[1]]))
      lmm.map[[names(lmm.map.tmp)[j]]][,iIn] <- lmm.map.tmp[[names(lmm.map.tmp)[j]]]
    }

    save(lmm.map , file = paste0(dataDir , '/map.row' , irow , '.RData'))

  }
  
  print(paste0('Mapped for row ' , irow))
  
} ### done

stop(paste0('Mapped for batch ' , iBatch))
