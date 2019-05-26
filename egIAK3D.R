##############################################################
### source('/export/home/ortont/scripts/soilConstraints2/mapSoilConstraints.R')
##############################################################
rm(list=ls())
assign("last.warning", NULL, envir = baseenv())


testCL <- FALSE

nRules <- 10

mapRegion <- 'NGR'

############################################################
### get the info from the bq arguments.
### the following must agree with what is in the loop for sending to bq.
############################################################
#  13917 x 7457
nrowsMap <- 13917

nRowsPerBatch <- 100
nBatches <- ceiling(nrowsMap / nRowsPerBatch)

##############################################################
### load libraries and source files...
##############################################################
library(sp)
library(raster)
rasterOptions(tmpdir="/scratch/tmp/")
library(rgdal)
library(Cubist)
library(mgcv)

library(Matrix)
library(MASS)
library(splines)
library(deldir)

wDir <- '/export/home/ortont/scripts/'
iakDir <- '/export/home/ortont/scripts/iakTests'
lmm2Dir <- '/export/home/ortont/scripts/fLMM2'

dataDir <- '/scratch/rsc6/ortont/Data/soilConstraints2/soilData/ecMap2000s'

setwd(wDir)
source(paste0(iakDir , '/fitIAK3D.R'))
source(paste0(iakDir , '/setInitsIAK3D.R'))
source(paste0(iakDir , '/iaCovMatern.R'))
source(paste0(iakDir , '/nrUpdatesIAK3D.R'))
source(paste0(iakDir , '/nrUpdatesIAK3DlnN.R'))
source(paste0(iakDir , '/makeXvX.R'))
source(paste0(iakDir , '/optifix.R'))
source(paste0(iakDir , '/cubist2XIAK3D.R')) 
source(paste0(iakDir , '/predictIAK3D.R')) 
source(paste0(iakDir , '/modelSelectIAK3D.R')) 
source(paste0(iakDir , '/compositeLikIAK3D.R')) 
source(paste0(lmm2Dir , '/fLMM2.R')) 

source(paste0(iakDir , '/optimIt_TEST.R')) 

printnllTime <<- F

crsAusAlbers <- CRS("+proj=aea +lat_1=-18 +lat_2=-36 +lat_0=0 +lon_0=132 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
crsLongLat <- CRS('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')

##########################################################
### set some options for iak...
##########################################################
prodSum <- TRUE

#nud <- NULL
nud <- 1.5

lnTfmdData <- FALSE
rqrBTfmdPreds <- FALSE

useReml <- TRUE
#useReml <- FALSE

#allKnotsd <- c()
allKnotsd <- c(0 , 0.3 , 1 , 6) # update the final value later in script to be the max lower depth.
if(length(allKnotsd) > 0){ incdSpline <- TRUE }else{ incdSpline <- FALSE }

sdfdTypeANDcmeInit <- c(-9 , -1 , -1 , 1)
selectCovModelNow <- FALSE

fitCubistModelNow <- TRUE

###########################################
### load fitting and validation datasets.
###########################################
load(paste0(dataDir , '/cFit.RData'))
load(file = paste0(dataDir , '/dIFit.RData'))
load(file = paste0(dataDir , '/zFit.RData'))
load(file = paste0(dataDir , '/covsFit.RData'))

load(file = paste0(dataDir , '/cVal.RData'))
load(file = paste0(dataDir , '/dIVal.RData'))
load(file = paste0(dataDir , '/zVal.RData'))
load(file = paste0(dataDir , '/covsVal.RData'))

if(fitCubistModelNow){
### fit cubist model...
  cmFit <- cubist(x = covsFit , y = zFit , committees = 1 , cubistControl(rules = nRules))

  print(summary(cmFit))

### convert to des mtx
  tmp <- cubist2X(cubistModel = cmFit, dataFit = covsFit , incdSpline = incdSpline)
  cmFit <- tmp$cubistModel
  XFit <- tmp$X
  matRulesFit <- tmp$matRuleData

  save(cmFit , file = paste0(dataDir , '/cmFit.RData'))

}else{

  load(file = paste0(dataDir , '/cmFit.RData'))

}

if(testCL){
### if fitting different models, save out the compLikMats so that it is the same each time.
#  compLikMats <- setVoronoiBlocksWrap(x = cFit , nPerBlock = 50 , plotVor = T , optnBalance = 2)
### also attach compLikOptn
#  compLikMats$compLikOptn <- 2
#  save(compLikMats , file = paste0(dataDir , '/compLikMats.RData'))

  load(file = paste0(dataDir , '/compLikMats.RData'))

}else{
  compLikMats <- list()
  compLikMats$compLikOptn <- 0
}

lmmFitFile <- paste0(dataDir , '/lmm.fit.selected.RData')

if(fitModelNow){
### refit cubist model as lmm...
  nmplt <- paste0(dataDir , '/plot.selected.pdf')
  
  start_time <- Sys.time()
  
  tmp <- fitIAK3D(x = cFit , dI = dIFit , z = zFit , covs = covsFit , modelX = cmFit , modelx = 'matern' , nud = nud , allKnotsd = allKnotsd , 
                    sdfdType_cd1 = sdfdTypeANDcmeInit[1] , sdfdType_cxd0 = sdfdTypeANDcmeInit[2] , sdfdType_cxd1 = sdfdTypeANDcmeInit[3] , 
                    cmeOpt = sdfdTypeANDcmeInit[4] , prodSum = prodSum , lnTfmdData = lnTfmdData , useReml = useReml , compLikMats = compLikMats ,
                    namePlot = nmplt , rqrBTfmdPreds = rqrBTfmdPreds) 

  end_time <- Sys.time()
  print('Time to fit was:')
  print(end_time - start_time)

  lmm.fit.selected <- tmp$lmmFit

  save(lmm.fit.selected , file = lmmFitFile)

}else{
  load(file = lmmFitFile)
}

#########################################################
### validation bit...
#########################################################
fnamezkVal <- paste0(dataDir , '/zkVal.RData')
fnamevkVal <- paste0(dataDir , '/vkVal.RData')
  
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

  namePlot = paste0(dataDir , '/plotVal.pdf')

  tmp <- plotProfilesIAK3D(namePlot = namePlot , x = cVal , dI = dIVal , z = zVal , 
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

#> sample(83 , 6)
#[1] 45  8 81 14 64 47

  rand6ForPlot <- c(45 , 8 , 81 , 14 , 64 , 47)

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

  plotOnOrigScale <- TRUE

  if(plotOnOrigScale){
    zVal4Plot_PLOT <- exp(zVal4Plot)

    zkVal4Plot_PLOT <- exp(zkVal4Plot)
    pi90LkVal4Plot_PLOT <- exp(zkVal4Plot - 1.64 * sqrt(vkVal4Plot))
    pi90UkVal4Plot_PLOT <- exp(zkVal4Plot + 1.64 * sqrt(vkVal4Plot))

    zkProfPred_PLOT <- exp(zkProfPred)
    pi90LkProfPred_PLOT <- exp(pi90LkProfPred)
    pi90UkProfPred_PLOT <- exp(pi90UkProfPred)
    
    xlab <- expression("EC, dS " ~ m^{-1}) 
    vecTmp <- c(zVal4Plot_PLOT , zkVal4Plot_PLOT , as.numeric(zkProfPred_PLOT))
    xlim <- c(0 , 2)
  }else{
    zVal4Plot_PLOT <- zVal4Plot

    zkVal4Plot_PLOT <- zkVal4Plot
    pi90LkVal4Plot_PLOT <- zkVal4Plot - 1.64 * sqrt(vkVal4Plot)
    pi90UkVal4Plot_PLOT <- zkVal4Plot + 1.64 * sqrt(vkVal4Plot)

    zkProfPred_PLOT <- zkProfPred
    pi90LkProfPred_PLOT <- pi90LkProfPred
    pi90UkProfPred_PLOT <- pi90UkProfPred

    xlab <- expression("log EC, log[dS " ~ m^{-1} ~ "]")
    vecTmp <- c(zVal4Plot_PLOT , zkVal4Plot_PLOT , as.numeric(zkProfPred_PLOT))
    xlim <- c(min(vecTmp) , max(vecTmp))
  }
  
  namePlot = paste0(dataDir , '/plotVal4Plot.pdf')

  tmp <- plotProfilesIAK3D(namePlot = namePlot , x = cVal4Plot , dI = dIVal4Plot , z = zVal4Plot_PLOT , 
                  xPred = cVal4PlotU , dIPred = dIPred , zPred = zkProfPred_PLOT , pi90LPred = pi90LkProfPred_PLOT , pi90UPred = pi90UkProfPred_PLOT , 
                  zhatxv = zkVal4Plot_PLOT , pi90Lxv = pi90LkVal4Plot_PLOT , pi90Uxv = pi90UkVal4Plot_PLOT , 
                  profNames = paste0('Profile ' , rand6ForPlot) , xlim = xlim , xlab = xlab) 
}
### plot pred vs obs for 6 gsm layers...

###########################################################################
### some plots of the fitted covariance model...
###########################################################################
plotVargiogramFit <- FALSE
if(plotVargiogramFit){
  source(paste0(iakDir , '/iaCovMatern.R'))

  dIPlot <- data.frame('dU' = c(0 , 20 , 50 , 80)/100 , 'dL' = c(10 , 30 , 60 , 90)/100)
  hx <- seq(0 , 50 , 10)
  pdf(file = paste0(dataDir , '/varioFitLargeScale.pdf'))
  tmp <- plotCovx(lmm.fit = lmm.fit.selected , hx = hx , dIPlot = dIPlot , addExpmntlV = TRUE , hzntlUnits = 'km')
  dev.off()
  
  hx <- seq(0 , 10 , 2)
  pdf(file = paste0(dataDir , '/varioFitSmallScale.pdf'))
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
  
### plot of variox and variod next to each ofther for paper...  
  pdf(file = paste0(dataDir , '/variosxd.pdf') , width = 12 , height = 6.5)
  par(mfrow = c(1,2))

  dIPlot <- data.frame('dU' = c(0 , 20 , 50 , 80)/100 , 'dL' = c(10 , 30 , 60 , 90)/100)
  hx <- seq(0 , 50 , 10)
  tmp <- plotCovx(lmm.fit = lmm.fit.selected , hx = hx , dIPlot = dIPlot , addExpmntlV = TRUE , hzntlUnits = 'km')

  dTmp <- seq(0 , 2 , 0.1)
  dIPlot <- data.frame('dU' = dTmp[-length(dTmp)] , 'dL' = dTmp[-1])
  qwe <- plotCovd(lmm.fit = lmm.fit.selected , dIPlot = dIPlot , vrtclUnits = 'm')
  dev.off()

  ### plot of the variances...
  pdf(file = paste0(dataDir , '/varComps.pdf'))
  dPlot <- seq(0 , 2 , 0.01)
  plotVarComps(lmm.fit = lmm.fit.selected , dPlot = dPlot , xlim = c(0 , 0.85))
  dev.off()
   

}else{}


########################################################################
### set up covaraites for mapping, based on the GSM mapping grid, taken from ...
########################################################################
if(mapNow){

#  dIMap <- data.frame('dU' = c(0 , 0.05 , 0.15 , 0.3 , 0.6 , 1.0) , 'dL' = c(0.05 , 0.15 , 0.3 , 0.6 , 1.0 , 2.0))
  dIMap <- data.frame('dU' = c(0 , 0.15 , 0.6) , 'dL' = c(0.05 , 0.3 , 1.0))

### get which rows are needed in this batch...
  rowsToDo <- firstRow:lastRow

  for(irow in rowsToDo){

### define cMap and extract covsMap for this row...

################################################
### calculate the predictions for the mapping depth...
################################################
    lmm.map.tmp <- predictIAK3D(xMap = cMap , covsMap = covsMap , dIMap = dIMap , lmmFit = lmm.fit.selected , rqrBTfmdPreds = rqrBTfmdPreds , constrainX4Pred = constrainX4Pred)

    save(lmm.map.tmp , file = paste0(dataDir , '/map.row' , irow , '.RData')))

  }
  
  print(paste0('Mapped for row ' , irow))
  
} ### done

stop(paste0('Mapped for batch ' , iBatch))
