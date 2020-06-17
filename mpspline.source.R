####################################################################################
### the following code is from the mpspline function source code in the GSIF package
### taken out some of the checks and uses a data.frame input to fit in with IAK code. 
### profiles with overlap are found and overlapping depth intervals removed. 
####################################################################################
mpspline.source <- function (dfHzns, ...) 
{
### dfHzns has ID, dU, dL, z...
### final column is target variable z
    .local <- function (obj, var.name, lam = 0.1, d = t(c(0, 
        5, 15, 30, 60, 100, 200)), vlow = 0, vhigh = 1000, show.progress = TRUE) 
    {
        dfHzns <- dfHzns[(!(dfHzns$dU < 0) & !(dfHzns$dL < 0)),,drop=FALSE]
        objd <- dfHzns
        if(ncol(dfHzns) > 4){ stop('dfHzns should have 4 cols, ID, dU, dL, and target variable called z.') }else{}
        var.name = 'z'
        if(!is.element('z' , names(dfHzns))){ stop('Error - for mpspline.source, dfHzns must have a column called z.') }else{}
        depthcols <- c('dU' , 'dL')
        if((!is.element('dU' , names(dfHzns))) | (!is.element('dL' , names(dfHzns)))){ stop('Error - for mpspline.source, dfHzns must have columns called dU and dL.') }else{}
        idcol <- 'ID'
        if(!is.element('ID' , names(dfHzns))){ stop('Error - for mpspline.source, dfHzns must have a column called ID.') }else{}

        # dfHznsTmp <- rmOverlapAll(dfHzns)
        # iTmp <- which(!duplicated(rbind(dfHznsTmp , dfHzns))[-seq(nrow(dfHznsTmp))])
        # dfHznsRmvd <- dfHzns[iTmp,c(idcol , depthcols , var.name)]
        # dfHzns <- dfHznsTmp[,c(idcol , depthcols , var.name)]

### make sure order is ID, dU, dL, and target variable called z
        dfHzns <- dfHzns[,c(idcol , depthcols , var.name)]
        
        IDU <- unique(dfHzns$ID)
        sel <- logical(length(IDU))
        mHzns <- NA * numeric(length(IDU))
        for(i in 1:length(IDU)){ 
          iThis <- which(dfHzns$ID == IDU[i])
          if(length(iThis) > 1){ sel[i] <- TRUE }else{ sel[i] <- FALSE }
          mHzns[i] <- length(iThis)
        }
        ndata <- length(IDU)
        np <- max(mHzns)

        mxd <- max(d)
        m_fyfit <- matrix(NA, ncol = length(c(1:mxd)), nrow = ndata)
        yave <- matrix(NA, ncol = length(d), nrow = ndata)
        sse <- matrix(NA, ncol = length(lam), nrow = 1)
        sset <- matrix(NA, ncol = length(lam), nrow = ndata)
        nl <- length(lam)
        s <- 0.05 * sd(dfHzns[[var.name]], na.rm = TRUE)
        s2 <- s * s
        dave <- matrix(NA, ncol = np, nrow = ndata)
        if (np < 2 | is.na(np)) {
            print("Submitted soil profiles all have 1 horizon")
        }

################################################
        for (st in as.vector(which(sel))) {
            iThis <- which(dfHzns$ID == IDU[st])
            subs <- dfHzns[iThis,,drop=FALSE]
          
            ir <- c(1:length(subs[, 1]))
            ir <- as.matrix(t(ir))
            u <- subs[ir, 2]
            u <- as.matrix(t(u))
            v <- subs[ir, 3]
            v <- as.matrix(t(v))
            y <- subs[ir, 4]
            y <- as.matrix(t(y))
            n <- length(y)
            if (n == 1) {
                message(paste("Spline not fitted to profile:", 
                              objd_m[st, 1], sep = " "))
                xfit <- as.matrix(t(c(1:mxd)))
                nj <- max(v)
                if (nj > mxd) {
                    nj <- mxd
                }
                yfit <- xfit
                yfit[, 1:nj] <- y
                if (nj < mxd) {
                    yfit[, (nj + 1):mxd] = NA
                }
                m_fyfit[st, ] <- yfit
                nd <- length(d) - 1
                dl <- d + 1
                for (cj in 1:nd) {
                    xd1 <- dl[cj]
                    xd2 <- dl[cj + 1] - 1
                    if (nj >= xd1 & nj <= xd2) {
                        xd2 <- nj - 1
                        yave[st, cj] <- mean(yfit[, xd1:xd2])
                    }
                    else {
                        yave[st, cj] <- mean(yfit[, xd1:xd2])
                    }
                    yave[st, cj + 1] <- max(v)
                }
            }
            else {
                np1 <- n + 1
                nm1 <- n - 1
                delta <- v - u
                del <- c(u[2:n], u[n]) - v
                r <- matrix(0, ncol = nm1, nrow = nm1)
                for (dig in 1:nm1) {
                    r[dig, dig] <- 1
                }
                for (udig in 1:nm1 - 1) {
                    r[udig, udig + 1] <- 1
                }
                d2 <- matrix(0, ncol = nm1, nrow = nm1)
                diag(d2) <- delta[2:n]
                r <- d2 %*% r
                r <- r + t(r)
                d1 <- matrix(0, ncol = nm1, nrow = nm1)
                diag(d1) <- delta[1:nm1]
                d3 <- matrix(0, ncol = nm1, nrow = nm1)
                diag(d3) <- del[1:nm1]
                r <- r + 2 * d1 + 6 * d3
                q <- matrix(0, ncol = n, nrow = n)
                for (dig in 1:n) {
                    q[dig, dig] <- -1
                }
                for (udig in 1:n - 1) {
                    q[udig, udig + 1] <- 1
                }
                q <- q[1:nm1, 1:n]
                dim.mat <- matrix(q[], ncol = length(1:n), nrow = length(1:nm1))
                rinv <- try(solve(r), TRUE)
                if (is.matrix(rinv)) {
                    ind <- diag(n)
                    pr.mat <- matrix(0, ncol = length(1:nm1), nrow = length(1:n))
                    pr.mat[] <- 6 * n * lam
                    fdub <- pr.mat * t(dim.mat) %*% rinv
                    z <- fdub %*% dim.mat + ind
                    sbar <- solve(z, t(y))
                    b <- 6 * rinv %*% dim.mat %*% sbar
                    b0 <- rbind(0, b)
                    b1 <- rbind(b, 0)
                    gamma <- (b1 - b0)/t(2 * delta)
                    alfa <- sbar - b0 * t(delta)/2 - gamma * t(delta)^2/3
                    xfit <- as.matrix(t(c(1:mxd)))
                    nj <- max(v)
                    if (nj > mxd) {
                        nj <- mxd
                    }
                    yfit <- xfit
                    for (k in 1:nj) {
                        xd <- xfit[k]
                        if (xd < u[1]) {
                            p <- alfa[1]
                        }
                        else {
                            for (its in 1:n) {
                                if (its < n) {
                                    tf2 = as.numeric(xd > v[its] & xd < 
                                                         u[its + 1])
                                }
                                else {
                                    tf2 <- 0
                                }
                                if (xd >= u[its] & xd <= v[its]) {
                                    p = alfa[its] + b0[its] * (xd - u[its]) + 
                                        gamma[its] * (xd - u[its])^2
                                }
                                else if (tf2) {
                                    phi = alfa[its + 1] - b1[its] * (u[its + 
                                                                           1] - v[its])
                                    p = phi + b1[its] * (xd - v[its])
                                }
                            }
                        }
                        yfit[k] = p
                    }
                    if (nj < mxd) {
                        yfit[, (nj + 1):mxd] = NA
                    }
                    yfit[which(yfit > vhigh)] <- vhigh
                    yfit[which(yfit < vlow)] <- vlow
                    m_fyfit[st, ] <- yfit
                    nd <- length(d) - 1
                    dl <- d + 1
                    for (cj in 1:nd) {
                        xd1 <- dl[cj]
                        xd2 <- dl[cj + 1] - 1
                        if (nj >= xd1 & nj <= xd2) {
                            xd2 <- nj - 1
                            yave[st, cj] <- mean(yfit[, xd1:xd2])
                        }
                        else {
                            yave[st, cj] <- mean(yfit[, xd1:xd2])
                        }
                        yave[st, cj + 1] <- max(v)
                    }
                    dave[st, 1:n] <- sbar
                    ssq <- sum((t(y) - sbar)^2)
                    g <- solve(z)
                    ei <- eigen(g)
                    ei <- ei$values
                    df <- n - sum(ei)
                    sig2w <- ssq/df
                    dfc <- n - 2 * sum(ei) + sum(ei^2)
                    sig2c <- ssq/dfc
                    tmse <- ssq/n - 2 * s2 * df/n + s2
                    sset[st] <- tmse
                }
            }
         }
        
        # ### a bit to predict the hzns that were removed because of overlap
        # ### should give some idea of measurement error.
        # if(nrow(dfHznsRmvd) > 0){
        #     dfHznsRmvd$predicted <- NA
        #     for(i in 1:nrow(dfHznsRmvd)){
        #         iThis <- which(IDU == dfHznsRmvd$ID[i])
        #         if(round(dfHznsRmvd$dL[i]) > ncol(m_fyfit)){
        #             dfHznsRmvd$predicted[i] <- NA
        #         }else{
        #             dfHznsRmvd$predicted[i] <- mean(m_fyfit[iThis,seq(round(dfHznsRmvd$dU[i]) , round(dfHznsRmvd$dL[i]))] , na.rm = TRUE)
        #         }
        #     }
        # }else{}

        yave <- as.data.frame(yave)
        jmat <- matrix(NA, ncol = 1, nrow = length(d))
        for (i in 1:length(d) - 1) {
            a1 <- paste(d[i], d[i + 1], sep = "-")
            a1 <- paste(a1, "cm", sep = " ")
            jmat[i] <- a1
        }
        jmat[length(d)] <- "soil depth"
        for (jj in 1:length(jmat)) {
            names(yave)[jj] <- jmat[jj]
        }
        retval <- list(idcol = IDU, var.fitted = dave, 
            var.std = yave, var.1cm = t(m_fyfit)) #  , dfHznsRmvd = dfHznsRmvd

        return(retval)
    }
    .local(obj, ...)
}

##################################################
### a wrapper that only fits eas if >1 depth and no overlaps
### else uses data from fitted splines to make regression for target depths
### with predictor variables = known depths in the target profile
### in paper, horizons were removed from profiles with overlap so that ea spline would work
### to redo this, put overlapsByRegression = FALSE; hzns get removed in mpspline fn.
### NOTE - dIData and dIStd are in m as inputs here,
###        but mpspline function uses cm.
###        Also, default for this fn is vlow = -1000, while default for mpspline.source fn is vlow = 0
##################################################
harmonizeMPS <- function(profIDData , dIData , zData , dIStd , vlow = -1000 , vhigh = 1000 , singlesByRegression = TRUE){

    profIDData <- as.character(profIDData)

### only include where we have some data...
    iok <- which(!is.na(zData))
    if(length(iok) > 0){
      profIDData <- profIDData[iok]
      dIData <- dIData[iok,,drop=FALSE]
      zData <- zData[iok]
    }else{
      stop('No data given to harmonizeMPS!')
    }
  
    profEAS <- data.frame('ID' = profIDData , 'dU' = 100 * dIData[,1] , 'dL' = 100 * dIData[,2] , 'z' = zData) 

    profEAS <- averageOverlapAll(profEAS , idVar = "ID")
    
    
    if(singlesByRegression){
      IDMT1 <- unique(profEAS$ID[which(duplicated(profEAS$ID))])
      iEAS <- which(is.element(profEAS$ID , IDMT1))
    }else{
      iEAS <- seq(nrow(profEAS))
    }

    maxdcm <- round(max(max(round(100*dIStd[,2])) , max(profEAS$dL)))
    
    z.sEAS <- mpspline.source(dfHzns = profEAS[iEAS,,drop=FALSE] , lam = 0 , vlow = vlow , vhigh = vhigh , d = t(c(0,maxdcm)))

    profIDDataH <- unique(profEAS$ID)
    nProfs <- length(profIDDataH)

    ### re order var.1cm according to profIDDataH...    
    var.1cm <- matrix(NA , maxdcm , nProfs)
    for(i in 1:length(z.sEAS$idcol)){
        var.1cm[,which(profIDDataH == z.sEAS$idcol[i])] <- z.sEAS$var.1cm[,i]
    }
    remove(z.sEAS)
    
### all profiles that had spline fitted used more than 1 depth interval; 
### for harmonizing, we will extend the prediction at deepest dL to cover the target depths
### in future, could be worth looking at doind these by regression, using the splined average 
### over the covered depths as single predictor.
    
    hrmnzdData <- matrix(NA , nProfs , nrow(dIStd))
    for (i in 1:nrow(dIStd)){
      ncmThis <- round((dIStd[i,2]-dIStd[i,1])*100)
      zTmp <- var.1cm[seq(round(dIStd[i,1]*100) + 1 , round(dIStd[i,2]*100)),,drop=FALSE]
### fill zTmp down if at least some values in a col of zTmp      
      getLastElement <- function(xVec){ if(all(is.na(xVec))){ return(NA) }else{ xVec <- xVec[!is.na(xVec)] ; return(xVec[length(xVec)]) } }
      zTmpFillVals = apply(zTmp , 2 , getLastElement)

      fillNAs <- function(xVec , fillVal){ iTmp <- which(is.na(xVec)) ; if(length(iTmp) == 0){ return(xVec) }else{ xVec[iTmp] <- fillVal ; return(xVec) }  }
      zTmp <- mapply(fillNAs , xVec = lapply(seq_len(ncol(zTmp)), function(i) zTmp[,i]) , fillVal = zTmpFillVals , SIMPLIFY = FALSE)
      zTmp <- matrix(unlist(zTmp) , ncmThis , nProfs)
      
      hrmnzdData[,i] <- colSums(zTmp) / ncmThis
    }
    hrmnzdDataEAS <- !is.na(hrmnzdData) # TRUE where EAS fitted, FALSE will be a regression or missing

    if(singlesByRegression){
      for(i in 1:nrow(hrmnzdData)){
        iThis <- which(profEAS$ID == profIDDataH[i])
        for(j in 1:ncol(hrmnzdData)){
            if(is.na(hrmnzdData[i,j])){
### what is missing?
                dITarget <- dIStd[j,,drop=FALSE] * 100 # this is now in cm
### what data do we have in this profile? 
                dIDataThis <- as.matrix(profEAS[iThis,c('dU','dL'),drop=FALSE]) # this is already in cm
### use fitted splines where the target and data depths are fully covered to fit regression...
                XDataThis <- matrix(1 , nProfs , length(iThis) + 1)
                XPredThis <- matrix(c(1 , profEAS$z[iThis]) , 1 , length(iThis) + 1)
                for(k in 1:length(iThis)){
                  XDataThis[,k+1] <- colSums(var.1cm[seq(round(dIDataThis[k,1]) + 1 , round(dIDataThis[k,2])),,drop=FALSE]) / (dIDataThis[k,2]-dIDataThis[k,1])
                }
                zDataThis <- colSums(var.1cm[seq(round(dITarget[1,1]) + 1 , round(dITarget[1,2])),,drop=FALSE]) / (dITarget[1,2]-dITarget[1,1])

### look for best regression here (fewer predictors could be more available data)
                nVec <- vPredVec <- zPredVec <- NA * numeric(length(iThis))
                for(k in 1:length(iThis)){
                    iOK <- which((!is.na(rowSums(XDataThis[,1:(k+1)]))) & !is.na(zDataThis))
                    if(length(iOK) >= (k+1)){
                      tmp <- nllLm(zData = zDataThis[iOK] , XData = XDataThis[iOK,1:(k+1),drop=FALSE] , REML = TRUE)
                      nll <- tmp$nll 
                      betahat <- tmp$betahat 
                      vbetahat <- tmp$vbetahat
                      sigma2hat <- tmp$sigma2hat
                      
                      nVec[k] <- length(iOK)
                      zPredVec[k] <- XPredThis[,1:(k+1),drop=FALSE] %*% betahat
                      vPredVec[k] <- sigma2hat + XPredThis[,1:(k+1),drop=FALSE] %*% vbetahat %*% t(XPredThis[,1:(k+1),drop=FALSE])
                    }else{}
                }
                if(all(is.na(zPredVec))){
                  hrmnzdData[i,j] <- NA
                }else{
                  iBest <- which.min(vPredVec)
                  hrmnzdData[i,j] <- zPredVec[iBest]
                }
                
            }else{}
        }
      }
    }else{}
    
    hrmnzdData[hrmnzdData < vlow] <- vlow
    hrmnzdData[hrmnzdData > vhigh] <- vhigh
    
    return(list('profIDDataH' = profIDDataH , 'hrmnzdData' = hrmnzdData , 'hrmnzdDataEAS' = hrmnzdDataEAS , 'var.1cm' = var.1cm))    
}

isOverlap <- function(profIn){
  dU <- profIn$dU
  dL <- profIn$dL
  
  io <- order(dU)
  dU <- dU[io]
  dL <- dL[io]
  
  if(min(dL - dU) <= 0){
    ### a profile with a point observation - define as overlap for now.
    overlap <- TRUE
  }else{
    if(length(dU) > 1){
      dU <- dU[-1]
      dL <- dL[-length(dL)]
      
      if(max(dL - dU) > 0){
        overlap <- TRUE
      }else{
        overlap <- FALSE
      }
    }else{
      overlap <- FALSE
    }
  }
  return(overlap)
}

whichProfilesOverlap <- function(dfAll , idVar = 'ID'){
  dfAll[[idVar]] <- as.character(dfAll[[idVar]])
  idsProfs <- unique(dfAll[[idVar]])
  lProfs <- list()
  for(i in 1:length(idsProfs)){
    iThis <- which(dfAll[[idVar]] == idsProfs[i])
    lProfs[[i]] <- dfAll[iThis,,drop=FALSE]
  }
  ### something funny going on with large datasets here - seems to complete loop but hang. not sure why. works ok when entire fn is run.

  iOVLP <- which(unlist(lapply(lProfs , isOverlap)))
  if(length(iOVLP) > 0){
    idsOVLP <- idsProfs[iOVLP]
    dfOVLP <- dfAll[which(is.element(as.character(dfAll[[idVar]]) , idsOVLP)),,drop=FALSE]
  }else{
    dfOVLP <- dfAll[c(),,drop=FALSE]
  }
  
  return(dfOVLP)
}

rmOverlapAll <- function(dfIn , idVar = 'ID'){
  dfIn[[idVar]] <- as.character(dfIn[[idVar]])

  idsProfs <- unique(dfIn[[idVar]])
  profList <- list()
  for(i in 1:length(idsProfs)){
    profList[[i]] <- dfIn[which(dfIn[[idVar]] == idsProfs[i]),,drop=FALSE]
  }
  
  profListOut <- lapply(profList , rmOverlap)
  profListOut <- do.call(rbind.data.frame, profListOut)
  
  return(profListOut)
}

rmOverlap <- function(profIn){
  profIn <- profIn[order(profIn$dU),,drop=FALSE]
  iOK <- which(((profIn$dL - profIn$dU) > 0) & (profIn$dU >= 0) & (profIn$dL > 0))
  profIn <- profIn[iOK,,drop=FALSE]
  
  if(nrow(profIn) > 1){
    overlap <- TRUE
    while(overlap){
      dUTmp <- profIn$dU[-1]
      dLTmp <- profIn$dL[-length(profIn$dL)]
      
      iBad <- which((dLTmp - dUTmp) > 0)
      if(length(iBad) == 0){
        overlap <- FALSE
      }else{
        iBad <- iBad[1]
        ### if one of the offending pair has dU = 0 and the other doesn't, keep the surface one.
        ### else remove the one of the offending pair with the larger depth interval...    
        if(profIn$dU[iBad] == 0 & profIn$dU[iBad+1] > 0){
          iRm <- iBad + 1
        }else if(profIn$dU[iBad] > 0 & profIn$dU[iBad+1] == 0){
          iRm <- iBad
        }else if((profIn$dL[iBad] - profIn$dU[iBad]) > (profIn$dL[iBad+1] - profIn$dU[iBad+1])){
          iRm <- iBad
        }else{
          iRm <- iBad + 1
        }
        profIn <- profIn[-iRm,,drop=FALSE]
        if(nrow(profIn) > 1){
          overlap <- TRUE
        }else{
          overlap <- FALSE
        }
      }
    }  
  }else{}
  
  return(profIn)
}

### averageOverlap for variable z in the input data.frame
### other info (apart from depth) will be copied down the data.frame.
### all depths will be rounded to 1 cm.
averageOverlap <- function(profIn){
  
  profIn$dU <- round(profIn$dU , digits = 2)
  profIn$dL <- round(profIn$dL , digits = 2)
  
  profIn <- profIn[order(profIn$dU),,drop=FALSE]
  iOK <- which(((profIn$dL - profIn$dU) > 0) & (profIn$dU >= 0) & (profIn$dL > 0))
  profIn <- profIn[iOK,,drop=FALSE]

  if (nrow(profIn)> 1){
    m1cm <- matrix(NA , round(max(profIn$dL) * 100) , nrow(profIn))
    for(i in 1:nrow(profIn)){
      m1cm[(round(100*profIn$dU[i])+1):round(100*profIn$dL[i]),i] <- profIn$z[i]
    }
    dUNew <- unique(c(profIn$dU, profIn$dL))
    dUNew <- dUNew[order(dUNew)]
    dLNew <- dUNew[-1]
    dUNew <- dUNew[-length(dUNew)]

    if(is.element('Observation_ID_LIST' , names(profIn))){
      oidThis <- unique(unlist(strsplit(as.character(profIn$Observation_ID_LIST) , '_AND_')))    
      profIn$Observation_ID_LIST <- paste(oidThis , collapse = '_AND_')
    }else{}
    
### samples combined - just join all profile sample ids here, as difficult to properly keep track anymore
    if(is.element('SampleID_LIST' , names(profIn))){
      sidThis <- unique(unlist(strsplit(as.character(profIn$SampleID_LIST) , '_AND_')))    
      profIn$SampleID_LIST <- paste(sidThis , collapse = '_AND_')
    }else{}
    
    profOut <- profIn[rep(1 , length(dUNew)),,drop=FALSE]
    profOut$dU <- dUNew
    profOut$dL <- dLNew
    rm(dUNew , dLNew)
 
   if(is.element('dIMidPts' , names(profOut))){
     profOut$dIMidPts <- 0.5 * (profOut$dL + profOut$dU)
    }else{}
    
    profOut$z <- NA

    for(i in 1:nrow(profOut)){
      mThis <- m1cm[(round(100*profOut$dU[i])+1):round(100*profOut$dL[i]),,drop=FALSE]
      if(!all(is.na(mThis))){
        profOut$z[i] <- mean(mThis , na.rm = T)
      }else{}
    }
    profOut <- profOut[which(!is.na(profOut$z)),,drop=FALSE]
  }else{
    profOut <- profIn
  }

  return(profOut)  
}

averageOverlapAll <- function(dfIn , idVar = 'ID'){
### split into profiles that overlap and ones that don't...
  dfIn[[idVar]] <- as.character(dfIn[[idVar]])

  dfOVLP <- whichProfilesOverlap(dfIn , idVar = idVar)
  if(nrow(dfOVLP) == 0){
    ### no overlap, so return inputted df...
    dfOut <- dfIn
  }else{
    ### start output with all profiles that didn't overlap...  
    if(nrow(dfOVLP) == nrow(dfIn)){
### every profile had overlap, so start with empty df...
      dfOut <- dfIn[c(),,drop=FALSE] 
    }else{
      idTmpOVLP <- unique(dfOVLP[[idVar]])
      iNoOvlp <- which(!is.element(dfIn[[idVar]] , idTmpOVLP))
      dfOut <- dfIn[iNoOvlp,,drop=FALSE]
    }
    
    idsU <- unique(dfOVLP[[idVar]])
    profList <- list()
    for(i in 1:length(idsU)){
      profList[[i]] <- dfOVLP[which(dfOVLP[[idVar]] == idsU[i]),,drop=FALSE]
    }

    for(i in 1:length(profList)){
      profListOut <- averageOverlap(profList[[i]])
    }

    profListOut <- lapply(profList , averageOverlap)
    profListOut <- do.call(rbind.data.frame, profListOut)
    dfOut <- rbind(dfOut , profListOut)
  }  

  return(dfOut)
}

