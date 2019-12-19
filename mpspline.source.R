####################################################################################
### the following code is from the mpspline function source code in the GSIF package
### taken out some of the checks and uses a data.frame input to fit in with IAK code. 
### profiles with overlap are found and overlapping depth intervals removed. 
####################################################################################
mpspline.source <- function (dfHzns, ...) 
{
### dfHzns has ID, Eastings, Northings, dU, dL, z...
### final column is target variable z
    .local <- function (obj, var.name, lam = 0.1, d = t(c(0, 
        5, 15, 30, 60, 100, 200)), vlow = 0, vhigh = 1000, show.progress = TRUE) 
    {
        dfHzns <- dfHzns[(!(dfHzns$dU < 0) & !(dfHzns$dL < 0)),,drop=FALSE]
        objd <- dfHzns
        if(ncol(dfHzns) > 6){ stop('dfHzns should have 6 cols, ID, Eastings, Northings, dU, dL, and target variable.') }else{}
        var.name <- names(dfHzns)[ncol(dfHzns)]
        depthcols <- c('dU' , 'dL')
        idcol <- 'ID'
        dfHznsTmp <- rmOverlapAll(dfHzns)
        iTmp <- which(!duplicated(rbind(dfHznsTmp , dfHzns))[-seq(nrow(dfHznsTmp))])
        dfHznsRmvd <- dfHzns[iTmp,c(idcol , depthcols , var.name)]
        dfHzns <- dfHznsTmp[,c(idcol , depthcols , var.name)]
        
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
            
            # ### a bit to predict the hzns that were removed because of overlap
            # ### should give some idea of measurement error.
            #             if(is.element(subs$ID[1] , dfHznsRmvd$ID)){
            #                 subsPred <- dfHznsRmvd[which(dfHznsRmvd$ID == subs$ID[1]),]
            #                 m_fyfit <- 
            #             }else{}
            
            
            
        }
        
        ### a bit to predict the hzns that were removed because of overlap
        ### should give some idea of measurement error.
        if(nrow(dfHznsRmvd) > 0){
            dfHznsRmvd$predicted <- NA
            for(i in 1:nrow(dfHznsRmvd)){
                iThis <- which(IDU == dfHznsRmvd$ID[i])
                if(round(dfHznsRmvd$dL[i]) > ncol(m_fyfit)){
                    dfHznsRmvd$predicted[i] <- NA
                }else{
                    dfHznsRmvd$predicted[i] <- mean(m_fyfit[iThis,seq(round(dfHznsRmvd$dU[i]) , round(dfHznsRmvd$dL[i]))] , na.rm = TRUE)
                }
            }
        }else{}

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
            var.std = yave, var.1cm = t(m_fyfit) , dfHznsRmvd = dfHznsRmvd)

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
##################################################
harmonizeMPS <- function(xData , dIData , zData , dIStd , overlapsByRegression = TRUE){
    profIDData <- makeProfID(xData)
    nProfs <- max(profIDData)
    profEAS <- data.frame('ID' = profIDData , 'Eastings' = xData[,1]  , 'Northings' = xData[,2] , 
                          'dU' = 100 * dIData[,1] , 'dL' = 100 * dIData[,2] , 'z' = zData) 

    if(overlapsByRegression){
      IDMT1 <- unique(profIDData[which(duplicated(profIDData))])
      ID1 <- setdiff(seq(nProfs) , IDMT1)
      IDOVLP <- unique(whichProfilesOverlap(profEAS)$ID) 
      
      iEAS <- which(!is.element(profIDData , c(ID1 , IDOVLP)))        
    }else{
      iEAS <- seq(length(zData))
    }

    maxd <- max(max(dIStd[,2]) , max(dIData[,2]))
    
    z.sEAS <- mpspline.source(dfHzns = profEAS[iEAS,,drop=FALSE] , lam = 0 , vlow = -1000 , vhigh = 1000 , d = t(c(0,100*maxd)))

    xDataH <- xData[!duplicated(profIDData),,drop=FALSE]
    var.1cm <- matrix(NA , 100*maxd , nProfs)
    for(i in 1:length(z.sEAS$idcol)){
        var.1cm[,z.sEAS$idcol[i]] <- z.sEAS$var.1cm[,i]
    }
    remove(z.sEAS)

    hrmnzdData <- matrix(NA , nProfs , nrow(dIStd))
    for (i in 1:nrow(dIStd)){
        hrmnzdData[,i] <- colSums(var.1cm[seq(round(dIStd[i,1]*100) + 1 , round(dIStd[i,2]*100)),,drop=FALSE]) / ((dIStd[i,2]-dIStd[i,1])*100)
    }
    hrmnzdDataEAS <- !is.na(hrmnzdData) # TRUE where EAS fitted, FALSE will be a regression or missing

    if(overlapsByRegression){
      for(i in 1:nrow(hrmnzdData)){
        iThis <- which(profIDData == i)
        for(j in 1:ncol(hrmnzdData)){
            if(is.na(hrmnzdData[i,j])){
### what is missing?
                dITarget <- dIStd[j,,drop=FALSE]
### what data do we have in this profile? 
                dIDataThis <- dIData[iThis,,drop=FALSE]
### use fitted splines where the target and data depths are covered to fit regression...
                XDataThis <- matrix(1 , nProfs , length(iThis) + 1)
                XPredThis <- matrix(c(1 , zData[iThis]) , 1 , length(iThis) + 1)
                for(k in 1:length(iThis)){
                    XDataThis[,k+1] <- colSums(var.1cm[seq(round(dIDataThis[k,1]*100) + 1 , round(dIDataThis[k,2]*100)),,drop=FALSE]) / ((dIDataThis[k,2]-dIDataThis[k,1])*100)
                }
                zDataThis <- colSums(var.1cm[seq(round(dITarget[1,1]*100) + 1 , round(dITarget[1,2]*100)),,drop=FALSE]) / ((dITarget[1,2]-dITarget[1,1])*100)

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

    return(list('xDataH' = xDataH , 'hrmnzdData' = hrmnzdData , 'hrmnzdDataEAS' = hrmnzdDataEAS))    
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

whichProfilesOverlap <- function(dfAll){
  lProfs <- list()
  cAll <- dfAll[,c('Eastings' , 'Northings')]
  cAll <- cAll[which(!duplicated(cAll)),,drop=FALSE]
  idsProfs <- NA * numeric(nrow(cAll))
  for(i in 1:nrow(cAll)){
    iThis <- which(dfAll$Eastings == cAll[i,1] & dfAll$Northings == cAll[i,2])
    lProfs[[i]] <- dfAll[iThis,,drop=FALSE]
    idsProfs[i] <- dfAll$ID[iThis[1]]
  }
  
  iOVLP <- which(unlist(lapply(lProfs , isOverlap)))
  if(length(iOVLP) > 0){
    idsOVLP <- idsProfs[iOVLP]
    dfOVLP <- dfAll[which(is.element(dfAll$ID , idsOVLP)),,drop=FALSE]
  }else{
    dfOVLP <- dfAll[c(),,drop=FALSE]
  }
  
  return(dfOVLP)
}

rmOverlapAll <- function(dfIn){
  profList <- list()
  cAll <- dfIn[,c('Eastings' , 'Northings')]
  cU <- cAll[!duplicated(cAll),,drop=FALSE]
  for(i in 1:nrow(cU)){
    profList[[i]] <- dfIn[which(dfIn$Eastings == cU[i,1] & dfIn$Northings == cU[i,2]),,drop=FALSE]
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

