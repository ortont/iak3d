cubist2X <- function(cubistModel , dataFit , zFit = NULL , profIDFit = NULL , allKnotsd = c() , removeColinCols = TRUE , refineCubistModel = FALSE){

### removeColinCols and refineCubistModel are used in the cubist2XSetup function. 
  
    if(length(allKnotsd) == 0){
      incdSpline <- FALSE
    }else if(length(allKnotsd) >= 2){
      incdSpline <- TRUE 
    }else{
      stop('Error - if entering knots, then enter at least 2 knots (2 boundary knots plus any internal knots)')
    }

  if(is.null(cubistModel$listRules) || (is.null(cubistModel$allKnotsd) & (length(allKnotsd) > 0))){    

      tmp <- cubist2XSetup(cubistModel = cubistModel , dataFit = dataFit , zFit = zFit , profIDFit = profIDFit , allKnotsd = allKnotsd , removeColinCols = removeColinCols , refineCubistModel = refineCubistModel)
      X <- tmp$X 
      matRuleData <- tmp$matRuleData
      cubistModel <- tmp$cubistModel
            
    }else{
### check order of columns...if different order, then resort...
      if((length(names(dataFit)) != length(cubistModel$namesDataFit)) || any(names(dataFit) != cubistModel$namesDataFit)){
        iOrdered <- NA * integer(length(cubistModel$namesDataFit))
        for (i in 1:length(cubistModel$namesDataFit)){
          iOrdered[i] <- which(names(dataFit) == cubistModel$namesDataFit[i])
        }
        dataFit <- dataFit[,iOrdered,drop=FALSE]
      }else{}
      
      tmp <- cubist2XGivenSetup(cubistModel , dataFit)
      X <- tmp$X 
      matRuleData <- tmp$matRuleData
    }

    iNA <- which(rowSums(is.na(X) | is.nan(X)) > 0)
    if(length(iNA) > 0){
      X[iNA,] <- NA
    }else{}

    return(list('X' = X , 'matRuleData' = matRuleData , 'cubistModel' = cubistModel))

}

cubist2XGivenSetup <- function(cubistModel , dataFit){ 

    nRules <- nrow(cubistModel$coefficients)
    n <- nrow(dataFit)

#################################################
### get the total number of linear parameters...
#################################################
    p <- length(cubistModel$names4XCubist)
    X <- matrix(0 , n , p)
    colnames(X) <- cubistModel$names4XCubist
    matRuleData <- matrix(0 , n , nRules)

    for (iR in 1:nRules){
      if(nRules > 1){
        
        nVThis <- nrow(cubistModel$listRules[[iR]])
        
### get all the data that fall into rule iR...        
        for (jS in 1:nVThis){

            dirThis <- cubistModel$listRules[[iR]]$dir[jS]
            valThis <- cubistModel$listRules[[iR]]$valUpdated[jS]
            ivariableThis <- cubistModel$listRules[[iR]]$ivariable[jS]

            if (dirThis == '<='){
                iThisIneq <- which(dataFit[,ivariableThis] <= valThis)
            }else if (dirThis == '>'){
                iThisIneq <- which(dataFit[,ivariableThis] > valThis)
            }else{
                vecCatsThis <- cubistModel$listSplitCats[[iR]][[jS]]
                iThisIneq <- which(is.element(dataFit[,ivariableThis] , vecCatsThis))
            }
            
            if (jS == 1){
                iThis <- iThisIneq
            }else{
                iThis <- intersect(iThis , iThisIneq)
            }
        }
                    
        matRuleData[iThis , iR] <- 1  
      
      }else{
### only one rule, so no condits in it!
          iThis <- seq(n)
          matRuleData[iThis,1] <- 1
      }
             
      if (iR == 1){
          jThis <- 1:length(which(!is.na(cubistModel$dfCoeffs[1:iR,])))
      }else{
          jThis <- (length(which(!is.na(cubistModel$dfCoeffs[1:iR-1,]))) + 1):length(which(!is.na(cubistModel$dfCoeffs[1:iR,])))
      }    

      ivThis <- which(!is.na(cubistModel$dfCoeffs[iR,]))
      namesThis <- names(cubistModel$dfCoeffs)[ivThis]
      pThis <- cubistModel$pVec[iR]

### put the values into X...
      X[iThis,jThis[1]] <- 1
      if(length(jThis) > 1){
        if(length(jThis) == (length(cubistModel$listCoeffs_iv[[iR]]) + 1)){
### first of these columns is const, not in listCoeffs_iv         
          X[iThis,jThis[-1]] <- as.matrix(dataFit[iThis,cubistModel$listCoeffs_iv[[iR]]])
        }else if(length(jThis) == length(cubistModel$listCoeffs_iv[[iR]])){
### this rule does not have a const, must have been removed because of colinearity between rules...
          X[iThis,jThis] <- as.matrix(dataFit[iThis,cubistModel$listCoeffs_iv[[iR]]])
        }else{
### something gone wrong.
          print(iR)
          
          print(jThis)
          print(length(jThis))
          
          print(cubistModel$listCoeffs_iv[[iR]])
          print(length(cubistModel$listCoeffs_iv[[iR]]))
          
          stop('Error - something has gone wrong with listCoeffs_iv!')
        }
      }else{}
    
    } # end of loop over rules.

############################################
### cubsti predictions are X %*% betaCbst...
############################################    
### divide by the number of rules applying for each row...    
    if(nRules > 1){

      normWithinComms <- TRUE

      if(!normWithinComms){
        nRulesPerRow <- rowSums(matRuleData)
        iTmp <- which(nRulesPerRow > 0)
        X[iTmp,] <- X[iTmp,,drop=FALSE] / matrix(nRulesPerRow[iTmp] , length(iTmp) , ncol(X))
      }else{
### within each committee...
        iRTmp <- 1
        for(iC in 1:cubistModel$committees){
          jThis <- which(cubistModel$comms4XCubist == iC)
          nRulesThis <- length(which(as.numeric(cubistModel$coefficients$committee) == iC))

          if(length(jThis) > 0 & nRulesThis > 0){
            matRuleDataThis <- matRuleData[,seq(iRTmp , iRTmp + nRulesThis - 1),drop=FALSE]
            nRulesPerRowThis <- rowSums(matRuleDataThis)
            iTmp <- which(nRulesPerRowThis > 0)
          
            X[iTmp,jThis] <- X[iTmp,jThis,drop=FALSE] / matrix(nRulesPerRowThis[iTmp] , length(iTmp) , length(jThis))

            iRTmp <- iRTmp + nRulesThis
          }else{}
        }   
        X <- X / length(unique(cubistModel$comms4XCubist))
      }
    }else{}    

###########################################################    
### if spline included, then cbind this...    
###########################################################    
    if(length(cubistModel$allKnotsd) > 0){
      
      XSpline <- allKnotsd2X(dIMidPts = dataFit$dIMidPts, allKnotsd = cubistModel$allKnotsd)
      X <- cbind(X , XSpline)

    }else{}
    
    return(list('X' = X , 'matRuleData' = matRuleData))

}

###############################################################
### fn to make X given knots and dIMidPts. 
### dIMidPts is a vector, evaluated as pt-support for these depths
### intercept is not included.
###############################################################
allKnotsd2X <- function(dIMidPts , allKnotsd){

  dIMidPts <- as.numeric(dIMidPts)
  if(length(allKnotsd) > 0){
    intKnots <- allKnotsd[-1]
    intKnots <- intKnots[-length(intKnots)]
    bdryKnots <- c(allKnotsd[1] , allKnotsd[length(allKnotsd)])
    
    XSpline <- bs(dIMidPts , knots = intKnots , degree = 3 , intercept = F , Boundary.knots = bdryKnots)
    colnames(XSpline) <- paste0('dSpline.' , seq(ncol(XSpline)))
  }else{
    XSpline <- matrix(NA , n , 0)
  }
  
  return(XSpline)
}

cubist2XSetup <- function(cubistModel , dataFit , zFit = NULL , profIDFit = NULL , allKnotsd = c() , removeColinCols = TRUE , refineCubistModel = FALSE){

  if(refineCubistModel & (!removeColinCols)){ stop('Error - if you want to refine the Cubist model, use removeColinCols = TRUE to make sure X is legal first.') }else{}
  
  if(length(allKnotsd) == 0){
    incdSpline <- FALSE
  }else if(length(allKnotsd) >= 2){
    incdSpline <- TRUE 
  }else{
    stop('Error - if entering knots, then enter at least 2 knots (2 boundary knots plus any internal knots)')
  }

##########################################################################
### check first row of data in cubistModel agrees with 1st row in dataFit
##########################################################################
    dataRow1 <- strsplit(cubistModel$data , split = '\n')
    nRowTmp <- length(dataRow1[[1]])
    dataRow1 <- strsplit(dataRow1[[1]][1] , split = ',')[[1]]
    dataRow1 <- dataRow1[-1] 
    
    if(length(dataRow1) != ncol(dataFit)){ 
      stop(paste0('Error - to run cubist2XSetup, dataFit must be the same as the data used to fit the cubist model! But they have different numbers of columns!' , 
                  '\nTo get prediction design matrix, first update cubistModel by running cubist2X with the fitting data!'))
    }else{}
    if(nRowTmp != nrow(dataFit)){ 
      stop(paste0('Error - to run cubist2XSetup, dataFit must be the same as the data used to fit the cubist model! But they have different numbers of rows!' , 
                  '\nTo get prediction design matrix, first update cubistModel by running cubist2X with the fitting data!'))
    }else{}
    for (j in 1:ncol(dataFit)){ 
      if(is.numeric(dataFit[1,j])){
#        checkThis <- dataFit[1,j] == as.numeric(dataRow1[j])
### only bother checking to 4 dp (rounding issues in above)...
### may not be sensible for some scales, but should be good enough to flag any big errors.
        checkThis <- abs(dataFit[1,j] - as.numeric(dataRow1[j])) < 1E-4
      }else{
        checkThis <- dataFit[1,j] == dataRow1[j]
      }      
      if(!checkThis){ 
        warning(paste0('Warning - to run cubist2XSetup, dataFit must be the same as the data used to fit the cubist model!' , 
                       '\nColumn ' , j , ' does not seem to agree, check it! (Though difference may be because of rounding, which would be ok.)')) 
      }else{} 
    } 
  
    nRules <- nrow(cubistModel$coefficients)
    dfCoeffs <- cubistModel$coefficients

#################################################
### check and remove the columns 'committee' and 'rule'
#################################################
    committeeCheck <- as.integer(dfCoeffs$committee)
    ruleCheck <- as.integer(dfCoeffs$rule)

    doChecks <- FALSE

    if(doChecks){
      if ((max(committeeCheck) == 1) & (min(committeeCheck) == 1)){
### all ok        
      }else{
        stop('More than 1 committees...not sure how to convert this to a design matrix')
      }
      if ((max(ruleCheck - seq(nRules)) == 0 ) & (min(ruleCheck - seq(nRules)) == 0 )){
### all ok        
      }else{
        stop('I think this should not happen - check that my code deals with rules not ordered!')
      }
    }else{}

    dfCoeffs$committee <- NULL
    dfCoeffs$rule <- NULL

    if(names(dfCoeffs)[1] != '(Intercept)'){ stop('Error - my coding assumes that intercept is the first column in dfCoeffs - check this!') }else{}

#################################################
### initialise listRules, listSplitCats and listCoeffs_iv...
### listCoeffs_iv is to quickly convert the variable names (col names of dfCoeffs without the intercept) into column numbers of the namesDataFit...
#################################################
    namesdataFit <- names(dataFit)
    listRules <- listSplitCats <- listCoeffs_iv <- vector("list" , length = nRules)

#################################################
### get the total number of linear parameters...
#################################################
    p <- length(which(!is.na(dfCoeffs)))
    pVec <- NA * integer(nRules)

    if (nRules > 1){
### a way of extracting splits, though still with rounding errors.
        vecSplitValues <- getCubistSplitsFromModel(cubistModel$model , dataFit)
    }else{
        vecSplitValues <- c()
    }
    
    for (iR in 1:nRules){

      iC <- as.numeric(cubistModel$coefficients$committee[iR])
      iRIniC <- as.numeric(cubistModel$coefficients$rule[iR])

      if(nRules > 1){
### get all the data that fall into rule iRIniC of this committee...        
### update for multiple committees, 28/02/19...        
        iTmp <- which(as.numeric(cubistModel$splits$rule) == iRIniC & as.numeric(cubistModel$splits$committee) == iC)
        
        splitsThis <- cubistModel$splits[iTmp,]
        splitValsThis <- vecSplitValues[iTmp]
        vNamesThis <- gsub("\"" , "" , splitsThis$variable) # gsub used to get rid of "", which somehow got added.

### number of variables conditions in conditions for this rule...
        listRules[[iR]] <- cubistModel$splits[iTmp,]
        nVThis <- nrow(listRules[[iR]])

        listRules[[iR]]$variable <- vNamesThis
### add which column number this is for quicker access...
        listRules[[iR]]$ivariable <- NA
        listRules[[iR]]$valUpdated <- splitValsThis
        
        listSplitCats[[iR]] <- vector("list" , length = nVThis)
        for (jS in 1:nVThis){

            ivThis <- which(namesdataFit == vNamesThis[jS])
            listRules[[iR]]$ivariable[jS] <- ivThis # order in the data frame.
            
            dirThis <- splitsThis$dir[jS]
            valThis <- splitValsThis[jS]
            
            if (dirThis == ''){
### categorical variable...so get all cats
                catsThis <- as.character(splitsThis$category[jS])
                tmp <-  gregexpr(',' , catsThis)
                startsTmp <- c(1 , tmp[[1]]+1)
                endsTmp <- c(tmp[[1]] - 1 , nchar(catsThis))
                nCatsThis <- length(tmp[[1]]) + 1
                vecCatsThis <- character(nCatsThis)
                for (iCat in 1:nCatsThis){
                    vecCatsThis[iCat] <- substr(catsThis , startsTmp[iCat] , endsTmp[iCat])
                }
                
                listSplitCats[[iR]][[jS]] <- vecCatsThis
            }else{}  
        }
        
      }else{
### only one rule, so no condits in it!
      }

      namesTmp <- names(dfCoeffs)[!is.na(dfCoeffs[iR,])]
      namesTmp <- namesTmp[-1]

      ivTmp <- NA * integer(length(namesTmp))
      for(j in 1:length(ivTmp)){
        ivTmp[j] <- which(namesdataFit == namesTmp[j])
      }      
      listCoeffs_iv[[iR]] <- ivTmp
      
      pThis <- length(which(!is.na(dfCoeffs[iR,])))
      pVec[iR] <- pThis
    
    } # end of loop over rules.

    if (nRules > 1){
        cubistModel$splits$valUpdated <- vecSplitValues
    }else{}

############################################
### cubsti predictions are X %*% betaCbst...
############################################
    betaCbst <- as.numeric(t(as.matrix(dfCoeffs)))
    iOK <- which(!is.na(betaCbst))
    betaCbst <- matrix(betaCbst[iOK] , ncol = 1)

    names4XCubist <- names(dfCoeffs)
    names4XCubist[1] <- 'const'
    names4XCubist <- rep(names4XCubist , nRules)
    names4XCubist <- paste0(names4XCubist , rep(paste0('_R' , seq(nRules)) , each = ncol(dfCoeffs)))

    comms4XCubist <- as.numeric(t(matrix(as.numeric(cubistModel$coefficients$committee) , nrow(dfCoeffs) , ncol(dfCoeffs))))
    
    names4XCubist <- names4XCubist[iOK]
    comms4XCubist <- comms4XCubist[iOK]

###########################################################################
#### if we have a d spline + d is a predictor in every rule, 
### then remove d as a predictor from the final rule of the cubist model...
### otherwise numberical errors (lin fn appears twice).
### CURRENTLY REMOVED FROM SETUP, BUT MAYBE EASIER TO TRACK IF STILL INCLUDED IN SETUP,
### BUT A LIST OF COLUMNS TO BE REMOVED FROM X IS SAVED.
###########################################################################
    if(incdSpline){
      if(is.element('dIMidPts' , names(cubistModel$coefficients)) && (length(which(is.na(cubistModel$coefficients$dIMidPts))) == 0)){

        dfCoeffs$dIMidPts[nRules] <- NA
        listCoeffs_iv[[nRules]] <- setdiff(listCoeffs_iv[[nRules]] , which(namesdataFit == 'dIMidPts'))

        ipRemoveFromCubist <- which(names4XCubist ==  paste0("dIMidPts_R" , nRules))    
        
        names4XCubist <- names4XCubist[-ipRemoveFromCubist]
        comms4XCubist <- comms4XCubist[-ipRemoveFromCubist] 
        betaCbst <- betaCbst[-ipRemoveFromCubist,,drop=FALSE]
        pVec[nRules] <- pVec[nRules] - 1
            
      }else{}         
    }else{}         

#######################################################
### add all useful info to the cubist model object...
#######################################################
    cubistModel$dfCoeffs <- dfCoeffs
    cubistModel$namesDataFit <- names(dataFit)
    cubistModel$listCoeffs_iv <- listCoeffs_iv
    cubistModel$listRules <- listRules
    cubistModel$listSplitCats <- listSplitCats
    cubistModel$pVec <- pVec 
    cubistModel$betaCbst <- betaCbst
    cubistModel$names4XCubist <- names4XCubist
    cubistModel$comms4XCubist <- comms4XCubist
    if(length(allKnotsd) > 0){
      cubistModel$allKnotsd <- allKnotsd
    }else{ # assigning as c() doesn't work, so: 
      cubistModel['allKnotsd'] <- list(NULL)
    }
    
#######################################################
### convert to X, and check for potential colinearity within rules...
### remove columns from X and the cubist set up variables if found.
#######################################################
    tmp <- cubist2XGivenSetup(cubistModel , dataFit)
    X <- tmp$X 
    matRuleData <- tmp$matRuleData

    if(any(rowSums(is.na(X)) > 0)){ stop('Some error in conversion of cubist model to X has produced NA values!') }else{}    

    if(removeColinCols){
      if(cubistModel$committees == 1){
        # ipRmvd <- getColinPreds(X = X)
        if(is.null(zFit)){ stop('Error - to use legalizeXIAK3D, must input zFit') }else{}
        ipRmvd <- legalizeXIAK3D(X = X , z = zFit)$ipRemove
        tmp <- removeColinPreds(X = X , cubistModel = cubistModel , dataFit = dataFit , ipRmvd = ipRmvd , reasonRmv = 0)
        X <- tmp$X 
        cubistModel <- tmp$cubistModel 

        if(refineCubistModel){
### further refine the model, by removing predictors with p > 0.15    
          if(is.null(profIDFit)){ stop('Error - to use refineXIAK3D, must input profIDFit') }else{}
          ipRmvd <- refineXIAK3D(X = X , z = zFit , profID = profIDFit , alpha = 0.15)$ipRemove
          tmp <- removeColinPreds(X = X , cubistModel = cubistModel , dataFit = dataFit , ipRmvd = ipRmvd , reasonRmv = 1)
          X <- tmp$X 
          cubistModel <- tmp$cubistModel 
        }else{}
        
      }else{

### look for any rule conditions that are duplicated...
#        listRulesTmp <- cubistModel$listRules
#        for (i in 1:length(listRulesTmp)){
#          listRulesTmp[[i]]$committee <- NULL
#          listRulesTmp[[i]]$rule <- NULL
#          listRulesTmp[[i]]$type <- NULL
#          listRulesTmp[[i]]$percentile <- NULL
#          listRulesTmp[[i]]$valUpdated <- NULL
#          oTmp <- order(listRulesTmp[[i]]$ivariable , listRulesTmp[[i]]$value)
#          listRulesTmp[[i]] <- listRulesTmp[[i]][oTmp,]
#        }
#        iDup <- which(duplicated(listRulesTmp))


### if multiple committees, remove colinearity from each committee first, then
### remove colinearity from combined model 
        iRTmp <- 1
        ipRmvd <- c()
        for(iC in 1:cubistModel$committees){
          jThis <- which(cubistModel$comms4XCubist == iC)
          nRulesThis <- length(which(as.numeric(cubistModel$coefficients$committee) == iC))
          iRThis <- seq(iRTmp , iRTmp + nRulesThis - 1)
        
          ipRmvdThis <- getColinPreds(X = X[,jThis,drop=FALSE])
        
          if(!is.null(ipRmvdThis)){ 
            ipRmvd <- c(ipRmvd , jThis[ipRmvdThis])
          }else{}
          
          iRTmp <- iRTmp + nRulesThis
        }

        tmp <- removeColinPreds(X = X , cubistModel = cubistModel , dataFit = dataFit , ipRmvd = ipRmvd)
        X <- tmp$X 
        cubistModel <- tmp$cubistModel 

        nRPerC <- NA * integer(cubistModel$committees)
        for(iC in 1:cubistModel$committees){
          nRPerC[iC] <- length(which(cubistModel$coefficients$committee == iC))
        }

### remove const from final rule of all comms > 1, due to perfect colinearity with sum of constants from comm = 1.
### though if final rule is const only, remove pre-ceding rule.
        ipRmvd <- c()
        for(iC in 2:cubistModel$committees){
          rTmp <- sum(nRPerC[1:iC])
          while(pVec[rTmp] <= 1){
            rTmp <- rTmp - 1
          }
          if(rTmp <= (sum(nRPerC[1:iC]) - nRPerC[iC])){ 
### all rules have just const. not sure if removing will work in my code.
            stop('All rules in this committee are constant only - not sure this will work!')           
          }else{}
          
          ipRmvd <- c(ipRmvd , which(cubistModel$names4XCubist == paste0('const_R' , rTmp)))
        }

        tmp <- removeColinPreds(X = X , cubistModel = cubistModel , dataFit = dataFit , ipRmvd = ipRmvd)
        X <- tmp$X 
        cubistModel <- tmp$cubistModel 

### now check full combined model, and if any colin exists, delete more...
        ipRmvd <- getColinPreds(X = X , comms4XCubist = cubistModel$comms4XCubist)
        tmp <- removeColinPreds(X = X , cubistModel = cubistModel , dataFit = dataFit , ipRmvd = ipRmvd)
        X <- tmp$X 
        cubistModel <- tmp$cubistModel 
      }

    }else{}

### add namesX to the cubistModel, which has the names4XCubist as well as the spline names.
    cubistModel$namesX <- colnames(X)    
      
    return(list('X' = X , 'cubistModel' = cubistModel , 'matRuleData' = matRuleData))

}

###########################################################################
### just get the indices of columns that can be removed from X to make XX invertible... 
###########################################################################
getColinPreds <- function(X , y = NULL , comms4XCubist = NULL){
    namesX <- colnames(X)
    ipRmvd <- namesXRmvd <- c()

    origColNumbers <- seq(ncol(X))
    XCurrent <- X
    XX <- t(XCurrent) %*% XCurrent
    eXX <- eigen(XX)
    iXX <- try(solve(XX) , silent = TRUE)

### overkill using ncol(X) here, but safe.    
#    isConstTmp <- is.element(namesX , paste0('const_R' , seq(ncol(X))))
#    isConstTmp <- is.element(namesX , paste0('constjibberjabber_R' , seq(ncol(X))))

### get which columns of X are the only column of that rule...
### and don't allow those columns to be removed.
    tmp <- strsplit(namesX , '_R')
    ruleNumbers <- NA * integer(length(tmp))
    for(i in 1:length(tmp)){
      ruleNumbers[i] <- as.numeric(tmp[[i]][length(tmp[[i]])])
    }
    
    tblTmp <- table(ruleNumbers)
    iTmp <- which(tblTmp == 1)
    rulesKeepTmp <- as.integer(rownames(tblTmp)[iTmp])
    colsKeepTmp <- which(is.element(ruleNumbers , rulesKeepTmp))
    
    while(is.character(iXX)){
      jMin <- which.min(eXX$value)
      evec <- eXX$vectors[,jMin]

      evec[colsKeepTmp] <- 0 # set these to 0 so they won't be removed.

      if(is.null(comms4XCubist)){
        ipRmvFromCurrent <- which.max(abs(evec))
      }else{
# take the top 5 contributing columns, get the subset of these which have the largest comm number
# then the max abs(evec) value of these
        iPoss <- order(-abs(evec))
        if(length(iPoss) > 5){
          iPoss <- iPoss[1:5]
        }else{}
        commsPoss <- comms4XCubist[iPoss]
        iPoss <- iPoss[which(commsPoss == max(commsPoss))]
        ipRmvFromCurrent <- iPoss[which.max(abs(evec[iPoss]))] 
      }    

      ipRmvThis <- origColNumbers[ipRmvFromCurrent]
      ipRmvd <- c(ipRmvd , ipRmvThis)
      namesXRmvd <- c(namesXRmvd , namesX[ipRmvThis])
      comms4XCubist <- comms4XCubist[-ipRmvFromCurrent]
      
      origColNumbers <- origColNumbers[-ipRmvFromCurrent]
      ruleNumbers <- ruleNumbers[-ipRmvFromCurrent]
      
      tblTmp <- table(ruleNumbers)
      iTmp <- which(tblTmp == 1)
      rulesKeepTmp <- as.integer(rownames(tblTmp)[iTmp])
      colsKeepTmp <- which(is.element(ruleNumbers , rulesKeepTmp))
      
      XCurrent <- XCurrent[,-ipRmvFromCurrent,drop=FALSE]
      XX <- XX[-ipRmvFromCurrent,-ipRmvFromCurrent,drop=FALSE]

      eXX <- eigen(XX)
      iXX <- try(solve(XX) , silent = TRUE)
    }
    
    return(ipRmvd)
}


removeColinPreds <- function(X , cubistModel , dataFit , ipRmvd , reasonRmv = 0){

    n <- nrow(X)
    
### if there are any NAs, remove...
    nRules <- length(cubistModel$listRules)

    namesX <- colnames(X)

#####################################################################
### update the things in the cubistModel in light of the removed parameters...
#####################################################################
    if(length(ipRmvd) > 0){

      cubistModel$names4XCubistRmvd <- c(cubistModel$names4XCubistRmvd , namesX[ipRmvd])
      cubistModel$betaCbst <- cubistModel$betaCbst[-ipRmvd,,drop=FALSE] 
      cubistModel$names4XCubist <- cubistModel$names4XCubist[-ipRmvd] 
      cubistModel$comms4XCubist <- cubistModel$comms4XCubist[-ipRmvd] 

      if(reasonRmv == 0){
        print('Some colinear columns found in X - legalizing X by removing the following columns:')
      }else if(reasonRmv == 1){
        print('Some redundant columns found in X - refining X by removing the following columns:')
      }else{
        stop('Enter valid reasonRmv value to removeColinPreds (0 or 1)!')
      }
      for(i in 1:length(ipRmvd)){ print(paste0(i , '.  ' , namesX[ipRmvd[i]])) }

      X <- X[,-ipRmvd,drop=FALSE]
    
      rulesBad <- colsBad <- colsBad2 <- NA * integer(length(ipRmvd))
      for (i in 1:length(ipRmvd)){
        badTmp <- strsplit(namesX[ipRmvd[i]] , '_R')
        rulesBad[i] <- as.numeric(badTmp[[1]][length(badTmp[[1]])]) # rule numbers in my notation (cts over all committees)
        cubistModel$pVec[rulesBad[i]] <- cubistModel$pVec[rulesBad[i]] - 1
        if(rulesBad[i] < 10){ 
          nameBad <- substr(namesX[ipRmvd[i]] , 1 , nchar(namesX[ipRmvd[i]]) - 3)
        }else if(rulesBad[i] < 100){ 
          nameBad <- substr(namesX[ipRmvd[i]] , 1 , nchar(namesX[ipRmvd[i]]) - 4)
        }else if(rulesBad[i] < 1000){ 
          nameBad <- substr(namesX[ipRmvd[i]] , 1 , nchar(namesX[ipRmvd[i]]) - 5)
        }else{        
          nameBad <- substr(namesX[ipRmvd[i]] , 1 , nchar(namesX[ipRmvd[i]]) - 6)
        }
### first within dfCoeffs...
        if(nameBad == 'const'){ nameBad <- '(Intercept)' }else{}
        jBad <- which(names(cubistModel$dfCoeffs) == nameBad)
        if(length(jBad) == 0){
          stop(paste0('Error - the name of the bad predictor (' , nameBad , ') has not been found in dfCoeffs!'))
        }else if(length(jBad) > 1){
          stop(paste0('Error - multiple names of the bad predictor (' , nameBad , ') have been found in dfCoeffs!'))
        }else{}
        colsBad[i] <- jBad
### second within dataFit...
        if(nameBad == '(Intercept)'){
          jBad <- NA
        }else{
          jBad <- which(names(dataFit) == nameBad)
        }
        if(length(jBad) == 0){
          stop(paste0('Error - the name of the bad predictor (' , nameBad , ') has not been found in dataFit!'))
        }else if(length(jBad) > 1){
          stop(paste0('Error - multiple names of the bad predictor (' , nameBad , ') have been found in dataFit!'))
        }else{}
        colsBad2[i] <- jBad
      }
      for(i in 1:length(rulesBad)){ 
        cubistModel$dfCoeffs[rulesBad[i] , colsBad[i]] <- NA 
      }
      for(i in 1:length(rulesBad)){ 
        if(!is.na(colsBad2[i])){
          tmp <- cubistModel$listCoeffs_iv[[rulesBad[i]]]
          itmp <- which(tmp == colsBad2[i])
          if(length(itmp) == 1){
            tmp <- tmp[-itmp]
            cubistModel$listCoeffs_iv[[rulesBad[i]]] <- tmp
          }else if(length(itmp) == 0){
            stop('Error - not found a removed variable in listCoeffs_iv!')
          }else{
            stop('Error - found multiple occurences of a removed variable in listCoeffs_iv!')
          }
        }else{}
      }
      for(i in 1:nRules){ cubistModel$listCoeffs_iv[[i]] <- cubistModel$listCoeffs_iv[[i]][which(!is.na(cubistModel$listCoeffs_iv[[i]]))] }

    }else{}

    return(list('X' = X , 'cubistModel' = cubistModel))
}

getCubistSplitsFromModel <- function(x , dataFit = NULL)
  {

### argument x is cubistModel$model
### dataFit (if given) is the data that were used to fit model
### if given, snap any splits to nearest data point if one within 1E-4
  
    x <- strsplit(x, "\n")[[1]]
    comNum <- ruleNum <- condNum <- rep(NA, length(x))
    comIdx <- rIdx <- 0
    for(i in seq(along = x))
      {
        tt <- parserTmp(x[i])
        ## Start of a new rule
        if(names(tt)[1] == "rules")
          {
            comIdx <- comIdx + 1
            rIdx <- 0
          }
        comNum[i] <-comIdx
        ## Start of a new condition
        if(names(tt)[1] == "conds")
          {
            rIdx <- rIdx + 1
            cIdx <- 0
          }
        ruleNum[i] <-rIdx
        ## Within a rule, type designates the type of conditional statement
        ## type = 2 appears to be a simple split of a continuous predictor
        if(names(tt)[1] == "type")
          {
            cIdx <- cIdx + 1
            condNum[i] <- cIdx
          }
      }
    
    numCom <- sum(grepl("^rules=", x))
    rulesPerCom <- unlist(lapply(split(ruleNum, as.factor(comNum)), max))
    rulesPerCom <- rulesPerCom[rulesPerCom > 0]
    if (! is.null(rulesPerCom) && numCom > 0)
      names(rulesPerCom) <- paste("Com", 1:numCom)

    ## In object x, what element starts a new rule
    isNewRule <- ifelse(grepl("^conds=", x), TRUE, FALSE)   
    splitVar <- rep("", length(x))
    splitVal <- rep(NA, length(x))
    splitCats <- rep("", length(x))
    splitDir <- rep("", length(x))

    ## This is a simple continuous split, such as
    ##
    ##   nox > 0.668
    ##
    ## or
    ##
    ## type="2" att="nox" cut="0.66799998" result=">"
    ##    
    isType2 <- grepl("^type=\"2\"", x)
    
    if(any(isType2))
      {
        splitVar[isType2] <- type2Tmp(x[isType2])$var
        splitVar[isType2] <- gsub("\"", "", splitVar[isType2])
        splitDir[isType2] <- type2Tmp(x[isType2])$rslt
        splitVal[isType2] <- type2Tmp(x[isType2])$val
      }
    ## This is a split of categorical data such as 
    ##
    ##   X4 in {c, d}
    ##
    ## or
    ##
    ## type="3" att="X4" elts="c","d"
    ##

    isType3 <- grepl("^type=\"3\"", x)
    if(any(isType3))
      {
        splitVar[isType3] <- type3Tmp(x[isType3])$var
        splitCats[isType3] <- type3Tmp(x[isType3])$val
        splitCats[isType3] <- gsub("[{}]", "", splitCats[isType3])
        splitCats[isType3] <- gsub("\"", "", splitCats[isType3])
        splitCats[isType3] <- gsub(" ", "", splitCats[isType3])
      }
    if(!any(isType2) & !any(isType3)) return(NULL)
    splitData <- data.frame(committee = comNum,
                            rule = ruleNum,
                            variable = splitVar,
                            dir = splitDir,
                            value = as.numeric(splitVal),
                            category = splitCats)
    splitData$type <- ""
    if(any(isType2)) splitData$type[isType2] <- "type2"
    if(any(isType3)) splitData$type[isType3] <- "type3"    
    splitData <- splitData[splitData$variable != "" ,]
    splitData


    if(!is.null(dataFit)){
### snap cts splits to nearest datapoint if within 1E-4 (as some rounding errors in display/table)    
### better way might be to round appropriately before fitting cubist model
      for(i in 1:length(splitData$val)){
        if(splitData$type[i] == 'type2'){
          zThis <- dataFit[[as.character(splitData$variable[i])]]
          iTmp <- which.min(abs(zThis - splitData$val[i]))
          if(abs(zThis[iTmp] - splitData$val[i]) < 1E-4){
            splitData$val[i] <- zThis[iTmp]
          }else{}
        }else{}
      }
    }else{}      
        
    return(splitData$val)
  }

type3Tmp <- function(x)
  {
    aInd <- regexpr("att=", x)
    eInd <- regexpr("elts=", x)
    var <- substring(x, aInd + 4, eInd - 2)
    val <- substring(x, eInd + 5)
    multVals <- grepl(",", val)
    val <- gsub(",", ", ", val) 
    val <- ifelse(multVals, paste("{", val, "}", sep = ""), val)
    txt <- ifelse(multVals,  paste(var, "in", val),  paste(var, "=", val))
 
    list(var = var, val = val, text = txt)
  }

type2Tmp <- function(x, dig = 16)
  {
    x <- gsub("\"", "", x)
    aInd <- regexpr("att=", x)
    cInd <- regexpr("cut=", x)
    rInd <- regexpr("result=", x)
    vInd <- regexpr("val=", x)

    var <- val <- rslt <- rep("", length(x))
    
    missingRule <- cInd < 1 & vInd > 0
    
    if(any(missingRule))
      {
        var[missingRule] <- substring(x[missingRule], aInd[missingRule] + 4, vInd[missingRule] - 2)
        val[missingRule] <- "NA"
        rslt[missingRule] <- "="  
      }
    if(any(!missingRule))
      {
        var[!missingRule] <- substring(x[!missingRule], aInd[!missingRule] + 4, cInd[!missingRule] - 2)        
        val[!missingRule] <- substring(x[!missingRule], cInd[!missingRule] + 4, rInd[!missingRule] - 1)
        val[!missingRule] <- format(as.numeric(val[!missingRule]), digits = dig)
        rslt[!missingRule] <- substring(x[!missingRule], rInd[!missingRule] + 7)
      }

    list(var = var, val = as.numeric(val), rslt = rslt,
         text = paste(var, rslt, val))
  }


parserTmp <- function(x)
  {
    x <- strsplit(x, " ")
    x <- lapply(x,
                function(y)
                {
                  y <- strsplit(y, "=")
                  nms <- unlist(lapply(y, function(z) z[1]))
                  val <- unlist(lapply(y, function(z) z[2]))
                  names(val) <- nms
                  val
                })
    if(length(x) == 1) x <- x[[1]]
    x
  }


##############################################################
### function to get all d Breaks from cubist model (which has already had cubist2X run on it).
##############################################################
getAlldBreaks <- function(cubistModel){
  if(!is.element('listRules' , names(cubistModel))){ stop('Error - run cubist2X function before getting all dBreaks!') }else{}
  dfdBreaks <- lapply(cubistModel$listRules , function(dfIn){ dfIn[which(dfIn$variable == 'dIMidPts'),,drop=FALSE] })
  dfdBreaks <- do.call(rbind , dfdBreaks)
  if(!is.null(dfdBreaks)){
    dBreaks <- unique(dfdBreaks$valUpdated)
    dBreaks <- dBreaks[order(dBreaks)]
  }else{
    dBreaks <- c()
  }
  return(dBreaks)
}

##############################################################
### function to get the indices of the rules with d in the rule conditions from cubist model (which has already had cubist2X run on it).
##############################################################
getRulesWithdInCondits <- function(cubistModel){
  if(!is.element('listRules' , names(cubistModel))){ stop('Error - run cubist2X function before getting all dBreaks!') }else{}
  
  ndBreaksPerRule <- unlist(lapply(cubistModel$listRules , function(dfIn){ length(which(dfIn$variable == 'dIMidPts')) }))
  
  return(which(ndBreaksPerRule > 0))
}

##############################################################
### function to rule numbers from X column names from cubist model (which has already had cubist2X run on it).
##############################################################
getRuleNumbersFromColNames <- function(cubistModel){
  if(!is.element('listRules' , names(cubistModel))){ stop('Error - run cubist2X function before getting all dBreaks!') }else{}
  
### split by '_'
  tmp <- strsplit(cubistModel$names4XCubist , '_')
### get final elements...  
  tmp <- unlist(lapply(tmp, tail , n = 1L))
### remove 'R' and convert to integer...    
  jRules <- as.integer(substr(tmp , 2 , nchar(tmp)))

  return(jRules)
}

##############################################################
### function to remove variables from X to make inv(XX) exist.
### does this by fitting models if possible with 1 column removed
### if any are legal, select the one with smallest nll
### else continue looking at all models with 2 columns removed
### if 2 doesn't work, function will stop and return error
### could continue but will get more time consuming - 
### so either continue or write something different.
##############################################################
legalizeXIAK3D <- function(X , z){
  
  XX <- t(X) %*% X
  Xz <- t(X) %*% z 
  
  n <- length(z) 
  
#  tmp0 <- try(solve(XX , Xz) , silent = TRUE)
#  tmp <- lndetANDinvCb(XX, Xz)

  iXXXz <- try(solve(XX , Xz) , silent = TRUE)

  nCols2Rmv <- 0
#  continueRemoving <- (is.na(tmp$lndet) | is.character(tmp0))
  continueRemoving <- (is.character(iXXXz))
  ipRemove <- c()
  while(continueRemoving){
    
    nCols2Rmv <- nCols2Rmv + 1
    
    cands <- setdiff(seq(ncol(X)) , which(is.element(colnames(X) , paste0('const_R' , seq(100))))) 
    cands <- setdiff(cands , which(colnames(X) == 'dSpline'))
    
    if (nCols2Rmv == 1){
      nllTest <- NA * numeric(ncol(X))
      for (j in cands){
        nllTest[j] <- nllLm(z = z , X = X[,-j,drop=FALSE] , REML = F)$nll
      }
      if(all(is.na(nllTest))){
        continueRemoving <- TRUE
      }else{
        ipRemove <- which.min(nllTest)
        namesRemove <- colnames(X)[ipRemove]
        continueRemoving <- FALSE
      }
    }else if (nCols2Rmv == 2){
      nllTest <- matrix(NA , ncol(X) , ncol(X))
      for (j in cands){
        for(k in cands){
          if(k>j){
            nllTest[j,k] <- nllLm(z = z , X = X[,-c(j,k),drop=FALSE] , REML = F)$nll
          }else{}
        }
      }
      if(all(is.na(nllTest))){
        continueRemoving <- TRUE
#        stop('Generalise this algorithm or write another to better remove colinear predictors when >2 need removing!')
      }else{
        ipRemove <- which(nllTest == min(nllTest , na.rm = TRUE) , arr.ind = TRUE)
        ipRemove <- as.numeric(ipRemove[1,])
        namesRemove <- colnames(X)[ipRemove]
        continueRemoving <- FALSE
      }
      
    }else if (nCols2Rmv == 3){
      nllTest <- array(NA , dim = c(ncol(X) , ncol(X) , ncol(X)))
      for (j in cands){
        for(k in cands){
          for(l in cands){
            if(k>j & l>k){
              nllTest[j,k,l] <- nllLm(z = z , X = X[,-c(j,k,l),drop=FALSE] , REML = F)$nll
            }else{}
          }
        }
      }
      if(all(is.na(nllTest))){
        continueRemoving <- TRUE
        stop('Generalise this algorithm or write another to better remove colinear predictors when >3 need removing!')
      }else{
        ipRemove <- which(nllTest == min(nllTest , na.rm = TRUE) , arr.ind = TRUE)
        ipRemove <- as.numeric(ipRemove[1,])
        namesRemove <- colnames(X)[ipRemove]
        continueRemoving <- FALSE
      }
      
    }else{
      stop('Generalise this algorithm or write another to better remove colinear predictors when >3 need removing!')
    }
  }
  
  if(length(ipRemove) > 0){
    namesRemove <- colnames(X)[ipRemove]
    X <- X[,-ipRemove,drop=FALSE]
  }else{
    namesRemove <- c()
  }
  
  return(list('X' = X , 'namesRemove' = namesRemove , 'ipRemove' = ipRemove))
}

#######################################################################
### function to remove variables from X by Wald tests. Assumes inv(XX) exists.
### fits models with a profile-specific random effect
### uses benj-hoch with all p values and the given alpha to see if all variables are significant
### if not, the least significant is removed
### continues until all are significant
### alpha = 0.5 intended to give a model that one would expect to be better at each stage.
#######################################################################
refineXIAK3D <- function(X , z , profID , alpha = 0.5){

  continueRefining <- TRUE
  n <- nrow(X)
  p <- ncol(X)
  colnamesXIn <- colnames(X)

  XIn <- X
### rescale the candidate columns...
  cands <- setdiff(seq(ncol(X)) , which(is.element(colnames(X) , paste0('const_R' , seq(100))))) 
  cands <- setdiff(cands , which(substr(colnames(X) , 1 , 8) == 'dSpline.'))

  sdX <- sqrt(apply(X,2,var))
  X[,cands] <- X[,cands,drop=FALSE]  / matrix(sdX[cands] , nrow(X) , length(cands) , byrow = TRUE)
  
  
  namesRemove <- c()
  while(continueRefining){

    formulaTmp <- paste0('z ~ 0 + ' , paste(colnames(X) , collapse = ' + ') , ' + (1|profID)')
    df4lmer <- data.frame(X , 'z' = z , 'profID' = profID)
    ft <- lmer(as.formula(formulaTmp) , data = df4lmer)

    betahat <- summary(ft)$coefficients[,1]
    vbetahat <- summary(ft)$vcov

    cands <- setdiff(seq(ncol(X)) , which(is.element(colnames(X) , paste0('const_R' , seq(100))))) 
    cands <- setdiff(cands , which(substr(colnames(X) , 1 , 8) == 'dSpline.'))

    wStats <- NA * numeric(p)
    for(i in cands){
      wStats[i] <- waldTest(betahat = betahat , vbetahat = vbetahat , ip0 = i)
    }

### are all of the variables in the mlrs significant?
### if not, remove the one with the largest p value.
    iOrder <- order(wStats)
    nTests <- length(which(!is.na(wStats)))
    wStatsOrdered <- wStats[iOrder[1:nTests]]
    alphaMHT <- seq(alpha/nTests , alpha , alpha/nTests) 
  
    if(max(wStatsOrdered - alphaMHT) >= 0){
### something not significant, so drop the least significant variable...
      iRemoveThis <- which.max(wStats)
      namesRemove <- c(namesRemove , colnames(X)[iRemoveThis])
      X <- X[,-iRemoveThis,drop=FALSE]
      p <- ncol(X)
    }else{
### everything significant, so stop dropping.
      continueRefining <- FALSE
    }
  }
  
  if(length(namesRemove) > 0){
    ipRemove <- NA * integer(length(namesRemove))
    for(j in 1:length(namesRemove)){
      ipRemove[j] <- which(colnamesXIn == namesRemove[j])
    }
### X was stdzd, so returning to the unstandardized version and remoivng cols...
    X <- XIn[,-ipRemove,drop=FALSE]
  }else{
    ipRemove <- c()
  }

  return(list('X' = X , 'namesRemove' = namesRemove , 'ipRemove' = ipRemove))
}

#############################################################################
### function to make profID from coordinates...
#############################################################################
makeProfID <- function(cAll , useOldVersion = TRUE){
  if(is.null(cAll)){
    return(NULL)
  }else{
    # if(!is.element(class(cAll) , c('data.frame' , 'matrix'))){ stop('Error - makeProfID assumes cAll entered as matrix or data.frame.') }else{}
    if(!(is(cAll ,'data.frame') | is(cAll , 'matrix'))){ stop('Error - makeProfID assumes cAll entered as matrix or data.frame.') }else{}
    
    ndim <- ncol(cAll)
    
    if(useOldVersion){
      cU <- cAll[which(!duplicated(cAll)),,drop=FALSE]
      profID <- NA * numeric(nrow(cAll))
      for (i in 1:nrow(cU)){
        if(ndim == 1){
          iThis <- which(cAll[,1] == cU[i,1])
        }else if(ndim == 2){
          iThis <- which(cAll[,1] == cU[i,1] & cAll[,2] == cU[i,2])
        }else if(ndim == 3){
          iThis <- which(cAll[,1] == cU[i,1] & cAll[,2] == cU[i,2] & cAll[,3] == cU[i,3])
        }else{
          stop('Really - more than 3 dims for making ids?')
        }
        profID[iThis] <- i
      }
      
    }else{
      ### or would be better...
      if(ndim == 1){
        profID <- as.character(cAll[,1])
      }else if(ndim == 2){
        profID <- paste0(as.character(cAll[,1]) , '_' , as.character(cAll[,2]))
      }else if(ndim == 3){
        profID <- paste0(as.character(cAll[,1]) , '_' , as.character(cAll[,2]) , '_' , as.character(cAll[,3]))
      }else{
        stop('Really - more than 3 dims for making ids?')
      }
### could convert to numeric id if rqd...      
      # profID <- as.numeric(factor(profID))
      
    }

    return(profID)
  }
}

#############################################################################
### function to do a xv routine with cFit to select best nRules and whether refineCubistModel = T/F is best
### xv is with full profiles removed/kept.
#############################################################################
selectCubistOptsXV <- function(cFit , zFit , covsFit , allKnotsd = c() , nRulesVec = seq(10) , refineCubistModelVec = c(FALSE , TRUE)){

### proportion of cFit profiles used to fit, the rest used as 'val' data
  prop4XVFit <- 0.7

### number of reps of this xval.  
  nXVReps <- 10

### use cFit to define profIDFit...    
  profIDFit <- makeProfID(cFit)
  
  profIDUFit <- unique(profIDFit)
  nProfsFit <- length(profIDUFit)
  nProfsXVFit <- ceiling(prop4XVFit * nProfsFit)
  nProfsXVVal <- nProfsFit - nProfsXVFit
  
  rmseMatList <- list()
  for(i in 1:length(refineCubistModelVec)){
    rmseMatList[[i]] <- matrix(NA , length(nRulesVec) , nXVReps)
  }

### for calculating a weighted average sqd err at some stage     
  # Kd <- sparseMatrix(i = seq(nrow(cFit)) , j = as.integer(factor(profIDFit)) , x = 1)
  # KdKd <- t(Kd) %*% Kd
  # iKdKd <- sparseMatrix(i = seq(nrow(KdKd)) , j = seq(nrow(KdKd)) , x = (1/diag(KdKd)))
  
### fit and store all initial cubist models using all Fit data to get actual predictors with each setting...
  XFitList <- list()
  for(jrc in 1:length(refineCubistModelVec)){
    refineCubistModel <- refineCubistModelVec[jrc]
    XFitList[[jrc]] <- list()
    for(inRules in 1:length(nRulesVec)){
      nRules <- nRulesVec[inRules]
      
      cmFit <- cubist(x = covsFit , y = zFit , committees = 1 , cubistControl(rules = nRules))
      
      ### convert to des mtx
      tmp <- cubist2X(cubistModel = cmFit, dataFit = covsFit , zFit = zFit , profIDFit = profIDFit , allKnotsd = allKnotsd , refineCubistModel = refineCubistModel)
      cmFit <- tmp$cubistModel
      XFitThis <- tmp$X
      
      cands <- setdiff(seq(ncol(XFitThis)) , which(is.element(colnames(XFitThis) , paste0('const_R' , seq(100))))) 
      cands <- setdiff(cands , which(substr(colnames(XFitThis) , 1 , 8) == 'dSpline.'))
      
      sdX <- sqrt(apply(XFitThis,2,var))
      XFitThis[,cands] <- XFitThis[,cands,drop=FALSE]  / matrix(sdX[cands] , nrow(XFitThis) , length(cands) , byrow = TRUE)
      
      XFitList[[jrc]][[inRules]] <- XFitThis
    }
  }
      
### do xv...  
  warningFlagFitList <- list()
  for(jrc in 1:length(refineCubistModelVec)){
    warningFlagFitList[[jrc]] <- matrix(0 ,  length(nRulesVec) , nXVReps)
  }
  for (ixv in 1:nXVReps){
    
    ### split profiles into fit/val...      
    iFitThis <- sample(profIDUFit , nProfsXVFit)
    iFitThis <- iFitThis[order(iFitThis)]
    iValThis <- setdiff(profIDUFit , iFitThis)
    
    ### get indices for horizon data...
    iFitThis <- which(is.element(profIDFit , iFitThis))
    iValThis <- which(is.element(profIDFit , iValThis))
    
    ### split data...                                                         
    zFitThis <- zFit[iFitThis]
    profIDFitThis <- profIDFit[iFitThis]
    
    zValThis <- zFit[iValThis]
    profIDValThis <- profIDFit[iValThis]
    
    for(jrc in 1:length(refineCubistModelVec)){
      refineCubistModel <- refineCubistModelVec[jrc]
      for(inRules in 1:length(nRulesVec)){
        XFitThis <- XFitList[[jrc]][[inRules]][iFitThis,,drop=FALSE]
        XValThis <- XFitList[[jrc]][[inRules]][iValThis,,drop=FALSE]
        
        formulaTmp <- paste0('z ~ 0 + ' , paste(colnames(XFitThis) , collapse = ' + ') , ' + (1|profID)')
        df4lmerFit <- data.frame(XFitThis , 'z' = zFitThis , 'profID' = profIDFitThis)
        lmerFt <- lmer(as.formula(formulaTmp) , data = df4lmerFit , REML = TRUE)
        
        betahatThis <- summary(lmerFt)$coefficients[,1]
        namesInc <- rownames(summary(lmerFt)$coefficients)
        if(length(namesInc) < ncol(XValThis)){
          iDrop <- which(!is.element(colnames(XValThis) , namesInc))
          XValThis <- XValThis[,-iDrop,drop=FALSE]
          warningFlagFitList[[jrc]][inRules,ixv] <- 1
        }else{}
        
        zkValThis <- XValThis %*% betahatThis
        
        errsThis <- zkValThis - zValThis
        rmseMatList[[jrc]][inRules , ixv] <- sqrt(mean(errsThis ^ 2))
      }
    }
  }
  
  if(length(refineCubistModelVec) == 1){
    refineCubistModel <- refineCubistModelVec[1]
    rmrmse <- rowMeans(rmseMatList[[1]])
    iFlagged <- which(rowSums(warningFlagFitList[[1]]) > 0)
    rmrmse[iFlagged] <- Inf
    rmrmseBest <- min(rmrmse)
    
    if(is.infinite(rmrmseBest)){
      nRules <- NA      
    }else{
      nRules <- nRulesVec[which.min(rmrmse)]
    }

  }else{
    rmrmse1 <- rowMeans(rmseMatList[[1]])
    rmrmse2 <- rowMeans(rmseMatList[[2]])
    
### only consider the ones that didn't ever need any variables removing...    
    iFlagged1 <- which(rowSums(warningFlagFitList[[1]]) > 0)
    iFlagged2 <- which(rowSums(warningFlagFitList[[2]]) > 0)
    rmrmse1[iFlagged1] <- Inf
    rmrmse2[iFlagged2] <- Inf

    rmrmse1Best <- min(rmrmse1)
    rmrmse2Best <- min(rmrmse2)
    
    if(is.infinite(rmrmse1Best) & is.infinite(rmrmse2Best)){
      refineCubistModel <- NA
      nRules <- NA
    }else{
      if(rmrmse1Best < rmrmse2Best){
        refineCubistModel <- refineCubistModelVec[1]
        nRules <- nRulesVec[which.min(rmrmse1)]
      }else{
        refineCubistModel <- refineCubistModelVec[2]
        nRules <- nRulesVec[which.min(rmrmse2)]
      }
    }
    
  }

  return(list('nRules' = nRules , 'refineCubistModel' = refineCubistModel , 'rmseMatList' = rmseMatList , 'warningFlagFitList' = warningFlagFitList))
  
}



