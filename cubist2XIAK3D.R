cubist2X <- function(cubistModel , dataFit , incdSpline = FALSE , removeColinCols = TRUE){

    if(is.null(cubistModel$listRules) || is.null(cubistModel$incdSpline) || (cubistModel$incdSpline != incdSpline)){
    
      tmp <- cubist2XSetup(cubistModel = cubistModel , dataFit = dataFit , incdSpline = incdSpline , removeColinCols = removeColinCols)
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
      
      tmp <- cubist2XGivenSetup(cubistModel , dataFit , incdSpline = FALSE)
      X <- tmp$X 
      matRuleData <- tmp$matRuleData
    }

    iNA <- which(rowSums(is.na(X) | is.nan(X)) > 0)
    if(length(iNA) > 0){
      X[iNA,] <- NA
    }else{}

    return(list('X' = X , 'matRuleData' = matRuleData , 'cubistModel' = cubistModel))

}

cubist2XGivenSetup <- function(cubistModel , dataFit , incdSpline = FALSE){

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

#      print(head(matRuleData))

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

    return(list('X' = X , 'matRuleData' = matRuleData))

}

cubist2XSetup <- function(cubistModel , dataFit , incdSpline = FALSE , removeColinCols = TRUE){

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
    cubistModel$incdSpline <- incdSpline
    
#######################################################
### convert to X, and check for potential colinearity within rules...
### remove columns from X and the cubist set up variables if found.
#######################################################
    tmp <- cubist2XGivenSetup(cubistModel , dataFit , incdSpline = FALSE)
    X <- tmp$X 
    matRuleData <- tmp$matRuleData

    if(any(rowSums(is.na(X)) > 0)){ stop('Some error in conversion of cubist model to X has produced NA values!') }else{}    

    if(removeColinCols){
      if(cubistModel$committees == 1){
        ipRmvd <- getColinPreds(X = X)
        tmp <- removeColinPreds(X = X , cubistModel = cubistModel , dataFit = dataFit , ipRmvd = ipRmvd)
        X <- tmp$X 
        cubistModel <- tmp$cubistModel 
    
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


removeColinPreds <- function(X , cubistModel , dataFit , ipRmvd){

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

      print('Some colinear columns found in X - removing the following columns:')
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


