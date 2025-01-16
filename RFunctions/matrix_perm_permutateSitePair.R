####################################
#
# PERMUTATE GDM TABLE
# in mixed-effects model gdm
#
####################################
# 02.02.2022 Noelle
# script contains 2 functions to perfrom matrix permutation by group
#    - permutateSitePair_edit()
#    - permutateVarSitePair_edit()
# analogously to the functions permutateSitePair() and permutateVarSitePair()
#    in the gdm package.


####################################
# Edit permutatesitepair
####################################
# spTab <- currSitePair
# # spTab <- currSitePair[sample(1:nrow(currSitePair), 50),]
# siteVarTab <- siteData
# vNames <- varNames
# indexTab
# IndPermGroup

permutateSitePair_edit <- function(spTab, siteVarTab, indexTab, vNames, IndPermGroup, includeYear, 
                                   spacefortimedataset, nameToRemove){
  # print("Using the edited matrix permutation.")
  # permute the siteVarTable
  # need id of permgroup
  names(IndPermGroup) <- sub("s1.", "", names(IndPermGroup))
  if(spacefortimedataset == "temporal"){
    torand <- data.table(merge(siteVarTab, IndPermGroup, by = c("xCoord", "yCoord")))}
  if(spacefortimedataset == "spatial"){
    torand <- data.table(merge(siteVarTab, IndPermGroup, by = c("Year")))}
  
  for(g in unique(torand$permgroup)){
    # permute every group individually
    temp <- torand[permgroup == g, ]
    temp <- temp[sample(nrow(temp), nrow(temp), replace = F)]
    torand[permgroup == g, ] <- temp
    rm(temp)
  }
  torand[, permgroup := NULL]
  randVarTab <- data.frame(torand)
  rm(torand)
  class(randVarTab) <- class(siteVarTab)
  
  ##### next paragraph : fill permuted data back into spTable
  if(includeYear == F){
    # fill permuted data back into spTable "spTab"
    s1xCoord <- sapply(1:nrow(spTab), function(i) {
      randVarTab[indexTab[i, 1], 1]
    })
    s1yCoord <- sapply(1:nrow(spTab), function(i) {
      randVarTab[indexTab[i, 1], 2]
    })
    s2xCoord <- sapply(1:nrow(spTab), function(i) {
      randVarTab[indexTab[i, 2], 1]
    })
    s2yCoord <- sapply(1:nrow(spTab), function(i) {
      randVarTab[indexTab[i, 2], 2]
    })
    varLists <- lapply(vNames, function(vn, rvTab, spt, inT) {
      if (vn != "Geographic") {
        randCols <- grep(paste("^", vn, "$", sep = ""), colnames(rvTab))
        spCols <- grep(vn, colnames(spt))
        s1var <- sapply(1:nrow(spt), function(i) {
          rvTab[inT[i, 1], randCols]
        })
        s2var <- sapply(1:nrow(spt), function(i) {
          rvTab[inT[i, 2], randCols]
        })
        return(list(s1var, s2var))
      }
    }, rvTab = randVarTab, spt = spTab, inT = indexTab)
    bySite <- lapply(1:2, function(i, vlist) {
      sapply(vlist, function(vl, k) {
        vl[[k]]
      }, k = i)
    }, vlist = varLists)
    if (is(bySite[[1]], "list")) {
      site1Vars <- do.call("cbind", bySite[[1]])
      site2Vars <- do.call("cbind", bySite[[2]])
    } else {
      site1Vars <- bySite[[1]]
      site2Vars <- bySite[[2]]
    }
    newSP <- as.data.frame(cbind(spTab$distance, spTab$weights, 
                                 s1xCoord, s1yCoord, s2xCoord, s2yCoord, site1Vars, site2Vars))
    colnames(newSP) <- colnames(spTab)
  }
  if(includeYear == T){
    #    new code : for includeYear = T
    # idea : indexTab contains the order of coordinates for the pairs. the pairs were
    #    constructed respecting that only few pairs are allowed to be built (only within group)
    #    shuffle around data within group
    #    fill back to spTable based on the indices saved in indexTab.
    newSP <- spTab[, grep("dist|weights", names(spTab), value = T, perl = T)]
    s1.randvartab <- data.table::copy(randVarTab)
    names(s1.randvartab) <- paste("s1.", names(s1.randvartab), sep = "")
    s2.randvartab <- data.table::copy(randVarTab)
    names(s2.randvartab) <- paste("s2.", names(s2.randvartab), sep = "")
    newSP <- cbind(newSP, s1.randvartab[indexTab[, 1]], s2.randvartab[indexTab[, 2]])
    # adjust column order
    setcolorder(newSP, neworder = colnames(spTab))
    # check number of columns
    if(ncol(newSP) != ncol(spTab)){
      # print("adjusting number of columns for permutation in reduced dataset")
      # print(paste("removing column", nameToRemove))
      newSP <- newSP[, grep(nameToRemove, colnames(newSP), invert = T)]
    }
    #####
  }
  class(newSP) <- c(class(spTab))
  return(newSP)
}

# a <- permutateSitePair_edit(currSitePair, siteData, indexTab, varNames, IndPermGroup)
# par(mfrow = c(2, 1))
# hist(c(a$s1.LU, a$s2.LU))
# hist(c(currSitePair$s1.LU, currSitePair$s2.LU))





####################################
# Edit permutateVARsitepair
####################################

# spTab <- currSitePair
# siteVarTab <- siteData
# vNames <- varNames
# indexTab
# IndPermGroup
# vName <- "LU"

permutateVarSitePair_edit <- function (spTab, siteVarTab, indexTab, vName, IndPermGroup, includeYear, spacefortimedataset){
  # print("Using the edited matrix permutation, restricted to 1 variable only.")
  ### EDIT START
  # REPLACED THE LINE BELOW
  # randVarTab <- siteVarTab[sample(nrow(siteVarTab), nrow(siteVarTab)), ]
  # instead of permuting all over the place, permutation is only allowed within groups.
  
  # get permgroup ID and add to siteVarTab
  check_coln <- ncol(spTab)
  names(IndPermGroup) <- sub("s1.", "", names(IndPermGroup))
  if(spacefortimedataset == "temporal"){
    torand <- data.table(merge(siteVarTab, IndPermGroup, by = c("xCoord", "yCoord")))}
  if(spacefortimedataset == "spatial"){
    torand <- data.table(merge(siteVarTab, IndPermGroup, by = c("Year")))}
  
  for(g in unique(torand$permgroup)){
    # permute every group individually
    temp <- torand[permgroup == g, ]
    temp <- temp[sample(nrow(temp), nrow(temp))]
    torand[permgroup == g, ] <- temp
    rm(temp)
  }
  torand[, permgroup := NULL]
  randVarTab <- data.frame(torand)
  rm(torand)
  class(randVarTab) <- class(siteVarTab)
  ### 
  
  # add back the permuted column to spTab
  if(includeYear == F){
    randCols <- grep(paste("^", vName, "$", sep = ""), colnames(randVarTab))
    spCols1 <- grep(paste("^s1.", vName, "$", sep = ""), colnames(spTab))
    spCols2 <- grep(paste("^s2.", vName, "$", sep = ""), colnames(spTab))
    s1var <- sapply(1:nrow(spTab), function(i) {
      randVarTab[indexTab[i, 1], randCols]
    })
    s2var <- sapply(1:nrow(spTab), function(i) {
      randVarTab[indexTab[i, 2], randCols]
    })
    spTab[, spCols1] <- s1var
    spTab[, spCols2] <- s2var
  }
  if(includeYear == T){
    # fill in permuted values of vName (name of columns to be permuted) to spTab
    # note that here we are working with data.frames, not like in the usual case with data.tables
    newSP <- spTab[, grep(vName, colnames(spTab), value = T, invert = T)] # all columns but the permuted ones
    s1.randvartab <- data.frame(randVarTab)
    names(s1.randvartab) <- paste("s1.", names(s1.randvartab), sep = "")
    s1.randvartab <- s1.randvartab[, grep(vName, colnames(s1.randvartab), value = T), drop = F] # chose the columns to be permuted
    s2.randvartab <- data.frame(randVarTab)
    names(s2.randvartab) <- paste("s2.", names(s2.randvartab), sep = "")
    s2.randvartab <- s2.randvartab[, grep(vName, colnames(s2.randvartab), value = T), drop = F] # chose the columns to be permuted
    newSP <- cbind(newSP, s1.randvartab[indexTab[, 1],, drop = F], s2.randvartab[indexTab[, 2],, drop = F])
    # adjust column order
    setcolorder(newSP, neworder = colnames(spTab))
    spTab <- newSP
    # convert to gdm class
    spTab <- formatsitepair(spTab, 4, predData = spTab)
    # check number of columns
    if(check_coln != ncol(spTab)){
      stop("need to write this. Most most probably this never happens because permutateSitePair, which is called before, already fixes it.")
      # below is some code used in permutatesitevar to fix this error.
      # print("adjusting number of columns for permutation in reduced dataset")
      # print(paste("removing column", nameToRemove))
      # newSP <- newSP[, grep(nameToRemove, colnames(newSP), invert = T)]
    }
  }
  return(spTab)
}

# # example use of function
# a <- permutateVarSitePair_edit(currSitePair, siteData, indexTab, vName, IndPermGroup)
# plot(currSitePair$s1.LU, a$s1.LU)
# plot(currSitePair$s2.LU, a$s2.LU)
# # plot(sort(unique(currSitePair$s1.LU)), sort(unique(a$s1.LU))) 
