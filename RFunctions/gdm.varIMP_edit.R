# EDIT GDM.VARIPM FUNCTION
# Caterina Penone, Noelle Schenk 26.03.21
# last edit : 03.02.22 by N.S.


# # by hand definition of function parameters.
# #    can be used to run sub-parts of the function.
# #TODO Noelle Main edit :
# # - be able to run varimp on dataset with less than 3 varaibles
# # - be able to run on data with groups : mixed effects GDM
# #     temporal dataset : plot is random effect
# #     spatial dataset : year is random effect
# # - taken out parallel --> not adjusted for it - better not keep this
# # - Austauschen der 2 permutate functions : restrict to within groups
# #     - unnecessary additional thing : wrote perm function to permute rows (non-matrix permutation)
# #   adjust for data with year : added includeYear == T parameter

#' QUANTIFY MODEL SIGNIFICANCE AND VARIABLE IMPORTANCE IN GDM USING MATRIX PERMUTATION
#' 
#' This function contains a few modifications from the original `gdm.varImp`. It allows (1) datasets with
#' less than 3 predictor variables, (2) allows duplicate coordinates and (3) performs matrix-permutation
#' within groups. More specifically : when running on the "time dataset", plot is treated as random
#' effect and when running on the "space dataset", year is treated as random effect. No parallel
#' computation is allowed any more. The two permutation functions called by `gdm.varImp` were modified
#' as well (`permutateSitePair` and `permutateVarSitePair`).
#' 
#' Only parameters with deviating functionality from the original function are described here.
#' For a full description of the function and its purpose, please read `help("gdm.varImp")`.
#' About the use of this function : Parameters includeYear and spacefortimedataset need to be given.
#' The parameter `fullModelOnly` needs to be set to `T`, because backward elimination is depreciated
#' in this function.
#' The packages `data.table` and `gdm` need to be loaded (this is done within this script), as well 
#' as the modified permutation functions (from script "matrix_perm_permutateSitePair.R").
#' 
#' *Side note* : this function was developed using the script "gdm_demo_script.R" in Summer 2021 and modified later
#' using the script "test_varimp.edit_function_february.R" in February 2022.
#' 
#' If you encounter any bugs/ errors/ unclarities, please report them back to us, e.g. 
#' via https://github.com/biodiversity-exploratories-synthesis .
#' 
#' It has to be noted that the versatility of the gdm.varImp function is strongly restricted by this 
#' modifications. With the modified function, only full models can be assessed, thus no backward 
#' elimination can be performed any more. We strongly recommend to use the original function in all 
#' cases but the specific case presented in the space-for time manuscript (#TODO add reference to manuscritp).
#' Please also consider the additional documentation in "XXXX.md" (#TODO add documentation)
#' 
#' @param includeYear either TRUE or FALSE. If true, requires columns in spTable "s1.Year" and
#' "s2.Year", no other naming is accepted (case sensitive).
#' @param spacefortimedataset either FALSE, "spatial" or "temporal".
#' 
require(data.table)
require(gdm)

####
# # FOR DEVELOPERS
# # please uncomment the following lines to get input variables
# source("GDM/matrix_perm_permutateSitePair.R")
# geo = T
# splines = NULL
# knots = NULL
# fullModelOnly = T # This is not the default setting
# #    changed to T because for the edited function,
# #    no backwards selection is needed.
# #    running code with fullModelOnly = F will raise a warning
# nPerm = 4
# sampleSites = 1
# sampleSitePairs = 1
# outFile = NULL
# includeYear = TRUE # put this to F if space model is runned
# # spacefortimedataset <- "temporal"
# spacefortimedataset <- "spatial"
#     # is either "temporal", "spatial" or FALSE (FALSE needs to be given as well.)
# # read in test data
# # spTable <- readRDS(file = "GDM/testdataset_time.Rds")
# # spTable <- readRDS(file = "GDM/testdataset_space.Rds")
####

gdm.varImp.synthesis.edit <- function(spTable,
                                      geo,
                                      splines = NULL,
                                      knots = NULL,
                                      fullModelOnly = TRUE,
                                      nPerm = 50,
                                      includeYear = FALSE,
                                      spacefortimedataset = FALSE,
                                      sampleSites = 1,
                                      sampleSitePairs = 1,
                                      outFile = NULL)
{
  k <- NULL
  # first part : testing input data
  if (class(spTable)[1] != "gdmData") {
    warning("spTable class does not include type 'gdmData'. Make sure your data is in site-pair format or the gdm model will not fit.")
  }
  if (!(class(spTable)[1] == "gdmData" | class(spTable)[1] == 
        "matrix" | class(spTable)[1] == "data.frame")) {
    stop("spTable argument needs to be gdmData, a matrix, or a data frame")
  }
  if (ncol(spTable) < 6) {
    stop("Not enough columns in data. (Minimum need: Observed, weights, X0, Y0, X1, Y1)")
  }
  if (nrow(spTable) < 1) {
    stop("Not enough rows in data")
  }
  if (!(geo == TRUE | geo == FALSE)) {
    stop("geo argument must be either TRUE or FALSE")
  }
  if (is.null(splines) == FALSE & class(splines) != "numeric") {
    stop("argument splines needs to be a numeric data type")
  }
  if (is.null(knots) == FALSE & class(knots) != "numeric") {
    stop("argument knots needs to be a numeric data type")
  }
  if (!(fullModelOnly == TRUE | fullModelOnly == FALSE)) {
    stop("fullModelOnly argument must be either TRUE or FALSE")
  }
  if ((is.null(nPerm) == FALSE & is.numeric(nPerm) == FALSE) | 
      nPerm < 1) {
    stop("argument nPerm needs to be a positive integer")
  }
  # if (parallel == TRUE) {
  #   stop("parallel argument was taken out, please use original function for running in parallel.")
  # }
  if (exists("cores") == TRUE) {
    stop("parallel running was depreciated in this function. Please do not specify cores argument. Use the original function for running in parallel.")
  }
  if (!spacefortimedataset %in% c("temporal", "spatial", F)){
    stop("argument spacefortimedataset needs to be specified, either 'temporal', 'spatial' or 'F'.")
  }
  if(!exists("includeYear")){
    stop("argument includeYear needs to be either 'TRUE' or 'FALSE'.")
  }
  if(!includeYear %in% c(T, F)){
    stop("argument includeYear needs to be either 'TRUE' or 'FALSE'.")
  }
  if (is.numeric(sampleSites) == FALSE | sampleSites < 0 | 
      sampleSites > 1) {
    stop("argument sampleSites needs to be a positive number between 0 and 1")
  }
  if (is.numeric(sampleSitePairs) == FALSE | sampleSitePairs < 
      0 | sampleSitePairs > 1) {
    stop("argument sampleSitePairs needs to be a positive number between 0 and 1")
  }
  if (sampleSites == 0) {
    stop("a sampleSites value of 0 will remove all sites from the analysis")
  }
  if (sampleSitePairs == 0) {
    stop("a sampleSitePairs value of 0 will remove all sites from the analysis")
  }
  if (is.null(outFile) == FALSE) {
    if (is.character(outFile) == FALSE) {
      stop("argument outFile needs to be a character string of the directory and file name you wish the tables to be written to")
    }
    outFileChar <- nchar(outFile)
    if (substr(outFile, outFileChar - 5, outFileChar) != 
        ".RData") {
      outFile <- paste(outFile, ".RData", sep = "")
    }
    if (length(strsplit(outFile, "/")[[1]]) > 1) {
      splitOutFile <- strsplit(outFile, "/")[[1]][-length(strsplit(outFile, 
                                                                   "/")[[1]])]
      dir.create(paste(splitOutFile, collapse = "/"))
    }
    else {
      outFile <- paste("./", outFile, sep = "")
    }
  }
  nPerm <- as.integer(nPerm)
  if (sampleSites < 1) {
    spTable <- removeSitesFromSitePair(spTable, sampleSites = sampleSites)
    if (sampleSitePairs < 1) {
      warning("You have selected to randomly remove sites and/or site-pairs.")
    }
  }
  if (sampleSitePairs < 1) {
    numRm <- sample(1:nrow(spTable), round(nrow(spTable) * 
                                             (1 - sampleSitePairs)))
    spTable <- spTable[-c(numRm), ]
  }
  rtmp <- spTable[, 1] # get distance as vector
  if (length(rtmp[rtmp < 0]) > 0) { # is there any negative distance?
    stop("Response spTable has negative values. Must be between 0 - 1.")
  }
  if (length(rtmp[rtmp > 1]) > 0) { # is there any distance > 1?
    stop("Response spTable has values greater than 1. Must be between 0 - 1.")
  }
  #### here, testing is over
  
  nVars <- (ncol(spTable) - 6)/2 # take away distance, weights, s1.xCoord, ..., s2.yCoord (6 columns)
  if(nVars == 0){
    stop("GDM with only 1 variable given. No p Value calculation can be done on that.")
  }
  varNames <- colnames(spTable[c(7:(6 + nVars))])
  varNames <- sapply(strsplit(varNames, "s1."), "[[", 2)
  if (geo == TRUE) {
    nVars <- nVars + 1
    varNames <- c("Geographic", varNames)
  }
  
  #EDIT : the following was put into an if statement to have the year-option separately
  # Indices are constructed based on unique coordinates.
  if(includeYear == F){
    sortMatX <- sapply(1:nrow(spTable), function(i, spTab) {
      c(spTab[i, 3], spTab[i, 5]) # get x Coordinates (for s1 and s2)
    }, spTab = spTable)
    sortMatY <- sapply(1:nrow(spTable), function(i, spTab) {
      c(spTab[i, 4], spTab[i, 6])# get y Coordinates (for s1 and s2)
    }, spTab = spTable)
    sortMatNum <- sapply(1:nrow(spTable), function(i) {
      c(1, 2) # create matrix with 2 rows and nrow columns. First row = 1, second one = 2
    })
    sortMatRow <- sapply(1:nrow(spTable), function(i) {
      c(i, i)
    })
    fullSortMat <- cbind(as.vector(sortMatX), as.vector(sortMatY),
                         as.vector(sortMatNum), as.vector(sortMatRow), rep(NA,
                                                                           length(sortMatX)))
    # matrix with nrow(spTable)*2 rows and 5 columns
    #   X, Y, Num, Row and NA
    siteByCoords <- as.data.frame(unique(fullSortMat[, 1:2]))
    # here, the identifiers come in! (X and Y coordinates)
    numSites <- nrow(siteByCoords)
    for (i in 1:numSites) { # This for loop creates the site identifiers. Could be changed here
      fullSortMat[which(fullSortMat[, 1] == siteByCoords[i,
                                                         1] & fullSortMat[, 2] == siteByCoords[i, 2]), 5] <- i
      # which row corresponds to x coordinate of plot1 AND the Y coordinate of plot 1? --> set to plot id "i"
    }
    indexTab <- matrix(NA, nrow(spTable), 2)
    for (iRow in 1:nrow(fullSortMat)) {
      indexTab[fullSortMat[iRow, 4], fullSortMat[iRow, 3]] <- fullSortMat[iRow,
                                                                          5]
    } # Table of plot comparison indices (e.g. plot 1 versus plot 2 is first row, ...)
    rm(fullSortMat)
    rm(sortMatX)
    rm(sortMatY)
    rm(sortMatNum)
    rm(sortMatRow)
    rm(siteByCoords)
  }
  
  ### SYNTHESIS EDIT FOR YEARS
  # additional if clause for an option to include measurements over various years.
  # Indices (unique identifiers) are constructed based on unique coordinates AND years, not unique 
  # coordinates as for includeYear = F and the original function.
  # This means, each plot in each year gets a unique identifier, which can be shuffled.
  # Main creation of this piece of code is indexTab, a data frame with indices in the same order
  #    as they appear in spTable.
  if(includeYear == T){
    print("You are working with time series data, the edited function is running now.")
    spHelper <- data.table::copy(spTable)
    spHelper <- data.table(spHelper)
    # spHelper[, .(s1.xCoord, s1.yCoord, s2.xCoord, s2.yCoord, s1.Year, s2.Year)]
    # get all included unique plot-year combinations
    #    need to stack on top of each other s1. and s2. coordinates + year
    IndPlotYear <- unique(rbind(spHelper[, .(s1.xCoord, s1.yCoord, s1.Year)],
                                spHelper[, .(s2.xCoord, s2.yCoord, s2.Year)], use.names = F))
    IndPlotYear[, index := seq(1, nrow(IndPlotYear), 1)] # create unique identifier (index) for each plot in each year
    
    if(spacefortimedataset == "temporal"){
      print("You are working in the temporal dataset. Plot is treated as random variable.")
      # create the permutation groups here. Permutation will only happen WITHIN those groups,
      #    in order to respect the hierarchical design.
      # if spacefortimedataset is spatial, the group is PLOT
      IndPermGroup <- unique(IndPlotYear[, .(s1.xCoord, s1.yCoord)])
      # create unique group identifiers : for temporal dataset : each plot gets a unique
      #    identifiers. Rows WITHIN the same group (within the same plot) will be
      #    permuted.
      IndPermGroup[, permgroup := seq(1, nrow(IndPermGroup))]
      IndPlotYear <- merge(IndPlotYear, IndPermGroup, by = c("s1.xCoord", "s1.yCoord"))
    }
    if(spacefortimedataset == "spatial"){
      print("You are working in the spatial dataset. Year is treated as random variable.")
      # create the permutation groups here. Permutation will only happen WITHIN those groups,
      #    in order to respect the hierarchical design.
      # if spacefortimedataset is spatial, the group is YEAR
      IndPermGroup <- unique(IndPlotYear[, .(s1.Year)])
      IndPermGroup[, permgroup := seq(1, nrow(IndPermGroup))]
      IndPlotYear <- merge(IndPlotYear, IndPermGroup, by = c("s1.Year"))
    }
    
    # read out the plot comparisons from spHelper and add to the dataset
    #    first add identifiers to s1.
    spHelper <- merge(spHelper, IndPlotYear, by = c("s1.xCoord", "s1.yCoord", "s1.Year"), all.x = T) # which rows in s1. correspond to which index?
    setnames(spHelper, old = "index", new = "index1")
    setnames(spHelper, old = "permgroup", new = "permgroup1")
    #    change name of columns to s2. and merge identifiers to s2.
    colnames(IndPlotYear) <- sub("s1.", "s2.", colnames(IndPlotYear))
    spHelper <- merge(spHelper, IndPlotYear, by = c("s2.xCoord", "s2.yCoord", "s2.Year"), all.x = T) # which rows in s2. correspond to which index?
    setnames(spHelper, old = "index", new = "index2")
    setnames(spHelper, old = "permgroup", new = "permgroup2")
    # graphical test : permgroup of all pairs should be the same (always permgroup1 == permgroup2)
    #    because in the dataset, only pairs within the same plot do exist.
    # plot(spHelper[, .(permgroup1, permgroup2)])
    if(!isFALSE(spacefortimedataset) && !all(spHelper$permgroup1 == spHelper$permgroup2)){
      stop("The dataset contains pairs which are not of the same permutation group. \n
           please consider your choice of the parameter spacefortimedataset.")
    }
    
    # create matrix indexTab containing indices where each row contains the index of the corresponding
    #    row of spTable.
    indexTab <- as.matrix(spHelper[, .(index1, index2)])
    colnames(indexTab) <- NULL
    if (nrow(indexTab) != nrow(spHelper)) {
      stop("Error in index creation. indexTab and spHelper do not have the same length, check code.")
    }
    spHelper[, index1 := NULL]
    spHelper[, index2 := NULL]
    
    # create data.table indexPermgroup containing group indices where each row contains the group ID of the
    #    corresponding row of spTable
    #    delete permgroup1 and permgroup2 from spHelper
    indexPermgroup <- spHelper[, .(permgroup1)]
    spHelper[, permgroup1 := NULL]
    spHelper[, permgroup2 := NULL]
    
    numSites <- length(unique(c(indexTab[,1], indexTab[,2])))
    spHelper <- as.data.frame(spHelper)
  }
  ### Synthesis EDIT END
  
  ###
  # goal of next code section : 
  #    create siteData : a data.frame (gdmData) with rows containing : 
  #    coordinates of each site and the variables for each site. (not in pairwise matrix format any more)
  # Construction for includeyear = F is different from includeyear = T
  #    includeyear = F : original function is used : looks up coordinates one by one in original dataset 
  #        and constructs list
  #    includeyear = T : stacks datasets on top of each other and uses unique rows.
  # Problem in both: if 2 different roundings/ variables are taken, duplicated rows with differing
  #    response variable are done. is ok for permutation. 
  if(includeYear == F){
    # prepare a list, where each element corresponds to all variables of one single plot
    #    and store it to siteData
    exBySite <- lapply(1:numSites, function(i, index, tab) { # start to put input variables of each site to separate element in list
      rowSites <- which(index[, 1] %in% i)
      if (length(rowSites) < 1) {
        rowSites <- which(index[, 2] %in% i)
      }
      exSiteData <- tab[rowSites[1], ]
      return(exSiteData)
    }, index = indexTab, tab = spTable)
    
    outSite <- which(!(1:numSites %in% indexTab[, 1])) # find plots which ARE only occurring in second column of indexTab
    # more than 1 if several years are given
    for (i in 1:length(exBySite)) {
      siteRow <- exBySite[[i]]
      if (i %in% outSite) { # for the last plot, the coordinates are in s2. for all others, in s1
        siteRow <- siteRow[grep("s2.", colnames(siteRow))]
        colnames(siteRow) <- sapply(strsplit(colnames(siteRow),
                                             "s2."), "[[", 2)
      }
      else {
        siteRow <- siteRow[grep("s1.", colnames(siteRow))]
        colnames(siteRow) <- sapply(strsplit(colnames(siteRow),
                                             "s1."), "[[", 2)
      }
      exBySite[[i]] <- siteRow
    }
    siteData <- unique(do.call("rbind", exBySite)) # Edit : took unique here, to avoid construction of duplicated rows
  }
  if(includeYear == T){
    ### edit start
    # new approach, without list.
    # create from spTable a new table "siteData" which contains all data in non-matrix format.
    # only the explanatory variables
    temp_SD <- data.table(data.table::copy(spTable))
    s1cols <- grep("s1.", names(temp_SD), value = T) # colnames of s1 to include in siteData
    s1temp <- unique(temp_SD[, ..s1cols])
    names(s1temp) <- sub("s1.", "", names(s1temp))
    s2cols <- grep("s2.", names(temp_SD), value = T) # colnames of s2 to include in siteData
    s2temp <- unique(temp_SD[, ..s2cols])
    names(s2temp) <- sub("s2.", "", names(s2temp))
    siteData <- unique(rbindlist(list(s1temp, s2temp), use.names = T))
    rm(temp_SD); rm(s1cols); rm(s2cols)
    ### edit end
  }
  
  
  ###
  # prepare return matrices
  #    here, the dimname error occurred. Solved with if clause.
  #    Source of error : Code stores the results of leave-one-out gdm (i.e. if one variable is 
  #       taken out) (backwards selection). 
  #       BUT we only have two variables, if <3 Variables are given! --> error when trying to 
  #       perform a GDM with <2 vars.
  #           In this case, just return NA.
  if(nVars > 1){
    modelTestValues <- matrix(NA, 4, nVars, dimnames = list(c("Model deviance", 
                                                              "Percent deviance explained", "Model p-value", "Fitted permutations"), 
                                                            c("fullModel", paste("fullModel-", seq(1, nVars - 1), 
                                                                                 sep = ""))))
  }
  if(nVars == 1){
    # modelTestValues <- matrix(NA, 4, nVars, dimnames = list(c("Model deviance", 
    #                                                           "Percent deviance explained", "Model p-value", "Fitted permutations"), 
    #                                                         c("fullModel")))
    stop("Special case: You have given less than 3 Variables, but the matrix devReductVars does not have 1 column. The function
           needs to be edited. Please report back :)")
  }
  
  devReductVars <- matrix(NA, nVars, nVars - 1) #takes GEO and all other variables 
  rownames(devReductVars) <- varNames
  # EDIT : the else case of the following statement was the original code.
  #    An if case and a test have been added.
  if(nVars < 3){
    print("You have given less than 3 Variables, the edited function is running now.")
    if(ncol(devReductVars) != 1){
      stop("Special case: You have given less than 3 Variables, but the matrix devReductVars does not have 1 column. The function
           needs to be edited. Please report back :)")
    }
    colnames(devReductVars) <- "fullModel"
  } else {
    colnames(devReductVars) <- c("fullModel", paste("fullModel-", 
                                                    seq(1, nVars - 2), sep = "")) # tries to assign 2 names minimum to a matrix with one column
  }
  
  ###
  # next code section :
  # ev. calc models with -1 variable?
  # Permutation takes part here
  # definitely here : calc p-value of full model
  pValues <- numPermsFit <- devReductVars
  currSitePair <- spTable
  nullGDMFullFit <- 0
  ## 1.) create the full model gdm without any permutation
  ##### HERE, large loop starts
  # for (v in 1:(nVars-1)) { # dirty edit : only run v through 1 and 2 to avoid fullmod-2 problems
  for(v in 1:min(nVars-1, 2)){
    # EDIT : IMPORTANT : changed from 1:nVars to 1:(nVars-1) in order to not run the last round.
    #    The function was adapted to run on datasets with less than 3 variables, which causes problems here.
    #    For space-for-time dataset, this means that the output from varimp usually only contains the full model
    #    and the full model -1.
    print(paste("round", v, "of", max(min(nVars-1, 2)), sep = " "))
    fullGDM <- gdm(currSitePair, geo = geo, splines = splines, 
                   knots = knots)
    if (is.null(fullGDM) == TRUE) {
      warning(paste("The model did not fit when testing variable number ", 
                    v, ". Terminating analysis. Returning output objects filled up to the point of termination.", 
                    sep = ""))
      break
    }
    ## 2.) create permuted datasets (are stored in a list permSitePairs of length nPerm)
    if (0 == TRUE) {
      stop("parallel option was taken out from this function. Please use the original function for running in parallel.")
    } else {
      # this loop : fills modelTestValues, one of the return objects.
      # From package descript : The environmental data are permuted nPerm times (by randomizing the 
      #    order of the rows) and a GDM is fit to each permuted table. Model significance is determined 
      #    by comparing the deviance explained by GDM fit to the un-permuted table to the distribution 
      #    of deviance explained values from GDM fit to the nPerm permuted tables.
      # creation of permSitePairs, the list of permuted datasets
      # Edit : option 1 : use original function for permutating across all groups (in case there are no groups)
      if(isFALSE(spacefortimedataset)){
        permSitePairs <- lapply(1:nPerm, function(i, csp, 
                                                  # This line creates a list of all variables, with length nPerms 
                                                  #    which contains permuted rows
                                                  sd, it, vn) {
          permutateSitePair(csp, sd, it, vn) # function which randomizes the rows of a given site table 
          #                                     (before pairs are built)
        }, csp = currSitePair, sd = siteData, it = indexTab, 
        vn = varNames) # permutateSitePair(currSitePair, siteData, indexTab, varNames)
      } else if(spacefortimedataset %in% c("temporal", "spatial")){
        # EDIT Start with else if : use the newly created permutation function to create permuted datsets
        #    in case we are working with a spacefortime dataset.
        # print("Restricted permutation to within-groups.")
        permSitePairs <- lapply(1:nPerm, function(i, csp, sd, it, vn, ipg, iy, sftd, ntr) {
          permutateSitePair_edit(csp, sd, it, vn, ipg, iy, sftd, ntr) # function which randomizes the rows of a given site table 
          #                                     (before pairs are built)
        }, csp = currSitePair, sd = siteData, it = indexTab, vn = varNames, ipg = IndPermGroup, 
        iy = includeYear, sftd = spacefortimedataset, ntr = nameToRemove)
      }
      
      ## 3.) calculate gdm on each of the permuted datasets from the list
      #      i.e. run GDM on all permuted variants in the list permSitePairs
      permGDM <- lapply(permSitePairs, gdm, geo = geo, 
                        splines = NULL, knots = NULL)
    }
    permModelDev <- sapply(permGDM, function(mod) {
      mod$gdmdeviance
    })
    modPerms <- length(which(sapply(permModelDev, is.null) == 
                               TRUE))
    if (modPerms > 0) {
      permModelDev <- unlist(permModelDev[-(which(sapply(permModelDev, 
                                                         is.null) == T))])
    }
    modelTestValues[1, v] <- fullGDM$gdmdeviance
    modelTestValues[2, v] <- fullGDM$explained
    modelTestValues[3, v] <- sum(permModelDev <= fullGDM$gdmdeviance)/(nPerm - 
                                                                         modPerms)
    modelTestValues[4, v] <- nPerm - modPerms
    if (length(varNames) < 2) {
      print("less than 2 variables given, this raises an error and breaks the function.")
      print(geo)
      print(varNames)
      break
    }
    if (geo == TRUE) {
      # compare model with  model without GEO.
      noGeoGDM <- gdm(currSitePair, geo = FALSE, splines = NULL, 
                      knots = NULL) # run model without GEO
      if (0 == TRUE) {
        stop("parallel option was taken out from this function. Please use the original function for running in parallel.")
      }
      else {
        # print(c("Compare to model without geo. iteration v = ", v))
        #  create permuted datasets
        if(isFALSE(spacefortimedataset)){
          # same as above, prepare list of length nPerm with randomised rows
          permSitePairs <- lapply(1:nPerm, function(i, csp, sd, it, vn) {
            permutateSitePair(csp, sd, it, vn)
          }, csp = currSitePair, sd = siteData, it = indexTab, 
          vn = varNames)
        } else if(spacefortimedataset %in% c("temporal", "spatial")){
          # Edit is identical to above creation of permuted datasets
          # EDIT this else if loop : use the newly created permutation function to create permuted datsets
          #    in case we are working with a spacefortime dataset.
          # print("Restricted permutation to within-groups.")
          permSitePairs <- lapply(1:nPerm, function(i, csp, sd, it, vn, ipg, iy, sftd, ntr) {
            permutateSitePair_edit(csp, sd, it, vn, ipg, iy, sftd, ntr) # function which randomizes the rows of a given site table 
            #                                     (before pairs are built)
          }, csp = currSitePair, sd = siteData, it = indexTab, vn = varNames, ipg = IndPermGroup, 
          iy = includeYear, sftd = spacefortimedataset, ntr = nameToRemove)
        }
        
        # NOTE : depending on the permutation, the below model raises a warning. In this case, the dev is NULL.
        permGDM <- lapply(permSitePairs, gdm, geo = geo,  # run GDM on the prepared list of permuted variables
                          splines = NULL, knots = NULL)
      }
      permModelDev <- sapply(permGDM, function(mod) {
        mod$gdmdeviance
      })
      modPerms <- length(which(sapply(permModelDev, is.null) == 
                                 TRUE))
      if (modPerms > 0) {
        permModelDev <- unlist(permModelDev[-(which(sapply(permModelDev, 
                                                           is.null) == T))])
      }
      if (is.null(noGeoGDM$gdmdeviance) == TRUE) {
        permDevReduct <- -9999
        devReductVars[1, v] <- -9999
        pValues[1, v] <- -9999
      } else {
        permDevReduct <- noGeoGDM$gdmdeviance - permModelDev
        devReductVars[1, v] <- 100 * abs((noGeoGDM$explained - 
                                            fullGDM$explained)/fullGDM$explained)
        pValues[1, v] <- sum(permDevReduct >= (noGeoGDM$gdmdeviance - 
                                                 fullGDM$gdmdeviance))/(nPerm - modPerms)
      }
      numPermsFit[1, v] <- nPerm - modPerms
    }
    for (varChar in varNames) {
      # print("permuting single variables.")
      if (varChar != "Geographic") {
        # take away the given  variable and run model without it.
        testVarCols1 <- grep(paste("^s1.", varChar, "$",
                                   sep = ""), colnames(currSitePair))
        testVarCols2 <- grep(paste("^s2.", varChar, "$",
                                   sep = ""), colnames(currSitePair))
        testSitePair <- currSitePair[, -c(testVarCols1,
                                          testVarCols2)]
        noVarGDM <- gdm(testSitePair, geo = geo, splines = NULL,
                        knots = NULL)
        if (0 == TRUE) {
          stop("parallel option was taken out from this function. Please use the original function for running in parallel.")
        }
        else {
          # prepare list of length nPerm with randomised rows again
          
          # EDIT START :
          # Permutation is either done as in original function (upper if case)
          # or Permutatino is restricted to within groups (else if case)
          if(isFALSE(spacefortimedataset)){
            noVarSitePairs <- lapply(1:nPerm, function(i,
                                                       csp, sd, it, vn) {
              permutateVarSitePair(csp, sd, it, vn)
            }, csp = currSitePair, sd = siteData, it = indexTab,
            vn = varChar)
          } else if(spacefortimedataset %in% c("temporal", "spatial")){
            noVarSitePairs <- lapply(1:nPerm, function(i, csp, sd, it, vn, ipg, iy, sftd) {
              permutateVarSitePair_edit(csp, sd, it, vn, ipg, iy, sftd)
            }, csp = currSitePair, sd = siteData, it = indexTab, vn = varChar, ipg = IndPermGroup, iy = includeYear, sftd = spacefortimedataset)
          }
          # Edit end
          
          permGDM <- lapply(noVarSitePairs, gdm, geo = geo,
                            splines = NULL, knots = NULL)
        }
        permModelDev <- sapply(permGDM, function(mod) {
          mod$gdmdeviance
        })
        modPerms <- length(which(sapply(permModelDev,
                                        is.null) == TRUE))
        if (modPerms > 0) {
          permModelDev <- unlist(permModelDev[-(which(sapply(permModelDev,
                                                             is.null) == T))])
        }
        if (is.null(noVarGDM$gdmdeviance) == TRUE) {
          permDevReduct <- -9999
          ggg <- which(rownames(devReductVars) %in% varChar)
          devReductVars[ggg, v] <- rep(-9999, times = length(ggg))
          pValues[ggg, v] <- rep(-9999, times = length(ggg))
        } else {
          # calc dev of full model and the model without the variable chosen
          #    in order to fill devReductVars
          permDevReduct <- noVarGDM$gdmdeviance - permModelDev
          devReductVars[which(rownames(devReductVars) %in%
                                varChar), v] <- 100 * abs((noVarGDM$explained -
                                                             fullGDM$explained)/fullGDM$explained)
          pValues[which(rownames(pValues) %in% varChar),
                  v] <- sum(permDevReduct >= (noVarGDM$gdmdeviance -
                                                fullGDM$gdmdeviance))/(nPerm - modPerms)
        }
        numPermsFit[which(rownames(numPermsFit) %in%
                            varChar), v] <- nPerm - modPerms
      }
    }
    if (fullModelOnly != TRUE) { # we don't need backward elimination
      print(fullModelOnly)
      warning("Please set the option fullModelOnly to TRUE. \nExplanation : At the moment, the function is runned with 
            fullModelOnly = F, the edited version of this function is not written
            for this parameter setting and will most probably cause problems.
            Please either use the original function or edit the function section for
            fullModelOnly = F.") # Note this warning has been added
      break
    }
    tempPVals <- as.numeric(pValues[c(1:nVars), v])
    tempDevs <- as.numeric(devReductVars[c(1:nVars), v])
    tempPVals <- tempPVals[!is.na(tempPVals)]
    tempDevs <- tempDevs[!is.na(tempDevs)]
    varToOmit <- which.max(tempPVals)
    for (iCheck in 1:length(varNames)) {
      if (tempPVals[iCheck] == tempPVals[varToOmit]) {
        if (tempDevs[iCheck] < tempDevs[varToOmit]) {
          varToOmit <- iCheck
        }
      }
    }
    if (varToOmit == 1 & geo == TRUE) {
      # check if a varaible can still be removed from currSitePair
      # removing geo as a last step
      geo <- FALSE
      varNames <- varNames[-1]
    } else {
      # take away a variable from currSitePair for next round.
      nameToRemove <- varNames[varToOmit]
      varNames <- varNames[-varToOmit]
      removeFromSitePs1 <- grep(paste("^s1.", nameToRemove, 
                                      "$", sep = ""), colnames(currSitePair))
      removeFromSitePs2 <- grep(paste("^s2.", nameToRemove, 
                                      "$", sep = ""), colnames(currSitePair))
      currSitePair <- currSitePair[, -c(removeFromSitePs1, 
                                        removeFromSitePs2)]
    }
  }
  ##### HERE, large loop ends
  outObject <- list(modelTestValues, devReductVars, pValues, 
                    numPermsFit)
  if (is.null(outFile) == FALSE) {
    save(outObject, file = outFile)
  }
  return(outObject)
}
NULL