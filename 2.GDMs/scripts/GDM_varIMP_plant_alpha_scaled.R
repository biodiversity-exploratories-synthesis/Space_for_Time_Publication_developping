######################################################
#       GDM p-values for plants alpha diversity (space and time)
#
######################################################
# set working directory to folder "Space_for_Time_Publication"
# upload data
# plants
load("./2.GDMs/data/InputData/pwise_time_plants.RData")
load("./2.GDMs/data/InputData/pwise_space_plants.RData")

pwise_time_ayrs<- pwise_time_plants
pwise_space<- pwise_space_plants

# if the uploaded, assembled data files are used, upload the complete insect files
# here.

# packages
library(tidyverse)
library(gdm)
library(data.table)
library(vegan)

# edited var.Imp function
source("./RFunctions/gdm.varIMP_edit.R")
source("./RFunctions/matrix_perm_permutateSitePair.R")

pwise_space$da0st <- decostand(pwise_space$da0, method = "range", na.rm = T)
pwise_space$da1st <- decostand(pwise_space$da1, method = "range", na.rm = T)
pwise_space$da2st <- decostand(pwise_space$da2, method = "range", na.rm = T)
pwise_space$da3st <- decostand(pwise_space$da3, method = "range", na.rm = T)
pwise_space$da4st <- decostand(pwise_space$da4, method = "range", na.rm = T)

pwise_time_ayrs$da0st <- decostand(pwise_time_ayrs$da0abs, method = "range", na.rm = T)
pwise_time_ayrs$da1st <- decostand(pwise_time_ayrs$da1abs, method = "range", na.rm = T)
pwise_time_ayrs$da2st <- decostand(pwise_time_ayrs$da2abs, method = "range", na.rm = T)
pwise_time_ayrs$da3st <- decostand(pwise_time_ayrs$da3abs, method = "range", na.rm = T)
pwise_time_ayrs$da4st <- decostand(pwise_time_ayrs$da4abs, method = "range", na.rm = T)


# GDMs space
LU1 <- c("LUI1res", "MOW1res", "GRA1res", "FER1res")
LU2 <- c("LUI2res", "MOW2res", "GRA2res", "FER2res")

LU <- c("LUI", "MOW", "GRA", "FER")

# a0
# list for results
var_plant_b0 <- list()

for (i in 1:4) {
  attach(pwise_space)

  # geo dist, LUI and dummy as predictors

  dat <- data.frame(
    da0st,
    rep(1, nrow(pwise_space)),
    HWG1, RWG1, HWG2, RWG2,
    YR,
    pwise_space[, colnames(pwise_space) == LU1[i]],
    YR,
    pwise_space[, colnames(pwise_space) == LU2[i]]
  )

  colnames(dat) <- c(
    "distance", "weights",
    "s1.xCoord", "s1.yCoord", "s2.xCoord", "s2.yCoord",
    "s1.Year", "s1.LU", "s2.Year", "s2.LU"
  )

  # # check dat
  # nplots <- nrow(unique(rbindlist(list(unique(dat[, c("s1.xCoord", "s1.yCoord")]),
  #                                      unique(dat[, c("s2.xCoord", "s2.yCoord")])),
  #                                 use.names = F))) # 150 plots
  # nyears <- length(unique(c(dat$s1.Year, dat$s2.Year)))
  # length(unique(c(dat$s1.LU, dat$s2.LU))) <= nyears * nplots
  # nrow(dat) == (choose((nplots + 2 - 1), 2) - nplots) * nyears # number of rows as expected
  # nrow(unique(rbindlist(list(unique(dat[which(dat$s1.Year == "2008"), c("s1.xCoord", "s1.yCoord", "s1.LU")]),
  #                            unique(dat[which(dat$s2.Year == "2008"), c("s2.xCoord", "s2.yCoord", "s2.LU")])),
  #                       use.names = F))) # 150 plots
  # test <- unique(rbindlist(list(unique(dat[which(dat$s1.Year == "2008"), c("s1.xCoord", "s1.yCoord", "s1.LU", "s1.Year")]),
  #                       unique(dat[which(dat$s2.Year == "2008"), c("s2.xCoord", "s2.yCoord", "s2.LU", "s2.Year")])),
  #                  use.names = F))
  # test <- test[with(test, order(s1.xCoord, s1.yCoord)),]
  # test # expected number of LUI values and rows present
  # #### testing over

  dat <- dat[complete.cases(dat), ]

  detach(pwise_space)

  gdmTab <- formatsitepair(dat, 4, predData = dat)
  var <- gdm.varImp.synthesis.edit(gdmTab,
    geo = T, includeYear = T,
    fullModelOnly = T, nPerm = 50, spacefortimedataset = "spatial"
  )

  var_plant_b0[[i]] <- var
}

names(var_plant_b0) <- LU
var_plant_a0s_sc <- var_plant_b0
save(var_plant_a0s_sc, file = "./2.GDMs/data/OutputData/var_plant_a0s_sc.RData")


# a1
# list for results
var_plant_b1 <- list()

for (i in 1:4) {
  attach(pwise_space)

  # geo dist, LUI and dummy as predictors


  dat <- data.frame(
    da1st,
    rep(1, nrow(pwise_space)),
    HWG1, RWG1, HWG2, RWG2,
    YR,
    pwise_space[, colnames(pwise_space) == LU1[i]],
    YR,
    pwise_space[, colnames(pwise_space) == LU2[i]]
  )

  colnames(dat) <- c(
    "distance", "weights",
    "s1.xCoord", "s1.yCoord", "s2.xCoord", "s2.yCoord",
    "s1.Year", "s1.LU", "s2.Year", "s2.LU"
  )

  # # check dat
  # nplots <- nrow(unique(rbindlist(list(unique(dat[, c("s1.xCoord", "s1.yCoord")]),
  #                                      unique(dat[, c("s2.xCoord", "s2.yCoord")])),
  #                                 use.names = F))) # 150 plots
  # nyears <- length(unique(c(dat$s1.Year, dat$s2.Year)))
  # length(unique(c(dat$s1.LU, dat$s2.LU))) <= nyears * nplots
  # nrow(dat) == (choose((nplots + 2 - 1), 2) - nplots) * nyears # number of rows as expected
  # nrow(unique(rbindlist(list(unique(dat[which(dat$s1.Year == "2008"), c("s1.xCoord", "s1.yCoord", "s1.LU")]),
  #                            unique(dat[which(dat$s2.Year == "2008"), c("s2.xCoord", "s2.yCoord", "s2.LU")])),
  #                       use.names = F))) # 150 plots
  # test <- unique(rbindlist(list(unique(dat[which(dat$s1.Year == "2008"), c("s1.xCoord", "s1.yCoord", "s1.LU", "s1.Year")]),
  #                               unique(dat[which(dat$s2.Year == "2008"), c("s2.xCoord", "s2.yCoord", "s2.LU", "s2.Year")])),
  #                          use.names = F))
  # test <- test[with(test, order(s1.xCoord, s1.yCoord)),]
  # test # expected number of LUI values and rows present
  # #### testing over

  dat <- dat[complete.cases(dat), ]

  detach(pwise_space)

  gdmTab <- formatsitepair(dat, 4, predData = dat)
  var <- gdm.varImp.synthesis.edit(gdmTab,
    geo = T, includeYear = T,
    fullModelOnly = T, nPerm = 50, spacefortimedataset = "spatial"
  )

  var_plant_b1[[i]] <- var
}

names(var_plant_b1) <- LU
var_plant_a1s_sc <- var_plant_b1
save(var_plant_a1s_sc, file = "./2.GDMs/data/OutputData/var_plant_a1s_sc.RData")



# a2
# list for results
var_plant_b2 <- list()

for (i in 1:4) {
  attach(pwise_space)

  # geo dist, LUI and dummy as predictors

  dat <- data.frame(
    da2st,
    rep(1, nrow(pwise_space)),
    HWG1, RWG1, HWG2, RWG2,
    YR,
    pwise_space[, colnames(pwise_space) == LU1[i]],
    YR,
    pwise_space[, colnames(pwise_space) == LU2[i]]
  )

  colnames(dat) <- c(
    "distance", "weights",
    "s1.xCoord", "s1.yCoord", "s2.xCoord", "s2.yCoord",
    "s1.Year", "s1.LU", "s2.Year", "s2.LU"
  )

  # # check dat
  # nplots <- nrow(unique(rbindlist(list(unique(dat[, c("s1.xCoord", "s1.yCoord")]),
  #                                      unique(dat[, c("s2.xCoord", "s2.yCoord")])),
  #                                 use.names = F))) # 150 plots
  # nyears <- length(unique(c(dat$s1.Year, dat$s2.Year)))
  # length(unique(c(dat$s1.LU, dat$s2.LU))) <= nyears * nplots
  # nrow(dat) == (choose((nplots + 2 - 1), 2) - nplots) * nyears # number of rows as expected
  # nrow(unique(rbindlist(list(unique(dat[which(dat$s1.Year == "2008"), c("s1.xCoord", "s1.yCoord", "s1.LU")]),
  #                            unique(dat[which(dat$s2.Year == "2008"), c("s2.xCoord", "s2.yCoord", "s2.LU")])),
  #                       use.names = F))) # 150 plots
  # test <- unique(rbindlist(list(unique(dat[which(dat$s1.Year == "2008"), c("s1.xCoord", "s1.yCoord", "s1.LU", "s1.Year")]),
  #                               unique(dat[which(dat$s2.Year == "2008"), c("s2.xCoord", "s2.yCoord", "s2.LU", "s2.Year")])),
  #                          use.names = F))
  # test <- test[with(test, order(s1.xCoord, s1.yCoord)),]
  # test # expected number of LUI values and rows present
  # #### testing over

  dat <- dat[complete.cases(dat), ]

  detach(pwise_space)

  gdmTab <- formatsitepair(dat, 4, predData = dat)
  var <- gdm.varImp.synthesis.edit(gdmTab,
    geo = T, includeYear = T,
    fullModelOnly = T, nPerm = 50, spacefortimedataset = "spatial"
  )

  var_plant_b2[[i]] <- var
}

names(var_plant_b2) <- LU
var_plant_a2s_sc <- var_plant_b1
save(var_plant_a2s_sc, file = "./2.GDMs/data/OutputData/var_plant_a2s_sc.RData")



# b3
# list for results
var_plant_b3 <- list()

for (i in 1:4) {
  attach(pwise_space)

  # geo dist, LUI and dummy as predictors

  dat <- data.frame(
    da3st,
    rep(1, nrow(pwise_space)),
    HWG1, RWG1, HWG2, RWG2,
    YR,
    pwise_space[, colnames(pwise_space) == LU1[i]],
    YR,
    pwise_space[, colnames(pwise_space) == LU2[i]]
  )

  colnames(dat) <- c(
    "distance", "weights",
    "s1.xCoord", "s1.yCoord", "s2.xCoord", "s2.yCoord",
    "s1.Year", "s1.LU", "s2.Year", "s2.LU"
  )

  # # check dat
  # nplots <- nrow(unique(rbindlist(list(unique(dat[, c("s1.xCoord", "s1.yCoord")]),
  #                                      unique(dat[, c("s2.xCoord", "s2.yCoord")])),
  #                                 use.names = F))) # 150 plots
  # nyears <- length(unique(c(dat$s1.Year, dat$s2.Year)))
  # length(unique(c(dat$s1.LU, dat$s2.LU))) <= nyears * nplots
  # nrow(dat) == (choose((nplots + 2 - 1), 2) - nplots) * nyears # number of rows as expected
  # nrow(unique(rbindlist(list(unique(dat[which(dat$s1.Year == "2008"), c("s1.xCoord", "s1.yCoord", "s1.LU")]),
  #                            unique(dat[which(dat$s2.Year == "2008"), c("s2.xCoord", "s2.yCoord", "s2.LU")])),
  #                       use.names = F))) # 150 plots
  # test <- unique(rbindlist(list(unique(dat[which(dat$s1.Year == "2008"), c("s1.xCoord", "s1.yCoord", "s1.LU", "s1.Year")]),
  #                               unique(dat[which(dat$s2.Year == "2008"), c("s2.xCoord", "s2.yCoord", "s2.LU", "s2.Year")])),
  #                          use.names = F))
  # test <- test[with(test, order(s1.xCoord, s1.yCoord)),]
  # test # expected number of LUI values and rows present
  # #### testing over

  dat <- dat[complete.cases(dat), ]

  detach(pwise_space)

  gdmTab <- formatsitepair(dat, 4, predData = dat)
  var <- gdm.varImp.synthesis.edit(gdmTab,
    geo = T, includeYear = T,
    fullModelOnly = T, nPerm = 50, spacefortimedataset = "spatial"
  )

  var_plant_b3[[i]] <- var
}

names(var_plant_b3) <- LU
var_plant_a3s_sc <- var_plant_b3
save(var_plant_a3s_sc, file = "./2.GDMs/data/OutputData/var_plant_a3s_sc.RData")

# a4
# list for results
var_plant_b4 <- list()

for (i in 1:4) {
  attach(pwise_space)

  # geo dist, LUI and dummy as predictors

  dat <- data.frame(
    da4st,
    rep(1, nrow(pwise_space)),
    HWG1, RWG1, HWG2, RWG2,
    YR,
    pwise_space[, colnames(pwise_space) == LU1[i]],
    YR,
    pwise_space[, colnames(pwise_space) == LU2[i]]
  )

  colnames(dat) <- c(
    "distance", "weights",
    "s1.xCoord", "s1.yCoord", "s2.xCoord", "s2.yCoord",
    "s1.Year", "s1.LU", "s2.Year", "s2.LU"
  )

  # # check dat
  # nplots <- nrow(unique(rbindlist(list(unique(dat[, c("s1.xCoord", "s1.yCoord")]),
  #                                      unique(dat[, c("s2.xCoord", "s2.yCoord")])),
  #                                 use.names = F))) # 150 plots
  # nyears <- length(unique(c(dat$s1.Year, dat$s2.Year)))
  # length(unique(c(dat$s1.LU, dat$s2.LU))) <= nyears * nplots
  # nrow(dat) == (choose((nplots + 2 - 1), 2) - nplots) * nyears # number of rows as expected
  # nrow(unique(rbindlist(list(unique(dat[which(dat$s1.Year == "2008"), c("s1.xCoord", "s1.yCoord", "s1.LU")]),
  #                            unique(dat[which(dat$s2.Year == "2008"), c("s2.xCoord", "s2.yCoord", "s2.LU")])),
  #                       use.names = F))) # 150 plots
  # test <- unique(rbindlist(list(unique(dat[which(dat$s1.Year == "2008"), c("s1.xCoord", "s1.yCoord", "s1.LU", "s1.Year")]),
  #                               unique(dat[which(dat$s2.Year == "2008"), c("s2.xCoord", "s2.yCoord", "s2.LU", "s2.Year")])),
  #                          use.names = F))
  # test <- test[with(test, order(s1.xCoord, s1.yCoord)),]
  # test # expected number of LUI values and rows present
  # #### testing over

  dat <- dat[complete.cases(dat), ]

  detach(pwise_space)

  gdmTab <- formatsitepair(dat, 4, predData = dat)
  var <- gdm.varImp.synthesis.edit(gdmTab,
    geo = T, includeYear = T,
    fullModelOnly = T, nPerm = 50, spacefortimedataset = "spatial"
  )

  var_plant_b4[[i]] <- var
}

names(var_plant_b4) <- LU

var_plant_a4s_sc <- var_plant_b4
save(var_plant_a4s_sc, file = "./2.GDMs/data/OutputData/var_plant_a4s_sc.RData")


#############################################################################
############################################################################
# GDMs time
###########

LU <- c("LUI", "MOW", "GRA", "FER")

LU1 <- c("LUI1res", "MOW1res", "GRA1res", "FER1res")
LU2 <- c("LUI2res", "MOW2res", "GRA2res", "FER2res")



# a0
# list for results
var_plant_b0t <- list()

for (i in 1:4) {
  attach(pwise_time_ayrs)

  # geo dist, LUI and dummy as predictors

  dat <- data.frame(
    da0st,
    rep(1, nrow(pwise_time_ayrs)),
    HWG1, RWG1, HWG2, RWG2,
    YR1,
    (pwise_time_ayrs[, colnames(pwise_time_ayrs) == LU1[i]]),
    YR2,
    (pwise_time_ayrs[, colnames(pwise_time_ayrs) == LU2[i]])
  )


  colnames(dat) <- c(
    "distance", "weights",
    "s1.xCoord", "s1.yCoord",
    "s2.xCoord", "s2.yCoord",
    "s1.Year", "s1.LU",
    "s2.Year", "s2.LU"
  )

  dat <- dat[complete.cases(dat), ]

  detach(pwise_time_ayrs)

  gdmTab <- formatsitepair(dat, 4, predData = dat)
  var <- gdm.varImp.synthesis.edit(gdmTab,
    geo = F, includeYear = T,
    fullModelOnly = T, nPerm = 50, spacefortimedataset = "temporal"
  )

  var_plant_b0t[[i]] <- var
}

names(var_plant_b0t) <- LU
var_plant_a0t_sc <- var_plant_b0t
save(var_plant_a0t_sc, file = "./2.GDMs/data/OutputData/var_plant_a0t_sc.RData")



# a1
# list for results
var_plant_b1t <- list()

for (i in 1:4) {
  attach(pwise_time_ayrs)

  # geo dist, LUI and dummy as predictors

  dat <- data.frame(
    da1st,
    rep(1, nrow(pwise_time_ayrs)),
    HWG1, RWG1, HWG2, RWG2,
    YR1,
    (pwise_time_ayrs[, colnames(pwise_time_ayrs) == LU1[i]]),
    YR2,
    (pwise_time_ayrs[, colnames(pwise_time_ayrs) == LU2[i]])
  )

  colnames(dat) <- c(
    "distance", "weights",
    "s1.xCoord", "s1.yCoord",
    "s2.xCoord", "s2.yCoord",
    "s1.Year", "s1.LU",
    "s2.Year", "s2.LU"
  )

  # ## check dataset
  # nplots <- nrow(unique(rbindlist(list(unique(dat[, c("s1.xCoord", "s1.yCoord")]),
  #                                      unique(dat[, c("s2.xCoord", "s2.yCoord")])),
  #                                 use.names = F))) # 150 plots
  # length(unique(c(paste(dat$s1.xCoord, dat$s1.yCoord), paste(dat$s2.xCoord, dat$s2.yCoord)))) # 150 plots (other method, same result)
  # nyears <- length(unique(c(dat$s1.Year, dat$s2.Year))) # 10 years
  # length(unique(c(dat$s1.LU, dat$s2.LU))) <= nyears * nplots # not too many LUI values
  # length(unique(c(dat$distance, dat$distance)))
  # nrow(dat) <= (choose((nplots + 2 - 1), 2) - nplots) * nyears # spatial dataset
  # ## check dataset end

  dat <- dat[complete.cases(dat), ]

  detach(pwise_time_ayrs)

  gdmTab <- formatsitepair(dat, 4, predData = dat)
  var <- gdm.varImp.synthesis.edit(gdmTab,
    geo = F, includeYear = T,
    fullModelOnly = T, nPerm = 50, spacefortimedataset = "temporal"
  )
  var_plant_b1t[[i]] <- var
}

names(var_plant_b1t) <- LU
var_plant_a1t_sc <- var_plant_b1t
save(var_plant_a1t_sc, file = "./2.GDMs/data/OutputData/var_plant_a1t_sc.RData")



# a2
# list for results
var_plant_b2t <- list()

for (i in 1:4) {
  attach(pwise_time_ayrs)

  # geo dist, LUI and dummy as predictors

  dat <- data.frame(
    da2st,
    rep(1, nrow(pwise_time_ayrs)),
    HWG1, RWG1, HWG2, RWG2,
    YR1,
    (pwise_time_ayrs[, colnames(pwise_time_ayrs) == LU1[i]]),
    YR2,
    (pwise_time_ayrs[, colnames(pwise_time_ayrs) == LU2[i]])
  )

  colnames(dat) <- c(
    "distance", "weights",
    "s1.xCoord", "s1.yCoord",
    "s2.xCoord", "s2.yCoord",
    "s1.Year", "s1.LU",
    "s2.Year", "s2.LU"
  )

  # ## check dataset
  # nplots <- nrow(unique(rbindlist(list(unique(dat[, c("s1.xCoord", "s1.yCoord")]),
  #                                      unique(dat[, c("s2.xCoord", "s2.yCoord")])),
  #                                 use.names = F))) # 150 plots
  # length(unique(c(paste(dat$s1.xCoord, dat$s1.yCoord), paste(dat$s2.xCoord, dat$s2.yCoord)))) # 150 plots (other method, same result)
  # nyears <- length(unique(c(dat$s1.Year, dat$s2.Year))) # 10 years
  # length(unique(c(dat$s1.LU, dat$s2.LU))) <= nyears * nplots # not too many LUI values
  # length(unique(c(dat$distance, dat$distance)))
  # nrow(dat) <= (choose((nplots + 2 - 1), 2) - nplots) * nyears # spatial dataset
  # ## check dataset end

  dat <- dat[complete.cases(dat), ]

  detach(pwise_time_ayrs)

  gdmTab <- formatsitepair(dat, 4, predData = dat)
  var <- gdm.varImp.synthesis.edit(gdmTab,
    geo = F, includeYear = T,
    fullModelOnly = T, nPerm = 50, spacefortimedataset = "temporal"
  )

  var_plant_b2t[[i]] <- var
}

names(var_plant_b2t) <- LU
var_plant_a2t_sc <- var_plant_b2t
save(var_plant_a2t_sc, file = "./2.GDMs/data/OutputData/var_plant_a2t_sc.RData")



# a3
# list for results
var_plant_b3t <- list()

for (i in 1:4) {
  attach(pwise_time_ayrs)

  # geo dist, LUI and dummy as predictors
  dat <- data.frame(
    da3st,
    rep(1, nrow(pwise_time_ayrs)),
    HWG1, RWG1, HWG2, RWG2,
    YR1,
    (pwise_time_ayrs[, colnames(pwise_time_ayrs) == LU1[i]]),
    YR2,
    (pwise_time_ayrs[, colnames(pwise_time_ayrs) == LU2[i]])
  )
  colnames(dat) <- c(
    "distance", "weights",
    "s1.xCoord", "s1.yCoord",
    "s2.xCoord", "s2.yCoord",
    "s1.Year", "s1.LU",
    "s2.Year", "s2.LU"
  )

  # ## check dataset
  # nplots <- nrow(unique(rbindlist(list(unique(dat[, c("s1.xCoord", "s1.yCoord")]),
  #                                      unique(dat[, c("s2.xCoord", "s2.yCoord")])),
  #                                 use.names = F))) # 150 plots
  # length(unique(c(paste(dat$s1.xCoord, dat$s1.yCoord), paste(dat$s2.xCoord, dat$s2.yCoord)))) # 150 plots (other method, same result)
  # nyears <- length(unique(c(dat$s1.Year, dat$s2.Year))) # 10 years
  # length(unique(c(dat$s1.LU, dat$s2.LU))) <= nyears * nplots # not too many LUI values
  # length(unique(c(dat$distance, dat$distance)))
  # nrow(dat) <= (choose((nplots + 2 - 1), 2) - nplots) * nyears # spatial dataset
  # ## check dataset end

  dat <- dat[complete.cases(dat), ]

  detach(pwise_time_ayrs)

  gdmTab <- formatsitepair(dat, 4, predData = dat)
  var <- gdm.varImp.synthesis.edit(gdmTab,
    geo = F, includeYear = T,
    fullModelOnly = T, nPerm = 50, spacefortimedataset = "temporal"
  )

  var_plant_b3t[[i]] <- var
}

names(var_plant_b3t) <- LU
var_plant_a3t_sc <- var_plant_b3t
save(var_plant_a3t_sc, file = "./2.GDMs/data/OutputData/var_plant_a3t_sc.RData")



# b4
# list for results
var_plant_b4t <- list()

for (i in 1:4) {
  attach(pwise_time_ayrs)

  # geo dist, LUI and dummy as predictors
  dat <- data.frame(
    da4st,
    rep(1, nrow(pwise_time_ayrs)),
    HWG1, RWG1, HWG2, RWG2,
    YR1,
    (pwise_time_ayrs[, colnames(pwise_time_ayrs) == LU1[i]]),
    YR2,
    (pwise_time_ayrs[, colnames(pwise_time_ayrs) == LU2[i]])
  )

  colnames(dat) <- c(
    "distance", "weights",
    "s1.xCoord", "s1.yCoord",
    "s2.xCoord", "s2.yCoord",
    "s1.Year", "s1.LU",
    "s2.Year", "s2.LU"
  )

  # ## check dataset
  # nplots <- nrow(unique(rbindlist(list(unique(dat[, c("s1.xCoord", "s1.yCoord")]),
  #                                      unique(dat[, c("s2.xCoord", "s2.yCoord")])),
  #                                 use.names = F))) # 150 plots
  # length(unique(c(paste(dat$s1.xCoord, dat$s1.yCoord), paste(dat$s2.xCoord, dat$s2.yCoord)))) # 150 plots (other method, same result)
  # nyears <- length(unique(c(dat$s1.Year, dat$s2.Year))) # 10 years
  # length(unique(c(dat$s1.LU, dat$s2.LU))) <= nyears * nplots # not too many LUI values
  # length(unique(c(dat$distance, dat$distance)))
  # nrow(dat) <= (choose((nplots + 2 - 1), 2) - nplots) * nyears # spatial dataset
  # ## check dataset end

  dat <- dat[complete.cases(dat), ]

  detach(pwise_time_ayrs)

  gdmTab <- formatsitepair(dat, 4, predData = dat)
  var <- gdm.varImp.synthesis.edit(gdmTab,
    geo = F, includeYear = T,
    fullModelOnly = T, nPerm = 50, spacefortimedataset = "temporal"
  )
  var_plant_b4t[[i]] <- var
}

names(var_plant_b4t) <- LU
var_plant_a4t_sc <- var_plant_b4t
save(var_plant_a4t_sc, file = "./2.GDMs/data/OutputData/var_plant_a4t_sc.RData")


remove(pwise_space, pwise_time_ayrs)