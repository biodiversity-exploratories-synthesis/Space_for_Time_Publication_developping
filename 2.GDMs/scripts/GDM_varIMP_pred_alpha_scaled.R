######################################################
#       GDM p-values for insect secondary consumers/predators (space and time)
######### differences in alpha diversity ########
######################################################

# set working directory to folder "Space_for_Time_Publication"

# upload data
# secondary consumers/predators
load("./2.GDMs/data/InputData/pwise_time_pred.RData")
load("./2.GDMs/data/InputData/pwise_space_pred.RData")

pwise_time_ayrs_in<- pwise_time_pred
pwise_space_in<- pwise_space_pred

# if the uploaded, assembled data files are used, upload the complete insect files
# here.


# packages
library(tidyverse)
library(gdm)
library(data.table)
library(vegan)

# edited var.Imp function
source("RFunctions/gdm.varIMP_edit.R")
source("RFunctions/matrix_perm_permutateSitePair.R")

pwise_space_in$dpa0st <- decostand(pwise_space_in$dpa0, method = "range", na.rm = T)
pwise_space_in$dpa1st <- decostand(pwise_space_in$dpa1, method = "range", na.rm = T)
pwise_space_in$dpa2st <- decostand(pwise_space_in$dpa2, method = "range", na.rm = T)
pwise_space_in$dpa3st <- decostand(pwise_space_in$dpa3, method = "range", na.rm = T)
pwise_space_in$dpa4st <- decostand(pwise_space_in$dpa4, method = "range", na.rm = T)

pwise_time_ayrs_in$dpa0st <- decostand(pwise_time_ayrs_in$dpa0abs, method = "range", na.rm = T)
pwise_time_ayrs_in$dpa1st <- decostand(pwise_time_ayrs_in$dpa1abs, method = "range", na.rm = T)
pwise_time_ayrs_in$dpa2st <- decostand(pwise_time_ayrs_in$dpa2abs, method = "range", na.rm = T)
pwise_time_ayrs_in$dpa3st <- decostand(pwise_time_ayrs_in$dpa3abs, method = "range", na.rm = T)
pwise_time_ayrs_in$dpa4st <- decostand(pwise_time_ayrs_in$dpa4abs, method = "range", na.rm = T)


# GDMs space
LU <- c("LUI", "MOW", "GRA", "FER")
LU1 <- c("LUI1res", "MOW1res", "GRA1res", "FER1res")
LU2 <- c("LUI2res", "MOW2res", "GRA2res", "FER2res")



# a0
# list for results
var_pred_b0 <- list()

for (i in 1:4) {
  attach(pwise_space_in)

  # geo dist, LUI and dummy as predictors
  dat <- data.frame(
    dpa0st,
    rep(1, nrow(pwise_space_in)),
    HWG1, RWG1, HWG2, RWG2,
    YR,
    pwise_space_in[, colnames(pwise_space_in) == LU1[i]],
    YR,
    pwise_space_in[, colnames(pwise_space_in) == LU2[i]]
  )

  colnames(dat) <- c(
    "distance", "weights",
    "s1.xCoord", "s1.yCoord", "s2.xCoord", "s2.yCoord",
    "s1.Year", "s1.LU", "s2.Year", "s2.LU"
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

  detach(pwise_space_in)

  gdmTab <- formatsitepair(dat, 4, predData = dat)
  var <- gdm.varImp.synthesis.edit(gdmTab,
    geo = T, includeYear = T,
    fullModelOnly = T, nPerm = 50, spacefortimedataset = "spatial"
  )

  var_pred_b0[[i]] <- var

  remove(var)
}

names(var_pred_b0) <- LU
var_pred_a0s_sc <- var_pred_b0
save(var_pred_a0s_sc, file = "./2.GDMs/data/OutputData/var_pred_a0s_sc.RData")


# a1
# list for results
var_pred_b1 <- list()

for (i in 1:4) {
  attach(pwise_space_in)

  # geo dist, LUI and dummy as predictors
  dat <- data.frame(
    dpa1st,
    rep(1, nrow(pwise_space_in)),
    HWG1, RWG1, HWG2, RWG2,
    YR,
    pwise_space_in[, colnames(pwise_space_in) == LU1[i]],
    YR,
    pwise_space_in[, colnames(pwise_space_in) == LU2[i]]
  )

  colnames(dat) <- c(
    "distance", "weights",
    "s1.xCoord", "s1.yCoord", "s2.xCoord", "s2.yCoord",
    "s1.Year", "s1.LU", "s2.Year", "s2.LU"
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

  detach(pwise_space_in)

  gdmTab <- formatsitepair(dat, 4, predData = dat)
  var <- gdm.varImp.synthesis.edit(gdmTab,
    geo = T, includeYear = T,
    fullModelOnly = T, nPerm = 50, spacefortimedataset = "spatial"
  )

  var_pred_b1[[i]] <- var

  remove(var)
}

names(var_pred_b1) <- LU
var_pred_a1s_sc <- var_pred_b1
save(var_pred_a1s_sc, file = "./2.GDMs/data/OutputData/var_pred_a1s_sc.RData")



# a2
# list for results
var_pred_b2 <- list()

for (i in 1:4) {
  attach(pwise_space_in)

  # geo dist, LUI and dummy as predictors
  dat <- data.frame(
    dpa2st,
    rep(1, nrow(pwise_space_in)),
    HWG1, RWG1, HWG2, RWG2,
    YR,
    pwise_space_in[, colnames(pwise_space_in) == LU1[i]],
    YR,
    pwise_space_in[, colnames(pwise_space_in) == LU2[i]]
  )

  colnames(dat) <- c(
    "distance", "weights",
    "s1.xCoord", "s1.yCoord", "s2.xCoord", "s2.yCoord",
    "s1.Year", "s1.LU", "s2.Year", "s2.LU"
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

  detach(pwise_space_in)

  gdmTab <- formatsitepair(dat, 4, predData = dat)
  var <- gdm.varImp.synthesis.edit(gdmTab,
    geo = T, includeYear = T,
    fullModelOnly = T, nPerm = 50, spacefortimedataset = "spatial"
  )

  var_pred_b2[[i]] <- var
  remove(var)
}

names(var_pred_b2) <- LU
var_pred_a2s_sc <- var_pred_b2
save(var_pred_a2s_sc, file = "./2.GDMs/data/OutputData/var_pred_a2s_sc.RData")



# a3
# list for results
var_pred_b3 <- list()

for (i in 1:4) {
  attach(pwise_space_in)

  # geo dist, LUI and dummy as predictors
  dat <- data.frame(
    dpa3st,
    rep(1, nrow(pwise_space_in)),
    HWG1, RWG1, HWG2, RWG2,
    YR,
    pwise_space_in[, colnames(pwise_space_in) == LU1[i]],
    YR,
    pwise_space_in[, colnames(pwise_space_in) == LU2[i]]
  )

  colnames(dat) <- c(
    "distance", "weights",
    "s1.xCoord", "s1.yCoord", "s2.xCoord", "s2.yCoord",
    "s1.Year", "s1.LU", "s2.Year", "s2.LU"
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

  detach(pwise_space_in)

  gdmTab <- formatsitepair(dat, 4, predData = dat)
  var <- gdm.varImp.synthesis.edit(gdmTab,
    geo = T, includeYear = T,
    fullModelOnly = T, nPerm = 50, spacefortimedataset = "spatial"
  )

  var_pred_b3[[i]] <- var
  remove(var)
}

names(var_pred_b3) <- LU
var_pred_a3s_sc <- var_pred_b3
save(var_pred_a3s_sc, file = "./2.GDMs/data/OutputData/var_pred_a3s_sc.RData")



# a4
# list for results
var_pred_b4 <- list()

for (i in 1:4) {
  attach(pwise_space_in)

  # geo dist, LUI and dummy as predictors
  dat <- data.frame(
    dpa4st,
    rep(1, nrow(pwise_space_in)),
    HWG1, RWG1, HWG2, RWG2,
    YR,
    pwise_space_in[, colnames(pwise_space_in) == LU1[i]],
    YR,
    pwise_space_in[, colnames(pwise_space_in) == LU2[i]]
  )

  colnames(dat) <- c(
    "distance", "weights",
    "s1.xCoord", "s1.yCoord", "s2.xCoord", "s2.yCoord",
    "s1.Year", "s1.LU", "s2.Year", "s2.LU"
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

  detach(pwise_space_in)

  gdmTab <- formatsitepair(dat, 4, predData = dat)
  var <- gdm.varImp.synthesis.edit(gdmTab,
    geo = T, includeYear = T,
    fullModelOnly = T, nPerm = 50, spacefortimedataset = "spatial"
  )

  var_pred_b4[[i]] <- var
  remove(var)
}

names(var_pred_b4) <- LU

var_pred_a4s_sc <- var_pred_b4
save(var_pred_a4s_sc, file = "./2.GDMs/data/OutputData/var_pred_a4s_sc.RData")


#############################################################################
############################################################################
# GDMs time
###########

LU <- c("LUI", "MOW", "GRA", "FER")
LU1 <- c("LUI1res", "MOW1res", "GRA1res", "FER1res")
LU2 <- c("LUI2res", "MOW2res", "GRA2res", "FER2res")



# a0
# list for results
var_pred_b0t <- list()

for (i in 1:4) {
  attach(pwise_time_ayrs_in)

  # geo dist, LUI and dummy as predictors
  dat <- data.frame(
    dpa0st,
    rep(1, nrow(pwise_time_ayrs_in)),
    HWG1, RWG1, HWG2, RWG2,
    YR1,
    (pwise_time_ayrs_in[, colnames(pwise_time_ayrs_in) == LU1[i]]),
    YR2,
    (pwise_time_ayrs_in[, colnames(pwise_time_ayrs_in) == LU2[i]])
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
  # nrow(dat) <= (choose((nyears + 2 - 1), 2) - nyears) * nplots # temporal dataset
  # ## check dataset end

  dat <- dat[complete.cases(dat), ]

  detach(pwise_time_ayrs_in)

  gdmTab <- formatsitepair(dat, 4, predData = dat)
  var <- gdm.varImp.synthesis.edit(gdmTab,
    geo = F, includeYear = T,
    fullModelOnly = T, nPerm = 50, spacefortimedataset = "temporal"
  )

  var_pred_b0t[[i]] <- var
  remove(var)
}

names(var_pred_b0t) <- LU
var_pred_a0t_sc <- var_pred_b0t
save(var_pred_a0t_sc, file = "./2.GDMs/data/OutputData/var_pred_a0t_sc.RData")



# a1
# list for results
var_pred_b1t <- list()

for (i in 1:4) {
  attach(pwise_time_ayrs_in)

  # geo dist, LUI and dummy as predictors
  dat <- data.frame(
    dpa1st,
    rep(1, nrow(pwise_time_ayrs_in)),
    HWG1, RWG1, HWG2, RWG2,
    YR1,
    (pwise_time_ayrs_in[, colnames(pwise_time_ayrs_in) == LU1[i]]),
    YR2,
    (pwise_time_ayrs_in[, colnames(pwise_time_ayrs_in) == LU2[i]])
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

  detach(pwise_time_ayrs_in)

  gdmTab <- formatsitepair(dat, 4, predData = dat)
  var <- gdm.varImp.synthesis.edit(gdmTab,
    geo = F, includeYear = T,
    fullModelOnly = T, nPerm = 50, spacefortimedataset = "temporal"
  )
  var_pred_b1t[[i]] <- var
  remove(var)
}

names(var_pred_b1t) <- LU
var_pred_a1t_sc <- var_pred_b1t
save(var_pred_a1t_sc, file = "./2.GDMs/data/OutputData/var_pred_a1t_sc.RData")



# a2
# list for results
var_pred_b2t <- list()

for (i in 1:4) {
  attach(pwise_time_ayrs_in)

  # geo dist, LUI and dummy as predictors
  dat <- data.frame(
    dpa2st,
    rep(1, nrow(pwise_time_ayrs_in)),
    HWG1, RWG1, HWG2, RWG2,
    YR1,
    (pwise_time_ayrs_in[, colnames(pwise_time_ayrs_in) == LU1[i]]),
    YR2,
    (pwise_time_ayrs_in[, colnames(pwise_time_ayrs_in) == LU2[i]])
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
  # nrow(dat) <= (choose((nplots + 2 - 1), 2) - nplots) * nyears # temporal dataset
  # ## check dataset end

  dat <- dat[complete.cases(dat), ]

  detach(pwise_time_ayrs_in)

  gdmTab <- formatsitepair(dat, 4, predData = dat)
  var <- gdm.varImp.synthesis.edit(gdmTab,
    geo = F, includeYear = T,
    fullModelOnly = T, nPerm = 50, spacefortimedataset = "temporal"
  )

  var_pred_b2t[[i]] <- var
  remove(var)
}

names(var_pred_b2t) <- LU
var_pred_a2t_sc <- var_pred_b2t
save(var_pred_a2t_sc, file = "./2.GDMs/data/OutputData/var_pred_a2t_sc.RData")


# a3
# list for results
var_pred_b3t <- list()

for (i in 1:4) {
  attach(pwise_time_ayrs_in)

  # geo dist, LUI and dummy as predictors
  dat <- data.frame(
    dpa3st,
    rep(1, nrow(pwise_time_ayrs_in)),
    HWG1, RWG1, HWG2, RWG2,
    YR1,
    (pwise_time_ayrs_in[, colnames(pwise_time_ayrs_in) == LU1[i]]),
    YR2,
    (pwise_time_ayrs_in[, colnames(pwise_time_ayrs_in) == LU2[i]])
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
  # nrow(dat) <= (choose((nplots + 2 - 1), 2) - nplots) * nyears # temporal dataset
  # ## check dataset end

  dat <- dat[complete.cases(dat), ]

  detach(pwise_time_ayrs_in)

  gdmTab <- formatsitepair(dat, 4, predData = dat)
  var <- gdm.varImp.synthesis.edit(gdmTab,
    geo = F, includeYear = T,
    fullModelOnly = T, nPerm = 50, spacefortimedataset = "temporal"
  )

  var_pred_b3t[[i]] <- var
  remove(var)
}

names(var_pred_b3t) <- LU
var_pred_a3t_sc <- var_pred_b3t
save(var_pred_a3t_sc, file = "./2.GDMs/data/OutputData/var_pred_a3t_sc.RData")


# a4
# list for results
var_pred_b4t <- list()

for (i in 1:4) {
  attach(pwise_time_ayrs_in)

  # geo dist, LUI and dummy as predictors
  dat <- data.frame(
    dpa4st,
    rep(1, nrow(pwise_time_ayrs_in)),
    HWG1, RWG1, HWG2, RWG2,
    YR1,
    (pwise_time_ayrs_in[, colnames(pwise_time_ayrs_in) == LU1[i]]),
    YR2,
    (pwise_time_ayrs_in[, colnames(pwise_time_ayrs_in) == LU2[i]])
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
  # nrow(dat) <= (choose((nplots + 2 - 1), 2) - nplots) * nyears # temporal dataset
  # ## check dataset end

  dat <- dat[complete.cases(dat), ]

  detach(pwise_time_ayrs_in)

  gdmTab <- formatsitepair(dat, 4, predData = dat)
  var <- gdm.varImp.synthesis.edit(gdmTab,
    geo = F, includeYear = T,
    fullModelOnly = T, nPerm = 50, spacefortimedataset = "temporal"
  )

  var_pred_b4t[[i]] <- var
  remove(var)
}

names(var_pred_b4t) <- LU
var_pred_a4t_sc <- var_pred_b4t
save(var_pred_a4t_sc, file = "./2.GDMs/data/OutputData/var_pred_a4t_sc.RData")

remove(pwise_space_in, pwise_time_ayrs_in)

