######################################################
#       GDM p-values for plants beta diversity (space and time)
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

# edited var.Imp function
source("./RFunctions/gdm.varIMP_edit.R")
source("RFunctions/matrix_perm_permutateSitePair.R")


# GDMs space
LU <- c("LUI", "MOW", "GRA", "FER")

# list for results
var_plants_b0 <- list()
var_plants_b0_sub <- list() # only within region plot pairs


# b0
# within region
for (i in 1:4) {
  LU1 <- str_c(LU[i], 1)
  LU2 <- str_c(LU[i], 2)

  attach(pwise_space)

  # geo dist, LUI

  dat <- data.frame(
    cqn0dis,
    rep(1, nrow(pwise_space)),
    HWG1, RWG1, HWG2, RWG2,
    YR,
    pwise_space[, colnames(pwise_space) == LU1],
    YR,
    pwise_space[, colnames(pwise_space) == LU2]
  )

  colnames(dat) <- c(
    "distance", "weights",
    "s1.xCoord", "s1.yCoord", "s2.xCoord", "s2.yCoord",
    "s1.Year", "s1.LU", "s2.Year", "s2.LU"
  )

  dat <- dat[complete.cases(dat), ]

  # ## check dataset
  # nplots <- nrow(unique(rbindlist(list(unique(dat[, c("s1.xCoord", "s1.yCoord")]),
  #                                      unique(dat[, c("s2.xCoord", "s2.yCoord")])),
  #                                 use.names = F))) # 150 plots
  # length(unique(c(paste(dat$s1.xCoord, dat$s1.yCoord), paste(dat$s2.xCoord, dat$s2.yCoord)))) # 150 plots (other method, same result)
  # nyears <- length(unique(c(dat$s1.Year, dat$s2.Year))) # 11 years
  # length(unique(c(dat$s1.LU, dat$s2.LU))) <= nyears * nplots # expected number of LUI values
  # length(unique(c(dat$distance, dat$distance)))
  # nrow(dat) <= (choose((nplots + 2 - 1), 2) - nplots) * nyears # spatial dataset : expected number of rows present
  # ## check dataset end

  detach(pwise_space)

  gdmTab <- formatsitepair(dat, 4, predData = dat)
  var <- gdm.varImp.synthesis.edit(gdmTab,
    geo = T, includeYear = T, nPerm = 50,
    fullModelOnly = T, spacefortimedataset = "spatial"
  )

  var_plants_b0_sub[[i]] <- var
}

names(var_plants_b0_sub) <- LU
save(var_plants_b0_sub, file = "./2.GDMs/data/OutputData/var_plants_b0_sub.RData")



# b1
# all regions
var_plants_b1 <- list()

for (i in 1:4) {
  LU1 <- str_c(LU[i], 1)
  LU2 <- str_c(LU[i], 2)

  attach(pwise_space)

  # geo dist, LUI and dummy as predictors

  dat <- data.frame(
    cqn1dis - 0.00001,
    rep(1, nrow(pwise_space)),
    HWG1, RWG1, HWG2, RWG2,
    YR,
    pwise_space[, colnames(pwise_space) == LU1],
    YR,
    pwise_space[, colnames(pwise_space) == LU2]
  )

  colnames(dat) <- c(
    "distance", "weights",
    "s1.xCoord", "s1.yCoord", "s2.xCoord", "s2.yCoord",
    "s1.Year", "s1.LU", "s2.Year", "s2.LU"
  )
  # # check dat
  # nplots <- nrow(unique(rbindlist(list(unique(dat[, c("s1.xCoord", "s1.yCoord")]),
  #                            unique(dat[, c("s2.xCoord", "s2.yCoord")])),
  #                       use.names = F))) # 150 plots
  # length(unique(c(paste(dat$s1.xCoord, dat$s1.yCoord), paste(dat$s2.xCoord, dat$s2.yCoord)))) # 150 plots (other method, same result)
  # nyears <- length(unique(c(dat$s1.Year, dat$s2.Year)))
  # length(unique(c(dat$s1.LU, dat$s2.LU))) <= nyears * nplots # expected number of LUI values
  # nrow(dat) == (choose((nplots + 2 - 1), 2) -nplots) * nyears # number of rows as expected in spatial dataset
  # ## check dataset end

  dat <- dat[complete.cases(dat), ]

  detach(pwise_space)

  gdmTab <- formatsitepair(dat, 4, predData = dat)
  var <- gdm.varImp.synthesis.edit(gdmTab,
    geo = T, includeYear = T, nPerm = 50,
    fullModelOnly = T, spacefortimedataset = "spatial"
  )

  var_plants_b1[[i]] <- var
}

names(var_plants_b1) <- LU

save(var_plants_b1, file = "./2.GDMs/data/OutputData/var_plants_b1.RData")


# b2
var_plants_b2 <- list()

for (i in 1:4) {
  LU1 <- str_c(LU[i], 1)
  LU2 <- str_c(LU[i], 2)

  attach(pwise_space)

  # geo dist, LUI and dummy as predictors

  dat <- data.frame(
    cqn2dis,
    rep(1, nrow(pwise_space)),
    HWG1, RWG1, HWG2, RWG2,
    YR,
    pwise_space[, colnames(pwise_space) == LU1],
    YR,
    pwise_space[, colnames(pwise_space) == LU2]
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
  # nrow(dat) <= (choose((nplots + 2 - 1), 2) - nplots) * nyears # spatial dataset : expected number of plots presen
  # ## check dataset end

  dat <- dat[complete.cases(dat), ]

  detach(pwise_space)

  gdmTab <- formatsitepair(dat, 4, predData = dat)
  var <- gdm.varImp.synthesis.edit(gdmTab,
    geo = T, includeYear = T, nPerm = 50,
    fullModelOnly = T, spacefortimedataset = "spatial"
  )

  var_plants_b2[[i]] <- var
}

names(var_plants_b2) <- LU

save(var_plants_b2, file = "./2.GDMs/data/OutputData/var_plants_b2.RData")

# b3
var_plants_b3 <- list()

for (i in 1:4) {
  LU1 <- str_c(LU[i], 1)
  LU2 <- str_c(LU[i], 2)

  attach(pwise_space)

  # geo dist, LUI

  dat <- data.frame(
    cqn3dis,
    rep(1, nrow(pwise_space)),
    HWG1, RWG1, HWG2, RWG2,
    YR,
    pwise_space[, colnames(pwise_space) == LU1],
    YR,
    pwise_space[, colnames(pwise_space) == LU2]
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
  # nrow(dat) <= (choose((nplots + 2 - 1), 2) - nplots) * nyears # spatial dataset : expected number of plots presen
  # ## check dataset end

  dat <- dat[complete.cases(dat), ]

  detach(pwise_space)

  gdmTab <- formatsitepair(dat, 4, predData = dat)
  var <- gdm.varImp.synthesis.edit(gdmTab,
    geo = T, includeYear = T, nPerm = 50,
    fullModelOnly = T, spacefortimedataset = "spatial"
  )

  var_plants_b3[[i]] <- var
}

names(var_plants_b3) <- LU
save(var_plants_b3, file = "./2.GDMs/data/OutputData/var_plants_b3.RData")


# b4

var_plants_b4 <- list()


for (i in 1:4) {
  LU1 <- str_c(LU[i], 1)
  LU2 <- str_c(LU[i], 2)

  attach(pwise_space)

  # geo dist, LUI

  dat <- data.frame(
    cqn4dis,
    rep(1, nrow(pwise_space)),
    HWG1, RWG1, HWG2, RWG2,
    YR,
    pwise_space[, colnames(pwise_space) == LU1],
    YR,
    pwise_space[, colnames(pwise_space) == LU2]
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
  # nrow(dat) <= (choose((nplots + 2 - 1), 2) - nplots) * nyears # spatial dataset : expected number of plots presen
  # ## check dataset end

  dat <- dat[complete.cases(dat), ]

  detach(pwise_space)

  gdmTab <- formatsitepair(dat, 4, predData = dat)
  var <- gdm.varImp.synthesis.edit(gdmTab,
    geo = T, includeYear = T, nPerm = 50,
    fullModelOnly = T, spacefortimedataset = "spatial"
  )

  var_plants_b4[[i]] <- var
}

names(var_plants_b4) <- LU
save(var_plants_b4, file = "./2.GDMs/data/OutputData/var_plants_b4.RData")


# bsim
var_plants_bsim <- list()

start.time <- Sys.time()
for (i in 1:4) {
  LU1 <- str_c(LU[i], 1)
  LU2 <- str_c(LU[i], 2)

  attach(pwise_space)

  # geo dist, LUI and dummy as predictors

  dat <- data.frame(
    bsim,
    rep(1, nrow(pwise_space)),
    HWG1, RWG1, HWG2, RWG2,
    YR,
    pwise_space[, colnames(pwise_space) == LU1],
    YR,
    pwise_space[, colnames(pwise_space) == LU2]
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
  # nrow(dat) <= (choose((nplots + 2 - 1), 2) - nplots) * nyears # spatial dataset : expected number of plots presen
  # ## check dataset end

  dat <- dat[complete.cases(dat), ]

  detach(pwise_space)

  gdmTab <- formatsitepair(dat, 4, predData = dat)
  var <- gdm.varImp.synthesis.edit(gdmTab,
    geo = T, includeYear = T, nPerm = 50,
    fullModelOnly = T, spacefortimedataset = "spatial"
  )

  var_plants_bsim[[i]] <- var
}
end.time <- Sys.time()

names(var_plants_bsim) <- LU

save(var_plants_bsim, file = "./2.GDMs/data/OutputData/var_plants_bsim.RData")


#############################################################################
############################################################################
# GDMs time
###########

LU <- c("LUI", "MOW", "GRA", "FER")

# list for results
var_plants_b0t <- list()


# b0
for (i in 1:4) {
  LU1 <- str_c(LU[i], 1)
  LU2 <- str_c(LU[i], 2)

  attach(pwise_time_ayrs)

  # geo dist, LUI and dummy as predictors

  dat <- data.frame(
    cqn0dis,
    rep(1, nrow(pwise_time_ayrs)),
    HWG1, RWG1, HWG2, RWG2,
    YR1,
    pwise_time_ayrs[, colnames(pwise_time_ayrs) == LU1],
    YR2,
    pwise_time_ayrs[, colnames(pwise_time_ayrs) == LU2]
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
  # nrow(dat) <= (choose((nyears + 2 - 1), 2) - nyears) * nplots # temporal dataset : expected number of plots presen
  # ## check dataset end

  dat <- dat[complete.cases(dat), ]

  detach(pwise_time_ayrs)

  gdmTab <- formatsitepair(dat, 4, predData = dat)
  var <- gdm.varImp.synthesis.edit(gdmTab,
    geo = F, includeYear = T, nPerm = 50,
    fullModelOnly = T, spacefortimedataset = "temporal"
  )

  var_plants_b0t[[i]] <- var
}

names(var_plants_b0t) <- LU
save(var_plants_b0t, file = "./2.GDMs/data/OutputData/var_plants_b0t.RData")



# b1
var_plants_b1t <- list()

for (i in 1:4) {
  LU1 <- str_c(LU[i], 1)
  LU2 <- str_c(LU[i], 2)

  attach(pwise_time_ayrs)

  # geo dist, LUI and dummy as predictors

  dat <- data.frame(
    cqn1dis - 0.00001,
    rep(1, nrow(pwise_time_ayrs)),
    HWG1, RWG1, HWG2, RWG2,
    YR1,
    pwise_time_ayrs[, colnames(pwise_time_ayrs) == LU1],
    YR2,
    pwise_time_ayrs[, colnames(pwise_time_ayrs) == LU2]
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
  # nrow(dat) <= (choose((nyears + 2 - 1), 2) - nyears) * nplots # temporal dataset : expected number of plots presen
  # ## check dataset end

  dat <- dat[complete.cases(dat), ]

  detach(pwise_time_ayrs)

  gdmTab <- formatsitepair(dat, 4, predData = dat)
  var <- gdm.varImp.synthesis.edit(gdmTab,
    geo = F, includeYear = T, nPerm = 50,
    fullModelOnly = T, spacefortimedataset = "temporal"
  )

  var_plants_b1t[[i]] <- var
}

names(var_plants_b1t) <- LU

save(var_plants_b1t, file = "./2.GDMs/data/OutputData/var_plants_b1t.RData")


# b2
var_plants_b2t <- list()

for (i in 1:4) {
  LU1 <- str_c(LU[i], 1)
  LU2 <- str_c(LU[i], 2)

  attach(pwise_time_ayrs)

  # geo dist, LUI and dummy as predictors

  dat <- data.frame(
    cqn2dis,
    rep(1, nrow(pwise_time_ayrs)),
    HWG1, RWG1, HWG2, RWG2,
    YR1,
    pwise_time_ayrs[, colnames(pwise_time_ayrs) == LU1],
    YR2,
    pwise_time_ayrs[, colnames(pwise_time_ayrs) == LU2]
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
  # nrow(dat) <= (choose((nyears + 2 - 1), 2) - nyears) * nplots # temporal dataset : expected number of plots presen
  # ## check dataset end

  dat <- dat[complete.cases(dat), ]

  detach(pwise_time_ayrs)

  gdmTab <- formatsitepair(dat, 4, predData = dat)
  var <- gdm.varImp.synthesis.edit(gdmTab,
    geo = F, includeYear = T, nPerm = 50,
    fullModelOnly = T, spacefortimedataset = "temporal"
  )

  var_plants_b2t[[i]] <- var
}


names(var_plants_b2t) <- LU

save(var_plants_b2t, file = "./2.GDMs/data/OutputData/var_plants_b2t.RData")

# b3
var_plants_b3t <- list()

for (i in 1:4) {
  LU1 <- str_c(LU[i], 1)
  LU2 <- str_c(LU[i], 2)

  attach(pwise_time_ayrs)

  # geo dist, LUI and dummy as predictors

  dat <- data.frame(
    cqn3dis,
    rep(1, nrow(pwise_time_ayrs)),
    HWG1, RWG1, HWG2, RWG2,
    YR1,
    pwise_time_ayrs[, colnames(pwise_time_ayrs) == LU1],
    YR2,
    pwise_time_ayrs[, colnames(pwise_time_ayrs) == LU2]
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
  # nrow(dat) <= (choose((nyears + 2 - 1), 2) - nyears) * nplots # temporal dataset : expected number of plots presen
  # ## check dataset end

  dat <- dat[complete.cases(dat), ]

  detach(pwise_time_ayrs)

  gdmTab <- formatsitepair(dat, 4, predData = dat)
  var <- gdm.varImp.synthesis.edit(gdmTab,
    geo = F, includeYear = T, nPerm = 50,
    fullModelOnly = T, spacefortimedataset = "temporal"
  )

  var_plants_b3t[[i]] <- var
}


names(var_plants_b3t) <- LU

save(var_plants_b3t, file = "./2.GDMs/data/OutputData/var_plants_b3t.RData")

# b4
var_plants_b4t <- list()

for (i in 1:4) {
  LU1 <- str_c(LU[i], 1)
  LU2 <- str_c(LU[i], 2)

  attach(pwise_time_ayrs)

  # geo dist, LUI and dummy as predictors

  dat <- data.frame(
    cqn4dis,
    rep(1, nrow(pwise_time_ayrs)),
    HWG1, RWG1, HWG2, RWG2,
    YR1,
    pwise_time_ayrs[, colnames(pwise_time_ayrs) == LU1],
    YR2,
    pwise_time_ayrs[, colnames(pwise_time_ayrs) == LU2]
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
  # nrow(dat) <= (choose((nyears + 2 - 1), 2) - nyears) * nplots # temporal dataset : expected number of plots presen
  # ## check dataset end

  dat <- dat[complete.cases(dat), ]

  detach(pwise_time_ayrs)

  gdmTab <- formatsitepair(dat, 4, predData = dat)
  var <- gdm.varImp.synthesis.edit(gdmTab,
    geo = F, includeYear = T, nPerm = 50,
    fullModelOnly = T, spacefortimedataset = "temporal"
  )

  var_plants_b4t[[i]] <- var
}


names(var_plants_b4t) <- LU

save(var_plants_b4t, file = "./2.GDMs/data/OutputData/var_plants_b4t.RData")

# bsim
var_plants_bsimt <- list()

start.time <- Sys.time()
for (i in 1:4) {
  LU1 <- str_c(LU[i], 1)
  LU2 <- str_c(LU[i], 2)

  attach(pwise_time_ayrs)

  # geo dist, LUI and dummy as predictors

  dat <- data.frame(
    bsim,
    rep(1, nrow(pwise_time_ayrs)),
    HWG1, RWG1, HWG2, RWG2,
    YR1,
    pwise_time_ayrs[, colnames(pwise_time_ayrs) == LU1],
    YR2,
    pwise_time_ayrs[, colnames(pwise_time_ayrs) == LU2]
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
  # nrow(dat) <= (choose((nyears + 2 - 1), 2) - nyears) * nplots # temporal dataset : expected number of plots presen
  # ## check dataset end

  dat <- dat[complete.cases(dat), ]

  detach(pwise_time_ayrs)

  gdmTab <- formatsitepair(dat, 4, predData = dat)
  var <- gdm.varImp.synthesis.edit(gdmTab,
    geo = T, includeYear = T, nPerm = 50,
    fullModelOnly = T, spacefortimedataset = "temporal"
  )

  var_plants_bsimt[[i]] <- var
}
end.time <- Sys.time()

names(var_plants_bsimt) <- LU

save(var_plants_bsimt, file = "./2.GDMs/data/OutputData/var_plants_bsimt.RData")

remove(pwise_space, pwise_time_ayrs)