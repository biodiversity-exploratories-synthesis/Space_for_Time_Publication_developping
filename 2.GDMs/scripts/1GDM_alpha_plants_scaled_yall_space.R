# # # # # # # # # # # # # # #  # # # #
#                                    #
#             FIT GDMS               # 
#      ALPHA DIVERSITY PLANTS        # ----
#                                    #
# # # # # # # # # # # # # # #  # # # #
#
# spatial and temporal dataset within this script

# # # # # # # # # # # # # #
# CONTENT                  ----
# AIM : Fit GDM models of plants alpha diversity, for time and space datasets.
# OUTPUT : 
#TODO

# # # # # # # # # # # # # #
# 0 - REQUIREMENTS              ----
#
# packages
# install.packages("gdm")
library(gdm)
library(vegan)
library(tidyverse)
library(ggplot2)
library(patchwork)


# # # # #
# 0.a. - DATA ----
#
# plants
#set working directory to folder "2.GDM"
#plants
load("./data/InputData/pwise_time_plants.RData")
load("./data/InputData/pwise_space_plants.RData")

# if the uploaded, assembled data files are used, upload those plant files
# here.


# Preparing datasets
#
# make sure distances are not <0, >1 --> restrict to a tiny bit >0 and <1 --> DO NOT SAVE this
# vegan::decostand()
#TODO is this only for alphadiversity? not for beta plants, 
pwise_space_plants$da0st<- decostand(pwise_space_plants$da0, method="range", na.rm=T)
pwise_space_plants$da1st<- decostand(pwise_space_plants$da1, method="range", na.rm=T)
pwise_space_plants$da2st<- decostand(pwise_space_plants$da2, method="range", na.rm=T)
pwise_space_plants$da3st<- decostand(pwise_space_plants$da3, method="range", na.rm=T)
pwise_space_plants$da4st<- decostand(pwise_space_plants$da4, method="range", na.rm=T)

pwise_time_plants$da0st<- decostand(pwise_time_plants$da0abs, method="range", na.rm=T)
pwise_time_plants$da1st<- decostand(pwise_time_plants$da1abs, method="range", na.rm=T)
pwise_time_plants$da2st<- decostand(pwise_time_plants$da2abs, method="range", na.rm=T)
pwise_time_plants$da3st<- decostand(pwise_time_plants$da3abs, method="range", na.rm=T)
pwise_time_plants$da4st<- decostand(pwise_time_plants$da4abs, method="range", na.rm=T)


pwise_space<- pwise_space_plants
pwise_time<- pwise_time_plants

# # # # # # # # # # # # # #
# 1 - SPATIAL DATASET              ----
#
# create a data frame for the splines

splines_alpha_plant<- as.data.frame(matrix(nrow=1000, ncol=16))
colnames(splines_alpha_plant)<- c("LUIx", "MOWx", "GRAx", "FERx", 
                          "COvLUIx","COvMOWx", "COvGRAx", "COvFERx",
                          "LUIy", "MOWy", "GRAy", "FERy", 
                          "COvLUIy","COvMOWy", "COvGRAy", "COvFERy")
splines_alpha_plant$alpha_type<- rep(c("a0", "a1", "a2", "a3", "a4"), each=200)



LU1<- c("LUI1res", "MOW1res", "GRA1res", "FER1res")
LU2<- c("LUI2res", "MOW2res", "GRA2res", "FER2res")

#randomized data frames - input for GDMs
#plants
subs<- matrix(ncol=1000, nrow=8250)
set.seed(123)

for (m in 1:1000){
  sub_row<- sample(c(1:122925), 8250, replace =F)#number of rows in space dat object - complete cases (all div)/time data
  subs[,m]<- sub_row}


# # # # #
# 1.0. -  a0 ----
# 
# a0 : alpha diversity with Chao = 0

for(i in 1:4){
  attach(pwise_space)
  
  #LUI and GEO as predictors
  dat<- data.frame(da0st, 
                   rep(1, nrow(pwise_space)), 
                   HWG1, RWG1, HWG2, RWG2, 
                   pwise_space[,colnames(pwise_space)== LU1[i]], 
                   pwise_space[,colnames(pwise_space)== LU2[i]])
  
  
  colnames(dat)<- c("distance", "weights", 
                    "s1.xCoord", "s1.yCoord", "s2.xCoord", "s2.yCoord", 
                    "s1.LU", "s2.LU")
  
  dat<- dat[complete.cases(dat),]
  
    detach(pwise_space)
  
  LUsplx<- matrix(ncol=1000, nrow=200)
  LUsply<- matrix(ncol=1000, nrow=200)
  GEOsplx<- matrix(ncol=1000, nrow=200)
  GEOsply<- matrix(ncol=1000, nrow=200)
  
  
  for (m in 1:1000){
    
    sub<- dat[subs[,m],]
    sub<- sub[complete.cases(sub),]
    gdmTab<- formatsitepair(sub, 4, predData=sub)
    gdm_sub<- gdm(gdmTab, geo=T)
    
    
    spline<- isplineExtract(gdm_sub)
    LUsplx[,m]<- spline$x[,2]
    LUsply[,m]<- spline$y[,2]
    GEOsplx[,m]<- spline$x[,1]
    GEOsply[,m]<- spline$y[,1]
    
  
  }
  
  splines_alpha_plant[1:200,i]<- apply(LUsplx, 1, mean, na.rm=T)
  splines_alpha_plant[1:200,i+4]<- apply(GEOsplx, 1, mean, na.rm=T)
  
  splines_alpha_plant[1:200,i+8]<- apply(LUsply, 1, mean, na.rm=T)
  splines_alpha_plant[1:200,i+12]<- apply(GEOsply, 1, mean, na.rm=T)
}


# # # # #
# 1.1. -  a1 ----
# 
# a1 : alpha diversity with Chao = 1

for(i in 1:4){
  
  attach(pwise_space)
  
  #LUI and GEO as predictors
  dat<- data.frame(da1st, 
                   rep(1, nrow(pwise_space)), 
                   HWG1, RWG1, HWG2, RWG2, 
                   pwise_space[,colnames(pwise_space)== LU1[i]], 
                   pwise_space[,colnames(pwise_space)== LU2[i]])
  
  
  colnames(dat)<- c("distance", "weights", 
                    "s1.xCoord", "s1.yCoord", "s2.xCoord", "s2.yCoord", 
                    "s1.LU", "s2.LU")
  
  dat<- dat[complete.cases(dat),]
  
  detach(pwise_space)
  
  LUsplx<- matrix(ncol=1000, nrow=200)
  LUsply<- matrix(ncol=1000, nrow=200)
  GEOsplx<- matrix(ncol=1000, nrow=200)
  GEOsply<- matrix(ncol=1000, nrow=200)
  
  
  for (m in 1:1000){
    
    sub<- dat[subs[,m],]
    sub<- sub[complete.cases(sub),]
    gdmTab<- formatsitepair(sub, 4, predData=sub)
    gdm_sub<- gdm(gdmTab, geo=T)
    
    spline<- isplineExtract(gdm_sub)
    LUsplx[,m]<- spline$x[,2]
    LUsply[,m]<- spline$y[,2]
    GEOsplx[,m]<- spline$x[,1]
    GEOsply[,m]<- spline$y[,1]
    
  }
  
  splines_alpha_plant[201:400,i]<- apply(LUsplx, 1, mean, na.rm=T)
  splines_alpha_plant[201:400,i+4]<- apply(GEOsplx, 1, mean, na.rm=T)
  
  splines_alpha_plant[201:400,i+8]<- apply(LUsply, 1, mean, na.rm=T)
  splines_alpha_plant[201:400,i+12]<- apply(GEOsply, 1, mean, na.rm=T)
}


# # # # #
# 1.2. -  a2 ----
# 
# a2 : alpha diversity with Chao = 2

for(i in 1:4){
  attach(pwise_space)
  
  #LUI and GEO as predictors
  dat<- data.frame(da2st, 
                   rep(1, nrow(pwise_space)), 
                   HWG1, RWG1, HWG2, RWG2, 
                   pwise_space[,colnames(pwise_space)== LU1[i]], 
                   pwise_space[,colnames(pwise_space)== LU2[i]])
  
  
  colnames(dat)<- c("distance", "weights", 
                    "s1.xCoord", "s1.yCoord", "s2.xCoord", "s2.yCoord", 
                    "s1.LU", "s2.LU")
  
  dat<- dat[complete.cases(dat),]
  
  detach(pwise_space)
  
  LUsplx<- matrix(ncol=1000, nrow=200)
  LUsply<- matrix(ncol=1000, nrow=200)
  GEOsplx<- matrix(ncol=1000, nrow=200)
  GEOsply<- matrix(ncol=1000, nrow=200)
  
  
  for (m in 1:1000){
    
    sub<- dat[subs[,m],]
    sub<- sub[complete.cases(sub),]
    gdmTab<- formatsitepair(sub, 4, predData=sub)
    gdm_sub<- gdm(gdmTab, geo=T)
    
    spline<- isplineExtract(gdm_sub)
    LUsplx[,m]<- spline$x[,2]
    LUsply[,m]<- spline$y[,2]
    GEOsplx[,m]<- spline$x[,1]
    GEOsply[,m]<- spline$y[,1]
    
  }
  
  splines_alpha_plant[401:600,i]<- apply(LUsplx, 1, mean, na.rm=T)
  splines_alpha_plant[401:600,i+4]<- apply(GEOsplx, 1, mean, na.rm=T)
  
  splines_alpha_plant[401:600,i+8]<- apply(LUsply, 1, mean, na.rm=T)
  splines_alpha_plant[401:600,i+12]<- apply(GEOsply, 1, mean, na.rm=T)
}


# # # # #
# 1.3. -  a3 ----
# 
# a3 : alpha diversity with Chao = 3

for(i in 1:4){
  
  attach(pwise_space)
  
  #LUI and GEO as predictors
  dat<- data.frame(da3st, 
                   rep(1, nrow(pwise_space)), 
                   HWG1, RWG1, HWG2, RWG2, 
                   pwise_space[,colnames(pwise_space)== LU1[i]], 
                   pwise_space[,colnames(pwise_space)== LU2[i]])
  
  
  colnames(dat)<- c("distance", "weights", 
                    "s1.xCoord", "s1.yCoord", "s2.xCoord", "s2.yCoord", 
                    "s1.LU", "s2.LU")
  
  dat<- dat[complete.cases(dat),]
  
  detach(pwise_space)
  
  LUsplx<- matrix(ncol=1000, nrow=200)
  LUsply<- matrix(ncol=1000, nrow=200)
  GEOsplx<- matrix(ncol=1000, nrow=200)
  GEOsply<- matrix(ncol=1000, nrow=200)
  
  
  for (m in 1:1000){
    
    sub<- dat[subs[,m],]
    sub<- sub[complete.cases(sub),]
    gdmTab<- formatsitepair(sub, 4, predData=sub)
    gdm_sub<- gdm(gdmTab, geo=T)
    
    spline<- isplineExtract(gdm_sub)
    LUsplx[,m]<- spline$x[,2]
    LUsply[,m]<- spline$y[,2]
    GEOsplx[,m]<- spline$x[,1]
    GEOsply[,m]<- spline$y[,1]
    
  }
  
  splines_alpha_plant[601:800,i]<- apply(LUsplx, 1, mean, na.rm=T)
  splines_alpha_plant[601:800,i+4]<- apply(GEOsplx, 1, mean, na.rm=T)
  
  splines_alpha_plant[601:800,i+8]<- apply(LUsply, 1, mean, na.rm=T)
  splines_alpha_plant[601:800,i+12]<- apply(GEOsply, 1, mean, na.rm=T)
}



# # # # #
# 1.4. -  a4 ----
# 
# a4 : alpha diversity with Chao = 4

for(i in 1:4){
  
  attach(pwise_space)
  
  #LUI and GEO as predictors
  dat<- data.frame(da4st, 
                   rep(1, nrow(pwise_space)), 
                   HWG1, RWG1, HWG2, RWG2, 
                   pwise_space[,colnames(pwise_space)== LU1[i]], 
                   pwise_space[,colnames(pwise_space)== LU2[i]])
  
  
  colnames(dat)<- c("distance", "weights", 
                    "s1.xCoord", "s1.yCoord", "s2.xCoord", "s2.yCoord", 
                    "s1.LU", "s2.LU")
  
  dat<- dat[complete.cases(dat),]
  
  detach(pwise_space)
  
  LUsplx<- matrix(ncol=1000, nrow=200)
  LUsply<- matrix(ncol=1000, nrow=200)
  GEOsplx<- matrix(ncol=1000, nrow=200)
  GEOsply<- matrix(ncol=1000, nrow=200)
  
  
  for (m in 1:1000){
    
    sub<- dat[subs[,m],]
    sub<- sub[complete.cases(sub),]
    gdmTab<- formatsitepair(sub, 4, predData=sub)
    gdm_sub<- gdm(gdmTab, geo=T)
    
    spline<- isplineExtract(gdm_sub)
    LUsplx[,m]<- spline$x[,2]
    LUsply[,m]<- spline$y[,2]
    GEOsplx[,m]<- spline$x[,1]
    GEOsply[,m]<- spline$y[,1]
    
  }
  
  splines_alpha_plant[801:1000,i]<- apply(LUsplx, 1, mean, na.rm=T)
  splines_alpha_plant[801:1000,i+4]<- apply(GEOsplx, 1, mean, na.rm=T)
  
  splines_alpha_plant[801:1000,i+8]<- apply(LUsply, 1, mean, na.rm=T)
  splines_alpha_plant[801:1000,i+12]<- apply(GEOsply, 1, mean, na.rm=T)
}

splines_alpha_plant_space_sc<- splines_alpha_plant

save(splines_alpha_plant_space_sc, file="./data/OutputData/splines_alpha_plant_space_sc.RData")


# # # # # # # # # # # # # #
# 2 - TEMPORAL DATASET              ----
#
# create a data frame for the splines
#
splines_alpha_plant<- as.data.frame(matrix(nrow=1000, ncol=16))
colnames(splines_alpha_plant)<- c("LUIx", "MOWx", "GRAx", "FERx", 
                                  "COvLUIx","COvMOWx", "COvGRAx", "COvFERx",
                                  "LUIy", "MOWy", "GRAy", "FERy", 
                                  "COvLUIy","COvMOWy", "COvGRAy", "COvFERy")
splines_alpha_plant$alpha_type<- rep(c("a0", "a1", "a2", "a3", "a4"), each=200)



LU1<- c("LUI1res", "MOW1res", "GRA1res", "FER1res")
LU2<- c("LUI2res", "MOW2res", "GRA2res", "FER2res")


# # # # #
# 2.0. -  a0 ----
# 
# a0 : alpha diversity with Chao = 0
# plants

for(i in 1:4){

  #mLU<- mwLUI[i]
  attach(pwise_time)
  
  #LUI and GEO as predictors
  dat<- data.frame(da0st, 
                   rep(1, nrow(pwise_time)), 
                   HWG1, RWG1, HWG2, RWG2, 
                   (pwise_time[,colnames(pwise_time)== LU1[i]]), 
                   YR1,
                   (pwise_time[,colnames(pwise_time)== LU2[i]]),
                   YR2) 
  
  
  colnames(dat)<- c("distance", "weights", 
                    "s1.xCoord", "s1.yCoord", "s2.xCoord", "s2.yCoord", 
                    "s1.LU","s1.Year", "s2.LU", "s2.Year")
  
  dat<- dat[complete.cases(dat),]
  
  detach(pwise_time)
  
  gdmTab<- formatsitepair(dat, 4, predData=dat)
  gdm_sub<- gdm(gdmTab, geo=F)
    
    
  spline<- isplineExtract(gdm_sub)
  splines_alpha_plant[1:200,i]<- spline$x[,1]#LUI
  splines_alpha_plant[1:200,i+4]<- spline$x[,2]#year
  
  splines_alpha_plant[1:200,i+8]<- spline$y[,1]
  splines_alpha_plant[1:200,i+12]<- spline$y[,2]

}

remove(dat, gdmTab)

# # # # #
# 2.1. -  a1 ----
# 
# a1 : alpha diversity with Chao = 1

for(i in 1:4){
  #mLU<- mwLUI[i]
  attach(pwise_time)
  
  #LUI and GEO as predictors
  dat<- data.frame(da1st, 
                   rep(1, nrow(pwise_time)), 
                   HWG1, RWG1, HWG2, RWG2, 
                   (pwise_time[,colnames(pwise_time)== LU1[i]]), 
                   YR1,
                   (pwise_time[,colnames(pwise_time)== LU2[i]]),
                   YR2) 
  
  colnames(dat)<- c("distance", "weights", 
                    "s1.xCoord", "s1.yCoord", "s2.xCoord", "s2.yCoord", 
                    "s1.LU","s1.Year", "s2.LU", "s2.Year")
  
  dat<- dat[complete.cases(dat),]
  
  detach(pwise_time)
  
  gdmTab<- formatsitepair(dat, 4, predData=dat)
  gdm_sub<- gdm(gdmTab, geo=F)
  
  spline<- isplineExtract(gdm_sub)
  splines_alpha_plant[201:400,i]<- spline$x[,1]
  splines_alpha_plant[201:400,i+4]<- spline$x[,2]
  
  splines_alpha_plant[201:400,i+8]<- spline$y[,1]
  splines_alpha_plant[201:400,i+12]<- spline$y[,2]
  
}

remove(dat, gdmTab)

# # # # #
# 2.2. -  a2 ----
# 
# a2 : alpha diversity with Chao = 2

for(i in 1:4){
  #mLU<- mwLUI[i]
  attach(pwise_time)
  
  #LUI and GEO as predictors
  dat<- data.frame(da2st, 
                   rep(1, nrow(pwise_time)), 
                   HWG1, RWG1, HWG2, RWG2, 
                   (pwise_time[,colnames(pwise_time)== LU1[i]]), 
                   YR1,
                   (pwise_time[,colnames(pwise_time)== LU2[i]]),
                   YR2) 
  
  
  colnames(dat)<- c("distance", "weights", 
                    "s1.xCoord", "s1.yCoord", "s2.xCoord", "s2.yCoord", 
                    "s1.LU","s1.Year", "s2.LU", "s2.Year")
  
  dat<- dat[complete.cases(dat),]
  
  detach(pwise_time)
  

  gdmTab<- formatsitepair(dat, 4, predData=dat)
  gdm_sub<- gdm(gdmTab, geo=F)
  
  
  spline<- isplineExtract(gdm_sub)
  splines_alpha_plant[401:600,i]<- spline$x[,1]
  splines_alpha_plant[401:600,i+4]<- spline$x[,2]
  
  splines_alpha_plant[401:600,i+8]<- spline$y[,1]
  splines_alpha_plant[401:600,i+12]<- spline$y[,2]
  
}

remove(dat, gdmTab)


# # # # #
# 2.3. -  a3 ----
# 
# a3 : alpha diversity with Chao = 3

for(i in 1:4){
  #mLU<- mwLUI[i]
  attach(pwise_time)
  
  #LUI and GEO as predictors
  dat<- data.frame(da3st, 
                   rep(1, nrow(pwise_time)), 
                   HWG1, RWG1, HWG2, RWG2, 
                   (pwise_time[,colnames(pwise_time)== LU1[i]]), 
                   YR1,
                   (pwise_time[,colnames(pwise_time)== LU2[i]]),
                   YR2) 
  
  
  colnames(dat)<- c("distance", "weights", 
                    "s1.xCoord", "s1.yCoord", "s2.xCoord", "s2.yCoord", 
                    "s1.LU","s1.Year", "s2.LU", "s2.Year")
  
  dat<- dat[complete.cases(dat),]
  
  detach(pwise_time)

  
  gdmTab<- formatsitepair(dat, 4, predData=dat)
  gdm_sub<- gdm(gdmTab, geo=F)
  
  spline<- isplineExtract(gdm_sub)
  splines_alpha_plant[601:800,i]<- spline$x[,1]
  splines_alpha_plant[601:800,i+4]<- spline$x[,2]
  
  splines_alpha_plant[601:800,i+8]<- spline$y[,1]
  splines_alpha_plant[601:800,i+12]<- spline$y[,2]
  
}

remove(dat, gdmTab)

# # # # #
# 2.4. -  a4 ----
# 
# a4 : alpha diversity with Chao = 4

for(i in 1:4){

  #mLU<- mwLUI[i]
  attach(pwise_time)
  
  #LUI and GEO as predictors
  dat<- data.frame(da4st, 
                   rep(1, nrow(pwise_time)), 
                   HWG1, RWG1, HWG2, RWG2, 
                   (pwise_time[,colnames(pwise_time)== LU1[i]]), 
                   YR1,
                   (pwise_time[,colnames(pwise_time)== LU2[i]]),
                   YR2) 
  
  colnames(dat)<- c("distance", "weights", 
                    "s1.xCoord", "s1.yCoord", "s2.xCoord", "s2.yCoord", 
                    "s1.LU","s1.Year", "s2.LU", "s2.Year")
  
  dat<- dat[complete.cases(dat),]
  
  detach(pwise_time)
  
  gdmTab<- formatsitepair(dat, 4, predData=dat)
  gdm_sub<- gdm(gdmTab, geo=F)
  
  spline<- isplineExtract(gdm_sub)
  splines_alpha_plant[801:1000,i]<- spline$x[,1]
  splines_alpha_plant[801:1000,i+4]<- spline$x[,2]
  
  splines_alpha_plant[801:1000,i+8]<- spline$y[,1]
  splines_alpha_plant[801:1000,i+12]<- spline$y[,2]
  
}

remove(dat, gdmTab)

splines_alpha_plant_time_sc<- splines_alpha_plant

save(splines_alpha_plant_time_sc, file="./data/OutputData/splines_alpha_plant_time_sc.RData")

#set working directory to folder "4. Plotting results"
save(splines_alpha_plant_sc_time, file="./data/InputData/splines_alpha_plant_sc_time.RData")

remove(pwise_time, pwise_space)
