# # # # # # # # # # # # # # #  # # # #
#                                    #
#             FIT GDMS               # 
#   BETA DIVERSITY PLANTS TEMPORAL   # ----
#                                    #
# # # # # # # # # # # # # # #  # # # #

# # # # # # # # # # # # # #
# CONTENT                  ----
# AIM : Fit GDM models of plants beta diversity, for temporal datasets.
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
#set working directory to folder "2.GDM"
#plants
load("./data/InputData/pwise_time_plants.RData")

pwise_time -> pwise_time_plants
# if the uploaded, assembled data files are used, upload those files
# here.


# # # # # # # # # # # # # #
# 2 - TEMPORAL DATASET              ----
#
# create a data frame for the splines
splines_plants_yall_sc<- as.data.frame(matrix(nrow=1200, ncol=16))
colnames(splines_plants_yall_sc)<- c("LUIx", "MOWx", "GRAx", "FERx", 
                                  "COvLUIx","COvMOWx", "COvGRAx", "COvFERx",
                                  "LUIy", "MOWy", "GRAy", "FERy", 
                                  "COvLUIy","COvMOWy", "COvGRAy", "COvFERy")
splines_plants_yall_sc$beta_type<- rep(c("cqn0", "cqn1", "cqn2", "cqn3", "cqn4", "bsim"), each=200)

#residuals of lm(LUI1/2 - EP)
LU1<- c("LUI1res", "MOW1res", "GRA1res", "FER1res")
LU2<- c("LUI2res", "MOW2res", "GRA2res", "FER2res")


# # # # #
# 2.0. -  b0 ----
# 
# b0 : beta diversity with Chao = 0
for(i in 1:4){
  
  attach(pwise_time)
  
  #year and LU as predictors
  dat<- data.frame(cqn0dis, 
                   rep(1, nrow(pwise_time)), 
                   HWG1, RWG1, HWG2, RWG2, 
                   pwise_time[,colnames(pwise_time)== LU1[i]], 
                   YR1, 
                   pwise_time[,colnames(pwise_time)== LU2[i]], 
                   YR2)
  
  colnames(dat)<- c("distance", "weights", 
                    "s1.xCoord", "s1.yCoord", "s2.xCoord", "s2.yCoord", 
                    "s1.LU","s1.Year", "s2.LU", "s2.Year")
  
  dat<- dat[complete.cases(dat),]
  
  detach(pwise_time)
  
  gdmTab<- formatsitepair(dat, 4, predData=dat)
  GDM<- gdm(gdmTab, geo=F)
  
  spline<- isplineExtract(GDM)
  splines_plants_yall_sc[,i]<- spline$x[,1]
  splines_plants_yall_sc[,i+4]<- spline$x[,2]
  
  splines_plants_yall_sc[,i+8]<- spline$y[,1]
  splines_plants_yall_sc[,i+12]<- spline$y[,2]
}

remove(dat, gdmTab, GDM)
  

# # # # #
# 2.1. -  b1 ----
# 
# b1 : beta diversity with Chao = 1
for(i in 1:4){
   
  attach(pwise_time)
  
  #year and LU as predictors
  
  dat<- data.frame(cqn1dis-0.0001, 
                   rep(1, nrow(pwise_time)), 
                   HWG1, HWG1, HWG2, HWG2, 
                   pwise_time[,colnames(pwise_time)== LU1[i]], 
                   YR1, 
                   pwise_time[,colnames(pwise_time)== LU2[i]], 
                   YR2)
  
  colnames(dat)<- c("distance", "weights", 
                    "s1.xCoord", "s1.yCoord", "s2.xCoord", "s2.yCoord", 
                    "s1.LU","s1.Year", "s2.LU", "s2.Year")
  
  dat<- dat[complete.cases(dat),]
  
  detach(pwise_time)
  
  gdmTab<- formatsitepair(dat, 4, predData=dat)
  GDM<- gdm(gdmTab, geo=F)
  
  spline<- isplineExtract(GDM)
  splines_plants_yall_sc[201:400,i]<- spline$x[,1]
  splines_plants_yall_sc[201:400,i+4]<- spline$x[,2]
  
  splines_plants_yall_sc[201:400,i+8]<- spline$y[,1]
  splines_plants_yall_sc[201:400,i+12]<- spline$y[,2]
}

remove(dat, gdmTab, GDM)


# # # # #
# 2.2. -  b2 ----
# 
# b2 : beta diversity with Chao = 2
for(i in 1:4){
  
  attach(pwise_time)
  
  #year and LU as predictors
  dat<- data.frame(cqn2dis, 
                   rep(1, nrow(pwise_time)), 
                   HWG1, HWG1, HWG2, HWG2, 
                   pwise_time[,colnames(pwise_time)== LU1[i]], 
                   YR1, 
                   pwise_time[,colnames(pwise_time)== LU2[i]], 
                   YR2)
  
  colnames(dat)<- c("distance", "weights", 
                    "s1.xCoord", "s1.yCoord", "s2.xCoord", "s2.yCoord", 
                    "s1.LU","s1.Year", "s2.LU", "s2.Year")
  
  dat<- dat[complete.cases(dat),]
  
  detach(pwise_time)
  
  gdmTab<- formatsitepair(dat, 4, predData=dat)
  GDM<- gdm(gdmTab, geo=F)
  
  spline<- isplineExtract(GDM)
  splines_plants_yall_sc[401:600,i]<- spline$x[,1]
  splines_plants_yall_sc[401:600,i+4]<- spline$x[,2]
  
  splines_plants_yall_sc[401:600,i+8]<- spline$y[,1]
  splines_plants_yall_sc[401:600,i+12]<- spline$y[,2]
}

remove(dat,gdmTab, GDM)


# # # # #
# 2.3. -  b3 ----
# 
# b3 : beta diversity with Chao = 3
for(i in 1:4){
  attach(pwise_time)
  
  #year and LU as predictors
  dat<- data.frame(cqn3dis, 
                   rep(1, nrow(pwise_time)), 
                   HWG1, HWG1, HWG2, HWG2, 
                   pwise_time[,colnames(pwise_time)== LU1[i]], 
                   YR1, 
                   pwise_time[,colnames(pwise_time)== LU2[i]], 
                   YR2)
  
  colnames(dat)<- c("distance", "weights", 
                    "s1.xCoord", "s1.yCoord", "s2.xCoord", "s2.yCoord", 
                    "s1.LU","s1.Year", "s2.LU", "s2.Year")
  
  dat<- dat[complete.cases(dat),]
  
  detach(pwise_time)
  
  gdmTab<- formatsitepair(dat, 4, predData=dat)
  GDM<- gdm(gdmTab, geo=F)
  
  spline<- isplineExtract(GDM)
  splines_plants_yall_sc[601:800,i]<- spline$x[,1]
  splines_plants_yall_sc[601:800,i+4]<- spline$x[,2]
  
  splines_plants_yall_sc[601:800,i+8]<- spline$y[,1]
  splines_plants_yall_sc[601:800,i+12]<- spline$y[,2]
}

remove(dat, gdmTab,GDM)


# # # # #
# 2.4. -  b4 ----
# 
# b4 : beta diversity with Chao = 4
for(i in 1:4){
 
  attach(pwise_time)
  
  #year and LU as predictors
  dat<- data.frame(cqn4dis, 
                   rep(1, nrow(pwise_time)), 
                   HWG1, HWG1, HWG2, HWG2, 
                   pwise_time[,colnames(pwise_time)== LU1[i]], 
                   YR1, 
                   pwise_time[,colnames(pwise_time)== LU2[i]], 
                   YR2)
  
  colnames(dat)<- c("distance", "weights", 
                    "s1.xCoord", "s1.yCoord", "s2.xCoord", "s2.yCoord", 
                    "s1.LU","s1.Year", "s2.LU", "s2.Year")
  
  dat<- dat[complete.cases(dat),]
  
  detach(pwise_time)
  
  gdmTab<- formatsitepair(dat, 4, predData=dat)
  GDM<- gdm(gdmTab, geo=F)
  
  spline<- isplineExtract(GDM)
  splines_plants_yall_sc[801:1000,i]<- spline$x[,1]
  splines_plants_yall_sc[801:1000,i+4]<- spline$x[,2]
  
  splines_plants_yall_sc[801:1000,i+8]<- spline$y[,1]
  splines_plants_yall_sc[801:1000,i+12]<- spline$y[,2]
}

remove(dat,gdmTab, GDM)


# # # # #
# 2.5. -  bsim ----
# 
# bsim : beta diversity simpson
for(i in 1:4){
 
  attach(pwise_time)
  
  #year and LUI as predictors
  dat<- data.frame(bsim, 
                   rep(1, nrow(pwise_time)), 
                   HWG1, HWG1, HWG2, HWG2, 
                   pwise_time[,colnames(pwise_time)== LU1[i]], 
                   YR1, 
                   pwise_time[,colnames(pwise_time)== LU2[i]], 
                   YR2)
  
  colnames(dat)<- c("distance", "weights", 
                    "s1.xCoord", "s1.yCoord", "s2.xCoord", "s2.yCoord", 
                    "s1.LU","s1.Year", "s2.LU", "s2.Year")
  
  dat<- dat[complete.cases(dat),]
  
  detach(pwise_time)
  
  gdmTab<- formatsitepair(dat, 4, predData=dat)
  GDM<- gdm(gdmTab, geo=F)
  
  spline<- isplineExtract(GDM)
  splines_plants_yall_sc[1001:1200,i]<- spline$x[,1]
  splines_plants_yall_sc[1001:1200,i+4]<- spline$x[,2]
  
  splines_plants_yall_sc[1001:1200,i+8]<- spline$y[,1]
  splines_plants_yall_sc[1001:1200,i+12]<- spline$y[,2]
}

remove(dat, gdmTab, GDM)

save(splines_plants_yall_sc, file="data/OutputData/splines_plants_yall_sc_new.RData")

#set working directory to folder "4. Plotting results"
save(splines_plants_yall_sc, file="./data/InputData/splines_pred_yall_sc.RData")


remove(pwise_time)