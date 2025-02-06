# # # # # # # # # # # # # # #  # # # #
#                                    #
#             FIT GDMS               # 
#  BETA DIVERSITY PREDATORS SPACE    # ----
#                                    #
# # # # # # # # # # # # # # #  # # # #

# # # # # # # # # # # # # #
# CONTENT                  ----
# AIM : Fit GDM models of predators beta diversity, for space datasets
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
#secondary consumers/predators
load("./data/InputData/pwise_space_pred.RData")

pwise_space -> pwise_space_pred
# if the uploaded, assembled data files are used, upload the complete insect files
# here.


# Preparing datasets
#

#make sure distances are not <0, >1 --> restrict to a tiny bit >0 and <1 --> DO NOT SAVE this
pwise_space$pcqn0dis[pwise_space$pcqn0dis<0.001]<- 0.001
pwise_space$pcqn0dis[pwise_space$pcqn0dis>0.999]<- 0.999

pwise_space$pcqn1dis[pwise_space$pcqn1dis<0.001]<- 0.001
pwise_space$pcqn1dis[pwise_space$pcqn1dis>0.999]<- 0.999

pwise_space$pcqn2dis[pwise_space$pcqn2dis<0.001]<- 0.001
pwise_space$pcqn2dis[pwise_space$pcqn2dis>0.999]<- 0.999

pwise_space$pcqn3dis[pwise_space$pcqn3dis<0.001]<- 0.001
pwise_space$pcqn3dis[pwise_space$pcqn3dis>0.999]<- 0.999

pwise_space$pcqn4dis[pwise_space$pcqn4dis<0.001]<- 0.001
pwise_space$pcqn4dis[pwise_space$pcqn4dis>0.999]<- 0.999



# # # # # # # # # # # # # #
# 1 - SPATIAL DATASET              ----
#
# create a data frame for the splines
splines_pred_space_sc<- as.data.frame(matrix(nrow=1200, ncol=16))
colnames(splines_pred_space_sc)<- c("LUIx", "MOWx", "GRAx", "FERx", 
                                "COvLUIx","COvMOWx", "COvGRAx", "COvFERx",
                                "LUIy", "MOWy", "GRAy", "FERy", 
                                "COvLUIy","COvMOWy", "COvGRAy", "COvFERy")
splines_pred_space_sc$beta_type<- rep(c("pcqn0", "pcqn1", "pcqn2", "pcqn3", "pcqn4", "pbsim"), each=200)

LU1<- c("LUI1res", "MOW1res", "GRA1res", "FER1res")
LU2<- c("LUI2res", "MOW2res", "GRA2res", "FER2res")

LU<- c("LUI", "MOW", "GRA", "FER")

#randomized data frames - input for GDMs

subs<- matrix(ncol=1000, nrow=5943)
set.seed(123)
for (m in 1:1000){
  sub_row<- sample(c(1:98694), 5943, replace =F)#number of rows in dat object - complete cases (all div)
  subs[,m]<- sub_row}


# # # # #
# 1.0. -  b0 ----
# 
# b0 : beta diversity with Chao = 0
for(i in 1:4){
   
  attach(pwise_space)
  
  #LUI and GEO as predictors
  dat<- data.frame(pcqn0dis, 
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
    gdmTab<- formatsitepair(sub, 4, predData=sub)
    gdm_sub<- gdm(gdmTab, geo=T)
    
    spline<- isplineExtract(gdm_sub)
    LUsplx[,i]<- spline$x[,2]
    LUsply[,i]<- spline$y[,2]
    GEOsplx[,i]<- spline$x[,1]
    GEOsply[,i]<- spline$y[,1]
    

  }
  
  splines_pred_space_sc[,i]<- apply(LUsplx, 1, mean, na.rm=T)
  splines_pred_space_sc[,i+4]<- apply(GEOsplx, 1, mean, na.rm=T)
  
  splines_pred_space_sc[,i+8]<- apply(LUsply, 1, mean, na.rm=T)
  splines_pred_space_sc[,i+12]<- apply(GEOsply, 1, mean, na.rm=T)
}



# # # # #
# 1.1. -  b1 ----
# 
# b1 : beta diversity with Chao = 1
for(i in 1:4){
  attach(pwise_space)
  
  #LUI and GEO as predictors
  dat<- data.frame(pcqn1dis, 
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
    gdmTab<- formatsitepair(sub, 4, predData=sub)
    gdm_sub<- gdm(gdmTab, geo=T)
    
    
    spline<- isplineExtract(gdm_sub)
    LUsplx[,i]<- spline$x[,2]
    LUsply[,i]<- spline$y[,2]
    GEOsplx[,i]<- spline$x[,1]
    GEOsply[,i]<- spline$y[,1]
    
    
    
  }
  
  splines_pred_space_sc[201:400,i]<- apply(LUsplx, 1, mean, na.rm=T)
  splines_pred_space_sc[201:400,i+4]<- apply(GEOsplx, 1, mean, na.rm=T)
  
  splines_pred_space_sc[201:400,i+8]<- apply(LUsply, 1, mean, na.rm=T)
  splines_pred_space_sc[201:400,i+12]<- apply(GEOsply, 1, mean, na.rm=T)
}



# # # # #
# 1.2. -  b2 ----
# 
# b2 : beta diversity with Chao = 2
for(i in 1:4){
  
  attach(pwise_space)
  
  #LUI and GEO as predictors
  dat<- data.frame(pcqn2dis, 
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
    gdmTab<- formatsitepair(sub, 4, predData=sub)
    gdm_sub<- gdm(gdmTab, geo=T)
    
    spline<- isplineExtract(gdm_sub)
    LUsplx[,i]<- spline$x[,2]
    LUsply[,i]<- spline$y[,2]
    GEOsplx[,i]<- spline$x[,1]
    GEOsply[,i]<- spline$y[,1]
    
  }
  
  splines_pred_space_sc[401:600,i]<- apply(LUsplx, 1, mean, na.rm=T)
  splines_pred_space_sc[401:600,i+4]<- apply(GEOsplx, 1, mean, na.rm=T)
  
  splines_pred_space_sc[401:600,i+8]<- apply(LUsply, 1, mean, na.rm=T)
  splines_pred_space_sc[401:600,i+12]<- apply(GEOsply, 1, mean, na.rm=T)
}



# # # # #
# 1.3. -  b3 ----
# 
# b3 : beta diversity with Chao = 3
for(i in 1:4){
  
  attach(pwise_space)
  
  #LUI and GEO as predictors
  dat<- data.frame(pcqn3dis, 
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
    gdmTab<- formatsitepair(sub, 4, predData=sub)
    gdm_sub<- gdm(gdmTab, geo=T)
    
    spline<- isplineExtract(gdm_sub)
    LUsplx[,i]<- spline$x[,2]
    LUsply[,i]<- spline$y[,2]
    GEOsplx[,i]<- spline$x[,1]
    GEOsply[,i]<- spline$y[,1]
    
  
  }
  
  
  splines_pred_space_sc[601:800,i]<- apply(LUsplx, 1, mean, na.rm=T)
  splines_pred_space_sc[601:800,i+4]<- apply(GEOsplx, 1, mean, na.rm=T)
  
  splines_pred_space_sc[601:800,i+8]<- apply(LUsply, 1, mean, na.rm=T)
  splines_pred_space_sc[601:800,i+12]<- apply(GEOsply, 1, mean, na.rm=T)
}


# # # # #
# 1.4. -  b4 ----
# 
# b4 : beta diversity with Chao = 4
for(i in 1:4){
  
  attach(pwise_space)
  
  #LUI and GEO as predictors
  dat<- data.frame(pcqn4dis, 
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
    gdmTab<- formatsitepair(sub, 4, predData=sub)
    gdm_sub<- gdm(gdmTab, geo=T)
    
    spline<- isplineExtract(gdm_sub)
    LUsplx[,i]<- spline$x[,2]
    LUsply[,i]<- spline$y[,2]
    GEOsplx[,i]<- spline$x[,1]
    GEOsply[,i]<- spline$y[,1]
    
  }
  
  

  splines_pred_space_sc[801:1000,i]<- apply(LUsplx, 1, mean, na.rm=T)
  splines_pred_space_sc[801:1000,i+4]<- apply(GEOsplx, 1, mean, na.rm=T)
  
  splines_pred_space_sc[801:1000,i+8]<- apply(LUsply, 1, mean, na.rm=T)
  splines_pred_space_sc[801:1000,i+12]<- apply(GEOsply, 1, mean, na.rm=T)
}

# # # # #
# 1.5. -  bsim ----
# 
# beta sim for predators : pbsim

for(i in 1:4){
 
  attach(pwise_space)
  
  #LUI and GEO as predictors
  dat<- data.frame(pbsim, 
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
    LUsplx[,i]<- spline$x[,2]
    LUsply[,i]<- spline$y[,2]
    GEOsplx[,i]<- spline$x[,1]
    GEOsply[,i]<- spline$y[,1]
    
    
  }
  
  
  splines_pred_space_sc[1001:1200,i]<- apply(LUsplx, 1, mean, na.rm=T)
  splines_pred_space_sc[1001:1200,i+4]<- apply(GEOsplx, 1, mean, na.rm=T)
  
  splines_pred_space_sc[1001:1200,i+8]<- apply(LUsply, 1, mean, na.rm=T)
  splines_pred_space_sc[1001:1200,i+12]<- apply(GEOsply, 1, mean, na.rm=T)
}


save(splines_pred_space_sc, file="./data/Output/splines_pred_space_sc_new.RData")

#set working directory to folder "4. Plotting results"
save(splines_pred_space_sc, file="./data/InputData/splines_pred_space_sc.RData")

remove(pwise_space)
