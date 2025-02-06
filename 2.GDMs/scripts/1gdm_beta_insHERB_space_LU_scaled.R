# # # # # # # # # # # # # # #  # # # #
#                                    #
#             FIT GDMS               # 
#  BETA DIVERSITY HERBIVORES SPACE   # ----
#                                    #
# # # # # # # # # # # # # # #  # # # #

# # # # # # # # # # # # # #
# CONTENT                  ----
# AIM : Fit GDM models of herbivores beta diversity, for space datasets
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
load("./data/InputData/pwise_space_herb.RData")

pwise_space -> pwise_space_herb
# if the uploaded, assembled data files are used, upload the complete insect files
# here.


# Preparing datasets
#
# make sure distances are not <0, >1 --> restrict to a tiny bit >0 and <1 --> DO NOT SAVE this
pwise_space$hcqn0dis[pwise_space$hcqn0dis<0.001]<- 0.001
pwise_space$hcqn0dis[pwise_space$hcqn0dis>0.999]<- 0.999

pwise_space$hcqn1dis[pwise_space$hcqn1dis<0.001]<- 0.001
pwise_space$hcqn1dis[pwise_space$hcqn1dis>0.999]<- 0.999

pwise_space$hcqn2dis[pwise_space$hcqn2dis<0.001]<- 0.001
pwise_space$hcqn2dis[pwise_space$hcqn2dis>0.999]<- 0.999

pwise_space$hcqn3dis[pwise_space$hcqn3dis<0.001]<- 0.001
pwise_space$hcqn3dis[pwise_space$hcqn3dis>0.999]<- 0.999

pwise_space$hcqn4dis[pwise_space$hcqn4dis<0.001]<- 0.001
pwise_space$hcqn4dis[pwise_space$hcqn4dis>0.999]<- 0.999

pwise_space$hbsim[pwise_space$hbsim<0.001]<- 0.001
pwise_space$hbsim[pwise_space$hbsim>0.999]<- 0.999
#TODO is this in the methods?


# # # # # # # # # # # # # #
# 1 - SPATIAL DATASET              ----
#
# create a data frame for the splines
splines_her_space_sc<- as.data.frame(matrix(nrow=1200, ncol=16))
colnames(splines_her_space_sc)<- c("LUIx", "MOWx", "GRAx", "FERx", 
                                "COvLUIx","COvMOWx", "COvGRAx", "COvFERx",
                                "LUIy", "MOWy", "GRAy", "FERy", 
                                "COvLUIy","COvMOWy", "COvGRAy", "COvFERy")
splines_her_space_sc$beta_type<- rep(c("hcqn0", "hcqn1", "hcqn2", "hcqn3", "hcqn4", "hbsim"), each=200)

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
  dat<- data.frame(hcqn0dis, 
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
  
  splines_her_space_sc[,i]<- apply(LUsplx, 1, mean, na.rm=T)
  splines_her_space_sc[,i+4]<- apply(GEOsplx, 1, mean, na.rm=T)
  
  splines_her_space_sc[,i+8]<- apply(LUsply, 1, mean, na.rm=T)
  splines_her_space_sc[,i+12]<- apply(GEOsply, 1, mean, na.rm=T)
}



# # # # #
# 1.1. -  b1 ----
# 
# b1 : beta diversity with Chao = 1
for(i in 1:4){
  attach(pwise_space)
  
  #LUI and GEO as predictors
  dat<- data.frame(hcqn1dis, 
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
  
  splines_her_space_sc[201:400,i]<- apply(LUsplx, 1, mean, na.rm=T)
  splines_her_space_sc[201:400,i+4]<- apply(GEOsplx, 1, mean, na.rm=T)
  
  splines_her_space_sc[201:400,i+8]<- apply(LUsply, 1, mean, na.rm=T)
  splines_her_space_sc[201:400,i+12]<- apply(GEOsply, 1, mean, na.rm=T)
}



# # # # #
# 1.2. -  b2 ----
# 
# b2 : beta diversity with Chao = 2
for(i in 1:4){
  
  attach(pwise_space)
  
  #LUI and GEO as predictors
  dat<- data.frame(hcqn2dis, 
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
  
  splines_her_space_sc[401:600,i]<- apply(LUsplx, 1, mean, na.rm=T)
  splines_her_space_sc[401:600,i+4]<- apply(GEOsplx, 1, mean, na.rm=T)
  
  splines_her_space_sc[401:600,i+8]<- apply(LUsply, 1, mean, na.rm=T)
  splines_her_space_sc[401:600,i+12]<- apply(GEOsply, 1, mean, na.rm=T)
}



# # # # #
# 1.3. -  b3 ----
# 
# b3 : beta diversity with Chao = 3
for(i in 1:4){
  
  attach(pwise_space)
  
  #LUI and GEO as predictors
  dat<- data.frame(hcqn3dis, 
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
  
  
  splines_her_space_sc[601:800,i]<- apply(LUsplx, 1, mean, na.rm=T)
  splines_her_space_sc[601:800,i+4]<- apply(GEOsplx, 1, mean, na.rm=T)
  
  splines_her_space_sc[601:800,i+8]<- apply(LUsply, 1, mean, na.rm=T)
  splines_her_space_sc[601:800,i+12]<- apply(GEOsply, 1, mean, na.rm=T)
}


# # # # #
# 1.4. -  b4 ----
# 
# b4 : beta diversity with Chao = 4
for(i in 1:4){
  
  attach(pwise_space)
  
  #LUI and GEO as predictors
  dat<- data.frame(hcqn4dis, 
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
  
  

  splines_her_space_sc[801:1000,i]<- apply(LUsplx, 1, mean, na.rm=T)
  splines_her_space_sc[801:1000,i+4]<- apply(GEOsplx, 1, mean, na.rm=T)
  
  splines_her_space_sc[801:1000,i+8]<- apply(LUsply, 1, mean, na.rm=T)
  splines_her_space_sc[801:1000,i+12]<- apply(GEOsply, 1, mean, na.rm=T)
}

# # # # #
# 1.5. -  bsim ----
# 
# beta sim for herbovires : hbsim

for(i in 1:4){
 
  attach(pwise_space)
  
  #LUI and GEO as predictors
  dat<- data.frame(hbsim, 
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
  
  
  splines_her_space_sc[1001:1200,i]<- apply(LUsplx, 1, mean, na.rm=T)
  splines_her_space_sc[1001:1200,i+4]<- apply(GEOsplx, 1, mean, na.rm=T)
  
  splines_her_space_sc[1001:1200,i+8]<- apply(LUsply, 1, mean, na.rm=T)
  splines_her_space_sc[1001:1200,i+12]<- apply(GEOsply, 1, mean, na.rm=T)
}


save(splines_her_space_sc, file="./data/OutputData/splines_her_space_sc_new.RData")

#set working directory to folder "4. Plotting results"
save(splines_herb_space_sc, file="./data/InputData/splines_herb_space_sc.RData")

remove(pwise_space)
remove(pwise_space)
