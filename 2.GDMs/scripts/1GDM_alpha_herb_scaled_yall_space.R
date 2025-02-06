# # # # # # # # # # # # # # #  # # # #
#                                    #
#             FIT GDMS               # 
#   ALPHA DIVERSITY HERBIVORES       # ----
#                                    #
# # # # # # # # # # # # # # #  # # # #
#
# spatial and temporal dataset within this script

# # # # # # # # # # # # # #
# CONTENT                  ----
# AIM : Fit GDM models of herbivores alpha diversity, for time and space datasets.
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
#herbivores
load("./data/InputData/pwise_time_herb.RData")
load("./data/InputData/pwise_space_herb.RData")

# if the uploaded, assembled data files are used, upload the complete insect files
# here.

# # # # # # # # # # # # # #
# 1 - SPATIAL DATASET              ----
#
# create a data frame for the splines
splines_alpha_herb_sc<- as.data.frame(matrix(nrow=1000, ncol=16))
colnames(splines_alpha_herb_sc)<- c("LUIx", "MOWx", "GRAx", "FERx", 
                                  "COvLUIx","COvMOWx", "COvGRAx", "COvFERx",
                                  "LUIy", "MOWy", "GRAy", "FERy", 
                                  "COvLUIy","COvMOWy", "COvGRAy", "COvFERy")
splines_alpha_herb_sc$alpha_type<- rep(c("a0", "a1", "a2", "a3", "a4"), each=200)

LU1<- c("LUI1res", "MOW1res", "GRA1res", "FER1res")
LU2<- c("LUI2res", "MOW2res", "GRA2res", "FER2res")

#LU<- c("LUI", "MOW", "GRA", "FER")
#TODO delete above?

# randomized data frames - input for GDMs
# herbs
subs<- matrix(ncol=1000, nrow=5943)
set.seed(123)

for (m in 1:1000){
  # number of rows in space dat object - complete cases (all div)/time data
  sub_row<- sample(c(1:98946), 5943, replace =F)
  subs[,m]<- sub_row}




# # # # #
# 1.0. -  a0 ----
# 
# a0 : alpha diversity with Chao = 0

for(i in 1:4){
  
  attach(pwise_space)
  
  #LUI and GEO as predictors
  dat<- data.frame(dha0st, 
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
    
    if(sum(gdm_sub$coefficients)==0){
      LUsplx[1:200,m]<- 0
      LUsply[1:200,m]<- 0
      GEOsplx[1:200,m]<- 0
      GEOsply[1:200,m]<- 0
      
    }else{
      spline<- isplineExtract(gdm_sub)
      
      LUsplx[,m]<- spline$x[,2]
      LUsply[,m]<- spline$y[,2]
      GEOsplx[,m]<- spline$x[,1]
      GEOsply[,m]<- spline$y[,1] }
    
    
    
  }
  
  splines_alpha_herb_sc[1:200,i]<- apply(LUsplx, 1, mean, na.rm=T)
  splines_alpha_herb_sc[1:200,i+4]<- apply(GEOsplx, 1, mean, na.rm=T)
  
  splines_alpha_herb_sc[1:200,i+8]<- apply(LUsply, 1, mean, na.rm=T)
  splines_alpha_herb_sc[1:200,i+12]<- apply(GEOsply, 1, mean, na.rm=T)
}


# # # # #
# 1.1. -  a1 ----
# 
# a1 : alpha diversity with Chao = 1

for(i in 1:4){
  
  attach(pwise_space)
  
  #LUI and GEO as predictors
  dat<- data.frame(dha1st, 
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
    
    if(sum(gdm_sub$coefficients)==0){
      LUsplx[1:200,m]<- 0
      LUsply[1:200,m]<- 0
      GEOsplx[1:200,m]<- 0
      GEOsply[1:200,m]<- 0
      
    }else{
      spline<- isplineExtract(gdm_sub)
      
      LUsplx[,m]<- spline$x[,2]
      LUsply[,m]<- spline$y[,2]
      GEOsplx[,m]<- spline$x[,1]
      GEOsply[,m]<- spline$y[,1] }
    
  }
  
  splines_alpha_herb_sc[201:400,i]<- apply(LUsplx, 1, mean, na.rm=T)
  splines_alpha_herb_sc[201:400,i+4]<- apply(GEOsplx, 1, mean, na.rm=T)
  
  splines_alpha_herb_sc[201:400,i+8]<- apply(LUsply, 1, mean, na.rm=T)
  splines_alpha_herb_sc[201:400,i+12]<- apply(GEOsply, 1, mean, na.rm=T)
}

# # # # #
# 1.2. -  a2 ----
# 
# a2 : alpha diversity with Chao = 2

for(i in 1:4){
 
  attach(pwise_space)
  
  #LUI and GEO as predictors
  dat<- data.frame(dha2st, 
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
    
    if(sum(gdm_sub$coefficients)==0){
      LUsplx[1:200,m]<- 0
      LUsply[1:200,m]<- 0
      GEOsplx[1:200,m]<- 0
      GEOsply[1:200,m]<- 0
      
    }else{
      spline<- isplineExtract(gdm_sub)
      
      LUsplx[,m]<- spline$x[,2]
      LUsply[,m]<- spline$y[,2]
      GEOsplx[,m]<- spline$x[,1]
      GEOsply[,m]<- spline$y[,1] }
    
  }
  
  splines_alpha_herb_sc[401:600,i]<- apply(LUsplx, 1, mean, na.rm=T)
  splines_alpha_herb_sc[401:600,i+4]<- apply(GEOsplx, 1, mean, na.rm=T)
  
  splines_alpha_herb_sc[401:600,i+8]<- apply(LUsply, 1, mean, na.rm=T)
  splines_alpha_herb_sc[401:600,i+12]<- apply(GEOsply, 1, mean, na.rm=T)
}


# # # # #
# 1.3. -  a3 ----
# 
# a3 : alpha diversity with Chao = 3

for(i in 1:4){
  
  attach(pwise_space)
  
  #LUI and GEO as predictors
  dat<- data.frame(dha3st, 
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
    
    if(sum(gdm_sub$coefficients)==0){
      LUsplx[1:200,m]<- 0
      LUsply[1:200,m]<- 0
      GEOsplx[1:200,m]<- 0
      GEOsply[1:200,m]<- 0
      
    }else{
      spline<- isplineExtract(gdm_sub)
      
      LUsplx[,m]<- spline$x[,2]
      LUsply[,m]<- spline$y[,2]
      GEOsplx[,m]<- spline$x[,1]
      GEOsply[,m]<- spline$y[,1] }
    
  }
  
  splines_alpha_herb_sc[601:800,i]<- apply(LUsplx, 1, mean, na.rm=T)
  splines_alpha_herb_sc[601:800,i+4]<- apply(GEOsplx, 1, mean, na.rm=T)
  
  splines_alpha_herb_sc[601:800,i+8]<- apply(LUsply, 1, mean, na.rm=T)
  splines_alpha_herb_sc[601:800,i+12]<- apply(GEOsply, 1, mean, na.rm=T)
}


# # # # #
# 1.4. -  a4 ----
# 
# a4

for(i in 1:4){
 
  attach(pwise_space)
  
  #LUI and GEO as predictors
  dat<- data.frame(dha4st, 
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
    
    if(sum(gdm_sub$coefficients)==0){
      LUsplx[1:200,m]<- 0
      LUsply[1:200,m]<- 0
      GEOsplx[1:200,m]<- 0
      GEOsply[1:200,m]<- 0
      
    }else{
      spline<- isplineExtract(gdm_sub)
      
      LUsplx[,m]<- spline$x[,2]
      LUsply[,m]<- spline$y[,2]
      GEOsplx[,m]<- spline$x[,1]
      GEOsply[,m]<- spline$y[,1] }
    
  }
  
  splines_alpha_herb_sc[801:1000,i]<- apply(LUsplx, 1, mean, na.rm=T)
  splines_alpha_herb_sc[801:1000,i+4]<- apply(GEOsplx, 1, mean, na.rm=T)
  
  splines_alpha_herb_sc[801:1000,i+8]<- apply(LUsply, 1, mean, na.rm=T)
  splines_alpha_herb_sc[801:1000,i+12]<- apply(GEOsply, 1, mean, na.rm=T)
}

splines_alpha_herb_sc_space<- splines_alpha_herb_sc

save(splines_alpha_herb_sc_space, file="./data/OutputData/splines_alpha_herb_sc_space.RData")


# # # # # # # # # # # # # #
# 2 - TEMPORAL DATASET              ----
#

#create a data frame for the splines
splines_alpha_herb_sc<- as.data.frame(matrix(nrow=1000, ncol=16))
colnames(splines_alpha_herb_sc)<- c("LUIx", "MOWx", "GRAx", "FERx", 
                                  "COvLUIx","COvMOWx", "COvGRAx", "COvFERx",
                                  "LUIy", "MOWy", "GRAy", "FERy", 
                                  "COvLUIy","COvMOWy", "COvGRAy", "COvFERy")
splines_alpha_herb_sc$alpha_type<- rep(c("a0", "a1", "a2", "a3", "a4"), each=200)

LU1<- c("LUI1res", "MOW1res", "GRA1res", "FER1res")
LU2<- c("LUI2res", "MOW2res", "GRA2res", "FER2res")
#LU<- c("LUI", "MOW", "GRA", "FER")


# # # # #
# 2.0. -  a0 ----
# 
# a0 : alpha diversity with Chao = 0
# herbivores

for(i in 1:4){
  
  attach(pwise_time)
  
  #LUI and GEO as predictors
  dat<- data.frame(dha0st, 
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
  gdm_sub<- gdm(gdmTab, geo=F)
  
  
  if(sum(gdm_sub$coefficients)==0){
    splines_alpha_herb_sc[1:200,i]<- 0
    splines_alpha_herb_sc[1:200,i+4]<- 0
    
    splines_alpha_herb_sc[1:200,i+8]<- 0
    splines_alpha_herb_sc[1:200,i+12]<- 0
    
  }else{
    spline<- isplineExtract(gdm_sub)
    
    spline<- isplineExtract(gdm_sub)
    splines_alpha_herb_sc[1:200,i]<- spline$x[,1]#LUI
    splines_alpha_herb_sc[1:200,i+4]<- spline$x[,2]#year
    
    splines_alpha_herb_sc[1:200,i+8]<- spline$y[,1]
    splines_alpha_herb_sc[1:200,i+12]<- spline$y[,2] }
}

remove(dat, gdmTab)



# # # # #
# 2.1. -  a1 ----
# 
# a1 : alpha diversity with Chao = 1

for(i in 1:4){
  
  attach(pwise_time)
  
  #LUI and GEO as predictors
  dat<- data.frame(dha1st, 
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
  gdm_sub<- gdm(gdmTab, geo=F)
  
  if(sum(gdm_sub$coefficients)==0){
    splines_alpha_herb_sc[201:400,i]<- 0
    splines_alpha_herb_sc[201:400,i+4]<- 0
    
    splines_alpha_herb_sc[201:400,i+8]<- 0
    splines_alpha_herb_sc[201:400,i+12]<- 0
    
  } else {
    spline<- isplineExtract(gdm_sub)
    
    spline<- isplineExtract(gdm_sub)
    splines_alpha_herb_sc[201:400,i]<- spline$x[,1]#LUI
    splines_alpha_herb_sc[201:400,i+4]<- spline$x[,2]#year
    
    splines_alpha_herb_sc[201:400,i+8]<- spline$y[,1]
    splines_alpha_herb_sc[201:400,i+12]<- spline$y[,2] }
}


remove(dat, gdmTab)

# # # # #
# 2.2. -  a2 ----
# 
# a2 : alpha diversity with Chao = 2

for(i in 1:4){

  attach(pwise_time)
  
  #LUI and GEO as predictors
  dat<- data.frame(dha2st, 
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
  gdm_sub<- gdm(gdmTab, geo=F)
  
  
  if(sum(gdm_sub$coefficients)==0){
    splines_alpha_herb_sc[401:600,i]<- 0
    splines_alpha_herb_sc[401:600,i+4]<- 0
    
    splines_alpha_herb_sc[401:600,i+8]<- 0
    splines_alpha_herb_sc[401:600,i+12]<- 0
    
  }else{
    spline<- isplineExtract(gdm_sub)
    
    spline<- isplineExtract(gdm_sub)
    splines_alpha_herb_sc[401:600,i]<- spline$x[,1]#LUI
    splines_alpha_herb_sc[401:600,i+4]<- spline$x[,2]#year
    
    splines_alpha_herb_sc[401:600,i+8]<- spline$y[,1]
    splines_alpha_herb_sc[401:600,i+12]<- spline$y[,2] }
}

remove(dat, gdmTab)



# # # # #
# 2.3. -  a3 ----
# 
# a3 : alpha diversity with Chao = 3

for(i in 1:4){
 
  attach(pwise_time)
  
  #LUI and GEO as predictors
  dat<- data.frame(dha3st, 
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
  gdm_sub<- gdm(gdmTab, geo=F)
  
  if(sum(gdm_sub$coefficients)==0){
    splines_alpha_herb_sc[601:800,i]<- 0
    splines_alpha_herb_sc[601:800,i+4]<- 0
    
    splines_alpha_herb_sc[601:800,i+8]<- 0
    splines_alpha_herb_sc[601:800,i+12]<- 0
    
  }else{
    spline<- isplineExtract(gdm_sub)
    
    spline<- isplineExtract(gdm_sub)
    splines_alpha_herb_sc[601:800,i]<- spline$x[,1]#LUI
    splines_alpha_herb_sc[601:800,i+4]<- spline$x[,2]#year
    
    splines_alpha_herb_sc[601:800,i+8]<- spline$y[,1]
    splines_alpha_herb_sc[601:800,i+12]<- spline$y[,2] }
}


remove(dat, gdmTab)



# # # # #
# 2.4. -  a4 ----
# 
# a4

for(i in 1:4){
 
  attach(pwise_time)
  
  #LUI and GEO as predictors
  dat<- data.frame(dha4st, 
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
  gdm_sub<- gdm(gdmTab, geo=F)
  
  if(sum(gdm_sub$coefficients)==0){
    splines_alpha_herb_sc[801:1000,i]<- 0
    splines_alpha_herb_sc[801:1000,i+4]<- 0
    
    splines_alpha_herb_sc[801:1000,i+8]<- 0
    splines_alpha_herb_sc[801:1000,i+12]<- 0
    
  }else{
    spline<- isplineExtract(gdm_sub)
    
    spline<- isplineExtract(gdm_sub)
    splines_alpha_herb_sc[801:1000,i]<- spline$x[,1]#LUI
    splines_alpha_herb_sc[801:1000,i+4]<- spline$x[,2]#year
    
    splines_alpha_herb_sc[801:1000,i+8]<- spline$y[,1]
    splines_alpha_herb_sc[801:1000,i+12]<- spline$y[,2] }
}


remove(dat, gdmTab)

splines_alpha_herb_sc_time<- splines_alpha_herb_sc

save(splines_alpha_herb_sc_time, file="./data/OutputData/splines_alpha_herb_sc_time.RData")


#saving in the "plotting results folder"

#set working directory to folder "4. Plotting results"
save(splines_alpha_herb_sc_time, file="./data/InputData/splines_alpha_herb_sc_time.RData")


remove(pwise_space, pwise_time)
