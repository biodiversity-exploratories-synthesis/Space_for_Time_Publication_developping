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
# insects
#load("./pwise_space_new.RData")
#TODO clean above
load("./sitepred_pwise_space_step3.RData")#new pwise_space with updated LUI residuals


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


save(splines_her_space_sc, file="InterimDatasets/splines_her_space_sc_new.RData")


# # # # #
# 1.5. -  PLOT ----
# 
splines_her_space_sc$beta_type<- as.factor(splines_her_space_sc$beta_type)
levels(splines_her_space_sc$beta_type)

load("./splines_her_space_sc_new.RData")

# LUI
p1<- ggplot(splines_her_space_sc, aes(x=LUIx, y=LUIy, colour=beta_type))+
  geom_line(size=1.5)+
  scale_y_sqrt(limits = c(0, 3), expand = c(0, 0))+
  scale_colour_hue(name="Type of Beta Diversity", labels=c("Turnover", "Chao 0",
                                                           "Chao 1", "Chao 2",
                                                           "Chao 3", "Chao 4"),l=50)+
  xlab("LUI (scaled to year average)")+
  ylab("Effect of LUI differences on beta-diversity")+
  ggtitle("a - Land use intensity (LUI)")+
  #geom_text(aes(x = 3.9, y = 1.5, label = "EV% = 2.8"), color="black") + 
  #geom_text(aes(x = 3.9, y = 1.1, label = "EV% = 4.0"),color="black") + 
  #geom_text(aes(x = 3.9, y = 0.8, label = "EV% = 3.5"),color="black") + 
  #geom_text(aes(x = 3.9, y = 0.5, label = "EV% = 3.1"),color="black") + 
  theme_classic()
p1

# Mowing
p2<- ggplot(splines_her_space_sc, aes(x=MOWx, y=MOWy, colour=beta_type))+
  geom_line(size=1.5)+
  scale_y_sqrt(limits = c(0, 3), expand = c(0, 0))+
  scale_colour_hue(name="Type of Beta Diversity", labels=c("Turnover", "Chao 0",
                                                           "Chao 1", "Chao 2",
                                                           "Chao 3", "Chao 4"),l=50)+
  xlab("MOW (scaled to year average)")+
  ylab("Effect of MOW differences on beta-diversity")+
  ggtitle("b - Mowing intensity (MOW)")+
  #geom_text(aes(x = 3.5, y = 0.19, label = "EV% = 0.6"), color="black") + 
  #geom_text(aes(x = 3.5, y = 0.145, label = "EV% = 0.9"),color="black") + 
  #geom_text(aes(x = 3.5, y = 0.123, label = "EV% = 0.6"),color="black") + 
  #geom_text(aes(x = 3.5, y = 0.09, label = "EV% = 0.5"),color="black") + 
  theme_classic()
p2

# Grazing
p3<- ggplot(splines_her_space_sc, aes(x=GRAx, y=GRAy, colour=beta_type))+
  geom_line(size=1.5)+
  scale_y_sqrt(limits = c(0, 3), expand = c(0, 0))+
  scale_colour_hue(name="Type of Beta Diversity", labels=c("Turnover", "Chao 0",
                                                           "Chao 1", "Chao 2",
                                                           "Chao 3", "Chao 4"),l=50)+
  xlab("GRA (scaled to year average)")+
  ylab("Effect of GRA differences on beta-diversity")+
  ggtitle("c - Grazing intensity (GRA)")+
  #geom_text(aes(x = 11.5, y = 0.01, label = "EV% < 0.1"), color="black") + 
  theme_classic()

p3

# Fertilisation
p4<- ggplot(splines_her_space_sc, aes(x=FERx, y=FERy, colour=beta_type))+
  geom_line(size=1.5)+
  scale_y_sqrt(limits = c(0, 3), expand = c(0, 0))+
  scale_colour_hue(name="Type of Beta Diversity", labels=c("Turnover", "Chao 0",
                                                           "Chao 1", "Chao 2",
                                                           "Chao 3", "Chao 4"),l=50)+
  xlab("FER (scaled to year average)")+
  ylab("Effect of FER differences on beta-diversity")+
  ggtitle("d - Fertilisation intensity (FER)")+
  #geom_text(aes(x = 13.5, y = 0.055, label = "EV% = 0.1"), color="black") + 
  #geom_text(aes(x = 13.5, y = 0.032, label = "EV% = 0.1"),color="black") + 
  #geom_text(aes(x = 13.5, y = 0.02, label = "EV% = 0.1"),color="black") + 
  #geom_text(aes(x = 13.5, y = 0.007, label = "EV% < 0.1"),color="black") + 
  theme_classic()

p4


p1 + p2 + p3 + p4 + plot_layout(guides="collect")
ggsave("ResultsPlots/herb_space_LUsplines_scaled.png")
