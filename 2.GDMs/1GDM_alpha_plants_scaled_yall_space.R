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
#load("./pwise_space_new.RData")
#load("./pwise_time_new.RData")
#TODO above is deprecated, right?
load("./sitepred_pwise_space_step3.RData")#new pwise_space with updated LUI residuals
load("./sitepred_pwise_time_step3.RData")#new pwise_space with updated LUI residuals


# Preparing datasets
#
# make sure distances are not <0, >1 --> restrict to a tiny bit >0 and <1 --> DO NOT SAVE this
# vegan::decostand()
#TODO is this only for alphadiversity? not for beta plants, 
pwise_space$da0st<- decostand(pwise_space$da0, method="range", na.rm=T)
pwise_space$da1st<- decostand(pwise_space$da1, method="range", na.rm=T)
pwise_space$da2st<- decostand(pwise_space$da2, method="range", na.rm=T)
pwise_space$da3st<- decostand(pwise_space$da3, method="range", na.rm=T)
pwise_space$da4st<- decostand(pwise_space$da4, method="range", na.rm=T)

pwise_time$da0st<- decostand(pwise_time$da0abs, method="range", na.rm=T)
pwise_time$da1st<- decostand(pwise_time$da1abs, method="range", na.rm=T)
pwise_time$da2st<- decostand(pwise_time$da2abs, method="range", na.rm=T)
pwise_time$da3st<- decostand(pwise_time$da3abs, method="range", na.rm=T)
pwise_time$da4st<- decostand(pwise_time$da4abs, method="range", na.rm=T)




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

save(splines_alpha_plant_space_sc, file="splines_alpha_plant_space_sc.RData")

# write.table(splines_alpha_plant_space_sc, file="splines.txt")



# # # # #
# 1.5. -  PLOT ----
# 
load("./splines_alpha_plant_space_sc.RData")

# LUI
p1<- ggplot(splines_alpha_plant_space_sc, aes(x=LUIx, y=LUIy, colour=alpha_type))+
  geom_line(size=1.5)+
  scale_y_sqrt(limits = c(0, 0.5), expand = c(0, 0))+
  scale_colour_hue(name="Type of Alpha Diversity", l=50)+
  xlab("LUI (scaled to year average")+
  ylab("Effect of LUI differences on alpha-diversity")+
  ggtitle("a - Land use intensity (LUI)")+
  #geom_text(aes(x = 3.9, y = 0.45, label = "EV% = 0.1"), color="black") + 
  theme_classic()
p1

# Mowing
p2<- ggplot(splines_alpha_plant_space_sc, aes(x=MOWx, y=MOWy, colour=alpha_type))+
  geom_line(size=1.5)+
  scale_y_sqrt(limits = c(0, 0.5), expand = c(0, 0))+
  scale_colour_hue(name="Type of Alpha Diversity", l=50)+
  xlab("MOW (scaled to year average)")+
  ylab("Effect of MOW differences on alpha-diversity")+
  ggtitle("b - Mowing frequency (MOW)")+
  #geom_text(aes(x = 3.9, y = 0.45, label = "EV% = 0.1"), color="black") + 
  theme_classic()
p2

# Grazing
p3<- ggplot(splines_alpha_plant_space_sc, aes(x=GRAx, y=GRAy, colour=alpha_type))+
  geom_line(size=1.5)+
  scale_y_sqrt(limits = c(0, 0.5), expand = c(0, 0))+
  scale_colour_hue(name="Type of Alpha Diversity", l=50)+
  xlab("GRA")+
  ylab("Effect of GRA differences on alpha-diversity")+
  ggtitle("c - Grazing intensity (GRA)")+
  #geom_text(aes(x = 3.9, y = 0.45, label = "EV% = 0.1"), color="black") + 
  theme_classic()
p3

# Fertilisation
p4<- ggplot(splines_alpha_plant_space_sc, aes(x=FERx, y=FERy, colour=alpha_type))+
  geom_line(size=1.5)+
  scale_y_sqrt(limits = c(0, 0.5), expand = c(0, 0))+
  scale_colour_hue(name="Type of Alpha Diversity", l=50)+
  xlab("FER (scaled to year average)")+
  ylab("Effect of FER differences on alpha-diversity")+
  ggtitle("d - Fertilisation intensity (FER)")+
  #geom_text(aes(x = 3.9, y = 0.45, label = "EV% = 0.1"), color="black") + 
  theme_classic()
p4

p1 + p2 + p3 + p4 + plot_layout(guides="collect")
ggsave("ResultsPlots/plant_space_LUsplines_alpha_scaled.png")




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

save(splines_alpha_plant_time_sc, file="splines_alpha_plant_time_sc.RData")

#write.table(splines_alpha_plant_space_sc, file="splines.txt")


# # # # #
# 2.5. -  PLOT ----
# 

# LUI index
p1<- ggplot(splines_alpha_plant_time_sc, aes(x=LUIx, y=LUIy, colour=alpha_type))+
  geom_line(size=1.5)+
  scale_y_sqrt(limits = c(0, 0.5), expand = c(0, 0))+
  scale_colour_hue(name="Type of Alpha Diversity", l=50)+
  xlab("LUI (scaled to EP average)")+
  ylab("Effect of LUI differences on alpha-diversity")+
  ggtitle("a - Land use intensity (LUI)")+
  #geom_text(aes(x = 3.9, y = 0.45, label = "EV% = 0.1"), color="black") + 
  theme_classic()
p1

# Mowing
p2<- ggplot(splines_alpha_plant_time_sc, aes(x=MOWx, y=MOWy, colour=alpha_type))+
  geom_line(size=1.5)+
  scale_y_sqrt(limits = c(0, 0.5), expand = c(0, 0))+
  scale_colour_hue(name="Type of Alpha Diversity", l=50)+
  xlab("MOW (scaled to EP average)")+
  ylab("Effect of MOW differences on alpha-diversity")+
  ggtitle("b - Mowing frequency (MOW)")+
  #geom_text(aes(x = 3.9, y = 0.45, label = "EV% = 0.1"), color="black") + 
  theme_classic()
p2

# Grazing
p3<- ggplot(splines_alpha_plant_time_sc, aes(x=GRAx, y=GRAy, colour=alpha_type))+
  geom_line(size=1.5)+
  scale_y_sqrt(limits = c(0, 0.5), expand = c(0, 0))+
  scale_colour_hue(name="Type of Alpha Diversity", l=50)+
  xlab("GRA (scaled to EP average)")+
  ylab("Effect of GRA differences on alpha-diversity")+
  ggtitle("c - Grazing intensity (GRA)")+
  #geom_text(aes(x = 3.9, y = 0.45, label = "EV% = 0.1"), color="black") + 
  theme_classic()
p3

# Fertilisation
p4<- ggplot(splines_alpha_plant_time_sc, aes(x=FERx, y=FERy, colour=alpha_type))+
  geom_line(size=1.5)+
  scale_y_sqrt(limits = c(0, 0.5), expand = c(0, 0))+
  scale_colour_hue(name="Type of Alpha Diversity", l=50)+
  xlab("FER (scaled to EP average)")+
  ylab("Effect of FER differences on alpha-diversity")+
  ggtitle("d - Fertilisation intensity (FER)")+
  #geom_text(aes(x = 3.9, y = 0.45, label = "EV% = 0.1"), color="black") + 
  theme_classic()
p4

p1 + p2 + p3 + p4 + plot_layout(guides="collect")
ggsave("ResultsPlots/plant_time_LUsplines_alpha_scaled.png")
