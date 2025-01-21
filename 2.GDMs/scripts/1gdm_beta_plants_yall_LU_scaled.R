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
#plants
#load("./pwise_time_new.RData")
load("./sitepred_pwise_time_step3.RData")#new pwise_time with updated LUI residuals



#TODO is this the range part?
#    saving this will not help here, absolutely need to update!

#dataprep - now done and saved
pwise_time$HWG1<- pwise_time$HWG

pwise_time$RWG1<- pwise_time$RWG
pwise_time$RWG2<- pwise_time$RWG

YR1<- rep(rep(c(2008:2017), times=(10:1)), times=150)
pwise_time$YR1<- YR1

YR2<- c(c(2009:2018), c(2010:2018), c(2011:2018),
        c(2012:2018), c(2013:2018), c(2014:2018),
        c(2015:2018), c(2016:2018), c(2017:2018),
        2018)
pwise_time$YR2<- rep(YR2, times=150)
save(pwise_time, file="pwise_time_new.RData")




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

save(splines_plants_yall_sc, file="splines_plants_yall_sc_new.RData")


# # # # #
# 2.5. -  PLOT ----
# 
splines_plants_yall_sc$beta_type<- as.factor(splines_plants_yall_sc$beta_type)
levels(splines_plants_yall_sc$beta_type)

load("./splines_plants_yall_sc_new.RData")

# LUI
p1<- ggplot(splines_plants_yall_sc, aes(x=LUIx, y=LUIy, colour=beta_type))+
  geom_line(size=1.5)+
  scale_y_sqrt(limits = c(0, 3), expand = c(0, 0))+
  scale_colour_hue(name="Type of Beta Diversity",
                    labels=c("Turnover",  "Chao 0",
                              "Chao 1", "Chao 2",
                             "Chao 3", "Chao 4"),
                   l=50)+
  xlab("LUI (scaled to EP average)")+
  ylab("Effect of LUI differences on beta-diversity")+
  ggtitle("a - Land use intensity (LUI)")+
  #geom_text(aes(x = 3.8, y = 0.25, label = "EV% = 1.5"), color="black") + 
  #geom_text(aes(x = 3.8, y = 0.19, label = "EV% = 2"),color="black") + 
  #geom_text(aes(x = 3.8, y = 0.15, label = "EV% = 1.3"),color="black") + 
  #geom_text(aes(x = 3.8, y = 0.07, label = "EV% = 1.3"),color="black") + 
  theme_classic()

p1

# Mowing
p2<- ggplot(splines_plants_yall_sc, aes(x=MOWx, y=MOWy, colour=beta_type))+
  geom_line(size=1.5)+
  scale_y_sqrt(limits = c(0, 3), expand = c(0, 0))+
  scale_colour_hue(name="Type of Beta Diversity", labels=c("Turnover", "Chao 0",
                                                           "Chao 1", "Chao 2",
                                                           "Chao 3", "Chao 4"),l=50)+
  xlab("MOW (scaled to EP average)")+
  ylab("Effect of MOW differences on beta-diversity")+
  ggtitle("b - Mowing intensity (MOW)")+
  #geom_text(aes(x = 4.5, y = 0.19, label = "EV% = 0.6"), color="black") + 
  #geom_text(aes(x = 4.5, y = 0.10, label = "EV% = 1.1"),color="black") + 
  #geom_text(aes(x = 4.5, y = 0.08, label = "EV% = 0.5"),color="black") + 
  #geom_text(aes(x = 4.5, y = 0.05, label = "EV% = 1.1"),color="black") + 
  theme_classic()
p2

# Grazing
p3<- ggplot(splines_plants_yall_sc, aes(x=GRAx, y=GRAy, colour=beta_type))+
  geom_line(size=1.5)+
  scale_y_sqrt(limits = c(0, 3), expand = c(0, 0))+
  scale_colour_hue(name="Type of Beta Diversity", labels=c("Turnover", "Chao 0",
                                                           "Chao 1", "Chao 2",
                                                           "Chao 3", "Chao 4"),l=50)+
  xlab("GRA (scaled to EP average)")+
  ylab("Effect of GRA differences on beta-diversity")+
  ggtitle("c - Grazing intensity (GRA)")+
  #geom_text(aes(x = 11.5, y = 0.16, label = "EV% = 0.7"), color="black") + 
  #geom_text(aes(x = 11.5, y = 0.11, label = "EV% = 1.4"),color="black") + 
  #geom_text(aes(x = 11.5, y = 0.09, label = "EV% = 1.8"),color="black") + 
  #geom_text(aes(x = 11.5, y = 0.055, label = "EV% = 1.6"),color="black") + 
  theme_classic()

p3

# Fertilisation
p4<- ggplot(splines_plants_yall_sc, aes(x=FERx, y=FERy, colour=beta_type))+
  geom_line(size=1.5)+
  scale_y_sqrt(limits = c(0, 3), expand = c(0, 0))+
  scale_colour_hue(name="Type of Beta Diversity", labels=c("Turnover", "Chao 0",
                                                           "Chao 1", "Chao 2",
                                                           "Chao 3", "Chao 4"),l=50)+
  xlab("FER (scaled to EP average)")+
  ylab("Effect of FER differences on beta-diversity")+
  ggtitle("d - Fertilisation intensity (FER)")+
  #geom_text(aes(x = 12.5, y = 0.19, label = "EV% = 1.5"), color="black") + 
  #geom_text(aes(x = 12.5, y = 0.16, label = "EV% = 0.9"),color="black") + 
  #geom_text(aes(x = 12.5, y = 0.14, label = "EV% = 0.9"),color="black") + 
  #geom_text(aes(x = 12.5, y = 0.11, label = "EV% = 1.1"),color="black") + 
  theme_classic()

p4

p1 + p2 + p3 + p4 + plot_layout(guides="collect")
ggsave("ResultsPlots/plant_yall_LUsplines_sc.png")
