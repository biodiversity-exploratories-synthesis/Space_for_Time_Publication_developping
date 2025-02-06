#load data - too large to load all at once
# # # # # # # # # # # # # # #  # # # #
#                                    #
#             FIT LMs                # 
#         Beta DIVERSITY             #
#     land use = Grazing             # 
#                                    #
# # # # # # # # # # # # # # #  # # # #
#
# spatial and temporal dataset within this script

# # # # # # # # # # # # # #
# CONTENT                  ----
# AIM : Fit LM models of beta diversity, for time and space datasets.


# # # # # # # # # # # # # #
# 0 - REQUIREMENTS              ----
#
# packages
# install.packages("gdm")
library(lme4)
library(nlme)
library(lmPerm)
library(MuMIn)#r.quared.GLMM
library(car)#Anova(type="II")

tlorder <- c("plants", "herbivores", "secondary.consumers")

# # # # #
# 0.a. - DATA ----
#
#set working directory to folder "3.LMs"
#herbivores
load("./data/InputData/pwise_time_herb.RData")
load("./data/InputData/pwise_space_herb.RData")

#secondary consumers/predators
load("./data/InputData/pwise_time_pred.RData")
load("./data/InputData/pwise_space_pred.RData")

#plants
load("./data/InputData/pwise_time_plants.RData")
load("./data/InputData/pwise_space_plants.RData")

# if the uploaded, assembled data files are used, upload the complete insect 
# files and plant files here.

#################PREPARATION####################
###############################################
#scaling all differences and means to max
####space
#plants
vars<- c("cqn0dis", "cqn1dis", "cqn2dis", "cqn3dis", "cqn4dis", "bsim", 
         "mGRA", "dGRA", "geo.dist", "YR", 
         "msur", "dsur","mph", "dph","msoilPC1", "dsoilPC1")
sub_plants_space<- pwise_space_plants[,colnames(pwise_space_plants) %in% vars]
sub_plants_space<- sub_plants_space[complete.cases(sub_plants_space),]

colnames(sub_plants_space)
preds<- apply(sub_plants_space[,c(1,2, 7:14)], 2, function(x)x/max(x, na.rm=T))
sub_plants_space<- cbind(sub_plants_space[c(3:6, 15,16)], preds)
colnames(sub_plants_space)

#sec. consumers
vars<- c("pcqn0dis", "pcqn1dis", "pcqn2dis", "pcqn3dis","pcqn4dis","pbsim", 
         "mGRA", "dGRA", "geo.dist", "YR",
         "msur", "dsur","mph", "dph","msoilPC1", "dsoilPC1")
sub_pred_space<- pwise_space_insects[,colnames(pwise_space_insects) %in% vars]
sub_pred_space<- sub_pred_space[complete.cases(sub_pred_space),]

colnames(sub_pred_space)
preds<- apply(sub_pred_space[,c(1,2,6:13)], 2, function(x)x/max(x, na.rm=T))
sub_pred_space<- cbind(sub_pred_space[c(14,3:5, 15, 16)], preds)
colnames(sub_pred_space)

#herbivores
vars<- c("hcqn0dis", "hcqn1dis", "hcqn2dis", "hcqn3dis","hcqn4dis","hbsim", 
         "mGRA", "dGRA", "geo.dist", "YR",
         "msur", "dsur","mph", "dph","msoilPC1", "dsoilPC1")
sub_her_space<- pwise_space_insects[,colnames(pwise_space_insects) %in% vars]
sub_her_space<- sub_her_space[complete.cases(sub_her_space),]

colnames(sub_her_space)
preds<- apply(sub_her_space[,c(1,2,7:14)], 2, function(x)x/max(x, na.rm=T))
sub_her_space<- cbind(sub_her_space[c(3:6, 15, 16)], preds)
colnames(sub_her_space)

scaled_space<- list()
scaled_space[[1]]<- sub_plants_space
scaled_space[[2]]<- sub_her_space
scaled_space[[3]]<- sub_pred_space
names(scaled_space)<- tlorder

save(scaled_space, file="./data/InputData/scaled_beta_space_GRA.RData")


####time####
###########
#plants
vars<- c("cqn0dis", "cqn1dis", "cqn2dis", "cqn3dis","cqn4dis","bsim", 
         "mGRA", "dGRA", "dYR", "HWG", "RWG",
         "G500", "pH","soilPCA")
sub_plants_time<- pwise_time_plants[,colnames(pwise_time_plants) %in% vars]
sub_plants_time<- sub_plants_time[complete.cases(sub_plants_time),]

colnames(sub_plants_time)
preds<- apply(sub_plants_time[,c(1:3, 8:12)], 2, function(x)x/max(x, na.rm=T))
sub_plants_time<- cbind(sub_plants_time[c(4:7, 13,14)], preds)
colnames(sub_plants_time)

#sec. consumers
vars<- c("pcqn0dis", "pcqn1dis", "pcqn2dis", "pcqn3dis","pcqn4dis","pbsim", 
         "mGRA", "dGRA", "dYR", "HWG", "RWG",
         "G500", "pH","soilPCA")
sub_pred_time<- pwise_time_insects[,colnames(pwise_time_insects) %in% vars]
sub_pred_time<- sub_pred_time[complete.cases(sub_pred_time),]

colnames(sub_pred_time)
preds<- apply(sub_pred_time[,c(1:4,6:8,12)], 2, function(x)x/max(x, na.rm=T))
sub_pred_time<- cbind(sub_pred_time[c(5,9:11, 13,14)], preds)
colnames(sub_pred_time)


#herbivores
vars<- c("hcqn0dis", "hcqn1dis", "hcqn2dis","hcqn3dis","hcqn4dis", "hbsim", 
         "mGRA", "dGRA", "dYR", "HWG", "RWG",
         "G500", "pH","soilPCA")
sub_her_time<- pwise_time_insects[,colnames(pwise_time_insects) %in% vars]
sub_her_time<- sub_her_time[complete.cases(sub_her_time),]

colnames(sub_her_time)
preds<- apply(sub_her_time[,c(1:4,6:8,12)], 2, function(x)x/max(x, na.rm=T))
sub_her_time<- cbind(sub_her_time[c(5,9:11, 13,14)], preds)
colnames(sub_her_time)

scaled_time<- list()
scaled_time[[1]]<- sub_plants_time
scaled_time[[2]]<- sub_her_time
scaled_time[[3]]<- sub_pred_time
names(scaled_time)<- tlorder

save(scaled_time, file="./data/InputData/scaled_beta_time_GRA.RData")

########################
###LINEAR MODELS
#######################
load("./data/InputData/scaled_beta_space_GRA.RData")
load("./data/InputData/scaled_beta_time_GRA.RData")

#space
EV_space<-list()
for(k in 1:3){
  all<- c()
  mGRApp<- c()
  dGRApp<- c()
  sharedpp<- c()
  cofm<- c()
  cofd<- c()
  int<- c()
  pvm<- c()
  pvd<-c()
  
  
  for(l in 1:6){
    sub<- scaled_space[[k]]
    lm<- lm(as.formula(paste(colnames(sub[l]),"~mGRA+dGRA+dph+mph+msur+dsur+dsoilPC1+msoilPC1+geo.dist+YR",sep="")),data=sub)
    lm1<- update(lm,~.-mGRA)
    lm2<- update(lm,~.-dGRA)
    lm3<- update(lm1,~.-dGRA)
    
    
    a<- summary(lm)$r.squared
    m<- (summary(lm)$r.squared-summary(lm1)$r.squared)#variance expl by mGRA
    d<- (summary(lm)$r.squared-summary(lm2)$r.squared)#variance expl by dGRA
    sh<-(summary(lm)$r.squared-summary(lm3)$r.squared) - m - d
    
    pp <- 100*(c(m, d, sh)/a)
    
    
    
    all<- c(all, a)
    mGRApp<- c(mGRApp,as.numeric(pp[1]))
    dGRApp<- c(dGRApp,as.numeric(pp[2]))
    sharedpp<- c(sharedpp, as.numeric(pp[3]))
    cofm<- c(cofm,lm$coefficients[2])
    cofd<- c(cofd,lm$coefficients[3])
    int<- c(int, lm$coefficients[1])
    pvm<- c(pvm, round(anova(lm)$`Pr(>F)`[1], digits=3))
    pvd<- c(pvd, round(anova(lm)$`Pr(>F)`[2], digits=3))
  }
  EV_space[[k]]<- data.frame(all, mGRApp, dGRApp, sharedpp, cofm, cofd, int, pvm, pvd)
}

names(EV_space)<- tlorder

save(EV_space, file="./data/OutputData/EV_space_GRA.RData")



#time
EV_time<-list()
for(k in 1:3){
  all<- c()
  mGRApp<- c()
  dGRApp<- c()
  sharedpp<- c()
  cofm<- c()
  cofd<- c()
  int<- c()
  pvm<- c()
  pvd<- c()
  
  #beta 0,1,2,sim
  for(l in 1:6){
    sub<- scaled_time[[k]]
    lm<- lm(as.formula(paste(colnames(sub[l]),"~mGRA+dGRA+pH+G500+soilPCA+HWG+RWG+dYR",sep="")),data=sub)
    lm1<- update(lm,~.-mGRA)
    lm2<- update(lm,~.-dGRA)
    lm3<- update(lm1,~.-dGRA)
    
    
    a<- summary(lm)$r.squared
    m<- (summary(lm)$r.squared-summary(lm1)$r.squared)
    d<- (summary(lm)$r.squared-summary(lm2)$r.squared)
    sh<-(summary(lm)$r.squared-summary(lm3)$r.squared) -m - d
    
    pp <- 100*(c(m, d, sh)/a)
    
    
    
    all<- c(all, a)
    mGRApp<- c(mGRApp,as.numeric(pp[1]))
    dGRApp<- c(dGRApp,as.numeric(pp[2]))
    sharedpp<- c(sharedpp, as.numeric(pp[3]))
    cofm<- c(cofm,lm$coefficients[2])
    cofd<- c(cofd,lm$coefficients[3])
    int<- c(int, lm$coefficients[1])
    pvm<- c(pvm, round(anova(lm)$`Pr(>F)`[1], digits=3))
    pvd<- c(pvd, round(anova(lm)$`Pr(>F)`[2], digits=3))
  }
  EV_time[[k]]<- data.frame(all, mGRApp, dGRApp, sharedpp, cofm, cofd, int, pvm, pvd)
}

names(EV_time)<- tlorder

save(EV_time, file="./data/OutputData/EV_time_GRA_beta.RData")

#set working directory to folder "4. Plotting results"
save(EV_time, file="./data/InputData/EV_time_GRA_beta.RData")
