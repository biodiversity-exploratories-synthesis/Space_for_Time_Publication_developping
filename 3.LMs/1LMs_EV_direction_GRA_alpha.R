# # # # # # # # # # # # # # #  # # # #
#                                    #
#             FIT LMs                # 
#         ALPHA DIVERSITY            #
#     land use = Grazing             # 
#                                    #
# # # # # # # # # # # # # # #  # # # #
#
# spatial and temporal dataset within this script

# # # # # # # # # # # # # #
# CONTENT                  ----
# AIM : Fit LM models of alpha diversity, for time and space datasets.


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


# # # # # # # # # # # # # #
# 0 - CALC SCALING              ----
# scaling all differences and means to max

####space
#plants
vars<- c("da0st", "da1st", "da2st", "da3st", "da4st", 
         "mGRA", "dGRA", "geo.dist", "YR", 
         "msur", "dsur","mph", "dph","msoilPC1", "dsoilPC1")
sub_plants_space<- pwise_space_plants[,colnames(pwise_space_plants) %in% vars]
sub_plants_space<- sub_plants_space[complete.cases(sub_plants_space),]

colnames(sub_plants_space)
preds<- apply(sub_plants_space[,c(2:10)], 2, function(x)x/max(x, na.rm=T))
sub_plants_space<- cbind(sub_plants_space[c(11:15, 1)], preds)
colnames(sub_plants_space)


#sec. consumers
vars<- c("dpa0st", "dpa1st", "dpa2st", "dpa3st", "dpa4st",
         "mGRA", "dGRA", "geo.dist", "YR",
         "msur", "dsur","mph", "dph","msoilPC1", "dsoilPC1")
sub_pred_space<- pwise_space_insects[,colnames(pwise_space_insects) %in% vars]
sub_pred_space<- sub_pred_space[complete.cases(sub_pred_space),]

colnames(sub_pred_space)
preds<- apply(sub_pred_space[,c(2:10)], 2, function(x)x/max(x, na.rm=T))
sub_pred_space<- cbind(sub_pred_space[c(11:15, 1)], preds)
colnames(sub_pred_space)

#herbivores
vars<- c("dha0st", "dha1st", "dha2st", "dha3st","dha4st", 
         "mGRA", "dGRA", "geo.dist", "YR",
         "msur", "dsur","mph", "dph","msoilPC1", "dsoilPC1")
sub_her_space<- pwise_space_insects[,colnames(pwise_space_insects) %in% vars]
sub_her_space<- sub_her_space[complete.cases(sub_her_space),]

colnames(sub_her_space)
preds<- apply(sub_her_space[,c(2:10)], 2, function(x)x/max(x, na.rm=T))
sub_her_space<- cbind(sub_her_space[c(11:15,1)], preds)
colnames(sub_her_space)

scaled_space_alpha<- list()
scaled_space_alpha[[1]]<- sub_plants_space
scaled_space_alpha[[2]]<- sub_her_space
scaled_space_alpha[[3]]<- sub_pred_space
names(scaled_space_alpha)<- tlorder

save(scaled_space_alpha, file="./data/InputData/scaled_alpha_space_GRA.RData")



####time
#plants
vars<- c("da0st", "da1st", "da2st", "da3st","da4st", 
         "mGRA", "dGRA", "dYR", "HWG", "RWG",
         "G500", "pH","soilPCA", "EP")
sub_plants_time<- pwise_time_plants[,colnames(pwise_time_plants) %in% vars]
sub_plants_time<- sub_plants_time[complete.cases(sub_plants_time),]

colnames(sub_plants_time)
preds<- apply(sub_plants_time[,c(2:9)], 2, function(x)x/max(x, na.rm=T))
sub_plants_time<- cbind(sub_plants_time[c(10:14,1)], preds)
colnames(sub_plants_time)

#sec. consumers
vars<- c("dpa0st", "dpa1st", "dpa2st", "dpa3st", "dpa4st",
         "mGRA", "dGRA", "dYR", "HWG", "RWG",
         "G500", "pH","soilPCA", "EP")
sub_pred_time<- pwise_time_insects[,colnames(pwise_time_insects) %in% vars]
sub_pred_time<- sub_pred_time[complete.cases(sub_pred_time),]

colnames(sub_pred_time)
preds<- apply(sub_pred_time[,c(2:9)], 2, function(x)x/max(x, na.rm=T))
sub_pred_time<- cbind(sub_pred_time[c(10:14,1)], preds)
colnames(sub_pred_time)


#herbivores
vars<- c("dha0st", "dha1st", "dha2st", "dha3st","dha4st", 
         "mGRA", "dGRA", "dYR", "HWG", "RWG",
         "G500", "pH","soilPCA", "EP")
sub_her_time<- pwise_time_insects[,colnames(pwise_time_insects) %in% vars]
sub_her_time<- sub_her_time[complete.cases(sub_her_time),]

colnames(sub_her_time)
preds<- apply(sub_her_time[,c(2:9)], 2, function(x)x/max(x, na.rm=T))
sub_her_time<- cbind(sub_her_time[c(10:14,1)], preds)
colnames(sub_her_time)

scaled_time_alpha<- list()
scaled_time_alpha[[1]]<- sub_plants_time
scaled_time_alpha[[2]]<- sub_her_time
scaled_time_alpha[[3]]<- sub_pred_time
names(scaled_time_alpha)<- tlorder

save(scaled_time_alpha, file="./data/InputData/scaled_alpha_time_GRA.RData")





# # # # # # # # # # # # # #
# 1 - LINEAR MODELS              ----
#
load("./data/InputData/scaled_alpha_space_GRA.RData")
load("./data/InputData/scaled_alpha_time_GRA.RData")

# # # # #
# 1.a. - SPATIAL DATASET ----
# 
EV_space_alpha<-list()
lm_space_alpha<- list()
#I ran this loop manually since sometimes convergence problems make the loop stop
for(k in 1:3){ 
  lm_space1<- list()
  
  sub<- scaled_space_alpha[[k]]
  sub$YR<- as.factor(sub$YR)
  
  names(sub)[1]<- "a0"
  names(sub)[2]<- "a1"
  names(sub)[3]<- "a2"
  names(sub)[4]<- "a3"
  names(sub)[5]<- "a4"
  
  
  lm11<- lme(a0~mGRA+dGRA+dph+mph+msur+dsur+dsoilPC1+msoilPC1+geo.dist,
             data=sub, random=~1|YR)
  lm22<- lme(a1~mGRA+dGRA+dph+mph+msur+dsur+dsoilPC1+msoilPC1+geo.dist,
             data=sub, random=~1|YR)
  lm33<- lme(a2~mGRA+dGRA+dph+mph+msur+dsur+dsoilPC1+msoilPC1+geo.dist,
             data=sub, random=~1|YR)
  lm44<- lme(a3~mGRA+dGRA+dph+mph+msur+dsur+dsoilPC1+msoilPC1+geo.dist,
             data=sub, random=~1|YR)
  lm55<- lme(a4~mGRA+dGRA+dph+mph+msur+dsur+dsoilPC1+msoilPC1+geo.dist,
             data=sub, random=~1|YR)
  
  lm_space1<- list(lm11,lm22, lm33, lm44, lm55)#if one model does not converge, set to list[0]
  
  all<- c()
  mGRApp<- c()
  dGRApp<- c()
  sharedpp<- c()
  cofm<- c()
  cofd<- c()
  int<- c()
  pvm<- c()
  pvd<- c()
  
  for(l in 1:5){
    
    lm<- lm_space1[[l]]
    
    if(length(lm)<10){
      all<- c(all,0)
      mGRApp<- c(mGRApp,0)
      dGRApp<- c(dGRApp,0)
      sharedpp<- c(sharedpp,0)
      cofm<- c(cofm,0)
      cofd<- c(cofd,0)
      int<- c(int,0)
      pvm<- c(pvm,1)
      pvd<- c(pvd,1)
      
    }else{
      
    suppressWarnings(lm1<- update(lm,~.-mGRA))
    suppressWarnings(lm2<- update(lm,~.-dGRA))
    suppressWarnings(lm3<- update(lm1,~.-dGRA))
      
    lm1<- update(lm,~.-mLUI)
    lm2<- update(lm,~.-dLUI)
    lm3<- update(lm1,~.-dLUI)
    
    a<- r.squaredGLMM(lm)[1]#marginal R2 of the fixed effects
    m<- (r.squaredGLMM(lm)[1]-r.squaredGLMM(lm1)[1])#variance expl by mLUI
    d<- (r.squaredGLMM(lm)[1]-r.squaredGLMM(lm2)[1])#variance expl by dLUI
    sh<-(r.squaredGLMM(lm)[1]-r.squaredGLMM(lm3)[1]) - m - d
    
    pp <- 100*(c(m, d, sh)/a)
    
        }
      
      all<- c(all, a)
      mGRApp<- c(mGRApp,as.numeric(pp[1]))
      dGRApp<- c(dGRApp,as.numeric(pp[2]))
      sharedpp<- c(sharedpp, as.numeric(pp[3]))
      cofm<- c(cofm,lm$coefficients$fixed[2])
      cofd<- c(cofd,lm$coefficients$fixed[3])
      int<- c(int, lm$coefficients$fixed[1])
      pvm<- c(pvm, round(Anova(lm, type="II")$`Pr(>Chisq)`[1], digits=3))
      pvd<- c(pvd, round(Anova(lm, type="II")$`Pr(>Chisq)`[2], digits=3))
      
      rm(lm1, lm2, lm3)
    }
    
  EV_space_alpha[[k]]<- data.frame(all, mGRApp, dGRApp, sharedpp, cofm, cofd, int, pvm, pvd)
  lm_space_alpha[[k]]<- lm_space1
}

all<- c(all, NA)
mGRApp<- c(mGRApp,NA)
dGRApp<- c(dGRApp,NA)
sharedpp<- c(sharedpp, NA)
cofm<- c(cofm,NA)
cofd<- c(cofd,NA)
int<- c(int, NA)
pvm<- c(pvm, NA)
pvd<- c(pvd, NA)

names(EV_space_alpha)<- tlorder
names(lm_space_alpha)<- tlorder

save(EV_space_alpha, file="./data/OutputData/EV_space_GRA_alpha.RData")
save(lm_space_alpha, file="./data/OutputData/lm_space_GRA_alpha.RData")


# # # # #
# 1.b. - TEMPORAL DATASET ----
# 
lm_time_alpha<- list()
EV_time_alpha<-list()
for(k in 1:3){ #I ran this loop manually since sometimes convergence problems make the loop stop
  all<- c()
  mGRApp<- c()
  dGRApp<- c()
  sharedpp<- c()
  cofm<- c()
  cofd<- c()
  int<- c()
  pvm<- c()
  pvd<- c()
  lm_time1<- list()
  
  sub<- scaled_time_alpha[[k]]
  names(sub)[1]<- "a0"
  names(sub)[2]<- "a1"
  names(sub)[3]<- "a2"
  names(sub)[4]<- "a3"
  names(sub)[5]<- "a4"
  
  
  lm11<- lme(a0~mGRA+dGRA+pH+G500+soilPCA+dYR,
             data=sub, random=~1|EP)
  lm22<- lme(a1~mGRA+dGRA+pH+G500+soilPCA+dYR,
             data=sub, random=~1|EP)
  lm33<- lme(a2~mGRA+dGRA+pH+G500+soilPCA+dYR,
             data=sub, random=~1|EP)
  lm44<- lme(a3~mGRA+dGRA+pH+G500+soilPCA+dYR,
             data=sub, random=~1|EP)
  lm55<- lme(a4~mGRA+dGRA+pH+G500+soilPCA+dYR,
             data=sub, random=~1|EP)
  lm_time1<- list(lm11, lm22, lm33, lm44, lm55)
  
  #beta 0,1,2,3,4, sim
  for(l in 1:5){
    lm<- lm_time1[[l]]
    lm1<- update(lm,~.-mGRA)
    lm2<- update(lm,~.-dGRA)
    lm3<- update(lm1,~.-dGRA)
    
    
    a<- r.squaredGLMM(lm)[1]#marginal R2 of the fixed effects
    m<- (r.squaredGLMM(lm)[1]-r.squaredGLMM(lm1)[1])#variance expl by mGRA
    d<- (r.squaredGLMM(lm)[1]-r.squaredGLMM(lm2)[1])#variance expl by dGRA
    sh<-(r.squaredGLMM(lm)[1]-r.squaredGLMM(lm3)[1]) - m - d
    
    pp <- 100*(c(m, d, sh)/a)
    
    
    
    all<- c(all, a)
    mGRApp<- c(mGRApp,as.numeric(pp[1]))
    dGRApp<- c(dGRApp,as.numeric(pp[2]))
    sharedpp<- c(sharedpp, as.numeric(pp[3]))
    cofm<- c(cofm,lm$coefficients$fixed[2])
    cofd<- c(cofd,lm$coefficients$fixed[3])
    int<- c(int, lm$coefficients$fixed[1])
    pvm<- c(pvm, round(Anova(lm, type="II")$`Pr(>Chisq)`[1], digits=3))
    pvd<- c(pvd, round(Anova(lm, type="II")$`Pr(>Chisq)`[2], digits=3))
  }
  
  EV_time_alpha[[k]]<- data.frame(all, mGRApp, dGRApp, sharedpp, cofm, cofd, int, pvm, pvd)
  lm_time_alpha[[k]]<- lm_time1
}

names(EV_time_alpha)<- tlorder
names(lm_time_alpha)<- tlorder

save(EV_time_alpha, file="./data/OutputData/EV_time_GRA_alpha.RData")
save(lm_time_alpha, file="./data/OutputData/lm_time_GRA_alpha.RData")

#set working directory to folder "4. Plotting results"
save(EV_time_alpha, file="./data/InputData/EV_time_GRA_alpha.RData")
