#load data - too large to load all at once
# # # # # # # # # # # # # # #  # # # #
#                                    #
#             FIT LMs                # 
#         alpha DIVERSITY           #
#     land use = LUI index           # 
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
vars<- c("da0st", "da1st", "da2st", "da3st", "da4st", 
         "mLUI", "dLUI", "geo.dist", "YR", 
         "msur", "dsur","mph", "dph","msoilPC1", "dsoilPC1")
sub_plants_space<- pwise_space_plants[,colnames(pwise_space_plants) %in% vars]
sub_plants_space<- sub_plants_space[complete.cases(sub_plants_space),]

colnames(sub_plants_space)
preds<- apply(sub_plants_space[,c(2:10)], 2, function(x)x/max(x, na.rm=T))
sub_plants_space<- cbind(sub_plants_space[c(11:15, 1)], preds)
colnames(sub_plants_space)


#sec. consumers
vars<- c("dpa0st", "dpa1st", "dpa2st", "dpa3st", "dpa4st",
         "mLUI", "dLUI", "geo.dist", "YR",
         "msur", "dsur","mph", "dph","msoilPC1", "dsoilPC1")
sub_pred_space<- pwise_space_insects[,colnames(pwise_space_insects) %in% vars]
sub_pred_space<- sub_pred_space[complete.cases(sub_pred_space),]

colnames(sub_pred_space)
preds<- apply(sub_pred_space[,c(2:10)], 2, function(x)x/max(x, na.rm=T))
sub_pred_space<- cbind(sub_pred_space[c(11:15, 1)], preds)
colnames(sub_pred_space)

#herbivores
vars<- c("dha0st", "dha1st", "dha2st", "dha3st","dha4st", 
         "mLUI", "dLUI", "geo.dist", "YR",
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
#scaled_space[[3]]<- sub_pol_space
#scaled_space[[4]]<- sub_dec_space
scaled_space_alpha[[3]]<- sub_pred_space
names(scaled_space_alpha)<- tlorder

save(scaled_space_alpha, file="./data/InputData/scaled_alpha_space_LUI.RData")


####time####
###########
#plants
vars<- c("da0st", "da1st", "da2st", "da3st","da4st", 
         "mLUI", "dLUI", "dYR", "HWG", "RWG",
         "G500", "pH","soilPCA", "EP")
sub_plants_time<- pwise_time_plants[,colnames(pwise_time_plants) %in% vars]
sub_plants_time<- sub_plants_time[complete.cases(sub_plants_time),]

colnames(sub_plants_time)
preds<- apply(sub_plants_time[,c(2:9)], 2, function(x)x/max(x, na.rm=T))
sub_plants_time<- cbind(sub_plants_time[c(10:14,1)], preds)
colnames(sub_plants_time)

#sec. consumers
vars<- c("dpa0st", "dpa1st", "dpa2st", "dpa3st", "dpa4st",
         "mLUI", "dLUI", "dYR", "HWG", "RWG",
         "G500", "pH","soilPCA", "EP")
sub_pred_time<- pwise_time_insects[,colnames(pwise_time_insects) %in% vars]
sub_pred_time<- sub_pred_time[complete.cases(sub_pred_time),]

colnames(sub_pred_time)
preds<- apply(sub_pred_time[,c(2:9)], 2, function(x)x/max(x, na.rm=T))
sub_pred_time<- cbind(sub_pred_time[c(10:14,1)], preds)
colnames(sub_pred_time)


#herbivores
vars<- c("dha0st", "dha1st", "dha2st", "dha3st","dha4st", 
         "mLUI", "dLUI", "dYR", "HWG", "RWG",
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
#scaled_time[[3]]<- sub_pol_time
#scaled_time[[4]]<- sub_dec_time
scaled_time_alpha[[3]]<- sub_pred_time
names(scaled_time_alpha)<- tlorder

save(scaled_time_alpha, file="./data/InputData/scaled_alpha_time_LUI.RData")

########################
###LINEAR MODELS
#######################
load("./data/InputData/scaled_alpha_space_LUI.RData")
load("./data/InputData/scaled_alpha_time_LUI.RData")

#space
EV_space_alpha<-list()
lm_space_alpha<- list()

#for(k in 1:3){ I ran this loop manually since sometimes convergence problems make the loop stop
  all<- c()
  mLUIpp<- c()
  dLUIpp<- c()
 sharedpp<- c()
 cofm<- c()
 cofd<- c()
 int<- c()
 pvm<- c()
 pvd<- c()
 lm_space1<- list()
 
 sub<- scaled_space_alpha[[k]]
 sub$YR<- as.factor(sub$YR)
 levels(sub$YR)
 levels(sub$YR) <- c("1","2","3", "4", "5", "6", "7", "8", "9", "10")
 
 names(sub)[1]<- "a0"
 names(sub)[2]<- "a1"
 names(sub)[3]<- "a2"
 names(sub)[4]<- "a3"
 names(sub)[5]<- "a4"
 

 lm11<- lme(a0~mLUI+dLUI+dph+mph+msur+dsur+dsoilPC1+msoilPC1+geo.dist,
           data=sub, random=~1|YR)
 lm22<- lme(a1~mLUI+dLUI+dph+mph+msur+dsur+dsoilPC1+msoilPC1+geo.dist,
           data=sub, random=~1|YR)
 lm33<- lme(a2~mLUI+dLUI+dph+mph+msur+dsur+dsoilPC1+msoilPC1+geo.dist,
           data=sub, random=~1|YR)
 lm44<- lme(a3~mLUI+dLUI+dph+mph+msur+dsur+dsoilPC1+msoilPC1+geo.dist,
            data=sub, random=~1|YR)
 lm55<- lme(a4~mLUI+dLUI+dph+mph+msur+dsur+dsoilPC1+msoilPC1+geo.dist,
            data=sub, random=~1|YR)
 
 lm_space1<- list(lm11,lm22, lm33, lm44, lm55)#if one model does not converge, set to list[0]
 
 for(l in 1:5){
    
   lm<- lm_space1[[l]]
   
   if(length(lm)<10){
     all<- c(all,0)
     mLUIpp<- c(mLUIpp,0)
     dLUIpp<- c(dLUIpp,0)
     sharedpp<- c(sharedpp,0)
     cofm<- c(cofm,0)
     cofd<- c(cofd,0)
     int<- c(int,0)
     pvm<- c(pvm,1)
     pvd<- c(pvd,1)
     
     }else{
     suppressWarnings(lm1<- update(lm,~.-mLUI))#missing convergence
     suppressWarnings(lm2<- update(lm,~.-dLUI))
     suppressWarnings(lm3<- update(lm1,~.-dLUI))
     
     
     a<- r.squaredGLMM(lm)[1]#marginal R2 of the fixed effects
     m<- (r.squaredGLMM(lm)[1]-r.squaredGLMM(lm1)[1])#variance expl by mLUI
     d<- (r.squaredGLMM(lm)[1]-r.squaredGLMM(lm2)[1])#variance expl by dLUI
     sh<-(r.squaredGLMM(lm)[1]-r.squaredGLMM(lm3)[1]) - m - d
     
     pp <- 100*(c(m, d, sh)/a)
     
     all<- c(all, a)
     mLUIpp<- c(mLUIpp,as.numeric(pp[1]))
     dLUIpp<- c(dLUIpp,as.numeric(pp[2]))
     sharedpp<- c(sharedpp, as.numeric(pp[3]))
     cofm<- c(cofm,lm$coefficients$fixed[2])
     cofd<- c(cofd,lm$coefficients$fixed[3])
     int<- c(int, lm$coefficients$fixed[1])
     pvm<- c(pvm, round(Anova(lm, type="II")$`Pr(>Chisq)`[1], digits=3))
     pvd<- c(pvd, round(Anova(lm, type="II")$`Pr(>Chisq)`[2], digits=3))
     }
   
  }
EV_space_alpha[[k]]<- data.frame(all, mLUIpp, dLUIpp, sharedpp, cofm, cofd, int, pvm, pvd)
lm_space_alpha[[k]]<- lm_space1
 # }

names(EV_space_alpha)<- tlorder
names(lm_space_alpha)<- tlorder

save(EV_space_alpha, file="./data/OutputData/EV_space_LUI_alpha.RData")
save(lm_space_alpha, file="./data/OutputData/lm_space_LUI_alpha.RData")


#time
lm_time_alpha<- list()
EV_time_alpha<-list()
for(k in 1:3){ #I ran this loop manually since sometimes convergence problems make the loop stop
  all<- c()
  mLUIpp<- c()
  dLUIpp<- c()
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
  

  lm11<- lme(a0~mLUI+dLUI+pH+G500+soilPCA+dYR,
             data=sub, random=~1|EP)
  lm22<- lme(a1~mLUI+dLUI+pH+G500+soilPCA+dYR,
             data=sub, random=~1|EP)
  lm33<- lme(a2~mLUI+dLUI+pH+G500+soilPCA+dYR,
             data=sub, random=~1|EP)
  lm44<- lme(a3~mLUI+dLUI+pH+G500+soilPCA+dYR,
             data=sub, random=~1|EP)
  lm55<- lme(a4~mLUI+dLUI+pH+G500+soilPCA+dYR,
             data=sub, random=~1|EP)
  lm_time1<- list(lm11, lm22, lm33, lm44, lm55)
  
  #beta 0,1,2,3,4, sim
  for(l in 1:5){
    lm<- lm_time1[[l]]
    lm1<- update(lm,~.-mLUI)
    lm2<- update(lm,~.-dLUI)
    lm3<- update(lm1,~.-dLUI)
    
    
    a<- r.squaredGLMM(lm)[1]#marginal R2 of the fixed effects
    m<- (r.squaredGLMM(lm)[1]-r.squaredGLMM(lm1)[1])#variance expl by mLUI
    d<- (r.squaredGLMM(lm)[1]-r.squaredGLMM(lm2)[1])#variance expl by dLUI
    sh<-(r.squaredGLMM(lm)[1]-r.squaredGLMM(lm3)[1]) - m - d
    
    pp <- 100*(c(m, d, sh)/a)
    
    
    
    all<- c(all, a)
    mLUIpp<- c(mLUIpp,as.numeric(pp[1]))
    dLUIpp<- c(dLUIpp,as.numeric(pp[2]))
    sharedpp<- c(sharedpp, as.numeric(pp[3]))
    cofm<- c(cofm,lm$coefficients$fixed[2])
    cofd<- c(cofd,lm$coefficients$fixed[3])
    int<- c(int, lm$coefficients$fixed[1])
    pvm<- c(pvm, round(Anova(lm, type="II")$`Pr(>Chisq)`[1], digits=3))
    pvd<- c(pvd, round(Anova(lm, type="II")$`Pr(>Chisq)`[2], digits=3))
  }
  
  EV_time_alpha[[k]]<- data.frame(all, mLUIpp, dLUIpp, sharedpp, cofm, cofd, int, pvm, pvd)
  lm_time_alpha[[k]]<- lm_time1
}

names(EV_time_alpha)<- tlorder
names(lm_time_alpha)<- tlorder

save(EV_time_alpha, file="./data/OutputData/EV_time_LUI_alpha.RData")
save(lm_time_alpha, file="./data/OutputData/lm_time_LUI_alpha.RData")


#set working directory to folder "4. Plotting results"
save(EV_time_alpha, file="./data/InputData/EV_time_LUI_alpha.RData")
