#load data - too large to load all at once
#load("./pwise_space_plants_new.RData")
#load("./pwise_time_ayrs_new.RData")
load("./sitepred_pwise_space_step3.RData")#new pwise_plants_space with updated GRA residuals
load("./sitepred_pwise_time_step3.RData")#new pwise_plants_space with updated GRA residuals

pwise_time_plants<-as.data.frame(pwise_time) #new names as RData files
pwise_space_plants<-as.data.frame(pwise_space) #new names as RData files

rm(pwise_space)
rm(pwise_time)

#insects
#NB - both plant and insects files have the same names--> add org name to files
#load("./pwise_plants_space_in_new.RData")
#load("./pwise_time_ayrs_in_new.RData")
load("./sitepred_pwise_space_in_step3.RData")#new pwise_plants_space with updated GRA residuals
load("./sitepred_pwise_time_in_step3.RData")#new pwise_plants_space with updated GRA residuals

pwise_time_insects<-as.data.frame(pwise_time) #new names as RData files
pwise_space_insects<-as.data.frame(pwise_space) #new names as RData files

rm(pwise_space)
rm(pwise_time)

# # # # # # # # # # # # # #
# 0 - CALC GEODIST              ----

#geo dist for the insect datasets - DONE
a2<- ((pwise_space_insects$HWG2-pwise_space_insects$HWG1))^2
b2<- ((pwise_space_insects$RWG2-pwise_space_insects$RWG1))^2

pwise_space_insects$geo.dist<- sqrt((a2+b2))

a2<- ((pwise_time_insects$HWG2-pwise_time_insects$HWG1))^2
b2<- ((pwise_time_insects$RWG2-pwise_time_insects$RWG1))^2

pwise_time_insects$geo.dist<- sqrt((a2+b2))

save(pwise_space_insects, file="pwise_space_insects_new.RData")
save(pwise_time_insects, file="pwise_time_insects_new.RData")


#coordinates -not needed anymore
pwise_time_insects$HWG1<- pwise_time_insects$HWG
pwise_time_insects$RWG1<- pwise_time_insects$RWG

colnames(pwise_time_insects)[8]<- "HWG"
colnames(pwise_time_insects)[8]

colnames(pwise_time_insects)[10]<- "RWG"
colnames(pwise_time_insects)[10]

RWG<- c()

for(i in 1:8250) {
  rwg<-pwise_space_plants$RWG1[pwise_space_plants$EP1 %in% pwise_time_plants$EP[i]]
  RWG<- c(RWG, rwg[1])}
RWG[1:57]
pwise_time_plants$RWG<- RWG
colnames(pwise_time_plants)[3]<- "HWG"


# # # # # # # # # # # # # #
# 0 - CALC PAIRWISE MEAN              ----
#meanGRA as residuals normalised for EP and year differences
#temporal data - levelled for plot differences
#plants

pwise_time_plants$mGRA<- (pwise_time_plants$GRA1res+pwise_time_plants$GRA2res)/2 

#insects
pwise_time_insects$mGRA<- (pwise_time_insects$GRA1res+pwise_time_insects$GRA2res)/2 


#spatial data - levelled for year differences
#plants
pwise_space_plants$mGRA<- (pwise_space_plants$GRA1res+pwise_space_plants$GRA2res)/2 

#insects
pwise_space_insects$mGRA<- (pwise_space_insects$GRA1res+pwise_space_insects$GRA2res)/2 


save(pwise_time_insects, file="pwise_time_insects_new.RData")
save(pwise_time_plants, file="pwise_time_plants_new.RData")


#data prep - check whether all parameters are there
#soil parameters, temporal/geographic distance and isolation

colnames(pwise_space_plants)#G500_1, G500_2, pH_1, pH_2, soilPCA_1, soilPCA_2, msoil_PC1, dsoil_PC1, msur, dsur, mph, dph
colnames(pwise_space_insects)#G500_1, G500_2, pH_1, pH_2, soilPCA_1, soilPCA_2, msoil_PC1, dsoil_PC1, msur, dsur, mph, dph


colnames(pwise_time_plants)#"pH" "soilPCA" "G500"  
colnames(pwise_time_insects)#"pH" "soilPCA"   "G500"    


# # # # # # # # # # # # # #
# 0 - CALC STANDARDISING              ----
#standardize alpha diversities to range between 0 and 1, like for GDM analysis
#check before if not already saved in the file!!! 
#Should be corrected and saved 06.01.2022
#TODO check this! needs clean

library(vegan)
#plants
pwise_space_plants$da0st<- decostand(pwise_space_plants$da0, method="range", na.rm=T)
pwise_space_plants$da1st<- decostand(pwise_space_plants$da1, method="range", na.rm=T)
pwise_space_plants$da2st<- decostand(pwise_space_plants$da2, method="range", na.rm=T)
pwise_space_plants$da3st<- decostand(pwise_space_plants$da3, method="range", na.rm=T)
pwise_space_plants$da4st<- decostand(pwise_space_plants$da4, method="range", na.rm=T)

pwise_time_plants$da0st<- decostand(pwise_time_plants$da0abs, method="range", na.rm=T)
pwise_time_plants$da1st<- decostand(pwise_time_plants$da1abs, method="range", na.rm=T)
pwise_time_plants$da2st<- decostand(pwise_time_plants$da2abs, method="range", na.rm=T)
pwise_time_plants$da3st<- decostand(pwise_time_plants$da3abs, method="range", na.rm=T)
pwise_time_plants$da4st<- decostand(pwise_time_plants$da4abs, method="range", na.rm=T)

#herbivores
pwise_space_insects$dha0st<- decostand(pwise_space_insects$dha0, method="range", na.rm=T)
pwise_space_insects$dha1st<- decostand(pwise_space_insects$dha1, method="range", na.rm=T)
pwise_space_insects$dha2st<- decostand(pwise_space_insects$dha2, method="range", na.rm=T)
pwise_space_insects$dha3st<- decostand(pwise_space_insects$dha3, method="range", na.rm=T)
pwise_space_insects$dha4st<- decostand(pwise_space_insects$dha4, method="range", na.rm=T)

pwise_time_insects$dha0st<- decostand(pwise_time_insects$dha0abs, method="range", na.rm=T)
pwise_time_insects$dha1st<- decostand(pwise_time_insects$dha1abs, method="range", na.rm=T)
pwise_time_insects$dha2st<- decostand(pwise_time_insects$dha2abs, method="range", na.rm=T)
pwise_time_insects$dha3st<- decostand(pwise_time_insects$dha3abs, method="range", na.rm=T)
pwise_time_insects$dha4st<- decostand(pwise_time_insects$dha4abs, method="range", na.rm=T)

#sec. consumers
pwise_space_insects$dpa0st<- decostand(pwise_space_insects$dpa0, method="range", na.rm=T)
pwise_space_insects$dpa1st<- decostand(pwise_space_insects$dpa1, method="range", na.rm=T)
pwise_space_insects$dpa2st<- decostand(pwise_space_insects$dpa2, method="range", na.rm=T)
pwise_space_insects$dpa3st<- decostand(pwise_space_insects$dpa3, method="range", na.rm=T)
pwise_space_insects$dpa4st<- decostand(pwise_space_insects$dpa4, method="range", na.rm=T)

pwise_time_insects$dpa0st<- decostand(pwise_time_insects$dpa0abs, method="range", na.rm=T)
pwise_time_insects$dpa1st<- decostand(pwise_time_insects$dpa1abs, method="range", na.rm=T)
pwise_time_insects$dpa2st<- decostand(pwise_time_insects$dpa2abs, method="range", na.rm=T)
pwise_time_insects$dpa3st<- decostand(pwise_time_insects$dpa3abs, method="range", na.rm=T)
pwise_time_insects$dpa4st<- decostand(pwise_time_insects$dpa4abs, method="range", na.rm=T)

#Packages
library(lme4)
library(nlme)
library(lmPerm)
library(MuMIn)#r.quared.GLMM
library(car)#Anova(type="II")

tlorder <- c("plants", "herbivores", "secondary.consumers")



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

save(scaled_space_alpha, file="scaled_alpha_space_GRA.RData")



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

save(scaled_time_alpha, file="scaled_alpha_time_GRA.RData")





# # # # # # # # # # # # # #
# 1 - LINEAR MODELS              ----
#
load("./scaled_alpha_space_GRA.RData")
load("scaled_alpha_time_GRA.RData")

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

save(EV_space_alpha, file="EV_space_GRA_alpha.RData")
save(lm_space_alpha, file="lm_space_GRA_alpha.RData")


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

save(EV_time_alpha, file="EV_time_GRA_alpha.RData")
save(lm_time_alpha, file="lm_time_GRA_alpha.RData")



# # # # #
# 2 - DEPREC PLOTTING  ----
# 
#use Figure.R script for plotting
#TODO seems that below part is done in other script.

load("./EV_time_GRA_alpha.RData")
load("./EV_space_GRA_alpha.RData")

EV_time<- EV_time_alpha
EV_space<- EV_space_alpha

# # # # #
# 2.a. - DEPREC time plotting ----
# 
#EV
a0_time_mGRA_EV<-c(EV_time[[1]]$mGRApp[1],EV_time[[2]]$mGRApp[1],EV_time[[3]]$mGRApp[1])
a0_time_dGRA_EV<-c(EV_time[[1]]$dGRApp[1],EV_time[[2]]$dGRApp[1],EV_time[[3]]$dGRApp[1])

a1_time_mGRA_EV<-c(EV_time[[1]]$mGRApp[2],EV_time[[2]]$mGRApp[2],EV_time[[3]]$mGRApp[2])
a1_time_dGRA_EV<-c(EV_time[[1]]$dGRApp[2],EV_time[[2]]$dGRApp[2],EV_time[[3]]$dGRApp[2])

a2_time_mGRA_EV<-c(EV_time[[1]]$mGRApp[3],EV_time[[2]]$mGRApp[3],EV_time[[3]]$mGRApp[3])
a2_time_dGRA_EV<-c(EV_time[[1]]$dGRApp[3],EV_time[[2]]$dGRApp[3],EV_time[[3]]$dGRApp[3])

a3_time_mGRA_EV<-c(EV_time[[1]]$mGRApp[4],EV_time[[2]]$mGRApp[4],EV_time[[3]]$mGRApp[4])
a3_time_dGRA_EV<-c(EV_time[[1]]$dGRApp[4],EV_time[[2]]$dGRApp[4],EV_time[[3]]$dGRApp[4])

a4_time_mGRA_EV<-c(EV_time[[1]]$mGRApp[5],EV_time[[2]]$mGRApp[5],EV_time[[3]]$mGRApp[5])
a4_time_dGRA_EV<-c(EV_time[[1]]$dGRApp[5],EV_time[[2]]$dGRApp[5],EV_time[[3]]$dGRApp[5])

#coefficients
a0_time_mGRA<-c(EV_time[[1]]$cofm[1],EV_time[[2]]$cofm[1],EV_time[[3]]$cofm[1])
a0_time_dGRA<-c(EV_time[[1]]$cofd[1],EV_time[[2]]$cofd[1],EV_time[[3]]$cofd[1])

a1_time_mGRA<-c(EV_time[[1]]$cofm[2],EV_time[[2]]$cofm[2],EV_time[[3]]$cofm[2])
a1_time_dGRA<-c(EV_time[[1]]$cofd[2],EV_time[[2]]$cofd[2],EV_time[[3]]$cofd[2])

a2_time_mGRA<-c(EV_time[[1]]$cofm[3],EV_time[[2]]$cofm[3],EV_time[[3]]$cofm[3])
a2_time_dGRA<-c(EV_time[[1]]$cofd[3],EV_time[[2]]$cofd[3],EV_time[[3]]$cofd[3])

a3_time_mGRA<-c(EV_time[[1]]$cofm[4],EV_time[[2]]$cofm[4],EV_time[[3]]$cofm[4])
a3_time_dGRA<-c(EV_time[[1]]$cofd[4],EV_time[[2]]$cofd[4],EV_time[[3]]$cofd[4])

a4_time_mGRA<-c(EV_time[[1]]$cofm[5],EV_time[[2]]$cofm[5],EV_time[[3]]$cofm[5])
a4_time_dGRA<-c(EV_time[[1]]$cofd[5],EV_time[[2]]$cofd[5],EV_time[[3]]$cofd[5])

# # # # #
# 2.b. - DEPREC space plotting ----
# 
#EV
a0_space_mGRA_EV<-c(EV_space[[1]]$mGRApp[1],EV_space[[2]]$mGRApp[1],EV_space[[3]]$mGRApp[1])
a0_space_dGRA_EV<-c(EV_space[[1]]$dGRApp[1],EV_space[[2]]$dGRApp[1],EV_space[[3]]$dGRApp[1])

a1_space_mGRA_EV<-c(EV_space[[1]]$mGRApp[2],EV_space[[2]]$mGRApp[2],EV_space[[3]]$mGRApp[2])
a1_space_dGRA_EV<-c(EV_space[[1]]$dGRApp[2],EV_space[[2]]$dGRApp[2],EV_space[[3]]$dGRApp[2])

a2_space_mGRA_EV<-c(EV_space[[1]]$mGRApp[3],EV_space[[2]]$mGRApp[3],EV_space[[3]]$mGRApp[3])
a2_space_dGRA_EV<-c(EV_space[[1]]$dGRApp[3],EV_space[[2]]$dGRApp[3],EV_space[[3]]$dGRApp[3])

a3_space_mGRA_EV<-c(EV_space[[1]]$mGRApp[4],EV_space[[2]]$mGRApp[4],EV_space[[3]]$mGRApp[4])
a3_space_dGRA_EV<-c(EV_space[[1]]$dGRApp[4],EV_space[[2]]$dGRApp[4],EV_space[[3]]$dGRApp[4])

a4_space_mGRA_EV<-c(EV_space[[1]]$mGRApp[5],EV_space[[2]]$mGRApp[5],EV_space[[3]]$mGRApp[5])
a4_space_dGRA_EV<-c(EV_space[[1]]$dGRApp[5],EV_space[[2]]$dGRApp[5],EV_space[[3]]$dGRApp[5])

#coefficients
a0_space_mGRA<-c(EV_space[[1]]$cofm[1],EV_space[[2]]$cofm[1],EV_space[[3]]$cofm[1])
a0_space_dGRA<-c(EV_space[[1]]$cofd[1],EV_space[[2]]$cofd[1],EV_space[[3]]$cofd[1])

a1_space_mGRA<-c(EV_space[[1]]$cofm[2],EV_space[[2]]$cofm[2],EV_space[[3]]$cofm[2])
a1_space_dGRA<-c(EV_space[[1]]$cofd[2],EV_space[[2]]$cofd[2],EV_space[[3]]$cofd[2])

a2_space_mGRA<-c(EV_space[[1]]$cofm[3],EV_space[[2]]$cofm[3],EV_space[[3]]$cofm[3])
a2_space_dGRA<-c(EV_space[[1]]$cofd[3],EV_space[[2]]$cofd[3],EV_space[[3]]$cofd[3])

a3_space_mGRA<-c(EV_space[[1]]$cofm[4],EV_space[[2]]$cofm[4],EV_space[[3]]$cofm[4])
a3_space_dGRA<-c(EV_space[[1]]$cofd[4],EV_space[[2]]$cofd[4],EV_space[[3]]$cofd[4])

a4_space_mGRA<-c(EV_space[[1]]$cofm[5],EV_space[[2]]$cofm[5],EV_space[[3]]$cofm[5])
a4_space_dGRA<-c(EV_space[[1]]$cofd[5],EV_space[[2]]$cofd[5],EV_space[[3]]$cofd[5])

tlnames <- c("plants","herbivores","secondary consumers")

cols <- c("#1b9e77", "#d95f02", "#7570b3")
letmat <- matrix(letters[1:9], nrow=3,byrow=T)



# # # # #
# 2.x. - DEPREC a0 time ----
# 
#a0
#time
par(mfrow=c(1,2))
par(mar=c(2,15,2,2))
par(oma=c(8,6,3,0))

y <- barplot(rev(a0_time_mGRA), col = rev(cols),las=1,
             horiz=T, cex.main=2.5, cex.axis = 2.5, cex.names=2.5, 
             axisnames = T, names=rev(tlnames),
             xlim=c(-0.4, 0.4),
             main="mean GRA")
text(x=0.3, y=y, cex=2.5,
     labels= paste(round(rev(abs(a0_time_mGRA_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

z <- barplot(rev(a0_time_dGRA), col = rev(cols),las=1,
             horiz=T, cex.main=2.5, cex.axis = 2.5, cex.names=2.5, 
             axisnames = T,
             xlim=c(-0.4, 0.4),
             main="delta GRA")

text(x=-0.3, y=z,cex=2.5, 
     labels= paste(round(rev(abs(a0_time_dGRA_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

## combined x axis
text(-0.7,-1, expression(paste("Effect on temporal differences in ", 
                                alpha," diversity (richness, q = 0)",sep="")),
     xpd=NA,cex=2.5)


# # # # #
# 2.x. - DEPREC a0 space ----
# 
#space
par(mfrow=c(1,2))
par(mar=c(2,15,2,2))
par(oma=c(8,6,3,0))


y <- barplot(rev(a0_space_mGRA), col = rev(cols),las=1,
             horiz=T, cex.main=2.5, cex.axis = 2.5, cex.names=2.5,
             axisnames = T,names=rev(tlnames),
             xlim=c(-0.4,0.4),
             main="mean GRA")
text(x=0.3, y=y, cex=2.5,
     labels= paste(round(rev(abs(a0_space_mGRA_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))


z <- barplot(rev(a0_space_dGRA), col = rev(cols),las=1,
             horiz=T, cex.main=2.5,cex.axis = 2.5, cex.names=2.5,
             axisnames = T,
             xlim=c(-0.4,0.4),
             main="delta GRA")
text(x=-0.3, y=y, cex=2.5,
     labels= paste(round(rev(abs(a0_space_dGRA_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

## combined x axis
text(-0.65,-1,expression(paste("Effect on spatial differences in ",
                               alpha, " diversity (richness, q = 0)",sep="")),
     xpd=NA,cex=2.5)


# # # # #
# 2.x. - DEPREC a1 time ----
# 
#time
par(mfrow=c(1,2))
par(mar=c(2,15,2,2))
par(oma=c(8,6,3,0))

y <- barplot(rev(a1_time_mGRA), col = rev(cols),las=1,
             horiz=T, cex.main=2.5, cex.axis = 2.5, cex.names=2.5, 
             axisnames = T, names=rev(tlnames),
             xlim=c(-0.4, 0.4),
             main="mean GRA")
text(x=0.3, y=y, cex=2.5,
     labels= paste(round(rev(abs(a1_time_mGRA_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

z <- barplot(rev(a1_time_dGRA), col = rev(cols),las=1,
             horiz=T, cex.main=2.5, cex.axis = 2.5, cex.names=2.5, 
             axisnames = T,
             xlim=c(-0.4,0.4),
             main="delta GRA")

text(x=-0.3, y=z,cex=2.5, 
     labels= paste(round(rev(abs(a1_time_dGRA_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

## combined x axis
text(-0.65,-1, expression(paste("Effect on temporal differences in ", 
                                alpha," diversity (q = 1)",sep="")),
     xpd=NA,cex=2.5)


# # # # #
# 2.x. - DEPREC a1 space ----
# 
par(mfrow=c(1,2))
par(mar=c(2,15,2,2))
par(oma=c(8,6,3,0))


y <- barplot(rev(a1_space_mGRA), col = rev(cols),las=1,
             horiz=T, cex.main=2.5, cex.axis = 2.5, cex.names=2.5,
             axisnames = T, names=rev(tlnames),
             xlim=c(-0.4,0.4),
             main="mean GRA")
text(x=0.3, y=y, cex=2.5,
     labels= paste(round(rev(abs(a1_space_mGRA_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))


z <- barplot(rev(a1_space_dGRA), col = rev(cols),las=1,
             horiz=T, cex.main=2.5,cex.axis = 2.5, cex.names=2.5,
             axisnames = T,
             xlim=c(-0.4,0.4),
             main="delta GRA")
text(x=-0.3, y=y, cex=2.5,
     labels= paste(round(rev(abs(a1_space_dGRA_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

## combined x axis
text(-0.65,-1.0,expression(paste("Effect on spatial differences in ",
                                 alpha, " diversity (q = 1)",sep="")),
     xpd=NA,cex=2.5)


# # # # #
# 2.x. - DEPREC a2 time ----
# 
# a2
#time
par(mfrow=c(1,2))
par(mar=c(2,15,2,2))
par(oma=c(8,6,3,0))

y <- barplot(rev(a2_time_mGRA), col = rev(cols),las=1,
             horiz=T, cex.main=2.5, cex.axis = 2.5, cex.names=2.5, 
             axisnames = T, names=rev(tlnames),
             xlim=c(-0.4, 0.4),
             main="mean GRA")
text(x=0.3, y=y, cex=2.5,
     labels= paste(round(rev(abs(a2_time_mGRA_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

z <- barplot(rev(a2_time_dGRA), col = rev(cols),las=1,
             horiz=T, cex.main=2.5, cex.axis = 2.5, cex.names=2.5, 
             axisnames = T,
             xlim=c(-0.4,0.4),
             main="delta GRA")

text(x=-0.3, y=z,cex=2.5, 
     labels= paste(round(rev(abs(a2_time_dGRA_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

## combined x axis
text(-0.65,-1, expression(paste("Effect on temporal differences in ", 
                                alpha," diversity (q = 2)",sep="")),
     xpd=NA,cex=2.5)

# # # # #
# 2.x. - DEPREC a2 space ----
# 
#space
par(mfrow=c(1,2))
par(mar=c(2,15,2,2))
par(oma=c(8,6,3,0))


y <- barplot(rev(a2_space_mGRA), col = rev(cols),las=1,
             horiz=T, cex.main=2.5, cex.axis = 2.5, cex.names=2.5,
             axisnames = T, names=rev(tlnames),
             xlim=c(-0.4,0.4),
             main="mean GRA")
text(x=0.3, y=y, cex=2.5,
     labels= paste(round(rev(abs(a2_space_mGRA_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))


z <- barplot(rev(a2_space_dGRA), col = rev(cols),las=1,
             horiz=T, cex.main=2.5,cex.axis = 2.5, cex.names=2.5,
             axisnames = T,
             xlim=c(-0.4,0.4),
             main="delta GRA")
text(x=-0.3, y=y, cex=2.5,
     labels= paste(round(rev(abs(a2_space_dGRA_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

## combined x axis
text(-0.65,-1.0,expression(paste("Effect on spatial differences in ",
                                 alpha, " diversity (q = 2)",sep="")),
     xpd=NA,cex=2.5)


# # # # #
# 2.x. - DEPREC a3 time ----
# 
# a3

#time
par(mfrow=c(1,2))
par(mar=c(2,15,2,2))
par(oma=c(8,6,3,0))

y <- barplot(rev(a3_time_mGRA), col = rev(cols),las=1,
             horiz=T, cex.main=2.5, cex.axis = 2.5, cex.names=2.5, 
             axisnames = T, names=rev(tlnames),
             xlim=c(-0.4, 0.4),
             main="mean GRA")
text(x=0.3, y=y, cex=2.5,
     labels= paste(round(rev(abs(a3_time_mGRA_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

z <- barplot(rev(a3_time_dGRA), col = rev(cols),las=1,
             horiz=T, cex.main=2.5, cex.axis = 2.5, cex.names=2.5, 
             axisnames = T,
             xlim=c(-0.4,0.4),
             main="delta GRA")

text(x=-0.3, y=z,cex=2.5, 
     labels= paste(round(rev(abs(a3_time_dGRA_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

## combined x axis
text(-0.65,-1, expression(paste("Effect on temporal differences in ", 
                                alpha," diversity (q = 3)",sep="")),
     xpd=NA,cex=2.5)

# # # # #
# 2.x. - DEPREC a3 space ----
# 
#space
par(mfrow=c(1,2))
par(mar=c(2,15,2,2))
par(oma=c(8,6,3,0))


y <- barplot(rev(a3_space_mGRA), col = rev(cols),las=1,
             horiz=T, cex.main=2.5, cex.axis = 2.5, cex.names=2.5,
             axisnames = T,names=rev(tlnames),
             xlim=c(-0.4,0.4),
             main="mean GRA")
text(x=0.3, y=y, cex=2.5,
     labels= paste(round(rev(abs(a3_space_mGRA_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))


z <- barplot(rev(a3_space_dGRA), col = rev(cols),las=1,
             horiz=T, cex.main=2.5,cex.axis = 2.5, cex.names=2.5,
             axisnames = T,
             xlim=c(-0.4,0.4),
             main="delta GRA")
text(x=-0.3, y=y, cex=2.5,
     labels= paste(round(rev(abs(a3_space_dGRA_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

## combined x axis
text(-0.65,-1.0,expression(paste("Effect on spatial differences in ",
                                 alpha, " diversity (q = 3)",sep="")),
     xpd=NA,cex=2.5)


# # # # #
# 2.x. - DEPREC a4 time ----
# 
# a4

#time
par(mfrow=c(1,2))
par(mar=c(2,15,2,2))
par(oma=c(8,6,3,0))

y <- barplot(rev(a4_time_mGRA), col = rev(cols),las=1,
             horiz=T, cex.main=2.5, cex.axis = 2.5, cex.names=2.5, 
             axisnames = T,names=rev(tlnames),
             xlim=c(-0.4, 0.4),
             main="mean GRA")
text(x=0.3, y=y, cex=2.5,
     labels= paste(round(rev(abs(a4_time_mGRA_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

z <- barplot(rev(a4_time_dGRA), col = rev(cols),las=1,
             horiz=T, cex.main=2.5, cex.axis = 2.5, cex.names=2.5, 
             axisnames = T,
             xlim=c(-0.4,0.4),
             main="delta GRA")

text(x=-0.3, y=z,cex=2.5, 
     labels= paste(round(rev(abs(a4_time_dGRA_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

## combined x axis
text(-0.65,-1, expression(paste("Effect on temporal differences in ", 
                                alpha," diversity (q = 4)",sep="")),
     xpd=NA,cex=2.5)

# # # # #
# 2.x. - DEPREC a4 space ----
# 
#space
par(mfrow=c(1,2))
par(mar=c(2,15,2,2))
par(oma=c(8,6,3,0))


y <- barplot(rev(a4_space_mGRA), col = rev(cols),las=1,
             horiz=T, cex.main=2.5, cex.axis = 2.5, cex.names=2.5,
             axisnames = T,names=rev(tlnames),
             xlim=c(-0.4,0.4),
             main="mean GRA")
text(x=0.3, y=y, cex=2.5,
     labels= paste(round(rev(abs(a4_space_mGRA_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))


z <- barplot(rev(a4_space_dGRA), col = rev(cols),las=1,
             horiz=T, cex.main=2.5,cex.axis = 2.5, cex.names=2.5,
             axisnames = T,
             xlim=c(-0.4,0.4),
             main="delta GRA")
text(x=-0.3, y=y, cex=2.5,
     labels= paste(round(rev(abs(a4_space_dGRA_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

## combined x axis
text(-0.65,-1.0,expression(paste("Effect on spatial differences in ",
                                 alpha, " diversity (q = 4)",sep="")),
     xpd=NA,cex=2.5)

