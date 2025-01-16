#load data - too large to load all at once
#load("./pwise_space_plants_new.RData")
#load("./pwise_time_ayrs_new.RData")
load("./sitepred_pwise_space_step3.RData")#new pwise_plants_space with updated MOW residuals
load("./sitepred_pwise_time_step3.RData")#new pwise_plants_space with updated MOW residuals

pwise_time_plants<-as.data.frame(pwise_time) #new names as RData files
pwise_space_plants<-as.data.frame(pwise_space) #new names as RData files

rm(pwise_space)
rm(pwise_time)

#insects
#NB - both plant and insects files have the same names--> add org name to files
#load("./pwise_plants_space_in_new.RData")
#load("./pwise_time_ayrs_in_new.RData")
load("./sitepred_pwise_space_in_step3.RData")#new pwise_plants_space with updated MOW residuals
load("./sitepred_pwise_time_in_step3.RData")#new pwise_plants_space with updated MOW residuals

pwise_time_insects<-as.data.frame(pwise_time) #new names as RData files
pwise_space_insects<-as.data.frame(pwise_space) #new names as RData files

rm(pwise_space)
rm(pwise_time)

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

#meanMOW as residuals normalised for EP and year differences
#temporal data - levelled for plot differences
#plants

pwise_time_plants$mMOW<- (pwise_time_plants$MOW1res+pwise_time_plants$MOW2res)/2 

#insects
pwise_time_insects$mMOW<- (pwise_time_insects$MOW1res+pwise_time_insects$MOW2res)/2 


#spatial data - levelled for year differences
#plants
pwise_space_plants$mMOW<- (pwise_space_plants$MOW1res+pwise_space_plants$MOW2res)/2 

#insects
pwise_space_insects$mMOW<- (pwise_space_insects$MOW1res+pwise_space_insects$MOW2res)/2 


save(pwise_time_insects, file="pwise_time_insects_new.RData")
save(pwise_time_plants, file="pwise_time_plants_new.RData")
#data prep - check whether all parameters are there
#soil parameters, temporal/geographic distance and isolation

colnames(pwise_space_plants)#G500_1, G500_2, pH_1, pH_2, soilPCA_1, soilPCA_2, msoil_PC1, dsoil_PC1, msur, dsur, mph, dph
colnames(pwise_space_insects)#G500_1, G500_2, pH_1, pH_2, soilPCA_1, soilPCA_2, msoil_PC1, dsoil_PC1, msur, dsur, mph, dph


colnames(pwise_time_plants)#"pH" "soilPCA" "G500"  
colnames(pwise_time_insects)#"pH" "soilPCA"   "G500"    

#standardize alpha diversities to range between 0 and 1, like for GDM analysis
#check before if not already saved in the file!!! 
#Should be corrected and saved 06.01.2022

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
#################PREPARATION####################
###############################################
#scaling all differences and means to max
####space
#plants
vars<- c("da0st", "da1st", "da2st", "da3st", "da4st", 
         "mMOW", "dMOW", "geo.dist", "YR", 
         "msur", "dsur","mph", "dph","msoilPC1", "dsoilPC1")
sub_plants_space<- pwise_space_plants[,colnames(pwise_space_plants) %in% vars]
sub_plants_space<- sub_plants_space[complete.cases(sub_plants_space),]

colnames(sub_plants_space)
preds<- apply(sub_plants_space[,c(2:10)], 2, function(x)x/max(x, na.rm=T))
sub_plants_space<- cbind(sub_plants_space[c(11:15, 1)], preds)
colnames(sub_plants_space)


#sec. consumers
vars<- c("dpa0st", "dpa1st", "dpa2st", "dpa3st", "dpa4st",
         "mMOW", "dMOW", "geo.dist", "YR",
         "msur", "dsur","mph", "dph","msoilPC1", "dsoilPC1")
sub_pred_space<- pwise_space_insects[,colnames(pwise_space_insects) %in% vars]
sub_pred_space<- sub_pred_space[complete.cases(sub_pred_space),]

colnames(sub_pred_space)
preds<- apply(sub_pred_space[,c(2:10)], 2, function(x)x/max(x, na.rm=T))
sub_pred_space<- cbind(sub_pred_space[c(11:15, 1)], preds)
colnames(sub_pred_space)

#herbivores
vars<- c("dha0st", "dha1st", "dha2st", "dha3st","dha4st", 
         "mMOW", "dMOW", "geo.dist", "YR",
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

save(scaled_space_alpha, file="scaled_alpha_space_MOW.RData")


####time####
###########
#plants
vars<- c("da0st", "da1st", "da2st", "da3st","da4st", 
         "mMOW", "dMOW", "dYR", "HWG", "RWG",
         "G500", "pH","soilPCA", "EP")
sub_plants_time<- pwise_time_plants[,colnames(pwise_time_plants) %in% vars]
sub_plants_time<- sub_plants_time[complete.cases(sub_plants_time),]

colnames(sub_plants_time)
preds<- apply(sub_plants_time[,c(2:9)], 2, function(x)x/max(x, na.rm=T))
sub_plants_time<- cbind(sub_plants_time[c(10:14,1)], preds)
colnames(sub_plants_time)

#sec. consumers
vars<- c("dpa0st", "dpa1st", "dpa2st", "dpa3st", "dpa4st",
         "mMOW", "dMOW", "dYR", "HWG", "RWG",
         "G500", "pH","soilPCA", "EP")
sub_pred_time<- pwise_time_insects[,colnames(pwise_time_insects) %in% vars]
sub_pred_time<- sub_pred_time[complete.cases(sub_pred_time),]

colnames(sub_pred_time)
preds<- apply(sub_pred_time[,c(2:9)], 2, function(x)x/max(x, na.rm=T))
sub_pred_time<- cbind(sub_pred_time[c(10:14,1)], preds)
colnames(sub_pred_time)


#herbivores
vars<- c("dha0st", "dha1st", "dha2st", "dha3st","dha4st", 
         "mMOW", "dMOW", "dYR", "HWG", "RWG",
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

save(scaled_time_alpha, file="scaled_alpha_time_MOW.RData")

########################
###LINEAR MODELS
#######################
load("./scaled_alpha_space_MOW.RData")
load("scaled_alpha_time_MOW.RData")

#space
EV_space_alpha<-list()
lm_space_alpha<- list()
#I ran this loop manually since sometimes convergence problems make the loop stop
for(k in 1:3){ 
all<- c()
mMOWpp<- c()
dMOWpp<- c()
sharedpp<- c()
cofm<- c()
cofd<- c()
int<- c()
pvm<- c()
pvd<- c()
lm_space1<- list()

sub<- scaled_space_alpha[[k]]
sub$YR<- as.factor(sub$YR)

names(sub)[1]<- "a0"
names(sub)[2]<- "a1"
names(sub)[3]<- "a2"
names(sub)[4]<- "a3"
names(sub)[5]<- "a4"


lm11<- lme(a0~mMOW+dMOW+dph+mph+msur+dsur+dsoilPC1+msoilPC1+geo.dist,
           data=sub, random=~1|YR)
lm22<- lme(a1~mMOW+dMOW+dph+mph+msur+dsur+dsoilPC1+msoilPC1+geo.dist,
           data=sub, random=~1|YR)
lm33<- lme(a2~mMOW+dMOW+dph+mph+msur+dsur+dsoilPC1+msoilPC1+geo.dist,
           data=sub, random=~1|YR)
lm44<- lme(a3~mMOW+dMOW+dph+mph+msur+dsur+dsoilPC1+msoilPC1+geo.dist,
           data=sub, random=~1|YR)
lm55<- lme(a4~mMOW+dMOW+dph+mph+msur+dsur+dsoilPC1+msoilPC1+geo.dist,
           data=sub, random=~1|YR)

lm_space1<- list(lm11,lm22, lm33, lm44, lm55)#if one model does not converge, set to list[0]

for(l in 1:5){
  
  lm<- lm_space1[[l]]
  
  if(length(lm)<10){
    all<- c(all,0)
    mMOWpp<- c(mMOWpp,0)
    dMOWpp<- c(dMOWpp,0)
    sharedpp<- c(sharedpp,0)
    cofm<- c(cofm,0)
    cofd<- c(cofd,0)
    int<- c(int,0)
    pvm<- c(pvm,1)
    pvd<- c(pvd,1)
    
  }else{
    suppressWarnings(lm1<- update(lm,~.-mMOW))#missing convergence
    suppressWarnings(lm2<- update(lm,~.-dMOW))
    suppressWarnings(lm3<- update(lm1,~.-dMOW))
    
    
    a<- r.squaredGLMM(lm)[1]#marginal R2 of the fixed effects
    m<- (r.squaredGLMM(lm)[1]-r.squaredGLMM(lm1)[1])#variance expl by mMOW
    d<- (r.squaredGLMM(lm)[1]-r.squaredGLMM(lm2)[1])#variance expl by dMOW
    sh<-(r.squaredGLMM(lm)[1]-r.squaredGLMM(lm3)[1]) - m - d
    
    pp <- 100*(c(m, d, sh)/a)
    
    all<- c(all, a)
    mMOWpp<- c(mMOWpp,as.numeric(pp[1]))
    dMOWpp<- c(dMOWpp,as.numeric(pp[2]))
    sharedpp<- c(sharedpp, as.numeric(pp[3]))
    cofm<- c(cofm,lm$coefficients$fixed[2])
    cofd<- c(cofd,lm$coefficients$fixed[3])
    int<- c(int, lm$coefficients$fixed[1])
    pvm<- c(pvm, round(Anova(lm, type="II")$`Pr(>Chisq)`[1], digits=3))
    pvd<- c(pvd, round(Anova(lm, type="II")$`Pr(>Chisq)`[2], digits=3))
  }
  
}
EV_space_alpha[[k]]<- data.frame(all, mMOWpp, dMOWpp, sharedpp, cofm, cofd, int, pvm, pvd)
lm_space_alpha[[k]]<- lm_space1
}

names(EV_space_alpha)<- tlorder
names(lm_space_alpha)<- tlorder

save(EV_space_alpha, file="EV_space_MOW_alpha.RData")
save(lm_space_alpha, file="lm_space_MOW_alpha.RData")


#time
lm_time_alpha<- list()
EV_time_alpha<-list()
for(k in 1:3){ #I ran this loop manually since sometimes convergence problems make the loop stop
  all<- c()
  mMOWpp<- c()
  dMOWpp<- c()
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
  
  
  lm11<- lme(a0~mMOW+dMOW+pH+G500+soilPCA+dYR,
             data=sub, random=~1|EP)
  lm22<- lme(a1~mMOW+dMOW+pH+G500+soilPCA+dYR,
             data=sub, random=~1|EP)
  lm33<- lme(a2~mMOW+dMOW+pH+G500+soilPCA+dYR,
             data=sub, random=~1|EP)
  lm44<- lme(a3~mMOW+dMOW+pH+G500+soilPCA+dYR,
             data=sub, random=~1|EP)
  lm55<- lme(a4~mMOW+dMOW+pH+G500+soilPCA+dYR,
             data=sub, random=~1|EP)
  lm_time1<- list(lm11, lm22, lm33, lm44, lm55)
  
  #beta 0,1,2,3,4, sim
  for(l in 1:5){
    lm<- lm_time1[[l]]
    lm1<- update(lm,~.-mMOW)
    lm2<- update(lm,~.-dMOW)
    lm3<- update(lm1,~.-dMOW)
    
    
    a<- r.squaredGLMM(lm)[1]#marginal R2 of the fixed effects
    m<- (r.squaredGLMM(lm)[1]-r.squaredGLMM(lm1)[1])#variance expl by mMOW
    d<- (r.squaredGLMM(lm)[1]-r.squaredGLMM(lm2)[1])#variance expl by dMOW
    sh<-(r.squaredGLMM(lm)[1]-r.squaredGLMM(lm3)[1]) - m - d
    
    pp <- 100*(c(m, d, sh)/a)
    
    
    
    all<- c(all, a)
    mMOWpp<- c(mMOWpp,as.numeric(pp[1]))
    dMOWpp<- c(dMOWpp,as.numeric(pp[2]))
    sharedpp<- c(sharedpp, as.numeric(pp[3]))
    cofm<- c(cofm,lm$coefficients$fixed[2])
    cofd<- c(cofd,lm$coefficients$fixed[3])
    int<- c(int, lm$coefficients$fixed[1])
    pvm<- c(pvm, round(Anova(lm, type="II")$`Pr(>Chisq)`[1], digits=3))
    pvd<- c(pvd, round(Anova(lm, type="II")$`Pr(>Chisq)`[2], digits=3))
  }
  
  EV_time_alpha[[k]]<- data.frame(all, mMOWpp, dMOWpp, sharedpp, cofm, cofd, int, pvm, pvd)
  lm_time_alpha[[k]]<- lm_time1
}

names(EV_time_alpha)<- tlorder
names(lm_time_alpha)<- tlorder

save(EV_time_alpha, file="EV_time_MOW_alpha.RData")
save(lm_time_alpha, file="lm_time_MOW_alpha.RData")

####PLOTTING
#use Figure.R script for plotting

load("./EV_time_MOW_alpha.RData")
load("./EV_space_MOW_alpha.RData")

EV_time<- EV_time_alpha
EV_space<- EV_space_alpha

#time
#EV
a0_time_mMOW_EV<-c(EV_time[[1]]$mMOWpp[1],EV_time[[2]]$mMOWpp[1],EV_time[[3]]$mMOWpp[1])
a0_time_dMOW_EV<-c(EV_time[[1]]$dMOWpp[1],EV_time[[2]]$dMOWpp[1],EV_time[[3]]$dMOWpp[1])


a1_time_mMOW_EV<-c(EV_time[[1]]$mMOWpp[2],EV_time[[2]]$mMOWpp[2],EV_time[[3]]$mMOWpp[2])
a1_time_dMOW_EV<-c(EV_time[[1]]$dMOWpp[2],EV_time[[2]]$dMOWpp[2],EV_time[[3]]$dMOWpp[2])

a2_time_mMOW_EV<-c(EV_time[[1]]$mMOWpp[3],EV_time[[2]]$mMOWpp[3],EV_time[[3]]$mMOWpp[3])
a2_time_dMOW_EV<-c(EV_time[[1]]$dMOWpp[3],EV_time[[2]]$dMOWpp[3],EV_time[[3]]$dMOWpp[3])

a3_time_mMOW_EV<-c(EV_time[[1]]$mMOWpp[4],EV_time[[2]]$mMOWpp[4],EV_time[[3]]$mMOWpp[4])
a3_time_dMOW_EV<-c(EV_time[[1]]$dMOWpp[4],EV_time[[2]]$dMOWpp[4],EV_time[[3]]$dMOWpp[4])

a4_time_mMOW_EV<-c(EV_time[[1]]$mMOWpp[5],EV_time[[2]]$mMOWpp[5],EV_time[[3]]$mMOWpp[5])
a4_time_dMOW_EV<-c(EV_time[[1]]$dMOWpp[5],EV_time[[2]]$dMOWpp[5],EV_time[[3]]$dMOWpp[5])

#coefficients
a0_time_mMOW<-c(EV_time[[1]]$cofm[1],EV_time[[2]]$cofm[1],EV_time[[3]]$cofm[1])
a0_time_dMOW<-c(EV_time[[1]]$cofd[1],EV_time[[2]]$cofd[1],EV_time[[3]]$cofd[1])

a1_time_mMOW<-c(EV_time[[1]]$cofm[2],EV_time[[2]]$cofm[2],EV_time[[3]]$cofm[2])
a1_time_dMOW<-c(EV_time[[1]]$cofd[2],EV_time[[2]]$cofd[2],EV_time[[3]]$cofd[2])

a2_time_mMOW<-c(EV_time[[1]]$cofm[3],EV_time[[2]]$cofm[3],EV_time[[3]]$cofm[3])
a2_time_dMOW<-c(EV_time[[1]]$cofd[3],EV_time[[2]]$cofd[3],EV_time[[3]]$cofd[3])

a3_time_mMOW<-c(EV_time[[1]]$cofm[4],EV_time[[2]]$cofm[4],EV_time[[3]]$cofm[4])
a3_time_dMOW<-c(EV_time[[1]]$cofd[4],EV_time[[2]]$cofd[4],EV_time[[3]]$cofd[4])

a4_time_mMOW<-c(EV_time[[1]]$cofm[5],EV_time[[2]]$cofm[5],EV_time[[3]]$cofm[5])
a4_time_dMOW<-c(EV_time[[1]]$cofd[5],EV_time[[2]]$cofd[5],EV_time[[3]]$cofd[5])

#space
#EV
a0_space_mMOW_EV<-c(EV_space[[1]]$mMOWpp[1],EV_space[[2]]$mMOWpp[1],EV_space[[3]]$mMOWpp[1])
a0_space_dMOW_EV<-c(EV_space[[1]]$dMOWpp[1],EV_space[[2]]$dMOWpp[1],EV_space[[3]]$dMOWpp[1])

a1_space_mMOW_EV<-c(EV_space[[1]]$mMOWpp[2],EV_space[[2]]$mMOWpp[2],EV_space[[3]]$mMOWpp[2])
a1_space_dMOW_EV<-c(EV_space[[1]]$dMOWpp[2],EV_space[[2]]$dMOWpp[2],EV_space[[3]]$dMOWpp[2])

a2_space_mMOW_EV<-c(EV_space[[1]]$mMOWpp[3],EV_space[[2]]$mMOWpp[3],EV_space[[3]]$mMOWpp[3])
a2_space_dMOW_EV<-c(EV_space[[1]]$dMOWpp[3],EV_space[[2]]$dMOWpp[3],EV_space[[3]]$dMOWpp[3])

a3_space_mMOW_EV<-c(EV_space[[1]]$mMOWpp[4],EV_space[[2]]$mMOWpp[4],EV_space[[3]]$mMOWpp[4])
a3_space_dMOW_EV<-c(EV_space[[1]]$dMOWpp[4],EV_space[[2]]$dMOWpp[4],EV_space[[3]]$dMOWpp[4])

a4_space_mMOW_EV<-c(EV_space[[1]]$mMOWpp[5],EV_space[[2]]$mMOWpp[5],EV_space[[3]]$mMOWpp[5])
a4_space_dMOW_EV<-c(EV_space[[1]]$dMOWpp[5],EV_space[[2]]$dMOWpp[5],EV_space[[3]]$dMOWpp[5])

#coefficients
a0_space_mMOW<-c(EV_space[[1]]$cofm[1],EV_space[[2]]$cofm[1],EV_space[[3]]$cofm[1])
a0_space_dMOW<-c(EV_space[[1]]$cofd[1],EV_space[[2]]$cofd[1],EV_space[[3]]$cofd[1])

a1_space_mMOW<-c(EV_space[[1]]$cofm[2],EV_space[[2]]$cofm[2],EV_space[[3]]$cofm[2])
a1_space_dMOW<-c(EV_space[[1]]$cofd[2],EV_space[[2]]$cofd[2],EV_space[[3]]$cofd[2])

a2_space_mMOW<-c(EV_space[[1]]$cofm[3],EV_space[[2]]$cofm[3],EV_space[[3]]$cofm[3])
a2_space_dMOW<-c(EV_space[[1]]$cofd[3],EV_space[[2]]$cofd[3],EV_space[[3]]$cofd[3])

a3_space_mMOW<-c(EV_space[[1]]$cofm[4],EV_space[[2]]$cofm[4],EV_space[[3]]$cofm[4])
a3_space_dMOW<-c(EV_space[[1]]$cofd[4],EV_space[[2]]$cofd[4],EV_space[[3]]$cofd[4])

a4_space_mMOW<-c(EV_space[[1]]$cofm[5],EV_space[[2]]$cofm[5],EV_space[[3]]$cofm[5])
a4_space_dMOW<-c(EV_space[[1]]$cofd[5],EV_space[[2]]$cofd[5],EV_space[[3]]$cofd[5])

tlnames <- c("plants","herbivores","secondary consumers")

cols <- c("#1b9e77", "#d95f02", "#7570b3")
letmat <- matrix(letters[1:9], nrow=3,byrow=T)

#a0
#time
par(mfrow=c(1,2))
par(mar=c(2,15,2,2))
par(oma=c(8,6,3,0))

y <- barplot(rev(a0_time_mMOW), col = rev(cols),las=1,
             horiz=T, cex.main=2.5, cex.axis = 2.5, cex.names=2.5, 
             axisnames = T, names=rev(tlnames),
             xlim=c(-0.07, 0.07),
             main="mean MOW")
text(x=0.06, y=y, cex=2.5,
     labels= paste(round(rev(abs(a0_time_mMOW_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

z <- barplot(rev(a0_time_dMOW), col = rev(cols),las=1,
             horiz=T, cex.main=2.5, cex.axis = 2.5, cex.names=2.5, 
             axisnames = T,
             xlim=c(-0.07,0.07),
             main="delta MOW")

text(x=-0.06, y=z,cex=2.5, 
     labels= paste(round(rev(abs(a0_time_dMOW_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

## combined x axis
text(-0.12,-1, expression(paste("Effect on temporal differences in ", 
                                alpha," diversity (richness, q = 0)",sep="")),
     xpd=NA,cex=2.5)

#space
par(mfrow=c(1,2))
par(mar=c(2,15,2,2))
par(oma=c(8,6,3,0))


y <- barplot(rev(a0_space_mMOW), col = rev(cols),las=1,
             horiz=T, cex.main=2.5, cex.axis = 2.5, cex.names=2.5,
             axisnames = T,names=rev(tlnames),
             xlim=c(-0.4,0.4),
             main="mean MOW")
text(x=0.3, y=y, cex=2.5,
     labels= paste(round(rev(abs(a0_space_mMOW_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))


z <- barplot(rev(a0_space_dMOW), col = rev(cols),las=1,
             horiz=T, cex.main=2.5,cex.axis = 2.5, cex.names=2.5,
             axisnames = T,
             xlim=c(-0.4,0.4),
             main="delta MOW")
text(x=-0.3, y=y, cex=2.5,
     labels= paste(round(rev(abs(a0_space_dMOW_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

## combined x axis
text(-0.65,-1,expression(paste("Effect on spatial differences in ",
                               alpha, " diversity (richness, q = 0)",sep="")),
     xpd=NA,cex=2.5)


#time
par(mfrow=c(1,2))
par(mar=c(2,15,2,2))
par(oma=c(8,6,3,0))

y <- barplot(rev(a1_time_mMOW), col = rev(cols),las=1,
             horiz=T, cex.main=2.5, cex.axis = 2.5, cex.names=2.5, 
             axisnames = T, names=rev(tlnames),
             xlim=c(-0.07, 0.07),
             main="mean MOW")
text(x=0.06, y=y, cex=2.5,
     labels= paste(round(rev(abs(a1_time_mMOW_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

z <- barplot(rev(a1_time_dMOW), col = rev(cols),las=1,
             horiz=T, cex.main=2.5, cex.axis = 2.5, cex.names=2.5, 
             axisnames = T,
             xlim=c(-0.07,0.07),
             main="delta MOW")

text(x=-0.06, y=z,cex=2.5, 
     labels= paste(round(rev(abs(a1_time_dMOW_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

## combined x axis
text(-0.12,-1, expression(paste("Effect on temporal differences in ", 
                                alpha," diversity (q = 1)",sep="")),
     xpd=NA,cex=2.5)

#space
par(mfrow=c(1,2))
par(mar=c(2,15,2,2))
par(oma=c(8,6,3,0))


y <- barplot(rev(a1_space_mMOW), col = rev(cols),las=1,
             horiz=T, cex.main=2.5, cex.axis = 2.5, cex.names=2.5,
             axisnames = T, names=rev(tlnames),
             xlim=c(-0.4,0.4),
             main="mean MOW")
text(x=0.3, y=y, cex=2.5,
     labels= paste(round(rev(abs(a1_space_mMOW_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))


z <- barplot(rev(a1_space_dMOW), col = rev(cols),las=1,
             horiz=T, cex.main=2.5,cex.axis = 2.5, cex.names=2.5,
             axisnames = T,
             xlim=c(-0.4,0.4),
             main="delta MOW")
text(x=-0.3, y=y, cex=2.5,
     labels= paste(round(rev(abs(a1_space_dMOW_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

## combined x axis
text(-0.65,-1.0,expression(paste("Effect on spatial differences in ",
                                 alpha, " diversity (q = 1)",sep="")),
     xpd=NA,cex=2.5)

# a2
#####

#time
par(mfrow=c(1,2))
par(mar=c(2,15,2,2))
par(oma=c(8,6,3,0))

y <- barplot(rev(a2_time_mMOW), col = rev(cols),las=1,
             horiz=T, cex.main=2.5, cex.axis = 2.5, cex.names=2.5, 
             axisnames = T, names=rev(tlnames),
             xlim=c(-0.07, 0.07),
             main="mean MOW")
text(x=0.06, y=y, cex=2.5,
     labels= paste(round(rev(abs(a2_time_mMOW_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

z <- barplot(rev(a2_time_dMOW), col = rev(cols),las=1,
             horiz=T, cex.main=2.5, cex.axis = 2.5, cex.names=2.5, 
             axisnames = T,
             xlim=c(-0.07,0.06),
             main="delta MOW")

text(x=-0.06, y=z,cex=2.5, 
     labels= paste(round(rev(abs(a2_time_dMOW_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

## combined x axis
text(-0.12,-1, expression(paste("Effect on temporal differences in ", 
                                alpha," diversity (q = 2)",sep="")),
     xpd=NA,cex=2.5)

#space
par(mfrow=c(1,2))
par(mar=c(2,15,2,2))
par(oma=c(8,6,3,0))


y <- barplot(rev(a2_space_mMOW), col = rev(cols),las=1,
             horiz=T, cex.main=2.5, cex.axis = 2.5, cex.names=2.5,
             axisnames = T, names=rev(tlnames),
             xlim=c(-0.4,0.4),
             main="mean MOW")
text(x=0.3, y=y, cex=2.5,
     labels= paste(round(rev(abs(a2_space_mMOW_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))


z <- barplot(rev(a2_space_dMOW), col = rev(cols),las=1,
             horiz=T, cex.main=2.5,cex.axis = 2.5, cex.names=2.5,
             axisnames = T,
             xlim=c(-0.4,0.4),
             main="delta MOW")
text(x=-0.3, y=y, cex=2.5,
     labels= paste(round(rev(abs(a2_space_dMOW_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

## combined x axis
text(-0.65,-1.0,expression(paste("Effect on spatial differences in ",
                                 alpha, " diversity (q = 2)",sep="")),
     xpd=NA,cex=2.5)

# a3
#####

#time
par(mfrow=c(1,2))
par(mar=c(2,15,2,2))
par(oma=c(8,6,3,0))

y <- barplot(rev(a3_time_mMOW), col = rev(cols),las=1,
             horiz=T, cex.main=2.5, cex.axis = 2.5, cex.names=2.5, 
             axisnames = T, names=rev(tlnames),
             xlim=c(-0.07, 0.07),
             main="mean MOW")
text(x=0.06, y=y, cex=2.5,
     labels= paste(round(rev(abs(a3_time_mMOW_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

z <- barplot(rev(a3_time_dMOW), col = rev(cols),las=1,
             horiz=T, cex.main=2.5, cex.axis = 2.5, cex.names=2.5, 
             axisnames = T,
             xlim=c(-0.07,0.07),
             main="delta MOW")

text(x=-0.06, y=z,cex=2.5, 
     labels= paste(round(rev(abs(a3_time_dMOW_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

## combined x axis
text(-0.12,-1, expression(paste("Effect on temporal differences in ", 
                                alpha," diversity (q = 3)",sep="")),
     xpd=NA,cex=2.5)

#space
par(mfrow=c(1,2))
par(mar=c(2,15,2,2))
par(oma=c(8,6,3,0))


y <- barplot(rev(a3_space_mMOW), col = rev(cols),las=1,
             horiz=T, cex.main=2.5, cex.axis = 2.5, cex.names=2.5,
             axisnames = T,names=rev(tlnames),
             xlim=c(-0.4,0.4),
             main="mean MOW")
text(x=0.3, y=y, cex=2.5,
     labels= paste(round(rev(abs(a3_space_mMOW_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))


z <- barplot(rev(a3_space_dMOW), col = rev(cols),las=1,
             horiz=T, cex.main=2.5,cex.axis = 2.5, cex.names=2.5,
             axisnames = T,
             xlim=c(-0.4,0.4),
             main="delta MOW")
text(x=-0.3, y=y, cex=2.5,
     labels= paste(round(rev(abs(a3_space_dMOW_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

## combined x axis
text(-0.65,-1.0,expression(paste("Effect on spatial differences in ",
                                 alpha, " diversity (q = 3)",sep="")),
     xpd=NA,cex=2.5)

# a4
#####

#time
par(mfrow=c(1,2))
par(mar=c(2,15,2,2))
par(oma=c(8,6,3,0))

y <- barplot(rev(a4_time_mMOW), col = rev(cols),las=1,
             horiz=T, cex.main=2.5, cex.axis = 2.5, cex.names=2.5, 
             axisnames = T,names=rev(tlnames),
             xlim=c(-0.07, 0.07),
             main="mean MOW")
text(x=0.06, y=y, cex=2.5,
     labels= paste(round(rev(abs(a4_time_mMOW_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

z <- barplot(rev(a4_time_dMOW), col = rev(cols),las=1,
             horiz=T, cex.main=2.5, cex.axis = 2.5, cex.names=2.5, 
             axisnames = T,
             xlim=c(-0.07,0.07),
             main="delta MOW")

text(x=-0.06, y=z,cex=2.5, 
     labels= paste(round(rev(abs(a4_time_dMOW_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

## combined x axis
text(-0.12,-1, expression(paste("Effect on temporal differences in ", 
                                alpha," diversity (q = 4)",sep="")),
     xpd=NA,cex=2.5)

#space
par(mfrow=c(1,2))
par(mar=c(2,15,2,2))
par(oma=c(8,6,3,0))


y <- barplot(rev(a4_space_mMOW), col = rev(cols),las=1,
             horiz=T, cex.main=2.5, cex.axis = 2.5, cex.names=2.5,
             axisnames = T,names=rev(tlnames),
             xlim=c(-0.4,0.4),
             main="mean MOW")
text(x=0.3, y=y, cex=2.5,
     labels= paste(round(rev(abs(a4_space_mMOW_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))


z <- barplot(rev(a4_space_dMOW), col = rev(cols),las=1,
             horiz=T, cex.main=2.5,cex.axis = 2.5, cex.names=2.5,
             axisnames = T,
             xlim=c(-0.4,0.4),
             main="delta MOW")
text(x=-0.3, y=y, cex=2.5,
     labels= paste(round(rev(abs(a4_space_dMOW_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

## combined x axis
text(-0.65,-1.0,expression(paste("Effect on spatial differences in ",
                                 alpha, " diversity (q = 4)",sep="")),
     xpd=NA,cex=2.5)

