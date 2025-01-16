#load data - too large to load all at once
#plants
#load("./pwise_space_new.RData")
#load("./pwise_time_plants_new.RData")
load("./sitepred_pwise_space_step3.RData")#new pwise_space with updated LUI residuals
load("./sitepred_pwise_time_step3.RData")#new pwise_space with updated LUI residuals

pwise_time_plants<-as.data.frame(pwise_time) #new names as RData files
pwise_space_plants<-as.data.frame(pwise_space) #new names as RData files

rm(pwise_space)
rm(pwise_time)

#insects
#NB - both plant and insects files have the same names--> add org name to files
#load("./pwise_space_insects_new.RData")
#load("./pwise_time_insects_new.RData")
load("./sitepred_pwise_space_in_step3.RData")#new pwise_space with updated LUI residuals
load("./sitepred_pwise_time_in_step3.RData")#new pwise_space with updated LUI residuals

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


#coordinates -not needed anymore
pwise_time_insects$HWG1<- pwise_time_insects$HWG
pwise_time_insects$RWG1<- pwise_time_insects$RWG

colnames(pwise_time_insects)[8]<- "HWG"
colnames(pwise_time_insects)[8]

colnames(pwise_time_insects)[10]<- "RWG"
colnames(pwise_time_insects)[10]

RWG<- c()
  
  for(i in 1:8250) {
  rwg<-pwise_space$RWG1[pwise_space$EP1 %in% pwise_time_plants$EP[i]]
  RWG<- c(RWG, rwg[1])}
RWG[1:57]
pwise_time_plants$RWG<- RWG
colnames(pwise_time_plants)[3]<- "HWG"

#meanLUI as average of residuals normalised for EP and year differences
#temporal data
#plants

pwise_time_plants$mLUI<- (pwise_time_plants$LUI1res+pwise_time_plants$LUI2res)/2 

#insects
pwise_time_insects$mLUI<- (pwise_time_insects$LUI1res+pwise_time_insects$LUI2res)/2 


#spatial data - levelled for year differences
#plants
pwise_space_plants$mLUI<- (pwise_space_plants$LUI1res+pwise_space_plants$LUI2res)/2 

#insects
pwise_space_insects$mLUI<- (pwise_space_insects$LUI1res+pwise_space_insects$LUI2res)/2 


####SAVE DATA
save(pwise_time_insects, file="pwise_time_insects.RData")
save(pwise_time_plants, file="pwise_time_plants.RData")

save(pwise_space_insects, file="pwise_space_insects.RData")
save(pwise_space_plants, file="pwise_space_plants.RData")

#data prep - check whether all parameters are there
#soil parameters, temporal/geographic distance and isolation

colnames(pwise_space)#G500_1, G500_2, pH_1, pH_2, soilPCA_1, soilPCA_2, msoil_PC1, dsoil_PC1, msur, dsur, mph, dph
colnames(pwise_space_insects)#G500_1, G500_2, pH_1, pH_2, soilPCA_1, soilPCA_2, msoil_PC1, dsoil_PC1, msur, dsur, mph, dph


colnames(pwise_time_plants)#"pH" "soilPCA" "G500"  
colnames(pwise_time_insects)#"pH" "soilPCA"   "G500"    

#Packages
library(lme4)
library(nlme)
library(lmPerm)
library(MuMIn)#r.quared.GLMM
library(car)#Anova(type="II")

#tlorder <- c("plants", "herbivore", "pollinator", "decomposer", "secondary.consumer")
tlorder <- c("plants", "herbivore", "secondary.consumer")
#################PREPARATION####################
###############################################
#scaling all differences and means to max
####space
#plants
vars<- c("cqn0dis", "cqn1dis", "cqn2dis", "cqn3dis", "cqn4dis", "bsim", 
         "mLUI", "dLUI", "geo.dist", "YR", 
         "msur", "dsur","mph", "dph","msoilPC1", "dsoilPC1")

sub_plants_space<-pwise_space_plants[,colnames(pwise_space_plants) %in% vars]
sub_plants_space<- sub_plants_space[complete.cases(sub_plants_space),]

colnames(sub_plants_space)
preds<- apply(sub_plants_space[,c(2, 7:14)], 2, function(x)x/max(x, na.rm=T))
sub_plants_space<- cbind(sub_plants_space[c(3:6, 15,16, 1)], preds)
colnames(sub_plants_space)


#sec. consumers
vars<- c("pcqn0dis", "pcqn1dis", "pcqn2dis", "pcqn3dis", "pcqn4dis","pbsim", 
         "mLUI", "dLUI", "geo.dist", "YR",
         "msur", "dsur","mph", "dph","msoilPC1", "dsoilPC1")
sub_pred_space<- data.frame(pwise_space_insects[,colnames(pwise_space_insects) %in% vars])
sub_pred_space<- sub_pred_space[complete.cases(sub_pred_space),]

colnames(sub_pred_space)
preds<- apply(sub_pred_space[,c(2,6:13)], 2, function(x)x/max(x, na.rm=T))
sub_pred_space<- cbind(sub_pred_space[c(14, 3:5, 15,16, 1)], preds)
colnames(sub_pred_space)

#herbivores
vars<- c("hcqn0dis", "hcqn1dis", "hcqn2dis", "hcqn3dis","hcqn4dis","hbsim", 
         "mLUI", "dLUI", "geo.dist", "YR",
         "msur", "dsur","mph", "dph","msoilPC1", "dsoilPC1")
sub_her_space<- pwise_space_insects[,colnames(pwise_space_insects) %in% vars]
sub_her_space<- sub_her_space[complete.cases(sub_her_space),]

colnames(sub_her_space)
preds<- apply(sub_her_space[,c(2,7:14)], 2, function(x)x/max(x, na.rm=T))
sub_her_space<- cbind(sub_her_space[c(3:6, 15, 16,1)], preds)
colnames(sub_her_space)

scaled_space<- list()
scaled_space[[1]]<- sub_plants_space
scaled_space[[2]]<- sub_her_space
#scaled_space[[3]]<- sub_pol_space
#scaled_space[[4]]<- sub_dec_space
scaled_space[[3]]<- sub_pred_space
names(scaled_space)<- tlorder

save(scaled_space, file="scaled_beta_space_LUI.RData")


####time####
###########
#plants
vars<- c("cqn0dis", "cqn1dis", "cqn2dis", "cqn3dis","cqn4dis","bsim", 
         "mLUI", "dLUI", "dYR", "HWG", "RWG",
         "G500", "pH","soilPCA", "EP")
sub_plants_time<- pwise_time_plants[,colnames(pwise_time_plants) %in% vars]
sub_plants_time<- sub_plants_time[complete.cases(sub_plants_time),]

colnames(sub_plants_time)
preds<- apply(sub_plants_time[,c(2:4, 9:13)], 2, function(x)x/max(x, na.rm=T))
sub_plants_time<- cbind(sub_plants_time[c(5:8,14,15,1)], preds)
colnames(sub_plants_time)

#sec. consumers
vars<- c("pcqn0dis", "pcqn1dis", "pcqn2dis", "pcqn3dis", "pcqn4dis", "pbsim", 
         "mLUI", "dLUI", "dYR", "HWG", "RWG",
         "G500", "pH","soilPCA", "EP")
sub_pred_time<- pwise_time_insects[,colnames(pwise_time_insects) %in% vars]
sub_pred_time<- sub_pred_time[complete.cases(sub_pred_time),]

colnames(sub_pred_time)
preds<- apply(sub_pred_time[,c(2:5,7:9,13)], 2, function(x)x/max(x, na.rm=T))
sub_pred_time<- cbind(sub_pred_time[c(6,10:12,14,15,1)], preds)
colnames(sub_pred_time)


#herbivores
vars<- c("hcqn0dis", "hcqn1dis", "hcqn2dis", "hcqn3dis","hcqn4dis","hbsim", 
         "mLUI", "dLUI", "dYR", "HWG", "RWG",
         "G500", "pH","soilPCA", "EP")
sub_her_time<- pwise_time_insects[,colnames(pwise_time_insects) %in% vars]
sub_her_time<- sub_her_time[complete.cases(sub_her_time),]

colnames(sub_her_time)
preds<- apply(sub_her_time[,c(2:5,7:9,13)], 2, function(x)x/max(x, na.rm=T))
sub_her_time<- cbind(sub_her_time[c(6,10:12,14,15,1)], preds)
colnames(sub_her_time)

scaled_time<- list()
scaled_time[[1]]<- sub_plants_time
scaled_time[[2]]<- sub_her_time
#scaled_time[[3]]<- sub_pol_time
#scaled_time[[4]]<- sub_dec_time
scaled_time[[3]]<- sub_pred_time
names(scaled_time)<- tlorder

save(scaled_time, file="scaled_beta_time_LUI.RData")

########################
###LINEAR MODELS
#######################
load("./scaled_beta_space_LUI.RData")
load("scaled_beta_time_LUI.RData")

#space
EV_space<-list()
lm_space<- list()

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
 
 sub<- scaled_space[[k]]
 sub$YR<- as.factor(sub$YR)
 levels(sub$YR)
 levels(sub$YR) <- c("1","2","3", "4", "5", "6", "7", "8", "9", "10")
 
 names(sub)[1]<- "bsim"
 names(sub)[2]<- "b0"
 names(sub)[3]<- "b1"
 names(sub)[4]<- "b2"
 names(sub)[5]<- "b3"
 names(sub)[6]<- "b4"
 
 lm11<- lme(bsim~mLUI+dLUI+dph+mph+msur+dsur+dsoilPC1+msoilPC1+geo.dist,
         data=sub, random=~1|YR)
 lm22<- lme(b0~mLUI+dLUI+dph+mph+msur+dsur+dsoilPC1+msoilPC1+geo.dist,
           data=sub, random=~1|YR)
 lm33<- lme(b1~mLUI+dLUI+dph+mph+msur+dsur+dsoilPC1+msoilPC1+geo.dist,
           data=sub, random=~1|YR)
 lm44<- lme(b2~mLUI+dLUI+dph+mph+msur+dsur+dsoilPC1+msoilPC1+geo.dist,
           data=sub, random=~1|YR)
 lm55<- lme(b3~mLUI+dLUI+dph+mph+msur+dsur+dsoilPC1+msoilPC1+geo.dist,
            data=sub, random=~1|YR)
 lm66<- lme(b4~mLUI+dLUI+dph+mph+msur+dsur+dsoilPC1+msoilPC1+geo.dist,
            data=sub, random=~1|YR)
 
 lm_space1<- list(lm11,lm22, lm33, lm44, lm55,  lm66)
 
 for(l in 1:6){
  
    
   lm<- lm_space1[[l]]
   
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
EV_space[[k]]<- data.frame(all, mLUIpp, dLUIpp, sharedpp, cofm, cofd, int, pvm, pvd)
lm_space[[k]]<- lm_space1
#}

names(EV_space)<- tlorder
names(lm_space)<- tlorder

save(EV_space, file="EV_space_LUI_beta.RData")
save(lm_space, file="lm_space_LUI_beta.RData")


#time
lm_time<- list()
EV_time<-list()
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
  
  sub<- scaled_time[[k]]
  names(sub)[1]<- "bsim"
  names(sub)[2]<- "b0"
  names(sub)[3]<- "b1"
  names(sub)[4]<- "b2"
  names(sub)[5]<- "b3"
  names(sub)[6]<- "b4"
  
  lm11<- lme(bsim~mLUI+dLUI+pH+G500+soilPCA+dYR,
             data=sub, random=~1|EP)
  lm22<- lme(b0~mLUI+dLUI+pH+G500+soilPCA+dYR,
             data=sub, random=~1|EP)
  lm33<- lme(b1~mLUI+dLUI+pH+G500+soilPCA+dYR,
             data=sub, random=~1|EP)
  lm44<- lme(b2~mLUI+dLUI+pH+G500+soilPCA+dYR,
             data=sub, random=~1|EP)
  lm55<- lme(b3~mLUI+dLUI+pH+G500+soilPCA+dYR,
             data=sub, random=~1|EP)
  lm66<- lme(b4~mLUI+dLUI+pH+G500+soilPCA+dYR,
             data=sub, random=~1|EP)
  lm_time1<- list(lm11, lm22, lm33, lm44, lm55, lm66)
  
  #beta 0,1,2,3,4, sim
  for(l in 1:6){
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
  
  EV_time[[k]]<- data.frame(all, mLUIpp, dLUIpp, sharedpp, cofm, cofd, int, pvm, pvd)
  lm_time[[k]]<- lm_time1
}

names(EV_time)<- tlorder
names(lm_time)<- tlorder

save(EV_time, file="EV_time_LUI_beta.RData")
save(lm_time, file="lm_time_LUI_beta.RData")

####PLOTTING
#see/Use script Figures.R as the uptodate script for the figures

load("./EV_time_LUI_beta.RData")
load("./EV_space_LUI_beta.RData")

EV_time<- EV_time_LUI
EV_space<- EV_space_LUI
#time
#EV
bsim_time_mLUI_EV<-c(EV_time[[1]]$mLUIpp[1],EV_time[[2]]$mLUIpp[1],EV_time[[3]]$mLUIpp[1])
bsim_time_dLUI_EV<-c(EV_time[[1]]$dLUIpp[1],EV_time[[2]]$dLUIpp[1],EV_time[[3]]$dLUIpp[1])
                

b0_time_mLUI_EV<-c(EV_time[[1]]$mLUIpp[2],EV_time[[2]]$mLUIpp[2],EV_time[[3]]$mLUIpp[2])
b0_time_dLUI_EV<-c(EV_time[[1]]$dLUIpp[2],EV_time[[2]]$dLUIpp[2],EV_time[[3]]$dLUIpp[2])

b1_time_mLUI_EV<-c(EV_time[[1]]$mLUIpp[3],EV_time[[2]]$mLUIpp[3],EV_time[[3]]$mLUIpp[3])
b1_time_dLUI_EV<-c(EV_time[[1]]$dLUIpp[3],EV_time[[2]]$dLUIpp[3],EV_time[[3]]$dLUIpp[3])

b2_time_mLUI_EV<-c(EV_time[[1]]$mLUIpp[4],EV_time[[2]]$mLUIpp[4],EV_time[[3]]$mLUIpp[4])
b2_time_dLUI_EV<-c(EV_time[[1]]$dLUIpp[4],EV_time[[2]]$dLUIpp[4],EV_time[[3]]$dLUIpp[4])

b3_time_mLUI_EV<-c(EV_time[[1]]$mLUIpp[5],EV_time[[2]]$mLUIpp[5],EV_time[[3]]$mLUIpp[5])
b3_time_dLUI_EV<-c(EV_time[[1]]$dLUIpp[5],EV_time[[2]]$dLUIpp[5],EV_time[[3]]$dLUIpp[5])

b4_time_mLUI_EV<-c(EV_time[[1]]$mLUIpp[6],EV_time[[2]]$mLUIpp[6],EV_time[[3]]$mLUIpp[6])
b4_time_dLUI_EV<-c(EV_time[[1]]$dLUIpp[6],EV_time[[2]]$dLUIpp[6],EV_time[[3]]$dLUIpp[6])

#coefficients
bsim_time_mLUI<-c(EV_time[[1]]$cofm[1],EV_time[[2]]$cofm[1],EV_time[[3]]$cofm[1])
bsim_time_dLUI<-c(EV_time[[1]]$cofd[1],EV_time[[2]]$cofd[1],EV_time[[3]]$cofd[1])

b0_time_mLUI<-c(EV_time[[1]]$cofm[2],EV_time[[2]]$cofm[2],EV_time[[3]]$cofm[2])
b0_time_dLUI<-c(EV_time[[1]]$cofd[2],EV_time[[2]]$cofd[2],EV_time[[3]]$cofd[2])

b1_time_mLUI<-c(EV_time[[1]]$cofm[3],EV_time[[2]]$cofm[3],EV_time[[3]]$cofm[3])
b1_time_dLUI<-c(EV_time[[1]]$cofd[3],EV_time[[2]]$cofd[3],EV_time[[3]]$cofd[3])

b2_time_mLUI<-c(EV_time[[1]]$cofm[4],EV_time[[2]]$cofm[4],EV_time[[3]]$cofm[4])
b2_time_dLUI<-c(EV_time[[1]]$cofd[4],EV_time[[2]]$cofd[4],EV_time[[3]]$cofd[4])

b3_time_mLUI<-c(EV_time[[1]]$cofm[5],EV_time[[2]]$cofm[5],EV_time[[3]]$cofm[5])
b3_time_dLUI<-c(EV_time[[1]]$cofd[5],EV_time[[2]]$cofd[5],EV_time[[3]]$cofd[5])

b4_time_mLUI<-c(EV_time[[1]]$cofm[6],EV_time[[2]]$cofm[6],EV_time[[3]]$cofm[6])
b4_time_dLUI<-c(EV_time[[1]]$cofd[6],EV_time[[2]]$cofd[6],EV_time[[3]]$cofd[6])

#space
#EV
bsim_space_mLUI_EV<-c(EV_space[[1]]$mLUIpp[1],EV_space[[2]]$mLUIpp[1],EV_space[[3]]$mLUIpp[1])
bsim_space_dLUI_EV<-c(EV_space[[1]]$dLUIpp[1],EV_space[[2]]$dLUIpp[1],EV_space[[3]]$dLUIpp[1])

b0_space_mLUI_EV<-c(EV_space[[1]]$mLUIpp[2],EV_space[[2]]$mLUIpp[2],EV_space[[3]]$mLUIpp[2])
b0_space_dLUI_EV<-c(EV_space[[1]]$dLUIpp[2],EV_space[[2]]$dLUIpp[2],EV_space[[3]]$dLUIpp[2])

b1_space_mLUI_EV<-c(EV_space[[1]]$mLUIpp[3],EV_space[[2]]$mLUIpp[3],EV_space[[3]]$mLUIpp[3])
b1_space_dLUI_EV<-c(EV_space[[1]]$dLUIpp[3],EV_space[[2]]$dLUIpp[3],EV_space[[3]]$dLUIpp[3])

b2_space_mLUI_EV<-c(EV_space[[1]]$mLUIpp[4],EV_space[[2]]$mLUIpp[4],EV_space[[3]]$mLUIpp[4])
b2_space_dLUI_EV<-c(EV_space[[1]]$dLUIpp[4],EV_space[[2]]$dLUIpp[4],EV_space[[3]]$dLUIpp[4])

b3_space_mLUI_EV<-c(EV_space[[1]]$mLUIpp[5],EV_space[[2]]$mLUIpp[5],EV_space[[3]]$mLUIpp[5])
b3_space_dLUI_EV<-c(EV_space[[1]]$dLUIpp[5],EV_space[[2]]$dLUIpp[5],EV_space[[3]]$dLUIpp[5])

b4_space_mLUI_EV<-c(EV_space[[1]]$mLUIpp[6],EV_space[[2]]$mLUIpp[6],EV_space[[3]]$mLUIpp[6])
b4_space_dLUI_EV<-c(EV_space[[1]]$dLUIpp[6],EV_space[[2]]$dLUIpp[6],EV_space[[3]]$dLUIpp[6])

#coefficients
bsim_space_mLUI<-c(EV_space[[1]]$cofm[1],EV_space[[2]]$cofm[1],EV_space[[3]]$cofm[1])
bsim_space_dLUI<-c(EV_space[[1]]$cofd[1],EV_space[[2]]$cofd[1],EV_space[[3]]$cofd[1])

b0_space_mLUI<-c(EV_space[[1]]$cofm[2],EV_space[[2]]$cofm[2],EV_space[[3]]$cofm[2])
b0_space_dLUI<-c(EV_space[[1]]$cofd[2],EV_space[[2]]$cofd[2],EV_space[[3]]$cofd[2])

b1_space_mLUI<-c(EV_space[[1]]$cofm[3],EV_space[[2]]$cofm[3],EV_space[[3]]$cofm[3])
b1_space_dLUI<-c(EV_space[[1]]$cofd[3],EV_space[[2]]$cofd[3],EV_space[[3]]$cofd[3])

b2_space_mLUI<-c(EV_space[[1]]$cofm[4],EV_space[[2]]$cofm[4],EV_space[[3]]$cofm[4])
b2_space_dLUI<-c(EV_space[[1]]$cofd[4],EV_space[[2]]$cofd[4],EV_space[[3]]$cofd[4])

b3_space_mLUI<-c(EV_space[[1]]$cofm[5],EV_space[[2]]$cofm[5],EV_space[[3]]$cofm[5])
b3_space_dLUI<-c(EV_space[[1]]$cofd[5],EV_space[[2]]$cofd[5],EV_space[[3]]$cofd[5])

b4_space_mLUI<-c(EV_space[[1]]$cofm[6],EV_space[[2]]$cofm[6],EV_space[[3]]$cofm[6])
b4_space_dLUI<-c(EV_space[[1]]$cofd[6],EV_space[[2]]$cofd[6],EV_space[[3]]$cofd[6])

tlnames <- c("plants","herbivores","secondary consumers")

cols <- c("#1b9e77", "#d95f02", "#7570b3")
letmat <- matrix(letters[1:9], nrow=3,byrow=T)

#Turnover
#time
par(mfrow=c(1,2))
par(mar=c(2,2,1,2))
par(oma=c(6,8,3,0))

y <- barplot(rev(bsim_time_mLUI), col = rev(cols),las=1,
             horiz=T, cex.main=1.2, axisnames = T,
             xlim=c(-0.02, 0.03), names=rev(tlnames), 
             main="mean LUI")
  text(x=0.02, y=y, 
       labels= paste(round(rev(abs(bsim_time_mLUI_EV)), 0), "%", sep=""),
       adj=c(0.5,0.5))


z <- barplot(rev(bsim_time_dLUI), col = rev(cols),las=1,
             horiz=T, cex.main=1.2, axisnames = T,
             xlim=c(-0.02,0.03),
             main="delta LUI")

text(x=-0.01, y=z, 
     labels= paste(round(rev(abs(bsim_time_dLUI_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

## combined x axis
text(-0.025,-1,expression(paste("Effect on temporal ",beta, " diversity (turnover)",sep="")),xpd=NA,cex=1.2)

#space
par(mfrow=c(1,2))
par(mar=c(2,2,1,2))
par(oma=c(6,8,3,0))

y <- barplot(rev(bsim_space_mLUI), col = rev(cols),las=1,
             horiz=T, cex.main=1.2, axisnames = T,
             xlim=c(-0.3,0.15), names=rev(tlnames), 
             main="mean LUI")
text(x=0.05, y=y, 
     labels= paste(round(rev(abs(bsim_space_mLUI_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

z <- barplot(rev(bsim_space_dLUI), col = rev(cols),las=1,
             horiz=T, cex.main=1.2, axisnames = T,
             xlim=c(-0.1,0.2),
             main="delta LUI")
text(x=-0.05, y=y, 
     labels= paste(round(rev(abs(bsim_space_dLUI_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

## combined x axis
text(-0.1,-1.0,expression(paste("Effect on spatial ",beta, " diversity (turnover)",sep="")),xpd=NA,cex=1.2)



#Chao 0
#time
par(mfrow=c(1,2))
par(mar=c(2,2,1,2))
par(oma=c(6,8,3,0))

y <- barplot(rev(b0_time_mLUI), col = rev(cols),las=1,
             horiz=T, cex.main=1.2, axisnames = T,
             xlim=c(-0.05,0.05), names=rev(tlnames), 
             main="mean LUI")
text(x=0.03, y=y, 
     labels= paste(round(rev(abs(b0_time_mLUI_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))


z <- barplot(rev(b0_time_dLUI), col = rev(cols),las=1,
             horiz=T, cex.main=1.2, axisnames = T,
             xlim=c(-0.02,0.02),
             main="delta LUI")

text(x=-0.01, y=z, 
     labels= paste(round(rev(abs(b0_time_dLUI_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

## combined x axis
text(-0.03,-1.0,expression(paste("Effect on temporal ",beta, " diversity (Chao 0)",sep="")),xpd=NA,cex=1.2)

#space
par(mfrow=c(1,2))
par(mar=c(2,2,1,2))
par(oma=c(6,8,3,0))

y <- barplot(rev(b0_space_mLUI), col = rev(cols),las=1,
             horiz=T, cex.main=1.2, axisnames = T,
             xlim=c(-0.2,0.2), names=rev(tlnames), 
             main="mean LUI")

text(x=0.1, y=y, 
     labels= paste(round(rev(abs(b0_space_mLUI_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

z <- barplot(rev(b0_space_dLUI), col = rev(cols),las=1,
             horiz=T, cex.main=1.2, axisnames = T,
             xlim=c(-0.1,0.2),
             main="delta LUI")

text(x=-0.05, y=y, 
     labels= paste(round(rev(abs(b0_space_dLUI_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

## combined x axis
text(-0.15,-1,expression(paste("Effect on spatial ",beta, " diversity (Chao 0)",sep="")),xpd=NA,cex=1.2)


#Chao 1
#time
par(mfrow=c(1,2))
par(mar=c(2,2,1,2))
par(oma=c(6,8,3,0))

y <- barplot(rev(b1_time_mLUI), col = rev(cols),las=1,
             horiz=T, cex.main=1.2, axisnames = T,
             xlim=c(-0.05,0.05), names=rev(tlnames), 
             main="mean LUI")

text(x=0.03, y=y, 
     labels= paste(round(rev(abs(b1_time_mLUI_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))


z <- barplot(rev(b1_time_dLUI), col = rev(cols),las=1,
             horiz=T, cex.main=1.2, axisnames = T,
             xlim=c(-0.05,0.05),
             main="delta LUI")

text(x=-0.03, y=z, 
     labels= paste(round(rev(abs(b1_time_dLUI_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

## combined x axis
text(-0.06,-1,expression(paste("Effect on temporal ",beta, " diversity (Chao 1)",sep="")),xpd=NA,cex=1.2)

#space
par(mfrow=c(1,2))
par(mar=c(2,2,1,2))
par(oma=c(6,8,3,0))

y <- barplot(rev(b1_space_mLUI), col = rev(cols),las=1,
             horiz=T, cex.main=1.2, axisnames = T,
             xlim=c(-0.25,0.1), names=rev(tlnames), 
             main="mean LUI")
text(x=0.05, y=y, 
     labels= paste(round(rev(abs(b1_space_mLUI_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

z <- barplot(rev(b1_space_dLUI), col = rev(cols),las=1,
             horiz=T, cex.main=1.2, axisnames = T,
             xlim=c(-0.1,0.3),
             main="delta LUI")

text(x=-0.05, y=y, 
     labels= paste(round(rev(abs(b1_space_dLUI_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

## combined x axis
text(-0.1,-1,expression(paste("Effect on spatial ",beta, " diversity (Chao 1)",sep="")),xpd=NA,cex=1.2)


#Chao 2
#time
par(mfrow=c(1,2))
par(mar=c(2,2,1,2))
par(oma=c(6,8,3,0))

y <- barplot(rev(b2_time_mLUI), col = rev(cols),las=1,
             horiz=T, cex.main=1.2, axisnames = T,
             xlim=c(-0.02,0.02), names=rev(tlnames), 
             main="mean LUI")

text(x=0.01, y=y, 
     labels= paste(round(rev(abs(b2_time_mLUI_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))


z <- barplot(rev(b2_time_dLUI), col = rev(cols),las=1,
             horiz=T, cex.main=1.2, axisnames = T,
             xlim=c(-0.05,0.05),
             main="delta LUI")
text(x=-0.03, y=z, 
     labels= paste(round(rev(abs(b2_time_dLUI_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))


## combined x axis
text(-0.05,-1,expression(paste("Effects on temporal ",beta, " diversity (Chao 2)",sep="")),xpd=NA,cex=1.2)

#space
par(mfrow=c(1,2))
par(mar=c(2,2,1,2))
par(oma=c(6,8,3,0))

y <- barplot(rev(b2_space_mLUI), col = rev(cols),las=1,
             horiz=T, cex.main=1.2, axisnames = T,
             xlim=c(-0.2,0.15), names=rev(tlnames), 
             main="mean LUI")

text(x=0.08, y=y, 
     labels= paste(round(rev(abs(b2_space_mLUI_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))


z <- barplot(rev(b2_space_dLUI), col = rev(cols),las=1,
             horiz=T, cex.main=1.2, axisnames = T,
             xlim=c(-0.1,0.3),
             main="delta LUI")

text(x=-0.08, y=y, 
     labels= paste(round(rev(abs(b2_space_dLUI_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

## combined x axis
text(-0.2,-1,expression(paste("Effect on spatial ",beta, " diversity (Chao 2)",sep="")),xpd=NA,cex=1.2)


#Chao 3
#time
par(mfrow=c(1,2))
par(mar=c(2,2,1,2))
par(oma=c(6,8,3,0))

y <- barplot(rev(b3_time_mLUI), col = rev(cols),las=1,
             horiz=T, cex.main=1.2, axisnames = T,
             xlim=c(-0.02,0.02), names=rev(tlnames), 
             main="mean LUI")

text(x=0.01, y=y, 
     labels= paste(round(rev(abs(b3_time_mLUI_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))


z <- barplot(rev(b3_time_dLUI), col = rev(cols),las=1,
             horiz=T, cex.main=1.2, axisnames = T,
             xlim=c(-0.05,0.05),
             main="delta LUI")
text(x=-0.03, y=z, 
     labels= paste(round(rev(abs(b3_time_dLUI_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))


## combined x axis
text(-0.05,-1,expression(paste("Effects on temporal ",beta, " diversity (Chao 3)",sep="")),xpd=NA,cex=1.2)

#space
par(mfrow=c(1,2))
par(mar=c(2,2,1,2))
par(oma=c(6,8,3,0))

y <- barplot(rev(b3_space_mLUI), col = rev(cols),las=1,
             horiz=T, cex.main=1.2, axisnames = T,
             xlim=c(-0.2,0.15), names=rev(tlnames), 
             main="mean LUI")

text(x=0.08, y=y, 
     labels= paste(round(rev(abs(b3_space_mLUI_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))


z <- barplot(rev(b3_space_dLUI), col = rev(cols),las=1,
             horiz=T, cex.main=1.2, axisnames = T,
             xlim=c(-0.1,0.3),
             main="delta LUI")

text(x=-0.08, y=y, 
     labels= paste(round(rev(abs(b3_space_dLUI_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

## combined x axis
text(-0.2,-1,expression(paste("Effect on spatial ",beta, " diversity (Chao 3)",sep="")),xpd=NA,cex=1.2)

#Chao 4
#time
par(mfrow=c(1,2))
par(mar=c(2,2,1,2))
par(oma=c(6,8,3,0))

y <- barplot(rev(b4_time_mLUI), col = rev(cols),las=1,
             horiz=T, cex.main=1.2, axisnames = T,
             xlim=c(-0.02,0.02), names=rev(tlnames), 
             main="mean LUI")

text(x=0.01, y=y, 
     labels= paste(round(rev(abs(b4_time_mLUI_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))


z <- barplot(rev(b4_time_dLUI), col = rev(cols),las=1,
             horiz=T, cex.main=1.2, axisnames = T,
             xlim=c(-0.05,0.05),
             main="delta LUI")
text(x=-0.03, y=z, 
     labels= paste(round(rev(abs(b4_time_dLUI_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))


## combined x axis
text(-0.05,-1,expression(paste("Effects on temporal ",beta, " diversity (Chao 4)",sep="")),xpd=NA,cex=1.2)

#space
par(mfrow=c(1,2))
par(mar=c(2,2,1,2))
par(oma=c(6,8,3,0))

y <- barplot(rev(b4_space_mLUI), col = rev(cols),las=1,
             horiz=T, cex.main=1.2, axisnames = T,
             xlim=c(-0.2,0.15), names=rev(tlnames), 
             main="mean LUI")

text(x=0.08, y=y, 
     labels= paste(round(rev(abs(b4_space_mLUI_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))


z <- barplot(rev(b4_space_dLUI), col = rev(cols),las=1,
             horiz=T, cex.main=1.2, axisnames = T,
             xlim=c(-0.1,0.3),
             main="delta LUI")

text(x=-0.08, y=y, 
     labels= paste(round(rev(abs(b4_space_dLUI_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

## combined x axis
text(-0.2,-1,expression(paste("Effect on spatial ",beta, " diversity (Chao 4)",sep="")),xpd=NA,cex=1.2)
