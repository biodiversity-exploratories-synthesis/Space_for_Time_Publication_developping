#load data - too large to load all at once
#pairwise differences
#plants
#load("./pwise_space_new.RData")
#load("./pwise_time_plants_new.RData")
load("./sitepred_pwise_space_step3.RData")#new pwise_space with updated LUI residuals
load("./sitepred_pwise_time_step3.RData")#new pwise_space with updated LUI residuals

pwise_time_plants<-pwise_time
pwise_space_plants<-pwise_space


#insects
#NB - both plant and insects files have the same names--> add org name to files
#load("./pwise_space_insects_new.RData")
#load("./pwise_time_insects_new.RData")
load("./sitepred_pwise_space_in_step3.RData")#new pwise_space with updated LUI residuals
load("./sitepred_pwise_time_in_step3.RData")#new pwise_space with updated LUI residuals

pwise_time_insects<-pwise_time
pwise_space_insects<-pwise_space

#remove unused files to prevent confusion
remove(pwise_space)
remove(pwise_time)

#Packages
library(lme4)
library(nlme)
library(lmPerm)
tlorder <- c("plants", "herbivore","secondary.consumer")
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

save(scaled_space, file="scaled_beta_space_GRA.RData")


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

save(scaled_time, file="scaled_beta_time_GRA.RData")

########################
###LINEAR MODELS
#######################
load("./scaled_beta_space_GRA.RData")
load("scaled_beta_time_GRA.RData")

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

save(EV_space, file="EV_space_GRA.RData")



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

save(EV_time, file="EV_time_GRA.RData")


####PLOTTING

load("./EV_time_GRA.RData")
load("./EV_space_GRA.RData")
#time
#EV
bsim_time_mGRA_EV<-c(EV_time[[1]]$mGRApp[1],EV_time[[2]]$mGRApp[1],EV_time[[3]]$mGRApp[1])
bsim_time_dGRA_EV<-c(EV_time[[1]]$dGRApp[1],EV_time[[2]]$dGRApp[1],EV_time[[3]]$dGRApp[1])

b0_time_mGRA_EV<-c(EV_time[[1]]$mGRApp[2],EV_time[[2]]$mGRApp[2],EV_time[[3]]$mGRApp[2])
b0_time_dGRA_EV<-c(EV_time[[1]]$dGRApp[2],EV_time[[2]]$dGRApp[2],EV_time[[3]]$dGRApp[2])

b1_time_mGRA_EV<-c(EV_time[[1]]$mGRApp[3],EV_time[[2]]$mGRApp[3],EV_time[[3]]$mGRApp[3])
b1_time_dGRA_EV<-c(EV_time[[1]]$dGRApp[3],EV_time[[2]]$dGRApp[3],EV_time[[3]]$dGRApp[3])

b2_time_mGRA_EV<-c(EV_time[[1]]$mGRApp[4],EV_time[[2]]$mGRApp[4],EV_time[[3]]$mGRApp[4])
b2_time_dGRA_EV<-c(EV_time[[1]]$dGRApp[4],EV_time[[2]]$dGRApp[4],EV_time[[3]]$dGRApp[4])

b3_time_mGRA_EV<-c(EV_time[[1]]$mGRApp[5],EV_time[[2]]$mGRApp[5],EV_time[[3]]$mGRApp[5])
b3_time_dGRA_EV<-c(EV_time[[1]]$dGRApp[5],EV_time[[2]]$dGRApp[5],EV_time[[3]]$dGRApp[5])

b4_time_mGRA_EV<-c(EV_time[[1]]$mGRApp[6],EV_time[[2]]$mGRApp[6],EV_time[[3]]$mGRApp[6])
b4_time_dGRA_EV<-c(EV_time[[1]]$dGRApp[6],EV_time[[2]]$dGRApp[6],EV_time[[3]]$dGRApp[6])

#coefficients
bsim_time_mGRA<-c(EV_time[[1]]$cofm[1],EV_time[[2]]$cofm[1],EV_time[[3]]$cofm[1])
bsim_time_dGRA<-c(EV_time[[1]]$cofd[1],EV_time[[2]]$cofd[1],EV_time[[3]]$cofd[1])

b0_time_mGRA<-c(EV_time[[1]]$cofm[2],EV_time[[2]]$cofm[2],EV_time[[3]]$cofm[2])
b0_time_dGRA<-c(EV_time[[1]]$cofd[2],EV_time[[2]]$cofd[2],EV_time[[3]]$cofd[2])

b1_time_mGRA<-c(EV_time[[1]]$cofm[3],EV_time[[2]]$cofm[3],EV_time[[3]]$cofm[3])
b1_time_dGRA<-c(EV_time[[1]]$cofd[3],EV_time[[2]]$cofd[3],EV_time[[3]]$cofd[3])

b2_time_mGRA<-c(EV_time[[1]]$cofm[4],EV_time[[2]]$cofm[4],EV_time[[3]]$cofm[4])
b2_time_dGRA<-c(EV_time[[1]]$cofd[4],EV_time[[2]]$cofd[4],EV_time[[3]]$cofd[4])

b3_time_mGRA<-c(EV_time[[1]]$cofm[5],EV_time[[2]]$cofm[5],EV_time[[3]]$cofm[5])
b3_time_dGRA<-c(EV_time[[1]]$cofd[5],EV_time[[2]]$cofd[5],EV_time[[3]]$cofd[5])

b4_time_mGRA<-c(EV_time[[1]]$cofm[6],EV_time[[2]]$cofm[6],EV_time[[3]]$cofm[6])
b4_time_dGRA<-c(EV_time[[1]]$cofd[6],EV_time[[2]]$cofd[6],EV_time[[3]]$cofd[6])


#space
#EV
bsim_space_mGRA_EV<-c(EV_space[[1]]$mGRApp[1],EV_space[[2]]$mGRApp[1],EV_space[[3]]$mGRApp[1])
bsim_space_dGRA_EV<-c(EV_space[[1]]$dGRApp[1],EV_space[[2]]$dGRApp[1],EV_space[[3]]$dGRApp[1])

b0_space_mGRA_EV<-c(EV_space[[1]]$mGRApp[2],EV_space[[2]]$mGRApp[2],EV_space[[3]]$mGRApp[2])
b0_space_dGRA_EV<-c(EV_space[[1]]$dGRApp[2],EV_space[[2]]$dGRApp[2],EV_space[[3]]$dGRApp[2])

b1_space_mGRA_EV<-c(EV_space[[1]]$mGRApp[3],EV_space[[2]]$mGRApp[3],EV_space[[3]]$mGRApp[3])
b1_space_dGRA_EV<-c(EV_space[[1]]$dGRApp[3],EV_space[[2]]$dGRApp[3],EV_space[[3]]$dGRApp[3])

b2_space_mGRA_EV<-c(EV_space[[1]]$mGRApp[4],EV_space[[2]]$mGRApp[4],EV_space[[3]]$mGRApp[4])
b2_space_dGRA_EV<-c(EV_space[[1]]$dGRApp[4],EV_space[[2]]$dGRApp[4],EV_space[[3]]$dGRApp[4])

b3_space_mGRA_EV<-c(EV_space[[1]]$mGRApp[5],EV_space[[2]]$mGRApp[5],EV_space[[3]]$mGRApp[5])
b3_space_dGRA_EV<-c(EV_space[[1]]$dGRApp[5],EV_space[[2]]$dGRApp[5],EV_space[[3]]$dGRApp[5])

b4_space_mGRA_EV<-c(EV_space[[1]]$mGRApp[6],EV_space[[2]]$mGRApp[6],EV_space[[3]]$mGRApp[6])
b4_space_dGRA_EV<-c(EV_space[[1]]$dGRApp[6],EV_space[[2]]$dGRApp[6],EV_space[[3]]$dGRApp[6])

#coefficients
bsim_space_mGRA<-c(EV_space[[1]]$cofm[1],EV_space[[2]]$cofm[1],EV_space[[3]]$cofm[1])
bsim_space_dGRA<-c(EV_space[[1]]$cofd[1],EV_space[[2]]$cofd[1],EV_space[[3]]$cofd[1])

b0_space_mGRA<-c(EV_space[[1]]$cofm[2],EV_space[[2]]$cofm[2],EV_space[[3]]$cofm[2])
b0_space_dGRA<-c(EV_space[[1]]$cofd[2],EV_space[[2]]$cofd[2],EV_space[[3]]$cofd[2])

b1_space_mGRA<-c(EV_space[[1]]$cofm[3],EV_space[[2]]$cofm[3],EV_space[[3]]$cofm[3])
b1_space_dGRA<-c(EV_space[[1]]$cofd[3],EV_space[[2]]$cofd[3],EV_space[[3]]$cofd[3])

b2_space_mGRA<-c(EV_space[[1]]$cofm[4],EV_space[[2]]$cofm[4],EV_space[[3]]$cofm[4])
b2_space_dGRA<-c(EV_space[[1]]$cofd[4],EV_space[[2]]$cofd[4],EV_space[[3]]$cofd[4])

b3_space_mGRA<-c(EV_space[[1]]$cofm[5],EV_space[[2]]$cofm[5],EV_space[[3]]$cofm[5])
b3_space_dGRA<-c(EV_space[[1]]$cofd[5],EV_space[[2]]$cofd[5],EV_space[[3]]$cofd[5])

b4_space_mGRA<-c(EV_space[[1]]$cofm[6],EV_space[[2]]$cofm[6],EV_space[[3]]$cofm[6])
b4_space_dGRA<-c(EV_space[[1]]$cofd[6],EV_space[[2]]$cofd[6],EV_space[[3]]$cofd[6])

tlnames <- c("plants","herbivores","secondary consumers")

cols <- c("#1b9e77", "#d95f02", "#7570b3")
letmat <- matrix(letters[1:9], nrow=3,byrow=T)

#Turnover
#time
par(mfrow=c(1,2))
par(mar=c(2,15,2,2))
par(oma=c(8,6,3,0))

y <- barplot(rev(bsim_time_mGRA), col = rev(cols),las=1,
             horiz=T, cex.main=2.5, cex.axis = 2.5, cex.names=2.5,
             axisnames = T,
             xlim=c(-0.4, 0.4), names=rev(tlnames), 
             main="mean GRA")
text(x=0.3, y=y, ,cex=2.5,
     labels= paste(round(rev(abs(bsim_time_mGRA_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))


z <- barplot(rev(bsim_time_dGRA), col = rev(cols),las=1,
             horiz=T, cex.main=2.5, cex.axis = 2.5, cex.names=2.5,
             axisnames = T,
             xlim=c(-0.4,0.4),
             main="delta GRA")

text(x=-0.3, y=z, ,cex=2.5,
     labels= paste(round(rev(abs(bsim_time_dGRA_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

## combined x axis
text(-0.75,-1,expression(paste("Effect on temporal ",beta, " diversity (turnover)",
                               sep="")),xpd=NA,cex=2.5)

#space
par(mfrow=c(1,2))
par(mar=c(2,15,2,2))
par(oma=c(8,6,3,0))

y <- barplot(rev(bsim_space_mGRA), col = rev(cols),las=1,
             horiz=T, cex.main=2.5, cex.axis = 2.5, cex.names=2.5,
             xlim=c(-0.4,0.4), names=rev(tlnames), axisnames = T,
             main="mean GRA")
text(x=0.3, y=y, cex=2.5,
     labels= paste(round(rev(abs(bsim_space_mGRA_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

z <- barplot(rev(bsim_space_dGRA), col = rev(cols),las=1,
             horiz=T, cex.main=2.5, cex.axis = 2.5, cex.names=2.5, axisnames = T,
             xlim=c(-0.4,0.4),
             main="delta GRA")
text(x=-0.3, y=y, cex=2.5,
     labels= paste(round(rev(abs(bsim_space_dGRA_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

## combined x axis
text(-0.75,-1.0,expression(paste("Effect on spatial ",beta, " diversity (turnover)",
                                 sep="")),xpd=NA,cex=2.5,)



#q 0
#time
par(mfrow=c(1,2))
par(mar=c(2,15,2,2))
par(oma=c(8,6,3,0))

y <- barplot(rev(b0_time_mGRA), col = rev(cols),las=1,
             horiz=T, cex.main=2.5, cex.axis = 2.5, cex.names=2.5, axisnames = T,
             xlim=c(-0.4,0.4), names=rev(tlnames), 
             main="mean GRA")
text(x=0.3, y=y, cex=2.5,
     labels= paste(round(rev(abs(b0_time_mGRA_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))


z <- barplot(rev(b0_time_dGRA), col = rev(cols),las=1,
             horiz=T, cex.main=2.5, cex.axis = 2.5, cex.names=2.5, axisnames = T,
             xlim=c(-0.4,0.4),
             main="delta GRA")

text(x=-0.3, y=z, cex=2.5,
     labels= paste(round(rev(abs(b0_time_dGRA_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

## combined x axis
text(-0.75,-1.0,expression(paste("Effect on temporal ",beta, " diversity (q = 0)",sep="")),xpd=NA,cex=2.5)

#space
par(mfrow=c(1,2))
par(mar=c(2,15,2,2))
par(oma=c(8,6,3,0))

y <- barplot(rev(b0_space_mGRA), col = rev(cols),las=1,
             horiz=T, cex.main=2.5, cex.axis = 2.5, cex.names=2.5, axisnames = T,
             xlim=c(-0.4,0.4), names=rev(tlnames), 
             main="mean GRA")

text(x=0.3, y=y, cex=2.5,
     labels= paste(round(rev(abs(b0_space_mGRA_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

z <- barplot(rev(b0_space_dGRA), col = rev(cols),las=1,
             horiz=T, cex.main=2.5, cex.axis = 2.5, cex.names=2.5, axisnames = T,
             xlim=c(-0.4,0.4),
             main="delta GRA")

text(x=-0.3, y=y, cex=2.5,
     labels= paste(round(rev(abs(b0_space_dGRA_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

## combined x axis
text(-0.75,-1.0,expression(paste("Effect on spatial ",beta, " diversity (q = 0)",sep="")),xpd=NA,cex=2.5)


#q 1
#time
par(mfrow=c(1,2))
par(mar=c(2,15,2,2))
par(oma=c(8,6,3,0))

y <- barplot(rev(b1_time_mGRA), col = rev(cols),las=1,
             horiz=T, cex.main=2.5, cex.axis = 2.5, cex.names=2.5, axisnames = T,
             xlim=c(-0.4,0.4), names=rev(tlnames), 
             main="mean GRA")

text(x=0.3, y=y, cex=2.5,
     labels= paste(round(rev(abs(b1_time_mGRA_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))


z <- barplot(rev(b1_time_dGRA), col = rev(cols),las=1,
             horiz=T, cex.main=2.5, cex.axis = 2.5, cex.names=2.5, axisnames = T,
             xlim=c(-0.4,0.4),
             main="delta GRA")

text(x=-0.3, y=z, cex=2.5,
     labels= paste(round(rev(abs(b1_time_dGRA_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

## combined x axis
text(-0.75,-1.0,expression(paste("Effect on temporal ",beta, " diversity (q = 1)",sep="")),xpd=NA,cex=2.5)

#space
par(mfrow=c(1,2))
par(mar=c(2,15,2,2))
par(oma=c(8,6,3,0))

y <- barplot(rev(b1_space_mGRA), col = rev(cols),las=1,
             horiz=T, cex.main=2.5, cex.axis = 2.5, cex.names=2.5, axisnames = T,
             xlim=c(-0.4,0.4), names=rev(tlnames), 
             main="mean GRA")
text(x=0.3, y=y, cex=2.5,
     labels= paste(round(rev(abs(b1_space_mGRA_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

z <- barplot(rev(b1_space_dGRA), col = rev(cols),las=1,
             horiz=T, cex.main=2.5, cex.axis = 2.5, cex.names=2.5, axisnames = T,
             xlim=c(-0.4,0.4),
             main="delta GRA")

text(x=-0.3, y=y, cex=2.5,
     labels= paste(round(rev(abs(b1_space_dGRA_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

## combined x axis
text(-0.75,-1.0,expression(paste("Effect on spatial ",beta, " diversity (q = 1)",sep="")),xpd=NA,cex=2.5)


#q 2
#time
par(mfrow=c(1,2))
par(mar=c(2,15,2,2))
par(oma=c(8,6,3,0))

y <- barplot(rev(b2_time_mGRA), col = rev(cols),las=1,
             horiz=T, cex.main=2.5, cex.axis = 2.5, cex.names=2.5, axisnames = T,
             xlim=c(-0.4,0.4), names=rev(tlnames), 
             main="mean GRA")

text(x=0.3, y=y, cex=2.5,
     labels= paste(round(rev(abs(b2_time_mGRA_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))


z <- barplot(rev(b2_time_dGRA), col = rev(cols),las=1,
             horiz=T, cex.main=2.5, cex.axis = 2.5, cex.names=2.5, axisnames = T,
             xlim=c(-0.4,0.4),
             main="delta GRA")
text(x=-0.3, y=z, cex=2.5,
     labels= paste(round(rev(abs(b2_time_dGRA_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))


## combined x axis
text(-0.75,-1.0,expression(paste("Effects on temporal ",beta, " diversity (q = 2)",sep="")),xpd=NA,cex=2.5)

#space
par(mfrow=c(1,2))
par(mar=c(2,15,2,2))
par(oma=c(8,6,3,0))

y <- barplot(rev(b2_space_mGRA), col = rev(cols),las=1,
             horiz=T, cex.main=2.5, cex.axis = 2.5, cex.names=2.5, axisnames = T,
             xlim=c(-0.4,0.4), names=rev(tlnames), 
             main="mean GRA")

text(x=0.3, y=y, cex=2.5,
     labels= paste(round(rev(abs(b2_space_mGRA_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))


z <- barplot(rev(b2_space_dGRA), col = rev(cols),las=1,
             horiz=T, cex.main=2.5, cex.axis = 2.5, cex.names=2.5, axisnames = T,
             xlim=c(-0.4,0.4),
             main="delta GRA")

text(x=-0.3, y=y, cex=2.5,
     labels= paste(round(rev(abs(b2_space_dGRA_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

## combined x axis
text(-0.75,-1.0,expression(paste("Effect on spatial ",beta, " diversity (q = 2)",
                                 sep="")),
     xpd=NA,cex=2.5)

#q 3
#time
par(mfrow=c(1,2))
par(mar=c(2,15,2,2))
par(oma=c(8,6,3,0))

y <- barplot(rev(b3_time_mGRA), col = rev(cols),las=1,
             horiz=T, cex.main=2.5, cex.axis = 2.5, cex.names=2.5, axisnames = T,
             xlim=c(-0.4,0.4), names=rev(tlnames), 
             main="mean GRA")

text(x=0.3, y=y, cex=2.5,
     labels= paste(round(rev(abs(b3_time_mGRA_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))


z <- barplot(rev(b3_time_dGRA), col = rev(cols),las=1,
             horiz=T, cex.main=2.5, cex.axis = 2.5, cex.names=2.5, axisnames = T,
             xlim=c(-0.4,0.4),
             main="delta GRA")
text(x=-0.3, y=z, cex=2.5,
     labels= paste(round(rev(abs(b3_time_dGRA_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))


## combined x axis
text(-0.75,-1.0,expression(paste("Effects on temporal ",beta, " diversity (q = 3)",sep="")),xpd=NA,cex=2.5)

#space
par(mfrow=c(1,2))
par(mar=c(2,15,2,2))
par(oma=c(8,6,3,0))

y <- barplot(rev(b3_space_mGRA), col = rev(cols),las=1,
             horiz=T, cex.main=2.5, cex.axis = 2.5, cex.names=2.5, axisnames = T,
             xlim=c(-0.4,0.4), names=rev(tlnames), 
             main="mean GRA")

text(x=0.3, y=y, cex=2.5,
     labels= paste(round(rev(abs(b3_space_mGRA_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))


z <- barplot(rev(b3_space_dGRA), col = rev(cols),las=1,
             horiz=T, cex.main=2.5, cex.axis = 2.5, cex.names=2.5, axisnames = T,
             xlim=c(-0.4,0.4),
             main="delta GRA")

text(x=-0.3, y=y, cex=2.5,
     labels= paste(round(rev(abs(b3_space_dGRA_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

## combined x axis
text(-0.75,-1.0,expression(paste("Effect on spatial ",beta, " diversity (q = 3)",
                                 sep="")),
     xpd=NA,cex=2.5)

#q 4
#time
par(mfrow=c(1,2))
par(mar=c(2,15,2,2))
par(oma=c(8,6,3,0))

y <- barplot(rev(b4_time_mGRA), col = rev(cols),las=1,
             horiz=T, cex.main=2.5, cex.axis = 2.5, cex.names=2.5, axisnames = T,
             xlim=c(-0.4,0.4), names=rev(tlnames), 
             main="mean GRA")

text(x=0.3, y=y, cex=2.5,
     labels= paste(round(rev(abs(b4_time_mGRA_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))


z <- barplot(rev(b4_time_dGRA), col = rev(cols),las=1,
             horiz=T, cex.main=2.5, cex.axis = 2.5, cex.names=2.5, axisnames = T,
             xlim=c(-0.4,0.4),
             main="delta GRA")
text(x=-0.3, y=z, cex=2.5,
     labels= paste(round(rev(abs(b4_time_dGRA_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))


## combined x axis
text(-0.75,-1.0,expression(paste("Effects on temporal ",beta, " diversity (q = 4)",sep="")),xpd=NA,cex=2.5)

#space
par(mfrow=c(1,2))
par(mar=c(2,15,2,2))
par(oma=c(8,6,3,0))

y <- barplot(rev(b4_space_mGRA), col = rev(cols),las=1,
             horiz=T, cex.main=2.5, cex.axis = 2.5, cex.names=2.5, axisnames = T,
             xlim=c(-0.4,0.4), names=rev(tlnames), 
             main="mean GRA")

text(x=0.3, y=y, cex=2.5,
     labels= paste(round(rev(abs(b4_space_mGRA_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))


z <- barplot(rev(b4_space_dGRA), col = rev(cols),las=1,
             horiz=T, cex.main=2.5, cex.axis = 2.5, cex.names=2.5, axisnames = T,
             xlim=c(-0.4,0.4),
             main="delta GRA")

text(x=-0.3, y=y, cex=2.5,
     labels= paste(round(rev(abs(b4_space_dGRA_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

## combined x axis
text(-0.75,-1.0,expression(paste("Effect on spatial ",beta, " diversity (q = 4)",
                                 sep="")),
     xpd=NA,cex=2.5)

