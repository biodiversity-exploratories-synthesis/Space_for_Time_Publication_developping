# # # # # # # # # # # # # # #  # # # #
#                                    #
#        FIT LINEAR MODELS           # ----
#           ALL DATASETS             #
#                                    #
# # # # # # # # # # # # # # #  # # # #


# # # # # # # # # # # # # #
# CONTENT                  ----
# AIM : A function to fit linear models
# OUTPUT : 
#TODO

# # # # # # # # # # # # # #
# 0 - REQUIREMENTS              ----
library(vegan)
library(lme4)
library(nlme)
library(lmPerm)
library(MuMIn)#r.quared.GLMM
library(car)#Anova(type="II")


# # # # #
# 0.a. - CHOSE PARAMETERS ----
#
diversity_level <- c("alpha", "beta") # default is both alpha and beta
landuse_index <- c("LUI", "MOW", "GRA", "FER") # default is all indices
# a list of spatial and temporal dataset
datasets <- list("spatial" = NA, "temporal" = NA)
datasets[["spatial"]] <- matrix(data = 0, nrow = 2, ncol = 2)
datasets[["temporal"]] <- matrix(data = 0, nrow = 2, ncol = 2)

# a list of spatial and temporal dataset
datasets <- list("spatial" = NA, "temporal" = NA)
datasets[["spatial"]] <- matrix(data = 0, nrow = 2, ncol = 2)
datasets[["temporal"]] <- matrix(data = 0, nrow = 2, ncol = 2)


# setting parameters for LUI alpha
datasets[["spatial"]] <- matrix(data = 0, nrow = 2, ncol = 2)
datasets[["temporal"]] <- matrix(data = 0, nrow = 2, ncol = 2)






# # # # # # # # # # # # # #
# 1 - LINEAR MODELS              ----
#
load("./scaled_alpha_space_GRA.RData")
load("scaled_alpha_time_GRA.RData")



# # # # #
# 1.a. - SPATIAL DATASET ----
# 


# preparing results list
EV_space_alpha<-list()
lm_space_alpha<- list()
#I ran this loop manually since sometimes convergence problems make the loop stop
for(k in 1:3){ 
  lm_space1<- list()
  
  sub <- scaled_space_alpha[[k]]
  sub$YR<- as.factor(sub$YR) #TODO move this to data assembly part
  
  names(sub)[1]<- "a0"
  names(sub)[2]<- "a1"
  names(sub)[3]<- "a2"
  names(sub)[4]<- "a3"
  names(sub)[5]<- "a4"
  
  # SPACE
  # fit linear mixed effects models
  # explanatory variables are fixed
  # changing response variable Chao 0-4
  lm11<- lme(a0 ~ mGRA + dGRA + dph + mph + msur + dsur + dsoilPC1 + msoilPC1 + geo.dist,
             data = sub, random = ~1|YR)
  lm22<- lme(a1 ~ mGRA + dGRA + dph + mph + msur + dsur + dsoilPC1 + msoilPC1 + geo.dist,
             data = sub, random = ~1|YR)
  lm33<- lme(a2 ~ mGRA + dGRA + dph + mph + msur + dsur + dsoilPC1 + msoilPC1 + geo.dist,
             data = sub, random = ~1|YR)
  lm44<- lme(a3 ~ mGRA + dGRA + dph + mph + msur + dsur + dsoilPC1 + msoilPC1 + geo.dist,
             data = sub, random = ~1|YR)
  lm55<- lme(a4 ~ mGRA + dGRA + dph + mph + msur + dsur + dsoilPC1 + msoilPC1 + geo.dist,
             data = sub, random = ~1|YR)
  # save all models in a list
  lm_space1<- list(lm11,lm22, lm33, lm44, lm55)# if one model does not converge, set to list[0]
  
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
      
      a<- r.squaredGLMM(lm)[1]# marginal R2 of the fixed effects
      m<- (r.squaredGLMM(lm)[1]-r.squaredGLMM(lm1)[1])# variance expl by mLUI
      d<- (r.squaredGLMM(lm)[1]-r.squaredGLMM(lm2)[1])# variance expl by dLUI
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

#TODO adding NA to those vectors : why?
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
  
  # TIME
  # fit linear mixed effects models
  # explanatory variables are fixed
  # changing response variable Chao 0-4
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
  
  #TODO what does beta here? this is taken from GRA alpha!
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

