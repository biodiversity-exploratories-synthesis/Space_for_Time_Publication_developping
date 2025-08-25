#Sensitivity analysis for linear models and GDMs 

# README ####

# content of the script
# 1 GDM with randomised LUI
# 2 Linear models with overall alpha diversity
# 3 Congruence analysis
# 4 Congruence plots

# The datasets used are synthesised of modelling
# results. For reasons of simplicity
# the aggregated datasets are uploaded with
# this code, as they were partially assembled 
# in R and partially in Excel.
# See method section in the publication
# Neuenkamp et al. 2025. Can space replace time? – Insights from a test 
# of multitrophic responses to land-use change. Nature Ecology & Evolution
# for detailed steps of dataset assembly.


# 1 GDM with randomised LUI ####
# idea: testing if asymptotic shape only shows
#       up if effect present
#       and not with random relationships

#       for this: testing 50 temporal GDMs
#                 with random LUI 
#                 (same range, random values)
#                 thought as expl.--> only plants

## a packages ####

library(gdm)
library(ggplot2)

## b data ####

#set working directory to folder "5. Sensitivity Analysis"
load("./InputData/pwise_space_plants.RData")
load("./InputData/pwise_space_herb.RData")
load("./InputData/pwise_space_pred.RData")

load("./InputData/pwise_time_plants.RData")
load("./InputData/pwise_time_herb.RData")
load("./InputData/pwise_time_pred.RData")

## c analysis ####

### i random LUI ####
#plants
LUIrandom<- matrix(NA, 8250, 100)

for(i in 1:100){
  
  m<- runif(8250, 0, 3.5)
  LUIrandom[,i]<- m
}

LUIrandom<- as.data.frame(LUIrandom)

save(LUIrandom, file="./InputData/LUIrandom.RData")

#herbivores
LUIrandomh<- matrix(NA, 5943, 100)

for(i in 1:100){
  
  m<- runif(5943, 0, 3.5)
  LUIrandomh[,i]<- m
}

LUIrandomh<- as.data.frame(LUIrandomh)

save(LUIrandomh, file="./InputData/LUIrandomh.RData")

#secondary consumers
LUIrandomp<- matrix(NA, 5943, 100)

for(i in 1:100){
  
  m<- runif(5943, 0, 3.5)
  LUIrandomp[,i]<- m
}

LUIrandomp<- as.data.frame(LUIrandomp)

save(LUIrandomp, file="./InputData/LUIrandomp.RData")

### ii GDMs ####

#### 1 plants ####
# create a data frame for the splines
splines<- list()

splines_plants_yall<- as.data.frame(matrix(nrow=200, ncol=4))
colnames(splines_plants_yall)<- c("LUIx",
                                  "COvLUIx",
                                  "LUIy",  
                                  "COvLUIy")

for(i in 1:50){
  
  attach(pwise_time_plants)
  
  #year and LU as predictors
  dat<- data.frame(cqn0dis, 
                   rep(1, nrow(pwise_time_plants)), 
                   HWG1, RWG1, HWG2, RWG2, 
                   LUIrandom[,i], 
                   YR1, 
                   LUIrandom[,i+50], 
                   YR2)
  
  colnames(dat)<- c("distance", "weights", 
                    "s1.xCoord", "s1.yCoord", "s2.xCoord", "s2.yCoord", 
                    "s1.LU","s1.Year", "s2.LU", "s2.Year")
  
  dat<- dat[complete.cases(dat),]
  
  detach(pwise_time)
  
  gdmTab<- formatsitepair(dat, 4, predData=dat)
  GDM<- gdm(gdmTab, geo=F)
  
  spline<- isplineExtract(GDM)
  splines_plants_yall[,1]<- spline$x[,1]
  splines_plants_yall[,2]<- spline$x[,2]
  
  splines_plants_yall[,3]<- spline$y[,1]
  splines_plants_yall[,4]<- spline$y[,2]
  
  splines[[i]]<- splines_plants_yall
}

remove(dat, gdmTab, GDM)

save(splines, file="./InputData/splines_plants.RData")

rm(splines)
rm(pwise_time)

#### 2 ins. herbivores ####
# create a data frame for the splines
splines<- list()

splines_herb_yall<- as.data.frame(matrix(nrow=200, ncol=4))
colnames(splines_herb_yall)<- c("LUIx",
                                "COvLUIx",
                                "LUIy",  
                                "COvLUIy")


for(i in 1:50){
  
  attach(pwise_time_herb)
  
  #year and LU as predictors
  dat<- data.frame(hcqn0dis, 
                   rep(1, nrow(pwise_time)), 
                   HWG1, RWG1, HWG2, RWG2, 
                   LUIrandomh[,i], 
                   YR1, 
                   LUIrandomh[,i+50], 
                   YR2)
  
  colnames(dat)<- c("distance", "weights", 
                    "s1.xCoord", "s1.yCoord", "s2.xCoord", "s2.yCoord", 
                    "s1.LU","s1.Year", "s2.LU", "s2.Year")
  
  dat<- dat[complete.cases(dat),]
  
  detach(pwise_time_herb)
  
  gdmTab<- formatsitepair(dat, 4, predData=dat)
  GDM<- gdm(gdmTab, geo=F)
  
  spline<- isplineExtract(GDM)
  splines_herb_yall[,1]<- spline$x[,1]
  splines_herb_yall[,2]<- spline$x[,2]
  
  splines_herb_yall[,3]<- spline$y[,1]
  splines_herb_yall[,4]<- spline$y[,2]
  
  splines[[i]]<- splines_herb_yall
}

remove(dat, gdmTab, GDM)

save(splines, file="./InputData/splines_herb.RData")

rm(splines)
rm(pwise_time_herb)
#### 3 ins. sec. consumers ####
# create a data frame for the splines

splines<- list()

splines_pred_yall<- as.data.frame(matrix(nrow=200, ncol=4))
colnames(splines_pred_yall)<- c("LUIx",
                                "COvLUIx",
                                "LUIy",  
                                "COvLUIy")

for(i in 1:50){
  
  attach(pwise_time_pred)
  
  #year and LU as predictors
  dat<- data.frame(pcqn0dis, 
                   rep(1, nrow(pwise_time)), 
                   HWG1, RWG1, HWG2, RWG2, 
                   LUIrandomh[,i], 
                   YR1, 
                   LUIrandomh[,i+50], 
                   YR2)
  
  colnames(dat)<- c("distance", "weights", 
                    "s1.xCoord", "s1.yCoord", "s2.xCoord", "s2.yCoord", 
                    "s1.LU","s1.Year", "s2.LU", "s2.Year")
  
  dat<- dat[complete.cases(dat),]
  
  detach(pwise_time_pred)
  
  gdmTab<- formatsitepair(dat, 4, predData=dat)
  GDM<- gdm(gdmTab, geo=F)
  
  spline<- isplineExtract(GDM)
  splines_pred_yall[,1]<- spline$x[,1]
  splines_pred_yall[,2]<- spline$x[,2]
  
  splines_pred_yall[,3]<- spline$y[,1]
  splines_pred_yall[,4]<- spline$y[,2]
  
  splines[[i]]<- splines_pred_yall
}

remove(dat, gdmTab, GDM)

save(splines, file="./InputData/splines_pred.RData")

rm(pwise_time)
rm(splines)

## d plotting ####

#### 1 plants ####
load("./InputData/splines_plants.RData")
plots<- list()

for(j in 1:50){
  sp<- splines[[j]]
  p1<- ggplot(sp, aes(x=LUIx, y=LUIy))+
    geom_line(linewidth=1.5)+
    xlim(0,3.5)+
    ylim(0, max(sp$LUIy))+
    xlab("mean LUI")+
    ylab(expression(paste("Effect of LUI differences on", Delta ," ", "alpha diversity",sep=" ")))+
    ggtitle("a - temporal land use intensity (LUI)\neffects on alpha diversity")+
    theme_classic()
  
  p1
  plots[[j]]<- p1
  
}

save(plots, file="./OutputData/plots_plants.RData")

for (k in 1:50) {
  file_name = paste("randomGDM_plants", k, ".tiff", sep="")
  tiff(file_name)
  print(plots[[k]])
  dev.off()
}

# Another option: create pdf where each page is a separate plot.
pdf("plots.pdf")
for (i in 1:50) {
  print(plot_list[[i]])
}
dev.off()

rm(plots)
rm(splines)

#### 2 ins. herbivores ####

load("./InputData/splines_herb.RData")
plots<- list()

for(j in 1:50){
  sp<- splines[[j]]
  p1<- ggplot(sp, aes(x=LUIx, y=LUIy))+
    geom_line(linewidth=1.5)+
    xlim(0,3.5)+
    ylim(0, max(sp$LUIy))+
    xlab("mean LUI")+
    ylab(expression(paste("Effect of LUI differences on", Delta ," ", "alpha diversity",sep=" ")))+
    ggtitle("a - temporal land use intensity (LUI)\neffects on alpha diversity")+
    theme_classic()
  
  p1
  plots[[j]]<- p1
  
}

save(plots, file="./OutputData/plots_herb.RData")

for (k in 1:50) {
  file_name = paste("randomGDM_herb", k, ".tiff", sep="")
  tiff(file_name)
  print(plots[[k]])
  dev.off()
}

# Another option: create pdf where each page is a separate plot.
pdf("plots.pdf")
for (i in 1:50) {
  print(plot_list[[i]])
}
dev.off()

rm(plots)
rm(splines)


#### 3 ins. sec. consumers ####

load("./InputData/splines_pred.RData")
plots<- list()

for(j in 1:50){
  sp<- splines[[j]]
  p1<- ggplot(sp, aes(x=LUIx, y=LUIy))+
    geom_line(linewidth=1.5)+
    xlim(0,3.5)+
    ylim(0, max(sp$LUIy))+
    xlab("mean LUI")+
    ylab(expression(paste("Effect of LUI differences on", Delta ," ", "alpha diversity",sep=" ")))+
    ggtitle("a - temporal land use intensity (LUI)\neffects on alpha diversity")+
    theme_classic()
  
  p1
  plots[[j]]<- p1
  
}

save(plots, file="./OutputData/plots_pred.RData")

for (k in 1:50) {
  file_name = paste("randomGDM_pred", k, ".tiff", sep="")
  tiff(file_name)
  print(plots[[k]])
  dev.off()
}

# Another option: create pdf where each page is a separate plot.
pdf("plots.pdf")
for (i in 1:50) {
  print(plot_list[[i]])
}
dev.off()

rm(plots)
rm(splines)


# 2 LMs with alpha diversity ####
# idea: comparing LUI effects on alpha diversity
#       across space and time
#       without using delta
#       time: residualing out EPs for all temp. variables
#       space: no prep needed
#       model structure otherwise like LMs

## a packages ####

library(lm)
library(ggplot2)
library(car)
library(patchwork)
library(vegan)
## b data ####

# option 1: use prepared data and continue with c analysis

# data prepared under 2 i. & ii. 00. data preparation
# all predictors are standardized (z scores)
# temporal data: EP residuals
load("./InputData/tplants.RData")
load("./InputData/therb.RData")
load("./InputData/tpred.RData")

# spatial data: year means
load("./InputData/splants.RData")
load("./InputData/sherb.RData")
load("./InputData/spred.RData")

# option 2: upload raw data and run first data preparation (i./ii.)
### 0. upload raw data ####
# this are the pair-wise datasets just reduced 
# so that each plot occurrs only once
# with its plot data
# This modification was done in Excel, data was uploaded here again.
# For this analysis insects data (herb, pred) were combined in one
# dataset. But can be done with the pairwise data and only in R also for each 
# group separately.

# Adjust folder to where data ara stored. 

div_plants<- read.csv("./InputData/div_plants_edited.csv",
                      header=T,
                      sep=";",
                      dec=".")
div_insects<- read.csv("./InputData/div_insects_edited.csv",
                       header=T,
                       sep=";",
                       dec=".")

plants<-div_plants[,c(1,4:28)]
plants<-plants[complete.cases(plants),]

herb<- div_insects[, c(1:12,18:31)]
herb<- herb[complete.cases(herb),]

pred<- div_insects[, c(1:7,13:31)]
pred<- pred[complete.cases(pred),]

### 00. data preparation ####

#### i. temporal analysis (EP residuals) ####
#scaled predictors
##### 1 plants ####
names(div_plants)

tempmodpl<- matrix(NA, 1613, 19)

for (i in 1:19){
  response<- plants[,(i+7)]
  explanatory<- plants$EP
  res<- residuals(lm (response ~ explanatory))
  tempmodpl[,i]<- res
  
}

#preds<- apply(tempmodpl, 2, function(x)x/max(x, na.rm=T))
preds<- decostand(tempmodpl, "range", 2, na.rm=T)
tplants<- data.frame(plants[,1:7], preds)
colnames(tplants)<- colnames(plants)

##### 2 ins. herbivores #### 
tempmodh<- matrix(NA, 1148, 19)
for (i in 1:19){
  response<- herb[,(i+7)]
  explanatory<- herb$EP
  res<- residuals(lm (response ~ explanatory))
  tempmodh[,i]<- res
  
}
#preds<- apply(tempmodh, 2, function(x)x/max(x, na.rm=T))
preds<- decostand(tempmodh, "range", 2, na.rm=T)
therb<- data.frame(herb[,1:7], preds)
colnames(therb)<- colnames(herb)

##### 3 ins. sec. consumers #### 

tempmodsc<- matrix(NA, 1111, 19)
for (i in 1:19){
  response<- pred[,(i+7)]
  explanatory<- pred$EP
  res<- residuals(lm (response ~ explanatory))
  tempmodsc[,i]<- res
  
}
#preds<- apply(tempmodsc, 2, function(x)x/max(x, na.rm=T))
preds<- decostand(tempmodsc, "range", 2, na.rm=T)
tpred<- data.frame(pred[,1:7], preds)
colnames(tpred)<- colnames(pred)


save(tplants, file="./InputData/tplants.RData")
save(therb, file="./InputData/therb.RData")
save(tpred, file="./InputData/tpred.RData")

#### ii. spatial analysis (year means) ####
?aggregate
#scaled predictors
##### 1 plants ####
length(levels(as.factor(plants$EP)))

spatmodpl<- matrix(NA, 149, 19)

for (i in 1:19){
  response<- plants[,(i+7)]
  factor<- plants$EP
  data<- data.frame(response, factor)
  mean<- aggregate(data, response ~ factor, mean)
  spatmodpl[,i]<- mean[,2]
  
}

HWG<- aggregate(plants, HWG1 ~ EP, mean)
RWG<- aggregate(plants, RWG1 ~ EP, mean)

preds<- apply(spatmodpl, 2, function(x)x/max(x, na.rm=T))
splants<- data.frame(levels(as.factor(factor)), HWG$HWG1, RWG$RWG1, preds)
colnames(splants)<- c("EP","HWG", "RWG", colnames(plants)[8:26])


##### 2 ins. herbivores #### 
length(levels(as.factor(herb$EP)))

spatmodh<- matrix(NA, 150, 19)

for (i in 1:19){
  response<- herb[,(i+7)]
  factor<- herb$EP
  data<- data.frame(response, factor)
  mean<- aggregate(data, response ~ factor, mean)
  spatmodh[,i]<- mean[,2]
  
}

HWG<- aggregate(herb, HWG1 ~ EP, mean)
RWG<- aggregate(herb, RWG1 ~ EP, mean)

preds<- apply(spatmodh, 2, function(x)x/max(x, na.rm=T))
sherb<- data.frame(levels(as.factor(factor)), HWG$HWG1, RWG$RWG1, preds)
colnames(sherb)<- c("EP","HWG", "RWG", colnames(herb)[8:26])


##### 3 ins. sec. consumers #### 

length(levels(as.factor(pred$EP)))

spatmodsc<- matrix(NA, 150, 19)

for (i in 1:19){
  response<- pred[,(i+7)]
  factor<- pred$EP
  data<- data.frame(response, factor)
  mean<- aggregate(data, response ~ factor, mean)
  spatmodsc[,i]<- mean[,2]
  
}

HWG<- aggregate(pred, HWG1 ~ EP, mean)
RWG<- aggregate(pred, RWG1 ~ EP, mean)

preds<- apply(spatmodsc, 2, function(x)x/max(x, na.rm=T))
spred<- data.frame(levels(as.factor(factor)), HWG$HWG1, RWG$RWG1, preds)
colnames(spred)<- c("EP","HWG", "RWG", colnames(pred)[8:26])

save(splants, file="./InputData/splants.RData")
save(sherb, file="./InputData/sherb.RData")
save(spred, file="./InputData/spred.RData")


## c analysis ####

# model equation LM on deltas

lm11<- lme(a0~mLUI+dLUI+dph+mph+msur+dsur+dsoilPC1+msoilPC1+geo.dist+
             mgperm+dgperm+mtemp+dtemp+mrain+drain,
           data=sub, random=~1|YR)

# we need only the plot values, i.e. no deltas and means

#lm_time<- lm(a0~LUI1+(LUI1)^2+pH+G500+soilPC1+HWG+RWG+
#            gperm+Temp_2m_Sum1+precip_raindays1,
#           data=data)

### i. temporal analysis ####

#### 1 plants ####
a0<- list()
a1<- list()
a2<- list()
a3<- list()
a4<- list()

model_time_plants<-list(a0, a1, a2, a3, a4)
anova_time_plants<- list(a0, a1, a2, a3, a4)

names(model_time_plants)<- c("a0", "a1", "a2", "a3", "a4")
names(anova_time_plants)<- c("a0", "a1", "a2", "a3", "a4")


names(tplants)

for (i in 1:5){
  a<- tplants[,i+7]
  modl<- lm(a~LUI1+I(LUI1^2)+pH+G500+soilPCA+HWG1+RWG1+
              gperm+Temp_2m_Sum1+precip_raindays1,
            data=tplants)
  modm<- lm(a~MOW1+I(MOW1^2)+pH+G500+soilPCA+HWG1+RWG1+
              gperm+Temp_2m_Sum1+precip_raindays1,
            data=tplants)
  modg<- lm(a~GRA1+I(GRA1^2)+pH+G500+soilPCA+HWG1+RWG1+
              gperm+Temp_2m_Sum1+precip_raindays1,
            data=tplants)
  modf<- lm(a~FER1+I(FER1^2)+pH+G500+soilPCA+HWG1+RWG1+
              gperm+Temp_2m_Sum1+precip_raindays1,
            data=tplants)
  
  testl<- Anova(modl, type="II")
  testm<- Anova(modm, type="II")
  testg<- Anova(modg, type="II")
  testf<- Anova(modf, type="II")
  
  model_time_plants[[i]][[1]]<- modl
  model_time_plants[[i]][[2]]<- modm
  model_time_plants[[i]][[3]]<- modg
  model_time_plants[[i]][[4]]<- modf
  
  anova_time_plants[[i]][[1]]<- testl
  anova_time_plants[[i]][[2]]<- testm
  anova_time_plants[[i]][[3]]<- testg
  anova_time_plants[[i]][[4]]<- testf
  
  names(model_time_plants[[i]])<-c("LUI", "MOW", "GRA", "FER") 
  names(anova_time_plants[[i]])<-c("LUI", "MOW", "GRA", "FER") 
}

save(model_time_plants, file="model_time_plants.RData")
save(anova_time_plants, file="anova_time_plants.RData")

rm(modl, modm, modg, modf,
   testl, testm, testg, testf)

#### 2 ins. herbivores ####

a0<- list()
a1<- list()
a2<- list()
a3<- list()
a4<- list()

model_time_herb<-list(a0, a1, a2, a3, a4)
anova_time_herb<- list(a0, a1, a2, a3, a4)

names(model_time_herb)<- c("a0", "a1", "a2", "a3", "a4")
names(anova_time_herb)<- c("a0", "a1", "a2", "a3", "a4")


names(therb)

for (i in 1:5){
  a<- therb[,i+7]
  modl<- lm(a~LUI1+I(LUI1^2)+pH+G500+soilPCA+HWG1+RWG1+
              gperm+Temp_2m_Sum1+precip_raindays1,
            data=therb)
  modm<- lm(a~MOW1+I(MOW1^2)+pH+G500+soilPCA+HWG1+RWG1+
              gperm+Temp_2m_Sum1+precip_raindays1,
            data=therb)
  modg<- lm(a~GRA1+I(GRA1^2)+pH+G500+soilPCA+HWG1+RWG1+
              gperm+Temp_2m_Sum1+precip_raindays1,
            data=therb)
  modf<- lm(a~FER1+I(FER1^2)+pH+G500+soilPCA+HWG1+RWG1+
              gperm+Temp_2m_Sum1+precip_raindays1,
            data=therb)
  
  testl<- Anova(modl, type="II")
  testm<- Anova(modm, type="II")
  testg<- Anova(modg, type="II")
  testf<- Anova(modf, type="II")
  
  model_time_herb[[i]][[1]]<- modl
  model_time_herb[[i]][[2]]<- modm
  model_time_herb[[i]][[3]]<- modg
  model_time_herb[[i]][[4]]<- modf
  
  anova_time_herb[[i]][[1]]<- testl
  anova_time_herb[[i]][[2]]<- testm
  anova_time_herb[[i]][[3]]<- testg
  anova_time_herb[[i]][[4]]<- testf
  
  names(model_time_herb[[i]])<-c("LUI", "MOW", "GRA", "FER") 
  names(anova_time_herb[[i]])<-c("LUI", "MOW", "GRA", "FER") 
}


save(model_time_herb, file="model_time_herb.RData")
save(anova_time_herb, file="anova_time_herb.RData")

rm(modl, modm, modg, modf,
   testl, testm, testg, testf)
#### 3 ins. dec. consumers ####
a0<- list()
a1<- list()
a2<- list()
a3<- list()
a4<- list()

model_time_pred<-list(a0, a1, a2, a3, a4)
anova_time_pred<- list(a0, a1, a2, a3, a4)

names(model_time_pred)<- c("a0", "a1", "a2", "a3", "a4")
names(anova_time_pred)<- c("a0", "a1", "a2", "a3", "a4")


names(tpred)

for (i in 1:5){
  a<- tpred[,i+7]
  modl<- lm(a~LUI1+I(LUI1^2)+pH+G500+soilPCA+HWG1+RWG1+
              gperm+Temp_2m_Sum1+precip_raindays1,
            data=tpred)
  modm<- lm(a~MOW1+I(MOW1^2)+pH+G500+soilPCA+HWG1+RWG1+
              gperm+Temp_2m_Sum1+precip_raindays1,
            data=tpred)
  modg<- lm(a~GRA1+I(GRA1^2)+pH+G500+soilPCA+HWG1+RWG1+
              gperm+Temp_2m_Sum1+precip_raindays1,
            data=tpred)
  modf<- lm(a~FER1+I(FER1^2)+pH+G500+soilPCA+HWG1+RWG1+
              gperm+Temp_2m_Sum1+precip_raindays1,
            data=tpred)
  
  testl<- Anova(modl, type="II")
  testm<- Anova(modm, type="II")
  testg<- Anova(modg, type="II")
  testf<- Anova(modf, type="II")
  
  model_time_pred[[i]][[1]]<- modl
  model_time_pred[[i]][[2]]<- modm
  model_time_pred[[i]][[3]]<- modg
  model_time_pred[[i]][[4]]<- modf
  
  anova_time_pred[[i]][[1]]<- testl
  anova_time_pred[[i]][[2]]<- testm
  anova_time_pred[[i]][[3]]<- testg
  anova_time_pred[[i]][[4]]<- testf
  
  names(model_time_pred[[i]])<-c("LUI", "MOW", "GRA", "FER") 
  names(anova_time_pred[[i]])<-c("LUI", "MOW", "GRA", "FER") 
}


save(model_time_pred, file="./OutputData/model_time_pred.RData")
save(anova_time_pred, file="./OutputData/anova_time_pred.RData")

rm(modl, modm, modg, modf,
   testl, testm, testg, testf)

### ii. spatial analysis ####

#### 1 plants ####
a0<- list()
a1<- list()
a2<- list()
a3<- list()
a4<- list()

model_space_plants<-list(a0, a1, a2, a3, a4)
anova_space_plants<- list(a0, a1, a2, a3, a4)

names(model_space_plants)<- c("a0", "a1", "a2", "a3", "a4")
names(anova_space_plants)<- c("a0", "a1", "a2", "a3", "a4")


names(splants)

for (i in 1:5){
  a<- splants[,i+3]
  modl<- lm(a~LUI1+I(LUI1^2)+pH+G500+soilPCA+HWG+RWG+
              gperm+Temp_2m_Sum1+precip_raindays1,
            data=splants)
  modm<- lm(a~MOW1+I(MOW1^2)+pH+G500+soilPCA+HWG+RWG+
              gperm+Temp_2m_Sum1+precip_raindays1,
            data=splants)
  modg<- lm(a~GRA1+I(GRA1^2)+pH+G500+soilPCA+HWG+RWG+
              gperm+Temp_2m_Sum1+precip_raindays1,
            data=splants)
  modf<- lm(a~FER1+I(FER1^2)+pH+G500+soilPCA+HWG+RWG+
              gperm+Temp_2m_Sum1+precip_raindays1,
            data=splants)
  
  testl<- Anova(modl, type="II")
  testm<- Anova(modm, type="II")
  testg<- Anova(modg, type="II")
  testf<- Anova(modf, type="II")
  
  model_space_plants[[i]][[1]]<- modl
  model_space_plants[[i]][[2]]<- modm
  model_space_plants[[i]][[3]]<- modg
  model_space_plants[[i]][[4]]<- modf
  
  anova_space_plants[[i]][[1]]<- testl
  anova_space_plants[[i]][[2]]<- testm
  anova_space_plants[[i]][[3]]<- testg
  anova_space_plants[[i]][[4]]<- testf
  
  names(model_space_plants[[i]])<-c("LUI", "MOW", "GRA", "FER") 
  names(anova_space_plants[[i]])<-c("LUI", "MOW", "GRA", "FER") 
}


save(model_space_plants, file="./OutputData/model_space_plants.RData")
save(anova_space_plants, file="./OutputData/anova_space_plants.RData")

rm(modl, modm, modg, modf,
   testl, testm, testg, testf)

#### 2 ins. herbivores ####

a0<- list()
a1<- list()
a2<- list()
a3<- list()
a4<- list()

model_space_herb<-list(a0, a1, a2, a3, a4)
anova_space_herb<- list(a0, a1, a2, a3, a4)

names(model_space_herb)<- c("a0", "a1", "a2", "a3", "a4")
names(anova_space_herb)<- c("a0", "a1", "a2", "a3", "a4")


names(sherb)

for (i in 1:5){
  a<- sherb[,i+3]
  modl<- lm(a~LUI1+I(LUI1^2)+pH+G500+soilPCA+HWG+RWG+
              gperm+Temp_2m_Sum1+precip_raindays1,
            data=sherb)
  modm<- lm(a~MOW1+I(MOW1^2)+pH+G500+soilPCA+HWG+RWG+
              gperm+Temp_2m_Sum1+precip_raindays1,
            data=sherb)
  modg<- lm(a~GRA1+I(GRA1^2)+pH+G500+soilPCA+HWG+RWG+
              gperm+Temp_2m_Sum1+precip_raindays1,
            data=sherb)
  modf<- lm(a~FER1+I(FER1^2)+pH+G500+soilPCA+HWG+RWG+
              gperm+Temp_2m_Sum1+precip_raindays1,
            data=sherb)
  
  testl<- Anova(modl, type="II")
  testm<- Anova(modm, type="II")
  testg<- Anova(modg, type="II")
  testf<- Anova(modf, type="II")
  
  model_space_herb[[i]][[1]]<- modl
  model_space_herb[[i]][[2]]<- modm
  model_space_herb[[i]][[3]]<- modg
  model_space_herb[[i]][[4]]<- modf
  
  anova_space_herb[[i]][[1]]<- testl
  anova_space_herb[[i]][[2]]<- testm
  anova_space_herb[[i]][[3]]<- testg
  anova_space_herb[[i]][[4]]<- testf
  
  names(model_space_herb[[i]])<-c("LUI", "MOW", "GRA", "FER") 
  names(anova_space_herb[[i]])<-c("LUI", "MOW", "GRA", "FER") 
}


save(model_space_herb, file="./OutputData/model_space_herb.RData")
save(anova_space_herb, file="./OutputData/anova_space_herb.RData")

rm(modl, modm, modg, modf,
   testl, testm, testg, testf)

#### 3 ins. dec. consumers ####
a0<- list()
a1<- list()
a2<- list()
a3<- list()
a4<- list()

model_space_pred<-list(a0, a1, a2, a3, a4)
anova_space_pred<- list(a0, a1, a2, a3, a4)

names(model_space_pred)<- c("a0", "a1", "a2", "a3", "a4")
names(anova_space_pred)<- c("a0", "a1", "a2", "a3", "a4")


names(spred)

for (i in 1:5){
  a<- spred[,i+3]
  modl<- lm(a~LUI1+I(LUI1^2)+pH+G500+soilPCA+HWG+RWG+
              gperm+Temp_2m_Sum1+precip_raindays1,
            data=spred)
  modm<- lm(a~MOW1+I(MOW1^2)+pH+G500+soilPCA+HWG+RWG+
              gperm+Temp_2m_Sum1+precip_raindays1,
            data=spred)
  modg<- lm(a~GRA1+I(GRA1^2)+pH+G500+soilPCA+HWG+RWG+
              gperm+Temp_2m_Sum1+precip_raindays1,
            data=spred)
  modf<- lm(a~FER1+I(FER1^2)+pH+G500+soilPCA+HWG+RWG+
              gperm+Temp_2m_Sum1+precip_raindays1,
            data=spred)
  
  testl<- Anova(modl, type="II")
  testm<- Anova(modm, type="II")
  testg<- Anova(modg, type="II")
  testf<- Anova(modf, type="II")
  
  model_space_pred[[i]][[1]]<- modl
  model_space_pred[[i]][[2]]<- modm
  model_space_pred[[i]][[3]]<- modg
  model_space_pred[[i]][[4]]<- modf
  
  anova_space_pred[[i]][[1]]<- testl
  anova_space_pred[[i]][[2]]<- testm
  anova_space_pred[[i]][[3]]<- testg
  anova_space_pred[[i]][[4]]<- testf
  
  names(model_space_pred[[i]])<-c("LUI", "MOW", "GRA", "FER") 
  names(anova_space_pred[[i]])<-c("LUI", "MOW", "GRA", "FER") 
}


save(model_space_pred, file="./OutputData/model_space_pred.RData")
save(anova_space_pred, file="./OutputData/anova_space_pred.RData")

rm(modl, modm, modg, modf,
   testl, testm, testg, testf)

### iii. joined analysis ####

#### 1 plants ####
# mean LUI
response<- plants$LUI1
factor<- plants$EP
data<- data.frame(response, factor)
mean<- aggregate(data, response ~ factor, mean)

mLUI<- c()
for (i in 1:length(plants$EP)){
  ep<- plants$EP[i]
  m<- mean$response[mean$factor %in% ep]
  mLUI<- c(mLUI,m)
}

plants$mLUI<- mLUI
rm(mLUI)

#### 2 ins. herbivores ####
# mean LUI
response<- herb$LUI1
factor<- herb$EP
data<- data.frame(response, factor)
mean<- aggregate(data, response ~ factor, mean)

mLUI<- c()
for (i in 1:length(herb$EP)){
  ep<- herb$EP[i]
  m<- mean$response[mean$factor %in% ep]
  mLUI<- c(mLUI,m)
}

herb$mLUI<- mLUI
rm(mLUI)

#### 3 ins. sec. consumers ####
# mean LUI
response<- pred$LUI1
factor<- pred$EP
data<- data.frame(response, factor)
mean<- aggregate(data, response ~ factor, mean)

mLUI<- c()
for (i in 1:length(pred$EP)){
  ep<- pred$EP[i]
  m<- mean$response[mean$factor %in% ep]
  mLUI<- c(mLUI,m)
}

pred$mLUI<- mLUI
rm(mLUI)

par(mfrow=c(1,3))

plot(plants$LUI1, plants$mLUI)
plot(herb$LUI1, herb$mLUI)
plot(pred$LUI1, pred$mLUI)

cor(plants$LUI1, plants$mLUI)
cor(herb$LUI1, herb$mLUI)
cor(pred$LUI1, pred$mLUI)
## d model validation ####

# no comments means: all ok

### i. temporal analysis ####

#### 1 plants ####

load("./OutputData/model_time_plants.RData")

par(mfrow=c(2,2))
plot(model_time_plants$a0$LUI)
plot(model_time_plants$a1$LUI)
plot(model_time_plants$a2$LUI)
plot(model_time_plants$a3$LUI)
plot(model_time_plants$a4$LUI)

par(mfrow=c(2,2))
plot(model_time_plants$a0$MOW)
plot(model_time_plants$a1$MOW)
plot(model_time_plants$a2$MOW)
plot(model_time_plants$a3$MOW)
plot(model_time_plants$a4$MOW)

par(mfrow=c(2,2))
plot(model_time_plants$a0$GRA)
plot(model_time_plants$a1$GRA)
plot(model_time_plants$a2$GRA)
plot(model_time_plants$a3$GRA)
plot(model_time_plants$a4$GRA)

par(mfrow=c(2,2))
plot(model_time_plants$a0$FER)
plot(model_time_plants$a1$FER)
plot(model_time_plants$a2$FER)
plot(model_time_plants$a3$FER)
plot(model_time_plants$a4$FER)

#### 2 ins. herbivores ####
#AEG01 as outlier --> remove
load("./OutputData/model_time_herb.RData")

# two plots with high diversity, especially for a0,a1
# for now ok, see when plotting if problematic

par(mfrow=c(2,2))
plot(model_time_herb$a0$LUI)
plot(model_time_herb$a1$LUI)
plot(model_time_herb$a2$LUI)
plot(model_time_herb$a3$LUI)
plot(model_time_herb$a4$LUI)

par(mfrow=c(2,2))
plot(model_time_herb$a0$MOW)
plot(model_time_herb$a1$MOW)
plot(model_time_herb$a2$MOW)
plot(model_time_herb$a3$MOW)
plot(model_time_herb$a4$MOW)

par(mfrow=c(2,2))
plot(model_time_herb$a0$GRA)
plot(model_time_herb$a1$GRA)
plot(model_time_herb$a2$GRA)
plot(model_time_herb$a3$GRA)
plot(model_time_herb$a4$GRA)

par(mfrow=c(2,2))
plot(model_time_herb$a0$FER)
plot(model_time_herb$a1$FER)
plot(model_time_herb$a2$FER)
plot(model_time_herb$a3$FER)
plot(model_time_herb$a4$FER)

#### 2 ins. sec. consumers ####
#AEG01 as outlier
load("./OutputData/model_time_pred.RData")

par(mfrow=c(2,2))
plot(model_time_pred$a0$LUI)
plot(model_time_pred$a1$LUI)
plot(model_time_pred$a2$LUI)
plot(model_time_pred$a3$LUI)
plot(model_time_pred$a4$LUI)

par(mfrow=c(2,2))
plot(model_time_pred$a0$MOW)
plot(model_time_pred$a1$MOW)
plot(model_time_pred$a2$MOW)
plot(model_time_pred$a3$MOW)
plot(model_time_pred$a4$MOW)

par(mfrow=c(2,2))
plot(model_time_pred$a0$GRA)
plot(model_time_pred$a1$GRA)
plot(model_time_pred$a2$GRA)
plot(model_time_pred$a3$GRA)
plot(model_time_pred$a4$GRA)

par(mfrow=c(2,2))
plot(model_time_pred$a0$FER)
plot(model_time_pred$a1$FER)
plot(model_time_pred$a2$FER)
plot(model_time_pred$a3$FER)
plot(model_time_pred$a4$FER)

### ii. spatial analysis ####

#### 1 plants ####

load("./OutputData/model_space_plants.RData")

# plots 59,70,91 and 93 inflate the variance
# higher diversity than the others
# check results, how influencial they are


par(mfrow=c(2,2))
plot(model_space_plants$a0$LUI)
plot(model_space_plants$a1$LUI)
plot(model_space_plants$a2$LUI)
plot(model_space_plants$a3$LUI)
plot(model_space_plants$a4$LUI)

par(mfrow=c(2,2))
plot(model_space_plants$a0$MOW)
plot(model_space_plants$a1$MOW)
plot(model_space_plants$a2$MOW)
plot(model_space_plants$a3$MOW)
plot(model_space_plants$a4$MOW)

par(mfrow=c(2,2))
plot(model_space_plants$a0$GRA)
plot(model_space_plants$a1$GRA)
plot(model_space_plants$a2$GRA)
plot(model_space_plants$a3$GRA)
plot(model_space_plants$a4$GRA)

par(mfrow=c(2,2))
plot(model_space_plants$a0$FER)
plot(model_space_plants$a1$FER)
plot(model_space_plants$a2$FER)
plot(model_space_plants$a3$FER)
plot(model_space_plants$a4$FER)

#### 2 ins. herbivores ####

load("./model_space_herb.RData")

par(mfrow=c(2,2))
plot(model_space_herb$a0$LUI)
plot(model_space_herb$a1$LUI)
plot(model_space_herb$a2$LUI)
plot(model_space_herb$a3$LUI)
plot(model_space_herb$a4$LUI)

par(mfrow=c(2,2))
plot(model_space_herb$a0$MOW)
plot(model_space_herb$a1$MOW)
plot(model_space_herb$a2$MOW)
plot(model_space_herb$a3$MOW)
plot(model_space_herb$a4$MOW)

par(mfrow=c(2,2))
plot(model_space_herb$a0$GRA)
plot(model_space_herb$a1$GRA)
plot(model_space_herb$a2$GRA)
plot(model_space_herb$a3$GRA)
plot(model_space_herb$a4$GRA)

par(mfrow=c(2,2))
plot(model_space_herb$a0$FER)
plot(model_space_herb$a1$FER)
plot(model_space_herb$a2$FER)
plot(model_space_herb$a3$FER)
plot(model_space_herb$a4$FER)

#### 2 ins. sec. consumers ####

load("./OutputData/model_space_pred.RData")

par(mfrow=c(2,2))
plot(model_space_pred$a0$LUI)
plot(model_space_pred$a1$LUI)
plot(model_space_pred$a2$LUI)
plot(model_space_pred$a3$LUI)
plot(model_space_pred$a4$LUI)

par(mfrow=c(2,2))
plot(model_space_pred$a0$MOW)
plot(model_space_pred$a1$MOW)
plot(model_space_pred$a2$MOW)
plot(model_space_pred$a3$MOW)
plot(model_space_pred$a4$MOW)

par(mfrow=c(2,2))
plot(model_space_pred$a0$GRA)
plot(model_space_pred$a1$GRA)
plot(model_space_pred$a2$GRA)
plot(model_space_pred$a3$GRA)
plot(model_space_pred$a4$GRA)

par(mfrow=c(2,2))
plot(model_space_pred$a0$FER)
plot(model_space_pred$a1$FER)
plot(model_space_pred$a2$FER)
plot(model_space_pred$a3$FER)
plot(model_space_pred$a4$FER)

## d model results ####

### 0. result table preparation ####

anova_time_plants[[1]][[1]]$Df

model_time_plants[[1]][[1]]$coefficients[1] #intercept
model_time_plants[[1]][[1]]$coefficients[2] #LU
model_time_plants[[1]][[1]]$coefficients[3] #LU2

LMresults<- matrix(NA, nrow=(2*3*4*5*3), ncol=9)
colnames(LMresults)<- c("model type",
                        "land use component",
                        "organism group",
                        "diversity abundance weighting",
                        "predictor",
                        "estimate",
                        "df",
                        "F-value",
                        "p-value")
LMresults[,1]<- (rep(c("temporal", "spatial"), each=180))
LMresults[,2]<- (rep(c("LUI", "MOW", "GRA", "FER"), each=9, times=10))
LMresults[,3]<- (rep(c("plants", "herbivores", "sec.consumers"), each=3, times=40))
LMresults[,4]<- (rep(c("q0", "q1", "q2", "q3", "q4"), each=36,times=2))
LMresults[,5]<- (rep(c("intL", "LUI", "LUI2",
                       "intL", "LUI", "LUI2",
                       "intL", "LUI", "LUI2",
                       "intM","MOW", "MOW2",
                       "intM","MOW", "MOW2",
                       "intM","MOW", "MOW2",
                       "intG","GRA", "GRA2",
                       "intG","GRA", "GRA2",
                       "intG","GRA", "GRA2",
                       "intF","FER", "FER2",
                       "intF","FER", "FER2",
                       "intF","FER", "FER2"), times=10))
LMresults[,7]<- rep(1, times=360)

### i. temporal analysis ####

load("./OutputData/anova_time_plants.RData")
load("./OutputData/model_time_plants.RData")

load("./OutputData/anova_time_herb.RData")
load("./OutputData/model_time_herb.RData")

load("./OutputData/anova_time_pred.RData")
load("./OutputData/model_time_pred.RData")


estt<- c()
Fvalt<- c()
pvalt<- c()


for(i in 1:5){
  for(k in 1:4){
    
    #estimate
    estt<- c(estt,model_time_plants[[i]][[k]]$coefficients[1])
    estt<-c(estt, model_time_plants[[i]][[k]]$coefficients[2])
    estt<-c(estt, model_time_plants[[i]][[k]]$coefficients[3])
    
    estt<- c(estt,model_time_herb[[i]][[k]]$coefficients[1])
    estt<-c(estt, model_time_herb[[i]][[k]]$coefficients[2])
    estt<-c(estt, model_time_herb[[i]][[k]]$coefficients[3])
    
    estt<- c(estt,model_time_pred[[i]][[k]]$coefficients[1])
    estt<-c(estt, model_time_pred[[i]][[k]]$coefficients[2])
    estt<-c(estt, model_time_pred[[i]][[k]]$coefficients[3])
    
    #Fvalue
    Fvalt<- c(Fvalt,NA)
    Fvalt<-c(Fvalt, anova_time_plants[[i]][[k]]$`F value`[1])
    Fvalt<-c(Fvalt, anova_time_plants[[i]][[k]]$`F value`[2])
    
    Fvalt<- c(Fvalt,NA)
    Fvalt<-c(Fvalt, anova_time_herb[[i]][[k]]$`F value`[1])
    Fvalt<-c(Fvalt, anova_time_herb[[i]][[k]]$`F value`[2])
    
    Fvalt<- c(Fvalt,NA)
    Fvalt<-c(Fvalt, anova_time_pred[[i]][[k]]$`F value`[1])
    Fvalt<-c(Fvalt, anova_time_pred[[i]][[k]]$`F value`[2])
    
    #pvalue
    pvalt<- c(pvalt,NA)
    pvalt<-c(pvalt, anova_time_plants[[i]][[k]]$`Pr(>F)`[1])
    pvalt<-c(pvalt, anova_time_plants[[i]][[k]]$`Pr(>F)`[2])
    
    pvalt<- c(pvalt,NA)
    pvalt<-c(pvalt, anova_time_herb[[i]][[k]]$`Pr(>F)`[1])
    pvalt<-c(pvalt, anova_time_herb[[i]][[k]]$`Pr(>F)`[2])
    
    pvalt<- c(pvalt,NA)
    pvalt<-c(pvalt, anova_time_pred[[i]][[k]]$`Pr(>F)`[1])
    pvalt<-c(pvalt, anova_time_pred[[i]][[k]]$`Pr(>F)`[2])
    
    length(estt)
    length(Fvalt)
    length(pvalt)
  }
}


### ii. spatial analysis ####

load("./OutputData/anova_space_plants.RData")
load("./OutputData/model_space_plants.RData")

load("./OutputData/anova_space_herb.RData")
load("./OutputData/model_space_herb.RData")

load("./OutputData/anova_space_pred.RData")
load("./OutputData/model_space_pred.RData")

ests<- c()
Fvals<- c()
pvals<- c()


for(i in 1:5){
  for(k in 1:4){
    
    #estimate
    ests<- c(ests,model_space_plants[[i]][[k]]$coefficients[1])
    ests<-c(ests, model_space_plants[[i]][[k]]$coefficients[2])
    ests<-c(ests, model_space_plants[[i]][[k]]$coefficients[3])
    
    ests<- c(ests,model_space_herb[[i]][[k]]$coefficients[1])
    ests<-c(ests, model_space_herb[[i]][[k]]$coefficients[2])
    ests<-c(ests, model_space_herb[[i]][[k]]$coefficients[3])
    
    ests<- c(ests,model_space_pred[[i]][[k]]$coefficients[1])
    ests<-c(ests, model_space_pred[[i]][[k]]$coefficients[2])
    ests<-c(ests, model_space_pred[[i]][[k]]$coefficients[3])
    
    #Fvalue
    Fvals<- c(Fvals,NA)
    Fvals<-c(Fvals, anova_space_plants[[i]][[k]]$`F value`[1])
    Fvals<-c(Fvals, anova_space_plants[[i]][[k]]$`F value`[2])
    
    Fvals<- c(Fvals,NA)
    Fvals<-c(Fvals, anova_space_herb[[i]][[k]]$`F value`[1])
    Fvals<-c(Fvals, anova_space_herb[[i]][[k]]$`F value`[2])
    
    Fvals<- c(Fvals,NA)
    Fvals<-c(Fvals, anova_space_pred[[i]][[k]]$`F value`[1])
    Fvals<-c(Fvals, anova_space_pred[[i]][[k]]$`F value`[2])
    
    #pvalue
    pvals<- c(pvals,NA)
    pvals<-c(pvals, anova_space_plants[[i]][[k]]$`Pr(>F)`[1])
    pvals<-c(pvals, anova_space_plants[[i]][[k]]$`Pr(>F)`[2])
    
    pvals<- c(pvals,NA)
    pvals<-c(pvals, anova_space_herb[[i]][[k]]$`Pr(>F)`[1])
    pvals<-c(pvals, anova_space_herb[[i]][[k]]$`Pr(>F)`[2])
    
    pvals<- c(pvals,NA)
    pvals<-c(pvals, anova_space_pred[[i]][[k]]$`Pr(>F)`[1])
    pvals<-c(pvals, anova_space_pred[[i]][[k]]$`Pr(>F)`[2])
    
    length(ests)
    length(Fvals)
    length(pvals)
  }
}

#### iii. merging data ####
est<- c()
Fval<-c()
pval<- c()

est<-c(round(as.numeric(paste(estt)), 3),round(as.numeric(paste(ests)), 3))
Fval<- c(round(Fvalt, 3), round(Fvals,3))
pval<- c(round(pvalt, 3), round(pvals,3))

LMresults[,6]<- est
LMresults[,8]<- Fval
LMresults[,9]<- pval

LMresults<- as.data.frame(LMresults)

save(LMresults, file="./OutputData/LMresults.RData")
write.table(LMresults, file="./OutputData/LMresults.txt")

## e plotting ####
### 0. data ####

# option 1: use raw data and continue with data preparation
#raw data - NB: use prepared data
load("./OutputData/LMresults.RData")
load("./InputData/tplants.RData")
load("./InputData/splants.RData")
load("./InputData/therb.RData")
load("./InputData/sherb.RData")
load("./InputData/tpred.RData")
load("./InputData/spred.RData")

# option 2: use prepared data and continue with plotting
#prepared data
load("./OutputData/ptlong.RData")
load("./OutputData/pslong.RData")

load("./OutputData/htlong.RData")
load("./OutputData/hslong.RData")

load("./OutputData/sctlong.RData")
load("./OutputData/scslong.RData")

### 00. code examples ####
ggplot(df_long, aes(x1, value, color = variable)) +
  geom_point() +
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1)

p1bs<- ggplot(splines_plants_space_sc, aes(x=LUIx, y=LUIy, colour=beta_type))+
  geom_line(linewidth=1.5)+
  scale_y_sqrt(limits = c(0, 3), expand = c(0, 0))+
  scale_colour_hue(name="Type of Beta Diversity",                     
                   labels = c("bsim", "b0", "b1", "b2", "b3", "b4"), 
                   l=50)+
  xlab("LUI (scaled to year average)")+
  ylab("Effect of LUI differences on beta diversity")+
  ggtitle("c. Spatial LUI effects on beta diversity")+
  #geom_text(aes(x = 3.9, y = 0.45, label = "EV% = 0.1"), color="black") + 
  theme_classic()+
  theme(text=element_text(size = 20))

x <- "USA"
paste("location =", x)
#[1] "location = USA"

### 000. data prepatation ####

# long format
library(tidyr)

#plants
pslong <- gather(splants, "abundance_weighting", "diversity", 
                 a0.1, a1.1, a2.1, a3.1, a4.1)
pslong$abundance_weighting<- as.factor(pslong$abundance_weighting)

ptlong <- gather(tplants, "abundance_weighting", "diversity", 
                 a0.1, a1.1, a2.1, a3.1, a4.1)
ptlong$abundance_weighting<- as.factor(ptlong$abundance_weighting)

#ins. herbivores
names(sherb)
hslong <- gather(sherb, "abundance_weighting", "diversity", 
                 ha0.1, ha1.1, ha2.1, ha3.1, ha4.1)
hslong$abundance_weighting<- as.factor(hslong$abundance_weighting)

htlong <- gather(therb, "abundance_weighting", "diversity", 
                 ha0.1, ha1.1, ha2.1, ha3.1, ha4.1)
htlong$abundance_weighting<- as.factor(htlong$abundance_weighting)

#ins. sec. consumers
scslong <- gather(spred, "abundance_weighting", "diversity", 
                  pa0.1, pa1.1, pa2.1, pa3.1, pa4.1)
scslong$abundance_weighting<- as.factor(scslong$abundance_weighting)

sctlong <- gather(tpred, "abundance_weighting", "diversity", 
                  pa0.1, pa1.1, pa2.1, pa3.1, pa4.1)
sctlong$abundance_weighting<- as.factor(sctlong$abundance_weighting)

save(ptlong, file="./OutputData/ptlong.RData")
save(pslong, file="./OutputData/pslong.RData")

save(htlong, file="./OutputData/htlong.RData")
save(hslong, file="./OutputData/hslong.RData")

save(sctlong, file="./OutputData/sctlong.RData")
save(scslong, file="./OutputData/scslong.RData")

#layout
# 8 panel plot (cols: space, time; rows: LUI+LU components; col: hill numbers)

### 1 plants ####
names(pslong)
names(ptlong)

LU<- c("LUI", "MOW", "GRA", "FER")

plants_plot_list = list()

for(i in 1:4){
  
  x<- names(pslong[i+3])
  z<- LU[i]
  
  data<- data.frame(pslong[,c((3+i),18,19)])
  names(data)<- c("LU",colnames(pslong)[c(18,19)])
  
  pls<- ggplot(data, aes(x=LU, y=diversity, colour=abundance_weighting))+
    geom_point(shape=19, alpha=0.4)+
    stat_smooth(method = "lm", formula = y ~ x + I(x^2), 
                linewidth = 1, se=FALSE)+
    scale_colour_hue(name="Type of Diversity",                     
                     labels = c("a0", "a1", "a2", "a3", "a4"), 
                     l=50)+
    xlab(paste(z))+
    ylab("alpha diversity")+
    ggtitle(paste("Spatial ",z,"effects on alpha diversity"))+
    #geom_text(aes(x = 3.9, y = 0.45, label = "EV% = 0.1"), color="black") + 
    theme_classic()+
    theme(text=element_text(size = 15))
  
  plants_plot_list[[i]]<- pls
  
  rm(data)
}

for(i in 1:4){
  
  x<- names(ptlong[i+7])
  z<- LU[i]
  
  data<- data.frame(ptlong[,c((7+i),22,23)])
  names(data)<- c("LU",colnames(ptlong)[c(22,23)])
  
  plt<- ggplot(data, aes(x=LU, y=diversity, colour=abundance_weighting))+
    geom_point(shape=19, alpha=0.4)+
    stat_smooth(method = "lm", formula = y ~ x + I(x^2), 
                linewidth = 1, se=FALSE)+
    scale_colour_hue(name="Type of Diversity",                     
                     labels = c("a0", "a1", "a2", "a3", "a4"), 
                     l=50)+
    xlab(paste(z))+
    ylab("alpha diversity")+
    ggtitle(paste("Temporal ",z,"effects on alpha diversity"))+
    #geom_text(aes(x = 3.9, y = 0.45, label = "EV% = 0.1"), color="black") + 
    theme_classic()+
    theme(text=element_text(size = 15))
  
  plants_plot_list[[i+4]]<- plt
  rm(data)
}
### 2 ins. herbivores ####
names(hslong)
names(htlong)

LU<- c("LUI", "MOW", "GRA", "FER")

herb_plot_list = list()

for(i in 1:4){
  
  x<- names(hslong[i+3])
  z<- LU[i]
  
  data<- data.frame(hslong[,c((3+i),18,19)])
  names(data)<- c("LU",colnames(hslong)[c(18,19)])
  
  pls<- ggplot(data, aes(x=LU, y=diversity, colour=abundance_weighting))+
    geom_point(shape=19, alpha=0.4)+
    stat_smooth(method = "lm", formula = y ~ x + I(x^2), 
                linewidth = 1, se=FALSE)+
    scale_colour_hue(name="Type of Diversity",                     
                     labels = c("a0", "a1", "a2", "a3", "a4"), 
                     l=50)+
    xlab(paste(z))+
    ylab("alpha diversity")+
    ggtitle(paste("Spatial ",z,"effects on alpha diversity"))+
    #geom_text(aes(x = 3.9, y = 0.45, label = "EV% = 0.1"), color="black") + 
    theme_classic()+
    theme(text=element_text(size = 15))
  
  herb_plot_list[[i]]<- pls
  rm(data)
}

for(i in 1:4){
  
  x<- names(htlong[i+7])
  z<- LU[i]
  
  data<- data.frame(htlong[,c((7+i),22,23)])
  names(data)<- c("LU",colnames(htlong)[c(22,23)])
  
  plt<- ggplot(data, aes(x=LU, y=diversity, colour=abundance_weighting))+
    geom_point(shape=19, alpha=0.4)+
    stat_smooth(method = "lm", formula = y ~ x + I(x^2), 
                linewidth = 1, se=FALSE)+
    scale_colour_hue(name="Type of Diversity",                     
                     labels = c("a0", "a1", "a2", "a3", "a4"), 
                     l=50)+
    xlab(paste(z))+
    ylab("alpha diversity")+
    ggtitle(paste("Temporal ",z,"effects on alpha diversity"))+
    #geom_text(aes(x = 3.9, y = 0.45, label = "EV% = 0.1"), color="black") + 
    theme_classic()+
    theme(text=element_text(size = 15))
  
  herb_plot_list[[i+4]]<- plt
  rm(data)
}
### 3 ins. sec. consumers ####

names(scslong)
names(sctlong)

LU<- c("LUI", "MOW", "GRA", "FER")

pred_plot_list = list()

for(i in 1:4){
  
  x<- names(scslong[i+3])
  z<- LU[i]
  
  data<- data.frame(scslong[,c((3+i),18,19)])
  names(data)<- c("LU",colnames(scslong)[c(18,19)])
  
  pls<- ggplot(data, aes(x=LU, y=diversity, colour=abundance_weighting))+
    geom_point(shape=19, alpha=0.4)+
    stat_smooth(method = "lm", formula = y ~ x + I(x^2), 
                linewidth = 1, se=FALSE)+
    scale_colour_hue(name="Type of Diversity",                     
                     labels = c("a0", "a1", "a2", "a3", "a4"), 
                     l=50)+
    xlab(paste(z))+
    ylab("alpha diversity")+
    ggtitle(paste("Spatial ",z,"effects on alpha diversity"))+
    #geom_text(aes(x = 3.9, y = 0.45, label = "EV% = 0.1"), color="black") + 
    theme_classic()+
    theme(text=element_text(size = 15))
  
  pred_plot_list[[i]]<- pls
  rm(data)
}

for(i in 1:4){
  
  x<- names(sctlong[i+7])
  z<- LU[i]
  
  data<- data.frame(sctlong[,c((7+i),22,23)])
  names(data)<- c("LU",colnames(sctlong)[c(22,23)])
  
  plt<- ggplot(data, aes(x=LU, y=diversity, colour=abundance_weighting))+
    geom_point(shape=19, alpha=0.4)+
    stat_smooth(method = "lm", formula = y ~ x + I(x^2), 
                linewidth = 1, se=FALSE)+
    scale_colour_hue(name="Type of Diversity",                     
                     labels = c("a0", "a1", "a2", "a3", "a4"), 
                     l=50)+
    xlab(paste(z))+
    ylab("alpha diversity")+
    ggtitle(paste("Temporal ",z,"effects on alpha diversity"))+
    #geom_text(aes(x = 3.9, y = 0.45, label = "EV% = 0.1"), color="black") + 
    theme_classic()+
    theme(text=element_text(size = 15))
  
  pred_plot_list[[i+4]]<- plt
  rm(data)
}

names_list<- c("LUI.space", "MOW.space", "GRA.space", "FER.space",
               "LUI.time", "MOW.time", "GRA.time", "FER.time")

names(plants_plot_list)<- names_list
names(herb_plot_list)<- names_list
names(pred_plot_list)<- names_list

save(plants_plot_list, file="./OutputData/plants_plot_list.RData")
save(herb_plot_list, file="./OutputData/herb_plot_list.RData")
save(pred_plot_list, file="./OutputData/pred_plot_list.RData")

library(plyr)
library(ggplot2)

li = structure(plants_plot_list, class = c("gglist", "ggplot"))
print.gglist = function(x, ...) l_ply(x, print, ...)
ggsave(li, file = "plant_div.pdf")

li = structure(herb_plot_list, class = c("gglist", "ggplot"))
print.gglist = function(x, ...) l_ply(x, print, ...)
ggsave(li, file = "herb_div.pdf")

li = structure(pred_plot_list, class = c("gglist", "ggplot"))
print.gglist = function(x, ...) l_ply(x, print, ...)
ggsave(li, file = "pred_div.pdf")


# 3 Model significant results GDM ####

## a packages ####
library(car)
library(lme4)
library(emmeans)
library(stats)
library(Rmisc)
## b data ####
#table summarises model results and is provided prepared
tab<- read.csv("./InputData/SupMat_LMGDM_Summary.csv", 
               header=T, 
               sep=";",
               dec=".")


## c analysis ####

# idea: test hypotheses with beta GLM based on GDM results

### i. H1 ####

# same shape and direction = congruence for spatial and temporal models
# congruence == 1

#### shape ####
modH1s<- glm((congruence2) ~ trophic.group*diversity.response + trophic.group*measure.of.landuse.intensity,
             family="quasibinomial",
             data=tab[tab$relationship.aspect=="shape",])

Anova(modH1s, type="III")
summary(modH1s)

modH1s0<- glm((congruence2) ~ diversity.response,
              family="quasibinomial",
              data=tab[tab$relationship.aspect=="shape",])


emm0<- emmeans(modH1s0, specs = pairwise ~ diversity.response, type="response")
emm0$emmeans

emm1<- emmeans(modH1s, specs = pairwise ~ trophic.group|diversity.response, type="response")
emm1$emmeans

#### direction mean ####
modH1dm<- glm(congruence2 ~ trophic.group*diversity.response + trophic.group*measure.of.landuse.intensity,
              family="quasibinomial", data=tab[tab$relationship.aspect=="direction_mean",])

modH1dmx<- glm(congruence ~ trophic.group*diversity.response + trophic.group*measure.of.landuse.intensity,
               family="quasibinomial", data=tab[tab$relationship.aspect=="direction_mean",])

Anova(modH1dm, type="III")

modH1dm0<- glm(congruence2 ~ diversity.response,
               family="quasibinomial", data=tab[tab$relationship.aspect=="direction_mean",])


#summary(modH1dm)
emm0<- emmeans(modH1dm0, specs = pairwise ~ diversity.response, type="response")
emm0$emmeans


emm1 = emmeans(modH1dm, specs = pairwise ~ trophic.group|diversity.response, type="response")
emm1$emmeans

emm2<- emmeans(modH1dm, specs = pairwise ~ trophic.group|measure.of.landuse.intensity, type="response")
emm2$emmeans

emm3x<- emmeans(modH1dmx, specs = pairwise ~ trophic.group|diversity.response, type="response")
emm3x$emmeans

emm4<- emmeans(modH1dm, specs = pairwise ~ diversity.response, type="response")
emm4$emmeans

#### direction delta ####
modH1dd<- glm(congruence2 ~ trophic.group*diversity.response + trophic.group*measure.of.landuse.intensity,
              family="quasibinomial", data=tab[tab$relationship.aspect=="direction_delta",])

modH1ddx<- glm(congruence ~ trophic.group*diversity.response + trophic.group*measure.of.landuse.intensity,
               family="quasibinomial", data=tab[tab$relationship.aspect=="direction_delta",])

Anova(modH1dd, type="III")

modH1dd0<- glm(congruence2 ~ diversity.response,
               family="quasibinomial", data=tab[tab$relationship.aspect=="direction_delta",])



emm0 = emmeans(modH1dd0, specs = pairwise ~ diversity.response, type="response")
emm0$emmeans


emm1 = emmeans(modH1dd, specs = pairwise ~ trophic.group|diversity.response, type="response")
emm1$emmeans

emm2<- emmeans(modH1dd, specs = pairwise ~ trophic.group|measure.of.landuse.intensity, type="response")
emm2$emmeans

emm3x<- emmeans(modH1ddx, specs = pairwise ~ trophic.group|diversity.response, type="response")
emm3x$emmeans

emm4<- emmeans(modH1dd, specs = pairwise ~ diversity.response, type="response")
emm4$emmeans

### ii. H2/H3 ####

#### effect size ####
# lower effect size in temporal vs spatial models
# and for plants vs animals


modH23es<- lm(congruence2 ~ trophic.group*diversity.response + trophic.group*measure.of.landuse.intensity,
              data=tab[tab$relationship.aspect=="effectsize" & tab$congruence2>-10,])
Anova(modH23es, type="III")
summary(modH23es)

modH23es0<- lm(congruence2 ~ diversity.response,
               data=tab[tab$relationship.aspect=="effectsize" & tab$congruence2>-10,])

#summary
emm0<- emmeans(modH23es0, specs = pairwise ~ diversity.response, type="response")
emm0$emmeans

emm1 = emmeans(modH23es, specs = pairwise ~ trophic.group|diversity.response, type="response")
emm1$emmeans

emm2<- emmeans(modH23es, specs = pairwise ~ trophic.group|measure.of.landuse.intensity, type="response")
emm2$emmeans


#### significance ####
modH23sig<- glm(congruence2 ~ trophic.group*diversity.response + trophic.group*measure.of.landuse.intensity,
                family="quasibinomial", data=tab[tab$relationship.aspect=="significance",])
Anova(modH23sig, type="III")
summary(modH23sig)

modH23sig0<- glm(congruence2 ~ diversity.response,
                 family="quasibinomial", data=tab[tab$relationship.aspect=="significance",])

#summary

emm0<- emmeans(modH23sig0, specs = pairwise ~ diversity.response, type="response")
emm0$emmeans

emm1 = emmeans(modH23sig, specs = pairwise ~ trophic.group|diversity.response, type="response")
emm1$emmeans

emm2<- emmeans(modH23sig, specs = pairwise ~ diversity.response|measure.of.landuse.intensity, type="response")
emm2$emmeans


# 4 Congruence plots ####

## a data ####

# data upload ####
# option 1: use raw data and assemble datasets yourself
# continue with "data prep"

# option 2: use prepared data sets "tab_agg1/2/2a/3 
# and continue with uploading data below and then plotting

# tab_agg1
load("./InputData/tab_agg1.RData")

tab_agg1$id1<- factor(tab_agg1$id1, levels=c( "arthropod secondary consumers",
                                              "arthropod herbivores",
                                              "plants"))
tab_agg1$id4<- factor(tab_agg1$id4, levels=c("beta", "alpha"))
tab_agg1$id3<- factor(tab_agg1$id3, levels=c("FER", "GRA", "MOW", "LUI"))
tab_agg1$int<- as.factor(interaction(tab_agg1$id3,
                                     tab_agg1$id1,
                                     tab_agg1$id4))
# tab_agg2
load("./InputData/tab_agg2.RData")

tab_agg2$id2<- factor(tab_agg2$id2, levels=c( "significance",
                                              "effectsize",
                                              "direction_delta",
                                              "direction_mean",
                                              "shape"))
tab_agg2$id4<- factor(tab_agg2$id4, levels=c("beta", "alpha"))

tab_agg2$int<- as.factor(interaction(tab_agg2$id2, tab_agg2$id3))

# tab_agg2a
load("./InputData/tab_agg2a.RData")

tab_agg2a$id1<- factor(tab_agg2a$id1, levels=c( "arthropod secondary consumers",
                                                "arthropod herbivores",
                                                "plants"))
tab_agg2a$id4<- factor(tab_agg2a$id4, levels=c("FER", "GRA", "MOW", "LUI"))
tab_agg2a$int<- as.factor(interaction(tab_agg2a$id1,
                                      tab_agg2a$id4))

# tab_agg3
load("./InputData/tab_agg3.RData")

tab_agg3$id2<- factor(tab_agg3$id2, levels=c( "significance",
                                              "effectsize",
                                              "direction_delta",
                                              "direction_mean",
                                              "shape"))
tab_agg3$id3<- factor(tab_agg3$id3, levels=c("beta", "alpha"))

tab_agg3$int<- as.factor(interaction(tab_agg3$id2, tab_agg3$id3))

# data prep ####

# means over abundance weighting (conservative)
## all combinations ####
tab_agg1<- aggregate(tab$congruence2[tab$congruence2>-10], 
                     list(id1 = tab$trophic.group[tab$congruence2>-10],
                          id2 = tab$relationship.aspect[tab$congruence2>-10],
                          id3 = tab$measure.of.landuse.intensity[tab$congruence2>-10],
                          id4 = tab$diversity.response[tab$congruence2>-10]),
                     FUN=mean)

head(tab_agg1)

tab_agg1$id1<- factor(tab_agg1$id1, levels=c( "arthropod secondary consumers",
                                              "arthropod herbivores",
                                              "plants"))
tab_agg1$id4<- factor(tab_agg1$id4, levels=c("beta", "alpha"))
tab_agg1$id3<- factor(tab_agg1$id3, levels=c("FER", "GRA", "MOW", "LUI"))
tab_agg1$int<- as.factor(interaction(tab_agg1$id3,
                                     tab_agg1$id1,
                                     tab_agg1$id4))

# add manually the means and CIs from the emmeans objects as they are based on the
# data distribution

write.table(tab_agg1, file="./InputData/tab_agg1.txt")
tab_agg1<- read.csv("tab_agg1_ed.csv", 
                    header=T, 
                    sep=";",
                    dec=".")

save(tab_agg1, file="./InputData/tab_agg1.RData")
#less conservative congruence
# see dir-mean/d x in tab_agg


## across all LU Types ####
tab_agg2<- aggregate(tab$congruence2[tab$congruence2>-10], 
                     list(id1 = tab$trophic.group[tab$congruence2>-10],
                          id2 = tab$relationship.aspect[tab$congruence2>-10],
                          id4 = tab$diversity.response[tab$congruence2>-10]),
                     FUN=mean)

head(tab_agg2)

tab_agg2$id1<- factor(tab_agg2$id1, levels=c( "arthropod secondary consumers",
                                              "arthropod herbivores",
                                              "plants"))
tab_agg2$id4<- factor(tab_agg2$id4, levels=c("beta", "alpha"))
tab_agg2$int<- as.factor(interaction(tab_agg2$id1,
                                     tab_agg2$id4))

# add manually the means and CIs from the emmeans objects as they are based on the
# data distribution

write.table(tab_agg2, file="./InputData/tab_agg2.txt")
tab_agg2<- read.csv("tab_agg2_ed.csv", 
                    header=T, 
                    sep=";",
                    dec=".")

save(tab_agg2, file="./InputData/tab_agg2.RData")
#less conservative congruence
# see dir-mean/d x in tab_agg

## across all div types ####
tab_agg2a<- aggregate(tab$congruence2[tab$congruence2>-10], 
                      list(id1 = tab$trophic.group[tab$congruence2>-10],
                           id2 = tab$relationship.aspect[tab$congruence2>-10],
                           id4 = tab$measure.of.landuse.intensity[tab$congruence2>-10]),
                      FUN=mean)

head(tab_agg2a)

tab_agg2a$id1<- factor(tab_agg2a$id1, levels=c( "arthropod secondary consumers",
                                                "arthropod herbivores",
                                                "plants"))
tab_agg2a$id4<- factor(tab_agg2a$id4, levels=c("FER", "GRA", "MOW", "LUI"))
tab_agg2a$int<- as.factor(interaction(tab_agg2a$id1,
                                      tab_agg2a$id4))

# add manually the means and CIs from the emmeans objects as they are based on the
# data distribution

write.table(tab_agg2a, file="./InputData/tab_agg2a.txt")
tab_agg2a<- read.csv("tab_agg2a_ed.csv", 
                     header=T, 
                     sep=";",
                     dec=".")

save(tab_agg2a, file="./InputData/tab_agg2a.RData")

## accross all LU Types & organims ####
tab_agg3<- aggregate(tab$congruence2[tab$congruence2>-10], 
                     list(id2 = tab$relationship.aspect[tab$congruence2>-10],
                          id3 = tab$diversity.response[tab$congruence2>-10]))



tab_agg3$id2<- factor(tab_agg3$id2, levels=c( "significance",
                                              "effectsize",
                                              "direction_delta",
                                              "direction_mean",
                                              "shape"))
tab_agg3$id3<- factor(tab_agg3$id3, levels=c("beta", "alpha"))

tab_agg3$int<- as.factor(interaction(tab_agg3$id2, tab_agg3$id3))

# add manually the means and CIs from the emmeans objects as they are based on the
# data distribution

write.table(tab_agg3, file="./InputData/tab_agg3.txt")
tab_agg3<- read.csv("tab_agg3_ed.csv", 
                    header=T, 
                    sep=";",
                    dec=".")

save(tab_agg3, file="./InputData/tab_agg3.RData")

# plotting ####
# with data aggregated accross LU Types to have more
# replicates and more precise CIs

# packages ####
library(ggplot2)

## plot subsets ####
shape<- tab_agg2[tab_agg2$id2=="shape",]
dir_del<- tab_agg2[tab_agg2$id2=="direction_delta",]
dir_mea<- tab_agg2[tab_agg2$id2=="direction_mean",]
sig<- tab_agg2[tab_agg2$id2=="significance",]
es<- tab_agg2[tab_agg2$id2=="effectsize",]

dir_delx<- tab_agg2[tab_agg2$id2=="direction_deltax",]
dir_meax<- tab_agg2[tab_agg2$id2=="direction_meanx",]

lues<- tab_agg2a[tab_agg2a$id2=="effectsize",]
ludirmean<- tab_agg2a[tab_agg2a$id2=="direction_mean",]

yed2<- t(c("significance",
           "effect size",
           "direction delta",
           "direction mean",
           "shape",
           "significance",
           "effect size",
           "direction delta",
           "direction mean",
           "shape"))

og<- t(c("sec. consumers",
         "herbivores",
         "plants",
         "sec. consumers",
         "herbivores",
         "plants"))

## 0 - summary ####

summary_plot<- ggplot(tab_agg3, aes(prob, int, col=id3, size = id3, shape=id3))+
  geom_pointrange(aes(xmin=asyCI_lower, xmax=asyCI_upper)) +
  geom_point() +
  xlim(0,1)+
  xlab("average congruence %")+
  ylab("")+
  scale_color_manual(values=c("black","grey"))+
  scale_size_manual(values=c(1,1))+
  scale_shape_manual(values=c(15,17))+
  guides(alpha = "none")+
  scale_y_discrete(labels=yed2)+
  labs(title= "(a) overall congruence")+
  geom_hline(yintercept=c(5.5), col="grey") +
  geom_vline(xintercept=c(0,1), linetype = "dotted")+
  theme_classic()+
  theme(text=element_text(size = 15))
summary_plot


## 1 - shape ####
shape_plot<- ggplot(shape, aes(prob, int, col=id1, shape=id4, size=id4))+
  geom_pointrange(aes(xmin=asyCI_lower, xmax=asyCI_upper)) +
  geom_point() +
  xlim(0,1)+
  xlab("average congruence %")+
  ylab("")+
  scale_color_manual(values=c("#7570b3","#d95f02", "#1b9e77"))+
  guides(alpha = "none")+
  scale_y_discrete(labels=og)+
  scale_size_manual(values=c(1,1))+
  scale_shape_manual(values=c(15,17))+
  labs(title= "(b) relationship shape",
       subtitle="organism: ns - diversity type: ns - org x div:**")+
  geom_hline(yintercept=c(12.5), col="grey") +
  geom_hline(yintercept=c(3.5), col="grey", linetype="dotted") +
  geom_vline(xintercept=c(0,1), linetype = "dotted")+
  theme_classic()+
  theme(text=element_text(size = 15))
shape_plot


## 2 - direction,mean ####
drm_plot<- ggplot(dir_mea, aes(prob, int, col=id1, shape=id4, size=id4))+
  geom_pointrange(aes(xmin=asyCI_lower, xmax=asyCI_upper)) +
  geom_point() +
  xlim(0,1)+
  xlab("average congruence %")+
  ylab("")+
  scale_color_manual(values=c("#7570b3","#d95f02", "#1b9e77"))+
  guides(alpha = "none")+
  scale_y_discrete(labels=og)+
  scale_size_manual(values=c(1,1))+
  scale_shape_manual(values=c(15,17))+
  labs(title= "(c) relationship direction (mean)",
       subtitle="organism: ns - diversity type: ns - org x div: ns")+
  geom_hline(yintercept=c(12.5), col="grey") +
  geom_hline(yintercept=c(3.5), col="grey", linetype="dotted") +
  geom_vline(xintercept=c(0,1), linetype = "dotted")+
  theme_classic()+
  theme(text=element_text(size = 15))
drm_plot

## 3 - direction,delta ####
drd_plot<- ggplot(dir_del, aes(prob, int, col=id1, shape=id4, size=id4))+
  geom_pointrange(aes(xmin=asyCI_lower, xmax=asyCI_upper)) +
  geom_point() +
  xlim(0,1)+
  xlab("average congruence %")+
  ylab("")+
  scale_color_manual(values=c("#7570b3","#d95f02", "#1b9e77"))+
  guides(alpha = "none")+
  scale_y_discrete(labels=og)+
  scale_size_manual(values=c(1,1))+
  scale_shape_manual(values=c(15,17))+
  labs(title= "(d) relationship direction (delta)",
       subtitle="organism: ns - diversity type: ns - org x div:**")+
  geom_hline(yintercept=c(12.5), col="grey") +
  geom_hline(yintercept=c(3.5), col="grey", linetype="dotted") +
  geom_vline(xintercept=c(0,1), linetype = "dotted")+
  theme_classic()+
  theme(text=element_text(size = 15))
drd_plot


## 4 - significance ####
sig_plot<- ggplot(sig, aes(prob, int, col=id1, shape=id4, size=id4))+
  geom_pointrange(aes(xmin=asyCI_lower, xmax=asyCI_upper)) +
  geom_point() +
  xlim(0,1)+
  xlab("average congruence %")+
  ylab("")+
  scale_color_manual(values=c("#7570b3","#d95f02", "#1b9e77"))+
  scale_y_discrete(labels=og)+
  scale_size_manual(values=c(1,1))+
  scale_shape_manual(values=c(15,17))+
  labs(title= "(f) relationship significance",
       subtitle="organism: ns - diversity type: ns - org x div:**")+
  geom_hline(yintercept=c(12.5), col="grey") +
  geom_hline(yintercept=c(3.5), col="grey", linetype="dotted") +
  geom_vline(xintercept=c(0,1), linetype = "dotted")+
  theme_classic()+
  theme(text=element_text(size = 15))
sig_plot


## 5 - effect size ####
es_plot<- ggplot(es, aes(prob, int, col=id1, shape=id4, size=id4))+
  geom_pointrange(aes(xmin=asyCI_lower, xmax=asyCI_upper)) +
  geom_point() +
  xlim(-0.5,1)+
  xlab("average congruence %")+
  ylab("")+
  scale_color_manual(values=c("#7570b3","#d95f02", "#1b9e77"))+
  guides(alpha = "none")+
  scale_y_discrete(labels=og)+
  scale_size_manual(values=c(1,1))+
  scale_shape_manual(values=c(15,17))+
  labs(title= "(e) relationship strength",
       subtitle="organism:* - diversity type:* - org x div:***")+
  geom_hline(yintercept=c(12.5), col="grey") +
  geom_hline(yintercept=c(3.5), col="grey", linetype="dotted") +
  geom_vline(xintercept=c(0,1), linetype = "dotted")+
  theme_classic()+
  theme(text=element_text(size = 15))
es_plot


# less conservative

## 2* - direction,mean ####
drmx_plot<- ggplot(dir_meax, aes(prob, int, col=id1, shape=id4, size=id4))+
  geom_pointrange(aes(xmin=asyCI_lower, xmax=asyCI_upper)) +
  geom_point() +
  xlim(0,1)+
  xlab("average congruence %")+
  ylab("")+
  scale_color_manual(values=c("#7570b3","#d95f02", "#1b9e77"))+
  guides(alpha = "none")+
  scale_y_discrete(labels=og)+
  scale_size_manual(values=c(1,1))+
  scale_shape_manual(values=c(15,17))+
  labs(title= "(a) relationship direction (mean*)")+
  geom_hline(yintercept=c(12.5), col="grey") +
  geom_hline(yintercept=c(3.5), col="grey", linetype="dotted") +
  geom_vline(xintercept=c(0,1), linetype = "dotted")+
  theme_classic()+
  theme(text=element_text(size = 20))
drmx_plot

## 3* - direction,delta ####
drdx_plot<- ggplot(dir_delx, aes(prob, int, col=id1, shape=id4, size=id4))+
  geom_pointrange(aes(xmin=asyCI_lower, xmax=asyCI_upper)) +
  geom_point() +
  xlim(0,1)+
  xlab("average congruence %")+
  ylab("")+
  scale_color_manual(values=c("#7570b3","#d95f02", "#1b9e77"))+
  guides(alpha = "none")+
  scale_y_discrete(labels=og)+
  scale_size_manual(values=c(1,1))+
  scale_shape_manual(values=c(15,17))+
  labs(title= "(b) relationship direction (delta*)")+
  geom_hline(yintercept=c(12.5), col="grey") +
  geom_hline(yintercept=c(3.5), col="grey", linetype="dotted") +
  geom_vline(xintercept=c(0,1), linetype = "dotted")+
  theme_classic()+
  theme(text=element_text(size = 20))
drdx_plot


## 6 - LU (effect size) ####
lues_plot<- ggplot(lues, aes(prob, int, col=id1, shape=id4, size=id4))+
  geom_pointrange(aes(xmin=asyCI_lower, xmax=asyCI_upper)) +
  geom_point() +
  xlim(-0.75,1.5)+
  xlab("average deviation (spatial-temporal,%)")+
  ylab("")+
  scale_color_manual(values=c("#7570b3","#d95f02", "#1b9e77"))+
  guides(alpha = "none")+
  #scale_y_discrete(labels=og)+
  scale_size_manual(values=c(1,1,1,1))+
  labs(title= "(b) relationship strength")+
  geom_hline(yintercept=c(12.5), col="grey") +
  geom_hline(yintercept=c(3.5, 6.5, 9.5), col="grey", linetype="dotted") +
  geom_vline(xintercept=c(0,1), linetype = "dotted")+
  theme_classic()+
  theme(text=element_text(size = 20))
lues_plot

## 7 - LU (direction mean) ####
ludm_plot<- ggplot(ludirmean, aes(prob, int, col=id1, shape=id4, size=id4))+
  geom_pointrange(aes(xmin=asyCI_lower, xmax=asyCI_upper)) +
  geom_point() +
  xlim(0,1)+
  xlab("average congruence %")+
  ylab("")+
  scale_color_manual(values=c("#7570b3","#d95f02", "#1b9e77"))+
  guides(alpha = "none")+
  #scale_y_discrete(labels=og)+
  scale_size_manual(values=c(1,1,1,1))+
  labs(title= "(a) relationship direction (mean)")+
  geom_hline(yintercept=c(12.5), col="grey") +
  geom_hline(yintercept=c(3.5, 6.5, 9.5), col="grey", linetype="dotted") +
  geom_vline(xintercept=c(0,1), linetype = "dotted")+
  theme_classic()+
  theme(text=element_text(size = 20))
ludm_plot
