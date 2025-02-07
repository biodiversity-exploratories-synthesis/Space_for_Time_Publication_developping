# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
######## CHAO Numbers Alpha Beta                                                ----
########  This script calculates the Chao numbers both alpha and beta-diversity ###
########  as well as delta alpha


### It starts from the complete pairwise data sets and adds the species abundance tables
### TODO @ Noelle: how to reduce the final pairwise datasets to only those lines that contain
### diversity indices? E.g. the insect data set has less plots and years and thus
## comparisons.

#load package
# install.packages("hillR")
library(hillR)

# install.packages("betapart")
library(betapart)

library(vegan)

#load data

#pairwise datasets as created from github

#set working directory to folder "step 2- Alpha_Beta_calculations_Chao"
#pwise_space
load("data/RawData/pwise_space.RData")
load("data/RawData/pwise_time.RData")

#community data
load("data/InputData/plants_rd.RData")
load("data/InputData/data/herb.RData")
load("data/InputData/pred.RData")


# 1. Plants ----
# 1.1. space ---- 
pwisedata<- pwise_space
com<- plants_rd
com_sub<- com[, c(6:372)]#cut out all non-species columns

# important but default: relative abundances are used, so function input is simple
a0_EP1<- c()
a1_EP1<- c()
a2_EP1<- c()
a3_EP1<- c()
a4_EP1<- c()

a0_EP2<- c()
a1_EP2<- c()
a2_EP2<- c()
a3_EP2<- c()
a4_EP2<- c()

#delta of alpha diversity 
da0<- c()
da1<- c()
da2<- c()
da3<- c()
da4<- c()

da0st<- c()
da1st<- c()
da2st<- c()
da3st<- c()
da4st<- c()

bsim<-c()
b0<- c()
b1<- c()
b2<- c()
b3<- c()
b4<- c()

for(i in 1:nrow(pwisedata)){
  
  com1<- com_sub[com$EPy == pwisedata$EPy1[i],]
  com2<- com_sub[com$EPy == pwisedata$EPy2[i],]
  
  com12<- as.matrix(rbind(com1,com2))
  
  
  #how to make sure that empty rows in the abundance data sets are skipped
  
  if(nrow(com1)< 1 || nrow(com2)<1 || rowSums(com1)<1 || rowSums(com2)<1){
    a0_EP1<- c(a0_EP1, NA)
    a0_EP2<- c(a0_EP2, NA)
    a1_EP1<- c(a1_EP1, NA)
    a1_EP2<- c(a1_EP2, NA)
    a2_EP1<- c(a2_EP1, NA)
    a2_EP2<- c(a2_EP2, NA)
    a3_EP1<- c(a3_EP1, NA)
    a3_EP2<- c(a3_EP2, NA)
    a4_EP1<- c(a4_EP1, NA)
    a4_EP2<- c(a4_EP2, NA)
    
    da0<- c(da0, NA)
    da1<- c(da1, NA)
    da2<- c(da2, NA)
    da3<- c(da3, NA)
    da4<- c(da4, NA)
    
    da0st<- c(da0st, NA)
    da1st<- c(da1st, NA)
    da2st<- c(da2st, NA)
    da3st<- c(da3st, NA)
    da4st<- c(da4st, NA)
    
    bsim<- c(bsim,NA)
    b0<- c(b0,NA)
    b1<- c(b1, NA)
    b2<- c(b2, NA)
    b3<- c(b3, NA)
    b4<- c(b4, NA)
  
  }else {
    
    com12<- com12[,colSums(com12)>0]
    
    alpha0<- hill_taxa(comm=com12, q=0)
    alpha1<- hill_taxa(comm=com12, q=1)
    alpha2<- hill_taxa(comm=com12, q=2)
    alpha3<- hill_taxa(comm=com12, q=3)
    alpha4<- hill_taxa(comm=com12, q=4)
    
    com12pa<- decostand(com12, method="pa")
    
    betasim<- beta.pair(com12pa, index.family="sorensen")
    
    beta0<- hill_taxa_parti(comm=com12, q=0)
    beta1<- hill_taxa_parti(comm=com12, q=1)    
    beta2<- hill_taxa_parti(comm=com12, q=2)    
    beta3<- hill_taxa_parti(comm=com12, q=3)    
    beta4<- hill_taxa_parti(comm=com12, q=4)    
    
    a0_EP1<- c(a0_EP1, alpha0[1])
    a0_EP2<- c(a0_EP2, alpha0[2])
    a1_EP1<- c(a1_EP1, alpha1[1])
    a1_EP2<- c(a1_EP2, alpha1[2])
    a2_EP1<- c(a2_EP1, alpha2[1])
    a2_EP2<- c(a2_EP2, alpha2[2])
    a3_EP1<- c(a3_EP1, alpha3[1])
    a3_EP2<- c(a3_EP2, alpha3[2])
    a4_EP1<- c(a4_EP1, alpha4[1])
    a4_EP2<- c(a4_EP2, alpha4[2])
    
    da0<- c(da0, (abs(alpha0[1]-alpha0[2])))
    da1<- c(da1, (abs(alpha1[1]-alpha1[2])))
    da2<- c(da2, (abs(alpha2[1]-alpha2[2])))
    da3<- c(da3, (abs(alpha3[1]-alpha3[2])))
    da4<- c(da4, (abs(alpha4[1]-alpha4[2])))
    
    da0st<- decostand(da0, method="range", na.rm=T)
    da1st<- decostand(da1, method="range", na.rm=T)
    da2st<- decostand(da2, method="range", na.rm=T)
    da3st<- decostand(da3, method="range", na.rm=T)
    da4st<- decostand(da4, method="range", na.rm=T)
    
    bsim<- c(bsim,betasim$beta.sim)
    b0<- c(b0, beta0$TD_beta)
    b1<- c(b1, beta1$TD_beta)
    b2<- c(b2, beta2$TD_beta)
    b3<- c(b3, beta3$TD_beta)
    b4<- c(b4, beta4$TD_beta)
    
  }
  
}

pwisedata_alphabeta<- data.frame(pwisedata,
                                 a0_EP1, a0_EP2,
                                 a1_EP1, a1_EP2,
                                 a2_EP1, a2_EP2,
                                 a3_EP1, a3_EP2,
                                 a4_EP1, a4_EP2,
                                 da0, da1, da2, da3, da4,
                                 da0st, da1st, da2st, da3st, da4st,
                                 bsim, b0, b1, b2, b3, b4
)

#reduce dataset to only those combination for which diversity estimates are there
# some of the combinations do not exist, because no data was collected in some
#   of the years. Exclude those cases.
div_names <- grep("a[0-9]_|bsim|[ab][0-4]", names(pwisedata_alphabeta), value = T)
# find rows with 21 NA values (all diversity indices are NA)
exclude_rows <- apply(pwisedata_alphabeta[, which(colnames(pwisedata_alphabeta) %in% div_names)],
      1, function(x) sum(is.na(x)) == 21)
# remove rows with NA in all diversity columns
pwisedata_cleaned <- pwisedata_alphabeta[!exclude_rows, ]
#
# double check 1: are there any NA values left?
pwisedata_cleaned[is.na(pwisedata_cleaned$b1), ]
# double check 2 : are there all rows which are present in the pwise space template?
nrow(pwise_space) >= nrow(pwisedata_cleaned)


#save as
pwise_space_plants<- pwisedata_cleaned

save(pwise_space_plants, file="data/OutputData/pwise_space_plants.RData")
# backup save as csv (large file)
write.table(pwise_space_plants, 
            file="data/OutputData/pwise_space_plants.csv",
            sep = ",", quote = F, row.names = F)

rm(pwisedata_alphabeta)
rm(pwisedata)
rm(pwisedata_cleaned)


#1.2. time ----
pwisedata<- pwise_time
com<- plants_rd
com_sub<- com[, c(6:372)]

#important but default: relative abundances are used, so function input is simple
a0_EP1<- c()
a1_EP1<- c()
a2_EP1<- c()
a3_EP1<- c()
a4_EP1<- c()

a0_EP2<- c()
a1_EP2<- c()
a2_EP2<- c()
a3_EP2<- c()
a4_EP2<- c()

#delta of alpha diversity 
da0<- c()
da1<- c()
da2<- c()
da3<- c()
da4<- c()

da0st<- c()
da1st<- c()
da2st<- c()
da3st<- c()
da4st<- c()

bsim<-c()
b0<- c()
b1<- c()
b2<- c()
b3<- c()
b4<- c()

for(i in 1:nrow(pwisedata)){
  
  com1<- com_sub[com$EPy == pwisedata$EPy1[i],]
  com2<- com_sub[com$EPy == pwisedata$EPy2[i],]
  
  com12<- as.matrix(rbind(com1,com2))
  
  #how to make sure that empty rows in the abundance data sets are skipped
  
  if(nrow(com1)< 1 || nrow(com2)<1 || rowSums(com1)<1 || rowSums(com2)<1){
    a0_EP1<- c(a0_EP1, NA)
    a0_EP2<- c(a0_EP2, NA)
    a1_EP1<- c(a1_EP1, NA)
    a1_EP2<- c(a1_EP2, NA)
    a2_EP1<- c(a2_EP1, NA)
    a2_EP2<- c(a2_EP2, NA)
    a3_EP1<- c(a3_EP1, NA)
    a3_EP2<- c(a3_EP2, NA)
    a4_EP1<- c(a4_EP1, NA)
    a4_EP2<- c(a4_EP2, NA)
    
    da0<- c(da0, NA)
    da1<- c(da1, NA)
    da2<- c(da2, NA)
    da3<- c(da3, NA)
    da4<- c(da4, NA)
    
    da0st<- c(da0st, NA)
    da1st<- c(da1st, NA)
    da2st<- c(da2st, NA)
    da3st<- c(da3st, NA)
    da4st<- c(da4st, NA)
    
    bsim<- c(bsim,NA)
    b0<- c(b0,NA)
    b1<- c(b1, NA)
    b2<- c(b2, NA)
    b3<- c(b3, NA)
    b4<- c(b4, NA)
    
  }else {
    
    com12<- com12[,colSums(com12)>0]
    
    alpha0<- hill_taxa(comm=com12, q=0)
    alpha1<- hill_taxa(comm=com12, q=1)
    alpha2<- hill_taxa(comm=com12, q=2)
    alpha3<- hill_taxa(comm=com12, q=3)
    alpha4<- hill_taxa(comm=com12, q=4)
    
    com12pa<- decostand(com12, method="pa")
    
    betasim<- beta.pair(com12pa, index.family="sorensen")
    
    beta0<- hill_taxa_parti(comm=com12, q=0)
    beta1<- hill_taxa_parti(comm=com12, q=1)    
    beta2<- hill_taxa_parti(comm=com12, q=2)    
    beta3<- hill_taxa_parti(comm=com12, q=3)    
    beta4<- hill_taxa_parti(comm=com12, q=4)    
    
    a0_EP1<- c(a0_EP1, alpha0[1])
    a0_EP2<- c(a0_EP2, alpha0[2])
    a1_EP1<- c(a1_EP1, alpha1[1])
    a1_EP2<- c(a1_EP2, alpha1[2])
    a2_EP1<- c(a2_EP1, alpha2[1])
    a2_EP2<- c(a2_EP2, alpha2[2])
    a3_EP1<- c(a3_EP1, alpha3[1])
    a3_EP2<- c(a3_EP2, alpha3[2])
    a4_EP1<- c(a4_EP1, alpha4[1])
    a4_EP2<- c(a4_EP2, alpha4[2])
    
    da0<- c(da0, (abs(alpha0[1]-alpha0[2])))
    da1<- c(da1, (abs(alpha1[1]-alpha1[2])))
    da2<- c(da2, (abs(alpha2[1]-alpha2[2])))
    da3<- c(da3, (abs(alpha3[1]-alpha3[2])))
    da4<- c(da4, (abs(alpha4[1]-alpha4[2])))
    
    da0st<- decostand(da0, method="range", na.rm=T)
    da1st<- decostand(da1, method="range", na.rm=T)
    da2st<- decostand(da2, method="range", na.rm=T)
    da3st<- decostand(da3, method="range", na.rm=T)
    da4st<- decostand(da4, method="range", na.rm=T)
    
    bsim<- c(bsim,betasim$beta.sim)
    b0<- c(b0, beta0$TD_beta)
    b1<- c(b1, beta1$TD_beta)
    b2<- c(b2, beta2$TD_beta)
    b3<- c(b3, beta3$TD_beta)
    b4<- c(b4, beta4$TD_beta)
    
  }
  
}

pwisedata_alphabeta<- data.frame(pwisedata,
                                 a0_EP1, a0_EP2,
                                 a1_EP1, a1_EP2,
                                 a2_EP1, a2_EP2,
                                 a3_EP1, a3_EP2,
                                 a4_EP1, a4_EP2,
                                 da0, da1, da2, da3, da4,
                                 da0st, da1st, da2st, da3st, da4st,
                                 bsim, b0, b1, b2, b3, b4
)

#reduce dataset to only those combination for which diversity estimates are there

# some of the combinations do not exist, because no data was collected in some
#   of the years. Exclude those cases.
div_names <- grep("a[0-9]_|bsim|[ab][0-4]", names(pwisedata_alphabeta), value = T)
# find rows with 26 NA values (all diversity indices are NA)
exclude_rows <- apply(pwisedata_alphabeta[, which(colnames(pwisedata_alphabeta) %in% div_names)],
                      1, function(x) sum(is.na(x)) == length(div_names))
# remove rows with NA in all diversity columns
pwisedata_cleaned <- pwisedata_alphabeta[!exclude_rows, ]


# save as
pwise_time_plants<- pwisedata_cleaned

save(pwise_time_plants, file="data/OutputData/pwise_time_plants.RData")
# backup save as csv (large file)
write.table(pwise_time_plants, 
            file="data/OutputData/pwise_time_plants.csv",
            sep = ",", quote = F, row.names = F)

rm(pwisedata_alphabeta)
rm(pwisedata)
rm(pwisedata_cleaned)

# 2. Insects ----

# 2.1. space ----

pwisedata <- pwise_space
hcom<- herb
pcom<- pred

#important but default: relative abundances are used, so function input is simple
# vectors herbivores
ha0_EP1<- c()
ha1_EP1<- c()
ha2_EP1<- c()
ha3_EP1<- c()
ha4_EP1<- c()

ha0_EP2<- c()
ha1_EP2<- c()
ha2_EP2<- c()
ha3_EP2<- c()
ha4_EP2<- c()

#delta of alpha diversity 
dha0<- c()
dha1<- c()
dha2<- c()
dha3<- c()
dha4<- c()

dha0st<- c()
dha1st<- c()
dha2st<- c()
dha3st<- c()
dha4st<- c()

hbsim<-c()
hb0<- c()
hb1<- c()
hb2<- c()
hb3<- c()
hb4<- c()

# vectors predators
pa0_EP1<- c()
pa1_EP1<- c()
pa2_EP1<- c()
pa3_EP1<- c()
pa4_EP1<- c()

pa0_EP2<- c()
pa1_EP2<- c()
pa2_EP2<- c()
pa3_EP2<- c()
pa4_EP2<- c()

#delta of alpha diversity 
dpa0<- c()
dpa1<- c()
dpa2<- c()
dpa3<- c()
dpa4<- c()

dpa0st<- c()
dpa1st<- c()
dpa2st<- c()
dpa3st<- c()
dpa4st<- c()

pbsim<-c()
pb0<- c()
pb1<- c()
pb2<- c()
pb3<- c()
pb4<- c()

for(i in 1:nrow(pwisedata)){
  hcom1<- hcom[rownames(hcom) == pwisedata$EPy1[i],]
  hcom2<- hcom[rownames(hcom) == pwisedata$EPy2[i],]
  
  hcom12<- as.matrix(rbind(hcom1,hcom2))
  
  
  #how to make sure that empty rows in the abundance data sets are skipped
  
  if(nrow(hcom1)< 1 || nrow(hcom2)<1 || rowSums(hcom1)<1 || rowSums(hcom2)<1){
    ha0_EP1<- c(ha0_EP1, NA)
    ha0_EP2<- c(ha0_EP2, NA)
    ha1_EP1<- c(ha1_EP1, NA)
    ha1_EP2<- c(ha1_EP2, NA)
    ha2_EP1<- c(ha2_EP1, NA)
    ha2_EP2<- c(ha2_EP2, NA)
    ha3_EP1<- c(ha3_EP1, NA)
    ha3_EP2<- c(ha3_EP2, NA)
    ha4_EP1<- c(ha4_EP1, NA)
    ha4_EP2<- c(ha4_EP2, NA)
    
    dha0<- c(dha0, NA)
    dha1<- c(dha1, NA)
    dha2<- c(dha2, NA)
    dha3<- c(dha3, NA)
    dha4<- c(dha4, NA)
    
    dha0st<- c(dha0st, NA)
    dha1st<- c(dha1st, NA)
    dha2st<- c(dha2st, NA)
    dha3st<- c(dha3st, NA)
    dha4st<- c(dha4st, NA)
    
    hbsim<- c(hbsim,NA)
    hb0<- c(hb0,NA)
    hb1<- c(hb1, NA)
    hb2<- c(hb2, NA)
    hb3<- c(hb3, NA)
    hb4<- c(hb4, NA)
    
  }else {
    
    hcom12<- hcom12[,colSums(hcom12)>0]
    
    halpha0<- hill_taxa(comm=hcom12, q=0)
    halpha1<- hill_taxa(comm=hcom12, q=1)
    halpha2<- hill_taxa(comm=hcom12, q=2)
    halpha3<- hill_taxa(comm=hcom12, q=3)
    halpha4<- hill_taxa(comm=hcom12, q=4)
    
    hcom12pa<- decostand(hcom12, method="pa")
    
    hbetasim<- beta.pair(hcom12pa, index.family="sorensen")
    
    hbeta0<- hill_taxa_parti(comm=hcom12, q=0)
    hbeta1<- hill_taxa_parti(comm=hcom12, q=1)    
    hbeta2<- hill_taxa_parti(comm=hcom12, q=2)    
    hbeta3<- hill_taxa_parti(comm=hcom12, q=3)    
    hbeta4<- hill_taxa_parti(comm=hcom12, q=4)    
    
    ha0_EP1<- c(ha0_EP1, halpha0[1])
    ha0_EP2<- c(ha0_EP2, halpha0[2])
    ha1_EP1<- c(ha1_EP1, halpha1[1])
    ha1_EP2<- c(ha1_EP2, halpha1[2])
    ha2_EP1<- c(ha2_EP1, halpha2[1])
    ha2_EP2<- c(ha2_EP2, halpha2[2])
    ha3_EP1<- c(ha3_EP1, halpha3[1])
    ha3_EP2<- c(ha3_EP2, halpha3[2])
    ha4_EP1<- c(ha4_EP1, halpha4[1])
    ha4_EP2<- c(ha4_EP2, halpha4[2])
    
    dha0<- c(dha0, (abs(halpha0[1]-halpha0[2])))
    dha1<- c(dha1, (abs(halpha1[1]-halpha1[2])))
    dha2<- c(dha2, (abs(halpha2[1]-halpha2[2])))
    dha3<- c(dha3, (abs(halpha3[1]-halpha3[2])))
    dha4<- c(dha4, (abs(halpha4[1]-halpha4[2])))
    
    dha0st<- decostand(dha0, method="range", na.rm=T)
    dha1st<- decostand(dha1, method="range", na.rm=T)
    dha2st<- decostand(dha2, method="range", na.rm=T)
    dha3st<- decostand(dha3, method="range", na.rm=T)
    dha4st<- decostand(dha4, method="range", na.rm=T)
    
    hbsim<- c(hbsim,hbetasim$beta.sim)
    hb0<- c(hb0, hbeta0$TD_beta)
    hb1<- c(hb1, hbeta1$TD_beta)
    hb2<- c(hb2, hbeta2$TD_beta)
    hb3<- c(hb3, hbeta3$TD_beta)
    hb4<- c(hb4, hbeta4$TD_beta)
    
  }

  pcom1<- pcom[rownames(pcom) == pwisedata$EPy1[i],]
  pcom2<- pcom[rownames(pcom) == pwisedata$EPy2[i],]
  
  pcom12<- as.matrix(rbind(pcom1,pcom2))
  
  
  #how to make sure that empty rows in the abundance data sets are skipped
  
  if(nrow(pcom1)< 1 || nrow(pcom2)<1 || rowSums(pcom1)<1 || rowSums(pcom2)<1){
    pa0_EP1<- c(pa0_EP1, NA)
    pa0_EP2<- c(pa0_EP2, NA)
    pa1_EP1<- c(pa1_EP1, NA)
    pa1_EP2<- c(pa1_EP2, NA)
    pa2_EP1<- c(pa2_EP1, NA)
    pa2_EP2<- c(pa2_EP2, NA)
    pa3_EP1<- c(pa3_EP1, NA)
    pa3_EP2<- c(pa3_EP2, NA)
    pa4_EP1<- c(pa4_EP1, NA)
    pa4_EP2<- c(pa4_EP2, NA)
    
    dpa0<- c(dpa0, NA)
    dpa1<- c(dpa1, NA)
    dpa2<- c(dpa2, NA)
    dpa3<- c(dpa3, NA)
    dpa4<- c(dpa4, NA)
    
    dpa0st<- c(dpa0st, NA)
    dpa1st<- c(dpa1st, NA)
    dpa2st<- c(dpa2st, NA)
    dpa3st<- c(dpa3st, NA)
    dpa4st<- c(dpa4st, NA)
    
    pbsim<- c(pbsim,NA)
    pb0<- c(pb0,NA)
    pb1<- c(pb1, NA)
    pb2<- c(pb2, NA)
    pb3<- c(pb3, NA)
    pb4<- c(pb4, NA)
    
  }else {
    
    pcom12<- pcom12[,colSums(pcom12)>0, drop = F]
    
    palpha0<- hill_taxa(comm=pcom12, q=0)
    palpha1<- hill_taxa(comm=pcom12, q=1)
    palpha2<- hill_taxa(comm=pcom12, q=2)
    palpha3<- hill_taxa(comm=pcom12, q=3)
    palpha4<- hill_taxa(comm=pcom12, q=4)
    
    pcom12pa<- decostand(pcom12, method="pa")
    
    pbetasim<- beta.pair(pcom12pa, index.family="sorensen")
    
    pbeta0<- hill_taxa_parti(comm=pcom12, q=0)
    pbeta1<- hill_taxa_parti(comm=pcom12, q=1)    
    pbeta2<- hill_taxa_parti(comm=pcom12, q=2)    
    pbeta3<- hill_taxa_parti(comm=pcom12, q=3)    
    pbeta4<- hill_taxa_parti(comm=pcom12, q=4)    
    
    pa0_EP1<- c(pa0_EP1, palpha0[1])
    pa0_EP2<- c(pa0_EP2, palpha0[2])
    pa1_EP1<- c(pa1_EP1, palpha1[1])
    pa1_EP2<- c(pa1_EP2, palpha1[2])
    pa2_EP1<- c(pa2_EP1, palpha2[1])
    pa2_EP2<- c(pa2_EP2, palpha2[2])
    pa3_EP1<- c(pa3_EP1, palpha3[1])
    pa3_EP2<- c(pa3_EP2, palpha3[2])
    pa4_EP1<- c(pa4_EP1, palpha4[1])
    pa4_EP2<- c(pa4_EP2, palpha4[2])
    
    dpa0<- c(dpa0, (abs(palpha0[1]-palpha0[2])))
    dpa1<- c(dpa1, (abs(palpha1[1]-palpha1[2])))
    dpa2<- c(dpa2, (abs(palpha2[1]-palpha2[2])))
    dpa3<- c(dpa3, (abs(palpha3[1]-palpha3[2])))
    dpa4<- c(dpa4, (abs(palpha4[1]-palpha4[2])))
    
    dpa0st<- decostand(dpa0, method="range", na.rm=T)
    dpa1st<- decostand(dpa1, method="range", na.rm=T)
    dpa2st<- decostand(dpa2, method="range", na.rm=T)
    dpa3st<- decostand(dpa3, method="range", na.rm=T)
    dpa4st<- decostand(dpa4, method="range", na.rm=T)
    
    pbsim<- c(pbsim,pbetasim$beta.sim)
    pb0<- c(pb0, pbeta0$TD_beta)
    pb1<- c(pb1, pbeta1$TD_beta)
    pb2<- c(pb2, pbeta2$TD_beta)
    pb3<- c(pb3, pbeta3$TD_beta)
    pb4<- c(pb4, pbeta4$TD_beta)
    
  }
}

pwisedata_alphabeta<- data.frame(pwisedata,
                                 ha0_EP1, ha0_EP2,
                                 ha1_EP1, ha1_EP2,
                                 ha2_EP1, ha2_EP2,
                                 ha3_EP1, ha3_EP2,
                                 ha4_EP1, ha4_EP2,
                                 dha0, dha1, dha2, dha3, dha4,
                                 dha0st, dha1st, dha2st, dha3st, dha4st,
                                 hbsim, hb0, hb1, hb2, hb3, hb4,
                                 pa0_EP1, pa0_EP2,
                                 pa1_EP1, pa1_EP2,
                                 pa2_EP1, pa2_EP2,
                                 pa3_EP1, pa3_EP2,
                                 pa4_EP1, pa4_EP2,
                                 dpa0, dpa1, dpa2, dpa3, dpa4,
                                 dpa0st, dpa1st, dpa2st, dpa3st, dpa4st,
                                 pbsim, pb0, pb1, pb2, pb3, pb4
)

#reduce dataset to only those combination for which diversity estimates are there
# some of the combinations do not exist, because no data was collected in some
#   of the years. Exclude those cases.
# 42 length (21 * 2, because individually for herbivores and predators)
div_names <- grep("a[0-9]_|ha[0-9]|hb[0-9]|bsim|[ab][0-4]", names(pwisedata_alphabeta), value = T)
# find rows with 21 * 2 NA values (all diversity indices are NA)

# individually assessed predator and herbivore.
#       for predators : removed rows where all diversity indices of predators are NA (ignoring what happens with herbivores)
#       for herbivores : remove rows where all diversity indices of predators are NA (ignoring what happens with predators)
#       --> predators and herbivores may have different sets of comparisons. (Because data situation is different)
# That means the datasets have to be separated here. they could be kept together, but we would
#   have a similar situation as above : The plots present in herb, missing in pred would show as 
#   line of NA in pred diversity indices.

# individually for herbivores
exclude_rows_herb <- apply(pwisedata_alphabeta[, which(colnames(pwisedata_alphabeta) %in%
                                                    grep("^h|^dh", div_names, value = T))],
                      1, function(x) sum(is.na(x)) == length(div_names))
pwise_space_herb <- pwisedata_alphabeta[!exclude_rows_herb, ]
pwise_space_herb <- pwise_space_herb[, which(!colnames(pwise_space_herb) %in% grep("^p|^dp", div_names, value = T))] # exclude predators
# 24231 rows are removed, 98694 kept of a total of a total 122925
#
# individually for predators
exclude_rows_pred <- apply(pwisedata_alphabeta[, which(colnames(pwisedata_alphabeta) %in%
                                                         grep("^p|^dp", div_names, value = T))],
                           1, function(x) sum(is.na(x)) == length(div_names))
pwise_space_pred <- pwisedata_alphabeta[!exclude_rows_pred, ]
pwise_space_pred <- pwise_space_pred[, which(!colnames(pwise_space_pred) %in% grep("^h|^dh", div_names, value = T))] # exclude herbivores
# 30052 rows are removed, 92873 are kept of a total of 122925

#save
save(pwise_space_herb, file="data/OutputData/pwise_space_herb.RData")
save(pwise_space_pred, file="data/OutputData/pwise_space_pred.RData")


rm(pwisedata_alphabeta)
rm(pwisedata)




# 2.2. time ----
pwisedata <- pwise_time
hcom<- herb
pcom<- pred

#important but default: relative abundances are used, so function input is simple
# vectors herbivores
ha0_EP1<- c()
ha1_EP1<- c()
ha2_EP1<- c()
ha3_EP1<- c()
ha4_EP1<- c()

ha0_EP2<- c()
ha1_EP2<- c()
ha2_EP2<- c()
ha3_EP2<- c()
ha4_EP2<- c()

#delta of alpha diversity 
dha0<- c()
dha1<- c()
dha2<- c()
dha3<- c()
dha4<- c()

dha0st<- c()
dha1st<- c()
dha2st<- c()
dha3st<- c()
dha4st<- c()

hbsim<-c()
hb0<- c()
hb1<- c()
hb2<- c()
hb3<- c()
hb4<- c()

# vectors predators
pa0_EP1<- c()
pa1_EP1<- c()
pa2_EP1<- c()
pa3_EP1<- c()
pa4_EP1<- c()

pa0_EP2<- c()
pa1_EP2<- c()
pa2_EP2<- c()
pa3_EP2<- c()
pa4_EP2<- c()

#delta of alpha diversity 
dpa0<- c()
dpa1<- c()
dpa2<- c()
dpa3<- c()
dpa4<- c()

dpa0st<- c()
dpa1st<- c()
dpa2st<- c()
dpa3st<- c()
dpa4st<- c()

pbsim<-c()
pb0<- c()
pb1<- c()
pb2<- c()
pb3<- c()
pb4<- c()

for(i in 1:nrow(pwisedata)){
  
  hcom1<- hcom[rownames(hcom) == pwisedata$EPy1[i],]
  hcom2<- hcom[rownames(hcom) == pwisedata$EPy2[i],]
  
  hcom12<- as.matrix(rbind(hcom1,hcom2))
  hcom12<- hcom12[,colSums(hcom12)>0]
  
  #how to make sure that empty rows in the abundance data sets are skipped
  
  if(nrow(hcom1)<1 || nrow(hcom2)<1 || rowSums(hcom1)<1 || rowSums(hcom2)<1){
    ha0_EP1<- c(ha0_EP1, NA)
    ha0_EP2<- c(ha0_EP2, NA)
    ha1_EP1<- c(ha1_EP1, NA)
    ha1_EP2<- c(ha1_EP2, NA)
    ha2_EP1<- c(ha2_EP1, NA)
    ha2_EP2<- c(ha2_EP2, NA)
    ha3_EP1<- c(ha3_EP1, NA)
    ha3_EP2<- c(ha3_EP2, NA)
    ha4_EP1<- c(ha4_EP1, NA)
    ha4_EP2<- c(ha4_EP2, NA)
    
    dha0<- c(dha0, NA)
    dha1<- c(dha1, NA)
    dha2<- c(dha2, NA)
    dha3<- c(dha3, NA)
    dha4<- c(dha4, NA)
    
    dha0st<- c(dha0, NA)
    dha1st<- c(dha1, NA)
    dha2st<- c(dha2, NA)
    dha3st<- c(dha3, NA)
    dha4st<- c(dha4, NA)
    
    hbsim<- c(hbsim,NA)
    hb0<- c(hb0,NA)
    hb1<- c(hb1, NA)
    hb2<- c(hb2, NA)
    hb3<- c(hb3, NA)
    hb4<- c(hb4, NA)
    
  }else {
    
    halpha0<- hill_taxa(comm=hcom12, q=0)
    halpha1<- hill_taxa(comm=hcom12, q=1)
    halpha2<- hill_taxa(comm=hcom12, q=2)
    halpha3<- hill_taxa(comm=hcom12, q=3)
    halpha4<- hill_taxa(comm=hcom12, q=4)
    
    hcom12pa<- decostand(hcom12, method="pa")
    
    hbetasim<- beta.pair(hcom12pa, index.family="sorensen")
    
    hbeta0<- hill_taxa_parti(comm=hcom12, q=0)
    hbeta1<- hill_taxa_parti(comm=hcom12, q=1)    
    hbeta2<- hill_taxa_parti(comm=hcom12, q=2)    
    hbeta3<- hill_taxa_parti(comm=hcom12, q=3)    
    hbeta4<- hill_taxa_parti(comm=hcom12, q=4)    
    
    ha0_EP1<- c(ha0_EP1, halpha0[1])
    ha0_EP2<- c(ha0_EP2, halpha0[2])
    ha1_EP1<- c(ha1_EP1, halpha1[1])
    ha1_EP2<- c(ha1_EP2, halpha1[2])
    ha2_EP1<- c(ha2_EP1, halpha2[1])
    ha2_EP2<- c(ha2_EP2, halpha2[2])
    ha3_EP1<- c(ha3_EP1, halpha3[1])
    ha3_EP2<- c(ha3_EP2, halpha3[2])
    ha4_EP1<- c(ha4_EP1, halpha4[1])
    ha4_EP2<- c(ha4_EP2, halpha4[2])
    
    dha0<- c(dha0, (abs(halpha0[1]-halpha0[2])))
    dha1<- c(dha1, (abs(halpha1[1]-halpha1[2])))
    dha2<- c(dha2, (abs(halpha2[1]-halpha2[2])))
    dha3<- c(dha3, (abs(halpha3[1]-halpha3[2])))
    dha4<- c(dha4, (abs(halpha4[1]-halpha4[2])))
    
    dha0st<- decostand(dha0, method="range", na.rm=T)
    dha1st<- decostand(dha1, method="range", na.rm=T)
    dha2st<- decostand(dha2, method="range", na.rm=T)
    dha3st<- decostand(dha3, method="range", na.rm=T)
    dha4st<- decostand(dha4, method="range", na.rm=T)
    
    hbsim<- c(hbsim,hbetasim$beta.sim)
    hb0<- c(hb0, hbeta0$TD_beta)
    hb1<- c(hb1, hbeta1$TD_beta)
    hb2<- c(hb2, hbeta2$TD_beta)
    hb3<- c(hb3, hbeta3$TD_beta)
    hb4<- c(hb4, hbeta4$TD_beta)
    
  }
  
  pcom1<- pcom[rownames(pcom) == pwisedata$EPy1[i],]
  pcom2<- pcom[rownames(pcom) == pwisedata$EPy2[i],]
  
  pcom12<- as.matrix(rbind(pcom1,pcom2))
  pcom12<- pcom12[,colSums(pcom12)>0]
  
  #how to make sure that empty rows in the abundance data sets are skipped
  
  if(nrow(hcom1)<1 || nrow(hcom2)<1 || rowSums(hcom1)<1 || rowSums(hcom2)<1){
    pa0_EP1<- c(pa0_EP1, NA)
    pa0_EP2<- c(pa0_EP2, NA)
    pa1_EP1<- c(pa1_EP1, NA)
    pa1_EP2<- c(pa1_EP2, NA)
    pa2_EP1<- c(pa2_EP1, NA)
    pa2_EP2<- c(pa2_EP2, NA)
    pa3_EP1<- c(pa3_EP1, NA)
    pa3_EP2<- c(pa3_EP2, NA)
    pa4_EP1<- c(pa4_EP1, NA)
    pa4_EP2<- c(pa4_EP2, NA)
    
    dpa0<- c(dpa0, NA)
    dpa1<- c(dpa1, NA)
    dpa2<- c(dpa2, NA)
    dpa3<- c(dpa3, NA)
    dpa4<- c(dpa4, NA)
    
    dpa0st<- c(dpa0st, NA)
    dpa1st<- c(dpa1st, NA)
    dpa2st<- c(dpa2st, NA)
    dpa3st<- c(dpa3st, NA)
    dpa4st<- c(dpa4st, NA)
    
    pbsim<- c(pbsim,NA)
    pb0<- c(pb0,NA)
    pb1<- c(pb1, NA)
    pb2<- c(pb2, NA)
    pb3<- c(pb3, NA)
    pb4<- c(pb4, NA)
    
  }else {
    
    palpha0<- hill_taxa(comm=pcom12, q=0)
    palpha1<- hill_taxa(comm=pcom12, q=1)
    palpha2<- hill_taxa(comm=pcom12, q=2)
    palpha3<- hill_taxa(comm=pcom12, q=3)
    palpha4<- hill_taxa(comm=pcom12, q=4)
    
    pcom12pa<- decostand(pcom12, method="pa")
    
    pbetasim<- beta.pair(pcom12pa, index.family="sorensen")
    
    pbeta0<- hill_taxa_parti(comm=pcom12, q=0)
    pbeta1<- hill_taxa_parti(comm=pcom12, q=1)    
    pbeta2<- hill_taxa_parti(comm=pcom12, q=2)    
    pbeta3<- hill_taxa_parti(comm=pcom12, q=3)    
    pbeta4<- hill_taxa_parti(comm=pcom12, q=4)    
    
    pa0_EP1<- c(pa0_EP1, palpha0[1])
    pa0_EP2<- c(pa0_EP2, palpha0[2])
    pa1_EP1<- c(pa1_EP1, palpha1[1])
    pa1_EP2<- c(pa1_EP2, palpha1[2])
    pa2_EP1<- c(pa2_EP1, palpha2[1])
    pa2_EP2<- c(pa2_EP2, palpha2[2])
    pa3_EP1<- c(pa3_EP1, palpha3[1])
    pa3_EP2<- c(pa3_EP2, palpha3[2])
    pa4_EP1<- c(pa4_EP1, palpha4[1])
    pa4_EP2<- c(pa4_EP2, palpha4[2])
    
    dpa0<- c(dpa0, (abs(palpha0[1]-palpha0[2])))
    dpa1<- c(dpa1, (abs(palpha1[1]-palpha1[2])))
    dpa2<- c(dpa2, (abs(palpha2[1]-palpha2[2])))
    dpa3<- c(dpa3, (abs(palpha3[1]-palpha3[2])))
    dpa4<- c(dpa4, (abs(palpha4[1]-palpha4[2])))
    
    dpa0st<- decostand(dpa0, method="range", na.rm=T)
    dpa1st<- decostand(dpa1, method="range", na.rm=T)
    dpa2st<- decostand(dpa2, method="range", na.rm=T)
    dpa3st<- decostand(dpa3, method="range", na.rm=T)
    dpa4st<- decostand(dpa4, method="range", na.rm=T)
    
    pbsim<- c(pbsim,pbetasim$beta.sim)
    pb0<- c(pb0, pbeta0$TD_beta)
    pb1<- c(pb1, pbeta1$TD_beta)
    pb2<- c(pb2, pbeta2$TD_beta)
    pb3<- c(pb3, pbeta3$TD_beta)
    pb4<- c(pb4, pbeta4$TD_beta)
    
  }
}

# some of the diversity vectors do not have the expected number of rows
# fill the missing ones with NA

pwisedata_alphabeta<- data.frame(pwisedata,
                                 ha0_EP1, ha0_EP2,
                                 ha1_EP1, ha1_EP2,
                                 ha2_EP1, ha2_EP2,
                                 ha3_EP1, ha3_EP2,
                                 ha4_EP1, ha4_EP2,
                                 dha0, dha1, dha2, dha3, dha4,
                                 dha0st, dha1st, dha2st, dha3st, dha4st,
                                 hbsim, hb0, hb1, hb2, hb3, hb4,
                                 pa0_EP1, pa0_EP2,
                                 pa1_EP1, pa1_EP2,
                                 pa2_EP1, pa2_EP2,
                                 pa3_EP1, pa3_EP2,
                                 pa4_EP1, pa4_EP2,
                                 dpa0, dpa1, dpa2, dpa3, dpa4,
                                 dpa0st, dpa1st, dpa2st, dpa3st, dpa4st,
                                 pbsim, pb0, pb1, pb2, pb3, pb4
)


#reduce dataset to only those combination for which diversity estimates are there

# -- DONE Noelle
# 42 length (21 * 2, because individually for herbivores and predators)
div_names <- grep("a[0-9]_|ha[0-9]|hb[0-9]|bsim|[ab][0-4]", names(pwisedata_alphabeta), value = T)
# find rows with 21 * 2 NA values (all diversity indices are NA)

# individually assessed predator and herbivore.
#       for predators : removed rows where all diversity indices of predators are NA (ignoring what happens with herbivores)
#       for herbivores : remove rows where all diversity indices of predators are NA (ignoring what happens with predators)
#       --> predators and herbivores may have different sets of comparisons. (Because data situation is different)
# That means the datasets have to be separated here. they could be kept together, but we would
#   have a similar situation as above : The plots present in herb, missing in pred would show as 
#   line of NA in pred diversity indices.

# individually for herbivores
exclude_rows_herb <- apply(pwisedata_alphabeta[, which(colnames(pwisedata_alphabeta) %in%
                                                         grep("^h|^dh", div_names, value = T))],
                           1, function(x) sum(is.na(x)) == length(div_names))
pwise_time_herb <- pwisedata_alphabeta[!exclude_rows_herb, ]
pwise_time_herb <- pwise_time_herb[, which(!colnames(pwise_time_herb) %in% grep("^p|^dp", div_names, value = T))] # exclude predators
# XY rows are removed, XY kept of a total of a total XY
#TODO @Lena fill in true values
#
# individually for predators
exclude_rows_pred <- apply(pwisedata_alphabeta[, which(colnames(pwisedata_alphabeta) %in%
                                                         grep("^p|^dp", div_names, value = T))],
                           1, function(x) sum(is.na(x)) == length(div_names))
pwise_time_pred <- pwisedata_alphabeta[!exclude_rows_pred, ]
pwise_time_pred <- pwise_time_pred[, which(!colnames(pwise_time_pred) %in% grep("^h|^dh", div_names, value = T))] # exclude herbivores
# XY rows are removed, XY are kept of a total of XY
#TODO @Lena fill in true values

#save
save(pwise_time_herb, file="data/OutputData/pwise_time_herb.RData")
save(pwise_time_pred, file="data/OutputData/pwise_time_pred.RData")

#################################
#save all data to analysis and figures folders

#GDM Analysis - set first working directory to folder"2.GDM" then save
save(pwise_time_plants, file="data/InputData/pwise_time_plants.RData")
save(pwise_time_herb, file="data/InputData/pwise_time_herb.RData")
save(pwise_time_pred, file="data/InputData/pwise_time_pred.RData")

save(pwise_space_plants, file="data/InputData/pwise_space_plants.RData")
save(pwise_space_herb, file="data/InputData/pwise_space_herb.RData")
save(pwise_space_pred, file="data/InputData/pwise_space_pred.RData")

#LM Analysis - set first working directory to folder"3.LM" then save
save(pwise_time_plants, file="data/InputData/pwise_time_plants.RData")
save(pwise_time_herb, file="data/InputData/pwise_time_herb.RData")
save(pwise_time_pred, file="data/InputData/pwise_time_pred.RData")

save(pwise_space_plants, file="data/InputData/pwise_space_plants.RData")
save(pwise_space_herb, file="data/InputData/pwise_space_herb.RData")
save(pwise_space_pred, file="data/InputData/pwise_space_pred.RData")


#Figures - set first working directory to folder"4.Plotting Results" then save
save(pwise_time_plants, file="data/InputData/pwise_time_plants.RData")
save(pwise_time_herb, file="data/InputData/pwise_time_herb.RData")
save(pwise_time_pred, file="data/InputData/pwise_time_pred.RData")

save(pwise_space_plants, file="data/InputData/pwise_space_plants.RData")
save(pwise_space_herb, file="data/InputData/pwise_space_herb.RData")
save(pwise_space_pred, file="data/InputData/pwise_space_pred.RData")



rm(pwisedata_alphabeta)
rm(pwisedata)
rm(pwise_time_herb)
rm(pwise_time_pred)
