# # # # # # # # # # # # # # #  # # # #
#                                    #
# calculate input datasets for       # 
#     ALPHA BETA CALCULATIONS        # ----
#                                    #
# # # # # # # # # # # # # # #  # # # #
#
# Aim : prepare specific input data for all 4 files in Alpha_Beta_Calcluations_Chao/


# # # # # # # # # # # # # #
# 1 - PLANTS              ----
#

#load plant data from Bexis: dataset 24247 (plants_long.txt) and 20766 (the script to transform the dataset to a short format)
# script below is the one from 20766

install.packages("reshape2") 
library(reshape2) 

# set working directory to folder "step 2 - Alpha_Beta_Calculations_Chao"

# replace path name to the folder "step 2 - Alpha_Beta_Calculations_Chao" 
mypath<-"C:/Users/Lenovo/Documents/R/space_time_pairwise_HOME" 
	
# open your data
# replace with the correct file name and delimiter 	
plant<-read.table(paste(mypath,"data/RawData/plants_long.txt",sep="/"),h=T, sep=";") 
	
# build the crossed table
# replace the parameters with the appropriate variables regarding your file
# in this example the regular table has the following variables: PlotID, Year, EP_PlotID, Useful_EP_PlotID, Species, Cover
# the crossed table will have: PlotID, Year, EP_PlotID, Useful_EP_PlotID, all other columns will be plant species names ("Species"), values will be based on entries in column "Cover"
plants_short<-dcast(plant,PlotID+Year+EP_PlotID+Useful_EP_PlotID~Species,value.var="Cover") 

# save table
save(plants_short, file="data/InputData/plants_short.RData")

#create column with combined plot year pairs (should be the same in environmental dataplot), ideally usefulEPPlotId.year
#should work with the function (interaction)

#this vector can later be used to select the needed plots from the dataframe during the diversity calculations
rownames(plants_short)<- interaction(plants_short$Useful_EP_PlotID, plants_short$Year)

#remove plot-year combination that have no values/were not sampled/have only NAs

plants_complete<- plants_short[complete.cases(plants_short),c(5:372)]


#replace cover values of 0.1 and 0.5 by 1 in the dataset
plants_complete[plants_complete == 0.1]<- 1
plants_complete[plants_complete == 0.5]<- 1

plants_rd<- plants_complete						      
                                                          
#save data
save(plants_rd, file="data/InputData/plants_rd.RData")


# # # # # # # # # # # # # #
# 1 - INSECTS

# Code preparing the abundance dataset for insects for pairwise diversity calculations

# replace path name to the folder "step 2 - Alpha_Beta_Calculations_Chao" 
mypath<-"C:/Users/Lenovo/Documents/R/space_time_pairwise_HOME" 

# open your data
# replace with the correct file name and delimiter 

# insect dataset (abudances): dataset 21969 (this has 2 data files the abundance data file, inabun_long and the file which plots to explude (insect_exclude)
							      
# opening the txt files does not work, because R does not detect the sep = "" correctly
# alternative: transform to csv in excel and upload csv

insect<-read.csv(paste(mypath,"data/RawData/inabun_long.txt",sep="/"),h=T, sep=";")							
							      
#dataprep - long to short format
#data in long format - transformation into short format

# dataprep
#remove "spec." species

inabun_long<- inabun_long[is.na(inabun_long$spec),]
droplevels(inabun_long$Species)
droplevels(inabun_long$Family)

length(levels(inabun_long$Species))
#1382

length(levels(inabun_long$Family))
#94
							      
# build the crossed table
# replace the parameters with the appropriate variables regarding your file (like plants)
# the crossed table will have: Exploratory, EP, CollectionYear, EPYear, all other columns will be insect species names ("Species"), values will be based on entries in column "NumberAdults"

#transform to short format
inabun_short<- reshape2::dcast(inabun_long, Exploratory + EP + CollectionYear + EPYear ~ Species, 
                     value.var = "NumberAdults", fun=sum)

#all plot*year combinations for which data on certain species are missing
# need to be excluded from the dataset
							      
inex<- read.csv("insect_exclude.csv", header=T, sep=";") 
inabun_short<- inabun_short[!inabun_short$EPYear %in% inex$EPYear,]

save(inabun_short, file="data/InputData/inabun_short.RData")
							      
# 1.1. HERBIVORES              ----
#
# Code preparing the input dataset for pairwise_beta_in_her_Chao.R

#load insect trait data that allow us to have separate abundance tables for herbivores (herbovores+pollinators) and sec. consumers
#data set in Bexis: synthesis dataset
#Assembled species information from grassland EPs (2008-2020) for multidiversity synthesis - November 2020
#Dataset id: 27706, Caterina Penone

#Comment by CP: This is an updated and more complete version of datasets 21726 and 25646.
#This dataset is connected to dataset 27707 (raw abundance/presence-absence) by the "Species" column.


# adjust to file name of the above mentioned traits data set (Bexis id 27706)
insect_traits<-read.csv(paste(mypath,"data/RawData/insect_traits.txt",sep="/"),h=T, sep=";")		

insect_traits <- insect_traits[insect_traits$Species %in% inabun_long$Species,]

#only herbivore abundance (here herbivores and pollinators, which are in our case also herbivores)

herbsp <- insect_traits[insect_traits$Trophic_level=="herbivore",]

# drop unused levels from factors (now all rows have "herbivore")
herbsp$Trophic_level<- droplevels(herbsp$Trophic_level) 
herbsp$Species<- droplevels(herbsp$Species)

# prepare herb dataset in required plot x species format
herb <- inabun_short[,colnames(inabun_short) %in% herbsp$Species]
rownames(herb)<- inabun_short$EPYear
herb<- herb[, colSums(herb)>0] # only include columns with data

save(herb, file="data/InputData//herb.RData")

rm(herb, herbsp, inabun_long, inabun_short, insect_traits)


# # # # # # # # # # # # # #
# 1 - INSECT PREDATORS              ----
#
# Code leading to the prepared pred.RData file - i.e. final rounded/complete compositional data 
# of insect predators only

#only predator abundance
predsp<- insect_traits[insect_traits$Trophic_level=="secondary.consumer",]
predsp$Trophic_level<- droplevels(predsp$Trophic_level)
predsp$Species<- droplevels(predsp$Species)
pred<- inabun_short[,colnames(inabun_short) %in% predsp$Species]

rownames(pred)<- inabun_short$EPYear

pred <- pred[, colSums(pred)>0]# only include columns with data

save(pred, file="data/InputData/pred.RData")
