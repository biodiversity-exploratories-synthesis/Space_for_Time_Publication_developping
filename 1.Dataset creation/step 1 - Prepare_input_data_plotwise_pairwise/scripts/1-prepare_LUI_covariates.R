# # # # # # # # # # # # # # #  # # # #
#                                    #
#        prepare datasets            # ----
#              LUI                   #
#           Covariates               #
#                                    #
# # # # # # # # # # # # # # #  # # # #

# Aim : prepare LUI and Covariates data from BExIS raw data.
#       values are continuously added to the growing variable lui_covariates.
# Output : 1.Dataset creation/step 1 - Prepare_input_data_plotwise_pairwise/Outputdata/lui_covariates.RData
#       which is a large table of all lui and all covariate values across plots and years.
#       note that not all variables were measured in every year. therefore, this output
#       table contains quite some duplicated values (e.g. pH is the same across years)
# Note : in the next step, pairwise datasets will be created from the output of this script

# # # # # #
# 0 - Requirements ----
# # # # # #
# PACKAGES
library(data.table)
library(naniar)
library(ade4) # for pca
library(factoextra) # visualise pca
library(corrplot)
#

# set working directory to folder "Space_for_Time_Publication".

# FUNCTIONS
source("RFunctions/BEplotZeros.R")
# impute by mean (soil)
f1 <- function(vec) {
  m <- mean(vec, na.rm = TRUE)
  vec[is.na(vec)] <- m
  return(vec)
}
#
# DATA


# data upload is with exemplary file names, please insert the respective
# file names for the data sets downloaded from BeXis.

# Links and IDs can be found in the document "Supp_Dataset_description_GitHub.docx".

# LUI data as downloaded from BExIS
lui <- fread("1.Dataset creation/step 1 - Prepare_input_data_plotwise_pairwise/data/Rawdata/LUI_default components set__2023-03-13.txt")

# Coordinates
plot.info <- fread("1.Dataset creation/step 1 - Prepare_input_data_plotwise_pairwise/data/Rawdata/20907_7_Dataset/20907_7_data.csv")
plot.names <- fread("1.Dataset creation/step 1 - Prepare_input_data_plotwise_pairwise/data/Rawdata/20826_7_Dataset/20826_7_data.csv")
# Use Longitude_Dec_ and Latitude_Dec (as in Gossner 2016 analysis)

# load climate data as downloaded from BExIS
clim <- fread("1.Dataset creation/step 1 - Prepare_input_data_plotwise_pairwise/data/Rawdata/Ta_200_Ta_200_DTR_Ta_200_gruenlandtemperatursumme_2008_2018_ef4298c30197645c/plots.csv")
prec <- fread("1.Dataset creation/step 1 - Prepare_input_data_plotwise_pairwise/data/Rawdata/result_7726019124401210128/plots.csv")
# soil
soil1 <- read.table("1.Dataset creation/step 1 - Prepare_input_data_plotwise_pairwise/data/Rawdata/14446_soil_grassland.txt", h = T)
soil2 <- read.table("1.Dataset creation/step 1 - Prepare_input_data_plotwise_pairwise/data/Rawdata/17086_pools_Grasslands.txt", h = T)
pH <- read.table("1.Dataset creation/step 1 - Prepare_input_data_plotwise_pairwise/data/Rawdata/14447_pH_grassland.txt", h=T)
# landscape measures
landscape <- fread("1.Dataset creation/step 1 - Prepare_input_data_plotwise_pairwise/data/Rawdata/30318_NEE2020_Aggregated environmental and land use covariates of grassland EPs.csv")
landscape <- data.table(BEplotZeros(landscape, column = "EP_PlotID", plotnam = "EP"))
setnames(landscape, old = c("EP_PlotID", "EP"), new = c("oldEP", "EP_PlotID"))
plt.sur <- fread("1.Dataset creation/step 1 - Prepare_input_data_plotwise_pairwise/data/Rawdata/18148_landscape.txt")
#
# # # # # #
# DOCUMENTATION
# # # # # #


# # # # # # # # # # # #
#
# 1 - Prepare LUI data ----
#
# # # # # # # # # # # #

# # # # # #
# 1.1 - wrangle data ----
# # # # # #
lui$Year <- as.numeric(sub(" 12:00:00 AM$", "", sub("^[0-9]\\/[0-9]\\/", "", lui$Year))) # convert datetime to year
# add useful plotIDs
lui <- data.table(BEplotZeros(lui, column = "EP_PlotID", plotnam = "EP0"))
lui[, EP_PlotID := NULL]
setnames(lui, old = "EP0", new = "EP_PlotID")
# rename Year to LUIyear
setnames(lui, old = "Year", new = "LUIyear")
# drop unused columns
lui[, c("isVIP", "isMIP") := NULL]
# use better plotIDs : those with 0!

# # # # # #
# select years of interest
# # # # # #
# We use the LUI of the previous year to affect the next year.
# Because LUI of 2008 includes LUI action until autumn. However, the arthropods
# and plants are rather measured in spring, and for measurement reasons 
# before the mowing. Therefore, if we want to see what drives biotic communities
# in 2008, we need measures of LUI from 2007.
lui <- lui[LUIyear %in% seq(2007, 2017, 1), ]
lui[, Year := LUIyear + 1 ]
# create EPy column
lui$EPy <- paste(lui$EP_PlotID, lui$Year, sep = ".")
setnames(lui, old = "EP_PlotID", new = "EP")

# # # # # #
# 1.2 - calculate LUI manually ----
# # # # # #
# global standardisation (over all regions)
# global standardisation (over all years)
# --> make LUI comparable across years and regions

# Source : BExIS LUI tool manual : 
# https://www.bexis.uni-jena.de/LUI/LUICalculation/DownloadPDF?fileName=HowToPlusDiagram.pdf
# global standardization (FG, MG, GG) across regions takes care that LUI
# values are directly comparable, but the regions may then cover different sections of the
# LUI gradient. After careful examination of the properties, we now recommend global
# standardization for most analyses, which makes the interpretation of the land use effects
# easier.
# long-term standardization of Gmean,
# Mmean, and Fmean for a LUI for one specific year removes the equal weight of Fi, Mi, and Gi
# within that year, but allows for direct comparisons of LUI values across different years. If
# you compare effects across different years, we recommend that you use the mean in that
# time span for standardization.

# calculate the global means
lui$MeanF <- rep(mean(lui$TotalFertilization, na.rm = T), times = nrow(lui))
lui$MeanG <- rep(mean(lui$TotalGrazing,na.rm = T), times = nrow(lui))
lui$MeanM <- rep(mean(lui$TotalMowing,na.rm = T), times = nrow(lui))
# calculate the LUI and its components for each year and each plot, using the global mean
lui$LUI<- sqrt((lui$TotalFertilization/lui$MeanF) + (lui$TotalMowing/lui$MeanM) + (lui$TotalGrazing / lui$MeanG))
lui$Fstd<- lui$TotalFertilization/lui$MeanF
lui$Mstd<- lui$TotalMowing/lui$MeanM
lui$Gstd<- lui$TotalGrazing/lui$MeanG

# # # # # #
# clean    
# # # # # #
lui_covariates <- lui[, .(EPy, EP, Year, LUIyear, Exploratory, LUI, Fstd, Mstd, Gstd)]
rm(lui)

# # # # # # # # # # # #
#
# 2 - Prepare Covariates ----
#
# # # # # # # # # # # #
#

# # # # # #
# 2.1 - soil ----
# # # # # #
# approach from Gossner et al. 2016
# Below code is exactly the same as in Gossner 2016 --> comparable
# Note : code occurs in script "plot_heterogeneity.R" AND in script "alphadiv_LAST.R"
# Source was not entirely clear --> using variant from alphadiv_LAST.R 
#    excluding "Total_C" (because = Inorganic_C + Organic_C)

# # # # # #
# soil nutrients
#
# select only variables included in Gossner et al and only grasslands
soil1 <- soil1[soil1$Type=="G", c("EP_Plotid", "Inorganic_C", "Organic_C", "Total_N", "CN_ratio")]
soil2 <- soil2[soil2$Type=="G", c("EP_Plotid", "OC_stock", "N_stock")]

# merge datasets and replace NA's with mean values (no influence in PCA)
soil <- merge(soil1, soil2, by = "EP_Plotid")
rm(soil1, soil2)
apply(soil, 2, function(x) sum(is.na(x)))
soil[is.na(soil$OC_stock), ]
# 4 missing values each in OC_stock and N_stock, at the same place.
# imputing by the mean

soil <- data.frame(EP_Plotid = soil$EP_Plotid, 
                   apply(soil[, -which(names(soil) == "EP_Plotid")], 2, f1))
apply(soil, 2, function(x) sum(is.na(x))) # missing values are imputed

# perform PCA
PCAsoil <- ade4::dudi.pca(soil[, -which(names(soil) == "EP_Plotid")], 
                          center = TRUE, scale = TRUE, scann = F)
screeplot(PCAsoil)
factoextra::fviz_pca_biplot(PCAsoil)

# extract principal component 1 and add plot name
PC1soil <- data.frame(EP_PlotID = soil$EP_Plotid, 
                      soilPCA = PCAsoil$li[,1])

# reverse soil nutrients so fertile is high values
PC1soil$soilPCA <- -1 * PC1soil$soilPCA
rm(soil, PCAsoil, f1)

# merge PC1soil to lui_covariates
PC1soil <- BEplotZeros(PC1soil, column = "EP_PlotID", plotnam = "EP")
PC1soil <- PC1soil[, c("EP", "soilPCA")]

# # # # # #
# soil pH
pH <- pH[pH$Type=="G", c("EP_Plotid","pH_1","pH_2")]

pH$pH <- rowMeans(pH[, c("pH_1","pH_2")], na.rm = T)
pH <- pH[,c("EP_Plotid","pH")]
pH <- BEplotZeros(pH, column = "EP_Plotid", plotnam = "EP")
pH <- pH[, c("EP", "pH")]

# merge soil and pH
env <- merge(pH, PC1soil, by="EP", all = T)
nrow(env) == 150
length(unique(env$EP)) == 150

# add to lui_covariates dataset
lui_covariates <- merge(lui_covariates, env, by = "EP", all = T)

rm(pH, PC1soil, env)


# # # # # #
# 2.2 - Landscape measures ----
# # # # # #
# Several landscape measures could be used. Goal : Stick to Gossner 2014 analysis,
# therefore use plot isolation with 500m radius.
# New : add grassland permanency

# # # # # #
# 2.2.1. plot coordinates ----
# 
plot.info <- plot.info[,c("Plot_ID", "Longitude_Dec_Plotcenter", "Latitude_Dec_Plotcenter"), with=F]
# get plot names
setnames(plot.info,"Plot_ID", "PlotID")
setkey(plot.info, PlotID)
setkey(plot.names, PlotID)
plot.info <- merge(plot.info, plot.names, by="PlotID")
plot.info <- plot.info[Landuse == "G" & ActivePlot == "yes", .(EP_PlotID, Longitude_Dec_Plotcenter, Latitude_Dec_Plotcenter)]
plot.info <- BEplotZeros(plot.info, column = "EP_PlotID", plotnam = "EP")
plot.info <- plot.info[, c("EP", "Longitude_Dec_Plotcenter", "Latitude_Dec_Plotcenter")]
nrow(plot.info) == 150
length(unique(plot.info$EP)) == 150
# Rename to HWG and RWG : Hochwert (HWG) and Rechtswert (RWG)
#  Rechtswert (RWG) corresponds to Latitude
#  Hochwert (HWG) corresponds to Longitude
setnames(plot.info, old = c("Longitude_Dec_Plotcenter", "Latitude_Dec_Plotcenter"),
         new = c("HWG", "RWG"))

# merge to growing dataset
lui_covariates <- merge(lui_covariates, plot.info, by = "EP", all = T)

# # # # # #
# 2.2.1. grl perm + plot isol. ----
# 
# grassland permanency and plot isolation
# plot isolation from dataset 18148
# grassland permanency from dataset 30318

# note the current plot isolation values are 0.99 correlated with Gossner 2016 data.
#   newer dataset used here.

# the below outcommented code is an assessment. not needed for variable cleaning, but
# was used to assess from which dataset values were used.
#
# # # # # # #
# # ASSESSMENT which radius?
# #
# # radii 500 and 1km are correlated
# plot(100 - landscape$Grassland.500 * 100, 100 - landscape$Grassland.1000 * 100)
# plot(landscape$PLAND_G_500, 100 - landscape$Grassland.1000 * 100)
# plot(landscape$Arable.500, landscape$Arable.1000) # correlated % crop
# # Arable 500m and 1km are slightly correlated with landscape land-use
# # note landscape-land-use is at 1km
# plot(landscape$Arable.1000, 100 - landscape$LII * 100)
# plot(landscape$Arable.500, 100 - landscape$LII * 100)
# m <- cor(landscape[, .(Grassland.500, Grassland.1000, Arable.500, Arable.1000, LII,
#                        Grassland.perm.500, Grassland.perm.1000, G500)],
#          method = "spearman")
# corrplot(m, diag = F, method = "number", type = "lower")
# rm(m)
# #
# # Conclusion : we take Grassland 500 radius, because close to Gossner analysis.

# merge both datasets
# clean plot isolation data
plt.sur <- unique(plt.sur) # for some reason each row has two values
plt.sur <- plt.sur[,c("EP_PLOTID","G500")]
plt.sur <- BEplotZeros(plt.sur, column = "EP_PLOTID", plotnam = "EP") 
plt.sur$PLAND_G_500 <- 100-plt.sur$G500*100    # change so more isolated plots have higher values
# clean landscape data
landscape <- BEplotZeros(landscape, column = "EP_PlotID", plotnam = "EP")
landscape <- data.table(merge(landscape, plt.sur, by = "EP", all = T))

landscape <- landscape[, .(EP, PLAND_G_500, Grassland.perm.500)]
setnames(landscape, old = c("PLAND_G_500", "Grassland.perm.500"), 
         new = c("G500", "grlperm500"))

# add to lui_covariates
lui_covariates <- merge(lui_covariates, landscape, by = "EP", all = T)
rm(landscape, plt.sur)




# # # # # #
# 2.3 - climate ----
# # # # # #
# Reference : https://www.sciencedirect.com/science/article/pii/S009830042030618X?via%3Dihub
#        Wöllauer et al. 2021 "TubeDB: An on-demand processing database system for climate station data"
#
# Temperature : Ta_200_gruenlandtemperatursumme
#    reason : Gruenlandtemperatursumme (sum of temperature in grassland is a derived
#       measure, recommended by climate reasearch team). We preferred the 
#       sum over the "DTR" measure, because the sum is closer to the actual
#       measurement (in terms of correlation.)
#     In theory, we would prefer to use temperature 10 cm aboveground, but the years 2008
#       and 2009 are quite gappy. Temp at 2 m is not gappy. And Temp at 2 m is correlated
#       with Temperature at 10 cm.
# Precipitation : precipitation_radolan_rain_days : only derived measure, recommended
#    by climate group.

# basic clean
clim$plotID <- sub("f", "", clim$plotID)
allclim <- merge(clim, prec, by = c("plotID", "datetime"))
# select indices for dataset
allclim <- allclim[, .(plotID, datetime, Ta_200_gruenlandtemperatursumme, 
                       precipitation_radolan_rain_days)]
# give nicer names
setnames(allclim, old = c("plotID", "datetime", "Ta_200_gruenlandtemperatursumme", "precipitation_radolan_rain_days"),
         new = c("EP_PlotID", "Year", "Temp_twom_Sum", "precip_raindays"))
# delete unused datasets
rm(clim, prec)
#
# check
nrow(allclim) == nrow(allclim[Year %in% seq(2008, 2018), ])

# add to lui_covariates dataset
setnames(allclim, old = "EP_PlotID", new = "EP")
lui_covariates <- merge(lui_covariates, allclim, by = c("EP", "Year"), all = T)
rm(allclim)

# # # # # # # # # # # #
#
# 3 - Save script output ----
#
# # # # # # # # # # # #

save(lui_covariates, file="1.Dataset creation/step 1 - Prepare_input_data_plotwise_pairwise/data/InputData/lui_covariates.RData")


# Next script : 
# 2-create_temporal_spatial_dataset_template.R 
