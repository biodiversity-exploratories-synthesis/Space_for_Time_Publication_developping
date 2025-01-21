###################################################################
#### Adding LUI and component as well as covariate ################
#### means and deltas                          ####################
#### to the pairwise space and time dataset #######################
###################################################################

#########
# 0. Load Data
########

load("1.Dataset creation/step 1 - Prepare_input_data_plotwise_pairwise/data/pwise_space.RData")
load("1.Dataset creation/step 1 - Prepare_input_data_plotwise_pairwise/data/pwise_time.RData")

#######
## 1. LUI and components
#######
# 1.1. LUI and component means
# The means need to be adjusted/standardized across years and EPs 
# for the spatial and temporal datasets, respectively
# Therefore, means are calculated from LUI residuals. 

## 1.1.1.SPATIAL DATA SET
# Calculating means
pwise_space$mLUI<- (pwise_space$LUI1res + pwise_space$LUI2res)/2 
pwise_space$mMOW<- (pwise_space$MOW1res + pwise_space$MOW2res)/2 
pwise_space$mGRA<- (pwise_space$GRA1res + pwise_space$GRA2res)/2 
pwise_space$mFER<- (pwise_space$FER1res + pwise_space$FER2res)/2 


## 1.1.2.TEMPORAL DATA SET
# Calculating means 
pwise_time$mLUI<- (pwise_time$LUI1res + pwise_time$LUI2res)/2 
pwise_time$mMOW<- (pwise_time$MOW1res + pwise_time$MOW2res)/2 
pwise_time$mGRA<- (pwise_time$GRA1res + pwise_time$GRA2res)/2 
pwise_time$mFER<- (pwise_time$FER1res + pwise_time$FER2res)/2 

##########
## 1.2. LUI and component differences/deltas between plot pairs
# Here we use the LUI and component values that are not standardized to YR and EP.
# As the difference in LUI should not be affected by this. 
# Yet, we calculate only absolute deltas as we are only interested in the
# raw difference more than the direction.


## 1.2.1. SPATIAL DATA SET
# Calculating absolute deltas
pwise_space$dLUI<- abs(pwise_space$LUI1-pwise_space$LUI2)
pwise_space$dMOW<- abs(pwise_space$Mstd1-pwise_space$Mstd2)
pwise_space$dGRA<- abs(pwise_space$Gstd1-pwise_space$Gstd2)
pwise_space$dFER<- abs(pwise_space$Fstd1-pwise_space$Fstd2)

## 1.2.2. TEMPORAL DATA SET
# Calculating absolute deltas
pwise_time$dLUI<- abs(pwise_time$LUI1-pwise_time$LUI2)
pwise_time$dMOW<- abs(pwise_time$Mstd1-pwise_time$Mstd2)
pwise_time$dGRA<- abs(pwise_time$Gstd1-pwise_time$Gstd2)
pwise_time$dFER<- abs(pwise_time$Fstd1-pwise_time$Fstd2)


##########
## 2. Covariates

##########
## 2.1. Covariate means
# Here we calculate all means for the spatial dataset and those from which we have temporal replicates (climate) for the temporal dataset

## 2.1.1. SPATIAL DATA SET

### SOIL VARIABLES
# # # # # #
# soil nutrients
pwise_space$mPC1soil<- (pwise_space$PC1soil1 + pwise_space$PC1soil2)/2

# # # # # #
# soil pH
pwise_space$mph<- (pwise_space$pH1 + pwise_space$pH2)/2


### LANDSCAPE MEASURES (radius 500m around each plot)
# # # # # #
# plot isolation
pwise_space$msur<- (pwise_space$G500_1 + pwise_space$G500_2)/2

# # # # # #
# grassland permanency
pwise_space$mgperm<- (pwise_space$grlperm5001 + pwise_space$grlperm5002)/2


### CLIMATE VARIABLES
# # # # # # 
# Temperature sum (Grünlandtemperatursumme), 2m, pear year
pwise_space$mtemp<- (pwise_space$Temp_2m_Sum1 + pwise_space$Temp_2m_Sum2)2

# # # # # # 
# Number of rainy days per year
pwise_space$mrain<- (pwise_space$precip_raindays1 + pwise_space$precip_raindays2)2


## 2.1.2. TEMPORAL DATA SET

### LANDSCAPE MEASURES (radius 500m around each plot)

# # # # # #
# grassland permanency
pwise_space$mgperm<- (pwise_space$grlperm5001 + pwise_space$grlperm5002)/2


### CLIMATE VARIABLES
# # # # # # 
# Temperature sum (Grünlandtemperatursumme), 2m, pear year
pwise_space$mtemp<- (pwise_space$Temp_2m_Sum1 + pwise_space$Temp_2m_Sum2)2

# # # # # # 
# Number of rainy days per year
pwise_space$mrain<- (pwise_space$precip_raindays1 + pwise_space$precip_raindays2)2


##########
## 2.2. Covariate deltas
# Here we calculate all differences for the spatial dataset and those from which we have temporal replicates (climate) for the temporal dataset
# As for land use data, we only calculate absolute differences.

## 2.2.1. SPATIAL DATA SET

### PLOT DISTANCE

# # # # # #
# Euclidean distance between plots - ok as a proxy

a2<- ((pwise_space$HWG2-pwise_space$HWG1))^2
b2<- ((pwise_space$RWG2-pwise_space$RWG1))^2

pwise_space$geo.dist<- sqrt((a2+b2))


### SOIL VARIABLES
# # # # # #
# soil nutrients
pwise_space$dPC1soil<- abs(pwise_space$PC1soil1 - pwise_space$PC1soil2)

# # # # # #
# soil pH
pwise_space$dph<- abs(pwise_space$pH1 - pwise_space$pH2)


### LANDSCAPE MEASURES (radius 500m around each plot)
# # # # # #
# plot isolation
pwise_space$dsur<- abs(pwise_space$G500_1 - pwise_space$G500_2)

# # # # # #
# grassland permanency
pwise_space$dgperm<- abs(pwise_space$grlperm5001 - pwise_space$grlperm5002)


### CLIMATE VARIABLES
# # # # # # 
# Temperature sum (Grünlandtemperatursumme), 2m, pear year
pwise_space$dtemp<- abs(pwise_space$Temp_2m_Sum1 - pwise_space$Temp_2m_Sum2)

# # # # # # 
# Number of rainy days per year
pwise_space$drain<- abs(pwise_space$precip_raindays1 - pwise_space$precip_raindays2)


## 2.1.2. TEMPORAL DATA SET

### LANDSCAPE MEASURES (radius 500m around each plot)

# # # # # #
# grassland permanency
pwise_space$dgperm<- abs(pwise_space$grlperm5001 - pwise_space$grlperm5002)


### CLIMATE VARIABLES
# # # # # # 
# Temperature sum (Grünlandtemperatursumme), 2m, pear year
pwise_space$dtemp<- (pwise_space$Temp_2m_Sum1 - pwise_space$Temp_2m_Sum2)

# # # # # # 
# Number of rainy days per year
pwise_space$drain<- (pwise_space$precip_raindays1 - pwise_space$precip_raindays2)


##########
## 3. Saving data

save(pwise_space, file="1.Dataset creation/step 1 - Prepare_input_data_plotwise_pairwise/output/pwise_space.RData")
save(pwise_time, file="1.Dataset creation/step 1 - Prepare_input_data_plotwise_pairwise/output/pwise_time.RData")

save(pwise_space, file="1.Dataset creation/step 2 - Alpha_Beta_Calculations_Chao/RawData/pwise_space.RData")
save(pwise_time, file="1.Dataset creation/step 2 - Alpha_Beta_Calculations_Chao/RawData/pwise_time.RData")
