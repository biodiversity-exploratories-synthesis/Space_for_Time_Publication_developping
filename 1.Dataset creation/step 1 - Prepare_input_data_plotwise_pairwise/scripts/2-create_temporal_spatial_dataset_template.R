# # # # # # # # # # # # # # #  # # # #
#                                    #
#   Skript 2                         #
#   convert to pairwise distances    # ----
#                                    #
#                                    #
# # # # # # # # # # # # # # #  # # # #

# AIM : convert covariates to pairwise distances
#   based on script : pairwise_differences_plants.R
#   contains code snipped from noelle to convert to pairwise --> use after calculating alpha and beta diversity

# # # # # # # # # # # # # #
# 0 - REQUIREMENTS              ----
library(data.table)
library(vegan)

# # # # #
# 0.a. - DATA ----
# set working directory to folder "Space_for_Time_Publication".
load("./1.Dataset creation/step 1 - Prepare_input_data_plotwise_pairwise/data/InputData/lui_covariates.RData") #called: lui_covariates

# # # # #
# 0.b. - FUNCTIONS ----
#
source("R/BEcreateBEplotnames.R")



# # # # #
# 1 - TEMPORAL AYRS ----
#
# find all year combinations

# # #
# CREATE DIST OBJECT
# create all combinations by using a dist function from vegdist, and filter out the unwanted comparisons
EPnames <- BEcreateBEplotnames(habitat = c("G"))
Years <- seq(2008, 2018)
EPYear <- as.vector(outer(EPnames, Years, paste, sep=".")) # 150 * 11 = 1650 length
temp_test <- data.frame("EPyear" = as.numeric(sub("^.....\\.", "", EPYear)))
rownames(temp_test) <- EPYear

m <- vegdist(temp_test, method = "euclidean") #
# reformatting to pariwise table
m <- as.matrix(m)
m[1:10, 1:10] # visualise
m[!lower.tri(m)] <- 1000 # only select upper triangle --> no double comparisons wanted.
# mark values with value = 1000 and filter out later.
m <- reshape2::melt(m, value.name = "dYR")
temp_test <- m[which(m[,3] != 1000), ] # filter out lower triangle values
temp_test <- data.table(temp_test)

# filter out comparisons across plots: only WITHIN plot is allowed
temp_test[, EP1 := sub("*\\.20[0-9][0-9]", "", temp_test$Var1)]
temp_test[, EP2 := sub("*\\.20[0-9][0-9]", "", temp_test$Var2)]
temp_test <- temp_test[EP1 == EP2, ]
# EP is the same --> only one column needed
temp_test[, EP := EP1]
temp_test[, EP1 := NULL]
temp_test[, EP2 := NULL]

# rename columns
data.table::setnames(temp_test, old = c("Var1", "Var2"), new = c("EPy1", "EPy2"))

# expected number of rows : 
# choose((11), 2) = 55
# 55 * 150 = 8250
nrow(temp_test) == 8250

# can be used to create "pwise_time_ayrs.RData" files
# use as formatvorlage für temporal dataset --> merge everything in there.

pwise_time <- temp_test
save(pwise_time, file = "1.Prepare_input_data_plotwise_pairwise/data/InputData/pwise_time_template.RData")

# # #
# 1.1 FILL IN DATA ----
#
# fill in data to be converted
# for filling in : merge column once with EPy1 and once with EPy2
# work with dataset: lui_covariates
# 
# Naming : LUI - LUI1, LUI2 ; M/G/Fstd - M/G/Fstd1, M/G/Fstd1
#   coordinates as Hochwert (HWG) and Rechtswert (RWG)  - RWG1, RWG2, and HWG1, HWG2

# Step 1 : rename names to <columnname>1 for merging
names(lui_covariates) <- paste0(names(lui_covariates), "1")
setnames(lui_covariates, old = "EP1", new = "EP")
pwise_time <- merge(pwise_time, lui_covariates, by = c("EPy1", "EP"), all.x = T)

# Step 2 : rename all colnames in lui_covariates to <colnames>2
names(lui_covariates) <- sub("1", "2", names(lui_covariates))
pwise_time <- merge(pwise_time, lui_covariates, by = c("EPy2", "EP"), all.x = T)
# rename Year
setnames(pwise_time, old = c("Year1", "Year2"), new = c("YR1", "YR2"))

# plot-specific measurements do not need double mentioning.
pwise_time[, test := HWG1 == HWG2]
unique(pwise_time$test) == TRUE
pwise_time[, test := RWG1 == RWG2]
unique(pwise_time$test) == TRUE
pwise_time[, HWG := HWG1]
pwise_time[, RWG := RWG1]
pwise_time[, c("HWG1", "HWG2", "RWG1", "RWG2", "test") := NULL]

# save data
save(pwise_time, file = "1.Prepare_input_data_plotwise_pairwise/data/InputData/pwise_time_template.RData")

# clean up
rm(m, temp_test, EPnames, EPYear, Years)
names(lui_covariates) <- sub("2", "", names(lui_covariates))





# # # # #
# 2 - SPATIA DATA ----
#
# find all plot combinations

# # #
# CREATE DIST OBJECT
# create all combinations by using a dist function from vegdist, and filter out the unwanted comparisons
EPnames <- BEcreateBEplotnames(habitat = c("G"))
Years <- seq(2008, 2018)
EPYear <- as.vector(outer(EPnames, Years, paste, sep=".")) # 150 * 11 = 1650 length
spatial_test <- data.frame(EPy = EPYear,
                           YR = as.numeric(sub("[A,H,S]EG[0-9][0-9]\\.", "", EPYear)))
rownames(spatial_test) <- spatial_test$EPy
spatial_test$EPy <- NULL
nrow(spatial_test) == 11 * 150

m <- vegdist(spatial_test, method = "euclidean", diag = F, upper = F)
# reformatting to pariwise table
m <- as.matrix(m)
m[1:10, 1:10] # visualise
m[!lower.tri(m, diag = F)] <- 1000 # only select upper triangle --> no double comparisons wanted.
#    mark values with value = 1000 and filter out later.
m <- reshape2::melt(m, value.name = "dYR")

spatial_test <- m[which(m[,3] != 1000), ] # filter out lower triangle values (keep rows with across-years comparisons)
spatial_test <- unique(spatial_test)
spatial_test <- data.table(spatial_test)

# filter out comparisons across plots : only WITHIN year is allowed
spatial_test[, Year1 := as.numeric(sub("[A,H,S]EG[0-9][0-9]\\.", "", spatial_test$Var1))]
spatial_test[, Year2 := as.numeric(sub("[A,H,S]EG[0-9][0-9]\\.", "", spatial_test$Var2))]
spatial_test <- spatial_test[Year1 == Year2, ]
# years are the same, therefore only have one column with Year
spatial_test[, YR := Year1]
spatial_test[, Year1 := NULL]
spatial_test[, Year2 := NULL]

# rename columns
data.table::setnames(spatial_test, old = c("Var1", "Var2"), new = c("EPy1", "EPy2"))

# expected number of rows : 
# choose((150), 2) = 11175
# 11175*11 = 122925
nrow(spatial_test) == 122925

# can be used to create "pwise_space.RData" files
# use as formatvorlage für spatial dataset --> merge everything in there.

pwise_space <- spatial_test
save(pwise_space, file = "1.Prepare_input_data_plotwise_pairwise/InputData/pwise_space_template.RData")

# clean up
rm(m, spatial_test, EPnames, EPYear, Years)


# # #
# 2.1 FILL IN DATA ----

# Step 1 : rename names to <columnname>1 for merging
names(lui_covariates) <- paste0(names(lui_covariates), "1")
setnames(lui_covariates, old = "Year1", new = "YR")
pwise_space <- merge(pwise_space, lui_covariates, by = c("EPy1", "YR"), all.x = T)

# Step 2 : rename all colnames in lui_covariates to <colnames>2
names(lui_covariates) <- sub("1", "2", names(lui_covariates))
pwise_space <- merge(pwise_space, lui_covariates, by = c("EPy2", "YR"), all.x = T)

#Save data
nrow(pwise_space) == 122925
save(pwise_space, file = "1.Dataset creation/step 1 - Prepare_input_data_plotwise_pairwise/InputData/pwise_space_template.RData")

# Next script
# 3-LUIresid_calculation.R

