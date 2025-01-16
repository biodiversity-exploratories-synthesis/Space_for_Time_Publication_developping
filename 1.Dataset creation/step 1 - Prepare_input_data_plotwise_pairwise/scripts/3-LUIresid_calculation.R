# # # # # # # # # # # # # # #  #
#                              #
# Re-calculation of LUIresids  #
#                              #
# # # # # # # # # # # # # # #  #
# script initialised on 09.2.22 by Lena Neuenkamp, Noelle Schenk
# last edit : 04.09.2023
# This script belongs to the manuscript : can space replace time?

# # # # # # # # # # # # # #
# CONTENT                  ----
#
## GOAL
# Re-calculating LUIresids from table pwise_space and pwise_time which has been created
#   in 2-create_temporal_spatial_dataset
# step 1 (Noelle) : produce site-predictor tables from `pwise_space` and `pwise_time datasets`
# step 2 (Lena)   : Calculate LUIresid values from the site-predictor tables
# step 3 (Noelle) : re-convert to site-pair-table suitable for GDM analysis
#
## REQUIREMENTS
# - input data :
#    - pwise space
#    - pwise time
# in a stage BEFORE biotic data was added.
#
## OUTPUT
# 2 new site-pair-tables analogous to `pwise_space` and `pwise_time` with updated
#    LUIresid values.
# Note that running this script on a dataset with already re-calculated LUI residuals
#    would not cause so much harm.

# # # # # # # # # # # # # #
# REQUIREMENTS              ----
#
# load packages
library(data.table)
#
save_interim_output <- T # should the interim datasets be saved? (between step 1 and 3?)
#                            note that the final dataset is saved anyway

# data

# set working directory to folder "1.Dataset creation".
load("1.Dataset creation/step 1 - Prepare_input_data_plotwise_pairwise/data/InputData/pwise_space_template.RData")
load("1.Dataset creation/step 1 - Prepare_input_data_plotwise_pairwise/data/InputData/pwise_time_template.RData") 


# # # # # # # # # # # # # # # # # # # # # #
# STEP 1  : PRODUCE SITE-PREDICTOR TABLE  # ----
# # # # # # # # # # # # # # # # # # # # # #
# convert site-pair table to site-predictor table

# # # # # # # # # # #
# step1.1. : SPATIAL DATASET ----
# included columns are : (all column names which are there after script 2)
# "EP1", "EP2", "EPy1", "EPy2", "YR", "dYR", "HWG1", "HWG2", "RWG1", "RWG2",
# "LUI1", "LUI2", "LUIYear1", "LUIYear2", "Mstd1", "Mstd2", "Gstd1", "Gstd2",
# "Fstd1", "Fstd2", 
# "pH1", "pH2", "soilPCA1", "soilPCA2",
# "G5001", "G5002", "grlperm5001", "grlperm5002",
# "Temp_twom_Sum1", "Temp_twom_Sum1", "precip_raindays1", "precip_raindays2"
# filter columns
colnames_include_site_pair <- c("YR", "EP1", "EP2", "EPy1", "EPy2", "dYR", "HWG1", "HWG2", "RWG1", "RWG2",
                                "LUI1", "LUI2", "LUIyear1", "LUIyear2", "Mstd1", "Mstd2", "Gstd1", "Gstd2",
                                "Fstd1", "Fstd2", 
                                "pH1", "pH2", "soilPCA1", "soilPCA2",
                                "G5001", "G5002", "grlperm5001", "grlperm5002",
                                "Temp_twom_Sum1", "Temp_twom_Sum2", "precip_raindays1", "precip_raindays2")
sitepred_pwise_space <- pwise_space[, ..colnames_include_site_pair] # note : only exclude Exploratory1 and Exploratory2 columns
# stack value 1 and 2 on top of each other (e.g. EP1 and EP2 to same column "EP")
s1cols <- grep("2", names(sitepred_pwise_space), value = T, invert = T) # invert grep : chose all names that do not include "1"
s2cols <- grep("1", names(sitepred_pwise_space), value = T, invert = T)
# check if order of s1cols and s2cols is the same
all(sub("1", "", s1cols) == sub("2", "", s2cols))
sitepred_pwise_space <- unique(rbindlist(list(unique(sitepred_pwise_space[, ..s1cols]),
                                              unique(sitepred_pwise_space[, ..s2cols])), 
                                         use.names = F))
nrow(sitepred_pwise_space) == 150*11 # check : expected number of rows present
length(unique(sitepred_pwise_space$EP1)) == 150 # check : expected number of plots
length(unique(sitepred_pwise_space$LUI1)) <= 150*11 # check : not too many LUI values
nrow(unique(sitepred_pwise_space[, .(HWG1, RWG1)])) # check : 150 unique coordinates there

if(save_interim_output == T){
  save(sitepred_pwise_space, file="1.Dataset creation/step 1 - Prepare_input_data_plotwise_pairwise/data/InputData/sitepred_pwise_space_step1.RData")
  }
rm(sitepred_pwise_space) # remove so it's not mixed up with pwise_sitepred_time
rm(pwise_space)
rm(colnames_include_site_pair); rm(s1cols); rm(s2cols)


# # # # # # # # # # #
# step1.2. : TEMPORAL DATASET ----
all(pwise_time$HWG == pwise_time$HWG1)
all(pwise_time$HWG1 == pwise_time$HWG2)
all(na.omit(pwise_time$RWG) == na.omit(pwise_time$RWG1))
all(na.omit(pwise_time$RWG1) == na.omit(pwise_time$RWG2))
colnames_include_site_pair <- c("EP", "EPy1", "EPy2", "YR1", "YR2", "HWG", "RWG",
                                "LUI1", "LUI2", "LUIyear1", "LUIyear2", "Mstd1", "Mstd2", "Gstd1", "Gstd2",
                                "Fstd1", "Fstd2", 
                                "pH1", "pH2", "soilPCA1", "soilPCA2",
                                "G5001", "G5002", "grlperm5001", "grlperm5002",
                                "Temp_twom_Sum1", "Temp_twom_Sum2", "precip_raindays1", "precip_raindays2")
sitepred_pwise_time <- pwise_time[, ..colnames_include_site_pair] # only exclude Exploratory1 and Exploratory2
# stack value 1 and 2 on top of each other (e.g. EP1 and EP2 to same column "EP")
s1cols <- grep("2", colnames_include_site_pair, value = T, invert = T) # invert grep : chose all names that do not include "1"
s2cols <- grep("1", colnames_include_site_pair, value = T, invert = T)
# check if order of s1cols and s2cols is the same
all(sub("1", "", s1cols) == sub("2", "", s2cols))
sitepred_pwise_time <- unique(rbindlist(list(unique(sitepred_pwise_time[, ..s1cols]),
                                             unique(sitepred_pwise_time[, ..s2cols])), 
                                        use.names = F))
nrow(sitepred_pwise_time) == 150*11 # check : expected number of rows present
length(unique(sitepred_pwise_time$EP)) == 150 # check : expected number of plots
length(unique(sitepred_pwise_time$LUI1)) <= 150*11 # check : not too many LUI values
nrow(unique(sitepred_pwise_time[, .(HWG, RWG)])) # check : 150 unique coordinates there
nrow(unique(sitepred_pwise_time[, .(EP, YR1)])) == 150 * 11

if(save_interim_output == T){
  save(sitepred_pwise_time, file="1.Dataset creation/step 1 - Prepare_input_data_plotwise_pairwise/data/InputData/sitepred_pwise_time_step1.RData")}
rm(sitepred_pwise_time)
rm(pwise_time)
# # # # # # # # # # # # # #

# # # # # # # # # # # # # #
# STEP 2 - CALC RESIDUALS # ----
# # # # # # # # # # # # # #
#load data
if(save_interim_output == T){
  load("1.Dataset creation/step 1 - Prepare_input_data_plotwise_pairwise/data/InputData/sitepred_pwise_space_step1.RData")
  load("1.Dataset creation/step 1 - Prepare_input_data_plotwise_pairwise/data/InputData/sitepred_pwise_time_step1.RData")
}


# # # # # # # # # # #
# step2.1. : SPATIAL DATASET ----
#
# residualling out year, calculate residuals from tables
#
sitepred_pwise_space$LUIres<- residuals(lm(LUI1 ~ YR, data=sitepred_pwise_space))
sitepred_pwise_space$MOWres<- residuals(lm(Mstd1 ~ YR, data=sitepred_pwise_space))
sitepred_pwise_space$GRAres<- residuals(lm(Gstd1 ~ YR, data=sitepred_pwise_space))
sitepred_pwise_space$FERres<- residuals(lm(Fstd1 ~ YR, data=sitepred_pwise_space))

if(save_interim_output == T){
  save(sitepred_pwise_space, file="1.Dataset creation/step 1 - Prepare_input_data_plotwise_pairwise/data/InputData/sitepred_pwise_space_step2.RData")
}



# # # # # # # # # # #
# step1.1. : TEMPORAL DATASET ----
#
# residualling out year
#
sitepred_pwise_time$LUIres<- residuals(lm(LUI1 ~ EP, data=sitepred_pwise_time))
sitepred_pwise_time$MOWres<- residuals(lm(Mstd1 ~ EP, data=sitepred_pwise_time))
sitepred_pwise_time$GRAres<- residuals(lm(Gstd1 ~ EP, data=sitepred_pwise_time))
sitepred_pwise_time$FERres<- residuals(lm(Fstd1 ~ EP, data=sitepred_pwise_time))

if(save_interim_output == T){
  save(sitepred_pwise_time, file="1.Dataset creation/step 1 - Prepare_input_data_plotwise_pairwise/data/InputData/sitepred_pwise_time_step2.RData")
}
rm(sitepred_pwise_time)




# # # # # # # # # # # # # # # # # # # # # #
# STEP 3  : PRODUCE SITE-PAIR TABLE       # ----
# # # # # # # # # # # # # # # # # # # # # #
#
# merge back the newly calculated rows of the site-predictor table to the site-pair table
#
if(save_interim_output == T){
  # load data
  load("1.Dataset creation/step 1 - Prepare_input_data_plotwise_pairwise/data/InputData/sitepred_pwise_space_step2.RData")
  load("1.Dataset creation/step 1 - Prepare_input_data_plotwise_pairwise/data/InputData/sitepred_pwise_time_step2.RData")
  # load site-pair table data
  load("1.Dataset creation/step 1 - Prepare_input_data_plotwise_pairwise/data/InputData/pwise_space_template.RData")
  load("1.Dataset creation/step 1 - Prepare_input_data_plotwise_pairwise/data/InputData/pwise_time_template.RData")
  pwise_space <- data.table(pwise_space)
  old_pwise_space <- data.table(pwise_space)
  pwise_time <- data.table(pwise_time)
  old_pwise_time <- data.table(pwise_time)
}

# remove LUI, MOW, GRA, FER columns from site-predictor table
rem_cols <- c("LUI1", "Mstd1", "Gstd1", "Fstd1")
sitepred_pwise_space[, (rem_cols) := NULL]
sitepred_pwise_time[, (rem_cols) := NULL]
rm(rem_cols)



# # # # # # # # # # #
# step3.1. : SPATIAL DATASET ----
#
# add the newly calculated residuals to the site-pair table
# adjust column names of sitepred_pwise_space (site-predictor) table to match pwise_space (site-pair)
#
nrow(pwise_space) == 122925
s1_match_cols <- c("EP1", "YR", "HWG1", "RWG1")
s2_match_cols <- c("EP2", "YR", "HWG2", "RWG2")

# rename sitepred_pwise_space residual columns for s1 and merge
colnames(sitepred_pwise_space)[grep("res", colnames(sitepred_pwise_space))] <- 
  paste(unlist(strsplit(grep("res", colnames(sitepred_pwise_space), value = T), split = "res")), "1", "res", sep = "")
pwise_space <- merge(pwise_space, 
                     sitepred_pwise_space[, .(EP1, YR, HWG1, RWG1, LUI1res, MOW1res, GRA1res, FER1res)], 
                     by = s1_match_cols, all.x = T)
nrow(pwise_space) == 11175 * 11 # expected number of columns
# rename sitepred_pwise_space residual columns for s2 and merge
colnames(sitepred_pwise_space) <- sub("1", "2", colnames(sitepred_pwise_space))
pwise_space <- merge(pwise_space, 
                     sitepred_pwise_space[, .(EP2, YR, HWG2, RWG2, LUI2res, MOW2res, GRA2res, FER2res)], 
                     by = s2_match_cols, all.x = T)
nrow(pwise_space) == 11175 * 11 

# checking new dataset
test <- unique(rbindlist(list(unique(pwise_space[, ..s1_match_cols]),unique(pwise_space[, ..s2_match_cols])), 
                         use.names = F))
nrow(test) == 150*11 # check : expected number of rows present
length(unique(test$EP1)) == 150 # check : expected number of plots
nrow(unique(test[, .(HWG1, RWG1)])) # check : 150 unique coordinates there
length(unique(pwise_space$LUI1res)) <= 150*11 # check : not too many LUI values
length(unique(pwise_space$GRA1res)) <= 150*11 # check : not too many LUI values
# compare old and new dataset , e.g. pwise_space == old_pwise_space
all(colnames(pwise_space) %in% colnames(old_pwise_space)) # FALSE
# above is false because residual columns were not there yet, and now added.
colnames(pwise_space)[which(!colnames(pwise_space) %in% colnames(old_pwise_space))]
all(colnames(old_pwise_space) %in% colnames(pwise_space))
setcolorder(pwise_space, neworder = colnames(old_pwise_space))
# order rows in same way
pwise_space <- pwise_space[order(EPy2, EPy1),]
old_pwise_space <- old_pwise_space[order(EPy2, EPy1)]
all.equal(data.frame(pwise_space), data.frame(old_pwise_space)) # mismatch
all.equal(data.frame(pwise_space[, 1:34]), data.frame(old_pwise_space)) # only the edited columns differ
all.equal(data.frame(pwise_space[, .(EPy1, EPy2, EP1, EP2, LUI1)]), data.frame(old_pwise_space[, .(EPy1, EPy2, EP1, EP2, LUI1)]))

# save
save(pwise_space, file="1.Dataset creation/step 1 - Prepare_input_data_plotwise_pairwise/data/InputData/sitepred_pwise_space_step3.RData")

rm(pwise_space); rm(old_pwise_space); rm(sitepred_pwise_space)
rm(test); rm(s1_match_cols); rm(s2_match_cols)



# # # # # # # # # # #
# step3.2. : TEMPORAL DATASET ----
#
# add the newly calculated residuals to the site-pair table
# adjust column names of sitepred_pwise_space (site-predictor) table to match pwise_space (site-pair)
#
nrow(pwise_time)
s1_match_cols <- c("EPy1", "EP", "YR1", "HWG", "RWG")
s2_match_cols <- c("EPy2", "EP", "YR2", "HWG", "RWG")

# Rename YR column in sitepred_pwise_time
setnames(sitepred_pwise_time, old = "YR1", new = "YR")
# rename sitepred_pwise_time residual/ year columns for s1 and merge
colnames(sitepred_pwise_time)[grep("res", colnames(sitepred_pwise_time))] <- 
  paste(unlist(strsplit(grep("res", colnames(sitepred_pwise_time), value = T), split = "res")), "1", "res", sep = "")
colnames(sitepred_pwise_time)[grep("YR", colnames(sitepred_pwise_time))]  <- 
  paste(grep("YR", colnames(sitepred_pwise_time), value = T), "1", sep = "")


pwise_time <- merge(pwise_time, 
                    sitepred_pwise_time[, .(EPy1, EP, YR1, HWG, RWG, LUI1res, MOW1res, GRA1res, FER1res)], 
                    by = s1_match_cols, all.x = T)
nrow(pwise_time) == 55 * 150 # expected number of columns
# rename sitepred_pwise_time residual/ year columns for s2 and merge
colnames(sitepred_pwise_time) <- sub("1", "2", colnames(sitepred_pwise_time))
pwise_time <- merge(pwise_time, 
                    sitepred_pwise_time[, .(EPy2, EP, YR2, HWG, RWG, LUI2res, MOW2res, GRA2res, FER2res)], 
                    by = s2_match_cols, all.x = T)
nrow(pwise_time) == 55 * 150

# checking new dataset
test <- unique(rbindlist(list(unique(pwise_time[, ..s1_match_cols]),unique(pwise_time[, ..s2_match_cols])), 
                         use.names = F))
nrow(test) == 150*11 # check : expected number of rows present
length(unique(test$EP)) == 150 # check : expected number of plots
nrow(unique(test[, .(HWG, RWG)])) # check : 150 unique coordinates there
nrow(unique(test[, .(YR1)]))

length(unique(pwise_time$LUI1res)) <= 150*11 # check : not too many LUI values
length(unique(pwise_time$GRA1res)) <= 150*11 # check : not too many LUI values
# compare old and new dataset , e.g. pwise_time == old_pwise_time
all(colnames(pwise_time) %in% colnames(old_pwise_time)) # fALSE because some new columns!
all(colnames(old_pwise_time) %in% colnames(pwise_time))
setcolorder(pwise_time, neworder = colnames(old_pwise_time))
# order rows in same way
pwise_time <- pwise_time[order(EPy2, EPy1),]
old_pwise_time <- old_pwise_time[order(EPy2, EPy1)]
all.equal(data.frame(pwise_time[, 1:32]), data.frame(old_pwise_time)) # only the edited columns differ
all.equal(data.frame(pwise_time[, .(EPy1, EPy2, EP, YR1, YR2,LUI1)]), data.frame(old_pwise_time[, .(EPy1, EPy2, EP, YR1, YR2, LUI1)]))

# save
save(pwise_time, file="1.Dataset creation/step 1 - Prepare_input_data_plotwise_pairwise/data/InputData/sitepred_pwise_time_in_step3.RData")

rm(pwise_time); rm(old_pwise_time); rm(sitepred_pwise_time)
rm(test); rm(s1_match_cols); rm(s2_match_cols)

