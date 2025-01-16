# # # # # # # # # # # # # # #  # # # #
#                                    #
#             SHORTCUT               # ----
#      add new data to step3         #
#                                    #
# # # # # # # # # # # # # # #  # # # #

# Aim : after June 2023, new data was added.
#    In order to not block the later analysis while including it into earlier processes,
#    the new data was added directly to the analysis input dataset, called "*step3*".
#  New data : 
#    - climate : Temp_2m_Sum, precip_raindays
#    - landscape : grlperm500
# Note : in the next step, pairwise datasets will be created from the output of this script

# # # # # #
# 0 - Requirements ----
# # # # # #
#
# PACKAGES
library(data.table)
#
# DATA
load("1.Prepare_input_data_plotwise_pairwise/Outputdata/lui_covariates.RData") # LUI and covariates data
load("OutputDatasets/sitepred_pwise_space_step3.RData") # pwise_space



# # # # # #
# 1 - Quality check ----
# # # # # #
#
# The variables which occur in step3 as well as in lui_covariates can be compared to see if they correspond
# --> quality check
# only compared to pwise_space, could do check with pwise_time (to exclude some wrangling errors) but not done
# at the moment

# get a dataset from the pwise to compare to lui_covariates
replotwised <- unique(pwise_space[, .(EP1, YR, pH_1, soilPCA_1, G500_1, LUI1, FER1, MOW1, GRA1)])
setnames(replotwised, old = c("EP1", "YR"), new = c("EP_PlotID", "Year"))

replotwised <- merge(lui_covariates, replotwised, by = c("EP_PlotID", "Year"))

# LUI
plot(replotwised$LUI, replotwised$LUI1) # corresponds with a few deviations
# show the deviating cases :
plot(replotwised$LUI - replotwised$LUI1)
replotwised[(replotwised$LUI - replotwised$LUI1) > 0, ]
# plots SEG42 and 43, years 2008 and 2018 (see below that this is due to deviations
# in fertilisation measure)
#TODO potentially check if this is an error in the dataset which has been corrected in the
# last few years.

# LUI components
plot(replotwised$Mstd - replotwised$MOW1) # perfectly the same
plot(replotwised$Gstd - replotwised$GRA1) # perfectly the same
plot(replotwised$Fstd, replotwised$FER1) # corresponds with a few deviations
# show the deviating cases :
plot(replotwised$Fstd - replotwised$FER1)
replotwised[(replotwised$Fstd - replotwised$FER1) > 0, ]
# plots SEG42 and 43, years 2008 and 2018


# plot isolation
plot(replotwised$G500 - replotwised$G500_1) # corresponds very well (minimal deviations)
plot(replotwised$G500 - replotwised$G500_1, ylim = c(-0.1, 0.1)) # range is 0 - 100 --> corresponds well


# clean
rm(replotwised)


# # # # # #
# 2 - Adding variables ----
# # # # # #
#
#  New data : 
#    - climate : Temp_2m_Sum, precip_raindays
#    - landscape : grlperm500
new_variable_names <- c("Temp_2m_Sum", "precip_raindays", "grlperm500")

# # # # # #
# 2.1 - pwise_space plants ----
#
load("OutputDatasets/sitepred_pwise_space_step3.RData") # pwise_space
# add to pwise_space
# rename columns in lui_covariates to fit names in pwise datasets
lui_covariates[, EPy1 := EPYear]
lui_covariates[, EPy2 := EPYear]

# merge to EPy1 and to EPy2
# merge to EPy1
pwise_space <- merge(pwise_space, lui_covariates[, .(EPy1, Temp_2m_Sum, precip_raindays, grlperm500)], 
      by = c("EPy1"), all.x = T)
# rename new columns  to plot 1
setnames(pwise_space, old = new_variable_names,
         new = paste0(new_variable_names, "1"))
# merge to EPy2
pwise_space <- merge(pwise_space, lui_covariates[, .(EPy2, Temp_2m_Sum, precip_raindays, grlperm500)], 
                     by = c("EPy2"), all.x = T)
setnames(pwise_space, old = new_variable_names,
         new = paste0(new_variable_names, "2"))

# quality checks
nrow(pwise_space) == 122925
ncol(pwise_space) == 120 + 3 + 3
length(unique(pwise_space$EP1)) == 149
length(unique(pwise_space$EP2)) == 149
length(unique(c(pwise_space$EP1, pwise_space$EP2))) == 150
plot(lui_covariates[EP_PlotID == "AEG01", Temp_2m_Sum],
     unique(pwise_space[EP1 == "AEG01", Temp_2m_Sum1]))
plot(unique(lui_covariates[EP_PlotID == "SEG43", precip_raindays]),
     unique(pwise_space[EP1 == "SEG43", precip_raindays1]))

# save
save(pwise_space, file = "OutputDatasets/sitepred_pwise_space_step3_ch.RData") # ch stands for : climate - history
rm(pwise_space)


# # # # # #
# 2.2 - pwise_space insects ----
#
load("OutputDatasets/sitepred_pwise_space_in_step3.RData") # pwise_space
# 
# merge to EPy1 and to EPy2
pwise_space <- merge(pwise_space, lui_covariates[, .(EPy1, Temp_2m_Sum, precip_raindays, grlperm500)], 
                     by = c("EPy1"), all.x = T)
# rename new columns  to plot 1
setnames(pwise_space, old = new_variable_names,
         new = paste0(new_variable_names, "1"))
# merge to EPy2
pwise_space <- merge(pwise_space, lui_covariates[, .(EPy2, Temp_2m_Sum, precip_raindays, grlperm500)], 
                     by = c("EPy2"), all.x = T)
setnames(pwise_space, old = new_variable_names,
         new = paste0(new_variable_names, "2"))

# quality checks
nrow(pwise_space) == 98946
ncol(pwise_space) == 176 + 3 + 3
length(unique(pwise_space$EP1)) == 149
length(unique(pwise_space$EP2)) == 149
length(unique(c(pwise_space$EP1, pwise_space$EP2))) == 150
length(unique(pwise_space$YR)) == 10
plot(unique(lui_covariates[EP_PlotID == "AEG01" & Year != "2018", Temp_2m_Sum]),
     unique(pwise_space[EP1 == "AEG01", Temp_2m_Sum1]))
plot(unique(lui_covariates[EP_PlotID == "SEG43" & Year != "2018" & Year != "2016", precip_raindays]),
     unique(pwise_space[EP1 == "SEG43", precip_raindays1]))

# There are not all the years in each plot.
# already in step 1 and step2.
# table(pwise_space[, .(EP1, YR)])

# save
save(pwise_space, file = "OutputDatasets/sitepred_pwise_space_in_step3_ch.RData") # ch stands for : climate - history
rm(pwise_space)



# # # # # #
# 2.3 - pwise_time plants ----
#
load("OutputDatasets/sitepred_pwise_time_step3.RData") # pwise_time, 8250 rows

# merge to EPy1 and to EPy2
pwise_time <- merge(pwise_time, lui_covariates[, .(EPy1, Temp_2m_Sum, precip_raindays, grlperm500)], 
                     by = c("EPy1"), all.x = T)
# rename new columns  to plot 1
setnames(pwise_time, old = new_variable_names,
         new = paste0(new_variable_names, "1"))
# merge to EPy2
pwise_time <- merge(pwise_time, lui_covariates[, .(EPy2, Temp_2m_Sum, precip_raindays, grlperm500)], 
                     by = c("EPy2"), all.x = T)
setnames(pwise_time, old = new_variable_names,
         new = paste0(new_variable_names, "2"))

# quality checks
nrow(pwise_time) == 8250
ncol(pwise_time) == 111 + 3 + 3
length(unique(pwise_time$EP)) == 150
length(unique(pwise_time$YR1)) == 10
length(unique(pwise_time$YR2)) == 10
length(unique(c(pwise_time$YR1, pwise_time$YR2))) == 11

# save
save(pwise_time, file = "OutputDatasets/sitepred_pwise_time_step3_ch.RData") # ch stands for : climate - history
rm(pwise_time)



# # # # # #
# 2.4 - pwise_time insects ----
#
#
load("OutputDatasets/sitepred_pwise_time_in_step3.RData") # pwise_time, 8250 rows

# merge to EPy1 and to EPy2
pwise_time <- merge(pwise_time, lui_covariates[, .(EPy1, Temp_2m_Sum, precip_raindays, grlperm500)], 
                    by = c("EPy1"), all.x = T)
# rename new columns  to plot 1
setnames(pwise_time, old = new_variable_names,
         new = paste0(new_variable_names, "1"))
# merge to EPy2
pwise_time <- merge(pwise_time, lui_covariates[, .(EPy2, Temp_2m_Sum, precip_raindays, grlperm500)], 
                    by = c("EPy2"), all.x = T)
setnames(pwise_time, old = new_variable_names,
         new = paste0(new_variable_names, "2"))

# quality checks
nrow(pwise_time) == 5943
ncol(pwise_time) == 181 + 3 + 3
length(unique(pwise_time$EP)) == 150
length(unique(pwise_time$YR1)) == 9
length(unique(pwise_time$YR2)) == 9
length(unique(c(pwise_time$YR1, pwise_time$YR2))) == 10

# save
save(pwise_time, file = "OutputDatasets/sitepred_pwise_time_in_step3_ch.RData") # ch stands for : climate - history
rm(pwise_time)
