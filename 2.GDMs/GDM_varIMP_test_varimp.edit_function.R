########
# Testing Varimp function February 2022
#
# this script aims to test the functionality of the varImp function. It was used during development of the varimp function.
# Script written by Noelle Schenk, 07.02.2022
# was previously called `22-02_test_varimp.edit_function.R`

# necessary functions
source("GDM/gdm.varIMP_edit.R")
source("GDM/matrix_perm_permutateSitePair.R")

# read in test data
spTable <- readRDS(file = "GDM/testdataset_time.Rds")
spTable <- readRDS(file = "GDM/testdataset_space.Rds")
spTable <- readRDS(file = "GDM/testdataset_space_allregions.rds")

#####
# check input dataset
# number of unique coordinates
nrow(unique(rbindlist(list(unique(spTable[, c("s1.xCoord", "s1.yCoord")]), unique(spTable[, c("s2.xCoord", "s2.yCoord")])), use.names = F)))
choose((150+2-1), 2) * 11 >= nrow(spTable)

# run gdmvarimp

# temporal
# start.time <- Sys.time()
out <- gdm.varImp.synthesis.edit(spTable = spTable, geo = T, fullModelOnly = T, nPerm = 20,
                          includeYear = T, spacefortimedataset = "temporal")
# end.time <- Sys.time()


# spatial

out <- gdm.varImp.synthesis.edit(spTable = spTable, geo = T, fullModelOnly = T, nPerm = 3,
                                 includeYear = T, spacefortimedataset = "spatial")
# geo = T does not always lead to an error
    # error found in cases
        # geo = T, no error
        # geo = T, no error

# ##### assessing runtime
# # runtime <- data.table(nperm = rep(1, 10), time = rep(1, 10), mod = rep("mod", 10))
# # saveRDS(runtime, file = "runtime_assessment.rds")
# runtime <- readRDS("runtime_assessment.rds")
# runtime <- rbindlist(list(runtime, data.table(nperm = 50,
#                                               time = as.numeric(end.time - start.time)*60/4,
#                                               mod = "temporal")))
# ggplot(runtime, aes(x = nperm, y = time, col = mod, pch = mod, lty = mod)) +
#   geom_point() +
#   geom_line() +
#   xlim(c(0, 55)) +
#   ylim(c(0, 360))
# ##### stop assessing runtime



# #####
# # TEST RESULTS with dummy variable
# # Add dummy column to see if results are the same
# # using time dataset
# spTable2 <- data.table(copy(spTable))
# spTable2$s1.dum <- runif(nrow(spTable), min=0, max=1)
# spTable2$s2.dum <- runif(nrow(spTable), min=0, max=1)
# colorder <- c("distance", "weights", "s1.xCoord", "s1.yCoord", "s2.xCoord", "s2.yCoord",
#               "s1.Year", "s1.LU", "s1.dum",
#               "s2.Year", "s2.LU", "s2.dum")
# spTable2 <- spTable2[, ..colorder]
# # remove duplicate rows for time dataset
# spTable2 <- unique(spTable2, by=c("s1.yCoord","s2.xCoord"))
# spTable2 <- formatsitepair(spTable2, bioFormat = 4, predData = 7:10)
# 
# spTable2 <- formatsitepair(spTable2, bioFormat = 4, predData = 7:10)
# set.seed(16)
# test_old_dummy <- gdm.varImp(spTable2, geo = T, nPerm = 3, fullModelOnly = T)
# set.seed(16)
# test_new_dummy <- gdm.varImp.synthesis.edit(spTable2, geo = T, nPerm = 3, fullModelOnly = T, includeYear = T, spacefortimedataset = "temporal")
# # comparison is not 1:1 possible, but model deviance, % dev expl and pvalues match ok
# #####
