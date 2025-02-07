# These function create the dataset and graphical output of
# scatterplots comparing mean and SD of spatial and temporal diversity measures 
# per trophic group and region.
# Functions correspond to step 3. in the Figures script
# and Fig. 4 in the associated publication.


# function 1 ####

# The following part contains 2 functions which calculate mean and sd 
# for each measure and return a table for plotting.

calc_temporal_regionwise_mean_sd <- function(div_colname, div_dataset, outname){
  # Temporal
  # mean
  temp <- aggregate(formula = formula(paste(div_colname, " ~ reg + group", sep = "")),
                    data = div_dataset,
                    FUN = function(x) mean(x, na.rm = T))
  names(temp)[grep(div_colname, names(temp))] <- paste("temp_mean", outname, sep = "_") # name new column
  # standard error
  div_dataset[, se := sd(get(div_colname), na.rm = T), by = .(reg, group)] # for "special" SE calc, go with data.table
  div_dataset[, sd := se] # save standard deviation for later
  div_dataset[, se := se / sqrt(N_unique_obs)]
  temp <- merge(temp, unique(div_dataset[, .(reg, group, se)]), by = c("reg", "group"))
  div_dataset[, se := NULL]
  names(temp)[grep("se", names(temp))] <- paste("temp_sd", outname, sep = "_")
  # standard deviation
  temp <- merge(temp, unique(div_dataset[, .(reg, group, sd)]), by = c("reg", "group")) # add standard deviation to table
  div_dataset[, sd := NULL]
  names(temp)[grep("^sd", names(temp))] <- paste("temp_standarddeviation", outname, sep = "_")
  return(temp)
}

calc_spatial_regionwise_mean_sd_and_add_to_temporal <- function(temp, sdiv_colname, sdiv_dataset, soutname){
  # Spatial
  # mean
  temp <- merge(temp, aggregate(formula = formula(paste(sdiv_colname, " ~ reg + group", sep = "")), 
                                data = sdiv_dataset, 
                                FUN = function(x) mean(x, na.rm = T)), by = c("reg", "group"))
  names(temp)[grep(sdiv_colname, names(temp))] <- paste("space_mean", soutname, sep = "_")
  # standard error
  sdiv_dataset[, se := sd(get(sdiv_colname), na.rm = T), by = .(reg, group)]
  sdiv_dataset[, sd := se]
  sdiv_dataset[, se := se / sqrt(N_unique_obs)]
  temp <- merge(temp, unique(sdiv_dataset[, .(reg, group, se)]), by = c("reg", "group"))
  sdiv_dataset[, se := NULL]
  names(temp)[grep("se", names(temp))] <- paste("space_sd", soutname, sep = "_")
  # standard deviation
  temp <- merge(temp, unique(sdiv_dataset[, .(reg, group, sd)]), by = c("reg", "group")) # add standard deviation to table
  sdiv_dataset[, sd := NULL]
  names(temp)[grep("^sd", names(temp))] <- paste("space_standarddeviation", soutname, sep = "_")
  return(temp)
}

# function 2 ####

# This function reproduces the same plot over and over
generate_time_space_comparison_scatterplot <- function(dat, dat_name, pretty_name, errorbartype = "standarderror"){
  # function arguments : 
  # dat <- bsim                 # dataset which contains the values to plot
  # dat_name <- "bsim"          # suffix in the dataset column names
  # pretty_name <- "beta sim"   # plot title
  # errorbartype = "standarddeviation"
  names(dat) <- sub(paste("_", dat_name, sep = ""), "", names(dat))  # make column names general
  p <- ggplot(dat, aes(x = temp_mean, y = space_mean, col = col_by_group)) +
    scale_colour_identity() +
    geom_point() + theme_half_open() +
    xlim(c(0, 1)) + ylim(c(0,1 )) +
    geom_abline(intercept = 0, slope =1, linetype = "dashed") +
    xlab("temporal beta-diversity") +
    ylab("spatial beta-diversity") +
    ggtitle(pretty_name)
  if(errorbartype == "standarderror"){
    p <- p + geom_errorbar(data = dat, mapping = aes(x = temp_mean, 
                                                     ymin = space_mean - space_sd,
                                                     ymax = space_mean + space_sd), width = errorbar_width) +
      geom_errorbar(data = dat, mapping = aes(y = space_mean, 
                                              xmin = temp_mean - temp_sd,
                                              xmax = temp_mean + temp_sd), width = errorbar_width)
  }
  if(errorbartype == "standarddeviation"){
    print("plotting standard deviation")
    p <- p + geom_errorbar(data = dat, mapping = aes(x = temp_mean, 
                                                     ymin = space_mean - space_standarddeviation,
                                                     ymax = space_mean + space_standarddeviation), width = errorbar_width) +
      geom_errorbar(data = dat, mapping = aes(y = space_mean, 
                                              xmin = temp_mean - temp_standarddeviation,
                                              xmax = temp_mean + temp_standarddeviation), width = errorbar_width) +
      xlim(c(0, 1.4)) + ylim(c(0, 1.4))
  }
  return(p)
}

# example
#without std
# p_bsim <- generate_time_space_comparison_scatterplot(dat = bsim, dat_name = "bsim", 
                                                     pretty_name = "Turnover, bsim")

#with std
# p_bsim_std <- generate_time_space_comparison_scatterplot(dat = bsim, dat_name = "bsim", 
#                                                         pretty_name = "Turnover (beta sim)",
#                                                         errorbartype = "standarddeviation")
