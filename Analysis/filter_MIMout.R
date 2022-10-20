library(tidyverse)
library(ggplot2)
library(Metrics)

setwd("C:/github/MIMICS_MSBio")


########################
# Bring in data
########################

# Get MIMout data and join it all together
MIMout_files <- list.files(path = "HPC_output/", pattern=".rds")

MIMout_all <- NULL
for(i in length(MIMout_files)){
  loop_df <- NULL
  loop_df <- readRDS(paste0("path", MIMout_files[i]))
  MIMout_all <- rbind(MIMout_all, loop_df)
}


#############################
# Filter out bad data
#############################

# Plausible filters
MCf <- MC_all %>% filter(between(MICpropSOC, 0.01, 0.08)) %>%
  filter(between(MIM_CO_mn, 0.01, 100)) #...etc....


# Remove pset if two rows (temp 15 and 25) don't exist
#...how to do this?


#######################################################
# Apply cost function and grab subset of best data
#######################################################

# Apply cost function
### BALLPARK GUESS EQN!!! 
### --> CHANGE LATER, USE ACTUAL SOIL AND LITTER CN CONTENT FOR THE INCUBATIONS
$rel_cmin <- cmin/30
$cost <- $rel_cmin - MIMICS_prop_cmin

# Combine cost function outcome for each pset
pset_cost <- df %>% group_by(pset_no) %>% summarize(cost_tot = abs sum of cost)

# Filter for top set of results...best 10%...?
#...how good are the results?


#############################
# Parameter space plots
#############################

# Parameter distribution plots
pairs(prm_spc[,1:4],
      col = alpha("blue", 0.3), # Change color
      pch = 16, # Change shape of points
      cex=1.5,
      labels = c("Vslope", "Vint", "Kslope", "Kint"), # Change labels of diagonal
      main = "MIMICS Monte Carlo Top 1% Lowest RMSE",
      upper.panel = NULL,
      xlim=c(0,1),
      ylim=c(0,1),
      breaks=c(0,1))

# Parameter change vs. site properties



