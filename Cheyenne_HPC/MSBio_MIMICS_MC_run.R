### MIMICS MC

#setwd("C:/github/MIMICS_MSBio")

########################################
# Load R packages
########################################
library(dplyr)
library(rootSolve)
library(purrr)
library(furrr)


########################################
# Load MIMICS data and ftns 
########################################
source("MIMICS_ftns/MIMICS_set_parameters.R") #Sets the initial parameters
source("MIMICS_ftns/MIMICS_INC_daily.R")
source("MIMICS_ftns/MIMICS_repeat_base.R")


########################################
# Load forcing data
########################################
data <- read.csv("Data/MIMICS_forcings/MSBio_forcing_temp_trts_only.csv", as.is=T)

####################################
# Use the brute force MIMICS ftn
####################################

# Set desired number of random parameter runs
MIM_runs <- 10000

### Create random parameter dataframe
## Parameter range informed by range observed over 10+ MCMC analysis results
rand_params <- data.frame(Vslope_x = runif(MIM_runs, 0.8, 1.5),  
                          Vint_x = runif(MIM_runs, 0.8, 1.5),  
                          Kslope_x = runif(MIM_runs, 0.8, 1.5),  
                          Kint_x = runif(MIM_runs, 0.8, 1.5)#,  
                          # Tau_x = runif(MIM_runs, 0.3, 3),  
                          # CUE_x = runif(MIM_runs, 0.5, 1.5),  
                          # desorb_x = runif(MIM_runs, 0.001, 0.3),  
                          # fPHYS_x = runif(MIM_runs, 0.4, 4)  
)

rand_params$run_num <- seq(1,MIM_runs,1)

#DEBUG: run default paramaters
#rand_params[1,1:4] = 1

# Set number of cores to use
no_cores <- availableCores() - 1
plan(multicore, gc = TRUE, workers = no_cores)

# Run MIMICS!

print(paste0("Starting ", MIM_runs, " runs"))
print(paste0("Start time: ", Sys.time()))

start_time <- Sys.time()
MC_MIMICS <- rand_params %>% split(1:nrow(rand_params)) %>% future_map(~MIMrepeat(forcing_df = data, rparams = .), .progress=TRUE) %>% bind_rows()

#check_df <- MC_MIMICS %>% filter(DAY == 200) %>% filter(SITE == "BART") %>% filter(run_num == 1)

wall_time <- Sys.time() - start_time
print(paste0("Task time: ", as.character(wall_time)))


# Release CPU cores
plan(sequential)
nbrOfWorkers()

# Clean up memory
gc()


## Join forcing data and parameters to MIMICS output table
MC_MIMICS <- MC_MIMICS %>% left_join(data %>% select(-SITE), by="ID")
MC_MIMICS <- MC_MIMICS %>% left_join(rand_params)


##########################################
# Save MC output data
##########################################
saveRDS(MC_MIMICS, paste0("Cheyenne_HPC/HPC_output/MSBio_MIM_MC_runs-", as.character(MIM_runs), "_", format(Sys.time(), "%Y%m%d_%H%M%S_"),  ".rds"))




