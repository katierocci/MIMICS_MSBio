## Set working drive to repo parent
#setwd("C:/github/MIMICS_MSBio")

# Load packages
library(dplyr) # used for pipes (%>%)
#library(purrr) # used for map() function
library(furrr) # used for map() function

# Load MIMICS model ftn script
source("MIMICS_ftns/MIMICS_INC_daily.R")

# Load parameters for MIMICS
source("MIMICS_set_parameters.R")


#########################
# Example MIMICS runs
#########################

### Single-point run
df <- data.frame(SITE = 'HARV',
                   ANPP = 750,
                   MAT = 25,
                   CLAY = 15,
                   LIG = 20,
                   N = 1,
                   CN = 49,
                   fW = 0.5,
                   soilGWC = 50)

# Run MIMICS and store output as MIMout 
MIMout <- MIMICS_INC_DAILY(df[1,], ndays=200) #hourly step 

# Save output 
write.csv(MIMout, "Cheyenne_HPC/HPC_output/MIMout.csv")
getwd()

###############################################################################
### Multi-point run (hourly)
df_n <- df[rep(seq_len(nrow(df)), each = 250),] # Create example dataframe with n replicate rows
df_n$SITE <-  sprintf("SITE%d",seq(1:250)) # Set unique SITE names

# Run MIMICS using each row of SITE data and row bind the model output into one dataframe
# Set number of cores to use
no_cores <- availableCores()/2 #Use half of availabel cores
print(paste0("Using ", as.character(no_cores), " cores"))
plan(multicore, gc = FALSE, workers = no_cores)

# Run MIMICS!
system.time(
MIMout_n <- df_n %>% split(1:nrow(df_n)) %>% future_map(~ MIMICS_INC_DAILY(df=., ndays=200)) %>% bind_rows()
)

# Release CPU cores
plan(sequential)
nbrOfWorkers()

#Save output
write.csv(MIMout_n, "Cheyenne_HPC/HPC_output/MIMout.csv")

# Clean up memory
gc()

### Map: 76 sec
### Multicore: 87 sec
### Multisession: 9.6 sec, 65.5 sec for 2500
