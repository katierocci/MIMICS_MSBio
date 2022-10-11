## Set working drive to repo parent
setwd("C:/github/MIMICS_MSBio")

# Load packages
library(dplyr) # used for pipes (%>%)
library(purrr) # used for map() function
library(furrr) # used for map() function

# Load MIMICS model ftn script
source("MIMICS_ftns/MIMICS_INC_hourly.R")

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
MIMout <- MIMICS_INC_HOUR(df[1,], ndays=200) #hourly step 

############
# Plots
###########
library(ggpubr)

# Litter mass
plot_LIT <- ggplot(MIMout, aes(y=LITs, x=DAY, color="Structural")) + geom_line(size=1) +
  geom_line(aes(y=LITm, x=DAY, color="Metabolic"), size=1) +
  theme_bw() +
  ylab("Litter mass remaining (%)") +
  xlab("Incubation Time (days)") +
  labs(color = "Litter Pool")

# SOM & MIC pools
plot_SOM_MIC <- ggplot(MIMout, aes(SOMc, x=DAY, color="SOMc")) + geom_line(size=1) +
  geom_line(aes(y=SOMa, x=DAY, color="SOMa"), size=1) +
  geom_line(aes(y=MICr, x=DAY, color="MIC-r"), size=1) +
  geom_line(aes(y=MICK, x=DAY, color="MIC-K"), size=1) +
  theme_bw() +
  ylab("Microbial and soil C") +
  xlab("Incubation Time (days)") +
  labs(color = "C Pool") +
  ylim(0, 3)

# CO2 fraction
plot_CO2 <- ggplot(MIMout, aes(y=rowSums(MIMout[,10:11])/rowSums(MIMout[,3:11]),
                               x=DAY, color="CO2-C")) + geom_line(size=1) +
  theme_bw() +
  ylab("CO2 (fraction of initial") +
  xlab("Incubation Time (days)") +
  labs(color = "C Pool")


# Build a panel plot
ggarrange(plot_LIT, plot_SOM_MIC, plot_CO2,
          nrow=3,
          ncol=1)

###############################################################################
### Multi-point run (hourly)
df50 <- df[rep(seq_len(nrow(df)), each = 250),] # Create example dataframe with 10 replicate rows
df50$SITE <-  sprintf("SITE%d",seq(1:250)) # Set unique SITE names

# Run MIMICS using each row of SITE data and row bind the model output into one dataframe
# Set number of cores to use
no_cores <- availableCores() - 1
plan(multicore, gc = TRUE, workers = no_cores)

# Run MIMICS!
system.time(
MIMout_50 <- df50 %>% split(1:nrow(df50)) %>% future_map(~ MIMICS_INC_HOUR(df=., ndays=200)) %>% bind_rows()
)

# Release CPU cores
plan(sequential)
nbrOfWorkers()

# Clean up memory
gc()



### Map: 76 sec
### Multicore: 87 sec
### Multisession: 9.6 sec, 65.5 sec for 2500
