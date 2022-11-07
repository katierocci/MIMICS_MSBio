### MIMICS MC

setwd("C:/github/MIMICS_MSBio")

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
data <- read.csv("Data/MIMICS_forcings/MSBio_forcing_INC_mock.csv", as.is=T)
colnames(data)[1] <- "ID"

####################################
# Use the brute force MIMICS ftn
####################################

# Set desired number of random parameter runs
MIM_runs <- 1000

### Create random parameter dataframe
## Parameter range informed by range observed over 10+ MCMC analysis results
rand_params <- data.frame(#Vslope_x = runif(MIM_runs, 0.5, 2),  
                          Vint_x = runif(MIM_runs, 0.8, 1.3)#,  
                          #Kslope_x = runif(MIM_runs, 0.5, 2),  
                          #Kint_x = runif(MIM_runs, 0.5, 2)#,  
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
plan(multisession, gc = TRUE, workers = no_cores)

# Run MIMICS!

print(paste0("Starting ", MIM_runs, " runs"))
print(paste0("Start time: ", Sys.time()))

start_time <- Sys.time()
MC_MIMICS <- rand_params %>% split(1:nrow(rand_params)) %>% future_map(~MIMrepeat(forcing_df = data, rparams = .), .progress=TRUE) %>% bind_rows()

#check_df <- MC_MIMICS %>% filter(DAY == 200) %>% filter(SITE == "BART") %>% filter(run_num == 1)

wall_time <- Sys.time() - start_time
print(paste0("Wall time: ", as.character(wall_time)))


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
saveRDS(MC_MIMICS, paste0("temp/MSBio_MC_", as.character(MIM_runs), "_", format(Sys.time(), "%Y%m%d_%H%M%S_"),  ".rds"))


######################
# Diagnostic Plots
######################
# library(ggplot2)
# library(ggpubr)
# 
# MC_200 <- MC_MIMICS %>% filter(DAY == 200)
# 
# MIMplot <- MC_MIMICS %>% 
#         filter(SITE == "A") %>%
#         filter(run_num == 1)
# 
# # Litter mass
# plot_LIT <- ggplot(MIMplot, aes(y=LITs, x=DAY, color="Structural")) + geom_line(size=1) +
#   geom_line(aes(y=LITm, x=DAY, color="Metabolic"), size=1) +
#   theme_bw() +
#   ylab("Litter mass remaining (%)") +
#   xlab("Incubation Time (days)") +
#   labs(color = "Litter Pool") +
#   facet_wrap(~MAT)
# 
# # SOM & MIC pools
# plot_SOM_MIC <- ggplot(MIMplot, aes(SOMc, x=DAY, color="SOMc")) + geom_line(size=1) +
#   geom_line(aes(y=SOMa, x=DAY, color="SOMa"), size=1) +
#   geom_line(aes(y=MICr, x=DAY, color="MIC-r"), size=1) +
#   geom_line(aes(y=MICK, x=DAY, color="MIC-K"), size=1) +
#   theme_bw() +
#   ylab("Microbial and soil C") +
#   xlab("Incubation Time (days)") +
#   labs(color = "C Pool") +
#   #ylim(0, 50) +
#   facet_wrap(~MAT)
# 
# # CO2 fraction
# plot_CO2_prop <- ggplot(MIMplot, aes(y=rowSums(MIMplot[,11:12])/rowSums(MIMplot[,4:12]),
#                    x=DAY, color="CO2-C")) + geom_line(size=1) + #, color="blue") +
#   theme_bw() +
#   ylab("CO2 \n(fraction of total)") +
#   xlab("Incubation Time (days)") +
#   labs(color = "C Pool") +
#   #ylim(0,1) +
#   facet_wrap(~MAT)
# 
# #CO2 total
# plot_CO2_tot <- ggplot(MIMplot, aes(y=rowSums(MIMplot[,11:12]),
#                                     x=DAY, color="CO2-C")) + geom_line(size=1) + #, color="light blue") +
#   theme_bw() +
#   ylab("Cumulative \nRespiration (unit C)") +
#   xlab("Incubation Time (days)") +
#   labs(color = "C Pool") +
#   facet_wrap(~MAT)
# 
# ptbl <- ggtexttable(round(MIMplot[1, 21:24], 3), rows = NULL,
#                     theme = ttheme("light"))
# 
# # Build a panel plot
# ggarrange(plot_LIT, plot_SOM_MIC, plot_CO2_prop, ptbl,
#           nrow=4,
#           ncol=1)



