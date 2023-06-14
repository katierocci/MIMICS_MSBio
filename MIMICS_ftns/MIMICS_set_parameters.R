########################################
# Set MIMICS default parameters
########################################
#why so many 6 repeats? Maybe since the max number of dif values below is 6?
Vslope  <- rep(0.063, 6) #regression coefficient for Vmax equation
Vint    <- rep(5.47, 6) #regression intercept for Vmax equation
aV      <- rep(0.000000075, 6) #tuning coefficient for Vmax equation
Kslope  <- rep(0.02, 6) #regression coefficient for Km equation
Kint    <- rep(3.19, 6) #regression intercept for Km equation
aK      <- rep(0.15625, 6) #tuning coefficient for Km equation
vMOD    <- c(2, 0.4, 2, 0.6, 0.6, 0.4) #modifiers for Vmax for LITm, LITs, SOMa entering MICr and then MICk
kMOD    <- c(8, 2, 4, 2, 4, 6)#modifiers for Km for LITm, LITs, SOMa entering MICr and then MICk
KO      <- c(6, 6) #modifies Km for oxidation of SOMc
CUE     <- c(0.5, 0.25, 0.7, 0.35) #initial values: c(0.5, 0.25, 0.7, 0.35)
#Microbial growth efficiency (r-MET, r-STRUC, K-MET, K-STRUC) - not varying with climate or anything
tau_r   <- c(0.00052, 0.3) #initially 0.00052, only changing first one #microbial biomass turnover rate for r-selected, second number is some sort of modifier?
tau_K   <- c(0.00024, 0.1) #initially 0.00024, only changing first one #microbial biomass turnover rate for K-selected, second number is some sort of modifier?
Tau_MOD <- c(100, 0.6, 1.3, 3.5) #modifier for microbial turnover rate - implemented in a more complicated way in MIMICS function
Tau_MULT <- 1 #modifer for microbial turnover rate (MULT usually indicates Derek modifier)
fPHYS_r <- c(0, 0) #fraction of turned over r-selected microbial biomass partitioned to SOMp - zero here bc of incubation setup
fPHYS_K <- c(0, 0) #fraction of turned over K-selected microbial biomass partitioned to SOMp - zero here bc of incubation setup
fCHEM_r <- c(0.1, -3, 1) #fraction of turned over r-selected microbial biomass partitioned to SOMc
fCHEM_K <- c(0.3, -3, 1) #fraction of turned over K-selected microbial biomass partitioned to SOMc
fSOM_p  <- c(0.000015, -1.5) #fraction of SOMp that can be desorbed to available pool?
PHYS_scalar <- c(2, -2, NA, NA, NA, NA) #Scalar for texture effects on SOMp
FI      <- c(0.05, 0.05) #fraction of litter inputs directly entering SOMp and SOMc
fmet_p <- c(1, 0.85, 0.013) #partitioning of litter into the metabolic pool - inputs for Fm equation, very similar to Parton DayCent one
depth <- 5 # Set soil depth

#Set required multipliers to 1 (e.g. use default)
Tau_MULT = 1
desorb_MULT = 1
fPHYS_MULT = 1
VMAX_MULT = 1
KM_MULT = 1

# Set defualts for MIMrepeat use

Vslope_default <- rep(0.063, 6)
Vint_default <- rep(5.47, 6)
Kslope_default <- rep(0.02, 6)
Kint_default <- rep(3.19, 6)
CUE_default <- c(0.5, 0.25, 0.7, 0.35)
Tau_MULT_default <- 1
vMOD_default    <- c(2, 0.4, 2, 0.6, 0.6, 0.4)
kMOD_default    <- c(8, 2, 4, 2, 4, 6)

########################################
# Apply parameter multipliers
########################################
# Vslope = Vslope * 1
# Vint = Vint * 1
# Kslope = Kslope * 1
# Kint = Kint * 1
# CUE = CUE * 1
# Tau_MULT = 1
# desorb_MULT = 1
# fPHYS_MULT = 1
# VMAX_MULT = 1
# KM_MULT = 1
