
# Bring in MIMICS RXEQ function
source("MIMICS_ftns/RXEQ_day_ftn.R")

# Data table input to the MIMICS_24 ftn requires the following data 
# ..with matching column names

# "SITE"  Identifier for the sampling location
# "MAT"   Mean annual temperature (or temp of incubation) in deg C
# "CLAY"  Soil clay content in percent
# "LIG"   Lignin content of the litter 
# "N"     N content of the litter
# "CN"    C:N ratio of the litter 

###########################################
# MIMICS single point function
###########################################
MIMICS_INC_DAILY <- function(df, ndays=200, parm_mult){
  
  #DEBUG
  #df <- data[1,]
  #ndays <- 200
  
  ### note var may be used to collect run notes
  note <- ""
  
  ###Bring in forcing data
  ANPP       <- df$ANPP/2 # Estimate SOM inputs to be 50% of ANPP
  fCLAY      <- df$CLAY/100  # Convert clay from percent to fraction
  TSOI       <- df$MAT # Use MAT as proxy for soil temperature
  #lig_N      <- df$lig_N 
  fW         <- df$fW # Soil moisture scalar
  soilGWC    <- df$soilGWC/100 # For use as scalar on decomposition rate
  
  ### Set fMET (=partioning coefficient for metabolic vs. structural litter pools)
  
  ## Option A: Defualt fMET equation using lig:N values
  #fMET <- fmet_p[1] * (fmet_p[2] - fmet_p[3] * lig_N) 
  
  ## Option B: LTER "SHORTCUT" fMET value (average from LiDET)
  #fMET <- 0.3846423
  
  ## Option C: 
  lig    <- df$LIG/100
  Nnew   <- 1/df$CN/2.5                  	                       # N in litter additions
  fMET  <- fmet_p[1] * (fmet_p[2] - fmet_p[3] * lig / Nnew) 
  
  
  ############################################################
  # MIMICS MODEL PARAMETERIZATION STARTS HERE
  ############################################################
  
  # Calc litter input rate
  EST_LIT <- (ANPP / (365*24)) #* 1e3 / 1e4
  #print(EST_LIT)# gC/m2/h (from gC/m2/y) then mgC/cm2/h(from gC/m2/h) 
  
  # ------------ caclulate parameters ---------------
  Vmax     <- exp(TSOI * Vslope + Vint) * aV 
  Km       <- exp(TSOI * Kslope + Kint) * aK
  
  #ANPP strongly correlated with MAP
  Tau_MOD1 <- sqrt(ANPP/Tau_MOD[1])         
  Tau_MOD2 <- Tau_MOD[4]                        
  Tau_MOD1[Tau_MOD1 < Tau_MOD[2]] <- Tau_MOD[2]
  Tau_MOD1[Tau_MOD1 > Tau_MOD[3]] <- Tau_MOD[3] 
  tau <- c(tau_r[1]*exp(tau_r[2]*fMET), 
           tau_K[1]*exp(tau_K[2]*fMET))   
  tau <- tau * Tau_MOD1 * Tau_MOD2 * Tau_MULT 
  
  fPHYS    <- c(fPHYS_r[1] * exp(fPHYS_r[2]*fCLAY), 
                fPHYS_K[1] * exp(fPHYS_K[2]*fCLAY)) 	            
  fCHEM    <- c(fCHEM_r[1] * exp(fCHEM_r[2]*fMET) * fCHEM_r[3], 
                fCHEM_K[1] * exp(fCHEM_K[2]*fMET) * fCHEM_K[3]) 	
  fAVAI    <- 1 - (fPHYS + fCHEM)
  
  desorb   <- fSOM_p[1] * exp(fSOM_p[2]*(fCLAY))                  
  desorb <- desorb * desorb_MULT
  
  fPHYS <- fPHYS * fPHYS_MULT
  
  pSCALAR  <- PHYS_scalar[1] * exp(PHYS_scalar[2]*(sqrt(fCLAY)))  #Scalar for texture effects on SOMp
  v_MOD    <- vMOD  
  k_MOD    <- kMOD 
  k_MOD[3] <- k_MOD[3] * pSCALAR    
  k_MOD[6] <- k_MOD[6] * pSCALAR    
  
  VMAX     <- Vmax * v_MOD 
  KM       <- Km / k_MOD

  #MC scalars on VMAX and Km
  VMAX <- VMAX * VMAX_MULT
  KM <- KM * KM_MULT
  
  ############################################################
  # SITE CONDITIONS SCALARS
  ############################################################
  
  # Soil moisture scalar
  #VMAX <- VMAX #* (0.5 + (0.5 * fW))  # Revamp with a more complex function
    
  # Site moisture condition scalars on VMAX and Km
  #VMAX <- VMAX * (0.7 + (0.3 * soilGWC)) 
  #KM <- KM * (0.7 + (0.3 * soilGWC)) 
  
  ############################################################
  # MIMICS MODEL RUN STARTS HERE
  ############################################################
 
  # Set run length
  nday   <- ndays # number of days
   
  # Open matrices to store model output
  LITmin  <- rep(NA, dim=4)
  MICtrn  <- rep(NA, dim=6)
  SOMmin  <- rep(NA, dim=2)
  DEsorb  <- rep(NA, dim=1)
  OXIDAT  <- rep(NA, dim=1)
  
  # Create dataframe to store model output
  MIMout <- data.frame(SITE = rep(df$SITE, nday),
                       DAY = rep(NA, nday),
                       LITm = rep(NA, nday),
                       LITs = rep(NA, nday),
                       MICr = rep(NA, nday),
                       MICK = rep(NA, nday),
                       SOMp = rep(NA, nday),
                       SOMc = rep(NA, nday),
                       SOMa = rep(NA, nday),
                       CO2_MICr = rep(NA, nday),
                       CO2_MICK = rep(NA, nday))
  
  # Initialize model pools and fluxes
  I        <- rep(0,2)
  LIT_1    <- 100   
  LIT_2    <- 100
  MIC_1    <- 0.01
  MIC_2    <- 0.01
  SOM_1    <- 0
  SOM_2    <- 0
  SOM_3    <- 0
  CO2_1    <- 0
  CO2_2    <- 0
  
  ## Set global parameters to pass to RXEQ function
#DEBUG: Double check this works correctly when run outside of the loop
  .GlobalEnv$VMAX <- VMAX
  .GlobalEnv$KM <- KM
  .GlobalEnv$fPHYS <- fPHYS
  .GlobalEnv$fCHEM <- fCHEM
  .GlobalEnv$fAVAI <- fAVAI
  .GlobalEnv$I <- I
  .GlobalEnv$tau <- tau
  .GlobalEnv$LITmin <- LITmin
  .GlobalEnv$SOMmin <- SOMmin
  .GlobalEnv$MICtrn <- MICtrn
  .GlobalEnv$desorb <- desorb
  .GlobalEnv$DEsorb <- DEsorb
  .GlobalEnv$OXIDAT <- OXIDAT
  
  # Create vector of parameter values
  tpars <- c(I = I, VMAX = VMAX, KM = KM, CUE = CUE, 
             fPHYS = fPHYS, fCHEM = fCHEM, fAVAI = fAVAI, FI = FI, 
             tau   = tau, LITmin = LITmin, SOMmin = SOMmin, MICtrn = MICtrn, 
             desorb= desorb, DEsorb = DEsorb, OXIDAT = OXIDAT, KO = KO)
  
  ### BEGIN MODEL LOOP ###
  # Loop over 24 hours for specified number of days 
  for (d in 1:nday)  {
      
    # Get model output from RXEQ ftn
    update <- RXEQ_day(y = c(LIT_1 = LIT_1, LIT_2 = LIT_2, 
                             MIC_1 = MIC_1, MIC_2 = MIC_2, 
                             SOM_1 = SOM_1, SOM_2 = SOM_2, 
                             SOM_3 = SOM_3),
                       pars = tpars)
    
    # Update C pools
    LIT_1  <- LIT_1 + update[[1]][1]
    LIT_2  <- LIT_2 + update[[1]][2]
    MIC_1  <- MIC_1 + update[[1]][3]
    MIC_2  <- MIC_2 + update[[1]][4]
    SOM_1  <- SOM_1 + update[[1]][5]
    SOM_2  <- SOM_2 + update[[1]][6]
    SOM_3  <- SOM_3 + update[[1]][7]
    CO2_1  <- CO2_1 + update[[1]][8]
    CO2_2  <- CO2_2 + update[[1]][9]
    #remove(UPpars, UPy, update)
    
    # Store daily output
    MIMout$DAY[d] <- d
    MIMout$LITm[d] <- LIT_1
    MIMout$LITs[d] <- LIT_2
    MIMout$MICr[d] <- MIC_1
    MIMout$MICK[d] <- MIC_2
    MIMout$SOMp[d] <- SOM_1
    MIMout$SOMc[d] <- SOM_2
    MIMout$SOMa[d] <- SOM_3
    MIMout$CO2_MICr[d] <- CO2_1
    MIMout$CO2_MICK[d] <- CO2_2
  }	#close daily loop	

  # Return daily model output as datatable
  return(MIMout) #Tbl with each of the 200 days
  #return(MIMout %>% filter(DAY == 200)) #Only day 200
  
} #close ftn

#####################
# Example use of 
#####################

# ##############################################
# #single point run
# ##############################################
# df <- data.frame(SITE = 'TEST',
#                    ANPP = 750,
#                    MAT = 15,
#                    CLAY = 15,
#                    LIG = 21,
#                    N = 1,
#                    CN = 49,
#                    fw=0.6)
# 
# source("MIMICS_ftns/MIMICS_set_parameters.R")
# Vslope = Vslope_default * 1.4425676
# Vint = Vint_default * 0.9413847
# Kslope = Kslope_default * 1.4416737 
# Kint = Kint_default * 0.4436813 
# 
# MIMout <- MIMICS_INC_DAILY(df[1,])


# ##############################################
# # Full forcing dataset run
# ##############################################
# data <- data <- read.csv("RCrk_Modelling_Data/LTER_SITE_1.csv", as.is=T)

# MIMrun <- data %>% split(1:nrow(data)) %>% map(~ MIMICS1(df=.)) %>% bind_rows()
# MIMrun <- data %>% cbind(MIMrun %>% select(-Site, -TSOI))


############
# Plots
###########
# library(ggplot2)
# library(ggpubr)

# Litter mass
# plot_LIT <- ggplot(MIMout, aes(y=LITs, x=DAY, color="Structural")) + geom_line(size=1) +
#   geom_line(aes(y=LITm, x=DAY, color="Metabolic"), size=1) +
#   theme_bw() +
#   ylab("Litter mass remaining (%)") +
#   xlab("Incubation Time (days)") +
#   labs(color = "Litter Pool")
# 
# # SOM & MIC pools
# plot_SOM_MIC <- ggplot(MIMout, aes(SOMc, x=DAY, color="SOMc")) + geom_line(size=1) +
#   geom_line(aes(y=SOMa, x=DAY, color="SOMa"), size=1) +
#   geom_line(aes(y=MICr, x=DAY, color="MIC-r"), size=1) +
#   geom_line(aes(y=MICK, x=DAY, color="MIC-K"), size=1) +
#   theme_bw() +
#   ylab("Microbial and soil C") +
#   xlab("Incubation Time (days)") +
#   labs(color = "C Pool") +
#   ylim(0, 3)
# 
# # CO2 fraction
# plot_CO2 <- ggplot(MIMout, aes(y=rowSums(MIMout[,10:11])/rowSums(MIMout[,3:11]),
#                    x=DAY, color="CO2-C")) + geom_line(size=1) +
#   theme_bw() +
#   ylab("CO2 (fraction of initial") +
#   xlab("Incubation Time (days)") +
#   labs(color = "C Pool")
# 
# # Build a panel plot
# ggarrange(plot_LIT, plot_SOM_MIC, plot_CO2,
#           nrow=3,
#           ncol=1)

