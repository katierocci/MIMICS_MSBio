## Set working drive
setwd("C:/github/MIMICS_MSBio/LIDET")

#Libraries
library(rootSolve)
library(dplyr)
library(ggplot2)

#bring in RXEQ function
source("MIMICS_ftns.R")

########################################
# Set MIMICS default parameters
########################################
Vslope  <- rep(0.063, 6)
Vint    <- rep(5.47, 6)
aV      <- rep(0.000008, 6)  
#Kslope  <- rep(c(0.025, 0.035, 0.025),2)
Kslope  <- rep(c(0.017, 0.027, 0.017),2) #Different from MIMICS-HiRes
Kint    <- rep(3.19, 6)
aK      <- rep(10, 6)
vMOD    <- c(10, 2, 10, 3, 3, 2)
kMOD    <- c(8, 2, 4, 2, 4, 6)
#KO      <- c(6, 6)
KO      <- c(4, 4) #Different from MIMICS-HiRes
CUE     <- c(0.55, 0.25, 0.75, 0.35) #for LITm and LITs entering MICr and MICK, respectively
tau_r   <- c(0.00052, 0.3)
tau_K   <- c(0.00024, 0.1)
Tau_MOD <- c(100, 0.8, 1.2, 2)
Tau_MULT <- 1
fPHYS_r <- c(0.3, 1.3)
fPHYS_K <- c(0.2, 0.8)
fCHEM_r <- c(0.1, -3, 1)
fCHEM_K <- c(0.3, -3, 1)
fSOM_p  <- c(0.000015, -1.5)
PHYS_scalar <- c(2, -2, NA, NA, NA, NA)
FI      <- c(0.05, 0.05)
fmet_p <- c(1, 0.85, 0.013)
depth       <- 30
h2y         <- 24*365
MICROtoECO  <- depth * 1e4 * 1e6 / 1e6   	#mgC/cm3 to kgC/km2

#Set default multipliers
Tau_MULT = 1
desorb_MULT = 1
fPHYS_MULT = 1

########################################
# Apply parameter multipliers
########################################
# Vslope = Vslope * 3.85
# Vint = Vint * 1.16
# Kslope = Kslope * 0.84
# Kint = Kint * 1.55
# CUE = CUE * 0.46
# Tau_MULT = 0.54
# desorb_MULT = 0.08
# fPHYS_MULT = 0.34

#Litter bag data
LITtype  <- c('TRAEf', 'PIREf','THPLf','ACSAf','QUPRf','DRGLf')
bagMET   <- c(10.6, 36.2, 37.4, 56.8, 37.1, 49.3) #from Gordon's LitterCharacteristics.txt
bagLIG   <- c(16.2, 19.2, 26.7, 15.9, 23.5, 10.9) # % from Gordon's LitterCharacteristics.txt
bagN     <- c(0.38, 0.59, 0.62, 0.81, 1.03, 1.97) # %N 
bagCN    <- c(133.3,92.7, 83.1, 61.8, 50.5, 24.2)
calcN    <- (1 / bagCN) / 2.5 * 100    
calcMET  <- 0.85 - 0.013 * bagLIG/calcN 			#as calculated in DAYCENT
bagMET   <- bagMET / 100
bagMET   <- calcMET
BAG_LIT  <- 100  * 1e3 / 1e4

#Litter bag partitioned to layers
BAG      <- array(NA, dim=c(6,2))  #litter BAG inputs to MET/STR
for (i in 1:6) {
  BAG[i,1]   <- (BAG_LIT / depth) * bagMET[i]     
  BAG[i,2]   <- (BAG_LIT / depth) * (1-bagMET[i])
}


###########################################
# MIMICS single point function
###########################################
MIMICS1 <- function(df){

  #DEBUG
  #df <- data[1,]
  
  note = ''
  
  fMET        <- mean(calcMET)
  TSOI        <- df$tsoi  
  ANPP        <- df$ANPP/2
  fCLAY       <- df$clay/100
  LIG_N       <- df$lig_N
  EST_LIT_in  <- ANPP/(365*24)   		#gC/m2/h (from g/m2/y, Knapp et al. Science 2001)
  BAG_LIT_in  <- 100      					#gC/m2/h
  EST_LIT     <- EST_LIT_in  * 1e3 / 1e4    #mgC/cm2/h 
  BAG_LIT     <- BAG_LIT_in  * 1e3 / 1e4    
  
###ADJUSTING VSLOPE AND KINT by MAT
  Vslope_ADJ = Vslope * ((0.0145*TSOI) + 0.5972)
  Kint_ADJ = Kint * ((-0.0255*TSOI) + 1.096)
  
  
  #-----------------caclulate parameters---------------------------
  #Calculate Vmax & (using parameters from German 2012, as in Wieder et al. 2013 Nature Climate Change)
  Vmax     <- exp(TSOI * Vslope_ADJ + Vint) * aV
  Km       <- exp(Kslope * TSOI + Kint_ADJ) * aK  
  
  Tao_MOD1 <- sqrt(ANPP/100)  #basicaily standardize against NWT
  tao      <- c(5.2e-4*exp(0.3*fMET), 2.4e-4*exp(0.1*fMET))	
  tau      <- tao * Tao_MOD1 * Tau_MULT
  
  fPHYS    <- c(fPHYS_r[1] * exp(fPHYS_r[2]*fCLAY), 
                fPHYS_K[1] * exp(fPHYS_K[2]*fCLAY)) 
  fCHEM    <- c(fCHEM_r[1] * exp(fCHEM_r[2]*fMET) * fCHEM_r[3], 
                fCHEM_K[1] * exp(fCHEM_K[2]*fMET) * fCHEM_K[3]) 	
  fAVAI    <- 1 - (fPHYS + fCHEM)
  
  desorb   <- fSOM_p[1] * exp(fSOM_p[2]*(fCLAY))       
  
  desorb <- desorb * desorb_MULT
  fPHYS <- fPHYS * fPHYS_MULT
  
  #cMAX     <- 1.4                    #ORIG 1.4 Maximum CHEM SOM scalar w/   0% Clay 
  #cMIN     <- 1.2                    #ORIG 1.4 Minimum CHEM SOM scalar w/ 100% Clay 
  #cSLOPE   <- cMIN - cMAX            #Slope of linear function of cSCALAR for CHEM SOM  
  
  pSCALAR  <- PHYS_scalar[1] * exp(PHYS_scalar[2]*(sqrt(fCLAY)))  #Scalar for texture effects on SOMp
  
  ###########################################################
  
  #------------!!MODIFIERS AS IN MIMICS2_b!!---------------
  v_MOD    <- vMOD  
  k_MOD    <- kMOD 
  k_MOD[3] <- k_MOD[3] * pSCALAR    
  k_MOD[6] <- k_MOD[6] * pSCALAR    
  
  VMAX     <- Vmax * v_MOD 
  KM       <- Km / k_MOD
  
  
  #Litter bag partitioned to layers
  BAG      <- array(NA, dim=c(6,2))  #litter BAG inputs to MET/STR
  for (i in 1:6) {
    BAG[i,1]   <- (BAG_LIT / depth) * bagMET[i]     
    BAG[i,2]   <- (BAG_LIT / depth) * (1-bagMET[i])
  }

  
  #initialize pools
  I       <- array(NA, dim=2)             
  I[1]    <- (EST_LIT / depth) * fMET     
  I[2]    <- (EST_LIT / depth) * (1-fMET)
  LIT     <- I   
  MIC     <- I  
  SOM     <- rep(NA, 3) 
  SOM[1]  <- I[1]
  SOM[2]  <- I[2]
  SOM[3]  <- I[1] 
  LITmin  <- rep(NA, dim=4)
  MICtrn  <- rep(NA, dim=6)
  SOMmin  <- rep(NA, dim=2)
  DEsorb  <- rep(NA, dim=1)
  OXIDAT  <- rep(NA, dim=1)
  
  
  #Calculate RXEQ pools  
  Tpars <- c( I = I, VMAX = VMAX, KM = KM, CUE = CUE, 
              fPHYS = fPHYS, fCHEM = fCHEM, fAVAI = fAVAI, FI = FI, 
              tao = tau, LITmin = LITmin, SOMmin = SOMmin, MICtrn = MICtrn, 
              desorb = desorb, DEsorb = DEsorb, OXIDAT = OXIDAT, KO = KO)
  Ty    <- c( LIT_1 = LIT[1], LIT_2 = LIT[2], 
              MIC_1 = MIC[1], MIC_2 = MIC[2], 
              SOM_1 = SOM[1], SOM_2 = SOM[2], SOM_3 = SOM[3])
  
  ## Set global parameters to allow pass to stode function
  .GlobalEnv$VMAX <- VMAX
  .GlobalEnv$KM <- KM
  .GlobalEnv$fPHYS <- fPHYS
  .GlobalEnv$fCHEM <- fCHEM
  .GlobalEnv$fAVAI <- fAVAI
  .GlobalEnv$I <- I
  .GlobalEnv$tao <- tau
  .GlobalEnv$LITmin <- LITmin
  .GlobalEnv$SOMmin <- SOMmin
  .GlobalEnv$MICtrn <- MICtrn
  .GlobalEnv$desorb <- desorb
  .GlobalEnv$DEsorb <- DEsorb
  .GlobalEnv$OXIDAT <- OXIDAT
  
  #Default stode
  test  <- stode(y = Ty, time = 1e6, fun = XEQ, parms = Tpars, positive = TRUE)

  ### Calc and get MIMICS output 
  MIMLIT    <- (test[[1]][[1]]+test[[1]][[2]])  * depth *1e4 / 1e6 #convert kgC/m2 from mgC/cm3 (0-30 cm) 
  MIMMIC    <- (test[[1]][[3]]+test[[1]][[4]])  * depth *1e4 / 1e6
  MIM_CO    <-  test[[1]][[3]]/test[[1]][[4]]
  MIMSOC    <- sum(test[[1]])  * depth *1e4 / 1e6   
  
  MIMout <- list()
  MIMout[[1]] <- as.numeric(test[[1]])
  MIMout[[2]] <- Tpars
  MIMout[[3]] <- df
  MIMout[[4]] <- BAG
  
  return(MIMout)
}


# ##############################################
#   Spin steady state forward ftn
# ##############################################

MIMsteady <- function(MIMstode, nday, make_plot=T){
  
  #DEBUG
  #MIMstode = MIMout_single
  #nday = 365 * 30
  #make_plot=T
  
  #Grab forcings
  Site   <- MIMstode[[3]][1]
  ANPP   <- MIMstode[[3]][2]
  tsoi   <- MIMstode[[3]][3]
  clay   <- MIMstode[[3]][4]/100
  lig_N  <- MIMstode[[3]][5]
  
  # Set run time vars
  day    <- seq(1,nday,1)
  year   <- day/365
  doy    <- 1
  
  #Init pool arrays
  LIT    <- array(NA, dim = c(2,nday))
  MIC    <- array(NA, dim = c(2,nday))
  SOM    <- array(NA, dim = c(3,nday))
  
  #Init pools
  LIT_1    <- MIMstode[[1]][1]   
  LIT_2    <- MIMstode[[1]][2]
  MIC_1    <- MIMstode[[1]][3]
  MIC_2    <- MIMstode[[1]][4]
  SOM_1    <- MIMstode[[1]][5]
  SOM_2    <- MIMstode[[1]][6]
  SOM_3    <- MIMstode[[1]][7]
  
  I        <- MIMstode[[2]][1:2]
  VMAX     <- MIMstode[[2]][3:8]
  KM       <- MIMstode[[2]][9:14]
  CUE      <- MIMstode[[2]][15:18]
  fPHYS    <-  MIMstode[[2]][19:20]
  fCHEM    <-  MIMstode[[2]][21:22]
  fAVAI    <-  MIMstode[[2]][23:24]
  fI       <-  MIMstode[[2]][25:26]
  tao      <- MIMstode[[2]][27:28]
  desorb   <- MIMstode[[2]][32]
  DEsorb   <- MIMstode[[2]][33]
  OXIDAT   <- MIMstode[[2]][34]
  KO       <- MIMstode[[2]][35:36]
  
  # Run MIMICS hourly for ndays
  for (d in 1:nday)  {
    for (h in 1:24)   {
      
      #Fluxes at each time step
      LITmin  <- rep(NA, dim=4)
      MICtrn  <- rep(NA, dim=6)
      SOMmin  <- rep(NA, dim=2)
      DEsorb  <- rep(NA, dim=1)
      OXIDAT  <- rep(NA, dim=1)
      
      #Flows to and from MIC_1
      LITmin[1] = MIC_1 * VMAX[1] * LIT_1 / (KM[1] + LIT_1)   #MIC_1 decomp of MET lit
      LITmin[2] = MIC_1 * VMAX[2] * LIT_2 / (KM[2] + LIT_2)   #MIC_1 decomp of STRUC lit
      SOMmin[1] = MIC_1 * VMAX[3] * SOM_3 / (KM[3] + SOM_3)   #Decomp of SOMa by MIC_1
      
      MICtrn[1] = MIC_1 * tao[1]  * fPHYS[1]                  #MIC_1 turnover to SOMp
      MICtrn[2] = MIC_1 * tao[1]  * fCHEM[1]                  #MIC_1 turnover to SOMc  
      MICtrn[3] = MIC_1 * tao[1]  * fAVAI[1]                  #MIC_1 turnover to SOMa 
      
      #Flows to and from MIC_2
      LITmin[3] = MIC_2 * VMAX[4] * LIT_1 / (KM[4] + LIT_1)   #decomp of MET litter
      LITmin[4] = MIC_2 * VMAX[5] * LIT_2 / (KM[5] + LIT_2)   #decomp of SRUCTURAL litter
      SOMmin[2] = MIC_2 * VMAX[6] * SOM_3 / (KM[6] + SOM_3)   #decomp of PHYSICAL SOM by MIC_1
      
      MICtrn[4] = MIC_2 * tao[2]  * fPHYS[2]                  #MIC_2 turnover to SOMp 
      MICtrn[5] = MIC_2 * tao[2]  * fCHEM[2]                  #MIC_2 turnover to SOMc 
      MICtrn[6] = MIC_2 * tao[2]  * fAVAI[2]                  #MIC_2 turnover to SOMa  
      
      
      DEsorb    = SOM_1 * desorb  #* (MIC_1 + MIC_2)  #desorbtion of PHYS to AVAIL (function of fCLAY)
      OXIDAT    = ((MIC_2 * VMAX[5] * SOM_2 / (KO[2]*KM[5] + SOM_2)) +
                     (MIC_1 * VMAX[2] * SOM_2 / (KO[1]*KM[2] + SOM_2)))  #oxidation of C to A
      
      
      LIT_1 = LIT_1 + I[1]*(1-FI[1]) - LITmin[1] - LITmin[3]
      MIC_1 = MIC_1 + CUE[1]*(LITmin[1]+ SOMmin[1]) + CUE[2]*(LITmin[2]) - sum(MICtrn[1:3])
      SOM_1 = SOM_1 + I[1]*FI[1] + MICtrn[1] + MICtrn[4]- DEsorb 
      
      LIT_2 = LIT_2 + I[2] * (1-FI[2]) - LITmin[2] - LITmin[4]
      MIC_2 = MIC_2 + CUE[3]*(LITmin[3]+ SOMmin[2]) + CUE[4]*(LITmin[4]) - sum(MICtrn[4:6])  
      SOM_2 = SOM_2 + I[2]*FI[2] + MICtrn[2] + MICtrn[5] - OXIDAT
      
      SOM_3 = SOM_3 + MICtrn[3] + MICtrn[6] + DEsorb + OXIDAT - SOMmin[1] - SOMmin[2]	
  
      #write out daily results
      if (h == 24) {
        LIT[1,d] <- LIT_1
        LIT[2,d] <- LIT_2
        MIC[1,d] <- MIC_1
        MIC[2,d] <- MIC_2
        SOM[1,d] <- SOM_1
        SOM[2,d] <- SOM_2
        SOM[3,d] <- SOM_3
        
      #advancy day of year counter
      if (doy == 365) {
        doy <- 1 
        print(paste0(c(Site, "- Finished initial year:", year[d])))
      } else {
        doy <- doy + 1
        }               #close day of year counter
      }	   						#close daily results counter
      
      #    remove(Vmax, VMAX, Km, KM)
    }							#close hour loop
  }								#close daily loop
  
  #Output plot
  if(make_plot) {
    par(mfrow=c(3,1), mar=c(4,4,1,1))
    test <- MIMstode[[1]]
    
    plot(year,  LIT[1,], lwd=3, ylim=c(min(LIT)*0.7, max(LIT))*1.15, type ="l", xlab="", main = paste(Site))
    lines(year, LIT[2,], lwd=3, col = 2)
    abline(h=test[1], col=1, lty=2)
    abline(h=test[2], col=2, lty=2)
    legend("topright", legend=c("Met","Struc"), col=c(1,2), lty = 1, 
           lwd = 3, cex=1.3, bty="n")
    
    plot(year,  MIC[1,], lwd=3, type ="l", xlab="", ylim=c(min(MIC)*0.7, max(MIC))*1.15, )
    lines(year, MIC[2,], lwd=3, col = 2) 
    abline(h=test[3], col=1, lty=2)
    abline(h=test[4], col=2, lty=2)
    legend("topright", legend=c("Mic_r","Mic_K"), col=c(1,2), lty = 1, 
           lwd = 3, cex=1.3, bty="n")
    
    plot(year,  SOM[1,], lwd=3, type ="l", ylim=c(min(SOM)*0.7, max(SOM))*1.15, )
    lines(year, SOM[2,], lwd=3, col = 2)
    lines(year, SOM[3,], lwd=3, col = 4)
    abline(h=test[5], col=1, lty=2)
    abline(h=test[6], col=2, lty=2)
    abline(h=test[7], col=4, lty=2)
    legend("topright", legend=c("Phys","Chem","Avail"), col=c(1,2,4), lty = 1, 
           lwd = 3, cex=1.3, bty="n")
  }
  
  return(list(LIT, MIC, SOM))
}


# ----------------------------------------------------------
# (D)     start litter bag experiment
#            add litter Oct 1, d=144
# ----------------------------------------------------------

litbag_decomp <- function(MIMss, litBAG, nspin_yrs=10, nspin_days=200, litadd_day=143){

  #DEBUG
  #MIMss <- MIMout_SS
  #litBAG <- BAG[1,]
  #nspin_yrs <- 10
  #nspin_days <- 200
  #litadd_day <- 143
  
  nday   <- 365 * nspin_yrs + nspin_days
  day    <- seq(1,nday,1)
  year   <- (day-litadd_day)/365
  doy    <- 1
  
  #Init arrays to store daily output data
  LIT    <- array(NA, dim = c(2,nday))#, dimnames = list(c("LITpool", "cumDAY")))
  LITBAG <- array(NA, dim = c(2,nday))#, dimnames = list(rep(c("BAGpool", "cumDAY"),6), c("A"), c("B"))) #for litter bag study
  MIC    <- array(NA, dim = c(2,nday))#, dimnames = list(c("MICpool", "cumDAY")))
  SOM    <- array(NA, dim = c(3,nday))#, dimnames = list(c("SOMpool", "SOMpool2", "cumDAY")))
  
  #Init MIMICS pools
  MIMout_SS_n <- length(MIMout_SS[[1]][1,])
  
  LIT_1    <- MIMout_SS[[1]][1,MIMout_SS_n]   
  LIT_2    <- MIMout_SS[[1]][2,MIMout_SS_n] 
  LITbag_1 <- MIMout_SS[[1]][1,MIMout_SS_n] 
  LITbag_2 <- MIMout_SS[[1]][2,MIMout_SS_n] 
  MIC_1    <- MIMout_SS[[2]][1,MIMout_SS_n] 
  MIC_2    <- MIMout_SS[[2]][2,MIMout_SS_n] 
  MICbag_1 <- MIMout_SS[[2]][1,MIMout_SS_n]
  MICbag_2 <- MIMout_SS[[2]][2,MIMout_SS_n]
  SOM_1    <- MIMout_SS[[3]][1,MIMout_SS_n]
  SOM_2    <- MIMout_SS[[3]][2,MIMout_SS_n]
  SOM_3    <- MIMout_SS[[3]][3,MIMout_SS_n]
  
  for (d in 1:nday)  {
    for (h in 1:24)   {
      #Fluxes at each time step
      LITmin  <- rep(NA, dim=4)
      LITbag  <- rep(NA, dim=4)
      MICtrn  <- rep(NA, dim=6)
      SOMmin  <- rep(NA, dim=2)
      DEsorb  <- rep(NA, dim=1)
      OXIDAT  <- rep(NA, dim=1)
      
      #Flows to and from MIC_1
      LITmin[1] = MIC_1 * VMAX[1] * LIT_1 / (KM[1] + LIT_1)   #MIC_1 decomp of MET lit
      LITmin[2] = MIC_1 * VMAX[2] * LIT_2 / (KM[2] + LIT_2)   #MIC_1 decomp of STRUC lit
      LITbag[1] = MIC_1 * VMAX[1] * LITbag_1 / (KM[1] + LITbag_1)   #MIC_1 mineralization of METABOLIC litter
      LITbag[2] = MIC_1 * VMAX[2] * LITbag_2 / (KM[2] + LITbag_2)   #MIC_1 mineralization of STRUC litter
      SOMmin[1] = MIC_1 * VMAX[3] * SOM_3 / (KM[3] + SOM_3)   #Decomp of SOMa by MIC_1
      
      MICtrn[1] = MIC_1 * tao[1]  * fPHYS[1]                  #MIC_1 turnover to SOMp
      MICtrn[2] = MIC_1 * tao[1]  * fCHEM[1]                  #MIC_1 turnover to SOMc  
      MICtrn[3] = MIC_1 * tao[1]  * fAVAI[1]                  #MIC_1 turnover to SOMa 
      
      #Flows to and from MIC_2
      LITmin[3] = MIC_2 * VMAX[4] * LIT_1 / (KM[4] + LIT_1)   #decomp of MET litter
      LITmin[4] = MIC_2 * VMAX[5] * LIT_2 / (KM[5] + LIT_2)   #decomp of SRUCTURAL litter
      LITbag[3] = MIC_2 * VMAX[4] * LITbag_1 / (KM[4] + LITbag_1)   #mineralization of MET litter
      LITbag[4] = MIC_2 * VMAX[5] * LITbag_2 / (KM[5] + LITbag_2)   #mineralization of SRUCTURAL litter
      SOMmin[2] = MIC_2 * VMAX[6] * SOM_3 / (KM[6] + SOM_3)   #decomp of PHYSICAL SOM by MIC_1
      
      MICtrn[4] = MIC_2 * tao[2]  * fPHYS[2]                  #MIC_2 turnover to SOMp 
      MICtrn[5] = MIC_2 * tao[2]  * fCHEM[2]                  #MIC_2 turnover to SOMc 
      MICtrn[6] = MIC_2 * tao[2]  * fAVAI[2]                  #MIC_2 turnover to SOMa  
      
      
      DEsorb    = SOM_1 * desorb  #* (MIC_1 + MIC_2)  #desorbtion of PHYS to AVAIL (function of fCLAY)
      OXIDAT    = (MIC_2 * VMAX[5] * SOM_2 / (KO[2]*KM[5] + SOM_2)) 
      + (MIC_1 * VMAX[2] * SOM_2 / (KO[1]*KM[2] + SOM_2))  #oxidation of C to A
      
      
      LIT_1 = LIT_1 + I[1]*(1-FI[1]) - LITmin[1] - LITmin[3]
      LITbag_1 <- LITbag_1 + I[1]*(1-FI[1]) - LITbag[1] - LITbag[3]
      MIC_1 = MIC_1 + CUE[1]*(LITmin[1]+ SOMmin[1]) + CUE[2]*(LITmin[2]) - sum(MICtrn[1:3])
      SOM_1 = SOM_1 + I[1]*FI[1] + MICtrn[1] + MICtrn[4]- DEsorb 
      
      LIT_2 = LIT_2 + I[2] * (1-FI[2]) - LITmin[2] - LITmin[4]
      LITbag_2 <- LITbag_2 + I[2] * (1-FI[2]) - LITbag[2] - LITbag[4]
      MIC_2 = MIC_2 + CUE[3]*(LITmin[3]+ SOMmin[2]) + CUE[4]*(LITmin[4]) - sum(MICtrn[4:6])  
      SOM_2 = SOM_2 + I[2]*FI[2] + MICtrn[2] + MICtrn[5] - OXIDAT
      
      SOM_3 = SOM_3 + MICtrn[3] + MICtrn[6] + DEsorb + OXIDAT - SOMmin[1] - SOMmin[2]	
      #add litter bag on Oct 1
      if (d == 143)   {
        if (h == 24)  {
          LITbag_1 <- LITbag_1 + litBAG[1]
          LITbag_2 <- LITbag_2 + litBAG[2]
          #print(paste("------added litter",LITtype[i]))
        }
      }
      
      #write out daily results
      if (h == 24) {
        LIT[1,d] <- LIT_1
        LIT[2,d] <- LIT_2
        LITBAG[1,d]  <- LITbag_1
        LITBAG[2,d]  <- LITbag_2
        MIC[1,d] <- MIC_1
        MIC[2,d] <- MIC_2
        SOM[1,d] <- SOM_1
        SOM[2,d] <- SOM_2
        SOM[3,d] <- SOM_3
        
        #advancy day of year counter
        if (doy == 365) {
          doy <- 1 
          #	    print(paste(c(strSite[s], "finished year", year[d])))
        } else {
          doy <- doy + 1
        }                         #close day of year counter
      }	   						#close daily results counter
      
      #    remove(Vmax, VMAX, Km, KM)
    }							#close hour loop
  }								#close daily loop
  
  print(paste(c("finished litter", LITtype[i])))
  return(list(LIT, LITBAG, MIC, SOM))
}


#--------------------------------------------------
# RUN LITBAG DECOMP FOR A SINGLE SITE
#--------------------------------------------------

LTER <- read.csv("LTER_SITE_1.csv")
site_no <- 7
Site <- LTER$Site[site_no]
data <- data.frame(Site = LTER$Site[site_no],
                   ANPP = LTER$ANPP[site_no],
                   tsoi = LTER$MAT[site_no],
                   clay = LTER$CLAY2[site_no],
                   lig_N = 0)

MIMout_single <- MIMICS1(data[1,])
MIMout_SS <- MIMsteady(MIMout_single, nday=365*30)
singleBag_decomp <- litbag_decomp(MIMout_SS, BAG[1,])

# Run for all bags
All_bagsDecomp <- list()
for(i in 1:nrow(BAG)){
  All_bagsDecomp[[i]] <- litbag_decomp(MIMout_SS, BAG[i,])
}


#--------------------------------------------------
# Calculate average litter decomp
#--------------------------------------------------

# Get LIT totals
LIT <- as.data.frame(All_bagsDecomp[[1]][1])
allLIT  <- colSums(LIT, dims = 1)

# Init arrays
nday   <- 365 * 10 + 200
day    <- seq(1,nday,1)
year   <- (day-143)/365
doy    <- 1

allBAG  <- array(NA, dim=c(6,nday))
difBAG  <- array(NA, dim=c(6,nday))
maxBAG  <- array(NA, dim=c(6,nday))
BAGleft <- array(NA, dim=c(6,nday), dimnames=list(LITtype,c(as.character(year))))

# Calc how much of the LIT bag remains
for (i in 1:length(All_bagsDecomp)) {
  allBAG[i,]  <- colSums(as.data.frame(All_bagsDecomp[[i]][2]))  #SPEED UP ANALYSIS change LITBAG[i ] to 1
  difBAG[i,]  <- allBAG[i,] - allLIT
  maxBAG      <- max(difBAG[i,])
  BAGleft[i,] <-  100* difBAG[i,] / maxBAG        #by taking LIT + BAG mineralization at each step
}

# Remove days before lit decomp start
BAGleft[1:6,1:143] <- NA

# Store mean decomp results in dataframe
mean_bag_decomp <- data.frame(yr = year,
                              meanBAG = colMeans(BAGleft),
                              sdBAG <- apply(BAGleft, 2, sd))


#--------------------------------------------------
# Plot it!
#--------------------------------------------------

ggplot(mean_bag_decomp, aes(x=yr, y=meanBAG)) + 
  geom_line(size=1.5) +
  ggtitle(paste("LTER SITE: ", Site, " - Leaf litter decomp")) +
  ylab("Mass remaining (%)") + xlab("Time (y)") +
  theme_minimal()


### Plot each litter type decomp
bag_decomp <- NULL
for(i in 1:length(row.names(BAGleft))){
  df <- data.frame(yr = year,
                   type = row.names(BAGleft)[i],
                   bag = BAGleft[i,])
  bag_decomp <- rbind(bag_decomp, df)
}

#png("LTER_site_LIT_decomp.png", width = 1600, height = 1200, res = 300)
ggplot(bag_decomp, aes(x=yr, y=bag, color=type)) + 
  geom_line(size=1.5, alpha=1) +
  geom_line(data=mean_bag_decomp, aes(x=yr, y=meanBAG, color=' MEAN'), linetype = "dashed", size=1) +
  scale_color_manual(values=c("#000032", "#137177", "#188977", "#1D9A6C", "#39A96B", "#74C67A", "#99D492")) +
  #ggtitle(paste(Site, " - Leaf litter decomp")) +
  ggtitle(paste("H.J. Andrews EF", " - Leaf litter decomp")) +
  ylab("Mass remaining (%)") + xlab("Time (y)") +
  labs(color = "Litter Type") +
  theme_minimal()
#dev.off()


#--------------------------------------------------------------
# Collect data for all litter type bags across all LTER sites
#--------------------------------------------------------------
LTER <- read.csv("LTER_SITE_1.csv")
LTER_bag_decomp <- NULL
for(s in 1:nrow(LTER)) {
  site_no <- s
  Site <- LTER$Site[site_no]
  data <- data.frame(Site = LTER$Site[site_no],
                     ANPP = LTER$ANPP[site_no],
                     tsoi = LTER$MAT[site_no],
                     clay = LTER$CLAY2[site_no],
                     lig_N = 0)

  MIMout_single <- MIMICS1(data[1,])
  MIMout_SS <- MIMsteady(MIMout_single, nday=365*30)


  # Run for all bags
  All_bagsDecomp <- list()
  for(i in 1:nrow(BAG)){
    All_bagsDecomp[[i]] <- litbag_decomp(MIMout_SS, BAG[i,])
  }

  #--------------------------------------------------
  # Calculate average litter decomp
  #--------------------------------------------------

  # Get LIT totals
  LIT <- as.data.frame(All_bagsDecomp[[1]][1])
  allLIT  <- colSums(LIT, dims = 1)

  # Init arrays
  nday   <- 365 * 10 + 200
  day    <- seq(1,nday,1)
  year   <- (day-143)/365
  doy    <- 1

  allBAG  <- array(NA, dim=c(6,nday))
  difBAG  <- array(NA, dim=c(6,nday))
  maxBAG  <- array(NA, dim=c(6,nday))
  BAGleft <- array(NA, dim=c(6,nday), dimnames=list(LITtype,c(as.character(year))))

  # Calc how much of the LIT bag remains
  for (i in 1:length(All_bagsDecomp)) {
    allBAG[i,]  <- colSums(as.data.frame(All_bagsDecomp[[i]][2]))  #SPEED UP ANALYSIS change LITBAG[i ] to 1
    difBAG[i,]  <- allBAG[i,] - allLIT
    maxBAG      <- max(difBAG[i,])
    BAGleft[i,] <-  100* difBAG[i,] / maxBAG        #by taking LIT + BAG mineralization at each step
  }

  # Remove days before lit decomp start
  BAGleft[1:6,1:143] <- NA
  
  # Store mean decomp results in dataframe
  # mean_bag_decomp <- data.frame(yr = year,
  #                               meanBAG = colMeans(BAGleft),
  #                               sdBAG <- apply(BAGleft, 2, sd))
  
  bag_decomp <- NULL
  for(i in 1:length(row.names(BAGleft))){
    df <- data.frame(yr = year,
                     site = Site,
                     type = row.names(BAGleft)[i],
                     bag = BAGleft[i,])
    bag_decomp <- rbind(bag_decomp, df)
  }
  
  LTER_bag_decomp <- rbind(LTER_bag_decomp, bag_decomp)
}

#Save the output
write.csv(LTER_bag_decomp, "output/LTER_MIMICS_Litter_Decomp_output_ADJ.csv")  
#(LTER_bag_decomp, "output/LTER_MIMICS_Litter_Decomp_output.rds")    


#--------------------------------------------------
# Plot mean bag decomp across sites
#--------------------------------------------------
meanBAGS <- LTER_bag_decomp %>% group_by(site, yr) %>%
              summarize(n = n(),
                        BagLeft_mn = mean(bag),
                        BagLeft_sd = sd(bag))

ggplot(meanBAGS, aes(x=yr, y=BagLeft_mn, color=site)) + 
  geom_line(size=1, alpha=0.8) +
  ggtitle("LTER SITES: MIMICS mean leaf litter decompostion") +
  ylab("Mass remaining (%)") + xlab("Time (y)") +
  labs(color = "Litter Type") +
  theme_minimal()  
ggsave('output/LTER_mean_decomp_ADJ.png', dpi=300, width = 5, height = 5)

ggplot(meanBAGS, aes(x=yr, y=BagLeft_sd, color=site)) + 
  geom_line(size=1, alpha=0.8) +
  ggtitle("LTER SITES: MIMICS mean leaf litter decompostion") +
  ylab("Mass remaining Std Dev (%)") + xlab("Time (y)") +
  labs(color = "Litter Type") +
  theme_minimal() 
#ggsave('output/LTER_stdev_decomp.tiff', dpi=300, width = 8, height = 5)
