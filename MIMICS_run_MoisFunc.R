#this script compares initial moisture function with moisture function from parameter calibration exercises

#FOR VSLOPE COMPARISONS
#initial function: Vmax * (0.7+ (0.3* soilGWC)); Km / (0.7+ (0.3* soilGWC))
#add Vslope multiplier that depends on GWC with equation: multiplier = 1.17 + 0.056*MS_GWC

#FOR VINT AND KINT COMPARISONS
#initial function: Vmax * (0.7+ (0.3* soilGWC)); Km / (0.7+ (0.3* soilGWC))
#add Vint multiplier: Vint multiplier = 1.1 + 0.04*MS_GWC

#FOR TAU COMPARISONS
#initial function: tau <- tau * Tau_MOD1 * Tau_MOD2 * Tau_MULT where Tau_MULT = 1
#add Tau_MULT multiplier above initial function: Tau_MULT = (1.63 - 0.29*MS_GWC) * Tau_MULT

#FOR TAU + vMOD COMPARISONS
#initial function: Tau_MULT and vMOD both = 1
#add Tau_MULT multiplier above initial function: Tau_MULT = (1.63 - 0.29*MS_GWC) * Tau_MULT
#add vMOD multiplier above initial function: vMOD <- vMOD * (1.6+0.15*MS_GWC)

# Load packages
library(dplyr) # used for pipes (%>%)
library(purrr) # used for map() function
library(furrr) # used for map() function

# Load MIMICS model ftn script
source("MIMICS_ftns/MIMICS_INC_daily.R")

# Load parameters for MIMICS
source("MIMICS_ftns/MIMICS_set_parameters.R")

#Vslope comparisons

#########################

####
# Example MIMICS runs
###

### Single-point run for different moisture functions

#for loop for default moisture functions
Sites <- c(rep('2', 5),rep('3', 5),rep('8', 5),rep('11',5),rep('17',5),rep('20', 5),rep('21', 5),rep('22', 5),rep('23', 5),
           rep('26', 5),rep('30', 5),rep('32', 5))
Incub_GWC <- rep(c(24.5, 8.5, 16.5, 32.5, 40.5), 12)
GWC <- c(rep(62.3, 5),rep(67.4, 5),rep(65.8, 5),rep(47.0, 5),rep(62.1, 5),rep(49.1, 5),rep(54.4, 5), rep(46.6, 5),
         rep(59.1, 5), rep(65.6, 5), rep(51.0, 5), rep(37.6, 5))
datalist = list()
for(i in 1:length(GWC)) {
df <- data.frame(ID = "BART",
                 SITE = Sites[i],
                 ANPP = 0,
                 MAT = 15,
                 CLAY = 3.3,
                 LIG = 10,
                 N = 0.8,
                 CN = 61,
                 fW = 0.5, #not used so doesn't matter
                 soilGWC = Incub_GWC[i],
                 MS_GWC = GWC[i])

# Run MIMICS and store output as MIMout
MIMout <- MIMICS_INC_DAILY(df[1,], ndays=3600)
MIMout$Site <- Sites[i]
MIMout$Inst_GWC <- Incub_GWC[i]
MIMout$hist_GWC <- GWC[i]
#hourly step

datalist[[i]] <- MIMout

}

MIMout_default_10 <- do.call(rbind, datalist) #with initial moisture function
MIMout_moist_10 <- do.call(rbind, datalist) #modified MIMICS_INC_daily to have moisture multiplier for Vslope only!
                                         #Maintained moisture function bc I think that's what I am supposed to do

###
# Plots
###
library(ggpubr)

# Litter mass loss
#can either limit to one moisture or create panels (currently moisutre = 60 WHC)
MD10_60 <- filter(MIMout_default_10, Inst_GWC == 24.5)
MM10_60 <- filter(MIMout_moist_10, Inst_GWC == 24.5)
plot_LIT_def <- ggplot(MD10_60, aes(y=LITs+LITm, x=DAY, color=Site, group=Site)) + geom_line(linewidth=1) +
  theme_bw() +
  ylab("Litter mass remaining (%)") +
  xlab("Incubation Time (days)") +
  labs(color = "Litter Pool")
  #facet_grid(.~Inst_GWC)
plot_LIT_def
#modified soil moisture function
plot_LIT_alt <- ggplot(MM10_60, aes(y=LITs+LITm, x=DAY, color=Site, group=Site)) + geom_line(linewidth=1) +
  theme_bw() +
  ylab("Litter mass remaining (%)") +
  xlab("Incubation Time (days)") +
  labs(color = "Litter Pool")
  #facet_grid(.~Inst_GWC)
plot_LIT_alt

#differences between mass loss at each site
MIM_def_105 <- MIMout_default_10 %>% filter(DAY==3600) %>% mutate(LIT = LITm+LITs) %>% mutate(site.ISM = paste(Site, Inst_GWC, sep = "."))
MIM_moi_105 <- MIMout_moist_10 %>% filter(DAY==3600) %>% mutate(LIT = LITm+LITs) %>% mutate(site.ISM = paste(Site, Inst_GWC, sep = "."))
MIMout_diff <- left_join(MIM_def_105, MIM_moi_105, by = "site.ISM")
MIMout_diff$LIT_dif <- (MIMout_diff$LIT.y - MIMout_diff$LIT.x)/MIMout_diff$LIT.x
#Site_moisture <- data.frame(Sites, GWC, Incub_GWC)
#MIMout_diff$Sites <- MIMout_diff$SITE.x
#MIMout_diff <- left_join(MIMout_diff, Site_moisture, by = "Sites")
MIMout_diff <- transform(MIMout_diff, SITE.x = reorder(SITE.x, LIT_dif))
MIMout_diff_60 <- filter(MIMout_diff, Inst_GWC.x == 24.5)
ggplot(MIMout_diff, aes(x=SITE.x, y=LIT_dif, fill=hist_GWC.x)) + geom_bar(stat = "identity") +
  ylab("Percent difference in litter mass loss (%)") + xlab("Microsite") + facet_grid(.~Inst_GWC.x)

#CO2 vs soil moisture
MIMout_diff$CO2_def <- MIMout_diff$CO2_MICr.x + MIMout_diff$CO2_MICK.x
MIMout_diff$CO2_moi <- MIMout_diff$CO2_MICr.y + MIMout_diff$CO2_MICK.y
ggplot(MIMout_diff, aes(x=Inst_GWC.x, y=CO2_def, shape = "Default", color=hist_GWC.x)) + geom_point(size=4) +geom_point(size =4, aes(x=Inst_GWC.x, y=CO2_moi, shape = "Alternate", color=hist_GWC.x)) +
  ylab("Total CO2 emited at day 3600") + theme_bw()

#######################################################################################


#Vint and Kint

##################################
#for loop for default moisture functions
Sites <- c(rep('2', 5),rep('3', 5),rep('8', 5),rep('11',5),rep('17',5),rep('20', 5),rep('21', 5),rep('22', 5),rep('23', 5),
           rep('26', 5),rep('30', 5),rep('32', 5))
Incub_GWC <- rep(c(24.5, 8.5, 16.5, 32.5, 40.5), 12)
GWC <- c(rep(62.3, 5),rep(67.4, 5),rep(65.8, 5),rep(47.0, 5),rep(62.1, 5),rep(49.1, 5),rep(54.4, 5), rep(46.6, 5),
         rep(59.1, 5), rep(65.6, 5), rep(51.0, 5), rep(37.6, 5))
datalist = list()
for(i in 1:length(GWC)) {
  df <- data.frame(ID = "BART",
                   SITE = Sites[i],
                   ANPP = 0,
                   MAT = 15,
                   CLAY = 3.3,
                   LIG = 10,
                   N = 0.8,
                   CN = 61,
                   fW = 0.5, #not used so doesn't matter
                   soilGWC = Incub_GWC[i],
                   MS_GWC = GWC[i])

  # Run MIMICS and store output as MIMout
  MIMout <- MIMICS_INC_DAILY(df[1,], ndays=3600)
  MIMout$Site <- Sites[i]
  MIMout$Inst_GWC <- Incub_GWC[i]
  MIMout$hist_GWC <- GWC[i]
  #hourly step

  datalist[[i]] <- MIMout

}

MIMout_default_10 <- do.call(rbind, datalist) #with initial moisture function
MIMout_moist_10 <- do.call(rbind, datalist) #modified MIMICS_INC_daily to have moisture multiplier for Vslope only!
#Maintained moisture function bc I think that's what I am supposed to do

###
# Plots
###
library(ggpubr)

# Litter mass loss
#can either limit to one moisture or create panels (currently moisutre = 60 WHC)
MD10_60 <- filter(MIMout_default_10, Inst_GWC == 24.5)
MM10_60 <- filter(MIMout_moist_10, Inst_GWC == 24.5)
plot_LIT_def <- ggplot(MD10_60, aes(y=LITs+LITm, x=DAY, color=Site, group=Site)) + geom_line(linewidth=1) +
  theme_bw() +
  ylab("Litter mass remaining (%)") +
  xlab("Incubation Time (days)") +
  labs(color = "Litter Pool")
#facet_grid(.~Inst_GWC)
plot_LIT_def
#modified soil moisture function
plot_LIT_alt <- ggplot(MM10_60, aes(y=LITs+LITm, x=DAY, color=Site, group=Site)) + geom_line(linewidth=1) +
  theme_bw() +
  ylab("Litter mass remaining (%)") +
  xlab("Incubation Time (days)") +
  labs(color = "Litter Pool")
#facet_grid(.~Inst_GWC)
plot_LIT_alt

#differences between mass loss at each site
MIM_def_105 <- MIMout_default_10 %>% filter(DAY==3600) %>% mutate(LIT = LITm+LITs) %>% mutate(site.ISM = paste(Site, Inst_GWC, sep = "."))
MIM_moi_105 <- MIMout_moist_10 %>% filter(DAY==3600) %>% mutate(LIT = LITm+LITs) %>% mutate(site.ISM = paste(Site, Inst_GWC, sep = "."))
MIMout_diff <- left_join(MIM_def_105, MIM_moi_105, by = "site.ISM")
MIMout_diff$LIT_dif <- (MIMout_diff$LIT.y - MIMout_diff$LIT.x)/MIMout_diff$LIT.x
#Site_moisture <- data.frame(Sites, GWC, Incub_GWC)
#MIMout_diff$Sites <- MIMout_diff$SITE.x
#MIMout_diff <- left_join(MIMout_diff, Site_moisture, by = "Sites")
MIMout_diff <- transform(MIMout_diff, SITE.x = reorder(SITE.x, LIT_dif))
MIMout_diff_60 <- filter(MIMout_diff, Inst_GWC.x == 24.5)
ggplot(MIMout_diff, aes(x=SITE.x, y=LIT_dif, fill=hist_GWC.x)) + geom_bar(stat = "identity") +
  ylab("Percent difference in litter mass loss (%)") + xlab("Microsite") + facet_grid(.~Inst_GWC.x)

#CO2 vs soil moisture
MIMout_diff$CO2_def <- MIMout_diff$CO2_MICr.x + MIMout_diff$CO2_MICK.x
MIMout_diff$CO2_moi <- MIMout_diff$CO2_MICr.y + MIMout_diff$CO2_MICK.y
ggplot(MIMout_diff, aes(x=Inst_GWC.x, y=CO2_def, shape = "Default", color=hist_GWC.x)) + geom_point(size=4) +geom_point(size =4, aes(x=Inst_GWC.x, y=CO2_moi, shape = "Alternate", color=hist_GWC.x)) +
  ylab("Total CO2 emited at day 3600") + theme_bw()
ggplot(MIMout_diff, aes(x=hist_GWC.x, y=CO2_def, shape = "Default")) + geom_point(size=4) +geom_point(size =4, aes(x=hist_GWC.x, y=CO2_moi, shape = "Alternate")) +
  ylab("Total CO2 emited at day 3600") + theme_bw()
######################################

#vMOD and kMOD
###################################

#for loop for default moisture functions
Sites <- c(rep('2', 5),rep('3', 5),rep('8', 5),rep('11',5),rep('17',5),rep('20', 5),rep('21', 5),rep('22', 5),rep('23', 5),
           rep('26', 5),rep('30', 5),rep('32', 5))
Incub_GWC <- rep(c(24.5, 8.5, 16.5, 32.5, 40.5), 12)
GWC <- c(rep(62.3, 5),rep(67.4, 5),rep(65.8, 5),rep(47.0, 5),rep(62.1, 5),rep(49.1, 5),rep(54.4, 5), rep(46.6, 5),
         rep(59.1, 5), rep(65.6, 5), rep(51.0, 5), rep(37.6, 5))
datalist = list()
for(i in 1:length(GWC)) {
  df <- data.frame(ID = "BART",
                   SITE = Sites[i],
                   ANPP = 0,
                   MAT = 15,
                   CLAY = 3.3,
                   LIG = 10,
                   N = 0.8,
                   CN = 61,
                   fW = 0.5, #not used so doesn't matter
                   soilGWC = Incub_GWC[i],
                   MS_GWC = GWC[i])

  # Run MIMICS and store output as MIMout
  MIMout <- MIMICS_INC_DAILY(df[1,], ndays=3600)
  MIMout$Site <- Sites[i]
  MIMout$Inst_GWC <- Incub_GWC[i]
  MIMout$hist_GWC <- GWC[i]
  #hourly step

  datalist[[i]] <- MIMout

}

MIMout_default_10 <- do.call(rbind, datalist) #with initial moisture function
MIMout_moist_10 <- do.call(rbind, datalist) #modified MIMICS_INC_daily to have moisture multiplier for Vslope only!
#Maintained moisture function bc I think that's what I am supposed to do
########################################################################################

#tau or tau+vMOD or tau+CUE
#############################################

###
#BART
###

#for loop for default moisture functions
Sites <- c(rep('2', 5),rep('3', 5),rep('8', 5),rep('11',5),rep('17',5),rep('20', 5),rep('21', 5),rep('22', 5),rep('23', 5),
           rep('26', 5),rep('30', 5),rep('32', 5))
Incub_GWC <- rep(c(24.5, 8.5, 16.5, 32.5, 40.5), 12)
GWC <- c(rep(62.3, 5),rep(67.4, 5),rep(65.8, 5),rep(47.0, 5),rep(62.1, 5),rep(49.1, 5),rep(54.4, 5), rep(46.6, 5),
         rep(59.1, 5), rep(65.6, 5), rep(51.0, 5), rep(37.6, 5))
datalist = list()
for(i in 1:length(GWC)) {
  df <- data.frame(ID = "BART",
                   SITE = Sites[i],
                   ANPP = 0,
                   MAT = 15,
                   CLAY = 3.3,
                   LIG = 10,
                   N = 0.8,
                   CN = 61,
                   fW = 0.5, #not used so doesn't matter
                   soilGWC = Incub_GWC[i],
                   MS_GWC = GWC[i])

  # Run MIMICS and store output as MIMout
  MIMout <- MIMICS_INC_DAILY(df[1,], ndays=3600)
  MIMout$Site <- Sites[i]
  MIMout$Inst_GWC <- Incub_GWC[i]
  MIMout$hist_GWC <- GWC[i]
  #hourly step

  datalist[[i]] <- MIMout

}

#need to run the sourcing code & loop above with and without new equation and then choose the right one below
MIMout_default_10 <- do.call(rbind, datalist) #with initial moisture function
MIMout_moist_10 <- do.call(rbind, datalist) #with tau multiplier or tau+vMOD multiplier


###
#GRSM
###

#for loop for default moisture functions
Sites <- c(rep('1', 5),rep('3', 5),rep('5', 5),rep('8',5),rep('14',5),rep('18', 5),rep('19', 5),rep('20', 5),rep('22', 5),
           rep('25', 5),rep('26', 5),rep('27', 5))
Incub_GWC <- rep(c(24.5, 8.5, 16.5, 32.5, 40.5), 12)
GWC <- c(rep(43, 5),rep(48, 5),rep(28, 5),rep(44, 5),rep(22, 5),rep(35, 5),rep(33, 5), rep(33, 5),
         rep(37, 5), rep(40, 5), rep(54, 5), rep(29, 5))
datalist = list()
for(i in 1:length(GWC)) {
  df <- data.frame(ID = "GRSM",
                   SITE = Sites[i],
                   ANPP = 0,
                   MAT = 15,
                   CLAY = 20,
                   LIG = 10,
                   N = 0.8,
                   CN = 61,
                   fW = 0.5, #not used so doesn't matter
                   soilGWC = Incub_GWC[i],
                   MS_GWC = GWC[i])

  # Run MIMICS and store output as MIMout
  MIMout <- MIMICS_INC_DAILY(df[1,], ndays=3600)
  MIMout$Site <- Sites[i]
  MIMout$Inst_GWC <- Incub_GWC[i]
  MIMout$hist_GWC <- GWC[i]
  #hourly step

  datalist[[i]] <- MIMout

}

#need to run the sourcing code & loop above with and without new equation and then choose the right one below
MIMout_default_10 <- do.call(rbind, datalist) #with initial moisture function
MIMout_moist_10 <- do.call(rbind, datalist) #with tau multiplier or tau+vMOD multiplier


###
#TREE
###

#tau and CUE together doesn't work - creating crazy high CUE even though mathematically it shouldn't be... must be bug in code...

#for loop for default moisture functions
Sites <- c(rep('3', 5),rep('4', 5),rep('8', 5),rep('11',5),rep('13',5),rep('14', 5),rep('15', 5),rep('19', 5),rep('21', 5),
           rep('31', 5),rep('34', 5),rep('35', 5))
Incub_GWC <- rep(c(24.5, 8.5, 16.5, 32.5, 40.5), 12)
GWC <- c(rep(54.9, 5),rep(32.5, 5),rep(79.4, 5),rep(47.5, 5),rep(36.9, 5),rep(29.5, 5),rep(20.4, 5), rep(23.2, 5),
         rep(18.0, 5), rep(17.9, 5), rep(15.8, 5), rep(16.6, 5))
datalist = list()
for(i in 1:length(GWC)) {
  df <- data.frame(ID = "TREE",
                   SITE = Sites[i],
                   ANPP = 0,
                   MAT = 15,
                   CLAY = 4.1,
                   LIG = 10,
                   N = 0.8,
                   CN = 61,
                   fW = 0.5, #not used so doesn't matter
                   soilGWC = Incub_GWC[i],
                   MS_GWC = GWC[i])

  # Run MIMICS and store output as MIMout
  MIMout <- MIMICS_INC_DAILY(df[1,], ndays=3600)
  MIMout$Site <- Sites[i]
  MIMout$Inst_GWC <- Incub_GWC[i]
  MIMout$hist_GWC <- GWC[i]
  #hourly step

  datalist[[i]] <- MIMout

}

#need to run the sourcing code & loop above with and without new equation and then choose the right one below
MIMout_default_10 <- do.call(rbind, datalist) #with initial moisture function
MIMout_moist_10 <- do.call(rbind, datalist) #with tau multiplier or tau+vMOD multiplier

######################################################################################################


###
# Plots
###
library(ggpubr)

# Litter mass loss
#can either limit to one moisture or create panels (currently moisutre = 60 WHC)
MD10_60 <- filter(MIMout_default_10, Inst_GWC == 24.5)
MM10_60 <- filter(MIMout_moist_10, Inst_GWC == 24.5)
plot_LIT_def <- ggplot(MD10_60, aes(y=LITs+LITm, x=DAY, color=Site, group=Site)) + geom_line(linewidth=1) +
  theme_bw() +
  ylab("Litter mass remaining (%)") +
  xlab("Incubation Time (days)") +
  labs(color = "Litter Pool")
#facet_grid(.~Inst_GWC)
plot_LIT_def
#modified soil moisture function
plot_LIT_alt <- ggplot(MM10_60, aes(y=LITs+LITm, x=DAY, color=Site, group=Site)) + geom_line(linewidth=1) +
  theme_bw() +
  ylab("Litter mass remaining (%)") +
  xlab("Incubation Time (days)") +
  labs(color = "Litter Pool")
#facet_grid(.~Inst_GWC)
plot_LIT_alt

#differences between mass loss at each site
#3600 is 10 yeas
#1000 for TREE since litter decomposes too quickly
MIM_def_105 <- MIMout_default_10 %>% filter(DAY==1000) %>% mutate(LIT = LITm+LITs) %>% mutate(site.ISM = paste(Site, Inst_GWC, sep = "."))
MIM_moi_105 <- MIMout_moist_10 %>% filter(DAY==1000) %>% mutate(LIT = LITm+LITs) %>% mutate(site.ISM = paste(Site, Inst_GWC, sep = "."))
MIMout_diff <- left_join(MIM_def_105, MIM_moi_105, by = "site.ISM")
MIMout_diff$LIT_dif <- (MIMout_diff$LIT.y - MIMout_diff$LIT.x)/MIMout_diff$LIT.x
#Site_moisture <- data.frame(Sites, GWC, Incub_GWC)
#MIMout_diff$Sites <- MIMout_diff$SITE.x
#MIMout_diff <- left_join(MIMout_diff, Site_moisture, by = "Sites")
MIMout_diff <- transform(MIMout_diff, SITE.x = reorder(SITE.x, LIT_dif))
MIMout_diff$CO2_dif <- ((MIMout_diff$CO2_MICr.y + MIMout_diff$CO2_MICK.y)
                        - (MIMout_diff$CO2_MICr.x + MIMout_diff$CO2_MICK.x))/(MIMout_diff$CO2_MICr.x + MIMout_diff$CO2_MICK.x)
MIMout_diff_60 <- filter(MIMout_diff, Inst_GWC.x == 24.5)
ggplot(MIMout_diff_60, aes(x=SITE.x, y=LIT_dif, fill=hist_GWC.x)) + geom_bar(stat = "identity") +
  ylab("Percent difference in litter mass remaining (%)") + xlab("Microsite") + theme_bw(base_size = 16) #+ facet_grid(.~Inst_GWC.x)

#is max CO2 related to soil moisture?
  #by default right now, max CO2 is last day of incubation and there are no differeces at instantaneous moisture
MIMout_diff$CO2.x <- MIMout_diff$CO2_MICr.x + MIMout_diff$CO2_MICK.x
MIMout_diff$CO2.y <- MIMout_diff$CO2_MICr.y + MIMout_diff$CO2_MICK.y
ggplot(MIMout_diff %>% filter(CO2.y > 0), aes(x=hist_GWC.y, y=CO2.y, color=SITE.y)) + geom_point(size = 3) + ylab("Max modeled CO2") + xlab("Historcal GWC (%)")

#CO2 vs soil moisture
MIMout_diff$CO2_def <- MIMout_diff$CO2_MICr.x + MIMout_diff$CO2_MICK.x
MIMout_diff$CO2_moi <- MIMout_diff$CO2_MICr.y + MIMout_diff$CO2_MICK.y
ggplot(MIMout_diff, aes(x=Inst_GWC.x, y=CO2_def, shape = "Default", color=hist_GWC.x)) + geom_point(size=4) +geom_point(size =4, aes(x=Inst_GWC.x, y=CO2_moi, shape = "Alternate", color=hist_GWC.x)) +
  ylab("Total CO2 emited at day 3600") + theme_bw()
ggplot(MIMout_diff, aes(x=hist_GWC.x, y=CO2_def, shape = "Default")) + geom_point(size=4) +geom_point(size =4, aes(x=hist_GWC.x, y=CO2_moi, shape = "Alternate")) +
  ylab("Total CO2 emited at day 3600") + theme_bw()
#######################################

#moisture function from testbed

################################################################################

#description of function is in CASA testbed manual found here: https://github.com/wwieder/biogeochem_testbed/tree/master

####
#this is the description from the manual
####

#the soil moisture multiplier on Vmax is given as f(theta) = (0.05, P(theta.l/theta.sat)^3 * (1-theta.l/theta.sat - theta.f/theta.sat)^gas_diff)
#where all fractions are calculated first and theta.l, theta.f, and theta.sat are the liquid, frozen, and saturated volumetric water content, respectively, gas_diff=2.5, and P=44.247 (a scalar that normalizes the function to a max value of 1)
#currently unclear what the "0.05," at the beginning means, maybe if you are above that this equation applies?
#also unclear how to determine theta.sat and theta.f for a given site - these are presumably global input values


###
#this is the code that matched the manual description above (in some non R language)
###

#R version of equation below and full code below that
#Read in soil moisture data as in CORPSE
#note: not sure what moistavg, frznmoistavg, and ssat represent;
       #I think min and max are checks to prevent values less than 0 and greater than 1!
theta_liq  = min(1.0, moistavg/soil$ssat)   # fraction of liquid water-filled pore space (0.0 - 1.0)
theta_frzn = min(1.0, frznmoistavg/ssat) # fraction of frozen water-filled pore space (0.0 - 1.0)
air_filled_porosity = max(0.0, 1.0-theta_liq-theta_frzn)
#CORPSE water scalar, adjusted to give maximum values of 1
fW = (theta_liq**3 * air_filled_porosity**2.5)/0.022600567942709
fW = max(0.05, fW)
#fW multiplies Vmax directly to down regulate decomposition

#
# ! Read in soil moisture data as in CORPSE
# theta_liq  = min(1.0, casamet%moistavg(npt)/soil%ssat(npt))     ! fraction of liquid water-filled pore space (0.0 - 1.0)
# theta_frzn = min(1.0, casamet%frznmoistavg(npt)/soil%ssat(npt)) ! fraction of frozen water-filled pore space (0.0 - 1.0)
# air_filled_porosity = max(0.0, 1.0-theta_liq-theta_frzn)
#
# if (mimicsbiome%fWFunction .eq. CORPSE) then
# ! CORPSE water scalar, adjusted to give maximum values of 1
# fW = (theta_liq**3 * air_filled_porosity**2.5)/0.022600567942709
# fW = max(0.05, fW)
# elseif (mimicsbiome%fWFunction .eq. CASACNP) then
# ! CASA water scalar, does not use frozen water in the calculation!
#   ! local variables
# fW = ((theta_liq-wfpscoefb)/(wfpscoefa-wfpscoefb))**wfpscoefe &
#   * ((theta_liq-wfpscoefc)/(wfpscoefa-wfpscoefc))**wfpscoefd
# fW = min(fW, 1.0)
# fW = max(0.01, fW)
# else
#   fW = 1.0
# endif
#
# mimicspool%fW(npt) =  fW
# mimicspool%thetaLiq(npt)  =  theta_liq
# mimicspool%thetaFrzn(npt) =  theta_frzn


####
#below here is the wrong function but keeping in case its useful!!
##

#this is a function in the testbed (seem to be different)


#soil moisture function in testbed isfrom Kelly et al., 2000 Fig 2b which uses a series of coefficients to create a relationship with soil water filled pore space
xkwater(npt) = (((fwps(npt)-wfpscoefb)/(wfpscoefa-wfpscoefb))^wfpscoefe)*(((fwps(npt)-wfpscoefc)/(wfpscoefa-wfpscoefc))^wfpscoefd)
#where:
 #wfpscoefa = 0.55 (optimal wfps), wfpscoefb = 1.70, wfpscoefc = -0.007, wfpscoefd = 3.22, wfpscoefe = 6.6481 (e is equivalent to wfpscoefd*(wfpscoefb-wfpscoefa)/(wfpscoefa-wfpscoefc))
 #npt is index to tell the code to do an equation for each grid cell (can be removed for MIMICS purposes)
#fwps is presumably the water filled pore space of the site - unclear how calculated - presumably nobody provides that outright...
   #equation: WFPS = (VWC/(1-(BD/PD)))*100 where VWC = volumetric water content (cm/cm), BD = bulk density (g/cm3), PD = particle density (2.65 g/cm3); eq from: https://doi.org/10.3390/app9030496
 #above coefficients technically for coarse textured soils

#actual function from testbed in below (in some other coding language)

# SUBROUTINE mimics_xratesoil(veg,soil,casamet,casabiome)
# !  to account for effects of T and W on litter decomposition: xklitter
# !  inputs:
#   !     ivt(mp)  :       biome type
# !     tsoilavg(mp):    soil temperature in K
# !     moistavg(mp):    volumetric soil moisture
#
# IMPLICIT NONE
# TYPE (veg_parameter_type),    INTENT(INOUT) :: veg  ! vegetation parameters
# TYPE (soil_parameter_type),   INTENT(INOUT) :: soil ! soil parameters
# TYPE (casa_met),              INTENT(INOUT) :: casamet
# TYPE (casa_biome),            INTENT(INOUT) :: casabiome
#
# ! local variables
# REAL(r_2), parameter :: wfpscoefa=0.55   ! Kelly et al. (2000) JGR, Figure 2b), optimal wfps
# REAL(r_2), parameter :: wfpscoefb=1.70   ! Kelly et al. (2000) JGR, Figure 2b)
# REAL(r_2), parameter :: wfpscoefc=-0.007 ! Kelly et al. (2000) JGR, Figure 2b)
# REAL(r_2), parameter :: wfpscoefd=3.22   ! Kelly et al. (2000) JGR, Figure 2b)
# REAL(r_2), parameter :: wfpscoefe=6.6481 ! =wfpscoefd*(wfpscoefb-wfpscoefa)/(wfpscoefa-wfpscoefc)
# ! Kirschbaum function parameters
# REAL(r_2), parameter :: xkalpha=-3.764   ! Kirschbaum (1995, SBB)
# REAL(r_2), parameter :: xkbeta=0.204
# REAL(r_2), parameter :: xktoptc=36.9
# REAL(r_2), DIMENSION(mp)       :: xkwater,xktemp,xklitter
# REAL(r_2), DIMENSION(mp)       :: fwps,tsavg
# INTEGER :: npt
#
# xklitter(:) = 1.0
# casaflux%klitter(:,:) = 0.0        ! initialize klitter for all 3 casa litter pool.  Only cwd will be reset.
# casaflux%fromLtoCO2(:,:) = 0.0     ! flow from L to CO2
#
# fwps(:)     =  min(1.0, casamet%moistavg(:)/soil%ssat(:))
# tsavg(:)    =  casamet%tsoilavg(:)
#
# DO npt=1,mp
# IF(casamet%iveg2(npt) /= icewater) THEN
# xktemp(npt)  = casabiome%q10soil(veg%iveg(npt))**(0.1*(tsavg(npt)-TKzeroC-35.0))
# xkwater(npt) = ((fwps(npt)-wfpscoefb)/(wfpscoefa-wfpscoefb))**wfpscoefe    &
#   * ((fwps(npt)-wfpscoefc)/(wfpscoefa-wfpscoefc))**wfpscoefd
# IF (veg%iveg(npt) == cropland .OR. veg%iveg(npt) == croplnd2) &
#   xkwater(npt)=1.0
# xklitter(npt) = casabiome%xkoptlitter(veg%iveg(npt)) * xktemp(npt) * xkwater(npt)
# casaflux%klitter(npt,cwd) = xklitter(npt) * casabiome%litterrate(veg%iveg(npt),cwd)
# !         write(*,*)
# !         write(*,*) 'mimics_xratesoil:'
# !         write(*,'(a42,f10.6)') 'casabiome%xkoptlitter(veg%iveg(npt)) =', casabiome%xkoptlitter(veg%iveg(npt))
# !         write(*,'(a42,f10.6)') 'xktemp(npt) =', xktemp(npt)
# !         write(*,'(a42,f10.6)') 'xkwater(npt) =', xkwater(npt)
# !         write(*,'(a42,f10.6)') 'casabiome%litterrate(veg%iveg(npt),cwd) =', casabiome%litterrate(veg%iveg(npt),cwd)
# !         write(*,'(a42,f10.6)') 'casaflux%klitter(npt,cwd) =', casaflux%klitter(npt,cwd)
#
#
# !! from casa_coeffsoil for reference
# !! casaflux%fromLtoS(:,mic,cwd)   = 0.40*(1.0 - casabiome%fracLigninplant(veg%iveg(:),wood)) ! CWD -> fmic
# !! casaflux%fromLtoS(:,slow,cwd)  = 0.7 * casabiome%fracLigninplant(veg%iveg(:),wood)        ! CWD -> slow
#
# ! Fraction of cwd decomposition that goes to heterotrophic respiration
# casaflux%fromLtoCO2(npt,cwd) = 1.0 - 0.40*(1.0 - casabiome%fracLigninplant(veg%iveg(npt),wood)) &
#   - 0.7 * casabiome%fracLigninplant(veg%iveg(npt),wood)
# casaflux%fromLtoCO2(npt,cwd) = MAX(0.0, casaflux%fromLtoCO2(npt,cwd))
# casaflux%fromLtoCO2(npt,cwd) = MIN(1.0, casaflux%fromLtoCO2(npt,cwd))
#
# END IF
# END DO
#
# END SUBROUTINE mimics_xratesoil

#################################################
