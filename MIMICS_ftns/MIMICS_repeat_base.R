###########################################
# MIMICS repeat run function
###########################################

source("MIMICS_ftns/MIMICS_INC_daily.R")

MIMrepeat <- function(forcing_df, rparams) {

  #debug
  #forcing_df = data
  #rparams = rand_params

  # Set global model parameters
  #Vslope <<- Vslope_default * rparams$Vslope_x[1]
  #Vint <<- Vint_default * rparams$Vint_x[1]
  #Kslope <<- Kslope_default * rparams$Kslope_x[1]
  #Kint <<- Kint_default * rparams$Kint_x[1]

  # Tau_MULT <<- Tau_MULT_default * rparams$Tau_x[1]
   CUE <<- CUE_default * rparams$CUE_x[1]
  # desorb_MULT <<- desorb_MULT_default * rparams$desorb_x[1]
  # fPHYS_MULT <<- fPHYS_MULT_default * rparams$fPHYS_x[1]
  # avGWC <<- rparams$avGWC[1]
  # akGWC <<- rparams$akGWC[1]

  #vMOD <<- vMOD_default * rparams$vMOD_x[1]
  #kMOD <<- kMOD_default * rparams$kMOD_x[1]

  #full run of forcing data csv
  MIMrun <- forcing_df %>% split(1:nrow(forcing_df)) %>% map(MIMICS_INC_DAILY) %>% bind_rows()

  #Optional combine MIMout with forcing data
  #MIMrun <- MIMrun %>% left_join(forcing_df %>% select(-SITE), by="ID")

  #add run number
  MIMrun$run_num <- rparams$run_num[1]

  # return MIMrun
  return(MIMrun)
}


# DIAGNOSTIC CODE BITS
#output <- MIMrepeat(forcing_df = data, rparams = rand_params) %>% filter(DAY==200)
