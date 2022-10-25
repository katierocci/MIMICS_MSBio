library(tidyverse)
library(ggplot2)
library(Metrics)

df <- mcdf %>% filter (DAY==200)
df$CO2_of_tot <- rowSums(df[,11:12])/rowSums(df[,4:12])
df$LIT_tot <- df$LITm + df$LITs


dfA15 <- df %>% filter(SITE == "A") %>% 
                filter(MAT == 25) %>% 
                filter(CLAY==10) %>%
                filter(LIG == 10)

pairs(dfA15[,c(25, 11,12, 10, 9, 21:24)],
      col = alpha("red", 0.05), # Change color
      pch = 16, # Change shape of points
      cex=1.5,
      labels = c("CO2_frac_of_tot", "CO2-r", "CO2-K", "SOMa", "SOMc", "Vslope", "Vint", "Kslope", "Kint"), # Change labels of diagonal
      main = "n=1000 Random Parameter Sets",
      upper.panel = NULL,
      breaks=c(0,0.2))


dfA15_lit_limit <- df %>% filter(SITE == "A") %>% 
  filter(MAT == 25) %>% 
  filter(CLAY==10) %>%
  filter(LIG == 10) %>%
  filter(LIT_tot > 50) %>%
  filter(LIT_tot < 90)

pairs(dfA15_lit_limit[,c(25, 11,12, 26,6,7,10,9, 21:24)],
      col = alpha("black", 0.2), # Change color
      pch = 16, # Change shape of points
      cex=1.5,
      labels = c("CO2_frac_of_tot", "CO2-r", "CO2-K", "LIT", "MICr", "MICK", "SOMa", "SOMc", "Vslope", "Vint", "Kslope", "Kint"), # Change labels of diagonal
      main = "n=1000 Random Parameter Set \nLimited to runs that consumed 10-40% of litter",
      upper.panel = NULL,
      breaks=c(0,0.2))
