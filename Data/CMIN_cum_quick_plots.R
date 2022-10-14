library(dplyr)
library(ggplot2)

setwd("C:/GitHub/MIMICS_MSBio/Data")
df <- read.csv("CMIN_cum.csv", as.is=T)

df2 <- df %>% filter(site == "BART")

ggplot(df, aes(x=as.character(temp.trt), y=cumulativeCO2Flux, fill = as.factor(moisture.trt))) +
  geom_boxplot() + 
  theme_bw() +
  ylim(0,6) +
  xlab("Temperature Group (deg C)") +
  ylab("Cumulative CO2 Flux (mg C/ g soil)") +
  geom_vline(xintercept = 1.5) +
  facet_wrap(~ site)


ggplot(df, aes(x=as.factor(site), y=cumulativeCO2Flux, fill = as.factor(moisture.trt))) +
  geom_boxplot() + 
  theme_bw() +
  ylim(0,6) +
  xlab("Temperature Group (deg C)") +
  ylab("Cumulative CO2 Flux (mg C/ g soil)") +
  facet_wrap(~ as.factor(temp.trt))

ggplot(df %>% filter(moisture.trt == 60), aes(x=as.character(temp.trt), y=cumulativeCO2Flux, fill = as.factor(moisture.trt))) +
  geom_boxplot() + 
  theme_bw() +
  ylim(0,6) +
  xlab("Temperature Group (deg C)") +
  ylab("Cumulative CO2 Flux (mg C/ g soil)") +
  geom_vline(xintercept = 1.5) +
  facet_wrap(~ site)


