library(dplyr)
library(ggplot2)
library(forcats)

setwd("C:/github/MIMICS_MSBio/Analysis/Incubation")

cmin <- read.csv("CMIN_cum.csv", as.is=T)

cmin <- cmin %>%
  mutate(site = fct_relevel(site, levels = "TREE", "BART", "GRSM", "SERC", "TALL", "LENO"))

# Site x MAT
#------------------------------------------------------------------------------------
png("site_MAT.png", width = 1600, height = 1200, res = 300)
ggplot(cmin, aes(x=site, y=MAT_c, color=MAT_c)) + geom_point(size=6) +
  scale_color_gradient2(midpoint=mean(cmin$MAT_c), low="blue", mid="yellow",
                        high="red" ) +
  ggtitle("NEON Study Sites") +
  xlab("Site") +
  ylab("Mean Annual Temperature (°C)") +
  labs(color="MAT (deg C)") +
  ylim(0,20) +
  theme_minimal() #+ theme(axis.text = element_text(size = 16))  
dev.off()

png("site_MAT_MAP.png", width = 1600, height = 1200, res = 300)
ggplot(cmin, aes(x=MAT_c, y=MAP_mm, color=MAT_c)) + geom_point(size=6) +
  scale_color_gradient2(midpoint=mean(cmin$MAT_c), low="blue", mid="yellow",
                        high="red" ) +
  ggtitle("NEON Study Sites") +
  ylab("Mean annual precipitation (mm)") +
  xlab("Mean annual temperature (°C)") +
  labs(color="MAT (deg C)") +
  xlim(4,20) +
  ylim(600, 1600) +
  theme_minimal() #+ theme(axis.text = element_text(size = 16))  
dev.off()


# CO2 Flux x Site MAT @ 15
#------------------------------------------------------------------------------------
png("Cflux_15.png", width = 1800, height = 1200, res = 300)
ggplot(cmin %>% filter(temp.trt == 15), aes(y=cumulativeCO2Flux, x=moisture.trt, group=unique.id, color=MAT_c)) +
  geom_line(size=0.5, alpha=0.8) +
  scale_color_gradient2(midpoint=mean(cmin$MAT_c), low="blue", mid="yellow",
                        high="red" ) +
  ggtitle("Soil Incubation @ 15 °C") +
  ylab(paste0("CO2 Flux (mg C/g soil)")) +
  xlab("Soil moisture (%)") +
  labs(color="MAT (deg C)") +
  ylim(0,75) +
  theme_minimal()
dev.off()

png("Cflux_25.png", width = 1800, height = 1200, res = 300)
ggplot(cmin %>% filter(temp.trt == 25), aes(y=cumulativeCO2Flux, x=moisture.trt, group=unique.id, color=MAT_c)) +
  geom_line(size=0.5, alpha=0.8) +
  scale_color_gradient2(midpoint=mean(cmin$MAT_c), low="blue", mid="yellow",
                        high="red" ) +
  ggtitle("Soil Incubation @ 25 °C") +
  ylab(paste0("CO2 Flux (mg C/g soil)")) +
  xlab("Soil moisture (%)") +
  labs(color="MAT (deg C)") +
  ylim(0,75) +
  theme_minimal()
dev.off()




# CO2 Flux x Site MAT @ 15
#------------------------------------------------------------------------------------
cmin_box <- cmin %>%
  mutate(site = fct_relevel(site, levels = "TREE", "SERC", "TALL", "BART", "GRSM", "LENO"))

png("Cflux_15_box.png", width = 1800, height = 1200, res = 300)
ggplot(cmin_box %>% filter(temp.trt == 15), aes(x=moisture.trt, y=cumulativeCO2Flux, group=moisture.trt, fill = factor(site), alpha = moisture.trt)) +
  scale_fill_manual(values=c("#0000cc", "#ccb909",  "#cb6208", "#6c34b3", "#ccb007", "#cc4a05"), guide = 'none') + 
  geom_boxplot() + 
  ggtitle("Soil Incubation @ 15 °C") +
  ylab(paste0("CO2 Flux (mg C/g soil)")) +
  xlab("Soil moisture (%)") +
  scale_alpha(guide = 'none') +
  theme_bw() +
  facet_wrap(~site)  
dev.off()
