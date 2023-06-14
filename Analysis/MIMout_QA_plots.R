library(tidyverse)
library(ggplot2)
library(ggridges)
library(ggpubr)
library(gridExtra)

#KT edited this to make it work with microsite data! Commented out lines are from Derek
#first tried with 100 runs and only top 10 values but now trying with 1000 runs and top 100 values (density plots not very helpful rn)

#df <- readRDS("C:/github/MIMICS_MSBio/Cheyenne_HPC/HPC_output/MSBio_MIM_MC_runs-1e+05_20221028_112647_.rds")
#MC_MIMICS <- readRDS("C:/Users/katie/OneDrive - Colostate/Documents/Postdoc/CU/Macrosystems project/Project Code/MIMICS_MSBio_KR/temp/MSBio_MC_10000_20230524_123409_.rds")
df <- left_join(MC_MIMICS, data, by = "ID") #MC output has to be first for CO2 rows below to be right
df <- left_join(df, rand_params, by = "run_num")
df$CO2_of_tot <- rowSums(df[,11:12])/rowSums(df[,4:12])
df$LIT_tot <- df$LITm + df$LITs
df$MIM_CO <- df$MICr/df$MICK
#df$cost <- abs(df$CO2_of_tot - df$CO2C_prop)
#alt cost function - trying to match max CO2 prop across WHCs
#Will suggests the alt cost function!
Max_CO2 <- data %>% group_by(SITE) %>% summarise(max=max(CO2C_prop))
Max_CO2$GWC <- c(0.623, 0.674, 0.658, 0.47, 0.621, 0.491, 0.544, 0.466, 0.591, 0.656, 0.51, 0.376)
ggplot(Max_CO2, aes(x=GWC, y=max)) +geom_point() + geom_smooth() #vaguely positive
#cost function
df_cost <- df %>% filter(DAY == 105) %>% group_by(SITE.x, run_num) %>% summarise(max_model= max(CO2_of_tot))
df_cost$max_incub <- c(rep(0.0798,1000), rep(0.0846, 1000), rep(0.0645, 1000), rep(0.0698, 1000),
                       rep(0.0699, 1000), rep(0.0866, 1000), rep(0.0774, 1000), rep(0.0676, 1000),
                       rep(0.0742, 1000), rep(0.0790, 1000), rep(0.0750, 1000), rep(0.0724, 1000))
df_cost <-left_join(df_cost, rand_params, by = "run_num")
df_MIM <- df %>% filter(DAY == 105) %>% group_by(SITE.x, run_num) %>% summarise(MIM_CO= mean(MIM_CO))
df_cost$Site_RN <- paste(df_cost$SITE.x, df_cost$run_num, sep = "_")
df_MIM$Site_RN <- paste(df_MIM$SITE.x, df_MIM$run_num, sep = "_")
df_cost <-left_join(df_cost, df_MIM, by = "Site_RN")
df_cost$cost <- abs(df_cost$max_model - df_cost$max_incub)
#debugging
ggplot(df_cost, aes(x=SITE.x.x, y=cost, group=SITE.x.x)) + geom_boxplot()
#cost is different between sites but distributed the same
ggplot(df_cost, aes(x=CUE_x, y=cost, group=SITE.x.x, color=SITE.x.x)) + geom_point()
#CUE so tightly tied to CO2 that the lowest cost is always the same parameter
#tying vMOD below to see if also strongly constrained like CUE
ggplot(df_cost, aes(x=vMOD_x, y=cost, group=SITE.x.x, color=SITE.x.x)) + geom_point()
#tau also looks the same...
ggplot(df_cost, aes(x=Tau_x, y=cost, group=SITE.x.x, color=SITE.x.x)) + geom_point()

### Ftn to find best psets (this is an objective function?)

best_psets <- function(df, best_n) {
  cost_df <- df %>% group_by(run_num) %>% summarize(SITE = unique(SITE.x.x), n = n(), cost = sum(cost))
  return(cost_df %>% slice_min(cost, n=best_n))
}

best_data <- function(df, best_n) {
  cost_df <- df %>% group_by(run_num) %>% summarize(SITE = unique(SITE.x.x),
                                                    n = n(),
                                                    cost = sum(cost),
                                                    min_MIM_CO = min(MIM_CO),
                                                    max_MIM_CO = max(MIM_CO))
  low_cost <- cost_df %>%
                filter(min_MIM_CO > 0.01) %>%
                filter(max_MIM_CO < 100) %>%
                slice_min(cost, n=best_n)

  out <- df %>% filter(run_num %in% low_cost$run_num)
  return(out)
}


### Plot best psets
df$SITE<- df$SITE.x
BART.2 <- best_psets(df = df %>% filter(SITE == "2"), best_n = 100) #originally BART and 1000 for this and below
BART.2_data <- best_data(df = df %>% filter(SITE == "2"), best_n = 100)

plot_df <- df %>% filter(SITE == "2") %>%
                  #filter(MAT == 15) %>%
                  filter(run_num %in% BART.2$run_num)

#not sure how Derek's data table is set up so not sure this code works for my purposes
pairs(plot_df[,c(30,28,29, 11,12,6,7,10,9,24:27)],
      col = alpha("black", 0.2), # Change color
      pch = 16, # Change shape of points
      cex=1.5,
      labels = c("cost", "CO2_frac_of_tot", "LIT", "CO2-r", "CO2-K", "MICr", "MICK", "SOMa", "SOMc", "Vslope", "Vint", "Kslope", "Kint"), # Change labels of diagonal
      main = "n=30 Lowest Cost Parameter Sets",
      upper.panel = NULL,
      breaks=c(0,0.2))


### 3D plot example
library(plotly)
plot_ly(plot_df, x = ~Vint_x, y = ~Vslope_x, z = ~Kslope_x, color = ~cost) #, colors = c('#BF382A', '#0C4B8E'))
plot_ly(plot_df, x = ~Vint_x, y = ~Kint_x, z = ~-cost, color = ~cost)


####
#this is most helpful code
####

### Collect all the best psets by site
df_cost$run_num<- df_cost$run_num.x
sites <- unique(df_cost$SITE.x.x)
top_df <- NULL
for(i in 1:length(sites)) {
  df1 <- best_data(df = df_cost %>% filter(SITE.x.x == sites[i]), best_n = 100)
  top_df <- rbind(top_df, df1)
}

### Ridge plots - had to add group to aesthetics and comment out fill to get to work
      #likely because site is number, if make factor would likely work
pCost <- ggplot(top_df, aes(x = (cost/CO2C_prop)*100, y = SITE.x.x, fill=SITE.x.x, group=SITE.x.x)) +
  geom_density_ridges(scale = 2) +
  scale_y_discrete(expand = c(0, 0)) +     # will generally have to set the `expand` option
  scale_x_continuous(expand = c(0, 0)) +   # for both axes to remove unneeded padding
  coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
  theme_ridges() +
  #scale_fill_brewer(palette = "Greys") +
  labs(title="Cumulative %diff from incubation CO2-C target",
       subtitle="n = 100 lowest cost parameter sets") +
  theme(legend.position = "none") #removes the legend
pCost

pVint <- ggplot(top_df, aes(x = Vint_x, y = SITE.x.x, fill=SITE.x.x, group=SITE.x.x)) +
  geom_density_ridges(scale = 2) +
  scale_y_discrete(expand = c(0, 0)) +     # will generally have to set the `expand` option
  scale_x_continuous(expand = c(0, 0)) +   # for both axes to remove unneeded padding
  coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
  theme_ridges() +
  #scale_fill_brewer(palette = "Blues") +
  labs(title="Vint Multiplier") +#,
          #subtitle="n=200 lowest cost parameter sets") +
  theme(legend.position = "none") #removes the legend
pVint

pVslope <- ggplot(top_df, aes(x = Vslope_x, y = SITE, fill=SITE, group=SITE)) +
  geom_density_ridges(scale = 2) +
  scale_y_discrete(expand = c(0, 0)) +     # will generally have to set the `expand` option
  scale_x_continuous(expand = c(0, 0)) +   # for both axes to remove unneeded padding
  coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
  theme_ridges() +
  #scale_fill_brewer(palette = "Blues") +
  labs(title="Vslope Multiplier") +#,
       #subtitle="n=200 lowest cost parameter sets") +
  theme(legend.position = "none") #removes the legend
pVslope

pKslope <- ggplot(top_df, aes(x = Kslope_x, y = SITE, fill=SITE, group=SITE)) +
  geom_density_ridges(scale = 2) +
  scale_y_discrete(expand = c(0, 0)) +     # will generally have to set the `expand` option
  scale_x_continuous(expand = c(0, 0)) +   # for both axes to remove unneeded padding
  coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
  theme_ridges() +
  #scale_fill_brewer(palette = "Oranges") +
  labs(title="Kslope Multiplier") +#,
       #subtitle="n=200 lowest cost parameter sets") +
  theme(legend.position = "none") #removes the legend
pKslope

pKint <- ggplot(top_df, aes(x = Kint_x, y = SITE.x.x, fill=SITE.x.x, group=SITE.x.x)) +
  geom_density_ridges(scale = 2) +
  scale_y_discrete(expand = c(0, 0)) +     # will generally have to set the `expand` option
  scale_x_continuous(expand = c(0, 0)) +   # for both axes to remove unneeded padding
  coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
  theme_ridges() +
  #scale_fill_brewer(palette = "Oranges") +
  labs(title="Kint Multiplier") +#,
       #subtitle="n=200 lowest cost parameter sets") +
  theme(legend.position = "none") #removes the legend
pKint


pVmod <- ggplot(top_df, aes(x = vMOD_x, y = SITE.x.x, fill=SITE.x.x, group=SITE.x.x)) +
  geom_density_ridges(scale = 2) +
  scale_y_discrete(expand = c(0, 0)) +     # will generally have to set the `expand` option
  scale_x_continuous(expand = c(0, 0)) +   # for both axes to remove unneeded padding
  coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
  theme_ridges() +
  #scale_fill_brewer(palette = "Oranges") +
  labs(title="vMOD Multiplier") +#,
  #subtitle="n=200 lowest cost parameter sets") +
  theme(legend.position = "none") #removes the legend
pVmod


pKmod <- ggplot(top_df, aes(x = kMOD_x, y = SITE.x.x, fill=SITE.x.x, group=SITE.x.x)) +
  geom_density_ridges(scale = 2) +
  scale_y_discrete(expand = c(0, 0)) +     # will generally have to set the `expand` option
  scale_x_continuous(expand = c(0, 0)) +   # for both axes to remove unneeded padding
  coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
  theme_ridges() +
  #scale_fill_brewer(palette = "Oranges") +
  labs(title="kMOD Multiplier") +#,
  #subtitle="n=200 lowest cost parameter sets") +
  theme(legend.position = "none") #removes the legend
pKmod

pCUE <- ggplot(top_df, aes(x = CUE_x, y = SITE.x.x, fill=SITE.x.x, group=SITE.x.x)) +
  geom_density_ridges(scale = 2) +
  scale_y_discrete(expand = c(0, 0)) +     # will generally have to set the `expand` option
  scale_x_continuous(expand = c(0, 0)) +   # for both axes to remove unneeded padding
  coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
  theme_ridges() +
  #scale_fill_brewer(palette = "Oranges") +
  labs(title="CUE Multiplier") +#,
  #subtitle="n=200 lowest cost parameter sets") +
  theme(legend.position = "none") #removes the legend
pCUE
#all the same for all 1000 parameter multipliers as well so not a top_df issue

pTAU <- ggplot(top_df, aes(x = Tau_x, y = SITE.x.x, fill=SITE.x.x, group=SITE.x.x)) +
  geom_density_ridges(scale = 2) +
  scale_y_discrete(expand = c(0, 0)) +     # will generally have to set the `expand` option
  scale_x_continuous(expand = c(0, 0)) +   # for both axes to remove unneeded padding
  coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
  theme_ridges() +
  #scale_fill_brewer(palette = "Oranges") +
  labs(title="Tau Multiplier") +#,
  #subtitle="n=200 lowest cost parameter sets") +
  theme(legend.position = "none") #removes the legend
pTAU

#don't have VMAX automatically in my df - would need to calculate for each multiplier
pVMAX <- ggplot(top_df %>% filter(VMAX < 2e-04), aes(x = VMAX, y = SITE, fill=SITE, group=SITE)) +
  geom_density_ridges(scale = 2) +
  scale_y_discrete(expand = c(0, 0)) +     # will generally have to set the `expand` option
  scale_x_continuous(expand = c(0, 0)) +   # for both axes to remove unneeded padding
  coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
  theme_ridges() +
  #scale_fill_brewer(palette = "Greens") +
  labs(title="VMAX[1]",
       subtitle="Clipped: VMAX < 2e-04") +
  theme(legend.position = "none") #removes the legend
pVMAX

#don't have VMAX automatically in my df - would need to calculate for each multiplier
pKM <- ggplot(top_df %>% filter(KM < 2), aes(x = KM, y = SITE, fill=SITE, group=SITE)) +
  geom_density_ridges(scale = 2) +
  scale_y_discrete(expand = c(0, 0)) +     # will generally have to set the `expand` option
  scale_x_continuous(expand = c(0, 0)) +   # for both axes to remove unneeded padding
  coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
  theme_ridges() +
  #scale_fill_brewer(palette = "Greens") +
  labs(title="KM[1]",
       subtitle="Clipped: KM < 2") +
  theme(legend.position = "none") #removes the legend
pKM

top_row = ggarrange(pCost, ncol = 1, labels = c("A"))
bottom_row = ggarrange(
          pVint, pVslope,
          pKint, pKslope,
          pVMAX, pKM,
          nrow=3,
          ncol=2,
          labels = c("B", "C", "D", "E"))

png(file = "C:/github/MIMICS_MSBio/Cheyenne_HPC/HPC_output/best_pset_ridge_plot.png", width = 1000, height = 1200, res=108)
ggarrange(top_row, bottom_row,
          nrow=2,
          ncol=1,
          heights = c(1,3))
dev.off()

######
#Katie's code down here - looking for relationship between multipliers and soil moisture
##########
#summarize by mean and median
top_df$SITE <- as.character(top_df$SITE.x.x)
data$SITE <- as.character(data$SITE)
#Vint_sum <- top_df %>% group_by(SITE) %>%
#  summarise(
#                     mean.Vi=mean(Vint_x),
#                     med.Vi = median(Vint_x))
#Vslope_sum <- top_df %>% group_by(SITE) %>%
#  summarise(
#    mean.Vs=mean(Vslope_x),
#    med.Vs = median(Vslope_x))
#Kint_sum <- top_df %>% group_by(SITE) %>%
#  summarise(
#    mean.Ki=mean(Kint_x),
#    med.Ki = median(Kint_x))
#Kslope_sum <- top_df %>% group_by(SITE) %>%
#  summarise(
#    mean.Ks=mean(Kslope_x),
#    med.Ks = median(Kslope_x))
# vMOD_sum <- top_df %>% group_by(SITE) %>%
#   summarise(
#     mean.vMOD=mean(vMOD_x),
#     med.vMOD = median(vMOD_x))
# kMOD_sum <- top_df %>% group_by(SITE) %>%
#   summarise(
#     mean.kMOD=mean(kMOD_x),
#     med.kMOD = median(kMOD_x))
##all the same for CUE and Tau
CUE_sum <- top_df %>% group_by(SITE) %>%
     summarise(
       mean.CUE=mean(CUE_x),
       med.CUE = median(CUE_x))
Tau_sum <- top_df %>% group_by(SITE) %>%
  summarise(
    mean.tau=mean(Tau_x),
    med.tau = median(Tau_x))
MC_sum <- vMOD_sum %>% #left_join(Vslope_sum, by = "SITE") %>%
  #left_join(Kint_sum, by = "SITE") %>%
  #left_join(Kslope_sum, by = "SITE") %>%
  #left_join(data, by = "SITE") %>%
  left_join(kMOD_sum, by = "SITE") %>%
  left_join(data, by = "SITE") %>%
  pivot_longer(cols=c(2,4), values_to = "mean", names_to = "Parameter") %>%
  pivot_longer(cols=c(2,3), values_to = "median", names_to = "Parameter2")

ggplot(MC_sum, aes(x = GWC_site, y=mean, group=Parameter, color=Parameter)) + geom_point() +geom_smooth(method = "lm")
ggplot(MC_sum, aes(x = GWC_site, y=median, group=Parameter2, color=Parameter2)) + geom_point() +geom_smooth(method = "lm")
MC_sum_Vs.mean <- filter(MC_sum, Parameter == "mean.Vs")
MC_sum_Vs.med <- filter(MC_sum, Parameter2 == "med.Vs")
Vs.mean_GWC <- lm(mean ~ GWC_site, data=MC_sum_Vs.mean)
summary(Vs.mean_GWC)
Vs.med_GWC <- lm(median ~ GWC_site, data=MC_sum_Vs.med)
summary(Vs.med_GWC)
MC_sum_Ks.mean <- filter(MC_sum, Parameter == "mean.Ks")
MC_sum_Ks.med <- filter(MC_sum, Parameter2 == "med.Ks")
Ks.mean_GWC <- lm(GWC ~ mean, data=MC_sum_Ks.mean)
summary(Ks.mean_GWC)
Ks.med_GWC <- lm(GWC ~ median, data=MC_sum_Ks.med)
summary(Ks.med_GWC)
MC_sum_Vi.mean <- filter(MC_sum, Parameter == "mean.Vi")
MC_sum_Vi.med <- filter(MC_sum, Parameter2 == "med.Vi")
Vi.mean_GWC <- lm(mean ~ GWC_site, data=MC_sum_Vi.mean)
summary(Vi.mean_GWC)
Vi.med_GWC <- lm(median ~ GWC_site, data=MC_sum_Vi.med)
summary(Vi.med_GWC)
MC_sum_Ki.mean <- filter(MC_sum, Parameter == "mean.Ki")
MC_sum_Ki.med <- filter(MC_sum, Parameter2 == "med.Ki")
Ki.mean_GWC <- lm(mean ~ GWC_site, data=MC_sum_Ki.mean)
summary(Ki.mean_GWC)
Ki.med_GWC <- lm(median ~ GWC_site, data=MC_sum_Ki.med)
summary(Ki.med_GWC)
MC_sum_vMOD.mean <- filter(MC_sum, Parameter == "mean.vMOD")
MC_sum_vMOD.med <- filter(MC_sum, Parameter2 == "med.vMOD")
vMOD.mean_GWC <- lm(mean ~ GWC_site, data=MC_sum_vMOD.mean)
summary(vMOD.mean_GWC)
vMOD.med_GWC <- lm(median ~ GWC_site, data=MC_sum_vMOD.med)
summary(vMOD.med_GWC)
MC_sum_kMOD.mean <- filter(MC_sum, Parameter == "mean.kMOD")
MC_sum_kMOD.med <- filter(MC_sum, Parameter2 == "med.kMOD")
kMOD.mean_GWC <- lm(mean ~ GWC_site, data=MC_sum_kMOD.mean)
summary(kMOD.mean_GWC)
kMOD.med_GWC <- lm(median ~ GWC_site, data=MC_sum_kMOD.med)
summary(kMOD.med_GWC)
#plotting each separately
mean.Vs_plot <- MC_sum %>% filter(Parameter == "mean.Vs") %>% ggplot(aes(x = GWC_site, y=mean, group=Parameter)) +
  geom_point(aes(color = Parameter)) + geom_smooth(method = "lm", aes(color=Parameter))
mean.Vi_plot <- MC_sum %>% filter(Parameter == "mean.Vi") %>% ggplot(aes(x = GWC_site, y=mean, group=Parameter)) +
  geom_point(aes(color = Parameter)) + geom_smooth(method = "lm", aes(color=Parameter))
mean.Ks_plot <- MC_sum %>% filter(Parameter == "mean.Ks") %>% ggplot(aes(x = GWC_site, y=mean, group=Parameter)) +
  geom_point(aes(color = Parameter)) + geom_smooth(method = "lm", aes(color=Parameter))
mean.Ki_plot <- MC_sum %>% filter(Parameter == "mean.Ki") %>% ggplot(aes(x = GWC_site, y=mean, group=Parameter)) +
  geom_point(aes(color = Parameter)) + geom_smooth(method = "lm", aes(color=Parameter))
grid.arrange(mean.Vs_plot, mean.Vi_plot, mean.Ks_plot, mean.Ki_plot, ncol = 2, nrow = 2)
