library(tidyverse)
library(ggplot2)
library(ggridges)
library(ggpubr)

df <- readRDS("C:/github/MIMICS_MSBio/Cheyenne_HPC/HPC_output/MSBio_MIM_MC_runs-1e+05_20221028_112647_.rds")
df$CO2_of_tot <- rowSums(df[,11:12])/rowSums(df[,4:12])
df$LIT_tot <- df$LITm + df$LITs
df$cost <- abs(df$CO2_of_tot - df$CO2C_prop)
df$MIM_CO <- df$MICr/df$MICK


### Ftn to find best psets

best_psets <- function(df, best_n) {
  cost_df <- df %>% group_by(run_num) %>% summarize(SITE = unique(SITE), n = n(), cost = sum(cost))
  return(cost_df %>% slice_min(cost, n=best_n))
}

best_data <- function(df, best_n) {
  cost_df <- df %>% group_by(run_num) %>% summarize(SITE = unique(SITE), 
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
BART <- best_psets(df = df %>% filter(SITE == "BART"), best_n = 1000)
BART_data <- best_data(df = df %>% filter(SITE == "BART"), best_n = 1000)

plot_df <- df %>% filter(SITE == "BART") %>%
                  filter(MAT == 15) %>%
                  filter(run_num %in% BART$run_num)

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


### Collect all the best psets by site
sites <- unique(df$SITE)
top_df <- NULL
for(i in 1:length(sites)) {
  df1 <- best_data(df = df %>% filter(SITE == sites[i]), best_n = 100)
  top_df <- rbind(top_df, df1)
}

### Ridge plots
pCost <- ggplot(top_df, aes(x = (cost/CO2C_prop)*100, y = SITE, fill=SITE)) +
  geom_density_ridges(scale = 2) + 
  scale_y_discrete(expand = c(0, 0)) +     # will generally have to set the `expand` option
  scale_x_continuous(expand = c(0, 0)) +   # for both axes to remove unneeded padding
  coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
  theme_ridges() +
  scale_fill_brewer(palette = "Greys") +
  labs(title="Cumulative %diff from incubation CO2-C target",
       subtitle="n = 100 lowest cost parameter sets") +
  theme(legend.position = "none") #removes the legend
pCost

pVint <- ggplot(top_df, aes(x = Vint_x, y = SITE, fill=SITE)) +
  geom_density_ridges(scale = 2) + 
  scale_y_discrete(expand = c(0, 0)) +     # will generally have to set the `expand` option
  scale_x_continuous(expand = c(0, 0)) +   # for both axes to remove unneeded padding
  coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
  theme_ridges() +
  scale_fill_brewer(palette = "Blues") +
  labs(title="Vint Multiplier") +#,
          #subtitle="n=200 lowest cost parameter sets") +
  theme(legend.position = "none") #removes the legend
pVint

pVslope <- ggplot(top_df, aes(x = Vslope_x, y = SITE, fill=SITE)) +
  geom_density_ridges(scale = 2) + 
  scale_y_discrete(expand = c(0, 0)) +     # will generally have to set the `expand` option
  scale_x_continuous(expand = c(0, 0)) +   # for both axes to remove unneeded padding
  coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
  theme_ridges() +
  scale_fill_brewer(palette = "Blues") +
  labs(title="Vslope Multiplier") +#,
       #subtitle="n=200 lowest cost parameter sets") +
  theme(legend.position = "none") #removes the legend
pVslope

pKslope <- ggplot(top_df, aes(x = Kslope_x, y = SITE, fill=SITE)) +
  geom_density_ridges(scale = 2) + 
  scale_y_discrete(expand = c(0, 0)) +     # will generally have to set the `expand` option
  scale_x_continuous(expand = c(0, 0)) +   # for both axes to remove unneeded padding
  coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
  theme_ridges() +
  scale_fill_brewer(palette = "Oranges") +
  labs(title="Kslope Multiplier") +#,
       #subtitle="n=200 lowest cost parameter sets") +
  theme(legend.position = "none") #removes the legend
pKslope

pKint <- ggplot(top_df, aes(x = Kint_x, y = SITE, fill=SITE)) +
  geom_density_ridges(scale = 2) + 
  scale_y_discrete(expand = c(0, 0)) +     # will generally have to set the `expand` option
  scale_x_continuous(expand = c(0, 0)) +   # for both axes to remove unneeded padding
  coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
  theme_ridges() +
  scale_fill_brewer(palette = "Oranges") +
  labs(title="Kint Multiplier") +#,
       #subtitle="n=200 lowest cost parameter sets") +
  theme(legend.position = "none") #removes the legend
pKint

pVMAX <- ggplot(top_df %>% filter(VMAX < 2e-04), aes(x = VMAX, y = SITE, fill=SITE)) +
  geom_density_ridges(scale = 2) + 
  scale_y_discrete(expand = c(0, 0)) +     # will generally have to set the `expand` option
  scale_x_continuous(expand = c(0, 0)) +   # for both axes to remove unneeded padding
  coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
  theme_ridges() +
  scale_fill_brewer(palette = "Greens") +
  labs(title="VMAX[1]",
       subtitle="Clipped: VMAX < 2e-04") +
  theme(legend.position = "none") #removes the legend
pVMAX

pKM <- ggplot(top_df %>% filter(KM < 2), aes(x = KM, y = SITE, fill=SITE)) +
  geom_density_ridges(scale = 2) + 
  scale_y_discrete(expand = c(0, 0)) +     # will generally have to set the `expand` option
  scale_x_continuous(expand = c(0, 0)) +   # for both axes to remove unneeded padding
  coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
  theme_ridges() +
  scale_fill_brewer(palette = "Greens") +
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
