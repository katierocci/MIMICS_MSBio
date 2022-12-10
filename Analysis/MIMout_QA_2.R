suppressMessages(library(dplyr))
library(ggplot2)
library(ggridges)
library(ggpubr)
library(DT)
suppressMessages(library(plotly))
library(forcats)

setwd("C:/GitHub/MIMICS_MSBio")

### LOAD MC RUN DATA
#For a single file
data_dir <- 'C:/github/MIMICS_MSBio/Cheyenne_HPC/HPC_output/'
run_files <- list.files(path=data_dir, pattern="20221201", all.files=TRUE, full.names=TRUE)
MIMout <- NULL
for(i in 1:length(run_files)) {
  run_rds <- readRDS(run_files[i])
  run_num_add <- 1200000*(i-1)
  run_rds$run_num <- run_rds$run_num + run_num_add 
  MIMout <- rbind(MIMout, run_rds)
  rm(run_rds)
}


### ADD CALCULATIONS
MIMout$CO2_of_tot <- rowSums(MIMout[,11:12])/rowSums(MIMout[,4:12])
MIMout$LIT_tot <- MIMout$LITm + MIMout$LITs
MIMout$MIM_CO <- MIMout$MICr/MIMout$MICK
MIMout$cost <- abs(MIMout$CO2_of_tot - MIMout$CO2C_prop)

### FTN TO GRAB BEST n PARAMETERIZATIONS
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


### Collect all the best psets by site
sites <- unique(MIMout$SITE)
top_df <- NULL
best_n_sets <- 100
for(i in 1:length(sites)) {
  df1 <- best_data(df = MIMout %>% filter(SITE == sites[i]), best_n = best_n_sets)
  top_df <- rbind(top_df, df1)
}

# Set SITE order
top_df <- top_df %>%
  mutate(SITE = fct_relevel(SITE, levels = "TREE", "BART", "GRSM", "SERC", "TALL", "LENO"))


### RIDGE PLOT FTN
ridge_plot <- function(plot_df, x_var, y_var, fill_var, color_pal, plot_title, plot_subtitle){
  ggplot(plot_df, aes(x = x_var, y = y_var, fill=fill_var)) +
    geom_density_ridges(scale = 2) + 
    scale_y_discrete(expand = c(0, 0)) +     # will generally have to set the `expand` option
    scale_x_continuous(expand = c(0, 0)) +   # for both axes to remove unneeded padding
    coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
    theme_ridges() +
    scale_fill_brewer(palette = color_pal) +
    labs(title=plot_title,
         subtitle=plot_subtitle) +
    xlab("") + ylab("") +
    theme(legend.position = "none") #removes the legend
}


# HOW TO GET DENSITY PEAK VALUE

sites <- c("TREE",
           "BART",
           "GRSM",
           "SERC",
           "TALL",
           "LENO")

for(i in 1:6) {
density_df <- top_df %>% filter(SITE == sites[i]) 
d_max_loc <- which.max(density(density_df$Kslope_x)$y)
print(density(density_df$Kslope_x)$x[d_max_loc])
}


ggplot(density_df, aes(Vslope_x)) + geom_density() + geom_vline(xintercept = density(density_df$Vslope_x)$x[d_max_loc])




# COST RIDGEPLOT
pCost <- ridge_plot(plot_df = top_df,
                    x_var = (top_df$cost/top_df$CO2C_prop)*100, 
                    y_var = top_df$SITE, 
                    fill_var = top_df$SITE, 
                    color_pal = "Greys", 
                    plot_title = "Cumulative %diff from incubation CO2-C target", 
                    plot_subtitle = "n = 30 lowest cost parameter sets")

suppressMessages(print(pCost))

#ggplot(top_df, aes((top_df$cost/top_df$CO2C_prop)*100)) + geom_histogram(binwidth = 0.3) 

# Vslope ridgeplot
pVslope <- ridge_plot(plot_df = top_df,
                      x_var = top_df$Vslope_x, 
                      y_var = top_df$SITE, 
                      fill_var = top_df$SITE, 
                      color_pal = "Blues", 
                      plot_title = "Vslope Multiplier", 
                      plot_subtitle = "")

suppressMessages(print(pVslope))

# Vint ridgeplot
pVint <- ridge_plot(plot_df = top_df,
                    x_var = top_df$Vint_x, 
                    y_var = top_df$SITE, 
                    fill_var = top_df$SITE, 
                    color_pal = "Blues", 
                    plot_title = "Vint Multiplier", 
                    plot_subtitle = "")

suppressMessages(print(pVint))

# Kslope ridgeplot
pKslope <- ridge_plot(plot_df = top_df,
                      x_var = top_df$Kslope_x, 
                      y_var = top_df$SITE, 
                      fill_var = top_df$SITE, 
                      color_pal = "Oranges", 
                      plot_title = "Kslope Multiplier", 
                      plot_subtitle = "")

suppressMessages(print(pKslope))

# Kint ridgeplot
pKint <- ridge_plot(plot_df = top_df,
                    x_var = top_df$Kint_x, 
                    y_var = top_df$SITE, 
                    fill_var = top_df$SITE, 
                    color_pal = "Oranges", 
                    plot_title = "Kint Multiplier", 
                    plot_subtitle = "")

suppressMessages(print(pKint))


# OUTPUT COMBO PLOT .png
top_row = ggarrange(pCost, ncol = 1, labels = c("A"))

bottom_row = ggarrange(
  pVint, pVslope,
  pKint, pKslope,
  nrow=2,
  ncol=2,
  labels = c("B", "C", "D", "E"))

png(file = "best_pset_ridge_plots.png", width = 1000, height = 1200, res=100)
ggarrange(top_row, bottom_row,
          nrow=2,
          ncol=1,
          heights = c(2,3))
dev.off()

