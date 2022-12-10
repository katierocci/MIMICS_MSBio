library(dplyr)
library(ggplot2)
library(forcats)

setwd("C:/github/MIMICS_MSBio/LIDET")

base <- read.csv("output/NEON_MIMICS_Litter_Decomp_output.csv")

# Get mean years to decomp 50% of litter
base_year_50 <- base %>% group_by(site, type) %>% slice(which.min(abs(bag - 50))) %>% 
            group_by(site) %>% summarize(year50 = mean(yr))


adj <- read.csv("output/NEON_MIMICS_Litter_Decomp_output_ADJ_3.csv")

# Get mean years to decomp 50% of litter
adj_year_50 <- adj %>% group_by(site, type) %>% slice(which.min(abs(bag - 50))) %>% 
  group_by(site) %>% summarize(year50 = mean(yr))


comp_df <- left_join(base_year_50, adj_year_50, by="site")
comp_df$diff <- comp_df$year50.x - comp_df$year50.y 
comp_df$p_diff <- (comp_df$year50.x/comp_df$year50.y)-1

#write.csv(comp_df, "comp_df.csv")

NEON_MAT <- read.csv("NEON_MAT.csv")
colnames(NEON_MAT)[1] <- "site"
plot_df <- left_join(comp_df, NEON_MAT, by=)

# Set SITE order
plot_df <- plot_df %>% 
  mutate(site = fct_relevel(site, levels = "TREE","BART","SERC","GRSM","TALL","LENO"))#,"LUQ"))

#png("NEON6_LIT_Bag_diff.png", width = 2200, height = 1200, res = 300)
ggplot(plot_df, aes(x=site, y=p_diff*100, fill=MAT)) + geom_bar(stat="identity")+
  scale_fill_gradient2(midpoint=mean(plot_df$MAT), low="blue", mid="orange",
                        high="red" ) +
  labs(title="Parameter adjustment effect on litter bag decomposition rates",
          subtitle ="NEON study sites (n = 6)") +
  xlab("Site") +
  ylab("Change in deocmposition rate (%)") +
  labs(fill="MAT (°C)") +
  #scale_y_continuous(limits=c(-10,30), breaks = seq(-10, 30, 10)) +
  theme_minimal() #+ theme(axis.text = element_text(size = 16))  
#dev.off()


