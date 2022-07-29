# Code from the manuscript: Insect herbivores drive sex allocation in angiosperm flowers

# Script to run sensitivity analysis on the pgls regression between flower 
# maleness against the number of insect species
# Excluding pollinators
# Last update: 2022.07.27
# Author: Gustavo Brant Paterno (paternogbc@gmail.com)

# Load packages -----------------------------------------------------------
library(tidyverse) 
library(sensiPhy)  
library(phylolm)  
library(patchwork) 
library(broom)     
library(writexl)   
library(ggtext)    
library(ggExtra)   
library(ggtext)    
library(cowplot)   
library(factoextra)
library(caper)     
library(rr2)

# source local functions--------------------------------------------------------
source("scripts/zzz_functions.R")
set.seed(1234522)

# Load data ---------------------------------------------------------------
d <- read.csv("data/full_data.csv")
t <- read.tree(file = "data/phylogentic_tree_s3_tre")
rownames(d) <- d$tip_name

# only herbivores that feed on 1spe
d1spe <- read.csv("data/diver_processed_exclude_pollinators.csv")
d1spe <- left_join(d1spe, d %>% dplyr::select(tip_name, maleness))
rownames(d1spe) <- d1spe$tip_name

# Pure regression---------------------------------------------------------------
# Excluding pottential pollinators
# maleness ~ nspe---------------------------------------------------------------
m1spe <- phylolm(maleness ~ log2(nspe+1), data = d1spe, phy = t, model = "lambda",
                 boot = 1000)
summary(m1spe)
R2.pred(m1spe, phy = t)

# save table
write_xlsx(tidy(m1spe), path = "output/supp/STable_sensi_exclude_pollinators.xlsx")

# plot regression
bs <- 14
pcolor = "steelblue"
lcolor = "gray"
psize = 3
cc   <- coef(m1spe)
ci   <- m1spe$bootstrap
x    <- d1spe$nspe +1
y    <- pred_log2x(a = cc[1], x = x, b = cc[2])
r2   <- round(m1spe$r.squared, digits = 3)
pva  <- round(summary(m1spe)[[2]][2,6], digits = 5)

pred1 <- data.frame(nspe = d1spe$nspe, maleness = d1spe$maleness, y, x)

g1 <-
  ggplot(pred1, aes(y = maleness, x = nspe+1)) +
  scale_x_continuous(trans = "log2", 
                     breaks = c(1, 2, 4, 8, 16, 32, 64, 128, 256, 512),
                     limits = c(1, 512)) +
  scale_y_continuous(breaks = seq(0,.8,.1)) +
  geom_point(fill = pcolor, color = "black", shape = 21, size = psize, alpha = .75) +
  geom_abline(intercept = ci[,1], slope = ci[,2], color = lcolor, alpha = .01, size = 1) +
  geom_line(aes(y = y, x = nspe+1), size = 1.5) +
  theme_classic(base_size = bs) +
  labs(y = "Flower maleness", 
       x = "Number of insect species",
       subtitle = "Excluding potential pollinators") +
  theme(plot.caption = element_markdown(),
        axis.text = element_text(size = 9),
        plot.subtitle = element_text(size = 10)) +
  annotate(geom = "text",  x = 100, y = 0.1, label = paste("Y = ",  round(cc[1], digits = 3),
                                                           " + ", round(cc[2], digits = 3), "X",
                                                           " | p = ", pva,
                                                           sep = ""), size = 3)

g1

# Save plots
ggsave(g1, filename = "output/supp/SFigure_sensitivity_exclude_pollinators.png",
       width = 6, height = 5)
saveRDS(g1, file =  "output/temp/sensi_exclude_polli.Rds")
