# Code from the manuscript: Insect herbivores drive sex allocation in angiosperm flowers

# Script to run PGLS analysis between outcrossing rate and flower maleness
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
library(psych)       
library(GGally)    
library(ggtree)

# source local functions--------------------------------------------------------
source("scripts/zzz_functions.R")

# Load data ---------------------------------------------------------------
tm_data <- read.csv("data/data_outcrossing_rate.csv")
rownames(tm_data) <- tm_data$Accepted_species
tm_tree <-  read.tree("data/phylogentic_tree_s3_tm_analysis.tre")

# Match data and phylogeny ------------------------------------------------
dc <- match_dataphy(tm ~ maleness, data = tm_data, phy = tm_tree)

# PGLS regression ---------------------------------------------------------
modphy <- phylolm(tm ~ maleness, data = dc$data, phy = dc$phy, model = "lambda", boot = 1000)
summary(modphy)
mspe_slope_CI <- modphy$bootconfint95[,2]

# save table
write_xlsx(tidy(modphy), path = "output/supp/STable_outcrossing_maleness.xlsx")

# Plot regression-----------------------------------------------------------
bs <- 14
pcolor = "steelblue"
lcolor = "gray"
psize = 3
cc   <- coef(modphy)
ci   <- modphy$bootstrap
x    <- tm_data$maleness
y    <- pred_x(a = cc[1], x = x, b = cc[2])
r2   <- round(modphy$r.squared, digits = 3)
pva  <- round(summary(modphy)[[2]][2,6], digits = 5)

pred1 <- data.frame(tm = tm_data$tm, maleness = tm_data$maleness, y, x)

#PGLS
g1 <-
  ggplot(pred1, aes(y = tm, x = maleness)) +
  geom_point(fill = pcolor, color = "black", shape = 21, size = psize, alpha = .75) +
  geom_abline(intercept = ci[,1], slope = ci[,2], color = lcolor, alpha = .01, size = 1) +
  geom_line(aes(y = y, x = maleness), size = 1.5) +
  scale_x_continuous(limits = c(0.1,.9)) +
  theme_classic(base_size = bs) +
  labs(y = "Outcrossing rate [tm]", 
       x = "Flower maleness",
       subtitle = "") +
  theme(plot.caption = element_markdown(),
        axis.text = element_text(size = 9),
        plot.subtitle = element_text(size = 10)) +
  annotate(geom = "text",  x = .1, hjust = -.15, y = 0, 
           label = paste("Y = ", round(cc[2], digits = 3), "X",
                         " + ",  round(cc[1], digits = 3),
                         " | p = ", pva,
                         sep = ""),
           size = 3) +
  scale_y_continuous(breaks = seq(0,1,.2)) +
  scale_x_continuous(breaks = seq(0,1,.2))

# Phylogenetic tree
ptree <- 
  ggtree(dc$phy) +
  geom_tiplab(size = 2) +
  xlim(0, 200)

ptree / g1 + plot_annotation(tag_levels = "a")

# Save figure
ggsave(filename = 'output/supp/SFigure_pgls_outcrossing_maleness.png', 
       width = 93, height = 186, units = "mm")
