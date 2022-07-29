# Code from the manuscript: Insect herbivores drive sex allocation in angiosperm flowers

# Script to run perform sensitivity analysis controling for life span
# Last update: 2022.07.27
# Author: Gustavo Brant Paterno (paternogbc@gmail.com)

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
library(TNRS)
library(inspectdf)

# source local functions--------------------------------------------------------
source("scripts/zzz_functions.R")

# Load data ---------------------------------------------------------------
set.seed(1234522)
d <- read_csv("data/full_data.csv")
t <- read.tree(file = "data/phylogentic_tree_s3_tre")
rownames(d) <- d$tip_name

m1 <- phylolm(maleness ~ log2(nspe+1), data = d, phy = t, model = "lambda")
summary(m1)

# Load life span data Biolflor ---------------------------------------------
d3 <- openxlsx::read.xlsx("data/data_lifespan.xlsx")
rownames(d3) <- d3$tip_name

inspectdf::inspect_num(d3) %>% inspectdf::show_plot()
inspectdf::inspect_cat(d3) %>% inspectdf::show_plot()

# 1. Adding life span in the model as a covariate
m11 <- phylolm(maleness ~ log2(nspe+1) * life_span, data = d3, phy = t, model = "lambda")
summary(m11)
tidy(m11)

# 2. maleness against life span
m11 <- phylolm(maleness ~ life_span, data = d3, phy = t, model = "lambda")
summary(m11)
tidy(m11)

# Plot overview
ggplot() +
  scale_x_continuous(trans = "log2", 
                     breaks = c(1, 2, 4, 8, 16, 32, 64, 128, 256, 512),
                     limits = c(1, 512)) +
  scale_y_continuous(breaks = seq(0,.8,.1)) +
  geom_point(data =  d3 %>% drop_na(life_span), alpha = .75, size = 2, aes(y = maleness, x = nspe+1, shape = life_span, color = life_span)) +
  geom_abline(intercept = coef(m1)[1], slope = coef(m1)[2]) +
  theme_classic(base_size = 14) +
  labs(y = "Flower maleness", 
       x = "Number of insect species",
       subtitle = "",
       color = "Life span",
       shape = "Life span") +
  scale_color_discrete() +
  theme(axis.text = element_text(size = 9),
        plot.subtitle = element_text(size = 10))
ggsave(filename = "output/supp/SFigure_sensitivity_life_span.png",
       width = 6, height = 5)
