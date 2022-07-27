# Code from the manuscript: Insect herbivores drive sex allocation in angiosperm flowers

# Check if pgls residuals of maless ~ herbivores is explained by number of family/diversity
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

# source local functions--------------------------------------------------------
source("scripts/zzz_functions.R")

# Load data ---------------------------------------------------------------
set.seed(1234522)
d <- read_csv("data/full_data.csv")
t <- read.tree(file = "data/phylogentic_tree_s3_tre")
rownames(d) <- d$tip_name

m1 <- phylolm(maleness ~ log2(nspe+1), data = d, phy = t, model = "lambda")

# Add to the data.frame
res <- data.frame(
  tip_name = names(residuals(m1)),
  res      = as.numeric(residuals(m1))
)
d2 <- left_join(d, res)
rownames(d2) <- d2$tip_name

# Residuals ~ number of families ------------------------------------------
m2 <- phylolm(res ~ log(nfam + 2), data = d2, phy = t, model = "lambda")
summary(m2)

# Residuals ~ number of guild ------------------------------------------
m3 <- phylolm(res ~ log(ngui + 2), data = d2, phy = t, model = "lambda")
summary(m3)

# Residuals ~ number of guild ------------------------------------------
m4 <- phylolm(res ~ shan, data = d2, phy = t, model = "lambda")
summary(m4)

tab_res <- bind_rows(tidy(m2), tidy(m3), tidy(m4))
write.csv(x = tab_res, file = "output/supp/STable_pgls_res_vs_tax.csv")
