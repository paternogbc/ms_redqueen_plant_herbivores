# Code from the manuscript: Insect herbivores drive sex allocation in angiosperm flowers

# Script to run main pgls between flower maleness and 4 alternative 
# proxies of herbivore evolutionary pressure.
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
library(rr2)

# source local functions--------------------------------------------------------
source("scripts/zzz_functions.R")

# Load data ---------------------------------------------------------------
set.seed(1234522)
d <- read_csv("data/full_data.csv")
t <- read.tree(file = "data/phylogentic_tree_s3_tre")
rownames(d) <- d$tip_name

# Fit PGLS----------------------------------------------------------------------
# 1. maleness ~ nspe------------------------------------------------------------
m1 <- phylolm(maleness ~ log2(nspe+1), data = d, phy = t, model = "lambda", boot = 1000)
summary(m1)
R2(m1, phy = t)

# Model Diagnostics (USING CAPER)-----
cdat <- comparative.data(phy = t, data = d,
                            names.col = tip_name, vcv = TRUE,
                            na.omit = FALSE, warn.dropped = TRUE)
model.pgls <- pgls(maleness ~ log2(nspe+1),
                   data = cdat, lambda = "ML")
summary(model.pgls)
par(mfrow = c(2, 2))
plot(model.pgls)

# withouth phylogeny
summary(lm(maleness ~ log2(nspe+1), data = d))

# save model
writexl::write_xlsx(x = tidy(m1), path = "output/supp/STable_pgls_maleness_vs_nspe.xlsx")
saveRDS(object = m1, file = 'output/temp/pgls_maleness_vs_nspe.RDs')

# 2. maleness ~ nfam------------------------------------------------------------
m2 <- phylolm(maleness ~ log2(nfam+1), data = d, phy = t, model = "lambda", boot = 1000)
summary(m2) 
R2(m2, phy = t)

# withouth phylogeny
summary(lm(maleness ~ log2(nfam+1), data = d))

# save model
writexl::write_xlsx(x = tidy(m2), path = "output/supp/STable_pgls_maleness_vs_nfam.xlsx")
saveRDS(object = m2, file = 'output/temp/pgls_maleness_vs_nfam.RDs')

# 3. maleness ~ ngui----------------------------------------------------------------
m3 <- phylolm(maleness ~ log2(ngui+1), data = d, phy = t, model = "lambda", boot = 1000)
summary(m3) 
R2(m3, phy = t)

# withouth phylogeny
summary(lm(maleness ~ log2(ngui+1), data = d))

# save model
writexl::write_xlsx(x = tidy(m3), path = "output/supp/STable_pgls_maleness_vs_ngui.xlsx")
saveRDS(object = m3, file = 'output/temp/pgls_maleness_vs_ngui.RDs')

# 4. maleness ~ shan----------------------------------------------------------------
m4 <- phylolm(maleness ~ shan, data = d, phy = t, model = "lambda", boot = 1000)
summary(m4) 
R2(m4, phy = t)

# withouth phylogeny
summary(lm(maleness ~ shan, data = d))

# save model
writexl::write_xlsx(x = tidy(m4), path = "output/supp/STable_pgls_maleness_vs_shan.xlsx")
saveRDS(object = m4, file = 'output/temp/pgls_maleness_vs_shan.RDs')

tab_r2 <- data.frame(
  model = c("N. of sp", "N. of fam", "N. of gui", "Shannon"),
  r2pred = c(R2(m1, phy = t)[[3]],
             R2(m2, phy = t)[[3]],
             R2(m3, phy = t)[[3]],
             R2(m4, phy = t)[[3]])
)
writexl::write_xlsx(x = tab_r2, path = "output/supp/STable_pgls_r2_pred.xlsx")
tidy(m1)
