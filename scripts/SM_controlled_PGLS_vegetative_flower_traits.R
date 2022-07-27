# Code from the manuscript: Insect herbivores drive sex allocation in angiosperm flowers

# Script to run pgls regressions between flower maleness against number 
# of insect species after controlling for vegetative and reproductive traits.
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
library(caper)     

# source local functions--------------------------------------------------------
source("scripts/zzz_functions.R")
set.seed(1234522)

# Load data ---------------------------------------------------------------
d <- read.csv("data/full_data.csv")
t <- read.tree(file = "data/phylogentic_tree_s3_tre")
rownames(d) <- d$tip_name

# prepare subsets of data-------------------------------------------------------
# Vegetative traits
dveg  <- match_dataphy(maleness ~ log2(nspe+1) + sla_imp + height, data = d, phy = t)
dsla  <- match_dataphy(maleness ~ log2(nspe+1) + sla_imp, data = d, phy = t)
dhei  <- match_dataphy(maleness ~ log2(nspe+1) + height, data = d, phy = t)

# flower traits
dflo <- match_dataphy(maleness ~ log2(nspe+1) + nectar + color + polli + shape, data = d, phy = t)
dnec <- match_dataphy(maleness ~ log2(nspe+1) + nectar, data = d, phy = t)
dcol <- match_dataphy(maleness ~ log2(nspe+1) + color, data = d, phy = t)
dsha <- match_dataphy(maleness ~ log2(nspe+1) + shape, data = d, phy = t)
dpol <- match_dataphy(maleness ~ log2(nspe+1) + polli, data = d, phy = t)

# 1. Vegetative traits------------------------------------------
# sla_imp-------------------------------------------------------------------------
mv1 <- phylolm(maleness ~ log2(nspe+1) + sla_imp, 
               data = dveg$data, phy = dveg$phy, 
               model = "lambda", boot = 1000)
summary(mv1)
writexl::write_xlsx(x = tidy(mv1), path = "output/supp/STable_controlled_pgls_sla.xlsx")

# fit with caper
mv1.2 <- pgls(maleness ~ log2(nspe+1) + sla_imp, 
              data = comparative.data(t, dsla$data %>% dplyr::select(maleness, tip_name, nspe, sla_imp),
                                      tip_name, vcv=TRUE, vcv.dim=3),
              lambda='ML')
anova(mv1.2)
# save anova table
writexl::write_xlsx(x = tidy(anova(mv1.2)), 
                    path = "output/supp/STable_controlled_pgls_sla_anova.xlsx")

# height-------------------------------------------------------------------------
mv2 <- phylolm(maleness ~ log2(nspe+1) + height, 
               data = dveg$data, phy = dveg$phy, 
               model = "lambda", boot = 1000)
summary(mv2)
writexl::write_xlsx(x = tidy(mv2), path = "output/supp/STable_controlled_pgls_height.xlsx")

# fit with caper
mv2.2 <- pgls(maleness ~ log2(nspe+1) + height, 
              data = comparative.data(t, dveg$data %>% dplyr::select(maleness, tip_name, nspe, height),
                                      tip_name, vcv=TRUE, vcv.dim=3),
              lambda='ML')
anova(mv2.2)
# save anova table
writexl::write_xlsx(x = tidy(anova(mv2.2)), 
                    path = "output/supp/STable_controlled_pgls_height_anova.xlsx")

# all vegetative----------------------------------------------------------------
mv3 <- phylolm(maleness ~ log2(nspe+1) + height + sla_imp, 
               data = dveg$data, phy = dveg$phy, 
               model = "lambda", boot = 1000)
summary(mv3)
writexl::write_xlsx(x = tidy(mv3), path = "output/supp/STable_controlled_pgls_all_vegetative.xlsx")

# fit with caper
mv3.2 <- pgls(maleness ~ log2(nspe+1) + height + sla_imp, 
              data = comparative.data(t, dveg$data %>% dplyr::select(maleness, tip_name, nspe, height, sla_imp),
                                      tip_name, vcv=TRUE, vcv.dim=3),
              lambda='ML')
anova(mv3.2)
# save anova table
writexl::write_xlsx(x = tidy(anova(mv3.2)), 
                    path = "output/supp/STable_controlled_pgls_all_vegetative_anova.xlsx")

# 2. Flower traits------------------------------------------
# nectar-------------------------------------------------------------------------
mf1 <- phylolm(maleness ~ log2(nspe+1) + nectar, 
               data = dnec$data, phy = dnec$phy, 
               model = "lambda", boot = 1000)
summary(mf1)
writexl::write_xlsx(x = tidy(mf1), path = "output/supp/STable_controlled_pgls_nectar.xlsx")

# fit with caper
mf1.2 <- pgls(maleness ~ log2(nspe+1) + nectar, 
              data = comparative.data(t, dnec$data %>% dplyr::select(maleness, tip_name, nspe, nectar),
                                      tip_name, vcv=TRUE, vcv.dim=3),
              lambda='ML')
anova(mf1.2)
# save anova table
writexl::write_xlsx(x = tidy(anova(mf1.2)), 
                    path = "output/supp/STable_controlled_pgls_nectar_anova.xlsx")

# shape-------------------------------------------------------------------------
mf2 <- phylolm(maleness ~ log2(nspe+1) + shape, 
               data = dsha$data, phy = dsha$phy, 
               model = "lambda", boot = 1000)
summary(mf2)
writexl::write_xlsx(x = tidy(mf2), path = "output/supp/STable_controlled_pgls_shape.xlsx")

# fit with caper
mf2.2 <- pgls(maleness ~ log2(nspe+1) + shape, 
              data = comparative.data(t, dsha$data %>% dplyr::select(maleness, tip_name, nspe, shape),
                                      tip_name, vcv=TRUE, vcv.dim=3),
              lambda='ML')
summary(mf2.2)

# save anova table
writexl::write_xlsx(x = tidy(anova(mf2.2)), 
                    path = "output/supp/STable_controlled_pgls_shape_anova.xlsx")

# polli-------------------------------------------------------------------------
mf3 <- phylolm(maleness ~ log2(nspe+1) + polli, 
               data = dpol$data, phy = dpol$phy, 
               model = "lambda", boot = 1000)
summary(mf3)

writexl::write_xlsx(x = tidy(mf3), path = "output/supp/STable_controlled_pgls_polli.xlsx")

# fit with caper
mf3.2 <- pgls(maleness ~ log2(nspe+1) + polli, 
              data = comparative.data(t, dpol$data %>% dplyr::select(maleness, tip_name, nspe, polli),
                                      tip_name, vcv=TRUE, vcv.dim=3),
              lambda='ML')
summary(mf3.2)
anova(mf3.2)

# save anova table
writexl::write_xlsx(x = tidy(anova(mf3.2)), 
                    path = "output/supp/STable_controlled_pgls_polli_anova.xlsx")

# color-------------------------------------------------------------------------
mf4 <- phylolm(maleness ~ log2(nspe+1) + color, 
               data = dcol$data, phy = dcol$phy, 
               model = "lambda", boot = 1000)
summary(mf4)
writexl::write_xlsx(x = tidy(mf4), path = "output/supp/STable_controlled_pgls_color.xlsx")

# fit with caper
mf4.2 <- pgls(maleness ~ log2(nspe+1) + color, 
              data = comparative.data(t, dcol$data %>% dplyr::select(maleness, tip_name, nspe, color),
                                      tip_name, vcv=TRUE, vcv.dim=3),
              lambda='ML')
summary(mf4.2)
anova(mf4.2)

# save anova table
writexl::write_xlsx(x = tidy(anova(mf4.2)), 
                    path = "output/supp/STable_controlled_pgls_color_anova.xlsx")

# all reproductive-------------------------------------------------------------------------
mf5 <- phylolm(maleness ~ log2(nspe+1) + nectar + shape + polli + color, 
               data = dflo$data, phy = dflo$phy, 
               model = "lambda", boot = 1000)
summary(mf5)
writexl::write_xlsx(x = tidy(mf5), path = "output/supp/STable_controlled_pgls_all_reproductive.xlsx")

# Table 1-----------------------------------------------------------------------
# sample size
N <- c(mv2$n, mv1$n, mf3$n, mf4$n, mf2$n, mf1$n)

ssla <- summary(mv1)[[2]][2, ]
shei <- summary(mv2)[[2]][2, ]

snec <- summary(mf1)[[2]][2, ]
ssha <- summary(mf2)[[2]][2, ]
spol <- summary(mf3)[[2]][2, ]
scol <- summary(mf4)[[2]][2, ]


sall <- rbind(shei, ssla, spol, scol, ssha, snec)

tibble(
  covariate = 
    c("Height", "SLA", 
      "Pollinator group", "Flower color", "Flower shape", "Nectar offer"
    ),
  N,
) %>% 
  cbind(sall) %>% 
  as_tibble() -> tab_cont_pgls

tab_cont_pgls
write_xlsx(x = tab_cont_pgls, path = "output/supp/STable_summary_controlled_pgls_regressions.xlsx")

# Coefficients plot--------------------------------------------------------------
g1 <- 
  ggplot(tab_cont_pgls, aes(Estimate, x = covariate, label = paste("N = ", N,  "| p-value =",round(p.value, digits = 5)))) +
  geom_pointrange(aes(ymin = lowerbootCI, max = upperbootCI), color = "steelblue", size = 2) +
  geom_text(nudge_y = .001, vjust = -2) +
  geom_hline(yintercept = 0, lty = 2) +
  scale_x_discrete(limits = c("Height", "SLA", "Pollinator group", "Flower color", "Flower shape", "Nectar offer")[6:1]) +
  coord_flip() +
  theme_bw(base_size = 18) +
  labs(y = "Estimate [number of insect species]",
       x = "Covariate");g1
ggsave(plot = g1, filename = "output/supp/SFigure_controlled_pgls_estimates.png")
