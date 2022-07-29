# Code from the manuscript: Insect herbivores drive sex allocation in angiosperm flowers

# Script to perform sensitivity analysis for intraspecific variation
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
d <- read.csv("data/full_data.csv")
t <- read.tree(file = "data/phylogentic_tree_s3_tre")
rownames(d) <- d$tip_name

# Full model coef
m1 <- readxl::read_excel("output/supp/STable_pgls_maleness_vs_nspe.xlsx")

# std_names biomass -------------------------------------------------------
ttintra <- read_csv("data/data_maleness_sd.csv") %>% as.data.frame()
row.names(ttintra) <- ttintra$tip_name

# Variancd 1X -------------------------------------------------------------
sensi_intra_1 <- sensiPhy::intra_phylm(male_avg ~ log_nspe, data = ttintra, Vy = "male_sd_1",
                                       phy = t, model = "lambda", n.intra = 10000)
saveRDS(sensi_intra_1, file = "output/temp/sensi_intra_pop_var1.Rds")

summary(sensi_intra_1)
sensi_plot(sensi_intra_1)

# Variancd 1.5X -------------------------------------------------------------
sensi_intra_1.5 <- sensiPhy::intra_phylm(male_avg ~ log_nspe, data = ttintra, Vy = "male_sd_1.5",
                                         phy = t, model = "lambda", n.intra = 10000)
saveRDS(sensi_intra_1.5, file = "output/temp/sensi_intra_pop_var1-5.Rds")

summary(sensi_intra_1.5)
sensi_plot(sensi_intra_1.5)

# Variancd 2X -------------------------------------------------------------
sensi_intra_2 <- sensiPhy::intra_phylm(male_avg ~ log_nspe, data = ttintra, Vy = "male_sd_2",
                                       phy = t, model = "lambda", n.intra = 10000)
saveRDS(sensi_intra_2, file = "output/temp/sensi_intra_pop_var2.Rds")

summary(sensi_intra_2)
sensi_plot(sensi_intra_2)


# Plot estimates ----------------------------------------------------------
# read simulations --------------------------------------------------------
intra1 <- readRDS(file = "output/temp/sensi_intra_pop_var1.Rds")
intra1.5 <- readRDS(file = "output/temp/sensi_intra_pop_var1-5.Rds")
intra2 <- readRDS(file = "output/temp/sensi_intra_pop_var2.Rds")

# Full model --------------------------------------------------------------
# 1. maleness ~ nspe------------------------------------------------------------
m1 <- phylolm(maleness ~ log2(nspe+1), data = d, phy = t, model = "lambda")
summary(m1)

## 1x variance -------------------------------------------------------------
med <- round(mean(intra1$sensi.estimates$estimate), digits = 3)
cil <- round(quantile(intra1$sensi.estimates$estimate, probs = c(.025)), digits = 3) %>% as.numeric()
cih <- round(quantile(intra1$sensi.estimates$estimate, probs = c(.975)), digits = 3) %>% as.numeric()
subtit <- paste("Average Slope = ", med,  " | [", cil, "-", cih, "]", sep = "")

g11 <-  
  ggplot(intra1$sensi.estimates, aes(x = estimate)) +
  geom_histogram(fill = "gray") +
  theme_classic(base_size = 10) +
  theme(plot.caption = element_markdown(),
        axis.text = element_text(size = 9),
        plot.subtitle = element_text(size = 10)) +
  geom_vline(xintercept = c(med, cil, cih), lty = c(1,2,2)) +
  labs(y = "Frequency", x = "Slope",
       title = "1x the intraspecific variance",
       subtitle = subtit) +
  theme(axis.title = element_text(size = 14),
        axis.text =  element_text(size = 12)) +
  xlim(-0,0.045)
g11

# intercetp
med <- round(mean(intra1$sensi.estimates$intercept), digits = 3)
cil <- round(quantile(intra1$sensi.estimates$intercept, probs = c(.025)), digits = 3) %>% as.numeric()
cih <- round(quantile(intra1$sensi.estimates$intercept, probs = c(.975)), digits = 3) %>% as.numeric()
subtit <- paste("Average Intercept = ", med,  " | [", cil, "-", cih, "]", sep = "")

g12 <-  
  ggplot(intra1$sensi.estimates, aes(x = intercept)) +
  geom_histogram(fill = "gray") +
  theme_classic(base_size = 10) +
  theme(plot.caption = element_markdown(),
        axis.text = element_text(size = 9),
        plot.subtitle = element_text(size = 10)) +
  geom_vline(xintercept = c(med, cil, cih), lty = c(1,2,2)) +
  labs(y = "Frequency", x = "Intercept",
       title = "1x the intraspecific variance",
       subtitle = subtit) +
  theme(axis.title = element_text(size = 14),
        axis.text =  element_text(size = 12)) +
  xlim(0.23,.45)
g12

nrow(intra1$sensi.estimates)

pval1 <-  sum(intra1$sensi.estimate$pval.estimate < 0.05)/10000
g13 <-  
  ggplot(intra1$sensi.estimates, aes(x = pval.estimate)) +
  geom_density(fill = "gray") +
  theme_classic(base_size = 10) +
  theme(plot.caption = element_markdown(),
        axis.text = element_text(size = 9),
        plot.subtitle = element_text(size = 10)) +
  geom_vline(xintercept = 0.05, color = "red", lty = 2) +
  labs(y = "Frequency", x = "P-value", 
       title = "1x the intraspecific variance",
       subtitle = paste("Percentage of significant regressions: ", (pval1)*100, "%")) +
  theme(axis.title = element_text(size = 14),
        axis.text =  element_text(size = 12)) +
  xlim(0.0001, 0.5)
g13

range(intra2$sensi.estimate$pval.estimate)

# 1.5x variance -------------------------------------------------------------
med <- round(mean(intra1.5$sensi.estimates$estimate), digits = 3)
cil <- round(quantile(intra1.5$sensi.estimates$estimate, probs = c(.025)), digits = 3) %>% as.numeric()
cih <- round(quantile(intra1.5$sensi.estimates$estimate, probs = c(.975)), digits = 3) %>% as.numeric()
subtit <- paste("Average Slope = ", med,  " | [", cil, "-", cih, "]", sep = "")

g21 <-  
  ggplot(intra1.5$sensi.estimates, aes(x = estimate)) +
  geom_histogram(fill = "gray") +
  theme_classic(base_size = 10) +
  theme(plot.caption = element_markdown(),
        axis.text = element_text(size = 9),
        plot.subtitle = element_text(size = 10)) +
  geom_vline(xintercept = c(med, cil, cih), lty = c(1,2,2)) +
  labs(y = "Frequency", x = "Slope",
       title = "1.5x the intraspecific variance",
       subtitle = subtit) +
  theme(axis.title = element_text(size = 14),
        axis.text =  element_text(size = 12)) +
  xlim(-0,0.045)
g21

# intercetp
med <- round(mean(intra1.5$sensi.estimates$intercept), digits = 3)
cil <- round(quantile(intra1.5$sensi.estimates$intercept, probs = c(.025)), digits = 3) %>% as.numeric()
cih <- round(quantile(intra1.5$sensi.estimates$intercept, probs = c(.975)), digits = 3) %>% as.numeric()
subtit <- paste("Average Intercept = ", med,  " | [", cil, "-", cih, "]", sep = "")

g22 <-  
  ggplot(intra1.5$sensi.estimates, aes(x = intercept)) +
  geom_histogram(fill = "gray") +
  theme_classic(base_size = 10) +
  theme(plot.caption = element_markdown(),
        axis.text = element_text(size = 9),
        plot.subtitle = element_text(size = 10)) +
  geom_vline(xintercept = c(med, cil, cih), lty = c(1,2,2)) +
  labs(y = "Frequency", x = "Intercept",
       title = "1.5x the intraspecific variance",
       subtitle = subtit) +
  theme(axis.title = element_text(size = 14),
        axis.text =  element_text(size = 12)) +
  xlim(0.23,.45)
g22


pval1.5 <-  sum(intra1.5$sensi.estimate$pval.estimate < 0.05)/10000
g23 <-  
  ggplot(intra1.5$sensi.estimates, aes(x = pval.estimate)) +
  geom_density(fill = "gray") +
  theme_classic(base_size = 10) +
  theme(plot.caption = element_markdown(),
        axis.text = element_text(size = 9),
        plot.subtitle = element_text(size = 10)) +
  geom_vline(xintercept = 0.05, color = "red", lty = 2) +
  labs(y = "Frequency", x = "P-value", 
       title = "1.5x the intraspecific variance",
       subtitle = paste("Percentage of significant regressions: ", (pval1.5)*100, "%")) +
  theme(axis.title = element_text(size = 14),
        axis.text =  element_text(size = 12)) +
  xlim(0.0001, 0.5)
g23

# 2 variance -------------------------------------------------------------
med <- round(mean(intra2$sensi.estimates$estimate), digits = 3)
cil <- round(quantile(intra2$sensi.estimates$estimate, probs = c(.025)), digits = 3) %>% as.numeric()
cih <- round(quantile(intra2$sensi.estimates$estimate, probs = c(.975)), digits = 3) %>% as.numeric()
subtit <- paste("Average Slope = ", med,  " | [", cil, "-", cih, "]", sep = "")

g31 <-  
  ggplot(intra2$sensi.estimates, aes(x = estimate)) +
  geom_histogram(fill = "gray") +
  theme_classic(base_size = 10) +
  theme(plot.caption = element_markdown(),
        axis.text = element_text(size = 9),
        plot.subtitle = element_text(size = 10)) +
  geom_vline(xintercept = c(med, cil, cih), lty = c(1,2,2)) +
  labs(y = "Frequency", x = "Slope",
       title = "2x the intraspecific variance",
       subtitle = subtit) +
  theme(axis.title = element_text(size = 14),
        axis.text =  element_text(size = 12)) +
  xlim(-0,0.045)
g31
range(intra2$sensi.estimates$estimate)


# intercetp
med <- round(mean(intra2$sensi.estimates$intercept), digits = 3)
cil <- round(quantile(intra2$sensi.estimates$intercept, probs = c(.025)), digits = 3) %>% as.numeric()
cih <- round(quantile(intra2$sensi.estimates$intercept, probs = c(.975)), digits = 3) %>% as.numeric()
subtit <- paste("Average Intercept = ", med,  " | [", cil, "-", cih, "]", sep = "")

g32 <-  
  ggplot(intra2$sensi.estimates, aes(x = intercept)) +
  geom_histogram(fill = "gray") +
  theme_classic(base_size = 10) +
  theme(plot.caption = element_markdown(),
        axis.text = element_text(size = 9),
        plot.subtitle = element_text(size = 10)) +
  geom_vline(xintercept = c(med, cil, cih), lty = c(1,2,2)) +
  labs(y = "Frequency", x = "Intercept",
       title = "2x the intraspecific variance",
       subtitle = subtit) +
  theme(axis.title = element_text(size = 14),
        axis.text =  element_text(size = 12)) +
  xlim(0.23,.45)
g32

range(intra2$sensi.estimates$intercept)
pval2 <-  sum(intra2$sensi.estimate$pval.estimate < 0.05)/10000
g33 <-  
  ggplot(intra2$sensi.estimates, aes(x = pval.estimate)) +
  geom_density(fill = "gray") +
  theme_classic(base_size = 10) +
  theme(plot.caption = element_markdown(),
        axis.text = element_text(size = 9),
        plot.subtitle = element_text(size = 10)) +
  geom_vline(xintercept = 0.05, color = "red", lty = 2) +
  labs(y = "Frequency", x = "P-value", 
       title = "2x the intraspecific variance",
       subtitle = paste("Percentage of significant regressions: ", (pval2)*100, "%")) +
  theme(axis.title = element_text(size = 14),
        axis.text =  element_text(size = 12)) +
  xlim(0.0001, 0.5)
g33

g11 + g12 + g13 +
  g21 + g22 + g23 +
  g31 + g32 + g33 + plot_annotation(tag_level = "a")

ggsave("output/supp/SFigure_intra_specific_variation_10000.pdf", width = 12, height = 9)
