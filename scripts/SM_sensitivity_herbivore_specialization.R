# Code from the manuscript: Insect herbivores drive sex allocation in angiosperm flowers

# Script to run sensitivity analysis on the pgls regression between flower 
# maleness against the number of insect species
# Herbivore specialization
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
d1spe <- read.csv("data/diver_processed_1spe.csv")
d1spe <- left_join(d1spe, d %>% dplyr::select(tip_name, maleness))
rownames(d1spe) <- d1spe$tip_name

# only herbivores that feed on 1gen
d1gen <- read.csv("data/diver_processed_1gen.csv")
d1gen <- left_join(d1gen, d %>% dplyr::select(tip_name, maleness))
rownames(d1gen) <- d1gen$tip_name

# only herbivores that feed on 1fam
d1fam <- read.csv("data/diver_processed_1fam.csv")
d1fam <- left_join(d1fam, d %>% dplyr::select(tip_name, maleness))
rownames(d1fam) <- d1fam$tip_name

# only herbivores that feed on mfam
dmfam <- read.csv("data/diver_processed_mfam.csv")
dmfam <- left_join(dmfam, d %>% dplyr::select(tip_name, maleness))
rownames(dmfam) <- dmfam$tip_name

# Calculate specialists and generalists
# quick check
d1spe$tip_name == d1gen$tip_name
d1spe$tip_name == d1fam$tip_name
d1spe$tip_name == dmfam$tip_name

# Define specialization groups
monophagous     <- d1spe$nspe 
famspecialist   <- d1fam$nspe + d1gen$nspe
generalists     <- dmfam$nspe

d1spe$monophagous    <- monophagous
d1spe$famspecialist  <- famspecialist
d1spe$generalists    <- generalists

# Check numbers
d <- d[match(d1spe$tip_name, d$tip_name), ]
any(!d$nspe == c(d1spe$monophagous + d1spe$famspecialist + d1spe$generalists))

# Add area
d1spe$range <- d$range

# Pure regression---------------------------------------------------------------
# maleness ~ nspe---------------------------------------------------------------
m1 <- phylolm(maleness ~ log2(nspe+1), data = d, phy = t, model = "lambda")
summary(m1)

# Figure 3  ---------------------------------------------------------------
## Monophagous -------------------------------------------------------------
mspe <- phylolm(maleness ~ log2(monophagous+1), data = d1spe, phy = t,
                model = "lambda", boot = 1000)
summary(mspe)
mspe_slope_CI <- mspe$bootconfint95[,2]
R2.pred(mspe, phy = t)

# save table
write_xlsx(tidy(mspe), path = "output/supp/STable_sensi_taxa_monophagous_only.xlsx")

# plot regression
bs <- 12
pcolor = "steelblue"
lcolor = "gray"
psize = 3
cc   <- coef(mspe)
ci   <- mspe$bootstrap
x    <- d1spe$monophagous +1
y    <- pred_log2x(a = cc[1], x = x, b = cc[2])
r2   <- round(mspe$r.squared, digits = 3)
pva  <- round(summary(mspe)[[2]][2,6], digits = 5)

pred1 <- data.frame(nspe = d1spe$monophagous, maleness = d1spe$maleness, y, x)

g5 <-
  ggplot(pred1, aes(y = maleness, x = monophagous+1)) +
  scale_x_continuous(trans = "log2", 
                     breaks = c(1, 2, 4, 8, 16, 32, 64, 128, 256),
                     limits = c(1, 300)) +
  scale_y_continuous(breaks = seq(0,.8,.1)) +
  geom_point(fill = pcolor, color = "black", shape = 21, size = psize, alpha = .5) +
  geom_abline(intercept = ci[,1], slope = ci[,2], color = lcolor, alpha = .01, size = .25) +
  geom_line(aes(y = y, x = nspe+1), size = 1.5) +
  theme_classic(base_size = bs) +
  labs(y = "Flower maleness", 
       x = "N. of monophagous herbivores",
       subtitle = "") +
  theme(plot.caption = element_markdown(),
        axis.text = element_text(size = 9)) +
  annotate(geom = "text", x = 0, hjust = -.2,
           y = 0, label = paste("Y = ", round(cc[2], digits = 3), "X",
                                " + ",  round(cc[1], digits = 3),
                                " | p = ", pva,
                                sep = ""), size = 3) 

g5

## Family specialists -------------------------------------------------------------
mfam <- phylolm(maleness ~ log2(famspecialist + 1), data = d1spe, phy = t,
                model = "lambda", boot = 1000)
summary(mfam)
mgen_slope_CI <- mfam$bootconfint95[,2]
R2.pred(mfam, phy = t)

# save table
write_xlsx(tidy(mfam), path = "output/supp/STable_sensi_taxa_mfamily_specialist_only.xlsx")

# plot regression
bs <- 12
pcolor = "steelblue"
lcolor = "gray"
psize = 3
cc   <- coef(mfam)
ci   <- mfam$bootstrap
x    <- d1spe$famspecialist +1
y    <- pred_log2x(a = cc[1], x = x, b = cc[2])
r2   <- round(mfam$r.squared, digits = 3)
pva  <- round(summary(mfam)[[2]][2,6], digits = 5)

pred1 <- data.frame(nspe = d1spe$famspecialist, maleness = d1spe$maleness, y, x)

g6 <-
  ggplot(pred1, aes(y = maleness, x = famspecialist+1)) +
  scale_x_continuous(trans = "log2", 
                     breaks = c(1, 2, 4, 8, 16, 32, 64),
                     limits = c(1, 80)) +
  scale_y_continuous(breaks = seq(0,.8,.1)) +
  geom_point(fill = pcolor, color = "black", shape = 21, size = psize, alpha = .5) +
  geom_abline(intercept = ci[,1], slope = ci[,2], color = lcolor, alpha = .01, size = .25) +
  geom_line(aes(y = y, x = nspe+1), size = 1.5) +
  theme_classic(base_size = bs) +
  labs(y = "Flower maleness", 
       x = "N. of family specialist herbivores",
       subtitle = "") +
  theme(plot.caption = element_markdown(),
        axis.text = element_text(size = 9)) +
  annotate(geom = "text", x = 0, hjust = -.2,
           y = 0, label = paste("Y = ", round(cc[2], digits = 3), "X",
                                " + ",  round(cc[1], digits = 3),
                                " | p = ", pva,
                                sep = ""), size = 3) 


g6

## Generalists only-------------------------------------------------------------
mgen <- phylolm(maleness ~ log2(generalists + 1), data = d1gen, phy = t,
                model = "lambda", boot = 1000)
summary(mgen)
mgen_slope_CI <- mgen$bootconfint95[,2]
R2.pred(mgen, phy = t)

# save table
write_xlsx(tidy(mgen), path = "output/supp/STable_sensi_taxa_mgeneralists_only.xlsx")

# plot regression
bs <- 12
pcolor = "steelblue"
lcolor = "gray"
psize = 3
cc   <- coef(mgen)
ci   <- mgen$bootstrap
x    <- d1spe$generalists +1
y    <- pred_log2x(a = cc[1], x = x, b = cc[2])
r2   <- round(mgen$r.squared, digits = 3)
pva  <- round(summary(mgen)[[2]][2,6], digits = 5)

pred1 <- data.frame(nspe = d1spe$generalists, maleness = d1spe$maleness, y, x)

g7 <-
  ggplot(pred1, aes(y = maleness, x = generalists+1)) +
  scale_x_continuous(trans = "log2", 
                     breaks = c(1, 2, 4, 8, 16, 32, 64, 128, 256),
                     limits = c(1, 256)) +
  scale_y_continuous(breaks = seq(0,.8,.1)) +
  geom_point(fill = pcolor, color = "black", shape = 21, size = psize, alpha = .5) +
  geom_abline(intercept = ci[,1], slope = ci[,2], color = lcolor, alpha = .01, size = .25) +
  geom_line(aes(y = y, x = nspe+1), size = 1.5) +
  theme_classic(base_size = bs) +
  labs(y = "Flower maleness", 
       x = "N. of generalist herbivores",
       subtitle = "") +
  theme(plot.caption = element_markdown(),
        axis.text = element_text(size = 9)) +
  annotate(geom = "text", x = 0, hjust = -.2,
           y = 0, label = paste("Y = ", round(cc[2], digits = 3), "X",
                                " + ",  round(cc[1], digits = 3),
                                " | p = ", pva,
                                sep = ""), size = 3) 


g7

# merge plots
gtaxa <- g5 / g6 / g7 + plot_annotation(tag_levels = "a")
gtaxa + theme(axis.text = element_text(size = 10))

ggsave(gtaxa, filename = 'output/Figure_3.pdf', device = cairo_pdf,
       width = 82, height = 246, units = "mm")


# Slope comparison between specialists and generalists
mspe_slope_CI 
mspe$coefficients
mgen_slope_CI
mgen$coefficients

# Controling for specialists and range area----------------------------------
modspe <- lm(log2(monophagous+1) ~ (range), d1spe)
modgen <- lm(log2(monophagous+1) ~ (range), d1spe)

mod_res <- data.frame(tip_name = names(residuals(modspe)), 
                      specialists_res = residuals(modspe),
                      generalists_res = residuals(modgen))

dres <- left_join(d1spe, mod_res)
nrow(dres)
nrow(d1spe)
rownames(dres) <- dres$tip_name

# fit model with residual specialists
mres_specialists <- phylolm(maleness ~ specialists_res, 
                            data = dres, phy = t, model = "lambda", boot = 1000)
summary(mres_specialists)

# fit model with residual generalists
mres_generalists <- phylolm(maleness ~ generalists_res, 
                            data = dres, phy = t, model = "lambda", boot = 1000)
summary(mres_generalists)
# save table
write_xlsx(tidy(mres_specialists), path = "output/supp/STable_controlled_pgls_species_range_residuals_specialists.xlsx")
write_xlsx(tidy(mres_generalists), path = "output/supp/STable_controlled_pgls_species_range_residuals_generalists.xlsx")

# plot regression specialiosts
cc   <- coef(mres_specialists)
ci   <- mres_specialists$bootstrap
x    <- dres$specialists_res
y    <- pred_x(a = cc[1], x = x, b = cc[2])
r2   <- round(mres_specialists$r.squared, digits = 3)
pva  <- round(summary(mres_specialists)[[2]][2,6], digits = 5)

pred1 <- data.frame(nspe_res = dres$specialists_res, maleness = d1spe$maleness, y, x)

g1a <-
  ggplot(pred1, aes(y = maleness, x = nspe_res)) +
  scale_y_continuous(breaks = seq(0,.8,.1)) +
  geom_abline(intercept = ci[,1], slope = ci[,2], color = lcolor, alpha = .01, size = 1) +
  geom_point(fill = pcolor, color = "black", shape = 21, size = psize, alpha = .75) +
  geom_line(aes(y = y, x = nspe_res), size = 1.5) +
  theme_classic(base_size = 14) +
  labs(y = "Flower maleness", 
       x = "Residual number of insect species",
       subtitle = "[specialists: feed on 1 genera only]") +
  annotate(geom = "text", x = 1.5, y = 0.1, label = paste("Y = ", round(cc[1], digits = 3),
                                                          " + ", round(cc[2], digits = 3), "X",
                                                          " | p = ", pva,
                                                          sep = ""), size = 3)


# plot regression generalists
cc   <- coef(mres_generalists)
ci   <- mres_generalists$bootstrap
x    <- dres$generalists_res
y    <- pred_x(a = cc[1], x = x, b = cc[2])
r2   <- round(mres_generalists$r.squared, digits = 3)
pva  <- round(summary(mres_generalists)[[2]][2,6], digits = 5)

pred1 <- data.frame(nspe_res = dres$generalists_res, maleness = d1spe$maleness, y, x)

g1b <-
  ggplot(pred1, aes(y = maleness, x = nspe_res)) +
  scale_y_continuous(breaks = seq(0,.8,.1)) +
  geom_abline(intercept = ci[,1], slope = ci[,2], color = lcolor, alpha = .01, size = 1) +
  geom_point(fill = pcolor, color = "black", shape = 21, size = psize, alpha = .75) +
  geom_line(aes(y = y, x = nspe_res), size = 1.5) +
  theme_classic(base_size = 14) +
  labs(y = "Flower maleness", 
       x = "Residual number of insect species",
       subtitle = "[generalists: feed on more than 1 genera]") +
  annotate(geom = "text", x = 1.5, y = 0.1, label = paste("Y = ", round(cc[1], digits = 3),
                                                          " + ", round(cc[2], digits = 3), "X",
                                                          " | p = ", pva,
                                                          sep = ""), size = 3)

grange <- g1a + g1b + plot_annotation(tag_levels = "a") 
ggsave(grange, filename = "output/supp/SFigure_sensitivity_occupancy_specialists_generalists.png",
       width = 9, height = 5)
