# Code from the manuscript: Insect herbivores drive sex allocation in angiosperm flowers

# Script to run sensitivity analysis for sampling, taxonomic, phylogenetic and also 
# regressions controlling for species distribution rage. 
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
library(ggtree)
library(ggtext)    
library(cowplot)   
library(factoextra)
library(caper)     

# source local functions--------------------------------------------------------
source("scripts/zzz_functions.R")
set.seed(1234522)

# Load data ---------------------------------------------------------------
d <- read.csv("data/full_data.csv")
t <- read.tree(file = "data/phylogentic_tree_s3_tre")
trees <- read.tree(file = "data/phylogentic_trees_300_s2.tre")
rownames(d) <- d$tip_name

# data including non-insect taxa
dfull <- read.csv("data/data_full_records.csv")
rownames(dfull) <- dfull$tip_name

# data excluding records on reproductive guilds.
drepr <- read.csv("data/data_exclude_repro.csv")
rownames(drepr) <- drepr$tip_name

# Pure regression---------------------------------------------------------------
# maleness ~ nspe---------------------------------------------------------------
m1 <- phylolm(maleness ~ log2(nspe+1), data = d, phy = t, model = "lambda")
summary(m1)

# 1. Sampling uncertainty-------------------------------------------------------
# 300 repetitions for each interval (5%, 10%, 20%, 30%, 40%, 50%)
samp <- samp_phylm(maleness ~ log2(nspe+1), data = d, phy = t, model = "lambda", 
                   n.sim = 1000, 
                   breaks = c(.05,.1,.2,.3,.4,.5))
summary(samp)
saveRDS(object = samp, file = "output/temp/sensitivity_sampling_uncertainty_estimates.RDs")

# Make raw data plots
# distribution of estimates
samp <- readRDS("output/temp/sensitivity_sampling_uncertainty_estimates.RDs")
gs1 <- ggplot(samp$sensi.estimates, aes(y = estimate, x = as.factor(n.percent))) +
  geom_jitter(width = .1, color = "gray") +
  geom_boxplot(fill = "white", alpha = .55, outlier.colour = NA) +
  geom_hline(yintercept = samp$full.model.estimates$coef[[2]], lty = 2, color = 'red') +
  theme_classic(base_size = 18) +
  labs(x = "Percentage of species removed (%)",
       y = "Estimated slope")
gs1  

# distribution of p.values
gs2 <- ggplot(samp$sensi.estimates, aes(y = pval.estimate , x = as.factor(n.percent))) +
  geom_jitter(width = .1, color = "gray") +
  geom_boxplot(fill = "white", alpha = .55, outlier.colour = NA) +
  geom_hline(yintercept = 0.05, lty = 2, color = "red") +
  scale_y_log10(breaks = c(0.000001,0.000001,0.00001, 0.0001, 0.001, 0.01, 0.1, 1)) +
  theme_classic(base_size = 18) +
  labs(x = "Percentage of species removed (%)",
       y = "P-value")
gs2  

gs <- gs1 + gs2 + plot_annotation(tag_levels = "a")
ggsave(gs, filename = "output/supp/SFigure_sensitivity_sampling.pdf", height =5, width = 11)

# save table
write_xlsx(round(samp$sign.analysis, digits = 6), path = "output/supp/STable_sensitivity_sampling.xlsx")

# 2. Taxonomic uncertainty-------------------------------------------------------
# 2.1 Orders removal----------
# 300 repetitions for each order removed
# orders: 
cla_ord <- clade_phylm(maleness ~ log2(nspe+1), data = d, phy = t, model = "lambda", 
                       clade.col = "order", n.sim = 1000, n.species = 10)
summary(cla_ord)

# save table
write_xlsx(x = summary(cla_ord)[[1]], path = "output/supp/STable_sensi_clade_ord.xlsx")

# 2.2 families removal-----------------
# 300 repetitions for each family removed
# families
cla_fam <- clade_phylm(maleness ~ log2(nspe+1), data = d, phy = t, model = "lambda", 
                       clade.col = "family", n.sim = 1000, n.species = 10)
summary(cla_fam)

sensi_plot(cla_fam, clade = "Lamiaceae") # Manual save png (1000 x 500)
write_xlsx(x = summary(cla_fam)[[1]], path = "output/supp/STable_sensi_clade_fam.xlsx")

# 3. Phylogenetic uncertainty-------------------------------------------------------
stree <- tree_phylm(maleness ~ log2(nspe+1), data = d, phy = trees, model = "lambda",
                    n.tree = 300) 
#saveRDS(stree, file = "output/temp/sensitivity_tree_uncertainty_estimates.RDs")
#stree <- readRDS("output/temp/sensitivity_tree_uncertainty_estimates.RDs")

range(stree$sensi.estimates$estimate)
summary(stree)

# slope variation due to tree uncertainty
range(stree$sensi.estimates$estimate)

# intercept variation due to tree uncertainty
range(stree$sensi.estimates$intercept)

# p-value (slope) variation due to tree uncertainty
range(stree$sensi.estimates$pval.estimate)

# save plots
# distribution of estimates
gtreea <- 
  ggplot(stree$sensi.estimates, aes(x = estimate)) +
  geom_histogram(color = "white", fill = gray(.6)) +
  geom_vline(xintercept = mean(stree$sensi.estimates$estimate), lty = 2, color = "red") +
  theme_classic(base_size = 14) +
  labs(y = "Frequency", x = "Slope") 
gtreea

# distribution of lambda
gtreeb <- 
  ggplot(stree$sensi.estimates, aes(x = optpar)) +
  geom_histogram(color = "white", fill = gray(.6)) +
  theme_classic(base_size = 14) +
  labs(y = "Frequency", x = "Estimated Lambda") 
gtreeb

# distribution of p-values
gtreec <- 
  ggplot(stree$sensi.estimates, aes(x = pval.estimate)) +
  geom_histogram(color = "white", fill = gray(.6)) +
  theme_classic(base_size = 14) +
  labs(y = "Frequency", x = "P-value")
gtreeb
gtree <- gtreeb + gtreeb + gtreec 
gtree

# Phylogenetic uncertainty on tree----------------------------------------------
guncert <- ggdensitree(trees[sample(x = 1:100, size = 100)], alpha = .3,
                       colour= "darkgreen", size = .5, layout = "circular") +
  theme(plot.margin=unit(c(0,0,0,0),"mm"))
ggsave(guncert, filename = "output/supp/SFigure_phylogenetic_tree_uncertainty.png",
       width = 10, height = 10)

guncert <- ggdraw() + draw_image("output/supp/SFigure_phylogenetic_tree_uncertainty.png")

gsensitree <- plot_grid(guncert, gtreea, gtreec, ncol = 1, labels = c("a", "b", "c"), 
                        label_fontface = "plain", label_size = 18)
ggsave(plot = gsensitree, filename = "output/supp/SFigure_sensitivity_tree_uncertainty.pdf", 
       height = 11, width = 5)


# 4. Feeding guilds-------------------------------------------------------------
# Excluding reproductive feeding guilds
mguil <- phylolm(maleness ~ log2(nspe+1), data = drepr, phy = t, model = "lambda", boot = 1000)
summary(mguil)

# save table
write_xlsx(tidy(mguil), path = "output/supp/STable_sensitivity_taxa_repro_guilds.xlsx")

# plot regression
bs <- 14
pcolor = "steelblue"
lcolor = "gray"
psize = 3
cc   <- coef(mguil)
ci   <- mguil$bootstrap
x    <- drepr$nspe +1
y    <- pred_log2x(a = cc[1], x = x, b = cc[2])
r2   <- round(mguil$r.squared, digits = 3)
pva  <- round(summary(mguil)[[2]][2,6], digits = 5)

pred1 <- data.frame(nspe = drepr$nspe, maleness = drepr$maleness, y, x)

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
       subtitle = "Excluding insects feeding on reproductive structures") +
  theme(plot.caption = element_markdown(),
        axis.text = element_text(size = 9),
        plot.subtitle = element_text(size = 10)) +
  annotate(geom = "text", x = 100, y = 0.1, label = paste("Y = ", round(cc[1], digits = 3),
                                                          " + ", round(cc[2], digits = 3), "X",
                                                          " | p = ", pva,
                                                          sep = ""), size = 3)

g1
saveRDS(g1, file = "output/temp/sensi_exclude_repro.Rds")

# 5. Non-insect enemies--------------------------------------------------------
# Including the full enemy matrix (nematodes, etc)
mall <- phylolm(maleness ~ log(nspe+1), data = dfull, phy = t, model = "lambda", boot = 1000)
summary(mall)

# save table
write_xlsx(tidy(mall), path = "output/supp/STable_sensi_taxa_all.xlsx")

# plot regression
cc   <- coef(mall)
ci   <- mall$bootstrap
x    <- dfull$nspe +1
y    <- pred_log2x(a = cc[1], x = x, b = cc[2])
r2   <- round(mall$r.squared, digits = 3)
pva  <- round(summary(mall)[[2]][2,6], digits = 5)

pred1 <- data.frame(nspe = dfull$nspe, maleness = dfull$maleness, y, x)

g2 <-
  ggplot(pred1, aes(y = maleness, x = nspe+1)) +
  scale_x_continuous(trans = "log2", 
                     breaks = c(1, 2, 4, 8, 16, 32, 64, 128, 256, 512),
                     limits = c(1, 512)) +
  scale_y_continuous(breaks = seq(0,.8,.1)) +
  geom_abline(intercept = ci[,1], slope = ci[,2], color = lcolor, alpha = .01, size = 1) +
  geom_point(fill = pcolor, color = "black", shape = 21, size = psize, alpha = .75) +
  geom_line(aes(y = y, x = nspe+1), size = 1.5) +
  theme_classic(base_size = bs) +
  labs(y = "Flower maleness", 
       x = "Number of insect species",
       subtitle = "Including non-insect invertebrate herbivores") +
  theme(plot.caption = element_markdown(),
        axis.text = element_text(size = 9),
        plot.subtitle = element_text(size = 10)) +
  annotate(geom = "text", x = 100, y = 0.1, label = paste("Y = ", round(cc[1], digits = 3),
                                                          " + ", round(cc[2], digits = 3), "X",
                                                          " | p = ", pva,
                                                          sep = ""), size = 3)
g2
saveRDS(g2, file = "output/temp/sensi_add_non_insect.Rds")


# Ploting all together
ggrepro <- readRDS("output/temp/sensi_exclude_repro.Rds")
ggpolli <- readRDS("output/temp/sensi_exclude_polli.Rds")
ggfull  <- readRDS("output/temp/sensi_add_non_insect.Rds")

gall <- ggrepro + ggpolli + ggfull + plot_annotation(tag_levels = "a")
gall
ggsave(gall, filename = "output/supp/SFigure_polli_repro_non_insect.png", 
       width = 12, height = 5)

# 6. Controlling for species distribution range---------------------------------
# 6.1 Using distribution as a covariate-----------------------------------------
mdist <- phylolm(maleness ~ log2(nspe+1) + range,
                 data = d, phy = t, model = "lambda", boot = 1000)
summary(mdist)
write_xlsx(tidy(mdist), path = "output/supp/STable_controlled_pgls_species_range.xlsx")

# fit with caper
mdist.2 <- pgls(maleness ~ log2(nspe+1) + range,
                data = comparative.data(t, d %>% dplyr::select(maleness, tip_name, nspe, range),
                                        tip_name, vcv=TRUE, vcv.dim=3),
                lambda='ML')
anova(mdist.2)
summary(mdist.2)

# save anova table
writexl::write_xlsx(x = tidy(anova(mdist.2)), 
                    path = "output/supp/STable_controlled_pgls_species_range_anova.xlsx")

# 6.2 Using the residuals-------------------------------------------------------
# Distribution of area of occupancy
ga1 <- 
  ggplot(d, aes(x = range)) +
  geom_histogram(color = "white", fill = gray(.4)) +
  theme_classic(base_size = 14) +
  labs(x = "Area of occupancy",
       y = "Frequency")
ga1

# regression between number of isect species and range
pcolor = "steelblue"
lcolor = "gray"
psize = 3

mod <- lm(log2(nspe+1) ~ (range), d)
summary(mod)

# plot regression 
cc   <- coef(mod)
pva  <- round(summary(mod)[[4]][8], digits = 5)
ga2 <- 
  ggplot(d, aes(y= nspe+1, x = (range))) +
  geom_point(fill = pcolor, color = "black", shape = 21, size = psize, alpha = .75) +
  scale_y_continuous(trans = "log2") +
  geom_smooth(method = lm, color = "black") +
  theme_classic(base_size = 14) +
  labs(y = "Number of insect species",
       x = "Area of occupancy") +
  annotate(geom = "text", x = 2000, y = 1.5, 
           label = paste("Y = ", round(cc[1], digits = 3),
                         " + ", round(cc[2], digits = 3), "X",
                         " | p < 0.00001",
                         sep = ""), size = 3)
ga2
mod_res <- data.frame(tip_name = names(residuals(mod)), nspe_res = residuals(mod))

d <- left_join(d, mod_res)
rownames(d) <- d$tip_name

mres <- phylolm(maleness ~ nspe_res, data = d, phy = t, model = "lambda", boot = 1000)
summary(mres)

# save table
write_xlsx(tidy(mres), path = "output/supp/STable_controlled_pgls_species_range_residuals.xlsx")

# plot regression 
cc   <- coef(mres)
ci   <- mres$bootstrap
x    <- d$nspe_res
y    <- pred_x(a = cc[1], x = x, b = cc[2])
r2   <- round(mres$r.squared, digits = 3)
pva  <- round(summary(mres)[[2]][2,6], digits = 5)

pred1 <- data.frame(nspe_res = d$nspe_res, maleness = d$maleness, y, x)

ga3 <-
  ggplot(pred1, aes(y = maleness, x = nspe_res)) +
  scale_y_continuous(breaks = seq(0,.8,.1)) +
  geom_abline(intercept = ci[,1], slope = ci[,2], color = lcolor, alpha = .01, size = 1) +
  geom_point(fill = pcolor, color = "black", shape = 21, size = psize, alpha = .75) +
  geom_line(aes(y = y, x = nspe_res), size = 1.5) +
  theme_classic(base_size = 14) +
  labs(y = "Flower maleness", 
       x = "Residual number of insect species",
       subtitle = "") +
  annotate(geom = "text", x = 1.5, y = 0.1, label = paste("Y = ", round(cc[1], digits = 3),
                                                          " + ", round(cc[2], digits = 3), "X",
                                                          " | p = ", pva,
                                                          sep = ""), size = 3)

grange <- ga1 / ga2 / ga3 + plot_annotation(tag_levels = "a") 
ggsave(grange, filename = "output/supp/SFigure_controlling_area_occupancy.pdf",
       width = 4.5, height = 10)

# 7. Integrating all predictors (PCA)---------------------------------
dpred <- d %>% 
  mutate(nspe = log2(nspe + 1),
         nfam = log2(nfam + 1),
         ngui = log2(nfam + 1),
         shan = shan)

dpred <- dpred[, c("nspe", "nfam", "ngui", "shan")]
pca_out <- prcomp(dpred, scale. = TRUE)

# PCA output
fviz_eig(pca_out)
fviz_pca_biplot(pca_out, col.ind = gray(.5), col.var = "tomato",
                geom.ind = "point") +
  theme_bw(base_size = 18) 
pca_out 

d$pc1 <- pca_out$x[,1]
d$pc2 <- pca_out$x[,2]

mpca1 <- phylolm(maleness ~ pc1, data = d, phy = t, model = "lambda", boot = 1000)
mpca2 <- phylolm(maleness ~ pc2, data = d, phy = t, model = "lambda", boot = 1000)
summary(mpca1)
summary(mpca2)

# save table
write_xlsx(tidy(mpca1), path = "output/supp/STable_controlled_pgls_pca1_predictors.xlsx")

# plot regression 
pcolor = "steelblue"
lcolor = "gray"
psize = 3
cc   <- coef(mpca1)
ci   <- mpca1$bootstrap
x    <- d$pc1
y    <- pred_x(a = cc[1], x = x, b = cc[2])
r2   <- round(mpca1$r.squared, digits = 3)
pva  <- round(summary(mpca1)[[2]][2,6], digits = 5)

pred1 <- data.frame(nspe_res = d$pc1, maleness = d$maleness, y, x)

g7 <-
  ggplot(pred1, aes(y = maleness, x = nspe_res)) +
  scale_y_continuous(breaks = seq(0,.8,.1)) +
  geom_abline(intercept = ci[,1], slope = ci[,2], color = lcolor, alpha = .01, size = 1) +
  geom_point(fill = pcolor, color = "black", shape = 21, size = psize, alpha = .75) +
  geom_line(aes(y = y, x = nspe_res), size = 1.5) +
  theme_classic(base_size = 18) +
  labs(y = "Flower maleness", 
       x = "Principal component (PC1)",
       subtitle = "") +
  annotate(geom = "text", x = 2, y = 0.1, label = paste("Y = ", round(cc[1], digits = 3),
                                                        " + ", round(cc[2], digits = 3), "X",
                                                        " | p = ", pva,
                                                        sep = ""), size = 4)
g7
ggsave(g7, filename = "output/supp/SFigure_sensi_pgls_pca_predictors.png", width = 6, height = 5)
