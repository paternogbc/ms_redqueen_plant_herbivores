# Code from the manuscript: Insect herbivores drive sex allocation in angiosperm flowers

# Script to plot main figures: the phylogenetic tree and PGLS regressions
# Last update: 2022.07.27
# Author: Gustavo Brant Paterno (paternogbc@gmail.com)

# Load packages -----------------------------------------------------------
library(tidyverse)  # CRAN v1.3.0
library(phylolm)    # CRAN v2.6.2
library(ggtext)     # [github::wilkelab/ggtext] v0.1.0.9000
library(ggtree)     # Bioconductor v2.2.1
library(patchwork)  # CRAN v1.0.1
library(ggnewscale) # CRAN v0.4.2
library(phytools)   # CRAN v0.7-47
library(cowplot)    # CRAN v1.0.0

# source local functions--------------------------------------------------------
source("scripts/zzz_functions.R")

# Load data ---------------------------------------------------------------
d  <- read.csv("data/full_data.csv") 
m1 <- readRDS("output/temp/pgls_maleness_vs_nspe.RDs")
m2 <- readRDS("output/temp/pgls_maleness_vs_nfam.RDs")
m4 <- readRDS("output/temp/pgls_maleness_vs_ngui.RDs")
m5 <- readRDS("output/temp/pgls_maleness_vs_shan.RDs")

tax <- read.csv("data/study_species_list_taxonomy.csv")
trees <- read.tree("data/phylogentic_trees_300_s2.tre")
t <- read.tree(file = "data/phylogentic_tree_s3_tre")
rownames(d) <- d$tip_name

# PGLS regressions--------------------------------------------------------------

# Plot parameters
bs <- 12
pcolor = "steelblue"
psize = 2.5
palpha = .5
lcolor = gray(.5)
lsize  = .25

# 1. Number of insect species------------------------------------------------------
cc   <- coef(m1)
ci   <- m1$bootstrap
x    <- d$nspe +1
y    <- pred_log2x(a = cc[1], x = x, b = cc[2])
r2   <- round(m1$r.squared, digits = 3)
pva  <- round(summary(m1)[[2]][2,6], digits = 5)

pred1 <- data.frame(nspe = d$nspe, maleness = d$maleness, y, x)

g1 <-
  ggplot(pred1, aes(y = maleness, x = nspe+1)) +
  scale_x_continuous(trans = "log2", 
                     breaks = c(1, 2, 4, 8, 16, 32, 64, 128, 256, 512),
                     limits = c(1, 512)) +
  scale_y_continuous(breaks = seq(0,.8,.1)) +
  geom_abline(intercept = ci[,1], slope = ci[,2], color = lcolor, alpha = .01, size = lsize) +
  geom_point(fill = pcolor, color = "black", shape = 21, size = psize, alpha = palpha) +
  geom_line(aes(y = y, x = nspe+1), size = 1.5) +
  theme_classic(base_size = bs) +
  labs(y = "Flower maleness", x = "Number of insect species") +
  theme(plot.caption = element_markdown(),
        axis.text = element_text(size = 9)) +
  annotate(geom = "text", x = 0, hjust = -.2,
           y = 0, label = paste("Y = ", round(cc[2], digits = 3), "X",
                                " + ",  round(cc[1], digits = 3),
                                " | p = ", pva,
                                sep = ""), size = 3) 

# 2. Number of insect families------------------------------------------------------
cc   <- coef(m2)
ci   <- m2$bootstrap
x    <- d$nfam +1
y    <- pred_log2x(a = cc[1], b = cc[2], x = x)
pva  <- round(summary(m2)[[2]][2,6], digits = 5)

pred2 <- data.frame(nfam = d$nfam, maleness = d$maleness, y, x)

g2 <-
  ggplot(pred2, aes(y = maleness, x = nfam+1)) +
  scale_x_continuous(trans = "log2", breaks = c(1, 2, 4, 8, 16, 32, 64), limits = c(1,64)) +
  scale_y_continuous(breaks = seq(0,.8,.1)) +
  geom_abline(intercept = ci[,1], slope = ci[,2], color = lcolor, alpha = .01, size = lsize) +
  geom_point(fill = pcolor, color = "black", shape = 21, size = psize, alpha = palpha) +
  geom_line(aes(y = y, x = nfam + 1), size = 1.5) +
  theme_classic(base_size = bs) +
  labs(y = "Flower maleness", x = "Number of insect families") +
  theme(plot.caption = element_markdown(),
        axis.text = element_text(size = 9)) +
  annotate(geom = "text", x = 0, hjust = -.2,
           y = 0, label = paste("Y = ", round(cc[2], digits = 3), "X",
                                " + ",  round(cc[1], digits = 3),
                                " | p = ", pva,
                                sep = ""), size = 3) 
g2

# 3, Number of feeding guilds------------------------------------------------------
cc   <- coef(m4)
ci   <- m4$bootstrap
x    <- d$ngui + 1
y    <- pred_log2x(a = cc[1], b = cc[2], x = x)
pva  <- round(summary(m4)[[2]][2,6], digits = 5)

pred3 <- data.frame(ngui = d$ngui, maleness = d$maleness, y, x)

g3 <-
  ggplot(pred3, aes(y = maleness, x = ngui + 1)) +
  scale_x_continuous(trans = "log2", breaks = c(1, 2, 4, 8, 16, 32), limits = c(1,32)) +
  scale_y_continuous(breaks = seq(0,.8,.1)) +
  geom_abline(intercept = ci[,1], slope = ci[,2], color = lcolor, alpha = .01, size = lsize) +
  geom_point(fill = pcolor, color = "black", shape = 21, size = psize, alpha = palpha) +
  geom_line(aes(y = y, x = ngui + 1), size = 1.5) +
  theme_classic(base_size = bs) +
  labs(y = "Flower maleness", x = "Number of feeding guilds") +
  theme(plot.caption = element_markdown(),
        axis.text = element_text(size = 9)) +
  annotate(geom = "text", x = 0, hjust = -.2,
           y = 0, label = paste("Y = ", round(cc[2], digits = 3), "X",
                                " + ",  round(cc[1], digits = 3),
                                " | p = ", pva,
                                sep = ""), size = 3) 
g3

# 4. Shannon------------------------------------------------------
cc   <- coef(m5)
ci   <- m5$bootstrap
x    <- d$shan
y    <- pred_x(a = cc[1], b = cc[2], x = x)
pva  <- round(summary(m5)[[2]][2,6], digits = 5)

pred4 <- data.frame(shan = d$shan, maleness = d$maleness, y, x)

g4 <-
  ggplot(pred4, aes(y = maleness, x = shan)) +
  scale_x_continuous(breaks = seq(0,3,.5), limits = c(0,2.8)) +
  scale_y_continuous(breaks = seq(0,.8,.1)) +
  geom_abline(intercept = ci[,1], slope = ci[,2], color = lcolor, alpha = .01, size = lsize) +
  geom_point(fill = pcolor, color = "black", shape = 21, size = psize, alpha = palpha) +
  geom_line(aes(y = y, x = shan), size = 1.5) +
  theme_classic(base_size = bs) +
  labs(y = "Flower maleness", x = "Diversity of feeding guilds") +
  theme(plot.caption = element_markdown(),
        axis.text = element_text(size = 9)) +
  annotate(geom = "text", x = 0, hjust = -.2,
           y = 0, label = paste("Y = ", round(cc[2], digits = 3), "X",
                                " + ",  round(cc[1], digits = 3),
                                " | p = ", pva,
                                sep = ""), size = 3) 
g4

reg_plots <- list(g1, g2, g3, g4)
saveRDS(object = reg_plots, file = "output/temp/main_regression_plots_revision.RDs")

# 5. Phylogentic tree--------------------------------------------------------------
circ <- ggtree(t, layout = "fan", open.angle = 25, color = gray(.5), size = .25)
circ

d1 <- d %>% dplyr::select(Maleness = maleness)
d2 <- d %>% dplyr::mutate(
  "N insects" = log2(nspe + 1)
) %>% 
  dplyr::select(`N insects`) %>% 
  as.data.frame()
rownames(d2) <- d$tip_name

# colls <- c('#9e0142','#d53e4f','#f46d43','#fdae61','#fee08b','#ffffbf','#e6f598','#abdda4','#66c2a5','#3288bd','#5e4fa2')
# colls <- c('#d73027','#f46d43','#fdae61','#fee090','white','#e0f3f8','#abd9e9','#74add1','#4575b4')
colls <- c('#67001f','#b2182b','#d6604d','#f4a582','#fddbc7','#f7f7f7','#d1e5f0','#92c5de','#4393c3','#2166ac','#053061')

colls <- colls[10:2]
p1 <- 
  gheatmap(circ, d1, offset=.5, width =.1, 
           colnames = T,font.size = 2,
           colnames_angle = 80,
           colnames_offset_y = 5,
           colnames_position = "top") +
  scale_fill_gradientn(colours = colls,
                       breaks = seq(0,.8,.2), limits = c(0,.8), 
                       name = "Flower\nMaleness")
p1
p2 <- p1 + new_scale_fill()
p2
pout <- 
  gheatmap(p2, d2, colnames = T, font.size = 2,
           colnames_position = "top",
           colnames_angle = 80,
           colnames_offset_y = 5,
           offset=15, width = .1) +
  scale_fill_gradientn(colours = colls,
                       breaks = seq(0,8,2),
                       name = "N insects\n(log2)")
pout <- 
  pout + theme(plot.margin=unit(c(0,0,0,0),"mm")) +
  theme(
    # legend.position = c(-0,.5),
    legend.key.width = unit(.3, "cm"),
    legend.key.height = unit(.4, "cm"),
    legend.title =  element_text(size = 7),
    legend.text =  element_text(size = 6))

# Add clade labels-
### find species by Order
sp.ord <-
  split(as.character(d$tip_name),
        as.character(d$order))
ord.na <- names(sp.ord)

d %>% group_by(order) %>% 
  tally(sort = T) %>% 
  slice(1:10) %>% 
  pull(order) -> top_ord

### Find RMCA for each order:
ord.na <- top_ord
nods <- list()

for (i in ord.na) {
  nod <- findMRCA(tree = t, sp.ord[[i]])
  if (is.null(nod)) {
    nod <- NA
  }
  nods[i] <- nod
}

pout + 
  geom_cladelabel(node = nods[[1]], label = names(nods[1]), angle = 23, color = gray(.3), 
                  offset=38, offset.text = 9, hjust = .5, fontsize = 2.7) +
  geom_cladelabel(node = nods[[2]], label = names(nods[2]), angle = 59, color = gray(.3), 
                  offset=38, offset.text = 9, hjust = .5, fontsize = 2.7) +
  geom_cladelabel(node = nods[[3]], label = names(nods[3]), angle = -59, color = gray(.3), 
                  offset=38, offset.text = 9, hjust = .5, fontsize = 2.7) +
  geom_cladelabel(node = nods[[4]], label = names(nods[4]), angle = 98, color = gray(.3), 
                  offset=38, offset.text = 9, hjust = .5, fontsize = 2.7) +
  geom_cladelabel(node = nods[[5]], label = names(nods[5]), angle = 29, color = gray(.3), 
                  offset=38, offset.text = 9, hjust = .5, fontsize = 2.7) +
  geom_cladelabel(node = nods[[6]], label = names(nods[6]), angle = -32, color = gray(.3),
                  offset=38, offset.text = 9, hjust = .5, fontsize = 2.7) +
  geom_cladelabel(node = nods[[7]], label = names(nods[7]), angle = 0, color = gray(.3), 
                  offset=38, offset.text = 9, hjust = .5, fontsize = 2.7) -> g5

# # Figure 1----------------------------------------------------------------------
# Figure 1
fig1 <- g1 + g5 + plot_annotation(tag_levels = "a")

# save plot
ggsave(filename = "output/Figure_1.pdf",  device = cairo_pdf,
       plot = fig1, width = 173, height = 87, units = "mm")

# Figure 2 ----------------------------------------------------------------
fig2 <- g2 / g3 / g4 + plot_annotation(tag_levels = "a") 
fig2 + theme(axis.text = element_text(size = 10))

ggsave(fig2, filename = 'output/Figure_2.pdf', device = cairo_pdf,
       width = 82, height = 246, units = "mm")

