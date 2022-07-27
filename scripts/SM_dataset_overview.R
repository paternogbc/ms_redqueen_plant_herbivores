# Code from the manuscript: Insect herbivores drive sex allocation in angiosperm flowers

# Script to check and illustrate the distribution of variables
# Last update: 2022.07.27
# Author: Gustavo brant Paterno (paternogbc@gmail.com)

# Load packages ----------------------------------------------------------------
library(tidyverse) 
library(phylolm)   
library(patchwork) 
library(broom)     
library(writexl)   
library(ggtext)    
library(ggExtra)   
library(tidyverse) 
library(sensiPhy)  
library(phytools)  
library(ggtree)    
library(cowplot)   
library(GGally)

# source local functions--------------------------------------------------------
source("scripts/zzz_functions.R")
set.seed(1234522)

# Load data ---------------------------------------------------------------
d <- read.csv("data/full_data.csv")
tax <- read.csv("data/study_species_list_taxonomy.csv")
trees <- read.tree("data/phylogentic_trees_300_s2.tre")
t <- read.tree(file = "data/phylogentic_tree_s3_tre")
rownames(d) <- d$tip_name

# Taxonomic and phylogenetic distribution---------------------------------------
# Phylogenetic tree-------------------------------------------------------------
rectree <- 
  ggtree(t, ladderize = TRUE, color = gray(.4), size = .5) +
  geom_tiplab(size = 1.5) +
  xlim(0, 200)
rectree

# Distribution of orders
tax %>% 
  group_by(order) %>% 
  tally() %>% 
  arrange(-n) %>% 
  ggplot(aes(x = reorder(order, n), y = n, label = n)) +
  geom_col(fill = "darkgreen", alpha = .5,  width = .75) +
  coord_flip() +
  scale_y_continuous(breaks = seq(0,35,5)) +
  theme_classic(base_size = 12) +
  labs(y = "Number of species", x = "",
       subtitle = "Distribution of orders") -> g1
g1

# Distribution of familes
tax %>% 
  group_by(family) %>% 
  tally() %>% 
  arrange(-n) %>% 
  ggplot(aes(x = reorder(family, n), y = n, label = n)) +
  geom_col(fill = "darkgreen", alpha = .5, width = .75) +
  scale_y_continuous(breaks = seq(0,35,2)) +
  coord_flip() +
  theme_classic(base_size = 12) +
  labs(y = "Number of species", x = "",
       subtitle = "Distribution of families") -> g2
g2

n_distinct(d$order)
n_distinct(d$family)

gtaxtree <- rectree + (g1 / g2) + plot_annotation(tag_levels = "a")
ggsave(plot = gtaxtree, 
       filename = "output/supp/SFigure_phylogenetic_tree.pdf",
       width = 9.5, height = 11)
ggsave(plot = gtaxtree, 
       filename = "output/supp/SFigure_phylogenetic_tree.png",
       width = 9.5, height = 11)

# Plant Variables distribution--------------------------------------------------------
range(d$maleness)
mean(d$maleness)
sd(d$maleness)
gmale <-
  ggplot(d, aes(x = maleness, fill = maleness)) +
  geom_histogram(fill = "darkgreen", alpha = .5, color = "white", bins = 30) +
  scale_x_continuous(breaks = seq(0,1,.1)) +
  scale_y_continuous(breaks = seq(0,12,2), limits = c(0, 12)) +
  theme_classic(base_size = 12) +
  labs(y= "Frequency",
       x = "Flower maleness") +
  theme(axis.title = element_text(size = 16))
gmale

# Vegetative traits
gsla <-
  ggplot(d, aes(x = sla_imp)) +
  geom_histogram(fill = "darkgreen", alpha = .5, color = "white", bins = 20) +
  scale_x_continuous(trans = "log10" )+
  #scale_y_continuous(breaks = seq(0,10,2)) +
  theme_classic(base_size = 12) +
  labs(y= "Frequency",
       x = expression("Specific leaf area" ~ (mm^2 * mg^-1))) +
  theme(axis.title = element_text(size = 16))

gsla

# Height
ghei <-
  ggplot(d, aes(x = height)) +
  geom_histogram(fill = "darkgreen", alpha = .5, color = "white", bins = 20) +
  scale_x_continuous(trans = "log10" )+
  #scale_y_continuous(breaks = seq(0,10,2)) +
  theme_classic(base_size = 12) +
  labs(y= "Frequency",
       x = "Maximum plant height (m)") +
  theme(axis.title = element_text(size = 16))
ghei

# Reproductive traits
# pollinators
gpoll <-
  ggplot(d, aes(x = polli)) +
  geom_bar(fill = "darkgreen", alpha = .5, width = .75) +
  scale_x_discrete(limits = c("A", "B", "F", "H", "P", "other", NA))+
  #scale_y_continuous(breaks = seq(0,10,2)) +
  theme_classic(base_size = 12) +
  labs(y= "Frequency",
       x = "Pollinator guild") +
  theme(axis.title = element_text(size = 16))

gpoll

# flower shape
gshape <-
  ggplot(d, aes(x = shape)) +
  geom_bar(fill = "darkgreen", alpha = .5, width = .75) +
  scale_y_continuous(breaks = seq(0,60, 20), limits = c(0, 60)) +
  theme_classic(base_size = 12) +
  labs(y= "Frequency",
       x = "Flower shape") +
  theme(axis.title = element_text(size = 16))

gshape

# flower color
gcolor <-
  ggplot(d, aes(x = color)) +
  geom_bar(fill = "darkgreen", alpha = .5, width = .75) +
  #scale_y_continuous(breaks = seq(0,10,2)) +
  theme_classic(base_size = 12) +
  labs(y= "Frequency",
       x = "Flower colour")+
  theme(axis.title = element_text(size = 16))

gcolor

# Nectar offer
gnectar <-
  ggplot(d, aes(x = nectar)) +
  geom_bar(fill = "darkgreen", alpha = .5, width = .75) +
  scale_y_continuous(breaks = seq(0,100,25), limits = c(0, 100)) +
  theme_classic(base_size = 12) +
  labs(y= "Frequency",
       x = "Nectar amount") +
  theme(axis.title = element_text(size = 16))
gnectar

gvar <- (gmale + (ghei / gsla )) / (gpoll + gshape + gcolor + gnectar) + 
  plot_annotation(tag_levels = "a") 
gvar
ggsave(gvar, filename = "output/supp/SFigure_plant_variables.pdf", width = 8.5, height = 10)
ggsave(gvar, filename = "output/supp/SFigure_plant_variables.png", width = 8.5, height = 10)

# Insect Variables distribution--------------------------------------------------------
# Correlation between predictors------------------------------------------------
dcorr <- select(d, nspe, nfam, ngui, shan)
dcorr %>% 
  mutate(
    nspecies  = log2(nspe+1),
    nfamilies =  log2(nfam+1),
    nguilds   = log2(ngui+1),
  ) %>% 
  select(-nspe, -nfam, -ngui) %>% 
  select(nspecies, nfamilies, nguilds, divguilds = shan) -> dcorr

# save correlation matrix
write.csv(x = cor(dcorr), file =  "output/supp/STable_correlation_matrix.csv")

# correlation matrix 
gcor <- 
  ggpairs(dcorr,
          lower = list(continuous = wrap("smooth", fill = "steelblue", method = "loess", color = "tomato", alpha = .5)),
          diag  = list(continuous = wrap("densityDiag", alpha=0.25, fill = "orange")),
          upper = list(continuous = wrap("cor", size = 9))) +
  theme_classic(base_size = 18)
gcor

ggsave(gcor, filename = "output/supp/SFigure_correlation_alternative_herbivore_predictors.png",
       width = 9, height = 8)
saveRDS(object = gcor, file = "output/temp/Figure_correlation_matrix_insect_predictors.RDs")

# Numbe rof inscet species
d %>% arrange(-nspe) %>% 
  select(nspe, nfam, ngui, shan) 
range(d$nspe)
mean(d$nspe)
median(d$nspe)

gnsp <-
  ggplot(d, aes(x = log2(nspe+1))) +
  geom_histogram(fill = "orange", alpha = .7, color = "white", bins = 10) +
  scale_x_continuous(breaks = c(0,2,4,6,8), 
                     labels = 2^c(0,2,4,6,8)) +
  theme_classic(base_size = 12) +
  labs(y= "Frequency",
       x = expression("Number of insect species"))
gnsp

# Number of insect families
range(d$nfam)
mean(d$nfam)
median(d$nfam)

gfam <-
  ggplot(d, aes(x = log2(nfam + 1))) +
  geom_histogram(fill = "orange", alpha = .7, color = "white", bins = 10) +
  scale_x_continuous(breaks = c(0,2,4,6), 
                     labels = 2^c(0,2,4,6)) +
  theme_classic(base_size = 12) +
  labs(y= "Frequency",
       x = expression("Number of insect families"))
gfam

# Number of feeding guilds
range(d$ngui)
mean(d$nfam)
median(d$ngui)
ggui <-
  ggplot(d, aes(x = log2(ngui + 1))) +
  geom_histogram(fill = "orange", alpha = .7, color = "white", bins = 10) +
  scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5),
                     labels = 2^c(0, 1, 2, 3, 4, 5)) + 
  theme_classic(base_size = 12) +
  labs(y= "Frequency",
       x = expression("Number of feeding guilds"))
ggui

# Diversity of feeding guilds
range(d$shan)
median(d$shan)

gsha <-
  ggplot(d, aes(x = shan )) +
  geom_histogram(fill = "orange", alpha = .7, color = "white", bins = 10) +
  theme_classic(base_size = 12) +
  labs(y= "Frequency",
       x = expression("Diversity of feeding guilds"))
gsha

gins <- gnsp + gfam + ggui + gsha

# load correaltion matrix
gcor <- ggdraw() + draw_image("output/supp/SFigure_correlation_alternative_herbivore_predictors.png")

ginsout <-  (gnsp + gfam + ggui + gsha) / gcor + plot_annotation(tag_levels = "a")

ggsave(ginsout, filename = "output/supp/SFigure_insect_variables.pdf",
       width = 6, height = 10)
ggsave(ginsout, filename = "output/supp/SFigure_insect_variables.png",
       width = 6, height = 10)
