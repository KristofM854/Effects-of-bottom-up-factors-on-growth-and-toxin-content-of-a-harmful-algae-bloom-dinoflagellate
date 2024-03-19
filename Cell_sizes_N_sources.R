##########################################
## AP1 N-growth Experiment; statistics of A pseudogonyaulax cell sizes of strains L2-D2 and L4-B1
## Published in Limnology & Oceanography: "Effects of bottom-up factors on growth and toxin content of a harmful algae bloom dinoflagellate"
## All raw-data available on PANGAEA: https://doi.pangaea.de/10.1594/PANGAEA.965195 
## Questions to: kristof-moeller@outlook.de
## Kristof Möller 09/22
## Alfred-Wegener-Institute Bremerhaven
##########################################

# Rows 1-2 exponential/stationary C; 3-4 NO3, 5-6 NH4
# Column 1: Treatment; Column 2: Mean; Column 3: SD
# Date: 25.11.21

# INSTALL AND LOAD PACKAGES ################################
library(datasets)  # Load base packages manually

# Installs pacman ("package manager") if needed
if (!require("pacman"))
  install.packages("pacman")

# Use pacman to load add-on packages as desired
pacman::p_load(pacman,
               psych,
               rio,
               tidyverse,
               data.table,
               flextable,
               officer,
               ggpubr,
               rstatix,
               conover.test)

# LOAD DATA strain L2-D2 ################################################

Stat_2D2 = read.csv(
  "C:\\Users\\krist\\OneDrive\\Dokumente\\AWI\\Promotion\\AP1\\N-AP1a\\R_txt_files\\cell_size_L2D2_all.txt",
  sep = "",
  header = T
)

# General data transformations
Stat_2D2 <- Stat_2D2 %>% t() %>% data.frame()

Stat_2D2 <-
  Stat_2D2 %>% dplyr::rename("ex" = "X1", "st" = "X2") %>% # ex = exponential and st = stationary phase
  mutate(group = as.factor(rep(c( # introduce treatments
    "C", "NO3", "NH4", "Urea"
  ), each = 3))) %>%
  pivot_longer(cols = c(1:2)) %>%
  filter(group != "Urea" | name != "ex") %>%
  mutate(treat2 = paste(group, name, sep = "_")) %>%
  dplyr::rename(treat = group, gp = name, data = value) %>%
  convert_as_factor(treat, treat2, gp)

# Statistics L2-D2 ################################################
# Perform Kruskal-Wallis Test
Krus1 <- kruskal_test(data = Stat_2D2, data ~ treat2)

# Calculate Cohen's effect size
Stat_2D2 %>% group_by(gp) %>% cohens_d(data ~ treat)

# Perform Conover-Iman Post Hoc test
Con_2D2 <-
  conover.test(Stat_2D2$data,
               Stat_2D2$treat2,
               method = "BH",
               altp = T)

# LOAD DATA strain L4-B1 ################################################

Stat_4B1 = read.csv(
  "C:\\Users\\krist\\OneDrive\\Dokumente\\AWI\\Promotion\\AP1\\N-AP1a\\R_txt_files\\cell_size_L4B1_all.txt",
  sep = "",
  header = T
)

# General data transformations
Stat_4B1 <- Stat_4B1 %>% t() %>% data.frame()

Stat_4B1 <-
  Stat_4B1 %>% dplyr::rename("ex" = "X1", "st" = "X2") %>%
  mutate(group = as.factor(rep(c(
    "C", "NO3", "NH4", "Urea"
  ), each = 3))) %>%
  pivot_longer(cols = c(1:2)) %>%
  filter(group != "Urea" | name != "ex") %>%
  #filter(group != "NH4" | name != "st") %>%
  mutate(treat2 = paste(group, name, sep = "_")) %>%
  dplyr::rename(treat = group, gp = name, data = value) %>%
  convert_as_factor(treat, treat2, gp)

# Statistics L4-B1 ################################################
# Perform Kruskal-Wallis Test
Stat_4B1 %>% group_by(gp) %>% kruskal_test(data ~ treat)
Stat_4B1 %>% kruskal_test(data ~ treat2)

# Calculate Cohen's effect size
Stat_4B1 %>% group_by(gp) %>% cohens_d(data ~ treat)

# Perform Conover Iman Post Hoc test
Con_4B1 <-
  conover.test(Stat_4B1$data,
               Stat_4B1$treat2,
               method = "BH",
               altp = T)

# Between-strains statistics ################################################

# Combine both strains in one data set
Stat_2D2 <- Stat_2D2 %>% mutate(strain = "L2-D2")
Stat_4B1 <- Stat_4B1 %>% mutate(strain = "L4-B1")

Stat_all <- rbind(Stat_2D2, Stat_4B1) %>% convert_as_factor(treat, strain, treat2)

# Perform Kruskal-Wallis Test for each treatment in between strains
cell_size_all <- Stat_all %>%
  group_by(treat2) %>%
  nest() %>%
  mutate(kruskal_results = map(data, ~kruskal_test(data ~ strain, data = .x))) %>%
  unnest(c(kruskal_results)) 

# Garbage collection: call after large objects have been removed  ################################

gc()

# Delete workspace, clean environment/variables and restarte R if used memory piles up

dev.off()

rm(list = ls())

.rs.restartR()
