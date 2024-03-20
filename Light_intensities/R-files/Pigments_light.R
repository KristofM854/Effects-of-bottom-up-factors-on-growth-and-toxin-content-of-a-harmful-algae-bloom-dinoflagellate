##########################################
## Pigment composition and cell quotas of A. pseudogonyaulax strains L2-D2, L4-B1 and L4-B9 exposed to different light intensities
## Published in Limnology & Oceanography: "Effects of bottom-up factors on growth and toxin content of a harmful algae bloom dinoflagellate"
## All raw-data available on PANGAEA: https://doi.pangaea.de/10.1594/PANGAEA.965195 
## Questions to: kristof-moeller@outlook.de
## Kristof Möller 09/22
## Alfred-Wegener-Institute Bremerhaven
##########################################

# INSTALL AND LOAD PACKAGES ################################
# Use pacman to load add-on packages as desired
pacman::p_load(
  pacman,
  tidyverse,
  ggplot2,
  dplyr,
  DescTools,
  rstatix,
  broom,
  conover.test,
  flextable,
  officer
)

# Load and transform data ########
# Extract first two lines for the header (Pigment name and unit)
headers <-
  read.table(
    "C:\\Users\\krist\\OneDrive\\Dokumente\\AWI\\Promotion\\AP1\\Light-AP1b\\raw\\Pigments.txt",
    nrows = 2,
    header = FALSE
  )

# Combine both header lines in a single line
headers_names <- sapply(headers, paste, collapse = "_")

data <-
  read.csv(
    "C:\\Users\\krist\\OneDrive\\Dokumente\\AWI\\Promotion\\AP1\\Light-AP1b\\raw\\Pigments.txt",
    skip = 2,
    header = FALSE,
    sep = ""
  )

# Assign extracted header to the pigment dataframe
names(data) <- headers_names

# replace all , with . + make columns numeric
data[, 2:17] <-
  lapply(data[, 2:17], function(x)
    as.numeric(gsub(",", ".", x)))

# Introduce groups according to triplicate treatments
data$group <- rep(c(1:9), each = 3)

# Change ng/L of pigments to pg/cell and change colnames
data[, 9:17] <- (data[, 9:17] / data$Cell_count_Cells_L) * 10 ^ 3
beta <- intToUtf8(946) # beta symbol for plot
colnames(data)[9:17] <-
  c(
    "chl_c1/2",
    "peridinin",
    "viola",
    "diadino",
    "dino",
    "diato",
    "zea",
    "chl-a",
    paste(beta, "-carotene")
  )

# Introduce light harvesting LH and light protecting pigment LP sum
data <- data %>%
  mutate(
    LH = peridinin + `chl-a` + `chl_c1/2`,
    LP = diadino + dino + diato + viola + zea + `β -carotene`
  )

## Export pigment data for PANGAEA upload
# data %>%
#   select(treat, Day_adjusted, data1, data, replicate) %>%
#   write.table("C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/AP1/Light-AP1b/PANGAEA/cell_counts_L2D2.txt", sep = "\t", row.names = FALSE)
#

# Reshape data for statistical analysis ####

data1 <-
  data %>% dplyr::select(Sample_unit, all_of(colnames(data[, 9:17]))) %>%
  pivot_longer(cols = -Sample_unit,
               names_to = "pigment",
               values_to = "data") %>%
  mutate(
    group = as.factor(rep(1:9, each = 27)),
    strain = as.factor(str_sub(Sample_unit, 1, 4)),
    light = as.factor(str_extract(Sample_unit, "(?<=_)[0-9]+(?=_)"))
  ) %>%
  convert_as_factor(pigment)

# Aggregate mean of each treatment (group) and each pigment
agg <- aggregate(data ~ group + pigment, data = data1, FUN = "mean")

agg$grouping <- rep(c("L2-D2", "L4-B1", "L4-B9"), each = 3)

xlabels = rep(c("20", "100", "200"), length.out = 9) # x-labels for plot

# In-between strains statistics ######
# Perform Kruskal-Wallis test comparing pigment quotas between strains
Pig_all <- data1 %>%
  group_by(pigment, light) %>%
  nest() %>%
  mutate(kruskal_results = map(data, ~ kruskal_test(data ~ strain, data = .x))) %>%
  unnest(c(kruskal_results))

# Perform Kruskal-Wallis test comparing LH:LP ratios between strains
# extract only sample name and LH, LP values 
LH_LP <-
  data %>% dplyr::select(Sample_unit, LH, LP) %>% mutate(
    strain = as.factor(str_sub(Sample_unit, 1, 4)),
    light = as.factor(str_extract(Sample_unit, "(?<=_)[0-9]+(?=_)")),
    ratio = LH / LP
  )

LH_LP %>%
  group_by(light) %>%
  nest() %>%
  mutate(kruskal_results = map(data, ~ kruskal_test(ratio ~ strain, data = .x))) %>%
  unnest(c(kruskal_results))  %>%
  mutate(conover_results = map(data, ~ conover.test(
    .x$ratio, .x$strain, altp = T, method = "BH"
  ))) %>%
  unnest(c(conover_results))

# Calculate mean ratios of pigments +- SD and export as table
ratios <-
  read.csv(
    "C:\\Users\\krist\\OneDrive\\Dokumente\\AWI\\Promotion\\AP1\\Light-AP1b\\raw\\Pigments_ratio.txt",
    header = T,
    sep = ""
  )

# introduce replicate grouping variable (triplicates)
ratios$group <- rep(c(1:9), each = 3)

# replace all , with . + make columns numeric
ratios[, 2:9] <-
  sapply(ratios[, 2:9], function(x)
    as.numeric(gsub(",", ".", x)))

ratios_stat <- ratios %>% dplyr::select(Sample, everything()) %>%
  pivot_longer(cols = -c(Sample, group), names_to = "ratio", values_to = "data") %>%
  mutate(group = as.factor(str_extract(Sample, "(?<=_)[0-9]+(?=_)")),
         grouping = as.factor(str_sub(Sample, 1, 4)),
         group2 = paste(group, grouping, sep = "_"),
         group2 = paste(group2, ratio, sep = "_"))

# Statistical analysis of pigment ratios ######
# Comparison between pigment ratios at different light intensities of the same strain
# Perform Kruskal-Wallis test
for (i in unique(ratios_stat$ratio)) {
  print(i)
  for (strain in unique(ratios_stat$grouping)) {
    sub <-
      subset(ratios_stat,
             ratios_stat$ratio == i &
               ratios_stat$grouping == strain)
    
    print(strain)
    print(kruskal.test(sub$data ~ sub$group2))
  }
}

# Perform Conover-Iman posthoc test
Conover <- data.frame()
for (i in unique(ratios_stat$ratio)) {
  for (strain in unique(ratios_stat$grouping)) {
    sub <-
      subset(ratios_stat,
             ratios_stat$ratio == i &
               ratios_stat$grouping == strain)
    if (kruskal.test(sub$data ~ sub$group2)$'p.value' < 0.05) {
      Conover <-
        rbind(Conover, as.data.frame(conover.test(
          sub$data,
          factor(sub$group2),
          method = "bh",
          wrap = T
        )))
    }
  }
}

# Perform Kruskal-Wallis and Conover-Iman posthoc test in the same tidyverse call 
LH_LP_stat <- LH_LP %>%
  group_by(strain) %>%
  nest() %>%
  mutate(kruskal_results = map(data, ~ kruskal_test(ratio ~ light, data = .x))) %>%
  unnest(c(kruskal_results)) %>%
  mutate(conover_results = map(data, ~ conover.test(
    .x$ratio, .x$light, altp = T, method = "BH"
  ))) %>%
  unnest(c(conover_results))

Conover$strain <- sub("^(\\d+)_([^_]+)_.+", "\\2", Conover$comparisons)
Conover$ratio <- sub("^[^_]+_[^_]+_(.*?) - .*", "\\1", Conover$comparisons)

# # aggregate mean and standard deviation of each pigment ratio of each strain of each light intensity and export as table ######
# ratio_agg <-
#   aggregate(data = ratios_stat, data ~ group + grouping + ratio, FUN = "mean")
# 
# ratio_agg$sd <-
#   aggregate(data = ratios_stat, data ~ group + grouping + ratio, FUN = "sd")$data
# 
# # combine mean and SD in one expression for table in manuscript
# ratio_agg$mean_sd <-
#   paste(format(round(ratio_agg$data, 2), nsmall = 2), format(round(ratio_agg$sd, 2), nsmall =
#                                                                2), sep = " +/- ")
# # pivot wider for easier transfer to manuscript
# ratio_agg_pivot <-
#   pivot_wider(ratio_agg[, c(1:3, 6)], names_from = ratio, values_from = mean_sd)
# 
# ratio_agg_pivot <- flextable(ratio_agg_pivot) %>% autofit %>%
#   save_as_docx(path = "C:\\Users\\krist\\OneDrive\\Dokumente\\AWI\\Promotion\\AP1\\Light-AP1b\\ratio_agg_pivot.docx")

# Garbage collection: call after large objects have been removed
gc()

# Delete workspace, clean environment/variables and restarte R if used memory piles up
dev.off()

rm(list = ls())

.rs.restartR()