##########################################
## Analysis of dark parameters of the FRRF measurements of A. pseudogonyaulax strains L2-D2, L4-B1 and L4-B9 exposed to different light intensities
## Published in Limnology & Oceanography: "Effects of bottom-up factors on growth and toxin content of a harmful algae bloom dinoflagellate"
## All raw-data available on PANGAEA: https://doi.pangaea.de/10.1594/PANGAEA.965195 
## Questions to: kristof-moeller@outlook.de
## Kristof Möller 09/22
## Alfred-Wegener-Institute Bremerhaven
##########################################

# INSTALL AND LOAD PACKAGES ################################
library(datasets, ggplot2, easyGgplot2)  # Load base packages manually

# Installs pacman ("package manager") if needed
if (!require("pacman"))
  install.packages("pacman")

# Use pacman to load add-on packages as desired
pacman::p_load(
  tidyverse,
  dplyr,
  DescTools,
  rstatix,
  conover.test,
  flextable,
  officer
)

# Load Windows Fonts and define Times font
loadfonts(device = "win")
windowsFonts(Times = windowsFont("Times"))

# Load and transform data ########
data <-
  read.csv(
    "C:\\Users\\krist\\OneDrive\\Dokumente\\AWI\\Promotion\\AP1\\Light-AP1b\\raw\\Dark_parameters.txt",
    header = TRUE,
    sep = ""
  )

# replace all , with . + make columns numeric
data[, 2:15] <-
  lapply(data[, 2:15], function(x)
    as.numeric(gsub(",", ".", x)))

# introduce strain grouping factor
data$group <-
  as.factor(rep(c("L2-D2", "L4-B1", "L4-B9"), each = length(data$Fm...) /
                  3))

# introduce light intensity grouping factor
data$grouping <-
  as.factor(rep(c("20", "100", "200"), each = length(data$Fm...) / (3 * 3)))

# reorder factor levels of light intensity grouping factor
data$grouping <- factor(data$grouping, levels = c("20", "100", "200"))

# extract names of dark parameters
names <- colnames(data[, c(2:7, 10, 12:14)])

# Reshape data
data1 <- data %>%
  dplyr::select(all_of(names), group, grouping) %>%  # Select relevant columns including group and grouping
  pivot_longer(
    cols = -c(group, grouping),  # Exclude group and grouping from reshaping
    names_to = "parameter", 
    values_to = "data"
  )  # Reshape into long format

# Aggregate data ########################################

# Aggregate mean and standard deviation
agg <-
  aggregate(data ~ parameter + group + grouping,
            data = data1,
            FUN = "mean")
agg$sd <-
  aggregate(data ~ group + parameter + grouping,
            data = data1,
            FUN = "sd")$data

# combine mean and sd for export
agg$total <-
  paste(formatC(agg$data, digits = 2, format = "f"),
        formatC(agg$sd, digits = 2, format = "f"),
        sep = "+/-")

# pivot wider for export
agg_wider <-
  pivot_wider(agg[, c(1:3, 6)], names_from = parameter, values_from = total)
agg_wider <- agg_wider %>% arrange(group)

# export dark parameters
# agg_wider <- flextable(agg_wider) %>% autofit %>%
#   save_as_docx(path="C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/AP1/Light-AP1b/Dark_parameters.docx")

# STATISTICS ############################################

data1$parameter <- as.factor(data1$parameter)

data1 <- subset(data1, parameter != "Fm..." & parameter != "Fv.Fq..")

# Group-wise statistical tests (Kruskal Wallis Test, followed by Conover Iman posthoc test)
dark_stats <- data1 %>%
  group_by(parameter, group) %>%
  nest() %>%
  mutate(kruskal_results = map(data, ~ kruskal_test(data ~ grouping, data = .x))) %>%
  unnest(c(kruskal_results))  %>%
  filter(p <= 0.1) %>%
  mutate(conover_results = map(data, ~ conover.test(
    .x$data, .x$grouping, altp = T, method = "BH"
  ))) %>%
  unnest(c(conover_results))

# extract group comparisons 
Con_groups <- data1 %>%
  group_split(parameter, group)

# perform conover iman test alone 
Con <- data1 %>%
  group_split(parameter, group) %>%
  map( ~ .x %>%
         ConoverTest(
           data ~ grouping,
           method = "BH",
           data = .,
           out.list = T
         ))

# extract and reorder test results 
x <- c(1:3)
Con2 <- as.data.frame(sapply(Con, function(x) {
  do.call(rbind, x)
}))


Con2 <- as.data.frame(Con2[5:7, ])

Con3 <-
  data.frame(
    data = unlist(Con2),
    sample = rep(c("L2-D2", "L4-B1", "L4-B9"), each = 3),
    parameter = rep(
      levels(data1$parameter),
      each = 9,
      length.out = length(unlist(Con2))
    )
  )

# Calculate effect size for ETR max

data1_ETR <- subset(data1, data1$parameter == "maxETR")

data1_ETR %>%
  group_split(group) %>%
  map( ~ .x %>%
         ConoverTest(
           data ~ grouping,
           method = "BH",
           data = .,
           out.list = T
         ))

data1_ETR %>%
  group_split(group) %>%
  map( ~ .x %>%
         cohens_d(data ~ grouping, data = .))

# Garbage collection: call after large objects have been removed

gc()

# Delete workspace, clean environment/variables and restarte R if used memory piles up

dev.off()

rm(list = ls())

.rs.restartR()
