##########################################
## AP1 N-growth Experiment; POC/PON quota of both A. pseudogonyaulax strains L2-D2 and L4-B1
## Published in Limnology & Oceanography: "Effects of bottom-up factors on growth and toxin content of a harmful algae bloom dinoflagellate"
## All raw-data available on PANGAEA: https://doi.pangaea.de/10.1594/PANGAEA.965195 
## Questions to: kristof-moeller@outlook.de
## Kristof Möller 09/22
## Alfred-Wegener-Institute Bremerhaven
##########################################

# Install and load packages ################################

# Installs pacman ("package manager") if needed
if (!require("pacman"))
  install.packages("pacman")

# Use pacman to load add-on packages as desired
pacman::p_load(pacman,
               ggplot2,
               scales,
               stats,
               rstatix,
               dplyr,
               extrafont,
               conover.test,
               flextable)

# Load Windows Fonts and define Times font
loadfonts(device = "win")
windowsFonts(Times = windowsFont("Times"))

# List used packages in the file
# used_packages <- list.functions.in.file("C:\\Users\\krist\\OneDrive\\Dokumente\\AWI\\Promotion\\AP1\\N-AP1a\\R_files\\Toxin_all.R")

# Import data and general data transformations ################################
# 1-9: exponential harvest without Urea;
# 10-21: stationary harvest

# Toxin abundance per Cell 
data <-
  read.csv(
    "C:\\Users\\krist\\OneDrive\\Dokumente\\AWI\\Promotion\\AP1\\N-AP1a\\R_txt_files\\GDA_quota_L2D2.txt",
    sep = "",
    header = T
  )

# POC/PON per Cell
data1 <-
  read.csv(
    "C:\\Users\\krist\\OneDrive\\Dokumente\\AWI\\Promotion\\AP1\\N-AP1a\\R_txt_files\\POC_PON_L2D2.txt",
    sep = "",
    header = TRUE
  )

# Clean data by replacing commas with periods and converting to numeric
data1[, 2:3] <-
  sapply(data1[, 2:3], function(x)
    as.numeric(gsub(",", ".", x)))

data[, 2:3] <-
  sapply(data[, 2:3], function(x)
    as.numeric(gsub(",", ".", x)))

# Create data frame for analysis
ex = data.frame(
  treat = rep(
    c("C_e", "C_s", "NH4_e", "NH4_s", "NO3_e", "NO3_s", "Urea_s"),
    each = 3
  ),
  GDA = data$GDA_cell,
  POC = data1$C_.ng.Zelle. * 10 ^ -9 / 12.001, # Convert mass per cell to molecular mass per cell 
  PON = data1$N_.ng.Zelle. * 10 ^ -9 / 14.0067, # Convert mass per cell to molecular mass per cell 
  "C/N" = (data1$C_.ng.Zelle. * 12.011) / (data1$N_.ng.Zelle. * 14.0067)
)

# POC / PON ################################
# statistics
# Prepare data
Stat = ex

# Color vector for exponential / stationary phase:
Stat$col <-
  c(rep(
    c("exponential", "stationary"),
    each = 3,
    length.out = 15
  ), rep(c("stationary"), length.out = 6))

# Perform Kruskal Wallis Test
kruskal.test(Stat$C.N ~ Stat$treat)

# Perform Conover Iman Post-Hoc test
Con <- conover.test(Stat$C.N, Stat$treat, method = "BH", altp = T)

# Calculate Cohen's d effect size
effsize <- cohens_d(Stat, C.N ~ treat, hedges.correction = T)
effsize_sub <-
  abs(as.numeric(format(round(
    rbind(effsize$effsize[c(2, 4, 8, 10, 11)]), 1
  ) * -1, nsmall = 1)))


# Make ggplot
# Transformation function for x-axis labels to have 1 digit ################################
scaleFUN <- function(x)
  sprintf("%.1f", x)

P1 <- ggplot(Stat, aes(x = treat, y = C.N, fill = col)) +
  geom_boxplot(show.legend = F) +
  geom_point(show.legend = F, position = position_nudge(x = -0.45)) +
  ggtitle("A)") +
  theme_classic() +
  scale_fill_grey(start = 0.4, end = 0.8) +
  theme(
    legend.title = element_blank(),
    panel.background = element_rect(fill = 'white'),
    legend.position = "bottom"
  ) +
  xlab(expression(paste(""))) +
  scale_x_discrete(labels = c("C", "C", "NH4", "NH4", "NO3", "NO3", "Urea")) +
  scale_y_continuous(expression("C:N ratio"), labels = scaleFUN) +
  xlab("")

## Export Mean POC/PON and standard deviation for table in manuscript ################################

agg <- aggregate(Stat, C.N ~ treat, FUN = "mean")
agg$sd <- aggregate(Stat, C.N ~ treat, FUN = "SD")$'C.N'

agg$POC <-
  aggregate(Stat, C_.ng.Zelle. ~ treat, FUN = "mean")$'C_.ng.Zelle.'
agg$POC_sd <-
  aggregate(Stat, C_.ng.Zelle. ~ treat, FUN = "SD")$'C_.ng.Zelle.'

agg$PON <-
  aggregate(Stat, N_.ng.Zelle. ~ treat, FUN = "mean")$'N_.ng.Zelle.'
agg$PON_sd <-
  aggregate(Stat, N_.ng.Zelle. ~ treat, FUN = "SD")$'N_.ng.Zelle.'

agg[, c(2:7)] <- format(round(agg[, c(2:7)], 2), nsmall = 2)
agg$data_CN <- paste(agg$C.N, agg$sd, sep = " +/- ")
agg$data_POC <- paste(agg$POC, agg$POC_sd, sep = " +/- ")
agg$data_PON <- paste(agg$PON, agg$PON_sd, sep = " +/- ")

#
# ## Mean and sd for table (GDA:C):
#
# agg <- aggregate(Stat, C.N ~ treat, FUN= "mean")
# agg$sd <- aggregate(Stat, C.N ~ treat, FUN= "sd")$'C.N'
#
# agg[,c(2:3)] <- format(round(agg[,c(2:3)], 2), nsmall=2)
# agg$data <- paste(agg$C.N, agg$sd, sep = " +/- ")

# add A. pseudogonyaulax strain L4-B1 as well  ################################
# Import data and general data transformations ################################

# Toxin abundance per Cell 
data <-
  read.csv(
    "C:\\Users\\krist\\OneDrive\\Dokumente\\AWI\\Promotion\\AP1\\N-AP1a\\R_txt_files\\GDA_quota_L4B1.txt",
    sep = "",
    header = T
  )

# POC/PON per Cell
data1 <-
  read.csv(
    "C:\\Users\\krist\\OneDrive\\Dokumente\\AWI\\Promotion\\AP1\\N-AP1a\\R_txt_files\\POC_PON_L4B1.txt",
    sep = "",
    header = TRUE
  )

# Clean data by replacing commas with periods and converting to numeric
data1[, 2:3] <-
  sapply(data1[, 2:3], function(x)
    as.numeric(gsub(",", ".", x)))

data[, 2:3] <-
  sapply(data[, 2:3], function(x)
    as.numeric(gsub(",", ".", x)))

# Create data frame for analysis
ex2 = data.frame(
  treat = rep(
    c("C_e", "C_s", "NH4_e", "NH4_s", "NO3_e", "NO3_s", "Urea_s"),
    each = 3
  ),
  GDA = data$GDA_cell,
  POC = data1$C_.ng.Zelle. * 10 ^ -9 / 12.001, # Convert mass per cell to molecular mass per cell 
  PON = data1$N_.ng.Zelle. * 10 ^ -9 / 14.0067, # Convert mass per cell to molecular mass per cell 
  "C/N" = (data1$C_.ng.Zelle. * 12.011) / (data1$N_.ng.Zelle. * 14.0067)
)

# POC / PON ################################
#statistics
# Prepare data
Stat2 = ex

#Color vector for exponential / Stat2ionary phase:
Stat2$col <-
  c(rep(
    c("exponential", "stationary"),
    each = 3,
    length.out = 15
  ), rep(c("stationary"), length.out = 6))

# Perform Kruskal Wallis Test
kruskal.test(Stat2$C.N ~ Stat2$treat)

# Perform Conover Iman Post-Hoc test
Con <-
  conover.test(Stat2$C.N, Stat2$treat, method = "BH", altp = T)

# Calculate Cohen's d effect size
effsize <- cohens_d(Stat2, C.N ~ treat, hedges.correction = T)
effsize_sub <-
  abs(as.numeric(format(round(
    rbind(effsize$effsize[c(2, 4, 8, 10, 11)]), 1
  ) * -1, nsmall = 1)))

# Make ggplot ################################

## Transformation function for x-axis labels to have 1 digit
scaleFUN <- function(x)
  sprintf("%.1f", x)

P2 <- ggplot(Stat2, aes(x = group2, y = C.N, fill = col)) +
  geom_boxplot(show.legend = T) +
  geom_point(show.legend = F, position = position_nudge(x = -0.45)) +
  ggtitle("B)") +
  theme_classic() +
  scale_fill_grey(start = 0.4, end = 0.8) +
  theme(
    legend.title = element_blank(),
    panel.background = element_rect(fill = 'white'),
    legend.position = "bottom"
  ) +
  xlab(expression(paste(""))) +
  scale_x_discrete(labels = c("C", "C", "NH4", "NH4", "NO3", "NO3", "Urea")) +
  scale_y_continuous(expression("C:N ratio"), labels = scaleFUN) +
  xlab("")

## Export Mean POC/PON and standard deviation for table in manuscript  ################################

agg2 <- aggregate(Stat2, C.N ~ treat, FUN = "mean")
agg2$sd <- aggregate(Stat2, C.N ~ treat, FUN = "SD")$'C.N'

agg2$POC <-
  aggregate(Stat2, C_.ng.Zelle. ~ treat, FUN = "mean")$'C_.ng.Zelle.'
agg2$POC_sd <-
  aggregate(Stat2, C_.ng.Zelle. ~ treat, FUN = "SD")$'C_.ng.Zelle.'

agg2$PON <-
  aggregate(Stat2, N_.ng.Zelle. ~ treat, FUN = "mean")$'N_.ng.Zelle.'
agg2$PON_sd <-
  aggregate(Stat2, N_.ng.Zelle. ~ treat, FUN = "SD")$'N_.ng.Zelle.'

agg2[, c(2:7)] <- format(round(agg2[, c(2:7)], 2), nsmall = 2)
agg2$data_CN <- paste(agg2$C.N, agg2$sd, sep = " +/- ")
agg2$data_POC <- paste(agg2$POC, agg2$POC_sd, sep = " +/- ")
agg2$data_PON <- paste(agg2$PON, agg2$PON_sd, sep = " +/- ")

#
# ## Combine both agg and agg2 in one data.frame and export as table
#
# agg_total <- rbind(agg, agg2)
# agg_total$strain <- rep(c("L2-D2", "L4-B1"), each=length(agg_total$'C.N')/2)
#
# agg_total <- flextable(agg_total) %>% autofit %>%
# save_as_docx(path="C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/AP1/N-AP1a/POC_PON.docx")

# Between-strains statistics ################################
# Growth rates
ex <- ex %>% mutate(strain = "L2-D2")
ex2 <- ex2 %>% mutate(strain = "L4-B1")

ex_all <- rbind(ex, ex2)

POC_all <- ex_all %>%
  group_by(treat) %>%
  nest() %>%
  mutate(kruskal_results = map(data, ~kruskal_test(POC ~ strain, data = .x))) %>%
  unnest(c(kruskal_results)) 

PON_all <- ex_all %>%
  group_by(treat) %>%
  nest() %>%
  mutate(kruskal_results = map(data, ~kruskal_test(PON ~ strain, data = .x))) %>%
  unnest(c(kruskal_results)) 

CN_all <- ex_all %>%
  group_by(treat) %>%
  nest() %>%
  mutate(kruskal_results = map(data, ~kruskal_test(C.N ~ strain, data = .x))) %>%
  unnest(c(kruskal_results)) 

# combine cell count plots of both strains in one plot  ################################
g <- grid.arrange(P1, P2, ncol = 1)

# save the combined plot as a png
ggsave(
  "POC_PON_all.png",
  g,
  path = "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/AP1/N-AP1a/figures",
  dpi = 300,
  width = 15,
  height = 20,
  units = "cm"
)

# Garbage collection: call after large objects have been removed  ################################

gc()

# Delete workspace, clean environment/variables and restarte R if used memory piles up

dev.off()

rm(list = ls())

.rs.restartR()
