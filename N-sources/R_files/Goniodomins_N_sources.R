##########################################
## AP1 N-growth Experiment; GDA plots and statistical analysis of A pseudogonyaulax strains L2-D2 and L4-B1
## Published in Limnology & Oceanography: "Effects of bottom-up factors on growth and toxin content of a harmful algae bloom dinoflagellate"
## All raw-data available on PANGAEA: https://doi.pangaea.de/10.1594/PANGAEA.965195 
## Questions to: kristof-moeller@outlook.de
## Kristof Möller 09/22
## Alfred-Wegener-Institute Bremerhaven
##########################################

# Install and load packages ################################
pacman::p_load(
  ggplot2,
  stats,
  rstatix,
  dplyr,
  tidyverse,
  outliers,
  NCmisc,
  patchwork,
  ggthemes,
  extrafont,
  conover.test,
  purrr,
  grid,
  flextable,
  ggpubr,
  ggprism,
  scales,
  glue,
  markdown
)

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
data[, 2:5] <-
  sapply(data[, 2:5], function(x)
    as.numeric(gsub(",", ".", x)))

data1[, 2:3] <-
  sapply(data1[, 2:3], function(x)
    as.numeric(gsub(",", ".", x)))

# Import cell sizes to normalize toxins to cell volume ################################
cell_size <- read.csv(
  "C:\\Users\\krist\\OneDrive\\Dokumente\\AWI\\Promotion\\AP1\\N-AP1a\\R_txt_files\\cell_size_L2D2_all.txt",
  sep = "",
  header = TRUE
)

# General data transformations
cell_size <- cell_size %>% t() %>% data.frame()

cell_size <-
  cell_size %>% dplyr::rename("ex" = "X1", "st" = "X2") %>%
  mutate(group = as.factor(rep(c(
    "C", "NO3", "NH4", "Urea"
  ), each = 3))) %>%
  pivot_longer(cols = c(1:2)) %>%
  filter(group != "Urea" | name != "ex") %>%
  mutate(treat2 = paste(group, name, sep = "_")) %>%
  dplyr::rename(treat = group, gp = name, data = value) %>%
  convert_as_factor(treat, treat2, gp)

cell_size$volume <- (cell_size$data / 2) ^ 3 * (4 / 3) * pi

# Combine data for analysis ################################
ex = data.frame(
  treat = rep(
    c("C_e", "C_s", "NH4_e", "NH4_s", "NO3_e", "NO3_s", "Urea_s"),
    each = 3
  ),
  GDA_C = as.numeric(data$GDA_C),
  GDA_cell = as.numeric(data$GDA_cell),
  POC = data1$C_.ng.Zelle.,
  PON = data1$N_.ng.Zelle.,
  "C/N" = data1$C_.ng.Zelle. / data1$N_.ng.Zelle.,
  volume = cell_size$volume,
  cells_mL = data$Zellen.mL,
  day = data$day
)

ex$GDA_vol <- ex$GDA_cell / ex$volume


# GDA_C ################################
# statistics
Stat_2D2 <- data
Stat_2D2$GDA_vol <- ex$GDA_vol

Stat_2D2$group <- as.factor(c(rep(c(1:7), each = 3)))
Stat_2D2$group2 <-
  c(rep(
    c("Exponential phase", "Stationary phase"),
    each = 3,
    length.out = 18
  ), rep("Stationary phase", each = 3))

# Statistics  ################################
# Stat1: exponential phase; Stat2: stationary phase
Stat1 <- subset(Stat_2D2, group2 == "Exponential phase")
Stat2 <- subset(Stat_2D2, group2 == "Stationary phase")

# Perform Kruskal Wallis Test
kruskal.test(GDA_C ~ group, data = Stat1)
kruskal.test(GDA_C ~ group, data = Stat2)
kruskal.test(GDA_C ~ group, data = Stat_2D2)

# Perform Conover Iman Post-Hoc test
Con_all <-
  conover.test(Stat_2D2$GDA_C,
               Stat_2D2$group,
               method = "BH",
               altp = T)

# extract p-values < 0.05 and order for the plot
df_p_val <-
  data.frame(
    p.adj = sprintf("%.3f", Con_all$altP.adjusted),
    group1 = str_sub(Con_all$comparisons, 0, 1),
    group2 = str_sub(Con_all$comparisons, 5, 6)
  )

df_p_val <- df_p_val %>% subset(p.adj < 0.05)

df_p_val$arrange <-
  as.numeric(df_p_val$group2) - as.numeric(df_p_val$group1)

df_p_val <- df_p_val %>% arrange(as.numeric(abs(arrange)))

# report p-values below 0.001 as < 0.001
df_p_val$p.adj3 <-
  ifelse(df_p_val$p.adj == "0.000", 0.001, df_p_val$p.adj)

# Include Cohen's d effect size
effsize <- cohens_d(Stat_2D2, GDA_C ~ group, hedges.correction = T)

effsize_sub <-
  abs(as.numeric(format(round(
    rbind(effsize$effsize[c(2, 8, 4, 10, 11)]), 1
  ) * -1, nsmall = 1)))

text1 <- data.frame(
  group = unique(Stat_2D2$group),
  'GDA_C' = aggregate(Stat_2D2, GDA_C ~ group, FUN = "max")$'GDA_C',
  group2 = c(rep(
    c("Exponential phase", "Stationary phase"),
    each = 1,
    length.out = 6
  ), "Stationary phase"),
  effsize = c(NA, NA, effsize_sub)
)

# GDA_cell ################################

# Perform Kruskal Wallis Test
res.aov1 <- kruskal_test(GDA_cell ~ group, data = Stat_2D2)

# Perform Conover Iman Post-Hoc test
Con_all <-
  conover.test(Stat_2D2$GDA_cell,
               Stat_2D2$group,
               method = "BH",
               altp = T)

# extract p-values < 0.05 and order for the plot

df_p_val2 <-
  data.frame(
    p.adj = sprintf("%.3f", Con_all$altP.adjusted),
    group1 = str_sub(Con_all$comparisons, 0, 1),
    group2 = str_sub(Con_all$comparisons, 5, 6)
  )

df_p_val2 <- df_p_val2 %>% subset(p.adj < 0.05)

df_p_val2$arrange <-
  as.numeric(df_p_val2$group2) - as.numeric(df_p_val2$group1)

df_p_val2 <- df_p_val2 %>% arrange(as.numeric(abs(arrange)))

# report p-values below 0.001 as < 0.001
df_p_val2$p.adj3 <-
  ifelse(df_p_val2$p.adj == "0.000", 0.001, df_p_val2$p.adj)

# Calculate Cohen's d effect size
effsize <-
  cohens_d(Stat_2D2, GDA_cell ~ group, hedges.correction = T)

effsize_sub <-
  abs(as.numeric(format(round(
    rbind(effsize$effsize[c(2, 8, 4, 10, 11)]), 1
  ) * -1, nsmall = 1)))

text1a <- data.frame(
  group = unique(Stat_2D2$group),
  'GDA_cell' = aggregate(Stat_2D2, GDA_cell ~ group, FUN =
                           "max")$'GDA_cell',
  group2 = c(rep(
    c("Exponential phase", "Stationary phase"),
    each = 1,
    length.out = 6
  ), "Stationary phase"),
  effsize = c(NA, NA, effsize_sub)
)

# Make GDA_cell ggplot  ################################
dodge <- position_dodge(width = 0.9)

scientific_10 <- function(x) {
  parse(text1a = gsub("e", " %*% 10^", scales::scientific_format()(x)))
}

# colorblind colors but skip black (first color) and switch blueish and greenish because NH4 / NO3 are close to each other
colors <- c(colorblind_pal()(3))[c(2:3)]

# Format the statistic with two digits after the decimal point
formatted_statistic <- sprintf("%.2f", res.aov1$statistic)

# Manually create the subtitle with markdown syntax using glue
statistic_label <-
  glue::glue(
    "<b>a) </b>Kruskal-Wallis, <i>H</i><sub>{res.aov1$df}</sub> = {formatted_statistic}, <i>p</i> = {res.aov1$p}"
  )

Stat_2D2 <- Stat_2D2 %>% mutate(treat = c(rep(c("N-deplete", "NH4", "NO3"), each = 6), rep("urea", each = 3)))

df_pval_test <- df_p_val2[c(1,3,5), ]
df_pval_test$group1 <- c("N-deplete", "NH4", "NO3")
df_pval_test$group2 <- c("N-deplete", "NH4", "NO3")

labels = c(bquote("N-deplete"),
           bquote("NH"[4] ^ + ~ ""),
           bquote("NO"[3] ^ - ~ ""),
           bquote("urea"))

P1b <- ggplot(Stat_2D2, aes(y = GDA_cell, x = treat, group = group2)) +
  # geom_point(position=position_nudge(x = -0.45), show.legend = F) +
  add_pvalue(
    df_pval_test[,],
    label = "p < {p.adj3}",
    tip.length = 0.03,
    fontface = "bold",
    y.position = 7.5,
    label.size = 2.75,
    step.increase = 0,
    bracket.shorten = -0.35,
    vjust = -0.5
  ) +
  geom_point(aes(col = group2), position = position_dodge(width = 0.7)) +
  scale_color_manual(values = colors) +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    panel.background = element_rect(fill = 'white'),
    legend.position = "none",
    plot.margin = unit(c(0, 0, -0.5, 0), 'cm'),
    plot.subtitle = element_markdown(),
    strip.text = element_text(size = 14)
  ) +
  # geom_text(data = text1a, aes(x = group, label = effsize, y = GDA_cell), nudge_y = 0.4, show.legend = F)+
  scale_x_discrete(labels = rep(labels, each = 1))  +
  scale_y_continuous(expression("Toxin content (pg GDA cell " ^ -1 * ")")) +
  xlab("") +
  coord_cartesian(ylim = c(0, 10)) +
  labs(
    subtitle = statistic_label)

# Toxin normalization comparison ################################

# colorblind colors but skip black (first color) and switch blueish and greenish because NH4 / NO3 are close to each other
colors <- c(colorblind_pal()(3))[c(2, 3)]

Pex1 <- ggplot(Stat_2D2, aes(y = GDA_cell, x = treat, group = group2)) +
  geom_point(aes(col = group2), position = position_dodge(width = 0.7)) +
  theme_classic() +
  scale_color_manual(values = colors) +
  theme(
    legend.title = element_blank(),
    panel.background = element_rect(fill = 'white'),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, vjust = 0.5),
    plot.margin = unit(c(0, 0, 0, 0), 'cm'),
    plot.subtitle = element_markdown()
  ) +
  scale_x_discrete(labels = rep(labels, each = 1))  +
  scale_y_continuous(expression("Toxin content (pg GDA cell " ^ -1 * ")")) +
  xlab("") +
  labs(subtitle = "<b>a)</b> per cell:")

Pex2 <- ggplot(Stat_2D2, aes(y = GDA_C, x = treat, group = group2)) +
  geom_point(aes(col = group2), position = position_dodge(width = 0.7)) +
  theme_classic() +
  scale_color_manual(values = colors) +
  theme(
    legend.title = element_blank(),
    panel.background = element_rect(fill = 'white'),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, vjust = 0.5),
    plot.margin = unit(c(0, 0, 0, 0), 'cm'),
    plot.subtitle = element_markdown()
  ) +
  scale_x_discrete(labels = rep(labels, each = 1))  +
  scale_y_continuous(expression("Toxin content (GDA C" ^ -1 * " mol:mol)"),
                     label = scientific_format(digits = 2)) +
  xlab("") +
  labs(subtitle = "<b>a)</b> per carbon:")

Pex3 <- ggplot(Stat_2D2, aes(y = GDA_vol, x = treat, group = group2)) +
  geom_point(aes(col = group2), position = position_dodge(width = 0.7)) +
  theme_classic() +
  scale_color_manual(values = colors) +
  theme(
    legend.title = element_blank(),
    panel.background = element_rect(fill = 'white'),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, vjust = 0.5),
    plot.margin = unit(c(0, 0, 0, 0), 'cm'),
    plot.subtitle = element_markdown()
  ) +
  scale_x_discrete(labels = rep(labels, each = 1))  +
  scale_y_continuous(expression(paste("Toxin content (pg GDA ", mu, "m" ^ -3 *
                                        ")")), label = scientific_format(digits = 2)) +
  xlab("") +
  labs(subtitle = "<b>a)</b> per volume:")

all_plots_2D2 <-
  Pex1 + Pex3 + Pex2 + plot_layout(guides = "collect", ncol = 3) &
  theme(legend.position = "none")
# ggsave("Norm_comparison.png", all_plots_2D2 ,path = "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/AP1/N-AP1a/figures", dpi = 300, width = 30, height = 15, units = "cm")

# add A. pseudogonyaulax strain L4-B1 ################################

# Import data and general data transformations ################################
# 1-9: exponential harvest without Urea;
# 10-21: stationary harvest

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
data[, 2:4] <-
  sapply(data[, 2:4], function(x)
    as.numeric(gsub(",", ".", x)))

data1[, 2:3] <-
  sapply(data1[, 2:3], function(x)
    as.numeric(gsub(",", ".", x)))

# Import cell sizes to normalize toxins to cell volume ################################
cell_size = read.csv(
  "C:\\Users\\krist\\OneDrive\\Dokumente\\AWI\\Promotion\\AP1\\N-AP1a\\R_txt_files\\cell_size_L4B1_all.txt",
  sep = "",
  header = T
)

# General data transformations
cell_size <- cell_size %>% t() %>% data.frame()

cell_size <-
  cell_size %>% dplyr::rename("ex" = "X1", "st" = "X2") %>%
  mutate(group = as.factor(rep(c(
    "C", "NO3", "NH4", "Urea"
  ), each = 3))) %>%
  pivot_longer(cols = c(1:2)) %>%
  filter(group != "Urea" | name != "ex") %>%
  mutate(treat2 = paste(group, name, sep = "_")) %>%
  dplyr::rename(treat = group, gp = name, data = value) %>%
  convert_as_factor(treat, treat2, gp)

# cell_size <- aggregate(data=cell_size, data ~ treat2, FUN="mean")
cell_size$volume <- (cell_size$data / 2) ^ 3 * (4 / 3) * pi

# Combine data for analysis ################################
ex = data.frame(
  treat = rep(
    c("C_e", "C_s", "NH4_e", "NH4_s", "NO3_e", "NO3_s", "Urea_s"),
    each = 3
  ),
  GDA_C = as.numeric(data$GDA_C),
  GDA_cell = as.numeric(data$GDA_cell),
  POC = data1$C_.ng.Zelle.,
  PON = data1$N_.ng.Zelle.,
  "C/N" = data1$C_.ng.Zelle. / data1$N_.ng.Zelle.,
  volume = cell_size$volume
)

ex$GDA_vol <- ex$GDA_cell / ex$volume

# GDA_C ################################
#statistics
Stat_4B1 <- data
Stat_4B1$GDA_vol <- ex$GDA_vol

Stat_4B1$group <- as.factor(c(rep(c(1:7), each = 3)))
Stat_4B1$group2 <-
  c(rep(
    c("Exponential phase", "Stationary phase"),
    each = 3,
    length.out = 18
  ), rep("Stationary phase", each = 3))

## Stat1: exponential phase; Stat2: stationary phase

Stat1 <- subset(Stat_4B1, group2 == "Exponential phase")
Stat2 <- subset(Stat_4B1, group2 == "Stationary phase")

# Perform Kruskal Wallis Test
res.aov2 <- kruskal_test(GDA_C ~ group, data = Stat_4B1)

# Perform Conover Iman Post-Hoc test
Con_all <-
  conover.test(Stat_4B1$GDA_C,
               Stat_4B1$group,
               method = "BH",
               altp = T)

# extract p-values < 0.05 and order for the plot
df_p_val3 <-
  data.frame(
    p.adj = sprintf("%.3f", Con_all$altP.adjusted),
    group1 = str_sub(Con_all$comparisons, 0, 1),
    group2 = str_sub(Con_all$comparisons, 5, 6)
  )
df_p_val3 <- df_p_val3 %>% subset(p.adj < 0.05)
df_p_val3$arrange <-
  as.numeric(df_p_val3$group2) - as.numeric(df_p_val3$group1)
df_p_val3 <- df_p_val3 %>% arrange(as.numeric(abs(arrange)))

# report p-values below 0.001 as < 0.001
df_p_val3$p.adj3 <-
  ifelse(df_p_val3$p.adj == "0.000", 0.001, df_p_val3$p.adj)

# Include Cohen's d effect size
effsize <- cohens_d(Stat_4B1, GDA_C ~ group, hedges.correction = T)
effsize_sub <-
  abs(as.numeric(format(round(
    rbind(effsize$effsize[c(2, 8, 4, 10, 11)]), 1
  ) * -1, nsmall = 1)))

text2 <- data.frame(
  group = unique(Stat_4B1$group),
  'GDA_C' = aggregate(Stat_4B1, GDA_C ~ group, FUN = "max")$'GDA_C',
  group2 = c(rep(
    c("Exponential phase", "Stationary phase"),
    each = 1,
    length.out = 6
  ), "Stationary phase"),
  effsize = c(NA, NA, effsize_sub)
)

# GDA_cell ################################

# Perform Kruskal Wallis Test
res.aov3 <- kruskal_test(GDA_cell ~ group, data = Stat_4B1)

# Conover-Iman PostHoc test following a significant kruskal test
Con_all <-
  conover.test(
    Stat_4B1$GDA_cell,
    Stat_4B1$group,
    method = "BH",
    altp = T,
    kw = T,
    alpha = 0.1
  )

# extract p-values < 0.05 and order for the plot
df_p_val4 <-
  data.frame(
    p.adj = sprintf("%.3f", Con_all$altP.adjusted),
    group1 = str_sub(Con_all$comparisons, 0, 1),
    group2 = str_sub(Con_all$comparisons, 5, 6)
  )
df_p_val4 <- df_p_val4 %>% subset(p.adj < 0.1)
df_p_val4$arrange <-
  as.numeric(df_p_val4$group2) - as.numeric(df_p_val4$group1)
df_p_val4 <- df_p_val4 %>% arrange(as.numeric(abs(arrange)))

# Calculate Cohen's d effect size
effsize <-
  cohens_d(Stat_4B1, GDA_cell ~ group, hedges.correction = T)
effsize_sub <-
  abs(as.numeric(format(round(
    rbind(effsize$effsize[c(2, 8, 4, 10, 11)]), 1
  ) * -1, nsmall = 1)))

text2a <- data.frame(
  group = unique(Stat_4B1$group),
  GDA_cell = aggregate(Stat_4B1, GDA_cell ~ group, FUN = "max")$GDA_cell,
  group2 = c(rep(
    c("Exponential phase", "Stationary phase"),
    each = 1,
    length.out = 6
  ), "Stationary phase"),
  effsize = c(NA, NA, effsize_sub)
)

# Make GDA_cell ggplot  ################################
dodge <- position_dodge(width = 0.9)

scientific_10 <- function(x) {
  parse(text1a = gsub("e", " %*% 10^", scales::scientific_format()(x)))
}

# Format the statistic with two digits after the decimal point
formatted_statistic <- sprintf("%.2f", res.aov3$statistic)

# Manually create the subtitle with markdown syntax using glue
statistic_label <-
  glue::glue(
    "<b>b) </b>Kruskal-Wallis, <i>H</i><sub>{res.aov3$df}</sub> = {formatted_statistic}, <i>p</i> = {res.aov3$p}"
  )

Stat_4B1 <- Stat_4B1 %>% mutate(treat = c(rep(c("N-deplete", "NH4", "NO3"), each = 6), rep("urea", each = 3)))

p_val_plot <- df_p_val4[1, ]
p_val_plot <- p_val_plot %>% mutate(group1 = "NH4", group2 = "NH4")

P2a <- ggplot(Stat_4B1, aes(y = GDA_cell, x = treat, group = group2)) +
  # geom_point(position = position_nudge(x = -0.45), show.legend = F) +
  add_pvalue(
    p_val_plot,
    label = "p = {p.adj}",
    tip.length = 0.03,
    fontface = "bold",
    y.position = 8.5,
    label.size = 2.75,
    step.increase = 0,
    bracket.shorten = -0.35,
    vjust = -0.5
  ) +
  geom_point(aes(col = group2), position = position_dodge(width = 0.7)) +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    panel.background = element_rect(fill = 'white'),
    legend.position = "none",
    plot.margin = unit(c(-0.5, 0, 0, 0), 'cm'),
    plot.subtitle = element_markdown(),
    legend.box.margin = margin(-10, -10, -10, -10),
    strip.text = element_text(size = 14)) +
  scale_color_manual(values = colors) +
  # geom_text(data = text2a, aes(x = group, label = effsize, y = GDA_cell), nudge_y = 0.4, show.legend = F) +
  scale_x_discrete(labels = rep(labels, each = 1))  +
  scale_y_continuous(expression("Toxin content (pg GDA cell " ^ -1 * ")")) +
  xlab("N-source") +
  coord_cartesian(ylim = c(0, 10)) +
  labs(subtitle = statistic_label)


# combine toxin plots of both strains in one plot  ################################
all_plots <-
  P1b + P2a + plot_layout(guides = "collect", ncol = 1, axis_titles = "collect_y") &
  theme(legend.position = "bottom",
        legend.box.margin = margin(-10, -10, -10, -10))

ggsave(
  "Toxins_all.png",
  all_plots ,
  path = "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/AP1/N-AP1a/figures",
  dpi = 300,
  width = 3.5,
  height = 4.9,
  units = "in"
)

# Toxin normalization comparison ################################
# colorblind colors but skip black (first color) and switch blueish and greenish because NH4 / NO3 are close to each other
colors <- c(colorblind_pal()(3))[c(2, 3)]

labels = c(bquote("N-deplete"),
           bquote("NH"[4] ^ + ~ ""),
           bquote("NO"[3] ^ - ~ ""),
           bquote("urea"))

Stat_4B1 <- Stat_4B1 %>% mutate(treat = c(rep(c("N-deplete", "NH4", "NO3"), each = 6), rep("urea", each = 3)))

Pex1_4B1 <- ggplot(Stat_4B1, aes(y = GDA_cell, x = treat, group = group2)) +
  geom_point(aes(col = group2), position = position_dodge(width = 0.7)) +
  theme_classic() +
  scale_color_manual(values = colors) +
  theme(
    legend.title = element_blank(),
    panel.background = element_rect(fill = 'white'),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, vjust = 0.5),
    plot.margin = unit(c(0, 0, 0, 0), 'cm'),
    plot.subtitle = element_markdown()
  ) +
  scale_x_discrete(labels = rep(labels, each = 1))  +
  scale_y_continuous(expression("Toxin content (pg GDA cell " ^ -1 * ")")) +
  xlab("") +
  labs(subtitle = "<b>b)</b> per cell:")

Pex2_4B1 <- ggplot(Stat_4B1, aes(y = GDA_C, x = treat, group = group2)) +
  geom_point(aes(col = group2), position = position_dodge(width = 0.7)) +
  theme_classic() +
  scale_color_manual(values = colors) +
  theme(
    legend.title = element_blank(),
    panel.background = element_rect(fill = 'white'),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, vjust = 0.5),
    plot.margin = unit(c(0, 0, -0, 0), 'cm'),
    plot.subtitle = element_markdown()
  ) +
  scale_x_discrete(labels = rep(labels, each = 1))  +
  scale_y_continuous(expression("Toxin content (GDA C" ^ -1 * " mol:mol)"),
                     label = scientific_format(digits = 2)) +
  xlab("") +
  labs(subtitle = "<b>b)</b> per carbon:")

Pex3_4B1 <- ggplot(Stat_4B1, aes(y = GDA_vol, x = treat, group = group2)) +
  geom_point(aes(col = group2), position = position_dodge(width = 0.7)) +
  theme_classic() +
  scale_color_manual(values = colors) +
  theme(
    legend.title = element_blank(),
    panel.background = element_rect(fill = 'white'),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, vjust = 0.5),
    plot.margin = unit(c(0, 0, 0, 0), 'cm'),
    plot.subtitle = element_markdown()
  ) +
  scale_x_discrete(labels = rep(labels, each = 1))  +
  scale_y_continuous(expression(paste("Toxin content (pg GDA ", mu, "m" ^ -3 *
                                        ")")), label = scientific_format(digits = 2)) +
  xlab("") +
  labs(subtitle = "<b>b)</b> per volume:")

# combine toxin normalization comparisons of both strains in one plot  ################################
# Combine Pex1 with Pex1_4B1
plot1_combined <- Pex1 + Pex1_4B1 + plot_layout(guides = "collect", ncol = 1, axis_titles = "collect_y") &
  theme(legend.position = "none",
        legend.box.margin = margin(-10, -10, -10, -10))

# Combine Pex3 with Pex3_4B1
plot2_combined <- Pex3 + Pex3_4B1 + plot_layout(guides = "collect", ncol = 1, axis_titles = "collect_y")  &
  theme(legend.position = "bottom",
        legend.box.margin = margin(-10, -10, -10, -10))

# Combine Pex2 with Pex2_4B1
plot3_combined <- Pex2 + Pex2_4B1 + plot_layout(guides = "collect", ncol = 1, axis_titles = "collect_y")  &
  theme(legend.position = "none",
        legend.box.margin = margin(-10, -10, -10, -10))

# Arrange all combined plots in a grid layout
all_plots_comparison <- wrap_plots(plot1_combined, plot2_combined, plot3_combined, ncol = 3)

ggsave(
  "all_plots_comparison.png",
  all_plots_comparison ,
  path = "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/AP1/N-AP1a/figures",
  dpi = 300,
  width = 30,
  height = 15,
  units = "cm"
)

# Export Mean toxin quotas and standard deviations for table in manuscript  ################################
agg <-
  format(aggregate(Stat_2D2, GDA_C ~ group, FUN = "mean"),
         digits = 3,
         format = "e")
agg$sd <-
  format(
    aggregate(Stat_2D2, GDA_C ~ group, FUN = "sd")$'GDA_C',
    digits = 3,
    format = "e"
  )

agg$GDA_cell <-
  format(aggregate(Stat_2D2, GDA_cell ~ group, FUN = "mean")$'GDA_cell',
         digits = 2)
agg$GDA_cell_sd <-
  format(aggregate(Stat_2D2, GDA_cell ~ group, FUN = "sd")$'GDA_cell',
         digits = 2)

agg$GDA_vol <-
  format(
    aggregate(Stat_2D2, GDA_vol ~ group, FUN = "mean")$'GDA_vol',
    digits = 3,
    format = "e"
  )
agg$GDA_vol_sd <-
  format(
    aggregate(Stat_2D2, GDA_vol ~ group, FUN = "sd")$'GDA_vol',
    digits = 3,
    format = "e"
  )

L2D2_means <-
  data.frame(GDA_C = paste(agg$GDA_C, agg$sd, sep = "+/-"))

agg2 <-
  format(aggregate(Stat_4B1, GDA_C ~ group, FUN = "mean"),
         digits = 3,
         format = "e")
agg2$sd <-
  format(
    aggregate(Stat_4B1, GDA_C ~ group, FUN = "sd")$'GDA_C',
    digits = 3,
    format = "e"
  )

agg2$GDA_cell <-
  format(aggregate(Stat_4B1, GDA_cell ~ group, FUN = "mean")$'GDA_cell',
         digits = 2)
agg2$GDA_cell_sd <-
  format(aggregate(Stat_4B1, GDA_cell ~ group, FUN = "sd")$'GDA_cell',
         digits = 2)

agg2$GDA_vol <-
  format(
    aggregate(Stat_4B1, GDA_vol ~ group, FUN = "mean")$'GDA_vol',
    digits = 3,
    format = "e"
  )
agg2$GDA_vol_sd <-
  format(
    aggregate(Stat_4B1, GDA_vol ~ group, FUN = "sd")$'GDA_vol',
    digits = 3,
    format = "e"
  )

agg_total <- rbind(agg, agg2)

agg_total$strain <-
  rep(c("L2-D2", "L4-B1"), each = length(agg_total$'GDA_C') / 2)

agg_total <- agg_total %>% arrange(GDA_vol)

agg_total %>% arrange(GDA_C)

agg_total <- flextable(agg_total) %>% autofit %>%
  save_as_docx(path = "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/AP1/N-AP1a/GDA_C.docx")

# Garbage collection: call after large objects have been removed  ################################

gc()

# Delete workspace, clean environment/variables and restarte R if used memory piles up

dev.off()

rm(list = ls())

.rs.restartR()
