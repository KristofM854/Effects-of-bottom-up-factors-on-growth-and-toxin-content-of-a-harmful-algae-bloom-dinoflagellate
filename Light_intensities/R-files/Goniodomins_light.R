##########################################
## AP1 N-light Experiment; GDA plots and statistical analysis of A pseudogonyaulax strains L2-D2, L4-B1 and L4-B9
## Published in Limnology & Oceanography: "Effects of bottom-up factors on growth and toxin content of a harmful algae bloom dinoflagellate"
## All raw-data available on PANGAEA: https://doi.pangaea.de/10.1594/PANGAEA.965195 
## Questions to: kristof-moeller@outlook.de
## Kristof Möller 09/22
## Alfred-Wegener-Institute Bremerhaven
##########################################

# Install and load packages ################################
pacman::p_load(
  ggplot2,
  scales,
  emmeans,
  ggprism,
  stats,
  rstatix,
  dplyr,
  tidyverse,
  ggpattern,
  ggtext,
  gridExtra,
  ggpubr,
  conover.test,
  grid,
  flextable,
  patchwork,
  ggthemes,
  PMCMRplus,
  extrafont
)

# Load Windows Fonts and define Times font
loadfonts(device = "win")
windowsFonts(Times = windowsFont("Times"))

# Import data and general data transformations ################################

# cell size data of all strains
cell_size <-
  read.csv(
    "C:\\Users\\krist\\OneDrive\\Dokumente\\AWI\\Promotion\\AP1\\Light-AP1b\\raw\\cell_sizes.txt",
    sep = "",
    header = T
  )

# Clean data by replacing commas with periods and converting to numeric
cell_size[] <-
  lapply(cell_size, function(x)
    as.numeric(gsub(",", ".", x)))

# read in POC/PON data of all strains for toxin normalization per carbon 
CN <-
  read.csv(
    "C:\\Users\\krist\\OneDrive\\Dokumente\\AWI\\Promotion\\AP1\\Light-AP1b\\raw\\POC_PON_AP1_light.txt",
    sep = "",
    header = T
  )

# introduce replicate grouping factor (all triplicates)
CN$replicate <- rep(c(1:3))

# import data A. pseudogonyaulax strain L2-D2 #####
data <-
  read.csv(
    "C:\\Users\\krist\\OneDrive\\Dokumente\\AWI\\Promotion\\AP1\\Light-AP1b\\raw\\data_L2D2.txt",
    sep = "",
    header = FALSE,
    skip = 14,   # Skip the first 13 rows
    nrows = 9    # Read 9 rows
  )

# put first column as row names and transpose the data frame
data <- data %>% remove_rownames %>% column_to_rownames(var = "V1")
data <- as.data.frame(t(data))

#replace all , with .
data[] <- lapply(data, function(x)
  as.numeric(gsub(",", ".", x)))

# reshape dataframe 
data1 <-
  data.frame(
    data1 = unlist(data),
    treat = rep(
      c("20 Einstein", "100 Einstein", "200 Einstein"), # three light intensity treatments
      each = 15,
      length.out = 45
    ),
    group = c(
      rep(c(1:5), length.out = 15), # up to 5 toxin sampling days
      rep(c(6:10), length.out = 15),
      rep(c(11:15), length.out = 15)
    )
  )

# add cell sizes column, ex = Exponential and st = Stationary phase 
data1$size <- c(
  rep(cell_size[1:3, 2], each = 5, length.out = 15),
  rep(cell_size[4, 2], length.out = 3),
  rep(cell_size[4, 5], length.out = 2),
  rep(cell_size[5, 2], length.out = 3),
  rep(cell_size[5, 5], length.out = 2),
  rep(cell_size[6, 2], length.out = 3),
  rep(cell_size[6, 5], length.out = 2),
  rep(cell_size[7, 2], length.out = 3),
  rep(cell_size[7, 5], length.out = 2),
  rep(cell_size[8, 2], length.out = 3),
  rep(cell_size[8, 5], length.out = 2),
  rep(cell_size[9, 2], length.out = 3),
  rep(cell_size[9, 5], length.out = 2)
)

# add toxin per volume
data1$tox_vol <- data1$data1/((4/3)*pi*(data1$size)^3)

# Statistical analysis

# Subset of first Exponential and Stationary phase sampling point for manuscript
Stat <- data1 %>% filter(group == 1 | group == 6 | group == 9 | group == 11 | group == 14)

# Perform Kruskal-Wallis Test
res.aov1 <- kruskal_test(data1 ~ group, data = Stat)

# reorder factor levels
Stat <- Stat[order(Stat$group), ]
Stat$group <- factor(Stat$group)

# Perform Conover-Iman posthoc test
Con_a <- conover.test(Stat$data1, Stat$group, method = "BH", altp = T)

# Calculate Cohen's effect size
effsize <- cohens_d(Stat, data1 ~ group, hedges.correction = T)
effsize <-
  format(round(effsize[1:4, 4], 1) * -1, nsmall = 1, digits = 2)
effsize <- as.vector(abs(as.numeric(effsize$effsize)))

# introduce factor for Exponential / Stationary phase:
Stat$col <- as.factor(c(rep("Exponential phase", each = 6), rep(
  c("Stationary phase", "Exponential phase"),
  each = 3,
  length.out = 9
)))

text <- data.frame(
  group = unique(Stat$group),
  effsize = c(NA, effsize),
  data = aggregate(Stat, data1 ~ group, FUN = "max")$data,
  col = c(rep("Exponential phase", each = 2), rep(
    c("Stationary phase", "Exponential phase"),
    each = 1,
    length.out = 3
  )),
  GDA_sd = aggregate(Stat, data1 ~ group, FUN = "sd")$data
)

# extract p-values < 0.05 and order for the plot
df_p_val <-
  data.frame(
    p.adj = sprintf("%.3f", Con_a$altP.adjusted),
    group1 = sub("^([0-9]+) .*", "\\1", Con_a$comparisons),
    group2 = sub(".* - ([0-9]+)", "\\1", Con_a$comparisons)
  )

df_p_val <- df_p_val %>% subset(p.adj < 0.05)

df_p_val[] <-
  sapply(df_p_val, function(x)
    gsub(pattern = " ", replacement = "", x))

df_p_val$arrange <-
  as.numeric(df_p_val$group2) - as.numeric(df_p_val$group1)

df_p_val <- df_p_val %>% arrange(as.numeric(abs(arrange)))

Stat$group <- as.numeric(as.character(Stat$group))

# Prepare ggplot parameters
dodge <- position_dodge(width = 0.9)

# colorblind colors but skip black (first color) and switch blueish and greenish because NH4 / NO3 are close to each other
colors <- c(colorblind_pal()(3))[c(2, 3)]

# Format the statistic with two digits after the decimal point
formatted_statistic <- sprintf("%.2f", res.aov1$statistic)

# Manually create the subtitle with markdown syntax using glue
statistic_label <-
  glue::glue(
    "<b>a) </b>Kruskal-Wallis, <i>H</i><sub>{res.aov1$df}</sub> = {formatted_statistic}, <i>p</i> = {res.aov1$p}"
  )

Stat$treat <- factor(as.numeric(gsub("[^0-9]", "", Stat$treat)), level = c(20, 100, 200))

pval_P1 <- df_p_val[1:2,]
pval_P1 <- pval_P1 %>% mutate(group1 = c("200", "100"), group2 = c("200", "100"))

P1 <-  ggplot(Stat, aes(x = treat, y = data1, group = col)) +
  geom_point(aes(col = col), position = position_dodge(width = 0.7)) +
  add_pvalue(
    pval_P1,
    label = "p = {p.adj}",
    tip.length = 0.03,
    fontface = "bold",
    y.position = 10,
    label.size = 2.75,
    step.increase = 0,
    bracket.shorten = -0.35,
    vjust = -0.5
  ) +
  xlab("") +
  scale_color_manual(values = colors) +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    panel.background = element_rect(fill = 'white'),
    legend.position = "none",
    plot.margin = unit(c(0, 0, -0.5, 0), 'cm'),
    plot.subtitle = element_markdown()
  ) +
  scale_y_continuous(expression("Toxin content (pg GDA cell " ^ -1 * ")"), breaks = seq(0, 12.5, 2.5)) +
  scale_x_discrete(labels = c("20", "100", "200")) +
  coord_cartesian(ylim = c(0, 14)) +
  labs(subtitle = statistic_label)

# Toxin normalization comparison: 

# add toxin per carbon, but only for Exponential phase
CN_2D2 <- CN %>% filter(strain == "L2-D2")
CN_2D2_vec <- c(rep(NA, each = 3), CN_2D2$C_nmol.cell[1:3], rep(NA, each = 3), CN_2D2$C_nmol.cell[4:6], rep(NA, each = 3))

Stat$C_mol_per_cell <- CN_2D2_vec * 10^-9

# calculate GDA:C M(GDA) = 768.941
Stat$GDA_C <- (Stat$data1 * 10^-12 / 768.941) / Stat$C_mol_per_cell

Pex1_2D2 <- ggplot(Stat, aes(y = data1, x =  factor(group))) +
  geom_point(aes(col = col)) +
  theme_classic() +
  scale_color_manual(values = colors)  +
  theme(
    legend.title = element_blank(),
    panel.background = element_rect(fill = 'white'),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, vjust = 0.5),
    plot.margin = unit(c(0, 0, 0, 0), 'cm'),
    plot.subtitle = element_markdown(),
    axis.title.x = element_markdown()
  ) +
  scale_x_discrete(labels = c(20, 100, 100, 200, 200)) +
  scale_y_continuous(expression("Toxin content (pg GDA cell " ^ -1 * ")")) +
  xlab("") +
  labs(subtitle = "<b>a)</b> per cell:")

Pex2_2D2 <- ggplot(Stat, aes(y = GDA_C, x = factor(group))) +
  geom_point(aes(col = col)) +
  theme_classic() +
  scale_color_manual(values = colors)  +
  theme(
    legend.title = element_blank(),
    panel.background = element_rect(fill = 'white'),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, vjust = 0.5),
    plot.margin = unit(c(0, 0, 0, 0), 'cm'),
    plot.subtitle = element_markdown(),
    axis.title.x = element_markdown()
  ) +
  scale_x_discrete(labels = c(20, 100, 100, 200, 200)) +
  scale_y_continuous(expression("Toxin content (GDA C" ^ -1 * " mol:mol)"),
                     label = scientific_format(digits = 2)) +
  xlab("") +
  labs(subtitle = "<b>a)</b> per carbon:")

Pex3_2D2 <- ggplot(Stat, aes(y = tox_vol, x = factor(group))) +
  geom_point(aes(col = col)) +
  theme_classic() +
  scale_color_manual(values = colors)  +
  theme(
    legend.title = element_blank(),
    panel.background = element_rect(fill = 'white'),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, vjust = 0.5),
    plot.margin = unit(c(0, 0, 0, 0), 'cm'),
    plot.subtitle = element_markdown(),
    axis.title.x = element_markdown()
  ) +
  scale_x_discrete(labels = c(20, 100, 100, 200, 200)) +
  scale_y_continuous(expression(paste("Toxin content (pg GDA ", mu, "m" ^ -3 *
                                        ")")), label = scientific_format(digits = 2)) +
  xlab("") +
  labs(subtitle = "<b>a)</b> per volume:")


# import data A. pseudogonyaulax strain L4-B1 #####
data <-
  read.csv(
    "C:\\Users\\krist\\OneDrive\\Dokumente\\AWI\\Promotion\\AP1\\Light-AP1b\\raw\\data_L4B1.txt",
    sep = "",
    header = FALSE
  )

data <- data[15:23, ]

# put first column as row names and transpose the data frame
data <- data %>% remove_rownames %>% column_to_rownames(var = "V1")
data <- as.data.frame(t(data))
data <- data[1:5, ]

#replace all , with .
data[] <- lapply(data, function(x)
  as.numeric(gsub(",", ".", x)))

# reshape dataframe 
data1 <-
  data.frame(
    data1 = unlist(data),
    treat = rep(
      c("20 Einstein", "100 Einstein", "200 Einstein"),
      each = 15,
      length.out = 45
    ),
    group = c(
      rep(c(1:5), length.out = 15),
      rep(c(6:10), length.out = 15),
      rep(c(11:15), length.out = 15)
    )
  )

data1$size <- c(
  rep(cell_size[1:3, 2], each = 5, length.out = 15),
  rep(cell_size[4, 3], length.out = 3),
  rep(cell_size[4, 6], length.out = 2),
  rep(cell_size[5, 3], length.out = 3),
  rep(cell_size[5, 6], length.out = 2),
  rep(cell_size[6, 3], length.out = 3),
  rep(cell_size[6, 6], length.out = 2),
  rep(cell_size[7, 3], length.out = 3),
  rep(cell_size[7, 6], length.out = 2),
  rep(cell_size[8, 3], length.out = 3),
  rep(cell_size[8, 6], length.out = 2),
  rep(cell_size[9, 3], length.out = 3),
  rep(cell_size[9, 6], length.out = 2)
)

# add toxin per volume
data1$tox_vol <- data1$data1/((4/3)*pi*(data1$size)^3)

# Statistical analysis
Stat2 <- data1 %>% filter(group == 1 | group == 6 | group == 9 | group == 11 | group == 14)

# Perform Kruskal-Wallis Test
res.aov2 <- kruskal_test(data1 ~ group, data = Stat2)

# Perform Conover-Iman posthoc test
Con_b <- conover.test(Stat2$data1, Stat2$group, method = "BH", altp = T)

# Calculate Cohen's effect size
effsize <- cohens_d(Stat2, data1 ~ group, hedges.correction = T)
effsize <-
  format(round(effsize[1:4, 4], 1) * -1, nsmall = 1, digits = 2)
effsize <- as.vector(abs(as.numeric(effsize$effsize)))

# introduce factor for Exponential / Stationary phase:
Stat2 <- Stat2 %>% arrange(group)
Stat2$col <-
  as.factor(c(rep("Exponential phase", each = 6), rep(
    c("Stationary phase", "Exponential phase"),
    each = 3,
    length.out = 9
  )))

text2 <- data.frame(
  group = unique(Stat2$group),
  effsize = c(NA, effsize),
  data = aggregate(Stat2, data1 ~ group, FUN = "max")$data,
  col = c(rep("Exponential phase", each = 2), rep(
    c("Stationary phase", "Exponential phase"),
    each = 1,
    length.out = 3
  )),
  GDA_sd = aggregate(Stat2, data1 ~ group, FUN = "sd")$data
)

# extract p-values < 0.05 and order for the plot
df_p_val2 <-
  data.frame(
    p.adj = sprintf("%.3f", Con_b$altP.adjusted),
    group1 = sub("^([0-9]+) .*", "\\1", Con_b$comparisons),
    group2 = sub(".* - ([0-9]+)", "\\1", Con_b$comparisons)
  )
df_p_val2 <- df_p_val2 %>% subset(p.adj < 0.05)

df_p_val2[] <-
  sapply(df_p_val2, function(x)
    gsub(pattern = " ", replacement = "", x))

df_p_val2$arrange <-
  as.numeric(df_p_val2$group2) - as.numeric(df_p_val2$group1)

df_p_val2 <- df_p_val2 %>% arrange(as.numeric(abs(arrange)))

df_p_val2$p.adj2 <-
  ifelse(df_p_val2$p.adj < 0.001, 0.001, df_p_val2$p.adj)

Stat2$group <- as.numeric(as.character(Stat2$group))

# Prepare ggplot parameters
dodge <- position_dodge(width = 0.9)

# Format the statistic with two digits after the decimal point
formatted_statistic <- sprintf("%.2f", res.aov2$statistic)

# Manually create the subtitle with markdown syntax using glue
statistic_label2 <- glue::glue("<b>b) </b>Kruskal-Wallis, <i>H</i><sub>{res.aov2$df}</sub> = {formatted_statistic}, <i>p</i> = {res.aov2$p}")

Stat2$treat <- factor(as.numeric(gsub("[^0-9]", "", Stat2$treat)), level = c(20, 100, 200))

pval_P2 <- df_p_val2[c(3, 2), ]
pval_P2 <- pval_P2 %>% mutate(group1 = c("100", "200"), group2 = c("100", "200"))

P2 <- ggplot(Stat2, aes(x = treat, y = data1, group = col)) +
  geom_point(aes(col = col), position = position_dodge(width = 0.7)) +
  add_pvalue(
    pval_P2[1, ],
    label = "p = {p.adj}",
    tip.length = 0.03,
    fontface = "bold",
    y.position = 9,
    label.size = 2.75,
    step.increase = 0,
    bracket.shorten = -0.35,
    vjust = -0.5
  ) +
  xlab("") +
  add_pvalue(
    pval_P2[2, ],
    label = "p < {p.adj2}",
    tip.length = 0.03,
    fontface = "bold",
    y.position = 9,
    label.size = 2.75,
    step.increase = 0,
    bracket.shorten = -0.35,
    vjust = -0.5
  ) +
  scale_color_manual(values = colors) +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    panel.background = element_rect(fill = 'white'),
    legend.position = "none",
    plot.margin = unit(c(0, 0, -0.5, 0), 'cm'),
    plot.subtitle = element_markdown()
  ) +
  scale_y_continuous(expression("Toxin content (pg GDA cell " ^ -1 * ")"), breaks = seq(0, 12.5, 2.5)) +
  scale_x_discrete(labels = c("20", "100", "200")) +
  coord_cartesian(ylim = c(0, 14)) +
  labs(subtitle = statistic_label2)

# Toxin comparison: 

# add toxin per carbon, but only for Exponential phase
CN_4B1 <- CN %>% filter(strain == "L4-B1")
CN_4B1_vec <- c(rep(NA, each = 3), CN_4B1$C_nmol.cell[1:3], rep(NA, each = 3), CN_4B1$C_nmol.cell[4:6], rep(NA, each = 3))

Stat2$C_mol_per_cell <- CN_4B1_vec * 10^-9

# calculate GDA:C M(GDA) = 768.941
Stat2$GDA_C <- (Stat2$data1 * 10^-12 / 768.941) / Stat2$C_mol_per_cell

Pex1_4B1 <- ggplot(Stat2, aes(y = data1, x =  factor(group))) +
  geom_point(aes(col = col)) +
  theme_classic() +
  scale_color_manual(values = colors) +
  theme(
    legend.title = element_blank(),
    panel.background = element_rect(fill = 'white'),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, vjust = 0.5),
    plot.margin = unit(c(0, 0, 0, 0), 'cm'),
    plot.subtitle = element_markdown(),
    axis.title.x = element_blank()
  ) +
  scale_x_discrete(labels = c(20, 100, 100, 200, 200)) +
  scale_y_continuous(expression("Toxin content (pg GDA cell " ^ -1 * ")")) +
  xlab("Photon flux densities (<i>&mu;</i>mol photons m<sup>-2</sup>s<sup>-1</sup>)") +
  labs(subtitle = "<b>b)</b> per cell:")

Pex2_4B1 <- ggplot(Stat2, aes(y = GDA_C, x = factor(group))) +
  geom_point(aes(col = col)) +
  theme_classic() +
  scale_color_manual(values = colors)+
  theme(
    legend.title = element_blank(),
    panel.background = element_rect(fill = 'white'),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, vjust = 0.5),
    plot.margin = unit(c(0, 0, 0, 0), 'cm'),
    plot.subtitle = element_markdown(),
    axis.title.x = element_blank()
  ) +
  scale_x_discrete(labels = c(20, 100, 100, 200, 200)) +
  scale_y_continuous(expression("Toxin content (GDA C" ^ -1 * " mol:mol)"),
                     label = scientific_format(digits = 2)) +
  xlab("Photon flux densities (<i>&mu;</i>mol photons m<sup>-2</sup>s<sup>-1</sup>)") +
  labs(subtitle = "<b>b)</b> per cell:")

Pex3_4B1 <- ggplot(Stat2, aes(y = tox_vol, x = factor(group))) +
  geom_point(aes(col = col)) +
  theme_classic() +
  scale_color_manual(values = colors) +
  theme(
    legend.title = element_blank(),
    panel.background = element_rect(fill = 'white'),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, vjust = 0.5),
    plot.margin = unit(c(0, 0, 0, 0), 'cm'),
    plot.subtitle = element_markdown(),
    axis.title.x = element_blank()
  ) +
  scale_x_discrete(labels = c(20, 100, 100, 200, 200)) +
  scale_y_continuous(expression(paste("Toxin content (pg GDA ", mu, "m" ^ -3 *
                                        ")")), label = scientific_format(digits = 2)) +
  xlab("Photon flux densities (<i>&mu;</i>mol photons m<sup>-2</sup>s<sup>-1</sup>)") +
  labs(subtitle = "<b>b)</b> per volume:")

# import data A. pseudogonyaulax strain L4-B9 #####
data <-
  read.csv(
    "C:\\Users\\krist\\OneDrive\\Dokumente\\AWI\\Promotion\\AP1\\Light-AP1b\\raw\\data_L4B9.txt",
    sep = "",
    header = FALSE,
    skip = 15,
    nrow = 9
  )

# put first column as row names and transpose the data frame
data <- data %>% remove_rownames %>% column_to_rownames(var = "V1")
data <- as.data.frame(t(data))
data <- data[1:4, ]

# replace all , with . and change to numeric
data[] <- lapply(data, function(x)
  as.numeric(gsub(",", ".", x)))

# reshape dataframe 
data1 <-
  data.frame(
    data1 = unlist(data),
    treat = rep(
      c("20 Einstein", "100 Einstein", "200 Einstein"),
      each = 12,
      length.out = 36
    ),
    group = c(
      rep(c(1:4), length.out = 12),
      rep(c(5:8), length.out = 12),
      rep(c(9:12), length.out = 12)
    )
  )

# add cell sizes column, ex = Exponential and st = Stationary phase 
data1$size <- c(
  rep(cell_size[1:3, 2], each = 4, length.out = 12),
  rep(cell_size[4, 4], length.out = 3),
  rep(cell_size[4, 7], length.out = 1),
  rep(cell_size[5, 4], length.out = 3),
  rep(cell_size[5, 7], length.out = 1),
  rep(cell_size[6, 4], length.out = 3),
  rep(cell_size[6, 7], length.out = 1),
  rep(cell_size[7, 4], length.out = 3),
  rep(cell_size[7, 7], length.out = 1),
  rep(cell_size[8, 4], length.out = 3),
  rep(cell_size[8, 7], length.out = 1),
  rep(cell_size[9, 4], length.out = 3),
  rep(cell_size[9, 7], length.out = 1)
)

# add toxin per volume
data1$tox_vol <- data1$data1/((4/3)*pi*(data1$size)^3)

# Statistical analysis
# Subset of first Exponential and Stationary phase sampling point for manuscript
Stat3 <- data1 %>% filter(group == 1 |
           group == 5 | group == 8 | group == 9 | group == 12)

# Perform Kruskal-Wallis Test
res.aov3 <- kruskal_test(data1 ~ group, data = Stat3)

# Perform Conover-Iman posthoc test
Con_c <- conover.test(Stat3$data1, Stat3$group, method = "BH", altp = T)

# Calculate Cohen's effect size
effsize <- cohens_d(Stat3, data1 ~ group, hedges.correction = T)
effsize <-
  format(round(effsize[1:4, 4], 1) * -1, nsmall = 1, digits = 2)
effsize <-
  sprintf("%.1f", as.vector(abs(as.numeric(effsize$effsize))))

# introduce factor for Exponential / Stationary phase:
Stat3 <- Stat3 %>% arrange(group)
Stat3$col <-
  c(rep("Exponential phase", each = 6), rep(
    c("Stationary phase", "Exponential phase"),
    each = 3,
    length.out = 9
  ))

text3 <- data.frame(
  group = unique(Stat3$group),
  effsize = c(NA, effsize),
  data = aggregate(Stat3, data1 ~ group, FUN = "max")$data,
  col = c(rep("Exponential", each = 2), rep(
    c("Stationary", "Exponential"),
    each = 1,
    length.out = 3
  )),
  GDA_sd = aggregate(Stat2, data1 ~ group, FUN = "sd")$data
)

# extract p-values < 0.05 and order for the plot
df_p_val3 <-
  data.frame(
    p.adj = sprintf("%.3f", Con_c$altP.adjusted),
    group1 = sub("^([0-9]+) .*", "\\1", Con_c$comparisons),
    group2 = sub(".* - ([0-9]+)", "\\1", Con_c$comparisons)
  )
df_p_val3 <- df_p_val3 %>% subset(p.adj < 0.05)

df_p_val3[] <-
  sapply(df_p_val3, function(x)
    gsub(pattern = " ", replacement = "", x))
df_p_val3$arrange <-
  as.numeric(df_p_val3$group2) - as.numeric(df_p_val3$group1)
df_p_val3 <- df_p_val3 %>% arrange(as.numeric(arrange))

Stat3$group <- as.numeric(as.character(Stat3$group))
Stat3$col <- factor(Stat3$col)

# Format the statistic with two digits after the decimal point
formatted_statistic <- sprintf("%.2f", res.aov3$statistic)

# Manually create the subtitle with markdown syntax using glue
statistic_label3 <-
  glue::glue(
    "<b>c) </b>Kruskal-Wallis, <i>H</i><sub>{res.aov3$df}</sub> = {formatted_statistic}, <i>p</i> = {res.aov3$p}"
  )

Stat3$treat <- factor(as.numeric(gsub("[^0-9]", "", Stat3$treat)), level = c(20, 100, 200))

pval_P3 <- df_p_val3[c(4, 2), ]
pval_P3 <- pval_P3 %>% mutate(group1 = c("100", "200"), group2 = c("100", "200"))
  
  P3 <-  ggplot(Stat3, aes(x = treat, y = data1, group = col)) +
  geom_point(aes(col = col), position = position_dodge(width = 0.7)) +
  add_pvalue(
    pval_P3[, ],
    label = "p = {p.adj}",
    tip.length = 0.03,
    fontface = "bold",
    y.position = 13,
    label.size = 2.75,
    step.increase = 0,
    bracket.shorten = -0.35,
    vjust = -0.5
  ) +
  scale_color_manual(values = colors) +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    panel.background = element_rect(fill = 'transparent'),
    legend.position = "none",
    plot.margin = unit(c(0, 0, 1, 0), 'lines'),
    plot.subtitle = element_markdown(),
    axis.title.x = element_markdown()
  ) +
    xlab("Photon flux densities (<i>&mu;</i>mol photons m<sup>-2</sup>s<sup>-1</sup>)") +
  scale_y_continuous(expression("Toxin content (pg GDA cell " ^ -1 * ")"), breaks = seq(0, 12.5, 2.5)) +
  scale_x_discrete(labels = c(20,  100, 200)) +
  coord_cartesian(ylim = c(0, 14), clip = 'off') +
  labs(subtitle = statistic_label3)
  

# Toxin comparison: 
# add toxin per carbon, but only for Exponential phase
CN_4B9 <- CN %>% filter(strain == "L4-B9")
CN_4B9_vec <- c(rep(NA, each = 3), CN_4B9$C_nmol.cell[1:3], rep(NA, each = 3), CN_4B9$C_nmol.cell[4:6], rep(NA, each = 3))

Stat3$C_mol_per_cell <- CN_4B9_vec * 10^-9

# calculate GDA:C M(GDA) = 768.941
Stat3$GDA_C <- (Stat3$data1 * 10^-12 / 768.941) / Stat3$C_mol_per_cell

Pex1_4B9 <- ggplot(Stat3, aes(y = data1, x =  factor(group))) +
  geom_point(aes(col = col)) +
  theme_classic() +
  scale_color_manual(values = colors)  +
  theme(
    legend.title = element_blank(),
    panel.background = element_rect(fill = 'white'),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, vjust = 0.5),
    plot.margin = unit(c(0, 0, -0, 0), 'cm'),
    plot.subtitle = element_markdown(),
    axis.title.x = element_markdown()
  ) +
  scale_x_discrete(labels = c(20, 100, 100, 200, 200)) +
  scale_y_continuous(expression("Toxin content (pg GDA cell " ^ -1 * ")")) +
  xlab("Photon flux densities (<i>&mu;</i>mol photons m<sup>-2</sup>s<sup>-1</sup>)") +
  labs(subtitle = "<b>c)</b> per cell:")

Pex2_4B9 <- ggplot(Stat3, aes(y = GDA_C, x = factor(group))) +
  geom_point(aes(col = col)) +
  theme_classic() +
  scale_color_manual(values = colors)  +
  theme(
    legend.title = element_blank(),
    panel.background = element_rect(fill = 'white'),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, vjust = 0.5),
    plot.margin = unit(c(0, 0, -0, 0), 'cm'),
    plot.subtitle = element_markdown(),
    axis.title.x = element_markdown()
  ) +
  scale_x_discrete(labels = c(20, 100, 100, 200, 200)) +
  scale_y_continuous(expression("Toxin content (GDA C" ^ -1 * " mol:mol)"),
                     label = scientific_format(digits = 2)) +
  xlab("Photon flux densities (<i>&mu;</i>mol photons m<sup>-2</sup>s<sup>-1</sup>)") +
  labs(subtitle = "<b>c)</b> per carbon:")

Pex3_4B9 <- ggplot(Stat3, aes(y = tox_vol, x = factor(group))) +
  geom_point(aes(col = col)) +
  theme_classic() +
  scale_color_manual(values = colors)  +
  theme(
    legend.title = element_blank(),
    panel.background = element_rect(fill = 'white'),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, vjust = 0.5),
    plot.margin = unit(c(0, 0, -0, 0), 'cm'),
    plot.subtitle = element_markdown(),
    axis.title.x = element_markdown()
  ) +
  scale_x_discrete(labels = c(20, 100, 100, 200, 200)) +
  scale_y_continuous(expression(paste("Toxin content (pg GDA ", mu, "m" ^ -3 *
                                        ")")), label = scientific_format(digits = 2)) +
  xlab("Photon flux densities (<i>&mu;</i>mol photons m<sup>-2</sup>s<sup>-1</sup>)") +
  labs(subtitle = "<b>c)</b> per volume:")

# combine toxin plots of all strains in one plot and export ################################
all_plots <-
  P1 + P2 + P3 + plot_layout(guides = "collect", ncol = 1, axis_titles = "collect_y") &
  theme(legend.position = "bottom",
        legend.box.margin = margin(-10, -10, -10, -10))

# ggsave(
#   "Toxins_all_normiert.png",
#   all_plots,
#   path = "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/AP1/Light-AP1b/figures",
#   dpi = 300,
#   width = 3.5,
#   height = 4.9,
#   units = "in"
# )

# combine toxin normalization comparisons of all strains in one plot and export ################################

# combine toxin normalization comparisons of both strains in one plot  ################################
# Combine Pex1 with Pex1_4B1
plot1_combined <- Pex1_2D2 + Pex1_4B1 + Pex1_4B9 + plot_layout(guides = "collect", ncol = 1, axis_titles = "collect_y") &
  theme(legend.position = "bottom",
        legend.box.margin = margin(-10, -10, -10, -10))

# Combine Pex3 with Pex3_4B1
plot2_combined <- Pex3_2D2 + Pex3_4B1 + Pex3_4B9 + plot_layout(guides = "collect", ncol = 1, axis_titles = "collect_y")  &
  theme(legend.position = "bottom",
        legend.box.margin = margin(-10, -10, -10, -10))

# Combine Pex2 with Pex2_4B1
# plot3_combined <- Pex2_2D2 + Pex2_4B1 + Pex2_4B9 + plot_layout(guides = "collect", ncol = 1, axis_titles = "collect_y")  &
#   theme(legend.position = "none",
#         axis.title.x = element_blank(),
#         legend.box.margin = margin(-10, -10, -10, -10))

# Arrange all combined plots in a grid layout
all_plots_comparison <- wrap_plots(plot1_combined, plot2_combined, ncol = 2)

ggsave(
  "Toxins_comparison.png",
  all_plots_comparison ,
  path = "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/AP1/Light-AP1b/figures",
  dpi = 300,
  width = 20,
  height = 20,
  units = "cm"
)

# Between-strains statistics ################################

# combine toxin quotas of all strains in one dataframe
Stat <- Stat %>% mutate(strain = "L2-D2") %>% convert_as_factor(treat, col, strain)
Stat2 <- Stat2 %>% mutate(strain = "L4-B1") %>% convert_as_factor(treat, col, strain)
Stat3 <- Stat3 %>% mutate(strain = "L4-B9") %>% convert_as_factor(treat, col, strain)

Stat_all <- rbind(Stat, Stat2, Stat3) 

Stat_all$treat <- factor(Stat_all$treat, levels = c("20", "100", "200"))

# Perform Kruskal-Wallis and Conover-Iman posthoc test on all light intensities between the three strains

Tox_all <- Stat_all %>%
  group_by(treat, col) %>%
  nest() %>%
  mutate(kruskal_results = map(data, ~kruskal_test(data1 ~ strain, data = .x))) %>%
  unnest(c(kruskal_results)) %>%
  mutate(conover_results = map(data, ~ conover.test(.x$data1, .x$strain, altp = T, method = "BH"))) %>%
  unnest(c(conover_results))

# statistical analysis of Cell size in between strains

cell_size2 <-
  pivot_longer(
    cell_size,
    cols = c(2:7),
    names_to = "strain_gp",
    values_to = "data"
  ) %>% mutate(strain = str_sub(strain_gp, 1, 4), gp = str_sub(strain_gp, 6, 7)) %>% convert_as_factor(strain, gp)

cell_size$strain <- as.factor(str_sub(cell_size$treat2, 1, 4))
cell_size$gp <- as.factor(str_extract(cell_size$treat2, "(?<=_)[^_]+(?=_[0-9]+)"))
cell_size$light <- factor(str_extract(cell_size$treat, "([0-9]+)(?=_\\d)"), levels = c("20", "100", "200"))

cell_size_all <- cell_size2 %>%
  drop_na(data) %>%
  filter(treat  !=  20 | gp != "st") %>%
  group_by(strain, gp) %>%
  nest() %>%
  mutate(kruskal_results = map(data, ~kruskal_test(data ~ treat, data = .x))) %>%
  unnest(c(kruskal_results)) %>%
  mutate(conover_results = map(data, ~ conover.test(.x$data, .x$treat, altp = T, method = "BH"))) %>%
  unnest(c(conover_results))

# # Molar ratio POC/PONand GDA:C calculation + export as table
# CN <-  read.csv("C:\\Users\\krist\\OneDrive\\Dokumente\\AWI\\Promotion\\AP1\\Light-AP1b\\raw\\POC_PON_AP1_light.txt", sep="", header=T)
#
# CN$ratio <- CN$C_nmol.cell/CN$N_nmol.cell
# CN$GDA_pg_cell <- c(Stat$data1[c(4:6, 10:12)], Stat2$data1[c(4:6, 10:12)], Stat3$data1[c(4:6, 10:12)])
# CN$GDA_fmol_cell <- format(round(CN$GDA_pg_cell / 768.9410 *1000, 2), nsmall=2, digits=2)
# CN$GDA_per_C <- signif((as.numeric(CN$GDA_fmol_cell) * 10^-15) / (as.numeric(CN$C_nmol.cell) * 10^-9), 3)
#
# CN_agg <- aggregate(data=CN, ratio ~ strain + treat, FUN="mean")
# CN_agg$sd <- aggregate(data=CN, ratio ~ strain + treat, FUN="sd")$ratio
#
# CN_agg$GDA_per_C <- signif(aggregate(data=CN, GDA_per_C ~ strain + treat, FUN="mean")$GDA_per_C, 3)
# CN_agg$GDA_per_C_sd <- signif(aggregate(data=CN, GDA_per_C ~ strain + treat, FUN="sd")$GDA_per_C, 3)
#
# CN_agg$ratio2 <- paste(format(round(CN_agg$ratio, 2), nsmall=2, digits=2), format(round(CN_agg$sd, 2), nsmall=2, digits=2), sep=" +/- ")
# CN_agg$GDA_per_C_total <- paste(signif(CN_agg$GDA_per_C, 3), signif(CN_agg$GDA_per_C_sd, 3), sep=" +/- ")
#
# # Export as table
# CN_agg <- flextable(CN_agg) %>% autofit() %>%  save_as_docx(path="C:\\Users\\krist\\OneDrive\\Dokumente\\AWI\\Promotion\\AP1\\Light-AP1b\\POC_PON.docx")

# Garbage collection: call after large objects have been removed

gc()

# Delete workspace, clean environment/variables and restarte R if used memory piles up

dev.off()

rm(list = ls())

.rs.restartR()
