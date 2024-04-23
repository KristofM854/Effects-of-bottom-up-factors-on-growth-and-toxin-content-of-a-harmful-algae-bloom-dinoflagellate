##########################################
## Growth rate calculations and plots of A. pseudogonyaulax strains L2-D2, L4-B1 and L4-B9 exposed to different light intensities
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
  ggplot2,
  dplyr,
  tidyverse,
  rstatix,
  ggprism,
  ggpubr,
  FME,
  stats,
  latex2exp,
  growthrates,
  gridExtra,
  ggpmisc,
  grid,
  conover.test,
  patchwork,
  ggthemes,
  ggtext
)

# Load and transform data ########
data_2D2 <-
  read.csv(
    "C:\\Users\\krist\\OneDrive\\Dokumente\\AWI\\Promotion\\AP1\\Light-AP1b\\raw\\data_L2D2.txt",
    sep = "",
    header = FALSE,
    nrows = 22
  )

# extract maximum cell densities for the plot
data2 <-
  data.frame(treat = as.factor(rep(c("20", "100", "200"), each = 3)),
             data = c(
               as.numeric(data_2D2[3:5, 13]),
               as.numeric(data_2D2[6:8, 11]),
               as.numeric(data_2D2[9:11, 11])
             ))

# restrict data.frame to cell counts; remaining data is toxin quota and cell densities
data_2D2 <- data_2D2[2:11, ]

# put first column as row names and transpose the data frame
data_2D2 <-
  data_2D2 %>% remove_rownames %>% column_to_rownames(var = "V1")

data_2D2 <- data_2D2 %>% t() %>% as.data.frame()

# replace all , with . + make columns numeric
data_2D2[] <- lapply(data_2D2, function(x)
  gsub("NaN", "NA", x))

data_2D2[] <- lapply(data_2D2, function(x)
  as.numeric(gsub(",", ".", x)))

# Growth rate calculations ######
# make linear fit for each bottle between 5 data points
# h: amount of data points used for calculation of growth rate!!
fit_20_1 <-
  fit_easylinear((data_2D2[c(1:10, 12), 1]), (data_2D2[c(1:10, 12), 2]), h = 5)
fit_20_2 <-
  fit_easylinear((data_2D2[c(1:10, 12), 1]), (data_2D2[c(1:10, 12), 3]), h = 5)
fit_20_3 <-
  fit_easylinear((data_2D2[c(1:10, 12), 1]), (data_2D2[c(1:10, 12), 4]), h = 5)
fit_100_1 <-
  fit_easylinear(((data_2D2[1:11, 1])), (data_2D2[1:11, 5]), h = 5)
fit_100_2 <-
  fit_easylinear(((data_2D2[1:11, 1])), (data_2D2[1:11, 6]), h = 5)
fit_100_3 <-
  fit_easylinear(((data_2D2[1:11, 1])), (data_2D2[1:11, 7]), h = 5)
fit_200_1 <-
  fit_easylinear(((data_2D2[1:11, 1])), (data_2D2[1:11, 8]), h = 5)
fit_200_2 <-
  fit_easylinear(((data_2D2[1:11, 1])), (data_2D2[1:11, 9]), h = 5)
fit_200_3 <-
  fit_easylinear(((data_2D2[1:11, 1])), (data_2D2[1:11, 10]), h = 5)

# extract coefficients of the linear fits 
matrix_coef <-
  as.data.frame(rbind(
    coef(fit_20_1),
    coef(fit_20_2),
    coef(fit_20_3),
    coef(fit_100_1),
    coef(fit_100_2),
    coef(fit_100_3),
    coef(fit_200_1),
    coef(fit_200_2),
    coef(fit_200_3)
  ))

# construct new dataframe 'ex' with maximal growth rate and treatments for following analysis
ex = data.frame(value = matrix_coef$mumax, group = rep(factor(
  c("20", "100", "200"), levels = c("20", "100", "200")
), each = 3))

# Statistical Analysis  ######
# Perform a non-parametric Kruskal test since three samples do not satisfy the normality assumption
# Null hypothesis (H0): The growth rate of A. pseudogonyaulax exposed to different nitrogen sources is equal
res.aov1 <- kruskal_test(value ~ group, data = ex)

# Conduct a Conover-Iman non-parametric PostHoc test if p < 0.05, with a liberal p-value adjustment using Benjamini & Hochberg (FDR)
Con_a <- conover.test(ex$value, ex$group, method = "BH", altp = T)

# Extract p-values and group comparisons from the Conover-Iman test output
df_p_val <-
  data.frame(
    p.adj = sprintf("%.3f", Con_a$altP.adjusted),
    group1 = as.factor(str_sub(Con_a$comparisons, 1, 3)),
    group2 = as.factor(str_sub(Con_a$comparisons, 6, 9))
  )

df_p_val <- df_p_val %>% subset(p.adj < 0.1)

df_p_val[] <-
  sapply(df_p_val, function(x)
    gsub(pattern = " ", replacement = "", x))

df_p_val$arrange <-
  as.numeric(df_p_val$group2) - as.numeric(df_p_val$group1)

# swap the first comparison because otherwise the bracket.shorten attribute from add_pvalue somehow does not work properly ...
df_p_val[1, 2] <- 20
df_p_val[1, 3] <- 100

# Calculate Cohen's d effect size: ~0.2 (small effect), ~0.5 (medium effect), ~0.8 (large effect)
effsize <- cohens_d(ex, value ~ group, hedges.correction = T)

effsize_sub <-
  abs(as.numeric(format(
    round(effsize$effsize[c(1:2)], digits = 1), nsmall = 1
  )))

# Check for significant differences in cell densities
# Null hypothesis (H0): The maximum cell densities of A. pseudogonyaulax exposed to different nitrogen sources are equal
res.aov2 <- kruskal_test(data ~ treat, data = data2)

# Perform a Conover-Iman non-parametric PostHoc test if p < 0.05, with a liberal p-value adjustment using Benjamini & Hochberg (FDR)
Con_a2 <-
  conover.test(data2$data, data2$treat, method = "BH", altp = T)

# Prepare ggplot parameters  ######
dodge <- position_dodge(width = 0.4)

text <-
  data.frame(
    value = aggregate(ex, value ~ group, FUN = "max")$value,
    group = unique(ex$group),
    effsize = c(NA, effsize_sub)
  )

max_dens <- aggregate(data2, data ~ treat, FUN = "max")
max_dens$sd <- aggregate(data2, data ~ treat, FUN = "sd")$data

# Format the statistic with two digits after the decimal point
formatted_statistic <- sprintf("%.2f", res.aov1$statistic)

# Create subtitle with markdown syntax
statistic_label <-
  glue::glue(
    "<b>a) </b>Kruskal-Wallis, <i>H</i><sub>{res.aov1$df}</sub> = {formatted_statistic}, <i>p</i> = {res.aov1$p}"
  )
# L2D2 ggplot  ######
P1 <- ggplot(ex, aes(y = value, x = group)) +
  geom_point() +
  add_pvalue(
    df_p_val[c(1, 2), ],
    label = "p = {p.adj}",
    tip.length = 0.03,
    fontface = "bold",
    y.position = 0.4,
    label.size = 2.75,
    step.increase = 0,
    bracket.shorten = 0.1,
    vjust = -0.3
  ) +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    panel.background = element_rect(fill = 'white'),
    legend.position = "bottom",
    plot.margin = unit(c(0, 0, -0.5, 0), 'cm'),
    plot.subtitle = element_markdown()
  ) +
  xlab(expression(paste(""))) +
  scale_y_continuous(expression("Growth rate (d" ^ -1 * ")")) +
  coord_cartesian(ylim = c(0, 0.5))  +
  labs(subtitle = statistic_label)

#### Add A. pseudogonyaulax strain L4-B1  ######

# Load and transform data ########
data_4B1 <-
  read.csv(
    "C:\\Users\\krist\\OneDrive\\Dokumente\\AWI\\Promotion\\AP1\\Light-AP1b\\raw\\data_L4B1.txt",
    sep = "",
    header = FALSE
  )

#replace all , with . + make columns numeric
data_4B1[, c(2:13)] <-
  lapply(data_4B1[, c(2:13)], function(x)
    as.numeric(gsub(",", ".", x)))

# extract maximum cell densities for the plot
data3 <-
  data.frame(treat = as.factor(rep(c("20", "100", "200"), each = 3)),
             data = c(
               as.numeric(data_4B1[4:6, 5]),
               as.numeric(data_4B1[7:9, 12]),
               as.numeric(data_4B1[10:12, 12])
             ))

# restrict data.frame to cell counts; remaining data is toxin quota and cell densities
data_4B1 <- data_4B1[2:12, ]

# put first column as row names and transpose the data frame
data_4B1 <-
  data_4B1 %>% remove_rownames %>% column_to_rownames(var = "V1")

data_4B1 <- as.data.frame(t(data_4B1))

## Growth rate calculations ######
# make linear fit for each bottle between, if possible, 5 data points
# h: amount of data points used for calculation of growth rate!!
fit_60_1 <-
  fit_easylinear((data_4B1$Day_20), (data_4B1$`20_umol_A`), h = 5)
fit_60_2 <-
  fit_easylinear((data_4B1$Day_20), (data_4B1$`20_umol_B`), h = 5)
fit_60_3 <-
  fit_easylinear((data_4B1$Day_20), (data_4B1$`20_umol_C`), h = 5)
fit_100_1 <-
  fit_easylinear(((data_4B1[1:11, 1])), (data_4B1[1:11, 6]), h = 5)
fit_100_2 <-
  fit_easylinear(((data_4B1[1:11, 1])), (data_4B1[1:11, 7]), h = 5)
fit_100_3 <-
  fit_easylinear(((data_4B1[1:11, 1])), (data_4B1[1:11, 8]), h = 5)
fit_200_1 <-
  fit_easylinear(((data_4B1[1:11, 1])), (data_4B1[1:11, 9]), h = 5)
fit_200_2 <-
  fit_easylinear(((data_4B1[1:11, 1])), (data_4B1[1:11, 10]), h = 5)
fit_200_3 <-
  fit_easylinear(((data_4B1[1:11, 1])), (data_4B1[1:11, 11]), h = 5)
fit_20_1 <-
  fit_easylinear((data_4B1[1:8, 2]), (data_4B1[1:8, 3]), h = 3)
fit_20_2 <-
  fit_easylinear((data_4B1[1:8, 2]), (data_4B1[1:8, 4]), h = 3)
fit_20_3 <-
  fit_easylinear((data_4B1[1:8, 2]), (data_4B1[1:8, 5]), h = 3)

# extract coefficients of the linear fits 
matrix_coef <-
  as.data.frame(rbind(
    coef(fit_20_1),
    coef(fit_20_2),
    coef(fit_20_3),
    coef(fit_60_1),
    coef(fit_60_2),
    coef(fit_60_3),
    coef(fit_100_1),
    coef(fit_100_2),
    coef(fit_100_3),
    coef(fit_200_1),
    coef(fit_200_2),
    coef(fit_200_3)
  ))

# construct new dataframe 'ex2' with maximal growth rate and treatments for following analysis
ex2 = data.frame(value = matrix_coef$mumax, group = rep(factor(
  c("20", "60", "100", "200"), levels = c("20", "60", "100", "200")
), each = 3))

# Statistical Analysis  ######
# Perform a non-parametric Kruskal test since three samples do not satisfy the normality assumption
# Null hypothesis (H0): The growth rate of A. pseudogonyaulax exposed to different nitrogen sources is equal
res.aov2 <- kruskal_test(value ~ group, data = ex2)

# Conduct a Conover-Iman non-parametric PostHoc test if p < 0.05, with a liberal p-value adjustment using Benjamini & Hochberg (FDR)
Con_b <- conover.test(ex2$value, ex2$group, altp = T, method = "BH")

# Extract p-values and group comparisons from the Conover-Iman test output
df_p_val2 <-
  data.frame(
    p.adj = sprintf("%.3f", Con_b$altP.adjusted),
    group1 = str_sub(Con_b$comparisons, 1, 3),
    group2 = str_sub(Con_b$comparisons, 6, 9)
  )

df_p_val2 <- df_p_val2 %>% subset(p.adj < 0.1)

df_p_val2[] <-
  sapply(df_p_val2, function(x)
    gsub(pattern = " ", replacement = "", x))
df_p_val2$arrange <-
  as.numeric(df_p_val2$group2) - as.numeric(df_p_val2$group1)
df_p_val2 <- df_p_val2 %>% arrange(as.numeric(arrange))

df_p_val2[3, c(2:3)] <- df_p_val2[3, c(3:2)]
df_p_val2[c(3:4), ] <- df_p_val2[c(4:3), ]

# Calculate Cohen's d effect size: ~0.2 (small effect), ~0.5 (medium effect), ~0.8 (large effect)
effsize2 <- cohens_d(ex2, value ~ group, hedges.correction = T)

effsize_sub <- effsize2[c(1:3), 4]

effsize_sub <-
  format(round(effsize_sub * -1, digits = 1), nsmall = 1)

ex2$effsize <-
  c(rep(NA, each = 5),
    effsize_sub[1, 1],
    NA,
    NA,
    effsize_sub[2, 1],
    NA,
    NA,
    effsize_sub[3, 1])

text2 <-
  data.frame(
    value = aggregate(ex2, value ~ group, FUN = "max")$value,
    group = unique(ex2$group),
    effsize = c(NA, effsize_sub$effsize)
  )

# Prepare ggplot parameters
dodge <- position_dodge(width = 0.4)

max_dens2 <- aggregate(data3, data ~ treat, FUN = "max")
max_dens2$sd <- aggregate(data3, data ~ treat, FUN = "sd")$data

# Format the statistic with two digits after the decimal point
formatted_statistic <- sprintf("%.2f", res.aov2$statistic)

# Create subtitle with markdown syntax
statistic_label2 <-
  glue::glue(
    "<b>b) </b>Kruskal-Wallis, <i>H</i><sub>{res.aov2$df}</sub> = {formatted_statistic}, <i>p</i> = {res.aov2$p}"
  )

# L4-B1 ggplot  ######
P2 <- ggplot(ex2, aes(y = value, x = group)) +
  geom_point() +
  add_pvalue(
    df_p_val2[3:5, ],
    label = "p = {p.adj}",
    tip.length = 0.03,
    fontface = "bold",
    y.position = 0.425,
    label.size = 2.75,
    step.increase = 0,
    bracket.shorten = 0.1,
    vjust = -0.3
  ) +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    panel.background = element_rect(fill = 'white'),
    legend.position = "bottom",
    plot.margin = unit(c(-0.5, 0, -0.5, 0), 'cm'),
    plot.subtitle = element_markdown()
  ) +
  xlab(expression(paste(""))) +
  scale_y_continuous(expression("Growth rate (d" ^ -1 * ")")) +
  coord_cartesian(ylim = c(0, 0.5)) +
  labs(subtitle = statistic_label2)

#### Add A. pseudogonyaulax strain L4-B9  ######

# Load and transform data ########
data_4B9 <-
  read.csv(
    "C:\\Users\\krist\\OneDrive\\Dokumente\\AWI\\Promotion\\AP1\\Light-AP1b\\raw\\data_L4B9.txt",
    sep = "",
    header = FALSE
  )

# extract maximum cell densities for the plot
data4 <-
  data.frame(treat = as.factor(rep(c("20", "100", "200"), each = 3)),
             data = c(
               as.numeric(data_4B9[4:6, 5]),
               as.numeric(data_4B9[7:9, 12]),
               as.numeric(data_4B9[10:12, 12])
             ))

# restrict data.frame to cell counts; remaining data is toxin quota and cell densities
data_4B9 <- data_4B9[2:12,]

# put first column as row names and transpose the data frame
data_4B9 <-
  data_4B9 %>% remove_rownames %>% column_to_rownames(var = "V1")

data_4B9 <- as.data.frame(t(data_4B9))

#replace all , with . + make columns numeric
data_4B9[] <- lapply(data_4B9, function(x)
  as.numeric(gsub(",", ".", x)))

## Growth rate calculations ######
# make linear fit for each bottle between 5 data points
# h: amount of data points used for calculation of growth rate!!
fit_60_1 <-
  fit_easylinear((data_4B9[1:12, 2]), (data_4B9[1:12, 3]), h = 5)
fit_60_2 <-
  fit_easylinear((data_4B9[1:12, 2]), (data_4B9[1:12, 4]), h = 5)
fit_60_3 <-
  fit_easylinear((data_4B9[1:12, 2]), (data_4B9[1:12, 5]), h = 5)
fit_100_1 <-
  fit_easylinear((data_4B9$Day_100_200), (data_4B9$`100_umol_A`), h = 5)
fit_100_2 <-
  fit_easylinear((data_4B9$Day_100_200), (data_4B9$`100_umol_B`), h = 5)
fit_100_3 <-
  fit_easylinear((data_4B9$Day_100_200), (data_4B9$`100_umol_C`), h = 5)
fit_200_1 <-
  fit_easylinear((data_4B9$Day_100_200), (data_4B9$`200_umol_A`), h = 5)
fit_200_2 <-
  fit_easylinear((data_4B9$Day_100_200), (data_4B9$`200_umol_B`), h = 5)
fit_200_3 <-
  fit_easylinear((data_4B9$Day_100_200), (data_4B9$`200_umol_C`), h = 5)
fit_20_1 <-
  fit_easylinear((data_4B9[1:8, 2]), (data_4B9[1:8, 3]), h = 5)
fit_20_2 <-
  fit_easylinear((data_4B9[1:8, 2]), (data_4B9[1:8, 4]), h = 5)
fit_20_3 <-
  fit_easylinear((data_4B9[1:8, 2]), (data_4B9[1:8, 5]), h = 5)

# extract coefficients of the linear fits 
matrix_coef <-
  as.data.frame(rbind(
    coef(fit_20_1),
    coef(fit_20_2),
    coef(fit_20_3),
    coef(fit_60_1),
    coef(fit_60_2),
    coef(fit_60_3),
    coef(fit_100_1),
    coef(fit_100_2),
    coef(fit_100_3),
    coef(fit_200_1),
    coef(fit_200_2),
    coef(fit_200_3)
  ))

# construct new dataframe 'ex' with maximal growth rate and treatments for following analysis
ex3 = data.frame(value = matrix_coef$mumax, group = rep(factor(
  c("20", "60", "100", "200"), levels = c("20", "60", "100", "200")
), each = 3))

# Statistical Analysis  ######
# Perform a non-parametric Kruskal test since three samples do not satisfy the normality assumption
# Null hypothesis (H0): The growth rate of A. pseudogonyaulax exposed to different nitrogen sources is equal
res.aov3 <- kruskal_test(value ~ group, data = ex3)

# Conduct a Conover-Iman non-parametric PostHoc test if p < 0.05, with a liberal p-value adjustment using Benjamini & Hochberg (FDR)
Con_c <- conover.test(ex3$value, ex3$group, altp = T, method = "BH")

# Extract p-values and group comparisons from the Conover-Iman test output
df_p_val3 <-
  data.frame(
    p.adj = sprintf("%.3f", Con_c$altP.adjusted),
    group1 = str_sub(Con_c$comparisons, 1, 3),
    group2 = str_sub(Con_c$comparisons, 6, 9)
  )
df_p_val3 <- df_p_val3 %>% subset(p.adj < 0.1)

df_p_val3[] <-
  sapply(df_p_val3, function(x)
    gsub(pattern = " ", replacement = "", x))
df_p_val3$arrange <-
  as.numeric(df_p_val3$group2) - as.numeric(df_p_val3$group1)
df_p_val3 <- df_p_val3 %>% arrange(as.numeric(arrange))

df_p_val3[3, c(2:3)] <- df_p_val3[3, c(3:2)]
df_p_val3[c(3:4),] <- df_p_val3[c(4:3),]

# Calculate Cohen's d effect size: ~0.2 (small effect), ~0.5 (medium effect), ~0.8 (large effect)
effsize3 <- cohens_d(ex3, value ~ group, hedges.correction = T)
effsize_sub <- effsize3[c(1:3), 4]

effsize_sub <-
  format(round(effsize_sub * -1, digits = 1), nsmall = 1)
ex3$effsize <-
  c(rep(NA, each = 5),
    effsize_sub[1, 1],
    NA,
    NA,
    effsize_sub[2, 1],
    NA,
    NA,
    effsize_sub[3, 1])

text3 <-
  data.frame(
    value = aggregate(ex3, value ~ group, FUN = "max")$value,
    group = unique(ex3$group),
    effsize = c(NA, effsize_sub$effsize)
  )


dodge <- position_dodge(width = 0.4)

max_dens3 <- aggregate(data4, data ~ treat, FUN = "max")
max_dens3$sd <- aggregate(data4, data ~ treat, FUN = "sd")$data

# Format the statistic with two digits after the decimal point
formatted_statistic <- sprintf("%.2f", res.aov3$statistic)

# Create subtitle with markdown syntax
statistic_label3 <-
  glue::glue(
    "<b>c) </b>Kruskal-Wallis, <i>H</i><sub>{res.aov3$df}</sub> = {formatted_statistic}, <i>p</i> = {res.aov3$p}"
  )

# L4-B9 ggplot  ######
P3 <- ggplot(ex3, aes(y = value, x = group)) +
  geom_point() +
  add_pvalue(
    df_p_val3[3:5, ],
    label = "p = {p.adj}",
    tip.length = 0.03,
    fontface = "bold",
    y.position = 0.4,
    label.size = 2.75,
    step.increase = 0,
    bracket.shorten = 0.1,
    vjust = -0.3
  ) +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    panel.background = element_rect(fill = 'white'),
    legend.position = "bottom",
    plot.margin = unit(c(-0.5, 0, 0, 0), 'cm'),
    plot.subtitle = element_markdown(),
    axis.title.x = element_markdown()
  ) +
  xlab(expression(paste(""))) +
  scale_y_continuous(expression("Growth rate (d" ^ -1 * ")")) +
  xlab("Photon flux densities (<i>&mu;</i>mol photons m<sup>-2</sup>s<sup>-1</sup>)") +
  coord_cartesian(ylim = c(0, 0.5)) +
  labs(subtitle = statistic_label3)

# Combine all plots ######
all_plots <-
  P1 + P2 + P3 + plot_layout(guides = "collect", ncol = 1, axis_titles = "collect_y") &
  theme(legend.position = "bottom")

ggsave(
  "Growth_all.png",
  all_plots ,
  path = "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/AP1/Light-AP1b/figures",
  dpi = 300,
  width = 3.5,
  height = 4.9,
  units = "in"
)

# In-between strain statistics ######
# combine all three strains
ex <- ex %>% mutate(strain = as.factor("L2-D2"))
ex2 <- ex2 %>% mutate(strain = as.factor("L4-B1"))
ex3 <- ex3 %>% mutate(strain = as.factor("L4-B9"))

ex_all <- rbind(ex, ex2[, c(1:2, 4)], ex3[, c(1:2, 4)])
ex_all$group <- factor(ex_all$group, levels = c(20, 60, 100, 200))

# Perform Kruskal-Wallis test for each light level in-between strains
gr_all <- ex_all %>%
  filter(group != 60) %>%
  group_by(group) %>%
  nest() %>%
  mutate(kruskal_results = map(data, ~ kruskal_test(value ~ strain, data = .x))) %>%
  unnest(c(kruskal_results)) %>%
  mutate(conover_results = map(data, ~ conover.test(
    .x$value, .x$strain, altp = T, method = "BH"
  ))) %>%
  unnest(c(conover_results))

# Correlation between growth rates and Photon flux densities of all strains  ######
# Aggregate data from each strain
agg_total <- rbind(
  aggregate(value ~ group, data = ex, FUN = "mean"),
  aggregate(value ~ group, data = ex2, FUN = "mean"),
  aggregate(value ~ group, data = ex3, FUN = "mean")
)

# Add light intensity values to aggregated data
agg_total$light <-
  c(c(20, 100, 200), rep(c(20, 60, 100, 200), each = 1, length.out = 8))

# Aggregate data by light intensity
light <- c(20, 60, 100, 200)
value <- aggregate(agg_total, value ~ light, FUN = "mean")$value

# Calculate standard deviation for aggregated data
agg_total2 <- aggregate(agg_total, value ~ light, FUN = "mean")
agg_total2$sd <-
  aggregate(agg_total, value ~ light, FUN = "sd")$value

# Define parameters for the model
parms_a <- c(a = 0.5, b = 0, c = 0.02)

# Define model function (regular asymptotic regression)
model_a <- function(parms_a, light)
  with(as.list(parms_a),
       return(a - (a - b) * exp(-c * light)))

# Define cost function for optimization
ModelCost_a <- function(P) {
  out <- model_a(P, light)
  return(value - out)  # residuals
}

# Fit the model to the data through a Nealder-Mead method
fit_a <-
  modFit(
    f = ModelCost_a,
    p = parms_a,
    method = "Nelder-Mead",
    upper = c(0.6, Inf, Inf)
  )
summary(fit_a)

# Generate model predictions
model_a <-
  data.frame(
    x = seq(0, 1000, 1),
    y = fit_a$par[1] - (fit_a$par[1] - fit_a$par[2]) * exp(-fit_a$par[3] * seq(0, 1000, 1))
  )

# Define limits for error bars
limits <-
  aes(ymax = (value + sd),
      ymin = (value - sd),
      width = 5) #Set up the error bars

# Define labels for annotations
label1 <- expression("x " [y / 2] * " = 48")
label2 <- expression("y"[max] * " = 0.37 d" ^ -1 * "")
label3 <- expression("Formula: y = a - (a-b) e" ^ (-cx) * "")
label4 <- expression("RMSE = 0.031")
label5 <- expression(I["k"] ~ "=" ~ 57 ~ italic("\u03BC") ~ "mol photons m"^-2 * "s"^-1)
label6 <- expression("x"[italic("y = 0")] ~ "=" ~ 9 ~ italic("\u03BC") ~ "mol photons m"^-2 * "s"^-1)

label5 <- expression(I["k"] ~ "=" ~ 57)
label6 <- expression("x"[italic("y = 0")] ~ "=" ~ 9)
# calculate minimum saturation irradiance Ik (tangente) at x = 48
deriv_model <-
  data.frame(
    x = seq(0, 1000, 1),
    y = fit_a$par[1] * fit_a$par[3] * exp(
      -fit_a$par[3] * seq(0, 1000, 1) - fit_a$par[2] * fit_a$par[3] * exp(-fit_a$par[3] *
                                                                            seq(0, 1000, 1))
    )
  )

m_tangente <-
  (
    fit_a$par[1] * fit_a$par[3] * exp(-fit_a$par[3] * 0) - fit_a$par[2] * fit_a$par[3] *
      exp(-fit_a$par[3] * 0)
  )
n_tangente <- model_a[which(model_a$x == 0), 2] - (m_tangente * 0)

deriv_model$tan <- m_tangente * deriv_model$x + n_tangente

# Plot the data and model  ######
ggplot(agg_total2, aes(x = light, y = value)) +
  geom_point() +
  geom_errorbar(limits) +
  geom_line(
    data = deriv_model[1:58,],
    aes(x = x, y = tan),
    col = "red",
    linetype = "dotted"
  ) +
  annotate(
    "text",
    label = as.character(label2),
    x = 250,
    y = 0.39,
    parse = T
  ) +
  annotate(
    "text",
    label = as.character(label5),
    x = 60,
    y = 0.39,
    parse = T
  ) +
  annotate(
    "text",
    label = as.character(label6),
    x = 50,
    y = 0,
    parse = T
  ) +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    panel.background = element_rect(fill = 'white'),
    legend.position = "bottom",
    axis.title.x = element_markdown()
  ) +
  ylab(expression("Growth rate (d" ^ -1 * ")")) +
  xlab("Photon flux densities (<i>&mu;</i>mol photons m<sup>-2</sup>s<sup>-1</sup>)") +
  geom_line(data = model_a, aes(x = x, y = y), col = "black") +
  coord_cartesian(xlim = c(0, 300), ylim = c(0, 0.4)) +
  geom_segment(
    aes(
      x = -5,
      xend = 500,
      y = fit_a$par[1],
      yend = fit_a$par[1]
    ),
    linetype = "dotted",
    col = "black"
  )

ggsave(
  "Growth_cor.png",
  path = "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/AP1/Light-AP1b/figures",
  dpi = 300,
  width = 3.5,
  height = 2.8,
  units = "in"
)

# export word document with package citations used in ingestion rate and GDA R-files ######
pacman::p_load(grateful)
cite_packages(
  out.format = "docx",
  out.dir = "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/AP1/",
  pkgs = "Session",
  out.file = "GR_packages"
)

# Garbage collection: call after large objects have been removed  ################################

gc()

# Delete workspace, clean environment/variables and restarte R if used memory piles up

dev.off()

rm(list = ls())

.rs.restartR()


# restart R

.rs.restartR()
