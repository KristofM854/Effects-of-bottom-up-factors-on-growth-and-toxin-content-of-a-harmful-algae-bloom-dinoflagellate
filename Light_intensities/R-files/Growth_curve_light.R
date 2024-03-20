##########################################
## Growth curves of A. pseudogonyaulax strains L2-D2, L4-B1 and L4-B9 exposed to different light intensities
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
  psych,
  rio,
  plyr,
  ggplot2,
  scales,
  dplyr,
  gridExtra,
  tibble,
  rstatix,
  PMCMRplus,
  tidyr,
  purrr,
  ggpp,
  superb,
  ggpubr,
  grid,
  janitor
)

# Load and transform data of L2-D2 ################################################

data <-
  read.csv(
    "C:\\Users\\krist\\OneDrive\\Dokumente\\AWI\\Promotion\\AP1\\Light-AP1b\\raw\\data_L2D2.txt",
    sep = "",
    header = FALSE,
    skip = 1,
    nrow = 10
  )

data <- data %>%
  remove_rownames %>%
  column_to_rownames(var = "V1") %>%
  t() %>%
  as.data.frame() %>%
  mutate_all( ~ as.numeric(gsub(",", ".", .))) %>%
  pivot_longer(cols = c(2:10),
               names_to = "treat",
               values_to = "data") %>%
  mutate(
    treat = as.factor(rep(
      c("20", "100", "200"), length.out = 108, each = 3
    )),
    data1 = data,
    data = log10(data),
    replicate = as.factor(rep(1:3, length.out = 108))
  ) %>%
  na.omit()


data <- dplyr::rename(data, time = "Day_adjusted")

# aggregate data for plot
agg <- data %>%
  aggregate(data ~ treat + time, FUN = "mean")
agg$sd <-
  aggregate(data = data, data ~ treat + time, FUN = "sd")$data

# Two-way repeated measures ANOVA #######
# Check assumptions for two-way repeated measures ANOVA first: Sphericity and Normality

# for(each_treat in unique(data$treat)) {
#   assumptions_check <- data %>% filter(treat == each_treat)
#   assumptions_check2 <-
#     pivot_wider(assumptions_check,
#                 names_from = replicate,
#                 values_from = data)
#   # Extract non-NA values from each column separately
#   non_na_values <-
#     lapply(assumptions_check2[, c(4:6)], function(x)
#       x[complete.cases(x)])
#   
#   # Combine the non-NA values into a data frame
#   non_na_df <- as.data.frame(do.call(cbind, non_na_values))
#   
#   print(MauchlySphericityTest(non_na_df)) # Checks for Sphericity
# }
# 
# data %>%
#   group_by(treat, replicate) %>% shapiro_test(data)

# filter sampling days < 20. repeated measures ANOVA requires equal sampling days and low light treatment was sampled longer
data1 <- data %>% filter(time < 20)

# Perform repeated measures ANOVA 
res.aov <-
  anova_test(
    data = data1,
    dv = data, # use log-transformed data
    wid = replicate,
    within = time, 
    between = treat
  )

# Display ANOVA table
get_anova_table(res.aov)

# Perform one-way ANOVAs and pairwise comparisons at each time point after significant interactions in repeated measures ANOVA
one.way <- data1 %>%
  group_by(time) %>%
  anova_test(dv = data, wid = replicate, within = treat) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "BH")

# Display adjusted p-values for one-way ANOVA
one.way

# Perform pairwise comparisons
pwc <- data1 %>%
  group_by(time) %>%
  pairwise_t_test(
    data ~ treat, paired = F,
    p.adjust.method = "BH"
  )

# Display pairwise comparisons
pwc

# Define parameters for ggplot
# linetypes
lines <- c('20' = "solid",
           '100' = "longdash",
           '200' = "dotted")

# point types
points <- c('20' = 15,
            '100' = 17,
            '200' = 19)

# breaks for x-axis and legend
breaks <- c("20", "100", "200")
breaks2 <- c("20 PFDs", "100 PFDs", "200 PFDs")

# Set position dodge width and jitter
dodge <- position_jitternudge(width = 0.2, x = 0)

# colorblind colors but skip black (first color) and switch blueish and greenish because NH4 / NO3 are close to each other
colors <- c(colorblind_pal()(4))[c(4, 2, 3)]

# prepare dataframe to bind aggregated results of each strain
agg_all <- agg %>% mutate(strain = "strain A")

# ANOVA test results as subtitle label
subtitle_label <- get_test_label(res.aov, detailed = T)

P1_c <- ggplot(agg, aes(x = time, y = data, col = treat)) +
  
  geom_pointrange(
    aes(ymin = data - sd, ymax = data + sd),
    position = dodge,
    size = 0.25,
    alpha = 1,
    show.legend = F,
    linetype = "solid"
  ) +
  
  geom_line(show.legend = F) +
  
  xlab("") +
  ylab(bquote("log"[10] * "(cells mL" ^ -1 * ")")) +
  theme_classic() +
  theme(
    panel.background = element_rect(fill = 'white'),
    strip.background = element_blank(),
    strip.text = element_text(size = 14),
    plot.margin = unit(c(0, 0, -0.5, 0), "cm")
  ) +
  scale_color_manual(
    name = "Photon flux densities",
    values = colors,
    breaks = breaks,
    labels = breaks2
  ) +
  coord_cartesian(ylim = c(1.5, 4.0), xlim = c(0, 25)) +
  labs(subtitle = bquote(paste(bold("A) "), .(subtitle_label))))


## add A. pseudogonyaulax strain L4-B1 as well

# Load and transform data of L4-B1 ################################################
data <-
  read.csv(
    "C:\\Users\\krist\\OneDrive\\Dokumente\\AWI\\Promotion\\AP1\\Light-AP1b\\raw\\data_L4B1.txt",
    sep = "",
    skip = 1,
    header = FALSE,
    nrow = 11
  )

data <- data %>%
  remove_rownames %>%
  column_to_rownames(var = "V1") %>%
  t() %>%
  as.data.frame() %>%
  mutate_all(~ as.numeric(gsub(",", ".", .))) %>%
  pivot_longer(cols = c(3:11),
               names_to = "treat",
               values_to = "data") %>%
  mutate(
    treat = factor(rep(
      c("20", "100", "200"), length.out = 108, each = 3
    ), levels = c(20, 60, 100, 200)),
    data1 = data,
    data = log10(data), # log-transform cell counts data
    replicate = rep(1:3, length.out = 108)
  ) %>%
  arrange(treat) %>%
  mutate(time = c(Day_20[1:36], Day_100_200[37:108])) %>%
  arrange(time)

# Export for PANGAEA data repository
# data %>%
#   dplyr::select(treat, time, data1, data, replicate) %>%
#   write.table(
#     "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/AP1/Light-AP1b/PANGAEA/cell_counts_L4B1.txt",
#     sep = "\t",
#     row.names = FALSE
#   )

# aggregate data for plot
agg <- data %>%
  aggregate(data ~ treat + time, FUN = "mean")
agg$sd <-
  aggregate(data = data, data ~ treat + time, FUN = "sd")$data

# substitute light intensity of low light (20) treatment to 60 after day 9.5
# light intensity was increased as algae population was stagnating
agg$treat2 <- agg$treat
for (i in 1:33) {
  if (agg$treat[i] == "20" & agg$time[i] > 9.5) {
    agg$treat2[i] = "60"
  }
}

# Two-way repeated measures ANOVA #######
# Check assumptions for two-way repeated measures ANOVA first: Sphericity and Normality
# for(each_treat in unique(data$treat)) {
#   assumptions_check <- data1 %>% filter(treat == each_treat)
#   assumptions_check2 <-
#     pivot_wider(assumptions_check,
#                 names_from = replicate,
#                 values_from = data1)
#   # Extract non-NA values from each column separately
#   non_na_values <-
#     lapply(assumptions_check2[, c(4:6)], function(x)
#       x[complete.cases(x)])
# 
#   # Combine the non-NA values into a data frame
#   non_na_df <- as.data.frame(do.call(cbind, non_na_values))
# 
#   print(MauchlySphericityTest(non_na_df)) # Checks for Sphericity
# }
# 
# data1 %>%
#   group_by(treat, replicate) %>% shapiro_test(data1)

data1 <- data %>% convert_as_factor(time, treat)

# remove sampling days which do not contain all treatments for repeated measures ANOVA
data2 <-
  data1 %>% subset(
    as.numeric(as.character(time)) != 9.88 &
      as.numeric(as.character(time)) != 10.83 &
      as.numeric(as.character(time)) != 11.88 &
      as.numeric(as.character(time)) != 12.83 &
      as.numeric(as.character(time)) != 13.83 &
      as.numeric(as.character(time)) != 20.92
  )

# perform repeated measures ANOVA
res.aov2 <- anova_test(
  data = data2,
  dv = data,
  wid = replicate,
  within = time,
  between = treat
)

# Display ANOVA table
get_anova_table(res.aov)

# Perform one-way ANOVAs and pairwise comparisons at each time point after significant interactions in repeated measures ANOVA
one.way <- data2 %>%
  group_by(time) %>%
  anova_test(dv = data, wid = replicate, within = treat) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "BH")

# Display adjusted p-values for one-way ANOVA
one.way

# Perform pairwise comparisons
pwc <- data2 %>%
  group_by(time) %>%
  pairwise_t_test(
    data ~ treat, paired = F,
    p.adjust.method = "BH"
  )

# Display pairwise comparisons
pwc

# bind aggregated results of each strain
agg_all <-
  rbind(agg_all,
        agg %>% mutate(strain = "strain B") %>% dplyr::select(-treat2))

# ANOVA test results as subtitle label
subtitle_label2 <- get_test_label(res.aov2, detailed = T)

P2_c <- ggplot(agg, aes(x = time, y = data, col = treat)) +
  
  geom_pointrange(
    aes(ymin = data - sd, ymax = data + sd),
    position = dodge,
    size = 0.25,
    alpha = 1,
    show.legend = F,
    linetype = "solid"
  ) +
  
  geom_line(show.legend = F) +
  
  xlab("") +
  ylab(bquote("log"[10] * "(cells mL" ^ -1 * ")")) +
  theme_classic() +
  theme(
    panel.background = element_rect(fill = 'white'),
    strip.background = element_blank(),
    strip.text = element_text(size = 14),
    plot.margin = unit(c(-0.5, 0, -0.5, 0), "cm")
  ) +
  scale_color_manual(
    name = "Photon flux densities",
    values = colors,
    breaks = breaks,
    labels = breaks2
  ) +
  geom_segment(
    aes(
      x = 9.9,
      xend = 9.9,
      y = 1.45,
      yend = 1.65
    ),
    linewidth = 0.5,
    arrow = arrow(length = unit(0.15, "cm")),
    colour = "black",
    show.legend = F
  ) +
  coord_cartesian(ylim = c(1.5, 4.0), xlim = c(0, 25)) +
  labs(subtitle = bquote(paste(bold("B) "), .(subtitle_label2))))


## add A. pseudogonyaulax strain L4-B9 as well

# Load and transform data of L4-B9 ################################################
data <-
  read.csv(
    "C:\\Users\\krist\\OneDrive\\Dokumente\\AWI\\Promotion\\AP1\\Light-AP1b\\raw\\data_L4B9.txt",
    sep = "",
    skip = 1,
    header = FALSE,
    nrow = 11
  )

data1 <- data %>%
  remove_rownames %>%
  column_to_rownames(var = "V1") %>%
  t() %>%
  as.data.frame() %>%
  mutate_all(~ as.numeric(gsub(",", ".", .))) %>%
  pivot_longer(cols = c(3:11),
               names_to = "treat",
               values_to = "data") %>%
  mutate(
    treat = factor(rep(
      c("20", "100", "200"), length.out = 108, each = 3
    ), levels = c(20, 60, 100, 200)),
    data1 = data,
    data = log10(data),
    replicate = rep(1:3, length.out = 108)
  ) %>%
  arrange(treat) %>%
  mutate(time = c(Day_20[1:36], Day_100_200[37:108])) %>%
  arrange(time)

# Export for PANGAEA data repository
data1 %>%
  dplyr::select(treat, time, data1, replicate) %>%
  write.table(
    "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/AP1/Light-AP1b/PANGAEA/cell_counts_L4B9.txt",
    sep = "\t",
    row.names = FALSE
  )

# aggregate data for plot
agg <- aggregate(data = data1, data ~ treat + time, FUN = "mean")
agg$sd <-
  aggregate(data = data1, data ~ treat + time, FUN = "sd")$data
agg$treat <- factor(agg$treat, levels = c("20", "100", "200"))
agg$treat2 <-
  factor(agg$treat, levels = c("20", "60", "100", "200"))

# Change light setting of all 20s after 9.5 days to 60
for (i in 1:34) {
  if (agg$treat[i] == "20" & agg$time[i] > 14) {
    agg$treat2[i] = "60"
  }
}
# Two-way repeated measures ANOVA #######
# Check assumptions for two-way repeated measures ANOVA first: Sphericity and Normality

# for(each_treat in unique(data$treat)) {
#   assumptions_check <- data1 %>% filter(treat == each_treat)
#   assumptions_check2 <-
#     pivot_wider(assumptions_check,
#                 names_from = replicate,
#                 values_from = data)
#   # Extract non-NA values from each column separately
#   non_na_values <-
#     lapply(assumptions_check2[, c(4:6)], function(x)
#       x[complete.cases(x)])
#   
#   # Combine the non-NA values into a data frame
#   non_na_df <- as.data.frame(do.call(cbind, non_na_values))
#   
#   print(MauchlySphericityTest(non_na_df)) # Checks for Sphericity
# }
# 
# data1 %>%
#   group_by(treat, replicate) %>% shapiro_test(data)


data1 <- data1 %>% convert_as_factor(time, treat)

# remove sampling days which do not contain all treatments for repeated measures ANOVA
data2 <-
  data1 %>% subset(
    as.numeric(as.character(time)) != 12.96 &
      as.numeric(as.character(time)) != 15.83 &
      as.numeric(as.character(time)) != 13.88 &
      as.numeric(as.character(time)) != 14.83
  )

# Perform repeated measures ANOVA 
res.aov3 <- anova_test(
  data = data2,
  dv = data,
  wid = replicate,
  within = time,
  between = treat
)

# bind aggregated results of each strain
agg_all <-
  rbind(agg_all,
        agg %>% mutate(strain = "strain C") %>% dplyr::select(-treat2))

# ANOVA test results as subtitle label
subtitle_label3 <- get_test_label(res.aov3, detailed = TRUE)

P3_c <- ggplot(agg, aes(x = time, y = data, col = treat)) +
  
  geom_pointrange(
    aes(ymin = data - sd, ymax = data + sd),
    position = dodge,
    size = 0.25,
    alpha = 1,
    show.legend = T,
    linetype = "solid"
  ) +
  
  geom_line(show.legend = T) +
  
  xlab("days") +
  ylab(bquote("log"[10] * "(cells mL" ^ -1 * ")")) +
  theme_classic() +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    panel.background = element_rect(fill = "white"),
    legend.spacing = unit(0, "cm"),
    legend.key.height = unit(1, "lines"),
    strip.text = element_text(size = 14),
    plot.margin = unit(c(-0.5, 0, 0, 0), "cm"),
    legend.box.margin = margin(-10, -10, -10, -10)
  ) +
  geom_segment(
    aes(
      x = 13.9,
      xend = 13.9,
      y = 1.7,
      yend = 1.9
    ),
    linewidth = 0.5,
    arrow = arrow(length = unit(0.15, "cm")),
    colour = "black"
  ) +
  scale_color_manual(
    name = "Photon flux densities",
    values = colors,
    breaks = breaks,
    labels = breaks2
  ) +
  coord_cartesian(ylim = c(1.5, 4.0), xlim = c(0, 25)) +
  labs(subtitle = bquote(paste(bold("C) "), .(subtitle_label3)))) +
  guides(col = guide_legend(nrow = 1))

# combine plots of all three strains of A. pseudogonyaulax and export ################
all_plots <-
  P1_c + P2_c + P3_c + plot_layout(guides = "collect", ncol = 1) &
  theme(legend.position = "bottom",
        legend.box.margin = margin(-10, 0, 0, 0))

ggsave(
  "Growth_counts_all_color.png",
  all_plots,
  path = "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/AP1/Light-AP1b/figures",
  dpi = 300,
  width = 11,
  height = 17.5,
  units = "cm"
)

# Garbage collection: call after large objects have been removed
gc()

# Delete workspace, clean environment/variables and restarte R if used memory piles up
dev.off()

rm(list = ls())

# restart R
.rs.restartR()
