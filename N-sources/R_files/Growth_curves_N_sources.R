## Growth curves of both A. pseudogonyaulax strains L2-D2 and L4-B1 exposed to different nitrogen sources
## Published in Limnology & Oceanography: "Effects of bottom-up factors on growth and toxin content of a harmful algae bloom dinoflagellate"
## All raw-data available on PANGAEA: https://doi.pangaea.de/10.1594/PANGAEA.965195 
## Questions to: kristof-moeller@outlook.de
## Kristof Möller 09/22
## Alfred-Wegener-Institute Bremerhaven

# Use pacman to load add-on packages as desired
pacman::p_load(
  pacman,
  superb,
  stats,
  ggplot2,
  scales,
  dplyr,
  ggpubr,
  tibble,
  rstatix,
  NCmisc,
  ggthemes,
  patchwork,
  tidyr,
  extrafont,
  ggtext
)

# Load Windows Fonts and define Times font
loadfonts(device = "win")
windowsFonts(Times = windowsFont("Times"))

# Load and transform data ########

data <-
  read.csv(
    "C:\\Users\\krist\\OneDrive\\Dokumente\\AWI\\Promotion\\AP1\\N-AP1a\\R_txt_files\\L2D2_growth.txt",
    sep = "",
    header = FALSE
  )

# put first column as row names and transpose the data frame
data <-
  data %>% remove_rownames %>% column_to_rownames(var = "V1") %>% t() %>% as.data.frame()

# change variable type of all columns to numeric
data <- as.data.frame(sapply(data[, ], as.numeric))

# list used packages in the R-file
# used_packages <- list.functions.in.file("C:\\Users\\krist\\OneDrive\\Dokumente\\AWI\\Promotion\\AP1\\N-AP1a\\R_files\\Growth_all_counts.R")

data1 <- data %>%
  # Pivot the data from wide to long format
  pivot_longer(
    cols = -c(time_rest, time_C), 
    names_to = "treat",
    values_to = "data"
  ) %>%
  # Convert the 'treat' column into a factor based on pattern matching
  mutate(
    treat = factor(
      case_when(
        grepl("C", treat) ~ "N-deplete",
        grepl("NO", treat) ~ "NO3",
        grepl("NH", treat) ~ "NH4",
        grepl("U", treat) ~ "Urea"
      )
    ),
    # Determine 'time' based on 'treat' -> N-deplete treatment had different sampling days
    time = ifelse(grepl("N-deplete", treat), time_C, time_rest)
  ) %>%
  # Remove original time columns and arrange the data
  select(-c(time_C, time_rest)) %>%
  arrange(treat, time) %>% 
  drop_na(data) # Drop rows with NA values in the 'data' column

data1 <- data1 %>%
# Apply log10 transformation to cell counts and create a 'replicate' column
  mutate(
    data_log = log10(data),
    replicate = rep(1:3, length.out = length(data1$data))
  ) %>%
  # Arrange the data by 'time'
  arrange(time)

# Two-way repeated measures ANOVA #######
# Check assumptions for two-way repeated measures ANOVA first: Normality and Sphericity
# 
# for (each_treat in unique(data1$treat)) {
#   assumptions_check <- subset(data1, data1$treat == each_treat)
#   assumptions_check2 <-
#     pivot_wider(assumptions_check,
#                 names_from = replicate,
#                 values_from = data_log)
# 
#   # Extract non-NA values from each column separately
#   non_na_values <- lapply(assumptions_check2[, c(4:6)], function(x) x[complete.cases(x)])
# 
#   # Combine the non-NA values into a data frame
#   non_na_df <- as.data.frame(do.call(cbind, non_na_values))
# 
#   print(MauchlySphericityTest(non_na_df)) # Checks for Sphericity 
# }
# 
# data1 %>% group_by(treat,replicate) %>% shapiro_test(data) # Checks for normality

# # Exporting data for Pangaea upload
# write.table(data1, "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/AP1/N-AP1a/PANGAEA/cell_counts_L2D2.txt", sep = "\t", row.names = FALSE)

# repeated measures ANOVA is only valid for matching sampling times, hence N-deplete treatment & last Urea sampling day
# are removed for this analysis 
# time as within factor (within the same treatment) and treatment as between factor!

# Convert 'treat' column to factor type
data1 <- data1 %>% convert_as_factor(treat) 

# Perform repeated measures ANOVA
res.aov <- anova_test(
  data = data1 %>% subset(treat != "N-deplete"  & as.numeric(as.character(time)) < 17),
  dv = data_log, # use log-transformed cell-counts
  wid = replicate,
  within = time,
  between = treat
)

# Display ANOVA table
get_anova_table(res.aov)

# Perform one-way ANOVAs and pairwise comparisons at each time point after significant interactions in repeated measures ANOVA
one.way <- data1 %>% 
  subset(treat != "N-deplete"  & as.numeric(as.character(time)) < 17) %>%
  group_by(time) %>%
  anova_test(dv = data_log, wid = replicate, within = treat) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "BH")

# Display adjusted p-values for one-way ANOVA
one.way

# Perform pairwise comparisons
pwc <- data1 %>% 
  subset(treat != "N-deplete"  & as.numeric(as.character(time)) < 17) %>%
  group_by(time) %>%
  pairwise_t_test(data ~ treat, paired = FALSE, p.adjust.method = "BH")

# Display pairwise comparisons
pwc

# calculate mean and standard deviation of cell counts and log-transformed cell counts at each time point
agg <- aggregate(data1, data ~ treat + time, FUN = "mean")
agg$sd <- aggregate(data1, data ~ treat + time, FUN = "sd")$data

agg2 <- aggregate(data1, data_log ~ treat + time, FUN = "mean")
agg2$sd <-
  aggregate(data1, data_log ~ treat + time, FUN = "sd")$data

# Define parameters for ggplot
# Set position dodge width
dodge <- position_dodge(width = 0.5)

# Define linetypes for each treatment
lines <- c(
  'N-deplete' = "longdash",
  'NH4' = "solid",
  'NO3' = "dotted",
  "Urea" = "twodash"
)

# Define point shapes for each treatment
points <- c(
  'N-deplete' = 15,
  'NH4' = 17,
  'NO3' = 1,
  "Urea" = 4
)

# HTML-style treatment names for labels
names <- c(
  bquote("N-deplete"),
  bquote("NH"[4] ^ + ~ ""),
  bquote("NO"[3] ^ - ~ ""),
  bquote("urea")
)

# Treatment names according to dataframe --> Also needed for labels
names2 <- c("N-deplete", "NH4", "NO3", "Urea")

# Spike labels
spike <- ("+N")
spike2 <- deparse(bquote(paste("NO"[3] ^ - ~ "", "/", "NH"[4] ^ + ~ "", " spike")))

# ANOVA test results as subtitle label

# Format the statistic with two digits after the decimal point
formatted_statistic <- sprintf("%.2f", res.aov$F[3])

# Create subtitle with markdown syntax
subtitle_label <-
  glue::glue(
    "<b>a)</b> ANOVA, <i>F</i><sub>{res.aov$DFn[3]},{res.aov$DFd[3]}</sub> = {formatted_statistic}, <i>p</i> < 0.0001"
  )

# Colorblind-friendly colors
colors <- c(colorblind_pal()(5))[c(2:3, 5, 4)]

# Create ggplot object P1_c
P1_c <- ggplot(agg2, aes(x = time, y = data_log, col = treat)) +  # Set up the ggplot with data and aesthetics
  geom_pointrange(  # Add pointrange geom for plotting data points with error bars
    aes(ymin = data_log - sd, ymax = data_log + sd),  # Define y values for error bars
    show.legend = FALSE,  # Hide legend for this layer
    size = 0.25,  # Set size of points
    alpha = 1,  # Set alpha (transparency) of points
    linetype = "solid"  # Set line type for points
  ) +
  geom_line(  # Add line geom for connecting data points
    data = subset(agg2, agg2$treat != "Urea"),  # Subset data to exclude "Urea" treatment
    linewidth = 0.5,  # Set width of lines
    show.legend = FALSE  # Hide legend for this layer
  ) +
  geom_line(  # Add line geom for "Urea" treatment
    data = agg2 %>% filter(treat == "Urea" &  # Filter data for "Urea" treatment and specific time points
                             agg2$time < 15 |  
                             agg2$time == 17.81),
    linewidth = 0.5,  # Set width of lines
    show.legend = FALSE  # Hide legend for this layer
  ) +
  geom_line(  # Add line geom for "Urea" treatment (continued)
    data = agg2 %>% filter(  # Filter data for "Urea" treatment and specific time points
      treat == "Urea" &
        agg2$time > 15 &  # Time points greater than 15, not equal to 17.81
        agg2$time != 17.81 |
        treat == "Urea" &
        agg2$time == 16.71 | treat == "Urea" & agg2$time == 14.65
    ),
    linewidth = 0.5,  # Set width of lines
    show.legend = FALSE  # Hide legend for this layer
  ) +
  xlab("") +  # Set empty x-axis label
  theme_classic() +  # Apply classic theme
  ylab(bquote("Cell density (log"[10] * "(cells mL" ^ -1 * "))")) +  # Set y-axis label with mathematical expression
  theme(  # Customize plot theme
    panel.background = element_rect(fill = 'white'),  # Set panel background color
    strip.text = element_text(size = 14),  # Set text size for strip labels
    plot.margin = unit(c(0, 0, -0.5, 0), "cm"),  # Set plot margins
    plot.subtitle = element_markdown() 
  ) +
  scale_color_manual(  # Define manual color scale for treatments
    name = "N-source",  # Set legend title
    breaks = names2,  # Define breaks (legend labels)
    values = colors,  # Assign colors to treatments
    labels = names  # Set legend labels
  ) +
  coord_cartesian(xlim = c(0, 30), ylim = c(1.5, 4)) +  # Set Cartesian coordinate system limits
  scale_x_continuous(breaks = seq(0, 30, 5)) +  # Set x-axis breaks
  labs(subtitle = subtitle_label)  # Add subtitle with expression label
  
## add A. pseudogonyaulax strain L4-B1 as well

# Load and transform data ########
data2 <-
  read.csv(
    "C:\\Users\\krist\\OneDrive\\Dokumente\\AWI\\Promotion\\AP1\\N-AP1a\\\\R_txt_files\\L4B1_growth.txt",
    sep = "",
    header = FALSE
  )

# put first column as row names and transpose the data2 frame
data2 <-
  data2 %>% remove_rownames %>% column_to_rownames(var = "V1") %>% t() %>% as.data.frame()

# change variable type of all columns to numeric
data2 <- as.data.frame(sapply(data2[, ], as.numeric))

data2 <- data2 %>%
  # Pivot the data from wide to long format
  pivot_longer(
    cols = -c(time_rest, time_C), 
    names_to = "treat",
    values_to = "data"
  ) %>%
  # Convert the 'treat' column into a factor based on pattern matching
  mutate(
    treat = factor(
      case_when(
        grepl("C", treat) ~ "N-deplete",
        grepl("NO", treat) ~ "NO3",
        grepl("NH", treat) ~ "NH4",
        grepl("U", treat) ~ "Urea"
      )
    ),
    # Determine 'time' based on 'treat' -> N-deplete treatment had different sampling days
    time = ifelse(grepl("N-deplete", treat), time_C, time_rest)
  ) %>%
  # Remove original time columns and arrange the data
  select(-c(time_C, time_rest)) %>%
  drop_na(data) # Drop rows with NA values in the 'data' column

data2 <- data2 %>% na.omit() %>%
  mutate(data = log10(data),
         replicate = c(rep(1:3, length.out = 108),
                        1, 2,
                        rep(1:3, length.out = 6))) %>% arrange(time)

# # Exporting data for Pangaea
# write.table(data2, "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/AP1/N-AP1a/PANGAEA/cell_counts_L4B1.txt", sep = "\t", row.names = FALSE)

## Two-way repeated measures ANOVA
## Check assumptions for two-way repeated measures ANOVA first: Normality and Sphericity


# Two-way repeated measures ANOVA #######
# Check assumptions for two-way repeated measures ANOVA first: Normality and Sphericity
for (each_treat in unique(data2$treat)) {
  assumptions_check <- subset(data2, data2$treat == each_treat)
  assumptions_check2 <-
    pivot_wider(assumptions_check,
                names_from = replicate,
                values_from = data)

  # Extract non-NA values from each column separately
  non_na_values <- lapply(assumptions_check2[, c(3:5)], function(x) x[complete.cases(x)])

  # Combine the non-NA values into a data frame
  non_na_df <- as.data.frame(do.call(cbind, non_na_values))

  print(MauchlySphericityTest(non_na_df)) # Checks for Sphericity
}

data2 %>% group_by(treat, replicate) %>% shapiro_test(data) # Checks for normality

# repeated measures ANOVA is only valid for matching sampling times, hence N-deplete treatment & last Urea sampling day
# are removed for this analysis 
# time as within factor (within the same treatment) and treatment as between factor!

# Convert 'treat' column to factor type
data2 <- data2 %>% convert_as_factor(treat) 

# Perform repeated measures ANOVA
res.aov2 <- anova_test(
  data = data2 %>% subset(
    treat != "N-deplete"  &
      as.numeric(as.character(time)) < 20 &
      as.numeric(as.character(time)) != 17.17
  ),
  dv = data,
  wid = replicate,
  within = time,
  between = treat,
  white.adjust = T
)

# Display ANOVA table
get_anova_table(res.aov2, correction = "auto")

# Perform one-way ANOVAs and pairwise comparisons at each time point after significant interactions in repeated measures ANOVA
one.way <- data2 %>% subset(
  treat != "N-deplete"  &
    as.numeric(as.character(time)) < 20 &
    as.numeric(as.character(time)) != 17.17
) %>%
  group_by(time) %>%
  anova_test(dv = data, wid = replicate, within = treat) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "BH")

# Display adjusted p-values for one-way ANOVA
one.way

# Perform pairwise comparisons
pwc <- data2 %>% subset(
  treat != "N-deplete"  &
    as.numeric(as.character(time)) < 20 &
    as.numeric(as.character(time)) != 17.17
) %>%
  group_by(time) %>%
  pairwise_t_test(data ~ treat, paired = F,
                  p.adjust.method = "BH")

# Display pairwise comparisons
pwc

# calculate mean and standard deviation of cell counts and log-transformed cell counts at each time point
agg3 <- aggregate(data2, data ~ treat + time, FUN = "mean")
agg3$sd <-
  aggregate(data = data2, data ~ treat + time, FUN = "sd")$data

agg3$time <- as.numeric(as.character(agg3$time))


# Define parameters for ggplot

# Set position dodge width
dodge <- position_dodge(width = 0.5)

# Define linetypes for each treatment
lines <- c(
  'N-deplete' = "longdash",
  'NH4' = "solid",
  'NO3' = "dotted",
  "Urea" = "twodash"
)

# Define point shapes for each treatment
points <- c(
  'N-deplete' = 15,
  'NH4' = 17,
  'NO3' = 1,
  "Urea" = 4
)

# HTML-style treatment names for labels
names <- c(
  bquote("N-deplete"),
  bquote("NH"[4] ^ + ~ ""),
  bquote("NO"[3] ^ - ~ ""),
  bquote("urea")
)

# Treatment names according to dataframe --> Also needed for labels
names2 <- c("N-deplete", "NH4", "NO3", "Urea")

# Spike labels
spike2 = deparse(bquote(paste("NO"[3] ^ - ~ "", "/", "NH"[4] ^ + ~ "", " spike")))

# ANOVA test results as subtitle label
subtitle_label2 <-
  get_test_label(res.aov2, detailed = TRUE, type = "expression", correction = "none")

# Format the statistic with two digits after the decimal point
formatted_statistic <- sprintf("%.2f", res.aov2$F[3])

# Create subtitle with markdown syntax
subtitle_label2 <-
  glue::glue(
    "<b>b)</b> ANOVA, <i>F</i><sub>{res.aov2$DFn[3]},{res.aov2$DFd[3]}</sub> = {formatted_statistic}, <i>p</i> < 0.0001"
  )

# Create ggplot object P2_c
P2_c <- ggplot(agg3, aes(x = time, y = data, col = treat)) +  # Set up the ggplot with data and aesthetics
  geom_pointrange(  # Add pointrange geom for plotting data points with error bars
    aes(ymin = data - sd, ymax = data + sd),  # Define y values for error bars
    size = 0.25,  # Set size of points
    alpha = 1,  # Set alpha (transparency) of points
    show.legend = TRUE,  # Show legend for this layer
    linetype = "solid"  # Set line type for points
  ) +
  geom_line(  # Add line geom for connecting data points
    data = subset(agg3, treat != "Urea"),  # Subset data to exclude "Urea" treatment
    linewidth = 0.5,  # Set width of lines
    show.legend = TRUE  # Show legend for this layer
  ) +
  geom_line(  # Add line geom for "Urea" treatment
    data = agg3 %>% filter(treat == "Urea" &  # Filter data for "Urea" treatment and specific time points
                             time < 15 | time == 17.17),
    linewidth = 0.5,  # Set width of lines
    show.legend = TRUE  # Show legend for this layer
  ) +
  geom_line(  # Add line geom for "Urea" treatment (continued)
    data = agg3 %>% filter(  # Filter data for "Urea" treatment and specific time points
      treat == "Urea" & time > 15 & time != 17.17 |  # Time points greater than 15, not equal to 17.17
        treat == "Urea" &
        time == 16.06 | treat == "Urea" & time == 14.08  # Specific time points
    ),
    linewidth = 0.5,  # Set width of lines
    show.legend = TRUE  # Show legend for this layer
  ) +
  xlab("Time (days)") +  # Set x-axis label
  theme_classic() +  # Apply classic theme
  ylab(bquote("Cell density (log"[10] * "(cells mL" ^ -1 * "))")) +  # Set y-axis label with mathematical expression
  theme(  # Customize plot theme
    legend.position = "right",  # Set legend position
    legend.title = element_blank(),  # Remove legend title
    panel.background = element_rect(fill = "white"),  # Set panel background color
    legend.spacing = unit(0, "cm"),  # Set legend spacing
    legend.key.height = unit(1, "lines"),  # Set legend key height
    strip.text = element_text(size = 14),  # Set text size for strip labels
    plot.margin = unit(c(-0.5, 0, 0, 0), "cm"),  # Set plot margins
    plot.subtitle = element_markdown()
  ) +
  coord_cartesian(xlim = c(0, 30), ylim = c(1.5, 4)) +  # Set Cartesian coordinate system limits
  scale_color_manual(  # Define manual color scale for treatments
    name = "N-source",  # Set legend title
    breaks = names2,  # Define breaks (legend labels)
    values = colors,  # Assign colors to treatments
    labels = names  # Set legend labels
  ) +
  scale_x_continuous(breaks = seq(0, 30, 5)) +  # Set x-axis breaks
  labs(subtitle = subtitle_label2)   # Add subtitle with expression label

# combine cell count plots of both strains in one plot
all_plots <-
  P1_c + P2_c + plot_layout(guides = "collect", ncol = 1, axis_titles = "collect_y") &
  theme(legend.position = "bottom",
        legend.box.margin = margin(-10, 0, 0, 0))

# save the combined plot as a png
ggsave(
  "Growth_counts_all_color.png",
  all_plots,
  path = "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/AP1/N-AP1a/figures",
  dpi = 300,
  width = 3.5,
  height = 6,
  units = "in"
)

# Garbage collection: call after large objects have been removed
gc()

# Delete workspace, clean environment/variables and restarte R if used memory piles up
dev.off()

rm(list = ls())

# restart R
.rs.restartR()
