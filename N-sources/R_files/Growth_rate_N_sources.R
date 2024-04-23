##########################################
## Growth rate calculations and plots of both A. pseudogonyaulax strains L2-D2 and L4-B1 exposed to different nitrogen sources
## Published in Limnology & Oceanography: "Effects of bottom-up factors on growth and toxin content of a harmful algae bloom dinoflagellate"
## All raw-data available on PANGAEA: https://doi.pangaea.de/10.1594/PANGAEA.965195 
## Questions to: kristof-moeller@outlook.de
## Kristof Möller 09/22
## Alfred-Wegener-Institute Bremerhaven
##########################################

# INSTALL AND LOAD PACKAGES ################################

# Installs pacman ("package manager") if needed
if (!require("pacman"))
  install.packages("pacman")

# Use pacman to load add-on packages as desired
pacman::p_load(
  pacman,
  ggplot2,
  DescTools,
  tibble,
  rstatix,
  stats,
  gridExtra,
  tidyr,
  dplyr,
  conover.test,
  growthrates,
  outliers,
  ggprism,
  NCmisc,
  patchwork,
  ggtext,
  glue,
  ragg,
  purrr
)

# Load and transform data ########
data1 <- # cell counts data
  read.csv(
    "C:\\Users\\krist\\OneDrive\\Dokumente\\AWI\\Promotion\\AP1\\N-AP1a\\R_txt_files\\L2D2_growth.txt",
    sep = "",
    header = FALSE
  )

data2_2D2 <- # max cell density
  as.data.frame(
    read.csv(
      "C:\\Users\\krist\\OneDrive\\Dokumente\\AWI\\Promotion\\AP1\\N-AP1a\\R_txt_files\\L2D2_max_density.txt",
      sep = "",
      header = FALSE
    )
  )

# introduce treatment column = different nitrogen sources with N-deplete being K-medium without nitrogen addition as control
data2_2D2$treat <- rep(c("N_deplete", "NO3", "NH4", "Urea"), each = 3)

# List used packages in the file
# used_packages <- list.functions.in.file("C:\\Users\\krist\\OneDrive\\Dokumente\\AWI\\Promotion\\AP1\\N-AP1a\\R_files\\Growth_all_rate.R")

# put first column as row names and transpose the data frame
data1 <-
  data1 %>% remove_rownames %>% column_to_rownames(var = "V1") %>% t() %>% as.data.frame()

# change NaN to NA; make all columns numeric
data1[, ] <- sapply(data1[, ], function(x)
  gsub(NaN, NA, x))
data1[, ] <- sapply(data1[, ], function(x)
  as.numeric(x))

## Growth rate calculations ######
# make linear fit for each bottle between 5 data points
# h: amount of data points used for calculation of growth rate!!
fit_C1 <- fit_easylinear(na.omit(data1$time_C), na.omit(data1$C1), h = 5)
fit_C2 <- fit_easylinear(na.omit(data1$time_C), na.omit(data1$C2), h = 5)
fit_C3 <- fit_easylinear(na.omit(data1$time_C), na.omit(data1$C3), h = 5)
fit_NO1 <- fit_easylinear(na.omit(data1$time_rest)[1:10], na.omit(data1$NO1), h = 5)
fit_NO2 <- fit_easylinear(na.omit(data1$time_rest)[1:10], na.omit(data1$NO2), h = 5)
fit_NO3 <- fit_easylinear(na.omit(data1$time_rest)[1:10], na.omit(data1$NO3), h = 5)
fit_NH1 <- fit_easylinear(na.omit(data1$time_rest)[1:10], na.omit(data1$NH1), h = 5)
fit_NH2 <- fit_easylinear(na.omit(data1$time_rest)[1:10], na.omit(data1$NH2), h = 5)
fit_NH3 <- fit_easylinear(na.omit(data1$time_rest)[1:10], na.omit(data1$NH3), h = 5)
fit_U1 <- fit_easylinear(na.omit(data1[1:8, 1]), na.omit(data1[1:8, 12]), h = 5)
fit_U2 <- fit_easylinear(na.omit(data1[1:8, 1]), na.omit(data1[1:8, 13]), h = 5)
fit_U3 <- fit_easylinear(na.omit(data1[1:8, 1]), na.omit(data1[1:8, 14]), h = 5)

# extract coefficients of the linear fits 
matrix_coef <-
  as.data.frame(rbind(
    coef(fit_C1),
    coef(fit_C2),
    coef(fit_C3),
    coef(fit_NO1),
    coef(fit_NO2),
    coef(fit_NO3),
    coef(fit_NH1),
    coef(fit_NH2),
    coef(fit_NH3),
    coef(fit_U1),
    coef(fit_U2),
    coef(fit_U3)
  ))

# Assign column names based on the original data
matrix_coef$treat <- colnames(data1[, 3:14])

# construct new dataframe 'ex' with maximal growth rate and treatments for following analysis
ex = data.frame(treat = rep(as.factor(c(
  "N_deplete", "NO3", "NH4", "Urea"
)), each = 3), data = matrix_coef$mumax)

## Check for potential outliers and remove them #####
# 
# outliers <- c()
# 
# for (each_group in unique(ex$treat)) {
#   dix1 <-
#     if (dixon.test(subset(ex, ex$treat == each_group)$data)$p.val < 0.05) {
#       parse_number(dixon.test(subset(ex, ex$treat == each_group)$data)$alternative)
#     }
#   if (!is.null(dix1)) {
#     outliers <- c(outliers, dix1)
#   }
#   
# }
# 
# ex <- ex[!(ex$data %in% outliers),]
## No outliers found!!

# Statistical Analysis 
# Perform a non-parametric Kruskal test since three samples do not satisfy the normality assumption
# Null hypothesis (H0): The growth rate of A. pseudogonyaulax exposed to different nitrogen sources is equal
res.aov1 <- kruskal_test(data ~ treat, data = ex)

# Conduct a Conover-Iman non-parametric PostHoc test if p < 0.05, with a liberal p-value adjustment using Benjamini & Hochberg (FDR)
Con_a <- conover.test(ex$data, ex$treat, method = "BH", altp = TRUE, alpha = 0.1)

# Calculate Cohen's d effect size: ~0.2 (small effect), ~0.5 (medium effect), ~0.8 (large effect)
effsize <- cohens_d(ex, data ~ treat, hedges.correction = TRUE)
effsize_sub <- abs(round(effsize[c(1:3), 4], digits = 1))

# Check for significant differences in cell densities
# Null hypothesis (H0): The maximum cell densities of A. pseudogonyaulax exposed to different nitrogen sources are equal
res.aov2 <- kruskal_test(V1 ~ treat, data = data2_2D2)

# Perform a Conover-Iman non-parametric PostHoc test if p < 0.05, with a liberal p-value adjustment using Benjamini & Hochberg (FDR)
Con_a2 <- as.data.frame(rbind(ConoverTest(V1 ~ treat, data = data2_2D2, method = "BH")[[1]]))

# Extract p-values and group comparisons from the Conover-Iman test output
df_p_val <- data.frame(
  p.adj = sprintf("%.3f", Con_a$altP.adjusted),
  group1 = sub("-.*", "", Con_a$comparisons),
  group2 = sub(".*-", "", Con_a$comparisons)
)

df_p_val <- df_p_val %>% subset(p.adj < 0.1)
df_p_val[] <- sapply(df_p_val, function(x) gsub(pattern = " ", replacement = "", x))

# Get maximum densities for plotting
max_dens_2D2_2 <- aggregate(data = data2_2D2, V1 ~ treat, FUN = "max")

# Prepare ggplot parameters
dodge <- position_dodge(width = 0.4)

# HTML-style treatment names for labels
labels <-
  c(
    glue("N-deplete"),
    glue("NH<sub>4</sub><sup>+</sup>"),
    glue("NO<sub>3</sub><sup>-</sup>"),
    glue("urea")
  )

names2 <- c("N_deplete", "NH4", "NO3", "Urea")

# Format the statistic with two digits after the decimal point
formatted_statistic <- sprintf("%.2f", res.aov1$statistic)

# Create subtitle with markdown syntax
statistic_label <-
  glue::glue(
    "<b>a)</b> Kruskal-Wallis, <i>H</i><sub>{res.aov1$df}</sub> = {formatted_statistic}, <i>p</i> = {res.aov1$p}"
  )

df_p_val$p <- c("0.01", "0.01", "0.05", "0.05", "0.05")

# Create ggplot with treatment (x) and data (y) aesthetics
P1 <- ggplot(ex, aes(x = treat, y = data)) +
  # Add points to the plot with size 2.0
  geom_point(size = 2.0) +
  # Add p-values annotations using add_pvalue function
  add_pvalue(
    df_p_val[1, ],  # p-value info for the first comparison
    label = "p < {p}",  # label format for the p-value
    tip.length = 0.03,  # length of the tooltip
    fontface = "bold",  # font style
    y.position = 0.485,  # y-position for the p-value
    label.size = 2.75,  # size of the label
    step.increase = 0,  # step increase for label positioning
    bracket.shorten = 0,  # shorten bracket length
    vjust = -0.3,  # vertical justification
    remove.bracket = TRUE  # remove bracket around the p-value
  ) + add_pvalue(
    df_p_val[c(2), ],  # p-value info for the second comparison
    label = "p < {p}",  # label format for the p-value
    tip.length = 0.01,  # length of the tooltip
    fontface = "bold",  # font style
    y.position = 0.45,  # y-position for the p-value
    label.size = 0,  # size of the label
    step.increase = 0.175,  # step increase for label positioning
    bracket.shorten = 0,  # shorten bracket length
    vjust = -0.3  # vertical justification
  ) +
  # Additional p-value annotations for other comparisons
  add_pvalue(
    df_p_val[c(4), ],  # p-value info for the third comparison
    label = "p < {p}",  # label format for the p-value
    tip.length = 0.03,  # length of the tooltip
    fontface = "bold",  # font style
    y.position = 0.55,  # y-position for the p-value
    label.size = 0,  # size of the label
    step.increase = 0.175,  # step increase for label positioning
    bracket.shorten = 0,  # shorten bracket length
    vjust = -0.3  # vertical justification
  ) +
  add_pvalue(
    df_p_val[c(5), ],  # p-value info for the fourth comparison
    label = "p < {p}",  # label format for the p-value
    tip.length = 0.01,  # length of the tooltip
    fontface = "bold",  # font style
    y.position = 0.585,  # y-position for the p-value
    label.size = 2.75,  # size of the label
    step.increase = 0,  # step increase for label positioning
    bracket.shorten = 0,  # shorten bracket length
    vjust = -0.3,  # vertical justification
    remove.bracket = TRUE  # remove bracket around the p-value
  ) +
  # Set y-axis limits
  coord_cartesian(ylim = c(0, 0.7)) +
  # Apply classic theme settings
  theme_classic() +
  theme(
    legend.title = element_blank(),  # remove legend title
    panel.background = element_rect(fill = 'white'),  # set panel background color
    legend.position = "none",  # remove legend
    plot.margin = unit(c(-0.5, 0, 0, 0), 'cm'),  # set plot margins
    axis.text.x = element_markdown(vjust = 0.5, debug = FALSE),  # adjust x-axis text
    plot.subtitle = element_markdown(),  # enable markdown for subtitle
    strip.text = element_text(size = 14)  # Set text size for strip labels
    
  ) +
  # Set x-axis labels using 'labels' vector
  scale_x_discrete("", labels = labels) +
  # Set y-axis label as 'growth rate (d ^ -1)'
  scale_y_continuous(expression("Growth rate (d " ^ -1 * ")")) +
  # Set subtitle using previously defined 'statistic_label'
  labs(subtitle = statistic_label)

#### Add A. pseudogonyaulax strain L4-B1

# Load and transform data ########
data1 <- # cell counts data
  read.csv(
    "C:\\Users\\krist\\OneDrive\\Dokumente\\AWI\\Promotion\\AP1\\N-AP1a\\\\R_txt_files\\L4B1_growth.txt",
    sep = "",
    header = FALSE
  )

data2_4B1 <- # max cell density
  as.data.frame(
    read.csv(
      "C:\\Users\\krist\\OneDrive\\Dokumente\\AWI\\Promotion\\AP1\\N-AP1a\\\\R_txt_files\\L4B1_max_density.txt",
      sep = "",
      header = FALSE
    )
  )

# introduce treatment column = different nitrogen sources with N-deplete being K-medium without nitrogen addition as control
data2_4B1$treat <- rep(c("N_deplete", "NO3", "NH4", "Urea"), each = 3)

# Transpose the data frame and handle missing values
data1 <-
  data1 %>% remove_rownames %>% column_to_rownames(var = "V1") %>% t() %>% as.data.frame()

data1 <- data1 %>% mutate_all(~ ifelse(is.nan(.), NA, .)) %>%
  mutate_all(as.numeric)

## Growth rate calculations ######
# make linear fit for each bottle between 5 data points
# h: amount of data points used for calculation of growth rate!!
fit_C1 <-
  fit_easylinear(na.omit(data1$time_C)[1:7], na.omit(data1$C1)[1:7], h =
                   5)
fit_C2 <-
  fit_easylinear(na.omit(data1$time_C)[1:7], na.omit(data1$C2)[1:7], h =
                   5)
fit_C3 <-
  fit_easylinear(na.omit(data1$time_C)[1:7], na.omit(data1$C3)[1:7], h =
                   5)
fit_NO1 <-
  fit_easylinear(na.omit(data1$time_rest)[1:10], na.omit(data1$NO1)[1:10], h =
                   5)
fit_NO2 <-
  fit_easylinear(na.omit(data1$time_rest)[1:10], na.omit(data1$NO2)[1:10], h =
                   5)
fit_NO3 <-
  fit_easylinear(na.omit(data1$time_rest)[1:10], na.omit(data1$NO3)[1:10], h =
                   5)
fit_NH1 <-
  fit_easylinear(na.omit(data1$time_rest)[1:10], na.omit(data1$NH1)[1:10], h =
                   5)
fit_NH2 <-
  fit_easylinear(na.omit(data1$time_rest)[1:10], na.omit(data1$NH2)[1:10], h =
                   5)
fit_NH3 <-
  fit_easylinear(na.omit(data1[1:9, 1]), na.omit(data1$NH3)[1:9], h = 5)
fit_U1 <-
  fit_easylinear(na.omit(data1[1:8, 1]), na.omit(data1[1:8, 12]), h = 5)
fit_U2 <-
  fit_easylinear(na.omit(data1[1:8, 1]), na.omit(data1[1:8, 13]), h = 5)
fit_U3 <-
  fit_easylinear(na.omit(data1[1:8, 1]), na.omit(data1[1:8, 14]), h = 5)

# extract coefficients of the linear fits 
matrix_coef <-
  as.data.frame(rbind(
    coef(fit_C1),
    coef(fit_C2),
    coef(fit_C3),
    coef(fit_NO1),
    coef(fit_NO2),
    coef(fit_NO3),
    coef(fit_NH1),
    coef(fit_NH2),
    coef(fit_NH3),
    coef(fit_U1),
    coef(fit_U2),
    coef(fit_U3)
  ))

# Assign column names based on the original data
matrix_coef$treat <- colnames(data1[, 3:14])

# construct new dataframe 'ex2' with maximal growth rate and treatments for following analysis
ex2 = data.frame(treat = rep(as.factor(c(
  "N_deplete", "NO3", "NH4", "Urea"
)), each = 3), data = matrix_coef$mumax)

## Check for potential outliers and remove them #####
# 
# outliers <- c()
# 
# for (each_group in unique(ex2$treat)) {
#   dix1 <-
#     if (dixon.test(subset(ex2, ex2$treat == each_group)$data)$p.val < 0.05) {
#       parse_number(dixon.test(subset(ex2, ex2$treat == each_group)$data)$alternative)
#     }
#   if (!is.null(dix1)) {
#     outliers <- c(outliers, dix1)
#   }
#   
# }
# 
# ex2 <- ex2[!(ex2$data %in% outliers), ]
## No outliers found!!

# Statistical Analysis 
# Perform a non-parametric Kruskal test since three samples do not satisfy the normality assumption
# Null hypothesis (H0): The growth rate of A. pseudogonyaulax exposed to different nitrogen sources is equal
res.aov2 <- kruskal_test(data ~ treat, data = ex2)

# Conduct a Conover-Iman non-parametric PostHoc test if p < 0.05, with a liberal p-value adjustment using Benjamini & Hochberg (FDR)
Con_b <-
  conover.test(
    ex2$data,
    ex2$treat,
    method = "BH",
    altp = T,
    alpha = 0.1
  )

# Calculate Cohen's d effect size: ~0.2 (small effect), ~0.5 (medium effect), ~0.8 (large effect)
effsize <- cohens_d(ex2, data ~ treat, hedges.correction = T)
effsize_sub <- abs(round(effsize[c(1:3), 4], digits = 1))

# Extract p-values and group comparisons from the Conover-Iman test output
df_p_val2 <-
  data.frame(
    p.adj = sprintf("%.3f", Con_b$altP.adjusted),
    group1 = sub("-.*", "", Con_b$comparisons),
    group2 = sub(".*-", "", Con_b$comparisons)
  )

df_p_val2 <- df_p_val2 %>% subset(p.adj < 0.1)

df_p_val2[] <-
  sapply(df_p_val2, function(x)
    gsub(pattern = " ", replacement = "", x))

# Check for significant differences in cell densities
# Null hypothesis (H0): The maximum cell densities of A. pseudogonyaulax exposed to different nitrogen sources are equal
Con_b2 <-
  as.data.frame(rbind(ConoverTest(
    V1 ~ treat, data = data2_4B1, method = "BH"
  )[[1]]))

# Prepare ggplot parameters
dodge <- position_dodge(width = 0.4)

# Format the statistic with two digits after the decimal point
formatted_statistic <- sprintf("%.2f", res.aov2$statistic)

# Create subtitle with markdown syntax
statistic_label <-
  glue::glue(
    "<b>b)</b> Kruskal-Wallis, <i>H</i><sub>{res.aov2$df}</sub> = {formatted_statistic}, <i>p</i> = {res.aov2$p}"
  )

# format p-values
df_p_val2$p <- ifelse(df_p_val2$p.adj < 0.01, 0.01, 0.05)

# Create ggplot for L4-B1 dataset
P2 <- ggplot(ex2, aes(x = treat, y = data)) +  # Define ggplot object with x = treatment and y = growth rate
  geom_point(size = 2.0) +  # Add points to the plot with size 2.0
  add_pvalue(  # Add p-value annotations
    df_p_val2[c(1, 5),],  # p-value info for the first comparison
    tip.length = 0.03,  # Set length of the tooltip
    fontface = "bold",  # Set font style
    y.position = 0.45,  # Set y-position for the p-value
    label.size = 0,  # Set size of the label
    step.increase = 0,  # Set step increase for label positioning
    bracket.shorten = 0,  # Shorten bracket length
    vjust = -0.3  # Set vertical justification
  ) + 
  add_pvalue(  # Add p-value annotations
    df_p_val2[2,],  # p-value info for the second comparison
    label = "p < {p}",  # Label format for the p-value
    tip.length = 0.01,  # Set length of the tooltip
    fontface = "bold",  # Set font style
    y.position = 0.485,  # Set y-position for the p-value
    label.size = 2.75,  # Set size of the label
    step.increase = 0.175,  # Set step increase for label positioning
    bracket.shorten = 0,  # Shorten bracket length
    vjust = -0.3  # Set vertical justification
  ) +
  add_pvalue(  # Add p-value annotations
    df_p_val2[3,],  # p-value info for the third comparison
    label = "p = {p.adj}",  # Label format for the p-value
    tip.length = 0.03,  # Set length of the tooltip
    fontface = "bold",  # Set font style
    y.position = 0.615,  # Set y-position for the p-value
    label.size = 0,  # Set size of the label
    step.increase = 0.175,  # Set step increase for label positioning
    bracket.shorten = 0,  # Shorten bracket length
    vjust = -0.3  # Set vertical justification
  ) +
  add_pvalue(  # Add p-value annotations
    df_p_val2[4,],  # p-value info for the fourth comparison
    label = "p < {p}",  # Label format for the p-value
    tip.length = 0.01,  # Set length of the tooltip
    fontface = "bold",  # Set font style
    y.position = 0.65,  # Set y-position for the p-value
    label.size = 2.75,  # Set size of the label
    step.increase = 0.175,  # Set step increase for label positioning
    bracket.shorten = 0,  # Shorten bracket length
    vjust = -0.3  # Set vertical justification
  ) +
  theme_classic() +  # Apply classic theme settings
  theme(  # Customize theme settings
    legend.title = element_blank(),  # Remove legend title
    legend.text = element_markdown(),  # Use markdown for legend text
    panel.background = element_rect(fill = 'white'),  # Set panel background color
    legend.position = "bottom",  # Set legend position to bottom
    strip.text = element_text(size = 14),  # Set text size for strip labels
    plot.margin = unit(c(-0.5, 0, 0, 0), 'cm'),  # Set plot margins
    axis.text.x = element_markdown(vjust = 0.5, debug = F),  # Customize x-axis text
    plot.subtitle = element_markdown()  # Enable markdown for subtitle
  )  +
  scale_x_discrete("N-source", labels = labels) +  # Set x-axis label and custom labels
  coord_cartesian(ylim = c(0, 0.7)) +  # Set y-axis limits
  scale_y_continuous(expression("Growth rate (d " ^ -1 * ")")) +  # Set y-axis label as 'growth rate (d^-1)'
  labs(subtitle = statistic_label)  # Set subtitle using previously defined 'statistic_label'

# combine cell count plots of both strains in one plot
all_plots <-
  P1 + P2 + plot_layout(guides = "collect", ncol = 1, axis_titles = "collect_y") &
  theme(legend.position = "bottom")

# save the combined plot as a png
ggsave(
  "Growth_rate_all.png",
  all_plots,
  path = "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/AP1/N-AP1a/figures",
  dpi = 300,
  width = 3.5,
  height = 4.9,
  units = "in"
)

# Intraspecific statistical analysis of growth rates

# Combine toxin quotas of both A. pseudogonyaulax strains in one dataframe for analysis
ex <- ex %>% mutate(strain = "L2-D2") 
ex2 <- ex2 %>% mutate(strain = "L4-B1")

ex_all <- rbind(ex, ex2)

# conduct kruskal-test on each treatment of each strain
GR_all <- ex_all %>%
  group_by(treat) %>%
  nest() %>%
  mutate(kruskal_results = map(data, ~kruskal_test(data ~ strain, data = .x))) %>%
  unnest(c(kruskal_results)) 

# export mean growth rates and standard deviations for table in paper
agg <- aggregate(ex_all, data ~ treat + strain, FUN = "mean")
agg$sd <- aggregate(ex_all, data ~ treat + strain, FUN = "sd")$data

agg$range <-
  paste(format(round(agg$data, 2), nsmall = 2, digits = 2),
        format(round(agg$sd, 2), nsmall = 2, digits = 2),
        sep = " +/- ")

# Garbage collection: call after large objects have been removed

gc()

# Delete workspace, clean environment/variables and restarte R if used memory piles up
dev.off()

rm(list = ls())

# restart R
.rs.restartR()
