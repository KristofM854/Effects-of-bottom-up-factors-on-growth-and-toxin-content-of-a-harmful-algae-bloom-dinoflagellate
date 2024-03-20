##########################################
## Modeling and plotting of photoirradiance curves of A. pseudogonyaulax strains L2-D2, L4-B1 and L4-B9 exposed to different light intensities
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
  tidyverse,
  plyr,
  ggplot2,
  conover.test,
  growthrates,
  scales,
  installr,
  agricolae,
  ggforce,
  grid,
  dplyr,
  DescTools,
  rstatix,
  ggpmisc,
  acepack,
  agricolae,
  remotes,
  FME,
  lava,
  ggpubr,
  gridExtra,
  flextable,
  ggthemes
)

# Load Windows Fonts and define Times font
loadfonts(device = "win")
windowsFonts(Times = windowsFont("Times"))

##############
# Load data
data <-
  read.csv(
    "C:\\Users\\krist\\OneDrive\\Dokumente\\AWI\\Promotion\\AP1\\Light-AP1b\\raw\\PI_2D2.txt",
    sep = "",
    header = TRUE
  )

# replace all , with . + make columns numeric
data[c(2:84), ] <-
  lapply(data[2:84, ], function(x)
    as.numeric(gsub("," , "." , x)))

## Include PSII concentrations of each replicate
RCII_2D2 <-
  c(rep(mean(c(
    18.9122807, 15.36616702,	12.23394004
  )), each = 3),	rep(mean(c(
    40.90833333,	36.98975045,	46.58045977
  )), each = 3),	rep(mean(c(
    34.75931677,	25.47293364,	32.72371464
  )), each = 3))

## Make data.frame and name each column; introduce grouping factor for each experimental light intensity 20, 100, 200
data1 <-
  data.frame(
    PAR = as.numeric(unlist(data[2:84, c(1, 4, 7, 10, 13, 16, 19, 22, 25)])),
    ETR = as.numeric(unlist(data[2:84, c(3, 6, 9, 12, 15, 18, 21, 24, 27)])),
    NPQ = as.numeric(unlist(data[2:84, c(2, 5, 8, 11, 14, 17, 20, 23, 26)])),
    group = as.factor(c(rep(
      c("20", "100", "200"), each = (length(data$X2) - 1) * 3
    ))),
    RCII = as.numeric(c(rep(
      c(RCII_2D2[1:9]), each = (length(data$X2) - 1)
    )))
  )

# Calculate ETR per cell as multiplication of the absolute ETR (calculated in Excel file) with the PSII concentration and convert to fmol ETR per cell
data1$ETR_per_cell <- data1$ETR * data1$RCII * 10^-3

# Calculate mean and SD for each experimental light intensity at each light level of the FRRf (PAR)
# both for ETR and ETR per cell 
data1$PAR <- round(data1$PAR, digits = 0)

agg <- aggregate(ETR ~ PAR + group, data = data1, FUN = "mean")
agg$group <- factor(agg$group, levels = c("20", "100", "200"))
agg$sd <- aggregate(ETR ~ PAR + group, data = data1, FUN = "sd")$ETR

agg$ETR_per_cell <-
  aggregate(ETR_per_cell ~ PAR + group, data = data1, FUN = "mean")$ETR_per_cell
agg$ETR_per_cell_sd <-
  aggregate(ETR_per_cell ~ PAR + group, data = data1, FUN = "sd")$ETR_per_cell

# Subset light intensity 20 for modelling - first model with average points for plot; after each replicate separately for statistics
data1_sub <- subset(data1, data1$group == "20")
PI <- aggregate(data1_sub, ETR ~ PAR, FUN = "mean")

## Model and Modfit  as described in Silsbe and Kromkamp 2012  https://doi.org/10.4319/lom.2012.10.645
## Platt model:
## get starting value for alpha --> slope of linear PAR vs ETR #########
# PI2 <- PI[0:4,]
# ggplot(PI2,aes(x=PAR,y=ETR)) +
#   geom_point() +
#   stat_poly_line( se=F, fullrange = T) +
#   stat_poly_eq(aes(label=paste(..eq.label..,..adj.rr.label..,
#                                sep= "*plain(\",\")~")))

# par = individual light levels of the FRRf; parms = starting values; alpha from slope of linear PAR vs ETR
# ps as maximum observed ETR; beta between 0 and 1
par <- PI$PAR
par2 <- seq(0, 2000, 0.1)
parms_20f = c(alpha = 2.8, beta = 0.35, ps = 2000)

# set up model according to Platt; ModelCost returns the residuals; modFit finds local residual minimum
model_20f <-
  function(parms_20f, par)
    with(as.list(parms_20f), return(ps * (1 - exp((
      -1 * alpha * par
    ) / ps)) * exp((-1 * beta * par) / ps)))

ModelCost_20f <- function(P) {
  out <- model_20f(P, par)
  return(PI$ETR - out)  # residuals
}

fit_20f <-
  modFit(
    f = ModelCost_20f,
    p = parms_20f,
    method = "Nelder-Mead",
    lower = c(0, 0, 0),
    upper = c(Inf, 1, Inf)
  )

## plot experimental data points and model as line
#
# plot(PI$PAR,PI$ETR)
# lines(par2, model_20f(fit_20f$par, par2), lwd = 2, col = "red")
#
# #model statistics
# summary(fit_20f)
# plot(fit_20f) #residual plot

## extract model parameters with confidence intervals and residual standard error
# parameters_20_platt <- as.data.frame(summary(fit_20f)$par)
# parameters_20_platt$model <- "platt_20"
# parameters_20_platt$ssr <- fit_20f$ssr

# Subset light intensity 20 for modelling
data1_sub <- subset(data1, data1$group == "20")
data1_sub$grouping <- rep(c(1:3), each = length(data1_sub$PAR) / 3)

# Aggregate ETR over the different light intensities; grouping = triplicates of L2-D2
PI <- aggregate(data1_sub, ETR ~ PAR + grouping, FUN = "mean")

# PI for Pangea export
PI_pan <-
  aggregate(data1_sub, ETR ~ PAR + grouping, FUN = "mean") %>% mutate(treat = 20, strain = "L2-D2")

# PI subset for each triplicate
PI <- subset(PI, PI$grouping == "1")

##Model and Modfit  as described in Silsbe and Kromkamp 2012  https://doi.org/10.4319/lom.2012.10.645
##Platt model:
##get starting value for alpha --> slope of linear PAR vs ETR #########
#
# PI2 <- PI[0:4,]
# ggplot(PI2,aes(x=PAR,y=ETR)) +
#   geom_point() +
#   stat_poly_line( se=F, fullrange = T) +
#   stat_poly_eq(aes(label=paste(..eq.label..,..adj.rr.label..,
#                                sep= "*plain(\",\")~")))

# par = individual light levels of the FRRf; parms = starting values; alpha from slope of linear PAR vs ETR
# ps as maximum observed ETR; beta between 0 and 1
par <- PI$PAR
par2 <- seq(0, 2000, 0.1)
parms_20 = c(alpha = 1.98, beta = 0.5, ps = 1400)

##set up model according to Platt; ModelCost returns the residuals; modFit finds local residual minimum
model_20 <-
  function(parms_20, par)
    with(as.list(parms_20), return(ps * (1 - exp((
      -1 * alpha * par
    ) / ps)) * exp((-1 * beta * par) / ps)))

ModelCost_20 <- function(P) {
  out <- model_20(P, par)
  return(PI$ETR - out)  # residuals
}

fit_20a1 <-
  modFit(
    f = ModelCost_20,
    p = parms_20,
    method = "Nelder-Mead",
    lower = c(0, 0, 0),
    upper = c(Inf, 1, Inf)
  )

##plot experimental data points and model as line
# plot(PI$PAR,PI$ETR)
# lines(par2, model_20(fit_20a1$par, par2), lwd = 2, col = "red")
#
# ##model statistics
# summary(fit_20a1)
# plot(fit_20a1) #residual plot
#
# ##extract model parameters with confidence intervals and residual standard error
# parameters_20_platt1 <- as.data.frame(summary(fit_20a1)$par)
# parameters_20_platt1$model <- "platt_20"
# parameters_20_platt1$ssr <- fit_20a1$ssr

# PI subset for each triplicate
PI <- aggregate(data1_sub, ETR ~ PAR + grouping, FUN = "mean")
PI <- subset(PI, PI$grouping == "2")

## Model and Modfit  as described in Silsbe and Kromkamp 2012  https://doi.org/10.4319/lom.2012.10.645
## Platt model:
## get starting value for alpha --> slope of linear PAR vs ETR #########
#
# PI2 <- PI[0:4,]
# ggplot(PI2,aes(x=PAR,y=ETR)) +
#   geom_point() +
#   stat_poly_line( se=F, fullrange = T) +
#   stat_poly_eq(aes(label=paste(..eq.label..,..adj.rr.label..,
#                                sep= "*plain(\",\")~")))

# par = individual light levels of the FRRf; parms = starting values; alpha from slope of linear PAR vs ETR
# ps as maximum observed ETR; beta between 0 and 1
par <- PI$PAR
par2 <- seq(0, 2000, 0.1)
parms_20 = c(alpha = 2.55, beta = 0.5, ps = 2500)

# set up model according to Platt; ModelCost returns the residuals; modFit finds local residual minimum
model_20 <-
  function(parms_20, par)
    with(as.list(parms_20), return(ps * (1 - exp((
      -1 * alpha * par
    ) / ps)) * exp((-1 * beta * par) / ps)))

ModelCost_20 <- function(P) {
  out <- model_20(P, par)
  return(PI$ETR - out)  # residuals
}

fit_20a2 <-
  modFit(
    f = ModelCost_20,
    p = parms_20,
    method = "Nelder-Mead",
    lower = c(0, 0, 0),
    upper = c(Inf, 1, Inf)
  )

##plot experimental data points and model as line
# plot(PI$PAR,PI$ETR)
# lines(par2, model_20(fit_20a2$par, par2), lwd = 2, col = "red")
#
# ##model statistics
# summary(fit_20a2)
# plot(fit_20a2) #residual plot
#
# ##extract model parameters with confidence intervals and residual standard error
# parameters_20_platt2 <- as.data.frame(summary(fit_20a2)$par)
# parameters_20_platt2$model <- "platt_20"
# parameters_20_platt2$ssr <- fit_20a2$ssr

# PI subset for each triplicate
PI <- aggregate(data1_sub, ETR ~ PAR + grouping, FUN = "mean")
PI <- subset(PI, PI$grouping == "3")

## Model and Modfit  as described in Silsbe and Kromkamp 2012  https://doi.org/10.4319/lom.2012.10.645
## Platt model:
## get starting value for alpha --> slope of linear PAR vs ETR #########
#
# PI2 <- PI[0:4,]
# ggplot(PI2,aes(x=PAR,y=ETR)) +
#   geom_point() +
#   stat_poly_line( se=F, fullrange = T) +
#   stat_poly_eq(aes(label=paste(..eq.label..,..adj.rr.label..,
#                                sep= "*plain(\",\")~")))

# par = individual light levels of the FRRf; parms = starting values; alpha from slope of linear PAR vs ETR
# ps as maximum observed ETR; beta between 0 and 1
par <- PI$PAR
par2 <- seq(0, 2000, 0.1)
parms_20 = c(alpha = 2.8, beta = 0.1, ps = 2000)

# set up model according to Platt; ModelCost returns the residuals; modFit finds local residual minimum
model_20 <-
  function(parms_20, par)
    with(as.list(parms_20), return(ps * (1 - exp((
      -1 * alpha * par
    ) / ps)) * exp((-1 * beta * par) / ps)))

ModelCost_20 <- function(P) {
  out <- model_20(P, par)
  return(PI$ETR - out)  # residuals
}

fit_20a3 <-
  modFit(
    f = ModelCost_20,
    p = parms_20,
    method = "Nelder-Mead",
    lower = c(0, 0, 0),
    upper = c(Inf, 1, Inf)
  )

##plot experimental data points and model as line
# plot(PI$PAR,PI$ETR)
# lines(par2, model_20(fit_20a3$par, par2), lwd = 2, col = "red")
#
# ##model statistics
# summary(fit_20a3)
# plot(fit_20a3) #residual plot
#
# ##extract model parameters with confidence intervals and residual standard error
# parameters_20_platt3 <- as.data.frame(summary(fit_20a3)$par)
# parameters_20_platt3$model <- "platt_20"
# parameters_20_platt3$ssr <- fit_20a3$ssr

## Webb model; starting parameters randomly
# par <- PI$PAR
# parms = c(alpha=3,Ek=470)
#
# ## set up model according to Webb; ModelCost returns the residuals; modFit finds local residual minimum
# model <- function(parms,par) with (as.list(parms), return(alpha*Ek*(1-exp(-1*par/Ek))))
#
# ModelCost <- function(P) {
#   out <- model(P, par)
#   return(PI$ETR-out)  # residuals
# }
#
# fit_20 <- modFit(f = ModelCost, p = parms,method = "Nelder-Mead")
#
# ##plot experimental data points and model as line
# plot(PI$PAR,PI$ETR)
# lines(par2, model(fit_20$par, par2), lwd = 2, col = "red")
#
# ## model statistics
#
# summary(fit_20)
# plot(fit_20)
#
# ## extract model parameters and confidence intervals
# parameters_20_webb <- as.data.frame(summary(fit_20)$par)
# parameters_20_webb$model <- "webb_20"
# parameters_20_webb$ssr <- fit_20$ssr
#
# ## Jasby and Platt model; starting parameters randomly ####
# par <- PI$PAR
# parms = c(alpha=2.3,Ek=590)
#
# ## set up model according to Jasby and Platt; ModelCost returns the residuals; modFit finds local residual minimum
# model <- function(parms,par) with (as.list(parms), return(alpha*Ek*tanh(par/Ek)))
#
# ModelCost <- function(P) {
#   out <- model(P, par)
#   return(PI$ETR-out)  # residuals
# }
#
# fit_20 <- modFit(f = ModelCost, p = parms,method = "Nelder-Mead")
#
# ## plot experimental data points and model as line
# plot(PI$PAR,PI$ETR)
# lines(par2, model(fit_20$par, par2), lwd = 2, col = "red")
#
# ## model statistics
# summary(fit_20)
# plot(fit_20)
#
# ## extract model parameters and confidence intervals
# parameters_20_jasby <- as.data.frame(summary(fit_20)$par)
# parameters_20_jasby$model <- "jasby_20"
# parameters_20_jasby$ssr <- fit_20$ssr

### 100er light treatment
data1_sub <- subset(data1, data1$group == "100")
PI <- aggregate(data1_sub, ETR ~ PAR, FUN = "mean")

## Model and Modfit  as described in Silsbe and Kromkamp 2012  https://doi.org/10.4319/lom.2012.10.645
## Platt model:
## get starting value for alpha --> slope of linear PAR vs ETR #########
#
# PI2 <- PI[0:4,]
# ggplot(PI2,aes(x=PAR,y=ETR)) +
#   geom_point() +
#   stat_poly_line( se=F, fullrange = T) +
#   stat_poly_eq(aes(label=paste(..eq.label..,..adj.rr.label..,
#                                sep= "*plain(\",\")~")))

## par = individual light levels of the FRRf; parms_100 = starting values; alpha from slope of linear PAR vs ETR
## ps as maximum observed ETR; beta between 0 and 1
par <- PI$PAR
par2 <- seq(0, 2000, 0.1)
parms_100f = c(alpha = 3.8, beta = 0.1, ps = 795)

# set up model_100 according to Platt; ModelCost_100 returns the residuals; modFit finds local residual minimum
model_100f <-
  function(parms_100f, par)
    with(as.list(parms_100f),  return(ps * ((1 - exp((-1 * alpha * par) /
                                                        ps
    )) * exp((
      -1 * beta * par
    ) / ps))))

ModelCost_100f <- function(P) {
  out <- model_100f(P, par)
  return(PI$ETR - out)  # residuals
}

fit_100f <-
  modFit(
    f = ModelCost_100f,
    p = parms_100f,
    method = "Nelder-Mead",
    lower = c(0, 0, 0),
    upper = c(Inf, 1, Inf)
  )

## plot experimental data points and model as line
#
# plot(PI$PAR,PI$ETR)
# lines(par2, model_100f(fit_100f$par, par2), lwd = 2, col = "red")
#
# ##model statistics
# summary(fit_100f)
# plot(fit_100f) #residual plot
#
# ##extract model parameters with confidence intervals and residual standard error
# parameters_100_platt <- as.data.frame(summary(fit_100f)$par)
# parameters_100_platt$model <- "platt_100"
# parameters_100_platt$ssr <- fit_100f$ssr

# PI subset for each triplicate
data1_sub$grouping <- rep(c(1:3), each = length(data1_sub$group) / 3)
data1_sub$PAR <- round(data1_sub$PAR, digits = 0)
PI <- aggregate(data1_sub, ETR ~ PAR + grouping, FUN = "mean")

# PI for Pangea export
PI_pan2 <-
  aggregate(data1_sub, ETR ~ PAR + grouping, FUN = "mean")  %>% mutate(treat = 100, strain = "L2-D2")

PI <- subset(PI, PI$grouping == "1")

## Model and Modfit  as described in Silsbe and Kromkamp 2012  https://doi.org/10.4319/lom.2012.10.645
## Platt model:
## get starting value for alpha --> slope of linear PAR vs ETR #########
#
# PI2 <- PI[0:4,]
# ggplot(PI2,aes(x=PAR,y=ETR)) +
#   geom_point() +
#   stat_poly_line( se=F, fullrange = T) +
#   stat_poly_eq(aes(label=paste(..eq.label..,..adj.rr.label..,
#                                sep= "*plain(\",\")~")))

# par = individual light levels of the FRRf; parms_100 = starting values; alpha from slope of linear PAR vs ETR
# ps as maximum observed ETR; beta between 0 and 1
par <- PI$PAR
par2 <- seq(0, 2000, 0.1)
parms_100 = c(alpha = 3.8, beta = 0.1, ps = 795)

# set up model_100 according to Platt; ModelCost_100 returns the residuals; modFit finds local residual minimum
model_100 <-
  function(parms_100, par)
    with(as.list(parms_100),  return(ps * ((1 - exp((-1 * alpha * par) / ps
    )) * exp((
      -1 * beta * par
    ) / ps))))

ModelCost_100 <- function(P) {
  out <- model_100(P, par)
  return(PI$ETR - out)  # residuals
}

fit_100a1 <-
  modFit(
    f = ModelCost_100,
    p = parms_100,
    method = "Nelder-Mead",
    lower = c(0, 0, 0),
    upper = c(Inf, 1, Inf)
  )

## plot experimental data points and model as line
#
# plot(PI$PAR,PI$ETR)
# lines(par2, model_100(fit_100a1$par, par2), lwd = 2, col = "red")
#
# ##model statistics
# summary(fit_100a1)
# plot(fit_100a1) #residual plot
#
# ##extract model parameters with confidence intervals and residual standard error
# parameters_100_platt1 <- as.data.frame(summary(fit_100a1)$par)
# parameters_100_platt1$model <- "platt_100"
# parameters_100_platt1$ssr <- fit_100a1$ssr

# PI subset for each triplicate
PI <- aggregate(data1_sub, ETR ~ PAR + grouping, FUN = "mean")
PI <- subset(PI, PI$grouping == "2")

## Model and Modfit  as described in Silsbe and Kromkamp 2012  https://doi.org/10.4319/lom.2012.10.645
## Platt model:
## get starting value for alpha --> slope of linear PAR vs ETR #########
#
# PI2 <- PI[0:4,]
# ggplot(PI2,aes(x=PAR,y=ETR)) +
#   geom_point() +
#   stat_poly_line( se=F, fullrange = T) +
#   stat_poly_eq(aes(label=paste(..eq.label..,..adj.rr.label..,
#                                sep= "*plain(\",\")~")))

## par = individual light levels of the FRRf; parms_100 = starting values; alpha from slope of linear PAR vs ETR
## ps as maximum observed ETR; beta between 0 and 1
par <- PI$PAR
par2 <- seq(0, 2000, 0.1)
parms_100 = c(alpha = 3.8, beta = 0.1, ps = 795)

##set up model_100 according to Platt; ModelCost_100 returns the residuals; modFit finds local residual minimum
model_100 <-
  function(parms_100, par)
    with(as.list(parms_100),  return(ps * ((1 - exp((-1 * alpha * par) / ps
    )) * exp((
      -1 * beta * par
    ) / ps))))

ModelCost_100 <- function(P) {
  out <- model_100(P, par)
  return(PI$ETR - out)  # residuals
}

fit_100a2 <-
  modFit(
    f = ModelCost_100,
    p = parms_100,
    method = "Nelder-Mead",
    lower = c(0, 0, 0),
    upper = c(Inf, 1, Inf)
  )

## plot experimental data points and model as line
#
# plot(PI$PAR,PI$ETR)
# lines(par2, model_100(fit_100a2$par, par2), lwd = 2, col = "red")
#
# ## model statistics
# summary(fit_100a2)
# plot(fit_100a2) #residual plot
#
# ## extract model parameters with confidence intervals and residual standard error
# parameters_100_platt2 <- as.data.frame(summary(fit_100a2)$par)
# parameters_100_platt2$model <- "platt_100"
# parameters_100_platt2$ssr <- fit_100a2$ssr

## PI subset for each triplicate
PI <- aggregate(data1_sub, ETR ~ PAR + grouping, FUN = "mean")
PI <- subset(PI, PI$grouping == "3")

## Model and Modfit  as described in Silsbe and Kromkamp 2012  https://doi.org/10.4319/lom.2012.10.645
## Platt model:
## get starting value for alpha --> slope of linear PAR vs ETR #########
#
# PI2 <- PI[0:4,]
# ggplot(PI2,aes(x=PAR,y=ETR)) +
#   geom_point() +
#   stat_poly_line( se=F, fullrange = T) +
#   stat_poly_eq(aes(label=paste(..eq.label..,..adj.rr.label..,
#                                sep= "*plain(\",\")~")))

## par = individual light levels of the FRRf; parms_100 = starting values; alpha from slope of linear PAR vs ETR
## ps as maximum observed ETR; beta between 0 and 1
par <- PI$PAR
par2 <- seq(0, 2000, 0.1)
parms_100 = c(alpha = 3.8, beta = 0.1, ps = 795)

##set up model_100 according to Platt; ModelCost_100 returns the residuals; modFit finds local residual minimum
model_100 <-
  function(parms_100, par)
    with(as.list(parms_100),  return(ps * ((1 - exp((-1 * alpha * par) / ps
    )) * exp((
      -1 * beta * par
    ) / ps))))

ModelCost_100 <- function(P) {
  out <- model_100(P, par)
  return(PI$ETR - out)  # residuals
}

fit_100a3 <-
  modFit(
    f = ModelCost_100,
    p = parms_100,
    method = "Nelder-Mead",
    lower = c(0, 0, 0),
    upper = c(Inf, 1, Inf)
  )

## plot experimental data points and model as line
#
# plot(PI$PAR,PI$ETR)
# lines(par2, model_100(fit_100a3$par, par2), lwd = 2, col = "red")
#
# ##model statistics
# summary(fit_100a3)
# plot(fit_100a3) #residual plot
#
# ##extract model parameters with confidence intervals and residual standard error
# parameters_100_platt2 <- as.data.frame(summary(fit_100a3)$par)
# parameters_100_platt2$model <- "platt_100"
# parameters_100_platt2$ssr <- fit_100a3$ssr

## Webb model; starting parameters randomly
# par <- PI$PAR
# parms = c(alpha=3.66,Ek=230)
#
# ## set up model according to Webb; ModelCost returns the residuals; modFit finds local residual minimum
# model <- function(parms,par) with (as.list(parms), return(alpha*Ek*(1-exp(-1*par/Ek))))
#
# ModelCost <- function(P) {
#   out <- model(P, par)
#   return(PI$ETR-out)  # residuals
# }
#
# fit_100 <- modFit(f = ModelCost, p = parms,method = "Nelder-Mead")
#
#
# ## plot experimental data points and model as line
# plot(PI$PAR,PI$ETR)
# lines(par, model(fit_100$par, par), lwd = 2, col = "red")
#
# ## model statistics
# summary(fit_100)
# plot(fit_100)
#
# ## extract model parameters and confidence intervals
# parameters_100_webb <- as.data.frame(summary(fit_100)$par)
# parameters_100_webb$model <- "webb_100"
# parameters_100_webb$ssr <- fit_100$ssr
#
# ## Jasby and Platt model; starting parameters randomly ####
# par <- PI$PAR
# parms = c(alpha=2.77,Ek=299)
#
# ## set up model according to Jasby and Platt; ModelCost returns the residuals; modFit finds local residual minimum
# model <- function(parms,par) with (as.list(parms), return(alpha*Ek*tanh(par/Ek)))
#
# ModelCost <- function(P) {
#   out <- model(P, par)
#   return(PI$ETR-out)  # residuals
# }
#
# fit_100 <- modFit(f = ModelCost, p = parms,method = "Nelder-Mead")
#
# ## plot experimental data points and model as line
# plot(PI$PAR,PI$ETR)
# lines(par, model(fit_100$par, par), lwd = 2, col = "red")

## model statistics
# summary(fit_100)
# plot(fit_100)

## extract model parameters and confidence intervals
# parameters_100_jasby <- as.data.frame(summary(fit_100)$par)
# parameters_100_jasby$model <- "jasby_100"
# parameters_100_jasby$ssr <- fit_100$ssr

### 200er light treatment

data1_sub <- subset(data1, data1$group == "200")

PI <- aggregate(data1_sub, ETR ~ PAR, FUN = "mean")

## Model and Modfit  as described in Silsbe and Kromkamp 2012  https://doi.org/10.4319/lom.2012.10.645
## Platt model:
## get starting value for alpha --> slope of linear PAR vs ETR #########
#
# PI2 <- PI[0:4,]
# ggplot(PI2,aes(x=PAR,y=ETR)) +
#   geom_point() +
#   stat_poly_line( se=F, fullrange = T) +
#   stat_poly_eq(aes(label=paste(..eq.label..,..adj.rr.label..,
#                                sep= "*plain(\",\")~")))

## par = individual light levels of the FRRf; parms_200 = starting values; alpha from slope of linear PAR vs ETR
## ps as maximum observed ETR; beta between 0 and 1
par <- PI$PAR
par2 <- seq(0, 2000, 0.1)
parms_200f = c(alpha = 1.96, beta = 0.42, ps = 2000)

## set up model_200 according to Platt; ModelCost_200 returns the residuals; modFit finds local residual minimum
model_200f <-
  function(parms_200f, par)
    with(as.list(parms_200f),  return(ps * ((1 - exp((-1 * alpha * par) /
                                                        ps
    )) * exp((
      -1 * beta * par
    ) / ps))))

ModelCost_200f <- function(P) {
  out <- model_200f(P, par)
  return(PI$ETR - out)  # residuals
}

fit_200f <-
  modFit(
    f = ModelCost_200f,
    p = parms_200f,
    method = "Nelder-Mead",
    lower = c(0, 0, 0),
    upper = c(Inf, 1, Inf)
  )

## plot experimental data points and model as line
#
# plot(PI$PAR,PI$ETR)
# lines(par2, model_200f(fit_200f$par, par2), lwd = 2, col = "red")
#
# ##model statistics
# summary(fit_200f)
# plot(fit_200f) #residual plot
#
# ## extract model parameters with confidence intervals and residual standard error
# parameters_200_platt <- as.data.frame(summary(fit_200f)$par)
# parameters_200_platt$model <- "platt_200_1"
# parameters_200_platt$ssr <- fit_200f$ssr

# introduce grouping variable for the triplicate. Modeling has to be done separately to get standard deviations of parameters
data1_sub$grouping <- rep(c(1:3), each = length(data1_sub$group) / 3)
data1_sub$PAR <- round(data1_sub$PAR, digits = 0)
PI <- aggregate(data1_sub, ETR ~ PAR + grouping, FUN = "mean")

# PI for Pangea export
PI_pan3 <-
  aggregate(data1_sub, ETR ~ PAR + grouping, FUN = "mean") %>% mutate(treat = 200, strain = "L2-D2")

# Assuming your dataframe is named data2
PI_pan_L2D2 <- rbind(PI_pan, PI_pan2, PI_pan3)

# PI subset for each triplicate
PI <- subset(PI, PI$grouping == "1")

## Model and Modfit  as described in Silsbe and Kromkamp 2012  https://doi.org/10.4319/lom.2012.10.645
## Platt model:
## get starting value for alpha --> slope of linear PAR vs ETR #########
#
# PI2 <- PI[0:4,]
# ggplot(PI2,aes(x=PAR,y=ETR)) +
#   geom_point() +
#   stat_poly_line( se=F, fullrange = T) +
#   stat_poly_eq(aes(label=paste(..eq.label..,..adj.rr.label..,
#                                sep= "*plain(\",\")~")))

## par = individual light levels of the FRRf; parms_200 = starting values; alpha from slope of linear PAR vs ETR
## ps as maximum observed ETR; beta between 0 and 1
par <- PI$PAR
par2 <- seq(0, 2000, 0.1)
parms_200 = c(alpha = 3.5, beta = 0.1, ps = 1200)

# set up model_200 according to Platt; ModelCost_200 returns the residuals; modFit finds local residual minimum
model_200 <-
  function(parms_200, par)
    with(as.list(parms_200),  return(ps * ((1 - exp((-1 * alpha * par) / ps
    )) * exp((
      -1 * beta * par
    ) / ps))))

ModelCost_200 <- function(P) {
  out <- model_200(P, par)
  return(PI[1:10, 3] - out)  # residuals
}

fit_200a1 <-
  modFit(
    f = ModelCost_200,
    p = parms_200,
    method = "Nelder-Mead",
    lower = c(0, 0, 0),
    upper = c(Inf, 1, Inf)
  )

## plot experimental data points and model as line
#
# plot(PI$PAR,PI$ETR)
# lines(par2, model_200(fit_200a1$par, par2), lwd = 2, col = "red")
#
# ## model statistics
# summary(fit_200a1)
# plot(fit_200a1) #residual plot
#
# ## extract model parameters with confidence intervals and residual standard error
# parameters_200_platt1 <- as.data.frame(summary(fit_200a1)$par)
# parameters_200_platt1$model <- "platt_200_1"
# parameters_200_platt1$ssr <- fit_200a1$ssr

# PI subset for each triplicate
PI <- aggregate(data1_sub, ETR ~ PAR + grouping, FUN = "mean")
PI <- subset(PI, PI$grouping == "2")

## Model and Modfit  as described in Silsbe and Kromkamp 2012  https://doi.org/10.4319/lom.2012.10.645
## Platt model:
## get starting value for alpha --> slope of linear PAR vs ETR #########
# PI2 <- PI[0:4,]
# ggplot(PI2,aes(x=PAR,y=ETR)) +
#   geom_point() +
#   stat_poly_line( se=F, fullrange = T) +
#   stat_poly_eq(aes(label=paste(..eq.label..,..adj.rr.label..,
#                                sep= "*plain(\",\")~")))

## par = individual light levels of the FRRf; parms_200 = starting values; alpha from slope of linear PAR vs ETR
## ps as maximum observed ETR; beta between 0 and 1
par <- PI$PAR
par2 <- seq(0, 2000, 0.1)
parms_200 = c(alpha = 3.5, beta = 0.9, ps = 2000)

# set up model_200 according to Platt; ModelCost_200 returns the residuals; modFit finds local residual minimum
model_200 <-
  function(parms_200, par)
    with(as.list(parms_200),  return(ps * ((1 - exp((-1 * alpha * par) / ps
    )) * exp((
      -1 * beta * par
    ) / ps))))

ModelCost_200 <- function(P) {
  out <- model_200(P, par)
  return(PI[1:10, 3] - out)  # residuals
}

fit_200a2 <-
  modFit(
    f = ModelCost_200,
    p = parms_200,
    method = "Nelder-Mead",
    lower = c(0, 0, 0),
    upper = c(Inf, 1, Inf)
  )

## plot experimental data points and model as line
#
# plot(PI$PAR,PI$ETR)
# lines(par2, model_200(fit_200a2$par, par2), lwd = 2, col = "red")
#
# ## model statistics
# summary(fit_200a2)
# plot(fit_200a2) #residual plot
#
# ## extract model parameters with confidence intervals and residual standard error
# parameters_200_platt2 <- as.data.frame(summary(fit_200a2)$par)
# parameters_200_platt2$model <- "platt_200_2"
# parameters_200_platt2$ssr <- fit_200a2$ssr

# PI subset for each triplicate
PI <- aggregate(data1_sub, ETR ~ PAR + grouping, FUN = "mean")
PI <- subset(PI, PI$grouping == "3")

## Model and Modfit  as described in Silsbe and Kromkamp 2012  https://doi.org/10.4319/lom.2012.10.645
## Platt model:
## get starting value for alpha --> slope of linear PAR vs ETR #########
# PI2 <- PI[0:4,]
# ggplot(PI2,aes(x=PAR,y=ETR)) +
#   geom_point() +
#   stat_poly_line( se=F, fullrange = T) +
#   stat_poly_eq(aes(label=paste(..eq.label..,..adj.rr.label..,
#                                sep= "*plain(\",\")~")))

# par = individual light levels of the FRRf; parms_200 = starting values; alpha from slope of linear PAR vs ETR
# ps as maximum observed ETR; beta between 0 and 1
par <- PI$PAR
par2 <- seq(0, 2000, 0.1)
parms_200 = c(alpha = 3.5, beta = 0.9, ps = 2500)

# set up model_200 according to Platt; ModelCost_200 returns the residuals; modFit finds local residual minimum
model_200 <-
  function(parms_200, par)
    with(as.list(parms_200),  return(ps * ((1 - exp((-1 * alpha * par) / ps
    )) * exp((
      -1 * beta * par
    ) / ps))))

ModelCost_200 <- function(P) {
  out <- model_200(P, par)
  return(PI[1:10, 3] - out)  # residuals
}

fit_200a3 <-
  modFit(
    f = ModelCost_200,
    p = parms_200,
    method = "Nelder-Mead",
    lower = c(0, 0, 0),
    upper = c(Inf, 1, Inf)
  )

## plot experimental data points and model as line
#
# plot(PI$PAR,PI$ETR)
# lines(par2, model_200(fit_200a3$par, par2), lwd = 2, col = "red")
#
# ##model statistics
# summary(fit_200a3)
# plot(fit_200a3) #residual plot
#
# ##extract model parameters with confidence intervals and residual standard error
# parameters_200_platt3 <- as.data.frame(summary(fit_200a3)$par)
# parameters_200_platt3$model <- "platt_200_3"
# parameters_200_platt3$ssr <- fit_200a3$ssr


## Webb model; starting parameters randomly
# par <- PI$PAR
# parms = c(alpha=3.9,Ek=356)
#
# ## set up model according to Webb; ModelCost returns the residuals; modFit finds local residual minimum
# model <- function(parms,par) with (as.list(parms), return(alpha*Ek*(1-exp(-1*par/Ek))))
#
# ModelCost <- function(P) {
#   out <- model(P, par)
#   return(PI$ETR-out)  # residuals
# }
#
# fit_200 <- modFit(f = ModelCost, p = parms,method = "Nelder-Mead")
#
#
# ## plot experimental data points and model as line
# plot(PI$PAR,PI$ETR)
# lines(par, model(fit_200$par, par), lwd = 2, col = "red")
#
# ## model statistics
# summary(fit_200)
# plot(fit_200)
#
# ## extract model parameters and confidence intervals
# parameters_200_webb <- as.data.frame(summary(fit_200)$par)
# parameters_200_webb$model <- "webb_200"
# parameters_200_webb$ssr <- fit_200$ssr
#
# ## Jasby and Platt model; starting parameters randomly ####
# par <- PI$PAR
# parms = c(alpha=3.05,Ek=410)
#
# ## set up model according to Jasby and Platt; ModelCost returns the residuals; modFit finds local residual minimum
# model <- function(parms,par) with (as.list(parms), return(alpha*Ek*tanh(par/Ek)))
#
# ModelCost <- function(P) {
#   out <- model(P, par)
#   return(PI$ETR-out)  # residuals
# }
#
# fit_200 <- modFit(f = ModelCost, p = parms,method = "Nelder-Mead")
#
# ## plot experimental data points and model as line
# plot(PI$PAR,PI$ETR)
# lines(par, model(fit_200$par, par), lwd = 2, col = "red")

## model statistics
# summary(fit_200)
# plot(fit_200)

## extract model parameters and confidence intervals
# parameters_200_jasby <- as.data.frame(summary(fit_200)$par)
# parameters_200_jasby$model <- "jasby_200"
# parameters_200_jasby$ssr <- fit_200$ssr

# calculate maximum electron transport rate and minimum light saturation according to Oxborough et al. 2012
# 20 treatment
# replicate 1
a_20_1 = round(as.numeric(fit_20a1$par[1]), digits = 3)
b1 = round(as.numeric(fit_20a1$par[2]), digits = 3)
ps1 = round(as.numeric(fit_20a1$par[3]), digits = 3)

absETRmax_20_1 = ps1 * (a_20_1 / (a_20_1 + b1)) * (b1 / (a_20_1 + b1)) ^
  (b1 / a_20_1)
I_20_1 = absETRmax_20_1 / a_20_1

# replicate 2
a_20_2 = round(as.numeric(fit_20a2$par[1]), digits = 3)
b1 = round(as.numeric(fit_20a2$par[2]), digits = 3)
ps1 = round(as.numeric(fit_20a2$par[3]), digits = 3)

absETRmax_20_2 = ps1 * (a_20_2 / (a_20_2 + b1)) * (b1 / (a_20_2 + b1)) ^
  (b1 / a_20_2)
I_20_2 = absETRmax_20_2 / a_20_2

# replicate 3
a_20_3 = round(as.numeric(fit_20a3$par[1]), digits = 3)
b1 = round(as.numeric(fit_20a3$par[2]), digits = 3)
ps1 = round(as.numeric(fit_20a3$par[3]), digits = 3)

absETRmax_20_3 = ps1 * (a_20_3 / (a_20_3 + b1)) * (b1 / (a_20_3 + b1)) ^
  (b1 / a_20_3)
I_20_3 = absETRmax_20_3 / a_20_3

absETR_20_mean <-
  mean(c(absETRmax_20_1, absETRmax_20_2, absETRmax_20_3))
absETR_20_sd <- sd(c(absETRmax_20_1, absETRmax_20_2, absETRmax_20_3))

I_20_mean <- mean(c(I_20_1, I_20_2, I_20_3))
I_20_sd <- sd(c(I_20_1, I_20_2, I_20_3))

a_20_mean <- mean(c(a_20_1, a_20_2, a_20_3))
a_20_sd <- sd(c(a_20_1, a_20_2, a_20_3))

# 100 treatment

a_100_1 = round(as.numeric(fit_100a1$par[1]), digits = 3)
b1 = round(as.numeric(fit_100a1$par[2]), digits = 3)
ps1 = round(as.numeric(fit_100a1$par[3]), digits = 3)

absETRmax_100_1 = ps1 * (a_100_1 / (a_100_1 + b1)) * (b1 / (a_100_1 + b1)) ^
  (b1 / a_100_1)
I_100_1 = absETRmax_100_1 / a_100_1

a_100_2 = round(as.numeric(fit_100a2$par[1]), digits = 3)
b1 = round(as.numeric(fit_100a2$par[2]), digits = 3)
ps1 = round(as.numeric(fit_100a2$par[3]), digits = 3)

absETRmax_100_2 = ps1 * (a_100_2 / (a_100_2 + b1)) * (b1 / (a_100_2 + b1)) ^
  (b1 / a_100_2)
I_100_2 = absETRmax_100_2 / a_100_2

a_100_3 = round(as.numeric(fit_100a3$par[1]), digits = 3)
b1 = round(as.numeric(fit_100a3$par[2]), digits = 3)
ps1 = round(as.numeric(fit_100a3$par[3]), digits = 3)

absETRmax_100_3 = ps1 * (a_100_3 / (a_100_3 + b1)) * (b1 / (a_100_3 + b1)) ^
  (b1 / a_100_3)
I_100_3 = absETRmax_100_3 / a_100_3

absETR_100_mean <-
  mean(c(absETRmax_100_1, absETRmax_100_2, absETRmax_100_3))
absETR_100_sd <-
  sd(c(absETRmax_100_1, absETRmax_100_2, absETRmax_100_3))

I_100_mean <- mean(c(I_100_1, I_100_2, I_100_3))
I_100_sd <- sd(c(I_100_1, I_100_2, I_100_3))

a_100_mean <- mean(c(a_100_1, a_100_2, a_100_3))
a_100_sd <- sd(c(a_100_1, a_100_2, a_100_3))

# 200 treatment

a_200_1 = round(as.numeric(fit_200a1$par[1]), digits = 3)
b1 = round(as.numeric(fit_200a1$par[2]), digits = 3)
ps1 = round(as.numeric(fit_200a1$par[3]), digits = 3)

absETRmax_200_1 = ps1 * (a_200_1 / (a_200_1 + b1)) * (b1 / (a_200_1 + b1)) ^
  (b1 / a_200_1)
I_200_1 = absETRmax_200_1 / a_200_1

a_200_2 = round(as.numeric(fit_200a2$par[1]), digits = 3)
b1 = round(as.numeric(fit_200a2$par[2]), digits = 3)
ps1 = round(as.numeric(fit_200a2$par[3]), digits = 3)

absETRmax_200_2 = ps1 * (a_200_2 / (a_200_2 + b1)) * (b1 / (a_200_2 + b1)) ^
  (b1 / a_200_2)
I_200_2 = absETRmax_200_2 / a_200_2

a_200_3 = round(as.numeric(fit_200a3$par[1]), digits = 3)
b1 = round(as.numeric(fit_200a3$par[2]), digits = 3)
ps1 = round(as.numeric(fit_200a3$par[3]), digits = 3)

absETRmax_200_3 = ps1 * (a_200_3 / (a_200_3 + b1)) * (b1 / (a_200_3 + b1)) ^
  (b1 / a_200_3)
I_200_3 = absETRmax_200_3 / a_200_3

absETR_200_mean <-
  mean(c(absETRmax_200_1, absETRmax_200_2, absETRmax_200_3))
absETR_200_sd <-
  sd(c(absETRmax_200_1, absETRmax_200_2, absETRmax_200_3))

I_200_mean <- mean(c(I_200_1, I_200_2, I_200_3))
I_200_sd <- sd(c(I_200_1, I_200_2, I_200_3))

a_200_mean <- mean(c(a_200_1, a_200_2, a_200_3))
a_200_sd <- sd(c(a_200_1, a_200_2, a_200_3))

# Combine all for overview
photo_parameters_all <-
  data.frame(
    absETR_mean = rbind(absETR_20_mean, absETR_100_mean, absETR_200_mean),
    absETR_sd = rbind(absETR_20_sd, absETR_100_sd, absETR_200_sd),
    I_mean = rbind(I_20_mean, I_100_mean, I_200_mean),
    I_sd = rbind(I_20_sd, I_100_sd, I_200_sd),
    group = factor(c("20", "100", "200"), levels =
                     c("20", "100", "200"))
  )

rownames(photo_parameters_all) <- c("20", "100", "200")

# Calculate ETR_cell_max by multiplying each ETR_max with the respective RCII concentration (in 10^-18 mol); convert to fmol after; add IK also
# prepare dataframes for export
ETR_cell_2D2 <-
  data.frame(
    ETR_max = c(
      absETRmax_20_1,
      absETRmax_20_2,
      absETRmax_20_3,
      absETRmax_100_1,
      absETRmax_100_2,
      absETRmax_100_3,
      absETRmax_200_1,
      absETRmax_200_2,
      absETRmax_200_3
    ),
    Ik = c(
      I_20_1,
      I_20_2,
      I_20_3,
      I_100_1,
      I_100_2,
      I_100_3,
      I_200_1,
      I_200_2,
      I_200_3
    ),
    a = c(
      a_20_1,
      a_20_2,
      a_20_3,
      a_100_1,
      a_100_2,
      a_100_3,
      a_200_1,
      a_200_2,
      a_200_3
    ),
    group = rep(c("20", "100", "200"), each = 3),
    strain = rep("L2D2")
  )

ETR_cell_2D2$ETR_cell <-
  as.numeric(format(
    ETR_cell_2D2$ETR_max * RCII_2D2,
    digits = 1,
    format = "f"
  ))

ETR_cell_2D2_mean <-
  format(
    aggregate(data = ETR_cell_2D2, ETR_cell ~ group, FUN = "mean"),
    digits = 1,
    format = "f"
  )
ETR_cell_2D2_mean$sd <-
  format(
    aggregate(data = ETR_cell_2D2, ETR_cell ~ group, FUN = "sd"),
    digits = 1,
    format = "f"
  )$'ETR_cell'

ETR_cell_2D2_mean$ETR_max <-
  format(
    aggregate(data = ETR_cell_2D2, ETR_max ~ group, FUN = "mean"),
    digits = 1,
    format = "f"
  )$'ETR_max'
ETR_cell_2D2_mean$ETR_max_sd <-
  format(
    aggregate(data = ETR_cell_2D2, ETR_max ~ group, FUN = "sd"),
    digits = 1,
    format = "f"
  )$'ETR_max'

ETR_cell_2D2_mean$Ik <-
  format(
    aggregate(data = ETR_cell_2D2, Ik ~ group, FUN = "mean"),
    digits = 1,
    format = "f"
  )$'Ik'
ETR_cell_2D2_mean$Ik_sd <-
  format(
    aggregate(data = ETR_cell_2D2, Ik ~ group, FUN = "sd"),
    digits = 1,
    format = "f"
  )$'Ik'

ETR_cell_2D2_mean$a <-
  format(
    aggregate(data = ETR_cell_2D2, a ~ group, FUN = "mean"),
    digits = 3,
    format = "f"
  )$'a'
ETR_cell_2D2_mean$a_sd <-
  format(
    aggregate(data = ETR_cell_2D2, a ~ group, FUN = "sd"),
    digits = 2,
    format = "f"
  )$'a'

ETR_cell_2D2_mean$total_ETR_cell <-
  paste(ETR_cell_2D2_mean$ETR_cell, ETR_cell_2D2_mean$sd, sep = "+/-")
ETR_cell_2D2_mean$total_ETR_max <-
  paste(ETR_cell_2D2_mean$ETR_max,
        ETR_cell_2D2_mean$ETR_max_sd,
        sep = "+/-")
ETR_cell_2D2_mean$total_IK <-
  paste(ETR_cell_2D2_mean$Ik, ETR_cell_2D2_mean$Ik_sd, sep = "+/-")
ETR_cell_2D2_mean$total_a <-
  paste(ETR_cell_2D2_mean$a, ETR_cell_2D2_mean$a_sd, sep = "+/-")

# Statistical analysis of ETR_max and Ik

kruskal_test(data = ETR_cell_2D2, ETR_max ~ group)

conover.test(ETR_cell_2D2$ETR_max,
             ETR_cell_2D2$group,
             method = "BH",
             altp = T)

kruskal_test(data = ETR_cell_2D2, Ik ~ group)

conover.test(ETR_cell_2D2$Ik,
             ETR_cell_2D2$group,
             method = "BH",
             altp = T)

cohens_d(data = ETR_cell_2D2, ETR_max ~ group)

cohens_d(data = ETR_cell_2D2, Ik ~ group)

# Plot all together
PI_model <-
  data.frame(
    PAR = rep(par2),
    ETR = c(
      model_20f(fit_20f$par, par2),
      model_100f(fit_100f$par, par2),
      model_200f(fit_200f$par, par2)
    ),
    group = rep(c("20", "100", "200"), each = 20001)
  )

PI_model$group <- factor(PI_model$group, levels = c(20, 100, 200))

# convert to fmol RCII per cell
PI_2D2_model_Silke$ETR_per_cell <-
  PI_2D2_model_Silke$ETR_per_cell * 10 ^ -3

# Include Silke Thoms photophysiology model
PI_2D2_model_Silke <-
  read.csv(
    "C:\\Users\\krist\\OneDrive\\Dokumente\\AWI\\Promotion\\AP1\\Light-AP1b\\raw\\PI_2D2_model.txt",
    sep = "",
    header = TRUE
  )

# Prepare ggplot parameters
# colorblind colors but skip black (first color) and switch blueish and greenish because NH4 / NO3 are close to each other
colors <- c(colorblind_pal()(4))[c(4, 2, 3)]

limits <- aes(ymax = (ETR + sd), ymin = (ETR - sd)) #Set up the error bars

dodge <- position_jitternudge(width = 0.2, x = 0)

lines <- c('20' = "solid",
           '100' = "longdash",
           '200' = "dotted")

points <- c('20' = 15,
            '100' = 17,
            '200' = 19)

breaks2 <-
  c(bquote(paste("20 ", mu, "mol photons ", m ^ -2, s ^ -1)), bquote(paste("100 ", mu, "mol photons ", m ^
                                                                             -2, s ^ -1)), bquote(paste("200 ", mu, "mol photons ", m ^ -2, s ^ -1)))
breaks <- c("20 PFDs", "100 PFDs", "200 PFDs")

PI_2D2 <- ggplot(data = agg, aes(x = PAR, y = ETR)) +
  geom_point(
    data = subset(agg, agg$group == '20'),
    aes(col = '20'),
    position = dodge,
    show.legend = T
  ) +
  geom_point(
    data = subset(agg, agg$group == '100'),
    aes(col = '100'),
    position = dodge,
    show.legend = T
  ) +
  geom_point(
    data = subset(agg, agg$group == '200'),
    aes(col = '200'),
    position = dodge,
    show.legend = T
  ) +
  geom_line(
    data = subset(PI_2D2_model_Silke, PI_2D2_model_Silke$group == '20'),
    aes(
      x = PAR,
      y = ETR_per_RCII,
      group = '20',
      col = '20'
    ),
    linewidth = 0.5,
    show.legend = T
  ) +
  geom_line(
    data = subset(PI_2D2_model_Silke, PI_2D2_model_Silke$group == '100'),
    aes(
      x = PAR,
      y = ETR_per_RCII,
      group = '100',
      col = '100'
    ),
    linewidth = 0.5,
    show.legend = T
  ) +
  geom_line(
    data = subset(PI_2D2_model_Silke, PI_2D2_model_Silke$group == '200'),
    aes(
      x = PAR,
      y = ETR_per_RCII,
      group = '200',
      col = '200'
    ),
    linewidth = 0.5,
    show.legend = T
  ) +
  theme(
    plot.title = element_blank(),
    strip.background = element_blank(),
    legend.title = element_blank(),
    strip.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    panel.background = element_rect(fill = 'white'),
    legend.key.size = unit(1, "lines"),
    legend.key = element_rect(fill = 'white')
  ) +
  scale_color_manual(
    name = "Photon flux densities",
    values = colors,
    breaks = c(20, 100, 200),
    labels = breaks
  ) +
  xlab("") +
  ylab(expression(paste("ETR (e" ^ '-' * "PSII" ^ -1 * "s" ^ -1 * ")"))) +
  labs(subtitle = "A)")

PI_2D2_2 <- ggplot(data = agg, aes(x = PAR, y = ETR)) +
  geom_pointrange(
    data = subset(agg, agg$group == '20'),
    aes(
      col = '20',
      ymin = ETR - 0,
      ymax = ETR + 0
    ),
    position = dodge,
    show.legend = F,
    size = 0.25
  ) +
  geom_pointrange(
    data = subset(agg, agg$group == '100'),
    aes(
      col = '100',
      ymin = ETR - 0,
      ymax = ETR + 0
    ),
    position = dodge,
    show.legend = F,
    size = 0.25
  ) +
  geom_pointrange(
    data = subset(agg, agg$group == '200'),
    aes(
      col = '200',
      ymin = ETR - 0,
      ymax = ETR + 0
    ),
    position = dodge,
    show.legend = F,
    size = 0.25
  ) +
  geom_line(
    data = subset(PI_2D2_model_Silke, PI_2D2_model_Silke$group == '20'),
    aes(
      x = PAR,
      y = ETR_per_RCII,
      group = '20',
      col = '20'
    ),
    linewidth = 0.5,
    show.legend = F
  ) +
  geom_line(
    data = subset(PI_2D2_model_Silke, PI_2D2_model_Silke$group == '100'),
    aes(
      x = PAR,
      y = ETR_per_RCII,
      group = '100',
      col = '100'
    ),
    linewidth = 0.5,
    show.legend = F
  ) +
  geom_line(
    data = subset(PI_2D2_model_Silke, PI_2D2_model_Silke$group == '200'),
    aes(
      x = PAR,
      y = ETR_per_RCII,
      group = '200',
      col = '200'
    ),
    linewidth = 0.5,
    show.legend = F
  ) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    legend.title = element_blank(),
    panel.background = element_rect(fill = 'white'),
    legend.key.size = unit(1, "lines"),
    legend.key = element_rect(fill = 'white'),
    plot.margin = margin(0, 0.5, 0, 0, "cm")
  ) +
  xlab(expression(paste("PFDs (", mu, "mol photons m" ^ -2 * "s" ^ -1 *
                          ")"))) +
  scale_y_continuous(expression(paste("ETR (e" ^ '-' * "PSII" ^ -1 * "s" ^
                                        -1 * ")")), breaks = seq(0, 1500, 250)) +
  scale_color_manual(
    name = "Photon flux densities",
    values = colors,
    breaks = c(20, 100, 200),
    labels = breaks
  ) +
  labs(subtitle = "A)")

limits <-
  aes(
    ymax = (ETR_per_cell + ETR_per_cell_sd),
    ymin = (ETR_per_cell - ETR_per_cell_sd)
  ) #Set up the error bars

PI_2D2_per_cell <- ggplot(data = agg, aes(x = PAR, y = ETR_per_cell)) +
    geom_pointrange(
    data = subset(agg, agg$group == '20'),
    aes(
      col = '20',
      ymin = ETR_per_cell - 0,
      ymax = ETR_per_cell + 0
    ),
    position = dodge,
    show.legend = T,
    size = 0.25
  ) +
  geom_pointrange(
    data = subset(agg, agg$group == '100'),
    aes(
      col = '100',
      ymin = ETR_per_cell - 0,
      ymax = ETR_per_cell + 0
    ),
    position = dodge,
    show.legend = T,
    size = 0.25
  ) +
  geom_pointrange(
    data = subset(agg, agg$group == '200'),
    aes(
      col = '200',
      ymin = ETR_per_cell - 0,
      ymax = ETR_per_cell + 0
    ),
    position = dodge,
    show.legend = T,
    size = 0.25
  ) +
  geom_line(
    data = subset(PI_2D2_model_Silke, PI_2D2_model_Silke$group == '20'),
    aes(
      x = PAR,
      y = ETR_per_cell,
      group = '20',
      col = '20'
    ),
    linewidth = 0.5,
    show.legend = T
  ) +
  geom_line(
    data = subset(PI_2D2_model_Silke, PI_2D2_model_Silke$group == '100'),
    aes(
      x = PAR,
      y = ETR_per_cell,
      group = '100',
      col = '100'
    ),
    linewidth = 0.5,
    show.legend = T
  ) +
  geom_line(
    data = subset(PI_2D2_model_Silke, PI_2D2_model_Silke$group == '200'),
    aes(
      x = PAR,
      y = ETR_per_cell,
      group = '200',
      col = '200'
    ),
    linewidth = 0.5,
    show.legend = T
  ) +
  theme_classic() +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    panel.background = element_rect(fill = "white"),
    legend.spacing = unit(0, "cm"),
    legend.key.height = unit(1, "lines"),
    plot.margin = margin(0, 0.5, 0, 0, "cm")
  ) +
  scale_color_manual(
    name = "Photon flux densities",
    values = colors,
    breaks = c(20, 100, 200),
    labels = breaks
  ) +
  guides(shape = guide_legend(override.aes = list(
    linetype = lines, linewidth = 0.5
  )),
  linetype = guide_legend(override.aes = list(shape = points, size = 0.75))) +
  xlab(expression(paste("PFDs (", mu, "mol photons m" ^ -2 * "s" ^ -1 *
                          ")"))) +
  ylab(expression(paste(
    "ETR"[cell] * " (e" ^ '-' * "PSII" ^ -1 * "s" ^ -1 * " fmol cell" ^ -1 *
      ")"
  ))) +
  labs(subtitle = "B)")

PI_2D2_per_cell2 <- ggplot(data = agg, aes(x = PAR, y = ETR_per_cell)) +
  geom_pointrange(
    data = subset(agg, agg$group == '20'),
    aes(
      col = '20',
      ymin = ETR_per_cell - 0,
      ymax = ETR_per_cell + 0
    ),
    position = dodge,
    show.legend = T
  ) +
  geom_pointrange(
    data = subset(agg, agg$group == '100'),
    aes(
      col = '100',
      ymin = ETR_per_cell - 0,
      ymax = ETR_per_cell + 0
    ),
    position = dodge,
    show.legend = T
  ) +
  geom_pointrange(
    data = subset(agg, agg$group == '200'),
    aes(
      col = '200',
      ymin = ETR_per_cell - 0,
      ymax = ETR_per_cell + 0
    ),
    position = dodge,
    show.legend = T
  ) +
  geom_line(
    data = subset(PI_2D2_model_Silke, PI_2D2_model_Silke$group == '20'),
    aes(
      x = PAR,
      y = ETR_per_cell,
      group = '20',
      col = '20'
    ),
    linewidth = 0.5,
    show.legend = T
  ) +
  geom_line(
    data = subset(PI_2D2_model_Silke, PI_2D2_model_Silke$group == '100'),
    aes(
      x = PAR,
      y = ETR_per_cell,
      group = '100',
      col = '100'
    ),
    linewidth = 0.5,
    show.legend = T
  ) +
  geom_line(
    data = subset(PI_2D2_model_Silke, PI_2D2_model_Silke$group == '200'),
    aes(
      x = PAR,
      y = ETR_per_cell,
      group = '200',
      col = '200'
    ),
    linewidth = 0.5,
    show.legend = T
  ) +
  theme_classic() +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    panel.background = element_rect(fill = "white"),
    legend.spacing = unit(0, "cm"),
    legend.key.height = unit(1, "lines")
  ) +
  
  scale_color_manual(
    name = "Photon flux densities",
    values = colors,
    breaks = c(20, 100, 200),
    labels = breaks
  ) +
  guides(shape = guide_legend(override.aes = list(
    linetype = lines, linewidth = 0.5
  )),
  linetype = guide_legend(override.aes = list(shape = points, size = 0.75))) +
  ylab(expression(paste(
    "ETR"[cell] * " (e" ^ '-' * "PSII" ^ -1 * "s" ^ -1 * " fmol cell" ^ -1 *
      ")"
  ))) +
  xlab("") +
  labs(subtitle = "A)")

# Plot L2-D2 with Silkes model, both per PSII and per cell
g4 <-
  ggarrange(
    PI_2D2_2,
    PI_2D2_per_cell,
    ncol = 2,
    common.legend = T,
    legend = "bottom"
  )

# export plot
ggsave(
  "PI_all_2D2.png",
  g4,
  path = "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/AP1/Light-AP1b/figures",
  dpi = 300,
  width = 12.5,
  height = 8,
  units = "cm"
)

#### Add A. pseudogonyaulax strain L4-B1
# Load data 
data <-
  read.csv(
    "C:\\Users\\krist\\OneDrive\\Dokumente\\AWI\\Promotion\\AP1\\Light-AP1b\\raw\\PI_4B1.txt",
    sep = "",
    header = TRUE
  )

# replace all , with . + make columns numeric
data[c(2:84), ] <-
  lapply(data[2:84, ], function(x)
    as.numeric(gsub(",", ".", x)))
data <- data[c(1:84), ]

# Include PSII concentrations of each replicate
RCII_4B1 <-
  c(rep(mean(c(
    52.10839161,	23.62183544,	20.00740741
  )), each = 3),	rep(mean(c(
    32.65915119, 28.92719956,	26.12396409
  )), each = 3),	rep(mean(c(
    29.95507246,	24.45754004,	23.63696269
  )), each = 3))

# Make data.frame and name each column; introduce grouping factor for each experimental light intensity 20,100,200
data1 <-
  data.frame(
    PAR = as.numeric(unlist(data[2:84, c(1, 4, 7, 10, 13, 16, 19, 22, 25)])),
    ETR = as.numeric(unlist(data[2:84, c(3, 6, 9, 12, 15, 18, 21, 24, 27)])),
    NPQ = as.numeric(unlist(data[2:84, c(2, 5, 8, 11, 14, 17, 20, 23, 26)])),
    group = as.factor(c(rep(
      c("20", "100", "200"), each = (length(data$X2) - 1) * 3
    ))),
    RCII = as.numeric(c(rep(
      c(RCII_4B1[1:9]), each = (length(data$X2) - 1)
    )))
  )

# Calculate ETR per cell as multiplication of the absolute ETR (calculated in Excel file) with the PSII concentration and conver to fmol ETR per cell
data1$ETR_per_cell <- data1$ETR * data1$RCII * 10 ^ -3

# Calculate mean and SD for each experimental light intensity at each light level of the FRRf (PAR)
# both for ETR and ETR per cell 
data1$PAR <- round(data1$PAR, digits = 0)

agg <- aggregate(ETR ~ PAR + group, data = data1, FUN = "mean")
agg$group <- factor(agg$group, levels = c("20", "100", "200"))
agg$sd <- aggregate(ETR ~ PAR + group, data = data1, FUN = "sd")$ETR

agg$ETR_per_cell <-
  aggregate(ETR_per_cell ~ PAR + group, data = data1, FUN = "mean")$ETR_per_cell
agg$ETR_per_cell_sd <-
  aggregate(ETR_per_cell ~ PAR + group, data = data1, FUN = "sd")$ETR_per_cell

agg2 <- aggregate(NPQ ~ PAR + group, data = data1, FUN = "mean")
agg2$group <- factor(agg2$group, levels = c("20", "100", "200"))
agg2$sd <- aggregate(NPQ ~ PAR + group, data = data1, FUN = "sd")$NPQ

# Subset light intensity 20 for modelling - first model with average points for plot; after each replicate separately for statistics
data1_sub <- subset(data1, data1$group == "20")
PI <- aggregate(data1_sub, ETR ~ PAR, FUN = "mean")

## Model and Modfit  as described in Silsbe and Kromkamp 2012  https://doi.org/10.4319/lom.2012.10.645
## Platt model:
## get starting value for alpha --> slope of linear PAR vs ETR #########
# PI2 <- PI[0:4,]
# ggplot(PI2,aes(x=PAR,y=ETR)) +
#   geom_point() +
#   stat_poly_line( se=F, fullrange = T) +
#   stat_poly_eq(aes(label=paste(..eq.label..,..adj.rr.label..,
#                                sep= "*plain(\",\")~")))

## par = individual light levels of the FRRf; parms = starting values; alpha from slope of linear PAR vs ETR
## ps as maximum observed ETR; beta between 0 and 1
par <- PI$PAR
par2 <- seq(0, 2000, 0.1)
parms_20f = c(alpha = 2.8, beta = 0.35, ps = 800)

# set up model according to Platt; ModelCost returns the residuals; modFit finds local residual minimum
model_20f <-
  function(parms_20f, par)
    with(as.list(parms_20f), return(ps * (1 - exp((
      -1 * alpha * par
    ) / ps)) * exp((-1 * beta * par) / ps)))

ModelCost_20f <- function(P) {
  out <- model_20f(P, par)
  return(PI$ETR - out)  # residuals
}

fit_20f <-
  modFit(
    f = ModelCost_20f,
    p = parms_20f,
    method = "Nelder-Mead",
    lower = c(0, 0, 0),
    upper = c(Inf, 1, Inf)
  )

## plot experimental data points and model as line
# plot(PI$PAR,PI$ETR)
# lines(par2, model_20f(fit_20f$par, par2), lwd = 2, col = "red")
#
# ## model statistics
# summary(fit_20f)
# plot(fit_20f) #residual plot
#
# ## extract model parameters with confidence intervals and residual standard error
# parameters_20_platt <- as.data.frame(summary(fit_20f)$par)
# parameters_20_platt$model <- "platt_20"
# parameters_20_platt$ssr <- fit_20f$ssr

# Subset light intensity 20 for modelling
data1_sub <- subset(data1, data1$group == "20")
data1_sub$grouping <- rep(c(1:3), each = length(data1_sub$PAR) / 3)
PI <- aggregate(data1_sub, ETR ~ PAR + grouping, FUN = "mean")

# PI for Pangea export
PI_pan <-
  aggregate(data1_sub, ETR ~ PAR + grouping, FUN = "mean") %>% mutate(treat = 20, strain = "L4-B1")

# PI subset for each triplicate
PI <- subset(PI, PI$grouping == "1")

## Model and Modfit  as described in Silsbe and Kromkamp 2012  https://doi.org/10.4319/lom.2012.10.645
## Platt model:
## get starting value for alpha --> slope of linear PAR vs ETR #########
# PI2 <- PI[0:4,]
# ggplot(PI2,aes(x=PAR,y=ETR)) +
#   geom_point() +
#   stat_poly_line( se=F, fullrange = T) +
#   stat_poly_eq(aes(label=paste(..eq.label..,..adj.rr.label..,
#                                sep= "*plain(\",\")~")))

## par = individual light levels of the FRRf; parms = starting values; alpha from slope of linear PAR vs ETR
## ps as maximum observed ETR; beta between 0 and 1
par <- PI$PAR
par2 <- seq(0, 2000, 0.1)
parms_20 = c(alpha = 2.8, beta = 0.1, ps = 1200)

# set up model according to Platt; ModelCost returns the residuals; modFit finds local residual minimum
model_20 <-
  function(parms_20, par)
    with(as.list(parms_20), return(ps * (1 - exp((
      -1 * alpha * par
    ) / ps)) * exp((-1 * beta * par) / ps)))

ModelCost_20 <- function(P) {
  out <- model_20(P, par)
  return(PI$ETR - out)  # residuals
}

fit_20a1 <-
  modFit(
    f = ModelCost_20,
    p = parms_20,
    method = "Nelder-Mead",
    lower = c(0, 0, 0),
    upper = c(Inf, 1, Inf)
  )

## plot experimental data points and model as line
# plot(PI$PAR,PI$ETR)
# lines(par2, model_20(fit_20a1$par, par2), lwd = 2, col = "red")
#
# ## model statistics
# summary(fit_20a1)
# plot(fit_20a1) #residual plot
#
# ## extract model parameters with confidence intervals and residual standard error
# parameters_20_platt1 <- as.data.frame(summary(fit_20a1)$par)
# parameters_20_platt1$model <- "platt_20"
# parameters_20_platt1$ssr <- fit_20a1$ssr


# PI subset for each triplicate
PI <- aggregate(data1_sub, ETR ~ PAR + grouping, FUN = "mean")
PI <- subset(PI, PI$grouping == "2")

## Model and Modfit  as described in Silsbe and Kromkamp 2012  https://doi.org/10.4319/lom.2012.10.645
## Platt model:
## get starting value for alpha --> slope of linear PAR vs ETR #########
# PI2 <- PI[0:4,]
# ggplot(PI2,aes(x=PAR,y=ETR)) +
#   geom_point() +
#   stat_poly_line( se=F, fullrange = T) +
#   stat_poly_eq(aes(label=paste(..eq.label..,..adj.rr.label..,
#                                sep= "*plain(\",\")~")))

## par = individual light levels of the FRRf; parms = starting values; alpha from slope of linear PAR vs ETR
## ps as maximum observed ETR; beta between 0 and 1
par <- PI$PAR
par2 <- seq(0, 2000, 0.1)
parms_20 = c(alpha = 2.8, beta = 0.1, ps = 1000)

# set up model according to Platt; ModelCost returns the residuals; modFit finds local residual minimum
model_20 <-
  function(parms_20, par)
    with(as.list(parms_20), return(ps * (1 - exp((
      -1 * alpha * par
    ) / ps)) * exp((-1 * beta * par) / ps)))

ModelCost_20 <- function(P) {
  out <- model_20(P, par)
  return(PI$ETR - out)  # residuals
}

fit_20a2 <-
  modFit(
    f = ModelCost_20,
    p = parms_20,
    method = "Nelder-Mead",
    lower = c(0, 0, 0),
    upper = c(Inf, 1, Inf)
  )

## plot experimental data points and model as line
# plot(PI$PAR,PI$ETR)
# lines(par2, model_20(fit_20a2$par, par2), lwd = 2, col = "red")
#
# ## model statistics
# summary(fit_20a2)
# plot(fit_20a2) #residual plot
#
# ## extract model parameters with confidence intervals and residual standard error
# parameters_20_platt2 <- as.data.frame(summary(fit_20a2)$par)
# parameters_20_platt2$model <- "platt_20"
# parameters_20_platt2$ssr <- fit_20a2$ssr

# PI subset for each triplicate
PI <- aggregate(data1_sub, ETR ~ PAR + grouping, FUN = "mean")
PI <- subset(PI, PI$grouping == "3")

## Model and Modfit  as described in Silsbe and Kromkamp 2012  https://doi.org/10.4319/lom.2012.10.645
## Platt model:
## get starting value for alpha --> slope of linear PAR vs ETR #########
# PI2 <- PI[0:4,]
# ggplot(PI2,aes(x=PAR,y=ETR)) +
#   geom_point() +
#   stat_poly_line( se=F, fullrange = T) +
#   stat_poly_eq(aes(label=paste(..eq.label..,..adj.rr.label..,
#                                sep= "*plain(\",\")~")))

# par = individual light levels of the FRRf; parms = starting values; alpha from slope of linear PAR vs ETR
# ps as maximum observed ETR; beta between 0 and 1
par <- PI$PAR
par2 <- seq(0, 2000, 0.1)
parms_20 = c(alpha = 2.8, beta = 0.1, ps = 1100)

# set up model according to Platt; ModelCost returns the residuals; modFit finds local residual minimum
model_20 <-
  function(parms_20, par)
    with(as.list(parms_20), return(ps * (1 - exp((
      -1 * alpha * par
    ) / ps)) * exp((-1 * beta * par) / ps)))

ModelCost_20 <- function(P) {
  out <- model_20(P, par)
  return(PI$ETR - out)  # residuals
}

fit_20a3 <-
  modFit(
    f = ModelCost_20,
    p = parms_20,
    method = "Nelder-Mead",
    lower = c(0, 0, 0),
    upper = c(Inf, 1, Inf)
  )

## plot experimental data points and model as line
# plot(PI$PAR,PI$ETR)
# lines(par2, model_20(fit_20a3$par, par2), lwd = 2, col = "red")
#
# ##model statistics
# summary(fit_20a3)
# plot(fit_20a3) #residual plot
#
# ## extract model parameters with confidence intervals and residual standard error
# parameters_20_platt3 <- as.data.frame(summary(fit_20a3)$par)
# parameters_20_platt3$model <- "platt_20"
# parameters_20_platt3$ssr <- fit_20a3$ssr

# ## Webb model; starting parameters randomly
# par <- PI$PAR
# parms = c(alpha=3,Ek=470)
#
# ## set up model according to Webb; ModelCost returns the residuals; modFit finds local residual minimum
# model <- function(parms,par) with (as.list(parms), return(alpha*Ek*(1-exp(-1*par/Ek))))
#
# ModelCost <- function(P) {
#   out <- model(P, par)
#   return(PI$ETR-out)  # residuals
# }
#
# fit_20 <- modFit(f = ModelCost, p = parms,method = "Nelder-Mead",lower=c(0,0))
#
#
# ## plot experimental data points and model as line
# plot(PI$PAR,PI$ETR)
# lines(par2, model(fit_20$par, par2), lwd = 2, col = "red")
#
# ## model statistics
# summary(fit_20)
# plot(fit_20)
#
# ## extract model parameters and confidence intervals
# parameters_20_webb <- as.data.frame(summary(fit_20)$par)
# parameters_20_webb$model <- "webb_20"
# parameters_20_webb$ssr <- fit_20$ssr
#
# ## Jasby and Platt model; starting parameters randomly ####
# par <- PI$PAR
# parms = c(alpha=2.3,Ek=590)
#
# ## set up model according to Jasby and Platt; ModelCost returns the residuals; modFit finds local residual minimum
# model <- function(parms,par) with (as.list(parms), return(alpha*Ek*tanh(par/Ek)))
#
# ModelCost <- function(P) {
#   out <- model(P, par)
#   return(PI$ETR-out)  # residuals
# }
#
# fit_20 <- modFit(f = ModelCost, p = parms,method = "Nelder-Mead",lower=c(0,0))
#
# ## plot experimental data points and model as line
# plot(PI$PAR,PI$ETR)
# lines(par2, model(fit_20$par, par2), lwd = 2, col = "red")
#
# ## model statistics
# summary(fit_20)
# plot(fit_20)
#
# ## extract model parameters and confidence intervals
# parameters_20_jasby <- as.data.frame(summary(fit_20)$par)
# parameters_20_jasby$model <- "jasby_20"
# parameters_20_jasby$ssr <- fit_20$ssr

### 100er light treatment
data1_sub <- subset(data1, data1$group == "100")
PI <- aggregate(data1_sub, ETR ~ PAR, FUN = "mean")

## Model and Modfit  as described in Silsbe and Kromkamp 2012  https://doi.org/10.4319/lom.2012.10.645
## Platt model:
## get starting value for alpha --> slope of linear PAR vs ETR #########

# PI2 <- PI[0:4,]
# ggplot(PI2,aes(x=PAR,y=ETR)) +
#   geom_point() +
#   stat_poly_line( se=F, fullrange = T) +
#   stat_poly_eq(aes(label=paste(..eq.label..,..adj.rr.label..,
#                                sep= "*plain(\",\")~")))

## par = individual light levels of the FRRf; parms_100 = starting values; alpha from slope of linear PAR vs ETR
## ps as maximum observed ETR; beta between 0 and 1
par <- PI$PAR
par2 <- seq(0, 2000, 0.1)
parms_100f = c(alpha = 3.8, beta = 0.1, ps = 795)

# set up model_100 according to Platt; ModelCost_100 returns the residuals; modFit finds local residual minimum
model_100f <-
  function(parms_100f, par)
    with(as.list(parms_100f), return(ps * (1 - exp((
      -1 * alpha * par
    ) / ps)) * exp((-1 * beta * par) / ps)))

ModelCost_100f <- function(P) {
  out <- model_100f(P, par)
  return(PI$ETR - out)  # residuals
}

fit_100f <-
  modFit(
    f = ModelCost_100f,
    p = parms_100f,
    method = "Nelder-Mead",
    lower = c(0, 0, 0),
    upper = c(Inf, 1, Inf)
  )

## plot experimental data points and model as line
# plot(PI$PAR,PI$ETR)
# lines(par2, model_100f(fit_100f$par, par2), lwd = 2, col = "red")
#
# ## model statistics
# summary(fit_100f)
# plot(fit_100f) #residual plot
#
# ##extract model parameters with confidence intervals and residual standard error
# parameters_100_platt <- as.data.frame(summary(fit_100f)$par)
# parameters_100_platt$model <- "platt_100"
# parameters_100_platt$ssr <- fit_100f$ssr

# PI subset for each triplicate
data1_sub$grouping <- rep(c(1:3), each = length(data1_sub$group) / 3)
data1_sub$PAR <- round(data1_sub$PAR, digits = 0)
PI <- aggregate(data1_sub, ETR ~ PAR + grouping, FUN = "mean")

# PI for Pangea export
PI_pan2 <-
  aggregate(data1_sub, ETR ~ PAR + grouping, FUN = "mean") %>% mutate(treat = 100, strain = "L4-B1")

PI <- subset(PI, PI$grouping == "1")

## Model and Modfit  as described in Silsbe and Kromkamp 2012  https://doi.org/10.4319/lom.2012.10.645
## Platt model:
## get starting value for alpha --> slope of linear PAR vs ETR #########
# PI2 <- PI[0:4,]
# ggplot(PI2,aes(x=PAR,y=ETR)) +
#   geom_point() +
#   stat_poly_line( se=F, fullrange = T) +
#   stat_poly_eq(aes(label=paste(..eq.label..,..adj.rr.label..,
#                                sep= "*plain(\",\")~")))

## par = individual light levels of the FRRf; parms_100 = starting values; alpha from slope of linear PAR vs ETR
## ps as maximum observed ETR; beta between 0 and 1
par <- PI$PAR
par2 <- seq(0, 2000, 0.1)
parms_100 = c(alpha = 4.9, beta = 0.1, ps = 1500)

# set up model_100 according to Platt; ModelCost_100 returns the residuals; modFit finds local residual minimum
model_100 <-
  function(parms_100, par)
    with(as.list(parms_100), return(ps * (1 - exp((
      -1 * alpha * par
    ) / ps)) * exp((-1 * beta * par) / ps)))

ModelCost_100 <- function(P) {
  out <- model_100(P, par)
  return(PI$ETR - out)  # residuals
}

fit_100a1 <-
  modFit(
    f = ModelCost_100,
    p = parms_100,
    method = "Nelder-Mead",
    lower = c(0, 0, 0),
    upper = c(Inf, 1, Inf)
  )

## plot experimental data points and model as line
# plot(PI$PAR,PI$ETR)
# lines(par2, model_100(fit_100a1$par, par2), lwd = 2, col = "red")
#
# ## model statistics
# summary(fit_100a1)
# plot(fit_100a1) #residual plot
# ##extract model parameters with confidence intervals and residual standard error
# parameters_100_platt1 <- as.data.frame(summary(fit_100a1)$par)
# parameters_100_platt1$model <- "platt_100"
# parameters_100_platt1$ssr <- fit_100a1$ssr

# PI subset for each triplicate
PI <- aggregate(data1_sub, ETR ~ PAR + grouping, FUN = "mean")
PI <- subset(PI, PI$grouping == "2")

## Model and Modfit  as described in Silsbe and Kromkamp 2012  https://doi.org/10.4319/lom.2012.10.645
## Platt model:
## get starting value for alpha --> slope of linear PAR vs ETR #########

# PI2 <- PI[0:4,]
# ggplot(PI2,aes(x=PAR,y=ETR)) +
#   geom_point() +
#   stat_poly_line( se=F, fullrange = T) +
#   stat_poly_eq(aes(label=paste(..eq.label..,..adj.rr.label..,
#                                sep= "*plain(\",\")~")))

# par = individual light levels of the FRRf; parms_100 = starting values; alpha from slope of linear PAR vs ETR
# ps as maximum observed ETR; beta between 0 and 1
par <- PI$PAR
par2 <- seq(0, 2000, 0.1)
parms_100 = c(alpha = 3.8, beta = 0.1, ps = 1000)

# set up model_100 according to Platt; ModelCost_100 returns the residuals; modFit finds local residual minimum
model_100 <-
  function(parms_100, par)
    with(as.list(parms_100), return(ps * (1 - exp((
      -1 * alpha * par
    ) / ps)) * exp((-1 * beta * par) / ps)))

ModelCost_100 <- function(P) {
  out <- model_100(P, par)
  return(PI$ETR - out)  # residuals
}

fit_100a2 <-
  modFit(
    f = ModelCost_100,
    p = parms_100,
    method = "Nelder-Mead",
    lower = c(0, 0, 0),
    upper = c(Inf, 1, Inf)
  )

## plot experimental data points and model as line
# plot(PI$PAR,PI$ETR)
# lines(par2, model_100(fit_100a2$par, par2), lwd = 2, col = "red")
#
# ## model statistics
# summary(fit_100a2)
# plot(fit_100a2) #residual plot
#
# ## extract model parameters with confidence intervals and residual standard error
# parameters_100_platt2 <- as.data.frame(summary(fit_100a2)$par)
# parameters_100_platt2$model <- "platt_100"
# parameters_100_platt2$ssr <- fit_100a2$ssr

# PI subset for each triplicate
PI <- aggregate(data1_sub, ETR ~ PAR + grouping, FUN = "mean")
PI <- subset(PI, PI$grouping == "3")

## Model and Modfit  as described in Silsbe and Kromkamp 2012  https://doi.org/10.4319/lom.2012.10.645
## Platt model:
## get starting value for alpha --> slope of linear PAR vs ETR #########
# PI2 <- PI[0:4,]
# ggplot(PI2,aes(x=PAR,y=ETR)) +
#   geom_point() +
#   stat_poly_line( se=F, fullrange = T) +
#   stat_poly_eq(aes(label=paste(..eq.label..,..adj.rr.label..,
#                                sep= "*plain(\",\")~")))

# par = individual light levels of the FRRf; parms_100 = starting values; alpha from slope of linear PAR vs ETR
# ps as maximum observed ETR; beta between 0 and 1
par <- PI$PAR
par2 <- seq(0, 2000, 0.1)
parms_100 = c(alpha = 3.8, beta = 0.1, ps = 1000)

# set up model_100 according to Platt; ModelCost_100 returns the residuals; modFit finds local residual minimum
model_100 <-
  function(parms_100, par)
    with(as.list(parms_100), return(ps * (1 - exp((
      -1 * alpha * par
    ) / ps)) * exp((-1 * beta * par) / ps)))

ModelCost_100 <- function(P) {
  out <- model_100(P, par)
  return(PI$ETR - out)  # residuals
}

fit_100a3 <-
  modFit(
    f = ModelCost_100,
    p = parms_100,
    method = "Nelder-Mead",
    lower = c(0, 0, 0),
    upper = c(Inf, 1, Inf)
  )

## plot experimental data points and model as line
# plot(PI$PAR,PI$ETR)
# lines(par2, model_100(fit_100a3$par, par2), lwd = 2, col = "red")
#
# ## model statistics
# summary(fit_100a3)
# plot(fit_100a3) #residual plot
#
# ##extract model parameters with confidence intervals and residual standard error
# parameters_100_platt2 <- as.data.frame(summary(fit_100a3)$par)
# parameters_100_platt2$model <- "platt_100"
# parameters_100_platt2$ssr <- fit_100a3$ssr
#
# ## Webb model; starting parameters randomly
# par <- PI$PAR
# parms = c(alpha=3.66,Ek=230)
#
# ## set up model according to Webb; ModelCost returns the residuals; modFit finds local residual minimum
# model <- function(parms,par) with (as.list(parms), return(alpha*Ek*(1-exp(-1*par/Ek))))
#
# ModelCost <- function(P) {
#   out <- model(P, par)
#   return(PI$ETR-out)  # residuals
# }
#
# fit_100 <- modFit(f = ModelCost, p = parms,method = "Nelder-Mead",lower=c(0,0))

## plot experimental data points and model as line
# plot(PI$PAR,PI$ETR)
# lines(par, model(fit_100$par, par), lwd = 2, col = "red")
#
# ## model statistics
# summary(fit_100)
# plot(fit_100)
#
# ## extract model parameters and confidence intervals
# parameters_100_webb <- as.data.frame(summary(fit_100)$par)
# parameters_100_webb$model <- "webb_100"
# parameters_100_webb$ssr <- fit_100$ssr
#
# ## Jasby and Platt model; starting parameters randomly ####
# par <- PI$PAR
# parms = c(alpha=2.77,Ek=299)
#
# ## set up model according to Jasby and Platt; ModelCost returns the residuals; modFit finds local residual minimum
# model <- function(parms,par) with (as.list(parms), return(alpha*Ek*tanh(par/Ek)))
#
# ModelCost <- function(P) {
#   out <- model(P, par)
#   return(PI$ETR-out)  # residuals
# }
#
# fit_100 <- modFit(f = ModelCost, p = parms,method = "Nelder-Mead",lower=c(0,0))
#
# ## plot experimental data points and model as line
# plot(PI$PAR,PI$ETR)
# lines(par, model(fit_100$par, par), lwd = 2, col = "red")
#
# ## model statistics
# summary(fit_100)
# plot(fit_100)
#
# ##extract model parameters and confidence intervals
# parameters_100_jasby <- as.data.frame(summary(fit_100)$par)
# parameters_100_jasby$model <- "jasby_100"
# parameters_100_jasby$ssr <- fit_100$ssr

### 200er light treatment

data1_sub <- subset(data1, data1$group == "200")

PI <- aggregate(data1_sub, ETR ~ PAR, FUN = "mean")

## Model and Modfit  as described in Silsbe and Kromkamp 2012  https://doi.org/10.4319/lom.2012.10.645
## Platt model:
## get starting value for alpha --> slope of linear PAR vs ETR #########
# PI2 <- PI[0:4,]
# ggplot(PI2,aes(x=PAR,y=ETR)) +
#   geom_point() +
#   stat_poly_line( se=F, fullrange = T) +
#   stat_poly_eq(aes(label=paste(..eq.label..,..adj.rr.label..,
#                                sep= "*plain(\",\")~")))

# par = individual light levels of the FRRf; parms_200 = starting values; alpha from slope of linear PAR vs ETR
# ps as maximum observed ETR; beta between 0 and 1
par <- PI$PAR
par2 <- seq(0, 2000, 0.1)
parms_200f = c(alpha = 3.9, beta = 0.1, ps = 1000)

# set up model_200 according to Platt; ModelCost_200 returns the residuals; modFit finds local residual minimum
model_200f <-
  function(parms_200f, par)
    with(as.list(parms_200f), return(ps * (1 - exp((
      -1 * alpha * par
    ) / ps)) * exp((-1 * beta * par) / ps)))

ModelCost_200f <- function(P) {
  out <- model_200f(P, par)
  return(PI$ETR - out)  # residuals
}

fit_200f <-
  modFit(
    f = ModelCost_200f,
    p = parms_200f,
    method = "Nelder-Mead",
    lower = c(0, 0, 0),
    upper = c(Inf, 1, Inf)
  )

## plot experimental data points and model as line
# plot(PI$PAR,PI$ETR)
# lines(par2, model_200f(fit_200f$par, par2), lwd = 2, col = "red")
#
# ## model statistics
# summary(fit_200f)
# plot(fit_200f) #residual plot
#
# ## extract model parameters with confidence intervals and residual standard error
# parameters_200_platt <- as.data.frame(summary(fit_200f)$par)
# parameters_200_platt$model <- "platt_200_1"
# parameters_200_platt$ssr <- fit_200f$ssr

# introduce grouping variable for the triplicate. Modeling has to be done separately to get standard deviations of parameters
data1_sub$grouping <- rep(c(1:3), each = length(data1_sub$group) / 3)
data1_sub$PAR <- round(data1_sub$PAR, digits = 0)
PI <- aggregate(data1_sub, ETR ~ PAR + grouping, FUN = "mean")

# PI subset for each triplicate
PI <- subset(PI, PI$grouping == "1")

## Model and Modfit  as described in Silsbe and Kromkamp 2012  https://doi.org/10.4319/lom.2012.10.645
## Platt model:
## get starting value for alpha --> slope of linear PAR vs ETR #########
# PI2 <- PI[0:4,]
# ggplot(PI2,aes(x=PAR,y=ETR)) +
#   geom_point() +
#   stat_poly_line( se=F, fullrange = T) +
#   stat_poly_eq(aes(label=paste(..eq.label..,..adj.rr.label..,
#                                sep= "*plain(\",\")~")))

# par = individual light levels of the FRRf; parms_200 = starting values; alpha from slope of linear PAR vs ETR
# ps as maximum observed ETR; beta between 0 and 1
par <- PI$PAR
par2 <- seq(0, 2000, 0.1)
parms_200 = c(alpha = 3.5, beta = 0.1, ps = 1000)

# set up model_200 according to Platt; ModelCost_200 returns the residuals; modFit finds local residual minimum
model_200 <-
  function(parms_200, par)
    with(as.list(parms_200), return(ps * (1 - exp((
      -1 * alpha * par
    ) / ps)) * exp((-1 * beta * par) / ps)))

ModelCost_200 <- function(P) {
  out <- model_200(P, par)
  return(PI[1:10, 3] - out)  # residuals
}

fit_200a1 <-
  modFit(
    f = ModelCost_200,
    p = parms_200,
    method = "Nelder-Mead",
    lower = c(0, 0, 0),
    upper = c(Inf, 1, Inf)
  )

## plot experimental data points and model as line
#
# plot(PI$PAR,PI$ETR)
# lines(par2, model_200(fit_200a1$par, par2), lwd = 2, col = "red")
#
# ##model statistics
# summary(fit_200a1)
# plot(fit_200a1) #residual plot
#
# ##extract model parameters with confidence intervals and residual standard error
# parameters_200_platt1 <- as.data.frame(summary(fit_200a1)$par)
# parameters_200_platt1$model <- "platt_200_1"
# parameters_200_platt1$ssr <- fit_200a1$ssr


# PI subset for each triplicate
PI <- aggregate(data1_sub, ETR ~ PAR + grouping, FUN = "mean")
PI <- subset(PI, PI$grouping == "2")

## Model and Modfit  as described in Silsbe and Kromkamp 2012  https://doi.org/10.4319/lom.2012.10.645
## Platt model:
## get starting value for alpha --> slope of linear PAR vs ETR #########
# PI2 <- PI[0:4,]
# ggplot(PI2,aes(x=PAR,y=ETR)) +
#   geom_point() +
#   stat_poly_line( se=F, fullrange = T) +
#   stat_poly_eq(aes(label=paste(..eq.label..,..adj.rr.label..,
#                                sep= "*plain(\",\")~")))

## par = individual light levels of the FRRf; parms_200 = starting values; alpha from slope of linear PAR vs ETR
## ps as maximum observed ETR; beta between 0 and 1
par <- PI$PAR
par2 <- seq(0, 2000, 0.1)
parms_200 = c(alpha = 3.5, beta = 0.1, ps = 1200)

# set up model_200 according to Platt; ModelCost_200 returns the residuals; modFit finds local residual minimum
model_200 <-
  function(parms_200, par)
    with(as.list(parms_200), return(ps * (1 - exp((
      -1 * alpha * par
    ) / ps)) * exp((-1 * beta * par) / ps)))

ModelCost_200 <- function(P) {
  out <- model_200(P, par)
  return(PI[1:10, 3] - out)  # residuals
}

fit_200a2 <-
  modFit(
    f = ModelCost_200,
    p = parms_200,
    method = "Nelder-Mead",
    lower = c(0, 0, 0),
    upper = c(Inf, 1, Inf)
  )

## plot experimental data points and model as line
# plot(PI$PAR,PI$ETR)
# lines(par2, model_200(fit_200a2$par, par2), lwd = 2, col = "red")
#
# ## model statistics
# summary(fit_200a2)
# plot(fit_200a2) #residual plot
#
# ## extract model parameters with confidence intervals and residual standard error
# parameters_200_platt2 <- as.data.frame(summary(fit_200a2)$par)
# parameters_200_platt2$model <- "platt_200_2"
# parameters_200_platt2$ssr <- fit_200a2$ssr
#

## PI subset for each triplicate
PI <- aggregate(data1_sub, ETR ~ PAR + grouping, FUN = "mean")

# PI for Pangea export
PI_pan3 <-
  aggregate(data1_sub, ETR ~ PAR + grouping, FUN = "mean") %>% mutate(treat = 200, strain = "L4-B1")

# Assuming your dataframe is named data2
PI_pan_L4B1 <- rbind(PI_pan, PI_pan2, PI_pan3)

PI <- subset(PI, PI$grouping == "3")

## Model and Modfit  as described in Silsbe and Kromkamp 2012  https://doi.org/10.4319/lom.2012.10.645
## Platt model:
## get starting value for alpha --> slope of linear PAR vs ETR #########
# PI2 <- PI[0:4,]
# ggplot(PI2,aes(x=PAR,y=ETR)) +
#   geom_point() +
#   stat_poly_line( se=F, fullrange = T) +
#   stat_poly_eq(aes(label=paste(..eq.label..,..adj.rr.label..,
#                                sep= "*plain(\",\")~")))

# par = individual light levels of the FRRf; parms_200 = starting values; alpha from slope of linear PAR vs ETR
# ps as maximum observed ETR; beta between 0 and 1
par <- PI$PAR
par2 <- seq(0, 2000, 0.1)
parms_200 = c(alpha = 3.5, beta = 0.1, ps = 1200)

# set up model_200 according to Platt; ModelCost_200 returns the residuals; modFit finds local residual minimum
model_200 <-
  function(parms_200, par)
    with(as.list(parms_200), return(ps * (1 - exp((
      -1 * alpha * par
    ) / ps)) * exp((-1 * beta * par) / ps)))

ModelCost_200 <- function(P) {
  out <- model_200(P, par)
  return(PI[1:10, 3] - out)  # residuals
}

fit_200a3 <-
  modFit(
    f = ModelCost_200,
    p = parms_200,
    method = "Nelder-Mead",
    lower = c(0, 0, 0),
    upper = c(Inf, 1, Inf)
  )

## plot experimental data points and model as line
# plot(PI$PAR,PI$ETR)
# lines(par2, model_200(fit_200a3$par, par2), lwd = 2, col = "red")
#
# ## model statistics
# summary(fit_200a3)
# plot(fit_200a3) #residual plot
#
# ## extract model parameters with confidence intervals and residual standard error
# parameters_200_platt3 <- as.data.frame(summary(fit_200a3)$par)
# parameters_200_platt3$model <- "platt_200_3"
# parameters_200_platt3$ssr <- fit_200a3$ssr
#
# ## Webb model; starting parameters randomly
# par <- PI$PAR
# parms = c(alpha=3.9,Ek=356)
#
# ## set up model according to Webb; ModelCost returns the residuals; modFit finds local residual minimum
# model <- function(parms,par) with (as.list(parms), return(alpha*Ek*(1-exp(-1*par/Ek))))
#
# ModelCost <- function(P) {
#   out <- model(P, par)
#   return(PI$ETR-out)  # residuals
# }
#
# fit_200 <- modFit(f = ModelCost, p = parms,method = "Nelder-Mead",lower=c(0,0))
#
# ## plot experimental data points and model as line
# plot(PI$PAR,PI$ETR)
# lines(par, model(fit_200$par, par), lwd = 2, col = "red")
#
# ## model statistics
# summary(fit_200)
# plot(fit_200)
#
# ## extract model parameters and confidence intervals
# parameters_200_webb <- as.data.frame(summary(fit_200)$par)
# parameters_200_webb$model <- "webb_200"
# parameters_200_webb$ssr <- fit_200$ssr
#
# ## Jasby and Platt model; starting parameters randomly ####
# par <- PI$PAR
# parms = c(alpha=3.05,Ek=410)
#
# ## set up model according to Jasby and Platt; ModelCost returns the residuals; modFit finds local residual minimum
# model <- function(parms,par) with (as.list(parms), return(alpha*Ek*tanh(par/Ek)))
#
# ModelCost <- function(P) {
#   out <- model(P, par)
#   return(PI$ETR-out)  # residuals
# }
#
# fit_200 <- modFit(f = ModelCost, p = parms,method = "Nelder-Mead",lower=c(0,0))
#
# ## plot experimental data points and model as line
# plot(PI$PAR,PI$ETR)
# lines(par, model(fit_200$par, par), lwd = 2, col = "red")
#
# ## model statistics
# summary(fit_200)
# plot(fit_200)
#
# ## extract model parameters and confidence intervals
# parameters_200_jasby <- as.data.frame(summary(fit_200)$par)
# parameters_200_jasby$model <- "jasby_200"
# parameters_200_jasby$ssr <- fit_200$ssr

# calculate maximum electron transport rate and minimum light saturation according to Oxborough et al. 2012
# 20 treatment

a_20_1 = round(as.numeric(fit_20a1$par[1]), digits = 3)
b1 = round(as.numeric(fit_20a1$par[2]), digits = 3)
ps1 = round(as.numeric(fit_20a1$par[3]), digits = 3)

absETRmax_20_1 = ps1 * (a_20_1 / (a_20_1 + b1)) * (b1 / (a_20_1 + b1)) ^
  (b1 / a_20_1)
I_20_1 = absETRmax_20_1 / a_20_1

a_20_2 = round(as.numeric(fit_20a2$par[1]), digits = 3)
b1 = round(as.numeric(fit_20a2$par[2]), digits = 3)
ps1 = round(as.numeric(fit_20a2$par[3]), digits = 3)

absETRmax_20_2 = ps1 * (a_20_2 / (a_20_2 + b1)) * (b1 / (a_20_2 + b1)) ^
  (b1 / a_20_2)
I_20_2 = absETRmax_20_2 / a_20_2

a_20_3 = round(as.numeric(fit_20a3$par[1]), digits = 3)
b1 = round(as.numeric(fit_20a3$par[2]), digits = 3)
ps1 = round(as.numeric(fit_20a3$par[3]), digits = 3)

absETRmax_20_3 = ps1 * (a_20_3 / (a_20_3 + b1)) * (b1 / (a_20_3 + b1)) ^
  (b1 / a_20_3)
I_20_3 = absETRmax_20_3 / a_20_3

absETR_20_mean <-
  mean(c(absETRmax_20_1, absETRmax_20_2, absETRmax_20_3))
absETR_20_sd <- sd(c(absETRmax_20_1, absETRmax_20_2, absETRmax_20_3))

I_20_mean <- mean(c(I_20_1, I_20_2, I_20_3))
I_20_sd <- sd(c(I_20_1, I_20_2, I_20_3))

a_20_mean <- mean(c(a_20_1, a_20_2, a_20_3))
a_20_sd <- sd(c(a_20_1, a_20_2, a_20_3))
## 100 treatment

a_100_1 = round(as.numeric(fit_100a1$par[1]), digits = 3)
b1 = round(as.numeric(fit_100a1$par[2]), digits = 3)
ps1 = round(as.numeric(fit_100a1$par[3]), digits = 3)

absETRmax_100_1 = ps1 * (a_100_1 / (a_100_1 + b1)) * (b1 / (a_100_1 + b1)) ^
  (b1 / a_100_1)
I_100_1 = absETRmax_100_1 / a_100_1

a_100_2 = round(as.numeric(fit_100a2$par[1]), digits = 3)
b1 = round(as.numeric(fit_100a2$par[2]), digits = 3)
ps1 = round(as.numeric(fit_100a2$par[3]), digits = 3)

absETRmax_100_2 = ps1 * (a_100_2 / (a_100_2 + b1)) * (b1 / (a_100_2 + b1)) ^
  (b1 / a_100_2)
I_100_2 = absETRmax_100_2 / a_100_2

a_100_3 = round(as.numeric(fit_100a3$par[1]), digits = 3)
b1 = round(as.numeric(fit_100a3$par[2]), digits = 3)
ps1 = round(as.numeric(fit_100a3$par[3]), digits = 3)

absETRmax_100_3 = ps1 * (a_100_3 / (a_100_3 + b1)) * (b1 / (a_100_3 + b1)) ^
  (b1 / a_100_3)
I_100_3 = absETRmax_100_3 / a_100_3

absETR_100_mean <-
  mean(c(absETRmax_100_1, absETRmax_100_2, absETRmax_100_3))
absETR_100_sd <-
  sd(c(absETRmax_100_1, absETRmax_100_2, absETRmax_100_3))

I_100_mean <- mean(c(I_100_1, I_100_2, I_100_3))
I_100_sd <- sd(c(I_100_1, I_100_2, I_100_3))

a_100_mean <- mean(c(a_100_1, a_100_2, a_100_3))
a_100_sd <- sd(c(a_100_1, a_100_2, a_100_3))

## 200 treatment

a_200_1 = round(as.numeric(fit_200a1$par[1]), digits = 3)
b1 = round(as.numeric(fit_200a1$par[2]), digits = 3)
ps1 = round(as.numeric(fit_200a1$par[3]), digits = 3)

absETRmax_200_1 = ps1 * (a_200_1 / (a_200_1 + b1)) * (b1 / (a_200_1 + b1)) ^
  (b1 / a_200_1)
I_200_1 = absETRmax_200_1 / a_200_1

a_200_2 = round(as.numeric(fit_200a2$par[1]), digits = 3)
b1 = round(as.numeric(fit_200a2$par[2]), digits = 3)
ps1 = round(as.numeric(fit_200a2$par[3]), digits = 3)

absETRmax_200_2 = ps1 * (a_200_2 / (a_200_2 + b1)) * (b1 / (a_200_2 + b1)) ^
  (b1 / a_200_2)
I_200_2 = absETRmax_200_2 / a_200_2

a_200_3 = round(as.numeric(fit_200a3$par[1]), digits = 3)
b1 = round(as.numeric(fit_200a3$par[2]), digits = 3)
ps1 = round(as.numeric(fit_200a3$par[3]), digits = 3)

absETRmax_200_3 = ps1 * (a_200_3 / (a_200_3 + b1)) * (b1 / (a_200_3 + b1)) ^
  (b1 / a_200_3)
I_200_3 = absETRmax_200_3 / a_200_3

absETR_200_mean <-
  mean(c(absETRmax_200_1, absETRmax_200_2, absETRmax_200_3))
absETR_200_sd <-
  sd(c(absETRmax_200_1, absETRmax_200_2, absETRmax_200_3))

I_200_mean <- mean(c(I_200_1, I_200_2, I_200_3))
I_200_sd <- sd(c(I_200_1, I_200_2, I_200_3))

a_200_mean <- mean(c(a_200_1, a_200_2, a_200_3))
a_200_sd <- sd(c(a_200_1, a_200_2, a_200_3))

# Combine all for overview
photo_parameters_all <-
  data.frame(
    absETR_mean = rbind(absETR_20_mean, absETR_100_mean, absETR_200_mean),
    absETR_sd = rbind(absETR_20_sd, absETR_100_sd, absETR_200_sd),
    I_mean = rbind(I_20_mean, I_100_mean, I_200_mean),
    I_sd = rbind(I_20_sd, I_100_sd, I_200_sd),
    group = factor(c("20", "100", "200"), levels =
                     c("20", "100", "200"))
  )

rownames(photo_parameters_all) <- c("20", "100", "200")

# Calculate ETR_cell_max by multiplying each ETR_max with the respective RCII concentration (in 10^-18 mol); convert to fmol after; add IK also
ETR_cell_4B1 <-
  data.frame(
    ETR_max = c(
      absETRmax_20_1,
      absETRmax_20_2,
      absETRmax_20_3,
      absETRmax_100_1,
      absETRmax_100_2,
      absETRmax_100_3,
      absETRmax_200_1,
      absETRmax_200_2,
      absETRmax_200_3
    ),
    Ik = c(
      I_20_1,
      I_20_2,
      I_20_3,
      I_100_1,
      I_100_2,
      I_100_3,
      I_200_1,
      I_200_2,
      I_200_3
    ),
    a = c(
      a_20_1,
      a_20_2,
      a_20_3,
      a_100_1,
      a_100_2,
      a_100_3,
      a_200_1,
      a_200_2,
      a_200_3
    ),
    group = rep(c("20", "100", "200"), each = 3),
    strain = rep("L4B1")
  )

ETR_cell_4B1$ETR_cell <-
  as.numeric(format(
    ETR_cell_4B1$ETR_max * RCII_4B1,
    digits = 1,
    format = "f"
  ))

ETR_cell_4B1_mean <-
  format(
    aggregate(data = ETR_cell_4B1, ETR_cell ~ group, FUN = "mean"),
    digits = 1,
    format = "f"
  )
ETR_cell_4B1_mean$sd <-
  format(
    aggregate(data = ETR_cell_4B1, ETR_cell ~ group, FUN = "sd"),
    digits = 1,
    format = "f"
  )$'ETR_cell'

ETR_cell_4B1_mean$ETR_max <-
  format(
    aggregate(data = ETR_cell_4B1, ETR_max ~ group, FUN = "mean"),
    digits = 1,
    format = "f"
  )$'ETR_max'
ETR_cell_4B1_mean$ETR_max_sd <-
  format(
    aggregate(data = ETR_cell_4B1, ETR_max ~ group, FUN = "sd"),
    digits = 1,
    format = "f"
  )$'ETR_max'

ETR_cell_4B1_mean$Ik <-
  format(
    aggregate(data = ETR_cell_4B1, Ik ~ group, FUN = "mean"),
    digits = 1,
    format = "f"
  )$'Ik'
ETR_cell_4B1_mean$Ik_sd <-
  format(
    aggregate(data = ETR_cell_4B1, Ik ~ group, FUN = "sd"),
    digits = 1,
    format = "f"
  )$'Ik'

ETR_cell_4B1_mean$a <-
  format(
    aggregate(data = ETR_cell_4B1, a ~ group, FUN = "mean"),
    digits = 3,
    format = "f"
  )$'a'
ETR_cell_4B1_mean$a_sd <-
  format(
    aggregate(data = ETR_cell_4B1, a ~ group, FUN = "sd"),
    digits = 2,
    format = "f"
  )$'a'

ETR_cell_4B1_mean$total_ETR_cell <-
  paste(ETR_cell_4B1_mean$ETR_cell, ETR_cell_4B1_mean$sd, sep = "+/-")

ETR_cell_4B1_mean$total_ETR_max <-
  paste(ETR_cell_4B1_mean$ETR_max,
        ETR_cell_4B1_mean$ETR_max_sd,
        sep = "+/-")

ETR_cell_4B1_mean$total_IK <-
  paste(ETR_cell_4B1_mean$Ik, ETR_cell_4B1_mean$Ik_sd, sep = "+/-")

ETR_cell_4B1_mean$total_a <-
  paste(ETR_cell_4B1_mean$a, ETR_cell_4B1_mean$a_sd, sep = "+/-")

# Statistics for ETR_max and Ik
conover.test(ETR_cell_4B1$ETR_max,
             ETR_cell_4B1$group,
             method = "BH",
             altp = T)

conover.test(ETR_cell_4B1$Ik,
             ETR_cell_4B1$group,
             method = "BH",
             altp = T)

cohens_d(data = ETR_cell_4B1, ETR_max ~ group)

cohens_d(data = ETR_cell_4B1, Ik ~ group)

# Plot all together
PI_model <-
  data.frame(
    PAR = rep(par2),
    ETR = c(
      model_20f(fit_20f$par, par2),
      model_100f(fit_100f$par, par2),
      model_200f(fit_200f$par, par2)
    ),
    group = rep(c("20", "100", "200"), each = 20001)
  )
PI_model$group <- factor(PI_model$group, levels = c(20, 100, 200))

PI_model$ETR_per_cell <-
  c(
    subset(PI_model, PI_model$group == 20)$ETR * mean(RCII_4B1[1:3]) * 10 ^
      -3,
    subset(PI_model, PI_model$group == 100)$ETR * mean(RCII_4B1[4:6]) * 10 ^
      -3,
    subset(PI_model, PI_model$group == 200)$ETR * mean(RCII_4B1[7:9]) * 10 ^
      -3
  )

limits <- aes(ymax = (ETR + sd), ymin = (ETR - sd)) #Set up the error bars

dodge <- position_jitternudge(width = 0.2, x = 0)
dodge2 <- position_jitternudge(width = 0.2, x = 25)
dodge3 <- position_jitternudge(width = 0.2, x = 50)

lines <- c('20' = "solid",
           '100' = "longdash",
           '200' = "dotted")

points <- c('20' = 15,
            '100' = 17,
            '200' = 19)

PI_4B1 <- ggplot(data = agg, aes(x = PAR, y = ETR)) +
  geom_point(
    data = subset(agg, agg$group == '20'),
    aes(shape = '20'),
    position = dodge,
    show.legend = T
  ) +
  geom_point(
    data = subset(agg, agg$group == '100'),
    aes(shape = '100'),
    position = dodge2,
    show.legend = T
  ) +
  geom_point(
    data = subset(agg, agg$group == '200'),
    aes(shape = '200'),
    position = dodge3,
    show.legend = T
  ) +
  geom_errorbar(
    data = subset(agg, agg$group == '20'),
    limits,
    position = dodge,
    show.legend = F,
    width = 0
  ) +
  geom_errorbar(
    data = subset(agg, agg$group == '100'),
    limits,
    position = dodge2,
    show.legend = F,
    width = 0
  ) +
  geom_errorbar(
    data = subset(agg, agg$group == '200'),
    limits,
    position = dodge3,
    show.legend = F,
    width = 0
  ) +
  geom_line(
    data = subset(PI_model, PI_model$group == '20'),
    aes(
      x = PAR,
      y = ETR,
      group = '20',
      linetype = '20'
    ),
    show.legend = T
  ) +
  geom_line(
    data = subset(PI_model, PI_model$group == '100'),
    aes(
      x = PAR,
      y = ETR,
      group = '100',
      linetype = '100'
    ),
    show.legend = T
  ) +
  geom_line(
    data = subset(PI_model, PI_model$group == '200'),
    aes(
      x = PAR,
      y = ETR,
      group = '200',
      linetype = '200'
    ),
    show.legend = T
  ) +
  theme(
    plot.title = element_blank(),
    strip.background = element_blank(),
    legend.title = element_blank(),
    strip.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    panel.background = element_rect(fill = 'white'),
    legend.key.size = unit(2, "lines"),
    legend.key = element_rect(fill = 'white')
  ) +
  scale_linetype_manual(
    name = "Photon flux densities",
    values = lines,
    breaks = c(20, 100, 200),
    labels = breaks2
  ) +
  scale_shape_manual(
    name = "Photon flux densities",
    values = points,
    breaks = c(20, 100, 200),
    labels = breaks2
  ) +
  xlab("") +
  ylab(expression(paste("ETR (e" ^ '-' * "PSII" ^ -1 * "s" ^ -1 * ")"))) +
  labs(subtitle = "B)")

limits <-
  aes(
    ymax = (ETR_per_cell + ETR_per_cell_sd),
    ymin = (ETR_per_cell - ETR_per_cell_sd)
  ) #Set up the error bars

PI_4B1_per_cell <- ggplot(data = agg, aes(x = PAR, y = ETR_per_cell)) +
  geom_pointrange(
    data = subset(agg, agg$group == '20'),
    aes(
      col = '20',
      ymin = ETR_per_cell - 0,
      ymax = ETR_per_cell + 0
    ),
    position = dodge,
    show.legend = T
  ) +
  geom_pointrange(
    data = subset(agg, agg$group == '100'),
    aes(
      col = '100',
      ymin = ETR_per_cell - 0,
      ymax = ETR_per_cell + 0
    ),
    position = dodge2,
    show.legend = T
  ) +
  geom_pointrange(
    data = subset(agg, agg$group == '200'),
    aes(
      col = '200',
      ymin = ETR_per_cell - 0,
      ymax = ETR_per_cell + 0
    ),
    position = dodge3,
    show.legend = T
  ) +
  geom_line(
    data = subset(PI_model, PI_model$group == '20'),
    aes(
      x = PAR,
      y = ETR_per_cell,
      group = '20',
      col = '20'
    ),
    linewidth = 0.75,
    show.legend = T
  ) +
  geom_line(
    data = subset(PI_model, PI_model$group == '100'),
    aes(
      x = PAR,
      y = ETR_per_cell,
      group = '100',
      col = '100'
    ),
    linewidth = 0.75,
    show.legend = T
  ) +
  geom_line(
    data = subset(PI_model, PI_model$group == '200'),
    aes(
      x = PAR,
      y = ETR_per_cell,
      group = '200',
      col = '200'
    ),
    linewidth = 0.75,
    show.legend = T
  ) +
  theme_classic() +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    panel.background = element_rect(fill = "white"),
    legend.spacing = unit(0, "cm"),
    legend.key.height = unit(1, "lines"),
    strip.text = element_text(size = 14),
    axis.title = element_text(size = 14)
  ) +
  ggtitle("B)") +
  scale_color_manual(
    name = "Photon flux densities",
    values = colors,
    breaks = c(20, 100, 200),
    labels = breaks2
  ) +
  guides(shape = guide_legend(override.aes = list(
    linetype = lines, linewidth = 0.5
  )),
  linetype = guide_legend(override.aes = list(shape = points, size = 0.75))) +
  # xlab(expression(paste("Photon flux densities [",mu,"mol photons m"^-2*"s"^-1*"]"))) +
  ylab(expression(paste(
    "ETR"[cell] * " (e" ^ '-' * "PSII" ^ -1 * "s" ^ -1 * " fmol cell" ^ -1 *
      ")"
  ))) +
  xlab("")

#### Add A. pseudogonyaulax strain L4-B9
# Load data 
data <-
  read.csv(
    "C:\\Users\\krist\\OneDrive\\Dokumente\\AWI\\Promotion\\AP1\\Light-AP1b\\raw\\PI_4B9.txt",
    sep = "",
    header = TRUE
  )

# replace all , with . + make columns numeric
data[c(2:84), ] <-
  lapply(data[2:84, ], function(x)
    as.numeric(gsub(",", ".", x)))

data <- data[c(1:84), ]

# Include PSII concentrations of each replicate
RCII_4B9 <-
  c(rep(mean(c(
    16.73913043,	21.73279352,	26.32834101
  )), each = 3),	rep(mean(c(
    25.13218391,	19.20020325,	11.31441618
  )), each = 3),	rep(mean(c(
    21.79831288,	19.5063846,	12.33572568
  )), each = 3))

# Make data.frame and name each column; introduce grouping factor for each experimental light intensity 20,100,200
data1 <-
  data.frame(
    PAR = as.numeric(unlist(data[2:84, c(1, 4, 7, 10, 13, 16, 19, 22, 25)])),
    ETR = as.numeric(unlist(data[2:84, c(3, 6, 9, 12, 15, 18, 21, 24, 27)])),
    NPQ = as.numeric(unlist(data[2:84, c(2, 5, 8, 11, 14, 17, 20, 23, 26)])),
    group = as.factor(c(rep(
      c("20", "100", "200"), each = (length(data$X2) - 1) * 3
    ))),
    RCII = as.numeric(c(rep(
      c(RCII_4B9[1:9]), each = (length(data$X2) - 1)
    )))
  )

# Calculate ETR per cell as multiplication of the absolute ETR (calculated in Excel file) with the PSII concentration and conver to fmol ETR per cell
data1$ETR_per_cell <- data1$ETR * data1$RCII * 10 ^ -3

# Calculate mean for each experimental light intensity at each light level of the FRRf (PAR)
data1$PAR <- round(data1$PAR, digits = 0)

agg <- aggregate(ETR ~ PAR + group, data = data1, FUN = "mean")
agg$group <- factor(agg$group, levels = c("20", "100", "200"))
agg$sd <- aggregate(ETR ~ PAR + group, data = data1, FUN = "sd")$ETR

agg$ETR_per_cell <-
  aggregate(ETR_per_cell ~ PAR + group, data = data1, FUN = "mean")$ETR_per_cell
agg$ETR_per_cell_sd <-
  aggregate(ETR_per_cell ~ PAR + group, data = data1, FUN = "sd")$ETR_per_cell

agg2 <- aggregate(NPQ ~ PAR + group, data = data1, FUN = "mean")
agg2$group <- factor(agg2$group, levels = c("20", "100", "200"))
agg2$sd <- aggregate(NPQ ~ PAR + group, data = data1, FUN = "sd")$NPQ

# Subset light intensity 20 for modelling - first model with average points for plot; after each replicate separately for statistics
data1_sub <- subset(data1, data1$group == "20")
PI <- aggregate(data1_sub, ETR ~ PAR, FUN = "mean")

## Model and Modfit  as described in Silsbe and Kromkamp 2012  https://doi.org/10.4319/lom.2012.10.645
## Platt model:
## get starting value for alpha --> slope of linear PAR vs ETR #########
# PI2 <- PI[0:4,]
# ggplot(PI2,aes(x=PAR,y=ETR)) +
#   geom_point() +
#   stat_poly_line( se=F, fullrange = T) +
#   stat_poly_eq(aes(label=paste(..eq.label..,..adj.rr.label..,
#                                sep= "*plain(\",\")~")))

# par = individual light levels of the FRRf; parms = starting values; alpha from slope of linear PAR vs ETR
# ps as maximum observed ETR; beta between 0 and 1
par <- PI$PAR
par2 <- seq(0, 2000, 0.1)
parms_20f = c(alpha = 2.8, beta = 0.35, ps = 2000)

# set up model according to Platt; ModelCost returns the residuals; modFit finds local residual minimum
model_20f <-
  function(parms_20f, par)
    with(as.list(parms_20f), return(ps * (1 - exp((
      -1 * alpha * par
    ) / ps)) * exp((-1 * beta * par) / ps)))

ModelCost_20f <- function(P) {
  out <- model_20f(P, par)
  return(PI$ETR - out)  # residuals
}

fit_20f <-
  modFit(
    f = ModelCost_20f,
    p = parms_20f,
    method = "Nelder-Mead",
    lower = c(0, 0, 0),
    upper = c(Inf, 1, Inf)
  )

## plot experimental data points and model as line
# plot(PI$PAR,PI$ETR)
# lines(par2, model_20f(fit_20f$par, par2), lwd = 2, col = "red")
#
# ## model statistics
# summary(fit_20f)
# plot(fit_20f) #residual plot
#
# ## extract model parameters with confidence intervals and residual standard error
# parameters_20_platt <- as.data.frame(summary(fit_20f)$par)
# parameters_20_platt$model <- "platt_20"
# parameters_20_platt$ssr <- fit_20f$ssr

# Subset light intensity 20 for modelling
data1_sub <- subset(data1, data1$group == "20")
data1_sub$grouping <- rep(c(1:3), each = length(data1_sub$PAR) / 3)
PI <- aggregate(data1_sub, ETR ~ PAR + grouping, FUN = "mean")

# PI for Pangea export
PI_pan <-
  aggregate(data1_sub, ETR ~ PAR + grouping, FUN = "mean") %>% mutate(treat = 20, strain = "L4-B9")

# PI subset for each triplicate
PI <- subset(PI, PI$grouping == "1")

## Model and Modfit  as described in Silsbe and Kromkamp 2012  https://doi.org/10.4319/lom.2012.10.645
## Platt model:
## get starting value for alpha --> slope of linear PAR vs ETR #########
# PI2 <- PI[0:4,]
# ggplot(PI2,aes(x=PAR,y=ETR)) +
#   geom_point() +
#   stat_poly_line( se=F, fullrange = T) +
#   stat_poly_eq(aes(label=paste(..eq.label..,..adj.rr.label..,
#                                sep= "*plain(\",\")~")))

# par = individual light levels of the FRRf; parms = starting values; alpha from slope of linear PAR vs ETR
# ps as maximum observed ETR; beta between 0 and 1
par <- PI$PAR
par2 <- seq(0, 2000, 0.1)
parms_20 = c(alpha = 2.8, beta = 0.1, ps = 1200)

# set up model according to Platt; ModelCost returns the residuals; modFit finds local residual minimum
model_20 <-
  function(parms_20, par)
    with(as.list(parms_20), return(ps * (1 - exp((
      -1 * alpha * par
    ) / ps)) * exp((-1 * beta * par) / ps)))

ModelCost_20 <- function(P) {
  out <- model_20(P, par)
  return(PI$ETR - out)  # residuals
}

fit_20a1 <-
  modFit(
    f = ModelCost_20,
    p = parms_20,
    method = "Nelder-Mead",
    lower = c(0, 0, 0),
    upper = c(Inf, 1, Inf)
  )

## plot experimental data points and model as line
# plot(PI$PAR,PI$ETR)
# lines(par2, model_20(fit_20a1$par, par2), lwd = 2, col = "red")
#
# ## model statistics
# summary(fit_20a1)
# plot(fit_20a1) #residual plot
#
# ## extract model parameters with confidence intervals and residual standard error
# parameters_20_platt1 <- as.data.frame(summary(fit_20a1)$par)
# parameters_20_platt1$model <- "platt_20"
# parameters_20_platt1$ssr <- fit_20a1$ssr

# PI subset for each triplicate
PI <- aggregate(data1_sub, ETR ~ PAR + grouping, FUN = "mean")
PI <- subset(PI, PI$grouping == "2")

## Model and Modfit  as described in Silsbe and Kromkamp 2012  https://doi.org/10.4319/lom.2012.10.645
## Platt model:
## get starting value for alpha --> slope of linear PAR vs ETR #########
# PI2 <- PI[0:4,]
# ggplot(PI2,aes(x=PAR,y=ETR)) +
#   geom_point() +
#   stat_poly_line( se=F, fullrange = T) +
#   stat_poly_eq(aes(label=paste(..eq.label..,..adj.rr.label..,
#                                sep= "*plain(\",\")~")))

# par = individual light levels of the FRRf; parms = starting values; alpha from slope of linear PAR vs ETR
# ps as maximum observed ETR; beta between 0 and 1
par <- PI$PAR
par2 <- seq(0, 2000, 0.1)
parms_20 = c(alpha = 2.8, beta = 0.1, ps = 1100)

# set up model according to Platt; ModelCost returns the residuals; modFit finds local residual minimum
model_20 <-
  function(parms_20, par)
    with(as.list(parms_20), return(ps * (1 - exp((
      -1 * alpha * par
    ) / ps)) * exp((-1 * beta * par) / ps)))

ModelCost_20 <- function(P) {
  out <- model_20(P, par)
  return(PI$ETR - out)  # residuals
}

fit_20a2 <-
  modFit(
    f = ModelCost_20,
    p = parms_20,
    method = "Nelder-Mead",
    lower = c(0, 0, 0),
    upper = c(Inf, 1, Inf)
  )

## plot experimental data points and model as line
# plot(PI$PAR,PI$ETR)
# lines(par2, model_20(fit_20a2$par, par2), lwd = 2, col = "red")
#
# ## model statistics
# summary(fit_20a2)
# plot(fit_20a2) #residual plot
#
# ## extract model parameters with confidence intervals and residual standard error
# parameters_20_platt2 <- as.data.frame(summary(fit_20a2)$par)
# parameters_20_platt2$model <- "platt_20"
# parameters_20_platt2$ssr <- fit_20a2$ssr

# PI subset for each triplicate
PI <- aggregate(data1_sub, ETR ~ PAR + grouping, FUN = "mean")
PI <- subset(PI, PI$grouping == "3")

## Model and Modfit  as described in Silsbe and Kromkamp 2012  https://doi.org/10.4319/lom.2012.10.645
## Platt model:
## get starting value for alpha --> slope of linear PAR vs ETR #########
# PI2 <- PI[0:4,]
# ggplot(PI2,aes(x=PAR,y=ETR)) +
#   geom_point() +
#   stat_poly_line( se=F, fullrange = T) +
#   stat_poly_eq(aes(label=paste(..eq.label..,..adj.rr.label..,
#                                sep= "*plain(\",\")~")))

# par = individual light levels of the FRRf; parms = starting values; alpha from slope of linear PAR vs ETR
# ps as maximum observed ETR; beta between 0 and 1
par <- PI$PAR
par2 <- seq(0, 2000, 0.1)
parms_20 = c(alpha = 2.8, beta = 0.35, ps = 2000)

# set up model according to Platt; ModelCost returns the residuals; modFit finds local residual minimum
model_20 <-
  function(parms_20, par)
    with(as.list(parms_20), return(ps * (1 - exp((
      -1 * alpha * par
    ) / ps)) * exp((-1 * beta * par) / ps)))

ModelCost_20 <- function(P) {
  out <- model_20(P, par)
  return(PI$ETR - out)  # residuals
}

fit_20a3 <-
  modFit(
    f = ModelCost_20,
    p = parms_20,
    method = "Nelder-Mead",
    lower = c(0, 0, 0),
    upper = c(Inf, 1, Inf)
  )

## plot experimental data points and model as line
# plot(PI$PAR,PI$ETR)
# lines(par2, model_20(fit_20a3$par, par2), lwd = 2, col = "red")
#
# ## model statistics
# summary(fit_20a3)
# plot(fit_20a3) #residual plot
#
# ## extract model parameters with confidence intervals and residual standard error
# parameters_20_platt3 <- as.data.frame(summary(fit_20a3)$par)
# parameters_20_platt3$model <- "platt_20"
# parameters_20_platt3$ssr <- fit_20a3$ssr

## Webb model; starting parameters randomly
# par <- PI$PAR
# parms = c(alpha=3,Ek=470)
#
# ## set up model according to Webb; ModelCost returns the residuals; modFit finds local residual minimum
# model <- function(parms,par) with (as.list(parms), return(alpha*Ek*(1-exp(-1*par/Ek))))
#
# ModelCost <- function(P) {
#   out <- model(P, par)
#   return(PI$ETR-out)  # residuals
# }
#
# fit_20 <- modFit(f = ModelCost, p = parms,method = "Nelder-Mead")
#
# ## plot experimental data points and model as line
# plot(PI$PAR,PI$ETR)
# lines(par2, model(fit_20$par, par2), lwd = 2, col = "red")
#
# ## model statistics
# summary(fit_20)
# plot(fit_20)
#
# ## extract model parameters and confidence intervals
# parameters_20_webb <- as.data.frame(summary(fit_20)$par)
# parameters_20_webb$model <- "webb_20"
# parameters_20_webb$ssr <- fit_20$ssr
#
# ## Jasby and Platt model; starting parameters randomly ####
# par <- PI$PAR
# parms = c(alpha=2.3,Ek=590)
#
# ## set up model according to Jasby and Platt; ModelCost returns the residuals; modFit finds local residual minimum
# model <- function(parms,par) with (as.list(parms), return(alpha*Ek*tanh(par/Ek)))
#
# ModelCost <- function(P) {
#   out <- model(P, par)
#   return(PI$ETR-out)  # residuals
# }
#
# fit_20 <- modFit(f = ModelCost, p = parms,method = "Nelder-Mead")
#
# ## plot experimental data points and model as line
# plot(PI$PAR,PI$ETR)
# lines(par2, model(fit_20$par, par2), lwd = 2, col = "red")
#
# ## model statistics
# summary(fit_20)
# plot(fit_20)
#
# ## extract model parameters and confidence intervals
# parameters_20_jasby <- as.data.frame(summary(fit_20)$par)
# parameters_20_jasby$model <- "jasby_20"
# parameters_20_jasby$ssr <- fit_20$ssr

### 100er light treatment
## remove third replicate because the PI measurement stopped around 400 PAR
data1_sub <- subset(data1, data1$group == "100")
PI <- aggregate(na.omit(data1_sub), ETR ~ PAR, FUN = "mean")

## Model and Modfit  as described in Silsbe and Kromkamp 2012  https://doi.org/10.4319/lom.2012.10.645
## Platt model:
## get starting value for alpha --> slope of linear PAR vs ETR #########
# PI2 <- PI[0:4,]
# ggplot(PI2,aes(x=PAR,y=ETR)) +
#   geom_point() +
#   stat_poly_line( se=F, fullrange = T) +
#   stat_poly_eq(aes(label=paste(..eq.label..,..adj.rr.label..,
#                                sep= "*plain(\",\")~")))

# par = individual light levels of the FRRf; parms_100 = starting values; alpha from slope of linear PAR vs ETR
# ps as maximum observed ETR; beta between 0 and 1
par <- PI$PAR
par2 <- seq(0, 2000, 0.1)
parms_100f = c(alpha = 2, beta = 0.1, ps = 1000)

# set up model_100 according to Platt; ModelCost_100 returns the residuals; modFit finds local residual minimum
model_100f <-
  function(parms_100f, par)
    with(as.list(parms_100f), return(ps * (1 - exp((
      -1 * alpha * par
    ) / ps)) * exp((-1 * beta * par) / ps)))

ModelCost_100f <- function(P) {
  out <- model_100f(P, par)
  return(PI$ETR - out)  # residuals
}

fit_100f <-
  modFit(
    f = ModelCost_100f,
    p = parms_100f,
    method = "Nelder-Mead",
    lower = c(0, 0, 0),
    upper = c(Inf, 1, Inf)
  )

## plot experimental data points and model as line
# plot(PI$PAR,PI$ETR)
# lines(par2, model_100f(fit_100f$par, par2), lwd = 2, col = "red")
#
# ## model statistics
# summary(fit_100f)
# plot(fit_100f) #residual plot
#
# ## extract model parameters with confidence intervals and residual standard error
# parameters_100_platt <- as.data.frame(summary(fit_100f)$par)
# parameters_100_platt$model <- "platt_100"
# parameters_100_platt$ssr <- fit_100f$ssr

# PI subset for each triplicate
data1_sub$grouping <- rep(c(1:3), each = length(data1_sub$group) / 3)
data1_sub$PAR <- round(data1_sub$PAR, digits = 0)
PI <- aggregate(data1_sub, ETR ~ PAR + grouping, FUN = "mean")

# PI for Pangea export
PI_pan2 <-
  aggregate(data1_sub, ETR ~ PAR + grouping, FUN = "mean") %>% mutate(treat = 100, strain = "L4-B9")

PI <- subset(PI, PI$grouping == "1")

## Model and Modfit  as described in Silsbe and Kromkamp 2012  https://doi.org/10.4319/lom.2012.10.645
## Platt model:
## get starting value for alpha --> slope of linear PAR vs ETR #########
# PI2 <- PI[0:4,]
# ggplot(PI2,aes(x=PAR,y=ETR)) +
#   geom_point() +
#   stat_poly_line( se=F, fullrange = T) +
#   stat_poly_eq(aes(label=paste(..eq.label..,..adj.rr.label..,
#                                sep= "*plain(\",\")~")))

# par = individual light levels of the FRRf; parms_100 = starting values; alpha from slope of linear PAR vs ETR
# ps as maximum observed ETR; beta between 0 and 1
par <- PI$PAR
par2 <- seq(0, 2000, 0.1)
parms_100 = c(alpha = 2.5, beta = 0.1, ps = 1000)

# set up model_100 according to Platt; ModelCost_100 returns the residuals; modFit finds local residual minimum
model_100 <-
  function(parms_100, par)
    with(as.list(parms_100), return(ps * (1 - exp((
      -1 * alpha * par
    ) / ps)) * exp((-1 * beta * par) / ps)))

ModelCost_100 <- function(P) {
  out <- model_100(P, par)
  return(PI$ETR - out)  # residuals
}

fit_100a1 <-
  modFit(
    f = ModelCost_100,
    p = parms_100,
    method = "Nelder-Mead",
    lower = c(0, 0, 0),
    upper = c(Inf, 1, Inf)
  )

## plot experimental data points and model as line
# plot(PI$PAR,PI$ETR)
# lines(par2, model_100(fit_100a1$par, par2), lwd = 2, col = "red")
#
# ## model statistics
# summary(fit_100a1)
# plot(fit_100a1) #residual plot
#
# ## extract model parameters with confidence intervals and residual standard error
# parameters_100_platt1 <- as.data.frame(summary(fit_100a1)$par)
# parameters_100_platt1$model <- "platt_100"
# parameters_100_platt1$ssr <- fit_100a1$ssr

## PI subset for each triplicate
PI <- aggregate(data1_sub, ETR ~ PAR + grouping, FUN = "mean")
PI <- subset(PI, PI$grouping == "2")

## Model and Modfit  as described in Silsbe and Kromkamp 2012  https://doi.org/10.4319/lom.2012.10.645
## Platt model:
## get starting value for alpha --> slope of linear PAR vs ETR #########
# PI2 <- PI[0:4,]
# ggplot(PI2,aes(x=PAR,y=ETR)) +
#   geom_point() +
#   stat_poly_line( se=F, fullrange = T) +
#   stat_poly_eq(aes(label=paste(..eq.label..,..adj.rr.label..,
#                                sep= "*plain(\",\")~")))

# par = individual light levels of the FRRf; parms_100 = starting values; alpha from slope of linear PAR vs ETR
# ps as maximum observed ETR; beta between 0 and 1
par <- PI$PAR
par2 <- seq(0, 2000, 0.1)
parms_100 = c(alpha = 2.5, beta = 0.1, ps = 1000)

# set up model_100 according to Platt; ModelCost_100 returns the residuals; modFit finds local residual minimum
model_100 <-
  function(parms_100, par)
    with(as.list(parms_100), return(ps * (1 - exp((
      -1 * alpha * par
    ) / ps)) * exp((-1 * beta * par) / ps)))

ModelCost_100 <- function(P) {
  out <- model_100(P, par)
  return(PI$ETR - out)  # residuals
}

fit_100a2 <-
  modFit(
    f = ModelCost_100,
    p = parms_100,
    method = "Nelder-Mead",
    lower = c(0, 0, 0),
    upper = c(Inf, 1, Inf)
  )

# # plot experimental data points and model as line
# plot(PI$PAR,PI$ETR)
# lines(par2, model_100(fit_100a2$par, par2), lwd = 2, col = "red")
#
# ## model statistics
# summary(fit_100a2)
# plot(fit_100a2) #residual plot
#
# ## extract model parameters with confidence intervals and residual standard error
# parameters_100_platt2 <- as.data.frame(summary(fit_100a2)$par)
# parameters_100_platt2$model <- "platt_100"
# parameters_100_platt2$ssr <- fit_100a2$ssr

### 200er light treatment
data1_sub <- subset(data1, data1$group == "200")

PI <- aggregate(data1_sub, ETR ~ PAR, FUN = "mean")

## Model and Modfit  as described in Silsbe and Kromkamp 2012  https://doi.org/10.4319/lom.2012.10.645
## Platt model:
## get starting value for alpha --> slope of linear PAR vs ETR #########
# PI2 <- PI[0:4,]
# ggplot(PI2,aes(x=PAR,y=ETR)) +
#   geom_point() +
#   stat_poly_line( se=F, fullrange = T) +
#   stat_poly_eq(aes(label=paste(..eq.label..,..adj.rr.label..,
#                                sep= "*plain(\",\")~")))

# par = individual light levels of the FRRf; parms_200 = starting values; alpha from slope of linear PAR vs ETR
# ps as maximum observed ETR; beta between 0 and 1
par <- PI$PAR
par2 <- seq(0, 2000, 0.1)
parms_200f = c(alpha = 3.83, beta = 0.29, ps = 1767)

# set up model_200 according to Platt; ModelCost_200 returns the residuals; modFit finds local residual minimum
model_200f <-
  function(parms_200f, par)
    with(as.list(parms_200f), return(ps * (1 - exp((
      -1 * alpha * par
    ) / ps)) * exp((-1 * beta * par) / ps)))

ModelCost_200f <- function(P) {
  out <- model_200f(P, par)
  return(PI$ETR - out)  # residuals
}

fit_200f <-
  modFit(
    f = ModelCost_200f,
    p = parms_200f,
    method = "Nelder-Mead",
    lower = c(0, 0, 0),
    upper = c(Inf, 1, Inf)
  )

## plot experimental data points and model as line
# plot(PI$PAR,PI$ETR)
# lines(par2, model_200f(fit_200f$par, par2), lwd = 2, col = "red")
#
# ## model statistics
# summary(fit_200f)
# plot(fit_200f) #residual plot
#
# ## extract model parameters with confidence intervals and residual standard error
# parameters_200_platt <- as.data.frame(summary(fit_200f)$par)
# parameters_200_platt$model <- "platt_200_1"
# parameters_200_platt$ssr <- fit_200f$ssr

# introduce grouping variable for the triplicate. Modeling has to be done separately to get standard deviations of parameters
data1_sub$grouping <- rep(c(1:3), each = length(data1_sub$group) / 3)
data1_sub$PAR <- round(data1_sub$PAR, digits = 0)
PI <- aggregate(data1_sub, ETR ~ PAR + grouping, FUN = "mean")

# PI for Pangea export
PI_pan3 <-
  aggregate(data1_sub, ETR ~ PAR + grouping, FUN = "mean") %>% mutate(treat = 200, strain = "L4-B9")

PI_pan_L4B9 <- rbind(PI_pan, PI_pan2, PI_pan3)

# PI subset for each triplicate
PI <- subset(PI, PI$grouping == "1")

## Model and Modfit  as described in Silsbe and Kromkamp 2012  https://doi.org/10.4319/lom.2012.10.645
## Platt model:
## get starting value for alpha --> slope of linear PAR vs ETR #########
# PI2 <- PI[0:4,]
# ggplot(PI2,aes(x=PAR,y=ETR)) +
#   geom_point() +
#   stat_poly_line( se=F, fullrange = T) +
#   stat_poly_eq(aes(label=paste(..eq.label..,..adj.rr.label..,
#                                sep= "*plain(\",\")~")))

# par = individual light levels of the FRRf; parms_200 = starting values; alpha from slope of linear PAR vs ETR
# ps as maximum observed ETR; beta between 0 and 1
par <- PI$PAR
par2 <- seq(0, 2000, 0.1)
parms_200 = c(alpha = 4, beta = 0.1, ps = 2000)

# set up model_200 according to Platt; ModelCost_200 returns the residuals; modFit finds local residual minimum
model_200 <-
  function(parms_200, par)
    with(as.list(parms_200), return(ps * (1 - exp((
      -1 * alpha * par
    ) / ps)) * exp((-1 * beta * par) / ps)))

ModelCost_200 <- function(P) {
  out <- model_200(P, par)
  return(PI[1:10, 3] - out)  # residuals
}

fit_200a1 <-
  modFit(
    f = ModelCost_200,
    p = parms_200,
    method = "Nelder-Mead",
    lower = c(0, 0, 0),
    upper = c(Inf, 1, Inf)
  )

## plot experimental data points and model as line
# plot(PI$PAR,PI$ETR)
# lines(par2, model_200(fit_200a1$par, par2), lwd = 2, col = "red")
#
# ## model statistics
# summary(fit_200a1)
# plot(fit_200a1) #residual plot
#
# ## extract model parameters with confidence intervals and residual standard error
# parameters_200_platt1 <- as.data.frame(summary(fit_200a1)$par)
# parameters_200_platt1$model <- "platt_200_1"
# parameters_200_platt1$ssr <- fit_200a1$ssr

# PI subset for each triplicate
PI <- aggregate(data1_sub, ETR ~ PAR + grouping, FUN = "mean")
PI <- subset(PI, PI$grouping == "2")

## Model and Modfit  as described in Silsbe and Kromkamp 2012  https://doi.org/10.4319/lom.2012.10.645
## Platt model:
## get starting value for alpha --> slope of linear PAR vs ETR #########
# PI2 <- PI[0:4,]
# ggplot(PI2,aes(x=PAR,y=ETR)) +
#   geom_point() +
#   stat_poly_line( se=F, fullrange = T) +
#   stat_poly_eq(aes(label=paste(..eq.label..,..adj.rr.label..,
#                                sep= "*plain(\",\")~")))

# par = individual light levels of the FRRf; parms_200 = starting values; alpha from slope of linear PAR vs ETR
# ps as maximum observed ETR; beta between 0 and 1
par <- PI$PAR
par2 <- seq(0, 2000, 0.1)
parms_200 = c(alpha = 4, beta = 0.1, ps = 1200)

# set up model_200 according to Platt; ModelCost_200 returns the residuals; modFit finds local residual minimum
model_200 <-
  function(parms_200, par)
    with(as.list(parms_200), return(ps * (1 - exp((
      -1 * alpha * par
    ) / ps)) * exp((-1 * beta * par) / ps)))

ModelCost_200 <- function(P) {
  out <- model_200(P, par)
  return(PI[1:10, 3] - out)  # residuals
}

fit_200a2 <-
  modFit(
    f = ModelCost_200,
    p = parms_200,
    method = "Nelder-Mead",
    lower = c(0, 0, 0),
    upper = c(Inf, 1, Inf)
  )

## plot experimental data points and model as line
# plot(PI$PAR,PI$ETR)
# lines(par2, model_200(fit_200a2$par, par2), lwd = 2, col = "red")
#
# ## model statistics
# summary(fit_200a2)
# plot(fit_200a2) #residual plot
#
# ## extract model parameters with confidence intervals and residual standard error
# parameters_200_platt2 <- as.data.frame(summary(fit_200a2)$par)
# parameters_200_platt2$model <- "platt_200_2"
# parameters_200_platt2$ssr <- fit_200a2$ssr

# PI subset for each triplicate
PI <- aggregate(data1_sub, ETR ~ PAR + grouping, FUN = "mean")
PI <- subset(PI, PI$grouping == "3")

## Model and Modfit  as described in Silsbe and Kromkamp 2012  https://doi.org/10.4319/lom.2012.10.645
## Platt model:
## get starting value for alpha --> slope of linear PAR vs ETR #########
# PI2 <- PI[0:4,]
# ggplot(PI2,aes(x=PAR,y=ETR)) +
#   geom_point() +
#   stat_poly_line( se=F, fullrange = T) +
#   stat_poly_eq(aes(label=paste(..eq.label..,..adj.rr.label..,
#                                sep= "*plain(\",\")~")))

# par = individual light levels of the FRRf; parms_200 = starting values; alpha from slope of linear PAR vs ETR
# ps as maximum observed ETR; beta between 0 and 1
par <- PI$PAR
par2 <- seq(0, 2000, 0.1)
parms_200 = c(alpha = 4, beta = 0.1, ps = 2000)

# set up model_200 according to Platt; ModelCost_200 returns the residuals; modFit finds local residual minimum
model_200 <-
  function(parms_200, par)
    with(as.list(parms_200), return(ps * (1 - exp((
      -1 * alpha * par
    ) / ps)) * exp((-1 * beta * par) / ps)))

ModelCost_200 <- function(P) {
  out <- model_200(P, par)
  return(PI[1:10, 3] - out)  # residuals
}

fit_200a3 <-
  modFit(
    f = ModelCost_200,
    p = parms_200,
    method = "Nelder-Mead",
    lower = c(0, 0, 0),
    upper = c(Inf, 1, Inf)
  )

## plot experimental data points and model as line
# plot(PI$PAR,PI$ETR)
# lines(par2, model_200(fit_200a3$par, par2), lwd = 2, col = "red")
#
# ## model statistics
# summary(fit_200a3)
# plot(fit_200a3) #residual plot
#
# ## extract model parameters with confidence intervals and residual standard error
# parameters_200_platt3 <- as.data.frame(summary(fit_200a3)$par)
# parameters_200_platt3$model <- "platt_200_3"
# parameters_200_platt3$ssr <- fit_200a3$ssr
#
# ## Webb model; starting parameters randomly
# par <- PI$PAR
# parms = c(alpha=3.9,Ek=356)
#
# ## set up model according to Webb; ModelCost returns the residuals; modFit finds local residual minimum
# model <- function(parms,par) with (as.list(parms), return(alpha*Ek*(1-exp(-1*par/Ek))))
#
# ModelCost <- function(P) {
#   out <- model(P, par)
#   return(PI$ETR-out)  # residuals
# }
#
# fit_200 <- modFit(f = ModelCost, p = parms,method = "Nelder-Mead")
#
# ## plot experimental data points and model as line
# plot(PI$PAR,PI$ETR)
# lines(par, model(fit_200$par, par), lwd = 2, col = "red")
#
# ## model statistics
# summary(fit_200)
# plot(fit_200)
#
# ## extract model parameters and confidence intervals
# parameters_200_webb <- as.data.frame(summary(fit_200)$par)
# parameters_200_webb$model <- "webb_200"
# parameters_200_webb$ssr <- fit_200$ssr
#
# ## Jasby and Platt model; starting parameters randomly ####
# par <- PI$PAR
# parms = c(alpha=3.05,Ek=410)
#
# ## set up model according to Jasby and Platt; ModelCost returns the residuals; modFit finds local residual minimum
# model <- function(parms,par) with (as.list(parms), return(alpha*Ek*tanh(par/Ek)))
#
# ModelCost <- function(P) {
#   out <- model(P, par)
#   return(PI$ETR-out)  # residuals
# }
#
# fit_200 <- modFit(f = ModelCost, p = parms,method = "Nelder-Mead")

## plot experimental data points and model as line
# plot(PI$PAR,PI$ETR)
# lines(par, model(fit_200$par, par), lwd = 2, col = "red")
#
# ## model statistics
# summary(fit_200)
# plot(fit_200)
#
# ## extract model parameters and confidence intervals
# parameters_200_jasby <- as.data.frame(summary(fit_200)$par)
# parameters_200_jasby$model <- "jasby_200"
# parameters_200_jasby$ssr <- fit_200$ssr

## calculate maximum electron transport rate and minimum light saturation according to Oxborough et al. 2012
## 20 treatment
a_20_1 = round(as.numeric(fit_20a1$par[1]), digits = 3)
b1 = round(as.numeric(fit_20a1$par[2]), digits = 3)
ps1 = round(as.numeric(fit_20a1$par[3]), digits = 3)

absETRmax_20_1 = ps1 * (a_20_1 / (a_20_1 + b1)) * (b1 / (a_20_1 + b1)) ^
  (b1 / a_20_1)
I_20_1 = absETRmax_20_1 / a_20_1

a_20_2 = round(as.numeric(fit_20a2$par[1]), digits = 3)
b1 = round(as.numeric(fit_20a2$par[2]), digits = 3)
ps1 = round(as.numeric(fit_20a2$par[3]), digits = 3)

absETRmax_20_2 = ps1 * (a_20_2 / (a_20_2 + b1)) * (b1 / (a_20_2 + b1)) ^
  (b1 / a_20_2)
I_20_2 = absETRmax_20_2 / a_20_2

a_20_3 = round(as.numeric(fit_20a3$par[1]), digits = 3)
b1 = round(as.numeric(fit_20a3$par[2]), digits = 3)
ps1 = round(as.numeric(fit_20a3$par[3]), digits = 3)

absETRmax_20_3 = ps1 * (a_20_3 / (a_20_3 + b1)) * (b1 / (a_20_3 + b1)) ^
  (b1 / a_20_3)
I_20_3 = absETRmax_20_3 / a_20_3

absETR_20_mean <-
  mean(c(absETRmax_20_1, absETRmax_20_2, absETRmax_20_3))
absETR_20_sd <- sd(c(absETRmax_20_1, absETRmax_20_2, absETRmax_20_3))

I_20_mean <- mean(c(I_20_1, I_20_2, I_20_3))
I_20_sd <- sd(c(I_20_1, I_20_2, I_20_3))

a_20_mean <- mean(c(a_20_1, a_20_2, a_20_3))
a_20_sd <- sd(c(a_20_1, a_20_2, a_20_3))

## 100 treatment
a_100_1 = round(as.numeric(fit_100a1$par[1]), digits = 3)
b1 = round(as.numeric(fit_100a1$par[2]), digits = 3)
ps1 = round(as.numeric(fit_100a1$par[3]), digits = 3)

absETRmax_100_1 = ps1 * (a_100_1 / (a_100_1 + b1)) * (b1 / (a_100_1 + b1)) ^
  (b1 / a_100_1)
I_100_1 = absETRmax_100_1 / a_100_1

a_100_2 = round(as.numeric(fit_100a2$par[1]), digits = 3)
b1 = round(as.numeric(fit_100a2$par[2]), digits = 3)
ps1 = round(as.numeric(fit_100a2$par[3]), digits = 3)

absETRmax_100_2 = ps1 * (a_100_2 / (a_100_2 + b1)) * (b1 / (a_100_2 + b1)) ^
  (b1 / a_100_2)
I_100_2 = absETRmax_100_2 / a_100_2

absETR_100_mean <- mean(c(absETRmax_100_1, absETRmax_100_2))
absETR_100_sd <- sd(c(absETRmax_100_1, absETRmax_100_2))

I_100_mean <- mean(c(I_100_1, I_100_2))
I_100_sd <- sd(c(I_100_1, I_100_2))

a_100_mean <- mean(c(a_100_1, a_100_2))
a_100_sd <- sd(c(a_100_1, a_100_2))

## 200 treatment
a_200_1 = round(as.numeric(fit_200a1$par[1]), digits = 3)
b1 = round(as.numeric(fit_200a1$par[2]), digits = 3)
ps1 = round(as.numeric(fit_200a1$par[3]), digits = 3)

absETRmax_200_1 = ps1 * (a_200_1 / (a_200_1 + b1)) * (b1 / (a_200_1 + b1)) ^
  (b1 / a_200_1)
I_200_1 = absETRmax_200_1 / a_200_1

a_200_2 = round(as.numeric(fit_200a2$par[1]), digits = 3)
b1 = round(as.numeric(fit_200a2$par[2]), digits = 3)
ps1 = round(as.numeric(fit_200a2$par[3]), digits = 3)

absETRmax_200_2 = ps1 * (a_200_2 / (a_200_2 + b1)) * (b1 / (a_200_2 + b1)) ^
  (b1 / a_200_2)
I_200_2 = absETRmax_200_2 / a_200_2

a_200_3 = round(as.numeric(fit_200a3$par[1]), digits = 3)
b1 = round(as.numeric(fit_200a3$par[2]), digits = 3)
ps1 = round(as.numeric(fit_200a3$par[3]), digits = 3)

absETRmax_200_3 = ps1 * (a_200_3 / (a_200_3 + b1)) * (b1 / (a_200_3 + b1)) ^
  (b1 / a_200_3)
I_200_3 = absETRmax_200_3 / a_200_3

absETR_200_mean <-
  mean(c(absETRmax_200_1, absETRmax_200_2, absETRmax_200_3))
absETR_200_sd <-
  sd(c(absETRmax_200_1, absETRmax_200_2, absETRmax_200_3))

I_200_mean <- mean(c(I_200_1, I_200_2, I_200_3))
I_200_sd <- sd(c(I_200_1, I_200_2, I_200_3))

a_200_mean <- mean(c(a_200_1, a_200_2, a_200_3))
a_200_sd <- sd(c(a_200_1, a_200_2, a_200_3))

# Combine all for overview
photo_parameters_all <-
  data.frame(
    absETR_mean = rbind(absETR_20_mean, absETR_100_mean, absETR_200_mean),
    absETR_sd = rbind(absETR_20_sd, absETR_100_sd, absETR_200_sd),
    I_mean = rbind(I_20_mean, I_100_mean, I_200_mean),
    I_sd = rbind(I_20_sd, I_100_sd, I_200_sd),
    group = factor(c("20", "100", "200"), levels =
                     c("20", "100", "200"))
  )

rownames(photo_parameters_all) <- c("20", "100", "200")

limits <-
  aes(ymax = (absETR_mean + absETR_sd),
      ymin = (absETR_mean - absETR_sd)) #Set up the error bars

# Calculate ETR_cell_max by multiplying each ETR_max with the respective RCII concentration (in 10^-18 mol); convert to fmol after; add IK also
ETR_cell_4B9 <-
  data.frame(
    ETR_max = c(
      absETRmax_20_1,
      absETRmax_20_2,
      absETRmax_20_3,
      absETRmax_100_1,
      absETRmax_100_2,
      absETRmax_100_3,
      absETRmax_200_1,
      absETRmax_200_2,
      absETRmax_200_3
    ),
    Ik = c(
      I_20_1,
      I_20_2,
      I_20_3,
      I_100_1,
      I_100_2,
      I_100_3,
      I_200_1,
      I_200_2,
      I_200_3
    ),
    a = c(
      a_20_1,
      a_20_2,
      a_20_3,
      a_100_1,
      a_100_2,
      a_100_3,
      a_200_1,
      a_200_2,
      a_200_3
    ),
    group = rep(c("20", "100", "200"), each = 3),
    strain = rep("L4B9")
  )

ETR_cell_4B9$ETR_cell <-
  as.numeric(format(
    ETR_cell_4B9$ETR_max * RCII_4B9,
    digits = 1,
    format = "f"
  ))

ETR_cell_4B9_mean <-
  format(
    aggregate(data = ETR_cell_4B9, ETR_cell ~ group, FUN = "mean"),
    digits = 1,
    format = "f"
  )
ETR_cell_4B9_mean$sd <-
  format(
    aggregate(data = ETR_cell_4B9, ETR_cell ~ group, FUN = "sd"),
    digits = 1,
    format = "f"
  )$'ETR_cell'

ETR_cell_4B9_mean$ETR_max <-
  format(
    aggregate(data = ETR_cell_4B9, ETR_max ~ group, FUN = "mean"),
    digits = 1,
    format = "f"
  )$'ETR_max'

ETR_cell_4B9_mean$ETR_max_sd <-
  format(
    aggregate(data = ETR_cell_4B9, ETR_max ~ group, FUN = "sd"),
    digits = 1,
    format = "f"
  )$'ETR_max'

ETR_cell_4B9_mean$Ik <-
  format(
    aggregate(data = ETR_cell_4B9, Ik ~ group, FUN = "mean"),
    digits = 1,
    format = "f"
  )$'Ik'

ETR_cell_4B9_mean$Ik_sd <-
  format(
    aggregate(data = ETR_cell_4B9, Ik ~ group, FUN = "sd"),
    digits = 1,
    format = "f"
  )$'Ik'

ETR_cell_4B9_mean$a <-
  format(
    aggregate(data = ETR_cell_4B9, a ~ group, FUN = "mean"),
    digits = 3,
    format = "f"
  )$'a'

ETR_cell_4B9_mean$a_sd <-
  format(
    aggregate(data = ETR_cell_4B9, a ~ group, FUN = "sd"),
    digits = 2,
    format = "f"
  )$'a'

ETR_cell_4B9_mean$total_ETR_cell <-
  paste(ETR_cell_4B9_mean$ETR_cell, ETR_cell_4B9_mean$sd, sep = "+/-")

ETR_cell_4B9_mean$total_ETR_max <-
  paste(ETR_cell_4B9_mean$ETR_max,
        ETR_cell_4B9_mean$ETR_max_sd,
        sep = "+/-")

ETR_cell_4B9_mean$total_IK <-
  paste(ETR_cell_4B9_mean$Ik, ETR_cell_4B9_mean$Ik_sd, sep = "+/-")

ETR_cell_4B9_mean$total_a <-
  paste(ETR_cell_4B9_mean$a, ETR_cell_4B9_mean$a_sd, sep = "+/-")

# Statistics for ETR_max and Ik
conover.test(
  ETR_cell_4B9$ETR_max,
  ETR_cell_4B9$group,
  method = "BH",
  altp = T,
  alpha = 0.05
)

conover.test(
  ETR_cell_4B9$Ik,
  ETR_cell_4B9$group,
  method = "BH",
  altp = T,
  alpha = 0.05
)

cohens_d(data = ETR_cell_4B9, ETR_max ~ group)

cohens_d(data = ETR_cell_4B9, Ik ~ group)

# export all PI Pangea data
# rbind(PI_pan_L2D2, PI_pan_L4B1, PI_pan_L4B9) %>% select(treat, PAR, grouping, ETR, strain) %>%
#   write.table("C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/AP1/Light-AP1b/PANGAEA/PI_all.txt", sep = "\t", row.names = FALSE)

# ## combine all and export
# photo_parameters_all <- rbind(ETR_cell_2D2_mean, ETR_cell_4B1_mean, ETR_cell_4B9_mean)
# photo_parameters_all$strain <- rep(c("L2-D2", "L4-B1", "L4-B9"), each = 3)
# photo_parameters_all$group <- factor(photo_parameters_all$group, levels = c(20, 100, 200))
#
# photo_parameters_all <- photo_parameters_all[,c(1,10:14)] %>% arrange(strain, group) %>% flextable() %>% autofit %>%
#     save_as_docx(path="C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/AP1/Light-AP1b/photo_parameters_all.docx")

### Plot all together
PI_model <-
  data.frame(
    PAR = rep(par2),
    ETR = c(
      model_20f(fit_20f$par, par2),
      model_100f(fit_100f$par, par2),
      model_200f(fit_200f$par, par2)
    ),
    group = rep(c("20", "100", "200"), each = 20001)
  )
PI_model$group <- factor(PI_model$group, levels = c(20, 100, 200))

PI_model$ETR_per_cell <-
  c(
    subset(PI_model, PI_model$group == 20)$ETR * mean(RCII_4B9[1:3]) * 10 ^
      -3,
    subset(PI_model, PI_model$group == 100)$ETR * mean(RCII_4B9[4:5]) * 10 ^
      -3,
    subset(PI_model, PI_model$group == 200)$ETR * mean(RCII_4B9[7:9]) * 10 ^
      -3
  )

limits <- aes(ymax = (ETR + sd), ymin = (ETR - sd)) #Set up the error bars

dodge <- position_jitternudge(width = 0.2, x = 0)
dodge2 <- position_jitternudge(width = 0.2, x = 25)
dodge3 <- position_jitternudge(width = 0.2, x = 50)

# Remove sd from 100er treatment because only two replicates are available!
agg[1:10, 4] <- NA

lines <- c('20' = "solid",
           '100' = "longdash",
           '200' = "dotted")

points <- c('20' = 15,
            '100' = 17,
            '200' = 19)

PI_4B9 <- ggplot(data = agg, aes(x = PAR, y = ETR)) +
  geom_point(
    data = subset(agg, agg$group == '20'),
    aes(shape = '20'),
    position = dodge,
    show.legend = T
  ) +
  geom_point(
    data = subset(agg, agg$group == '100'),
    aes(shape = '100'),
    position = dodge2,
    show.legend = T
  ) +
  geom_point(
    data = subset(agg, agg$group == '200'),
    aes(shape = '200'),
    position = dodge3,
    show.legend = T
  ) +
  geom_errorbar(
    data = subset(agg, agg$group == '20'),
    limits,
    position = dodge,
    show.legend = F,
    width = 0
  ) +
  geom_errorbar(
    data = subset(agg, agg$group == '100'),
    limits,
    position = dodge2,
    show.legend = F,
    width = 0
  ) +
  geom_errorbar(
    data = subset(agg, agg$group == '200'),
    limits,
    position = dodge3,
    show.legend = F,
    width = 0
  ) +
  geom_line(
    data = subset(PI_model, PI_model$group == '20'),
    aes(
      x = PAR,
      y = ETR,
      group = '20',
      linetype = '20'
    ),
    show.legend = T
  ) +
  geom_line(
    data = subset(PI_model, PI_model$group == '100'),
    aes(
      x = PAR,
      y = ETR,
      group = '100',
      linetype = '100'
    ),
    show.legend = T
  ) +
  geom_line(
    data = subset(PI_model, PI_model$group == '200'),
    aes(
      x = PAR,
      y = ETR,
      group = '200',
      linetype = '200'
    ),
    show.legend = T
  ) +
  theme(
    plot.title = element_blank(),
    strip.background = element_blank(),
    legend.title = element_blank(),
    strip.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    panel.background = element_rect(fill = 'white'),
    legend.key.size = unit(2, "lines"),
    legend.key = element_rect(fill = 'white')
  ) +
  scale_linetype_manual(
    name = "Photon flux densities",
    values = lines,
    breaks = c(20, 100, 200),
    labels = breaks2
  ) +
  scale_shape_manual(
    name = "Photon flux densities",
    values = points,
    breaks = c(20, 100, 200),
    labels = breaks2
  ) +
  
  xlab(expression(
    paste("Photon flux densities (", mu, "mol photons m" ^ -2 * "s" ^ -1 * ")")
  )) +
  ylab(expression(paste("ETR (e" ^ '-' * "PSII" ^ -1 * "s" ^ -1 * ")"))) +
  labs(subtitle = "C)")


limits <-
  aes(
    ymax = (ETR_per_cell + ETR_per_cell_sd),
    ymin = (ETR_per_cell - ETR_per_cell_sd)
  ) #Set up the error bars

PI_4B9_per_cell <- ggplot(data = agg, aes(x = PAR, y = ETR_per_cell)) +
  geom_pointrange(
    data = subset(agg, agg$group == '20'),
    aes(
      col = '20',
      ymin = ETR_per_cell - 0,
      ymax = ETR_per_cell + 0
    ),
    position = dodge,
    show.legend = T
  ) +
  geom_pointrange(
    data = subset(agg, agg$group == '100'),
    aes(
      col = '100',
      ymin = ETR_per_cell - 0,
      ymax = ETR_per_cell + 0
    ),
    position = dodge2,
    show.legend = T
  ) +
  geom_pointrange(
    data = subset(agg, agg$group == '200'),
    aes(
      col = '200',
      ymin = ETR_per_cell - 0,
      ymax = ETR_per_cell + 0
    ),
    position = dodge3,
    show.legend = T
  ) +
  geom_line(
    data = subset(PI_model, PI_model$group == '20'),
    aes(
      x = PAR,
      y = ETR_per_cell,
      group = '20',
      col = '20'
    ),
    linewidth = 0.75,
    show.legend = T
  ) +
  geom_line(
    data = subset(PI_model, PI_model$group == '100'),
    aes(
      x = PAR,
      y = ETR_per_cell,
      group = '100',
      col = '100'
    ),
    linewidth = 0.75,
    show.legend = T
  ) +
  geom_line(
    data = subset(PI_model, PI_model$group == '200'),
    aes(
      x = PAR,
      y = ETR_per_cell,
      group = '200',
      col = '200'
    ),
    linewidth = 0.75,
    show.legend = T
  ) +
  theme_classic() +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    panel.background = element_rect(fill = "white"),
    legend.spacing = unit(0, "cm"),
    legend.key.height = unit(1, "lines"),
    strip.text = element_text(size = 14),
    axis.title = element_text(size = 14)
  ) +
  ggtitle("C)") +
  scale_color_manual(
    name = "Photon flux densities",
    values = colors,
    breaks = c(20, 100, 200),
    labels = breaks2
  ) +
  guides(shape = guide_legend(override.aes = list(
    linetype = lines, linewidth = 0.5
  )),
  linetype = guide_legend(override.aes = list(shape = points, size = 0.75))) +
  xlab(expression(
    paste("Photon flux densities (", mu, "mol photons m" ^ -2 * "s" ^ -1 * ")")
  )) +
  ylab(expression(paste(
    "ETR"[cell] * " (e" ^ '-' * "PSII" ^ -1 * "s" ^ -1 * " fmol cell" ^ -1 *
      ")"
  )))

# Combine all models in one ggplot
g2 <-
  grid.arrange(rbind(ggplotGrob(PI_2D2), ggplotGrob(PI_4B1), ggplotGrob(PI_4B9)))

g3 <-
  ggarrange(
    PI_2D2_per_cell2,
    PI_4B1_per_cell,
    PI_4B9_per_cell,
    ncol = 1,
    common.legend = T,
    legend = "bottom"
  )

ragg::agg_png(
  "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/AP1/Light-AP1b/figures/PI_all2.png",
  width = 20,
  height = 30,
  units = "in",
  res = 300,
  scaling = 3
)

dev.off()

ggsave(
  "PI_all.png",
  g2,
  path = "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/AP1/Light-AP1b/figures",
  dpi = 300,
  width = 20,
  height = 30,
  units = "cm"
)

ggsave(
  "PI_all_per_cell.png",
  g3,
  path = "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/AP1/Light-AP1b/figures",
  dpi = 300,
  width = 20,
  height = 30,
  units = "cm"
)

# Garbage collection: call after large objects have been removed
gc()

# Delete workspace, clean environment/variables and restarte R if used memory piles up
dev.off()

rm(list = ls())

.rs.restartR()