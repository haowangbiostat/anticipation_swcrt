library(dplyr)
library(ggplot2)
library(gridExtra)
library(lme4)
library(reshape2)

# number of individuals per cluster-period
n_per <- 50

# number of time periods
t_max = 8

# Number of clusters that cross over at each time period
each <- 40

# number of clusters
i_max <- (t_max-1)*each

# maximum number of exposure times observed
expt_max <- t_max - 1

# sd of individual level heterogeneity, or sd of residual error terms
sd_epsilon = 0.5

# treatment effect
overall_eff = 2

# anticipation treatment effect
ant_eff = 0

# cluster-level heterogeneity alpha_i
sigma_alpha = 0.141

# true treatment effect
true_effect_constant <- rep(overall_eff, t_max - 1)
true_effect_lagged <- c(0, 0, rep(2, t_max-3))
true_effect_curved <- c(0.5, 1.0, 1.5, 1.75, 1.875, 1.9375, 2)
true_effect_partially_convex <- c(1/3, 2/3, 1, 1.5, 2, 2, 2)

true_effects <- list(
  'Constant' = true_effect_constant,
  'Lagged' = true_effect_lagged,
  'Curved' = true_effect_curved,
  'Partially convex' = true_effect_partially_convex
)

run_simulation <- function(true_effect) {
  alphas <- rnorm(i_max, 0, sigma_alpha)
  
  data = expand.grid(i = 1:i_max, t = 1:t_max, id = 1:n_per)
  data$sequence <- ceiling(data$i / each)
  data$trt <- ifelse(data$t < (data$sequence+1), 0, 1)
  data$ant <- ifelse(data$t == data$sequence, 1, 0)
  data$expt <- data$t - data$sequence
  data[data$expt < 0,]$expt = 0
  data$q <- ceiling(data$i / 2)
  
  t_fixef_osc <- 0.5 * sin(pi * 2 * (1:t_max - 1) / (t_max - 1))
  
  data$Y <- 14 +
    apply(data, 1, function(x) {
      alphas[as.numeric(x["i"])] + t_fixef_osc[as.numeric(x["t"])]
    })
  
  expt_eff = c(0, true_effect)
  
  data$Y <- data$Y + 
    apply(data, 1, function(x) {
      expt_eff[as.numeric(x["expt"]) + 1]
    }) + rnorm(nrow(data), 0, sd_epsilon)
  
  data$Y <- data$Y + data$ant * ant_eff
  
  return(data)
}

n_sim <- 1000
df_all <- data.frame()

set.seed(123)
for (pattern in names(true_effects)) {
  true_effect <- true_effects[[pattern]]
  trt_coefficients_model_1 <- numeric(n_sim)
  trt_coefficients_model_2 <- numeric(n_sim)
  ant_coefficients_model_2 <- numeric(n_sim)
  
  for (sim in 1:n_sim) {
    data <- run_simulation(true_effect)
    model_1 <- lmer(Y ~ trt + as.factor(t) + (1 | i), data = data)
    model_2 <- lmer(Y ~ trt + ant + as.factor(t) + (1 | i), data = data)
    trt_coefficients_model_1[sim] <- fixef(model_1)["trt"]
    trt_coefficients_model_2[sim] <- fixef(model_2)["trt"]
    ant_coefficients_model_2[sim] <- fixef(model_2)["ant"]
  }
  
  # average treatment effect over n_sim simulations
  mean_trt_coeff_model_1 <- mean(trt_coefficients_model_1)
  mean_trt_coeff_model_2 <- mean(trt_coefficients_model_2)
  mean_ant_coeff_model_2 <- mean(ant_coefficients_model_2)
  
  # expectation of the treatment effect under the HH-ANT model
  phi <- sigma_alpha*sigma_alpha/(sigma_alpha*sigma_alpha + sd_epsilon*sd_epsilon/n_per)
  Q <- t_max - 1
  theoretical_value <- 0
  for (j in 1:Q) {
    numerator <- 6 * (phi * Q^3 - 3 * phi * Q^2 * j + phi * Q^2 + Q^2 + 
                        2 * phi * Q * j^2 - 2 * phi * Q * j + phi * Q - 2 * Q * j + j^2)
    denominator <- Q * (phi * Q^3 - 3 * phi * Q^2 + 2 * Q^2 + 2 * phi * Q - 3 * Q + 1)
    
    theoretical_value <- theoretical_value + (numerator / denominator)* true_effect[j]
  }
  
  theoretical_value_ant <- ant_eff
  for (j in 1:Q) {
    numerator <- 2 * phi * Q^3 - 8 * phi * Q^2 * j + 5 * phi * Q^2 + Q^2 +
      6 * phi * Q * j^2 - 8 * phi * Q * j + 3 * phi * Q -
      4 * Q * j + Q + 3 * j^2 - j
    denominator <- Q * (Q - 1) * (phi * Q^2 - 2 * phi * Q + 2 * Q - 1)
    
    theoretical_value_ant <- theoretical_value_ant + (numerator / denominator) * true_effect[j]
  }
  
  time_points <- 0:(t_max - 1)
  df <- data.frame(
    Time = time_points,
    True = c(ant_eff, true_effect),
    Estimated_model_1 = c(mean_trt_coeff_model_1, rep(mean_trt_coeff_model_1, t_max - 1)),
    Estimated_model_2 = c(mean_ant_coeff_model_2, rep(mean_trt_coeff_model_2, t_max - 1)),
    Theoretical = c(theoretical_value_ant, rep(theoretical_value, t_max - 1)),
    Pattern = pattern
  )
  
  df_all <- rbind(df_all, df)
}

df_long <- melt(df_all, id.vars = c("Time", "Pattern"), 
                variable.name = "EffectType", value.name = "Effect")
df_long$EffectType <- factor(df_long$EffectType, 
                             levels = c("Estimated_model_1", "Estimated_model_2", "Theoretical", "True"), 
                             labels = c("Estimated (HH)", "Estimated (HH-ANT)", "Theoretical (HH-ANT)", "True"))

# true treatment effect
p1 <- ggplot(df_long, aes(x = Time, y = Effect, color = EffectType)) +
  geom_line(data = subset(df_long, EffectType == "Estimated (HH)" & Time >= 1), 
            linetype = "longdash") +
  geom_line(data = subset(df_long, EffectType == "Estimated (HH-ANT)"), 
            linetype = "dashed") +
  geom_line(data = subset(df_long, EffectType == "True" & Pattern == "Constant"), 
            alpha = 0.4) +
  geom_line(data = subset(df_long, EffectType == "True" & Pattern != "Constant")) +
  geom_point(data = subset(df_long, EffectType == "True")) +
  scale_y_continuous(breaks = c(-1, 0, 1, 2), limits = c(-1, 2.02)) +
  scale_x_continuous(breaks = seq(0, t_max - 1, 1), 
                     labels = function(x) ifelse(x == 0, "A", x), position = "top") +
  scale_color_manual(values = c("Estimated (HH)" = "#f4c100", 
                                "Estimated (HH-ANT)" = "#3fa535", 
                                "True" = "#009cc4")) +
  labs(x = "Exposure Time",
       y = "Trt\nEffect\n\nAnt\nEffect",
       color = "Effect curve") +
  facet_wrap(~Pattern, ncol = 4, strip.position = "bottom",
             labeller = labeller(Pattern = c("Constant" = "(a) Constant",
                                             "Curved" = "(b) Curved",
                                             "Lagged" = "(c) Lagged",
                                             "Partially convex" = "(d) Partially convex"))) +
  theme(legend.position = "bottom", 
        legend.direction = "horizontal",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12, angle = 0, vjust = 0.5),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x.top = element_text(size = 12),
        axis.text.x.top = element_text(size = 12),
        strip.text = element_text(size = 12),
        strip.placement = "outside")

# number of individuals per cluster-period
n_per <- 50

# number of time periods
t_max = 8

# number of clusters per sequence
each <- 40

# number of clusters
i_max <- (t_max-1)*each

# maximum number of exposure times
expt_max <- t_max - 1

# sd of individual level heterogeneity
sd_epsilon = 0.5

# overall treatment effect
overall_eff = 2

# anticipation treatment effect
ant_eff = 0.5

# cluster-level heterogeneity alpha_i
sigma_alpha = 0.141

# true effect
true_effect_constant <- rep(overall_eff, t_max - 1)
true_effect_lagged <- c(0.5, 0.5, rep(2, t_max-3))
true_effect_curved <- c(0.5, 1.0, 1.5, 1.75, 1.875, 1.9375, 2)
true_effect_partially_convex <- c(0.5, 0.75, 1, 1.5, 2, 2, 2)

true_effects <- list(
  'Constant' = true_effect_constant,
  'Lagged' = true_effect_lagged,
  'Curved' = true_effect_curved,
  'Partially convex' = true_effect_partially_convex
)

run_simulation <- function(true_effect) {
  alphas <- rnorm(i_max, 0, sigma_alpha)
  
  data = expand.grid(i = 1:i_max, t = 1:t_max, id = 1:n_per)
  data$sequence <- ceiling(data$i / each)
  data$trt <- ifelse(data$t < (data$sequence+1), 0, 1)
  data$ant <- ifelse(data$t == data$sequence, 1, 0)
  data$expt <- data$t - data$sequence
  data[data$expt < 0,]$expt = 0
  data$q <- ceiling(data$i / 2)
  
  t_fixef_osc <- 0.5 * sin(pi * 2 * (1:t_max - 1) / (t_max - 1))
  
  data$Y <- 14 +
    apply(data, 1, function(x) {
      alphas[as.numeric(x["i"])] + t_fixef_osc[as.numeric(x["t"])]
    })
  
  expt_eff = c(0, true_effect)
  
  data$Y <- data$Y + 
    apply(data, 1, function(x) {
      expt_eff[as.numeric(x["expt"]) + 1]
    }) + rnorm(nrow(data), 0, sd_epsilon)
  
  data$Y <- data$Y + data$ant * ant_eff
  
  return(data)
}

n_sim <- 1000
df_all <- data.frame()

set.seed(123)
for (pattern in names(true_effects)) {
  true_effect <- true_effects[[pattern]]
  trt_coefficients_model_1 <- numeric(n_sim)
  trt_coefficients_model_2 <- numeric(n_sim)
  ant_coefficients_model_2 <- numeric(n_sim)
  
  for (sim in 1:n_sim) {
    data <- run_simulation(true_effect)
    model_1 <- lmer(Y ~ trt + as.factor(t) + (1 | i), data = data)
    model_2 <- lmer(Y ~ trt + ant + as.factor(t) + (1 | i), data = data)
    trt_coefficients_model_1[sim] <- fixef(model_1)["trt"]
    trt_coefficients_model_2[sim] <- fixef(model_2)["trt"]
    ant_coefficients_model_2[sim] <- fixef(model_2)["ant"]
  }
  
  # average treatment effect over n_sim simulations
  mean_trt_coeff_model_1 <- mean(trt_coefficients_model_1)
  mean_trt_coeff_model_2 <- mean(trt_coefficients_model_2)
  mean_ant_coeff_model_2 <- mean(ant_coefficients_model_2)
  
  # expectation of the treatment effect under the HH-ANT model
  phi <- sigma_alpha*sigma_alpha/(sigma_alpha*sigma_alpha + sd_epsilon*sd_epsilon/n_per)
  Q <- t_max - 1
  theoretical_value <- 0
  for (j in 1:Q) {
    numerator <- 6 * (phi * Q^3 - 3 * phi * Q^2 * j + phi * Q^2 + Q^2 + 
                        2 * phi * Q * j^2 - 2 * phi * Q * j + phi * Q - 2 * Q * j + j^2)
    denominator <- Q * (phi * Q^3 - 3 * phi * Q^2 + 2 * Q^2 + 2 * phi * Q - 3 * Q + 1)
    
    theoretical_value <- theoretical_value + (numerator / denominator)* true_effect[j]
  }
  
  theoretical_value_ant <- ant_eff
  for (j in 1:Q) {
    numerator <- 2 * phi * Q^3 - 8 * phi * Q^2 * j + 5 * phi * Q^2 + Q^2 +
      6 * phi * Q * j^2 - 8 * phi * Q * j + 3 * phi * Q -
      4 * Q * j + Q + 3 * j^2 - j
    denominator <- Q * (Q - 1) * (phi * Q^2 - 2 * phi * Q + 2 * Q - 1)
    
    theoretical_value_ant <- theoretical_value_ant + (numerator / denominator) * true_effect[j]
  }
  
  time_points <- 0:(t_max - 1)
  df <- data.frame(
    Time = time_points,
    True = c(ant_eff, true_effect),
    Estimated_model_1 = c(mean_trt_coeff_model_1, rep(mean_trt_coeff_model_1, t_max - 1)),
    Estimated_model_2 = c(mean_ant_coeff_model_2, rep(mean_trt_coeff_model_2, t_max - 1)),
    Theoretical = c(theoretical_value_ant, rep(theoretical_value, t_max - 1)),
    Pattern = pattern
  )
  
  df_all <- rbind(df_all, df)
}

df_long <- melt(df_all, id.vars = c("Time", "Pattern"), variable.name = "EffectType", 
                value.name = "Effect")
df_long$EffectType <- factor(df_long$EffectType, 
                             levels = c("Estimated_model_1", "Estimated_model_2", "Theoretical", "True"), 
                             labels = c("Estimated (HH)", "Estimated (HH-ANT)", "Theoretical (HH-ANT)", "True"))

p2 <- ggplot(df_long, aes(x = Time, y = Effect, color = EffectType)) +
  geom_line(data = subset(df_long, EffectType == "Estimated (HH)" & Time >= 1), 
            linetype = "longdash") +
  geom_line(data = subset(df_long, EffectType == "Estimated (HH-ANT)"), 
            linetype = "dashed") +
  geom_line(data = subset(df_long, EffectType == "True" & Pattern == "Constant"), 
            alpha = 0.4) +
  geom_line(data = subset(df_long, EffectType == "True" & Pattern != "Constant")) +
  geom_point(data = subset(df_long, EffectType == "True")) +
  scale_y_continuous(breaks = c(-1, 0, 1, 2), limits = c(-1, 2.02)) +
  scale_x_continuous(breaks = seq(0, t_max - 1, 1), 
                     labels = function(x) ifelse(x == 0, "A", x), position = "top") +
  scale_color_manual(values = c("Estimated (HH)" = "#f4c100", 
                                "Estimated (HH-ANT)" = "#3fa535",
                                "True" = "#009cc4")) +
  labs(x = "Exposure Time",
       y = "Trt\nEffect\n\nAnt\nEffect",
       color = "Effect curve") +
  facet_wrap(~Pattern, ncol = 4, strip.position = "bottom",
             labeller = labeller(Pattern = c("Constant" = "(e) Constant",
                                             "Curved" = "(f) Curved",
                                             "Lagged" = "(g) Lagged",
                                             "Partially convex" = "(h) Partially convex"))) +
  theme(legend.position = "bottom", 
        legend.direction = "horizontal",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12, angle = 0, vjust = 0.5),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x.top = element_text(size = 12),
        axis.text.x.top = element_text(size = 12),
        strip.text = element_text(size = 12),
        strip.placement = "outside")

p3 <- grid.arrange(p1, p2, nrow = 2)

ggsave(filename = "../figures/figure_HH-ANT.pdf", plot = p3, width = 12, height = 5)