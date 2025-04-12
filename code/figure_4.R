library(dplyr)
library(ggplot2)
library(gridExtra)
library(lme4)
library(reshape2)

# number of individuals per cluster-period
n_per <- 50

# number of time periods
t_max = 8

# maximum number of exposure times
expt_max <- t_max - 1

# sd of individual level heterogeneity
sd_epsilon = 0.5

# anticipation treatment effect
ant_eff = 0.5

# cluster-level heterogeneity
sigma_alpha = 0.141

# treatment effect
overall_eff = 2

# simulation
run_simulation <- function(true_effect, t_max, i_max) {
  alphas <- rnorm(i_max, 0, sigma_alpha)
  
  data <- expand.grid(i = 1:i_max, t = 1:t_max, id = 1:n_per)
  data$sequence <- ceiling(data$i / each)
  data$trt <- ifelse(data$t < (data$sequence + 1), 0, 1)
  data$ant <- ifelse(data$t == data$sequence, 1, 0)
  data$expt <- data$t - data$sequence
  data[data$expt < 0, ]$expt <- 0
  data$q <- ceiling(data$i / 2)
  
  t_fixef_osc <- 0.5 * sin(pi * 2 * (1:t_max - 1) / (t_max - 1))
  
  data$Y <- 14 +
    apply(data, 1, function(x) {
      alphas[as.numeric(x["i"])] + t_fixef_osc[as.numeric(x["t"])]
    })
  
  expt_eff <- c(0, true_effect)
  
  data$Y <- data$Y + 
    apply(data, 1, function(x) {
      expt_eff[as.numeric(x["expt"]) + 1]
    }) + rnorm(nrow(data), 0, sd_epsilon)
  
  data$Y <- data$Y + data$ant * ant_eff
  
  return(data)
}

# true effect
true_effect_constant <- rep(overall_eff, 2)
true_effect_lagged <- c(0.5, 2)
true_effect_curved <- c(1.75, 2)
true_effect_partially_convex <- c(1, 2)

true_effects_J3 <- list(
  'Constant' = true_effect_constant,
  'Lagged' = true_effect_lagged,
  'Curved' = true_effect_curved,
  'Partially convex' = true_effect_partially_convex
)

# true treatment effect
true_effect_constant <- rep(overall_eff, 7)
true_effect_lagged <- c(0.5, 0.5, rep(2, 5))
true_effect_curved <- c(0.5, 1.0, 1.5, 1.75, 1.875, 1.9375, 2)
true_effect_partially_convex <- c(0.5, 0.75, 1, 1.5, 2, 2, 2)

true_effects_J8 <- list(
  'Constant' = true_effect_constant,
  'Lagged' = true_effect_lagged,
  'Curved' = true_effect_curved,
  'Partially convex' = true_effect_partially_convex
)

n_sim <- 1000
df_all <- data.frame()

set.seed(123)
for (t_max in c(3, 8)) {
  if (t_max == 3) {
    true_effects = true_effects_J3
    # number of clusters per sequence
    each = 140
  }
  if (t_max == 8) {
    true_effects = true_effects_J8
    # number of clusters per sequence
    each = 40
  }
  i_max <- (t_max - 1) * each
  for (pattern in names(true_effects)) {
    true_effect <- true_effects[[pattern]]
    expt_coefficients <- matrix(NA, nrow = n_sim, ncol = t_max - 1)
    
    for (sim in 1:n_sim) {
      data <- run_simulation(true_effect, t_max, i_max)
      model_1 <- lmer(Y ~ as.factor(expt) + as.factor(t) + (1 | i), data = data)
      
      for (expt in 1:(t_max - 1)) {
        expt_coefficients[sim, expt] <- fixef(model_1)[paste0("as.factor(expt)", expt)]
      }
    }
    
    mean_expt_coeff <- colMeans(expt_coefficients, na.rm = TRUE)
    
    if (t_max == 3) {
      phi <- sigma_alpha^2 / (sigma_alpha^2 + sd_epsilon^2 / n_per)
      theoretical_value_1 <- true_effect[1] - ant_eff - phi * ant_eff
      theoretical_value_2 <- true_effect[2] - ant_eff - 2 * phi * ant_eff
      theoretical_values <- c(0, theoretical_value_1, theoretical_value_2)
    } else {
      theoretical_values <- rep(NA, t_max)
    }
    
    time_points <- 0:(t_max - 1)
    true_values <- c(ant_eff, true_effect)
    
    df <- data.frame(
      Time = time_points,
      True = c(ant_eff, true_effect),
      Estimated = c(0, mean_expt_coeff),
      Theoretical = theoretical_values,
      Pattern = pattern,
      t_max = t_max
    )
    
    df_all <- rbind(df_all, df)
  }
}

df_long <- melt(df_all, id.vars = c("Time", "Pattern", "t_max"), 
                         variable.name = "EffectType", value.name = "Effect")
df_long$EffectType <- factor(df_long$EffectType, 
                                      levels = c("Estimated", "Theoretical", "True"), 
                                      labels = c("Estimated (ETI)", "Theoretical (ETI)", "True"))

df_long$t_max <- factor(df_long$t_max, levels = c(3, 8), labels = c("J = 3", "J = 8"))

p1 <- ggplot(subset(df_long, t_max == "J = 3"), 
             aes(x = Time, y = Effect, color = EffectType)) +
  geom_line(data = subset(df_long, 
                          t_max == "J = 3" & EffectType == "Estimated (ETI)" & Time >= 1), 
            linetype = "dashed") +
  geom_line(data = subset(df_long, 
                          t_max == "J = 3" & EffectType == "True")) +
  geom_point(data = subset(df_long, 
                           t_max == "J = 3" & EffectType == "Estimated (ETI)" & Time >= 1)) +
  geom_point(data = subset(df_long, 
                           t_max == "J = 3" & EffectType == "True")) +
  scale_x_continuous(limits = c(0, 2), breaks = seq(0, 2, 1), 
                     labels = function(x) ifelse(x == 0, "A", x), position = "top") +
  labs(x = "Exposure Time", y = "Trt\nEffect\n\nAnt\nEffect", 
       color = "Effect curve (J = 3)") +
  scale_color_manual(values = c("Estimated (ETI)" = "#3fa535", 
                                "True" = "#009cc4")) +
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

p2 <- ggplot(subset(df_long, t_max == "J = 8"), 
             aes(x = Time, y = Effect, color = EffectType)) +
  geom_line(data = subset(df_long, 
                          t_max == "J = 8" & EffectType == "Estimated (ETI)" & Time >= 1), 
            linetype = "dashed") +
  geom_line(data = subset(df_long, 
                          t_max == "J = 8" & EffectType == "True")) +
  geom_point(data = subset(df_long, 
                           t_max == "J = 8" & EffectType == "Estimated (ETI)" & Time >= 1)) +
  geom_point(data = subset(df_long, 
                           t_max == "J = 8" & EffectType == "True")) +
  scale_x_continuous(limits = c(0, 7), breaks = seq(0, 7, 1), 
                     labels = function(x) ifelse(x == 0, "A", x),
                     position = "top") +
  labs(x = "Exposure Time", y = "Trt\nEffect\n\nAnt\nEffect", 
       color = "Effect curve (J = 8)") +
  scale_color_manual(values = c("Estimated (ETI)" = "#3fa535",
                                "True" = "#009cc4")) +
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

ggsave(filename = "../figures/figure_ETI.pdf", plot = p3, width = 12, height = 5)