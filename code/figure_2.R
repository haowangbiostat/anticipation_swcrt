library(dplyr)
library(ggplot2)
library(lme4)
library(reshape2)

# number of individuals per cluster-period
n_per <- 50

# number of time periods
t_max = 8

# number of clusters per sequence
each <- 40

# number of clusters
i_max <- (t_max-1)*each

# maximum number of exposure times observed
expt_max <- t_max - 1

# sd of individual level heterogeneity
sd_epsilon = 0.5

# anticipation treatment effect
ant_eff = 0.5

# cluster-level heterogeneity
sigma_alpha = 0.141

# treatment effect
overall_eff = 2

# simulate the data
run_simulation <- function(true_effect) {
  alphas <- rnorm(i_max, 0, sigma_alpha)
  
  data = expand.grid(i = 1:i_max, t = 1:t_max, id = 1:n_per)
  data$sequence <- ceiling(data$i / each)
  data$trt <- ifelse(data$t < (data$sequence + 1), 0, 1)
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

# true effect
true_effect_constant <- rep(2, t_max - 1)
true_effect_lagged <- c(0.5, 0.5, rep(2, t_max-3))
true_effect_curved <- c(0.5, 1.0, 1.5, 1.75, 1.875, 1.9375, 2)
true_effect_partially_convex <- c(0.5, 0.75, 1, 1.5, 2, 2, 2)

true_effects <- list(
  'Constant' = true_effect_constant,
  'Lagged' = true_effect_lagged,
  'Curved' = true_effect_curved,
  'Partially convex' = true_effect_partially_convex
)

n_sim <- 1000
df_all <- data.frame()

set.seed(123)
for (pattern in names(true_effects)) {
  true_effect <- true_effects[[pattern]]
  trt_coefficients <- numeric(n_sim)
  
  for (sim in 1:n_sim) {
    data <- run_simulation(true_effect)
    model_1 <- lmer(Y ~ trt + as.factor(t) + (1 | i), data = data)
    trt_coefficients[sim] <- fixef(model_1)["trt"]
  }
  
  # average treatment effect over n_sim simulations
  mean_trt_coeff <- mean(trt_coefficients)
  
  # expectation of the treatment effect under the HH model
  phi <- sigma_alpha* sigma_alpha / 
    (sigma_alpha * sigma_alpha + sd_epsilon*sd_epsilon/n_per)
  if (pattern == "Constant") {
    theoretical_value <- overall_eff - 6 * (1 + phi * expt_max) / 
      (t_max * (2 + phi * expt_max)) * ant_eff
  } else {
    Q <- t_max - 1
    theoretical_value <- 0
    for (j in 1:Q) {
      theoretical_value <- theoretical_value + 
        (6 * (j - Q - 1) * ((1 + 2 * phi * Q) * j - (1 + phi + phi * Q) * Q) /
           (Q * (Q + 1) * (phi * Q^2 + 2 * Q - phi * Q - 2))) * true_effect[j] -
        (6 * (phi * Q^2 - phi * Q - 2 + 2 * j) /
           (Q * (phi * Q^3 - phi * Q + 2 * Q^2 - 2))) * ant_eff
    }
  }
  
  time_points <- 0:(t_max - 1)
  df <- data.frame(
    Time = time_points,
    True = c(ant_eff, true_effect),
    Estimated = c(mean_trt_coeff, rep(mean_trt_coeff, t_max - 1)),
    Theoretical = c(theoretical_value, rep(theoretical_value, t_max - 1)),
    Pattern = pattern
  )
  
  df_all <- rbind(df_all, df)
}

df_long <- melt(df_all, id.vars = c("Time", "Pattern"), 
                variable.name = "EffectType", value.name = "Effect")
df_long$EffectType <- factor(df_long$EffectType, 
                             levels = c("Estimated", "Theoretical", "True"), 
                             labels = c("Estimated (HH)", "Theoretical (HH)", "True"))

p <- ggplot(df_long, aes(x = Time, y = Effect, color = EffectType)) +
  geom_line(data = subset(df_long, EffectType == "Estimated (HH)" & Time >= 1), linetype = "dashed") +
  geom_line(data = subset(df_long, EffectType == "True")) +
  geom_point(data = subset(df_long, EffectType == "True")) +
  scale_y_continuous(limits = c(-0.05, 2.02)) +
  scale_x_continuous(breaks = seq(0, t_max - 1, 1), 
                     labels = function(x) ifelse(x == 0, "A", x),
                     position = "top") +
  scale_color_manual(values = c("Estimated (HH)" = "#3fa535",
                                "True" = "#009cc4")) +
  labs(x = "Exposure Time",
       y = "Effect",
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
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x.top = element_text(size = 12),
        axis.text.x.top = element_text(size = 12),
        strip.text = element_text(size = 12),
        strip.placement = "outside")

ggsave(filename = "../figures/figure_HH.pdf", plot = p, width = 12, height = 2.5)