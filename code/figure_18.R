library(dplyr)
library(ggplot2)
library(gridExtra)
library(lme4)
library(reshape2)

# number of individuals per cluster-period
n_per <- 50

# number of time periods
t_max <- 8

# number of clusters per sequence
each <- 40

# number of clusters
i_max <- (t_max - 1) * each

# maximum number of exposure times
expt_max <- t_max - 1

# sd of individual level heterogeneity
sd_epsilon <- 0.5

# anticipation treatment effect
ant_eff <- 0.5

# cluster-level heterogeneity
sigma_alpha <- 0.141

# treatment effect
overall_eff <- 2

# simulation
run_simulation <- function(true_effect, ell) {
  alphas <- rnorm(i_max, 0, sigma_alpha)
  
  data <- expand.grid(i = 1:i_max, t = 1:t_max, id = 1:n_per)
  data$sequence <- ceiling(data$i / each)
  data$trt <- ifelse(data$t < (data$sequence + 1), 0, 1)
  
  # anticipation indicator based on ell
  data$ant <- ifelse(data$t >= pmax(1, data$sequence - (ell - 1)) & data$t <= data$sequence, 1, 0)
  
  data$expt <- data$t - data$sequence
  data[data$expt < 0, "expt"] <- 0
  data$q <- ceiling(data$i / 2)
  
  t_fixef_osc <- 0.5 * sin(pi * 2 * ((1:t_max) - 1) / (t_max - 1))
  
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

# true treatment effect
true_effect_constant <- rep(overall_eff, t_max - 1)

true_effects <- list(
  'First Order' = true_effect_constant,
  'Second Order' = true_effect_constant,
  'Third Order' = true_effect_constant
)

n_sim <- 1000
df_all <- data.frame()

set.seed(123)
for (pattern in names(true_effects)) {
  true_effect <- true_effects[[pattern]]
  
  current_ell <- switch(pattern,
                        "First Order" = 1,
                        "Second Order" = 2,
                        "Third Order" = 3)
  
  trt_coefficients <- numeric(n_sim)
  
  for (sim in 1:n_sim) {
    data <- run_simulation(true_effect, current_ell)
    model <- lmer(Y ~ trt + as.factor(t) + (1 | i), data = data)
    trt_coefficients[sim] <- fixef(model)["trt"]
  }
  
  # average treatment effect
  mean_trt_coeff <- mean(trt_coefficients)
  
  # theoretical value
  phi <- sigma_alpha^2 / (sigma_alpha^2 + sd_epsilon^2 / n_per)
  Q <- expt_max
  ell <- current_ell
  
  num <- -ell * (
    6  * phi * Q^3
    - 9 * phi * ell * Q^2
    + 3 * phi * Q^2
    + 6 * Q^2
    + 4 * phi * ell^2 * Q
    - 3 * phi * ell * Q
    - 6 * ell * Q
    - phi * Q
    + 2 * ell^2
    - 2
  )
  
  denom <- Q * (Q + 1) * (phi * Q^2 + 2 * Q - phi * Q - 2)
  
  theoretical_value <- overall_eff + num / denom * ant_eff
  
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

df_all$Pattern <- factor(df_all$Pattern, 
                         levels = c("First Order", "Second Order", 
                                    "Third Order", "Fourth Order"))

df_long <- melt(df_all, id.vars = c("Time", "Pattern"),
                variable.name = "EffectType", value.name = "Effect")
df_long$EffectType <- factor(df_long$EffectType,
                             levels = c("Estimated", "Theoretical", "True"),
                             labels = c("Estimated (HH)", "Theoretical (HH)", "True"))

df_long <- rbind(df_long, list(-1, "Second Order", "True", 0.5))
df_long <- rbind(df_long, list(-1, "Third Order", "True", 0.5))
df_long <- rbind(df_long, list(-2, "Third Order", "True", 0.5))

df_long$Pattern <- factor(df_long$Pattern, 
                          levels = c("First Order", "Second Order", "Third Order"))

p <- ggplot(df_long, aes(x = Time, y = Effect, color = EffectType)) +
  geom_line(data = subset(df_long, EffectType == "Estimated (HH)" & Time >= 1), 
            linetype = "dashed") +
  geom_line(data = subset(df_long, EffectType == "True")) +
  geom_point(data = subset(df_long, EffectType == "True")) +
  scale_y_continuous(limits = c(-0.05, 2.02)) +
  scale_x_continuous(
    breaks = c(0, -1, -2, 1, 2, 3, 4, 5, 6, 7),
    labels = c(expression(A[1]), expression(A[2]), expression(A[3]), 
               1, 2, 3, 4, 5, 6, 7),
    position = "top"
  ) +
  scale_color_manual(values = c("Estimated (HH)" = "#3fa535", 
                                "True" = "#009cc4")) +
  labs(x = "Exposure Time", y = "Effect", color = "Effect curve") +
  facet_wrap(~Pattern, ncol = 3, strip.position = "bottom",
             labeller = labeller(Pattern = c(
               "First Order" = "(a) First Order",
               "Second Order" = "(b) Second Order",
               "Third Order" = "(c) Third Order"
             ))) +
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

ggsave(filename = "figure_HH_higher_order.pdf", plot = p, width = 12, height = 2.5)