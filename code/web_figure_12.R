library(ggplot2)
library(reshape2)

# parameters
t_max <- 8

# define five functional forms
true_effects <- list(
  'Constant' = rep(2, t_max - 1),
  'Curved' = c(0.5, 2 - 0.5^(0:4), 2),
  'Lagged' = c(rep(0.5, 2), rep(2, t_max - 3)),
  'Partially convex' = c(0.25 + 0.25 * (1:3), 1.5, rep(2, 3)),
  'Sinusoidal' = -1.41 * sin(2 * pi * (0:6) / 7) + 0.12
)

# create data frame for plotting
df_all <- data.frame()
for (pattern in names(true_effects)) {
  true_effect <- true_effects[[pattern]]
  
  df <- data.frame(
    Time = 1:(t_max - 1),
    True = true_effect,
    Pattern = pattern
  )
  df_all <- rbind(df_all, df)
}

df_long <- melt(df_all, id.vars = c("Time", "Pattern"), 
                variable.name = "EffectType", value.name = "Effect")

# create plot matching your style
p <- ggplot(df_long, aes(x = Time, y = Effect)) +
  geom_line(color = "#009cc4") +
  geom_point(color = "#009cc4") +
  scale_y_continuous(limits = c(-2, 3)) +
  scale_x_continuous(breaks = 1:7, position = "top") +
  labs(x = "Exposure Time",
       y = "Treatment Effect") +
  facet_wrap(~Pattern, ncol = 5, strip.position = "bottom",
             labeller = labeller(Pattern = c(
               "Constant" = "(a) Constant",
               "Curved" = "(b) Curved",
               "Lagged" = "(c) Lagged",
               "Partially convex" = "(d) Partially convex",
               "Sinusoidal" = "(e) Sinusoidal"))) +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x.top = element_text(size = 12),
        axis.text.x.top = element_text(size = 12),
        strip.text = element_text(size = 12),
        strip.placement = "outside")

ggsave("../figures/figure_functional_forms.pdf", plot = p, width = 12, height = 2.5)