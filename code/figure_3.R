library(dplyr)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)

coeff_delta <- function(phi, Q, j) {
  numerator <- 6 * (j - Q - 1) * ((1 + 2 * phi * Q) * j - (1 + phi + phi * Q) * Q)
  denominator <- Q * (Q + 1) * (phi * Q^2 + 2 * Q - phi * Q - 2)
  return(numerator / denominator)
}

coeff_gamma <- function(phi, Q, j) {
  numerator <- 6 * (phi * Q^2 - phi * Q - 2 + 2 * j)
  denominator <- Q * (phi * Q^3 + 2 * Q^2 - phi * Q - 2)
  return(-numerator / denominator)
}

run_simulation <- function(Q, expression_func) {
  phi_values <- seq(0, 1, length.out = 200)
  data <- expand.grid(phi = phi_values, j = 1:Q) %>%
          mutate(value = mapply(expression_func, phi, Q, j),
                     j = as.factor(j),
                     Q = as.factor(Q))
  return(data)
}

Qs <- c(3, 5, 7)
df_delta <- do.call(rbind, lapply(Qs, run_simulation, coeff_delta))
df_gamma <- do.call(rbind, lapply(Qs, run_simulation, coeff_gamma))

colors <- brewer.pal(n = 7, name = "Set1")

generate_plot <- function(data, y_label, y_limits, facet_labels) {
  ggplot(data, aes(x = phi, y = value, color = j)) +
    geom_line() +
    scale_y_continuous(limits = y_limits, n.breaks = 4) +
    scale_x_continuous(breaks = seq(0, 1, by = 0.2), position = "top") +
    scale_color_manual(values = colors) +
    labs(x = expression(phi), y = y_label, color = "j") +
    facet_wrap(~ Q, ncol = 3, strip.position = "bottom", 
               labeller = labeller(Q = facet_labels)) +
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
          strip.placement = "outside") +
    guides(color = guide_legend(nrow = 1))
}

p1 <- generate_plot(df_delta, y_label = expression("Coefficients of " * delta), 
                    y_limits = c(-0.3, 1.5),
                    facet_labels = c("3" = "(a) Q = 3", 
                                     "5" = "(b) Q = 5", "7" = "(c) Q = 7"))

p2 <- generate_plot(df_gamma, y_label = expression("Coefficients of " * gamma), 
                    y_limits = c(-0.6, 0),
                    facet_labels = c("3" = "(d) Q = 3", 
                                     "5" = "(e) Q = 5", "7" = "(f) Q = 7"))

p3 <- grid.arrange(p1, p2, nrow = 2)

ggsave(filename = "../figures/figure_coeff_HH.pdf", plot = p3, width = 12, height = 5)