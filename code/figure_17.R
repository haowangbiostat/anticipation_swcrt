library(geomtextpath)
library(ggplot2)
library(gridExtra)
library(showtext)
showtext_auto()

coeff_gamma <- function(phi, Q, ell) {
  num <- -ell * (
    6*phi*Q^3 - 9*phi*ell*Q^2 + 3*phi*Q^2 + 6*Q^2 +
      4*phi*ell^2*Q - 3*phi*ell*Q - 6*ell*Q - phi*Q +
      2*ell^2 - 2
  )
  
  den <- Q * (Q + 1) * (phi*Q^2 + 2*Q - phi*Q - 2)
  
  return(num / den)
}

phi_values <- seq(0, 1, length.out = 5000)

# Q = 4
df_Q4 <- expand.grid(
  phi = phi_values,
  ell = 1:3
)
df_Q4$value <- with(df_Q4, coeff_gamma(phi, 4, ell))
df_Q4$Q <- 4

# Q = 8
df_Q8 <- expand.grid(
  phi = phi_values,
  ell = 1:7
)
df_Q8$value <- with(df_Q8, coeff_gamma(phi, 8, ell))
df_Q8$Q <- 8

df_all <- rbind(df_Q4, df_Q8)

p1 <- ggplot(subset(df_all, Q == 4), aes(x = phi, y = ell, z = value)) +
  geom_textcontour(
    color = "blue",
    size = 6,
    breaks = c(-0.7, -0.8, -0.9, -1, -1.1, -1.2),
    linewidth = 0.8
  ) +
  labs(
    x = expression(phi),
    y = expression("\u2113"),
    title = expression(paste("(a) ", omega[list(HH)]^{"HH-ANT, \u2113"}, " (Q = 4)"))
  ) +
  scale_y_continuous(breaks = 1:4) +
  theme_bw() +
  theme(
    text = element_text(size = 15),
    plot.title = element_text(hjust = 0.5),
    axis.title.y = element_text(angle = 0, vjust = 0.5)
  )

p2 <- ggplot(subset(df_all, Q == 8), aes(x = phi, y = ell, z = value)) +
  geom_textcontour(
    color = "blue",
    size = 6,
    breaks = c(-0.5, -0.6, -0.7, -0.8, -0.9, -1, -1.1, -1.2),
    linewidth = 0.8
  ) +
  labs(
    x = expression(phi),
    y = expression("\u2113"),
    title = expression(paste("(b) ", omega[list(HH)]^{"HH-ANT, \u2113"}, " (Q = 8)"))
  ) +
  scale_y_continuous(breaks = 1:8) +
  theme_bw() +
  theme(
    text = element_text(size = 15),
    plot.title = element_text(hjust = 0.5),
    axis.title.y = element_text(angle = 0, vjust = 0.5)
  )

pdf("../figures/figure_coeff_HH_higher_order.pdf", width = 12.5, height = 6, paper = "special")
grid.arrange(p1, p2, ncol = 2)
dev.off()