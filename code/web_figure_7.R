library(ggplot2)
library(geomtextpath)
library(gridExtra)
library(scales)
source("utils_HH_power_std.R")

# set up
I <- 32
K <- 100
sigma <- 1
alpha <- 0.05
gamma_over_delta <- 0.4

Js <- c(5, 9, 17, 33)
rhos <- seq(0.001, 0.25, length.out = 399)
delta_s <- seq(0, 0.1, length.out = 401)

df_all <- expand.grid(
  J = Js,
  rho = rhos,
  delta = delta_s
)

df_all$ratio <- mapply(function(J, rho, delta) {
  Q <- J - 1
  tau2 <- rho / (1 - rho)
  tau <- sqrt(tau2)
  phi <- tau2 / (tau2 + 1 / K)
  
  var_HH <- standard_HH(tau, sigma, I, J, K)
  var_HH_ANT <- standard_HH_ANT(tau, sigma, I, J, K)
  
  gamma <- gamma_over_delta * delta
  delta_HH <- delta - 6 * (1 + phi * Q) / (J * (2 + phi * Q)) * gamma
  
  power_HH <- power_calculation(delta_HH, var_HH, alpha)
  power_HH_ANT <- power_calculation(delta, var_HH_ANT, alpha)
  
  ratio_value <- power_HH_ANT / power_HH
  return(ratio_value)
},
df_all$J, df_all$rho, df_all$delta)

panel_labels <- c("(a)", "(b)", "(c)", "(d)")
labels_map   <- setNames(panel_labels, Js)
df_all$panel_id <- labels_map[as.character(df_all$J)]

plots_list <- list()

label_fmt <- label_number(accuracy = 0.01, trim = FALSE)

for (id in unique(df_all$panel_id)) {
  
  df_sub <- subset(df_all, panel_id == id)
  
  if(id == "(a)"){
    custom_breaks <- c(1.05, 1.1, 1.2, 1.3)
  } else if(id == "(b)"){
    custom_breaks <- c(1.05, 1.1, 1.15)
  } else if(id == "(c)"){
    custom_breaks <- c(1.02, 1.06, 1.08)
  } else if(id == "(d)"){
    custom_breaks <- c(1.01, 1.04)
  }
  
  J_current <- df_sub$J[1]
  panel_title <- bquote(
    .(id) ~ " " ~ Power^{"HH-ANT"} / Power^{"HH"} ~
      "(" ~ J == .(J_current) ~ ", " ~ gamma^{"HH-ANT"} / delta^{"HH-ANT"} ~ " = " ~ .(gamma_over_delta) ~ ")"
  )
  
  p <- ggplot(df_sub, aes(x = rho, y = delta, z = ratio)) +
    geom_textcontour(
      aes(label = label_fmt(after_stat(level))),
      color = "blue",
      size = 5,
      breaks = custom_breaks,
      linewidth = 0.8
    ) +
    scale_y_continuous(limits = c(0, 0.1),  breaks = c(0.00, 0.02, 0.04, 0.06, 0.08, 0.10)) +
    labs(
      x = expression(rho),
      y = expression(delta^{"HH-ANT"}),
      title = panel_title
    ) +
    theme_bw() +
    theme(
      text = element_text(size = 12),
      plot.title = element_text(hjust = 0.5),
      axis.title.y = element_text(angle = 0, vjust = 0.5)
    )
  
  plots_list[[id]] <- p
}

pdf("../figures/figure_power_ratio_0.4.pdf", width = 14, height = 14, paper = "special")
grid.arrange(grobs = plots_list, ncol = 2)
dev.off()