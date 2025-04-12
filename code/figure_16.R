library(geomtextpath)
library(ggplot2)
library(gridExtra)
library(scales)
source("utils_ETI_power_std.R")

# set up
I <- 32
K <- 100
sigma <- 1
alpha <- 0.05
sd_expt <- 1
mu <- 0

gamma_over_Delta_vec <- c(0.2, 0.3, 0.4)
J_values <- c(5, 9)
rho_values <- seq(0.001, 0.25, length.out = 99)
Delta_values <- seq(0, 0.20, length.out = 101)

panel_setup <- data.frame(
  gamma_over_Delta = c(0.2, 0.2, 0.3, 0.3, 0.4, 0.4),
  J = c(5, 9, 5, 9, 5, 9),
  panel_label = c("(a)","(b)","(c)","(d)","(e)","(f)"),
  stringsAsFactors = FALSE
)

df_all <- data.frame()

for (i in seq_len(nrow(panel_setup))) {
  
  gamma_ratio <- panel_setup$gamma_over_Delta[i]
  J_current <- panel_setup$J[i]
  label_i <- panel_setup$panel_label[i]
  
  Z_ETI <- build_Z_ETI(I, J_current)
  A_ETI <- build_A_ETI(I, J_current)
  
  tmp_grid <- expand.grid(rho = rho_values, Delta = Delta_values)
  
  tmp_grid$ratio <- mapply(function(rho, Delta) {
    tau2 <- rho / (1 - rho)
    tau <- sqrt(tau2)
    
    var_ETI <- general_ETI(tau, sigma, I, J_current, K, Z_ETI)
    var_ETI_ANT <- general_ETI_ANT(tau, sigma, I, J_current, K, Z_ETI, A_ETI)
    
    gamma_val <- gamma_ratio * Delta
    
    Delta_ETI <- ETI_est(
      I, J_current, K,
      tau, sigma,
      Delta, gamma_val, sd_expt,
      mu, beta = rep(0, J_current),
      func_form = "sine"
    )
    
    power_ETI <- power_calculation(Delta_ETI, var_ETI, alpha)
    power_ETI_ANT <- power_calculation(Delta, var_ETI_ANT, alpha)
    
    power_ETI_ANT / power_ETI
  }, tmp_grid$rho, tmp_grid$Delta)
  
  tmp_grid$panel_id <- label_i
  tmp_grid$J <- J_current
  tmp_grid$gamma_over_Delta <- gamma_ratio
  
  df_all <- rbind(df_all, tmp_grid)
}

df_all$panel_id <- factor(df_all$panel_id, levels = panel_setup$panel_label)

plots_list <- list()

label_fmt <- label_number(accuracy = 0.01, trim = FALSE)

for (panel_lab in levels(df_all$panel_id)) {
  
  df_sub <- subset(df_all, panel_id == panel_lab)
  
  J_cur <- df_sub$J[1]
  gamma_ratio <- df_sub$gamma_over_Delta[1]
  
  if(panel_lab == "(a)"){
    custom_breaks <- c(0.75, 0.80, 0.85, 0.90, 0.95)
  } else if(panel_lab == "(b)"){
    custom_breaks <- c(0.99, 0.96, 0.95)
  } else if(panel_lab == "(c)"){
    custom_breaks <- c(1.05, 1.15, 1.25, 1.30)
  } else if(panel_lab == "(d)"){
    custom_breaks <- c(1.05, 1.10,  1.20)
  } else if(panel_lab == "(e)"){
    custom_breaks <- c(1.2, 1.6, 2.0, 2.4, 2.8, 3)
  } else if(panel_lab == "(f)"){
    custom_breaks <- c(1.1, 1.4, 1.6)
  }
  
  panel_title <- bquote(
    .(panel_lab) ~ " " ~ Power[ETI-ANT] / Power[ETI] ~
      "(" ~ J == .(J_cur) ~ ", " ~ gamma / Delta ~ " = " ~ .(gamma_ratio) ~ ")"
  )
  
  p <- ggplot(df_sub, aes(x = rho, y = Delta, z = ratio)) +
    geom_textcontour(
      aes(label = label_fmt(after_stat(level))),
      color = "blue",
      size = 5,
      breaks = custom_breaks,
      linewidth = 0.8
    ) +
    scale_x_continuous(
      limits = c(0, 0.25),
      breaks = c(0.0, 0.1, 0.2, 0.25)
    ) +
    scale_y_continuous(
      limits = c(0, 0.20),
      breaks = c(0.00, 0.05, 0.10, 0.15, 0.20)
    ) +
    labs(
      x = expression(rho),
      y = expression(Delta),
      title = panel_title
    ) +
    theme_bw() +
    theme(
      text = element_text(size = 12),
      plot.title = element_text(hjust = 0.5),
      axis.title.y = element_text(angle = 0, vjust = 0.5)
    )
  
  plots_list[[panel_lab]] <- p
}

pdf("figure_power_ratio_fixed_ratio.pdf", width = 12, height = 15, paper = "special")
grid.arrange(
  grobs = plots_list,
  nrow = 3, ncol = 2
)
dev.off()