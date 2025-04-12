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
Delta_vec <- c(0.01, 0.04, 0.07)
mu <- 0

J_values <- c(5, 9)
rho_values <- seq(0.001, 0.25, length.out = 99)
gamma_over_Delta_values <- seq(0, 1, length.out = 101)

panel_setup <- data.frame(
  Delta = rep(Delta_vec, each = length(J_values)),
  J = rep(J_values, times = length(Delta_vec)),
  panel_label = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)"),
  stringsAsFactors = FALSE
)

df_all <- data.frame()

for (i in seq_len(nrow(panel_setup))) {
  
  Delta_cur <- panel_setup$Delta[i]
  J_cur <- panel_setup$J[i]
  label_i <- panel_setup$panel_label[i]
  
  Z_ETI <- build_Z_ETI(I, J_cur)
  A_ETI <- build_A_ETI(I, J_cur)
  
  tmp_grid <- expand.grid(
    rho = rho_values,
    gamma_over_Delta = gamma_over_Delta_values
  )
  
  tmp_grid$ratio <- mapply(function(rho, gamma_ratio) {
    tau2 <- rho / (1 - rho)
    tau <- sqrt(tau2)
    
    var_ETI <- general_ETI(tau, sigma, I, J_cur, K, Z_ETI)
    var_ETI_ANT <- general_ETI_ANT(tau, sigma, I, J_cur, K, Z_ETI, A_ETI)
    
    gamma_val <- gamma_ratio * Delta_cur
    
    Delta_ETI <- ETI_est(
      I, J_cur, K,
      tau, sigma,
      Delta_cur, gamma_val, sd_expt,
      mu, beta = rep(0, J_cur),
      func_form = "sine"
    )
    
    power_ETI <- power_calculation(Delta_ETI, var_ETI, alpha)
    power_ETI_ANT <- power_calculation(Delta_cur, var_ETI_ANT, alpha)
    
    ratio_val <- power_ETI_ANT / power_ETI
    return(ratio_val)
    
  }, tmp_grid$rho, tmp_grid$gamma_over_Delta)
  
  tmp_grid$panel_id <- label_i
  tmp_grid$J <- J_cur
  tmp_grid$Delta <- Delta_cur
  
  df_all <- rbind(df_all, tmp_grid)
}

df_all$panel_id <- factor(df_all$panel_id, levels = panel_setup$panel_label)

plots_list <- list()

label_fmt <- label_number(accuracy = 0.01, trim = FALSE)

for (panel_lab in levels(df_all$panel_id)) {
  
  df_sub <- subset(df_all, panel_id == panel_lab)
  J_cur <- df_sub$J[1]
  Delta_cur <- df_sub$Delta[1]
  
  if(panel_lab == "(a)"){
    custom_breaks <- c(0.95, 1, 1.05, 1.10)
  } else if(panel_lab == "(b)"){
    custom_breaks <- c(0.95, 1, 1.05, 1.10, 1.15)
  } else if(panel_lab == "(c)"){
    custom_breaks <- c(0.75, 1, 1.5, 2, 3)
  } else if(panel_lab == "(d)"){
    custom_breaks <- c(0.75, 1, 1.50, 2, 2.5, 3)
  } else if(panel_lab == "(e)"){
    custom_breaks <- c(0.75, 1, 1.5, 2, 3)
  } else if(panel_lab == "(f)"){
    custom_breaks <- c(0.75, 1, 2, 3, 4, 5, 6, 7)
  }
  
  panel_title <- bquote(
    .(panel_lab) ~ " " ~ Power[ETI-ANT] / Power[ETI] ~
      "(" ~ J == .(J_cur) ~ ", " ~ Delta == .(Delta_cur) ~ ")"
  )
  
  p <- ggplot(df_sub, aes(x = rho, y = gamma_over_Delta, z = ratio)) +
    geom_textcontour(
      aes(label = label_fmt(after_stat(level))),
      color = "blue",
      size = 5,
      linewidth = 0.8,
      breaks = custom_breaks
    ) +
    scale_x_continuous(
      limits = c(0, 0.25),
      breaks = c(0.0, 0.05, 0.1, 0.15, 0.2, 0.25)
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.2)
    ) +
    labs(
      x = expression(rho),
      y = expression(frac(gamma, Delta)),
      title = panel_title
    ) +
    theme_bw() +
    theme(
      text = element_text(size = 12),
      plot.title = element_text(hjust = 0.5),
      axis.title.y  = element_text(angle = 0, vjust = 0.5)
    )
  
  plots_list[[panel_lab]] <- p
}

pdf("../figures/figure_power_ratio_fixed_Delta.pdf", width = 12, height = 15, paper = "special")
grid.arrange(
  grobs = plots_list,
  nrow  = 3,
  ncol  = 2
)
dev.off()