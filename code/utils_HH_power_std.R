standard_HH <- function(tau, sigma, I, J, K) {
  Q = J - 1
  
  rho <- tau^2 / (tau^2 + sigma^2)
  lambda1 <- 1 - rho
  lambda2 <- 1 + (J * K - 1) * rho
  sigma_t2 <- tau^2 + sigma^2
  
  numerator <- 12 * Q * sigma_t2 * lambda1 * lambda2
  denominator <- K * I * (Q - 1) * (Q * lambda1 + (Q+2) * lambda2)
  
  result <- numerator / denominator
  return(result)
}

standard_HH_ANT <- function(tau, sigma, I, J, K) {
  Q = J - 1
  
  rho <- tau^2 / (tau^2 + sigma^2)
  lambda1 <- 1 - rho
  lambda2 <- 1 + (J * K - 1) * rho
  sigma_t2 <- tau^2 + sigma^2
  
  numerator <- 12 * Q * sigma_t2 * lambda1 * lambda2
  denominator <- K * I * (Q - 1) * (Q * lambda1 + (Q-1) * lambda2)
  
  result <- numerator / denominator
  return(result)
}

power_calculation <- function(delta, var_delta, alpha=0.05) {
  # z value
  Z_alpha <- qnorm(1 - alpha / 2)
  # two-sided power
  power <- pnorm(delta / sqrt(var_delta) - Z_alpha) + 
    pnorm(-delta / sqrt(var_delta) - Z_alpha)
  return(power)
}