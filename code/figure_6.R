library(geomtextpath)
library(ggplot2)
library(gridExtra)
library(scales)

# variance functions (HH, HH-ANT, ETI, ETI-ANT)
general_HH <- function(tau, sigma, I, J, K, Z) {
  rho <- tau^2 / (tau^2 + sigma^2)
  lambda1 <- 1 - rho
  lambda2 <- 1 + (J * K - 1) * rho
  sigma_t2 <- tau^2 + sigma^2
  
  U <- sum(Z)
  W1 <- sum((colSums(Z))^2)
  W2 <- sum((rowSums(Z))^2)
  
  numerator <- (sigma_t2 / K) * I * J * lambda1 * lambda2
  denominator <- ((U^2 + I * J * U - J * W1 - I * W2) * lambda2 -
                  (U^2 - I * W2) * lambda1)
  result <- numerator / denominator
  return(result)
}

general_HH_ANT <- function(tau, sigma, I, J, K, Z, A){
  rho <- tau^2 / (tau^2 + sigma^2)
  lambda1 <- 1 - rho
  lambda2 <- 1 + (J * K - 1) * rho
  sigma_t2 <- tau^2 + sigma^2
  
  U <- sum(Z)
  W1 <- sum((colSums(Z))^2)
  W2 <- sum((rowSums(Z))^2)
  
  W3 <- sum((colSums(A))^2)
  W5 <- sum(colSums(Z) * colSums(A))
  
  numerator <- I * J * lambda1 * lambda2 * sigma_t2 / K
  denominator <- lambda2 * (U^2 + I * J * U - J*W1 - I * W2 - 
                 (J * W5^2)/(I^2 - W3))  lambda1 * (I * W2 - U^2)
  result <- numerator / denominator
  return(result)
}

general_ETI <- function(tau, sigma, I, J, K, Z){
  rho <- tau^2 / (tau^2 + sigma^2)
  lambda1 <- 1 - rho
  lambda2 <- 1 + (J * K - 1) * rho
  sigma_t2 <- tau^2 + sigma^2
  S <- J - 1
  
  U1 <- matrix(0, nrow = S, ncol = 1)
  for (i in 1:I) {
    U1 = U1 + t(Z[((i-1)*J+1):(i*J), ]) %*% matrix(1, nrow = J, ncol = 1)
  }
  
  U2 <- matrix(0, nrow = S, ncol = S)
  for (i in 1:I) {
    U2 = U2 + t(Z[((i-1)*J+1):(i*J), ]) %*% Z[((i-1)*J+1):(i*J), ]
  }
  
  W1_sqrt <- matrix(0, nrow = S, ncol = J)
  for (i in 1:I) {
    W1_sqrt = W1_sqrt + t(Z[((i-1)*J+1):(i*J), ])
  }
  W1 <- W1_sqrt %*% t(W1_sqrt)
  
  W2 <- matrix(0, nrow = S, ncol = S)
  for (i in 1:I) {
    W2 = W2 + t(Z[((i-1)*J+1):(i*J), ]) %*% matrix(1, nrow = J, ncol = 1) %*% 
      matrix(1, nrow = 1, ncol = J) %*% Z[((i-1)*J+1):(i*J), ]
  }
  
  coeff <- I * J * lambda1 * lambda2 * sigma_t2 / (K * S^2)
  
  m <- lambda2 * (U1 %*% t(U1) + I * J * U2 - J * W1 - I * W2) + 
    lambda1 * (I * W2 - U1 %*% t(U1))
  
  result <- coeff * matrix(1, nrow = 1, ncol = S) %*% solve(m) %*% matrix(1, nrow = S, ncol = 1)
  return(result[1,1])
}

general_ETI_ANT <- function(tau, sigma, I, J, K, Z, A){
  rho <- tau^2 / (tau^2 + sigma^2)
  lambda1 <- 1 - rho
  lambda2 <- 1 + (J * K - 1) * rho
  sigma_t2 <- tau^2 + sigma^2
  S <- J - 1
  
  U1 <- matrix(0, nrow = S, ncol = 1)
  for (i in 1:I) {
    U1 = U1 + t(Z[((i-1)*J+1):(i*J), ]) %*% matrix(1, nrow = J, ncol = 1)
  }
  
  U2 <- matrix(0, nrow = S, ncol = S)
  for (i in 1:I) {
    U2 = U2 + t(Z[((i-1)*J+1):(i*J), ]) %*% Z[((i-1)*J+1):(i*J), ]
  }
  
  W1_sqrt <- matrix(0, nrow = S, ncol = J)
  for (i in 1:I) {
    W1_sqrt = W1_sqrt + t(Z[((i-1)*J+1):(i*J), ])
  }
  W1 <- W1_sqrt %*% t(W1_sqrt)
  
  W2 <- matrix(0, nrow = S, ncol = S)
  for (i in 1:I) {
    W2 = W2 + t(Z[((i-1)*J+1):(i*J), ]) %*% matrix(1, nrow = J, ncol = 1) %*% 
      matrix(1, nrow = 1, ncol = J) %*% Z[((i-1)*J+1):(i*J), ]
  }
  
  W3_sqrt <- matrix(0, nrow = 1, ncol = J)
  for (i in 1:I) {
    W3_sqrt = W3_sqrt + A[i, ]
  }
  W3 <- W3_sqrt %*% t(W3_sqrt)
  
  W5 <- W1_sqrt %*% t(W3_sqrt)
  
  coeff <- I * J * lambda1 * lambda2 * sigma_t2 / (K * S^2)
  
  m <- lambda2 * (U1 %*% t(U1) + I * J * U2 - J * W1 - I * W2 
                  - (J * W5 %*% t(W5)) / (I^2 - W3[1,1])) + 
    lambda1 * (I * W2 - U1 %*% t(U1))
  
  result <- coeff * matrix(1, nrow = 1, ncol = S) %*% solve(m) %*% matrix(1, nrow = S, ncol = 1)
  return(result[1,1])
}

# vary rho and J
I <- 32
K <- 100
sigma <- 1

Js <- c(5, 9, 17, 33)
rhos <- seq(0, 0.25, length.out = 2000)

df_all <- data.frame()

# design matrix
build_Z_HH <- function(I, J){
  nseq <- J - 1
  nclus_seq <- I / nseq
  Z_tmp <- matrix(0, nrow = nseq, ncol = J)
  for(i in 1:nseq){
    Z_tmp[i, (i + 1):J] <- 1
  }
  Z_HH <- Z_tmp %x% rep(1, nclus_seq)
  return(Z_HH)
}

# design matrix
build_A_HH <- function(I, J){
  nseq <- J - 1
  nclus_seq <- I / nseq
  A_tmp <- matrix(0, nrow = nseq, ncol = J)
  for(i in 1:nseq){
    A_tmp[i, i] <- 1
  }
  A_HH <- A_tmp %x% rep(1, nclus_seq)
  return(A_HH)
}

# design matrix
build_Z_ETI <- function(I, J){
  nseq <- J - 1
  nclus_seq <- I / nseq
  
  Z_ETI <- NULL
  for (i in 1:nseq) {
    Z_0 <- matrix(0, nrow = i, ncol = nseq)
    Z_1 <- matrix(0, nrow = J - i, ncol = nseq)
    for (j in 1:(J - i)) {
      Z_1[j, j] <- 1
    }
    Z_i <- rbind(Z_0, Z_1)
    
    Z_i_replicated <- do.call(rbind, replicate(nclus_seq, Z_i, simplify = FALSE))
    Z_ETI <- rbind(Z_ETI, Z_i_replicated)
  }
  return(Z_ETI)
}

# design matrix
build_A_ETI <- function(I, J){
  nseq <- J - 1
  nclus_seq <- I / nseq
  A_tmp <- matrix(0, nrow = nseq, ncol = J)
  for(i in 1:nseq){
    A_tmp[i, i] <- 1
  }
  A_ETI <- A_tmp %x% rep(1, nclus_seq)
  return(A_ETI)
}

for(J in Js){
  Z_HH <- build_Z_HH(I, J)
  A_HH <- build_A_HH(I, J)
  Z_ETI <- build_Z_ETI(I, J)
  A_ETI <- build_A_ETI(I, J)
  
  for(rho in rhos){
    tau2 <- rho/(1 - rho)
    tau <- sqrt(tau2)
    
    var_HH <- general_HH(tau, sigma, I, J, K, Z_HH)
    var_HH_ANT <- general_HH_ANT(tau, sigma, I, J, K, Z_HH, A_HH)
    var_ETI <- general_ETI(tau, sigma, I, J, K, Z_ETI)
    var_ETI_ANT <- general_ETI_ANT(tau, sigma, I, J, K, Z_ETI, A_ETI)
    
    df_all <- rbind(df_all,
                    data.frame(rho = rho,
                               J = J,
                               HH = var_HH,
                               HH_ANT = var_HH_ANT,
                               ETI = var_ETI,
                               ETI_ANT = var_ETI_ANT))
  }
}

df_all$ratio_HH <- df_all$HH_ANT/df_all$HH
df_all$ratio_ETI <- df_all$ETI_ANT/df_all$ETI

js <- sort(unique(df_all$J))
rhos <- sort(unique(df_all$rho))

ratioHH_mat <- matrix(NA, nrow = length(rhos), ncol = length(js))
ratioETI_mat <- matrix(NA, nrow = length(rhos), ncol = length(js))

for (i in seq_along(rhos)) {
  for (j in seq_along(js)) {
    ratioHH_mat[i, j] <- df_all$ratio_HH[
      df_all$rho == rhos[i] & df_all$J == js[j]
    ]
    
    ratioETI_mat[i, j] <- df_all$ratio_ETI[
      df_all$rho == rhos[i] & df_all$J == js[j]
    ]
  }
}

p1 <- ggplot(df_all, aes(x = rho, y = J, z = ratio_HH)) +
  geom_textcontour(
    color = "blue",
    size = 6,
    breaks = c(1.9, 1.4, 1.2, 1.1),
    linewidth = 0.8
  ) +
  labs(
    x = expression(rho),
    y = expression(J),
    title = expression(paste("(a) Variance Inflation (HH-ANT vs HH)"))
  ) +
  scale_y_continuous(limits = c(5, 33), breaks = Js) + 
  theme_bw() +
  theme(
    text = element_text(size = 15),
    plot.title = element_text(hjust = 0.5),
    axis.title.y = element_text(angle = 0, vjust = 0.5)
  )

p2 <- ggplot(df_all, aes(x = rho, y = J, z = ratio_ETI)) +
  geom_textcontour(
    color = "blue",
    size = 6,
    breaks = c(3.6, 1.9, 1.4, 1.2),
    linewidth = 0.8
  ) +
  labs(
    x = expression(rho),
    y = expression(J),
    title = expression(paste("(b) Variance Inflation (ETI-ANT vs ETI)"))
  ) +
  scale_y_continuous(limits = c(5, 33), breaks = Js) + 
  theme_bw() +
  theme(
    text = element_text(size = 15),
    plot.title = element_text(hjust = 0.5),
    axis.title.y = element_text(angle = 0, vjust = 0.5)
  )

pdf("../figures/figure_variance_inflation.pdf", width = 12.5, height = 6, paper = "special")
grid.arrange(p1, p2, ncol = 2)
dev.off()