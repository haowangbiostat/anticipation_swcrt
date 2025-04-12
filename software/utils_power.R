# variance of HH
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

# variance of HH-ANT
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
                              (J * W5^2)/(I^2 - W3)) + lambda1 * (I * W2 - U^2)
  result <- numerator / denominator
  return(result)
}

# variance of ETI
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

# variance of ETI-ANT
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

# create design matrix (HH)
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

# create design matrix (HH-ANT)
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

# create design matrix (ETI)
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

# create design matrix (ETI-ANT)
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

get_tau <- function(rho, sigma_sq){
  return(sqrt(rho * sigma_sq / (1 - rho)))
}

power_calculation <- function(est, var, alpha=0.05) {
  # z value
  Z_alpha <- qnorm(1 - alpha / 2)
  # two-sided power
  power <- pnorm(est / sqrt(var) - Z_alpha) + 
    pnorm(-est / sqrt(var) - Z_alpha)
  return(power)
}