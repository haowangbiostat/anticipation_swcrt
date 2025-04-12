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

# create design matrix
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

# create design matrix
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

power_calculation <- function(Delta, var_Delta, alpha=0.05) {
  # z value
  Z_alpha <- qnorm(1 - alpha / 2)
  # two-sided power
  power <- pnorm(Delta / sqrt(var_Delta) - Z_alpha) + 
    pnorm(-Delta / sqrt(var_Delta) - Z_alpha)
  return(power)
}

scale2 <- function(x, Delta, sd_expt) {
  out <- sd_expt * (x - mean(x)) / sd(x) + Delta
  return(out)
}

# E[Y_j(q)] = mu + beta[j] + (gamma if j == q) + (delta if j == q+s)
E_Y_jq <- function(j, q, s, mu, beta, gamma, delta) {
  val <- mu + beta[j]
  if (j == q) {
    val <- val + gamma
  }
  if (j == q + s) {
    val <- val + delta
  }
  return(val)
}

# Q x J matrix
create_Z <- function(Q, J) {
  Z <- matrix(0, nrow = Q, ncol = J)
  for (i in 1:Q) {
    for (j in 1:J) {
      if (j > i) {
        Z[i, j] <- 1
      }
    }
  }
  return(Z)
}

# Q x 1 vector
create_D2 <- function(Z) {
  sums <- rowSums(Z)
  D2   <- matrix(sums, nrow = nrow(Z), ncol = 1)
  return(D2)
}

# Q x Q matrix
create_Xi <- function(Q, x, y) {
  Xi <- matrix(0, nrow = Q, ncol = Q)
  for (a in 1:Q) {
    for (b in 1:Q) {
      if (a == b) {
        # diagonal
        Xi[a, b] <- (Q + 1 - a) * (x - y)
      } else {
        # off-diagonal
        Xi[a, b] <- -((Q + 1) - max(a, b)) * y
      }
    }
  }
  return(Xi)
}

# (J+Q) x (J+Q) matrix
create_big_matrix <- function(J, Q, x, y, Xi, Z, D2) {
  # top-left block
  I_J   <- diag(1, nrow = J, ncol = J)
  oneJ  <- matrix(1, nrow = J, ncol = 1)
  top_left <- Q*x*I_J - Q*y*(oneJ %*% t(oneJ))
  
  # top-right block
  top_right <- x * t(Z) - y * (oneJ %*% t(D2))
  
  # bottom-left block
  bottom_left <- x * Z - y * (D2 %*% t(oneJ))
  
  # bottom-right block
  bottom_right <- Xi
  
  # bind them all
  top <- cbind(top_left, top_right)
  bot <- cbind(bottom_left, bottom_right)
  M   <- rbind(top, bot)
  return(M)
}

# (J+Q) x 1 vector
create_big_vector <- function(J, Q, s, x, y, mu, beta, gamma, delta) {
  EY <- matrix(0, nrow = J, ncol = Q)
  for (j in 1:J) {
    for (qq in 1:Q) {  # use 'qq' instead of 'q' to avoid confusion
      EY[j, qq] <- E_Y_jq(j, qq, s, mu, beta, gamma, delta)
    }
  }
  
  # sum_q E[Y_j(q)] for each j
  sum_q_EYj <- rowSums(EY)
  
  # sum_q sum_j E[Y_j(q)]
  sum_q_sum_j_EY <- sum(EY)
  
  # top_J: length J
  top_J <- x * sum_q_EYj - y * sum_q_sum_j_EY
  
  # bottom_Q: length Q
  bot_Q <- numeric(Q)
  for (k in 1:Q) {
    # x part
    temp_sum_x <- 0
    for (qq in 1:Q) {
      for (j in 1:J) {
        if (j == qq + k) {
          temp_sum_x <- temp_sum_x + EY[j, qq]
        }
      }
    }
    # y part
    temp_sum_y <- 0
    for (qq in 1:(Q - k + 1)) {
      for (j in 1:J) {
        temp_sum_y <- temp_sum_y + EY[j, qq]
      }
    }
    bot_Q[k] <- x * temp_sum_x - y * temp_sum_y
  }
  
  return(c(top_J, bot_Q))
}

# (J+Q) x 1 vector
e_s_vector <- function(J, Q, s) {
  e_s <- rep(0, J + Q)
  e_s[J + s] <- 1
  return(e_s)
}

# E_hat_delta_s for a given s
E_hat_delta_s <- function(J, Q, s,
                          mu, beta, gamma, delta,
                          x, y,
                          Z, D2){
  Xi <- create_Xi(Q = Q, x = x, y = y)
  M  <- create_big_matrix(J = J, Q = Q, x = x, y = y, Xi = Xi, Z = Z, D2 = D2)
  v  <- create_big_vector(J = J, Q = Q, s = s, x = x, y = y,
                          mu = mu, beta = beta, gamma = gamma, delta = delta)
  
  invM <- solve(M)
  e_s  <- e_s_vector(J = J, Q = Q, s = s)
  
  val <- as.numeric(e_s %*% invM %*% v)
  return(val)
}

# phi = tau^2 / (tau^2 + sigma^2 / K)
# eta^2 = tau^2 / phi
# x = 1 / [eta^2 (1 - phi)]
# y = phi / [eta^2 (1 - phi)(1 + phi J - phi)]
calc_x_y <- function(tau2, sigma2, K, J) {
  phi  <- tau2 / (tau2 + sigma2 / K)
  eta2 <- tau2 / phi
  
  x <- 1 / (eta2 * (1 - phi))
  y <- phi / (eta2 * (1 - phi) * (1 + phi * J - phi))
  
  list(phi = phi, eta2 = eta2, x = x, y = y)
}

get_delta_s <- function(delta_true, s) {
  delta_true[s]
}

# expectation of the estimated TATE under the ETI model
ETI_est <- function(I, J, K,
                    tau, sigma,
                    Delta, gamma, sd_expt, 
                    mu, beta, func_form = "sine") {
  Q <- J - 1
  tau2 <- tau * tau
  sigma2 <- sigma * sigma
  store <- calc_x_y(tau2, sigma2, K, J)
  x <- store$x
  y <- store$y
  
  e <- 1:Q
  if (func_form == "sine") {
      delta_raw <- -0.5 * sin(2*pi*(e-1)/(J-2)) + log(1.2)
      delta_true <- scale2(delta_raw, Delta, sd_expt)
  } else if (func_form == "lagged") {
      delta_raw <- ifelse(e <= Q/2, 0.2 * Delta, Delta)
      delta_true <- scale2(delta_raw, Delta, sd_expt)  
  } else if (func_form == "curved") {
      t <- (e - 1) / (Q - 1)
      a <- 2.5 
      delta_raw <- Delta * (1 - exp(-a * t)) / (1 - exp(-a))
      delta_true <- scale2(delta_raw, Delta, sd_expt)  
  } else if (func_form == "partial_convex") {
      p <- 0.8
      delta_raw <- ifelse(e <= Q/2,
                          p * Delta * (e - 1) / ((Q/2) - 1), 
                          p * Delta + (Delta - p * Delta) * (e - Q/2) / (Q/2))
      delta_true <- scale2(delta_raw, Delta, sd_expt)
  } else if (func_form == "const") {
      delta_true <- rep(Delta, Q)
  }
  Z_ETI <- build_Z_ETI(I, J)
  Z <- create_Z(Q, J)
  D2 <- create_D2(Z)
  
  vals_hat <- numeric(Q)
  for (s in 1:Q) {
    vals_hat[s] <- E_hat_delta_s(
      J = J, Q = Q, s = s,
      mu = mu, beta = beta, gamma = gamma,
      delta = get_delta_s(delta_true, s),
      x = x, y = y,
      Z = Z, D2 = D2
    )
  }
  Delta_ETI <- mean(vals_hat)
  
  return(Delta_ETI)
}