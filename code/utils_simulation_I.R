library(tidyverse)
library(MASS)
library(lmerTest)

gendata <- function(
    nclus,    # number of clusters
    nperiod,  # number of periods
    nsubject, # number of individuals per cluster-period
    alpha,    # cluster-level heterogeneity
    sigma,    # individual error standard deviation
    delta,    # treatment effect
    gamma     # anticipation effect
){
  nseq = nperiod - 1
  nclus_seq = nclus / nseq
  
  # treatment effect indicator
  X <- matrix(0, nrow = nseq, ncol = nperiod)
  for(i in 1:nseq){
    X[i, (i + 1):nperiod] <- 1
  }
  X <- X %x% rep(1, nclus_seq)
  
  # anticipation effect indicator
  A <- matrix(0, nrow = nseq, ncol = nperiod)
  for(i in 1:nseq){
    A[i, i] <- 1
  }
  A <- A %x% rep(1, nclus_seq)
  
  # cluster random effects
  C <- rnorm(n = nclus, mean = 0, sd = alpha)
  
  # outcome vector
  Yijk <- numeric(nclus * nperiod * nsubject)
  
  # individual errors
  e <- rnorm(n = length(Yijk), mean = 0, sd = sigma)
  
  for (i in 1:nclus) {
    for (j in 1:nperiod) {
      for (k in 1:nsubject) {
        # index for the outcome vector
        l <- (i - 1) * nperiod * nsubject + (j - 1) * nsubject + k
        Yijk[l] <- j + X[i, j] * delta + A[i, j] * gamma + C[i] + e[l]
      }
    }
  }
  
  Xvec <- rep(as.vector(t(X)), each = nsubject)
  Avec <- rep(as.vector(t(A)), each = nsubject)
  Cvec <- rep(1:nclus, each = nperiod * nsubject)
  Pvec <- rep(rep(1:nperiod, each = nsubject), times = nclus)
  
  dat <- data.frame(
    response.var = Yijk,
    tx.var = Xvec,
    ant.var = Avec,
    cluster.var = Cvec,
    time.var = Pvec
  )
  
  first_time_assignment <- dat %>%
    filter(tx.var == 1) %>%
    group_by(cluster.var) %>%
    summarize(first_time = min(time.var))
  
  dat <- dat %>%
    left_join(first_time_assignment, by = "cluster.var") %>%
    mutate(expt.var = ifelse(tx.var == 0, 0, time.var - first_time + 1)) %>%
    dplyr::select(-first_time)
  
  dat$cluster.var <- as.factor(dat$cluster.var)
  dat$time.var <- as.factor(dat$time.var)
  dat$expt.var <- as.factor(dat$expt.var)
  
  return(dat)
}

# naive variance estimators of treatment effect
s1_get_nve_t <- function(m1, m2){
  nai_HH <- sqrt(vcov(m1)[2, 2])
  nai_HH_ANT <- sqrt(vcov(m2)[2, 2])
  return(c(nai_HH, nai_HH_ANT))
}

# naive variance estimators of anticipation effect
s1_get_nve_a <- function(m1, m2){
  nai_HH <- NA
  nai_HH_ANT <- sqrt(vcov(m2)[3, 3])
  return(c(nai_HH, nai_HH_ANT))
}

# estimate treatment effect
s1_get_est_t <- function(m1, m2){
  est_HH <- unname(fixef(m1)[2])
  est_HH_ANT <- unname(fixef(m2)[2])
  return(c(est_HH, est_HH_ANT))
}

# estimate anticipation effect
s1_get_est_a <- function(m1, m2){
  est_HH <- NA
  est_HH_ANT <- unname(fixef(m2)[3])
  return(c(est_HH, est_HH_ANT))
}

# estimate cluster random effect
s1_get_alpha <- function(m1, m2){
  est_HH <- as.data.frame(VarCorr(m1))[1,5]
  est_HH_ANT <- as.data.frame(VarCorr(m2))[1,5]
  return(c(est_HH, est_HH_ANT))
}

# check convergence
s1_get_convergence <- function(m1, m2){
  con_HH <- ifelse(is.null(m1@optinfo$conv$lme4$code), 1, 0)
  con_HH_ANT <- ifelse(is.null(m2@optinfo$conv$lme4$code), 1, 0)
  return(c(con_HH, con_HH_ANT))
}

# naive variance estimators of treatment effect
s2_get_nve_t <- function(m1, m2, m3, m4, m5, m6, nperiod){
  nai_HH <- sqrt(vcov(m1)[2, 2])
  nai_HH_ANT <- sqrt(vcov(m2)[2, 2])
  
  C <- rep(1,(nperiod-1))/(nperiod-1)
  vcov_params_ETI <- vcov(m3)[2:nperiod, 2:nperiod]
  nai_ETI <- as.numeric(sqrt(C %*% vcov_params_ETI %*% C))
  vcov_params_ETI_ANT <- vcov(m4)[2:nperiod, 2:nperiod]
  nai_ETI_ANT <- as.numeric(sqrt(C %*% vcov_params_ETI_ANT %*% C))
  
  
  nai_TEH <- sqrt(vcov(m5)[2, 2])
  nai_TEH_ANT <- sqrt(vcov(m6)[2, 2])
  return(c(nai_HH, nai_HH_ANT, nai_ETI, nai_ETI_ANT, nai_TEH, nai_TEH_ANT))
}

# naive variance estimators of anticipation effect
s2_get_nve_a <- function(m1, m2, m3, m4, m5, m6, nperiod){
  nai_HH <- NA
  nai_HH_ANT <- sqrt(vcov(m2)[3, 3])
  nai_ETI <- NA
  nai_ETI_ANT <- sqrt(vcov(m4)[3, 3])
  nai_TEH <- NA
  nai_TEH_ANT <- sqrt(vcov(m6)[3, 3])
  return(c(nai_HH, nai_HH_ANT, nai_ETI, nai_ETI_ANT, nai_TEH, nai_TEH_ANT))
}

# estimate treatment effect
s2_get_est_t <- function(m1, m2, m3, m4, m5, m6, nperiod){
  est_HH <- unname(fixef(m1)[2])
  est_HH_ANT <- unname(fixef(m2)[2])
  est_ETI <- mean(fixef(m3)[2:nperiod])
  est_ETI_ANT <- mean(fixef(m4)[2:nperiod])
  est_TEH <- unname(fixef(m5)[2])
  est_TEH_ANT <- unname(fixef(m6)[2])
  return(c(est_HH, est_HH_ANT, est_ETI, est_ETI_ANT, est_TEH, est_TEH_ANT))
}

# estimate anticipation effect
s2_get_est_a <- function(m1, m2, m3, m4, m5, m6, nperiod){
  est_HH <- NA
  est_HH_ANT <- unname(fixef(m2)[3])
  est_ETI <- NA
  est_ETI_ANT <- unname(fixef(m4)[nperiod+1])
  est_TEH <- NA
  est_TEH_ANT <- unname(fixef(m6)[3])
  return(c(est_HH, est_HH_ANT, est_ETI, est_ETI_ANT, est_TEH, est_TEH_ANT))
}

# estimate cluster random effect
s2_get_alpha <- function(m1, m2, m3, m4, m5, m6){
  est_HH <- as.data.frame(VarCorr(m1))[1,5]
  est_HH_ANT <- as.data.frame(VarCorr(m2))[1,5]
  est_ETI <- as.data.frame(VarCorr(m3))[1,5]
  est_ETI_ANT <- as.data.frame(VarCorr(m4))[1,5]
  
  df_varcorr_TEH = as.data.frame(VarCorr(m5))
  est_TEH <- df_varcorr_TEH[df_varcorr_TEH$grp == "cluster.var",]$sdcor
  
  df_varcorr_TEH_ANT = as.data.frame(VarCorr(m6))
  est_TEH_ANT <- df_varcorr_TEH_ANT[df_varcorr_TEH_ANT$grp == "cluster.var",]$sdcor
  return(c(est_HH, est_HH_ANT, est_ETI, est_ETI_ANT, est_TEH, est_TEH_ANT))
}

# estimate random treatment effect
s2_get_sigma_delta <- function(m1, m2, m3, m4, m5, m6){
  est_HH <- NA
  est_HH_ANT <- NA
  est_ETI <- NA
  est_ETI_ANT <- NA
  
  df_varcorr_TEH = as.data.frame(VarCorr(m5))
  est_TEH <- df_varcorr_TEH[df_varcorr_TEH$grp == "expt.var",]$sdcor
  
  df_varcorr_TEH_ANT = as.data.frame(VarCorr(m6))
  est_TEH_ANT <- df_varcorr_TEH_ANT[df_varcorr_TEH_ANT$grp == "expt.var",]$sdcor
  return(c(est_HH, est_HH_ANT, est_ETI, est_ETI_ANT, est_TEH, est_TEH_ANT))
}

# check convergence
s2_get_convergence <- function(m1, m2, m3, m4, m5, m6){
  con_HH <- ifelse(is.null(m1@optinfo$conv$lme4$code), 1, 0)
  con_HH_ANT <- ifelse(is.null(m2@optinfo$conv$lme4$code), 1, 0)
  con_ETI <- ifelse(is.null(m3@optinfo$conv$lme4$code), 1, 0)
  con_ETI_ANT <- ifelse(is.null(m4@optinfo$conv$lme4$code), 1, 0)
  con_TEH <- ifelse(is.null(m5@optinfo$conv$lme4$code), 1, 0)
  con_TEH_ANT <- ifelse(is.null(m6@optinfo$conv$lme4$code), 1, 0)
  return(c(con_HH, con_HH_ANT, con_ETI, con_ETI_ANT, con_TEH, con_TEH_ANT))
}