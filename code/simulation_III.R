source('utils_simulation_III.R')

# set up
nclus = 32 # number of cluster
nperiod = 9 # number of period
nsubject = 100 # number of individuals per period
alpha = 0.141 # cluster-level heterogeneity
sigma = 1 # individual error
delta = 0.12 # treatment effect
gamma = 0 # anticipation effect

n_sim <- 2000
set.seed(123)
initial_seed <- sample(1:10000, 1)
seeds <- initial_seed:(initial_seed + 2*n_sim - 1)

temp <- data.frame()

treatment_effects <- matrix(NA, nrow = n_sim, ncol = 6)
anticipation_effects <- matrix(NA, nrow = n_sim, ncol = 6)
naive_ses_t <- matrix(NA, nrow = n_sim, ncol = 6)
naive_ses_a <- matrix(NA, nrow = n_sim, ncol = 6)
naive_coverage_t <- matrix(NA, nrow = n_sim, ncol = 6)
naive_power_t <- matrix(NA, nrow = n_sim, ncol = 6)
naive_coverage_a <- matrix(NA, nrow = n_sim, ncol = 6)
naive_power_a <- matrix(NA, nrow = n_sim, ncol = 6)
est_alphas <- matrix(NA, nrow = n_sim, ncol = 6)
est_sigma_deltas <- matrix(NA, nrow = n_sim, ncol = 6)
convergence <- matrix(NA, nrow = n_sim, ncol = 6)

colnames(treatment_effects) <- c("HH", "HH-ANT", "ETI", "ETI-ANT", "TEH", "TEH-ANT")
colnames(anticipation_effects) <- c("HH", "HH-ANT", "ETI", "ETI-ANT", "TEH", "TEH-ANT")
colnames(naive_ses_t) <- c("HH", "HH-ANT", "ETI", "ETI-ANT", "TEH", "TEH-ANT")
colnames(naive_ses_a) <- c("HH", "HH-ANT", "ETI", "ETI-ANT", "TEH", "TEH-ANT")
colnames(naive_coverage_t) <- c("HH", "HH-ANT", "ETI", "ETI-ANT", "TEH", "TEH-ANT")
colnames(naive_power_t) <- c("HH", "HH-ANT", "ETI", "ETI-ANT", "TEH", "TEH-ANT")
colnames(naive_coverage_a) <- c("HH", "HH-ANT", "ETI", "ETI-ANT", "TEH", "TEH-ANT")
colnames(naive_power_a) <- c("HH", "HH-ANT", "ETI", "ETI-ANT", "TEH", "TEH-ANT")
colnames(est_alphas) <- c("HH", "HH-ANT", "ETI", "ETI-ANT", "TEH", "TEH-ANT")
colnames(est_sigma_deltas) <- c("HH", "HH-ANT", "ETI", "ETI-ANT", "TEH", "TEH-ANT")
colnames(convergence) <- c("HH", "HH-ANT", "ETI", "ETI-ANT", "TEH", "TEH-ANT")

sim_count <- 0
seed_count <- 1
batch_size <- 200


while (sim_count < n_sim) {
  set.seed(seeds[seed_count])
  sim <- sim_count + 1
  print(paste("Running simulation", sim, "with seed", seeds[seed_count]))
  
  result_sim <- try({
    data <- gendata(nclus = nclus,
                    nperiod = nperiod,
                    nsubject = nsubject,
                    alpha = alpha,
                    sigma = sigma,
                    delta = delta,
                    gamma = gamma
    )
    
    # fit HH model
    fit_HH <- lmerTest::lmer(response.var ~  tx.var + time.var + (1|cluster.var), data = data)
    
    # fit HH-ANT model
    fit_HH_ANT <- lmerTest::lmer(response.var ~  tx.var + ant.var + time.var + (1|cluster.var), data = data)
    
    # fit ETI model
    fit_ETI <- lmerTest::lmer(response.var ~  expt.var + time.var + (1|cluster.var), data = data)
    
    # fit ETI-ANT model
    fit_ETI_ANT <- lmerTest::lmer(response.var ~  expt.var + ant.var + time.var + (1|cluster.var), data = data)
    
    # fit TEH
    fit_TEH <- lmerTest::lmer(response.var ~  tx.var + time.var + (0 + tx.var | expt.var) + (1|cluster.var), data = data)
    
    # fit TEH-ANT
    fit_TEH_ANT <- lmerTest::lmer(response.var ~  tx.var + ant.var + time.var + (0 + tx.var | expt.var) + (1|cluster.var), data = data)
    
    # model-based se
    naive_ses_t[sim, ] <- s2_get_nve_t(fit_HH, fit_HH_ANT, fit_ETI, fit_ETI_ANT, fit_TEH, fit_TEH_ANT, nperiod)
    naive_ses_a[sim, ] <- s2_get_nve_a(fit_HH, fit_HH_ANT, fit_ETI, fit_ETI_ANT, fit_TEH, fit_TEH_ANT, nperiod)
    
    # est of treatment effects and anticipation effects
    treatment_effects[sim, ] <- s2_get_est_t(fit_HH, fit_HH_ANT, fit_ETI, fit_ETI_ANT, fit_TEH, fit_TEH_ANT, nperiod)
    anticipation_effects[sim, ] <- s2_get_est_a(fit_HH, fit_HH_ANT, fit_ETI, fit_ETI_ANT, fit_TEH, fit_TEH_ANT, nperiod)
    
    # est of cluster-level random effects
    est_alphas[sim, ] <- s2_get_alpha(fit_HH, fit_HH_ANT, fit_ETI, fit_ETI_ANT, fit_TEH, fit_TEH_ANT)
    
    # est of sigma delta
    est_sigma_deltas[sim, ] <- s2_get_sigma_delta(fit_HH, fit_HH_ANT, fit_ETI, fit_ETI_ANT, fit_TEH, fit_TEH_ANT)
    
    # convergence
    convergence[sim, ] <- s2_get_convergence(fit_HH, fit_HH_ANT, fit_ETI, fit_ETI_ANT, fit_TEH, fit_TEH_ANT)
    
    
    for (model_idx in 1:6) {
      naive_lower_t = treatment_effects[sim, model_idx] - qnorm(0.975) * naive_ses_t[sim, model_idx]
      naive_upper_t = treatment_effects[sim, model_idx] + qnorm(0.975) * naive_ses_t[sim, model_idx]
      naive_coverage_t[sim, model_idx] <- ifelse(delta > naive_lower_t & delta < naive_upper_t, 1, 0) * 100
      naive_power_t[sim, model_idx] <- ifelse(0 > naive_lower_t & 0 < naive_upper_t, 1, 0) * 100
      
      naive_lower_a = anticipation_effects[sim, model_idx] - qnorm(0.975) * naive_ses_a[sim, model_idx]
      naive_upper_a = anticipation_effects[sim, model_idx] + qnorm(0.975) * naive_ses_a[sim, model_idx]
      naive_coverage_a[sim, model_idx] <- ifelse(gamma > naive_lower_a & gamma < naive_upper_a, 1, 0) * 100
      naive_power_a[sim, model_idx] <- ifelse(0 > naive_lower_a & 0 < naive_upper_a, 1, 0) * 100
    }
    
    for (model_type in c("HH", "HH-ANT", "ETI", "ETI-ANT", "TEH", "TEH-ANT")) {
      temp <- rbind(temp, data.frame(
        simulation = sim,
        model = model_type,
        treatment_effect = treatment_effects[sim, model_type],
        anticipation_effect = anticipation_effects[sim, model_type],
        naive_se_t = naive_ses_t[sim, model_type],
        naive_se_a = naive_ses_a[sim, model_type],
        naive_coverage_t = naive_coverage_t[sim, model_type],
        naive_power_t = naive_power_t[sim, model_type],
        naive_coverage_a = naive_coverage_a[sim, model_type],
        naive_power_a = naive_power_a[sim, model_type],
        alpha_sd = est_alphas[sim, model_type],
        delta_sd = est_sigma_deltas[sim, model_type],
        convergence = convergence[sim, model_type]
      ))
    }
    
    if (sim %% batch_size == 0) {
      file_name <- sprintf("simulation_results_%04d_to_%04d.csv", sim - batch_size + 1, sim)
      write.csv(temp, file_name, row.names = FALSE)
      temp <- data.frame(simulation = integer(),
                         model = character(),
                         treatment_effect = numeric(),
                         anticipation_effect = numeric(),
                         naive_se_t = numeric(),
                         naive_se_a = numeric(),
                         naive_coverage_t = numeric(),
                         naive_power_t = numeric(),
                         naive_coverage_a = numeric(),
                         naive_power_a = numeric(),
                         alpha_sd = numeric(),
                         convergence = logical(),
                         stringsAsFactors = FALSE)
    }
    
    sim_count <- sim_count + 1
    seed_count <- seed_count + 1
    
  }, silent = FALSE)
  
  if (inherits(result_sim, "try-error")) {
    cat("An error occurred in simulation", sim, "with seed", seeds[seed_count], ". Skipping to the next iteration.\n")
    seed_count <- seed_count + 1
  }
  
}