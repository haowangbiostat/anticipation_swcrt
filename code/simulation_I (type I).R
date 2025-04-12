source('utils_simulation_I.R')

# set up
nclus = 32 # number of cluster
nperiod = 9 # number of period
nsubject = 100 # number of individuals per period
alpha = 0.141 # cluster-level heterogeneity
sigma = 1 # individual error
delta = 0 # treatment effect
gamma = 0 # anticipation effect

n_sim <- 2000
set.seed(123)
initial_seed <- sample(1:10000, 1)
seeds <- initial_seed:(initial_seed + 2*n_sim - 1)

temp <- data.frame()

treatment_effects <- matrix(NA, nrow = n_sim, ncol = 2)
anticipation_effects <- matrix(NA, nrow = n_sim, ncol = 2)
naive_ses_t <- matrix(NA, nrow = n_sim, ncol = 2)
naive_ses_a <- matrix(NA, nrow = n_sim, ncol = 2)
naive_coverage_t <- matrix(NA, nrow = n_sim, ncol = 2)
naive_coverage_a <- matrix(NA, nrow = n_sim, ncol = 2)
est_alphas <- matrix(NA, nrow = n_sim, ncol = 2)
convergence <- matrix(NA, nrow = n_sim, ncol = 2)

colnames(treatment_effects) <- c("HH", "HH-ANT")
colnames(anticipation_effects) <- c("HH", "HH-ANT")
colnames(naive_ses_t) <- c("HH", "HH-ANT")
colnames(naive_ses_a) <- c("HH", "HH-ANT")
colnames(naive_coverage_t) <- c("HH", "HH-ANT")
colnames(naive_coverage_a) <- c("HH", "HH-ANT")
colnames(est_alphas) <- c("HH", "HH-ANT")
colnames(convergence) <- c("HH", "HH-ANT")

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
    
    # model-based se
    naive_ses_t[sim, ] <- s1_get_nve_t(fit_HH, fit_HH_ANT)
    naive_ses_a[sim, ] <- s1_get_nve_a(fit_HH, fit_HH_ANT)
    
    # est of treatment effects and anticipation effects
    treatment_effects[sim, ] <- s1_get_est_t(fit_HH, fit_HH_ANT)
    anticipation_effects[sim, ] <- s1_get_est_a(fit_HH, fit_HH_ANT)
    
    # est of cluster-level random effects
    est_alphas[sim, ] <- s1_get_alpha(fit_HH, fit_HH_ANT)
    
    # convergence
    convergence[sim, ] <- s1_get_convergence(fit_HH, fit_HH_ANT)
    
    
    for (model_idx in 1:2) {
      naive_lower_t = treatment_effects[sim, model_idx] - qnorm(0.975) * naive_ses_t[sim, model_idx]
      naive_upper_t = treatment_effects[sim, model_idx] + qnorm(0.975) * naive_ses_t[sim, model_idx]
      naive_coverage_t[sim, model_idx] <- ifelse(delta > naive_lower_t & delta < naive_upper_t, 1, 0) * 100
      
      naive_lower_a = anticipation_effects[sim, model_idx] - qnorm(0.975) * naive_ses_a[sim, model_idx]
      naive_upper_a = anticipation_effects[sim, model_idx] + qnorm(0.975) * naive_ses_a[sim, model_idx]
      naive_coverage_a[sim, model_idx] <- ifelse(gamma > naive_lower_a & gamma < naive_upper_a, 1, 0) * 100
    }
    
    for (model_type in c("HH", "HH-ANT")) {
      temp <- rbind(temp, data.frame(
        simulation = sim,
        model = model_type,
        treatment_effect = treatment_effects[sim, model_type],
        anticipation_effect = anticipation_effects[sim, model_type],
        naive_se_t = naive_ses_t[sim, model_type],
        naive_se_a = naive_ses_a[sim, model_type],
        naive_coverage_t = naive_coverage_t[sim, model_type],
        naive_coverage_a = naive_coverage_a[sim, model_type],
        alpha_sd = est_alphas[sim, model_type],
        convergence = convergence[sim, model_type]
      ))
    }
    
    if (sim %% batch_size == 0) {
      file_name <- sprintf("../result/simulation_I/simulation_results_%04d_to_%04d.csv", sim - batch_size + 1, sim)
      write.csv(temp, file_name, row.names = FALSE)
      temp <- data.frame(simulation = integer(),
                         model = character(),
                         treatment_effect = numeric(),
                         anticipation_effect = numeric(),
                         naive_se_t = numeric(),
                         naive_se_a = numeric(),
                         naive_coverage_t = numeric(),
                         naive_coverage_a = numeric(),
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