library(dplyr)
library(xtable)
library(tidyverse)

cont_all <- list.files(pattern = "simulation_results", full.names = TRUE) %>%
  lapply(read_csv) %>%
  bind_rows()
write.csv(cont_all, "cont_all.csv", row.names = FALSE)

df <- read.csv("cont_all.csv")

result <- df %>%
  group_by(model = factor(model, levels = c("HH", "HH-ANT", "ETI", "ETI-ANT", "TEH", "TEH-ANT"))) %>%
  summarise(
    est_t = mean(treatment_effect),
    est_a = mean(anticipation_effect),
    alpha_sd = mean(alpha_sd),
    empirical_t = sd(treatment_effect),
    naive_se_t = mean(naive_se_t),
    naive_coverage_t = mean(naive_coverage_t),
    naive_power_t = 100 - mean(naive_power_t),
    empirical_a = sd(anticipation_effect),
    naive_se_a = mean(naive_se_a),
    naive_coverage_a = mean(naive_coverage_a),
    naive_power_a = 100-mean(naive_power_a),
    convergence_rate = mean(convergence) * 100
  ) %>%
  arrange(model)

print(xtable(result, digits = c(0, 4, 4, 4, 4, 4, 4, 2, 2, 4, 4, 2, 2, 1)), include.rownames = FALSE)