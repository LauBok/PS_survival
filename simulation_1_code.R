source("Simulation util.R")
source("gen_datasets.R")

# Simulation 1
truth_ER_fit_ER_result <- simulation(
  data_1, truth_ER = TRUE, fitted_ER = TRUE, misspecification = FALSE,
  fitted_type = 'Cox', estimand = 'prob'
)

truth_ER_fit_non_ER_result <- simulation(
  data_1, truth_ER = TRUE, fitted_ER = FALSE, misspecification = FALSE,
  fitted_type = 'Cox', estimand = 'prob'
)

truth_non_ER_fit_ER_result <- simulation(
  data_2, truth_ER = FALSE, fitted_ER = TRUE, misspecification = FALSE,
  fitted_type = 'Cox', estimand = 'prob'
)

truth_non_ER_fit_non_ER_result <- simulation(
  data_2, truth_ER = FALSE, fitted_ER = FALSE, misspecification = FALSE,
  fitted_type = 'Cox', estimand = 'prob'
)

# -------------------------- PLOT ------------------------------

sim1_df <- data.frame(
  true = numeric(0), mean = numeric(0), lwr = numeric(0), upr = numeric(0), time = numeric(0), 
  ER_type = character(0), strata_treatment = character(0)
)

sim1_df <- sim1_df %>% 
  add_to(truth_ER_fit_ER_result, "Truth: ER, Fit: ER") %>%
  add_to(truth_ER_fit_non_ER_result, "Truth: ER, Fit: Non-ER") %>%
  add_to(truth_non_ER_fit_ER_result, "Truth: Non-ER, Fit: ER") %>%
  add_to(truth_non_ER_fit_non_ER_result, "Truth: Non-ER, Fit: Non-ER")

sim1_df %>% ggplot() + 
  geom_line(aes(time, true), linetype = "dashed") + 
  geom_line(aes(time, mean), linetype = "solid") + 
  geom_ribbon(aes(x = time, ymin = lwr, ymax = upr), alpha = 0.3) + 
  facet_grid(ER_type ~ strata_treatment) + 
  theme_bw() + 
  xlab("Time") + ylab("Survival Probability")

# -------------------- For Supplement ----------------------------

truth_ER_fit_ER_RACE_result <- simulation(
  data_1, truth_ER = TRUE, fitted_ER = TRUE, misspecification = FALSE,
  fitted_type = 'Cox', estimand = 'RACE'
)


truth_ER_fit_non_ER_RACE_result <- simulation(
  data_1, truth_ER = TRUE, fitted_ER = FALSE, misspecification = FALSE,
  fitted_type = 'Cox', estimand = 'RACE'
)

truth_non_ER_fit_ER_RACE_result <- simulation(
  data_2, truth_ER = FALSE, fitted_ER = TRUE, misspecification = FALSE,
  fitted_type = 'Cox', estimand = 'RACE'
)

truth_non_ER_fit_non_ER_RACE_result <- simulation(
  data_2, truth_ER = FALSE, fitted_ER = FALSE, misspecification = FALSE,
  fitted_type = 'Cox', estimand = 'RACE'
)

# ----------------------- PLOT ---------------------------

sim1_RACE_df <- data.frame(
  true = numeric(0), mean = numeric(0), lwr = numeric(0), upr = numeric(0), time = numeric(0), 
  ER_type = character(0), strata_treatment = character(0)
)
sim1_RACE_df <- sim1_RACE_df %>% 
  add_to(truth_ER_fit_ER_RACE_result, "Truth: ER, Fit: ER") %>%
  add_to(truth_ER_fit_non_ER_RACE_result, "Truth: ER, Fit: Non-ER") %>%
  add_to(truth_non_ER_fit_ER_result, "Truth: Non-ER, Fit: ER") %>%
  add_to(truth_non_ER_fit_non_ER_RACE_result, "Truth: Non-ER, Fit: Non-ER")
sim1_RACE_df %>% ggplot() + 
  geom_line(aes(time, true), linetype = "dashed") + 
  geom_line(aes(time, mean), linetype = "solid") + 
  geom_ribbon(aes(x = time, ymin = lwr, ymax = upr), alpha = 0.3) + 
  facet_grid(ER_type ~ strata_treatment) + 
  theme_bw() + 
  xlab("Time") + ylab("RACE")    