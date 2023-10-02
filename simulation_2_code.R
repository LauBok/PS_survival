source("Simulation util.R")
source("gen_datasets.R")

# Simulation 2

sim2_misS_misY_result <- simulation(
  data_1, truth_ER = TRUE, fitted_ER = TRUE, misspecification = TRUE,
  fitted_type = 'Cox', estimand = 'prob'
)

sim2_AFT_result <- simulation(
  data_1, truth_ER = TRUE, fitted_ER = TRUE, misspecification = FALSE,
  fitted_type = 'AFT', estimand = 'prob'
)

sim2_misS_misY_nonER_result <- simulation(
  data_1, truth_ER = TRUE, fitted_ER = FALSE, misspecification = TRUE,
  fitted_type = 'Cox', estimand = 'prob'
)

sim2_AFT_nonER_result <- simulation(
  data_1, truth_ER = TRUE, fitted_ER = FALSE, misspecification = FALSE,
  fitted_type = 'AFT', estimand = 'prob'
)

# -------------------- PLOT ------------------------

sim2_truth_ER_df <- data.frame(
  true = numeric(0), mean = numeric(0), lwr = numeric(0), upr = numeric(0), time = numeric(0), 
  ER_type = character(0), strata_treatment = character(0)
)

sim2_truth_ER_df <- sim2_truth_ER_df %>% 
  add_to(sim2_misS_misY_result, "S: mis, T: mis (ER)") %>%
  add_to(sim2_AFT_result, "S: true, T: AFT (ER)") %>%
  add_to(sim2_misS_misY_nonER_result, "S: mis, T: mis (No ER)") %>%
  add_to(sim2_AFT_nonER_result, "S: true, T: AFT (No ER)") %>%
  mutate(ER_type = factor(ER_type,
                          levels = c(
                            "S: mis, T: mis (ER)", "S: true, T: AFT (ER)",
                            "S: mis, T: mis (No ER)", "S: true, T: AFT (No ER)"
                          )))
sim2_truth_ER_df %>% ggplot() + 
  geom_line(aes(time, true), linetype = "dashed") + 
  geom_line(aes(time, mean), linetype = "solid") + 
  geom_ribbon(aes(x = time, ymin = lwr, ymax = upr), alpha = 0.3) + 
  facet_grid(ER_type ~ strata_treatment) + 
  theme_bw() + 
  xlab("Time") + ylab("Survival Probability")

# ------------------------- For supplement ---------------------------

sim2_non_ER_misS_misY_result <- simulation(
  data_2, truth_ER = FALSE, fitted_ER = TRUE, misspecification = TRUE,
  fitted_type = 'Cox', estimand = 'prob'
)

sim2_non_ER_AFT_result <- simulation(
  data_2, truth_ER = FALSE, fitted_ER = TRUE, misspecification = FALSE,
  fitted_type = 'AFT', estimand = 'prob'
)

sim2_non_ER_misS_misY_nonER_result <- simulation(
  data_2, truth_ER = FALSE, fitted_ER = FALSE, misspecification = TRUE,
  fitted_type = 'Cox', estimand = 'prob'
)

sim2_non_ER_AFT_nonER_result <- simulation(
  data_2, truth_ER = FALSE, fitted_ER = FALSE, misspecification = FALSE,
  fitted_type = 'AFT', estimand = 'prob'
)

# ---------------------------- PLOT --------------------------------

sim2_truth_non_ER_df <- data.frame(
  true = numeric(0), mean = numeric(0), lwr = numeric(0), upr = numeric(0), time = numeric(0), 
  ER_type = character(0), strata_treatment = character(0)
)

sim2_truth_non_ER_df <- sim2_truth_non_ER_df %>% 
  add_to(sim2_non_ER_misS_misY_result, "S: mis, T: mis (ER)") %>%
  add_to(sim2_non_ER_AFT_result, "S: true, T: AFT (ER)") %>%
  add_to(sim2_non_ER_misS_misY_nonER_result, "S: mis, T: mis (No ER)") %>%
  add_to(sim2_non_ER_AFT_nonER_result, "S: true, T: AFT (No ER)") %>%
  mutate(ER_type = factor(ER_type,
                          levels = c(
                            "S: mis, T: mis (ER)", "S: true, T: AFT (ER)",
                            "S: mis, T: mis (No ER)", "S: true, T: AFT (No ER)"
                          )))
sim2_truth_non_ER_df %>% ggplot() + 
  geom_line(aes(time, true), linetype = "dashed") + 
  geom_line(aes(time, mean), linetype = "solid") + 
  geom_ribbon(aes(x = time, ymin = lwr, ymax = upr), alpha = 0.3) + 
  facet_grid(ER_type ~ strata_treatment) + 
  theme_bw() + 
  xlab("Time") + ylab("Survival Probability")