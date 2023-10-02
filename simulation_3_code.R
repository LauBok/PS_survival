source("Simulation util.R")
source("gen_datasets.R")

# Simulation 3

sim3_ER_1_result <- simulation(
  data_3[[1]], truth_ER = TRUE, fitted_ER = TRUE, misspecification = FALSE,
  fitted_type = 'Cox', estimand = 'prob'
)

sim3_ER_2_result <- simulation(
  data_3[[2]], truth_ER = TRUE, fitted_ER = TRUE, misspecification = FALSE,
  fitted_type = 'Cox', estimand = 'prob'
)

sim3_ER_3_result <- simulation(
  data_3[[3]], truth_ER = TRUE, fitted_ER = TRUE, misspecification = FALSE,
  fitted_type = 'Cox', estimand = 'prob'
)

sim3_ER_4_result <- simulation(
  data_3[[4]], truth_ER = TRUE, fitted_ER = TRUE, misspecification = FALSE,
  fitted_type = 'Cox', estimand = 'prob'
)

sim3_ER_5_result <- simulation(
  data_3[[5]], truth_ER = TRUE, fitted_ER = TRUE, misspecification = FALSE,
  fitted_type = 'Cox', estimand = 'prob'
)


sim3_df <- data.frame(
  true = numeric(0), mean = numeric(0), lwr = numeric(0), upr = numeric(0), time = numeric(0), 
  ER_type = character(0), strata_treatment = character(0)
)
sim3_df <- sim3_df %>% 
  add_to(sim3_ER_1_result, "Complier: 78%") %>%
  add_to(sim3_ER_2_result, "Complier: 60%") %>%
  add_to(sim3_ER_3_result, "Complier: 37%") %>%
  add_to(sim3_ER_4_result, "Complier: 21%") %>%
  add_to(sim3_ER_5_result, "Complier: 9%") %>%
  mutate(ER_type = factor(ER_type,
                          levels = c(
                            "Complier: 78%", "Complier: 60%", "Complier: 37%",
                            "Complier: 21%", "Complier: 9%"
                          )))
sim3_df %>% ggplot() + 
  geom_line(aes(time, true), linetype = "dashed") + 
  geom_line(aes(time, mean), linetype = "solid") + 
  geom_ribbon(aes(x = time, ymin = lwr, ymax = upr), alpha = 0.3) + 
  facet_grid(ER_type ~ strata_treatment) + 
  theme_bw() + 
  xlab("Time") + ylab("Survival Probability")