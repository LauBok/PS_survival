## Util functions

library(devtools)
library(tidyverse)
library(pbapply)
# devtools::install_github("LauBok/PStrata")
library(PStrata) 

# Sample -----------------------------
## Calculate stratum probability
calc_stratum_probability <- function(strata_coef_matrix, data_matrix) {
  t(apply(
    exp(data_matrix %*% strata_coef_matrix), 1, function(x) x / sum(x)
  ))
}

## Calculate survival probability
calc_surv_probability_AFT <- function(outcome_coef_matrix, sigma_vector, data_matrix, time_sequence) {
  mu <- data_matrix %*% outcome_coef_matrix
  lapply(1:ncol(mu), function(j) pnorm((rep(1, nrow(mu)) %*% t(log(time_sequence)) - mu[, j]) / sigma_vector[j], lower.tail = FALSE))
}

calc_surv_probability <- function(outcome_coef_matrix, phi_vector, data_matrix, time_sequence) {
  exp.mu <- exp(data_matrix %*% outcome_coef_matrix)
  lapply(1:ncol(exp.mu), function(j) exp(- (exp.mu[, j]) %*% t(time_sequence^phi_vector[j]) / phi_vector[j]))
}

calc_RACE <- function(outcome_coef_matrix, phi_vector, data_matrix, time_sequence) {
  exp.mu <- exp(data_matrix %*% outcome_coef_matrix)
  c = t(t(exp.mu) / phi_vector)
  lapply(1:ncol(exp.mu), function(j) (phi_vector[j] * c[, j]^(1/phi_vector[j]))^(-1) * pgamma(c[, j] %*% t(time_sequence^phi_vector[j]), 1/phi_vector[j]) * gamma(1/phi_vector[j]))
}

## generate survival time from Weibull Cox model
get_survival_cox <- function(n, mu, phi) {
  # inverse CDF sampling
  u <- runif(n, 0, 1)
  return ((-log(u) * phi * exp(-mu))^(1 / phi))
}


## generate survival time from normal AFT model
get_survival_AFT <- function(n, mu, sigma) {
  return (exp(rnorm(n, mu, sigma)))
}

## sample S from multinomial model
sample_S <- function(strata_coef_matrix, data_matrix_S, strata) {
  strata_prob <- calc_stratum_probability(strata_coef_matrix, data_matrix_S)
  apply(strata_prob, 1, function(x)
    sample(c(0, 1, 3), size = 1, replace = TRUE, prob = x)
  )
}

## sample Y from Weibull Cox model
sample_Y_cox <- function(outcome_coef_matrix, phi_vector, data_matrix_Y) {
  outcome_matrix <- data_matrix_Y %*% outcome_coef_matrix
  do.call(
    cbind, 
    map(
      1:length(phi_vector),
      ~get_survival_cox(nrow(data_matrix_Y), outcome_matrix[, .x], phi_vector[.x])
    )
  )
}

## sample Y from Weibull Cox model
sample_Y_AFT <- function(outcome_coef_matrix, sigma_vector, data_matrix_Y) {
  outcome_matrix <- data_matrix_Y %*% outcome_coef_matrix
  do.call(
    cbind, 
    map(
      1:length(sigma_vector),
      ~get_survival_AFT(nrow(data_matrix_Y), outcome_matrix[, .x], sigma_vector[.x])
    )
  )
}

# Analysis ---------------------------

stratum_surv_probability <- function(strata_coef_matrix, outcome_coef_matrix, phi_vector,
                                     data_matrix_S, data_matrix_G, time_sequence, strata_list) {
  stratum_prob <- calc_stratum_probability(strata_coef_matrix, data_matrix_S)
  surv_prob <- calc_surv_probability(outcome_coef_matrix, phi_vector, data_matrix_G, time_sequence)
  sapply(1:length(surv_prob), 
         function(j) t(surv_prob[[j]]) %*% stratum_prob[, strata_list[j]] / sum(stratum_prob[, strata_list[j]])
  )
}

stratum_surv_probability_AFT <- function(strata_coef_matrix, outcome_coef_matrix, sigma_vector,
                                     data_matrix_S, data_matrix_G, time_sequence, strata_list) {
  stratum_prob <- calc_stratum_probability(strata_coef_matrix, data_matrix_S)
  surv_prob <- calc_surv_probability_AFT(outcome_coef_matrix, sigma_vector, data_matrix_G, time_sequence)
  sapply(1:length(surv_prob), 
         function(j) t(surv_prob[[j]]) %*% stratum_prob[, strata_list[j]] / sum(stratum_prob[, strata_list[j]])
  )
}

stratum_RACE <- function(strata_coef_matrix, outcome_coef_matrix, phi_vector,
                                     data_matrix_S, data_matrix_G, time_sequence, strata_list) {
  stratum_prob <- calc_stratum_probability(strata_coef_matrix, data_matrix_S)
  RACE <- calc_RACE(outcome_coef_matrix, phi_vector, data_matrix_G, time_sequence)
  sapply(1:length(RACE), 
         function(j) t(RACE[[j]]) %*% stratum_prob[, strata_list[j]] / sum(stratum_prob[, strata_list[j]])
  )
}

## Reparameterization from (theta, beta_G) given in PStrata to (phi, alpha) in the paper.
## phi = exp(theta), alpha = beta_G + theta
get_surv_prob_matrices <- function(result, time_sequence, strata_list) {
  beta_S <- rstan::extract(result$post_samples, par = "beta_S")[[1]]
  beta_G <- rstan::extract(result$post_samples, par = "beta_G")[[1]]
  theta <- rstan::extract(result$post_samples, par = "theta")[[1]]
  tmp_adjust <- t(c(1, rep(0, ncol(beta_G[1, , ]) - 1)))
  data_matrix_S = result$PSobject$S.formula$fixed_eff_matrix
  data_matrix_G = result$PSobject$Y.formula$fixed_eff_matrix
  results <- pblapply(1:nrow(theta), function(i) {
    stratum_surv_probability(
      strata_coef_matrix = t(rbind(0, beta_S[i, ,])),
      outcome_coef_matrix = t(beta_G[i, ,] + theta[i, ] %*% tmp_adjust),
      phi_vector = exp(theta[i, ]),
      data_matrix_S = data_matrix_S,
      data_matrix_G = data_matrix_G,
      time_sequence = time_sequence,
      strata_list = strata_list
    )
  })
  return (results)
}

get_RACE_matrices <- function(result, time_sequence, strata_list) {
  beta_S <- rstan::extract(result$post_samples, par = "beta_S")[[1]]
  beta_G <- rstan::extract(result$post_samples, par = "beta_G")[[1]]
  theta <- rstan::extract(result$post_samples, par = "theta")[[1]]
  tmp_adjust <- t(c(1, rep(0, ncol(beta_G[1, , ]) - 1)))
  data_matrix_S = result$PSobject$S.formula$fixed_eff_matrix
  data_matrix_G = result$PSobject$Y.formula$fixed_eff_matrix
  results <- pblapply(1:nrow(theta), function(i) {
    stratum_RACE(
      strata_coef_matrix = t(rbind(0, beta_S[i, ,])),
      outcome_coef_matrix = t(beta_G[i, ,] + theta[i, ] %*% tmp_adjust),
      phi_vector = exp(theta[i, ]),
      data_matrix_S = data_matrix_S,
      data_matrix_G = data_matrix_G,
      time_sequence = time_sequence,
      strata_list = strata_list
    )
  })
  return (results)
}

## phi = exp(theta), alpha = beta_G + theta
get_surv_prob_matrices_AFT <- function(result, time_sequence, strata_list) {
  beta_S <- rstan::extract(result$post_samples, par = "beta_S")[[1]]
  beta_G <- rstan::extract(result$post_samples, par = "beta_G")[[1]]
  sigma <- rstan::extract(result$post_samples, par = "sigma")[[1]]
  data_matrix_S = result$PSobject$S.formula$fixed_eff_matrix
  data_matrix_G = result$PSobject$Y.formula$fixed_eff_matrix
  results <- pblapply(1:nrow(sigma), function(i) {
    stratum_surv_probability_AFT(
      strata_coef_matrix = t(rbind(0, beta_S[i, ,])),
      outcome_coef_matrix = t(beta_G[i, ,]),
      sigma_vector = sigma[i, ],
      data_matrix_S = data_matrix_S,
      data_matrix_G = data_matrix_G,
      time_sequence = time_sequence,
      strata_list = strata_list
    )
  })
  return (results)
}

summarize_surv_prob <- function(surv_prob_matrices) {
  num_group <- ncol(surv_prob_matrices[[1]])
  lapply(1:num_group, function(g) {
    apply(
      sapply(1:length(surv_prob_matrices), function(b) surv_prob_matrices[[b]][, g]),
      1,
      function(x) c(mean = mean(x[!is.nan(x)]), lwr = unname(quantile(x, 0.025, na.rm = T)), upr = unname(quantile(x, 0.975, na.rm = T)))
    )
  })
}

combine_fitted_true_surv_prob <- function(true_surv_prob, summary_surv_prob, 
                                          true_id, summary_id) {
  lapply(1:6, function(g) {
    true_surv_prob_g = true_surv_prob[, true_id[g]]
    summarize_surv_prob_g = summary_surv_prob[[summary_id[g]]]
    return (rbind(true = c(true_surv_prob_g), summarize_surv_prob_g))
  })
}

make_plot_data <- function(summary_data, time_sequence) {
  # syntax: summary_data[[group]][type, time]
  # output: data.frame(S, Z, time, true, mean, lwr, upr)
  S = c("never-taker", "never-taker", "complier", "complier", "always-taker", "always-taker")
  Z = c("control", "treatment", "control", "treatment", "control", "treatment")
  lapply(1:length(summary_data), function(g) {
    data.frame(S = S[g], Z = Z[g], time = time_sequence, 
               true = summary_data[[g]]["true", ], mean = summary_data[[g]]["mean", ],
               lwr = summary_data[[g]]["lwr", ], upr = summary_data[[g]]['upr', ])
  }) %>% {do.call(bind_rows, .)}
}

plot_surv_prob <- function(plot_data) {
  ggplot(plot_data) +
    geom_line(aes(x = time, y = true, linetype = "true")) +
    geom_line(aes(x = time, y = mean, linetype = "estimate")) +
    geom_ribbon(aes(x = time, ymin = lwr, ymax = upr), alpha = 0.2) +
    facet_grid(S ~ Z) + 
    xlab("Time") + ylab("Survival probability") + theme_bw()
}

add_to <- function(data, result, ER_type) {
  strata_treatments <- c(
    "never-taker (control arm)", 
    "never-taker (treatment arm)", 
    "complier (control arm)", 
    "complier (treatment arm)", 
    "always-taker (control arm)", 
    "always-taker (treatment arm)" 
  )
  for (i in 1:6) {
    new_df <- data.frame(
      true = result$combined_data[[i]]["true", ],
      mean = result$combined_data[[i]]["mean", ],
      lwr = result$combined_data[[i]]["lwr", ],
      upr = result$combined_data[[i]]["upr", ],
      time = result$time_sequence,
      ER_type = ER_type,
      strata_treatment = strata_treatments[i]
    )
    data <- bind_rows(data, new_df)
  }
  return (data)
}

simulation <- function(
    sim_data, truth_ER, fitted_ER, misspecification = FALSE,
    fitted_type = 'Cox', estimand = 'prob'
) {
  set.seed(0)
  if (truth_ER) {
    true_strata_list = c(1, 2, 2, 3)
    true_model_id = c(1, 1, 2, 3, 4, 4)
  }
  else {
    true_strata_list = c(1, 1, 2, 2, 3, 3)
    true_model_id = c(1, 2, 3, 4, 5, 6)
  }
  if (fitted_ER) {
    fitted_strata_list = c(1, 2, 2, 3)
    fitted_model_id = c(1, 1, 2, 3, 4, 4)
    strata = c(n = "00*", c = "01", a = "11*")
  }
  else {
    fitted_strata_list = c(1, 1, 2, 2, 3, 3)
    fitted_model_id = c(1, 2, 3, 4, 5, 6)
    strata = c(n = "00", c = "01", a = "11")
  }
  
  if (misspecification) {
    S.formula = Z + D ~ X1
    Y.formula = Y + delta ~ X1
  }
  else {
    S.formula = Z + D ~ X1 + X2
    Y.formula = Y + delta ~ X1 + X2
  }
  
  time_sequence <- seq(0.01, 15, length.out = 100)
  
  result <- PStrata::PStrata(
    S.formula = S.formula,
    Y.formula = Y.formula,
    Y.family = survival(fitted_type),
    data = sim_data$data,
    strata = strata,
    warmup = 500, iter = 1000,
    chains = 6, core = 6,
    prior_intercept = prior_normal(0, 10),
    survival.time.points = seq(0.01, 15, length.out = 100)
  )
  
  data_matrix <- cbind(1, sim_data$data$X1, sim_data$data$X2, 0)
  
  if (estimand == 'prob') {
    true_surv_prob <- stratum_surv_probability(
      sim_data$strata_coef_matrix, sim_data$outcome_coef_matrix, sim_data$phi,
      data_matrix, data_matrix, time_sequence, true_strata_list
    )
    if (fitted_type == "Cox") {
      fitted_surv_prob <- get_surv_prob_matrices(result, time_sequence, fitted_strata_list)
    }
    else {
      fitted_surv_prob <- get_surv_prob_matrices_AFT(result, time_sequence, fitted_strata_list)
    }
    summary_surv_prob <- summarize_surv_prob(fitted_surv_prob)
    combined_data <- combine_fitted_true_surv_prob(
      true_surv_prob, summary_surv_prob, true_model_id, fitted_model_id
    )
  }
  else {
    true_RACE <- stratum_RACE(
      sim_data$strata_coef_matrix, sim_data$outcome_coef_matrix, sim_data$phi,
      data_matrix, data_matrix, time_sequence, fitted_strata_list
    )
    if (fitted_type == "Cox") {
      fitted_RACE <- get_RACE_matrices(result, time_sequence, fitted_strata_list)
    }
    else {
      fitted_RACE <- get_RACE_matrices_AFT(result, time_sequence, fitted_strata_list)
    }
    summary_RACE <- summarize_surv_prob(fitted_RACE)
    combined_data <- combine_fitted_true_surv_prob(
      true_RACE, summary_RACE, true_model_id, fitted_model_id
    )
  }
  
  gPlot <- plot_surv_prob(make_plot_data(combined_data, time_sequence))
  return (list(
    time_sequence = time_sequence, data = sim_data, PStrata_result = result, combined_data = combined_data, plot = gPlot
  ))
}