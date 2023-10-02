source("Simulation util.R")

gen_sim_data <- function(
    N, strata_coef_matrix, outcome_coef_matrix, phi, 
    true_model_id, beta_C
) {
  set.seed(0)
  X1 <- rbinom(N, 1, 0.5)
  X2 <- rnorm(N)
  data_matrix <-  cbind(1, X1, X2, 0)
  S <- sample_S(strata_coef_matrix, data_matrix, c(0, 1, 3))
  Z <- rbinom(N, 1, 0.5)
  D <- ifelse(bitwAnd(S, 2^(1-Z)), 1, 0)
  T_all <- sample_Y_cox(outcome_coef_matrix, phi, data_matrix)
  T <- T_all[, true_model_id[1]] * (S == 0 & Z == 0) + 
    T_all[, true_model_id[2]] * (S == 0 & Z == 1) + 
    T_all[, true_model_id[3]] * (S == 1 & Z == 0) +
    T_all[, true_model_id[4]] * (S == 1 & Z == 1) +
    T_all[, true_model_id[5]] * (S == 3 & Z == 0) + 
    T_all[, true_model_id[6]] * (S == 3 & Z == 1) 
  C <- rexp(N, exp(data_matrix %*% beta_C))
  Y <- pmin(T, C)
  delta <- ifelse(T < C, 1, 0)
  sim_data <- data.frame(
    S, Z, D, T, C, Y, X1, X2, delta
  )
  return (list(
    data = sim_data,
    strata_coef_matrix = strata_coef_matrix,
    outcome_coef_matrix = outcome_coef_matrix,
    phi = phi
  ))
}

data_1 <- gen_sim_data(
  N = 2000, 
  strata_coef_matrix = t(matrix(
    c(
      -1, 0.4, 0.2, 0, 
      0, 0.5, 0, 1, 
      -1, -0.2, -0.6, 0
    ), ncol = 4, byrow = TRUE
  )),
  outcome_coef_matrix = t(matrix(
    c(
      -3, 0.4, 0.1, 0, 
      -2.4, -0.1, 0.3, 0.15, 
      -1.8, 0.2, -0.2, -0.05, 
      -1.2, 0.1, 0, 0
    ), ncol = 4, byrow = TRUE
  )),
  phi = c(2, 1.5, 1.5, 1.0),
  true_model_id = c(1, 1, 2, 3, 4, 4),
  beta_C = c(-0.3, 0.1, -0.2, 0)
)

data_2 <- gen_sim_data(
  N = 2000, 
  strata_coef_matrix = t(matrix(
    c(
      -1, 0.4, 0.2, 0, 
      0, 0.5, 0, 1, 
      -1, -0.2, -0.6, 0
    ), ncol = 4, byrow = TRUE
  )),
  outcome_coef_matrix = t(matrix(
    c(
      -3, 0.4, 0.1, 0, 
      -2.4, 0.2, -0.3, 0,
      -2.4, -0.1, 0.3, 0.15, 
      -1.8, 0.2, -0.2, -0.05, 
      -1.2, 0.1, 0, 0,
      -0.6, 0, 0.2, 0
    ), ncol = 4, byrow = TRUE
  )),
  phi = c(2, 2, 1.5, 1.5, 1.0, 1.0),
  true_model_id = c(1, 2, 3, 4, 5, 6),
  beta_C = c(-0.3, 0.1, -0.2, 0)
)

data_3 <- list()
for (i in 1:5) {
  data_3[[i]] <- gen_sim_data(
    N = 2000, 
    strata_coef_matrix = t(matrix(
      c(
        -3 + i, 0.4, 0.2, 0, 
        0, 0.5, 0, 1, 
        -3 + i, -0.2, -0.6, 0
      ), ncol = 4, byrow = TRUE
    )),
    outcome_coef_matrix = t(matrix(
      c(
        -3, 0.4, 0.1, 0, 
        -2.4, -0.1, 0.3, 0.15, 
        -1.8, 0.2, -0.2, -0.05, 
        -1.2, 0.1, 0, 0
      ), ncol = 4, byrow = TRUE
    )),
    phi = c(2, 1.5, 1.5, 1.0),
    true_model_id = c(1, 1, 2, 3, 4, 4),
    beta_C = c(-0.3, 0.1, -0.2, 0)
  )
}
