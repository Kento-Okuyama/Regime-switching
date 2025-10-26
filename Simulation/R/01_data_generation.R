# R/01_data_generation.R

library(mvtnorm)
library(here)

#' Function to load true parameters from a CSV and return them as a list
load_true_parameters <- function(file_path) {
  params_df <- read.csv(file_path)
  params <- setNames(as.list(params_df$Value), params_df$Parameter)
  return(params)
}

#' Generate simulation data based on true parameters
generate_data <- function(seed, N, Nt, O1, O2, U1, true_params, N_IMPUTE, mice_maxit) {
  
  set.seed(seed)
  
  # --- 1. Correctly Extract Parameters Present in the CSV ---
  
  # ▼▼▼【Correction 1】Specify the 1st and 2nd factors to match the CSV ▼▼▼
  factor_map_old_to_new <- c(1, 2) 
  
  true_params[["gamma1"]] <- 3.5 # This value is overwritten from the workflow.R step
  B11           <- sapply(factor_map_old_to_new, function(i) true_params[[paste0("B11_", i)]])
  log_B12_delta <- sapply(factor_map_old_to_new, function(i) true_params[[paste0("log_B12_delta_", i)]])
  B21d          <- sapply(factor_map_old_to_new, function(i) true_params[[paste0("B21d_", i)]])
  B22d          <- sapply(factor_map_old_to_new, function(i) true_params[[paste0("B22d_", i)]])
  B31           <- sapply(factor_map_old_to_new, function(i) true_params[[paste0("B31_", i)]])
  B32           <- sapply(factor_map_old_to_new, function(i) true_params[[paste0("B32_", i)]])
  log_Qd        <- sapply(factor_map_old_to_new, function(i) true_params[[paste0("log_Qd_", i)]])
  gamma2        <- sapply(factor_map_old_to_new, function(i) true_params[[paste0("gamma2_", i)]])
  gamma4        <- sapply(factor_map_old_to_new, function(i) true_params[[paste0("gamma4_", i)]])
  gamma1        <- true_params[["gamma1"]]
  gamma3        <- true_params[["gamma3"]]
  
  # Specify Lmd1f numbers 1 through 4 to match the CSV 
  lmd1f_names <- c("Lmd1f1", "Lmd1f2", "Lmd1f3", "Lmd1f4")
  Lmd1f <- sapply(lmd1f_names, function(name) true_params[[name]])
  
  # Specify log_R1d numbers 1 through 6 to match the CSV
  r1d_indices <- 1:6
  log_R1d <- sapply(r1d_indices, function(i) true_params[[paste0("log_R1d_", i)]])
  
  Lmd2f <- sapply(1:2, function(i) true_params[[paste0("Lmd2f_", i)]])
  R2d <- sapply(1:O2, function(i) true_params[[paste0("R2d_", i)]])
  
  # --- 2. Construct Matrices from Parameters ---
  B12 <- B11 + exp(log_B12_delta)
  B1 <- rbind(B11, B12)
  B21 <- diag(B21d)
  B22 <- diag(B22d)
  B2 <- array(c(B21, B22), dim = c(U1, U1, 2))
  B3 <- rbind(B31, B32)
  Q <- diag(exp(log_Qd))
  R1 <- diag(exp(log_R1d))
  Lmd2 <- c(1, Lmd2f)
  R2 <- diag(R2d)
  
  # Completely rewrite Lmd1 matrix for a 2-factor structure
  Lmd1 <- matrix(0, nrow = O1, ncol = U1)
  Lmd1[1,1] <- 1
  Lmd1[2,1] <- Lmd1f[1] 
  Lmd1[3,1] <- Lmd1f[2] 
  Lmd1[4,2] <- 1
  Lmd1[5,2] <- Lmd1f[3] 
  Lmd1[6,2] <- Lmd1f[4] 
  
  # --- 3. Data Generation (Complete Data) ---
  y1_complete <- array(NA, dim = c(N, Nt, O1))
  y2_complete <- array(NA, dim = c(N, O2))
  eta1 <- array(NA, dim = c(N, Nt + 1, U1))
  eta2 <- array(NA, dim = c(N, 1))
  S <- matrix(NA, nrow = N, ncol = Nt + 1)
  DO <- matrix(0, nrow = N, ncol = Nt)
  
  eta1[, 1, ] <- rmvnorm(N, mean = rep(0, U1), sigma = diag(U1)/2)
  S[, 1] <- 1  
  eta2[, 1] <- rnorm(N, 0, 1)
  y2_complete <- t(Lmd2 %*% t(eta2) + t(rmvnorm(N, mean = rep(0, O2), sigma = R2)))
  
  for (t in 1:Nt) {
    B1_t <- B1[S[, t], ]; B3_t <- B3[S[, t], ]; eta1_t_minus_1 <- eta1[, t, ]
    autoreg_term <- t(sapply(1:N, function(i) eta1_t_minus_1[i, ] %*% B2[, , S[i, t]]))
    eta1_mean <- B1_t + autoreg_term + (eta2[, 1] * B3_t)
    eta1[, t + 1, ] <- eta1_mean + rmvnorm(N, mean = rep(0, U1), sigma = Q)
    y1_complete[, t, ] <- eta1[, t + 1, ] %*% t(Lmd1) + rmvnorm(N, mean = rep(0, O1), sigma = R1)
    
    eta1_pred_for_tpr <- eta1[, t + 1, ]
    interaction_term <- rowSums(eta1_pred_for_tpr * gamma4) * eta2[,1]
    prob_s1 <- plogis(gamma1 + eta1_pred_for_tpr %*% gamma2 + eta2[,1] * gamma3 + interaction_term)
    probs <- ifelse(S[, t] == 1, prob_s1, 0)
    S[, t + 1] <- rbinom(N, 1, probs) * (-1) + 2
    
    if (t > 1) {
      dropped_out_indices <- which(DO[, t - 1] == 1)
      if (length(dropped_out_indices) > 0) S[dropped_out_indices, t + 1] <- 2
    }
    
    prob_if_S2 <- 0
    dropout_prob <- ifelse(S[, t + 1] == 2, prob_if_S2, 0) 
    DO[, t] <- rbinom(N, 1, dropout_prob)
    if (t > 1) DO[, t][DO[, t - 1] == 1] <- 1
  }
  
  # --- 4. Structure and Save the Complete Dataset to a List ---
  cat("  Structuring and saving the complete dataset...\n")
  y1_t1_mean <- apply(y1_complete[, 1, ], 2, mean, na.rm=TRUE)
  y1_centered <- array(NA, dim=dim(y1_complete))
  for (var in 1:O1) y1_centered[,,var] <- y1_complete[,,var] - y1_t1_mean[var]
  y2_std <- apply(y2_complete, 2, scale)
  
  df <- list(y1 = y1_centered, y2 = y2_std, O1 = O1, O2 = O2, U1 = U1, 
             N = N, Nt = Nt, DO = DO, S_true = S, eta1_true = eta1, eta2_true = eta2)
  
  complete_datasets <- list()
  complete_datasets[[1]] <- df
  
  cat("All data generation processes are complete.\n")
  return(complete_datasets)
}