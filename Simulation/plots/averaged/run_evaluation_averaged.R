# plots/averaged/run_evaluation_averaged.R 

# --- 1. Initialization ---
cat("===== Starting aggregation of simulation results (Final Duplication Check) =====\n")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(dplyr, tidyr, ggplot2, here, progress, forcats)

source(here("config.R"))
param_file_path <- here("true_params_for_simulation.csv")
if (file.exists(param_file_path)) {
  TRUE_PARAMS_DF <- read.csv(param_file_path)
} else {
  TRUE_PARAMS_DF <- data.frame(Parameter = character(), Value = numeric())
  warning("True parameter file not found.")
}
dir.create(here("plots", "averaged"), showWarnings = FALSE, recursive = TRUE)
dir.create(here(OUTPUT_RESULTS_PATH, "averaged"), showWarnings = FALSE, recursive = TRUE)


# --- 2. Load All Simulation Results ---
cat(" -> Loading all simulation results (RDS)...\n")
all_results <- list()
all_source_data <- list()
pb <- progress_bar$new(format = "  simulation :current/:total [:bar] :percent ETA: :eta", total = N_SIM)

for (i in 1:N_SIM) {
  result_file <- here(OUTPUT_RESULTS_PATH, paste0("filter_sim_seed_", i, "_m_1.RDS"))
  data_file <- here(OUTPUT_SIM_DATA_PATH, paste0("sim_data_seed_", i, "_m_1.RDS"))
  
  if(file.exists(result_file) && file.exists(data_file)) {
    all_results[[i]] <- readRDS(result_file)
    all_source_data[[i]] <- readRDS(data_file)
  }
  pb$tick()
}
cat(" -> Loading of all simulation results complete.\n")


# --- 2.5 Aggregate Evaluation Metrics (Sensitivity/Specificity) and Time-Point Percentages ---
cat(" -> Aggregating evaluation metrics and time-point percentages...\n")
all_metrics <- list()
all_percent_train_pred <- list()
all_percent_final_pred <- list()
all_percent_train_actual <- list()
all_percent_final_actual <- list()

# Helper function (from run_evaluation.R)
calculate_metrics <- function(true_states, pred_states, positive_class = 1) {
  TP <- sum(true_states == positive_class & pred_states == positive_class, na.rm = TRUE)
  FN <- sum(true_states == positive_class & pred_states != positive_class, na.rm = TRUE)
  TN <- sum(true_states != positive_class & pred_states != positive_class, na.rm = TRUE)
  FP <- sum(true_states != positive_class & pred_states == positive_class, na.rm = TRUE)
  
  sensitivity <- if ((TP + FN) > 0) TP / (TP + FN) else 0
  specificity <- if ((TN + FP) > 0) TN / (TN + FP) else 0
  
  return(list(sensitivity = sensitivity, specificity = specificity))
}

pb_metrics <- progress_bar$new(format = "  metrics calc :current/:total [:bar] :percent", total = N_SIM)

for (i in 1:N_SIM) {
  res <- all_results[[i]]
  dat <- all_source_data[[i]]
  pb_metrics$tick()
  
  if (is.null(res) || is.null(dat) || is.null(res$mPr_pred_best) || is.null(dat$S_true)) next
  
  prob_S2_raw <- res$mPr_pred_best[, 2:(dat$Nt + 1), 2] # (N, Nt)
  prob_S2 <- array(prob_S2_raw, dim = c(dat$N, dat$Nt))
  true_states_S <- dat$S_true[, 2:(dat$Nt + 1)]
  
  N_total <- nrow(prob_S2)
  Nt_final <- dat$Nt
  
  # Time-point percentage calculations (Predicted and Actual)
  prob_at_T_train <- prob_S2[, NT_TRAIN]
  true_at_T_train <- true_states_S[, NT_TRAIN]
  N_pred_train <- sum(prob_at_T_train > 0.5, na.rm = TRUE)
  N_actual_train <- sum(true_at_T_train == 2, na.rm = TRUE)
  all_percent_train_pred[[i]] <- (N_pred_train / N_total) * 100
  all_percent_train_actual[[i]] <- (N_actual_train / N_total) * 100
  
  prob_at_T_final <- prob_S2[, Nt_final]
  true_at_T_final <- true_states_S[, Nt_final]
  N_pred_final <- sum(prob_at_T_final > 0.5, na.rm = TRUE)
  N_actual_final <- sum(true_at_T_final == 2, na.rm = TRUE)
  all_percent_final_pred[[i]] <- (N_pred_final / N_total) * 100
  all_percent_final_actual[[i]] <- (N_actual_final / N_total) * 100
  
  # Evaluation metrics calculations (Sensitivity/Specificity)
  positive_class_S <- 2
  pred_states_S <- ifelse(prob_S2 > 0.5, 2, 1)
  
  observed_period <- 1:NT_TRAIN
  forecast_period <- (NT_TRAIN + 1):dat$Nt
  overall_period <- 1:dat$Nt
  
  if (max(forecast_period) > dat$Nt) {
    forecast_period <- forecast_period[forecast_period <= dat$Nt]
  }
  
  has_forecast_period <- length(forecast_period) > 0 && min(forecast_period) <= dat$Nt
  
  if(has_forecast_period) {
    metrics_obs <- calculate_metrics(true_states_S[, observed_period], pred_states_S[, observed_period], positive_class = positive_class_S)
    metrics_fore <- calculate_metrics(true_states_S[, forecast_period], pred_states_S[, forecast_period], positive_class = positive_class_S)
    metrics_all <- calculate_metrics(true_states_S[, overall_period], pred_states_S[, overall_period], positive_class = positive_class_S)
    
    metrics_for_seed <- data.frame(
      Seed = i,
      Sensitivity_Observed = metrics_obs$sensitivity,
      Specificity_Observed = metrics_obs$specificity,
      Sensitivity_Forecast = metrics_fore$sensitivity,
      Specificity_Forecast = metrics_fore$specificity,
      Sensitivity_Overall = metrics_all$sensitivity,
      Specificity_Overall = metrics_all$specificity
    )
  } else {
    metrics_obs <- calculate_metrics(true_states_S[, observed_period], pred_states_S[, observed_period], positive_class = positive_class_S)
    metrics_all <- metrics_obs
    
    metrics_for_seed <- data.frame(
      Seed = i,
      Sensitivity_Observed = metrics_obs$sensitivity,
      Specificity_Observed = metrics_obs$specificity,
      Sensitivity_Forecast = NA,
      Specificity_Forecast = NA,
      Sensitivity_Overall = metrics_all$sensitivity,
      Specificity_Overall = metrics_all$specificity
    )
  }
  all_metrics[[i]] <- metrics_for_seed
}
cat(" -> Aggregation of evaluation metrics and time-point percentages complete.\n")


# --- 2.6. Calculate and Output Average Time-Point Percentages (Output to Console) ---
cat(" -> Calculating average dropout intention percentages (Predicted vs Actual)...\n")
mean_percent_train_pred <- mean(unlist(all_percent_train_pred), na.rm = TRUE)
mean_percent_final_pred <- mean(unlist(all_percent_final_pred), na.rm = TRUE)
mean_percent_train_actual <- mean(unlist(all_percent_train_actual), na.rm = TRUE)
mean_percent_final_actual <- mean(unlist(all_percent_final_actual), na.rm = TRUE)

n_sim_included <- length(all_percent_train_pred[!sapply(all_percent_train_pred, is.null)])
Nt_final_ref <- all_source_data[[1]]$Nt

# Console Output
cat(paste0("    -> At T=", NT_TRAIN, " (end of observation):\n"))
cat(paste0("        Predicted: ", round(mean_percent_train_pred, 1), "% of individuals were predicted to be in Regime 2 (Average of ", n_sim_included, " simulations).\n"))
cat(paste0("        Actual:    ", round(mean_percent_train_actual, 1), "% of individuals were *actually* in Regime 2 (Average of ", n_sim_included, " simulations).\n"))
cat(paste0("    -> At T=", Nt_final_ref, " (final time point):\n"))
cat(paste0("        Predicted: ", round(mean_percent_final_pred, 1), "% of individuals were predicted to be in Regime 2 (Average of ", n_sim_included, " simulations).\n"))
cat(paste0("        Actual:    ", round(mean_percent_final_actual, 1), "% of individuals were *actually* in Regime 2 (Average of ", n_sim_included, " simulations).\n"))


# --- 3. Parameter Estimate Evaluation (Integrated CSV Output) ---
cat(" -> Starting aggregation and plotting of parameter estimates...\n")

add_param <- function(param_name, values, sim_idx, type_name) {
  df <- data.frame(Value = as.numeric(values), Simulation = sim_idx, Type = type_name)
  if (length(values) > 1) {
    df$Parameter <- paste0(param_name, "_", 1:length(values))
  } else {
    df$Parameter <- param_name
  }
  return(df)
}

if (length(all_results) > 0) {
  params_list <- list()
  for (i in 1:length(all_results)) {
    res <- all_results[[i]]
    if (is.null(res) || is.null(res$theta_best)) next
    theta <- res$theta_best
    
    # 3.1: Transform parameters back to original scale and collect into a list
    params_list[[length(params_list) + 1]] <- add_param("B11", theta$B11, i, "Intercepts (B1, State 1)")
    params_list[[length(params_list) + 1]] <- add_param("B12", theta$B11 + exp(theta$log_B12_delta), i, "Intercepts (B1, State 2)")
    params_list[[length(params_list) + 1]] <- add_param("B21d", theta$B21d, i, "Autoregressive (B2)")
    params_list[[length(params_list) + 1]] <- add_param("B22d", theta$B22d, i, "Autoregressive (B2)")
    params_list[[length(params_list) + 1]] <- add_param("B31", theta$B31, i, "Covariate Effects (B3)")
    params_list[[length(params_list) + 1]] <- add_param("B32", theta$B32, i, "Covariate Effects (B3)")
    
    lmd1f_params <- names(theta)[startsWith(names(theta), "Lmd1f")]
    if(length(lmd1f_params) > 0){
      lmd1f_names_sorted <- lmd1f_params[order(as.numeric(sub("Lmd1f", "", lmd1f_params)))]
      lmd1f_values <- sapply(lmd1f_names_sorted, function(p) theta[[p]])
      params_list[[length(params_list) + 1]] <- add_param("Lmd1f", lmd1f_values, i, "Factor Loadings (Lmd1)")
    }
    params_list[[length(params_list) + 1]] <- add_param("Lmd2f", theta$Lmd2f, i, "Factor Loadings (Lmd2)")
    
    params_list[[length(params_list) + 1]] <- add_param("Qd", exp(theta$log_Qd), i, "Variances (Q, R)")
    params_list[[length(params_list) + 1]] <- add_param("R1d", exp(theta$log_R1d), i, "Variances (Q, R)")
    params_list[[length(params_list) + 1]] <- add_param("R2d", exp(theta$log_R2d), i, "Variances (Q, R)")
    
    params_list[[length(params_list) + 1]] <- add_param("gamma1", theta$gamma1, i, "Markov Switching (Gamma)")
    params_list[[length(params_list) + 1]] <- add_param("gamma2", theta$gamma2, i, "Markov Switching (Gamma)")
    params_list[[length(params_list) + 1]] <- add_param("gamma3", theta$gamma3, i, "Markov Switching (Gamma)")
    params_list[[length(params_list) + 1]] <- add_param("gamma4", theta$gamma4, i, "Markov Switching (Gamma)")
    params_list[[length(params_list) + 1]] <- add_param("P(S1|S2)", plogis(theta$logit_tP_SB) * 0.1, i, "Transition Probabilities")
  }
  params_df <- bind_rows(params_list)
  
  if (nrow(params_df) > 0 && nrow(TRUE_PARAMS_DF) > 0) {
    
    # 3.2: Aggregate statistics (Mean, SD, Bias).
    
    # 1. Join estimated values (Value.x) to the True Values (Value.y)
    params_comparison_df <- params_df %>% inner_join(TRUE_PARAMS_DF, by = "Parameter")
    
    # 2. Group by Type and Parameter ONLY, then summarize
    final_params_summary <- params_comparison_df %>%
      group_by(Type, Parameter) %>%
      reframe(
        # True Value is constant within this group, so use first()
        True_Value = first(Value.y),
        # Calculate statistics across ALL simulation results (Value.x)
        Mean_Estimated = mean(Value.x, na.rm = TRUE),
        SD_Estimated = sd(Value.x, na.rm = TRUE),
        Bias = Mean_Estimated - True_Value
      ) %>%
      ungroup() %>%
      # Select final columns, ensuring 1 row per unique (Type, Parameter)
      select(Type, Parameter, True_Value, Mean_Estimated, SD_Estimated, Bias)
    
    # Save the integrated and non-redundant averaged_parameter_estimates.csv
    write.csv(
      final_params_summary,
      here(OUTPUT_RESULTS_PATH, "averaged", "averaged_parameter_estimates.csv"),
      row.names = FALSE
    )
    cat("    -> Saved: averaged_parameter_estimates.csv (Integrated with True_Value, SD, Bias) - 1 row per parameter.\n")
    
    # 3.3: Create distribution plot
    param_dist_plot <- ggplot(params_comparison_df, aes(x = Value.x)) +
      geom_density(aes(fill = "Estimated"), alpha = 0.5) +
      geom_vline(aes(xintercept = Value.y, color = "True"), linetype = "dashed", linewidth = 1) +
      facet_wrap(~ fct_inorder(Parameter), scales = "free", ncol = 5) +
      scale_fill_manual(name = "", values = c("Estimated" = "blue")) +
      scale_color_manual(name = "", values = c("True" = "red")) +
      labs(title = "Distribution of Parameter Estimates and True Values", x = "Parameter Value", y = "Density") +
      theme_bw() + theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggsave(here("plots", "averaged", "agg_param_distribution.png"), param_dist_plot, width = 16, height = 18)
    cat("    -> Saved: agg_param_distribution.png\n")
    
  } else {
    cat("    -> Skipping parameter aggregation/plot due to missing data.\n")
  }
}


# --- 4. Aggregate Forecast Error Data and Model Variance (Plotting only) ---
cat(" -> Aggregating forecast error and model variance...\n")
all_forecast_data_list <- list()

for (i in 1:N_SIM) {
  res <- all_results[[i]]
  dat <- all_source_data[[i]]
  if (is.null(res) || is.null(dat) || is.null(dat$eta1_true) || is.null(res$mEta_best)) next
  
  forecast_indices <- (NT_TRAIN + 1):NT_TOTAL
  max_time_data <- dim(dat$eta1_true)[2]
  max_time_res <- dim(res$mEta_best)[2]
  valid_forecast_indices <- forecast_indices[forecast_indices <= max_time_data & forecast_indices <= max_time_res]
  if (length(valid_forecast_indices) == 0) next
  
  eta_true_forecast <- dat$eta1_true[, valid_forecast_indices, , drop = FALSE]
  eta_est_raw <- res$mEta_best[, valid_forecast_indices, 1, , drop = FALSE]
  eta_est_forecast <- array(eta_est_raw, dim = dim(eta_true_forecast))
  
  p_est_var_forecast_diag <- array(NA, dim = dim(eta_true_forecast))
  if(length(dim(res$mP_best)) == 5) {
    for(subj in 1:dim(p_est_var_forecast_diag)[1]) {
      for(t_idx in 1:length(valid_forecast_indices)) {
        t_actual <- valid_forecast_indices[t_idx]
        p_matrix_raw <- res$mP_best[subj, t_actual, 1, , ]
        p_matrix <- array(p_matrix_raw, dim = dim(res$mP_best)[4:5])
        p_est_var_forecast_diag[subj, t_idx, ] <- diag(p_matrix)
      }
    }
  }
  
  squared_error <- (eta_true_forecast - eta_est_forecast)^2
  dimnames(squared_error) <- list(Subject = 1:dim(squared_error)[1], Time = 1:dim(squared_error)[2], Dimension = paste0("Dim ", 1:dim(squared_error)[3]))
  dimnames(p_est_var_forecast_diag) <- dimnames(squared_error)
  
  scores_df_long <- as.data.frame.table(squared_error, responseName = "Score")
  p_vars_df_long <- as.data.frame.table(p_est_var_forecast_diag, responseName = "p_var")
  
  combined_df <- merge(scores_df_long, p_vars_df_long) %>%
    mutate(Seed = i, Forecast_Time = as.numeric(Time)) %>%
    select(Seed, Forecast_Time, Dimension, Score, p_var)
  
  all_forecast_data_list[[i]] <- combined_df
}
forecast_df <- bind_rows(all_forecast_data_list)

# --- 5. Visualize Forecast Error (Plotting only) ---
cat(" -> Creating and saving forecast error plots...\n")
if (nrow(forecast_df) > 0) {
  plot_data <- forecast_df %>%
    group_by(Forecast_Time, Dimension) %>%
    summarise(Mean_Score = mean(Score, na.rm = TRUE), Mean_Var_P = mean(p_var, na.rm = TRUE), .groups = 'drop')
  
  score_plot_by_dim <- ggplot(plot_data, aes(x = Forecast_Time, group = Dimension)) +
    geom_ribbon(aes(ymin = pmax(0, Mean_Score - Mean_Var_P), ymax = Mean_Score + Mean_Var_P, fill = Dimension), alpha = 0.3) +
    geom_line(aes(y = Mean_Score, color = Dimension), linewidth = 1) +
    geom_point(aes(y = Mean_Score, color = Dimension), size = 2.5) +
    scale_color_manual(name = "Latent Dimension", values = c("Dim 1" = "blue4", "Dim 2" = "red4")) +
    scale_fill_manual(name = "Latent Dimension", values = c("Dim 1" = "skyblue", "Dim 2" = "salmon")) +
    labs(title = "Forecast Error (Score) Trend by Dimension", subtitle = "Line: Empirical Score (Squared Error), Ribbon: ±1 Model Predictive Variance (P)", x = "Forecast Horizon (Time Steps)", y = "Score (Squared Error)") +
    scale_x_continuous(breaks = 1:max(plot_data$Forecast_Time)) +
    theme_bw(base_size = 14) + theme(legend.position = "bottom")
  
  ggsave(here("plots", "averaged", "agg_score_plot_by_dimension.png"), score_plot_by_dim, width = 10, height = 7)
  cat("   -> 'agg_score_plot_by_dimension.png' saved.\n")
} else { cat("   -> No data for score or variance P. Skipped.\n") }


# --- 6. Aggregate and Visualize Evaluation Metrics (CSV Output) ---
cat(" -> Aggregating and creating evaluation metrics CSV...\n")

if (length(all_metrics) > 0) {
  metrics_df <- bind_rows(all_metrics[!sapply(all_metrics, is.null)])
  
  if (nrow(metrics_df) > 0) {
    # 6.1: Calculate the average of evaluation metrics across all seeds
    metrics_avg_summary <- metrics_df %>%
      summarise(
        Sensitivity_Observed = mean(Sensitivity_Observed, na.rm = TRUE),
        Specificity_Observed = mean(Specificity_Observed, na.rm = TRUE),
        Sensitivity_Forecast = mean(Sensitivity_Forecast, na.rm = TRUE),
        Specificity_Forecast = mean(Specificity_Forecast, na.rm = TRUE),
        Sensitivity_Overall = mean(Sensitivity_Overall, na.rm = TRUE),
        Specificity_Overall = mean(Specificity_Overall, na.rm = TRUE)
      )
    
    # 6.2: Format as averaged_evaluation_metrics_detailed.csv
    averaged_metrics_detailed_df <- data.frame(
      Period = c("Observed", "Forecast", "Overall"),
      Sensitivity = c(metrics_avg_summary$Sensitivity_Observed, metrics_avg_summary$Sensitivity_Forecast, metrics_avg_summary$Sensitivity_Overall),
      Specificity = c(metrics_avg_summary$Specificity_Observed, metrics_avg_summary$Specificity_Forecast, metrics_avg_summary$Specificity_Overall)
    )
    
    write.csv(
      averaged_metrics_detailed_df,
      here(OUTPUT_RESULTS_PATH, "averaged", "averaged_evaluation_metrics_detailed.csv"),
      row.names = FALSE
    )
    cat("    -> Saved: averaged_evaluation_metrics_detailed.csv\n")
    
    # 6.3: Create bar chart (logic unchanged)
    metrics_plot_data <- averaged_metrics_detailed_df %>%
      pivot_longer(
        cols = c(Sensitivity, Specificity),
        names_to = "Metric",
        values_to = "Mean_Value"
      ) %>%
      mutate(Period = factor(Period, levels = c("Observed", "Forecast", "Overall")))
    
    metrics_barchart <- ggplot(metrics_plot_data, aes(x = Metric, y = Mean_Value, fill = Period)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
      geom_text(aes(label = round(Mean_Value, 3)), vjust = -0.3, position = position_dodge(width = 0.9), size = 3.2) +
      labs(title = "Averaged Evaluation Metrics by Period", subtitle = paste0("Average of ", nrow(metrics_df), " simulations (out of ", N_SIM, ")"), x = "Metric", y = "Average Value", fill = "Period") +
      scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.2)) +
      theme_minimal(base_size = 14) + theme(legend.position = "bottom")
    
    ggsave(filename = here("plots", "averaged", "agg_metrics_barchart.jpg"), plot = metrics_barchart, width = 8, height = 7, dpi = 300)
    cat(" -> Evaluation metrics plot saved.\n")
  } else {
    warning("Evaluation metrics data was empty. Plot not created.")
  }
} else {
  warning("Evaluation metrics files not found. Plot not created.")
}


cat("===== All aggregation processes complete (Created averaged_parameter_estimates.csv and averaged_evaluation_metrics_detailed.csv) =====\n")