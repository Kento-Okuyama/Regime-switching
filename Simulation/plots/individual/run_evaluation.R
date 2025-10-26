# plots/run_evaluation.R

# --- 1. Initialization ---
source(here::here("config.R"))
source(here::here("R", "00_load_libraries.R"))
load_libraries()

# Define plot path not defined in config.R
OUTPUT_PLOTS_PATH <- here::here("plots") 

# Create directories using the variable name defined in config.R (OUTPUT_RESULTS_PATH)
dir.create(here(OUTPUT_RESULTS_PATH, "individual"), recursive = TRUE, showWarnings = FALSE)
dir.create(here(OUTPUT_PLOTS_PATH, "individual"), recursive = TRUE, showWarnings = FALSE)


# Helper function to calculate evaluation metrics
calculate_metrics <- function(true_states, pred_states, positive_class = 1) {
  TP <- sum(true_states == positive_class & pred_states == positive_class, na.rm = TRUE)
  FN <- sum(true_states == positive_class & pred_states != positive_class, na.rm = TRUE)
  TN <- sum(true_states != positive_class & pred_states != positive_class, na.rm = TRUE)
  FP <- sum(true_states != positive_class & pred_states == positive_class, na.rm = TRUE)
  
  sensitivity <- if ((TP + FN) > 0) TP / (TP + FN) else 0
  specificity <- if ((TN + FP) > 0) TN / (TN + FP) else 0
  
  return(list(sensitivity = sensitivity, specificity = specificity))
}


# --- 2. Data Collection for Participant Selection ---
cat("===== Step 1: Collecting data from all seeds for participant selection =====\n")

result_files <- list.files(OUTPUT_RESULTS_PATH, pattern = "filter_sim_seed_.*_m_.*\\.RDS$", full.names = TRUE) 
if (length(result_files) == 0) stop(paste("Result files not found in:", OUTPUT_RESULTS_PATH, ". Please run main.R."))

prob_arrays_for_selection <- list()
best_model_files <- character(N_SIM) 

# For simulation data, m_current (Imputation) is fixed at 1
m_current <- 1 

# Loop through each Seed
for (seed_current in 1:N_SIM) {
  files_for_seed_m <- grep(paste0("filter_sim_seed_", seed_current, "_m_", m_current, "\\.RDS$"), result_files, value = TRUE)
  
  if (length(files_for_seed_m) == 0) next
  
  best_file_for_seed_m <- files_for_seed_m[1] 
  
  if (!is.null(best_file_for_seed_m)) {
    filtered_best <- readRDS(best_file_for_seed_m)
    
    imputed_data_path <- file.path(OUTPUT_SIM_DATA_PATH, paste0("sim_data_seed_", seed_current, "_m_", m_current, ".RDS"))
    if(file.exists(imputed_data_path)) {
      df_temp <- readRDS(imputed_data_path)
      if(!is.null(filtered_best$mPr_pred_best)){
        prob_arrays_for_selection[[seed_current]] <- filtered_best$mPr_pred_best[, 2:(df_temp$Nt + 1), 2] 
        best_model_files[seed_current] <- best_file_for_seed_m
      }
    }

  }
}

valid_indices <- which(!sapply(prob_arrays_for_selection, is.null))
if (length(valid_indices) == 0) stop("Could not collect valid prediction probability data from any seed.")

# Average probability across all Seeds (N x Nt)
prob_avg <- Reduce('+', prob_arrays_for_selection[valid_indices]) / length(valid_indices)

# (Load Seed 1 data as reference)
df_ref_path <- file.path(OUTPUT_SIM_DATA_PATH, "sim_data_seed_1_m_1.RDS")
if (!file.exists(df_ref_path)) stop("Reference data file (sim_data_seed_1_m_1.RDS) not found.")
df_ref <- readRDS(df_ref_path)


# --- 3. Select Common Participants ---
cat("\n===== Step 2: Selecting 3 common participants =====\n")
# Set as dummy values
person_indices <- c(1, 2, 3) 
person_indices <- person_indices[person_indices <= df_ref$N] # Handle small N
plot_subtitle <- paste("Selected Participants:", paste(person_indices, collapse = ", "))


# --- 4. Loop through each seed to evaluate and visualize ---
cat("\n===== Step 3: Creating plots for each seed (m=1) =====\n")
evaluation_metrics <- data.frame()

for (seed_current in 1:N_SIM) {
  cat(paste0("\n===== Starting evaluation for Seed #", seed_current, " (m=", m_current, ") =====\n"))
  
  best_file_for_m <- best_model_files[seed_current]
  if (is.na(best_file_for_m) || !file.exists(best_file_for_m)) {
    cat("  Best model not found. Skipping.\n")
    next
  }
  
  cat("  Best model:", basename(best_file_for_m), "\n")
  
  cat("  Loading corresponding simulation data...\n")
  imputed_data_path <- file.path(OUTPUT_SIM_DATA_PATH, paste0("sim_data_seed_", seed_current, "_m_", m_current, ".RDS"))
  if (!file.exists(imputed_data_path)) {
    cat("  Data file not found. Skipping.\n"); next
  }
  df <- readRDS(imputed_data_path)
  
  data_true_S_raw <- df$S_true # (N, Nt+1)
  data_true_S <- data_true_S_raw[, 2:(df$Nt + 1)] # (N, Nt)
  data_true_DO <- df$DO # (N, Nt)
  
  filtered_best <- readRDS(best_file_for_m)
  
  if (is.null(filtered_best$mPr_pred_best)) {
    cat("  Model results do not contain mPr_pred_best data for plotting. Skipping.\n")
    next
  }
  
  cat("  [1/4] Saving parameter estimates to file...\n")
  params_df <- data.frame(Parameter = character(), Value = numeric())
  if (!is.null(filtered_best$theta_best)) {
    for (param_name in names(filtered_best$theta_best)) {
      values <- as.vector(filtered_best$theta_best[[param_name]])
      for (i in 1:length(values)) {
        indexed_name <- if (length(values) > 1) paste0(param_name, "_", i) else param_name
        params_df <- rbind(params_df, data.frame(Parameter = indexed_name, Value = values[i]))
      }
    }
    param_filename <- paste0("parameter_estimates_seed_", seed_current, "_m_", m_current, ".csv")
    write.csv(params_df, here(OUTPUT_RESULTS_PATH, "individual", param_filename), row.names = FALSE)
    cat(paste0("    -> Saved: ", param_filename, "\n"))
  } else {
    cat("    -> No theta_best found. Skipping parameter save.\n")
  }
  
  
  cat("  [2/4] Calculating sensitivity and specificity for all periods...\n")
  
  true_states_S <- data_true_S 
  positive_class_S <- 2 # S=2
  
  prob_S2 <- filtered_best$mPr_pred_best[, 2:(df$Nt + 1), 2] # (N, Nt)
  pred_states_S <- ifelse(prob_S2 > 0.5, 2, 1) 
  
  observed_period <- 1:NT_TRAIN
  forecast_period <- (NT_TRAIN + 1):df$Nt
  overall_period <- 1:df$Nt
  
  if (max(forecast_period) > df$Nt) {
    forecast_period <- forecast_period[forecast_period <= df$Nt]
  }
  
  if(length(forecast_period) > 0) {
    metrics_obs <- calculate_metrics(true_states_S[, observed_period], pred_states_S[, observed_period], positive_class = positive_class_S)
    metrics_fore <- calculate_metrics(true_states_S[, forecast_period], pred_states_S[, forecast_period], positive_class = positive_class_S)
    metrics_all <- calculate_metrics(true_states_S[, overall_period], pred_states_S[, overall_period], positive_class = positive_class_S)
    
    metrics_for_m <- data.frame(
      Seed = seed_current,
      Imputation = m_current,
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
    
    metrics_for_m <- data.frame(
      Seed = seed_current,
      Imputation = m_current,
      Sensitivity_Observed = metrics_obs$sensitivity,
      Specificity_Observed = metrics_obs$specificity,
      Sensitivity_Forecast = NA,
      Specificity_Forecast = NA,
      Sensitivity_Overall = metrics_all$sensitivity,
      Specificity_Overall = metrics_all$specificity
    )
  }
  
  evaluation_metrics <- rbind(evaluation_metrics, metrics_for_m)
  cat("    -> Metrics calculated.\n")
  
  cat("  [3/4D] Creating regime probability plots...\n")
  
  prob_S2_all_t <- prob_S2 
  true_S_binary <- ifelse(data_true_S == 2, 1, 0)
  
  df_pred_long <- as.data.frame(prob_S2_all_t) %>% 
    dplyr::mutate(ID = dplyr::row_number()) %>%
    tidyr::pivot_longer(cols = -ID, names_to = "Time", values_to = "Predicted_Prob_S2", names_prefix = "V") %>% 
    dplyr::mutate(Time = as.numeric(Time))
  
  df_true_long <- as.data.frame(true_S_binary) %>% 
    dplyr::mutate(ID = dplyr::row_number()) %>%
    tidyr::pivot_longer(cols = -ID, names_to = "Time", values_to = "True_Value_S2", names_prefix = "V") %>% 
    dplyr::mutate(Time = as.numeric(Time))
  
  df_combined <- dplyr::inner_join(df_pred_long, df_true_long, by = c("ID", "Time"))
  
  plot3a <- ggplot(df_combined %>% dplyr::filter(ID %in% person_indices), aes(x = Time)) +
    geom_line(aes(y = Predicted_Prob_S2, color = "Predicted")) +
    geom_line(aes(y = True_Value_S2, color = "True"), linetype = "dashed") +
    facet_wrap(~ ID) +
    labs(title = paste("Predicted Regime Probability vs. True Regime (Seed #", seed_current, ", m=", m_current, ")"), 
         subtitle = plot_subtitle,
         x = "Time Point", y = "Regime 2 Prob / State (1=S2)") +
    scale_color_manual(name = "Legend", values = c("Predicted" = "black", "True" = "red"), labels = c("Prediction (Regime 2)", "Actual (Regime 2)")) +
    geom_hline(yintercept = 0.5, linetype="dotted", color = "grey40") +
    geom_vline(xintercept = NT_TRAIN + 0.5, linetype="dashed", color = "blue") + 
    theme_bw() + theme(legend.position = "bottom")
  
  ggsave(here(OUTPUT_PLOTS_PATH, paste0("individual/regime_prob_selected_seed_", seed_current, "_m_", m_current, ".png")), plot3a, width = 12, height = 5)
  
  # Add plot for all participants
  
  plot3b <- ggplot(df_combined, aes(x = Time)) +
    geom_line(aes(y = Predicted_Prob_S2, color = "Predicted"), alpha = 0.8) +
    geom_line(aes(y = True_Value_S2, color = "True"), linetype = "dashed") +
    facet_wrap(~ ID) +
    labs(title = paste("Predicted Regime Probability vs. True Regime (All Participants, Seed #", seed_current, ", m=", m_current, ")"),
         x = "Time Point", y = "Regime 2 Prob / State (1=S2)") +
    scale_color_manual(name = "Legend", values = c("Predicted" = "blue", "True" = "red"), labels = c("Prediction (Regime 2)", "Actual (Regime 2)")) +
    geom_hline(yintercept = 0.5, linetype="dotted", color = "grey40") +
    geom_vline(xintercept = NT_TRAIN + 0.5, linetype="dashed", color = "blue") +
    theme_bw() + theme(legend.position = "bottom")
  
  ggsave(here(OUTPUT_PLOTS_PATH, paste0("individual/regime_prob_all_seed_", seed_current, "_m_", m_current, ".png")), plot3b, width = 16, height = 12)
  
  cat("    -> Saved: 2 types of regime probability plots\n")
  
  cat("  [4/4] Creating latent variable trajectory plots...\n")
  
  series_labels <- c("Factor 1 (Afraid to Fail)", "Factor 2 (No Positive Affect)")
  
  plot_data_eta <- data.frame()
  
  if (is.null(filtered_best$mEta_best) || length(dim(filtered_best$mEta_best)) != 4) {
    cat("    -> mEta_best data is missing or has wrong dimensions. Skipping latent plot.\n")
    next
  }
  
  eta_data <- filtered_best$mEta_best[, 2:(df$Nt+1), 1, ] # (N, Nt, U1)
  p_data <- filtered_best$mP_best[, 2:(df$Nt+1), 1, , ] # (N, Nt, U1, U1)
  
  for (person in person_indices) {
    if (person > dim(eta_data)[1]) next
    for (l in 1:U1) {
      plot_data_eta <- rbind(plot_data_eta, data.frame(
        Time = 1:df$Nt, 
        Value = eta_data[person, , l], 
        StdDev = sqrt(p_data[person, , l, l]), 
        Series = factor(series_labels[l], levels = series_labels), 
        Person = factor(person)
      ))
    }
  }
  
  if (nrow(plot_data_eta) > 0) {
    plot4 <- ggplot(plot_data_eta, aes(x = Time, y = Value, color = as.factor(Person), fill = as.factor(Person))) +
      geom_line(size = 1) +
      geom_ribbon(aes(ymin = Value - StdDev, ymax = Value + StdDev), alpha = 0.2, linetype = 0) +
      facet_wrap(~ Series, scales = "free_y") +
      labs(title = paste("Latent Variable Trajectories (Regime 1, Seed #", seed_current, ", m=", m_current, ")"), 
           subtitle = plot_subtitle,
           x = "Time Point", y = "Value", color = "Participant ID", fill = "Participant ID") +
      geom_vline(xintercept = NT_TRAIN + 0.5, linetype="dashed", color = "blue") + 
      theme_bw() + theme(legend.position = "bottom")
    
    ggsave(here(OUTPUT_PLOTS_PATH, paste0("individual/latent_trajectory_selected_seed_", seed_current, "_m_", m_current, ".png")), plot4, width = 12, height = 6)
    cat("    -> Saved: Latent variable trajectory plot\n")
  } else {
    cat("    -> No data generated for latent plot. Skipping.\n")
  }
}

# --- 5. Final Summary of Evaluation Metrics ---
cat("\n===== Summary of Evaluation Metrics Across All Seeds (m=1) =====\n")

summary_metrics <- evaluation_metrics %>%
  dplyr::summarise(
    Mean_Sens_Obs = mean(Sensitivity_Observed, na.rm = TRUE),
    Mean_Spec_Obs = mean(Specificity_Observed, na.rm = TRUE),
    Mean_Sens_Fore = mean(Sensitivity_Forecast, na.rm = TRUE),
    Mean_Spec_Fore = mean(Specificity_Forecast, na.rm = TRUE),
    Mean_Sens_All = mean(Sensitivity_Overall, na.rm = TRUE),
    Mean_Spec_All = mean(Specificity_Overall, na.rm = TRUE)
  )

cat("--- Average Metrics Across All Seeds ---\n")
print(summary_metrics)

cat("\n--- Detailed Metrics per Seed ---\n")
print(evaluation_metrics)

write.csv(evaluation_metrics, here(OUTPUT_RESULTS_PATH, "individual", "evaluation_metrics_summary_per_seed.csv"), row.names = FALSE)
write.csv(summary_metrics, here(OUTPUT_RESULTS_PATH, "individual", "evaluation_metrics_summary_average.csv"), row.names = FALSE)

cat(paste0("\nSummary saved to '", here(OUTPUT_RESULTS_PATH, "individual"), "'.\n"))

cat("\nAll evaluation and visualization tasks are complete.\n")