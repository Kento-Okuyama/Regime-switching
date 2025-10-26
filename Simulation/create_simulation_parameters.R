# create_simulation_parameters.R
# Purpose: Read parameter estimates from the empirical study results, average them, and 
# perform necessary inverse transformations to use them as ground truth for the simulation study.

# --- 1. Initialization ---
cat("===== Creating Averaged and Transformed True Values for Simulation =====\n")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(dplyr, here, tidyr)

# --- 2. Load all parameter estimate files from the empirical study ---
cat(" -> Loading parameter estimates from all imputation runs of the empirical study...\n")
empirical_results_path <- here("..", "Empirical", "output", "results", "averaged")
param_files <- list.files(
  path = empirical_results_path,
  pattern = "averaged_parameter_estimates.csv",
  full.names = TRUE
)
if (length(param_files) == 0) {
  stop(paste("No parameter estimate file found for the empirical study at the specified path:", empirical_results_path, ". Please run the empirical study first."))
}
cat(paste0("  -> ", length(param_files), " parameter files found.\n"))
all_params_df <- lapply(param_files, read.csv) %>%
  bind_rows(.id = "Imputation")

# --- 3. Calculate the mean value for each parameter ---
cat(" -> Calculating the mean value for each parameter...\n")
if ("Mean_Value" %in% names(all_params_df)) {
  averaged_params_df <- all_params_df %>%
    group_by(Parameter) %>%
    summarise(Value = mean(Mean_Value, na.rm = TRUE)) %>%
    ungroup()
} else {
  averaged_params_df <- all_params_df %>%
    group_by(Parameter) %>%
    summarise(Value = mean(Value, na.rm = TRUE)) %>%
    ungroup()
}


# --- 4. Inverse Transformation and Parameter Renaming ---
cat(" -> Performing inverse transformation and parameter renaming...\n")

# Extract B11 and B12 and calculate log_B12_delta
b11_b12_wide <- averaged_params_df %>%
  filter(grepl("^(B11_|B12_)", Parameter)) %>%
  mutate(
    BaseParam = sub("(_[0-9]+)$", "", Parameter),
    Index = sub("^.*_([0-9]+)$", "\\1", Parameter)
  ) %>%
  select(-Parameter) %>%
  pivot_wider(names_from = BaseParam, values_from = Value)

log_b12_delta_df <- b11_b12_wide %>%
  filter(!is.na(B11) & !is.na(B12)) %>%
  mutate(
    Parameter = paste0("log_B12_delta_", Index),
    Value = log(pmax(B12 - B11, 1e-9))
  ) %>%
  select(Parameter, Value)

# Create log-scale versions of Q and R as separate data frames
cat(" -> Creating log-scale versions of Q and R...\n")
log_qr_df <- averaged_params_df %>%
  filter(grepl("^(Qd|R1d|R2d)_", Parameter)) %>%
  mutate(
    Value = log(pmax(Value, 1e-9)), 
    Parameter = case_when(
      grepl("^Qd_", Parameter) ~ sub("Qd_", "log_Qd_", Parameter),
      grepl("^R1d_", Parameter) ~ sub("R1d_", "log_R1d_", Parameter),
      grepl("^R2d_", Parameter) ~ sub("R2d_", "log_R2d_", Parameter)
    )
  )


# Create the final data frame
final_params_df <- averaged_params_df %>%
  bind_rows(log_b12_delta_df) %>%     
  bind_rows(log_qr_df) %>%             
  mutate(
    Parameter = case_when(
      Parameter == "Lmd2_2" ~ "Lmd2f_1",
      Parameter == "Lmd2_3" ~ "Lmd2f_2",
      TRUE ~ Parameter
    )
  )

# --- 5. Save the final parameters to a new file ---
output_path <- here("true_params_for_simulation.csv")
write.csv(final_params_df, output_path, row.names = FALSE, quote = FALSE)
cat(paste0(" -> Correctly transformed parameters saved to '", basename(output_path), "'.\n"))
cat(paste0("    (Includes both pre- and post-transformed values for B12, Qd, R1d, and R2d)\n"))
cat("=======================================================================\n\n")