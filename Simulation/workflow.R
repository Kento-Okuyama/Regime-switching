# workflow.R

# Set the working directory
setwd("C:/Users/kento/OneDrive - UT Cloud (1)/PhD/Regime-switching/Simulation")

# --- NEW STEP ---
# Run this first to create the averaged parameter file for the simulation
source("create_simulation_parameters.R") 

cat(" -> Modifying true_params_for_simulation.csv: Setting gamma1 to 3.5\n")
param_file_path <- "true_params_for_simulation.csv"

if (file.exists(param_file_path)) {
  # Read the file (specifying stringsAsFactors = FALSE)
  params_df <- read.csv(param_file_path, stringsAsFactors = FALSE)
  
  # Find the value of gamma1 and replace it with 3.5
  gamma1_index <- which(params_df$Parameter == "gamma1")
  
  if (length(gamma1_index) > 0) {
    original_value <- params_df$Value[gamma1_index[1]]
    params_df$Value[gamma1_index[1]] <- 3.5
    
    # Overwrite the file and save
    write.csv(params_df, param_file_path, row.names = FALSE, quote = FALSE)
    
    cat(paste0("  -> Successfully modified gamma1 from ", original_value, " to 3.5.\n"))
  } else {
    warning("  -> 'gamma1' not found in true_params_for_simulation.csv. File was not modified.\n")
  }
} else {
  warning(paste("  ->", param_file_path, "not found. Cannot modify gamma1.\n"))
}

# Run the simulation using the newly created parameter file
source("main.R")

# Run the evaluation scripts for simulation results
source("plots/individual/run_evaluation.R")
source("plots/averaged/run_evaluation_averaged.R")