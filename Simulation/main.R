# main.R

# --- 1. Initialization ---
source("config.R")
source(here("R", "00_load_libraries.R"))
source(here("R", "01_data_generation.R"))
source(here("R", "02_filtering.R")) 
load_libraries()

# --- 2. Load true parameters for simulation ---
cat("Loading true parameters for simulation...\n")
param_file_path <- "true_params_for_simulation.csv" # Load the averaged parameters
if (!file.exists(param_file_path)) {
  stop(paste("Parameter file not found:", param_file_path, "- Please run create_simulation_parameters.R first."))
}
TRUE_PARAMS <- load_true_parameters(param_file_path)
cat("Parameters loaded successfully.\n")


# --- 3. Wrapper function to run filtering with multiple initializations ---
run_filtering <- function(seed, n_init, ...){
  best_lik <- -Inf
  best_result <- NULL
  cat(paste0("    -> Estimating model with ", n_init, " initializations...\n"))
  pb <- progress_bar$new(
    format = "      init :current/:total [:bar] :percent | Lik: :lik | ETA: :eta",
    total = n_init, width = 80, clear = FALSE
  )
  for(i in 1:n_init){
    # Run filtering, suppress output, and catch errors
    sink(tempfile())
    result <- tryCatch({filtering(seed = seed, init = i, ...)}, 
                       error = function(e) {cat(paste("Error in Init", i, ":", e$message, "\n")); return(NULL)})
    sink()
    lik_val <- if (!is.null(result$sumLik_best)) sprintf("%.2f", result$sumLik_best) else "Error"
    pb$tick(tokens = list(lik = lik_val))
    if(!is.null(result) && result$sumLik_best > best_lik){
      best_lik <- result$sumLik_best
      best_result <- result
    }
  }
  cat(paste0("\n    -> Best log-likelihood: ", sprintf("%.2f", best_lik), "\n"))
  return(best_result)
}


# --- 4. Run the simulation (with skip functionality) ---
cat("Starting simulation...\n")
N_IMPUTE <- if (exists("N_IMPUTE")) N_IMPUTE else 1

for (i in 1:length(N_VEC)) {
  N <- N_VEC[i]
  Nt <- NT_TOTAL 
  
  cat(paste0("\nSetting: N=", N, ", Nt=", Nt, " (Train on first ", NT_TRAIN, " points)\n"))
  
  for (seed in 1:N_SIM) {
    cat(paste0("\n--- Processing Seed ", seed, " ---\n"))
    
    # --- 4.1. Data Generation or Loading ---
    list_of_dfs <- list()
    
    data_files_exist <- all(sapply(1:1, function(m) {
      file.exists(file.path(OUTPUT_SIM_DATA_PATH, paste0("sim_data_seed_", seed, "_m_", m, ".RDS")))
    }))
    
    if (data_files_exist) {
      cat("  Data generation: Already complete. Loading saved files.\n")
      for (m in 1:N_IMPUTE) {
        list_of_dfs[[m]] <- readRDS(file.path(OUTPUT_SIM_DATA_PATH, paste0("sim_data_seed_", seed, "_m_", m, ".RDS")))
      }
    } else {
      cat(paste0("  Data generation: Running (m=", N_IMPUTE, ")...\n"))
      list_of_dfs <- generate_data(
        seed = seed, N = N, Nt = Nt, 
        O1 = O1, O2 = O2, U1 = U1,
        true_params = TRUE_PARAMS,
        N_IMPUTE = N_IMPUTE
      )
      
      cat("  Data generation: Complete. Saving each imputation data to file.\n")
      for (m in 1:N_IMPUTE) {
        saveRDS(list_of_dfs[[m]], file.path(OUTPUT_SIM_DATA_PATH, paste0("sim_data_seed_", seed, "_m_", m, ".RDS")))
      }
    }
    
    # --- 4.2. Model Estimation (filtering) ---
    cat("  Starting model estimation.\n")
    for (m in 1:N_IMPUTE) {
      
      filter_path <- file.path(OUTPUT_RESULTS_PATH, paste0("filter_sim_seed_", seed, "_m_", m, ".RDS"))
      
      if (file.exists(filter_path)) {
        cat(paste0("    Imputation ", m, "/", N_IMPUTE, ": Already complete. Skipping.\n"))
        next
      }
      
      cat(paste0("    Imputation ", m, "/", N_IMPUTE, ": Running estimation...\n"))
      
      current_df <- list_of_dfs[[m]]
      
      filter_results <- run_filtering(
        seed = seed, n_init = N_INIT, N = N, Nt = Nt, O1 = O1, O2 = O2, U1 = U1,
        y1 = current_df$y1, y2 = current_df$y2, DO = current_df$DO,
        maxIter = MAX_ITER, patience = PATIENCE, n_train = NT_TRAIN
      )
      
      # Calculate prediction accuracy metrics (Sensitivity, Specificity, Accuracy)
      if (!is.null(filter_results) && !is.null(filter_results$mPr_pred_best)) {
        prob_S2_t <- filter_results$mPr_pred_best[,, 2]
        pred_binary_t <- ifelse(prob_S2_t[, (NT_TRAIN + 2):(current_df$Nt + 1)] > 0.5, 1, 0)
        true_binary_t <- current_df$S_true[, (NT_TRAIN + 2):(current_df$Nt + 1)] - 1
        
        TP <- sum(pred_binary_t == 1 & true_binary_t == 1, na.rm = TRUE)
        TN <- sum(pred_binary_t == 0 & true_binary_t == 0, na.rm = TRUE)
        FP <- sum(pred_binary_t == 1 & true_binary_t == 0, na.rm = TRUE)
        FN <- sum(pred_binary_t == 0 & true_binary_t == 1, na.rm = TRUE)
        
        filter_results$sensitivity <- if ((TP + FN) > 0) TP / (TP + FN) else 0
        filter_results$specificity <- if ((TN + FP) > 0) TN / (TN + FP) else 0
        filter_results$accuracy <- if ((TP + TN + FP + FN) > 0) (TP + TN) / (TP + TN + FP + FN) else 0
        
      } else {
        filter_results$sensitivity <- NA
        filter_results$specificity <- NA
        filter_results$accuracy <- NA
      }
      
      saveRDS(filter_results, filter_path)
      cat(paste0("    Imputation ", m, "/", N_IMPUTE, ": Results saved.\n"))
    }
    
    rm(list_of_dfs)
    gc()
  }
}

cat("\nSimulation complete.\n")