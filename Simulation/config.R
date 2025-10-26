# config.R

# Execution Settings
N_VEC <- 50
NT_TOTAL <- 30           
NT_TRAIN <- 15            
N_SIM <- 100
N_INIT <- 3
MAX_ITER <- 1000
PATIENCE <- 5

# Model Dimensions
O1 <- 6 
O2 <- 3
U1 <- 2  

# Path Settings
library(here)
OUTPUT_SIM_DATA_PATH <- here("output", "sim_data")
OUTPUT_RESULTS_PATH <- here("output", "results")

# Create directories if they do not exist
dir.create(OUTPUT_SIM_DATA_PATH, recursive = TRUE, showWarnings = FALSE)
dir.create(OUTPUT_RESULTS_PATH, recursive = TRUE, showWarnings = FALSE)