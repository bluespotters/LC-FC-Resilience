# Copy this file to config/config.R and edit the paths

# Folder containing the data
DATA_DIR <- "/path/to/derivatives"

# Folder for results
RESULTS_DIR <- "/path/to/results"

# Python interpreter with dependencies installed
PYTHON <- "python3"

# Parrallelization: number of cores to use for model fitting (1 = no parallelization)
N_CORES <- max(1, parallel::detectCores() - 1)
