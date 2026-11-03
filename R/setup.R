# ==============
# Shared setup 
# ==============

REPO_DIR <- normalizePath(Sys.getenv("LCFC_REPO", "."), mustWork = TRUE)
config_file <- file.path(REPO_DIR, "config", "config.R")
if (!file.exists(config_file)) {
  stop("Copy config/config.example.R to config/config.R and set your paths ",
       "(run scripts from the repository root, or set LCFC_REPO).")
}
source(config_file)                      # DATA_DIR, RESULTS_DIR, PYTHON, N_CORES

for (f in c("models.R", "data.R", "plots.R")) source(file.path(REPO_DIR, "R", f))

#' Parse --key=value command-line options over a list of defaults
#' e.g. Rscript script.R --seed=VTA --control=LC
parse_options <- function(defaults) {
  args <- commandArgs(trailingOnly = TRUE)
  for (a in args) {
    kv <- regmatches(a, regexec("^--([^=]+)=(.*)$", a))[[1]]
    if (length(kv) != 3 || !kv[2] %in% names(defaults)) stop("Unknown option: ", a)
    defaults[[kv[2]]] <- if (is.numeric(defaults[[kv[2]]])) as.numeric(kv[3]) else kv[3]
  }
  defaults
}

null_if_none <- function(x) if (is.null(x) || identical(x, "none") || identical(x, "")) NULL else x
