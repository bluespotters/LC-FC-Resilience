# =====================================
# Cognitive trajectories: parcelwise
# =====================================

source(file.path(Sys.getenv("LCFC_REPO", "."), "R", "setup.R"))

opt <- parse_options(list(seed = "LC", parcellation = 200, control = "none"))
control <- null_if_none(opt$control)

run_name <- if (is.null(control)) opt$seed else paste0(opt$seed, "_control-", control)
out_dir <- file.path(RESULTS_DIR, "trajectories", run_name, paste0("parcels_", opt$parcellation))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

data <- prepare_trajectory_data(opt$seed, "parcel", parcellation = opt$parcellation, control = control)
roi_info <- read_roi_info(opt$parcellation)
cat(sprintf("%s | %d parcels | control: %s | %d participants\n",
            opt$seed, n_distinct(data$region), opt$control, n_distinct(data$sub)))

control_measures <- if (is.null(control)) character(0) else "control_value"
runs <- lapply(names(ANALYSES), function(a) list(name = a, analysis = a, controls = control_measures))
if (is.null(control)) {
  runs <- c(runs, list(list(name = "apoe_amyloid_controlled", analysis = "apoe", controls = "centiloid")))
}

all_results <- list()
for (r in runs) {
  cat("Running", r$name, "\n")
  fit <- run_regions(data, r$analysis, r$controls, n_cores = N_CORES, keep_models = FALSE)
  res <- add_roi_info(fit$results, roi_info) %>% mutate(model = r$name)
  all_results[[r$name]] <- res

  adir <- file.path(out_dir, r$name)
  dir.create(adir, showWarnings = FALSE)

  # Surface + subcortex maps of signed -log10(p)
  csv <- file.path(adir, "parcel_results.csv")
  write.csv(res %>% arrange(p_value), csv, row.names = FALSE)
  plot_surfaces(csv, "signed_neglog10p", "p_value_fdr", parcellation = opt$parcellation)
  cat(sprintf("  %d parcels with pFDR < 0.05\n", sum(res$p_value_fdr < 0.05, na.rm = TRUE)))
}

results <- bind_rows(all_results) %>%
  mutate(seed = opt$seed, parcellation = opt$parcellation, control = opt$control)
write.csv(results, file.path(out_dir, "results.csv"), row.names = FALSE)

cat("\nParcels with pFDR < 0.05 per model:\n")
print(results %>% group_by(model) %>% summarise(n_significant = sum(p_value_fdr < 0.05, na.rm = TRUE)))
flagged <- results %>% filter(!converged_ok)
if (nrow(flagged)) cat(sprintf("\n%d model pairs with convergence flags (see results.csv)\n", nrow(flagged)))
cat("Results written to", out_dir, "\n")
