# =====================================================
# Cognitive trajectories: global and network scale
# =====================================================

source(file.path(Sys.getenv("LCFC_REPO", "."), "R", "setup.R"))
suppressPackageStartupMessages(library(patchwork))

opt <- parse_options(list(seed = "LC", networks = "yeo17", control = "none"))
control <- null_if_none(opt$control)
Y_LIMITS <- c(-13, 4)    # shared y-axis

run_name <- if (is.null(control)) opt$seed else paste0(opt$seed, "_control-", control)
out_dir <- file.path(RESULTS_DIR, "trajectories", run_name, opt$networks)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

data <- prepare_trajectory_data(opt$seed, "network", networks = opt$networks, control = control)
cat(sprintf("%s | %s | control: %s | %d participants, %d observations per region\n",
            opt$seed, opt$networks, opt$control, n_distinct(data$sub), nrow(data) / n_distinct(data$region)))

control_measures <- if (is.null(control)) character(0) else "control_value"
seed_label <- paste0(sub("_LC$", "", opt$seed), "-FC")

# Main models & the e4 model controlling for AB x time
runs <- lapply(names(ANALYSES), function(a) list(name = a, analysis = a, controls = control_measures))
if (is.null(control)) {
  runs <- c(runs, list(list(name = "apoe_amyloid_controlled", analysis = "apoe", controls = "centiloid")))
}

all_results <- list()
for (r in runs) {
  cat("Running", r$name, "\n")
  fit <- run_regions(data, r$analysis, r$controls, n_cores = N_CORES)
  res <- fit$results %>% mutate(model = r$name)
  all_results[[r$name]] <- res

  adir <- file.path(out_dir, r$name)
  dir.create(adir, showWarnings = FALSE)
  secondary <- ANALYSES[[r$analysis]]$secondary_var
  fdr_of <- function(region) res$p_value_fdr[res$region == region]

  if ("global" %in% names(fit$models)) {
    ggsave(file.path(adir, "global_trajectory.png"),
           trajectory_plot(fit$models[["global"]], seed_label, secondary, fdr_of("global"),
                           is_global = TRUE, y_limits = Y_LIMITS),
           width = 4, height = 4, dpi = 600)
  }

  ggsave(file.path(adir, "network_lollipop.png"), lollipop_plot(res, opt$networks),
         width = 5.5, height = 5.5, dpi = 600, bg = "white")

  networks <- setdiff(names(fit$models), "global")
  panels <- lapply(networks, function(n) {
    trajectory_plot(fit$models[[n]], n, secondary, fdr_of(n), y_limits = Y_LIMITS)
  })
  if (length(panels)) {
    ggsave(file.path(adir, "network_trajectories.png"),
           wrap_plots(panels, ncol = 4) & theme(legend.position = "bottom"),
           width = 15, height = 15, dpi = 600)
  }
}

results <- bind_rows(all_results) %>%
  mutate(seed = opt$seed, networks = opt$networks, control = opt$control)
write.csv(results, file.path(out_dir, "results.csv"), row.names = FALSE)

cat("\nSignificant (pFDR < 0.05):\n")
print(results %>% filter(p_value_fdr < 0.05) %>%
        select(model, region, chisq, df, p_value, p_value_fdr), row.names = FALSE)
flagged <- results %>% filter(!converged_ok)
if (nrow(flagged)) {
  cat("\nConvergence flags (gradient > 0.001, singular or non-converged):\n")
  print(flagged %>% select(model, region, full_gradient_max, reduced_gradient_max), row.names = FALSE)
}
cat("Results written to", out_dir, "\n")
