# ==============================================================
# Are LC-FC effects independent of local tissue properties?
# ==============================================================

source(file.path(Sys.getenv("LCFC_REPO", "."), "R", "setup.R"))
suppressPackageStartupMessages(library(ggforce))

controls <- c("ct_value", "wd_value")
out_dir <- file.path(RESULTS_DIR, "sensitivity", "tissue")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- cortical parcels only, LC-FC + CT + WD per participant x parcel ----
read_cortical <- function(file, value, offset = 0) {
  read.csv(file.path(DATA_DIR, file)) %>%
    select(sub, matches("^ROI_[0-9]+$")) %>%
    pivot_longer(-sub, names_to = "region", values_to = value) %>%
    mutate(ROI_Number = as.numeric(sub("ROI_", "", region)) + offset) %>%
    filter(ROI_Number > 32, ROI_Number <= 232, !is.na(.data[[value]])) %>%
    select(-region)
}
tissue <- read_cortical("200_LC_ROI_correlations.csv", "fc") %>%
  inner_join(read_cortical("cortical_thickness_200.csv", "ct_value", offset = 32), by = c("sub", "ROI_Number")) %>%
  inner_join(read_cortical("weighted_degree_200.csv", "wd_value", offset = 32), by = c("sub", "ROI_Number")) %>%
  mutate(region = paste0("ROI_", ROI_Number))

data <- read.csv(file.path(DATA_DIR, "longitudinal_data_processed.csv")) %>%
  inner_join(tissue, by = "sub", relationship = "many-to-many") %>%
  filter(amyloid_status == "amyloid+") %>%
  mutate(QSDY_years = QSDY / 365.25,
         sex = factor(sex, levels = c("F", "M")),
         APOE_status = factor(APOE_status, levels = c("E4-", "E4+")),
         demo_group = factor(paste(sex, APOE_status, sep = "_"), levels = c("F_E4-", "F_E4+", "M_E4-", "M_E4+")))
roi_info <- read_roi_info(200)
cat(sprintf("%d cortical parcels | %d participants | controls: %s\n",
            n_distinct(data$region), n_distinct(data$sub), paste(controls, collapse = ", ")))

# ---- Before vs after tissue control -------------------------------------------
comparison <- bind_rows(lapply(names(ANALYSES), function(m) {
  cat("Running", m, "\n")
  before <- run_regions(data, m, character(0), n_cores = N_CORES, keep_models = FALSE)$results
  after  <- run_regions(data, m, controls,     n_cores = N_CORES, keep_models = FALSE)$results
  compare_before_after(before, after) %>% mutate(model = m)
})) %>% add_roi_info(roi_info)
write.csv(comparison, file.path(out_dir, "comparison.csv"), row.names = FALSE)

summary <- comparison %>% group_by(model) %>% summarise(
  n_sig_before = sum(sig_before), n_sig_after = sum(sig_after),
  remained = sum(flow == "Remained Significant"), lost = sum(flow == "Lost Significance"),
  gained = sum(flow == "Gained Significance"),
  mean_chisq_pct_change = round(mean(chisq_pct_change[sig_before]), 1),
  mean_added_marginal_r2 = round(mean(marginal_r2_change[sig_before]), 4),   # variance added by CT + WD
  .groups = "drop")
write.csv(summary, file.path(out_dir, "flow_summary.csv"), row.names = FALSE)
print(as.data.frame(summary))

# ---- significance-flow maps and alluvial plots by Yeo-17 network --------------
flow_code <- c("Remained Significant" = 1, "Gained Significance" = 2, "Lost Significance" = 3)
network_colours <- c(VisCent = "#4B0082", VisPeri = "#FF0000", SomMotA = "#87CEEB", SomMotB = "#00CED1",
                     DorsAttnA = "#228B22", DorsAttnB = "#006400", SalVentAttnA = "#FF69B4",
                     SalVentAttnB = "#FFB6C1", LimbicA = "#C0C0C0", LimbicB = "#808080", ContA = "#DAA520",
                     ContB = "#8B4513", ContC = "#A0522D", DefaultA = "#FFD700", DefaultB = "#FF1493",
                     DefaultC = "#000080", TempPar = "#0000FF")

for (m in names(ANALYSES)) {
  mdir <- file.path(out_dir, m); dir.create(mdir, showWarnings = FALSE)
  d <- comparison %>% filter(model == m)
  csv <- file.path(mdir, "flow.csv")
  write.csv(d %>% mutate(flow_value = unname(flow_code[flow])) %>% arrange(ROI_Number), csv, row.names = FALSE)
  plot_surfaces(csv, "flow_value", vmin = 1, vmax = 3, colormap = "significance_flow")

  sets <- d %>%
    mutate(flow = sub("Remained Non-Significant", "Remained\nNon-Significant", flow)) %>%
    count(Network, flow, name = "n") %>%
    mutate(Network = factor(Network, levels = unique(Network)), flow = factor(flow, levels = unique(flow))) %>%
    gather_set_data(1:2)
  axes <- unique(as.character(sets$x))                 # axis ids differ between ggforce versions
  sets$x <- factor(as.character(sets$x), levels = axes)
  labs_pos <- sets %>% distinct(x, y) %>% mutate(h = as.numeric(x == axes[1]), nudge = ifelse(x == axes[1], -0.15, 0.05))
  p <- ggplot(sets, aes(x = x, id = id, value = n, split = y)) +
    geom_parallel_sets(aes(fill = Network), axis.width = -0.05, alpha = 0.5) +
    geom_parallel_sets_axes(axis.width = 0.05, fill = "transparent", colour = "black") +
    geom_parallel_sets_labels(angle = 0, size = 3, lineheight = 0.8, hjust = labs_pos$h,
                              position = position_nudge(x = labs_pos$nudge)) +
    scale_fill_manual(values = network_colours, guide = "none") +
    scale_x_discrete(labels = c("Brain Network", "Tissue Control Result")) +
    theme_void()
  ggsave(file.path(mdir, "alluvial.png"), p, width = 5, height = 7, dpi = 300, bg = "white")
}

# ---- Mean WD and CT per parcel ---------------------------------------------------
means <- data %>% distinct(sub, ROI_Number, ct_value, wd_value) %>%
  group_by(ROI_Number) %>% summarise(ct_mean = mean(ct_value), wd_mean = mean(wd_value), .groups = "drop")
mdir <- file.path(out_dir, "tissue_means"); dir.create(mdir, showWarnings = FALSE)
csv <- file.path(mdir, "tissue_means.csv"); write.csv(means, csv, row.names = FALSE)
plot_surfaces(csv, "ct_mean", vmin = min(means$ct_mean), vmax = max(means$ct_mean), colormap = "Purples")
plot_surfaces(csv, "wd_mean", vmin = min(means$wd_mean), vmax = max(means$wd_mean), colormap = "Greens")
cat("Results written to", out_dir, "\n")
