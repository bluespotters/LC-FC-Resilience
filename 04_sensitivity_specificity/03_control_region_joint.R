# ====================================================================
# LC-FC effects before vs after controlling for a control region
# ====================================================================

source(file.path(Sys.getenv("LCFC_REPO", "."), "R", "setup.R"))

opt <- parse_options(list(control = "VTA", parcellation = 200, axis_max = 50))
parcel_dir <- function(run) file.path(RESULTS_DIR, "trajectories", run, paste0("parcels_", opt$parcellation))
out_dir <- file.path(RESULTS_DIR, "specificity", opt$control, "joint_models")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
y_label <- c(shifted_LC = "with shifted-FC control", VTA = "with VTA-FC control")[[opt$control]]

before <- read.csv(file.path(parcel_dir("LC"), "results.csv"))
after  <- read.csv(file.path(parcel_dir(paste0("LC_control-", opt$control)), "results.csv"))
roi_info <- read_roi_info(opt$parcellation)

fdr_range <- function(x) if (length(x)) paste(sprintf("%.3f", range(x)), collapse = "-") else NA

comparison <- bind_rows(lapply(names(ANALYSES), function(m) {
  compare_before_after(before %>% filter(model == m), after %>% filter(model == m)) %>% mutate(model = m)
})) %>% add_roi_info(roi_info)
write.csv(comparison, file.path(out_dir, "comparison.csv"), row.names = FALSE)

summary <- comparison %>% group_by(model) %>% summarise(
  n_sig_before = sum(sig_before), n_sig_after = sum(sig_after),
  remained = sum(flow == "Remained Significant"), lost = sum(flow == "Lost Significance"),
  gained = sum(flow == "Gained Significance"),
  mean_chisq_pct_change = round(mean(chisq_pct_change[sig_before]), 1),
  lost_fdr_before = fdr_range(p_value_fdr_before[flow == "Lost Significance"]),
  lost_fdr_after  = fdr_range(p_value_fdr_after[flow == "Lost Significance"]),
  r_chisq = round(cor(chisq_before, chisq_after), 3), .groups = "drop")
write.csv(summary, file.path(out_dir, "flow_summary.csv"), row.names = FALSE)
print(as.data.frame(summary))

# ---- Significance-flow surface maps --------------------------------
flow_code <- c("Remained Significant" = 1, "Gained Significance" = 2, "Lost Significance" = 3)
for (m in names(ANALYSES)) {
  mdir <- file.path(out_dir, m); dir.create(mdir, showWarnings = FALSE)
  csv <- file.path(mdir, "flow.csv")
  write.csv(comparison %>% filter(model == m) %>% mutate(flow_value = unname(flow_code[flow])) %>% arrange(ROI_Number),
            csv, row.names = FALSE)
  plot_surfaces(csv, "flow_value", parcellation = opt$parcellation, vmin = 1, vmax = 3, colormap = "significance_flow")
}

# ---- Chi-sq scatter without vs with the control -------------------
fdr_threshold <- function(chisq, sig) if (!any(sig) || all(sig)) NA else (min(chisq[sig]) + max(chisq[!sig])) / 2
flow_colours <- c("Remained Non-Significant" = "#D9D9D9", "Remained Significant" = "#696969",
                  "Gained Significance" = "#008080", "Lost Significance" = "#FF6347")

for (m in names(ANALYSES)) {
  d <- comparison %>% filter(model == m) %>%
    mutate(flow = factor(flow, levels = names(flow_colours))) %>% arrange(flow)
  x_thr <- fdr_threshold(d$chisq_before, d$sig_before)
  y_thr <- fdr_threshold(d$chisq_after, d$sig_after)
  p <- ggplot(d, aes(chisq_before, chisq_after)) +
    geom_abline(slope = 1, intercept = 0, colour = "grey60", linewidth = 0.4) +
    { if (!is.na(x_thr)) geom_vline(xintercept = x_thr, linetype = "dashed", colour = "grey45", linewidth = 0.35) } +
    { if (!is.na(y_thr)) geom_hline(yintercept = y_thr, linetype = "dashed", colour = "grey45", linewidth = 0.35) } +
    geom_point(aes(fill = flow), shape = 21, colour = "white", stroke = 0.3, size = 2.4) +
    annotate("text", x = opt$axis_max, y = 0, hjust = 1, vjust = 0, size = 5.5,
             label = sprintf("r = %.2f", cor(d$chisq_before, d$chisq_after)), colour = "grey20") +
    scale_fill_manual(values = flow_colours, drop = FALSE, guide = "none") +
    scale_x_continuous(limits = c(0, opt$axis_max), expand = expansion(mult = 0.02)) +
    scale_y_continuous(limits = c(0, opt$axis_max), expand = expansion(mult = 0.02)) +
    coord_fixed() +
    labs(x = expression(chi^2 ~ "without control"), y = bquote(chi^2 ~ .(y_label))) +
    theme_classic(base_size = 10) +
    theme(axis.title = element_text(size = 16), axis.text = element_text(colour = "grey30", size = 12),
          axis.line = element_line(colour = "grey40", linewidth = 0.3),
          axis.ticks = element_line(colour = "grey40", linewidth = 0.3))
  ggsave(file.path(out_dir, paste0("chisq_scatter_", m, ".png")), p, width = 3.6, height = 3.6, dpi = 600, bg = "white")
}
cat("Results written to", out_dir, "\n")
