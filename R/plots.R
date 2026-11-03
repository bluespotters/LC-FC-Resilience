# ====================================================================
# Plotting helpers: trajectory, network lollipop, brain surfaces
# ====================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(sjPlot)
  library(dplyr)
})

COLOURS <- c(low = "#008080", mid = "darkgrey", high = "#FF6347")   # teal / grey / tomato

# --------------------------------------------------------------------
# Model-predicted PACC trajectories at -1SD, mean and +1SD of connectivity
# --------------------------------------------------------------------

#' @param model fitted full model (from run_regions)
#' @param secondary_var moderator to facet by ("centiloid", "sex", "APOE_status", "demo_group")
#' @param fdr_p FDR-corr pval
#' @param is_global if TRUE the title shows only  pval
trajectory_plot <- function(model, label, secondary_var = NULL, fdr_p = NULL,
                            is_global = FALSE, y_limits = NULL, x_breaks = c(0, 8)) {
  terms <- c("QSDY_years [all]", "fc [meansd]", secondary_var)

  p <- tryCatch(
    plot_model(model, type = "pred", terms = terms, axis.title = c("Time (years)", "PACC")),
    error = function(e) ggplot() + annotate("text", x = 0, y = 0, label = conditionMessage(e)) + theme_void()
  )

  title <- if (is.null(fdr_p)) label
           else if (is_global) bquote(p[FDR] == .(sprintf("%.3f", fdr_p)))
           else bquote(.(label) ~ ":" ~ p[FDR] == .(sprintf("%.3f", fdr_p)))

  p <- p +
    scale_color_manual(values = unname(COLOURS), labels = c("-1SD", "mean", "+1SD"), name = label) +
    scale_fill_manual(values = unname(COLOURS), labels = c("-1SD", "mean", "+1SD"), name = label) +
    scale_x_continuous(breaks = x_breaks) +
    labs(title = title) +
    theme_minimal() +
    theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
          axis.title = element_text(size = 18, face = "bold"),
          axis.text = element_text(size = 12),
          legend.title = element_text(size = 16, face = "bold"),
          legend.text = element_text(size = 14),
          legend.position = "bottom",
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black", linewidth = 0.8),
          axis.ticks = element_line(colour = "black", linewidth = 0.3),
          axis.ticks.length = unit(0.15, "cm"),
          plot.margin = margin(10, 10, 10, 10))

  if (!is.null(y_limits)) p <- p + scale_y_continuous(limits = y_limits)
  p
}

# --------------------------------------------------------------------
# Network lollipop (radial) plot of signed -log10(p)
# --------------------------------------------------------------------

YEO7_COLOURS <- c(Global = "#808080", Visual = "#781286", SomMot = "#4682B4", DorsAttn = "#00760E",
                  VentAttn = "#C43AFA", Limbic = "#DCF8A4", Cont = "#E69422", Default = "#CD3E4E")

LOLLIPOP_LAYOUT <- list(
  yeo17 = list(
    order = c("global", "VisCent", "VisPeri", "SomMotA", "SomMotB", "DorsAttnA", "DorsAttnB",
              "SalVentAttnA", "SalVentAttnB", "LimbicA", "LimbicB", "ContA", "ContB", "ContC",
              "DefaultA", "DefaultB", "DefaultC", "TempPar"),
    parent = c(global = "Global", VisCent = "Visual", VisPeri = "Visual", SomMotA = "SomMot",
               SomMotB = "SomMot", DorsAttnA = "DorsAttn", DorsAttnB = "DorsAttn",
               SalVentAttnA = "VentAttn", SalVentAttnB = "VentAttn", LimbicA = "Limbic",
               LimbicB = "Limbic", ContA = "Cont", ContB = "Cont", ContC = "Cont",
               DefaultA = "Default", DefaultB = "Default", DefaultC = "Default", TempPar = "Default"),
    label = c(global = "global", VisCent = "Cent", VisPeri = "Peri", SomMotA = "A", SomMotB = "B",
              DorsAttnA = "A", DorsAttnB = "B", SalVentAttnA = "A", SalVentAttnB = "B",
              LimbicA = "A", LimbicB = "B", ContA = "A", ContB = "B", ContC = "C",
              DefaultA = "A", DefaultB = "B", DefaultC = "C", TempPar = "TP")),
  yeo7 = list(
    order = c("global", "Vis", "SomMot", "DorsAttn", "SalVentAttn", "Limbic", "Cont", "Default"),
    parent = c(global = "Global", Vis = "Visual", SomMot = "SomMot", DorsAttn = "DorsAttn",
               SalVentAttn = "VentAttn", Limbic = "Limbic", Cont = "Cont", Default = "Default"),
    label = NULL)
)

#' @param results run_regions()$results for the global + network regions
#' @param value_range fixed reference range for the radial scale
lollipop_plot <- function(results, networks = "yeo17", value_range = c(-5, 5)) {
  lay <- LOLLIPOP_LAYOUT[[networks]]
  scaling <- max(abs(value_range))

  d <- results %>%
    transmute(network = factor(as.character(region), levels = lay$order),
              value = signed_neglog10p, p_value_fdr) %>%
    filter(!is.na(network), !is.na(value)) %>%
    arrange(network) %>%
    mutate(parent = lay$parent[as.character(network)],
           colour = YEO7_COLOURS[parent],
           position = 0.55 + value / scaling * 0.55 * 0.7,    # 0.55 = zero baseline
           significant = p_value_fdr < 0.05,
           fill = if_else(significant, colour, "#F0F0F0"),
           size = if_else(significant, 6.25, 5),
           marker = if_else(significant, "*", ""))

  p <- ggplot(d, aes(x = network)) +
    geom_hline(yintercept = 1.0, colour = "grey90", linewidth = 0.8) +
    geom_hline(yintercept = 0, colour = "black", linewidth = 0.8) +
    geom_hline(yintercept = 0.55, colour = "black", linewidth = 0.8) +
    geom_segment(aes(xend = network, y = 0, yend = 1.0, colour = parent), linewidth = 1, alpha = 0.15) +
    scale_colour_manual(values = YEO7_COLOURS) +
    geom_segment(aes(xend = network, y = 0.55, yend = position), colour = "black", linewidth = 1, lineend = "round") +
    geom_point(aes(y = position), size = d$size, fill = d$fill, shape = 21, colour = "black", stroke = 0.5) +
    geom_text(aes(y = position, label = marker), size = 4, fontface = "bold") +
    coord_polar() +
    scale_y_continuous(limits = c(0, 1.8), breaks = NULL) +
    theme_void() +
    theme(plot.background = element_rect(fill = "white", colour = NA), legend.position = "none")

  # coloured bands grouping networks by their Yeo-7 parent
  runs <- rle(d$parent)
  ends <- cumsum(runs$lengths); starts <- ends - runs$lengths + 1
  band <- if (networks == "yeo7") c(1.1, 1.2) else c(1.3, 1.45)
  for (i in seq_along(starts)) {
    p <- p + annotate("rect", xmin = starts[i] - 0.475, xmax = ends[i] + 0.475,
                      ymin = band[1], ymax = band[2], fill = YEO7_COLOURS[[runs$values[i]]],
                      colour = "black", linewidth = 0.25)
  }
  if (!is.null(lay$label)) {
    for (i in seq_len(nrow(d))) {
      p <- p + annotate("text", x = i, y = 1.15, label = lay$label[[as.character(d$network[i])]],
                        colour = d$colour[i], fontface = "bold", size = 4)
    }
  }
  for (lab in list(c(0.05, value_range[1]), c(0.5, 0), c(0.95, value_range[2]))) {
    p <- p + annotate("text", x = 0.5, y = lab[1], label = as.character(lab[2]), fontface = "bold", size = 3.5)
  }
  p
}

# --------------------------------------------------------------------
# Cortical + subcortical surface maps (python/plot_surfaces.py)
# -------------------------------------------------------------------

#' @param csv CSV with ROI_Number and the value column; images are written next to it
#' @param ... further plotter arguments, e.g. "--colormap", "custom"
plot_surfaces <- function(csv, value_col, fdr_col = NULL, parcellation = 200,
                          vmin = -5, vmax = 5, colormap = "teal_tomato", ...) {
  args <- c(file.path(REPO_DIR, "python", "plot_surfaces.py"), csv, value_col, fdr_col,
            "--parcellation", parcellation, "--vmin", vmin, "--vmax", vmax,
            "--colormap", colormap, ...)
  out <- system2(PYTHON, args = args, stdout = TRUE, stderr = TRUE)
  if (!is.null(attr(out, "status"))) warning("Surface plotting failed for ", csv, ":\n", paste(tail(out, 5), collapse = "\n"))
  invisible(out)
}
