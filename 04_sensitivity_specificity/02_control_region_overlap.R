# ==================================================
# Overlap between LC-FC and control-region effects
# ==================================================

source(file.path(Sys.getenv("LCFC_REPO", "."), "R", "setup.R"))
suppressPackageStartupMessages({ library(VennDiagram); library(grid) })
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

opt <- parse_options(list(control = "shifted_LC", parcellation = 200))
label <- c(shifted_LC = "Shifted", VTA = "VTA")[[opt$control]]
parcel_dir <- function(seed) file.path(RESULTS_DIR, "trajectories", seed, paste0("parcels_", opt$parcellation))
out_dir <- file.path(RESULTS_DIR, "specificity", opt$control, "overlap")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

lc  <- read.csv(file.path(parcel_dir("LC"), "results.csv"))
ctl <- read.csv(file.path(parcel_dir(opt$control), "results.csv"))

venn <- function(lc_sig, ctl_sig, n_total, file) {
  png(file, width = 6, height = 5.5, units = "in", res = 300)
  par(mar = c(1, 1, 2, 1)); plot.new(); plot.window(xlim = c(0, 10), ylim = c(0, 11))
  rect(0.5, 0.5, 9.5, 9.5, border = "black", lwd = 2)
  pushViewport(viewport(x = 0.5, y = 0.45, width = 0.8, height = 0.7))
  v <- venn.diagram(list(LC = lc_sig, Control = ctl_sig), filename = NULL, category.names = c("", ""),
                    fill = unname(COLOURS[c("low", "high")]), col = unname(COLOURS[c("low", "high")]),
                    alpha = 0.7, lwd = 3, cex = 1.8, fontfamily = "sans", fontface = "bold", cat.cex = 0,
                    margin = 0.05, rotation.degree = if (length(ctl_sig) > length(lc_sig)) 180 else 0)
  grid.draw(v)
  grid.text("LC", x = 0.15, y = 0.9, gp = gpar(col = COLOURS[["low"]], fontsize = 20, fontface = "bold"))
  grid.text(label, x = 0.85, y = 0.9, gp = gpar(col = COLOURS[["high"]], fontsize = 20, fontface = "bold"))
  popViewport()
  text(8.8, 1.2, n_total - length(union(lc_sig, ctl_sig)), cex = 1.2, font = 2)   # neither
  dev.off()
}

summary <- bind_rows(lapply(names(ANALYSES), function(m) {
  a <- lc  %>% filter(model == m)
  b <- ctl %>% filter(model == m)
  lc_sig  <- a$region[a$p_value_fdr < 0.05]
  ctl_sig <- b$region[b$p_value_fdr < 0.05]
  venn(lc_sig, ctl_sig, nrow(a), file.path(out_dir, paste0("venn_", m, ".png")))
  data.frame(model = m, n_parcels = nrow(a), lc_significant = length(lc_sig),
             control_significant = length(ctl_sig), overlap = length(intersect(lc_sig, ctl_sig)),
             neither = nrow(a) - length(union(lc_sig, ctl_sig)))
}))
write.csv(summary, file.path(out_dir, "overlap_summary.csv"), row.names = FALSE)
print(summary, row.names = FALSE)
