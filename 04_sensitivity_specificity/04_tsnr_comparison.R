# ==========================================
# tSNR in the LC, shifted LC and VTA masks
# ==========================================

source(file.path(Sys.getenv("LCFC_REPO", "."), "R", "setup.R"))
suppressPackageStartupMessages(library(rstatix))

out_dir <- file.path(RESULTS_DIR, "tsnr")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
regions <- c(LC = "LC", shifted_LC = "Shifted LC", VTA = "VTA")
colours <- c("LC" = "#008080", "Shifted LC" = "#FF6347", "VTA" = "#6A5ACD")

wide <- read.csv(file.path(DATA_DIR, "tsnr.csv"))
long <- wide %>% pivot_longer(all_of(names(regions)), names_to = "region", values_to = "tsnr") %>%
  mutate(region = factor(regions[region], levels = regions), sub = factor(sub))

# ---- Descriptives -------------------------------------------------------------
desc <- long %>% group_by(region) %>%
  summarise(n = n(), mean = mean(tsnr), sd = sd(tsnr), .groups = "drop")
write.csv(desc, file.path(out_dir, "tsnr_descriptives.csv"), row.names = FALSE)

# ---- Omnibus: repeated-measures ANOVA -----------------------------------------
aov <- anova_test(data = long, dv = tsnr, wid = sub, within = region, effect.size = "ges")
tab <- get_anova_table(aov, correction = "auto")
omnibus <- data.frame(F = tab$F, DFn = tab$DFn, DFd = tab$DFd, p = tab$p, ges = tab$ges,
                      mauchly_p = aov$`Mauchly's Test for Sphericity`$p,
                      GG_epsilon = aov$`Sphericity Corrections`$GGe, row.names = NULL)
write.csv(omnibus, file.path(out_dir, "tsnr_omnibus.csv"), row.names = FALSE)

# ---- Paired comparisons (Holm-corrected across the three pairs) ---------------
pairs <- combn(names(regions), 2, simplify = FALSE)
tests <- bind_rows(lapply(pairs, function(pr) {
  a <- wide[[pr[1]]]; b <- wide[[pr[2]]]; d <- a - b
  t <- t.test(a, b, paired = TRUE)
  data.frame(comparison = paste(regions[pr], collapse = " vs "),
             mean_diff = mean(d), ci_low = t$conf.int[1], ci_high = t$conf.int[2],
             t = unname(t$statistic), df = unname(t$parameter),
             p = t$p.value, cohens_dz = mean(d) / sd(d))
})) %>% mutate(p_holm = p.adjust(p, method = "holm"))
write.csv(tests, file.path(out_dir, "tsnr_paired_tests.csv"), row.names = FALSE)

print(desc); print(omnibus); print(tests)

# ---- Density plot ---------------------------------------------------------------
p <- ggplot(long, aes(tsnr, colour = region, fill = region)) +
  geom_density(alpha = 0.18, linewidth = 0.8) +
  geom_vline(data = desc, aes(xintercept = mean, colour = region), linetype = "dashed",
             linewidth = 0.5, show.legend = FALSE) +
  scale_colour_manual(values = colours, name = NULL) +
  scale_fill_manual(values = colours, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Voxelwise tSNR", y = "Density") +
  theme_classic(base_size = 14) +
  theme(axis.title = element_text(size = 16), axis.text = element_text(colour = "grey30", size = 12),
        legend.position = c(0.85, 0.85), legend.text = element_text(size = 13))
ggsave(file.path(out_dir, "tsnr_density.png"), p, width = 6, height = 4.2, dpi = 600, bg = "white")
