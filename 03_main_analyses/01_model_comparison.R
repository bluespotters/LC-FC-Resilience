# =================================================
# Linear vs quadratic vs natural spline (df = 2)
# =================================================

source(file.path(Sys.getenv("LCFC_REPO", "."), "R", "setup.R"))

pacc <- read.csv(file.path(DATA_DIR, "longitudinal_data_processed.csv")) %>%
  filter(amyloid_status == "amyloid+") %>%
  mutate(QSDY_years = QSDY / 365.25) %>%
  filter(!is.na(PACC), !is.na(QSDY_years))
cat(sprintf("Abeta+ participants: %d | observations: %d\n", n_distinct(pacc$sub), nrow(pacc)))

ctrl <- lmerControl(optimizer = "bobyqa")
models <- list(
  Linear    = lmer(PACC ~ QSDY_years + (1 + QSDY_years | sub), data = pacc, control = ctrl, REML = FALSE),
  Quadratic = lmer(PACC ~ QSDY_years + I(QSDY_years^2) + (1 + QSDY_years | sub), data = pacc, control = ctrl, REML = FALSE),
  Spline    = lmer(PACC ~ ns(QSDY_years, df = 2) + (1 + QSDY_years | sub), data = pacc, control = ctrl, REML = FALSE)
)

fit <- data.frame(model = names(models),
                  df  = sapply(models, function(m) attr(logLik(m), "df")),
                  AIC = sapply(models, AIC),
                  BIC = sapply(models, BIC), row.names = NULL) %>%
  mutate(delta_AIC = AIC - min(AIC), delta_BIC = BIC - min(BIC))

print(fit, row.names = FALSE, digits = 6)
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
write.csv(fit, file.path(RESULTS_DIR, "spline_model_comparison.csv"), row.names = FALSE)
