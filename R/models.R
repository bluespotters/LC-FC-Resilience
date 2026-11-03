# ====================================================
# Mixed-effects models for PACC trajectories (LRT)
# ====================================================

suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
  library(splines)
  library(dplyr)
  library(MuMIn)       
  library(parallel)
})

# Main models
ANALYSES <- list(
  time     = list(label = "FCxTIME",           secondary_var = NULL),
  amyloid  = list(label = "FCxAMYLOIDxTIME",   secondary_var = "centiloid"),
  sex      = list(label = "FCxSEXxTIME",       secondary_var = "sex"),
  apoe     = list(label = "FCxAPOExTIME",      secondary_var = "APOE_status"),
  sex_apoe = list(label = "FCxSEXxAPOExTIME",  secondary_var = "demo_group")
)

TIME_TERM <- "ns(QSDY_years, df = 2)"

# --------------------------------------------------------------------
# Model formulas
# --------------------------------------------------------------------

#' Full and reduced formulas for one analysis
#'
#' @param analysis one of names(ANALYSES)
#' @param control_measures character vector of columns entered as <measure> x time
#' @return list(full, reduced, interaction_pattern)
model_formulas <- function(analysis, control_measures = character(0)) {
  if (!analysis %in% names(ANALYSES)) {
    stop("analysis must be one of: ", paste(names(ANALYSES), collapse = ", "))
  }

  tt <- TIME_TERM
  random_effects <- "(1 + QSDY_years | sub)"
  base <- "YrsEd + fd + Cumulative_Dose_Scaled"
  controls <- if (length(control_measures)) paste0(paste(control_measures, "*", tt, collapse = " + "), " + ") else ""
  age_t  <- paste0("age * ", tt, " + ")
  sex_t  <- paste0("sex * ", tt, " + ")
  apoe_t <- paste0("APOE_status * ", tt, " + ")
  rest   <- paste0(base, " + ", random_effects)

  # Moderator for the 3-/4-way interactions
  spec <- switch(analysis,
    time     = list(mod = NULL,          covs = paste0(age_t, sex_t, apoe_t)),
    amyloid  = list(mod = "centiloid",   covs = paste0(age_t, sex_t, apoe_t)),
    sex      = list(mod = "sex",         covs = paste0(age_t, apoe_t)),
    apoe     = list(mod = "APOE_status", covs = paste0(age_t, sex_t)),
    sex_apoe = list(mod = "demo_group",  covs = age_t)
  )

  if (is.null(spec$mod)) {
    full    <- paste0("PACC ~ fc * ", tt, " + ", controls, spec$covs, rest)
    reduced <- paste0("PACC ~ fc + ", tt, " + ", controls, spec$covs, rest)
    pattern <- "^fc:ns\\(QSDY_years"
  } else {
    m <- spec$mod
    full    <- paste0("PACC ~ fc * ", tt, " * ", m, " + ", controls, spec$covs, rest)
    reduced <- paste0("PACC ~ fc * ", tt, " + fc * ", m, " + ", tt, " * ", m, " + ",
                      controls, spec$covs, rest)
    pattern <- paste0("^fc:ns\\(QSDY_years.*:", m)
  }

  list(full = as.formula(full), reduced = as.formula(reduced), interaction_pattern = pattern)
}

# --------------------------------------------------------------------
# Diagnostics
# --------------------------------------------------------------------

convergence_diagnostics <- function(model) {
  grad <- model@optinfo$derivs$gradient
  msgs <- model@optinfo$conv$lme4$messages
  list(
    converged    = model@optinfo$conv$opt == 0,
    singular     = isSingular(model),
    gradient_max = if (is.null(grad)) NA else max(abs(grad)),
    messages     = if (length(msgs)) paste(msgs, collapse = "; ") else "No messages"
  )
}

#' Marginal R2 of the full model
marginal_r_squared <- function(model) {
  tryCatch(r.squaredGLMM(model)[1, "R2m"], error = function(e) {
    warning("R2 calculation failed: ", conditionMessage(e))
    NA
  })
}

# --------------------------------------------------------------------
# One region
# --------------------------------------------------------------------

#' Fit full and reduced models for one region and run the LRT
#'
#' @return list(result = one-row data frame, model = fitted full model or NULL)
fit_region <- function(data, region, analysis, control_measures = character(0)) {
  f <- model_formulas(analysis, control_measures)
  environment(f$full) <- environment(f$reduced) <- environment()
  ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))

  base_row <- data.frame(region = region, analysis = analysis,
                         controls = paste(control_measures, collapse = ","),
                         stringsAsFactors = FALSE)

  tryCatch({
    full_model    <- lmer(f$full,    data = data, control = ctrl)
    reduced_model <- lmer(f$reduced, data = data, control = ctrl)
    lrt <- anova(reduced_model, full_model)        # refits both with ML
    p   <- lrt$`Pr(>Chisq)`[2]

    # Direction of the effect = sign of the second (later) spline coefficient
    coefs  <- fixef(full_model)
    terms  <- grep(f$interaction_pattern, names(coefs), value = TRUE)
    s2 <- if (length(terms) >= 2) coefs[[terms[2]]] else NA

    fc_conv <- convergence_diagnostics(full_model)
    rc_conv <- convergence_diagnostics(reduced_model)

    result <- base_row %>% mutate(
      chisq = round(lrt$Chisq[2], 2),
      df = lrt$Df[2],
      p_value = p,
      signed_neglog10p = -log10(p) * sign(s2),
      marginal_r2 = round(marginal_r_squared(full_model), 3),
      full_converged = fc_conv$converged, reduced_converged = rc_conv$converged,
      full_singular = fc_conv$singular,   reduced_singular = rc_conv$singular,
      full_gradient_max = round(fc_conv$gradient_max, 6),
      reduced_gradient_max = round(rc_conv$gradient_max, 6),
      full_messages = fc_conv$messages, reduced_messages = rc_conv$messages,
      converged_ok = fc_conv$converged & rc_conv$converged & !fc_conv$singular & !rc_conv$singular &
        fc_conv$gradient_max <= 0.001 & rc_conv$gradient_max <= 0.001,
      analysis_success = TRUE
    )
    list(result = result, model = full_model)

  }, error = function(e) {
    list(result = base_row %>% mutate(chisq = NA, df = NA, p_value = NA,
                                      signed_neglog10p = NA, converged_ok = FALSE,
                                      analysis_success = FALSE,
                                      full_messages = paste("ERROR:", conditionMessage(e))),
         model = NULL)
  })
}

# --------------------------------------------------------------------
# All regions of one model
# --------------------------------------------------------------------

#' Run one analysis across all regions with FDR correction across regions
#'
#' @param data long format df with one row per participant x visit x region
#' @param analysis one of names(ANALYSES)
#' @param control_measures columns entered as <measure> x time in full and reduced models
#' @param n_cores regions are fitted in parallel (forked)
#' @param keep_models return the fitted full models (needed for trajectory plots)
#' @return list(results = data frame, models = named list of full models)
run_regions <- function(data, analysis, control_measures = character(0),
                        n_cores = 1, keep_models = TRUE) {
  regions <- unique(as.character(data$region))
  fits <- mclapply(regions, function(r) {
    fit_region(data[data$region == r, ], r, analysis, control_measures)
  }, mc.cores = n_cores)
  names(fits) <- regions

  results <- bind_rows(lapply(fits, `[[`, "result"))
  results$p_value_fdr <- NA_real_
  ok <- !is.na(results$p_value)
  results$p_value_fdr[ok] <- p.adjust(results$p_value[ok], method = "BY")

  models <- if (keep_models) Filter(Negate(is.null), lapply(fits, `[[`, "model")) else NULL
  list(results = results, models = models)
}

# --------------------------------------------------------------------
# Before vs after a control
# --------------------------------------------------------------------

#' Classify each region 
significance_flow <- function(before_sig, after_sig) {
  case_when(
    before_sig &  after_sig ~ "Remained Significant",
    before_sig & !after_sig ~ "Lost Significance",
    !before_sig & after_sig ~ "Gained Significance",
    TRUE ~ "Remained Non-Significant"
  )
}

#' Join uncontrolled and controlled results region-by-region & flow labels
compare_before_after <- function(before, after) {
  keep <- c("region", "chisq", "p_value", "p_value_fdr", "signed_neglog10p", "marginal_r2")
  inner_join(before[keep], after[keep], by = "region", suffix = c("_before", "_after")) %>%
    mutate(
      sig_before = !is.na(p_value_fdr_before) & p_value_fdr_before < 0.05,
      sig_after  = !is.na(p_value_fdr_after)  & p_value_fdr_after  < 0.05,
      flow = significance_flow(sig_before, sig_after),
      chisq_pct_change = round((chisq_after - chisq_before) / chisq_before * 100, 1),
      marginal_r2_change = marginal_r2_after - marginal_r2_before
    )
}
