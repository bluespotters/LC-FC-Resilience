# ===========================
# Data loading and merging
# ==========================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

NETWORK_ORDER <- list(
  yeo17 = c("global", "ContA", "ContB", "ContC", "DefaultA", "DefaultB", "DefaultC",
            "DorsAttnA", "DorsAttnB", "LimbicA", "LimbicB", "SalVentAttnA", "SalVentAttnB",
            "SomMotA", "SomMotB", "TempPar", "VisCent", "VisPeri"),
  yeo7  = c("global", "Vis", "SomMot", "DorsAttn", "SalVentAttn", "Limbic", "Cont", "Default")
)

seed_file <- function(seed, scale, parcellation = 200, networks = "yeo17") {
  if (scale == "parcel") {
    file.path(DATA_DIR, sprintf("%d_%s_ROI_correlations.csv", parcellation, seed))
  } else {
    suffix <- if (networks == "yeo7") "_yeo7" else ""
    file.path(DATA_DIR, sprintf("%d_%s_network_connectivity%s.csv", parcellation, seed, suffix))
  }
}

#' Seed connectivity in long form
read_fc_long <- function(path, scale, networks = "yeo17") {
  d <- read.csv(path, check.names = FALSE)
  if (scale == "parcel") {
    cols <- grep("^ROI_[0-9]+$", names(d), value = TRUE)
    d %>% select(sub, all_of(cols)) %>%
      pivot_longer(-sub, names_to = "region", values_to = "value")
  } else {
    cols <- grep("_mean$", names(d), value = TRUE)
    d %>% select(sub, all_of(cols)) %>%
      pivot_longer(-sub, names_to = "region", values_to = "value") %>%
      mutate(region = sub("_mean$", "", region),
             region = factor(region, levels = NETWORK_ORDER[[networks]]))
  }
}

#' Analysis data set: AB+ participants, all visits, all regions
#'
#' @param seed seed region whose connectivity is predictor ("LC", "shifted_LC" or "VTA")
#' @param scale "network" (global mean + networks) or "parcel"
#' @param control optional second seed whose connectivity to the SAME region is
#'        added as control_value (joint models)
prepare_trajectory_data <- function(seed, scale, parcellation = 200, networks = "yeo17",
                                    control = NULL) {
  fc <- read_fc_long(seed_file(seed, scale, parcellation, networks), scale, networks) %>%
    rename(fc = value)

  if (!is.null(control)) {
    ctrl <- read_fc_long(seed_file(control, scale, parcellation, networks), scale, networks) %>%
      rename(control_value = value)
    fc <- inner_join(fc, ctrl, by = c("sub", "region"))
  }

  read.csv(file.path(DATA_DIR, "longitudinal_data_processed.csv")) %>%
    inner_join(fc, by = "sub", relationship = "many-to-many") %>%
    filter(amyloid_status == "amyloid+") %>%
    mutate(
      QSDY_years = QSDY / 365.25,
      sex = factor(sex, levels = c("F", "M")),
      APOE_status = factor(APOE_status, levels = c("E4-", "E4+")),
      # four-level sex x e4 composite
      demo_group = factor(paste(sex, APOE_status, sep = "_"),
                          levels = c("F_E4-", "F_E4+", "M_E4-", "M_E4+"))
    )
}

#' Parcel labels (ROI_Number, ROI_Name, Network) 
read_roi_info <- function(parcellation = 200) {
  path <- file.path(DATA_DIR, sprintf("roi_network_assignments_yeo17_%d.csv", parcellation))
  if (file.exists(path)) read.csv(path)[, c("ROI_Number", "ROI_Name", "Network")] else NULL
}

add_roi_info <- function(results, roi_info) {
  results$ROI_Number <- as.numeric(sub("ROI_", "", results$region))
  if (is.null(roi_info)) return(results)
  left_join(results, roi_info, by = "ROI_Number")
}
