# -----------------------------------------------------------------------------
# 08_fit_plot_figure6c_southT_gamm.R
#
# Purpose:
# Fit the South_T GAMM relating nightly maximum sunset chorus SPL to
# backscatter strength, day of year, and lunar phase, select the best
# backscatter depth bin by AIC, and generate the predicted smooths used for
# Figure 6C.
#
# In the paper:
# This script reproduces the GAMM described in Supplemental Table 5 and the
# predicted smooths shown in Figure 6C of Kim et al. (2026). The final model
# evaluates nightly maximum sunset chorus SPL as a function of mean volume
# backscattering strength, day of year, and lunar phase, with an AR(1)
# correlation structure.
#
# Inputs:
#   - nightlyDepth_metrics.csv
#   - South_T_master_dryad.csv
#
# Outputs:
#   - Figure6C_gamm_smooths_South_T.pdf
#   - Figure6C_gamm_smooths_South_T.png
#   - Figure6C_gamm_summary_South_T.txt
#   - Figure6C_candidate_depth_models.csv
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(lubridate)
  library(ggplot2)
  library(mgcv)
  library(nlme)
  library(tibble)
  library(patchwork)
  library(broom)
})

# =========================
# User settings
# =========================

data_dir <- "data"
output_dir <- "figures"

depth_path <- file.path(data_dir, "nightlyDepth_metrics.csv")
spl_path <- file.path(data_dir, "South_T_master_dryad.csv")

metric_family <- "SvMean"
tz_local <- "America/Los_Angeles"
zmin <- 100
zmax <- 200

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# =========================
# Check inputs
# =========================

if (!file.exists(depth_path)) stop("Depth file not found: ", depth_path)
if (!file.exists(spl_path)) stop("SPL file not found: ", spl_path)

# =========================
# Helper functions
# =========================

moon_phase_fraction <- function(night_date, tz_local = "America/Los_Angeles") {
  if (is.na(night_date)) return(NA_real_)
  
  epoch <- as.POSIXct("2000-01-06 18:14:00", tz = "UTC")
  local_noon <- as.POSIXct(paste0(night_date, " 12:00:00"), tz = tz_local)
  t_utc <- with_tz(local_noon, tzone = "UTC")
  synodic <- 29.53058867
  days_since <- as.numeric(difftime(t_utc, epoch, units = "days"))
  
  as.numeric((days_since %% synodic) / synodic)
}

safe_num <- function(x, default = NA_real_) {
  if (is.null(x) || length(x) == 0) return(default)
  x <- suppressWarnings(as.numeric(x[1]))
  if (!is.finite(x)) return(default)
  x
}

parse_timestamp_utc <- function(x) {
  x_chr <- as.character(x)
  
  out <- suppressWarnings(ymd_hms(x_chr, tz = "UTC", quiet = TRUE))
  bad <- is.na(out) & !is.na(x_chr)
  
  if (any(bad)) {
    out[bad] <- suppressWarnings(parse_date_time(
      x_chr[bad],
      orders = c("Y-m-d H:M:S", "Y/m/d H:M:S", "ymd HMS", "ymd HM"),
      tz = "UTC"
    ))
  }
  
  out
}

# =========================
# Read inputs
# =========================

message("Reading input files...")
depth <- read_csv(depth_path, show_col_types = FALSE)
spl <- read_csv(spl_path, show_col_types = FALSE)

if (!("Date" %in% names(depth))) {
  stop("nightlyDepth_metrics.csv must contain a Date column.")
}
if (!("timestamp_utc" %in% names(spl))) {
  stop("South_T_master_dryad.csv must contain timestamp_utc.")
}
if (!("UF440_sunset_SPL" %in% names(spl))) {
  stop("South_T_master_dryad.csv must contain UF440_sunset_SPL.")
}

depth2 <- depth %>%
  mutate(
    Date = dmy(Date),
    NightDate = Date
  )

spl2 <- spl %>%
  mutate(
    timestamp_utc = parse_timestamp_utc(timestamp_utc),
    tLocal = with_tz(timestamp_utc, tzone = tz_local),
    NightDate = as_date(tLocal) + days(1),
    UF440_sunset_SPL = as.numeric(UF440_sunset_SPL)
  )

# Nightly maximum sunset SPL, matching manuscript methods for modeling
spl_night <- spl2 %>%
  group_by(NightDate) %>%
  summarise(
    UF440_sunset_SPL = if (all(is.na(UF440_sunset_SPL))) NA_real_ else max(UF440_sunset_SPL, na.rm = TRUE),
    .groups = "drop"
  )

merged <- depth2 %>%
  select(-Date, -any_of("Sunset")) %>%
  left_join(spl_night, by = "NightDate") %>%
  arrange(NightDate) %>%
  mutate(
    doy = yday(NightDate),
    moon_frac = vapply(NightDate, moon_phase_fraction, numeric(1), tz_local = tz_local)
  )

# =========================
# Candidate depth-bin selection
# =========================

fam_pattern <- paste0("^", metric_family, "_[0-9]+m$")
cand_cols <- names(merged)[grepl(fam_pattern, names(merged))]
cand_depths <- as.numeric(sub(".*_([0-9]+)m$", "\\1", cand_cols))

keep <- is.finite(cand_depths) & cand_depths >= zmin & cand_depths <= zmax
cand_cols <- cand_cols[keep]
cand_depths <- cand_depths[keep]

if (length(cand_cols) == 0) {
  stop("No candidate depth-bin columns found for metric_family = ", metric_family)
}

fit_one <- function(pred_col, depth_m) {
  dmod <- merged %>%
    select(NightDate, UF440_sunset_SPL, doy, moon_frac, all_of(pred_col)) %>%
    filter(
      is.finite(UF440_sunset_SPL),
      is.finite(.data[[pred_col]]),
      is.finite(doy),
      is.finite(moon_frac)
    ) %>%
    arrange(NightDate) %>%
    mutate(obs_index = seq_len(n()))
  
  n <- nrow(dmod)
  if (n < 30) return(NULL)
  
  pred_sd <- sd(dmod[[pred_col]], na.rm = TRUE)
  if (!is.finite(pred_sd) || pred_sd == 0) return(NULL)
  
  form <- as.formula(
    paste0(
      "UF440_sunset_SPL ~ ",
      "s(", pred_col, ", bs = 'ts', k = 4) + ",
      "s(doy, bs = 'cc', k = 4) + ",
      "s(moon_frac, bs = 'cc', k = 4)"
    )
  )
  
  fit <- tryCatch(
    mgcv::gamm(
      formula = form,
      data = dmod,
      correlation = nlme::corAR1(form = ~ obs_index),
      method = "ML",
      knots = list(doy = c(0.5, 365.5), moon_frac = c(0, 1))
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit) || is.null(fit$gam) || is.null(fit$lme)) return(NULL)
  
  aic_gam <- safe_num(tryCatch(AIC(fit$gam), error = function(e) numeric(0)))
  aic_lme <- safe_num(tryCatch(AIC(fit$lme), error = function(e) numeric(0)))
  score_AIC <- if (is.finite(aic_gam)) aic_gam else aic_lme
  if (!is.finite(score_AIC)) return(NULL)
  
  tibble(
    predictor = pred_col,
    depth_m = depth_m,
    score_AIC = score_AIC,
    n = n
  )
}

message("Selecting best depth bin...")
fit_list <- lapply(seq_along(cand_cols), function(i) fit_one(cand_cols[i], cand_depths[i]))
fit_list <- fit_list[!vapply(fit_list, is.null, logical(1))]

if (length(fit_list) == 0) {
  stop("All candidate GAMM fits failed or were skipped.")
}

comp_tbl <- bind_rows(fit_list) %>%
  filter(is.finite(score_AIC)) %>%
  arrange(score_AIC)

if (nrow(comp_tbl) == 0) {
  stop("All candidate fits had non-finite AIC.")
}

best_pred <- comp_tbl$predictor[1]

# =========================
# Refit best model with REML
# =========================

best_dat <- merged %>%
  select(NightDate, UF440_sunset_SPL, doy, moon_frac, all_of(best_pred)) %>%
  filter(
    is.finite(UF440_sunset_SPL),
    is.finite(.data[[best_pred]]),
    is.finite(doy),
    is.finite(moon_frac)
  ) %>%
  arrange(NightDate) %>%
  mutate(obs_index = seq_len(n()))

best_form <- as.formula(
  paste0(
    "UF440_sunset_SPL ~ ",
    "s(", best_pred, ", bs = 'ts', k = 4) + ",
    "s(doy, bs = 'cc', k = 4) + ",
    "s(moon_frac, bs = 'cc', k = 4)"
  )
)

message("Refitting best model with REML...")
best_fit <- mgcv::gamm(
  formula = best_form,
  data = best_dat,
  correlation = nlme::corAR1(form = ~ obs_index),
  method = "REML",
  knots = list(doy = c(0.5, 365.5), moon_frac = c(0, 1))
)

# =========================
# Predicted smooths
# =========================

theme_chorus <- theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 12, face = "bold"),
    plot.margin = margin(5, 5, 5, 5),
    legend.position = "none"
  )

n_pred <- 200

# A) Seasonal effect
pred_doy <- seq(min(best_dat$doy, na.rm = TRUE), max(best_dat$doy, na.rm = TRUE), length.out = n_pred)

newdat_doy <- data.frame(
  doy = pred_doy,
  moon_frac = rep(median(best_dat$moon_frac, na.rm = TRUE), n_pred)
)
newdat_doy[[best_pred]] <- rep(median(best_dat[[best_pred]], na.rm = TRUE), n_pred)

pred_doy_fit <- predict(best_fit$gam, newdata = newdat_doy, se.fit = TRUE, type = "response")

plot_doy <- tibble(
  x = pred_doy,
  fit = as.numeric(pred_doy_fit$fit),
  lower = as.numeric(pred_doy_fit$fit - 2 * pred_doy_fit$se.fit),
  upper = as.numeric(pred_doy_fit$fit + 2 * pred_doy_fit$se.fit)
)

rug_doy <- tibble(x = best_dat$doy)

p_doy <- ggplot(plot_doy, aes(x = x, y = fit)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey70", alpha = 0.3) +
  geom_line(linewidth = 0.9) +
  geom_rug(data = rug_doy, aes(x = x), inherit.aes = FALSE, sides = "b", alpha = 0.5) +
  theme_chorus +
  labs(
    title = "A. Seasonal effect",
    x = "Day of year",
    y = "Predicted SPL (dB)"
  )

# B) Backscatter effect
pred_backscatter <- seq(
  min(best_dat[[best_pred]], na.rm = TRUE),
  max(best_dat[[best_pred]], na.rm = TRUE),
  length.out = n_pred
)

newdat_backscatter <- data.frame(
  doy = rep(median(best_dat$doy, na.rm = TRUE), n_pred),
  moon_frac = rep(median(best_dat$moon_frac, na.rm = TRUE), n_pred)
)
newdat_backscatter[[best_pred]] <- pred_backscatter

pred_bs <- predict(best_fit$gam, newdata = newdat_backscatter, se.fit = TRUE, type = "response")

plot_bs <- tibble(
  x = pred_backscatter,
  fit = as.numeric(pred_bs$fit),
  lower = as.numeric(pred_bs$fit - 2 * pred_bs$se.fit),
  upper = as.numeric(pred_bs$fit + 2 * pred_bs$se.fit)
)

rug_bs <- tibble(x = best_dat[[best_pred]])

p_bs <- ggplot(plot_bs, aes(x = x, y = fit)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey70", alpha = 0.3) +
  geom_line(linewidth = 0.9) +
  geom_rug(data = rug_bs, aes(x = x), inherit.aes = FALSE, sides = "b", alpha = 0.5) +
  theme_chorus +
  labs(
    title = "B. Backscatter effect",
    x = best_pred,
    y = "Predicted SPL (dB)"
  )

# C) Lunar effect
pred_moon <- seq(0, 1, length.out = n_pred)

newdat_moon <- data.frame(
  doy = rep(median(best_dat$doy, na.rm = TRUE), n_pred),
  moon_frac = pred_moon
)
newdat_moon[[best_pred]] <- rep(median(best_dat[[best_pred]], na.rm = TRUE), n_pred)

pred_moon_fit <- predict(best_fit$gam, newdata = newdat_moon, se.fit = TRUE, type = "response")

plot_moon <- tibble(
  x = pred_moon,
  fit = as.numeric(pred_moon_fit$fit),
  lower = as.numeric(pred_moon_fit$fit - 2 * pred_moon_fit$se.fit),
  upper = as.numeric(pred_moon_fit$fit + 2 * pred_moon_fit$se.fit)
)

rug_moon <- tibble(x = best_dat$moon_frac)

p_moon <- ggplot(plot_moon, aes(x = x, y = fit)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey70", alpha = 0.3) +
  geom_line(linewidth = 0.9) +
  geom_rug(data = rug_moon, aes(x = x), inherit.aes = FALSE, sides = "b", alpha = 0.5) +
  scale_x_continuous(
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    labels = c("New", "1st Q", "Full", "3rd Q", "New")
  ) +
  theme_chorus +
  labs(
    title = "C. Lunar effect",
    x = "Moon phase",
    y = "Predicted SPL (dB)"
  )

figure_6c <- (p_doy | p_bs | p_moon) +
  plot_annotation(
    title = paste0("South_T GAMM predicted smooths (best depth bin: ", best_pred, ")")
  )

# =========================
# Save outputs
# =========================

write_csv(comp_tbl, file.path(output_dir, "Figure6C_candidate_depth_models.csv"))

ggsave(
  file.path(output_dir, "Figure6C_gamm_smooths_South_T.pdf"),
  figure_6c,
  width = 14,
  height = 4.8
)

ggsave(
  file.path(output_dir, "Figure6C_gamm_smooths_South_T.png"),
  figure_6c,
  width = 14,
  height = 4.8,
  dpi = 600,
  bg = "white"
)

sink(file.path(output_dir, "Figure6C_gamm_summary_South_T.txt"))
cat("=== Requested metric family ===\n")
print(metric_family)
cat("\n=== Candidate depth range ===\n")
print(c(zmin, zmax))
cat("\n=== Best predictor ===\n")
print(best_pred)
cat("\n=== Candidate models ranked by AIC ===\n")
print(comp_tbl)
cat("\n=== GAM summary (best model) ===\n")
print(summary(best_fit$gam))
cat("\n=== LME summary (best model) ===\n")
print(summary(best_fit$lme))
sink()

message("Done. Outputs written to: ", normalizePath(output_dir, winslash = "/"))