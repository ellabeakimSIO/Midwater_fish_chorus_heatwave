# -----------------------------------------------------------------------------
# 04_fit_plot_figure4b_northB_seasonal_presence_gam
#
# Purpose:
# Fit the seasonal chorus presence GAM used for Figure 4B in Kim et al. (2026),
# using North_B chorus presence and seasonal mean sea surface temperature
# anomaly (SSTA).
#
# In the paper:
# This script reproduces the North_B chorus presence GAM described in Methods 6.1.
# Chorus presence is modeled at the seasonal scale as a function of:
#   1) season (cyclic smooth), and
#   2) seasonal mean SSTA (thin-plate shrinkage smooth).
#
# Response definition:
# - Days are classified as present when at least one chorus detection occurs in
#   any on-effort 20 min bin.
# - Seasons are Jan-Mar, Apr-Jun, Jul-Sep, and Oct-Dec.
# - Seasonal chorus presence is scored as present when chorus is detected on at
#   least 10% of recorded days in that season.
# - Seasons with no recorded days are excluded.
#
# Model:
# presence_bin ~ s(Season_num, bs = "cc", k = 4) + s(SSTA_season, bs = "ts", k = 4)
# family = binomial(link = "logit"), weights = recorded_days, method = "REML"
#
# Inputs:
#   - data/North_B_master_dryad.csv
#   - data/Climate_Timeseries.csv
#
# Outputs:
#   - figures/Figure4B_North_B_season_GAM.png
#   - figures/Figure4B_North_B_season_GAM.pdf
#   - figures/Figure4B_North_B_season_GAM_summary.txt
#   - figures/Figure4B_North_B_season_GAM_data.csv
#
# Notes:
# - This is the manuscript-matched seasonal model, not the earlier monthly
#   trial-based GAM variants.
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(lubridate)
  library(mgcv)
  library(ggplot2)
  library(tibble)
  library(patchwork)
  library(readr)
})

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

# ---------------------------
# Settings
# ---------------------------
data_dir <- "data"
output_dir <- "figures"

acoustic_file <- file.path(data_dir, "North_B_master_dryad.csv")
climate_file  <- file.path(data_dir, "Climate_Timeseries.csv")

chorus_prop_thresh_season <- 0.10

season_k <- 4
season_knots <- list(Season_num = c(0.5, 4.5))
season_labels <- c("Jan–Mar", "Apr–Jun", "Jul–Sep", "Oct–Dec")

sst_k <- 4

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

theme_clean <- theme_bw(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

# ---------------------------
# Helpers
# ---------------------------
read_delim_auto <- function(fp) {
  x <- tryCatch(
    read_csv(fp, show_col_types = FALSE),
    error = function(e) read_delim(fp, delim = "\t", show_col_types = FALSE)
  )
  if (ncol(x) == 1) {
    x <- read_delim(fp, delim = "\t", show_col_types = FALSE)
  }
  x
}

parse_datetime_flexibly <- function(x, tz = "UTC") {
  x_chr <- as.character(x)
  
  out <- suppressWarnings(ymd_hms(x_chr, tz = tz, quiet = TRUE))
  bad <- is.na(out) & !is.na(x_chr)
  
  if (any(bad)) {
    out[bad] <- suppressWarnings(parse_date_time(
      x_chr[bad],
      orders = c("Y-m-d H:M:S", "Y/m/d H:M:S", "ymd HMS", "ymd HM"),
      tz = tz
    ))
  }
  
  out
}

# ---------------------------
# Read data
# ---------------------------
all_raw <- read_delim_auto(acoustic_file) %>%
  mutate(
    Site = "North_B",
    timestamp_dt = parse_datetime_flexibly(timestamp_dt, tz = "UTC"),
    UF440_roving_PSD = as.numeric(UF440_roving_PSD),
    UF440_PA_no20iso_noday = as.numeric(UF440_PA_no20iso_noday)
  ) %>%
  filter(!is.na(timestamp_dt))

idx_beuti_mhw <- read_delim_auto(climate_file)

stopifnot(all(c("Site", "Date", "SST_anom_detrended") %in% names(idx_beuti_mhw)))

# ---------------------------
# 1) Daily chorus effort + daily chorus presence
# ---------------------------
daily_chorus <- all_raw %>%
  mutate(
    Date = as.Date(with_tz(timestamp_dt, "UTC")),
    on_effort_bin = !is.na(UF440_roving_PSD),
    det_bin = on_effort_bin & (coalesce(UF440_PA_no20iso_noday, 0) > 0)
  ) %>%
  group_by(Site, Date) %>%
  summarise(
    recorded_day = any(on_effort_bin, na.rm = TRUE),
    present_day  = recorded_day & any(det_bin, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Year = year(Date),
    Season_num = quarter(Date),
    Season = factor(season_labels[Season_num], levels = season_labels)
  )

# Seasonal response (>=10% of recorded days present)
dat_season_core <- daily_chorus %>%
  group_by(Site, Year, Season_num, Season) %>%
  summarise(
    recorded_days = sum(recorded_day, na.rm = TRUE),
    present_days  = sum(present_day, na.rm = TRUE),
    chorus_prop   = if_else(recorded_days > 0, present_days / recorded_days, NA_real_),
    presence_bin  = as.integer(!is.na(chorus_prop) & chorus_prop >= chorus_prop_thresh_season),
    .groups = "drop"
  ) %>%
  filter(recorded_days > 0) %>%
  arrange(Site, Year, Season_num)

# ---------------------------
# 2) Seasonal mean SSTA from monthly Climate_Timeseries.csv
# ---------------------------
ssta_season_all <- idx_beuti_mhw %>%
  mutate(
    Date = as.Date(Date),
    Year = year(Date),
    Season_num = quarter(Date)
  ) %>%
  filter(Site %in% c("North_B", "CINMS_B")) %>%
  transmute(
    Site = "North_B",
    Year,
    Season_num,
    SSTA = as.numeric(SST_anom_detrended)
  ) %>%
  filter(!is.na(SSTA)) %>%
  group_by(Site, Year, Season_num) %>%
  summarise(
    SSTA_season = mean(SSTA, na.rm = TRUE),
    n_months_ssta = n(),
    .groups = "drop"
  )

# Join seasonal SSTA to chorus dataset
dat_season <- dat_season_core %>%
  left_join(ssta_season_all, by = c("Site", "Year", "Season_num")) %>%
  filter(!is.na(SSTA_season)) %>%
  arrange(Site, Year, Season_num)

write_csv(dat_season, file.path(output_dir, "Figure4B_North_B_season_GAM_data.csv"))

# ---------------------------
# 3) Fit GAM
# ---------------------------
fit_season_ssta_gam <- function(site_name) {
  d <- dat_season %>%
    filter(Site == site_name) %>%
    arrange(Year, Season_num)
  
  m <- mgcv::gam(
    presence_bin ~
      s(Season_num, bs = "cc", k = season_k) +
      s(SSTA_season, bs = "ts", k = sst_k),
    family = binomial(link = "logit"),
    data = d,
    weights = recorded_days,
    method = "REML",
    knots = season_knots,
    na.action = na.exclude
  )
  
  list(model = m, data = d)
}

res_ssta_north <- fit_season_ssta_gam("North_B")

# ---------------------------
# 4) Plot smooths
# ---------------------------
plot_season_effect <- function(res_obj, n = 200, title = NULL) {
  m <- res_obj$model
  d <- res_obj$data
  hold_ssta <- median(d$SSTA_season, na.rm = TRUE)
  
  xseq <- seq(1, 4, length.out = n)
  newd <- tibble(Season_num = xseq, SSTA_season = hold_ssta)
  
  pr <- predict(m, newdata = newd, type = "link", se.fit = TRUE)
  dfp <- newd %>%
    mutate(
      p  = plogis(pr$fit),
      lo = plogis(pr$fit - 1.96 * pr$se.fit),
      hi = plogis(pr$fit + 1.96 * pr$se.fit)
    )
  
  rug_df <- d %>%
    transmute(Season_num = Season_num + runif(n(), -0.06, 0.06))
  
  ggplot(dfp, aes(x = Season_num, y = p)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.18) +
    geom_line(linewidth = 1.1) +
    geom_rug(
      data = rug_df,
      aes(x = Season_num),
      inherit.aes = FALSE,
      sides = "b",
      alpha = 0.25
    ) +
    scale_x_continuous(breaks = 1:4, labels = season_labels) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(
      x = NULL,
      y = "Predicted P(season present)",
      title = title %||% "Season smooth (SSTA held at median)"
    ) +
    theme_clean
}

plot_ssta_effect <- function(res_obj, n = 240, title = NULL) {
  m <- res_obj$model
  d <- res_obj$data
  hold_season <- 2.5
  
  rng <- range(d$SSTA_season, na.rm = TRUE)
  xseq <- seq(rng[1], rng[2], length.out = n)
  newd <- tibble(Season_num = hold_season, SSTA_season = xseq)
  
  pr <- predict(m, newdata = newd, type = "link", se.fit = TRUE)
  dfp <- newd %>%
    mutate(
      p  = plogis(pr$fit),
      lo = plogis(pr$fit - 1.96 * pr$se.fit),
      hi = plogis(pr$fit + 1.96 * pr$se.fit)
    )
  
  ggplot(dfp, aes(x = SSTA_season, y = p)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.18) +
    geom_line(linewidth = 1.1) +
    geom_rug(
      data = d,
      aes(x = SSTA_season),
      inherit.aes = FALSE,
      sides = "b",
      alpha = 0.25
    ) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(
      x = "Seasonal mean SST anomaly (SST_anom_detrended)",
      y = "Predicted P(season present)",
      title = title %||% "SSTA smooth (Season held)"
    ) +
    theme_clean
}

p_north <- plot_season_effect(res_ssta_north, title = "North_B — Season") +
  plot_ssta_effect(res_ssta_north, title = "North_B — SSTA") +
  plot_layout(widths = c(1.15, 1))

print(p_north)

# ---------------------------
# 5) Export + summary
# ---------------------------
ggsave(
  file.path(output_dir, "Figure4B_North_B_season_GAM.png"),
  plot = p_north,
  width = 10,
  height = 4.2,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(output_dir, "Figure4B_North_B_season_GAM.pdf"),
  plot = p_north,
  width = 10,
  height = 4.2
)

sink(file.path(output_dir, "Figure4B_North_B_season_GAM_summary.txt"))
cat("=== Model data summary ===\n")
print(summary(res_ssta_north$data))
cat("\n=== GAM summary ===\n")
print(summary(res_ssta_north$model))
sink()