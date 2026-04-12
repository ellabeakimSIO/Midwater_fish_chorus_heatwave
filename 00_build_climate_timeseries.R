# -----------------------------------------------------------------------------
# 00_build_climate_timeseries.R
#
# Purpose:
# Build the monthly environmental covariate table used in Kim et al. (2026)
# for UF440 interannual and climate analyses. This script combines:
#   1) monthly climate indices (PDO, ONI, NPGO),
#   2) monthly BEUTI and BEUTI anomaly by site-specific latitude bin, and
#   3) monthly sea surface temperature anomaly (SSTA) and marine heatwave
#      metrics derived from site-level daily SST/MHW time series.
#
# Output:
# A single monthly table for North_B and South_P containing:
# Site, Date, Year, Month, PDO, ONI, phase, ONI_month_window, NPGO,
# BEUTI, BEUTI_anom, BEUTI_latbin, SST_anom_detrended, MHW_days_pct_detrended
# Note this output is available on the dryad repository.
#
# In the paper:
# This script supports the environmental/climate covariates described in
# Methods section 4 and the monthly environmental time-series visualizations
# shown in Figure 4C and Supplementary Figure 4C.
#
# Notes:
# - BEUTI anomaly is calculated here from the full available BEUTI record,
#   by subtracting the monthly climatological mean for each latitude bin.
# - SSTA and MHW metrics are summarized to monthly values from daily site-level
#   time series generated upstream.
# - This script does not generate figures; plotting should be handled in a
#   separate script.
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(lubridate)
  library(purrr)
})

# =========================
# User settings
# =========================

# Set directories relative to the repository root when possible.
data_dir <- "data"
output_dir <- "derived_data"

# Input files
idx_file <- file.path(data_dir, "climate_indices_PDO_ONI_NPGO.csv")
beuti_file <- file.path(data_dir, "BEUTI_monthly.csv")

# Site-level daily SST / MHW time series
mhw_files <- tibble::tribble(
  ~Site,      ~mhw_file,                                  ~BEUTI_latbin,
  "North_B",  file.path(data_dir, "CINMSB_MHW_timeseries_90th.csv"), "34N",
  "South_P",  file.path(data_dir, "SOCALP_MHW_timeseries_90th.csv"), "33N"
)

# Date window used in the manuscript
date_start <- as.Date("2008-01-01")
date_end   <- as.Date("2026-01-01")

# Output file
out_file <- file.path(output_dir, "Climate_Timeseries.csv")

# =========================
# Helper functions
# =========================

read_delim_auto <- function(fp) {
  x <- tryCatch(
    read_csv(fp, show_col_types = FALSE),
    error = function(e) read_delim(fp, delim = "\t", show_col_types = FALSE)
  )
  
  # Handle cases where a tab-delimited file is read as one column
  if (ncol(x) == 1) {
    x <- read_delim(fp, delim = "\t", show_col_types = FALSE)
  }
  
  x
}

parse_date_flexibly <- function(x) {
  x_chr <- as.character(x)
  
  out <- suppressWarnings(as.Date(x_chr))
  bad <- is.na(out) & !is.na(x_chr)
  
  if (any(bad)) {
    out[bad] <- suppressWarnings(as.Date(parse_date_time(
      x_chr[bad],
      orders = c("ymd", "Y-m-d", "Y/m/d", "mdy", "m/d/Y"),
      tz = "UTC"
    )))
  }
  
  out
}

make_monthly_mhw <- function(fp, site_name) {
  df <- read_delim_auto(fp)
  
  required_cols <- c("date_utc", "is_mhw_detrended", "time_series_detrended")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns in MHW file for ", site_name, ": ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  df %>%
    mutate(
      date_utc = case_when(
        inherits(date_utc, "Date") ~ date_utc,
        inherits(date_utc, "POSIXct") ~ as.Date(date_utc),
        TRUE ~ parse_date_flexibly(date_utc)
      ),
      is_mhw_detrended = as.logical(is_mhw_detrended),
      time_series_detrended = as.numeric(time_series_detrended),
      Year = year(date_utc),
      Month = month(date_utc)
    ) %>%
    filter(!is.na(date_utc)) %>%
    group_by(Year, Month) %>%
    summarise(
      Date = make_date(first(Year), first(Month), 1),
      SST_anom_detrended = if (all(is.na(time_series_detrended))) {
        NA_real_
      } else {
        mean(time_series_detrended, na.rm = TRUE)
      },
      MHW_days_pct_detrended = if (all(is.na(is_mhw_detrended))) {
        NA_real_
      } else {
        mean(is_mhw_detrended, na.rm = TRUE) * 100
      },
      .groups = "drop"
    ) %>%
    mutate(Site = site_name) %>%
    select(Date, Year, Month, Site, SST_anom_detrended, MHW_days_pct_detrended)
}

# =========================
# Read monthly climate indices
# =========================

idx_monthly <- read_delim_auto(idx_file) %>%
  mutate(
    Date = parse_date_flexibly(Date),
    Year = as.integer(Year),
    Month = as.integer(Month)
  )

# Keep only the expected index columns if present
idx_keep <- c("Date", "Year", "Month", "PDO", "ONI", "phase", "ONI_month_window", "NPGO")
idx_monthly <- idx_monthly %>%
  select(any_of(idx_keep))

# =========================
# Read and process BEUTI
# =========================

beuti_raw <- read_delim_auto(beuti_file)

# Assumes first two columns are year and month, as in the working script
beuti_full <- beuti_raw %>%
  rename(year = 1, month = 2) %>%
  mutate(
    year = as.integer(year),
    month = as.integer(month),
    Date = as.Date(sprintf("%d-%02d-01", year, month))
  )

lat_cols <- setdiff(names(beuti_full), c("year", "month", "Date"))

beuti_long <- beuti_full %>%
  pivot_longer(
    cols = all_of(lat_cols),
    names_to = "latbin",
    values_to = "BEUTI"
  ) %>%
  mutate(BEUTI = as.numeric(BEUTI))

# Monthly climatology across full available BEUTI record
beuti_clim <- beuti_long %>%
  group_by(latbin, month) %>%
  summarise(
    BEUTI_clim = if (all(is.na(BEUTI))) NA_real_ else mean(BEUTI, na.rm = TRUE),
    .groups = "drop"
  )

beuti_long <- beuti_long %>%
  left_join(beuti_clim, by = c("latbin", "month")) %>%
  mutate(BEUTI_anom = BEUTI - BEUTI_clim)

beuti_long_window <- beuti_long %>%
  filter(Date >= date_start, Date <= date_end)

beuti_sites <- mhw_files %>%
  select(Site, BEUTI_latbin) %>%
  left_join(
    beuti_long_window,
    by = c("BEUTI_latbin" = "latbin")
  ) %>%
  transmute(
    Site,
    Date,
    Year = year,
    Month = month,
    BEUTI,
    BEUTI_anom,
    BEUTI_latbin
  )

# =========================
# Read and summarize site-level SSTA / MHW files
# =========================

mhw_monthly_sites <- purrr::map2_dfr(
  .x = mhw_files$mhw_file,
  .y = mhw_files$Site,
  .f = make_monthly_mhw
)

# =========================
# Build full monthly backbone and merge
# =========================

month_seq <- seq(date_start, date_end, by = "1 month")

backbone <- tidyr::crossing(
  Site = mhw_files$Site,
  Date = month_seq
) %>%
  mutate(
    Year = year(Date),
    Month = month(Date)
  )

climate_timeseries <- backbone %>%
  left_join(idx_monthly, by = c("Date", "Year", "Month")) %>%
  left_join(beuti_sites, by = c("Site", "Date", "Year", "Month")) %>%
  left_join(mhw_monthly_sites, by = c("Site", "Date", "Year", "Month")) %>%
  arrange(Site, Date)

# =========================
# Quality checks
# =========================

message("Date range in output: ", min(climate_timeseries$Date), " to ", max(climate_timeseries$Date))
message("Sites in output: ", paste(unique(climate_timeseries$Site), collapse = ", "))
print(dplyr::glimpse(climate_timeseries))

# =========================
# Save output
# =========================

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write_csv(climate_timeseries, out_file)

message("Wrote monthly climate table to: ", out_file)