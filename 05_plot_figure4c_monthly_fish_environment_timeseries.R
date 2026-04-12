# -----------------------------------------------------------------------------
# 05_plot_figure4c_monthly_fish_environment_timeseries.R
#
# Purpose:
# Read analysis-ready UF440 acoustic master tables and the monthly climate table,
# compute monthly chorus summaries, and generate stacked monthly fish +
# environmental time-series panels for focal sites.
#
# Inputs:
#   - North_B_master_dryad.csv
#   - South_P_master_dryad.csv
#   - Climate_Timeseries.csv
#
# Outputs:
#   - One stacked figure per site showing:
#       1) monthly UF440 chorus presence (mean weekly chorus hours per month)
#          with monthly sunset SPL context,
#       2) PDO,
#       3) NPGO,
#       4) ONI,
#       5) BEUTI anomaly,
#       6) SST anomaly with marine heatwave months highlighted
#
# In the paper:
# This script supports the monthly environmental time-series visualizations used
# for site-level climate context figures.
#
# Notes:
# - UF440 presence is assumed to be coded in 20 min bins using the
#   UF440_PA_no20iso_noday column.
# - Sunset SPL is assumed to be stored in UF440_sunset_SPL.
# - Climate covariates are assumed to be provided as monthly values in
#   Climate_Timeseries.csv.
# - This script is written to run directly from analysis-ready Dryad files.
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(lubridate)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(grid)
  library(stringr)
  library(purrr)
})

# =========================
# User settings
# =========================

data_dir <- "data"
output_dir <- "figures"

bin_mins <- 20
week_start_day <- 1  # 1 = Monday

site_files <- tibble::tribble(
  ~Site,      ~file,
  "North_B",  file.path(data_dir, "North_B_master_dryad.csv"),
  "South_P",  file.path(data_dir, "South_P_master_dryad.csv")
)

climate_file <- file.path(data_dir, "Climate_Timeseries.csv")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# =========================
# Helper functions
# =========================

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
      orders = c(
        "Y-m-d H:M:S",
        "Y/m/d H:M:S",
        "ymd HMS",
        "ymd HM",
        "mdy HMS",
        "mdy HM"
      ),
      tz = tz
    ))
  }
  
  out
}

year_breaks <- function(start_date, end_date) {
  seq(floor_date(start_date, "year"), floor_date(end_date, "year"), by = "1 year")
}

quarter_breaks <- function(start_date, end_date) {
  seq(floor_date(start_date, "quarter"), floor_date(end_date, "quarter"), by = "3 months")
}

scale_x_months <- function(global_start, global_end, qbrks, ybrks) {
  scale_x_date(
    breaks = ybrks,
    labels = label_date("%Y"),
    minor_breaks = qbrks,
    limits = c(global_start, global_end %m+% months(1)),
    expand = expansion(mult = c(0, 0))
  )
}

env_range <- function(x) {
  r <- suppressWarnings(range(as.numeric(x), na.rm = TRUE))
  if (!all(is.finite(r)) || r[1] == r[2]) r <- c(-1, 1)
  pad <- 0.03 * (r[2] - r[1])
  c(r[1] - pad, r[2] + pad)
}

map_dot_to_env <- function(dot_hours, pres_top, ymin, ymax) {
  ymin + (dot_hours / pres_top) * (ymax - ymin)
}

secaxis_env <- function(pres_top, ymin, ymax) {
  sec_axis(~ (. - ymin) / (ymax - ymin) * pres_top, name = "Presence (hrs/wk)")
}

secaxis_env_blank <- function(pres_top, ymin, ymax) {
  sec_axis(~ (. - ymin) / (ymax - ymin) * pres_top, name = NULL)
}

hide_x_labels <- theme(
  axis.text.x = element_blank(),
  axis.title.x = element_blank()
)

theme_panel <- theme_bw(base_size = 11) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 9),
    axis.ticks.length = unit(3, "pt"),
    plot.margin = margin(t = 6, r = 10, b = 6, l = 10, unit = "pt"),
    legend.position = "none"
  )

# =========================
# Read acoustic master tables
# =========================

read_master_file <- function(fp, site_name, bin_mins = 20) {
  df <- read_delim_auto(fp)
  
  required_cols <- c(
    "Site",
    "timestamp_utc",
    "timestamp_dt",
    "UF440_PA_no20iso_noday",
    "UF440_sunset_SPL"
  )
  
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns in ", basename(fp), ": ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  df %>%
    mutate(
      Site = site_name,
      timestamp_utc = parse_datetime_flexibly(timestamp_utc, tz = "UTC"),
      timestamp_dt = parse_datetime_flexibly(timestamp_dt, tz = "UTC"),
      timestamp = coalesce(timestamp_dt, timestamp_utc),
      UF440_PA_no20iso_noday = as.numeric(UF440_PA_no20iso_noday),
      UF440_sunset_SPL = as.numeric(UF440_sunset_SPL)
    ) %>%
    filter(!is.na(timestamp)) %>%
    mutate(
      date = as.Date(timestamp),
      month = floor_date(date, "month"),
      week = floor_date(date, unit = "week", week_start = week_start_day),
      bin_hours = bin_mins / 60
    ) %>%
    select(Site, timestamp, date, month, week,
           UF440_PA_no20iso_noday, UF440_sunset_SPL, bin_hours)
}

acoustic_raw <- purrr::map2_dfr(site_files$file, site_files$Site, read_master_file)

if (nrow(acoustic_raw) == 0) {
  stop("No acoustic data were read. Check input file paths.")
}

global_start <- floor_date(min(acoustic_raw$month, na.rm = TRUE), "month")
global_end <- floor_date(max(acoustic_raw$month, na.rm = TRUE), "month")
month_seq <- seq(global_start, global_end, by = "1 month")
qbrks <- quarter_breaks(global_start, global_end)
ybrks <- year_breaks(global_start, global_end)

# =========================
# Monthly fish summaries
# =========================

compute_monthly_fish_summary <- function(df_site, month_seq, bin_mins = 20) {
  expected_bins_per_day <- 24 * 60 / bin_mins
  
  weekly_summary <- df_site %>%
    group_by(Site, month, week) %>%
    summarise(
      weekly_hours = sum(UF440_PA_no20iso_noday == 1, na.rm = TRUE) * (bin_mins / 60),
      bins_with_effort = sum(!is.na(UF440_PA_no20iso_noday)),
      .groups = "drop"
    )
  
  monthly_effort <- df_site %>%
    group_by(Site, month) %>%
    summarise(
      bins_with_effort = sum(!is.na(UF440_PA_no20iso_noday)),
      days_in_month = days_in_month(first(month)),
      mean_effort_frac = bins_with_effort / (days_in_month * expected_bins_per_day),
      month_sunset_mean = if (all(is.na(UF440_sunset_SPL))) NA_real_ else mean(UF440_sunset_SPL, na.rm = TRUE),
      .groups = "drop"
    )
  
  monthly_presence <- weekly_summary %>%
    group_by(Site, month) %>%
    summarise(
      avg_weekly_hours = if (all(is.na(weekly_hours))) NA_real_ else mean(weekly_hours, na.rm = TRUE),
      se_weekly_hours = if (sum(!is.na(weekly_hours)) <= 1) NA_real_ else sd(weekly_hours, na.rm = TRUE) / sqrt(sum(!is.na(weekly_hours))),
      n_weeks = sum(!is.na(weekly_hours)),
      .groups = "drop"
    )
  
  tibble(Site = unique(df_site$Site), month = month_seq) %>%
    left_join(monthly_presence, by = c("Site", "month")) %>%
    left_join(monthly_effort, by = c("Site", "month")) %>%
    mutate(
      avg_weekly_hours = if_else(is.na(mean_effort_frac) | mean_effort_frac == 0, NA_real_, avg_weekly_hours),
      dot_hours = if_else(!is.na(mean_effort_frac) & mean_effort_frac > 0, avg_weekly_hours, NA_real_)
    )
}

mon_list <- acoustic_raw %>%
  split(.$Site) %>%
  purrr::map(compute_monthly_fish_summary, month_seq = month_seq, bin_mins = bin_mins)

# Global SPL bins across sites for consistent color scale in fish panels
all_month_spl <- bind_rows(mon_list) %>% pull(month_sunset_mean)
spl_rng <- suppressWarnings(range(all_month_spl, na.rm = TRUE))
if (!all(is.finite(spl_rng)) || diff(spl_rng) == 0) {
  spl_rng <- c(0, 1)
}

spl_breaks <- seq(spl_rng[1], spl_rng[2], length.out = 6)
spl_breaks[6] <- spl_breaks[6] + 1e-6
spl_labels <- sprintf("%.1f", spl_breaks[1:5])
spl_levels <- c(spl_labels, "No data")

spl_cols6 <- c("blue", "cyan", "yellow", "orange", "red", "grey30")
names(spl_cols6) <- spl_levels

assign_spl_bin <- function(x, breaks, labels) {
  out <- cut(x, breaks = breaks, labels = labels, include.lowest = TRUE, right = FALSE)
  out <- as.character(out)
  out[is.na(out)] <- "No data"
  factor(out, levels = c(labels, "No data"))
}

# =========================
# Read climate data
# =========================

climate_timeseries <- read_delim_auto(climate_file) %>%
  mutate(
    Date = as.Date(Date),
    Site = as.character(Site),
    PDO = as.numeric(PDO),
    ONI = as.numeric(ONI),
    NPGO = as.numeric(NPGO),
    BEUTI_anom = as.numeric(BEUTI_anom),
    SST_anom_detrended = as.numeric(SST_anom_detrended),
    MHW_days_pct_detrended = as.numeric(MHW_days_pct_detrended),
    phase = as.character(phase)
  )

prep_idx_site <- function(site_name, month_seq, fish_summary) {
  climate_timeseries %>%
    filter(Site == site_name) %>%
    transmute(
      Site,
      month = Date,
      PDO,
      ONI,
      NPGO,
      phase,
      BEUTI_anom,
      SST_anom_detrended,
      MHW_days_pct_detrended
    ) %>%
    right_join(
      tibble(Site = site_name, month = month_seq),
      by = c("Site", "month")
    ) %>%
    left_join(
      fish_summary %>% select(Site, month, dot_hours, avg_weekly_hours, se_weekly_hours, mean_effort_frac, month_sunset_mean),
      by = c("Site", "month")
    ) %>%
    mutate(
      month_end = month %m+% months(1)
    ) %>%
    arrange(month)
}

# =========================
# Plotting functions
# =========================

col_pos <- "firebrick2"
col_neg <- "steelblue3"
col_neutral <- "chartreuse4"
col_bar_grey <- "grey30"
col_mhw <- "purple4"

dot_col <- "darkorchid1"
dot_size <- 1.2
dot_alpha <- 0.55

plot_fish_panel <- function(df, site_name, global_start, global_end, qbrks, ybrks,
                            spl_breaks, spl_labels, spl_cols6) {
  d <- df %>%
    mutate(
      spl_bin = assign_spl_bin(month_sunset_mean, spl_breaks, spl_labels),
      no_effort = is.na(mean_effort_frac) | mean_effort_frac == 0
    )
  
  ggplot(d, aes(x = month, y = avg_weekly_hours)) +
    geom_rect(
      data = d %>% filter(no_effort),
      aes(xmin = month, xmax = month_end, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE,
      fill = "grey90",
      color = NA
    ) +
    geom_col(aes(fill = spl_bin), width = 25, color = NA, na.rm = TRUE) +
    geom_errorbar(
      aes(ymin = avg_weekly_hours - se_weekly_hours, ymax = avg_weekly_hours + se_weekly_hours),
      width = 10,
      linewidth = 0.3,
      na.rm = TRUE
    ) +
    scale_fill_manual(values = spl_cols6, drop = FALSE) +
    scale_x_months(global_start, global_end, qbrks, ybrks) +
    labs(y = "Presence\n(hrs/wk)", x = NULL, title = site_name) +
    theme_panel +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 11)
    )
}

plot_pdo <- function(df, global_start, global_end, qbrks, ybrks, pres_top, show_x = FALSE) {
  yr <- env_range(df$PDO); ymin <- yr[1]; ymax <- yr[2]
  d <- df %>%
    mutate(
      pdo_cat = case_when(is.na(PDO) ~ NA_character_, PDO < 0 ~ "neg", TRUE ~ "pos"),
      dot_y = map_dot_to_env(dot_hours, pres_top, ymin, ymax)
    )
  
  p <- ggplot(d, aes(x = month, y = PDO)) +
    geom_hline(yintercept = 0, linewidth = 0.4, color = "grey40") +
    geom_col(aes(fill = pdo_cat), width = 25, color = NA, na.rm = TRUE) +
    scale_fill_manual(values = c(neg = col_neg, pos = col_pos), na.translate = FALSE) +
    geom_point(aes(y = dot_y), color = dot_col, alpha = dot_alpha, size = dot_size, na.rm = TRUE) +
    scale_x_months(global_start, global_end, qbrks, ybrks) +
    scale_y_continuous(limits = c(ymin, ymax), sec.axis = secaxis_env(pres_top, ymin, ymax)) +
    labs(y = "PDO", x = NULL) +
    theme_panel
  
  if (!show_x) p <- p + hide_x_labels
  p
}

plot_npgo <- function(df, global_start, global_end, qbrks, ybrks, pres_top, show_x = FALSE) {
  yr <- env_range(df$NPGO); ymin <- yr[1]; ymax <- yr[2]
  d <- df %>%
    mutate(dot_y = map_dot_to_env(dot_hours, pres_top, ymin, ymax))
  
  p <- ggplot(d, aes(x = month, y = NPGO)) +
    geom_hline(yintercept = 0, linewidth = 0.4, color = "grey40") +
    geom_col(width = 25, fill = col_bar_grey, color = NA, na.rm = TRUE) +
    geom_point(aes(y = dot_y), color = dot_col, alpha = dot_alpha, size = dot_size, na.rm = TRUE) +
    scale_x_months(global_start, global_end, qbrks, ybrks) +
    scale_y_continuous(limits = c(ymin, ymax), sec.axis = secaxis_env_blank(pres_top, ymin, ymax)) +
    labs(y = "NPGO", x = NULL) +
    theme_panel
  
  if (!show_x) p <- p + hide_x_labels
  p
}

plot_oni <- function(df, global_start, global_end, qbrks, ybrks, pres_top, show_x = FALSE) {
  yr <- env_range(df$ONI); ymin <- yr[1]; ymax <- yr[2]
  d <- df %>%
    mutate(
      oni_cat = case_when(
        is.na(ONI) ~ NA_character_,
        str_detect(phase %||% "", regex("Cool", ignore_case = TRUE)) ~ "Cool Phase/La Nina",
        str_detect(phase %||% "", regex("Warm", ignore_case = TRUE)) ~ "Warm Phase/El Nino",
        TRUE ~ "Neutral Phase"
      ),
      dot_y = map_dot_to_env(dot_hours, pres_top, ymin, ymax)
    )
  
  p <- ggplot(d, aes(x = month, y = ONI)) +
    geom_hline(yintercept = 0, linewidth = 0.4, color = "grey40") +
    geom_col(aes(fill = oni_cat), width = 25, color = NA, na.rm = TRUE) +
    scale_fill_manual(
      values = c(
        "Cool Phase/La Nina" = col_neg,
        "Neutral Phase" = col_neutral,
        "Warm Phase/El Nino" = col_pos
      ),
      na.translate = FALSE
    ) +
    geom_point(aes(y = dot_y), color = dot_col, alpha = dot_alpha, size = dot_size, na.rm = TRUE) +
    scale_x_months(global_start, global_end, qbrks, ybrks) +
    scale_y_continuous(limits = c(ymin, ymax), sec.axis = secaxis_env_blank(pres_top, ymin, ymax)) +
    labs(y = "ONI", x = NULL) +
    theme_panel
  
  if (!show_x) p <- p + hide_x_labels
  p
}

plot_beuti <- function(df, global_start, global_end, qbrks, ybrks, pres_top, show_x = FALSE) {
  yr <- env_range(df$BEUTI_anom); ymin <- yr[1]; ymax <- yr[2]
  d <- df %>%
    mutate(
      beuti_cat = case_when(is.na(BEUTI_anom) ~ NA_character_, BEUTI_anom < 0 ~ "neg", TRUE ~ "pos"),
      dot_y = map_dot_to_env(dot_hours, pres_top, ymin, ymax)
    )
  
  p <- ggplot(d, aes(x = month, y = BEUTI_anom)) +
    geom_hline(yintercept = 0, linewidth = 0.4, color = "grey40") +
    geom_col(width = 25, fill = "grey40", color = NA, na.rm = TRUE) +
    geom_point(aes(y = dot_y), color = dot_col, alpha = dot_alpha, size = dot_size, na.rm = TRUE) +
    scale_x_months(global_start, global_end, qbrks, ybrks) +
    scale_y_continuous(limits = c(ymin, ymax), sec.axis = secaxis_env_blank(pres_top, ymin, ymax)) +
    labs(y = "BEUTI\nanom", x = NULL) +
    theme_panel
  
  if (!show_x) p <- p + hide_x_labels
  p
}

plot_sst <- function(df, global_start, global_end, qbrks, ybrks, pres_top, show_x = TRUE) {
  yr <- env_range(df$SST_anom_detrended); ymin <- yr[1]; ymax <- yr[2]
  d <- df %>%
    mutate(
      sst_cat = case_when(
        is.na(SST_anom_detrended) ~ NA_character_,
        !is.na(MHW_days_pct_detrended) & MHW_days_pct_detrended > 0 ~ "mhw",
        SST_anom_detrended < 0 ~ "neg",
        TRUE ~ "pos"
      ),
      dot_y = map_dot_to_env(dot_hours, pres_top, ymin, ymax)
    )
  
  p <- ggplot(d, aes(x = month, y = SST_anom_detrended)) +
    geom_hline(yintercept = 0, linewidth = 0.4, color = "grey40") +
    geom_col(aes(fill = sst_cat), width = 25, color = NA, na.rm = TRUE) +
    scale_fill_manual(
      values = c(neg = col_neg, pos = col_pos, mhw = col_mhw),
      na.translate = FALSE
    ) +
    geom_point(aes(y = dot_y), color = dot_col, alpha = dot_alpha, size = dot_size, na.rm = TRUE) +
    scale_x_months(global_start, global_end, qbrks, ybrks) +
    scale_y_continuous(limits = c(ymin, ymax), sec.axis = secaxis_env_blank(pres_top, ymin, ymax)) +
    labs(y = "SSTA", x = NULL) +
    theme_panel
  
  if (!show_x) p <- p + hide_x_labels
  p
}

build_stack_site <- function(site_name) {
  fish_df <- mon_list[[site_name]]
  idx_df <- prep_idx_site(site_name, month_seq, fish_df)
  
  pres_top <- suppressWarnings(max(fish_df$avg_weekly_hours + fish_df$se_weekly_hours, na.rm = TRUE))
  if (!is.finite(pres_top) || pres_top <= 0) pres_top <- 1
  
  p_fish <- plot_fish_panel(
    fish_df,
    site_name,
    global_start,
    global_end,
    qbrks,
    ybrks,
    spl_breaks,
    spl_labels,
    spl_cols6
  )
  
  p_pdo <- plot_pdo(idx_df, global_start, global_end, qbrks, ybrks, pres_top, show_x = FALSE)
  p_npgo <- plot_npgo(idx_df, global_start, global_end, qbrks, ybrks, pres_top, show_x = FALSE)
  p_oni <- plot_oni(idx_df, global_start, global_end, qbrks, ybrks, pres_top, show_x = FALSE)
  p_beuti <- plot_beuti(idx_df, global_start, global_end, qbrks, ybrks, pres_top, show_x = FALSE)
  p_sst <- plot_sst(idx_df, global_start, global_end, qbrks, ybrks, pres_top, show_x = TRUE)
  
  (p_fish / p_pdo / p_npgo / p_oni / p_beuti / p_sst) +
    plot_layout(ncol = 1, heights = c(1, 0.75, 0.75, 0.75, 0.75, 0.75))
}

# =========================
# Build and save figures
# =========================

for (site_name in names(mon_list)) {
  p <- build_stack_site(site_name)
  
  ggsave(
    filename = file.path(output_dir, paste0(site_name, "_monthly_fish_environment_stack.png")),
    plot = p,
    width = 8.6,
    height = 8.6,
    dpi = 600
  )
  
  ggsave(
    filename = file.path(output_dir, paste0(site_name, "_monthly_fish_environment_stack.pdf")),
    plot = p,
    width = 8.6,
    height = 8.6,
    device = cairo_pdf
  )
}