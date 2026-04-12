# -----------------------------------------------------------------------------
# 02_plot_figure3_monthly_focal_timeseries.R
#
# Purpose:
# Generate Figure 3: monthly UF440 chorus time series at the three multi-year
# focal sites (North_B, South_P, South_T).
#
# Figure elements:
#   - monthly chorus presence shown as the mean of weekly total chorus hours
#   - error bars showing ± 1 SE across weeks within each month
#   - bar color showing mean sunset chorus SPL for that month
#   - gray bars for months with chorus detections but no first-hour-after-sunset SPL
#   - black points showing monthly recording effort
#   - light gray shading for periods without recording effort
#
# In the paper:
# This script reproduces Figure 3 in Kim et al. (2026), "Time series of UF440
# chorus at multi-year focal sites."
#
# Notes:
# - Input files are analysis-ready master tables exported after removal of
#   isolated detections and daytime detections.
# - This script uses the current manuscript site names: North_B, South_P,
#   and South_T.
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
  library(purrr)
})

# =========================
# User settings
# =========================

data_dir <- "data"
output_dir <- "figures"

site_files <- tibble::tribble(
  ~Site,      ~file,
  "North_B",  file.path(data_dir, "North_B_master_dryad.csv"),
  "South_P",  file.path(data_dir, "South_P_master_dryad.csv"),
  "South_T",  file.path(data_dir, "South_T_master_dryad.csv")
)

bin_mins <- 20
week_start_day <- 1  # 1 = Monday

page_w_in <- 11
page_h_in <- 8.5

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
  
  out <- suppressWarnings(as.POSIXct(x_chr, tz = tz, format = "%Y-%m-%dT%H:%M:%OSZ"))
  bad <- is.na(out) & !is.na(x_chr)
  
  if (any(bad)) {
    out[bad] <- suppressWarnings(as.POSIXct(x_chr[bad], tz = tz, format = "%Y-%m-%dT%H:%M:%OS"))
  }
  
  bad <- is.na(out) & !is.na(x_chr)
  if (any(bad)) {
    out[bad] <- suppressWarnings(ymd_hms(x_chr[bad], tz = tz, quiet = TRUE))
  }
  
  bad <- is.na(out) & !is.na(x_chr)
  if (any(bad)) {
    out[bad] <- suppressWarnings(parse_date_time(
      x_chr[bad],
      orders = c("Y-m-d H:M:S", "Y/m/d H:M:S", "ymd HMS", "ymd HM", "mdy HMS", "mdy HM"),
      tz = tz
    ))
  }
  
  out
}

read_master_file <- function(fp, site_name) {
  df <- read_delim_auto(fp)
  
  required_cols <- c(
    "timestamp_dt",
    "UF440_roving_PSD",
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
      timestamp_dt = parse_datetime_flexibly(timestamp_dt, tz = "UTC"),
      UF440_roving_PSD = readr::parse_double(as.character(UF440_roving_PSD), na = c("", "NA", "NaN", "NAN", "Inf", "-Inf")),
      UF440_PA_no20iso_noday = readr::parse_double(as.character(UF440_PA_no20iso_noday), na = c("", "NA", "NaN", "NAN", "Inf", "-Inf")),
      UF440_sunset_SPL = readr::parse_double(as.character(UF440_sunset_SPL), na = c("", "NA", "NaN", "NAN", "Inf", "-Inf"))
    ) %>%
    mutate(across(where(is.double), ~ replace(.x, is.nan(.x) | is.infinite(.x), NA_real_))) %>%
    filter(!is.na(timestamp_dt)) %>%
    select(Site, timestamp_dt, UF440_roving_PSD, UF440_PA_no20iso_noday, UF440_sunset_SPL)
}

weekly_summary <- function(df_site, bin_mins = 20, week_start_day = 1) {
  expected_bins_per_week <- as.integer(7 * 24 * 60 / bin_mins)
  
  df_site %>%
    mutate(
      week = as.Date(floor_date(timestamp_dt, unit = "week", week_start = week_start_day)),
      on_effort = !is.na(UF440_roving_PSD),
      pa_eff = if_else(on_effort, UF440_PA_no20iso_noday, NA_real_)
    ) %>%
    group_by(Site, week) %>%
    summarise(
      hours_week = sum(replace_na(pa_eff, 0), na.rm = TRUE) * (bin_mins / 60),
      effort_bins = sum(on_effort, na.rm = TRUE),
      effort_frac = effort_bins / expected_bins_per_week,
      .groups = "drop"
    )
}

monthly_summary_complete <- function(df_site, month_seq, bin_mins = 20, week_start_day = 1) {
  wk <- weekly_summary(df_site, bin_mins = bin_mins, week_start_day = week_start_day) %>%
    mutate(month = as.Date(floor_date(week, "month")))
  
  mon_spl <- df_site %>%
    mutate(
      month = as.Date(floor_date(timestamp_dt, "month")),
      on_effort = !is.na(UF440_roving_PSD)
    ) %>%
    group_by(Site, month) %>%
    summarise(
      month_sunset_mean = if (all(is.na(UF440_sunset_SPL[on_effort]))) NA_real_ else mean(UF440_sunset_SPL[on_effort], na.rm = TRUE),
      month_effort_bins = sum(on_effort, na.rm = TRUE),
      bins_present = sum(UF440_PA_no20iso_noday == 1, na.rm = TRUE),
      .groups = "drop"
    )
  
  wk_mon <- wk %>%
    group_by(Site, month) %>%
    summarise(
      avg_weekly_hours = mean(hours_week, na.rm = TRUE),
      n_weeks = sum(!is.na(hours_week)),
      se_weekly_hours = {
        n <- sum(!is.na(hours_week))
        if (n <= 1) 0 else sd(hours_week, na.rm = TRUE) / sqrt(n)
      },
      mean_effort_frac = mean(effort_frac, na.rm = TRUE),
      .groups = "drop"
    )
  
  mon <- wk_mon %>%
    left_join(mon_spl, by = c("Site", "month")) %>%
    complete(month = month_seq) %>%
    mutate(Site = first(na.omit(df_site$Site)))
  
  effort_months <- mon_spl %>% filter(month_effort_bins > 0) %>% pull(month)
  if (length(effort_months) > 0) {
    first_eff <- min(effort_months, na.rm = TRUE)
    last_eff <- max(effort_months, na.rm = TRUE)
  } else {
    ts_months <- df_site %>%
      mutate(month = as.Date(floor_date(timestamp_dt, "month"))) %>%
      pull(month)
    first_eff <- min(ts_months, na.rm = TRUE)
    last_eff <- max(ts_months, na.rm = TRUE)
  }
  
  mon %>%
    mutate(
      month_end = month %m+% months(1),
      in_effort_span = month >= first_eff & month <= last_eff,
      month_effort_bins_plot = coalesce(month_effort_bins, 0L),
      mean_effort_frac_plot = coalesce(mean_effort_frac, 0),
      full_off = month_effort_bins_plot == 0,
      partial = in_effort_span & mean_effort_frac_plot > 0 & mean_effort_frac_plot < 1,
      avg_weekly_hours = if_else(in_effort_span, avg_weekly_hours, NA_real_),
      se_weekly_hours = if_else(in_effort_span, se_weekly_hours, NA_real_),
      mean_effort_frac = if_else(in_effort_span, mean_effort_frac_plot, NA_real_),
      # grey fill for months with detections but no sunset SPL
      grey_bar = !is.na(avg_weekly_hours) & avg_weekly_hours > 0 & is.na(month_sunset_mean)
    ) %>%
    select(-in_effort_span, -mean_effort_frac_plot)
}

quarter_breaks <- function(start_date, end_date) {
  seq(floor_date(start_date, "quarter"), ceiling_date(end_date, "quarter"), by = "3 months")
}

year_breaks <- function(start_date, end_date) {
  seq(floor_date(start_date, "year"), ceiling_date(end_date, "year"), by = "1 year")
}

make_spl_bin_fixed <- function(x, breaks, labels, levels_all) {
  b <- cut(x, breaks = breaks, include.lowest = TRUE, right = FALSE, labels = labels)
  b <- as.character(b)
  b[is.na(b)] <- "NA"
  factor(b, levels = levels_all)
}

plot_site_monthly <- function(mon, site_name,
                              global_start, global_end,
                              qbrks, ybrks,
                              breaks, labels, levels_all, cols_map) {
  if (nrow(mon) == 0) return(NULL)
  
  y_top <- suppressWarnings(max(mon$avg_weekly_hours + mon$se_weekly_hours, na.rm = TRUE))
  if (!is.finite(y_top) || y_top <= 0) y_top <- 1
  
  off_df <- mon %>% filter(full_off)
  dot_df <- mon %>%
    filter(!is.na(mean_effort_frac)) %>%
    mutate(eff_y = mean_effort_frac * y_top)
  
  mon <- mon %>%
    mutate(
      spl_bin = make_spl_bin_fixed(month_sunset_mean, breaks, labels, levels_all),
      spl_bin = if_else(grey_bar, "NA", as.character(spl_bin)),
      spl_bin = factor(spl_bin, levels = levels_all)
    )
  
  ggplot(mon, aes(x = month, y = avg_weekly_hours)) +
    geom_rect(
      data = off_df,
      aes(xmin = month, xmax = month_end, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE,
      fill = "grey90",
      alpha = 0.85
    ) +
    geom_col(aes(fill = spl_bin), width = 25, color = NA, linewidth = 0, na.rm = TRUE) +
    geom_errorbar(
      aes(
        ymin = pmax(avg_weekly_hours - se_weekly_hours, 0),
        ymax = avg_weekly_hours + se_weekly_hours
      ),
      width = 10,
      linewidth = 0.4,
      color = "grey25",
      na.rm = TRUE
    ) +
    geom_point(
      data = dot_df,
      aes(x = month, y = eff_y),
      inherit.aes = FALSE,
      color = "black",
      size = 0.6
    ) +
    scale_x_date(
      breaks = ybrks,
      labels = label_date("%Y"),
      minor_breaks = qbrks,
      limits = c(global_start, global_end %m+% months(1)),
      expand = expansion(mult = c(0, 0))
    ) +
    scale_fill_manual(values = cols_map, drop = FALSE, guide = "none") +
    scale_y_continuous(
      name = "Presence",
      expand = expansion(mult = c(0, 0.05)),
      sec.axis = sec_axis(~ . / y_top * 100, name = "Effort (%)")
    ) +
    labs(title = site_name, x = NULL) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0, size = 12),
      axis.text.x = element_text(size = 10),
      axis.ticks.length = unit(3, "pt"),
      axis.ticks.x = element_line(),
      axis.ticks.x.minor = element_line(),
      plot.margin = margin(t = 6, r = 10, b = 6, l = 10, unit = "pt")
    )
}

# =========================
# Read and combine data
# =========================

all_raw <- purrr::map2_dfr(site_files$file, site_files$Site, read_master_file) %>%
  distinct(Site, timestamp_dt, .keep_all = TRUE)

if (nrow(all_raw) == 0) {
  stop("No rows found after reading master tables.")
}

north_df <- all_raw %>% filter(Site == "North_B")
if (nrow(north_df) > 0) {
  global_start <- as.Date(floor_date(min(north_df$timestamp_dt), "month"))
  global_end <- as.Date(floor_date(max(north_df$timestamp_dt), "month"))
} else {
  global_start <- as.Date(floor_date(min(all_raw$timestamp_dt), "month"))
  global_end <- as.Date(floor_date(max(all_raw$timestamp_dt), "month"))
}

month_seq <- seq(global_start, global_end, by = "1 month")
qbrks <- quarter_breaks(global_start, global_end)
ybrks <- year_breaks(global_start, global_end)

mon_list <- setNames(vector("list", length(site_files$Site)), site_files$Site)

for (s in site_files$Site) {
  df_s <- all_raw %>% filter(Site == s)
  mon_list[[s]] <- monthly_summary_complete(
    df_s,
    month_seq = month_seq,
    bin_mins = bin_mins,
    week_start_day = week_start_day
  )
}

all_month_means <- bind_rows(mon_list) %>% pull(month_sunset_mean)
gmin <- suppressWarnings(min(all_month_means, na.rm = TRUE))
gmax <- suppressWarnings(max(all_month_means, na.rm = TRUE))

if (!is.finite(gmin) || !is.finite(gmax) || gmin == gmax) {
  stop("Bad global min/max for monthly sunset SPL bins.")
}

spl_breaks <- seq(gmin, gmax, length.out = 6)
spl_breaks[6] <- spl_breaks[6] + 1e-6
spl_labels <- sprintf("%.1f", spl_breaks[1:5])
spl_levels <- c(spl_labels, "NA")

# blue -> red gradient, plus gray for NA / no sunset SPL
spl_cols6 <- c("blue", "cyan", "yellow", "orange", "red", "grey50")
names(spl_cols6) <- spl_levels

# =========================
# Build figure
# =========================

p1 <- plot_site_monthly(
  mon_list[["North_B"]], "North_B",
  global_start, global_end, qbrks, ybrks,
  spl_breaks, spl_labels, spl_levels, spl_cols6
)

p2 <- plot_site_monthly(
  mon_list[["South_P"]], "South_P",
  global_start, global_end, qbrks, ybrks,
  spl_breaks, spl_labels, spl_levels, spl_cols6
)

p3 <- plot_site_monthly(
  mon_list[["South_T"]], "South_T",
  global_start, global_end, qbrks, ybrks,
  spl_breaks, spl_labels, spl_levels, spl_cols6
)

figure_3 <- (p1 / p2 / p3) +
  plot_layout(ncol = 1, heights = c(1, 1, 1))

# =========================
# Save outputs
# =========================

png_fp <- file.path(output_dir, "Figure3_UF440_monthly_timeseries_focal_sites.png")
pdf_fp <- file.path(output_dir, "Figure3_UF440_monthly_timeseries_focal_sites.pdf")

ggsave(png_fp, figure_3, width = page_w_in, height = page_h_in, units = "in", dpi = 600, bg = "white")
ggsave(pdf_fp, figure_3, width = page_w_in, height = page_h_in, units = "in")

message("Saved:\n  ", png_fp, "\n  ", pdf_fp)