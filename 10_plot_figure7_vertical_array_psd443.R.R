# -----------------------------------------------------------------------------
# 10_plot_figure7_vertical_array_psd443.R
#
# Purpose:
# Generate Figure 7 for the Vertical_DM array using 1 min-binned PSD at
# 443 Hz for hydrophones at 100, 200, and 600 m depth.
#
# In the paper:
# This script reproduces the Figure 7 workflow in Kim et al. (2026):
#   - Panel A: full 24 h cycle, minute-binned median SPL shown as points
#   - Panels B-C: expanded sunset and sunrise views, points only
#   - Panels D-E: expanded sunset and sunrise views, median line with IQR ribbon
#
# Notes:
# - Inputs are PSD CSV files exported from Triton / Soundscape Metrics.
# - This script uses 1 min bins and a fixed frequency of 443 Hz.
# - Depths are inferred from file names when possible; if not, use the
#   manual_depth_lookup object below.
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(readr)
  library(patchwork)
})

# =========================
# User settings
# =========================

data_dir <- "data"
output_dir <- "figures"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# If file names do not clearly contain depths, define them here:
# Example:
# manual_depth_lookup <- c(
#   "VerticalDM_100m" = 100,
#   "VerticalDM_200m" = 200,
#   "VerticalDM_600m" = 600
# )
manual_depth_lookup <- c()

# Figure windows (UTC), matching the final code branch
sunset_window <- c(1, 4)
sunrise_window <- c(10, 14)

# Y-limits from your late-stage working version
ylim_all <- c(58, 78)
ylim_sunset <- c(62, 80)
ylim_sunrise <- c(60, 75)

# Optional small padding for zoom windows
xpad_hours <- 0.08

# Colors matching your late-stage code
depth_cols <- c(
  "100" = "#F8766D",
  "200" = "#00BA38",
  "600" = "#619CFF"
)

# =========================
# Find files
# =========================

csv_files <- list.files(data_dir, pattern = "\\.csv$", full.names = TRUE)

if (length(csv_files) == 0) {
  stop("No CSV files found in: ", data_dir)
}

# =========================
# Helpers
# =========================

extract_depth_from_name <- function(file_stem) {
  if (length(manual_depth_lookup) > 0 && file_stem %in% names(manual_depth_lookup)) {
    return(as.numeric(manual_depth_lookup[[file_stem]]))
  }
  
  hit <- str_extract(file_stem, "(?i)(?<!\\d)(100|200|600)\\s*m?(?!\\d)")
  if (!is.na(hit)) {
    return(as.numeric(str_extract(hit, "\\d+")))
  }
  
  NA_real_
}

read_one_psd443 <- function(f) {
  message("Reading: ", basename(f))
  
  hdr <- names(fread(f, nrows = 0, showProgress = FALSE))
  hdr_trim <- trimws(hdr)
  
  time_col <- hdr[1]
  psd_idx  <- which(hdr_trim == "PSD_443")
  
  if (length(psd_idx) != 1) {
    stop("Could not uniquely find PSD_443 in file: ", basename(f))
  }
  
  psd_col <- hdr[psd_idx]
  
  dat <- fread(
    f,
    select = c(time_col, psd_col),
    showProgress = FALSE
  )
  
  setnames(dat, old = c(time_col, psd_col), new = c("datetime_utc", "PSD_443"))
  
  file_stem <- tools::file_path_sans_ext(basename(f))
  depth_m <- extract_depth_from_name(file_stem)
  
  dat[, datetime_utc := as.POSIXct(
    datetime_utc,
    format = "%Y-%m-%dT%H:%M:%OSZ",
    tz = "UTC"
  )]
  
  # fallback parse if needed
  bad <- is.na(dat$datetime_utc)
  if (any(bad)) {
    dat[bad, datetime_utc := as.POSIXct(datetime_utc, tz = "UTC")]
  }
  
  dat[, file_id := file_stem]
  dat[, depth_m := depth_m]
  dat[, depth_label := ifelse(is.na(depth_m), file_id, as.character(depth_m))]
  
  as_tibble(dat)
}

make_point_plot <- function(df, xlim_hours, title_text, ylim_vals, xpad_hours = 0) {
  x_range <- diff(xlim_hours)
  
  x_breaks <- if (x_range <= 4) {
    seq(xlim_hours[1], xlim_hours[2], by = 1)
  } else if (x_range <= 8) {
    seq(xlim_hours[1], xlim_hours[2], by = 1)
  } else {
    seq(0, 24, by = 2)
  }
  
  ggplot(
    df,
    aes(x = hour_of_day_utc_1min, y = median_psd, colour = depth_label)
  ) +
    geom_point(size = 1.2, alpha = 0.8) +
    scale_colour_manual(values = depth_cols, drop = FALSE) +
    scale_x_continuous(
      breaks = x_breaks,
      expand = c(0, 0)
    ) +
    coord_cartesian(
      xlim = c(xlim_hours[1] - xpad_hours, xlim_hours[2] + xpad_hours),
      ylim = ylim_vals
    ) +
    labs(
      x = "Hour of Day (UTC)",
      y = expression("PSD at 443 Hz (dB re 1 " * mu * "Pa"^2 * "/Hz)"),
      colour = "Hydrophone Depth (m)",
      title = title_text
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "right"
    )
}

make_ribbon_plot <- function(df, xlim_hours, title_text, ylim_vals, xpad_hours = 0) {
  x_range <- diff(xlim_hours)
  
  x_breaks <- if (x_range <= 4) {
    seq(xlim_hours[1], xlim_hours[2], by = 1)
  } else if (x_range <= 8) {
    seq(xlim_hours[1], xlim_hours[2], by = 1)
  } else {
    seq(0, 24, by = 2)
  }
  
  ggplot(
    df,
    aes(x = hour_of_day_utc_1min, group = depth_label)
  ) +
    geom_ribbon(
      aes(ymin = q25, ymax = q75, fill = depth_label),
      alpha = 0.18,
      colour = NA
    ) +
    geom_line(
      aes(y = median_psd, colour = depth_label),
      linewidth = 0.8
    ) +
    scale_colour_manual(values = depth_cols, drop = FALSE) +
    scale_fill_manual(values = depth_cols, drop = FALSE) +
    scale_x_continuous(
      breaks = x_breaks,
      expand = c(0, 0)
    ) +
    coord_cartesian(
      xlim = c(xlim_hours[1] - xpad_hours, xlim_hours[2] + xpad_hours),
      ylim = ylim_vals
    ) +
    labs(
      x = "Hour of Day (UTC)",
      y = expression("PSD at 443 Hz (dB re 1 " * mu * "Pa"^2 * "/Hz)"),
      colour = "Hydrophone Depth (m)",
      fill   = "Hydrophone Depth (m)",
      title  = title_text
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "right"
    )
}

# =========================
# Read files and summarize to 1 min bins
# =========================

psd443_list <- lapply(csv_files, read_one_psd443)
names(psd443_list) <- tools::file_path_sans_ext(basename(csv_files))

all_psd443 <- bind_rows(psd443_list) %>%
  filter(!is.na(datetime_utc), !is.na(PSD_443)) %>%
  mutate(
    hour   = as.integer(format(datetime_utc, "%H")),
    minute = as.integer(format(datetime_utc, "%M")),
    minute_of_day = hour * 60L + minute,
    hour_of_day_utc_1min = minute_of_day / 60
  )

depth_order <- all_psd443 %>%
  distinct(depth_label, depth_m) %>%
  mutate(depth_sort = ifelse(is.na(depth_m), Inf, depth_m)) %>%
  arrange(depth_sort, depth_label) %>%
  pull(depth_label)

all_psd443 <- all_psd443 %>%
  mutate(depth_label = factor(depth_label, levels = depth_order))

psd443_1min_summary <- all_psd443 %>%
  group_by(depth_label, depth_m, minute_of_day, hour_of_day_utc_1min) %>%
  summarise(
    n          = n(),
    median_psd = median(PSD_443, na.rm = TRUE),
    q25        = quantile(PSD_443, 0.25, na.rm = TRUE),
    q75        = quantile(PSD_443, 0.75, na.rm = TRUE),
    .groups    = "drop"
  )

write_csv(
  psd443_1min_summary,
  file.path(output_dir, "Figure7_psd443_1min_summary.csv")
)

# =========================
# Optional peak timing / depth differences
# =========================

calc_peak_by_window <- function(df, start_hour, end_hour, window_name) {
  df %>%
    filter(hour_of_day_utc_1min >= start_hour,
           hour_of_day_utc_1min <= end_hour) %>%
    group_by(depth_label, depth_m) %>%
    slice_max(order_by = median_psd, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(
      window = window_name,
      peak_hour = floor(hour_of_day_utc_1min),
      peak_min  = round((hour_of_day_utc_1min - peak_hour) * 60),
      peak_time_utc = sprintf("%02d:%02d", peak_hour, peak_min)
    ) %>%
    select(window, depth_label, depth_m, peak_time_utc, peak_median_psd = median_psd)
}

peak_sunset <- calc_peak_by_window(psd443_1min_summary, sunset_window[1], sunset_window[2], "sunset")
peak_sunrise <- calc_peak_by_window(psd443_1min_summary, sunrise_window[1], sunrise_window[2], "sunrise")

peak_summary <- bind_rows(peak_sunset, peak_sunrise) %>%
  arrange(window, depth_m)

write_csv(
  peak_summary,
  file.path(output_dir, "Figure7_peak_summary.csv")
)

# =========================
# Build figure panels
# =========================

p_A <- make_point_plot(
  psd443_1min_summary,
  xlim_hours = c(0, 24),
  title_text = "A",
  ylim_vals = ylim_all
)

p_B <- make_point_plot(
  filter(psd443_1min_summary,
         hour_of_day_utc_1min >= sunset_window[1],
         hour_of_day_utc_1min <= sunset_window[2]),
  xlim_hours = sunset_window,
  title_text = "B",
  ylim_vals = ylim_sunset,
  xpad_hours = xpad_hours
)

p_C <- make_point_plot(
  filter(psd443_1min_summary,
         hour_of_day_utc_1min >= sunrise_window[1],
         hour_of_day_utc_1min <= sunrise_window[2]),
  xlim_hours = sunrise_window,
  title_text = "C",
  ylim_vals = ylim_sunrise,
  xpad_hours = xpad_hours
)

p_D <- make_ribbon_plot(
  filter(psd443_1min_summary,
         hour_of_day_utc_1min >= sunset_window[1],
         hour_of_day_utc_1min <= sunset_window[2]),
  xlim_hours = sunset_window,
  title_text = "D",
  ylim_vals = ylim_sunset,
  xpad_hours = xpad_hours
)

p_E <- make_ribbon_plot(
  filter(psd443_1min_summary,
         hour_of_day_utc_1min >= sunrise_window[1],
         hour_of_day_utc_1min <= sunrise_window[2]),
  xlim_hours = sunrise_window,
  title_text = "E",
  ylim_vals = ylim_sunrise,
  xpad_hours = xpad_hours
)

figure_7 <- (p_A) / (p_B | p_C) / (p_D | p_E) +
  plot_layout(guides = "collect", heights = c(1.1, 1, 1)) &
  theme(legend.position = "right")

print(figure_7)

# =========================
# Export
# =========================

ggsave(
  filename = file.path(output_dir, "Figure7_vertical_array_psd443.png"),
  plot = figure_7,
  width = 10,
  height = 12,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = file.path(output_dir, "Figure7_vertical_array_psd443.svg"),
  plot = figure_7,
  width = 10,
  height = 12,
  units = "in",
  bg = "white"
)

message("Saved Figure 7 outputs to: ", output_dir)
