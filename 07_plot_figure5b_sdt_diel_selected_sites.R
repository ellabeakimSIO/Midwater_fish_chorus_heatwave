# -----------------------------------------------------------------------------
# 07_plot_figure5b_sdt_diel_selected_sites.R
#
# Purpose:
# Generate Figure 5B: diel patterns of UF440 chorus detections during the 2017
# chorusing season at SDT_WQ, SDT_PR, and SDT_HP.
#
# In the paper:
# This script reproduces the Figure 5B workflow in Kim et al. (2026), showing
# detections at three San Diego Trough sites in one row, with detections
# colored by chorus SPL.
#
# Notes:
# - Uses the final filtered detection column when available:
#     UF440_PA_no20iso_noday
#   and otherwise falls back to:
#     UF440_PA_no20iso
# - Night shading is included for visual context.
# - No off-effort shading is applied for these SDT panels.
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(tidyr)
  library(lubridate)
  library(ggplot2)
  library(suncalc)
  library(patchwork)
  library(scales)
})

# =========================
# User settings
# =========================

data_dir <- "data"
output_dir <- "figures"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

tz_use <- "UTC"
bin_len <- minutes(20)

# Selected Figure 5B sites, in manuscript order
selected_sites <- c("SDT_WQ", "SDT_PR", "SDT_HP")

# Fixed SPL display scale used across diel figures
spl_limits <- c(55, 85)

# Plot styling
night_fill <- "#E9E0FF"
night_alpha <- 0.65
na_spl_color <- "grey70"

axis_text_x_size <- 16
axis_text_y_size <- 16
title_size <- 16
legend_text_size <- 11
legend_title_size <- 12

# SDT site coordinates (decimal degrees)
sdt_coords <- tibble::tribble(
  ~site,      ~lat,              ~lon,
  "SDT_WQ",   32 + 46.301/60,    -(117 + 47.905/60),
  "SDT_PR",   32 + 54.885/60,    -(117 + 29.772/60),
  "SDT_HP",   32 + 45.597/60,    -(117 + 39.302/60)
)

# Jet-like palette used in your working script
jet_cols <- grDevices::colorRampPalette(
  c("#00007F", "#0000FF", "#007FFF", "#00FFFF", "#7FFF7F",
    "#FFFF00", "#FF7F00", "#FF0000", "#7F0000")
)(256)

# =========================
# Helpers
# =========================

parse_dt_iso <- function(x, tz = "UTC") {
  if (inherits(x, "POSIXt")) {
    return(lubridate::force_tz(as.POSIXct(x), tzone = tz))
  }
  
  x <- trimws(as.character(x))
  x[x == ""] <- NA_character_
  
  x <- gsub("T", " ", x, fixed = TRUE)
  x <- sub("Z$", "", x)
  x <- sub("\\.\\d+", "", x)
  x <- sub("([+-]\\d\\d:?\\d\\d)$", "", x)
  
  out <- lubridate::ymd_hms(x, tz = tz, quiet = TRUE)
  miss <- is.na(out)
  if (any(miss)) out[miss] <- lubridate::ymd_hm(x[miss], tz = tz, quiet = TRUE)
  miss <- is.na(out)
  if (any(miss)) out[miss] <- lubridate::ymd(x[miss], tz = tz, quiet = TRUE)
  
  out
}

normalize_sdt_site <- function(x) {
  x <- as.character(x)
  x <- sub("_\\d+$", "", x)
  x
}

pick_detection_col <- function(df) {
  if ("UF440_PA_no20iso_noday" %in% names(df)) return("UF440_PA_no20iso_noday")
  if ("UF440_PA_no20iso" %in% names(df)) return("UF440_PA_no20iso")
  stop("No filtered UF440 detection column found. Expected UF440_PA_no20iso_noday or UF440_PA_no20iso.")
}

read_master_dir <- function(data_dir, tz = "UTC") {
  fns <- list.files(data_dir, pattern = "\\.csv$", full.names = TRUE, ignore.case = TRUE)
  if (length(fns) == 0) stop("No CSV files found in: ", data_dir)
  
  dfs <- map(fns, function(f) {
    df <- read_csv(f, show_col_types = FALSE)
    names(df) <- trimws(names(df))
    
    if (!("Site" %in% names(df))) return(NULL)
    
    ts_col <- if ("timestamp_dt" %in% names(df)) {
      "timestamp_dt"
    } else if ("timestamp_utc" %in% names(df)) {
      "timestamp_utc"
    } else {
      NA_character_
    }
    
    if (is.na(ts_col)) return(NULL)
    
    det_col <- tryCatch(pick_detection_col(df), error = function(e) NA_character_)
    if (is.na(det_col)) return(NULL)
    
    if (!("UF440_SPL" %in% names(df))) return(NULL)
    
    df %>%
      transmute(
        Site = normalize_sdt_site(as.character(.data$Site)),
        timestamp_dt = parse_dt_iso(.data[[ts_col]], tz = tz),
        detect = as.integer(.data[[det_col]]),
        SPL = suppressWarnings(as.numeric(.data$UF440_SPL))
      ) %>%
      filter(!is.na(timestamp_dt))
  })
  
  bind_rows(dfs)
}

solar_night_rects <- function(date_min, date_max, lat, lon, tz = "UTC") {
  dates <- seq(date_min, date_max, by = "day")
  
  solar <- getSunlightTimes(
    date = dates, lat = lat, lon = lon, tz = tz,
    keep = c("sunrise", "sunset")
  ) %>%
    mutate(
      date_num = as.numeric(as.Date(date, tz = tz)),
      sunrise_sec = as.numeric(lubridate::hms(format(sunrise, "%H:%M:%S"))),
      sunset_sec = as.numeric(lubridate::hms(format(sunset, "%H:%M:%S")))
    ) %>%
    select(date_num, sunrise_sec, sunset_sec)
  
  pmap_dfr(solar, function(date_num, sunrise_sec, sunset_sec) {
    if (is.na(sunrise_sec) || is.na(sunset_sec)) return(tibble())
    
    ymin <- date_num - 0.5
    ymax <- date_num + 0.5
    
    if (sunrise_sec <= sunset_sec) {
      bind_rows(
        tibble(xmin = 0, xmax = sunrise_sec, ymin = ymin, ymax = ymax),
        tibble(xmin = sunset_sec, xmax = 86400, ymin = ymin, ymax = ymax)
      )
    } else {
      tibble(xmin = sunset_sec, xmax = sunrise_sec, ymin = ymin, ymax = ymax)
    }
  })
}

bins_from_master <- function(master_df, bin_len = minutes(20)) {
  master_df %>%
    filter(!is.na(detect), detect == 1) %>%
    transmute(
      Site,
      start_time = timestamp_dt,
      end_time = timestamp_dt + bin_len,
      SPL = SPL
    )
}

explode_to_daily_segments_spl <- function(df, tz = "UTC") {
  df %>%
    filter(!is.na(start_time), !is.na(end_time)) %>%
    mutate(
      start_date = as.Date(start_time, tz = tz),
      end_date = as.Date(end_time, tz = tz),
      start_sec = as.numeric(lubridate::hms(format(start_time, "%H:%M:%S"))),
      end_sec = as.numeric(lubridate::hms(format(end_time, "%H:%M:%S")))
    ) %>%
    mutate(
      seg = pmap(
        list(Site, start_date, end_date, start_sec, end_sec, SPL),
        function(site, sd, ed, ss, es, spl) {
          days <- seq(sd, ed, by = "day")
          tibble(
            Site = site,
            date = days,
            start_sec = ifelse(days == sd, ss, 0),
            end_sec = ifelse(days == ed, es, 86400),
            SPL = spl
          ) %>%
            filter(end_sec > start_sec)
        }
      )
    ) %>%
    select(seg) %>%
    unnest(seg)
}

make_diel_plot_monthly_spl <- function(df_segments, lat, lon, title_text,
                                       spl_limits = NULL) {
  date_min <- min(df_segments$date, na.rm = TRUE)
  date_max <- max(df_segments$date, na.rm = TRUE)
  
  night <- solar_night_rects(date_min, date_max, lat, lon, tz = tz_use)
  
  brks_dates <- seq(
    floor_date(date_min, "month"),
    ceiling_date(date_max, "month"),
    by = "1 month"
  )
  brks_num <- as.numeric(brks_dates)
  labs_y <- format(brks_dates, "%Y-%m")
  
  ggplot() +
    geom_rect(
      data = night,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = night_fill,
      alpha = night_alpha,
      inherit.aes = FALSE
    ) +
    geom_segment(
      data = df_segments,
      aes(
        x = start_sec, xend = end_sec,
        y = as.numeric(date), yend = as.numeric(date),
        color = SPL
      ),
      linewidth = 0.8
    ) +
    scale_color_gradientn(
      colours = jet_cols,
      na.value = na_spl_color,
      limits = spl_limits,
      oob = scales::squish,
      name = "UF440 SPL"
    ) +
    scale_x_continuous(
      limits = c(0, 86400),
      breaks = seq(0, 86400, by = 4 * 3600),
      labels = function(x) as.character(as.integer(round(x / 3600))),
      expand = c(0, 0)
    ) +
    scale_y_reverse(
      breaks = brks_num,
      labels = labs_y,
      expand = c(0, 0)
    ) +
    labs(x = "Hour of day (UTC)", y = NULL, title = title_text) +
    theme_bw(base_size = 14) +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = axis_text_x_size),
      axis.text.y = element_text(size = axis_text_y_size),
      axis.title.x = element_text(size = axis_text_x_size),
      plot.title = element_text(face = "bold", size = title_size, hjust = 0.5),
      legend.position = "right",
      legend.title = element_text(size = legend_title_size),
      legend.text = element_text(size = legend_text_size)
    )
}

# =========================
# Read data
# =========================

master_all <- read_master_dir(data_dir, tz = tz_use) %>%
  filter(Site %in% selected_sites)

if (nrow(master_all) == 0) {
  stop("No matching SDT records found for selected sites.")
}

# =========================
# Build selected-site plots
# =========================

selected_plots <- map(selected_sites, function(site_nm) {
  df_site <- master_all %>% filter(Site == site_nm)
  
  if (nrow(df_site) == 0) {
    warning("No data found for ", site_nm)
    return(NULL)
  }
  
  ll <- sdt_coords %>% filter(site == site_nm)
  bins <- bins_from_master(df_site, bin_len = bin_len)
  seg <- explode_to_daily_segments_spl(bins, tz = tz_use)
  
  make_diel_plot_monthly_spl(
    df_segments = seg,
    lat = ll$lat,
    lon = ll$lon,
    title_text = site_nm,
    spl_limits = spl_limits
  )
})

names(selected_plots) <- selected_sites
selected_plots <- selected_plots[!vapply(selected_plots, is.null, logical(1))]

figure_5b <- wrap_plots(selected_plots, nrow = 1, ncol = 3, guides = "collect") &
  theme(legend.position = "right")

print(figure_5b)

# =========================
# Export
# =========================

ggsave(
  filename = file.path(output_dir, "Figure5B_SDT_diel_selected_sites.png"),
  plot = figure_5b,
  width = 14,
  height = 5,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = file.path(output_dir, "Figure5B_SDT_diel_selected_sites.svg"),
  plot = figure_5b,
  width = 14,
  height = 5,
  units = "in",
  bg = "white"
)

message("Saved Figure 5B outputs to: ", output_dir)