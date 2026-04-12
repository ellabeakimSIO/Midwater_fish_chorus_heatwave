# -----------------------------------------------------------------------------
# S02_plot_supplemental_figure2_diel_focal_sites.R
#
# Purpose:
# Generate Supplemental Figure 2: diel patterns of UF440 chorus detections
# across all deployments at the three multi-year focal sites (North_B, South_P,
# South_T).
#
# Figure elements:
#   - one horizontal line per detected 20 min chorus bin
#   - detections colored by UF440 SPL
#   - light purple shading for night
#   - gray shading for periods without acoustic recording effort
#
# In the paper:
# This script reproduces Supplemental Figure 2 in Kim et al. (2026),
# "Diel patterns of UF440 chorus at multi-year focal sites."
#
# Notes:
# - This script assumes analysis-ready master tables with the columns:
#   Site, timestamp_dt, UF440_roving_PSD, UF440_PA_no20iso_noday, UF440_SPL
# - The final manuscript references Supplemental Figure 2 when discussing
#   nocturnal chorus timing and daytime false positives.
# - This script uses the final filtered detection column
#   UF440_PA_no20iso_noday.
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
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

site_files <- tibble::tribble(
  ~Site,      ~file,
  "North_B",  file.path(data_dir, "North_B_master_dryad.csv"),
  "South_P",  file.path(data_dir, "South_P_master_dryad.csv"),
  "South_T",  file.path(data_dir, "South_T_master_dryad.csv")
)

# Site coordinates (decimal degrees)
site_coords <- tibble::tribble(
  ~Site,      ~lat,        ~lon,
  "North_B",  34.250000,  -120.026000,
  "South_P",  32.893000,  -117.378000,
  "South_T",  32.887000,  -117.556000
)

tz_use <- "UTC"
bin_len_mins <- 20

night_fill <- "#E9E0FF"
night_alpha <- 0.65
offeffort_fill <- "grey85"
offeffort_alpha <- 1
na_spl_color <- "grey70"

# Fixed SPL range used for display consistency
spl_limits <- c(55, 85)

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

read_master_file <- function(fp, site_name, tz = "UTC") {
  df <- read_delim_auto(fp)
  
  required_cols <- c("timestamp_dt", "UF440_roving_PSD", "UF440_PA_no20iso_noday", "UF440_SPL")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns in ", basename(fp), ": ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  df %>%
    transmute(
      Site = site_name,
      timestamp_dt = parse_dt_iso(timestamp_dt, tz = tz),
      UF440_roving_PSD = suppressWarnings(as.numeric(UF440_roving_PSD)),
      UF440_PA_no20iso_noday = suppressWarnings(as.numeric(UF440_PA_no20iso_noday)),
      UF440_SPL = suppressWarnings(as.numeric(UF440_SPL))
    ) %>%
    mutate(across(where(is.numeric), ~ ifelse(is.nan(.) | is.infinite(.), NA_real_, .))) %>%
    filter(!is.na(timestamp_dt))
}

make_effort_gaps_from_master <- function(df_site, tz = "UTC") {
  df_days <- df_site %>%
    mutate(date = as.Date(timestamp_dt, tz = tz)) %>%
    group_by(date) %>%
    summarise(on_effort = any(!is.na(UF440_roving_PSD)), .groups = "drop")
  
  all_days <- tibble(date = seq(min(df_days$date), max(df_days$date), by = "day")) %>%
    left_join(df_days, by = "date") %>%
    mutate(on_effort = replace_na(on_effort, FALSE))
  
  all_days %>%
    filter(!on_effort) %>%
    mutate(
      xmin = 0,
      xmax = 86400,
      ymin = as.numeric(date) - 0.5,
      ymax = as.numeric(date) + 0.5
    )
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

bins_from_master <- function(master_df, pa_col = "UF440_PA_no20iso_noday", bin_len_mins = 20) {
  master_df %>%
    filter(!is.na(.data[[pa_col]]), .data[[pa_col]] == 1) %>%
    transmute(
      Site,
      start_time = timestamp_dt,
      end_time = timestamp_dt + minutes(bin_len_mins),
      SPL = UF440_SPL
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
    mutate(seg = pmap(
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
    )) %>%
    select(seg) %>%
    unnest(seg)
}

jet_cols <- grDevices::colorRampPalette(
  c("#00007F", "#0000FF", "#007FFF", "#00FFFF", "#7FFF7F",
    "#FFFF00", "#FF7F00", "#FF0000", "#7F0000")
)(256)

make_diel_plot_long_spl <- function(df_segments, lat, lon, title_text,
                                    plot_start, plot_end,
                                    offeffort_rects = NULL,
                                    spl_limits = NULL,
                                    label_april_only = FALSE) {
  df_segments_plot <- df_segments %>%
    filter(date >= plot_start, date <= plot_end)
  
  night <- solar_night_rects(plot_start, plot_end, lat, lon, tz = tz_use)
  
  if (!is.null(offeffort_rects) && nrow(offeffort_rects) > 0) {
    offeffort_rects_plot <- offeffort_rects %>%
      filter(date >= plot_start, date <= plot_end)
  } else {
    offeffort_rects_plot <- NULL
  }
  
  brks_dates <- seq(from = plot_start, to = plot_end, by = "6 months")
  brks_num <- as.numeric(brks_dates)
  
  labs_y <- if (label_april_only) {
    ifelse(lubridate::month(brks_dates) == 4, format(brks_dates, "%Y-%m"), "")
  } else {
    format(brks_dates, "%Y-%m")
  }
  
  y_limits <- c(as.numeric(plot_end) + 0.5, as.numeric(plot_start) - 0.5)
  
  ggplot() +
    geom_rect(
      data = night,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = night_fill,
      alpha = night_alpha,
      inherit.aes = FALSE
    ) +
    geom_segment(
      data = df_segments_plot,
      aes(x = start_sec, xend = end_sec,
          y = as.numeric(date), yend = as.numeric(date),
          color = SPL),
      linewidth = 0.8
    ) +
    {
      if (!is.null(offeffort_rects_plot) && nrow(offeffort_rects_plot) > 0) {
        geom_rect(
          data = offeffort_rects_plot,
          aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
          fill = offeffort_fill,
          alpha = offeffort_alpha,
          inherit.aes = FALSE
        )
      }
    } +
    scale_color_gradientn(
      colours = jet_cols,
      na.value = na_spl_color,
      limits = spl_limits,
      oob = scales::squish,
      name = "UF440 SPL"
    ) +
    scale_x_continuous(
      limits = c(0, 86400),
      breaks = seq(0, 86400, by = 3600),
      labels = function(x) {
        h <- as.integer(round(x / 3600))
        ifelse(h %% 2 == 0, as.character(h), "")
      },
      expand = c(0, 0)
    ) +
    scale_y_reverse(
      limits = y_limits,
      breaks = brks_num,
      labels = labs_y,
      expand = c(0, 0)
    ) +
    labs(x = NULL, y = NULL, title = title_text) +
    theme_bw(base_size = 14) +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 0),
      legend.position = "none",
      plot.title = element_text(face = "bold", hjust = 0, size = 12)
    )
}

# =========================
# Read and prepare data
# =========================

master_all <- purrr::map2_dfr(site_files$file, site_files$Site, read_master_file, tz = tz_use)

plot_list <- purrr::map(site_files$Site, function(site_name) {
  df_site <- master_all %>% filter(Site == site_name)
  coords <- site_coords %>% filter(Site == site_name)
  
  bins <- bins_from_master(df_site, pa_col = "UF440_PA_no20iso_noday", bin_len_mins = bin_len_mins)
  seg <- explode_to_daily_segments_spl(bins, tz = tz_use)
  
  effort_gaps <- make_effort_gaps_from_master(df_site, tz = tz_use)
  
  plot_start <- min(as.Date(df_site$timestamp_dt, tz = tz_use), na.rm = TRUE)
  plot_end <- max(as.Date(df_site$timestamp_dt, tz = tz_use), na.rm = TRUE)
  
  label_april_only <- identical(site_name, "North_B")
  
  make_diel_plot_long_spl(
    df_segments = seg,
    lat = coords$lat,
    lon = coords$lon,
    title_text = site_name,
    plot_start = plot_start,
    plot_end = plot_end,
    offeffort_rects = effort_gaps,
    spl_limits = spl_limits,
    label_april_only = label_april_only
  )
})

names(plot_list) <- site_files$Site

# Split North_B into two panels if desired for visual balance
# For GitHub simplicity, keep one panel per site in a vertical stack.
supp_figure_2 <- plot_list[["North_B"]] / plot_list[["South_T"]] / plot_list[["South_P"]]

# =========================
# Save outputs
# =========================

png_fp <- file.path(output_dir, "Supplemental_Figure2_diel_focal_sites.png")
pdf_fp <- file.path(output_dir, "Supplemental_Figure2_diel_focal_sites.pdf")

ggsave(png_fp, supp_figure_2, width = 10, height = 12, units = "in", dpi = 600, bg = "white")
ggsave(pdf_fp, supp_figure_2, width = 10, height = 12, units = "in")

message("Saved:\n  ", png_fp, "\n  ", pdf_fp)