# -----------------------------------------------------------------------------
# 08_fit_plot_figure5c_sdt_depth_distance_gams.R
#
# Purpose:
# Fit the two GAMs underlying Figure 5C in Kim et al. (2026), relating
# maximum daily sunset chorus SPL to:
#   1) distance from coast, and
#   2) instrument depth,
# with a smooth term for date in each model.
#
# In the paper:
# This script reproduces the analyses described in Supplemental Table 4.
# Maximum daily sunset SPL is modeled with a Gaussian family (identity link)
# as a function of either distance from coast or instrument depth, with date
# included as a separate smooth term in each model. All smooths use thin-plate
# shrinkage splines with k = 4. A 14-day systematic subsample is used to
# reduce temporal autocorrelation.
#
# Inputs:
#   - data/SDT_ALL_master_dryad.csv
#   - data/South_T_master_dryad.csv
#   - data/South_P_master_dryad.csv
#   - data/SDT_SiteMetadata.csv NOTE this is not in the dryad folder but is in the github repository
#
# Outputs:
#   - figures/Figure5C_depth_distance_gams.png
#   - figures/Figure5C_depth_distance_gams.pdf
#   - figures/Figure5C_depth_model_predictions.csv
#   - figures/Figure5C_distance_model_predictions.csv
#   - figures/Figure5C_analysis_data.csv
#   - figures/Figure5C_14day_subsample.csv
#   - figures/Figure5C_model_summary.txt
#   - figures/Figure5C_supplement_table4.csv
#
# Notes:
# - Site metadata must include depth and coordinates for all analyzed sites,
#   including South_T and South_P.
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(lubridate)
  library(stringr)
  library(ggplot2)
  library(mgcv)
  library(sf)
  library(rnaturalearth)
  library(patchwork)
  library(tibble)
  library(purrr)
  library(units)
})

# =========================
# User settings
# =========================

data_dir <- "data"
output_dir <- "figures"

sdt_file     <- file.path(data_dir, "SDT_ALL_master_dryad.csv")
south_t_file <- file.path(data_dir, "South_T_master_dryad.csv")
south_p_file <- file.path(data_dir, "South_P_master_dryad.csv")
meta_file    <- file.path(data_dir, "SDT_SiteMetadata.csv")

analysis_start <- as.Date("2017-08-01")
analysis_end   <- as.Date("2017-11-13")

set.seed(123)  # reproducible systematic 14-day subsample offset

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# =========================
# Helpers
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
    out[bad] <- suppressWarnings(as.POSIXct(x_chr[bad], tz = tz, format = "%Y-%m-%d %H:%M:%OS"))
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

parse_dms_to_dd <- function(x) {
  x <- trimws(as.character(x))
  deg <- suppressWarnings(as.numeric(str_extract(x, "^[0-9]+")))
  min <- suppressWarnings(as.numeric(str_extract(x, "(?<=°\\s*)[0-9.]+")))
  deg + min / 60
}

normalize_site <- function(x) {
  x <- as.character(x)
  x <- str_replace(x, "^SOCAL_T$", "South_T")
  x <- str_replace(x, "^SOCAL_P$", "South_P")
  x
}

find_first_col <- function(df, choices) {
  hit <- choices[choices %in% names(df)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

standardize_metadata <- function(meta_raw) {
  site_col  <- find_first_col(meta_raw, c("Site", "site"))
  depth_col <- find_first_col(meta_raw, c("Depth [m]", "Depth..m.", "Depth_m", "Depth.m.", "Depth"))
  lat_col   <- find_first_col(meta_raw, c("Latitude [N]", "Latitude..N.", "lat", "Lat", "Latitude"))
  lon_col   <- find_first_col(meta_raw, c("Longitude [W]", "Longitude..W.", "lon", "Long", "Longitude"))
  
  needed <- c(site_col, depth_col, lat_col, lon_col)
  if (any(is.na(needed))) {
    stop("Could not find required metadata columns for Site, Depth, Latitude, and Longitude.")
  }
  
  meta_raw %>%
    transmute(
      Site = normalize_site(.data[[site_col]]),
      depth_m = suppressWarnings(as.numeric(.data[[depth_col]])),
      lat = case_when(
        str_detect(as.character(.data[[lat_col]]), "°") ~ parse_dms_to_dd(.data[[lat_col]]),
        TRUE ~ suppressWarnings(as.numeric(.data[[lat_col]]))
      ),
      lon = case_when(
        str_detect(as.character(.data[[lon_col]]), "°") ~ -1 * parse_dms_to_dd(.data[[lon_col]]),
        TRUE ~ suppressWarnings(as.numeric(.data[[lon_col]]))
      )
    ) %>%
    distinct() %>%
    filter(!is.na(Site), !is.na(depth_m), !is.na(lat), !is.na(lon))
}

read_master_file <- function(fp, site_fallback = NULL) {
  df <- read_delim_auto(fp)
  
  ts_col <- find_first_col(df, c("timestamp_utc", "timestamp_dt"))
  if (is.na(ts_col)) stop("No timestamp_utc or timestamp_dt column found in ", basename(fp))
  if (!("UF440_sunset_SPL" %in% names(df))) stop("UF440_sunset_SPL column missing in ", basename(fp))
  
  site_col <- if ("Site" %in% names(df)) "Site" else NA_character_
  
  df %>%
    transmute(
      Site = if (!is.na(site_col)) normalize_site(.data[[site_col]]) else site_fallback,
      timestamp_utc = parse_datetime_flexibly(.data[[ts_col]], tz = "UTC"),
      UF440_sunset_SPL = suppressWarnings(as.numeric(UF440_sunset_SPL))
    ) %>%
    filter(!is.na(timestamp_utc))
}

systematic_subsample_14d <- function(df) {
  df %>%
    arrange(Site, date) %>%
    group_by(Site) %>%
    group_modify(~{
      offset <- sample(0:13, 1)
      start_date <- min(.x$date, na.rm = TRUE) + offset
      .x %>%
        mutate(
          subsample_offset_days = offset,
          keep_14d = (as.integer(date - start_date) %% 14) == 0
        ) %>%
        filter(keep_14d)
    }) %>%
    ungroup()
}

extract_smooth_row <- function(model, var_label) {
  s <- summary(model)
  stab <- as.data.frame(s$s.table)
  stab$term <- rownames(stab)
  rownames(stab) <- NULL
  names(stab)[names(stab) == "p-value"] <- "p_value"
  
  row <- stab %>% filter(str_detect(term, fixed(var_label)))
  if (nrow(row) == 0) {
    return(tibble(edf = NA_real_, p_value = NA_real_))
  }
  
  tibble(
    edf = as.numeric(row$edf[1]),
    p_value = as.numeric(row$p_value[1])
  )
}

make_prediction_df <- function(model, data, focal_var, n = 200) {
  focal_rng <- range(data[[focal_var]], na.rm = TRUE)
  pred_x <- seq(focal_rng[1], focal_rng[2], length.out = n)
  
  newdat <- tibble(
    date_num = median(data$date_num, na.rm = TRUE)
  )
  
  newdat <- newdat[rep(1, n), , drop = FALSE]
  newdat[[focal_var]] <- pred_x
  
  pr <- predict(model, newdata = newdat, type = "response", se.fit = TRUE)
  
  tibble(
    x = pred_x,
    fit = as.numeric(pr$fit),
    lower = as.numeric(pr$fit - 2 * pr$se.fit),
    upper = as.numeric(pr$fit + 2 * pr$se.fit)
  )
}

theme_clean <- theme_bw(base_size = 11) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    legend.position = "none"
  )

# =========================
# Read data
# =========================

sdt_master <- read_master_file(sdt_file)
south_t_master <- read_master_file(south_t_file, site_fallback = "South_T")
south_p_master <- read_master_file(south_p_file, site_fallback = "South_P")

meta_raw <- read_delim_auto(meta_file)
meta <- standardize_metadata(meta_raw)

all_master <- bind_rows(sdt_master, south_t_master, south_p_master) %>%
  filter(
    !is.na(UF440_sunset_SPL),
    !is.na(Site)
  ) %>%
  mutate(
    date = as.Date(timestamp_utc)
  ) %>%
  filter(date >= analysis_start, date <= analysis_end)

# =========================
# Distance from coast
# =========================

usa <- ne_countries(
  country = "United States of America",
  returnclass = "sf"
)

california_bbox <- st_bbox(
  c(xmin = -125, xmax = -114,
    ymin = 32, ymax = 43),
  crs = st_crs(usa)
)

california_coast <- st_crop(usa, california_bbox)
coastline <- st_boundary(california_coast)

sites_sf <- st_as_sf(meta, coords = c("lon", "lat"), crs = 4326)
dist_matrix <- st_distance(sites_sf, coastline)

meta$distance_km <- apply(drop_units(dist_matrix), 1, min) / 1000

# =========================
# Daily maximum sunset SPL
# =========================

daily_max_spl <- all_master %>%
  left_join(
    meta %>% select(Site, depth_m, distance_km),
    by = "Site"
  ) %>%
  group_by(Site, date, depth_m, distance_km) %>%
  summarise(
    max_SPL = max(UF440_sunset_SPL, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    date_num = as.numeric(date)
  ) %>%
  filter(
    is.finite(max_SPL),
    is.finite(depth_m),
    is.finite(distance_km),
    is.finite(date_num)
  ) %>%
  arrange(Site, date)

write_csv(daily_max_spl, file.path(output_dir, "Figure5C_analysis_data.csv"))

# =========================
# 14-day systematic subsample
# =========================

dat_14d <- systematic_subsample_14d(daily_max_spl) %>%
  select(-keep_14d)

write_csv(dat_14d, file.path(output_dir, "Figure5C_14day_subsample.csv"))

# =========================
# Fit GAMs
# =========================

gam_depth <- gam(
  max_SPL ~
    s(depth_m, bs = "ts", k = 4) +
    s(date_num, bs = "ts", k = 4),
  family = gaussian(link = "identity"),
  data = dat_14d,
  method = "REML"
)

gam_distance <- gam(
  max_SPL ~
    s(distance_km, bs = "ts", k = 4) +
    s(date_num, bs = "ts", k = 4),
  family = gaussian(link = "identity"),
  data = dat_14d,
  method = "REML"
)

# =========================
# Predictions for Figure 5C
# =========================

pred_depth <- make_prediction_df(gam_depth, dat_14d, "depth_m")
pred_distance <- make_prediction_df(gam_distance, dat_14d, "distance_km")

write_csv(pred_depth, file.path(output_dir, "Figure5C_depth_model_predictions.csv"))
write_csv(pred_distance, file.path(output_dir, "Figure5C_distance_model_predictions.csv"))

# =========================
# Figure 5C plots
# =========================

p_depth <- ggplot(pred_depth, aes(x = x, y = fit)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line(linewidth = 1.1) +
  geom_rug(
    data = dat_14d,
    aes(x = depth_m),
    inherit.aes = FALSE,
    sides = "b",
    alpha = 0.3
  ) +
  labs(
    x = "Depth (m)",
    y = "Predicted sunset SPL (dB re 1 µPa)",
    title = "A. Depth"
  ) +
  theme_clean

p_distance <- ggplot(pred_distance, aes(x = x, y = fit)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line(linewidth = 1.1) +
  geom_rug(
    data = dat_14d,
    aes(x = distance_km),
    inherit.aes = FALSE,
    sides = "b",
    alpha = 0.3
  ) +
  labs(
    x = "Distance from coast (km)",
    y = "Predicted sunset SPL (dB re 1 µPa)",
    title = "B. Distance from coast"
  ) +
  theme_clean

figure_5c <- p_depth | p_distance

ggsave(
  file.path(output_dir, "Figure5C_depth_distance_gams.png"),
  plot = figure_5c,
  width = 10,
  height = 4.5,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(output_dir, "Figure5C_depth_distance_gams.pdf"),
  plot = figure_5c,
  width = 10,
  height = 4.5
)

# =========================
# Supplemental Table 4-style summary
# =========================

depth_s1 <- extract_smooth_row(gam_depth, "depth_m")
depth_s2 <- extract_smooth_row(gam_depth, "date_num")

dist_s1 <- extract_smooth_row(gam_distance, "distance_km")
dist_s2 <- extract_smooth_row(gam_distance, "date_num")

supp_table4 <- bind_rows(
  tibble(
    Formula = "sunset_SPL ~ s(Depth) + s(Date)",
    n = nrow(dat_14d),
    `Deviance Explained (%)` = round(summary(gam_depth)$dev.expl * 100, 1),
    R2 = round(summary(gam_depth)$r.sq, 2),
    Covariate = c("Depth", "Date"),
    edf = c(depth_s1$edf, depth_s2$edf),
    `P-value` = c(depth_s1$p_value, depth_s2$p_value)
  ),
  tibble(
    Formula = "sunset_SPL ~ s(Distance) + s(Date)",
    n = nrow(dat_14d),
    `Deviance Explained (%)` = round(summary(gam_distance)$dev.expl * 100, 1),
    R2 = round(summary(gam_distance)$r.sq, 2),
    Covariate = c("Distance", "Date"),
    edf = c(dist_s1$edf, dist_s2$edf),
    `P-value` = c(dist_s1$p_value, dist_s2$p_value)
  )
)

write_csv(supp_table4, file.path(output_dir, "Figure5C_supplement_table4.csv"))

# =========================
# Text summary
# =========================

sink(file.path(output_dir, "Figure5C_model_summary.txt"))
cat("=== Figure 5C GAMs ===\n")
cat("Analysis window:", as.character(analysis_start), "to", as.character(analysis_end), "\n")
cat("14-day systematic subsample used to reduce temporal autocorrelation.\n\n")

cat("=== Depth model ===\n")
print(summary(gam_depth))
cat("\n")
print(gam.check(gam_depth))
cat("\n\n")

cat("=== Distance model ===\n")
print(summary(gam_distance))
cat("\n")
print(gam.check(gam_distance))
sink()
