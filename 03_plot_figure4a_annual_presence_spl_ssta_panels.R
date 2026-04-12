# -----------------------------------------------------------------------------
# 03_plot_figure4a_annual_presence_spl_ssta_panels
#
# Purpose:
# Generate annual UF440 chorus panels for North_B and South_P, showing:
#   1) total chorus hours per year (bars),
#   2) annual sunset SPL distributions (violins) with annual median points,
#   3) annual recording effort (gray points, secondary axis),
#   4) annual mean sea surface temperature anomaly (SSTA) tiles.
#
# In the paper:
# This script reproduces the core annual plotting workflow used for:
#   - Figure 4A (North_B annual chorus presence + SSTA context), and
#   - Supplemental Figure 4A-B / Supplemental Figure 5 annual panels,
#     excluding the later Illustrator-added BTA layer.
#
# Notes:
# - This script uses RAW annual chorus hours, not effort-adjusted values.
# - Recording effort is shown explicitly as gray points on a secondary axis.
# - Light gray shading is used only for years without recording effort.
# - SSTA is read from the analysis-ready Climate_Timeseries.csv table.
# - Site names use final manuscript names: North_B and South_P.
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

north_file   <- file.path(data_dir, "North_B_master_dryad.csv")
southp_file  <- file.path(data_dir, "South_P_master_dryad.csv")
climate_file <- file.path(data_dir, "Climate_Timeseries.csv")

years_all <- 2008:2026
bin_minutes <- 20

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

safe_mean <- function(x) {
  if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

# =========================
# Read annual SSTA
# =========================

make_yearly_sst <- function(climate_df, site_label) {
  climate_df %>%
    filter(Site == site_label) %>%
    mutate(
      Date = as.Date(Date),
      Year = year(Date),
      SST_anom_detrended = as.numeric(SST_anom_detrended)
    ) %>%
    group_by(Year) %>%
    summarise(
      sst_anom_year_mean = safe_mean(SST_anom_detrended),
      .groups = "drop"
    ) %>%
    complete(Year = years_all) %>%
    mutate(site = site_label)
}

climate_df <- read_delim_auto(climate_file)

sst_year <- bind_rows(
  make_yearly_sst(climate_df, "North_B"),
  make_yearly_sst(climate_df, "South_P")
)

# =========================
# Read acoustic master tables
# =========================

read_master <- function(fp, site_label) {
  df <- read_delim_auto(fp)
  
  required_cols <- c("timestamp_dt", "UF440_roving_PSD", "UF440_PA_no20iso_noday", "UF440_sunset_SPL")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns in ", basename(fp), ": ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  df %>%
    mutate(
      Site = site_label,
      timestamp_dt = parse_datetime_flexibly(timestamp_dt, tz = "UTC"),
      UF440_roving_PSD = as.numeric(UF440_roving_PSD),
      UF440_PA_no20iso_noday = as.numeric(UF440_PA_no20iso_noday),
      UF440_sunset_SPL = as.numeric(UF440_sunset_SPL)
    ) %>%
    mutate(across(where(is.numeric), ~ ifelse(is.nan(.) | is.infinite(.), NA_real_, .))) %>%
    filter(!is.na(timestamp_dt))
}

north_master <- read_master(north_file, "North_B")
south_master <- read_master(southp_file, "South_P")

# =========================
# Build annual fish summaries
# =========================

make_yearly_fish <- function(master_df, site_label, bin_minutes = 20) {
  master_df %>%
    mutate(
      Year = year(as.Date(timestamp_dt)),
      on_effort = !is.na(UF440_roving_PSD),
      pa_eff = if_else(on_effort, coalesce(UF440_PA_no20iso_noday, 0), 0),
      spl_eff = if_else(on_effort, UF440_sunset_SPL, NA_real_)
    ) %>%
    group_by(Year) %>%
    summarise(
      effort_bins = sum(on_effort, na.rm = TRUE),
      effort_hours = effort_bins * (bin_minutes / 60),
      presence_hours_raw = sum(pa_eff == 1, na.rm = TRUE) * (bin_minutes / 60),
      spl_mean_presenceonly = {
        y_pres <- spl_eff[pa_eff == 1]
        if (length(y_pres) == 0 || all(is.na(y_pres))) NA_real_ else mean(y_pres, na.rm = TRUE)
      },
      .groups = "drop"
    ) %>%
    mutate(
      hours_in_year = if_else(leap_year(Year), 366, 365) * 24,
      effort_frac = pmin(pmax(effort_hours / hours_in_year, 0), 1)
    ) %>%
    complete(Year = years_all) %>%
    mutate(site = site_label) %>%
    arrange(Year)
}

fish_year <- bind_rows(
  make_yearly_fish(north_master, "North_B", bin_minutes = bin_minutes),
  make_yearly_fish(south_master, "South_P", bin_minutes = bin_minutes)
)

# =========================
# SPL violin data (presence-only bins)
# =========================

make_spl_violin_df <- function(master_df, site_label) {
  master_df %>%
    mutate(
      Year = year(as.Date(timestamp_dt)),
      on_effort = !is.na(UF440_roving_PSD),
      pa = as.numeric(UF440_PA_no20iso_noday),
      spl = as.numeric(UF440_sunset_SPL)
    ) %>%
    filter(
      Year %in% years_all,
      on_effort,
      pa == 1,
      !is.na(spl),
      is.finite(spl)
    ) %>%
    transmute(site = site_label, Year, spl)
}

spl_violin_all <- bind_rows(
  make_spl_violin_df(north_master, "North_B"),
  make_spl_violin_df(south_master, "South_P")
)

# =========================
# Combine annual fish + annual SSTA
# =========================

year_all <- fish_year %>%
  left_join(sst_year, by = c("site", "Year"))

# Shared SSTA scale
sst_lim <- max(abs(year_all$sst_anom_year_mean), na.rm = TRUE)
if (!is.finite(sst_lim) || sst_lim == 0) sst_lim <- 1

sst_cols <- c("#08306B", "#DBEAFE", "#FEE2E2", "#7F0000")
sst_values <- scales::rescale(c(-sst_lim, -0.0001, 0.0001, sst_lim), to = c(0, 1))

# =========================
# Plot helper
# =========================

plot_year_panel <- function(df_site,
                            yvar,
                            ylab,
                            title_lab,
                            as_bar = TRUE,
                            as_violin = FALSE,
                            violin_df = NULL,
                            violin_min_n = 10,
                            sst_lim = 1,
                            sst_cols = c("#08306B", "#DBEAFE", "#FEE2E2", "#7F0000"),
                            sst_values = c(0, 0.49, 0.51, 1)) {
  
  df_site <- df_site %>%
    mutate(
      Year = as.integer(Year),
      no_effort = is.na(effort_frac) | effort_frac <= 0,
      effort_frac_dot = coalesce(effort_frac, 0),
      effort_frac_dot = pmin(pmax(effort_frac_dot, 0), 1),
      y_plot = .data[[yvar]]
    )
  
  noeff_rect <- df_site %>%
    filter(no_effort) %>%
    transmute(xmin = Year - 0.5, xmax = Year + 0.5)
  
  violin_plot <- NULL
  violin_meds <- NULL
  
  if (as_violin) {
    if (is.null(violin_df)) stop("as_violin = TRUE but violin_df is NULL")
    
    violin_plot <- violin_df %>%
      left_join(df_site %>% select(Year, effort_frac), by = "Year") %>%
      filter(coalesce(effort_frac, 0) > 0) %>%
      group_by(Year) %>%
      filter(n() >= violin_min_n) %>%
      ungroup()
    
    violin_meds <- violin_plot %>%
      group_by(Year) %>%
      summarise(med = median(spl, na.rm = TRUE), .groups = "drop")
    
    yvals <- violin_plot$spl
  } else {
    yvals <- df_site$y_plot
  }
  
  y_finite <- yvals[is.finite(yvals)]
  if (length(y_finite) == 0) {
    y_min <- 0
    y_max <- 1
  } else {
    y_min <- min(y_finite)
    y_max <- max(y_finite)
  }
  
  if (as_bar) {
    y_lower <- 0
    y_upper_data <- ifelse(is.finite(y_max) && y_max > 0, y_max, 1)
    pad <- 0.06 * y_upper_data
    if (!is.finite(pad) || pad == 0) pad <- 0.5
    y_upper <- y_upper_data + pad
  } else {
    rng <- y_max - y_min
    if (!is.finite(rng) || rng == 0) rng <- max(1, abs(y_max))
    pad <- 0.08 * rng
    if (!is.finite(pad) || pad == 0) pad <- 0.5
    y_lower <- y_min - pad
    y_upper <- y_max + pad
  }
  
  stripe_height <- max(0.12 * (y_upper - y_lower), 0.5)
  
  stripe_df <- df_site %>%
    transmute(
      Year,
      xmin = Year - 0.5,
      xmax = Year + 0.5,
      ymin = y_upper + 0.02 * stripe_height,
      ymax = y_upper + 1.02 * stripe_height,
      sst_anom_year_mean = sst_anom_year_mean
    )
  
  y_upper_plot <- max(stripe_df$ymax, na.rm = TRUE)
  
  stripe_bottom <- suppressWarnings(min(stripe_df$ymin, na.rm = TRUE))
  if (!is.finite(stripe_bottom)) stripe_bottom <- y_upper
  
  y_eff_base <- y_lower
  y_eff_span <- stripe_bottom - y_lower
  if (!is.finite(y_eff_span) || y_eff_span <= 0) {
    y_eff_span <- max(1, y_upper - y_lower)
  }
  
  df_site <- df_site %>%
    mutate(eff_y = y_eff_base + effort_frac_dot * y_eff_span)
  
  eff_breaks_sec <- c(0, 25, 50, 100)
  eff_labels_sec <- c("0", "25", "50", "100")
  
  p <- ggplot(df_site, aes(x = Year)) +
    geom_rect(
      data = noeff_rect,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE,
      fill = "grey92",
      alpha = 1
    ) +
    geom_rect(
      data = stripe_df,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = sst_anom_year_mean),
      inherit.aes = FALSE,
      alpha = 1
    ) +
    {
      if (as_violin) {
        list(
          geom_violin(
            data = violin_plot,
            aes(x = Year, y = spl, group = Year),
            inherit.aes = FALSE,
            width = 0.82,
            trim = TRUE,
            scale = "width",
            fill = "grey70",
            color = "grey25",
            linewidth = 0.25,
            na.rm = TRUE
          ),
          geom_point(
            data = violin_meds,
            aes(x = Year, y = med),
            inherit.aes = FALSE,
            size = 2.4,
            color = "black",
            na.rm = TRUE
          )
        )
      } else if (as_bar) {
        geom_col(aes(y = y_plot), width = 0.82, fill = "black", na.rm = TRUE)
      } else {
        list(
          geom_line(aes(y = y_plot), linewidth = 0.5, color = "black", na.rm = TRUE),
          geom_point(aes(y = y_plot), size = 2.6, color = "black", na.rm = TRUE)
        )
      }
    } +
    geom_point(aes(y = eff_y), color = "grey45", size = 1.6, na.rm = TRUE) +
    scale_x_continuous(
      breaks = years_all,
      limits = c(min(years_all) - 0.5, max(years_all) + 0.5),
      expand = expansion(mult = c(0, 0))
    ) +
    scale_y_continuous(
      name = ylab,
      limits = c(y_lower, y_upper_plot),
      expand = expansion(mult = c(0.02, 0)),
      sec.axis = sec_axis(
        trans = ~ (. - y_eff_base) / y_eff_span * 100,
        breaks = eff_breaks_sec,
        labels = eff_labels_sec,
        name = "Effort (%)"
      )
    ) +
    scale_fill_gradientn(
      colors = sst_cols,
      values = sst_values,
      limits = c(-sst_lim, sst_lim),
      oob = squish,
      na.value = "grey85",
      name = "SST anom (°C)"
    ) +
    guides(fill = guide_colorbar(
      barheight = unit(100, "pt"),
      barwidth = unit(10, "pt"),
      title.position = "right",
      title.theme = element_text(angle = 270, vjust = 0.5, hjust = 0.5),
      label.position = "right"
    )) +
    labs(title = title_lab, x = NULL) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      plot.title = element_text(face = "bold", hjust = 0),
      axis.title.y.right = element_text(margin = margin(l = 6)),
      axis.text.y.right = element_text(margin = margin(l = 3)),
      plot.margin = margin(t = 4, r = 6, b = 4, l = 4, unit = "pt"),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0),
      legend.box.spacing = unit(2, "pt")
    )
  
  p
}

# =========================
# Build final per-site plots
# =========================

north <- year_all %>% filter(site == "North_B")
south <- year_all %>% filter(site == "South_P")

vio_north <- spl_violin_all %>% filter(site == "North_B")
vio_south <- spl_violin_all %>% filter(site == "South_P")

p_hr_north <- plot_year_panel(
  north,
  yvar = "presence_hours_raw",
  ylab = "Ttl hrs/yr",
  title_lab = "North_B",
  as_bar = TRUE,
  as_violin = FALSE,
  sst_lim = sst_lim,
  sst_cols = sst_cols,
  sst_values = sst_values
)

p_spl_north <- plot_year_panel(
  north,
  yvar = "spl_mean_presenceonly",
  ylab = "Sunset SPL",
  title_lab = "North_B",
  as_bar = FALSE,
  as_violin = TRUE,
  violin_df = vio_north,
  violin_min_n = 10,
  sst_lim = sst_lim,
  sst_cols = sst_cols,
  sst_values = sst_values
)

p_hr_south <- plot_year_panel(
  south,
  yvar = "presence_hours_raw",
  ylab = "Ttl hrs/yr",
  title_lab = "South_P",
  as_bar = TRUE,
  as_violin = FALSE,
  sst_lim = sst_lim,
  sst_cols = sst_cols,
  sst_values = sst_values
)

p_spl_south <- plot_year_panel(
  south,
  yvar = "spl_mean_presenceonly",
  ylab = "Sunset SPL",
  title_lab = "South_P",
  as_bar = FALSE,
  as_violin = TRUE,
  violin_df = vio_south,
  violin_min_n = 10,
  sst_lim = sst_lim,
  sst_cols = sst_cols,
  sst_values = sst_values
)

plot_north <- (p_hr_north | p_spl_north) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

plot_south <- (p_hr_south | p_spl_south) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

print(plot_north)
print(plot_south)

# =========================
# Export
# =========================

w_in <- 7.5 * 1.15
h_in <- 3
dpi <- 300

ggsave(
  file.path(output_dir, "annual_North_B_presence_SPL_SSTA.svg"),
  plot = plot_north,
  width = w_in, height = h_in, units = "in", device = "svg"
)

ggsave(
  file.path(output_dir, "annual_South_P_presence_SPL_SSTA.svg"),
  plot = plot_south,
  width = w_in, height = h_in, units = "in", device = "svg"
)

ggsave(
  file.path(output_dir, "annual_North_B_presence_SPL_SSTA.png"),
  plot = plot_north,
  width = w_in, height = h_in, units = "in",
  dpi = dpi, device = ragg::agg_png
)

ggsave(
  file.path(output_dir, "annual_South_P_presence_SPL_SSTA.png"),
  plot = plot_south,
  width = w_in, height = h_in, units = "in",
  dpi = dpi, device = ragg::agg_png
)