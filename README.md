# UF440 climate-linked redistribution manuscript code

This repository contains code used for data processing, statistical analyses, and figure generation for the manuscript:

**Kim et al. (2026).** Climate-linked redistribution of a midwater fish chorus in the Southern California Bight.

## Repository contents

### Optional preprocessing
- `00_build_climate_timeseries.R`  
  Builds the monthly climate/environment table used in downstream analyses and figures.

### Main figure scripts
- `01_plot_figure1a_study_region_map.m`  
  Generates the main study-region map for Figure 1A.

- `02_plot_figure3_monthly_focal_timeseries.R`  
  Generates monthly UF440 chorus time series for the focal multi-year sites in Figure 3.

- `03_plot_figure4a_annual_presence_spl_ssta_panels.R`  
  Generates annual presence, SPL, and SST anomaly panels for Figure 4A and related supplemental panels.

- `04_fit_plot_figure4b_northB_seasonal_presence_gam.R`  
  Fits and plots the North_B seasonal chorus presence GAM used in Figure 4B.

- `05_plot_figure4c_monthly_fish_environment_timeseries.R`  
  Generates monthly fish/environment time series panels for Figure 4C.

- `06_plot_figure5a_sdt_array_bubble_map.m`  
  Generates the SDT array bubble map used for Figure 5A.

- `07_plot_figure5b_sdt_diel_selected_sites.R`  
  Generates diel plots for selected SDT sites used in Figure 5B.

- `08_fit_plot_figure6c_southT_gamm.R`  
  Fits and plots the South_T GAMM used for Figure 6C.

- `09_plot_figure7_vertical_array_psd443.R`  
  Generates the vertical array PSD plots used for Figure 7.

### Supplemental figure scripts
- `S02_plot_supplemental_figure2_diel_focal_sites.R`  
  Generates the focal-site diel plots used in Supplemental Figure 2.

## Data availability

Analysis-ready acoustic and environmental data are archived in the Dryad Digital Repository:  
(see Kim et al., 2026 for link)

Code associated with this manuscript is archived in this GitHub repository

## Software

Most analyses and figures were produced in **R**. Some map figures were produced in **MATLAB**.

### R packages used
Common packages used across scripts include:
- dplyr
- tidyr
- readr
- lubridate
- ggplot2
- patchwork
- mgcv
- nlme
- purrr
- scales
- suncalc
- data.table

Some scripts may require additional packages noted in the script headers.

## Expected input files

Scripts assume analysis-ready input files are stored in a local `data/` directory. These include:
- site master matrices
- climate time series tables
- nightly depth/backscatter summaries
- vertical array PSD files
- map input files for MATLAB figures

See the script headers for the expected filenames and required columns.

## Notes

- Scripts were cleaned from working analysis code used during manuscript preparation.
- Some figure layout refinements and panel assembly steps were completed outside these scripts in illustrator.
- `00_build_climate_timeseries.R` is included for provenance but may not be necessary to rerun the main figure scripts if the analysis-ready climate table is already available from Dryad.

## Contact

For questions about the code or data, contact:  
ebkim@ucsd.edu
