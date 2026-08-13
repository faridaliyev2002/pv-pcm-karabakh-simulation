# Fatty Acid PCM Cooling of Photovoltaic Panels — A Simulation Study for Karabakh, Azerbaijan

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21917430.svg)](https://doi.org/10.5281/zenodo.21917430)

Simulation source code and climate data for the MSc dissertation of **Farid Aliyev** (WMG, University of Warwick, August 2026).

The project models the cooling of a photovoltaic panel by a fatty acid phase change material (PCM) layer fitted to its back, for the climate of Jabrayil, Azerbaijan (39.330°N, 47.052°E). The simulation covers one summer of typical weather data (June to August, 2,208 hourly records) and compares panel temperature and electrical yield with and without the PCM layer, for two fatty acid mixtures: MA-SA 64:36 (melting point 44.13 °C) and LA-SA 75.5:24.5 (melting point 37.0 °C). The heat transfer model is described in Chapter 3 of the dissertation.

## Repository structure

| Path | Contents |
|---|---|
| `Python_Code/PVGIS_TMY_Processing_Script.py` | Cleans and checks the raw PVGIS weather data, extracts the summer months, writes the MATLAB input files, and generates Figures 3.2 and 3.3 of the dissertation. |
| `MATLAB_Codes/Code1_Summer_Simulation_FINAL.m` | Main summer simulation at the final configuration (30 kg PCM per panel). Produces the hourly results, summary statistics, and Figures 4.1 to 4.7. Uses the heat transfer model of Sections 3.6 and 3.7. |
| `MATLAB_Codes/Code2_MeltingPoint_Sensitivity.m` | Tests PCM melting temperatures from 30 °C to 62 °C in 2 °C steps at a fixed PCM mass. Writes `sensitivity_results.csv`. |
| `MATLAB_Codes/Code3_Mass_Sensitivity.m` | Tests PCM masses from 4 kg to 40 kg per panel at two melting temperatures. Writes `mass_sensitivity_results.csv`. |
| `MATLAB_Codes/sensitivity_results.csv` | Output of Code 2, kept here because Code 1 reads it when drawing Figure 4.6. |
| `Data/` | Processed climate inputs (see below). |

## Requirements

The Python script needs Python 3.9 or later with pandas, NumPy, SciPy and Matplotlib. The MATLAB scripts run on MATLAB R2020a or later and use no toolboxes beyond base MATLAB.

## How to run

The MATLAB simulations can be run directly from this repository: copy the three `.m` files and `Data/PVGIS_TMY_Summer_JJA.mat` into one working directory and run `Code1_Summer_Simulation_FINAL`. Code 1 takes one to two minutes on a desktop PC. Run `Code2_MeltingPoint_Sensitivity` (17 cases) and `Code3_Mass_Sensitivity` (14 cases) before Code 1 if you want Figures 4.6 and 4.7 regenerated from fresh results; otherwise the included `sensitivity_results.csv` is used.

The Python script is only needed to rebuild the `Data/` files from scratch. It expects the raw PVGIS file `tmy_39_330_47_052_2005_2023.csv` in its working directory, which is not included here (see the data section below) but can be downloaded in about a minute.

## Data source

The `Data/` folder holds the processed weather data for the site: a Typical Meteorological Year (an artificial year built from the most representative months of the period 2005 to 2023), in full-year and summer-only versions, each as a clean CSV and as a `.mat` file for MATLAB. The data come from the PVGIS-SARAH3 database (European Commission Joint Research Centre), obtained from the [PVGIS web tool](https://re.jrc.ec.europa.eu/pvg_tools/en/) for coordinates 39.330, 47.052. PVGIS data are free to reuse with attribution to the European Union/JRC. Timestamps in the PVGIS files are UTC; dissertation figures use local time (UTC+4).

To rebuild everything from the source, download the weather CSV for the coordinates above from the PVGIS tool and run the Python script on it.

## Relation to the dissertation

The model equations and parameter values are described in Chapter 3 of the dissertation (Sections 3.6 and 3.7 in particular); the scripts here are the implementation used to produce every result in Chapter 4. Each script carries a header comment stating its inputs, outputs and the figures it generates.

## Author and citation

Farid Aliyev, August 2026. If you use this code, please cite:

> Aliyev, F. (2026) *Fatty Acid PCM Cooling of Photovoltaic Panels: A Simulation Study for Karabakh, Azerbaijan* [MATLAB and Python source code]. Zenodo. Available from: https://doi.org/10.5281/zenodo.21917430

Repository: https://github.com/faridaliyev2002/pv-pcm-karabakh-simulation

Code is released under the MIT License (see `LICENSE`).
