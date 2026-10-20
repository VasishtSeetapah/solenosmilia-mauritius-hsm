# Habitat suitability model for *Solenosmilia variabilis* in the Mauritius region

R code for Seetapah V. and Jogee S.Y. "Predicting deep-sea coral habitat in a data-poor region: a globally calibrated habitat model for *Solenosmilia variabilis* in the Mauritius region", submitted to Deep-Sea Research Part II.

The model is a Maximum Entropy (MaxEnt) model fitted with the R package maxnet. It is calibrated on global occurrence records from the Ocean Biodiversity Information System (OBIS) and 10,000 random background points, tuned with ENMeval on spatial block folds (200 km blocks, five folds), and projected over the exclusive economic zone of the Republic of Mauritius, which includes the Chagos Archipelago and Tromelin, and over the Mauritius and Seychelles Joint Management Area.

## Contents

| File | What it holds |
| --- | --- |
| `R/solenosmilia_hsm.R` | The full workflow in twelve sections, from occurrence records to figures |
| `CITATION.cff` | How to cite this code |
| `LICENSE` | MIT licence |

## Workflow

| Section | Step | Paper |
| --- | --- | --- |
| 0 | Settings, package check and shared definitions (run first in every session) | |
| 1 | Presences, 10,000 random background points and the model table | Section 2.2 |
| 2 | Spatial block folds | Section 2.5 |
| 3 | Predictor layers cropped to the projection window, and a 5 km display stack | Section 2.3 |
| 4 | Tuning of 35 candidate models with ENMeval and selection by minimum AICc | Section 2.4, Table S2 |
| 5 | Spatial-fold evaluation, permutation importance, response curves, bivariate surfaces and fold uncertainty | Sections 3.1 to 3.4 and 3.6, Table 2, Figs 4 to 6 and 8 |
| 6 | Extrapolation diagnostics (MESS and ExDet) | Section 3.6, Fig. 9 |
| 7 | Projection at full resolution and the area of each suitability class | Section 3.5, Table 3 |
| 8 | Sensitivity analyses (thinning, Indian Ocean hold-out, target-group background) | Tables S3 and S4 |
| 9 | Predictor correlations and the thresholds for the suitability classes | Section 2.6, Table S1 |
| 10 | Occurrence summaries and nearest-neighbour distances | Section 2.2, Figs S1 and S2 |
| 11 | Predictor maps | Fig. 3 |
| 12 | Maps with the maritime boundaries and place names | Figs 1 and 7 |

## Getting started

1. Install R (version 4.4.0 was used) and the packages.

   ```r
   install.packages(c("terra", "sf", "maxnet", "ENMeval", "blockCV", "pROC", "ecospat", "predicts",
                      "ggplot2", "patchwork", "viridisLite", "ragg"))
   ```

2. Put the input files in a folder called `inputs` next to the `R` folder, laid out as in the table below. Another location can be set with `INPUT_DIR` at the top of the script.
3. Open `R/solenosmilia_hsm.R`, run Section 0, then run the other sections in order. Results go to `outputs` (tables, figures and rasters), and the R and package versions of the run are written to `outputs/session_info.txt`.

Section 0 stops with the exact `install.packages()` command if a package is missing. Sections 4, 7 and 8.2 take the longest. Section 7 predicts over about 24 million grid cells, so a machine with at least 16 GB of memory is advisable.

## Input data

The input data are not distributed with this repository. They can be obtained from their providers and prepared as described in Sections 2.2 and 2.3 of the paper.

| File in `inputs` | Content | Source |
| --- | --- | --- |
| `bathymetry.tif` | Depth (m, negative below sea level) | GEBCO 2023 grid (https://www.gebco.net) |
| `slope_degrees.tif` | Slope (degrees), from the bathymetry with `terra::terrain` | Derived |
| `aspect_sin.tif`, `aspect_cos.tif` | Sine and cosine of aspect, from the bathymetry with `terra::terrain` | Derived |
| `temperature.tif` | Near-bottom temperature, 2000 to 2020 mean (°C) | Bio-ORACLE v3 (https://www.bio-oracle.org) |
| `kinetic_energy.tif` | KE = ½ρU² in N m⁻², with U the near-bottom current speed (m s⁻¹) and ρ = 1,025 kg m⁻³ | Derived from Bio-ORACLE v3 |
| `geomorphic_class.tif` | Integer class codes of the Seafloor Geomorphic Features Map | Harris et al. (2014) |
| `Solenosmilia_variabilis_OBIS.csv` | OBIS download (8,469 records), with columns `decimalLongitude`, `decimalLatitude`, `depth` and `year` | OBIS (https://obis.org) |
| `Solenosmilia_variabilis_sample400.csv` | The 400-record training sample, with the same columns | Random sample of the downloaded records, stratified by ocean basin |
| `OBIS_other_scleractinians.csv` | Records of four other scleractinian species, with columns `species`, `lon` and `lat` (target-group test only) | OBIS |
| `boundaries/EEZ_MU.shp`, `Chagos.shp`, `Tromelin.shp`, `JMCA.shp` | Maritime zones of the Republic of Mauritius and the Joint Management Area | See the paper |

The seven rasters must share one grid in WGS 84 (EPSG:4326). The grid used in the paper has cells of about 0.0051° × 0.0041° (about 560 m × 450 m).

## Reproducibility

Random seeds are fixed at 42. The background sample, the spatial folds and the permutations depend on random draws, and the way these draws are made can change between versions of R and of the packages. Runs with other versions therefore give values close to, but not always identical with, those in the paper. The versions used for the paper are listed below, and older package versions can be installed with `remotes::install_version()`.

## Software

R 4.4.0 with terra 1.7.78, sf 1.0.20, maxnet 0.1.4, ENMeval 2.0.5.2, blockCV 3.2.0, pROC 1.18.5, ecospat 4.1.2, predicts 0.1.19, ggplot2 4.0.1 and patchwork.

## Licence and citation

The code is released under the MIT licence. Please cite the paper, and this repository as given in `CITATION.cff`.
