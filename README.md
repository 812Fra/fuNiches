
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fuNiches

<img src="man/figures/logo.png" align="right" width="200" />

<!-- badges: start -->

[![R-CMD-check](https://github.com/812Fra/fuNiches/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/812Fra/fuNiches/actions/workflows/R-CMD-check.yaml)
[![License:
MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

## Overview

fuNiches is a flexible and reproducible R package designed for
quantifying and comparing ecological niches between fungi, particularly
ectomycorrhizal (ECM) fungi. Developed to support in-depth ecological
analyses, fuNiches integrates functionality from fundamental R packages
such as hypervolume, vegan, Boruta, and the entire tidyverse suite,
offering a comprehensive workflow from data preparation to detailed
report generation.

## Installation

You can install the development version of fuNiches from GitHub with:

``` r
# If you don't have devtools installed, you can do so with:
# install.packages("devtools")

devtools::install_github("812Fra/fuNiches")
```

## Features

- Data Preparation: Integrates species occurrences (GPKG) and
  environmental rasters (GeoTIFF), handling standardization, filtering,
  and extraction of environmental data.
- Automated Variable Selection: Implements methods like Boruta and
  Random Forest Importance, with built-in filters for multicollinearity
  (Correlation, VIF).
- Niche Overlap Indices: Calculates multiple metrics including
  Schoener’s D, Pianka’s O, Czekanowski’s PSI, and Hurlbert’s L.
- Advanced Niche Comparison: A suite of multivariate analyses including
  PERMANOVA, Ecological Niche Factor Analysis (ENFA), Hypervolume
  construction, and Linear Discriminant Analysis (LDA).
- Statistical Testing: Provides robust tests for niche equivalency,
  similarity validation (using Random Forest cross-validation), and
  associations with categorical variables.
- Network Analysis: Identifies modules of ecologically similar species
  through niche network modularity analysis.
- Comprehensive Outputs: Generates numerical tables (CSV), customizable
  plots (ggplot2), and automated summary reports (Markdown) for each
  analysis.

## Usage

The core workflow of fuNiches follows three main steps:

1.  **Prepare Data**: Load and process occurrences and environmental
    layers.
2.  **Define Settings**: Configure the target species and the analyses
    to run.
3.  **Run Analysis**: Execute the workflow and generate outputs.

### 1. Prepare the niche data object

``` r
library(fuNiches)
prepared_data <- prepare_niche_data(
  occurrences_gpkg_path = "path/to/your/occurrences.gpkg",
  raster_dir_path = "path/to/your/rasters/",
  min_occurrences_per_species = 10
)
```

### 2. Configure the analysis settings

``` r
settings <- set_analysis_settings(
  target_species = "Tuber_aestivum",
  analyses_to_run = c("overlap_indices", "permanova_pairwise", "enfa")
)
```

### 3. Run the analysis workflow

``` r
run_niche_analysis(
  niche_data_object = prepared_data,
  analysis_settings = settings,
  output_dir = "path/to/your/output/"
)
```
