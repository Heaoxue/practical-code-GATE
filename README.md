# practical-code-GATE

## Overview
This repository contains the code for the GATE prediction model project, organized into three main components:

1. **Training GATE model** (`1_Model_construction.R`)
   Code for training and constructing the GATE prediction model.

2. **Reproducing performance** figures (`2_performance_figures.R`)
   Code for reproducing the main performance figures presented in the study.

3. **Running the interception model** (`3_interception_model.R`)
   Code for running the interception related analysis.

## File Structure
- `20200728_dwell_group_timing.tsv` / `20200728_dwell_time_groups.tsv`: Raw input dataset
- `1_Model_construction.R`: Script for model construction
- `2_performance_figures.R`: Script to generate model performance plots
- `3_interception_model.R`: Script for interception analysis
- `*.RData`: Saved intermediate results and fitted model objects

## Usage
Run the scripts in numerical order:
```r
source("1_Model_construction.R")
source("2_performance_figures.R")
source("3_interception_model.R")
