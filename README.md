# Longitudinal monitoring exposes correlated temporal protein variations in the female plasma proteome

This repository contains the R analysis scripts used for the study:

**Longitudinal monitoring exposes correlated temporal protein variations in the female plasma proteome**

The scripts process DIA-NN plasma proteomics data and perform the analyses used to investigate **sex-specific and longitudinal changes in plasma protein concentrations**.

## Overview

The workflow:

* Processes and quality-controls DIA-NN proteomics output.
* Generates protein-level LFQ measurements and estimated protein concentrations.
* Compares protein concentrations between female and male participants.
* Examines longitudinal protein variation across repeated sampling time points.
* Investigates correlations between selected proteins.
* Generates the figures used to visualize sex differences and longitudinal protein profiles.

The analysis focuses in particular on **PZP, SHBG, FETUB, AGT, SERPINA6, SERPINA7, CP, APOL1, and KNG1**, which showed notable temporal variation in the study.

## Data

The analysis uses data from the **TIMES** cohort, with additional analyses using the **AICOVI** and **AMSBIO** datasets.


## Requirements

The analysis is written in **R** and requires packages for data processing, DIA-NN analysis, statistical testing, and visualization, including `tidyverse`, `arrow`, `diann`, `effsize`, `FSA`, `ggh4x`, and `ggstatsplot`.

