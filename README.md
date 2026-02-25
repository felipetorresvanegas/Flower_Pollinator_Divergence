# Pollinator assemblage composition predicts trait divergence in a pollination-generalized plant

**Dataset DOI:** https://doi.org/10.5061/dryad.xwdbrv1tf  

---

## Description of the data and file structure

Data were collected from 27 populations of *Viscaria vulgaris* in southern Sweden between 2021 and 2024.  

At the inflorescence and flower level, we measured traits related to:

- **Pollinator attraction:** inflorescence length, flower number, corolla diameter  
- **Pollen transfer:** tube width, tube length, nectary–stigma distance  

In total, 3,878 inflorescences from 1,797 individuals were measured.

Concurrent repeated 10-minute censuses recorded the number and identity of flower visitors. Visitors were classified into functional groups based on morphology and foraging behaviour.

These data capture spatial, temporal, and developmental variation in both plant traits and pollinator activity, forming the basis for multivariate statistical analyses linking divergence in pollination traits to variation in local pollinator assemblages.

---

## Repository structure

The repository is organized into three main folders:


---

## 1. `code/`

This folder contains R scripts used to analyze the data.

### `Flower_Traits_Variance_Partitioning.Rmd`
- Performs variance partitioning of pollination traits  
- Uses multivariate generalized linear mixed-effects models (GLMMs) in **Hmsc**  
- Includes code to estimate proportional divergence (dP) of pollination traits  

### `Pollinators_Variance_Partitioning.Rmd`
- Performs variance partitioning of pollinator assemblages  
- Uses multivariate GLMMs in **Hmsc**  

### `RRR_Abundance_Population.Rmd`
- Conducts reduced-rank regression (RRR)  
- Relates divergence in pollination traits to the composition of local pollinator assemblages  
- Implemented in **Hmsc**

### `RRR_Abundance_Population_Graphs.Rmd`
- Contains code to reproduce the graphs of the reduced-rank regression (RRR) analysis.
- Quantifies the contribution of pollinator functional groups to the 'pollinator assemblage axis'.
- Includes the graphical outputs summarizing the relationships of the 'pollinator assemblage axis' and pollination traits as obtained from the RRR analysis.

### `RRR_Predictive_Explanatory_Power.Rmd`
- Evaluates predictive and explanatory performance of the RRR models  

Running these scripts reproduces the main statistical analyses described in the associated manuscript.

---

## 2. `data/`

This folder contains raw and processed datasets used in the main analyses.

### `Data_Pollination_Traits_Master.rds`
- Master file containing measurements of pollination traits  

### `Data_Pollinators_Master.rds`
- Master file containing pollinator census data  

### `RRR_Model_Abundance.rds`
- Output from reduced-rank regression analyses  
- Links pollinator assemblage composition to pollination trait divergence  

### `RRR_Pollinator_Assemblage_Axis.rds`
- Contains the two main reduced axes derived from the RRR analysis  
- Stores population-level scores for:
  - Pollinator assemblage axis 1  
  - Pollinator assemblage axis 2  
- These axes represent gradients of geographic variation in pollinator assemblage composition that best explain divergence in pollination traits  
- Values correspond to composite variables constructed as linear combinations of pollinator functional groups during the RRR analysis  

---

## 3. `R_Functions/`

This folder contains custom R functions required for specific components of the RRR analyses.

### `constructGradientRRR.R`
Defines the R function used to construct predictor gradients along the reduced-rank regression (RRR) axes.

---

## Missing values

Missing or unavailable measurements are indicated as `NA`.  

Blank cells are not used to represent missing values to ensure compatibility with R functions.

---

## Code and software

All data can be viewed and analyzed using the free and open-source software **R**.

Analyses and data handling rely primarily on the following R packages:

- **Data import and manipulation:** `readr`, `dplyr`, `tibble`, `tidyr`, `reshape2`, `stringr`  
- **Visualization:** `ggplot2`, `ggpubr`  
- **Statistical modeling:** `Hmsc`  
- **General data organization**

These packages enable full reproduction of the data cleaning, visualization, and multivariate analyses conducted in the study.
