<!-- README.md is generated from README.Rmd. Please edit that file -->

# SpaPheno: Linking Spatial Transcriptomics to Clinical Phenotypes with Interpretable Machine Learning

## Overview

**SpaPheno** is an R package designed to identify, visualize, and
interpret spatial phenotype associations from spatial transcriptomics
and simulated bulk data. Linking spatial transcriptomic patterns to
clinical phenotypes is essential for advancing precision oncology. We
introduce SpaPheno, an interpretable machine learning framework that
bridges spatial transcriptomic data with clinically annotated bulk
RNA-seq datasets. SpaPheno integrates Elastic Net regression with
SHAP-based interpretation to identify spatially localized features—such
as cell types and tissue regions—that are predictive of patient
survival, tumor stage, and immunotherapy response. Through comprehensive
simulations and applications to spatial datasets from primary liver
cancer, kidney renal clear cell carcinoma (KIRC), breast cancer (BRCA),
and melanoma, SpaPheno demonstrates superior performance and broad
applicability. By jointly optimizing predictive power and biological
interpretability, SpaPheno enables the discovery of clinically
meaningful spatial biomarkers, offering a generalizable approach for
spatially informed precision medicine.

<div class="figure" style="text-align: center">

<img src="./man/figures/workflow.jpg" alt="The Overview of SpaPheno" width="80%" height="80%" />
<p class="caption">

The Overview of SpaPheno
</p>

</div>

## :sunny: Key Features

- **Construct predictive models** from cell type compositions and
  phenotypic labels (e.g., disease presence).
- **Map spatial risk distributions** across tissues or organs using
  spatially-aware neighborhood features.
- **Assess statistical significance** through permutation-based tests.
- **Interpret feature contributions** using SHAP (SHapley Additive
  exPlanations) analysis at single-cell resolution.

## :arrow_double_down: Installation

``` r
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

## Install suggested packages
# BiocManager::install(c(
#   "glmnet",
#   "FNN",
#   "survival"
# ))

# install.packages("devtools")
# devtools::install_github("bm2-lab/SpaDo")

# SpaPheno installation
# devtools::install_github("DuanLab1/SpaPheno", dependencies = c("Depends", "Imports", "LinkingTo"))

library(SpaPheno)
library(tidyverse)
library(ggplot2)
library(reshape2)
library(stringr)
library(survival)
```

## 🚀 Quick Start

### Data availability

The data required for the test are all listed in the following google
cloud directory [SpaPheno Demo
Data](https://drive.google.com/drive/folders/1tiSgMjhzvIsirvJwFDIAQIEIhR7qixUW?usp=drive_link).

    ├── BRCAsurvival.RData
    ├── HCC_stage.RData
    ├── HCC_survival.RData
    ├── KIRC_survival.RData
    ├── Melanoma_ICB.RData
    ├── Simulation_osmFISH.RData
    └── Simulation_STARmap.RData

### Simulation osmFISH

This tutorial demonstrates the workflow of **SpaPheno** using simulated
osmFISH data, including:

1.  Loading and visualizing spatial cell annotations.
2.  Defining simulated phenotypes across spatial regions.
3.  Generating pseudo-bulk samples for phenotype modeling.
4.  Building a logistic regression model via automated regularization
    selection.
5.  Performing spatial phenotype risk prediction.
6.  Interpreting key contributing cell types using SHAP values.
7.  Exploring spatial patterns of model residuals for biological
    insight.

Together, this pipeline allows researchers to integrate spatial
structure, cell composition, and predictive modeling to understand how
local cell-type environments contribute to complex phenotypes.

#### load demo data

``` r
rm(list = ls())

load(system.file("extdata", "Simulation_osmFISH.RData", package = "SpaPheno"))

ggplot(test_coordinate, aes(x = X, y = Y, color = sample_information_cellType)) +
  geom_point(size = 1)
ggplot(test_coordinate, aes(x = X, y = Y, color = sample_information_region)) +
  geom_point(size = 1)
```

#### Choose simulated phenotypes

``` r
region_all <- names(table(sample_information_region))
phenotype_simu <- region_all[c(3, 4)]
sample_information_region_choose <- sample_information_region
sample_information_region_choose[!sample_information_region_choose %in% phenotype_simu] <- NA
```

#### Ground truth of simulated phenotypes

``` r
custom_colors <- c("red", "lightgray", "blue")
names(custom_colors) <- c(phenotype_simu[2], "Background", phenotype_simu[1])
sample_information_region_choose[!sample_information_region_choose %in% phenotype_simu] <- "Background"
Ground_truth <- factor(sample_information_region_choose[row.names(test_coordinate)], levels = c(phenotype_simu[1], "Background", phenotype_simu[2]))
ggplot(test_coordinate, aes(x = X, y = Y, color = Ground_truth)) +
  geom_point(size = 1) +
  scale_color_manual(values = custom_colors)
```

#### Simulating bulk data with phenotypes

``` r
pseudo_bulk_simi <- generate_simulated_bulk_data(
  input_data = sample_information_cellType, 
  region_labels = sample_information_region, 
  phenotypes = phenotype_simu, 
  perturbation_percent = 0.1, 
  num_samples = 50, 
  mode = "proportion")

pseudo_bulk_df1 <- pseudo_bulk_simi[[1]]
pseudo_bulk_df2 <- pseudo_bulk_simi[[2]]

bulk_decon <- t(as.matrix(cbind(pseudo_bulk_df1, pseudo_bulk_df2)))
bulk_pheno <- rep(c(0, 1), each = 50)
names(bulk_pheno) <- c(colnames(pseudo_bulk_df1), colnames(pseudo_bulk_df2))
family <- "binomial"
```

#### Obtaining prediction results

``` r
PhenoResult <- SpatialPhenoMap(
  bulk_decon = bulk_decon, 
  bulk_pheno = bulk_pheno, 
  family = family, 
  coord = test_coordinate, 
  resolution = "single_cell", 
  sample_information_cellType = sample_information_cellType, 
  n_perm = 1, 
  p = 0.001)
```

#### SHAP analysis

``` r
pred_result <- PhenoResult[[1]]
phenoPlus <- row.names(pred_result[pred_result$label %in% "phenotype+", ])

model <- PhenoResult[[3]]
X <- as.data.frame(PhenoResult[[2]][phenoPlus, ])
### Not run, this step may take a few mininutes
# shap_test_plus<-compute_shap_spatial(model,as.data.frame(bulk_decon),bulk_pheno,X)
head(shap_test_plus)
```

#### SHAP summary plot

``` r
SpaPheno_SHAP_summary_plot(shap_test_plus, top_n = 31)
```

#### SHAP residual analysis

``` r
resi_result <- SpaPheno_SHAP_residual_analysis(
  shap_df = shap_test_plus,
  feature_name = "Perivascular.Macrophages",
  coordinate_df = test_coordinate, size = 0.8
)
resi_hot <- resi_result$residual_table
head(resi_hot[order(abs(resi_hot$phi_resid_z), decreasing = T), ], 5)
SpaPheno_SHAP_waterfall_plot(shap_test_plus, "cell_5593", top_n = 10)
resi_result$dependence_plot
resi_result$spatial_plot
```

## :book: Vignette

Using the following command and Choosing the `html` for more details.

``` r
utils::browseVignettes(package = "SpaPheno")
```

## :sparkling_heart: Contributing

Welcome any contributions or comments, and you can file them
[here](https://github.com/DuanLab1/SpaPheno/issues).

## :trophy: Acknowledgement

Thanks all the developers of the methods integrated into **SpaPheno**.

## :eight_pointed_black_star: Citation

Kindly cite by using `citation("SpaPheno")` if you think **SpaPheno**
helps you. Alternative way is Bin Duan (2025). *SpaPheno:
Spatially-Informed Phenotype Prediction and Interpretation Using
Single-Cell Data*. R package version 0.0.1,
\<URL:<https://github.com/DuanLab1/SpaPheno/>\>.

## :writing_hand: Authors

- [Bin Duan](mailto:binduan@sjtu.edu.cn)
