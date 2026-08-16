# TriSmart

`TriSmart` is an R package for hybrid missing-value imputation in quantitative
proteomics data. Input matrices use samples in rows and proteins/features in
columns. The workflow assigns proteins to three missingness zones and then
uses a method tailored to each zone.

## Installation

From the package directory:

```r
install.packages(".", repos = NULL, type = "source")
```

The zone model always uses `mgcv` and `mclust`. The selected imputation
methods determine which additional packages are needed:

- `KNN`: `impute`
- `RandomForest`: `missForest`
- `BPCA`: `pcaMethods`
- mixture-zone modeling: `glmnet` when enough observations are available
- `QRILC`: `imputeLCMD`
- `GSimp`: a package or sourced implementation that provides `GS_impute()`
- `MinProb` and `MinDet`: no extra package

Several of these packages are distributed through Bioconductor. For example:

```r
install.packages(c("mclust", "glmnet", "missForest", "imputeLCMD"))
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("impute", "pcaMethods"))
```

## Worked example using the included data

The package includes `exampledata`, the supplied example matrix with 20
samples and 6,611 proteins. This illustration selects evenly spaced proteins
to keep the run relatively quick while retaining a range of missingness
patterns.

```r
library(TriSmart)
data(exampledata)

dim(exampledata)
mean(is.na(exampledata))

# Select 120 proteins spread across the example matrix.
keep <- unique(round(seq(1, ncol(exampledata), length.out = 120)))
example_subset <- exampledata[, keep, drop = FALSE]

# Inspect the zone assignments before imputation.
zones <- classify_missingness_zones(example_subset, seed = 2026)
table(zones$protein_zone_table$final_zone)

# KNN + MinDet is a useful, reproducible introductory configuration.
fit <- TriSmart(
  X_miss = example_subset,
  firstImpute = "KNN",
  secondImpute = "MinDet",
  workers = 1,
  seed = 2026,
  verbose = TRUE,
  return_details = TRUE
)

imputed <- fit$imputed_data
anyNA(imputed)
table(fit$zone_result$protein_zone_table$final_zone)
fit$mixture_diagnostics

# TriSmart never changes values that were originally observed.
observed <- !is.na(example_subset)
stopifnot(identical(imputed[observed], example_subset[observed]))
```

For a full analysis, pass the entire matrix after confirming that every
protein and every sample has at least one observed value:

```r
full_fit <- TriSmart(
  exampledata,
  firstImpute = "KNN",
  secondImpute = "MinDet",
  seed = 2026,
  return_details = TRUE
)
```

## Input requirements

- rows are samples and columns are proteins/features;
- all values are numeric;
- every feature has a unique, non-empty column name;
- no feature or sample is completely missing; and
- infinite values are not allowed.

The lower-level exported functions are available for researchers who want to
inspect or customize individual stages: `classify_missingness_zones()`,
`impute_high_abundance()`, `fast_mixture_impute()`, and
`impute_mnar_by_col()`.
