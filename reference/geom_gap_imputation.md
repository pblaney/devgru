# Diagnostic plot used for reviewing the imputation of gaps in coverage profile

Creates a simple ggpplot2 style plot for quickly reviewing how the
imputation of each gap identified in a coverage file

Note, this function was designed to be run as part of the
[`get_cov_gap_imputation()`](https://pblaney.github.io/devgru/reference/get_cov_gap_imputation.md)
workflow.

## Usage

``` r
geom_gap_imputation(original_gap, imputed_gap_gr)
```

## Arguments

- original_gap:

  GenomicRanges object of gap

- imputed_gap_gr:

  GenomicRanges object of imputed values across gap

## Value

ggplot object

## Examples

``` r
# Snippet from within the function `` that calls this geom
# geom_gap_imputation_diagnostic(
# original_gap = gap_of_interest,
# imputed_gap_gr = partition_regression_data)
```
