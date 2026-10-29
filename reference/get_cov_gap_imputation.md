# Generate gap imputed coverage from fragCounter output

Given fragCounter coverage with small, zero signal gaps, run
neighbrohood-based partition regression gap imputation in a parallel
manner to produce a gap imputed coverage file that is provides cleaner
foreground signal in JaBbA and GISTIC.

## Usage

``` r
get_cov_gap_imputation(
  path_to_fragcounter_cov,
  cpus = 1,
  make_diagnostic_plots = FALSE,
  path_to_diagnostic_plots_dir = NULL,
  exp_colnames = c("reads", "gc", "map", "reads.corrected"),
  verbose = T
)
```

## Arguments

- path_to_fragcounter_cov:

  Path to fragCounter coverage RDS

- cpus:

  Number of CPUs to use for parallel execution, default: 1

- make_diagnostic_plots:

  Generate diagnostic plots for each imputed gap, default: FALSE

- path_to_diagnostic_plots_dir:

  Path to directory to place diagnostic plot PNGs

- exp_colnames:

  Expected column names for coverage input GR

- verbose:

  Output CLI during workflow, default: TRUE

## Value

GenomicRanges object with gap imputed coverage
