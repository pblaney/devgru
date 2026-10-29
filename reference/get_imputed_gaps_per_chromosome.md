# Single chromosome run imputation of small, NA read gaps in fragCounter coverage

Run imputation of small, NA read gaps in fragCounter using partition
regression of signal in neighborhood adjacent to the gap. The gaps are
typically small, less than 500 bp on average and inter spaced among
large segments of continuous, homogeneous signal. The gap signal is far
less than expected for real deletions, even homozygous events Imputation
greatly improves the foreground signal in downstream tools like JaBbA
and GISTIC.

Note, this function was designed to be run as part of the
[`get_cov_gap_imputation()`](https://pblaney.github.io/devgru/reference/get_cov_gap_imputation.md)
workflow. Also, see
[`chrompar()`](https://pblaney.github.io/devgru/reference/chrompar.md)
for executing in parallel.

## Usage

``` r
get_imputed_gaps_per_chromosome(
  chrom_for_imputation,
  gapless_coverage_whitelist,
  coverage_gaps,
  sample_id,
  make_plots = F,
  path_for_plots = NULL
)
```

## Arguments

- chrom_for_imputation:

  The chromosome to run gap imputation algorithm on

- gapless_coverage_whitelist:

  GenomicRanges object of coverage after removing blacklist regions and
  gaps

- coverage_gaps:

  GenomicRanges object of gaps from the original coverage

- sample_id:

  Name of sample

- make_plots:

  Generate diagnostic plots of imputation results? default: FALSE

- path_for_plots:

  Path to directory to hold diagnostic plot small PNGs (~80 KB per plot)

## Value

data.table with gap imputed coverage
