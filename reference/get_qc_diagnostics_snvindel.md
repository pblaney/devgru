# Generate diagnostic plots for SNV & InDel variant calling QC checks

Single command to generate a comprehensive set of diagnostic plots that
assist in performing QC checks of SNVs & InDels called from MGP1000
using union consensus. Can be used as both a first pass of unfiltered
variants to see consensus skew as well as used after filtering to show
changes in calling metrics.

## Usage

``` r
get_qc_diagnostics_snvindel(
  path_to_snv_dir = NULL,
  path_to_indel_dir = NULL,
  snv_indel_obj = NULL,
  plot_sample_names = T,
  include_caveman = F
)
```

## Arguments

- path_to_snv_dir:

  Path to directory of MGP1000 union consensus SNVs

- path_to_indel_dir:

  Path to directory of MGP1000 union consensus InDels

- snv_indel_obj:

  Data.table object with merged SNVs & InDels, used in place of paths to
  directories

- plot_sample_names:

  Output plots to include sample names on y-axis, default: TRUE

- include_caveman:

  Add caveman to consensus list if used in SNV variant calling, default:
  FALSE

## Value

Patchwork 'quilt'-like plot of ggplots
