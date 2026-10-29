# Run diagnostic workflow of after JaBbA run to help evaluate quality of graph

A simple rework to the JaBbA:::QCStats function for more specific local
use cases Provides a text file and quilt-like plot with metrics on
convergence, RSME of copy number estimates, purity, ploidy, loose ends,
and junctions.

## Usage

``` r
get_jabba_qc_diagnostic(tumor_normal_id, jabba_workdir_path, verbose = T)
```

## Arguments

- tumor_normal_id:

  The name of the tumor/normal pair with completed `JaBbA` run

- jabba_workdir_path:

  Path to the completed `JaBbA` run work directory

- verbose:

  Output CLI during workflow, default: TRUE

## Value

data.table of QC metrics and quilt-like ggplot

## Examples

``` r
# get_jabba_qc_diagnostic(
# tumor_normal_id = "CG15-2647-T_vs_190585-N",
# jabba_workdir_path = "jAbBa/CG15-2647-T_vs_190585-N_jabba_workdir/")
```
