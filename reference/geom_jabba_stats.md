# Diagnostic plot used for evaluating the performance of a JaBbA run

Creates a ggplot quilt-like plot for quickly evaluating the run
performance and output graph characterisics of a JaBbA run.

Note, this function was designed to be run as part of the
[`get_jabba_qc_diagnostic()`](https://pblaney.github.io/devgru/reference/get_jabba_qc_diagnostic.md)
workflow.

## Usage

``` r
geom_jabba_stats(path_to_stats_dt, karyograph_dt)
```

## Arguments

- path_to_stats_dt:

  Path to the `qc_stats.txt` output from
  [`get_jabba_qc_diagnostic()`](https://pblaney.github.io/devgru/reference/get_jabba_qc_diagnostic.md)
  run

- karyograph_dt:

  Karyograph data.table from `JaBbA` run

## Value

ggplot object

## Examples

``` r
# Snippet from within the function `get_jabba_qc_diagnostic` that calls this geom
# diagnostic_plot <- geom_jabba_stats(
# path_to_stats_dt = paste0(jabba_workdir_path, "qc_stats.txt"),
# karyograph_dt = data.table::data.table("cn" = kar$segstats$cn,
#                                        "cnmle" = kar$segstats$cnmle))
```
