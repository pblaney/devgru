# Generate diagnostic plots for alignment QC checks using Alfred summary files

Single command to generate a comprehensive set of diagnostic plots that
assist in performing QC checks of tumor/normal read alignment of various
sequencing protocol flavors using Alfred summary files.

## Usage

``` r
get_qc_diagnostics_alignment(
  path_to_tumor_dir = NULL,
  path_to_normal_dir = NULL,
  dataset_title = "",
  seq_center = "",
  seq_protocol = "WGS",
  output_path = getwd(),
  verbose = T
)
```

## Arguments

- path_to_tumor_dir:

  Path to directory of tumor sample Alfred summary files

- path_to_normal_dir:

  Path to directory of normal sample Alfred summary files

- dataset_title:

  Name of the dataset

- seq_center:

  Institute/site where sequencing was performed

- seq_protocol:

  Type of sequencing protocol for display purposes, default: WGS

- output_path:

  Location to save report PDF, default:
  [`getwd()`](https://rdrr.io/r/base/getwd.html)

- verbose:

  Output CLI during workflow, default: TRUE

## Value

Patchwork 'quilt'-like plot of ggplots

## Examples

``` r
# After downloading `*.qc.summary.txt` Alfred output file from the MGP1000 run,
# place tumor and normal sample files into respective subdirectories
# get_qc_diagnostics_alignment(
# path_to_tumor_dir = "data/sanger/examplePDexport_999/preprocessing/alfred/tumor/",
# path_to_normal_dir = "data/sanger/examplePDexport_999/preprocessing/alfred/normal/",
# dataset_title = "examplePDexport_999",
# seq_center = "Wellcome Sanger Institute",
# output_path = "data/sanger/examplePDexport_999/preprocessing/")
```
