# Refactor seqinfo of GInteractions object for easy harmony

Single command to refactor all seq details of a GInteractions object to
easily harmonize with any other genomic object. This function updates
seqinfo for both anchor regions. By default, this package uses the
autosome (1-22) and sex chromosomes (X,Y) of hg38, see
[`gUtils::hg_seqlengths()`](https://rdrr.io/pkg/gUtils/man/hg_seqlengths.html)
Users can adjust this using the `new_levels` parameter.

## Usage

``` r
gi_refactor_seqs(
  input_gi,
  new_levels = gUtils::hg_seqlengths(),
  sort_gi = TRUE
)
```

## Arguments

- input_gi:

  GInteractions object to refactor

- new_levels:

  Named vector object used as the template for new seq details, see
  [`gUtils::hg_seqlengths()`](https://rdrr.io/pkg/gUtils/man/hg_seqlengths.html)
  for example

- sort_gi:

  Sort the output GInteractions object

## Value

GInteractions object with updated seqinfo for both anchors
