# Refactor seqinfo of GRangesList object for easy harmony

Single command to refactor all seq details of a GRangesList object to
easily harmonize with any other genomic object. This function updates
seqinfo for all GRanges elements in the list. By default, this package
uses the autosome (1-22) and sex chromosomes (X,Y) of hg38, see
[`gUtils::hg_seqlengths()`](https://rdrr.io/pkg/gUtils/man/hg_seqlengths.html)
Users can adjust this using the `new_levels` parameter.

## Usage

``` r
grl_refactor_seqs(input_grl, new_levels = gUtils::hg_seqlengths())
```

## Arguments

- input_grl:

  GRangesList object to refactor

- new_levels:

  Named vector object used as the template for new seq details, see
  [`gUtils::hg_seqlengths()`](https://rdrr.io/pkg/gUtils/man/hg_seqlengths.html)
  for example

## Value

GRangesList object with updated seqinfo for all elements
