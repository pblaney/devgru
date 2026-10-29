# Check if input is a GenomicRanges object, plus optional sanity check of column names

Simple check if the input is a GenomicRanges object, plus includes an
additional functionality to perform a smart, transparent sanity check of
expected column names.

## Usage

``` r
gr_sanitycheck(query_gr, expected_cols = NULL)
```

## Arguments

- query_gr:

  Suspected GenomicRanges-like object to test

- expected_cols:

  Vector of character strings to check for consistency in columns of
  query GenomicRanges, default: NULL

## Value

TRUE or FALSE

## Examples

``` r
# Not a GR
gr_sanitycheck(query_gr = braf_dnase_demo_dt_hg38)
#> [1] FALSE
# Just a GR
gr_sanitycheck(query_gr = kras_dnase_demo_gr_hg38)
#> [1] TRUE
# GR with correct columns
gr_sanitycheck(query_gr = kras_dnase_demo_gr_hg38,
               expected_cols = c("signalValue","pValue","biospecimen","gene"))
#> [1] TRUE
# GR with less columns than expected
gr_sanitycheck(query_gr = kras_dnase_demo_gr_hg38[,-3],
               expected_cols = c("signalValue","pValue","biospecimen","gene"))
#> # A tibble: 1 × 8
#>   seqnames    start      end strand width signalValue pValue gene 
#>   <chr>       <int>    <int> <chr>  <dbl>       <dbl>  <dbl> <chr>
#> 1 chr12    25241134 25241673 *        540        13.8   22.7 KRAS 
#> ✖ Input GenomicRanges object had less columns than expected - missing columns: biospecimen
#> [1] FALSE
```
