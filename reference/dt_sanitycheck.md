# Check if input is a data.table object, plus optional sanity check of column names

Simple check if the input is a data.table object, plus includes an
additional functionality to perform a smart, transparent sanity check of
expected column names.

## Usage

``` r
dt_sanitycheck(query_dt, expected_cols = NULL)
```

## Arguments

- query_dt:

  Suspected data.table-like object to test

- expected_cols:

  Vector of character strings to check for consistency in columns of
  query data.table, default: NULL

## Value

TRUE or FALSE

## Examples

``` r
# Not a DT
dt_sanitycheck(query_dt = kras_dnase_demo_gr_hg38)
#> [1] FALSE
# Just a DT
dt_sanitycheck(query_dt = braf_dnase_demo_dt_hg38)
#> [1] TRUE
# DT with correct columns
dt_sanitycheck(query_dt = braf_dnase_demo_dt_hg38,
               expected_cols = c("seqnames","start","end","strand","signalValue",
                                 "pValue","biospecimen","gene"))
#> [1] TRUE
# DT with less columns than expected
dt_sanitycheck(query_dt = braf_dnase_demo_dt_hg38[,-7],
               expected_cols = c("seqnames","start","end","strand","signalValue",
                                 "pValue","biospecimen","gene"))
#> # A tibble: 4 × 7
#>   seqnames     start       end strand signalValue pValue gene 
#>   <chr>        <int>     <int> <chr>        <dbl>  <dbl> <chr>
#> 1 chr7     140740305 140740552 *             2.17   1.59 KRAS 
#> 2 chr7     140797423 140798072 *             2.18   1.60 KRAS 
#> 3 chr7     140826433 140826768 *             2.59   1.99 KRAS 
#> 4 chr7     140842477 140843425 *             6.56   7.66 KRAS 
#> ✖ Input data.table object had less columns than expected - missing columns: biospecimen
#> [1] FALSE
```
