# Convert GenomicRangesList to GInteractions

Converts a BEDPE representation GenomicRangesList object into a
GInteractions object while preserving object metadata

## Usage

``` r
grl_to_gi(input_grl, keep_metadata = TRUE)
```

## Arguments

- input_grl:

  GenomicRangesList where each element contains exactly 2 ranges

- keep_metadata:

  Logical, whether to preserve metadata (default: TRUE)

## Value

GInteractions object
