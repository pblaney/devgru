# Convert GInteractions to GRangesList

Single command to smartly convert GIinteractions to a GRangesList object
while preserving all metadata.

## Usage

``` r
gi_to_grl(input_gi, keep_metadata = TRUE)
```

## Arguments

- input_gi:

  GInteractions object

- keep_metadata:

  Logical, whether to preserve metadata (default: TRUE)

## Value

GRangesList where each element contains a pair of breakpoints
