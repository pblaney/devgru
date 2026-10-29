# Convert BEDPE-like data.table to GInteractions

Single command to smartly convert a BEPDE-like data.table to a
GIinteractions object while preserving all metadata.

## Usage

``` r
dt_to_gi(input_dt, keep_metadata = TRUE)
```

## Arguments

- input_dt:

  data.table with BEDPE-like columns

- keep_metadata:

  Logical, whether to preserve metadata (default: TRUE)

## Value

GInteractions object
