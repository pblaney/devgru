# Read BEDPE File with Flexible Output Format

Read in a BEDPE file and convert it to a GInteractions object with
refactored seq details. Flexible to different naming conventions but
expects the 6 columns for breakpoint chr1, start1, end1, and breakpoint
chr2, start2, end2. Can also return data.table or GRangesList for
convenience.

## Usage

``` r
read_bedpe_file(
  bedpe_file,
  output_type = c("gi", "dt", "grl"),
  keep_metadata = TRUE,
  verbose = TRUE
)
```

## Arguments

- bedpe_file:

  Path to BEDPE file or data.table/data.frame

- output_type:

  Character string specifying output format: gi="GInteractions"
  (default), dt="data.table", or grl="GRangesList"

- keep_metadata:

  Logical, whether to preserve metadata (default: TRUE)

- verbose:

  Logical, whether to print progress messages. (default: TRUE)

## Value

GInteractions, data.table, or GRangesList object depending on
output_type
