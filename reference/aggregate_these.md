# Read and aggregate genomic or non-genomic data files with smart format detection

Collect all files that match a specific pattern at a specific path,
detect the file format, use the appropriate reader function, then
aggregate into a single data.table. Automatically detects and handles
BED, BEDPE, MAF, VCF, and generic tabular formats. All files are read as
data.tables and aggregated, ensuring maximum compatibility and
performance.

## Usage

``` r
aggregate_these(
  path_to_files,
  pattern_to_grab,
  file_format = c("auto", "bed", "bedpe", "maf", "vcf", "generic"),
  delim = "auto",
  has_header = TRUE,
  cpus = 1,
  add_uniq_id = FALSE,
  id_column_name = "sample_id",
  verbose = TRUE,
  ...
)
```

## Arguments

- path_to_files:

  Path to location of files to be aggregated

- pattern_to_grab:

  Pattern (regex) used to identify files

- file_format:

  Character string to force specific format detection: "auto" (default),
  "bed", "bedpe", "maf", "vcf", "generic". Default: "auto"

- delim:

  Delimiter for generic files. Default: "auto"

- has_header:

  Logical, if files have header line. Default: TRUE

- cpus:

  Number of CPUs for reading data. Default: 1

- add_uniq_id:

  Logical, add unique identifier column from filename. Default: FALSE

- id_column_name:

  Character, name for unique identifier column. Default: "sample_id"

- verbose:

  Logical, print progress messages. Default: TRUE

- ...:

  Additional arguments passed to format-specific readers

## Value

data.table object with aggregated data from all files
