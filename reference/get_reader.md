# Read file using appropriate format-specific reader

Detect and utilize the necessary pre-built reader for each file of a
specific format. Otherwise, use generic data.table reader

## Usage

``` r
get_reader(
  file_path,
  file_format,
  delim,
  has_header,
  cpus,
  seq_lengths,
  verbose,
  ...
)
```

## Arguments

- file_path:

  Full path to file

- file_format:

  Character string to force specific format detection: "auto" (default),
  "bed", "bedpe", "maf", "vcf", "generic". Default: "auto"

- delim:

  Delimiter for generic files. Default: "auto"

- has_header:

  Logical, if files have header line. Default: TRUE

- cpus:

  Number of CPUs for reading data. Default: 1

- seq_lengths:

  Named vector for seqinfo refactoring (genomic formats only). Default:
  gUtils::hg_seqlengths()

- verbose:

  Logical, print progress messages. Default: FALSE

- ...:

  Additional arguments passed to format-specific readers
