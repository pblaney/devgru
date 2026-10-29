# Read in a BED file and convert to GenomicRanges or data.table

Read in a BED file and convert it to a GenomicRanges object or
data.table with refactored seq details. Expects the first 3 columns as
chromosome, start, end; However column names are not necessary. The BED
file can be either zipped or unzipped.

## Usage

``` r
read_bed_file(
  bed_file,
  output_type = c("gr", "dt"),
  has_header = FALSE,
  additional_col_names = NULL,
  cpus = 1,
  seq_lengths = gUtils::hg_seqlengths(),
  allow_nonstandard_chr = FALSE,
  verbose = TRUE
)
```

## Arguments

- bed_file:

  Path to BED file (supports .bed, .bed.gz)

- output_type:

  Character string specifying output format: gr="GenomicRanges"
  (default) or dt="data.table"

- has_header:

  Logical, does BED file have header line. Default: FALSE

- additional_col_names:

  Character vector of names for additional columns in BED file, beyond
  first three. Default: NULL

- cpus:

  Integer, number of CPUs for reading data, used by
  [`data.table::fread()`](https://rdrr.io/pkg/data.table/man/fread.html).
  Default: 1

- seq_lengths:

  Named vector object used as template for new seq details. See
  [`gUtils::hg_seqlengths()`](https://rdrr.io/pkg/gUtils/man/hg_seqlengths.html)
  for example. Default: hg38 seqlengths. Only used when output_type =
  "GenomicRanges"

- allow_nonstandard_chr:

  Logical, whether to allow non-standard chromosomes (decoys, HLA, alt
  contigs). Default: FALSE

- verbose:

  Logical, whether to print progress messages. Default: TRUE

## Value

GenomicRanges object or data.table with BED data. If GenomicRanges,
includes updated seqinfo, seqnames, seqlengths, seqlevels
