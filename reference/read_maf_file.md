# Read MAF file and convert to GenomicRanges or data.table

Read in a MAF file which contains a number of columns and convert it to
a GenomicRanges object or data.table. The MAF file can be either zipped
or unzipped. This function performs validation and filtering to ensure
proper MAF format compliance. For more specific MAFtools operations, see
`maftools::read.maf()`

## Usage

``` r
read_maf_file(
  maf_file,
  output_type = c("gr", "dt"),
  cpus = 2,
  seq_lengths = gUtils::hg_seqlengths(),
  verbose = TRUE
)
```

## Arguments

- maf_file:

  Path to MAF file (supports .maf, .maf.gz, .maf.txt, .txt)

- output_type:

  Character string specifying output format: gr="GenomicRanges"
  (default) or dt="data.table"

- cpus:

  Number of CPUs for reading data, used by
  [`data.table::fread()`](https://rdrr.io/pkg/data.table/man/fread.html).
  Default: 2

- seq_lengths:

  Named vector object used as template for new seq details. See
  [`gUtils::hg_seqlengths()`](https://rdrr.io/pkg/gUtils/man/hg_seqlengths.html)
  for example. Default: hg38 seqlengths. Only used when output_type =
  "GenomicRanges"

- verbose:

  Logical, whether to print progress messages. Default: TRUE

## Value

GenomicRanges object or data.table with MAF data. If GenomicRanges,
includes updated seqinfo, seqnames, seqlengths, seqlevels
