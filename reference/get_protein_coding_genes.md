# Extract protein-coding genes from GTF file

Convenience function to extract protein-coding gene features from a GTF
file. Returns a GenomicRanges or data.table object with only
protein-coding genes. This is a wrapper around read_gtf_file with preset
filters.

## Usage

``` r
get_protein_coding_genes(
  gtf_file,
  output_type = c("gr", "dt"),
  seq_lengths = gUtils::hg_seqlengths(),
  verbose = TRUE
)
```

## Arguments

- gtf_file:

  Path to GTF file (supports .gtf, .gtf.gz)

- output_type:

  Character string specifying output format: gr="GenomicRanges"
  (default) or dt="data.table"

- seq_lengths:

  Named vector object used as template for new seq details. See
  [`gUtils::hg_seqlengths()`](https://rdrr.io/pkg/gUtils/man/hg_seqlengths.html)
  for example. Default: hg38 seqlengths. Only used when output_type =
  "GenomicRanges"

- verbose:

  Logical, whether to print progress messages. Default: TRUE

## Value

GenomicRanges object or data.table with protein-coding genes only
