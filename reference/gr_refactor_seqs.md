# Refactor seqinfo, seqnames, seqlengths, seqlevels of GenomicRanges object for easy harmony

Single command to refactor all seq details of a GenomicRanges object to
easily harmonize with any other GenomicRanges object. By default, this
package uses the autosome (1-22) and sex chromosomes (X,Y) of hg38, see
[`gUtils::hg_seqlengths()`](https://rdrr.io/pkg/gUtils/man/hg_seqlengths.html)
Users can adjust this using the `new_levels` parameter.

## Usage

``` r
gr_refactor_seqs(
  input_gr,
  new_levels = gUtils::hg_seqlengths(),
  sort_gr = TRUE
)
```

## Arguments

- input_gr:

  GenomicRanges object to refactor

- new_levels:

  Named vector object used as the template for new seq details, see
  [`gUtils::hg_seqlengths()`](https://rdrr.io/pkg/gUtils/man/hg_seqlengths.html)
  for example

- sort_gr:

  Sort the output GenomicRanges object

## Value

GenomicRanges object with updated seqinfo, seqnames, seqlengths,
seqlevels

## Examples

``` r
# Converting between reference genomes is complicated and proper seqinfo is hard
# to achieve without many lines of code, this function helps alleviate the issue

# Look at seqinfo from a hg19 dataset
GenomeInfoDb::seqinfo(kras_dnase_demo_gr_hg19)
#> Seqinfo object with 25 sequences from GRCh37 genome:
#>   seqnames seqlengths isCircular genome
#>   1         249250621      FALSE GRCh37
#>   2         243199373      FALSE GRCh37
#>   3         198295559      FALSE GRCh37
#>   4         191154276      FALSE GRCh37
#>   5         181538259      FALSE GRCh37
#>   ...             ...        ...    ...
#>   21         48129895      FALSE GRCh37
#>   22         51304566      FALSE GRCh37
#>   X         156040895      FALSE GRCh37
#>   Y          59373566      FALSE GRCh37
#>   M             16571      FALSE GRCh37

# Look at seqinfo after lift over to hg38
GenomeInfoDb::seqinfo(example_dnase_hg38)
#> Seqinfo object with 23 sequences from an unspecified genome; no seqlengths:
#>   seqnames seqlengths isCircular genome
#>   chr1           <NA>       <NA>   <NA>
#>   chr2           <NA>       <NA>   <NA>
#>   chr3           <NA>       <NA>   <NA>
#>   chr4           <NA>       <NA>   <NA>
#>   chr5           <NA>       <NA>   <NA>
#>   ...             ...        ...    ...
#>   chr19          <NA>       <NA>   <NA>
#>   chr20          <NA>       <NA>   <NA>
#>   chr21          <NA>       <NA>   <NA>
#>   chr22          <NA>       <NA>   <NA>
#>   chrX           <NA>       <NA>   <NA>

# Use gr_refactor_seqs to fix the inconsistencies
GenomeInfoDb::seqinfo(gr_refactor_seqs(example_dnase_hg38))
#> Seqinfo object with 24 sequences from GRCh38 genome:
#>   seqnames seqlengths isCircular genome
#>   chr1      248956422      FALSE GRCh38
#>   chr2      242193529      FALSE GRCh38
#>   chr3      198295559      FALSE GRCh38
#>   chr4      190214555      FALSE GRCh38
#>   chr5      181538259      FALSE GRCh38
#>   ...             ...        ...    ...
#>   chr20      64444167      FALSE GRCh38
#>   chr21      46709983      FALSE GRCh38
#>   chr22      50818468      FALSE GRCh38
#>   chrX      156040895      FALSE GRCh38
#>   chrY       57227415      FALSE GRCh38
```
