# Convert data.table to GenomicRanges Object

Single command to smartly convert a data.table object to a GenomicRanges
object by wrapping
[`gr_refactor_seqs()`](https://pblaney.github.io/devgru/reference/gr_refactor_seqs.md)
around the
[`gUtils::dt2gr()`](https://rdrr.io/pkg/gUtils/man/dt2gr.html) to avoid
seqinfo conflicts and apply proper sorting.

## Usage

``` r
dt_to_gr(input_dt)
```

## Arguments

- input_dt:

  data.table object that minimally contains columns like chromosome
  start, and end position that describe the genomic coordinates of the
  range

## Value

GenomicRanges object with input columns and updated seqinfo, seqnames,
seqlengths, seqlevels

## Examples

``` r
# Converting from data.table to GenomicRanges can be easily done with gUtils::dt2gr()
# however this single command conversion leaves some loose ends with the seqinfo

# Look at the original data.table
braf_dnase_demo_dt_hg38
#>    seqnames     start       end strand signalValue  pValue biospecimen   gene
#>      <char>     <int>     <int> <char>       <num>   <num>      <char> <char>
#> 1:     chr7 140740305 140740552      *     2.16924 1.59134        K562   KRAS
#> 2:     chr7 140797423 140798072      *     2.17796 1.59886        K562   KRAS
#> 3:     chr7 140826433 140826768      *     2.59020 1.98645        K562   KRAS
#> 4:     chr7 140842477 140843425      *     6.56378 7.66134        K562   KRAS
# Now convert to GRanges
dt_to_gr(braf_dnase_demo_dt_hg38)
#> GRanges object with 4 ranges and 4 metadata columns:
#>       seqnames              ranges strand | signalValue    pValue biospecimen
#>          <Rle>           <IRanges>  <Rle> |   <numeric> <numeric> <character>
#>   [1]     chr7 140740305-140740552      * |     2.16924   1.59134        K562
#>   [2]     chr7 140797423-140798072      * |     2.17796   1.59886        K562
#>   [3]     chr7 140826433-140826768      * |     2.59020   1.98645        K562
#>   [4]     chr7 140842477-140843425      * |     6.56378   7.66134        K562
#>              gene
#>       <character>
#>   [1]        KRAS
#>   [2]        KRAS
#>   [3]        KRAS
#>   [4]        KRAS
#>   -------
#>   seqinfo: 24 sequences from GRCh38 genome
# Check the seqinfo
GenomeInfoDb::seqinfo(dt_to_gr(braf_dnase_demo_dt_hg38))
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

# Compared to the old way
# Conversion looks fine
gUtils::dt2gr(braf_dnase_demo_dt_hg38)
#> GRanges object with 4 ranges and 4 metadata columns:
#>       seqnames              ranges strand | signalValue    pValue biospecimen
#>          <Rle>           <IRanges>  <Rle> |   <numeric> <numeric> <character>
#>   [1]     chr7 140740305-140740552      * |     2.16924   1.59134        K562
#>   [2]     chr7 140797423-140798072      * |     2.17796   1.59886        K562
#>   [3]     chr7 140826433-140826768      * |     2.59020   1.98645        K562
#>   [4]     chr7 140842477-140843425      * |     6.56378   7.66134        K562
#>              gene
#>       <character>
#>   [1]        KRAS
#>   [2]        KRAS
#>   [3]        KRAS
#>   [4]        KRAS
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
# However the seqinfo will conflict with other hg38 GRanges
GenomeInfoDb::seqinfo(gUtils::dt2gr(braf_dnase_demo_dt_hg38))
#> Seqinfo object with 1 sequence from an unspecified genome:
#>   seqnames seqlengths isCircular genome
#>   chr7      140843425         NA   <NA>
```
