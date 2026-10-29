# Precise addition of flanking to GenomicRanges with boundary-aware options

Precise addition of flanks to the start and/or end of a GenomicRanges
object with the option of setting a boundary for the flanks, helpful to
account for boundaries like chromosome arms.

Note, the boundaries are compared directly to the new flanked start/end
and do not currently take the chromosome into account automatically.
Thus, this function is best applied in on per chromosome basis.

If the flank exceeds the boundary, the new flanked start/end will set at
the boundary.

## Usage

``` r
gr_flank(
  input_gr,
  start_flank = NULL,
  end_flank = NULL,
  start_flank_boundary = NULL,
  end_flank_boundary = NULL
)
```

## Arguments

- input_gr:

  GenomicRanges object to add flanks to

- start_flank:

  Amount of flank (in bp) to be added to the start of the ranges,
  default: NULL

- end_flank:

  Amount of flank (in bp) to be added to the end of the ranges, default:
  NULL

- start_flank_boundary:

  The genomic coordinate the start of the flanked GenomicRanges should
  to not exceed, default: NULL

- end_flank_boundary:

  The genomic coordinate the end of the flanked GenomicRanges should to
  not exceed, default: NULL

## Value

GenomicRanges object with updated start/end coordinates

## Examples

``` r
# The gene NOTCH2 is close to the p-arm centromere on chr1
notch2 <- gUtils::parse.gr("chr1:119911553-120069662")
notch2
#> GRanges object with 1 range and 2 metadata columns:
#>       seqnames              ranges strand |    grl.ix   grl.iix
#>          <Rle>           <IRanges>  <Rle> | <integer> <integer>
#>   [1]     chr1 119911553-120069662      * |         1         1
#>   -------
#>   seqinfo: 25 sequences from an unspecified genome
chr1_p_arm <- gUtils::`%Q%`(chromosome_arms_hg38, seqnames == "chr1" & arm == "p")
chr1_p_arm
#> GRanges object with 1 range and 1 metadata column:
#>       seqnames      ranges strand |         arm
#>          <Rle>   <IRanges>  <Rle> | <character>
#>   [1]     chr1 1-122026459      * |           p
#>   -------
#>   seqinfo: 24 sequences from GRCh38 genome
GenomicRanges::end(chr1_p_arm)
#> [1] 122026459

# If sampling data around this locus in a 3 Mb window, it's important to not
# extend into/beyond the centromere so we need to set this boundary
gr_flank(input_gr = notch2,
         start_flank = 3e6,
         end_flank = 3e6,
         start_flank_boundary = GenomicRanges::start(chr1_p_arm),
         end_flank_boundary = GenomicRanges::end(chr1_p_arm))
#> GRanges object with 1 range and 2 metadata columns:
#>       seqnames              ranges strand |    grl.ix   grl.iix
#>          <Rle>           <IRanges>  <Rle> | <integer> <integer>
#>   [1]     chr1 116911553-122026459      * |         1         1
#>   -------
#>   seqinfo: 24 sequences from GRCh38 genome
```
