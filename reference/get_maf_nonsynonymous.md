# Extract Non-Synonymous Mutations from MAF

Filters a MAF GenomicRanges object to retain only non-synonymous
mutations, including both SNVs and InDels and splice site variantss.

Variant Class: SNV, insertion, deletion

Non-Synonymous Classifications: Missense_Mutation, Nonsense_Mutation,
Frame_Shift_Del, Frame_Shift_Ins, In_Frame_Del, In_Frame_Ins,
Nonstop_Mutation, Translation_Start_Site, Splice_Site, Splice_Region

Note, this function was designed to be run downstream of
[`read_maf_file()`](https://pblaney.github.io/devgru/reference/read_maf_file.md)

## Usage

``` r
get_maf_nonsynonymous(
  maf_gr,
  variant_class = c("SNV", "insertion", "deletion"),
  include_splice = TRUE
)
```

## Arguments

- maf_gr:

  GenomicRanges object containing MAF data with metadata columns

- variant_class:

  Character vector of variant class to include. Default: c("SNV",
  "insertion", "deletion")

- include_splice:

  Logical, whether to include splice site variants (default: TRUE)

## Value

A GenomicRanges object containing only non-synonymous mutations
