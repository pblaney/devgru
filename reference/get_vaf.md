# Quick pull or explicitly calculate the VAF for mutation records of various flavors

Given a data.table or GenomicRanges VCF object, quickly extract or
explicitly calculate the VAF for all mutations. For clarity, the read
support for the VAF will be extracted as well. Currently supports
somatic SNV/InDel VCFs from Mutect, Strelka, Varscan, SvABA, CaVEMan.

## Usage

``` r
get_vaf(vcf_obj, caller, mut_type)
```

## Arguments

- vcf_obj:

  VCF file in data.table or GenomicRanges format

- caller:

  Name of the caller that generated input VCF to be converted,
  supported: mutect, strelka, varscan, svaba, caveman

- mut_type:

  Type of mutations within VCF, supported: snv, indel

## Value

data.table object with read support and VAF per record
