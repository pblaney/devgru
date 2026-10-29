# Read in CNV profile data.table of various flavors and extract a corrected profile

Given a data.table or GenomicRanges CNV object, extract the corrected
CNV profile. Any segment with a non-rounded value within 0.2 of the next
integer value is rounded to that value. The output will data.table will
contain 6 columns: sample, seqnames, start, end, total, minor The
`sample` is either user-provided or row count placeholder Currently
supports CNV calls from Battenberg and FACETS

## Usage

``` r
get_corrected_cnv_profile(cnv_obj, caller, sample_id = NULL)
```

## Arguments

- cnv_obj:

  CNV file in data.table or GenomicRanges format

- caller:

  Name of caller that generated input CNV to be converted, supported:
  Battenberg and FACETS

- sample_id:

  Unique identifier to add to output, default: NULL

## Value

data.table object with corrected CNV segments
