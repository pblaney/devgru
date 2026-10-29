# Convert segmentation GenomicRanges object to a data.table with standardized columns and format

Converts a segmentation GenomicRanges object into a data.table object
then extracts the appropriate columns and renames the columns to match
the expected format.

The input GenomicRanges is typically generated from
[`get_fragcounter_segmentation()`](https://pblaney.github.io/devgru/reference/get_fragcounter_segmentation.md)
workflow or other `DNAcopy` CBS algorithm runs.

## Usage

``` r
gr_to_seg(input_gr, exp_colnames = c("ID", "num.mark", "seg.mean"))
```

## Arguments

- input_gr:

  Segmentation GenomicRanges object to convert to data.table

- exp_colnames:

  Expected column names for CBS segmentation input GenomicRanges

## Value

data.table

## Examples

``` r
# Convert output from `get_fragcounter_segmentation()` workflow
# frag_seg <- get_fragcounter_segmentation(
# path_to_fragcounter_profile = "read_depth_demo_hg38.rds",
# cpus = 4)

# gr_to_seg(dryclean_seg)
```
