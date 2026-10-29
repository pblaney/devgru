# Single chromosome run circular binary segmentation (CBS) algorithm on read depth profile

Run CBS algorithm using DNAcopy on a per-chromosome basis on read depth
profile. This allows for parallel computation of the segmentation,
rapidly reducing run time. There should be no `NA`s in the columns or
`0`s in the `signal` column.

Note, this function was designed to be run as part of the
[`get_fragcounter_segmentation()`](https://pblaney.github.io/devgru/reference/get_fragcounter_segmentation.md)
workflow. Also, see
[`chrompar()`](https://pblaney.github.io/devgru/reference/chrompar.md)
for executing in parallel.

## Usage

``` r
get_cbs_per_chromosome(
  chrom_for_cbs,
  chromosome_names,
  signal,
  position,
  sample_id
)
```

## Arguments

- chrom_for_cbs:

  The chromosome to run CBS algorithm on

- chromosome_names:

  A vector of chromosome strings, equivalent to the seqnames column.
  Must contain at least the `chrom_for_cbs`

- signal:

  A vector of read depth signal to be used as input. These values will
  be `log`ged before use in CBS

- position:

  A vector of the genomic position of the signal measurement, equivalent
  to the start column

- sample_id:

  A string used to populate the ID column in the data.table

## Value

data.table of CBS segmentation values per chromosome

## Examples

``` r
# For this example, there are no NAs or zeros in the `signal`
# Check for these and may need to filter them out
read_depth_demo_hg38
#> GRanges object with 20002 ranges and 8 metadata columns:
#>           seqnames            ranges strand | background.log foreground.log
#>              <Rle>         <IRanges>  <Rle> |      <numeric>      <numeric>
#>       [1]    chr19          23000000      * |       0.127594     -0.0394919
#>       [2]    chr19 23000001-23000200      * |      -0.311801     -0.0394919
#>       [3]    chr19 23000201-23000400      * |      -0.280075      0.2346718
#>       [4]    chr19 23000401-23000600      * |      -0.261436     -0.0394919
#>       [5]    chr19 23000601-23000800      * |      -0.261440     -0.0394919
#>       ...      ...               ...    ... .            ...            ...
#>   [19998]    chr21 24999001-24999200      * |     -0.3118023     -0.0559702
#>   [19999]    chr21 24999201-24999400      * |     -0.2614352     -0.0559702
#>   [20000]    chr21 24999401-24999600      * |     -0.3118012     -0.0559702
#>   [20001]    chr21 24999601-24999800      * |      0.1297325     -0.0559702
#>   [20002]    chr21 24999801-25000000      * |      0.0268664     -0.0559702
#>           input.read.counts median.chr foreground background  log.reads
#>                   <numeric>  <numeric>  <numeric>  <numeric>  <numeric>
#>       [1]          1.155556   0.961278   0.961278   1.136092  0.1445814
#>       [2]          0.386038   0.961278   0.961278   0.732127 -0.9518190
#>       [3]          0.955857   0.961278   1.264494   0.755727 -0.0451467
#>       [4]          1.447640   0.961278   0.961278   0.769945  0.3699345
#>       [5]          0.927922   0.961278   0.961278   0.769942 -0.0748075
#>       ...               ...        ...        ...        ...        ...
#>   [19998]          0.896007   0.945567   0.945567   0.732126 -0.1098068
#>   [19999]          1.567746   0.945567   0.945567   0.769946  0.4496388
#>   [20000]          0.883972   0.945567   0.945567   0.732127 -0.1233302
#>   [20001]          0.938086   0.945567   0.945567   1.138524 -0.0639141
#>   [20002]          1.217163   0.945567   0.945567   1.027231  0.1965224
#>           germline.status
#>                 <logical>
#>       [1]            TRUE
#>       [2]            TRUE
#>       [3]           FALSE
#>       [4]            TRUE
#>       [5]            TRUE
#>       ...             ...
#>   [19998]            TRUE
#>   [19999]            TRUE
#>   [20000]            TRUE
#>   [20001]            TRUE
#>   [20002]            TRUE
#>   -------
#>   seqinfo: 24 sequences from GRCh38 genome

# Example of what the input looks like for
# chromosome_names
# head(as.character(GenomicRanges::seqnames(read_depth_demo_hg38)))
# signal
# head(as.double(GenomicRanges::values(read_depth_demo_hg38)[, "reads.corrected"]))
# position
# head(GenomicRanges::start(read_depth_demo_hg38))
```
