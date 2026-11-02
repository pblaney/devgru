# Generate a segmentation file from a fragCounter read depth profile

Given fragCounter read depth profile, run DNAcopy's CBS algorithm in a
parallel manner to produce a segmentation file that is compatible with
GISTIC2.0.

## Usage

``` r
get_fragcounter_segmentation(
  path_to_fragcounter_profile,
  cpus = 1,
  random_seed = 999,
  verbose = T,
  exp_colnames = c("reads", "gc", "map", "reads.corrected")
)
```

## Arguments

- path_to_fragcounter_profile:

  Path to fragCounter read depth profile, .rds file

- cpus:

  Number of CPUs to use for parallel execution, default: 1

- random_seed:

  Numeric value to use as random seed for consistent results

- verbose:

  Output CLI during workflow, default: TRUE

- exp_colnames:

  Expected column names for fragCounter read depth profile input

## Value

data.table object with required GISTIC2.0 columns

## Examples

``` r
# Demo fragCounter read depth profile
read_depth_demo_hg38
#> GRanges object with 20002 ranges and 8 metadata columns:
#>           seqnames            ranges strand | background.log
#>              <Rle>         <IRanges>  <Rle> |      <numeric>
#>       [1]    chr19          23000000      * |       0.127594
#>       [2]    chr19 23000001-23000200      * |      -0.311801
#>       [3]    chr19 23000201-23000400      * |      -0.280075
#>       [4]    chr19 23000401-23000600      * |      -0.261436
#>       [5]    chr19 23000601-23000800      * |      -0.261440
#>       ...      ...               ...    ... .            ...
#>   [19998]    chr21 24999001-24999200      * |     -0.3118023
#>   [19999]    chr21 24999201-24999400      * |     -0.2614352
#>   [20000]    chr21 24999401-24999600      * |     -0.3118012
#>   [20001]    chr21 24999601-24999800      * |      0.1297325
#>   [20002]    chr21 24999801-25000000      * |      0.0268664
#>           foreground.log input.read.counts median.chr foreground
#>                <numeric>         <numeric>  <numeric>  <numeric>
#>       [1]     -0.0394919          1.155556   0.961278   0.961278
#>       [2]     -0.0394919          0.386038   0.961278   0.961278
#>       [3]      0.2346718          0.955857   0.961278   1.264494
#>       [4]     -0.0394919          1.447640   0.961278   0.961278
#>       [5]     -0.0394919          0.927922   0.961278   0.961278
#>       ...            ...               ...        ...        ...
#>   [19998]     -0.0559702          0.896007   0.945567   0.945567
#>   [19999]     -0.0559702          1.567746   0.945567   0.945567
#>   [20000]     -0.0559702          0.883972   0.945567   0.945567
#>   [20001]     -0.0559702          0.938086   0.945567   0.945567
#>   [20002]     -0.0559702          1.217163   0.945567   0.945567
#>           background  log.reads germline.status
#>            <numeric>  <numeric>       <logical>
#>       [1]   1.136092  0.1445814            TRUE
#>       [2]   0.732127 -0.9518190            TRUE
#>       [3]   0.755727 -0.0451467           FALSE
#>       [4]   0.769945  0.3699345            TRUE
#>       [5]   0.769942 -0.0748075            TRUE
#>       ...        ...        ...             ...
#>   [19998]   0.732126 -0.1098068            TRUE
#>   [19999]   0.769946  0.4496388            TRUE
#>   [20000]   0.732127 -0.1233302            TRUE
#>   [20001]   1.138524 -0.0639141            TRUE
#>   [20002]   1.027231  0.1965224            TRUE
#>   -------
#>   seqinfo: 24 sequences from GRCh38 genome

# Execute the full workflow
# fragcounter_seg <- get_fragcounter_segmentation(
# path_to_fragcounter_profile = "read_depth_demo_hg38.rds",
# cpus = 4)
```
