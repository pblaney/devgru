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

# Execute the full workflow
# fragcounter_seg <- get_fragcounter_segmentation(
# path_to_fragcounter_profile = "read_depth_demo_hg38.rds",
# cpus = 4)
```
