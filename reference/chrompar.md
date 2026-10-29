# Scaffold function around `easypar` to run workflows in parallel by chromosome

A scaffold function to be used within other complex workflow functions
that provides easy API to `easypar` functionality and executes a
function in parallel by chromosome

The scaffold takes a function where the input can be subset by
chromosome and the expected output is either a data.table or
GenomicRanges object

## Usage

``` r
chrompar(
  par_function,
  par_chromosomes,
  ...,
  run_in_par = T,
  par_cpus = NULL,
  par_packages = NULL
)
```

## Arguments

- par_function:

  The function that will be executed in parallel by chromosome, written
  as the function name without the `()` a the end

- par_chromosomes:

  The set of chromosomes to execute in parallel, non-duplicated values
  and can be any chromosomes as long as present in input data

- ...:

  The input parameters that will be passed to `par_function` for
  execution, these must be provided in the order expected by
  `par_function`

- run_in_par:

  Run the jobs in parallel and not sequentially, default: TRUE

- par_cpus:

  Number of CPUs to use for parallel execution, `easypar` will safeguard
  against using all available CPUs

- par_packages:

  Vector of character strings designating special packages needed for
  executing the `par_function`

## Value

GenomicRanges object with properly sorted genomic coordinates

## Examples

``` r
# Run as part of the `get_fragcounter_segmentation()` workflow using
# `get_cbs_per_chromosome()` as the par_function
# final_cbs_segmentation <- chrompar(
# par_function = get_cbs_per_chromosome,
# par_chromosomes = chrom_iter_list,
# par_cpus = threads,
# seqnames_nona, rdp_cov_nona, start_nona, sample_id)
```
