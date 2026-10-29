# Verbose CLI output to designate a workflow start

Display to the user the workflow has begun, the function that launched
it, and start the stopwatch to track execution time,

## Usage

``` r
cli_stopwatch_start(package = "devgru", function_name)
```

## Arguments

- package:

  What is the package of the function, default: devgru

- function_name:

  What is the function name

## Value

CLI message and start time of process run

## Examples

``` r
process_start <- cli_stopwatch_start(package = "devgru", function_name = "workflow_demo_func")
#> ── Start ───────────────────────────────────────────────────────────────────────
#> devgru :: workflow_demo_func() 2026-10-29 04:39:20.508491
process_start
#>    user  system elapsed 
#>  13.732   1.509  16.666 
```
