# Verbose CLI output to designate a workflow end

Display to the user the workflow has ended, stop the stopwatch, and
display execution time.

## Usage

``` r
cli_stopwatch_end(package = "devgru", function_name, stopwatch_start)
```

## Arguments

- package:

  What is the package of the function, default: devgru

- function_name:

  What is the function name

- stopwatch_start:

  What was the start time of the stopwatch, typically captured and
  passed along by
  [`cli_stopwatch_start()`](https://pblaney.github.io/devgru/reference/cli_stopwatch_start.md),
  see examples

## Examples

``` r
process_start <- cli_stopwatch_start(package = "devgru",
                                     function_name = "workflow_demo_func")
#> ── Start ───────────────────────────────────────────────────────────────────────
#> devgru :: workflow_demo_func() 2026-10-29 04:39:20.335516
process_end <- cli_stopwatch_end(package = "devgru",
                                 function_name = "workflow_demo_func",
                                 stopwatch_start = process_start)
#> ── Stop ────────────────────────────────────────────────────────────────────────
#> devgru :: workflow_demo_func() 2026-10-29 04:39:20.346553
#> Error in "cli::cli_alert_info(\"Duration {crayon::white(clisymbols::symbol$ellipsis)}\")": ! Could not evaluate cli `{}` expression: `crayon::white(clisy…`.
#> Caused by error in `loadNamespace(x)`:
#> ! there is no package called ‘clisymbols’
```
