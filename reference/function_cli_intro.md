# Verbose CLI output for extended and rich description of steps along a complex workflow, typically used in a pipeline

Display to the user the main package, the function being executed, and
the parameters passed to the function.

## Usage

``` r
function_cli_intro(package = "devgru", function_name, ...)
```

## Arguments

- package:

  What is the package of the function, default: devgru

- function_name:

  What is the function name

- ...:

  What parameters were used in the function

## Examples

``` r
# For demo purposes, set these variables
parameter_x <- 100
parameter_y <- "demo"

# In practice, these variables are set within the function being called as in below:
# workflow_demo_func <- function(parameter_x, parameter_y)
# then the function_cli_intro is called within the function

function_cli_intro(package = "devgru",
                   function_name = "workflow_demo_func",
                   parameter_x, parameter_y)
#>  [ devgru v0.99.0 ] 
#> 
#> ────────────────────────────────────────────────────────────────────────────────
#> 
#> Error in "cli::cli_alert(text = \"Launching {.emph {.pkg {function_name}()}} {crayon::white(clisymbols::symbol$ellipsis)}\")": ! Could not evaluate cli `{}` expression: `crayon::white(clisy…`.
#> Caused by error in `loadNamespace(x)`:
#> ! there is no package called ‘clisymbols’
```
