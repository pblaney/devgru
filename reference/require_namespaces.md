# Vectorised version of requireNamespace

Carbon copy of the helper function described by Raphael Sonabend-Friend
in the package `survivalmodels`. Vectorises the `requireNamespace`
function and returns `TRUE` if all packages, `x`, are available and
`FALSE` otherwise. For internal use only.
https://github.com/RaphaelS1/survivalmodels/blob/main/R/helpers.R

## Usage

``` r
require_namespaces(pkgs)
```

## Arguments

- pkgs:

  Vector of character strings with required package name(s)
