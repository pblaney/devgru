# Create color palette for integer copy number states

Given a set of discrete integer copy number states, produce a diverging
blue/red color palette that scales to any number.

## Usage

``` r
copynumber_palettier(cn_states)
```

## Arguments

- cn_states:

  Discrete integer values of copy number states, typically for plotting

## Value

data.table of states and assigned color

## Examples

``` r
# Most common set of states
copynumber_palettier(cn_states = seq(0,4))
#>    copynumber   color
#>        <fctr>  <char>
#> 1:          0 #264DFF
#> 2:          1 #72D9FF
#> 3:          2 #999999
#> 4:          3 #F4B86E
#> 5:          4 #A10E39

# high level jump
copynumber_palettier(cn_states = c(1,2,3,8))
#>    copynumber   color
#>        <fctr>  <char>
#> 1:          1 #72D9FF
#> 2:          2 #999999
#> 3:          3 #F4B86E
#> 4:          8 #A10E39
```
