# Pulls in package data for use in package functions

Carbon copy of the the helper function described by stackoverflow user
henfiber. Pulls data into a temporary new environment and accesses it
for assignment to a variable in a function. For internal use only.
https://stackoverflow.com/questions/30951204/

## Usage

``` r
get_data(name_of_data)
```

## Arguments

- name_of_data:

  Character string with name of data as it appears in package
