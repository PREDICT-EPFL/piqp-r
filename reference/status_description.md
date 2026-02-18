# Return the solver status description string

Return the solver status description string

## Usage

``` r
status_description(code)
```

## Arguments

- code:

  a valid solver return code

## Value

a status description string

## Examples

``` r
status_description(1) ## for solved problem
#> [1] "Solver solved problem up to given tolerance."
status_description(-1) ## for max iterations limit reached
#> [1] "Iteration limit was reached."
```
