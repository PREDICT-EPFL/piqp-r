# The PIQP Solver Model Class

An S7 class wrapping the PIQP C++ Solver. Users will never need to
directly create instances of this class and should use the more
user-friendly functions [`piqp()`](piqp.md) and
[`solve_piqp()`](solve_piqp.md).

## Usage

``` r
piqp_model(solver_ptr = NULL, dims = list(), dense_backend = logical(0))
```

## Arguments

- solver_ptr:

  external pointer to the C++ solver object

- dims:

  a named list with elements `n`, `p`, and `m`

- dense_backend:

  logical flag indicating if the dense solver is used

## Value

An S7 object of class `piqp_model` with properties `solver_ptr` (an
external pointer to the C++ solver), `dims` (a named list with elements
`n`, `p`, and `m` giving the number of variables, equality constraints,
and inequality constraints), and `dense_backend` (a logical flag
indicating whether the dense solver backend is in use). The object is
used with [`base::solve()`](https://rdrr.io/r/base/solve.html) to
compute a solution and with
[`stats::update()`](https://rdrr.io/r/stats/update.html) to modify
problem data in place.
