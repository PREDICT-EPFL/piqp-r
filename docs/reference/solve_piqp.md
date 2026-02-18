# PIQP Solver

Solves \$\$arg\min_x 0.5 x'P x + c'x\$\$ s.t. \$\$A x = b\$\$ \$\$h_l
\leq G x \leq h_u\$\$ \$\$x_l \leq x \leq x_u\$\$ for real matrices P
(nxn, positive semidefinite), A (pxn) with p number of equality
constraints, and G (mxn) with m number of inequality constraints

## Usage

``` r
solve_piqp(
  P = NULL,
  c = NULL,
  A = NULL,
  b = NULL,
  G = NULL,
  h_l = NULL,
  h_u = NULL,
  x_l = NULL,
  x_u = NULL,
  settings = list(),
  backend = c("auto", "sparse", "dense")
)
```

## Arguments

- P:

  dense or sparse matrix of class dgCMatrix or coercible into such, must
  be positive semidefinite

- c:

  numeric vector

- A:

  dense or sparse matrix of class dgCMatrix or coercible into such

- b:

  numeric vector

- G:

  dense or sparse matrix of class dgCMatrix or coercible into such

- h_l:

  numeric vector of lower inequality bounds, default `NULL` indicating
  `-Inf` for all inequality constraints

- h_u:

  numeric vector of upper inequality bounds, default `NULL` indicating
  `Inf` for all inequality constraints

- x_l:

  a numeric vector of lower variable bounds, default `NULL` indicating
  `-Inf` for all variables

- x_u:

  a numeric vector of upper variable bounds, default `NULL` indicating
  `Inf` for all variables

- settings:

  list with optimization parameters, empty by default; see
  [`piqp_settings()`](piqp_settings.md) for a comprehensive list of
  parameters that may be used

- backend:

  which backend to use, if auto and P, A or G are sparse then sparse
  backend is used (`"auto"`, `"sparse"` or `"dense"`) (`"auto"`)

## Value

A list with elements solution elements

## References

Schwan, R., Jiang, Y., Kuhn, D., Jones, C.N. (2023). “PIQP: A Proximal
Interior-Point Quadratic Programming Solver.”
<doi:10.48550/arXiv.2304.00290>

## See also

[`piqp()`](piqp.md), [`piqp_settings()`](piqp_settings.md) and the
underlying PIQP documentation: <https://predict-epfl.github.io/piqp/>

## Examples

``` r
## example, adapted from PIQP documentation
library(piqp)
library(Matrix)

P <- Matrix(c(6., 0.,
              0., 4.), 2, 2, sparse = TRUE)
c <- c(-1., -4.)
A <- Matrix(c(1., -2.), 1, 2, sparse = TRUE)
b <- c(1.)
G <- Matrix(c(1., 2., -1., 0.), 2, 2, sparse = TRUE)
h_u <- c(0.2, -1.)
x_l <- c(-1., -Inf)
x_u <- c(1., Inf)

settings <- list(verbose = TRUE)

# Solve with PIQP
res <- solve_piqp(P, c, A, b, G, h_u = h_u, x_l = x_l, x_u = x_u, settings = settings)
#> ----------------------------------------------------------
#>                         PIQP v0.6.2                       
#>                     (c) Roland Schwan                     
#>    Ecole Polytechnique Federale de Lausanne (EPFL) 2025   
#> ----------------------------------------------------------
#> sparse backend (sparse_ldlt)
#> variables n = 2, nzz(P upper triangular) = 2
#> equality constraints p = 1, nnz(A) = 2
#> inequality constraints m = 2, nnz(G) = 3
#> inequality lower bounds n_h_l = 0
#> inequality upper bounds n_h_u = 2
#> variable lower bounds n_x_l = 1
#> variable upper bounds n_x_u = 1
#> 
#> iter  prim_obj       dual_obj       duality_gap   prim_res      dual_res      rho         delta       mu          p_step   d_step
#>   0    2.43682e+00    3.58907e+00   1.15226e+00   1.32082e+00   1.70871e+00   1.000e-06   1.000e-04   5.720e-01   0.0000   0.0000
#>   1    5.07849e+00    5.73194e+00   6.53446e-01   1.70329e-01   4.06583e-01   1.742e-07   1.742e-05   9.965e-02   0.8711   0.9900
#>   2    7.11263e+00    5.32500e+00   1.78763e+00   1.57002e-03   8.11342e-01   1.742e-07   1.742e-05   2.475e-01   0.9900   0.2596
#>   3    6.16885e+00    6.13621e+00   3.26450e-02   3.14042e-05   1.49577e-02   3.525e-09   3.525e-07   5.007e-03   0.9798   0.9900
#>   4    6.16009e+00    6.15976e+00   3.26769e-04   3.12909e-07   1.49770e-04   1.000e-10   3.532e-09   5.017e-05   0.9900   0.9900
#>   5    6.16000e+00    6.16000e+00   3.26788e-06   3.12888e-09   1.49770e-06   1.000e-10   1.000e-10   5.018e-07   0.9900   0.9900
#>   6    6.16000e+00    6.16000e+00   3.26813e-08   3.12887e-11   1.49772e-08   1.000e-10   1.000e-10   5.018e-09   0.9900   0.9900
#>   7    6.16000e+00    6.16000e+00   3.26830e-10   3.12894e-13   1.49767e-10   1.000e-10   1.000e-10   5.019e-11   0.9900   0.9900
#> 
#> status:               solved
#> number of iterations: 7
#> objective:            6.16000e+00
res$x
#> [1] -0.6 -0.8
```
