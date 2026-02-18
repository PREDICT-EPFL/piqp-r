# PIQP Solver Interface

## 1. Introduction

PIQP solves quadratic programs of the form

$$\begin{aligned}
{\min\limits_{x}\quad} & {\frac{1}{2}x^{\top}Px + c^{\top}x} \\
{\text{s.t.}\quad} & {Ax = b,} \\
 & {h_{l} \leq Gx \leq h_{u},} \\
 & {x_{l} \leq x \leq x_{u}}
\end{aligned}$$

with primal decision variables $x \in {\mathbb{R}}^{n}$, matrices
$P \in {\mathbb{S}}_{+}^{n}$, $A \in {\mathbb{R}}^{p \times n}$,
$G \in {\mathbb{R}}^{m \times n}$, and vectors $c \in {\mathbb{R}}^{n}$,
$b \in {\mathbb{R}}^{p}$, $h_{l} \in {\mathbb{R}}^{m}$,
$h_{u} \in {\mathbb{R}}^{m}$, $x_{l} \in {\mathbb{R}}^{n}$, and
$x_{u} \in {\mathbb{R}}^{n}$.

## 2. The Problem Solver Interface

Consider:

$$\begin{aligned}
{\min\limits_{x}\quad} & {\frac{1}{2}x^{\top}\begin{bmatrix}
6 & 0 \\
0 & 4
\end{bmatrix}x + \begin{bmatrix}
{- 1} \\
{- 4}
\end{bmatrix}^{\top}x} \\
{\text{s.t.}\quad} & {\begin{bmatrix}
1 & {- 2}
\end{bmatrix}x = 1,} \\
 & {\begin{bmatrix}
1 & {- 1} \\
2 & 0
\end{bmatrix}x \leq \begin{bmatrix}
0.2 \\
{- 1}
\end{bmatrix},} \\
 & {- 1 \leq x_{1} \leq 1.}
\end{aligned}$$

The data for this problem can be specified as below.

``` r
P <- matrix(c(6, 0, 0, 4), nrow = 2)
c <- c(-1, -4)
A <- matrix(c(1, -2), nrow = 1)
b <- 1
G <- matrix(c(1, 2, -1, 0), nrow = 2)
h_u <- c(0.2, -1)
x_l <- c(-1, -Inf)  ## 2 variables
x_u <- c(1, Inf)    ## 2 variables
```

The problem can now be solved via a call to
[`solve_piqp()`](../reference/solve_piqp.md).

``` r
sol <- solve_piqp(P, c, A, b, G, h_u = h_u, x_l = x_l, x_u = x_u, backend = "auto")
cat(sprintf("(Solution status, description): = (%d, %s)\n",
            sol$status, sol$info$status_desc))
#> (Solution status, description): = (1, solved)
cat(sprintf("Objective: %f, solution: (x1, x2) = (%f, %f)\n", sol$info$primal_obj, sol$x[1], sol$x[2]))
#> Objective: 6.160000, solution: (x1, x2) = (-0.600000, -0.800000)
```

`sol` contains many components as `str(sol)` will display but the most
important ones are:

- `status` : 1 if all goes well (more below),
- `x` : solution vector
- `y` : dual solution for the equality constraints
- `z_l` : dual solution for the lower inequality constraints
- `z_u` : dual solution for the upper inequality constraints
- `z_bl` : dual solution of lower bound box constraints
- `z_bu` : dual solution of upper bound box constraints
- `info$status_desc`: a descriptive string of the status
- `info$primal_obj` : primal objective value
- `info$run_time` : total runtime, if asked for in settings (see below).

One can always construct the descriptive string for the status using:

``` r
status_description(sol$status)
#> [1] "Solver solved problem up to given tolerance."
```

Note that PIQP can handle infinite box constraints well, i.e. when
elements of $x_{l}$ or $x_{u}$ are $- \infty$ or $\infty$, respectively.
The inequality constraints now support double-sided bounds
$h_{l} \leq Gx \leq h_{u}$. For one-sided inequalities $Gx \leq h$,
simply pass `h_u = h` (the default `h_l` is `-Inf`).

## 3. The Solver Model Object

Users who wish to solve QP problems will mostly use the
[`solve_piqp()`](../reference/solve_piqp.md) function. Behind the
scenes, [`solve_piqp()`](../reference/solve_piqp.md) creates a solver
object and calls methods on the object to obtain the solution. The
solver object can be created explicitly using
[`piqp()`](../reference/piqp.md) and provides more elaborate facilities
for updating problem data and warm-starting subsequent solves, which can
be very efficient when one is solving a sequence of related problems.

The above problem could be solved using the solver model object thus:

``` r
model <- piqp(P, c, A, b, G, h_u = h_u, x_l = x_l, x_u = x_u)
sol2 <- solve(model)
identical(sol, sol2)
#> [1] FALSE
```

Indeed, this is exactly what
[`solve_piqp()`](../reference/solve_piqp.md) does. The generic functions
[`solve()`](https://rdrr.io/pkg/Matrix/man/solve-methods.html),
[`update()`](https://rdrr.io/r/stats/update.html),
[`get_settings()`](../reference/get_settings.md),
[`get_dims()`](../reference/get_dims.md), and
[`update_settings()`](../reference/update_settings.md) are available on
the model object.

``` r
get_dims(model)
#> $n
#> [1] 2
#> 
#> $p
#> [1] 1
#> 
#> $m
#> [1] 2
```

### Updating problem data and re-solving

The real advantage of the model object is the ability to update problem
data and re-solve. When
[`update()`](https://rdrr.io/r/stats/update.html) is called, only the
changed data is passed to the solver, and the solver can warm-start from
the previous solution. This is much more efficient than creating a new
model from scratch each time.

Suppose we want to tighten the inequality bounds and change the linear
cost. We update the model in place and re-solve:

``` r
update(model, c = c(-2, -3), h_u = c(0.1, -1.5))
sol3 <- solve(model)
cat(sprintf("Status: %s\n", sol3$info$status_desc))
#> Status: solved
cat(sprintf("Objective: %f, solution: (x1, x2) = (%f, %f)\n",
            sol3$info$primal_obj, sol3$x[1], sol3$x[2]))
#> Objective: 7.840000, solution: (x1, x2) = (-0.800000, -0.900000)
```

We can continue updating. Here we relax the variable bounds and change
the equality constraint:

``` r
update(model, A = matrix(c(1, -1), nrow = 1), b = 0,
       x_l = c(-2, -2), x_u = c(2, 2))
sol4 <- solve(model)
cat(sprintf("Status: %s\n", sol4$info$status_desc))
#> Status: solved
cat(sprintf("Objective: %f, solution: (x1, x2) = (%f, %f)\n",
            sol4$info$primal_obj, sol4$x[1], sol4$x[2]))
#> Objective: 6.562500, solution: (x1, x2) = (-0.750000, -0.750000)
```

Dimension mismatches are caught:

``` r
update(model, b = c(5, 2))
#> Error in `update.piqp::piqp_model`:
#> ! Update parameters do not match original problem dimensions
```

Settings can also be updated between solves:

``` r
update_settings(model, new_settings = list(verbose = FALSE, max_iter = 100L))
```

## 4. Dense and Sparse Interfaces

PIQP supports dense and sparse problem formulations. For small and dense
problems the dense interface is preferred since vectorized instructions
and cache locality can be exploited more efficiently, but for sparse
problems the sparse interface and result in significant speedups.

Either interface can be requested explicitly via the `backend` parameter
which can take on any value among `"dense"`, `"sparse"`, or `"auto"`,
the default. The last value will automatically switch to a sparse
interface if any of the supplied inputs ($A$, $P$, or $G$) is a sparse
matrix; otherwise it uses the dense interface.

``` r
sparse_sol <- solve_piqp(P, c, A, b, G, h_u = h_u, x_l = x_l, x_u = x_u, backend = "sparse")
str(sparse_sol)
#> List of 12
#>  $ status: int 1
#>  $ x     : num [1:2] -0.6 -0.8
#>  $ y     : num -11.8
#>  $ z_l   : num [1:2] 0 0
#>  $ z_u   : num [1:2] 1.64e+01 8.96e-11
#>  $ z_bl  : num [1:2] 2.18e-10 0.00
#>  $ z_bu  : num [1:2] 1.55e-12 0.00
#>  $ s_l   : num [1:2] 1e+30 1e+30
#>  $ s_u   : num [1:2] 5.67e-12 2.00e-01
#>  $ s_bl  : num [1:2] 4e-01 1e+30
#>  $ s_bu  : num [1:2] 1.6 1.0e+30
#>  $ info  :List of 32
#>   ..$ status_desc       : chr "solved"
#>   ..$ iter              : num 7
#>   ..$ rho               : num 1e-10
#>   ..$ delta             : num 1e-10
#>   ..$ mu                : num 5.02e-11
#>   ..$ sigma             : num 1e-06
#>   ..$ primal_step       : num 0.99
#>   ..$ dual_step         : num 0.99
#>   ..$ primal_res        : num 3.13e-13
#>   ..$ primal_res_rel    : num 1.96e-13
#>   ..$ dual_res          : num 1.5e-10
#>   ..$ dual_res_rel      : num 2.08e-11
#>   ..$ primal_res_reg    : num 3.13e-11
#>   ..$ primal_res_reg_rel: num 1.96e-11
#>   ..$ dual_res_reg      : num 1.5e-08
#>   ..$ dual_res_reg_rel  : num 2.08e-09
#>   ..$ primal_prox_inf   : num 0
#>   ..$ dual_prox_inf     : num 0
#>   ..$ primal_obj        : num 6.16
#>   ..$ dual_obj          : num 6.16
#>   ..$ duality_gap       : num 3.27e-10
#>   ..$ duality_gap_rel   : num 2.77e-11
#>   ..$ factor_retires    : num 0
#>   ..$ reg_limit         : num 1e-10
#>   ..$ no_primal_update  : num 1
#>   ..$ no_dual_update    : num 0
#>   ..$ setup_time        : num 0
#>   ..$ update_time       : num 0
#>   ..$ solve_time        : num 0
#>   ..$ kkt_factor_time   : num 0
#>   ..$ kkt_solve_time    : num 1479
#>   ..$ run_time          : num 0
```

## 5. Another Example

Suppose that we want to solve the following 2-dimensional quadratic
programming problem:

$$\begin{array}{ll}
\text{minimize} & {3x_{1}^{2} + 2x_{2}^{2} - x_{1} - 4x_{2}} \\
\text{subject to} & {- 1 \leq x \leq 1,\ x_{1} = 2x_{2}}
\end{array}$$

Since the solver expects the objective in the form
$\frac{1}{2}x^{\top}Px + c^{\top}x$, we define

$$P = 2 \cdot \begin{bmatrix}
3 & 0 \\
0 & 2
\end{bmatrix}{\mspace{6mu}\text{and}\mspace{6mu}}q = \begin{bmatrix}
{- 1} \\
{- 4}
\end{bmatrix}.$$

We have one equality constraint and box constraints. This leads to the
following straight-forward formulation.

``` r
P <- matrix(2 * c(3, 0, 0, 2), nrow = 2, ncol = 2)
c <- c(-1, -4)
A <- matrix(c(1, -2), ncol = 2)
b <- 0
x_l <- rep(-1.0, 2)
x_u <- rep(1.0, 2)
sol <- solve_piqp(P = P, c = c, A = A, b = b, x_l = x_l, x_u = x_u)
cat(sprintf("(Solution status, description): = (%d, %s)\n",
            sol$status, sol$info$status_desc))
#> (Solution status, description): = (1, solved)
cat(sprintf("Objective: %f, solution: (x1, x2) = (%f, %f)\n", sol$info$primal_obj, sol$x[1], sol$x[2]))
#> Objective: -0.642857, solution: (x1, x2) = (0.428571, 0.214286)
```

But we can also choose to move the upper box constraints into the
inequalities.

``` r
G <- diag(2)
h_u <- c(1, 1)
sol <- solve_piqp(P = P, c = c, A = A, b = b, G = G, h_u = h_u,
                  x_l = c(-1, -1), x_u = c(Inf, Inf))
cat(sprintf("(Solution status, description): = (%d, %s)\n",
            sol$status, sol$info$status_desc))
#> (Solution status, description): = (1, solved)
cat(sprintf("Objective: %f, solution: (x1, x2) = (%f, %f)\n", sol$info$primal_obj, sol$x[1], sol$x[2]))
#> Objective: -0.642857, solution: (x1, x2) = (0.428571, 0.214286)
```

Or we can move both of them into the inequalities.

``` r
G <- Matrix::Matrix(c(1, 0, -1, 0, 0, 1, 0, -1), byrow = TRUE,
                    nrow = 4, sparse = TRUE)
h_u <- rep(1, 4)

sol <- solve_piqp(A = A, b = b, c = c, P = P, G = G, h_u = h_u)
cat(sprintf("(Solution status, description): = (%d, %s)\n",
            sol$status, status_description(sol$status)))
#> (Solution status, description): = (1, Solver solved problem up to given tolerance.)
cat(sprintf("Objective: %f, solution: (x1, x2) = (%f, %f)\n", sol$info$primal_obj, sol$x[1], sol$x[2]))
#> Objective: -0.642857, solution: (x1, x2) = (0.428571, 0.214286)
```

All of them will yield the same result.

## 6. Solver parameters

PIQP has a number of parameters that control its behavior, including
verbosity, tolerances, etc.; see help on
[`piqp_settings()`](../reference/piqp_settings.md). As an example, in
the last problem, we can reduce the number of iterations.

``` r
s <- solve_piqp(P = P, c = c, A = A, b = b, G = G, h_u = h_u,
          settings = list(max_iter = 3)) ## Reduced number of iterations
cat(sprintf("(Solution status, description): = (%d, %s)\n",
            s$status, s$info$status_desc))
#> (Solution status, description): = (-1, max iterations reached)
cat(sprintf("Objective: %f, solution: (x1, x2) = (%f, %f)\n", s$info$primal_obj, s$x[1], s$x[2]))
#> Objective: -0.642857, solution: (x1, x2) = (0.428564, 0.214282)
```

Note the different status, which should always be checked in code.
