# This file is part of PIQP-R.
#
# Copyright (c) 2023 piqp-r authors
#
# This source code is licensed under the BSD 2-Clause License found in the
# LICENSE file in the root directory of this source tree.

#' PIQP Solver
#'
#' Solves \deqn{arg\min_x 0.5 x'P x + c'x}{argmin_x 0.5 x'P x + c'x}
#' s.t. \deqn{A x = b}{A x = b}
#'      \deqn{G x \leq h}{G x <= h}
#'      \deqn{x_{lb} \leq x \leq x_{ub}}{x_lb <= x <= x_ub}
#' for real matrices P (nxn, positive semidefinite), A (pxn) with m number of equality constraints, and G (mxn) with m number of inequality constraints
#' @param P dense or sparse matrix of class dgCMatrix or coercible into such, must be positive semidefinite
#' @param c numeric vector
#' @param A dense or sparse matrix of class dgCMatrix or coercible into such
#' @param b numeric vector
#' @param G dense or sparse matrix of class dgCMatrix or coercible into such
#' @param h numeric vector
#' @param x_lb a numeric vector of lower bounds, default `NULL`
#'   indicating `-Inf` for all variables, otherwise should be number
#'   of variables long
#' @param x_ub a numeric vector of upper bounds, default `NULL`
#'   indicating `Inf` for all variables, otherwise should be number of
#'   variables long
#' @param settings list with optimization parameters, empty by default; see [piqp_settings()] for a comprehensive list of parameters that may be used
#' @param backend which backend to use, if auto and P, A or G are sparse then sparse backend is used (`"auto"`, `"sparse"` or `"dense"`) (`"auto"`)
#' @return A list with elements solution elements
#' @seealso [piqp()], [piqp_settings()] and the underlying PIQP documentation: \url{https://predict-epfl.github.io/piqp/}
#' @examples
#' ## example, adapted from PIQP documentation
#' library(piqp)
#' library(Matrix)
#'
#' P <- Matrix(c(6., 0.,
#'               0., 4.), 2, 2, sparse = TRUE)
#' c <- c(-1., -4.)
#' A <- Matrix(c(1., -2.), 1, 2, sparse = TRUE)
#' b <- c(1.)
#' G <- Matrix(c(1., 2., -1., 0.), 2, 2, sparse = TRUE)
#' h <- c(0.2, -1.)
#' x_lb <- c(-1., -Inf)
#' x_ub <- c(1., Inf)
#'
#' settings <- list(verbose = TRUE)
#'
#' # Solve with PIQP
#' res <- solve_piqp(P, c, A, b, G, h, x_lb, x_ub, settings)
#' res$x
#'
#' @references{
#' Schwan, R., Jiang, Y., Kuhn, D., Jones, C.N. (2023).
#' ``PIQP: A Proximal Interior-Point Quadratic Programming Solver.''
#' <doi:10.48550/arXiv.2304.00290>
#' }
#' @export solve_piqp
solve_piqp <- function(P = NULL, c = NULL, A = NULL, b = NULL, G = NULL, h = NULL, x_lb = NULL, x_ub = NULL,
                       settings = list(), backend = c("auto", "sparse", "dense")) {
  backend <- match.arg(backend)
  model <- piqp(P, c, A, b, G, h, x_lb, x_ub, settings, backend)
  model$solve()
}

