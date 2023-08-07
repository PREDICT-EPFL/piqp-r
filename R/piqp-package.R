#' R Interface to PIQP Solver
#'
#' @description PIQP is an embedded Proximal Interior Point Quadratic
#'   Programming solver, which can solve dense and sparse quadratic
#'   programs <https://doi.org/10.48550/arXiv.2304.00290>. Combining an infeasible
#'   interior point method with the proximal method of multipliers,
#'   the algorithm can handle ill-conditioned convex QP problems
#'   without the need for linear independence of the constraints. The
#'   solver is written in header only C++ 14 leveraging the Eigen
#'   library for vectorized linear algebra. For small dense problems,
#'   vectorized instructions and cache locality can be exploited more
#'   efficiently. Allocation free problem updates and re-solves are
#'   also provided.
#'
#' @name piqp-package
#' @docType package
#' @useDynLib piqp 
#' @importFrom Rcpp evalCpp
#' @author Balasubramanian Narasimhan, Roland Schwan (C), Yuning Jiang, Daniel Kuhn, Colin N. Jones
#' @keywords package
NULL
