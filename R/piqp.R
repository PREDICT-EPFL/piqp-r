# This file is part of PIQP-R. It is partly based on osqp-r
# (https://github.com/osqp/osqp-r) which is licensed
# under Apache License 2.0
#
# Copyright (c) 2023 piqp-r authors
# Copyright (c) 2019 Paul Goulart, Bartolomeo Stellato
#
# This source code is licensed under the BSD 2-Clause License found in the
# LICENSE file in the root directory of this source tree.

#' @title The PIQP Solver Model Class
#'
#' @description This class wraps around the PIQP C++ Solver and
#'   exposes methods and fields of the C++ object. Users will never
#'   need to directly create instances this class and should use the
#'   more user-friendly functions [piqp()] and [solve_piqp()].
#'
#' @importFrom Rcpp evalCpp
#' @importFrom R6 R6Class
piqp_model <-
  R6::R6Class("piqp_model",
              private = list(.solver_ptr = NULL, .dims = NULL, .dense_backend = NULL),          
              public =
                list(
                  #' @description
                  #' Create a new piqp_model object
                  #' @param P dense or sparse matrix of class dgCMatrix or coercible into such, must be positive semidefinite
                  #' @param c numeric vector
                  #' @param A dense or sparse matrix of class dgCMatrix or coercible into such
                  #' @param b numeric vector
                  #' @param G dense or sparse matrix of class dgCMatrix or coercible into such
                  #' @param h numeric vector
                  #' @param x_lb a numeric vector of lower bounds
                  #' @param x_ub a numeric vector of upper bounds
                  #' @param settings list with optimization parameters
                  #' @param dense_backend a flag indicating if the dense solver is to be used
                  #' @param dims the dimensions of the problem, a named list containing `n`, `p` and `m`.
                  #' @return a piqp_model object that can be used to solve the QP
                  initialize = function(P, c, A, b, G, h, x_lb, x_ub, settings = list(), dense_backend, dims) {
                    private$.dense_backend <- dense_backend
                    private$.dims <- dims
                    if (dense_backend) {
                      private$.solver_ptr <-.Call('_piqp_piqp_dense_setup', PACKAGE = 'piqp', P, c, A, b, G, h, x_lb, x_ub, settings)
                    } else {
                      private$.solver_ptr <- .Call('_piqp_piqp_sparse_setup', PACKAGE = 'piqp', P, c, A, b, G, h, x_lb, x_ub, settings)
                    }
                  }
                 ,
                  #' @description
                  #' Solve the QP model
                  #' @return a list containing the solution
                  solve = function() {
                    .Call('_piqp_solve_model', PACKAGE = 'piqp', private$.solver_ptr, private$.dense_backend)
                  }
                 ,
                  #' @description
                  #' Update the current piqp_model with new data
                  #' @param P dense or sparse matrix of class dgCMatrix or coercible into such, must be positive semidefinite
                  #' @param c numeric vector
                  #' @param A dense or sparse matrix of class dgCMatrix or coercible into such
                  #' @param b numeric vector
                  #' @param G dense or sparse matrix of class dgCMatrix or coercible into such
                  #' @param h numeric vector
                  #' @param x_lb a numeric vector of lower bounds
                  #' @param x_ub a numeric vector of upper bounds
                  #' @param settings list with optimization parameters
                  #' @param dense_backend a flag indicating if the dense solver is to be used
                  #' @param dims the dimensions of the problem, a named list containing `n`, `p` and `m`.
                  update = function(P = NULL, c = NULL, A = NULL, b = NULL, G = NULL, h = NULL,
                                    x_lb = NULL, x_ub = NULL) {
                    dims <- private$.dims
                    n <- dims$n; p <- dims$p; m <- dims$m;
                    if (! (length(c) %in% c(0, n) && 
                             length(b) %in% c(0, p) &&
                             length(h) %in% c(0, m) &&
                             length(x_lb) %in% c(0, n) &&
                             length(x_ub) %in% c(0, n)) ) {
                      stop("Update parameters not match original problem dimensions")
                    }
                    
                    sparse_backend <- !private$.dense_backend
                    
                    if (!is.null(P)) {
                      if (NCOL(P) != n || NROW(P) != n) {
                        stop("P dimension not matching original problem")
                      }
                      if (sparse_backend) P <- ensure_dgc_matrix(P)
                    }
                    
                    if (!is.null(A)) {
                      if (NCOL(A) != n || NROW(A) != p) {
                        stop("A dimension not matching original problem")
                      }
                      if (sparse_backend) A <- ensure_dgc_matrix(A)                      
                    }
                    
                    if (!is.null(G)) {
                      if (NCOL(G) != n || NROW(G) != m) {
                        stop("G dimension not matching original problem")
                      }
                      if (sparse_backend) G <- ensure_dgc_matrix(G)                                            
                    }

                   if (sparse_backend) {
                      invisible(.Call('_piqp_piqp_update_sparse', PACKAGE = 'piqp', private$.solver_ptr, P, c, A, b, G, h, x_lb, x_ub))
                    } else {
                      invisible(.Call('_piqp_piqp_update_dense', PACKAGE = 'piqp', private$.solver_ptr, P, c, A, b, G, h, x_lb, x_ub))
                    }
                  }
                 ,
                  #' @description
                  #' Obtain the current settings for this model
                  get_settings = function() {
                    .Call('_piqp_get_settings', PACKAGE = 'piqp', private$.solver_ptr, private$.dense_backend)
                  }
                 ,
                  #' @description
                  #' Obtain the dimensions of this model
                  get_dims = function() {
                    private$.dims
                  }
                 ,
                  #' @description
                  #' Update the current settings with new values for this model
                  #' @param new_settings a list of named values for settings, default empty list; see [piqp_settings()] for a comprehensive list of defaults
                  update_settings = function(new_settings = list()) {
                    invisible(.Call('_piqp_update_settings', PACKAGE = 'piqp', private$.solver_ptr, private$.dense_backend, settings))
                  }
                )
              )

#' PIQP Solver object
#'
#' @importFrom Matrix sparseMatrix
#' @inheritParams solve_piqp
#' @return An R6-object of class "piqp_model" with methods defined which can be further
#'   used to solve the problem with updated settings / parameters.
#' @seealso [solve_piqp()],  [piqp_settings()]
#' @section Usage:
#' \preformatted{model = piqp(P = NULL, c = NULL, A = NULL, b = NULL, G = NULL, h = NULL, x_lb = NULL, x_ub = NULL, settings = piqp_settings(), backend = c("auto", "sparse", "dense"))
#'
#' model$solve()
#' model$update(P = NULL, c = NULL, A = NULL, b = NULL, G = NULL, h = NULL, x_lb = NULL, x_ub = NULL)
#' model$get_settings()
#' model$get_dims()
#' model$update_settings(new_settings = piqp_settings())
#'
#' print(model)
#' }
#' @details
#' Allows one to solve a parametric
#' problem with for example warm starts between updates of the parameter, c.f. the examples.
#' The object returned by \code{piqp} contains several methods which can be used to either update/get details of the
#' problem, modify the optimization settings or attempt to solve the problem.
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
#' model <- piqp(P, c, A, b, G, h, x_lb, x_ub, settings)
#'
#' # Solve
#' res <- model$solve()
#' res$x
#'
#' # Define new data
#' A_new <- Matrix(c(1., -3.), 1, 2, sparse = TRUE)
#' h_new <- c(2., 1.)
#'
#' # Update model and solve again
#' model$update(A = A_new, h = h_new)
#' res <- model$solve()
#' res$x
#'
#' @export piqp
piqp <- function(P = NULL, c = NULL, A = NULL, b = NULL, G = NULL, h = NULL, x_lb = NULL, x_ub = NULL,
                 settings = list(), backend = c("auto", "sparse", "dense")) {

  # match possible options
  backend <- match.arg(backend)

  sparse_backend <- (backend == "sparse") || inherits(P, "simple_triplet_matrix") ||
    inherits(A, "simple_triplet_matrix") || inherits(G, "simple_triplet_matrix") ||
    inherits(P, "sparseMatrix") || inherits(A, "sparseMatrix") || inherits(G, "sparseMatrix")

  if (is.null(P)) {
    n <- length(c)
  } else {
    n <- NCOL(P)
  }

  if (n == 0) {
    stop("At least one of P and c must be supplied")
  }
  
  if (!sparse_backend) { ## dense

    if (is.null(P)) {
      P <- matrix(0, n, n)
    }

    if (is.null(A)) {
      p <- 0
      A <- matrix(0, 0, n)
      b <- numeric(0)
    } else {
      p <- nrow(A)
      if (length(b) != p)
        stop(sprintf("b length %d must match A number of rows %d", length(b), p))        
      if (NCOL(A) != n)
        stop(sprintf("A should have %d columns", n))      
    }
    
    if (is.null(G)) {
      m <- 0
      G <- matrix(0, 0, n)
      h <- numeric(0)
    } else {
      m <- nrow(G)
      if (length(h) != m)
        stop(sprintf("h length %d must match G number of rows %d", length(h), m))
      if (NCOL(G) != n)
        stop(sprintf("G should have %d columns", n))        
    }
    
  } else { ## sparse

    if (is.null(P)) {
      P <- Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0), dims = c(n, n))
    } else {
      P <- ensure_dgc_matrix(P, n, n)                      
    }

    if (is.null(A)) {
      p <- 0
      A <- Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0), dims = c(0, n))
      b <- numeric(0)
    } else {
      p <- NROW(A)
      if (length(b) != p) stop(sprintf("b length %d must match A number of rows %d", length(b), p))
      if (NCOL(A) != n) stop(sprintf("A should have %d columns", n))
      A <- ensure_dgc_matrix(A, p, n)
    }

    if (is.null(G)) {
      m <- 0
      G <- Matrix::sparseMatrix(i = integer(), j = integer(0), x = numeric(0), dims = c(0, n))
      h <- numeric(0)
    } else {
      m <- NROW(G)
      if (length(h) != m) stop(sprintf("h length %d must match G number of rows %d", length(h), m))
      if (NCOL(G) != n) stop(sprintf("G should have %d columns", n))        
      G <- ensure_dgc_matrix(G, m, n)
    }
  }
  if (is.null(c)) {
    c <- numeric(n)
  } else {
    if (length(c) != n) {
      stop(sprintf("P and c must have the same dimension %d", n))
    }
  }

  if (is.null(x_lb)) {
    x_lb <- rep(-Inf, n)
  } else {
    if (length(x_lb) != n) stop(sprintf("x_lb length should be %d", n))
  }
  
  if (is.null(x_ub)) {
    x_ub <- rep(Inf, n)
  } else {
    if (length(x_ub) != n) stop(sprintf("x_ub length should be %d", n))
  }
  
    piqp_model$new(P, c, A, b, G, h, x_lb, x_ub, settings, dense_backend = !sparse_backend,
                   dims = list(n = n, p = p, m = m))
}

#' @method format pipq_model
format.piqp_model <- function(x, ...) {
  dims <- x$get_dims()
  sprintf("PIQP-model object\n\nNumber of variables: %i\nNumber of constraints: %i", dims$n, dims$p)
}

#' @method print pipq_model
print.piqp_model <- function(x, ...)
  cat(format(x))


#' Settings parameters with default values and types in parenthesis
#'
#' @param rho_init Initial value for the primal proximal penalty parameter rho (default = 1e-6)
#' @param delta_init Initial value for the augmented lagrangian penalty parameter delta (default = 1e-4)
#' @param eps_abs Absolute tolerance (default = 1e-8)
#' @param eps_rel Relative tolerance (default = 1e-9)
#' @param check_duality_gap Check terminal criterion on duality gap (default = TRUE)
#' @param eps_duality_gap_abs Absolute tolerance on duality gap (default = 1e-8)
#' @param eps_duality_gap_rel Relative tolerance on duality gap (default = 1e-9)
#' @param reg_lower_limit Lower limit for regularization (default = 1e-10)
#' @param reg_finetune_lower_limit Fine tune lower limit regularization (default = 1e-13)
#' @param reg_finetune_primal_update_threshold Threshold of number of no primal updates to transition to fine tune mode (default = 7)
#' @param reg_finetune_dual_update_threshold Threshold of number of no dual updates to transition to fine tune mode (default = 5)
#' @param max_iter Maximum number of iterations (default = 250)
#' @param max_factor_retires Maximum number of factorization retires before failure (default = 10)
#' @param preconditioner_scale_cost Scale cost in Ruiz preconditioner (default = FALSE)
#' @param preconditioner_iter Maximum of preconditioner iterations (default = 10)
#' @param tau Maximum interior point step length (default = 0.99)
#' @param iterative_refinement_always_enabled Always run iterative refinement and not only on factorization failure (default = FALSE)
#' @param iterative_refinement_eps_abs Iterative refinement absolute tolerance (default = 1e-12)
#' @param iterative_refinement_eps_rel Iterative refinement relative tolerance (default = 1e-12)
#' @param iterative_refinement_max_iter Maximum number of iterations for iterative refinement (default = 10)
#' @param iterative_refinement_min_improvement_rate Minimum improvement rate for iterative refinement (default = 5.0)
#' @param iterative_refinement_static_regularization_eps Static regularization for KKT system for iterative refinement (default = 1e-7)
#' @param iterative_refinement_static_regularization_rel Static regularization w.r.t. the maximum abs diagonal term of KKT system. (default = .Machine$double.eps^2)
#' @param verbose Verbose printing (default = FALSE)
#' @param compute_timings Measure timing information internally (default = FALSE)
#' @return a list containing the settings parameters.
#' @export piqp_settings
piqp_settings <- function(
                          ## Main algorithm settings
                          rho_init = 1e-6,
                          delta_init = 1e-4,
                          eps_abs = 1e-8,
                          eps_rel = 1e-9,
                          check_duality_gap = TRUE,
                          eps_duality_gap_abs = 1e-8,
                          eps_duality_gap_rel = 1e-9,
                          reg_lower_limit = 1e-10,
                          reg_finetune_lower_limit = 1e-13,
                          reg_finetune_primal_update_threshold = 7L,
                          reg_finetune_dual_update_threshold = 5L,
                          max_iter = 250L,
                          max_factor_retires = 10L,
                          preconditioner_scale_cost = FALSE,
                          preconditioner_iter = 10L,
                          tau = 0.99,
                          iterative_refinement_always_enabled = FALSE,
                          iterative_refinement_eps_abs = 1e-12,
                          iterative_refinement_eps_rel = 1e-12,
                          iterative_refinement_max_iter = 10L,
                          iterative_refinement_min_improvement_rate = 5.0,
                          iterative_refinement_static_regularization_eps = 1e-7,
                          iterative_refinement_static_regularization_rel = .Machine$double.eps^2,
                          verbose = FALSE,
                          compute_timings = FALSE) {

  params <- as.list(environment())
  
  bool_params <- c("check_duality_gap",
                   "preconditioner_scale_cost",
                   "iterative_refinement_always_enabled",
                   "verbose",
                   "compute_timings")
  
  int_params <- c("reg_finetune_primal_update_threshold",
                  "reg_finetune_dual_update_threshold",
                  "max_iter",
                  "max_factor_retires",
                  "preconditioner_iter",
                  "iterative_refinement_max_iter")

  if (any(sapply(params, length) != 1L)) stop("piqp_settings: arguments should be scalars!")
  if (any(unlist(params[int_params]) < 0)) stop("piqp_settings: integer arguments should be >= 0!")
  
  ## The rest
  float_params <- setdiff(names(params), c(bool_params, int_params))

  for (x in bool_params) {
    params[[x]] <- as.logical(params[[x]])
  }
  for (x in int_params) {
    params[[x]] <- as.integer(params[[x]])
  }

  for (x in float_params) {
    params[[x]] <- as.numeric(params[[x]])
  }
  params
}

#' Return the solver status description string
#' @param code a valid solver return code
#' @return a status description string
#' @examples
#' status_description(1) ## for solved problem
#' status_description(-1) ## for max iterations limit reached
#' @export status_description
status_description <- function(code) {
  ## Solver descriptions (using return codes in reverse order and offset by 11)
  ## https://predict-epfl.github.io/piqp/interfaces/c_cpp/getting_started_cpp
  desc <- c("Invalid settings were provided to the solver.", #(-10)
            "The problem is unsolved, i.e., `solve` was never called.", #(-9) 
            "Numerical error occurred during solving.", #(-8),
            NA,
            NA,
            NA,
            NA,
            "The problem is dual infeasible.", #(-3) 
            "The problem is primal infeasible.", #(-2) 
            "Iteration limit was reached.",#(-1),
            NA,
            "Solver solved problem up to given tolerance." #(1)
            )
  desc[code + 11]
}
