#' @title Interface to 'piqp', a proximal interior point quadratic solver
#'
#' @description
#'
#' PIQP solves linear programs (LPs) and quadratic programs (QPs). ,
#'
#' Minimize \deqn{\frac{1}{2}x^TPx + c^Tx} subject to \deqn{Ax = b},
#'   \deqn{Gx \leq h}, \deqn{x_{lb} \leq x \leq x_{ub}}, where \eqn{x
#'   \in R^n}, \eqn{h \in R^m}, \eqn{P = P^T} and
#'   nonnegative-definite, \eqn{c \in R^n}, \eqn{A \in R^{p\times n}},
#'   \eqn{G \in R^{m\times n}}, \eqn{b \in R^p}, \eqn{x{lb} \in R^n},
#'   \eqn{x_{ub} \in R^n}.
#' @details PIQP uses compressed column sparse (CSC) format
#'   internally, which is the [Matrix::dgCMatrix-class] in the Matrix
#'   package. So even P should be in CSC format, not taking advantage
#'   of symmetry for instance. This package accepts either
#'   [Matrix::dgcMatrix-class] class objects or the triplet format sparse
#'   matrices as in the [slam::simple_triplet_matrix()]
#' 
#' @param A a matrix of constraint coefficients
#' @param b a numeric vector giving the primal constraints
#' @param c a numeric vector giving the primal objective
#' @param P a symmetric positive semidefinite matrix
#' @param G a matrix of inequality coefficients
#' @param h a numeric vector of the inequality right hand side
#' @param x_lb a numeric vector of lower bounds, default `NULL`
#'   indicating `-Inf` for all variables, otherwise should be number
#'   of variables long
#' @param x_lb a numeric vector of upper bounds, default `NULL`
#'   indicating `Inf` for all variables, otherwise should be number of
#'   variables long
#' @param control a list giving specific control parameters to use in
#'   place of default values, with empty list indicating defaults
#' @param interface whether to use sparse, or dense interface. Default
#'   is "auto" mode will switch to sparse interface if any of the
#'   specified matrices are sparse as in the Matrix package.
#' @return named list
#'   - `status` (status code)
#'   - `x` (primal solution)
#'   - `y` (dual solution of equality constraints)
#'   - `z` (dual solution of inequality constraints)
#'   - `z_lb` (dual solution of lower bound box constraints)
#'   - `z_ub` (dual solution of upper bound box constraints)
#'   - `pobj` (primal objective value)
#'   - `dobj` (dual objective value)
#'   - `dgap` (absolute duality gap absolute)
#'   - `dgap_rel` (relative duality gap)
#'   - `run_time` (total runtime, if so specified in settings)
#'
#' @importFrom Matrix sparseMatrix
#' @seealso [piqp_control()], [status_descriptions()]
#' @export piqp
#'
#' 
#  ---------------------------------------------------------
piqp <- function(A = matrix(0, nrow = 0, ncol = 0), b = numeric(0), c = numeric(0),
                 P = matrix(0, nrow = 0, ncol = 0), G = matrix(0, nrow = 0, ncol = 0),
                 h = numeric(0),
                 x_lb = NULL, x_ub = NULL,
                 control = list(),
                 interface = c("auto", "dense", "sparse")) {
  interface <- match.arg(interface)
  # Sanitize control parameters
  control <- do.call(piqp_control, control)

  ncol_A <- NCOL(A)
  nrow_A <- NROW(A)
  ncol_P <- NCOL(P)
  ncol_G <- NCOL(G)
  nrow_G <- NROW(G)
  
  n <- max(ncol_A, ncol_P, ncol_G)
  if (n <= 0) {
    stop("piqp: Bad input matrices")
  }
         
  if (ncol_A == 0) {
    A <- matrix(0, nrow = 0, ncol = n)
  }
  if (ncol_P == 0) {
    P <- matrix(0, nrow = n, ncol = n)
  }
  if (ncol_G == 0) {
    G <- matrix(0, nrow = 0, ncol = n)
  }
  
  if (is.null(x_lb)) x_lb <- rep(-Inf, n)
  if (is.null(x_ub)) x_ub <- rep(Inf, n)

  if (interface == "auto") {
    if (inherits(A, "dgCMatrix") || inherits(A, "simple_triplet_matrix") ||
          inherits(P, "dgCMatrix") || inherits(P, "simple_triplet_matrix") ||
          inherits(G, "dgCMatrix") || inherits(G, "simple_triplet_matrix")) {
      use_sparse <- TRUE
    } else {
      use_sparse <- FALSE
    }
  } else {
    use_sparse <- (interface == "sparse")
  }
      
  ## Checks needed.
  if (use_sparse) {
    if (!inherits(A, "dgCMatrix") ) {
      csc <- make_csc_matrix(A)
      Ai <- csc[["matind"]]
      Ap <- csc[["matbeg"]]
      Ax <- csc[["values"]]
      A <- Matrix::sparseMatrix(i = Ai, p = Ap, x = Ax, dims = c(nrow_A, n), index1 = FALSE)
    }

    if (!inherits(P, "dgCMatrix") ) {
      csc  <- make_csc_matrix(P)
      Pi <- csc[["matind"]]
      Pp <- csc[["matbeg"]]
      Px <- csc[["values"]]
      P <- Matrix::sparseMatrix(i = Pi, p = Pp, x = Px, dims = c(n, n), index1 = FALSE)      
    }

    if (!inherits(G, "dgCMatrix") ) {
      csc  <- make_csc_matrix(G)
      Gi <- csc[["matind"]]
      Gp <- csc[["matbeg"]]
      Gx <- csc[["values"]]
      G <- Matrix::sparseMatrix(i = Gi, p = Gp, x = Gx, dims = c(nrow_G, n), index1 = FALSE)      
    }
    ## piqp_sparse_solve(P, c, A, b, G, h, x_lb, x_ub, control)
    .Call('_piqp_piqp_sparse_solve', PACKAGE = 'piqp', P, c, A, b, G, h, x_lb, x_ub, control)
  } else {
    ##piqp_dense_solve(P, c, A, b, G, h, x_lb, x_ub, control)
    .Call('_piqp_piqp_dense_solve', PACKAGE = 'piqp', P, c, A, b, G, h, x_lb, x_ub, control)
  }
}

#' Control parameters with default values and types in parenthesis
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
#' @return a list containing the control parameters.
#' @export piqp_control
piqp_control <- function(
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

  if (any(sapply(params, length) != 1L)) stop("piqp_control: arguments should be scalars!")
  if (any(unlist(params[int_params]) < 0)) stop("piqp_control: integer arguments should be >= 0!")
 
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
#' @param a valid solver return code
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

