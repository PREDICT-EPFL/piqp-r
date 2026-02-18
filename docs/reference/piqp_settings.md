# Settings parameters with default values and types in parenthesis

Settings parameters with default values and types in parenthesis

## Usage

``` r
piqp_settings(
  rho_init = 1e-06,
  delta_init = 1e-04,
  eps_abs = 1e-08,
  eps_rel = 1e-09,
  check_duality_gap = TRUE,
  eps_duality_gap_abs = 1e-08,
  eps_duality_gap_rel = 1e-09,
  infeasibility_threshold = 0.9,
  reg_lower_limit = 1e-10,
  reg_finetune_lower_limit = 1e-13,
  reg_finetune_primal_update_threshold = 7L,
  reg_finetune_dual_update_threshold = 7L,
  max_iter = 250L,
  max_factor_retires = 10L,
  preconditioner_scale_cost = FALSE,
  preconditioner_reuse_on_update = FALSE,
  preconditioner_iter = 10L,
  tau = 0.99,
  iterative_refinement_always_enabled = FALSE,
  iterative_refinement_eps_abs = 1e-12,
  iterative_refinement_eps_rel = 1e-12,
  iterative_refinement_max_iter = 10L,
  iterative_refinement_min_improvement_rate = 5,
  iterative_refinement_static_regularization_eps = 1e-08,
  iterative_refinement_static_regularization_rel = .Machine$double.eps^2,
  verbose = FALSE,
  compute_timings = FALSE
)
```

## Arguments

- rho_init:

  Initial value for the primal proximal penalty parameter rho (default =
  1e-6)

- delta_init:

  Initial value for the augmented lagrangian penalty parameter delta
  (default = 1e-4)

- eps_abs:

  Absolute tolerance (default = 1e-8)

- eps_rel:

  Relative tolerance (default = 1e-9)

- check_duality_gap:

  Check terminal criterion on duality gap (default = TRUE)

- eps_duality_gap_abs:

  Absolute tolerance on duality gap (default = 1e-8)

- eps_duality_gap_rel:

  Relative tolerance on duality gap (default = 1e-9)

- infeasibility_threshold:

  Threshold for infeasibility detection (default = 0.9)

- reg_lower_limit:

  Lower limit for regularization (default = 1e-10)

- reg_finetune_lower_limit:

  Fine tune lower limit regularization (default = 1e-13)

- reg_finetune_primal_update_threshold:

  Threshold of number of no primal updates to transition to fine tune
  mode (default = 7)

- reg_finetune_dual_update_threshold:

  Threshold of number of no dual updates to transition to fine tune mode
  (default = 7)

- max_iter:

  Maximum number of iterations (default = 250)

- max_factor_retires:

  Maximum number of factorization retires before failure (default = 10)

- preconditioner_scale_cost:

  Scale cost in Ruiz preconditioner (default = FALSE)

- preconditioner_reuse_on_update:

  Reuse preconditioner on problem update (default = FALSE)

- preconditioner_iter:

  Maximum of preconditioner iterations (default = 10)

- tau:

  Maximum interior point step length (default = 0.99)

- iterative_refinement_always_enabled:

  Always run iterative refinement and not only on factorization failure
  (default = FALSE)

- iterative_refinement_eps_abs:

  Iterative refinement absolute tolerance (default = 1e-12)

- iterative_refinement_eps_rel:

  Iterative refinement relative tolerance (default = 1e-12)

- iterative_refinement_max_iter:

  Maximum number of iterations for iterative refinement (default = 10)

- iterative_refinement_min_improvement_rate:

  Minimum improvement rate for iterative refinement (default = 5.0)

- iterative_refinement_static_regularization_eps:

  Static regularization for KKT system for iterative refinement (default
  = 1e-8)

- iterative_refinement_static_regularization_rel:

  Static regularization w.r.t. the maximum abs diagonal term of KKT
  system. (default = .Machine\$double.eps^2)

- verbose:

  Verbose printing (default = FALSE)

- compute_timings:

  Measure timing information internally (default = FALSE)

## Value

a list containing the settings parameters.
