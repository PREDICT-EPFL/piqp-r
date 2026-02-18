# Changelog

## piqp 0.6.2

- Update to v0.6.2 of the underlying PIQP library
- Migrate from R6 to S7 OOP system
- Inequality constraints are now double-sided: `h_l <= Gx <= h_u`
- Variable bounds renamed: `x_lb`/`x_ub` to `x_l`/`x_u`
- Result fields renamed: `z` to `z_l`/`z_u`, `z_lb`/`z_ub` to
  `z_bl`/`z_bu`, `s` to `s_l`/`s_u`, `s_lb`/`s_ub` to `s_bl`/`s_bu`
- Info fields renamed: `primal_inf` to `primal_res`, `dual_inf` to
  `dual_res`
- New settings: `infeasibility_threshold`,
  `preconditioner_reuse_on_update`
- Use generic functions
  [`solve()`](https://rdrr.io/pkg/Matrix/man/solve-methods.html),
  [`update()`](https://rdrr.io/r/stats/update.html),
  [`get_settings()`](../reference/get_settings.md),
  [`get_dims()`](../reference/get_dims.md), and
  [`update_settings()`](../reference/update_settings.md) instead of R6
  methods
- Support for problem data updates and warm starts via
  [`update()`](https://rdrr.io/r/stats/update.html)

## piqp 0.3.1

- Update to v0.3.1 of the underlying PIQP library

## piqp 0.2.2

CRAN release: 2023-08-14

- First CRAN release
