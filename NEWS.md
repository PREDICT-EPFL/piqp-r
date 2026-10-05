# piqp 0.6.4.1

* Fix undefined behavior reported by CRAN's gcc-UBSAN checks: the dense
  backend's KKT constructor move-assigned a freshly constructed `Eigen::LLT`,
  whose `m_info` member released Eigen (3.4.x and 5.0.x) leaves
  uninitialized. The Cholesky object is now initialized by factorizing a
  full-size zero matrix, which also keeps all memory allocation in setup.
  This is the fix merged upstream in PIQP (PREDICT-EPFL/piqp#45), applied
  to the vendored v0.6.4 sources through the package's patch file.
* `R CMD build` with R-devel no longer leaves the `src/.r_patched` marker in
  the tarball: a `clean` rule in `src/Makevars` now removes it.
* New on-demand GitHub Actions workflow `sanitizers.yaml` that checks the
  package in the r-hub gcc16, gcc-asan, and clang-ubsan containers.

# piqp 0.6.4

* Update to v0.6.4 of the underlying PIQP library
* Compatibility with Eigen 5 (replaces removed `EIGEN_NOEXCEPT` macro),
  required for upcoming RcppEigen release (#6)

# piqp 0.6.2

* Update to v0.6.2 of the underlying PIQP library
* Migrate from R6 to S7 OOP system
* Inequality constraints are now double-sided: `h_l <= Gx <= h_u`
* Variable bounds renamed: `x_lb`/`x_ub` to `x_l`/`x_u`
* Result fields renamed: `z` to `z_l`/`z_u`, `z_lb`/`z_ub` to `z_bl`/`z_bu`,
  `s` to `s_l`/`s_u`, `s_lb`/`s_ub` to `s_bl`/`s_bu`
* Info fields renamed: `primal_inf` to `primal_res`, `dual_inf` to `dual_res`
* New settings: `infeasibility_threshold`, `preconditioner_reuse_on_update`
* Use generic functions `solve()`, `update()`, `get_settings()`,
  `get_dims()`, and `update_settings()` instead of R6 methods
* Support for problem data updates and warm starts via `update()`

# piqp 0.3.1

* Update to v0.3.1 of the underlying PIQP library

# piqp 0.2.2

* First CRAN release

