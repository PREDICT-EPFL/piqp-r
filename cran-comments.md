## piqp 0.6.4

This release updates the vendored PIQP C++ library to v0.6.4, which removes
the use of the `EIGEN_NOEXCEPT` macro that no longer exists in Eigen 5. This
was requested by the RcppEigen maintainers (piqp-r issue #6) ahead of the
RcppEigen update that bundles Eigen 5. The package has been checked against
both the current CRAN RcppEigen (Eigen 3.4) and the development RcppEigen
(Eigen 5.0.1).

## R CMD check results

0 errors | 0 warnings | 0 notes

* Local: macOS arm64, R 4.6.1, `R CMD check --as-cran`
* GitHub Actions: macOS (release), Windows (release), Ubuntu (devel,
  release, oldrel-1)

## Reverse dependencies

There are no reverse dependencies on CRAN.
