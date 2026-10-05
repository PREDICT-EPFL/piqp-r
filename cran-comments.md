## piqp 0.6.4.1

This is a maintenance release addressing the gcc-UBSAN issue reported by
CRAN for piqp 0.6.4
(https://www.stats.ox.ac.uk/pub/bdr/memtests/gcc-UBSAN/piqp/):

    Eigen/src/Cholesky/LLT.h:66:49: runtime error: load of value 31578,
    which is not a valid value for type 'ComputationInfo'

The dense backend's KKT constructor move-assigned a freshly constructed
`Eigen::LLT`, whose `m_info` member is left uninitialized by the LLT
constructors in released Eigen (both 3.4.x, as bundled in the current CRAN
RcppEigen, and 5.0.x, as bundled in the forthcoming RcppEigen). The object
is now initialized by factorizing a full-size zero matrix. This is the fix
merged into the upstream PIQP library (PREDICT-EPFL/piqp#45), applied to
the vendored sources through the package's patch file.

The fix was verified with a GCC 16 UBSan check on Fedora 44 (R-devel,
r-hub gcc16 container), matching CRAN's gcc-UBSAN setup. That check
reproduces the error above on 0.6.4 and is clean on 0.6.4.1. The r-hub
gcc-asan and clang-ubsan containers are also clean on 0.6.4.1.

This release also removes a `src/.r_patched` build marker that R-devel's
`R CMD build` would otherwise include in the tarball.

## R CMD check results

0 errors | 0 warnings | 1 note

* The NOTE is "Days since last update": this release follows 0.6.4
  closely in order to fix the UBSAN issue reported by CRAN above.

* Local: macOS arm64, R 4.6.1, `R CMD check --as-cran`
* GitHub Actions: macOS (release), Windows (release), Ubuntu (devel,
  release, oldrel-1)

## Reverse dependencies

There are no reverse dependencies on CRAN.
