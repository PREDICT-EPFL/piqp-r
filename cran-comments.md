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
is now initialized via `compute()` on an empty matrix. The fix is applied to
the vendored PIQP sources through the package's patch file and has been
proposed to the upstream PIQP library.

The fix was verified locally by building the package with GCC 16 and
`-fsanitize=undefined -fsanitize-undefined-trap-on-error`: the previous
version traps in the dense KKT constructor, the fixed version runs clean.

## R CMD check results

0 errors | 0 warnings | 0 notes

* Local: macOS arm64, R 4.6.1, `R CMD check --as-cran`
* GitHub Actions: macOS (release), Windows (release), Ubuntu (devel,
  release, oldrel-1)

## Reverse dependencies

There are no reverse dependencies on CRAN.
