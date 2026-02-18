# R Interface to PIQP Solver

PIQP is an Proximal Interior Point Quadratic Programming solver, which
can solve dense and sparse quadratic programs described in described in
Schwan, Jiang, Kuhn, and Jones (2023)
(<https://arxiv.org/abs/2304.00290>). Combining an infeasible interior
point method with the proximal method of multipliers, the algorithm can
handle ill-conditioned convex QP problems without the need for linear
independence of the constraints. The solver is written in header only
'C++ 14' leveraging the Eigen library for vectorized linear algebra. For
small dense problems, vectorized instructions and cache locality can be
exploited more efficiently. Allocation free problem updates and
re-solves are also provided.

## See also

Useful links:

- <https://predict-epfl.github.io/piqp-r/>

- Report bugs at <https://github.com/PREDICT-EPFL/piqp-r/issues>

## Author

Balasubramanian Narasimhan, Roland Schwan (C), Yuning Jiang, Daniel
Kuhn, Colin N. Jones
