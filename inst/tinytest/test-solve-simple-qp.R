# This file is part of PIQP-R.
#
# Copyright (c) 2023 piqp-r authors
#
# This source code is licensed under the BSD 2-Clause License found in the
# LICENSE file in the root directory of this source tree.

TOL <- 1e-06

define_simple_qp <- function(backend = "sparse") {
    P <- Matrix(c(6., 0.,
                  0., 4.), 2, 2, sparse = TRUE)
    c <- c(-1., -4.)
    A <- Matrix(c(1., -2.), 1, 2, sparse = TRUE)
    b <- c(0.)
    G <- Matrix(c(1., -1., 0., 0.), 2, 2, sparse = TRUE)
    h <- c(1., 1.)
    x_lb <- c(-Inf, -1.)
    x_ub <- c(Inf, 1.)

    settings <- list(verbose = FALSE)

    # Create PIQP model
    if (backend == "sparse") {
        model <- piqp(P, c, A, b, G, h, x_lb, x_ub, settings)
    } else {
        model <- piqp(as.matrix(P), c, as.matrix(A), b, as.matrix(G), h, x_lb, x_ub, settings)
    }

    return(model)
}

# Solve dense QP

model <- define_simple_qp(backend = "dense")

# Solve
res <- model$solve()

expect_equal(res$x, c(0.4285714, 0.2142857), tolerance = TOL)
expect_equal(res$y, c(-1.5714286), 1e-05)
expect_equal(res$z, c(0., 0.), tolerance = TOL)
expect_equal(res$z_lb, c(0., 0.), tolerance = TOL)
expect_equal(res$z_ub, c(0., 0.), tolerance = TOL)

# Solve sparse QP

model <- define_simple_qp(backend = "sparse")

# Solve
res <- model$solve()

expect_equal(res$x, c(0.4285714, 0.2142857), tolerance = TOL)
expect_equal(res$y, c(-1.5714286), tolerance = TOL)
expect_equal(res$z, c(0., 0.), tolerance = TOL)
expect_equal(res$z_lb, c(0., 0.), tolerance = TOL)
expect_equal(res$z_ub, c(0., 0.), tolerance = TOL)

## Update dense QP

model <- define_simple_qp(backend = "dense")

# Solve
model$solve()

P_new <- Matrix(c(8., 0., 0., 4.), 2, 2, sparse = TRUE)
A_new <- Matrix(c(1., -3.), 1, 2, sparse = TRUE)
h_new <- c(2., 1.)
x_ub_new <- c(Inf, 2.)

# Update and solve
model$update(P = as.matrix(P_new), A = as.matrix(A_new), h = h_new, x_ub = x_ub_new)
res <- model$solve()

expect_equal(res$x, c(0.2763157, 0.0921052), tolerance = TOL)
expect_equal(res$y, c(-1.2105263), tolerance = TOL)
expect_equal(res$z, c(0., 0.), tolerance = TOL)
expect_equal(res$z_lb, c(0., 0.), tolerance = TOL)
expect_equal(res$z_ub, c(0., 0.), tolerance = TOL)

## Update sparse QP

model <- define_simple_qp(backend = "sparse")

# Solve
model$solve()

P_new <- Matrix(c(8., 0., 0., 4.), 2, 2, sparse = TRUE)
A_new <- Matrix(c(1., -3.), 1, 2, sparse = TRUE)
h_new <- c(2., 1.)
x_ub_new <- c(Inf, 2.)

# Update and solve
model$update(P = P_new, A = A_new, h = h_new, x_ub = x_ub_new)
res <- model$solve()

expect_equal(res$x, c(0.2763157, 0.0921052), tolerance = TOL)
expect_equal(res$y, c(-1.2105263), tolerance = TOL)
expect_equal(res$z, c(0., 0.), tolerance = TOL)
expect_equal(res$z_lb, c(0., 0.), tolerance = TOL)
expect_equal(res$z_ub, c(0., 0.), tolerance = TOL)
