# Shared fixtures for the v-transform tests.

# Every v-transform: V(0) = 1, V(delta) = 0, V(1) = 1, V-shaped about the
# fulcrum, which lies strictly inside (0, 1).
vtransform_list <- function() {
  list(
    vsymmetric = list(x = vsymmetric(), delta = 0.5),
    vlinear    = list(x = vlinear(delta = 0.4), delta = 0.4),
    v2p        = list(x = v2p(delta = 0.4, kappa = 1.3), delta = 0.4),
    v2b        = list(x = v2b(delta = 0.35, kappa = 1.2), delta = 0.35),
    v3p        = list(x = v3p(delta = 0.45, kappa = 0.8, xi = 1.2), delta = 0.45),
    v3b        = list(x = v3b(delta = 0.3, kappa = 1.4, xi = 1.1), delta = 0.3)
  )
}

# A grid of interior points, excluding the endpoints.
u_grid <- function(n = 21) seq(0.02, 0.98, length.out = n)
