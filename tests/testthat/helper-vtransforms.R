# Shared fixtures for the v-transform tests.

# "Proper" v-transforms: V(0) = 1, V(delta) = 0, V(1) = 1, V-shaped about the
# fulcrum. Vdegenerate() (the identity) is deliberately excluded.
vtransform_list <- function() {
  list(
    Vsymmetric = list(x = Vsymmetric(), delta = 0.5),
    Vlinear    = list(x = Vlinear(delta = 0.4), delta = 0.4),
    V2p        = list(x = V2p(delta = 0.4, kappa = 1.3), delta = 0.4),
    V2b        = list(x = V2b(delta = 0.35, kappa = 1.2), delta = 0.35),
    V3p        = list(x = V3p(delta = 0.45, kappa = 0.8, xi = 1.2), delta = 0.45),
    V3b        = list(x = V3b(delta = 0.3, kappa = 1.4, xi = 1.1), delta = 0.3)
  )
}

# A grid of interior points, excluding the endpoints.
u_grid <- function(n = 21) seq(0.02, 0.98, length.out = n)
