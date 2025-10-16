from rng import uniform_mod_q, uniform_small
from poly import Poly


def uniform_poly(mod: int, deg):
    coeffs = uniform_mod_q(deg, mod)
    return Poly.from_list(coeffs)


def uniform_pair(mod, deg):
    return uniform_poly(mod, deg), uniform_poly(mod, deg)


def crown_sample(mod, deg, radius):
    coeffs = uniform_small(deg, radius)
    return Poly.from_list([(x % mod) for x in coeffs])


def crown_pair(mod, deg, radius):
    return crown_sample(mod, deg, radius), crown_sample(mod, deg, radius)
