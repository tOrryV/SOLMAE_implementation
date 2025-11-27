from unifcrown import uniform_poly, crown_sample
from poly import Poly


def pairgen(mod, deg, noise_radius):
    s = crown_sample(mod, deg, noise_radius)
    e = crown_sample(mod, deg, noise_radius)
    a = uniform_poly(mod, deg)
    b = a * s + e
    return a, s, b


def pairgen_seeded(mod, deg, seed, noise_radius):
    from rng import expand_seed_to_mod_q, uniform_small
    coeffs_a = expand_seed_to_mod_q(seed + b"A", deg, mod)
    coeffs_s = uniform_small(deg, noise_radius)
    coeffs_e = uniform_small(deg, noise_radius)
    a = Poly.from_list(coeffs_a)
    s = Poly.from_list([(x % mod) for x in coeffs_s])
    e = Poly.from_list([(x % mod) for x in coeffs_e])
    b = a * s + e
    return a, s, b
