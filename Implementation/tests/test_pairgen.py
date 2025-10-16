import unittest
import random

from poly import Poly
from unifcrown import uniform_poly, crown_sample
from pairgen import pairgen, pairgen_seeded
from rng import expand_seed_to_mod_q, uniform_small

Q = 12289
D = 64
RADIUS = 3

def to_centered_list(coeffs, q):
    half = q // 2
    out = []
    for x in coeffs:
        x = int(x) % q
        out.append(x if x <= half else x - q)
    return out


class TestUnifCrownPairGen(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        random.seed(2025)

    def test_uniform_poly_basic(self):
        Poly.ring(Q, D)
        f = uniform_poly(Q, D)
        self.assertEqual(len(f.to_list()), D)
        self.assertTrue(all(0 <= c < Q for c in f.to_list()))

    def test_crown_sample_centered_bounds(self):
        Poly.ring(Q, D)
        f = crown_sample(Q, D, RADIUS)
        centered = to_centered_list(f.to_list(), Q)
        self.assertTrue(all(-RADIUS <= c <= RADIUS for c in centered))

    def test_pairgen_relation(self):
        Poly.ring(Q, D)
        a, s, b = pairgen(Q, D, RADIUS)

        e = b - (a * s)
        s_cent = to_centered_list(s.to_list(), Q)
        e_cent = to_centered_list(e.to_list(), Q)
        self.assertTrue(all(-RADIUS <= x <= RADIUS for x in s_cent))
        self.assertTrue(all(-RADIUS <= x <= RADIUS for x in e_cent))

        self.assertEqual(((a * s) + e).to_list(), b.to_list())

    def test_pairgen_seeded_a_and_relation(self):
        Poly.ring(Q, D)
        seed = b"unit-test-seed"

        a1, s1, b1 = pairgen_seeded(Q, D, seed, RADIUS)
        a2, s2, b2 = pairgen_seeded(Q, D, seed, RADIUS)

        a_expected = expand_seed_to_mod_q(seed + b"A", D, Q)
        self.assertEqual(a1.to_list(), a_expected)
        self.assertEqual(a2.to_list(), a_expected)

        e1 = b1 - (a1 * s1)
        e2 = b2 - (a2 * s2)
        e1_cent = to_centered_list(e1.to_list(), Q)
        e2_cent = to_centered_list(e2.to_list(), Q)
        self.assertTrue(all(-RADIUS <= x <= RADIUS for x in e1_cent))
        self.assertTrue(all(-RADIUS <= x <= RADIUS for x in e2_cent))

