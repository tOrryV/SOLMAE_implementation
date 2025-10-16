import math
import unittest
import random
import time

from poly import Poly
from unifcrown import uniform_poly, uniform_pair, crown_sample, crown_pair

Q = 12289
D = 64
RADIUS = 3

def to_centered_list(coeffs, q):
    half = q // 2
    out = []
    for x in coeffs:
        x %= q
        out.append(x if x <= half else x - q)
    return out


class TestUnifCrown(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        random.seed(2025)
        Poly.ring(Q, D)

    def test_uniform_poly_basic(self):
        f = uniform_poly(Q, D)
        coeffs = f.to_list()
        self.assertEqual(len(coeffs), D)
        self.assertTrue(all(0 <= c < Q for c in coeffs))

    def test_uniform_pair_independence(self):
        a, b = uniform_pair(Q, D)
        A, B = a.to_list(), b.to_list()
        self.assertNotEqual(A, B)

        equal_positions = sum(1 for i in range(D) if A[i] == B[i])
        self.assertLess(equal_positions, D // 4)

        Ac = to_centered_list(A, Q)
        Bc = to_centered_list(B, Q)
        meanA = sum(Ac) / D
        meanB = sum(Bc) / D
        cov = sum((Ac[i] - meanA) * (Bc[i] - meanB) for i in range(D))
        varA = sum((x - meanA) ** 2 for x in Ac)
        varB = sum((y - meanB) ** 2 for y in Bc)

        if varA > 0 and varB > 0:
            r = cov / math.sqrt(varA * varB)
            self.assertLess(abs(r), 0.15, f"unexpected correlation r={r:.3f}")

    def test_crown_sample_bounds(self):
        f = crown_sample(Q, D, RADIUS)
        centered = to_centered_list(f.to_list(), Q)
        self.assertTrue(all(-RADIUS <= c <= RADIUS for c in centered))

    def test_crown_pair_bounds(self):
        s, e = crown_pair(Q, D, RADIUS)
        sC = to_centered_list(s.to_list(), Q)
        eC = to_centered_list(e.to_list(), Q)
        self.assertTrue(all(-RADIUS <= c <= RADIUS for c in sC))
        self.assertTrue(all(-RADIUS <= c <= RADIUS for c in eC))

    def test_uniform_poly_lsb_balance(self):
        trials = 200
        ones = 0
        total = 0
        for _ in range(trials):
            f = uniform_poly(Q, D)
            for c in f.to_list():
                ones += (c & 1)
                total += 1
        ratio = ones / total
        self.assertTrue(0.40 <= ratio <= 0.60, f"LSB ratio out of range: {ratio:.3f}")

    def test_perf_uniform_512(self):
        Poly.ring(Q, 512)
        f = uniform_poly(Q, 512)
        g = crown_sample(Q, 512, RADIUS)
        self.assertEqual(len(f.to_list()), 512)
        self.assertEqual(len(g.to_list()), 512)
        Poly.ring(Q, D)
