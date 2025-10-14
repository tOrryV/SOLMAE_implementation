import unittest
import random

from poly import Poly

PRIME_61 = (1 << 61) - 1
PRIME_127 = (1 << 127) - 1
PRIME_25519 = (1 << 255) - 19

BIG_MODS = [PRIME_61, PRIME_127, PRIME_25519]

def rand_poly(deg, mod):
    return Poly.from_list([(random.getrandbits(64) - random.getrandbits(64)) % mod
                           for _ in range(deg)])

def ref_mul(a_list, b_list, mod, deg):
    acc = [0] * deg
    for i, ai in enumerate(a_list):
        if ai == 0:
            continue
        for j, bj in enumerate(b_list):
            if bj == 0:
                continue
            k = i + j
            prod = (ai * bj) % mod
            if k < deg:
                acc[k] = (acc[k] + prod) % mod
            else:
                acc[k - deg] = (acc[k - deg] - prod) % mod
    return [x % mod for x in acc]


class TestPoly(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        random.seed(2024)

    def test_add_sub_neg_identities(self):
        mod, deg = PRIME_61, 32
        Poly.ring(mod, deg)
        zero, one = Poly.zero(), Poly.one()

        f = rand_poly(deg, mod)
        g = rand_poly(deg, mod)

        self.assertEqual((f + zero).to_list(), f.to_list())
        self.assertEqual((f - zero).to_list(), f.to_list())
        self.assertEqual((one * f).to_list(), f.to_list())

        self.assertEqual((f + (-f)).to_list(), zero.to_list())
        self.assertEqual((f - f).to_list(), zero.to_list())

        self.assertEqual((f + g).to_list(), (g + f).to_list())

    def test_degree_folding_rule(self):
        mod, deg = PRIME_127, 16
        Poly.ring(mod, deg)

        for r in range(0, 6):
            coeffs = [0] * (deg + r + 1)
            coeffs[deg + r] = 1
            h = Poly.from_list(coeffs)
            e_r = Poly.from_list([0]*r + [1])
            self.assertEqual(h.to_list(), (-e_r).to_list())

    def test_scalar_mul(self):
        mod, deg = PRIME_61, 32
        Poly.ring(mod, deg)

        f = rand_poly(deg, mod)
        num1 = random.randrange(0, mod)
        num2 = random.randrange(0, mod)

        self.assertEqual((num1 * f).to_list(), (f * num1).to_list())
        left = ((num1 + num2) % mod) * f
        right = num1 * f + num2 * f
        right2 = (num1 % mod) * f + (num2 % mod) * f
        L, R, R2 = left.to_list(), right.to_list(), right2.to_list()

        self.assertEqual(L, R)
        self.assertEqual(L, R2)

    def test_mul_matches_reference(self):
        mod, deg = PRIME_61, 64
        Poly.ring(mod, deg)
        for _ in range(20):
            f = rand_poly(deg, mod)
            g = rand_poly(deg, mod)
            got = (f * g).to_list()
            ref = ref_mul(f.to_list(), g.to_list(), mod, deg)
            self.assertEqual(got, ref)

    def test_ring_axioms_random(self):
        mod, deg = PRIME_25519, 64
        Poly.ring(mod, deg)
        for _ in range(10):
            f, g, h = rand_poly(deg, mod), rand_poly(deg, mod), rand_poly(deg, mod)

            left = f * (g + h)
            right = f * g + f * h
            self.assertEqual(left.to_list(), right.to_list())

            self.assertEqual(((f * g) * h).to_list(), (f * (g * h)).to_list())

    def test_from_list_reduction(self):
        mod, deg = PRIME_61, 16
        Poly.ring(mod, deg)

        long = [0] * (deg + 3)
        long[deg + 0] = 7
        long[deg + 1] = -5
        long[deg + 2] = 11
        f = Poly.from_list(long)

        expect = [0] * deg
        expect[0] = (-7) % mod
        expect[1] = 5 % mod
        expect[2] = (-11) % mod
        self.assertEqual(f.to_list(), expect)

    def test_perf_mul_512(self):
        mod, deg = PRIME_61, 512
        Poly.ring(mod, deg)
        f = rand_poly(deg, mod)
        g = rand_poly(deg, mod)
        h = f * g
        self.assertEqual(len(h.to_list()), deg)
