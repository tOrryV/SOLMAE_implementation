import unittest
import random
import math
import ntrusolve as ns

def ref_negacyclic_mul(num1, num2):
    n = len(num1)
    assert len(num2) == n
    c = [0]*n
    for i, ai in enumerate(num1):
        if ai == 0: continue
        for j, bj in enumerate(num2):
            if bj == 0: continue
            k = i + j
            if k < n:
                c[k] += ai * bj
            else:
                c[k - n] -= ai * bj
    return c


def rand_vec(n, lo=-10, hi=10):
    return [random.randint(lo, hi) for _ in range(n)]


class TestZArithmetic(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        random.seed(2025)

    def test_z_negacyclic_mul_matches_reference(self):
        for n in (8, 16, 32):
            for _ in range(20):
                num1 = rand_vec(n, -5, 5)
                num2 = rand_vec(n, -5, 5)
                got = ns.z_negacyclic_mul(num1, num2)
                ref = ref_negacyclic_mul(num1, num2)
                self.assertEqual(got, ref)

    def test_z_add_sub_scalar_identities(self):
        n = 16
        num1 = rand_vec(n, -10, 10)
        num2 = rand_vec(n, -10, 10)
        k = random.randint(-7, 7) or 3

        self.assertEqual(ns.z_add(num1, num2), [num1[i]+num2[i] for i in range(n)])
        self.assertEqual(ns.z_sub(num1, num2), [num1[i]-num2[i] for i in range(n)])
        self.assertEqual(ns.z_scalar_mul(num1, k), [k*num1[i] for i in range(n)])

    def test_centered_mod_q_bounds(self):
        n, mod = 16, 12289
        num1 = rand_vec(n, -3*mod, 3*mod)
        c = ns.z_centered_mod_q(num1, mod)
        self.assertEqual(len(c), n)
        half = mod // 2
        self.assertTrue(all(-half <= x <= half for x in c))

    def test_round_div_correctness(self):
        n = 32
        num1 = rand_vec(n, -1000, 1000)
        k = random.randint(2, 17)
        r1 = ns.z_round_div(num1, k)
        r2 = [int(round(x / k)) for x in num1]
        for i in range(n):
            self.assertTrue(abs(r1[i] - r2[i]) <= 1)

    def test_closest_rounding_div_q(self):
        n, mod = 32, 257
        num1 = rand_vec(n, -2000, 2000)
        c = ns.closest_rounding_div_q(num1, mod)
        for i in range(n):
            self.assertTrue(abs(mod*c[i] - num1[i]) <= mod//2 + 1)


class TestEvenOdd(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        random.seed(2025)

    def test_split_merge_roundtrip(self):
        for n in (8, 16, 32):
            num1 = rand_vec(n, -50, 50)
            ae, ao = ns.split_even_odd(num1)
            self.assertEqual(len(ae), n//2)
            self.assertEqual(len(ao), n//2)
            back = ns.merge_even_odd(ae, ao)
            self.assertEqual(back, num1)

    def test_lift_up_matches_merge(self):
        for n in (8, 16, 32):
            num1 = rand_vec(n, -20, 20)
            ae, ao = ns.split_even_odd(num1)
            m1 = ns.merge_even_odd(ae, ao)
            m2 = ns.lift_up(ae, ao)
            self.assertEqual(m1, m2)


class TestReduceAndBase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        random.seed(2025)

    def test_basecase_solve_invariant(self):
        for _ in range(100):
            F0 = random.randint(-50, 50)
            G0 = random.randint(-50, 50)
            if F0 == 0 and G0 == 0:
                G0 = 1
            g = math.gcd(F0, G0)
            mult = random.randint(1, 7)
            q = g * mult if g != 0 else mult
            f, gcoef = ns.basecase_solve([F0], [G0], q)
            self.assertEqual(len(f), 1)
            self.assertEqual(len(gcoef), 1)
            lhs = f[0]*G0 - gcoef[0]*F0
            self.assertEqual(lhs, q)

    def test_reduce_target_invariance_and_norm(self):
        n, mod = 16, 12289
        F = rand_vec(n, -50, 50)
        G = rand_vec(n, -50, 50)
        f = rand_vec(n, -100, 100)
        g = rand_vec(n, -100, 100)

        H_before = ns.z_sub(ns.z_negacyclic_mul(f, G), ns.z_negacyclic_mul(g, F))
        norm_before = sum(x*x for x in f) + sum(x*x for x in g)

        fr, gr = ns.reduce_target(F, G, f, g, mod)

        H_after = ns.z_sub(ns.z_negacyclic_mul(fr, G), ns.z_negacyclic_mul(gr, F))

        for i in range(n):
            self.assertEqual((H_before[i] - H_after[i]) % mod, 0)

        norm_after = sum(x*x for x in fr) + sum(x*x for x in gr)
        self.assertLessEqual(norm_after, norm_before + 1e-9)

    def test_ring_axioms_small(self):
        n = 8
        num1 = rand_vec(n, -5, 5)
        num2 = rand_vec(n, -5, 5)
        c = rand_vec(n, -5, 5)

        left = ns.z_negacyclic_mul(num1, ns.z_add(num2, c))
        right = ns.z_add(ns.z_negacyclic_mul(num1, num2), ns.z_negacyclic_mul(num1, c))
        self.assertEqual(left, right)

        left = ns.z_negacyclic_mul(ns.z_negacyclic_mul(num1, num2), c)
        right = ns.z_negacyclic_mul(num1, ns.z_negacyclic_mul(num2, c))
        self.assertEqual(left, right)


