import unittest
import random
import ntt as ntt
from poly import Poly


Q_SMALL, D_SMALL = 769, 8
Q_MED,  D_MED = 12289, 256
Q_LARGE, D_LARGE = 12289, 1024


def rand_vec(n, mod):
    return [random.randrange(0, mod) for _ in range(n)]

def ref_negacyclic_naive(num1, num2, mod, deg):
    c = [0] * deg
    for i in range(deg):
        ai = num1[i]
        if ai == 0: continue
        for j in range(deg):
            bj = num2[j]
            if bj == 0: continue
            k = i + j
            prod = (ai * bj) % mod
            if k < deg:
                c[k] = (c[k] + prod) % mod
            else:
                c[k - deg] = (c[k - deg] - prod) % mod
    return [x % mod for x in c]


class TestNTT(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        random.seed(2025)

    def test_roundtrip_small(self):
        mod, deg = Q_SMALL, D_SMALL
        psi = ntt.primitive_root_for(mod, 2*deg)
        roots, roots_inv, bitrev = ntt.precompute_roots(mod, deg, psi)

        for _ in range(20):
            a = rand_vec(deg, mod)
            A = a[:]
            ntt.ntt_inplace(A, mod, roots, bitrev)
            ntt.intt_inplace(A, mod, roots_inv, bitrev)
            self.assertEqual([x % mod for x in A], a)

    def test_roundtrip_medium(self):
        mod, deg = Q_MED, D_MED
        psi = ntt.primitive_root_for(mod, 2*deg)
        roots, roots_inv, bitrev = ntt.precompute_roots(mod, deg, psi)
        a = rand_vec(deg, mod)
        A = a[:]
        ntt.ntt_inplace(A, mod, roots, bitrev)
        ntt.intt_inplace(A, mod, roots_inv, bitrev)
        self.assertEqual([x % mod for x in A], a)

    def test_negacyclic_matches_naive_small(self):
        mod, deg = Q_SMALL, D_SMALL
        psi = ntt.primitive_root_for(mod, 2*deg)
        roots, roots_inv, bitrev = ntt.precompute_roots(mod, deg, psi)
        tw_fwd, tw_inv = ntt.precompute_twists(mod, deg, psi)

        for _ in range(30):
            a = rand_vec(deg, mod)
            b = rand_vec(deg, mod)
            got = ntt.negacyclic_convolution(a, b, mod, roots, roots_inv, bitrev, tw_fwd, tw_inv)
            ref = ref_negacyclic_naive(a, b, mod, deg)
            self.assertEqual(got, ref)

    def test_negacyclic_matches_naive_medium(self):
        mod, deg = Q_MED, D_MED
        psi = ntt.primitive_root_for(mod, 2*deg)
        roots, roots_inv, bitrev = ntt.precompute_roots(mod, deg, psi)
        tw_fwd, tw_inv = ntt.precompute_twists(mod, deg, psi)

        a = rand_vec(deg, mod)
        b = rand_vec(deg, mod)
        got = ntt.negacyclic_convolution(a, b, mod, roots, roots_inv, bitrev, tw_fwd, tw_inv)
        ref = ref_negacyclic_naive(a, b, mod, deg)
        self.assertEqual(got, ref)

    def test_poly_mul_rq_ntt_wrapper_small(self):
        mod, deg = Q_SMALL, D_SMALL
        for _ in range(10):
            a = rand_vec(deg, mod)
            b = rand_vec(deg, mod)
            got = ntt.poly_mul_rq_ntt(a, b, mod, deg)
            ref = ref_negacyclic_naive(a, b, mod, deg)
            self.assertEqual(got, ref)

    def test_poly_ntt_vs_poly_naive_medium(self):
        mod, deg = Q_MED, D_MED
        Poly.ring(mod, deg)
        for _ in range(5):
            f = Poly.from_list(rand_vec(deg, mod))
            g = Poly.from_list(rand_vec(deg, mod))
            h_ntt = ntt.poly_mul_rq_ntt(f.to_list(), g.to_list(), mod, deg)
            h_poly = (f * g).to_list()
            self.assertEqual(h_poly, h_ntt)

    def test_perf_large_1024(self):
        mod, deg = Q_LARGE, D_LARGE
        a = rand_vec(deg, mod)
        b = rand_vec(deg, mod)
        c = ntt.poly_mul_rq_ntt(a, b, mod, deg)
        self.assertEqual(len(c), deg)
