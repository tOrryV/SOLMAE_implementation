import unittest
import time
import statistics

import sample_precomp as sp
from cfft import fft_inplace, ifft_inplace, precompute_twiddles


class TestNTTPlan(unittest.TestCase):
    Q = 12289        # класичний простий модуль
    N = 256          # 512 | 12288

    def test_make_ntt_plan_basic(self):
        plan = sp.make_ntt_plan(self.Q, self.N)
        self.assertEqual(plan.q, self.Q)
        self.assertEqual(plan.n, self.N)
        self.assertEqual(len(plan.roots), self.N - 1)
        self.assertEqual(len(plan.roots_inv), self.N - 1)
        self.assertEqual(len(plan.bitrev), self.N)
        self.assertEqual(len(plan.tw_fwd), self.N)
        self.assertEqual(len(plan.tw_inv), self.N)

        psi = plan.psi
        self.assertEqual(pow(psi, 2*self.N, self.Q), 1)
        self.assertNotEqual(pow(psi, self.N, self.Q), 1)

        for i in range(len(plan.roots)):
            r = plan.roots[i]
            rinv = plan.roots_inv[i]
            self.assertEqual((r * rinv) % self.Q, 1)

    def test_precompute_for_sample_flags(self):
        params = sp.SampleParams(n=self.N, q=self.Q, sigma=3.0, use_ntt=True, use_fft=False)
        state = sp.precompute_for_sample(params)
        self.assertIsNotNone(state.ntt)
        self.assertIsNone(state.cfft)

        params2 = sp.SampleParams(n=self.N, q=self.Q, sigma=3.0, use_ntt=False, use_fft=True)
        state2 = sp.precompute_for_sample(params2)
        self.assertIsNone(state2.ntt)
        self.assertIsNotNone(state2.cfft)


class TestCFFTPlan(unittest.TestCase):
    N = 128

    def test_make_cfft_plan_and_roundtrip(self):
        plan = sp.make_cfft_plan(self.N)
        self.assertEqual(plan.n, self.N)
        self.assertEqual(len(plan.W), self.N)
        self.assertEqual(len(plan.Winv), self.N)
        self.assertEqual(len(plan.bitrev), self.N)

        for k in range(self.N):
            self.assertAlmostEqual((plan.W[k] * plan.Winv[k]).real, 1.0, places=12)
            self.assertAlmostEqual((plan.W[k] * plan.Winv[k]).imag, 0.0, places=12)

        import random
        random.seed(2025)
        a = [complex((random.random()-0.5)*2.0, (random.random()-0.5)*2.0) for _ in range(self.N)]
        x = list(a)
        fft_inplace(x, plan.W, plan.bitrev)
        ifft_inplace(x, plan.Winv, plan.bitrev)
        for i in range(self.N):
            self.assertAlmostEqual(x[i].real, a[i].real, places=9)
            self.assertAlmostEqual(x[i].imag, a[i].imag, places=9)


class TestCDT(unittest.TestCase):
    def test_cdt_monotone_and_sigma(self):
        sigma = 3.0
        table = sp.gaussian_cdt_build(sigma, tailcut=8.0)
        self.assertGreater(table.cutoff, 0)
        self.assertEqual(len(table.cdf_scaled), table.cutoff + 1)

        c = table.cdf_scaled
        for i in range(1, len(c)):
            self.assertGreaterEqual(c[i], c[i-1])

        S = 8192
        xs = sp.gaussian_vec_from_rng(table, S, sp.rand64_os)
        mu = statistics.mean(xs)
        sd = statistics.pstdev(xs)
        self.assertAlmostEqual(mu, 0.0, delta=0.3)
        self.assertTrue(0.8*sigma <= sd <= 1.2*sigma, msg=f"sd={sd:.3f}, sigma={sigma}")

    def test_cdt_sample_range_and_symmetry(self):
        sigma = 2.5
        table = sp.gaussian_cdt_build(sigma, tailcut=10.0)
        M = 4096
        xs = sp.gaussian_vec_from_rng(table, M, sp.rand64_os)
        self.assertTrue(all(abs(x) <= table.cutoff for x in xs))
        positives = sum(1 for x in xs if x > 0)
        ratio = positives / M
        self.assertTrue(0.40 <= ratio <= 0.60, f"pos ratio {ratio:.3f} out of [0.40,0.60]")


class TestSamplerState(unittest.TestCase):
    Q = 12289
    N = 256

    def test_noise_poly_length_and_modq(self):
        params = sp.SampleParams(n=self.N, q=self.Q, sigma=3.2, use_ntt=True, use_fft=True)
        state = sp.precompute_for_sample(params)

        e = sp.sample_noise_poly(state)
        self.assertEqual(len(e), self.N)
        sd = statistics.pstdev(e)
        self.assertTrue(0.75*params.sigma <= sd <= 1.25*params.sigma, f"sd={sd:.3f}")

        em = sp.sample_noise_poly_modq(state)
        self.assertEqual(len(em), self.N)
        self.assertTrue(all(0 <= x < self.Q for x in em))

    def test_attach_drbg_reproducible(self):
        params = sp.SampleParams(n=64, q=7681, sigma=2.0, use_ntt=False, use_fft=False)
        state1 = sp.precompute_for_sample(params)
        state2 = sp.precompute_for_sample(params)

        seed = b"sample-drbg-seed"
        sp.attach_drbg(state1, seed)
        sp.attach_drbg(state2, seed)

        v1 = sp.sample_noise_poly(state1)
        v2 = sp.sample_noise_poly(state2)
        self.assertEqual(v1, v2)

        state3 = sp.precompute_for_sample(params)
        sp.attach_drbg(state3, b"another-seed")
        v3 = sp.sample_noise_poly(state3)
        self.assertNotEqual(v1, v3)

    def test_precompute_state_internals(self):
        params = sp.SampleParams(n=128, q=12289, sigma=3.0, use_ntt=True, use_fft=True)
        state = sp.precompute_for_sample(params)
        self.assertIsNotNone(state.ntt)
        self.assertEqual(state.ntt.n, params.n)
        self.assertEqual(state.ntt.q, params.q)
        self.assertIsNotNone(state.cfft)
        self.assertEqual(state.cfft.n, params.n)
        self.assertGreater(state.cdt.cutoff, 0)
        self.assertEqual(len(state.cdt.cdf_scaled), state.cdt.cutoff + 1)

