import math
import statistics as stats
import unittest
from unittest.mock import patch
import samplers
from sample_precomp import (
    SampleParams,
    precompute_for_sample,
    attach_drbg,
    next_u64_from_state,
)


class RNGPatchedTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        params = SampleParams(
            n=8,
            q=12289,
            sigma=2.0,
        )
        cls._rng_state = precompute_for_sample(params)

        attach_drbg(cls._rng_state, seed=b"unit-test-seed-for-samplers")

        def _rand_u64_deterministic():
            return next_u64_from_state(cls._rng_state)

        cls._patcher = patch("samplers._rand_u64", side_effect=_rand_u64_deterministic)
        cls._patcher.start()

    @classmethod
    def tearDownClass(cls):
        cls._patcher.stop()


class TestNSampler(RNGPatchedTestCase):
    def test_mean_and_variance(self):
        N = 20000
        vals = []
        for _ in range(N // 2):
            z0, z1 = samplers.n_sampler(2)
            vals.append(z0)
            vals.append(z1)

        mu = stats.fmean(vals)
        var = stats.pvariance(vals)

        self.assertLess(abs(mu), 0.03)
        self.assertGreater(var, 0.95)
        self.assertLess(var, 1.05)


class TestZSampler(RNGPatchedTestCase):
    def _moments(self, mu, sigma, N=25000):
        xs = [samplers.z_sampler(mu, sigma) for _ in range(N)]
        m = stats.fmean(xs)
        v = stats.pvariance(xs)
        return m, v

    def test_mu_0_sigma_125(self):
        m, v = self._moments(0.0, 1.25)
        self.assertLess(abs(m - 0.0), 0.06)
        self.assertGreater(v, (1.25**2) * 0.85)
        self.assertLess(v, (1.25**2) * 1.20)

    def test_mu_025_sigma_2(self):
        m, v = self._moments(0.25, 2.0)
        self.assertLess(abs(m - 0.25), 0.06)
        self.assertGreater(v, (2.0**2) * 0.85)
        self.assertLess(v, (2.0**2) * 1.20)

    def test_mu_05_sigma_2(self):
        m, v = self._moments(0.5, 2.0)
        self.assertLess(abs(m - 0.5), 0.06)
        self.assertGreater(v, (2.0**2) * 0.85)
        self.assertLess(v, (2.0**2) * 1.20)


class TestPeikertSampler(RNGPatchedTestCase):
    def test_shape_and_variance(self):
        d = 64
        Sigma = [1.1] * d
        t = [0.0] * d
        eta = 1.75

        N = 15000
        samples0 = []
        for _ in range(N):
            z = samplers.peikert_sampler(t, Sigma, eta, d)
            self.assertIsInstance(z, list)
            self.assertEqual(len(z), d)
            samples0.append(z[0])

        m = stats.fmean(samples0)
        v = stats.pvariance(samples0)

        target = eta**2 + Sigma[0]**2
        self.assertLess(abs(m), 0.08)
        self.assertGreater(v, target * 0.75)
        self.assertLess(v, target * 1.35)


class DummySK:
    def __init__(self, d, eta):
        self.d = d

        one = [1.0 + 0.0j] * d
        zero = [0.0 + 0.0j] * d

        self.b1_fft = (one[:], zero[:])
        self.b2_fft = (zero[:], one[:])
        self.b2_tilde_fft = (zero[:], one[:])

        self.Sigma1_fft = [0.9] * d
        self.Sigma2_fft = [1.1] * d

        self.beta1_fft = (one[:], one[:])
        self.beta2_fft = (one[:], one[:])

        self.q = 12289
        self.eta = eta


class DummyParams:
    def __init__(self, eta):
        self.eta = eta


class TestSampleIntegration(RNGPatchedTestCase):
    def test_sample_smoke(self):
        d = 32
        eta = 1.5
        sk = DummySK(d, eta)
        params = DummyParams(eta)

        c_fft = ([0.0 + 0.0j] * d, [0.0 + 0.0j] * d)

        v1, v2 = samplers.sample(c_fft, sk, params)
        self.assertIsInstance(v1, list)
        self.assertIsInstance(v2, list)
        self.assertEqual(len(v1), d)
        self.assertEqual(len(v2), d)
        self.assertTrue(all(isinstance(x, complex) for x in v1))
        self.assertTrue(all(isinstance(x, complex) for x in v2))

        norm_sq = sum((x.real**2 + x.imag**2) for x in v1) + \
                  sum((x.real**2 + x.imag**2) for x in v2)
        self.assertLess(norm_sq, 5_000.0)

    def test_sample_repeated_runs_distribution(self):
        d = 16
        eta = 1.25
        sk = DummySK(d, eta)
        params = DummyParams(eta)
        c_fft = ([0.0 + 0.0j] * d, [0.0 + 0.0j] * d)

        sums = []
        for _ in range(200):
            v1, v2 = samplers.sample(c_fft, sk, params)
            sums.append(sum(x.real for x in v1) + sum(x.real for x in v2))

        mu = stats.fmean(sums)
        sd = math.sqrt(stats.pvariance(sums))

        self.assertLess(abs(mu), 5.0)
        self.assertGreater(sd, 1.0)
