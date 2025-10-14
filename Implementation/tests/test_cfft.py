import unittest
import random
import math
import cmath

import cfft


EPS = 1e-9


def rand_complex_vec(n, scale=1.0):
    return [(random.random() - 0.5) * 2 * scale + 1j * ((random.random() - 0.5) * 2 * scale) for _ in range(n)]


def rand_real_vec(n, scale=1.0):
    return [(random.random() - 0.5) * 2 * scale for _ in range(n)]


def dft(x):
    n = len(x)
    out = []
    for k in range(n):
        s = 0j
        for t in range(n):
            s += x[t] * cmath.exp(-2j * math.pi * k * t / n)
        out.append(s)
    return out


def idft(X):
    n = len(X)
    out = []
    for t in range(n):
        s = 0j
        for k in range(n):
            s += X[k] * cmath.exp(+2j * math.pi * k * t / n)
        out.append(s / n)
    return out


def circular_conv(num1, num2):
    n = len(num1)
    out = [0j] * n
    for i in range(n):
        s = 0j
        for j in range(n):
            s += num1[j] * num2[(i - j) % n]
        out[i] = s
    return out


class TestCFFT(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        random.seed(2025)

    def test_roundtrip_small(self):
        n = 16
        self.assertTrue(cfft.is_power_of_two(n))
        x = rand_complex_vec(n, scale=5.0)
        X = cfft.fft(x)
        x2 = cfft.ifft(X)
        self.assertLessEqual(cfft.max_abs_diff(x, x2), 1e-9)

    def test_roundtrip_medium_inplace(self):
        n = 256
        self.assertTrue(cfft.is_power_of_two(n))
        x = rand_complex_vec(n, scale=1.0)
        a = list(x)
        W, Winv = cfft.precompute_twiddles(n)
        perm = cfft.bitrev_permutation(n)
        cfft.fft_inplace(a, W, perm)
        cfft.ifft_inplace(a, Winv, perm)
        self.assertLessEqual(cfft.max_abs_diff(x, a), 1e-9)

    def test_matches_naive_dft(self):
        for n in (8, 16, 32):
            x = rand_complex_vec(n, scale=2.0)
            X_ref = dft(x)
            X = cfft.fft(x)
            self.assertLessEqual(cfft.max_abs_diff(X, X_ref), 1e-8)

            x_ref = idft(X_ref)
            x_back = cfft.ifft(X)
            self.assertLessEqual(cfft.max_abs_diff(x_ref, x_back), 1e-9)

    def test_circular_convolution_via_fft(self):
        n = 64
        a = rand_complex_vec(n)
        b = rand_complex_vec(n)
        A = cfft.fft(a)
        B = cfft.fft(b)
        C = cfft.hadamard_product(A, B)
        c_fft = cfft.ifft(C)

        c_ref = circular_conv(a, b)
        self.assertLessEqual(cfft.max_abs_diff(c_fft, c_ref), 1e-8)

    def test_real_wrappers(self):
        n = 128
        x = rand_real_vec(n, scale=3.0)

        X1 = cfft.fft_real(x)
        X2 = cfft.fft([complex(t, 0.0) for t in x])
        self.assertLessEqual(cfft.max_abs_diff(X1, X2), 1e-12)

        x_back = cfft.ifft_to_real(X1)
        self.assertLessEqual(cfft.max_abs_diff([complex(t, 0) for t in x], [complex(t, 0) for t in x_back]), 1e-9)

    def test_pointwise_ops(self):
        n = 32
        A = rand_complex_vec(n)
        B = rand_complex_vec(n)
        alpha = 1.234 - 0.5j

        self.assertEqual(len(cfft.add_complex(A, B)), n)
        self.assertEqual(len(cfft.sub_complex(A, B)), n)
        self.assertEqual(len(cfft.hadamard_product(A, B)), n)
        self.assertEqual(len(cfft.scale_complex(A, alpha)), n)

        Z = [0j] * n
        self.assertTrue(cfft.is_close_vec(cfft.add_complex(A, Z), A, tol=EPS))
        self.assertTrue(cfft.is_close_vec(cfft.sub_complex(A, Z), A, tol=EPS))
        self.assertTrue(cfft.is_close_vec(cfft.scale_complex(A, 1+0j), A, tol=EPS))
