import unittest
import random
import math
from comp_decom import compress, decompress


def centered_mod_q(a, q):
    a %= q
    half_down = q // 2
    if a > half_down:
        a -= q
    return a


class TestCompressDecompress(unittest.TestCase):
    def setUp(self):
        random.seed(12345)

    def test_round_trip_small_params(self):
        q = 769
        d = 32
        b = 9
        slen = d * b

        vals = [random.randint(-2**(b-1)//2, 2**(b-1)//2 - 1) for _ in range(d)]
        s1 = [(v % q) for v in vals]

        blob = compress(s1, slen, q=q)
        self.assertEqual(len(blob), math.ceil(slen / 8), "bytes len must be ceil(slen/8)")

        s1_dec = decompress(blob, slen, d, q)
        self.assertIsInstance(s1_dec, list)
        self.assertEqual(len(s1_dec), d)

        got = [centered_mod_q(x, q) for x in s1_dec]
        want = [centered_mod_q(x, q) for x in s1]
        self.assertEqual(got, want)

    def test_saturation_at_bounds(self):
        q = 12289
        d = 16
        b = 8
        slen = d * b

        min_v = -(1 << (b - 1))
        max_v = (1 << (b - 1)) - 1

        orig_centered = [
            min_v - 5, min_v - 1, min_v, min_v + 1,
            -1, 0, 1, max_v - 1, max_v, max_v + 1, max_v + 7,
            -1234, 1234, -(1 << 15), (1 << 15) - 1, 0
        ]

        s1 = [(v % q) for v in orig_centered]

        blob = compress(s1, slen, q=q)
        s1_dec = decompress(blob, slen, d, q)
        self.assertIsNotNone(s1_dec)
        dec_centered = [centered_mod_q(x, q) for x in s1_dec]

        expected_centered = [
            min(max(centered_mod_q(v % q, q), min_v), max_v)
            for v in orig_centered
        ]

        self.assertEqual(dec_centered, expected_centered)

    def test_exact_byte_length(self):
        q = 12289
        d = 17
        b = 7
        slen = d * b

        s1 = [random.randrange(q) for _ in range(d)]
        blob = compress(s1, slen, q=q)
        self.assertEqual(len(blob), 15)
        s1_dec = decompress(blob, slen, d, q)
        self.assertEqual(len(s1_dec), d)

    def test_invalid_b_raises_in_compress(self):
        q = 12289
        d = 10
        slen = 0
        s1 = [0] * d
        with self.assertRaises(ValueError):
            _ = compress(s1, slen, q=q)

    def test_insufficient_bytes_returns_none(self):
        q = 12289
        d = 16
        b = 9
        slen = d * b
        blob = bytes(10)

        res = decompress(blob, slen, d, q)
        self.assertIsNone(res, "insufficient data length should return None")

    def test_sign_extension_edges(self):
        q = 12289
        d = 8
        b = 10
        slen = d * b

        min_v = -(1 << (b - 1))
        max_v = (1 << (b - 1)) - 1

        center_vals = [min_v, min_v + 1, -1, 0, 1, max_v - 1, max_v]
        while len(center_vals) < d:
            center_vals.append(0)

        s1 = [(x % q) for x in center_vals]
        blob = compress(s1, slen, q=q)
        s1_dec = decompress(blob, slen, d, q)

        got = [centered_mod_q(x, q) for x in s1_dec]
        want = [min(max(v, min_v), max_v) for v in center_vals]
        self.assertEqual(got, want)

    def test_non_divisible_slen_floor_bits(self):
        q = 12289
        d = 10
        b = 9
        slen = d * b + 5

        center_vals = [random.randint(-(1 << (b - 2)), (1 << (b - 2)) - 1) for _ in range(d)]
        s1 = [(x % q) for x in center_vals]

        blob = compress(s1, slen, q=q)
        self.assertEqual(len(blob), math.ceil(slen / 8))

        s1_dec = decompress(blob, slen, d, q)
        got = [centered_mod_q(x, q) for x in s1_dec]
        want = center_vals
        self.assertEqual(got, want)

