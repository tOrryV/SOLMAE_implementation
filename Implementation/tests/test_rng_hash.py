import unittest
import statistics

import rng
import hashing


Q = 2305843009213693951


class TestRNG(unittest.TestCase):

    def test_hmacdrbg_reproducible(self):
        seed = b"abc123"
        gen1 = rng.HMACDRBG(seed).generate(64)
        gen2 = rng.HMACDRBG(seed).generate(64)
        self.assertEqual(gen1, gen2)
        gen3 = rng.HMACDRBG(b"xyz").generate(64)
        self.assertNotEqual(gen1, gen3)

    def test_uniform_mod_q_range(self):
        v = rng.uniform_mod_q(100, Q)
        for x in v:
            self.assertTrue(0 <= x < Q)

    def test_uniform_small_distribution(self):
        bound = 3
        d = 1000
        v = rng.uniform_small(d, bound)
        self.assertTrue(all(-bound <= x <= bound for x in v))
        self.assertAlmostEqual(statistics.mean(v), 0, delta=0.5)

    def test_sample_cbd_basic(self):
        data = bytes([0b10101010] * 8)
        out = rng.sample_cbd(data, eta=2)
        self.assertIsInstance(out, list)
        self.assertTrue(all(isinstance(x, int) for x in out))
        v = rng.sample_cbd_random(1000, eta=2)
        self.assertAlmostEqual(statistics.mean(v), 0, delta=0.3)

    def test_expand_seed_to_mod_q(self):
        seed = b"seed_seed"
        a = rng.expand_seed_to_mod_q(seed, 10, Q)
        b = rng.expand_seed_to_mod_q(seed, 10, Q)
        self.assertEqual(a, b)
        self.assertTrue(all(0 <= x < Q for x in a))


class TestHashing(unittest.TestCase):

    def test_sha256_stability(self):
        h1 = hashing.H_sha256([b"abc"], domain=hashing.DOM_SEP["H_msg"])
        h2 = hashing.H_sha256([b"abc"], domain=hashing.DOM_SEP["H_msg"])
        self.assertEqual(h1, h2)
        h3 = hashing.H_sha256([b"abc"], domain=hashing.DOM_SEP["H_pk"])
        self.assertNotEqual(h1, h3)
        h4 = hashing.H_sha256([b"abd"], domain=hashing.DOM_SEP["H_msg"])
        self.assertNotEqual(h1, h4)

    def test_sha256_int_mod_q(self):
        val = hashing.H_sha256_int([b"message"], Q)
        self.assertTrue(0 <= val < Q)

    def test_xof_consistency(self):
        data = [b"hello", b"world"]

        out1 = hashing.XOF_shake128(data, 64)
        out2 = hashing.XOF_shake128(data, 64)
        self.assertEqual(out1, out2)
        out3 = hashing.XOF_shake128(data, 32)
        self.assertEqual(out1[:32], out3)
        self.assertEqual(len(out3), 32)
        self.assertEqual(len(out1), 64)

        out4 = hashing.XOF_shake128([b"hello", b"WORLD"], 64)
        self.assertNotEqual(out1, out4)

        out5 = hashing.XOF_shake256(data, 64)
        self.assertNotEqual(out1, out5)

    def test_h_functions_basic(self):
        msg = b"msg"
        pk = b"pk"
        mu = b"mu"
        transcript = b"transcript"

        val_msg = hashing.H_msg_to_int_mod_q(msg, Q)
        self.assertTrue(0 <= val_msg < Q)

        bind1 = hashing.H_pk_bind(pk)
        bind2 = hashing.H_pk_bind(pk)
        self.assertEqual(bind1, bind2)

        chal = hashing.H_challenge(mu, transcript, 64)
        self.assertEqual(len(chal), 64)

        seed_exp = hashing.H_seed_expand(b"seed", 64)
        self.assertEqual(len(seed_exp), 64)

    def test_h_to_small_poly(self):
        seed = b"12345"
        d = 64
        eta = 2
        poly1 = hashing.H_to_small_poly(seed, d, eta)
        poly2 = hashing.H_to_small_poly(seed, d, eta)
        self.assertEqual(poly1, poly2)
        self.assertEqual(len(poly1), d)
        self.assertTrue(all(isinstance(x, int) for x in poly1))
