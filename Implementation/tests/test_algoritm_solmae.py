import unittest

from algoritm_solmae import keygen_solmae, sign_solmae, verify_solmae
from params import n, k, q, d
from poly import Poly


class TestAlgoritmSolmae(unittest.TestCase):
    def setUp(self):
        self.msg = b"solmae e2e test message"

    def test_keygen_structure(self):
        pk, sk = keygen_solmae()

        self.assertIsInstance(pk, tuple)
        self.assertEqual(len(pk), 2)
        rho, t1 = pk
        self.assertIsInstance(rho, (bytes, bytearray))
        self.assertIsInstance(t1, list)
        self.assertEqual(len(t1), k)
        for t1i in t1:
            self.assertIsInstance(t1i, (bytes, bytearray))
            self.assertEqual(len(t1i), n)

        self.assertIsInstance(sk, tuple)
        self.assertEqual(len(sk), 5)
        s, e, t0, tr, pk2 = sk
        self.assertEqual(pk2, pk)
        self.assertIsInstance(tr, (bytes, bytearray))

        self.assertEqual(len(s), k)
        self.assertEqual(len(e), k)
        for pi in s + e:
            self.assertIsInstance(pi, Poly)
            self.assertEqual(len(pi.a), n)
            self.assertTrue(all(isinstance(x, int) and 0 <= x < q for x in pi.a))

        self.assertIsInstance(t0, list)
        self.assertEqual(len(t0), k)
        for t0i in t0:
            self.assertIsInstance(t0i, list)
            self.assertEqual(len(t0i), n)
            self.assertTrue(all(isinstance(x, int) for x in t0i))

    def test_sign_verify_positive(self):
        pk, sk = keygen_solmae()
        sig = sign_solmae(sk, self.msg)

        self.assertIsInstance(sig, tuple)
        self.assertEqual(len(sig), 3)
        z, c, w1 = sig

        self.assertIsInstance(z, list)
        self.assertEqual(len(z), k)
        for zi in z:
            self.assertIsInstance(zi, list)
            self.assertEqual(len(zi), n)
            self.assertTrue(all(isinstance(x, int) and 0 <= x < q for x in zi))

        self.assertIsInstance(c, (bytes, bytearray))
        self.assertIsInstance(w1, list)
        self.assertEqual(len(w1), k)
        for w1i in w1:
            self.assertIsInstance(w1i, (bytes, bytearray))
            self.assertEqual(len(w1i), n)

        self.assertTrue(verify_solmae(pk, self.msg, sig))

    def test_verify_rejects_changed_message(self):
        pk, sk = keygen_solmae()
        sig = sign_solmae(sk, self.msg)
        other = self.msg + b" (tampered)"
        self.assertFalse(verify_solmae(pk, other, sig))

    def test_verify_rejects_tampered_challenge(self):
        pk, sk = keygen_solmae()
        z, c, w1 = sign_solmae(sk, self.msg)
        c_bad = bytearray(c)
        c_bad[0] ^= 0x01
        self.assertFalse(verify_solmae(pk, self.msg, (z, bytes(c_bad), w1)))

    def test_verify_rejects_tampered_z(self):
        pk, sk = keygen_solmae()
        z, c, w1 = sign_solmae(sk, self.msg)

        z_bad = [[0 for _ in zi] for zi in z]

        self.assertFalse(verify_solmae(pk, self.msg, (z_bad, c, w1)))

    def test_verify_idempotent(self):
        pk, sk = keygen_solmae()
        sig = sign_solmae(sk, self.msg)
        for _ in range(3):
            self.assertTrue(verify_solmae(pk, self.msg, sig))

