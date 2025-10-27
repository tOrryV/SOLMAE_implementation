import unittest
import random
import math
import modular as mod

PRIME_61 = (1 << 61) - 1                   # 2^61 - 1
PRIME_127 = (1 << 127) - 1                 # 2^127 - 1 (Мерсена)
PRIME_25519 = (1 << 255) - 19              # 2^255 - 19 (відомий простий)

BIG_MODS = [PRIME_61, PRIME_127, PRIME_25519]

def rand_big(bits):
    return random.getrandbits(bits) - random.getrandbits(bits // 2)

class TestModularBig(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        random.seed(1337)

    def test_add_sub_mod_properties(self):
        for mod_ in BIG_MODS:
            for _ in range(200):
                num1 = rand_big(1024)
                num2 = rand_big(1024)
                ref_add = ((num1 % mod_) + (num2 % mod_)) % mod_
                ref_sub = ((num1 % mod_) - (num2 % mod_)) % mod_

                got_add = mod.add_mod(num1, num2, mod_)
                got_sub = mod.sub_mod(num1, num2, mod_)

                self.assertEqual(got_add, ref_add, f"add_mod failed for q={mod_}")
                self.assertEqual(got_sub, ref_sub, f"sub_mod failed for q={mod_}")

    def test_add_sub_with_multiples_of_q(self):
        for mod_ in BIG_MODS:
            for k in [0, 1, 2, 10, 123456]:
                num1 = rand_big(512) + k * mod_
                num2 = -rand_big(512) - k * mod_
                self.assertEqual(mod.add_mod(num1, num2, mod_), ((num1 % mod_) + (num2 % mod_)) % mod_)
                self.assertEqual(mod.sub_mod(num1, num2, mod_), ((num1 % mod_) - (num2 % mod_)) % mod_)

    def test_mul_mod_large(self):
        for mod_ in BIG_MODS:
            for _ in range(100):
                num1 = rand_big(1024)
                num2 = rand_big(1024)
                ref = ((num1 % mod_) * (num2 % mod_)) % mod_
                got = mod.mul_mod(num1, num2, mod_)
                self.assertEqual(got, ref, f"mul_mod failed for q={mod_}")

    def test_mul_with_multiples_and_negatives(self):
        for mod_ in BIG_MODS:
            for k in [0, 1, 2, 5, 1000]:
                num1 = (rand_big(512) - rand_big(256)) + k * mod_
                num2 = -(rand_big(512)) - k * mod_
                ref = ((num1 % mod_) * (num2 % mod_)) % mod_
                got = mod.mul_mod(num1, num2, mod_)
                self.assertEqual(got, ref)

    def test_inv_mod_matches_builtin_pow(self):
        for mod_ in BIG_MODS:
            for _ in range(50):
                num1 = rand_big(1024) % mod_
                if num1 == 0:
                    num1 = 1
                self.assertEqual(math.gcd(num1, mod_), 1)
                inv = mod.inv_mod(num1, mod_)
                ref = pow(num1, -1, mod_)
                self.assertEqual(inv, ref, f"inv_mod != pow for q={mod_}")
                self.assertEqual((num1 * inv) % mod_, 1)

    def test_inv_mod_non_invertible_raises(self):
        for mod_ in BIG_MODS:
            with self.assertRaises(Exception):
                mod.inv_mod(0, mod_)
            with self.assertRaises(Exception):
                mod.inv_mod(mod_ * 123456789, mod_)

    def test_distributivity(self):
        for mod_ in BIG_MODS:
            for _ in range(50):
                num1 = rand_big(512)
                num2 = rand_big(512)
                num3 = rand_big(512)
                left = mod.mul_mod(mod.add_mod(num1, num2, mod_), num3, mod_)
                right = mod.add_mod(mod.mul_mod(num1, num3, mod_), mod.mul_mod(num2, num3, mod_), mod_)
                self.assertEqual(left, right)

    def test_associativity_add_mul(self):
        for mod_ in BIG_MODS:
            for _ in range(50):
                num1 = rand_big(256)
                num2 = rand_big(256)
                num3 = rand_big(256)
                self.assertEqual(mod.add_mod(mod.add_mod(num1, num2, mod_), num3, mod_),
                                 mod.add_mod(num1, mod.add_mod(num2, num3, mod_), mod_))
                self.assertEqual(mod.mul_mod(mod.mul_mod(num1, num2, mod_), num3, mod_),
                                 mod.mul_mod(num1, mod.mul_mod(num2, num3, mod_), mod_))

if __name__ == "__main__":
    unittest.main(verbosity=2)
