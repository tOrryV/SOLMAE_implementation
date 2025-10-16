from __future__ import annotations
from typing import List
import os, hmac, hashlib


def random_bytes_sys(n):
    return os.urandom(n)


def random_uint_below(mod):
    if mod <= 0:
        raise ValueError("q must be > 0")
    nbytes = (mod.bit_length() + 7) // 8
    limit = (1 << (8*nbytes)) - ((1 << (8*nbytes)) % mod)
    while True:
        x = int.from_bytes(os.urandom(nbytes), "big")
        if x < limit:
            return x % mod


class HMACDRBG:
    def __init__(self, seed):
        self._K = b"\x00" * 32
        self._V = b"\x01" * 32
        self._update(seed)

    def _hmac(self, key, data):
        return hmac.new(key, data, hashlib.sha256).digest()

    def _update(self, provided_data):
        self._K = self._hmac(self._K, self._V + b"\x00" + (provided_data or b""))
        self._V = self._hmac(self._K, self._V)
        if provided_data is not None:
            self._K = self._hmac(self._K, self._V + b"\x01" + provided_data)
            self._V = self._hmac(self._K, self._V)

    def generate(self, n, additional_data=None):
        if additional_data:
            self._update(additional_data)
        out = bytearray()
        while len(out) < n:
            self._V = self._hmac(self._K, self._V)
            out += self._V
        self._update(additional_data)
        return bytes(out[:n])

def uniform_mod_q(deg, mod):
    return [random_uint_below(mod) for _ in range(deg)]

def uniform_small(deg, bound):
    if bound < 0:
        raise ValueError("bound must be >= 0")
    width = 2*bound + 1
    res = []
    while len(res) < deg:
        x = random_uint_below(width)
        res.append(x - bound)
    return res

def sample_cbd(bytes_in, eta):
    if eta <= 0:
        raise ValueError("eta must be > 0")
    out: List[int] = []
    bitbuf = int.from_bytes(bytes_in, "big")
    nbits = len(bytes_in) * 8
    pos = 0
    while pos + 2*eta <= nbits:
        a = 0
        b = 0
        for _ in range(eta):
            a += (bitbuf >> (nbits - 1 - pos)) & 1
            pos += 1
        for _ in range(eta):
            b += (bitbuf >> (nbits - 1 - pos)) & 1
            pos += 1
        out.append(a - b)
    return out

def sample_cbd_random(deg, eta):
    nbits = deg * 2 * eta
    nbytes = (nbits + 7) // 8
    buf = os.urandom(nbytes)
    coeffs = sample_cbd(buf, eta)
    while len(coeffs) < deg:
        buf = os.urandom(32)
        coeffs.extend(sample_cbd(buf, eta))
    return coeffs[:deg]

def expand_seed_to_mod_q(seed, deg, mod):
    drbg = HMACDRBG(seed)
    res = []
    nbytes = (mod.bit_length() + 7) // 8
    limit = (1 << (8*nbytes)) - ((1 << (8*nbytes)) % mod)
    while len(res) < deg:
        chunk = drbg.generate(nbytes)
        x = int.from_bytes(chunk, "big")
        if x < limit:
            res.append(x % mod)
    return res
