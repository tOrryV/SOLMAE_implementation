from __future__ import annotations
from typing import List

from modular import add_mod, sub_mod, mul_mod, inv_mod

def is_power_of_two(n):
    return n > 0 and (n & (n - 1)) == 0


def factorize(n):
    fs: List[int] = []
    d = 2
    x = n
    while d * d <= x:
        if x % d == 0:
            fs.append(d)
            while x % d == 0:
                x //= d
        d += 1 if d == 2 else 2
    if x > 1:
        fs.append(x)
    return fs

def primitive_root_for(mod, order):
    assert (mod - 1) % order == 0, "order must divide q-1"
    phi = mod - 1
    primes = factorize(phi)

    g = 2
    while g < mod:
        ok = True
        for p in primes:
            if pow(g, phi // p, mod) == 1:
                ok = False
                break
        if ok:
            break
        g += 1

    h = pow(g, phi // order, mod)
    for p in factorize(order):
        if pow(h, order // p, mod) == 1:
            k = 2
            while pow(h, order // p, mod) == 1:
                h = pow(g, (phi // order) * k, mod)
                k += 1
    assert pow(h, order, mod) == 1
    return h


def bitrev_permutation(n):
    assert is_power_of_two(n)
    m = n.bit_length() - 1
    perm = [0] * n
    for i in range(n):
        r = 0
        x = i
        for _ in range(m):
            r = (r << 1) | (x & 1)
            x >>= 1
        perm[i] = r
    return perm


def precompute_roots(mod, n, psi):
    assert is_power_of_two(n)
    omega = pow(psi, 2, mod)
    roots: List[int] = []
    roots_inv: List[int] = []

    length = 2
    while length <= n:
        half = length // 2
        w_len = pow(omega, n // length, mod)
        w = 1
        for _ in range(half):
            roots.append(w)
            w = mul_mod(w, w_len, mod)
        w_len_inv = inv_mod(w_len, mod)
        w = 1
        for _ in range(half):
            roots_inv.append(w)
            w = mul_mod(w, w_len_inv, mod)
        length <<= 1

    return roots, roots_inv, bitrev_permutation(n)


def precompute_twists(mod, n, psi):
    tw_fwd = [1] * n
    tw_inv = [1] * n
    inv_psi = inv_mod(psi, mod)
    for i in range(1, n):
        tw_fwd[i] = mul_mod(tw_fwd[i - 1], psi, mod)
        tw_inv[i] = mul_mod(tw_inv[i - 1], inv_psi, mod)
    return tw_fwd, tw_inv


def _bitrev_shuffle(a, perm):
    n = len(a)
    for i in range(n):
        j = perm[i]
        if j > i:
            a[i], a[j] = a[j], a[i]


def ntt_inplace(a, mod, roots, bitrev):
    n = len(a)
    _bitrev_shuffle(a, bitrev)
    stage = 0
    length = 2
    while length <= n:
        half = length // 2
        layer = roots[stage: stage + half]
        for i in range(0, n, length):
            for j in range(half):
                u = a[i + j]
                v = mul_mod(a[i + j + half], layer[j], mod)
                a[i + j] = add_mod(u, v, mod)
                a[i + j + half] = sub_mod(u, v, mod)
        stage += half
        length <<= 1


def intt_inplace(a, mod, roots_inv, bitrev):
    n = len(a)
    _bitrev_shuffle(a, bitrev)
    stage = 0
    length = 2
    while length <= n:
        half = length // 2
        layer = roots_inv[stage: stage + half]
        for i in range(0, n, length):
            for j in range(half):
                u = a[i + j]
                v = mul_mod(a[i + j + half], layer[j], mod)
                a[i + j] = add_mod(u, v, mod)
                a[i + j + half] = sub_mod(u, v, mod)
        stage += half
        length <<= 1
    inv_n = inv_mod(n % mod, mod)
    for i in range(n):
        a[i] = mul_mod(a[i], inv_n, mod)


def negacyclic_convolution(num1, num2, mod, roots, roots_inv, bitrev, tw_fwd, tw_inv):
    n = len(num1)
    A = [mul_mod(ai % mod, tw_fwd[i], mod) for i, ai in enumerate(num1)]
    B = [mul_mod(bi % mod, tw_fwd[i], mod) for i, bi in enumerate(num2)]
    ntt_inplace(A, mod, roots, bitrev)
    ntt_inplace(B, mod, roots, bitrev)
    C = [mul_mod(A[i], B[i], mod) for i in range(n)]
    intt_inplace(C, mod, roots_inv, bitrev)
    C = [mul_mod(C[i], tw_inv[i], mod) for i in range(n)]
    return C


def poly_mul_rq_ntt(num1, num2, mod, deg):
    assert len(num1) == deg and len(num2) == deg, "polys must have length d"
    assert is_power_of_two(deg), "d must be a power of two"
    assert (mod - 1) % (2 * deg) == 0, "2d must divide (q-1)"

    psi = primitive_root_for(mod, 2 * deg)
    roots, roots_inv, bitrev = precompute_roots(mod, deg, psi)
    tw_fwd, tw_inv = precompute_twists(mod, deg, psi)

    return negacyclic_convolution(num1, num2, mod, roots, roots_inv, bitrev, tw_fwd, tw_inv)
