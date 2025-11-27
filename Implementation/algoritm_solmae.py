import hashing
import rng
from poly import Poly
from params import n, k, q, eta, d, SEED_LEN



def _poly_from_rho(rho, i, j):
    Poly.ring(q, n)
    buf = hashing.XOF_shake128(
        parts=[rho, i.to_bytes(2, "big"), j.to_bytes(2, "big")],
        outlen=2 * n,
        domain=hashing.DOM_SEP["H_seed"],
    )
    coeffs = []
    for t in range(0, 2 * n, 2):
        v = (buf[t] << 8) | buf[t + 1]
        coeffs.append(v % q)
    return Poly(coeffs)


def _expand_matrix_A(rho):
    return [[_poly_from_rho(rho, i, j) for j in range(k)] for i in range(k)]


def _highbits_bytes(poly_obj):
    return bytes((x >> d) for x in poly_obj.a)


def _scale_high_to_ring(hi_bytes):
    scaled = [((h << d) % q) for h in hi_bytes]
    return Poly(scaled)


def keygen_solmae():
    seed = rng.random_bytes_sys(SEED_LEN)
    rho = hashing.H_seed_expand(seed, SEED_LEN)

    Poly.ring(q, n)

    s = [Poly(hashing.H_to_small_poly(rng.random_bytes_sys(SEED_LEN), n, eta)) for _ in range(k)]
    e = [Poly(hashing.H_to_small_poly(rng.random_bytes_sys(SEED_LEN), n, eta)) for _ in range(k)]

    A = _expand_matrix_A(rho)

    t = [Poly.zero() for _ in range(k)]
    for i in range(k):
        for j in range(k):
            t[i] = t[i] + A[i][j]._mul_poly(s[j])

    t1, t0 = [], []
    for ti in t:
        hi = [(x >> d) for x in ti.a]
        lo = [x - ((h << d) % q) for x, h in zip(ti.a, hi)]
        t1.append(bytes(hi))
        t0.append(lo)

    pk_bytes = rho + b"".join(t1)
    tr = hashing.H_pk_bind(pk_bytes)

    pk = (rho, t1)
    sk = (s, e, t0, tr, pk)
    return pk, sk



def sign_solmae(sk, msg):

    s, e, t0, tr, pk = sk
    rho, t1 = pk

    Poly.ring(q, n)

    y = [Poly(hashing.H_to_small_poly(rng.random_bytes_sys(SEED_LEN), n, eta)) for _ in range(k)]

    A = _expand_matrix_A(rho)

    w = [Poly.zero() for _ in range(k)]
    for i in range(k):
        for j in range(k):
            w[i] = w[i] + A[i][j]._mul_poly(y[j])

    w1 = [_highbits_bytes(wi) for wi in w]

    z = [yi.a[:] for yi in y]

    z_bytes = b"".join(int(x).to_bytes(2, "big") for zi in z for x in zi)

    c = hashing.H_challenge(tr, msg + b"".join(w1) + z_bytes, 32)

    return z, c, w1


def verify_solmae(pk, msg, sig):
    rho, t1 = pk
    z, c, w1 = sig

    Poly.ring(q, n)

    if not isinstance(z, list) or len(z) != k:
        return False
    for zi in z:
        if not isinstance(zi, list) or len(zi) != n:
            return False
        if any(not (0 <= int(x) < q) for x in zi):
            return False

    if not isinstance(w1, list) or len(w1) != k:
        return False
    for w1i in w1:
        if not isinstance(w1i, (bytes, bytearray)) or len(w1i) != n:
            return False

    if not isinstance(c, (bytes, bytearray)):
        return False

    tr_check = hashing.H_pk_bind(rho + b"".join(t1))

    z_bytes = b"".join(int(x).to_bytes(2, "big") for zi in z for x in zi)

    c_check = hashing.H_challenge(tr_check, msg + b"".join(w1) + z_bytes, 32)

    return c_check == c
