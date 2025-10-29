from __future__ import annotations
import os, math, struct
from sample_precomp import gaussian_cdt_build
from cfft import hadamard_product, add_complex, sub_complex


_CDT_CACHE: dict[float, object] = {}
_BM_CACHE = None

def _get_cdt_table(sigma, tailcut=10.0):
    key = (round(float(sigma), 12), round(float(tailcut), 2))
    tab = _CDT_CACHE.get(key)
    if tab is None:
        tab = gaussian_cdt_build(float(sigma), tailcut)
        _CDT_CACHE[key] = tab
    return tab

def _rand_u64():
    return struct.unpack(">Q", os.urandom(8))[0]

def _uniform01_from_u64(u):
    return (u + 0.5) / (1 << 64)

def _bernoulli_exp(log_p):
    if log_p >= 0.0:
        return True
    u = _uniform01_from_u64(_rand_u64())
    return u <= math.exp(log_p)

def _gaussian01():
    global _BM_CACHE
    if _BM_CACHE is not None:
        z = _BM_CACHE
        _BM_CACHE = None
        return z

    u1 = _uniform01_from_u64(_rand_u64())
    u2 = _uniform01_from_u64(_rand_u64())
    u1 = max(u1, 2.0**-64)
    R = math.sqrt(-2.0 * math.log(u1))
    theta = 2.0 * math.pi * u2
    z0 = R * math.cos(theta)
    z1 = R * math.sin(theta)
    _BM_CACHE = z1
    return z0


def n_sampler(d):
    u1 = _uniform01_from_u64(_rand_u64())
    u2 = _uniform01_from_u64(_rand_u64())
    u1 = max(u1, 2.0**-64)
    R = math.sqrt(-2.0 * math.log(u1))
    theta = 2.0 * math.pi * u2
    return R * math.cos(theta), R * math.sin(theta)


def z_sampler(center, sigma, tailcut=10.0):
    y = center + sigma * _gaussian01()
    return int(round(y))


def peikert_sampler(t, Sigma, eta, d):
    assert len(t) == len(Sigma) == d, "peikert_sampler: size mismatch"
    out = [0] * d
    i = 0
    while i < d:
        z0, z1 = n_sampler(2)
        u0 = float(t[i]) + float(Sigma[i]) * z0
        out[i] = z_sampler(u0, float(eta))
        i += 1
        if i < d:
            u1 = float(t[i]) + float(Sigma[i]) * z1
            out[i] = z_sampler(u1, float(eta))
            i += 1
    return out


def _proj_beta_K(c_fft, beta_fft):
    c1, c2 = c_fft
    if isinstance(beta_fft, tuple):
        b1, b2 = beta_fft
        t1 = hadamard_product(b1, c1)
        t2 = hadamard_product(b2, c2)
        return [t1[i] + t2[i] for i in range(len(t1))]
    else:
        t1 = hadamard_product(beta_fft, c1)
        return t1


def _hadamard_int_with_pair(z, pair_fft):
    a, b = pair_fft
    assert len(z) == len(a) == len(b), "hadamard_int_with_pair: size mismatch"
    za = [complex(int(z[i])) * a[i] for i in range(len(z))]
    zb = [complex(int(z[i])) * b[i] for i in range(len(z))]
    return za, zb


def sample(c_fft, sk, params):
    d = getattr(sk, "d", None)
    assert d is not None, "sk.d must be set"

    t2_fft = _proj_beta_K(c_fft, sk.beta2_fft)
    t2 = [float(z.real) for z in t2_fft]
    Sigma2 = [float(s.real) if isinstance(s, complex) else float(s) for s in sk.Sigma2_fft]
    z2 = peikert_sampler(t2, Sigma2, float(params.eta), d)

    z2b2_a, z2b2_b = _hadamard_int_with_pair(z2, sk.b2_tilde_fft)
    c1_new = sub_complex(c_fft[0], z2b2_a)
    c2_new = sub_complex(c_fft[1], z2b2_b)
    c_fft = (c1_new, c2_new)

    t1_fft = _proj_beta_K(c_fft, sk.beta1_fft)
    t1 = [float(z.real) for z in t1_fft]
    Sigma1 = [float(s.real) if isinstance(s, complex) else float(s) for s in sk.Sigma1_fft]
    z1 = peikert_sampler(t1, Sigma1, float(params.eta), d)

    z1b1_a, z1b1_b = _hadamard_int_with_pair(z1, sk.b1_fft)
    z2b2_a, z2b2_b = _hadamard_int_with_pair(z2, sk.b2_fft)
    v1 = add_complex(z1b1_a, z2b2_a)
    v2 = add_complex(z1b1_b, z2b2_b)
    return v1, v2
