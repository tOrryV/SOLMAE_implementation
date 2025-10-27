from __future__ import annotations
from dataclasses import dataclass
from typing import List, Optional, Dict
import math, os, struct

from rng import HMACDRBG, sample_cbd_random
from ntt import primitive_root_for, precompute_roots, precompute_twists
from cfft import precompute_twiddles, bitrev_permutation as cfft_bitrev


@dataclass
class SampleParams:
    n: int
    q: int
    sigma: float
    tailcut: float = 10.0
    eta: int = 2
    use_ntt: bool = True
    use_fft: bool = True


@dataclass
class NTTPlan:
    q: int
    n: int
    psi: int
    roots: List[int]
    roots_inv: List[int]
    bitrev: List[int]
    tw_fwd: List[int]
    tw_inv: List[int]


@dataclass
class CFFTPlan:
    n: int
    W: List[complex]
    Winv: List[complex]
    bitrev: List[int]


@dataclass
class CDTTable:
    sigma: float
    cutoff: int
    cdf_scaled: List[int]


@dataclass
class SamplerState:
    params: SampleParams
    ntt: Optional[NTTPlan]
    cfft: Optional[CFFTPlan]
    cdt: CDTTable
    cache: Dict[str, object]


def make_ntt_plan(q, n):
    assert (q - 1) % (2 * n) == 0, "q-1 must be divisible by 2*n for negacyclic NTT"
    psi = primitive_root_for(q, 2 * n)
    roots, roots_inv, bitrev = precompute_roots(q, n, psi)
    tw_fwd, tw_inv = precompute_twists(q, n, psi)
    return NTTPlan(q, n, psi, roots, roots_inv, bitrev, tw_fwd, tw_inv)


def make_cfft_plan(n):
    W, Winv = precompute_twiddles(n)
    br = cfft_bitrev(n)
    return CFFTPlan(n, W, Winv, br)


def gaussian_cdt_build(sigma, tailcut=10.0):
    cutoff = max(1, int(math.ceil(sigma * tailcut)))
    scale = 1 << 64

    w = [math.exp(-(x*x) / (2*sigma*sigma)) for x in range(cutoff+1)]

    Z = w[0] + 2 * sum(w[1:])

    cdf_scaled = []
    cum = 0.0
    for x in range(cutoff+1):
        p = (w[x] / Z) if x == 0 else (2*w[x] / Z)
        cum += p
        cdf_scaled.append(int(min(scale-1, max(0, round(cum * scale)))))

    for i in range(1, len(cdf_scaled)):
        cdf_scaled[i] = max(cdf_scaled[i], cdf_scaled[i-1])

    return CDTTable(sigma, cutoff, cdf_scaled)



def gaussian_cdt_sample_u64(table, u64):
    cutoff, c = table.cutoff, table.cdf_scaled
    lo, hi = 0, cutoff
    while lo < hi:
        mid = (lo + hi) // 2
        if c[mid] >= u64:
            hi = mid
        else:
            lo = mid + 1
    x = lo
    sign = u64 & 1
    return x if sign == 0 else -x



def rand64_os():
    return struct.unpack(">Q", os.urandom(8))[0]


def gaussian_vec_from_rng(table, n, rand64=rand64_os):
    return [gaussian_cdt_sample_u64(table, rand64()) for _ in range(n)]


def precompute_for_sample(params):
    plan_ntt = make_ntt_plan(params.q, params.n) if params.use_ntt else None
    plan_fft = make_cfft_plan(params.n) if params.use_fft else None
    cdt = gaussian_cdt_build(params.sigma, params.tailcut)
    return SamplerState(params, plan_ntt, plan_fft, cdt, cache={})


def next_u64_from_state(state):
    drbg = state.cache.get("drbg", None)
    if drbg is not None:
        return struct.unpack(">Q", drbg.generate(8))[0]
    return rand64_os()


def attach_drbg(state, seed):
    state.cache["drbg"] = HMACDRBG(seed)


def sample_noise_poly(state):
    return gaussian_vec_from_rng(state.cdt, state.params.n, lambda: next_u64_from_state(state))


def sample_noise_poly_modq(state):
    q = state.params.q
    return [(x % q) for x in sample_noise_poly(state)]


def sample_cbd_poly(n, eta):
    return sample_cbd_random(n, eta)
