import cmath
import math


def is_power_of_two(n):
    return n > 0 and (n & (n - 1)) == 0


def bitrev_permutation(n):
    assert is_power_of_two(n), "bitrev_permutation: n must be power of two"
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


def precompute_twiddles(n):
    assert is_power_of_two(n), "precompute_twiddles: n must be power of two"
    tau = 2.0 * math.pi
    W = [complex(1.0, 0.0)] * n
    Winv = [complex(1.0, 0.0)] * n
    for k in range(1, n):
        angle = -tau * k / n
        W[k] = cmath.cos(angle) + 1j * cmath.sin(angle)
        Winv[k] = W[k].conjugate()
    return W, Winv


def _bitrev_shuffle(a, perm):
    n = len(a)
    for i in range(n):
        j = perm[i]
        if j > i:
            a[i], a[j] = a[j], a[i]


def fft_inplace(a, W, bitrev):
    n = len(a)
    assert is_power_of_two(n), "fft_inplace: n must be power of two"
    _bitrev_shuffle(a, bitrev)

    length = 2
    while length <= n:
        half = length // 2
        step = n // length
        for start in range(0, n, length):
            w = 1+0j
            for j in range(half):
                u = a[start + j]
                v = a[start + j + half] * w
                a[start + j] = u + v
                a[start + j + half] = u - v
                w *= W[step]
        length <<= 1

def ifft_inplace(A, Winv, bitrev):
    n = len(A)
    assert is_power_of_two(n), "ifft_inplace: n must be power of two"
    _bitrev_shuffle(A, bitrev)

    length = 2
    while length <= n:
        half = length // 2
        step = n // length
        for start in range(0, n, length):
            w = 1+0j
            for j in range(half):
                u = A[start + j]
                v = A[start + j + half] * w
                A[start + j] = u + v
                A[start + j + half] = u - v
                w *= Winv[step]
        length <<= 1

    inv_n = 1.0 / n
    for i in range(n):
        A[i] *= inv_n


def fft(x):
    a = [complex(z) for z in x]
    n = len(a)
    assert is_power_of_two(n), "fft: length must be power of two"
    W, _Winv = precompute_twiddles(n)
    perm = bitrev_permutation(n)
    fft_inplace(a, W, perm)
    return a


def ifft(X):
    A = [complex(z) for z in X]
    n = len(A)
    assert is_power_of_two(n), "ifft: length must be power of two"
    _W, Winv = precompute_twiddles(n)
    perm = bitrev_permutation(n)
    ifft_inplace(A, Winv, perm)
    return A


def hadamard_product(A, B):
    A = list(A)
    B = list(B)
    assert len(A) == len(B), "hadamard_product: size mismatch"
    return [A[i] * B[i] for i in range(len(A))]


def scale_complex(X, alpha):
    return [alpha * z for z in X]


def add_complex(A, B):
    A = list(A)
    B = list(B)
    assert len(A) == len(B), "add_complex: size mismatch"
    return [A[i] + B[i] for i in range(len(A))]


def sub_complex(A, B):
    A = list(A)
    B = list(B)
    assert len(A) == len(B), "sub_complex: size mismatch"
    return [A[i] - B[i] for i in range(len(A))]


def fft_real(x):
    return fft([complex(float(t), 0.0) for t in x])


def ifft_to_real(X):
    y = ifft(X)
    return [float(z.real) for z in y]


def max_abs_diff(A, B):
    A = list(A)
    B = list(B)
    assert len(A) == len(B), "max_abs_diff: size mismatch"
    return max(abs(A[i] - B[i]) for i in range(len(A))) if A else 0.0


def is_close_vec(A, B, tol=1e-9):
    return max_abs_diff(A, B) <= tol
