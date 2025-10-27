from __future__ import annotations
import math


def z_negacyclic_mul(num1, num2):
    n = len(num1)
    assert len(num2) == n
    acc = [0] * n
    for i, ai in enumerate(num1):
        if ai == 0:
            continue
        for j, bj in enumerate(num2):
            if bj == 0:
                continue
            k = i + j
            prod = ai * bj
            if k < n:
                acc[k] += prod
            else:
                acc[k - n] -= prod
    return acc


def z_add(num1, num2):
    assert len(num1) == len(num2)
    return [num1[i] + num2[i] for i in range(len(num1))]


def z_sub(num1, num2):
    assert len(num1) == len(num2)
    return [num1[i] - num2[i] for i in range(len(num1))]


def z_scalar_mul(num, k):
    return [k * x for x in num]


def z_centered_mod_q(num, mod):
    half = mod // 2
    out = []
    for x in num:
        r = x % mod
        out.append(r if r <= half else r - mod)
    return out


def z_round_div(num, k):
    if k == 0:
        raise ZeroDivisionError("division by zero")
    res = []
    for x in num:
        if k > 0:
            r = math.floor((x + k/2) / k)
        else:
            r = math.floor((x + k/2) / k)
        res.append(int(r))
    return res


def deg(num):
    return len(num)


def split_even_odd(num):
    n = len(num)
    assert n % 2 == 0, "length must be even to split"
    ae = num[0::2]
    ao = num[1::2]
    return ae, ao


def merge_even_odd(ae, ao):
    assert len(ae) == len(ao)
    n2 = len(ae)
    n = 2 * n2
    out = [0] * n
    out[0::2] = ae
    out[1::2] = ao
    return out


def proj_down(num):
    ae, ao = split_even_odd(num)
    return [(ae[i] + ao[i]) for i in range(len(ae))]


def lift_up(u, v):
    assert len(u) == len(v)
    n2 = len(u)
    out = [0] * (2 * n2)
    out[0::2] = u
    out[1::2] = v
    return out


def closest_rounding_div_q(num, mod):
    if mod == 0:
        raise ZeroDivisionError("division by zero")
    return [int(math.floor((x + mod/2) / mod)) for x in num]


def reduce_target(F, G, f, g, mod):
    assert len(F) == len(G) == len(f) == len(g)
    num = 0
    den = 0
    for i in range(len(F)):
        num += f[i] * F[i] + g[i] * G[i]
        den += F[i] * F[i] + G[i] * G[i]
    if den == 0:
        return f, g
    k = int(round(num / den))
    if k != 0:
        f = [f[i] - k * F[i] for i in range(len(F))]
        g = [g[i] - k * G[i] for i in range(len(G))]
    return z_centered_mod_q(f, mod), z_centered_mod_q(g, mod)


def _xgcd(num1, num2):
    if num2 == 0:
        return (abs(num1), 1 if num1 > 0 else -1, 0)
    g, x1, y1 = _xgcd(num2, num1 % num2)

    return (g, y1, x1 - (num1 // num2) * y1)


def basecase_solve(F, G, mod):
    assert len(F) == len(G) == 1, "basecase only for n=1"
    F0, G0 = int(F[0]), int(G[0])
    if F0 == 0 and G0 == 0:
        raise ValueError("Trivial basecase: F=G=0")
    g, x, y = _xgcd(G0, F0)
    if mod % g != 0:
        raise ValueError("No integer solution to f*G - g*F = q (gcd does not divide q)")
    t = mod // g
    f = x * t
    g_coef = -y * t
    return [f], [g_coef]


def ntru_solve_pp(F, G, mod):
    n = len(F)
    assert n == len(G), "F and G must have same length"
    if n == 1:
        return basecase_solve(F, G, mod)

    Fe, Fo = split_even_odd(F)
    Ge, Go = split_even_odd(G)

    fe, ge = ntru_solve_pp(Fe, Ge, mod)

    f_lift = lift_up(fe, [0]*len(fe))
    g_lift = lift_up(ge, [0]*len(ge))

    f_red, g_red = reduce_target(F, G, f_lift, g_lift, mod)

    raise NotImplementedError("PP-recursion merge step is not implemented yet (skeleton returned).")


def poly_to_z(p_coeffs, n):
    c = list(int(x) for x in p_coeffs)
    if len(c) > n:
        base = c[:n]
        for i in range(n, len(c)):
            if c[i]:
                base[i - n] -= c[i]
        c = base
    if len(c) < n:
        c += [0] * (n - len(c))
    return c


def z_to_poly(a, mod):
    return [x % mod for x in a]


def ntru_solve_poly(F_poly, G_poly, mod):
    from poly import Poly
    n = Poly._deg
    F = poly_to_z(F_poly.to_list(), n)
    G = poly_to_z(G_poly.to_list(), n)
    f, g = ntru_solve_pp(F, G, mod)
    return Poly.from_list([x % mod for x in f]), Poly.from_list([x % mod for x in g])
