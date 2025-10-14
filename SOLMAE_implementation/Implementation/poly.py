from typing import Iterable, List
import random

from modular import add_mod, sub_mod, mul_mod


class Poly:
    _mod = None
    _deg = None

    @classmethod
    def ring(cls, mod, deg):
        assert mod > 1 and deg >= 1
        Poly._mod = int(mod)
        Poly._deg = int(deg)

    def __init__(self, coeffs: Iterable[int]):
        if Poly._mod is None or Poly._deg is None:
            raise RuntimeError("First, call Poly.ring(q, d)")
        mod, deg = Poly._mod, Poly._deg
        c = list(int(x) for x in coeffs)
        if len(c) > deg:
            base = c[:deg]
            for idx in range(deg, len(c)):
                if c[idx]:
                    base[idx - deg] -= c[idx]
            c = base
        if len(c) < deg:
            c += [0] * (deg - len(c))
        self.a = [ci % mod for ci in c]

    def __add__(self, other):
        if not isinstance(other, Poly):
            return NotImplemented
        if Poly._mod != other._mod or Poly._deg != other._deg:
            raise ValueError("Various ring parameters for operands")
        mod, deg = Poly._mod, Poly._deg
        return Poly([(self.a[i] + other.a[i]) % mod for i in range(deg)])

    def __sub__(self, other):
        if not isinstance(other, Poly):
            return NotImplemented
        if Poly._mod != other._mod or Poly._deg != other._deg:
            raise ValueError("Various ring parameters for operands")
        mod, deg = Poly._mod, Poly._deg
        return Poly([(self.a[i] - other.a[i]) % mod for i in range(deg)])

    @staticmethod
    def _reduce_to_length_and_mod(c, mod, deg):
        if len(c) > deg:
            base = c[:deg] + [0] * 0
            for idx in range(deg, len(c)):
                if c[idx] == 0:
                    continue
                j = idx - deg
                base[j] -= c[idx]
            c = base
        if len(c) < deg:
            c = c + [0] * (deg - len(c))
        for i in range(deg):
            c[i] %= mod
        return c

    @classmethod
    def zero(cls):
        return cls([0] * cls._deg)

    @classmethod
    def one(cls):
        a = [0] * cls._deg
        a[0] = 1 % cls._mod
        return cls(a)

    @classmethod
    def from_list(cls, coeffs):
        return cls(coeffs)

    @classmethod
    def random(cls):
        mod, deg = cls._mod, cls._deg
        return cls([random.randrange(0, mod) for _ in range(deg)])

    def to_list(self) -> List[int]:
        return list(self.a)

    def __repr__(self):
        return f"Poly({self.a!r}; q={Poly._mod}, d={Poly._deg})"

    def __str__(self):
        terms = []
        for i, c in enumerate(self.a):
            if c == 0:
                continue
            if i == 0:
                terms.append(f"{c}")
            elif i == 1:
                terms.append(f"{c}·X")
            else:
                terms.append(f"{c}·X^{i}")
        return "0" if not terms else " + ".join(terms)

    def _binop(self, other, op):
        if not isinstance(other, Poly):
            return NotImplemented
        if Poly._mod != other._mod or Poly._deg != other._deg:
            raise ValueError("Various ring parameters for operands")
        mod, deg = Poly._mod, Poly._deg
        res = [0] * deg
        for i in range(deg):
            res[i] = op(self.a[i], other.a[i], mod)
        return Poly(res)

    def __neg__(self):
        mod = Poly._mod
        return Poly([(-x) % mod for x in self.a])

    def _scalar_mul(self, c):
        mod = Poly._mod
        c %= mod
        return Poly([(ai % mod) * c % mod for ai in self.a])

    def __rmul__(self, c):
        if isinstance(c, int):
            return self._scalar_mul(c)
        return NotImplemented

    def __mul__(self, other):
        if isinstance(other, int):
            return self._scalar_mul(other)
        if not isinstance(other, Poly):
            return NotImplemented
        mod = Poly._mod
        deg = Poly._deg
        a = self.a
        b = other.a
        acc = [0] * deg
        for i, ai in enumerate(a):
            if ai == 0:
                continue
            for j, bj in enumerate(b):
                if bj == 0:
                    continue
                k = i + j
                prod = (ai * bj) % mod
                if k < deg:
                    acc[k] = (acc[k] + prod) % mod
                else:
                    acc[k - deg] = (acc[k - deg] - prod) % mod
        return Poly(acc)

    def _mul_poly(self, other: "Poly"):
        if Poly._mod != other._mod or Poly._deg != other._deg:
            raise ValueError("Various ring parameters for operands")
        mod, deg = Poly._mod, Poly._deg
        acc = [0] * deg
        A, B = self.a, other.a
        for i in range(deg):
            ai = A[i]
            if ai == 0:
                continue
            for j in range(deg):
                bj = B[j]
                if bj == 0:
                    continue
                k = i + j
                prod = mul_mod(ai, bj, mod)
                if k < deg:
                    acc[k] = add_mod(acc[k], prod, mod)
                else:
                    kd = k - deg
                    acc[kd] = sub_mod(acc[kd], prod, mod)
        return Poly(acc)

    def __eq__(self, other: object):
        if not isinstance(other, Poly):
            return False
        return Poly._mod == other._mod and Poly._deg == other._deg and self.a == other.a
