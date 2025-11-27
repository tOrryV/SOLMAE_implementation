def add_mod(num1, num2, mod):
    return ((num1 % mod) + (num2 % mod)) % mod


def sub_mod(num1, num2, mod):
    return ((num1 % mod) - (num2 % mod)) % mod


def mul_mod(num1, num2, mod):
    return ((num1 % mod) * (num2 % mod)) % mod


def inv_mod(num1, mod):
    num1 %= mod
    if num1 == 0:
        raise ValueError("inv_mod: a ≡ 0 (mod q)")
    t, new_t = 0, 1
    r, new_r = mod, num1
    while new_r:
        qout = r // new_r
        t, new_t = new_t, t - qout * new_t
        r, new_r = new_r, r - qout * new_r
    if r != 1:
        raise ValueError("inv_mod: a is not inverse")
    return t % mod