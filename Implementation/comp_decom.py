def _centered_mod_q(a, q):
    a = a % q
    half_down = q // 2
    if a > half_down:
        a -= q
    return a


def _pack_bits(values, bits_per_val, total_bits):
    acc = 0
    acc_bits = 0
    out = bytearray()
    mask = (1 << bits_per_val) - 1
    bits_written = 0

    for v in values:
        v &= mask
        acc |= (v << acc_bits)
        acc_bits += bits_per_val
        while acc_bits >= 8 and bits_written + 8 <= total_bits:
            out.append(acc & 0xFF)
            acc >>= 8
            acc_bits -= 8
            bits_written += 8

    if bits_written < total_bits and acc_bits > 0:
        remain = min(8, total_bits - bits_written)
        out.append(acc & ((1 << remain) - 1))
        bits_written += remain

    need_bytes = (total_bits + 7) // 8
    if len(out) < need_bytes:
        out.extend(b"\x00" * (need_bytes - len(out)))
    elif len(out) > need_bytes:
        out = out[:need_bytes]
    return bytes(out)


def _unpack_bits(data, bits_per_val, count):
    acc = 0
    acc_bits = 0
    idx = 0
    out = []
    mask = (1 << bits_per_val) - 1
    data_len = len(data)

    while len(out) < count:
        while acc_bits < bits_per_val and idx < data_len:
            acc |= data[idx] << acc_bits
            acc_bits += 8
            idx += 1
        if acc_bits < bits_per_val:
            v = acc & mask
            acc = 0
            acc_bits = 0
            out.append(v)
            continue
        v = acc & mask
        acc >>= bits_per_val
        acc_bits -= bits_per_val
        out.append(v)
    return out


def compress(s1, slen, q=12289):
    d = len(s1)
    if d <= 0:
        return b""
    b = slen // d
    if b <= 0:
        raise ValueError("compress: bits per coefficient b must be >= 1")
    min_v = -(1 << (b - 1))
    max_v = (1 << (b - 1)) - 1
    packed = []

    for a in s1:
        v = _centered_mod_q(a, q)
        if v < min_v:
            v = min_v
        elif v > max_v:
            v = max_v
        if v < 0:
            v = (1 << b) + v
        packed.append(v)

    return _pack_bits(packed, b, slen)


def decompress(data, slen, d, q):
    if d <= 0:
        return []
    b = slen // d
    if b <= 0:
        return None
    need_bytes = (slen + 7) // 8
    if len(data) < need_bytes:
        return None

    raw = _unpack_bits(data[:need_bytes], b, d)
    s1 = []
    sign_bit = 1 << (b - 1)
    full = 1 << b

    for v in raw:
        if v & sign_bit:
            v = v - full
        s1.append(v % q)

    return s1
