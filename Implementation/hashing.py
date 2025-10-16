import rng
import hashlib

DOM_SEP = {
    "H_msg":   b"HMSG",
    "H_pk":    b"HPK",
    "H_chal":  b"HCH",
    "H_seed":  b"HSEED",
    "H_poly":  b"HPOLY",
}

def _enc_len(x):
    return x.to_bytes(4, "big")


def H_sha256(parts, domain=b""):
    h = hashlib.sha256()
    if domain:
        h.update(_enc_len(len(domain))); h.update(domain)
    for p in parts:
        h.update(_enc_len(len(p))); h.update(p)
    return h.digest()


def H_sha256_int(parts, mod, domain=b""):
    digest = H_sha256(parts, domain)
    return int.from_bytes(digest, "big") % mod


def XOF_shake128(parts, outlen, domain=b""):
    x = hashlib.shake_128()
    if domain:
        x.update(_enc_len(len(domain))); x.update(domain)
    for p in parts:
        x.update(_enc_len(len(p))); x.update(p)
    return x.digest(outlen)


def XOF_shake256(parts, outlen, domain=b""):
    x = hashlib.shake_256()
    if domain:
        x.update(_enc_len(len(domain))); x.update(domain)
    for p in parts:
        x.update(_enc_len(len(p))); x.update(p)
    return x.digest(outlen)


def H_msg_to_int_mod_q(msg, q):
    return H_sha256_int([msg], q, DOM_SEP["H_msg"])


def H_pk_bind(pk_bytes):
    return H_sha256([pk_bytes], DOM_SEP["H_pk"])


def H_challenge(mu, transcript, out_bytes):
    return XOF_shake128([mu, transcript], out_bytes, DOM_SEP["H_chal"])


def H_seed_expand(seed, out_bytes):
    return XOF_shake128([seed], out_bytes, DOM_SEP["H_seed"])


def H_to_small_poly(seed, deg, eta):
    nbits = deg * 2 * eta
    nbytes = (nbits + 7) // 8
    buf = H_seed_expand(seed, nbytes)
    coeffs = rng.sample_cbd(buf, eta)
    while len(coeffs) < deg:
        buf = XOF_shake128([seed, len(coeffs).to_bytes(4,"big")], 32, DOM_SEP["H_poly"])
        coeffs.extend(rng.sample_cbd(buf, eta))
    return coeffs[:deg]
