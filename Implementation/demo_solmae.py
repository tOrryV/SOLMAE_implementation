from algoritm_solmae import keygen_solmae, sign_solmae, verify_solmae
from params import n, k, q, eta, d


def short_list(lst, limit=8):
    if len(lst) <= limit:
        return str(lst)
    return f"{lst[:limit]} ... (total {len(lst)} elements)"


def solmae_demo():
    print("=" * 70)
    print("                  SOLMAE DIGITAL SIGNATURE DEMO")
    print("=" * 70)

    msg = b"Hello, Mr. Anderson. Do u remember me?"
    print("\n[1] Message")
    try:
        print(f"  ASCII: {msg.decode('utf-8')}")
    except UnicodeDecodeError:
        print("  ASCII: <not printable>")
    print(f"  Bytes (hex): {msg.hex()}")

    print("\n[2] Scheme Parameters")
    print(f"  n   = {n}")
    print(f"  k   = {k}")
    print(f"  q   = {q}")
    print(f"  eta = {eta}")
    print(f"  d   = {d}")

    pk, sk = keygen_solmae()
    rho, t1 = pk
    print("\n[3] Key Generation")
    print("  Public Key:")
    print(f"    rho: length = {len(rho)} bytes")
    print(f"    rho (first 16 bytes hex): {rho[:16].hex()}")
    print(f"    t1: list of length k = {len(t1)}")
    print(f"    t1[0]: length = {len(t1[0])} bytes, first 16 bytes: {t1[0][:16].hex()}")

    print("  Secret Key:")
    s, e, t0, tr, pk2 = sk
    print(f"    s: {len(s)} polynomials, each of size {len(s[0].a)}")
    print(f"    e: {len(e)} polynomials, each of size {len(e[0].a)}")
    print(f"    t0: {len(t0)} low-bit vectors (first 8 of t0[0]):")
    print(f"      {short_list(t0[0], 8)}")
    print(f"    tr: length = {len(tr)} bytes, first 16 bytes: {tr[:16].hex()}")

    print("\n[4] Signature Generation")
    sig = sign_solmae(sk, msg)
    z, c, w1 = sig

    print("  Signature σ = (z, c, w1)")
    print(f"    c: length = {len(c)} bytes")
    print(f"    c (hex): {c.hex()}")

    print(f"    z: list of k = {len(z)} polynomials, each with {len(z[0])} coefficients")
    print(f"      z[0] (first 8 coefficients): {short_list(z[0], 8)}")

    print(f"    w1: list of k = {len(w1)} byte arrays of length {len(w1[0])}")
    print(f"      w1[0] (first 16 bytes hex): {w1[0][:16].hex()}")

    print("\n[5] Signature Verification")
    ok = verify_solmae(pk, msg, sig)
    print(f"  verify_solmae(pk, msg, sig) ⇒ {ok}")
    print(f"  Result: {'VALID SIGNATURE' if ok else 'INVALID SIGNATURE'}")

    print("\n[6] Tampering Demonstration")
    tampered_msg = msg + b" (That's Mr. Smith)"
    ok_tampered = verify_solmae(pk, tampered_msg, sig)
    print(f"  Tampered message: {tampered_msg.decode('utf-8')}")
    print(f"  verify_solmae(pk, tampered_msg, sig) ⇒ {ok_tampered}")
    print(f"  Expected: {'VALID' if ok_tampered else 'INVALID'} (signature must fail)")

    print("\n" + "=" * 70)
    print("                       END OF SOLMAE DEMO")
    print("=" * 70)


if __name__ == "__main__":
    solmae_demo()
