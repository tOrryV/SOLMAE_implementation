# SOLMAE Implementation

Post-quantum lattice-based digital signature scheme (SOLMAE) — full Python implementation  
This project provides a research-oriented implementation of the SOLMAE digital signature scheme based on NTRU lattices.

---

## Project Overview

This repository contains a full software implementation of the post-quantum digital signature scheme **SOLMAE**.  
The implementation is based on NTRU lattice constructions and Gaussian sampling techniques.

The purpose of this project is:
- educational exploration of lattice-based cryptography
- demonstration of modern post-quantum signature construction
- experimental verification of SOLMAE correctness and structure
- research-oriented implementation and testing framework

The implementation covers the complete pipeline from polynomial arithmetic to cryptographic signing and verification.

---

## Implemented Functionality

This repository implements:

- Polynomial arithmetic in ring:  (R_q = Z_q[X] / (X^n + 1))
- Modular arithmetic operations
- NTT (Number Theoretic Transform) and FFT acceleration
- Gaussian, Peikert and CBD samplers
- Secure random generation (HMAC-DRBG)
- Cryptographic hashing (SHA-256, SHAKE)
- Uniform polynomial generation (UnifCrown)
- NTRU-based public key generation
- Signature compression and decompression
- SOLMAE KeyGen / Sign / Verify
- Deterministic seeding and reproducibility
- Unit and integration tests
- End-to-end demonstration program

---

## Repository Structure

/
README.md               ← This file
params.py               ← Cryptographic parameters
hashing.py              ← Hashing and XOF
rng.py                  ← Random generators (HMAC-DRBG, uniform, CBD)
poly.py                 ← Polynomial ring arithmetic
modular.py              ← Modular arithmetic
ntt.py                  ← NTT / INTT and convolution
cfft.py                 ← FFT implementation
pairgen.py              ← Pair generation
unifcrown.py            ← Uniform polynomial sampling
ntrusolve.py            ← NTRU lattice solving utilities
samplers.py             ← Gaussian and Peikert samplers
sample_precomp.py       ← Precomputation tables
comp_decom.py           ← Compression / decompression
algoritm_solmae.py      ← SOLMAE core: KeyGen, Sign, Verify
demo_solmae.py          ← Demonstration script
tests/                  ← Test suite
MAC_lab.pdf             ← Report

---

## How to Run

1. Clone the repository

git clone https://github.com/tOrryV/SOLMAE_implementation.git
cd SOLMAE_implementation

2. (Optional) Create virtual environment

python -m venv venv
source venv/bin/activate   # Linux/macOS
venv\Scripts\activate      # Windows

3. Run the demo

python demo_solmae.py

4. Run the test suite

python main.py

---

## Testing Summary

More than **80 automatic tests** validate correctness and stability.

---

## Example Output

✔ Signature is valid for original message  
✖ Signature verification fails after message modification  

---

## References

- NTRU Cryptosystem
- Lattice-based cryptography
- GPV signatures
- SOLMAE paper

---

