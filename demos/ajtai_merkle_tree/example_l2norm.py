"""
Ajtai-Merkle tree well-formedness, proved with LaBRADOR  (L2-norm version).

We commit to a secret witness w of 2**k ring elements, build a t-ary Ajtai-hash
Merkle tree over the base-b gadget decomposition of every node, and produce a
single LaBRADOR proof that the public root is well formed.

The gadget uses base b = 8, so node digits live in [0, b) and every committed
digit vector is certified with the L2-squared bound  #coeffs * (b-1)^2  -- a
genuine short-vector (L2-norm) constraint, which is exactly the binding
requirement for the Ajtai/SIS hash.  A larger base means fewer digits per element
(h = ceil(log_b q)), which keeps the witnesses norm-rich and lets the tree scale.

Security: estimate_sis_params.sage shows base-8 / t=6 reaches ~147-bit at n=2,
which already clears 128-bit, so we use n = 2.  A height-2 Ajtai hash needs arity
t >= 2n = 4 to compress; we take n = 2, t = 6, base 8, k = 10 (1024 leaves).

Run:   conda activate sage && python3 example_l2norm.py
"""

from ajtai_merkle import run_demo
from labrador import LAB_RING_40

if __name__ == "__main__":
    ok = run_demo(
        k=14,                # 2**10 = 1024 leaves
        t=8,                 # arity (>= 2n so the height-2 hash compresses)
        n=2,                 # Ajtai hash height -- estimator's 128-bit choice (~147-bit)
        base=8,              # base-8 gadget: digits in [0,8), genuine L2 bound
        lab_ring=LAB_RING_40,
        seed=b"\xb2" * 32,
        title="EXAMPLE -- NON-BINARY Ajtai-Merkle tree (base 8, L2-norm, 128-bit)",
    )
    print()
    print("RESULT:", "PROOF VERIFIES " if ok else "PROOF FAILED ")
    raise SystemExit(0 if ok else 1)
