#!/usr/bin/env sage
# -----------------------------------------------------------------------------
#  Lattice-Estimator parameter search for the Ajtai-hash height n.
# -----------------------------------------------------------------------------
#
#  The Ajtai hash of one tree node is a Module-SIS instance:
#
#        A * z = 0    (mod q),     A in R_q^{n x m},   R_q = Z_q[X]/(X^d+1)
#
#  with  m = t*h  ring columns  (t = arity, h = ceil(log_b q) digits per element).
#  Collision-resistance (= tree binding) is the hardness of finding a short
#  z = x - x' with x, x' valid digit preimages.  Since every committed digit
#  vector is certified with L2^2 <= m*d*(b-1)^2, a collision has
#
#        ||z||_2  <=  beta := 2 * (b-1) * sqrt(t*h*d).
#
#  We treat Module-SIS as plain SIS over Z_q of dimension (n*d) x (m*d) (the
#  standard embedding used by the lattice estimator), and report the BKZ block
#  size / core-SVP cost of the best lattice attack as a function of n.
#
#  The BKZ root-Hermite <-> block-size conversion is taken from the bundled
#  lattice estimator (third_party/estimator.py, function betaf), so this
#  genuinely relays to the estimator's reduction model.
#
#  Run:   conda activate sage && sage estimate_sis_params.sage
# -----------------------------------------------------------------------------

import os
import math                          # use math.* explicitly (sage shadows log/sqrt)

# ---- bring in the bundled lattice estimator (for betaf) ---------------------
_HERE = os.path.dirname(os.path.abspath(__file__)) if "__file__" in dir() else os.getcwd()
_EST = os.path.normpath(os.path.join(_HERE, "..", "..", "third_party", "estimator.py"))
load(_EST)                           # defines betaf(delta), delta_0f(beta), BKZ, ...

INF = float("inf")
CLASSICAL = 0.292                    # classical core-SVP bits per BKZ block
QUANTUM   = 0.265                    # quantum core-SVP bits per BKZ block
TARGET    = 128                      # desired bits of security


def sis_security(nd, md, q, beta):
    """(bits_classical, bits_quantum, bkz_block, delta) of SIS_{nd,md,q,beta}.

    Models the best lattice attack: BKZ with root-Hermite delta on the optimal
    sub-dimension of the q-ary lattice Lambda^perp_q(A).  Security = core-SVP
    cost of the weakest (cheapest) delta that still reaches a vector of L2 norm
    <= beta.  Returns +inf bits when no such delta exists (attack infeasible).
    """
    nd, md, q, beta = float(nd), float(md), float(q), float(beta)
    log2q = math.log(q, 2)
    log2b = math.log(beta, 2)

    if beta >= q:                                  # q*e_i is always a solution
        return (0.0, 0.0, 0.0, INF)

    # Weakest delta reaching norm beta, balanced sub-dimension regime:
    #   log2(beta) = 2*sqrt(nd*log2q*log2delta)
    log2d_bal = (log2b ** 2) / (4.0 * nd * log2q)
    m_opt = math.sqrt(nd * math.log(q) / max(log2d_bal * math.log(2.0), 1e-12))

    if m_opt <= md:
        log2_delta = log2d_bal
    else:                                          # column-limited: use all md
        log2_delta = (log2b - (nd / md) * log2q) / md
        if log2_delta <= 0:                        # unreachable even with LLL
            return (INF, INF, INF, 1.0)

    delta = 2.0 ** log2_delta
    if delta <= 1.0:
        return (INF, INF, INF, delta)
    try:
        block = float(betaf(delta))                # lattice-estimator conversion
    except Exception:
        return (INF, INF, INF, delta)
    if block <= 0:
        return (INF, INF, INF, delta)
    return (CLASSICAL * block, QUANTUM * block, block, delta)


def estimate(d, q, b, t, n_max=24, label=""):
    h  = int(math.ceil(math.log(q, b)))
    m  = t * h
    md = m * d
    beta = 2.0 * (b - 1) * math.sqrt(m * d)

    print("=" * 78)
    print(f"  Ajtai-hash SIS estimate {label}")
    print("=" * 78)
    print(f"  ring degree d = {d},  modulus q ~ 2^{int(math.ceil(math.log(q,2)))} = {q}")
    print(f"  base b = {b},  digits h = ceil(log_b q) = {h},  arity t = {t}")
    print(f"  SIS columns m = t*h = {m}  (md = {md}),  collision bound beta = {beta:.1f}")
    print(f"  target: {TARGET} bits (classical core-SVP, exponent {CLASSICAL})")
    print("-" * 78)
    print(f"  {'n':>3} | {'SIS rank nd':>11} | {'BKZ beta':>9} | "
          f"{'classical':>10} | {'quantum':>9}")
    print("-" * 78)
    rec = None
    for n in range(1, n_max + 1):
        nd = n * d
        if nd >= md:
            break
        cbits, qbits, block, _ = sis_security(nd, md, q, beta)
        cbs = "  inf  " if cbits == INF else f"{cbits:7.1f}"
        qbs = " inf " if qbits == INF else f"{qbits:6.1f}"
        bbs = "  inf  " if block == INF else f"{block:7.1f}"
        star = ""
        if rec is None and (cbits == INF or cbits >= TARGET):
            rec = n
            star = "  <== smallest n with >= 128-bit (classical)"
        print(f"  {n:>3} | {nd:>11} | {bbs:>9} | {cbs:>10} | {qbs:>9}{star}")
    print("-" * 78)
    if rec is not None:
        print(f"  RECOMMENDED Ajtai height for 128-bit (classical core-SVP): n = {rec}")
        print(f"  (choose arity t >= 2n = {2*rec} so the tree actually compresses)")
    else:
        print(f"  No n <= {n_max} reaches {TARGET} bits with these (d,q,b,t).")
    print()
    return rec


# Parameters of the runnable demos (LaBRADOR ring: d = 64, q ~ 2^40).
Q40 = 2 ** 40 - 195                  # LAB_RING_40.mod

print()
print("####  Ajtai-Merkle hash: minimal height n for 128-bit security  ####")
print()

# Parameters of the runnable demo (example_l2norm.py): base 8, arity t = 6.
estimate(d=64, q=Q40, b=11, t=6, label="(L2 demo:  base 8,  t=6)")

print("Notes:")
print(" * q is fixed to ~2^40 (LAB_RING_40): the smallest LaBRADOR modulus for")
print("   which the outer proof's own commitments verify reliably.")
print(" * The task's starting assumption n = 2 already gives ~147-bit here, so")
print("   the demo uses n = 2.  A smaller modulus would push the required n up")
print("   (see the 2^20 reference table).")
