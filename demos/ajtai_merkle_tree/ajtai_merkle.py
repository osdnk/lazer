"""
Proof of well-formedness of an Ajtai-hash Merkle tree, using LaBRADOR.
=====================================================================

This module proves, in zero knowledge with a *single* succinct LaBRADOR proof,
that a public root is the correct Ajtai-hash Merkle root of a (secret) witness
vector ``w`` of ``2**k`` ring elements.

The tree
--------
Let ``R_q = Z_q[X]/(X^d + 1)`` be the LaBRADOR ring (degree ``d``, modulus ``q``).
Fix a decomposition base ``b`` and let ``h = ceil(log_b q)`` so that every element
of ``R_q`` is recovered from its ``h`` base-``b`` digits by the gadget vector
``g = (1, b, b^2, ..., b^{h-1})``  (``<g, G_b^{-1}(x)> == x``).

The Ajtai hash uses one public matrix ``A in R_q^{n x (t*h)}`` (height ``n``,
arity ``t``).  One hash step eats ``t`` ring elements ``x = (x_0..x_{t-1})``:

    digits  = G_b^{-1}(x)              in [0,b)^{t*h}          (small)
    hash(x) = A * digits              in R_q^n

i.e. ``t`` elements are compressed to ``n`` elements (needs ``t > n``).

The whole tree is just this hash applied level by level to a flat vector:

    level 0 :  w                       (2**k elements, secret)
    level l :  chop the digit vector of the current level into chunks of t
               elements, hash every chunk -> n elements, concatenate.
    last    :  one chunk of <= t elements hashes to the public root in R_q^n.

The number of levels is derived automatically; the final level may use a
narrower slice of ``A`` if the count is not a clean power of ``t/n``.

What is committed / proved
--------------------------
* We commit to ``G_b^{-1}(w)``  (the level-0 digits) and
* for every chunk ``w_i`` of every level we commit to ``v_i = G_b^{-1}(A * w_i)``
  and prove the recomposition relation, for each output coordinate ``j in [n]``:

        <a_j, w_i>  ==  <g, v_{i,j}>            (intermediate levels)
        <a_j, w_i>  ==  root_j                  (final level, public RHS)

  where ``a_j`` is row ``j`` of ``A`` and ``v_{i,j}`` is the ``j``-th block of
  ``h`` digits of ``v_i``.  These are *linear* constraints over ``R_q`` -- exactly
  what the LaBRADOR "simple statement" interface supports.

* Every committed digit vector additionally carries an L2 bound (``betasq``):
  with base ``b`` and a digit vector of ``c`` polynomials, ``betasq = c*d*(b-1)**2``
  (the worst-case L2-squared of digits in [0,b)).  Combined with the exact
  recomposition relations, the L2 bound on every level's preimage is what makes
  the underlying SIS/Ajtai hash binding.

We use LaBRADOR's L2-norm witness mode (``betasq > 0``).  Its alternative
``betasq = 0`` "binary" mode -- meant to force coefficients into {0,1} -- is
broken in the shipped ``liblabrador`` (``polyvec_isbinary`` mis-returns and the
prover rejects even honest binary witnesses), so it is not used here.
"""

import os
import sys
import math
import time

# Make the lazer python module (and its compiled _lazer_cffi) importable.
_HERE = os.path.dirname(os.path.abspath(__file__))
_PYDIR = os.path.normpath(os.path.join(_HERE, "..", "..", "python"))
if _PYDIR not in sys.path:
    sys.path.insert(0, _PYDIR)

from lazer import *          # poly_t, polyvec_t, polymat_t, polyring_t, ffi, lib, ...
from labrador import (        # the LaBRADOR "simple statement" prover/verifier
    proof_statement, pack_verify,
    LAB_RING_24, LAB_RING_32, LAB_RING_40, LAB_RING_48,
)


# --------------------------------------------------------------------------- #
#  Gadget helpers (base-b decomposition / recomposition over R_q)             #
# --------------------------------------------------------------------------- #

def gadget_len(q, base):
    """h = ceil(log_base q): number of base-`base` digits to cover [0, q)."""
    return math.ceil(math.log(q, base))


def gadget_vec(ring, base, h):
    """g = (1, base, base^2, ..., base^{h-1}) as a length-h polyvec over `ring`."""
    return polyvec_t(ring, h, [poly_t(ring, {0: base ** i}) for i in range(h)])


def decompose(pol, base, h):
    """G_b^{-1}: return the h base-`base` digits of `pol` (coeff-wise, in [0,base)).

    Satisfies <gadget_vec(ring,base,h), decompose(pol,base,h)> == pol  (mod q).
    """
    p = pol.copy()
    p.redp()                       # bring every coefficient into [0, q)
    digits = polyvec_t(pol.ring, h)
    work = p.to_list()             # d integers in [0, q)
    for i in range(h):
        digits[i] = poly_t(pol.ring, [c % base for c in work])
        work = [c // base for c in work]
    assert all(c == 0 for c in work), "gadget_len too small for this modulus"
    return digits


# --------------------------------------------------------------------------- #
#  Tree shape (automatic number of levels)                                    #
# --------------------------------------------------------------------------- #

def compute_tree_shape(num_leaves, t, n):
    """Return the per-level chunking of a t-ary, height-n Ajtai tree.

    Each entry is a dict with:
        L      : number of input ring elements at this level
        sizes  : list of chunk sizes (each <= t, last one may be smaller)
        C      : number of chunks ( == number of hashes at this level )
        final  : True for the last level (single chunk -> public root)
    """
    assert t > n, "need arity t > height n for the tree to compress"
    assert num_leaves > n, "need 2**k > n leaves to compress at all"
    levels = []
    L = num_leaves
    while L > n:
        C = math.ceil(L / t)
        sizes = [t] * (C - 1) + [L - t * (C - 1)]
        final = (C == 1)
        levels.append({"L": L, "sizes": sizes, "C": C, "final": final})
        if final:
            break
        L_next = C * n
        assert L_next < L, "tree is not compressing (check t, n)"
        L = L_next
    return levels


# --------------------------------------------------------------------------- #
#  The proof builder                                                          #
# --------------------------------------------------------------------------- #

class LabradorProof:
    """A produced proof: LaBRADOR's witness commitment + the succinct argument.

    `witness_commitment` and `argument` are the two C structs returned by
    `composite_prove_simple`; `verify()` consumes both plus the public statement.
    They are owned by the originating AjtaiMerkleProof (kept alive via `_amp`);
    this build of the wrapper does not serialise them to bytes, so keep the object
    around rather than expecting to pickle it.
    """

    def __init__(self, amp, witness_commitment, argument, prove_secs):
        self._amp = amp                       # keep owner (and its C memory) alive
        self.witness_commitment = witness_commitment
        self.argument = argument
        self.prove_secs = prove_secs
        try:
            self.size_kb = float(argument.size)   # composite.size field (KB)
        except Exception:
            self.size_kb = float("nan")

    def __repr__(self):
        return f"<LabradorProof argument≈{self.size_kb:.1f} KB + witness commitment>"


class AjtaiMerkleProof:
    """Build a LaBRADOR statement+witness proving Ajtai-Merkle well-formedness.

    Parameters
    ----------
    k     : log2 of the number of leaves (witness length 2**k).
    t     : tree arity  (elements hashed per node).
    n     : Ajtai hash height (rows of A) -- the digest size in ring elements.
    base  : gadget decomposition base b (digits live in [0,b)).
    lab_ring : one of LAB_RING_{24,32,40,48}; 40 is the smallest that verifies.
    seed  : bytes seed for the public matrix A and the (secret) leaves.
    """

    def __init__(self, k=6, t=6, n=2, base=8,
                 lab_ring=LAB_RING_40, seed=b"\x00" * 32):
        self.k, self.t, self.n, self.base = k, t, n, base
        self.seed = seed

        self.ring = polyring_t(64, lab_ring.mod)
        self.q = self.ring.mod
        self.d = self.ring.deg
        self.primesize = str(math.ceil(math.log2(self.q)))
        self.h = gadget_len(self.q, base)

        self.num_leaves = 2 ** k
        self.shape = compute_tree_shape(self.num_leaves, t, n)
        # witness index where each level's chunk-witnesses start
        self.base_idx = []
        acc = 0
        for lv in self.shape:
            self.base_idx.append(acc)
            acc += lv["C"]
        self.num_witnesses = acc
        self.num_constraints = sum(lv["C"] for lv in self.shape) * n

        self.g = gadget_vec(self.ring, base, self.h)
        self.neg_g = polyvec_t(self.ring, self.h, [-self.g[i] for i in range(self.h)])
        self.neg_g.redc()

        # filled by build()
        self.A = None
        self.WV = None          # list of witness polyvecs (one per chunk)
        self.WV_npols = None    # declared #polys per witness
        self.root = None
        self.PS = None

    # -- public matrix -------------------------------------------------------
    def _sample_A(self):
        """Public Ajtai matrix A in R_q^{n x (t*h)}, uniform mod q (centralised)."""
        cols = self.t * self.h
        A = polymat_t(self.ring, self.n, cols)
        A.urandom(self.q, self.seed, 0)
        A.redc()
        return A

    def _betasq(self, npols):
        """L2-squared bound for a digit witness of `npols` polynomials in [0,b)."""
        return npols * self.d * (self.base - 1) ** 2   # tight short-vector bound

    # -- build witnesses by actually running the tree ------------------------
    def _build_witnesses(self):
        """Compute every chunk-witness value (and the public root) by hashing."""
        t, n, h, base = self.t, self.n, self.h, self.base
        WV = [None] * self.num_witnesses
        npols = [0] * self.num_witnesses

        # level 0 : decompose the secret leaves w  ->  G_b^{-1}(w)
        w = polyvec_t.urandom_bnd_static(self.ring, self.num_leaves,
                                         0, self.q - 1, self.seed, 7)
        w.redc()
        flat0 = []                                     # flat list of level-0 digit polys
        for e in range(self.num_leaves):
            flat0.extend(decompose(w[e], base, h).to_pol_list())
        self._pack_into_witnesses(flat0, level=0, WV=WV, npols=npols)

        # subsequent levels : hash chunk by chunk
        for lvl, info in enumerate(self.shape):
            if info["final"]:
                # one chunk -> public root in R_q^n  (no output digits)
                chunk = WV[self.base_idx[lvl]]
                s = info["sizes"][0]
                A_s = self._A_slice(s)
                self.root = A_s * chunk                # polyvec dim n, PUBLIC
                self.root.redc()
                break

            out_flat = []                              # flat list of output digit polys
            off = 0
            for i, s in enumerate(info["sizes"]):
                chunk = WV[self.base_idx[lvl] + i]
                A_s = self._A_slice(s)
                out = A_s * chunk                      # polyvec dim n  (= A * w_i)
                for j in range(n):
                    out_flat.extend(decompose(out[j], base, h).to_pol_list())
            self._pack_into_witnesses(out_flat, level=lvl + 1, WV=WV, npols=npols)

        self.WV, self.WV_npols = WV, npols

    def _A_slice(self, s):
        """First s*h columns of A (for a possibly-smaller last/final chunk)."""
        if s == self.t:
            return self.A
        return self.A.get_col_list(list(range(s * self.h)))

    def _pack_into_witnesses(self, flat_digits, level, WV, npols):
        """Group a flat digit list into this level's chunk-witnesses (size <= t*h)."""
        info = self.shape[level]
        pos = 0
        for i, s in enumerate(info["sizes"]):
            cnt = s * self.h
            idx = self.base_idx[level] + i
            WV[idx] = polyvec_t(self.ring, cnt, flat_digits[pos:pos + cnt])
            npols[idx] = cnt
            pos += cnt
        assert pos == len(flat_digits), (pos, len(flat_digits), level)

    # -- assemble the LaBRADOR statement ------------------------------------
    def _output_coeff(self, out_npols, block):
        """Coeff vector for an output witness: -g placed at digit-block `block`."""
        lst = [poly_t(self.ring) for _ in range(out_npols)]
        for r in range(self.h):
            lst[block * self.h + r] = -self.g[r]
        v = polyvec_t(self.ring, out_npols, lst)
        v.redc()
        return v

    def _input_coeff(self, j, s):
        """Row j of A restricted to the first s*h columns (the chunk width)."""
        row = self.A.get_row(j).to_pol_list()[: s * self.h]
        v = polyvec_t(self.ring, s * self.h, row)
        v.redc()
        return v

    def build(self, tamper=False):
        """Sample A, build all witnesses, and assemble the proof_statement.

        If ``tamper`` is True, a single coefficient of the first committed leaf
        digit is flipped *after* the honest root is computed -- the committed
        tree then no longer hashes to the (still-correct) public root, so the
        statement becomes false and the proof must be rejected.
        """
        self.A = self._sample_A()
        self._build_witnesses()
        self._sanity_check()
        if tamper:
            bad = self.WV[0].to_pol_list()
            p0 = bad[0].to_list()
            p0[0] += 1                      # flip one digit -> breaks a relation
            bad[0] = poly_t(self.ring, p0)
            self.WV[0] = polyvec_t(self.ring, self.WV_npols[0], bad)
            self._tampered = True

        deg_list = [self.d] * self.num_witnesses
        npols_list = list(self.WV_npols)
        if getattr(self, "tight_norms", False):
            norm_list = [max(1, self.WV[i].l2sqr()) for i in range(self.num_witnesses)]
        else:
            norm_list = [self._betasq(c) for c in npols_list]

        PS = proof_statement(deg_list, npols_list, norm_list,
                             self.num_constraints, self.primesize)
        for idx in range(self.num_witnesses):
            PS.append_witness(self.WV[idx])

        zero = poly_t(self.ring)
        t, n, h = self.t, self.n, self.h
        for lvl, info in enumerate(self.shape):
            for i, s in enumerate(info["sizes"]):
                in_idx = self.base_idx[lvl] + i
                for j in range(n):
                    a_j = self._input_coeff(j, s)
                    if info["final"]:
                        PS.append_statement([a_j], [in_idx], self.root[j])
                    else:
                        oi = i * n + j
                        out_idx = self.base_idx[lvl + 1] + oi // t
                        block = oi % t
                        coeff_out = self._output_coeff(self.WV_npols[out_idx], block)
                        PS.append_statement([a_j, coeff_out], [in_idx, out_idx], zero)
        self.PS = PS
        return PS

    # -- sanity: the witnesses really satisfy the relations & norms ----------
    def _sanity_check(self):
        t, n, h = self.t, self.n, self.h
        for lvl, info in enumerate(self.shape):
            for i, s in enumerate(info["sizes"]):
                chunk = self.WV[self.base_idx[lvl] + i]
                A_s = self._A_slice(s)
                out = A_s * chunk
                for j in range(n):
                    if info["final"]:
                        diff = out[j] - self.root[j]
                    else:
                        oi = i * n + j
                        out_idx = self.base_idx[lvl + 1] + oi // t
                        block = oi % t
                        ow = self.WV[out_idx]
                        v_block = polyvec_t(self.ring, h,
                                            [ow[block * h + r] for r in range(h)])
                        diff = out[j] - self.g * v_block
                    diff.redp()
                    assert diff.linf() == 0, f"recomposition mismatch at level {lvl}"
        # norm bounds
        for idx in range(self.num_witnesses):
            bnd = self._betasq(self.WV_npols[idx])
            if bnd > 0:
                assert self.WV[idx].l2sqr() <= bnd, f"witness {idx} exceeds L2 bound"

    # -- commitment / prove / verify (three separable phases) ----------------
    def total_witness_polys(self):
        return sum(self.WV_npols)

    @property
    def commitment(self):
        """The public commitment to the 2**k leaves: the Ajtai-Merkle root in R_q^n.

        Produced by build() purely by hashing -- no proving involved.  Publishing
        this binds the prover to the whole tree; the proof below shows it is well
        formed (i.e. that a short witness hashing to this root is known).
        """
        assert self.root is not None, "call build() first"
        return self.root

    def _verify_fn(self):
        return {
            "24": lib.labrador24_composite_verify_simple,
            "32": lib.labrador32_composite_verify_simple,
            "40": lib.labrador40_composite_verify_simple,
            "48": lib.labrador48_composite_verify_simple,
        }[self.primesize]

    def prove(self, run_smpl_verify=False):
        """Run the LaBRADOR prover. Returns a LabradorProof (commitment + argument).

        The returned object's `witness_commitment` is LaBRADOR's binding
        commitment to the secret witness; `argument` is the succinct proof.  Both
        are needed to verify; neither reveals the witness.
        """
        assert self.PS is not None, "call build() first"
        if run_smpl_verify:
            self.PS.smpl_verify()
        t0 = time.perf_counter()
        err, comp, comm = self.PS.pack_prove()
        dt = time.perf_counter() - t0
        if err != 0:
            raise RuntimeError(f"pack_prove failed with code {err}")
        return LabradorProof(self, comm, comp, dt)

    def verify(self, proof):
        """Verify a LabradorProof against this statement (uses only public data).

        Reads the public statement (the matrix A, the gadget g and the root) plus
        the proof's commitment and argument -- never the witness.  Returns bool.
        """
        stmnt = self.PS.output_statement()
        return self._verify_fn()(proof.argument, proof.witness_commitment, stmnt) == 0

    def prove_and_verify(self, run_smpl_verify=True):
        """Convenience: prove() then verify(). Returns (ok, prove_secs, verify_secs)."""
        try:
            proof = self.prove(run_smpl_verify=run_smpl_verify)
        except RuntimeError as e:
            print(f"[ERR] {e}")
            return False, 0.0, 0.0
        t0 = time.perf_counter()
        ok = self.verify(proof)
        return ok, proof.prove_secs, time.perf_counter() - t0


# --------------------------------------------------------------------------- #
#  Convenience runner used by the example scripts                             #
# --------------------------------------------------------------------------- #

def run_demo(k=6, t=6, n=2, base=8, lab_ring=LAB_RING_40,
             seed=b"\x11" * 32, title=None):
    """Build, prove and verify one Ajtai-Merkle instance; print a report."""
    if title:
        print("=" * 70)
        print(title)
        print("=" * 70)
    amp = AjtaiMerkleProof(k=k, t=t, n=n, base=base, lab_ring=lab_ring, seed=seed)
    print(f"ring R_q : degree d = {amp.d}, modulus q ~ 2^{amp.primesize} "
          f"({amp.q})")
    print(f"params   : k = {k}  (#leaves = {amp.num_leaves}),  arity t = {t},  "
          f"Ajtai height n = {n}")
    print(f"gadget   : base b = {base},  digits h = ceil(log_b q) = {amp.h}")
    print(f"check    : L2-norm bound betasq = #coeffs*(b-1)^2  (digits in [0,{base}))")
    print(f"levels   : {len(amp.shape)}  (chunk counts per level: "
          f"{[lv['C'] for lv in amp.shape]})")
    for lvl, info in enumerate(amp.shape):
        tag = " (final -> public root)" if info["final"] else ""
        print(f"   level {lvl}: {info['L']:>6} elems -> {info['C']:>4} hashes"
              f"{tag}")

    t0 = time.perf_counter()
    amp.build()
    build_s = time.perf_counter() - t0
    print(f"witness  : {amp.num_witnesses} committed digit-vectors, "
          f"{amp.total_witness_polys()} ring elements; "
          f"{amp.num_constraints} linear constraints")

    # ---- phase 1: COMMIT -- just hashing the tree, no proving --------------
    root = amp.commitment
    print(f"[1] COMMIT : root = {amp.n} ring elements binding all "
          f"{amp.num_leaves} leaves  (built in {build_s:.3f}s, no proving)")

    # ---- phase 2: PROVE  -- witness commitment + succinct argument ---------
    proof = amp.prove(run_smpl_verify=True)
    print(f"[2] PROVE  : {proof}  in {proof.prove_secs:.3f}s")

    # ---- phase 3: VERIFY -- public statement + proof only (no witness) -----
    t0 = time.perf_counter()
    ok = amp.verify(proof)
    print(f"[3] VERIFY : {'PROOF VERIFIES' if ok else 'PROOF FAILED'}  "
          f"in {time.perf_counter() - t0:.3f}s")
    return ok


if __name__ == "__main__":
    run_demo(title="Ajtai-Merkle well-formedness (base-8 gadget, L2-norm check)")
