# Ajtai-hash Merkle tree well-formedness, proved with LaBRADOR

This demo proves, with a **single succinct LaBRADOR proof**, that a public root
is the correct **Ajtai-hash Merkle root** of a secret witness `w` of `2**k` ring
elements. It is the LaBRADOR version of the construction sketched in
`python/treethings/tree.py` and inspired by
[`aayux/lazer/python/ajtaitree`](https://github.com/aayux/lazer/tree/main/python/ajtaitree)
(which instead uses the per-node `lin_prover` system and notes that proving each
node separately is *not sound* — here the whole tree is one aggregated relation,
which is sound).

The single runnable variant uses a base-`b` gadget with a genuine **L2-norm
check** on every committed digit vector:

| script              | gadget | per-witness check          | params (estimator-backed)            |
|---------------------|--------|----------------------------|--------------------------------------|
| `example_l2norm.py` | base 8 | L2-norm bound `#coeffs·(b-1)²` | `n=2, t=6, k=10` → **~147-bit** |

Run it — it prints `PROOF VERIFIES`.

---

## The construction

Let `R_q = Z_q[X]/(X^d + 1)` be the LaBRADOR ring (`d = 64`). Fix a decomposition
base `b` and let `h = ceil(log_b q)`, with gadget vector
`g = (1, b, b², …, b^{h-1})`. Then `<g, G_b^{-1}(x)> == x` for every `x ∈ R_q`,
where `G_b^{-1}(x)` are the `h` base-`b` digits of `x` (each in `[0,b)`).

**Ajtai hash** (one public matrix `A ∈ R_q^{n × (t·h)}`, height `n`, arity `t`):
hashing `t` ring elements `x = (x₀,…,x_{t-1})` means
```
        digits = G_b^{-1}(x) ∈ [0,b)^{t·h}        (short)
        hash(x) = A · digits ∈ R_q^n
```
i.e. `t` elements are compressed to `n` (needs `t > n`; in fact `t ≥ 2n` so the
tree always compresses).

**The tree** is this hash applied level by level to a flat element-vector:
```
   level 0 :  w                          (2**k elements, SECRET)
   level ℓ :  cut the digit-vector of the level into chunks of t elements,
              hash each chunk -> n elements, concatenate (n·#chunks elements)
   final   :  one chunk of ≤ t elements hashes to the PUBLIC root ∈ R_q^n
```
The number of levels is derived automatically (`compute_tree_shape`). The last
chunk of a level — and the final level — may be smaller than `t`; then a
narrower column-slice of `A` is used ("the final level may be smaller if the
dims do not align").

**What is committed / proved.** We commit to `G_b^{-1}(w)` (the level-0 digits)
and, for every chunk `w_i` of every level, to `v_i = G_b^{-1}(A · w_i)`. For each
output coordinate `j ∈ [n]` we add one **linear** constraint over `R_q`:
```
        <a_j, w_i>  ==  <g, v_{i,j}>          (intermediate levels)
        <a_j, w_i>  ==  root_j                (final level; RHS is the public root)
```
where `a_j` is row `j` of `A` and `v_{i,j}` is the `j`-th `h`-digit block of `v_i`.
The digits flow straight through: `v_i` *is* the input of the next level, so each
constraint touches at most two committed vectors (one input chunk, one output
block) — that is what keeps the aggregated statement small and sound.

Every committed digit vector additionally carries an **L2 bound** (`betasq`):
for a vector of `c` polynomials with digits in `[0,b)`, `betasq = c·d·(b-1)²`.
Combined with the exact recomposition relations, the L2 bound on each level's
preimage is exactly the short-vector condition that makes the Ajtai/SIS hash
binding (collision-resistant), i.e. that makes the Merkle tree well-formed.

This maps onto LaBRADOR's "simple statement" interface (`python/labrador.py`):
linear constraints `<φ, witness> = b` plus a per-witness L2² bound. The whole
tree is **one** `proof_statement`, so the proof is succinct and the soundness is
that of a single relation `A·s = t` (no per-node "no-signaling" gap).

> We deliberately do **not** use LaBRADOR's `betasq = 0` "binary" mode (meant to
> force digits into {0,1}): it is broken in the shipped `liblabrador`
> (`polyvec_isbinary`, `src/labrador/poly.c:38`, mis-returns, and the prover
> rejects even honest binary witnesses). The L2-norm bound above is the working,
> sound mechanism — and it is what the SIS binding argument needs anyway.

---

## Parameters and the Lattice Estimator

The free Ajtai parameter is `n` (the height of `A`). As the task suggests, start
from `n = 2`, then run the estimator for the exact 128-bit value:

```
conda activate sage
sage estimate_sis_params.sage
```

`estimate_sis_params.sage` models the per-node hash as Module-SIS
`A·z = 0 (mod q)`, `||z||₂ ≤ β = 2·(b-1)·sqrt(t·h·d)` (a collision is the
difference of two short preimages), treats it as plain SIS of dimension
`(n·d) × (t·h·d)`, and uses the bundled lattice estimator
(`third_party/estimator.py`, `betaf`) for the BKZ root-Hermite ↔ block-size
conversion plus a core-SVP cost (`0.292·β_BKZ` classical, `0.265·β_BKZ` quantum).

Representative output (`d = 64`, `q ≈ 2^40`, base 8, `t = 6`):

```
n=1 →  51-bit
n=2 → 147-bit    ← smallest n with ≥128-bit (classical); used by the demo
n=3 → 255-bit
```

Because `q` is large (≈2^40) and the preimages are short, `n = 2` already clears
128-bit, so `example_l2norm.py` uses `n = 2` (with `t = 6`; a height-2 hash needs
`t ≥ 2n = 4` to compress). A smaller modulus would push `n` up (the script also
prints a `q ≈ 2^20` reference column).

> Why `q ≈ 2^40`? It is the **smallest** LaBRADOR modulus (`LAB_RING_40`) for
> which the *outer* proof's own SIS commitments verify reliably in this build:
> `q ≈ 2^24` fails (`init_proof`: "Cannot make inner commitments secure",
> `kappa > 32`) and `q ≈ 2^32` produces proofs that fail verification
> (error 12, "outer commitments not secure"). So "as small as permittable" = 40.

---

## Files

| file | purpose |
|------|---------|
| `ajtai_merkle.py`        | core library: `compute_tree_shape`, gadget helpers, `AjtaiMerkleProof` (build witnesses + assemble/prove/verify the LaBRADOR statement), `run_demo` |
| `example_l2norm.py`      | the demo — base-8 Ajtai-Merkle tree, L2-norm check, 128-bit params |
| `test_soundness.py`      | negative test: honest tree verifies, tampered tree is rejected |
| `estimate_sis_params.sage` | lattice-estimator search for the 128-bit Ajtai height `n` |

## How to run

```
conda activate sage
cd demos/ajtai_merkle_tree

python3 example_l2norm.py        # base-8 tree -> PROOF VERIFIES
python3 test_soundness.py        # honest verifies, tampered rejected
sage    estimate_sis_params.sage # 128-bit parameter tables

# or, with the convenience Makefile:
make            # runs the example + the soundness test
make estimate   # runs the lattice-estimator script
```

To explore other parameters, call `run_demo` directly, e.g.
```python
from ajtai_merkle import run_demo
from labrador import LAB_RING_40
run_demo(k=8, t=6, n=2, base=16, lab_ring=LAB_RING_40)   # 256-leaf base-16 tree
```
(`AjtaiMerkleProof` asserts `t > n` compression and verifies every
recomposition/norm relation in `_sanity_check` before proving.)
