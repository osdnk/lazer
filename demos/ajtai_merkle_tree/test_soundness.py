"""
Soundness sanity check: an honest tree VERIFIES, a tampered one is REJECTED.

We build the same instance twice.  The second time we flip one coefficient of
the first committed leaf-digit *after* fixing the (correct) public root, so the
committed witness no longer hashes to that root.  LaBRADOR must reject it:
  * smpl_verify returns non-zero (a linear recomposition constraint fails), and
  * the proof either fails to generate or fails to verify.

Run:   conda activate sage && python3 test_soundness.py
"""

from ajtai_merkle import AjtaiMerkleProof
from labrador import LAB_RING_40


def attempt(tamper):
    amp = AjtaiMerkleProof(k=4, t=6, n=2, base=8, lab_ring=LAB_RING_40,
                           seed=b"\xc3" * 32)
    amp.build(tamper=tamper)
    sv = amp.PS.simple_verify(amp.PS.smplstmnt_ptr, amp.PS.witness_ptr)
    ok, _, _ = amp.prove_and_verify(run_smpl_verify=False)
    return sv, ok


if __name__ == "__main__":
    print("=" * 70)
    print("SOUNDNESS CHECK")
    print("=" * 70)

    sv, ok = attempt(tamper=False)
    print(f"honest tree   : smpl_verify={sv} (0=ok), proof verifies={ok}")
    assert sv == 0 and ok, "honest tree should verify!"

    sv, ok = attempt(tamper=True)
    print(f"tampered tree : smpl_verify={sv} (non-0=caught), proof verifies={ok}")
    assert sv != 0, "tamper should break a linear constraint (smpl_verify != 0)"
    assert not ok, "tampered tree must NOT verify!"

    print()
    print("RESULT: honest tree verified, tampered tree rejected -- proof is sound. ")
