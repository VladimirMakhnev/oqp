"""L1 gate — finite-difference validation of the 2e spin-spin (SS) dipolar integral.

This is the FIRST regression test of the SSC/ZFS project (branch `ssc-zfs`). See `CLAUDE.md`
(§6 L1, §7) and `benchmarks.md` (L1). It is a SKELETON: the SS integral routine does not exist
yet, so the test is skipped until it does. Do NOT weaken the tolerances to make it pass — fix the
integral code (`CLAUDE.md §8`).

What L1 must verify, once the integral exists:

  The rank-2 dipolar kernel  T_kl(r12) = (3 r12,k r12,l - delta_kl r12^2) / r12^5
  is the Hessian of 1/r12 w.r.t. the interelectronic vector. Therefore each two-electron SS
  integral component <mu nu | T_kl | kappa tau> can be reproduced by FINITE DIFFERENCES of the
  existing ERI engine (the second mixed/again-diagonal derivative of the Coulomb operator),
  with NO external reference number.

  Checks:
    1. each of the 6 components (xx, yy, zz, xy, xz, yz) agrees with the ERI finite-difference
       value to 6-8 significant figures (h ~ 1e-4 a.u.; Richardson-extrapolate if needed);
    2. the trace is zero:  T_xx + T_yy + T_zz = 0  to <= 1e-10 for every shell quartet tested.

Test molecules/quartets: start with an (s,s,s,s) quartet (e.g. He2 or H2/STO-3G), then add p and
d shells once general angular momentum is implemented (P1.2).
"""

import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

# Tolerances (L1 gate; keep in sync with benchmarks.md L1 — do not loosen to pass).
SIGFIG_RTOL = 1e-6      # ~6-8 significant figures vs ERI finite difference
TRACE_ATOL = 1e-10     # Sum_k T_kk must vanish analytically


def _ssc_integrals_available() -> bool:
    """Return True once the native SS dipolar integral entry point is exposed.

    Update this probe when the integral routine lands (e.g. a `oqp.ssc_dipolar_int2` C-binding,
    or a Fortran self-test module wired like `hess1_selftest`). Until then L1 is skipped.
    """
    try:
        import oqp  # noqa: F401
    except Exception:
        return False
    # TODO(P1.x): replace with the real symbol once implemented, e.g.:
    #   return hasattr(oqp, "ssc_dipolar_int2")
    return False


@unittest.skipUnless(
    _ssc_integrals_available(),
    "SS dipolar 2e integral not implemented yet (L1 skeleton — see PROGRESS.md P1.1-P1.3).",
)
class TestSSCIntegralsFiniteDifference(unittest.TestCase):
    def test_components_match_eri_finite_difference(self):
        # P1.3: for each component kl in {xx,yy,zz,xy,xz,yz}, compare the analytic SS integral
        # to a central finite difference built from the existing ERI engine; assert relative
        # agreement <= SIGFIG_RTOL.
        self.fail("Not implemented: wire SS integral + ERI finite-difference reference (P1.3).")

    def test_trace_is_traceless(self):
        # P1.3: assert |T_xx + T_yy + T_zz| <= TRACE_ATOL for every shell quartet.
        self.fail("Not implemented: assert traceless kernel to TRACE_ATOL (P1.3).")


if __name__ == "__main__":
    unittest.main()
