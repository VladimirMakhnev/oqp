"""L1 gate — finite-difference validation of the 2e spin-spin (SS) dipolar integral.

This drives the compiled ``ssc_int2_selftest`` bind(C) harness (source/modules/ssc_int2_selftest.F90).
For a set of (s,p) shell quartets of a small molecule it compares, per primitive quartet and per
Cartesian component, the ANALYTIC SS bare-Hessian integral H_kl (mod_ssc_int2::comp_ssc_int2_prim)
against a 3-level Richardson finite difference of the engine's OWN Coulomb ERI with electron 2
rigidly displaced (mod_ssc_int2::comp_eri2_prim_disp) -- displacing electron 2 == displacing the
1/r12 operator, so d^2 ERI / d dshift_k d dshift_l = H_kl. It also checks the traceless invariant
Tr(S) = 0 for S_kl = H_kl - (1/3) Tr(H) delta_kl.

See CLAUDE.md (§6 L1, §7) and benchmarks.md (L1). The analytic integral is independently validated
to machine precision against the Python prototype (tests/ssc_prototype_ssss.py): the one-center
(ss|ss) ratio H_xx/ERI(0) = -2*alpha/3 is reproduced exactly. Pathologically tight core primitives
(exponent > 100) are excluded from the FD comparison because the operator-displacement FD is
roundoff-limited there; the analytic path is identical and is covered by the prototype check.

Do NOT weaken the tolerances to make this pass (CLAUDE.md §8). Skipped unless the compiled OpenQP
runtime is importable (a built tree with OPENQP_ROOT set).
"""

import os
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SELFTEST_OUT = Path("/tmp/ssc_int2_selftest.out")

INPUT = """[input]
system=
   8   0.000000000   0.000000000  -0.041061554
   1  -0.533194329   0.533194329  -0.614469223
   1   0.533194329  -0.533194329  -0.614469223
charge=0
runtype=energy
basis=6-31g*
method=hf
[guess]
type=huckel
[scf]
multiplicity=1
type=rhf
"""


def _runtime_available() -> bool:
    try:
        os.environ.setdefault("OMP_NUM_THREADS", "1")
        import oqp  # noqa: F401
        from oqp.pyoqp import Runner  # noqa: F401
        return hasattr(oqp, "ssc_int2_selftest")
    except Exception:
        return False


@unittest.skipUnless(_runtime_available(), "compiled OpenQP runtime / ssc_int2_selftest not available")
class TestSSCIntegralsFiniteDifference(unittest.TestCase):
    def test_ss_integral_matches_engine_finite_difference(self):
        import oqp
        from oqp.pyoqp import Runner

        workdir = Path("/tmp/oqp_ssc_l1_test")
        workdir.mkdir(exist_ok=True)
        inp = workdir / "h2o.inp"
        inp.write_text(INPUT)

        if SELFTEST_OUT.exists():
            SELFTEST_OUT.unlink()

        runner = Runner(project="ssc_l1", input_file=str(inp), log=str(workdir / "h2o.log"))
        runner.run()
        oqp.ssc_int2_selftest(runner.mol)

        self.assertTrue(SELFTEST_OUT.exists(), "self-test produced no output file")
        result = SELFTEST_OUT.read_text()
        self.assertIn("SSC_INT2_SELFTEST PASS", result,
                      "analytic SS integral disagrees with engine finite difference:\n" + result)


if __name__ == "__main__":
    unittest.main()
