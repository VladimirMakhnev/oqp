"""P2.1 — AO contraction of the first-order SS contribution to the ZFS D-tensor.

Drives the compiled ``ssc_dtensor_selftest`` bind(C) harness (source/modules/ssc_zfs.F90), which
contracts the validated native SS dipolar integral (mod_ssc_int2::comp_ssc_int2_prim) with the
ROHF spin density P^(alpha-beta) = DM_A - DM_B in the Coulomb-like (P_munu P_kaptau) and
exchange-like (P_mukap P_nutau) patterns, producing the six symmetric D-tensor components.

This test only checks the structural invariant available WITHOUT the L2 prefactor pin: Tr(D) = 0
(the dipolar kernel is traceless). The absolute scale/sign C and the a.u.->cm^-1 unit factor are
pinned numerically against the O2 benchmark at L2 (not here); the harness uses C = 1.

Skipped unless the compiled OpenQP runtime / ssc_dtensor_selftest is importable.
"""

import os
import unittest
from pathlib import Path

SELFTEST_OUT = Path("/tmp/ssc_dtensor_selftest.out")

# O2 3Sigma_g- at r(O-O) = 1.207 A, ROHF triplet, polarized (d-containing) basis.
INPUT = """[input]
system=
   8   0.000000000   0.000000000   0.000000000
   8   0.000000000   0.000000000   1.207000000
charge=0
runtype=energy
basis=6-31g*
method=hf
[guess]
type=huckel
[scf]
multiplicity=3
type=rohf
"""


def _runtime_available() -> bool:
    try:
        os.environ.setdefault("OMP_NUM_THREADS", "1")
        import oqp  # noqa: F401
        from oqp.pyoqp import Runner  # noqa: F401
        return hasattr(oqp, "ssc_dtensor_selftest")
    except Exception:
        return False


@unittest.skipUnless(_runtime_available(), "compiled OpenQP runtime / ssc_dtensor_selftest not available")
class TestSSCDTensorContraction(unittest.TestCase):
    def test_contraction_runs_and_is_traceless(self):
        import oqp
        from oqp.pyoqp import Runner

        workdir = Path("/tmp/oqp_ssc_dtensor_test")
        workdir.mkdir(exist_ok=True)
        inp = workdir / "o2.inp"
        inp.write_text(INPUT)

        if SELFTEST_OUT.exists():
            SELFTEST_OUT.unlink()

        runner = Runner(project="o2_ssc", input_file=str(inp), log=str(workdir / "o2.log"))
        runner.run()
        oqp.ssc_dtensor_selftest(runner.mol)

        self.assertTrue(SELFTEST_OUT.exists(), "self-test produced no output file")
        result = SELFTEST_OUT.read_text()
        self.assertIn("SSC_DTENSOR_SELFTEST PASS", result,
                      "SS D-tensor contraction failed the traceless invariant:\n" + result)


if __name__ == "__main__":
    unittest.main()
