import importlib.util
import re
import sys
import types
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ENERGY = ROOT / "source" / "modules" / "tdhf_mrsf_energy.F90"
LIB = ROOT / "source" / "tdhf_mrsf_lib.F90"
GRADIENT = ROOT / "source" / "modules" / "tdhf_mrsf_gradient.F90"
SINGLE_POINT = ROOT / "pyoqp" / "oqp" / "library" / "single_point.py"
OQPDATA = ROOT / "pyoqp" / "oqp" / "molecule" / "oqpdata.py"
INPUT_CHECKER = ROOT / "pyoqp" / "oqp" / "utils" / "input_checker.py"

# UMRSF implements only the energy path in production, plus an in-progress
# gradient gated behind OQP_UMRSF_GRAD_DEV that is enabled for runtype=grad
# ONLY. Every other non-energy runtype (Hessian/NAC/optimization, plus the
# gradient-driven prop/data) stays blocked even with the dev flag, because the
# developmental gradient must not silently drive numerical Hessians, couplings,
# or optimizations. ("thermo" is rejected earlier as an unknown runtype.)
UMRSF_ALWAYS_BLOCKED_RUNTYPES = (
    "prop", "data", "hess", "nac", "nacme",
    "optimize", "meci", "mecp", "mep", "ts", "irc", "neb",
)
# Enabled by the dev flag (and only by it).
UMRSF_DEV_GATED_RUNTYPES = ("grad",)


def compact(text: str) -> str:
    return re.sub(r"\s+", "", text.lower())


def _load_input_checker():
    """Import input_checker.py in isolation with a minimal MPI stub."""
    sys.modules.setdefault("oqp", types.ModuleType("oqp"))
    sys.modules.setdefault("oqp.utils", types.ModuleType("oqp.utils"))
    mpi_utils = types.ModuleType("oqp.utils.mpi_utils")

    class MPIManager:
        size = 1
        use_mpi = False

    mpi_utils.MPIManager = MPIManager
    sys.modules["oqp.utils.mpi_utils"] = mpi_utils

    spec = importlib.util.spec_from_file_location(
        "input_checker_umrsf_under_test", INPUT_CHECKER
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _umrsf_config(runtype):
    return {
        "input": {
            "runtype": runtype,
            "method": "tdhf",
            "basis": "6-31g",
            "system": "\nO 0.0 0.0 0.0\nH 0.0 0.0 0.95\nH 0.9 0.0 -0.3",
        },
        "guess": {},
        "scf": {"type": "uhf", "multiplicity": 3},
        "tdhf": {"type": "umrsf", "nstate": 2},
        "properties": {"grad": [1]},
        "optimize": {"lib": "geometric", "istate": 1, "jstate": 2},
        "nac": {"states": [[1, 2]]},
        "neb": {"product": "", "nimage": 3},
    }


def _umrsf_guard_errors(report):
    return [
        diag
        for diag in report.errors
        if diag.path == "tdhf.type" and "supports runtype=energy" in diag.message.lower()
    ]


class UMRSFEnergyRegressionTests(unittest.TestCase):
    def test_umrsf_mixed_exchange_channels_follow_mrsf_permutation_pattern(self):
        source = compact(LIB.read_text())
        self.assertIn(
            "f3(:nf,9:10,i,k)=f3(:nf,9:10,i,k)-xval*d3(:nf,9:10,j,l)",
            source,
        )
        self.assertIn(
            "f3(:nf,9:10,k,i)=f3(:nf,9:10,k,i)-xval*d3(:nf,9:10,l,j)",
            source,
        )
        self.assertIn(
            "f3(:nf,9:10,i,l)=f3(:nf,9:10,i,l)-xval*d3(:nf,9:10,j,k)",
            source,
        )
        self.assertIn(
            "f3(:nf,9:10,l,i)=f3(:nf,9:10,l,i)-xval*d3(:nf,9:10,k,j)",
            source,
        )

    def test_umrsf_flag_is_scoped_to_umrsf_entry_point(self):
        source = compact(ENERGY.read_text())
        self.assertIn("subroutinetdhf_mrsf_energy_c", source)
        self.assertIn("inf%tddft%umrsf=.false.", source)
        self.assertIn("logical::previous_umrsf", source)
        self.assertIn("previous_umrsf=inf%tddft%umrsf", source)
        self.assertIn("inf%tddft%umrsf=previous_umrsf", source)

    def test_umrsf_jacobi_rotation_intent_and_diagonal_are_consistent(self):
        lib = LIB.read_text().lower()
        energy = compact(ENERGY.read_text())
        self.assertRegex(
            lib,
            r"real\(kind=dp\),\s*intent\(inout\),\s*dimension\(:,:\)\s*::\s*mo_a,\s*mo_b",
        )
        self.assertIn("mo_energy_work_a", energy)
        self.assertIn("mo_energy_work_a(i)=fa(i,i)", energy)
        self.assertIn("mo_energy_work_b(i)=fb(i,i)", energy)
        self.assertIn("callmrinivec(infos,mo_energy_work_a,mo_energy_work_b", energy)

    def test_umrsf_is_registered_as_a_tdhf_type(self):
        oqpdata = compact(OQPDATA.read_text())
        single = SINGLE_POINT.read_text().lower()

        self.assertIn("'umrsf'", oqpdata)
        self.assertIn("self._data.tddft.umrsf=td_type=='umrsf'", oqpdata)
        self.assertIn("umrsf-tddft gradients are not implemented", single)

    @staticmethod
    def _set_grad_dev(value):
        import os
        if value is None:
            os.environ.pop("OQP_UMRSF_GRAD_DEV", None)
        else:
            os.environ["OQP_UMRSF_GRAD_DEV"] = value

    def _report(self, runtype):
        checker = _load_input_checker()
        return checker.check_input_values(
            _umrsf_config(runtype), raise_error=False, emit=False
        )

    def test_umrsf_energy_runtype_is_not_blocked(self):
        for dev in (None, "1"):
            with self.subTest(dev=dev):
                self._set_grad_dev(dev)
                try:
                    self.assertEqual(
                        _umrsf_guard_errors(self._report("energy")),
                        [],
                        "UMRSF energy must never be rejected by the runtype guard.",
                    )
                finally:
                    self._set_grad_dev(None)

    def test_umrsf_grad_is_gated_by_the_dev_flag(self):
        # Without the dev flag: grad is blocked. With it: grad is allowed.
        self._set_grad_dev(None)
        try:
            self.assertEqual(
                len(_umrsf_guard_errors(self._report("grad"))), 1,
                "Without OQP_UMRSF_GRAD_DEV, runtype=grad must be blocked.",
            )
            self._set_grad_dev("1")
            self.assertEqual(
                _umrsf_guard_errors(self._report("grad")), [],
                "With OQP_UMRSF_GRAD_DEV, runtype=grad must be allowed.",
            )
        finally:
            self._set_grad_dev(None)

    def test_umrsf_nongrad_runtypes_blocked_even_with_dev_flag(self):
        # Hessian/NAC/optimization etc. must stay blocked EVEN with the dev
        # flag: the developmental gradient must not drive them.
        for dev in (None, "1"):
            self._set_grad_dev(dev)
            try:
                for runtype in UMRSF_ALWAYS_BLOCKED_RUNTYPES:
                    with self.subTest(runtype=runtype, dev=dev):
                        guard_errors = _umrsf_guard_errors(self._report(runtype))
                        self.assertEqual(
                            len(guard_errors), 1,
                            f"runtype={runtype} (dev={dev}) must raise exactly "
                            f"one UMRSF guard error, got {len(guard_errors)}.",
                        )
                        self.assertIn(runtype, guard_errors[0].value)
            finally:
                self._set_grad_dev(None)

    def test_umrsf_intra_gamma_sp_is_not_transposed(self):
        # Two-set intra Gamma-SP: differentiating G_CO = sum (mu nu|ka la)
        # D^CO_{mu ka} D^CO_{nu la} gives the (ij|kl) cofactor co12(i,k)*
        # co12(j,l) -- the SECOND mixed-set density factor must NOT be the
        # transpose co12(l,j). The RO-copied transposed form was a bug
        # (invisible only when the density is symmetric, i.e. the RO limit).
        src = compact(GRADIENT.read_text())
        # corrected, un-transposed leading terms must be present (UMRSF path)
        self.assertIn("co12(i1,k1)*co12(j1,l1)", src)
        self.assertIn("o21v(i1,k1)*o21v(j1,l1)", src)
        # the RO-copied transposed leading term must NOT drive the UMRSF
        # intra-CO density (db1).  It legitimately remains in the RO get_density,
        # so we only require that the corrected term is the one feeding df1 via
        # umrsf_db1scale (i.e. the corrected pattern exists).
        self.assertIn("df1=df1+sgnk*qfspcp1*umrsf_db1scale()*db1", src)

    def test_umrsf_energy_does_not_use_mrsf_transition_density_output_path(self):
        source = compact(ENERGY.read_text())
        self.assertIn("if(umrsf)then", source)
        self.assertIn("trden=0.0_dp", source)
        self.assertIn("else", source)
        self.assertIn("callget_mrsf_transition_density", source)

    def test_spin_pair_scaling_avoids_hfscale_division_by_zero(self):
        source = compact(ENERGY.read_text())
        self.assertIn("if(abs(infos%tddft%hfscale)>epsilon(1.0_dp))then", source)
        self.assertIn("spc_scale_coco", source)
        self.assertIn("spc_scale_ovov", source)
        self.assertIn("spc_scale_coov", source)


if __name__ == "__main__":
    unittest.main()
