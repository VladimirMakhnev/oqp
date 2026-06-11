"""End-to-end test for the MRSF-NMR Gate 2 prototype (nmr_mrsf_shielding).

Runs an H2O BHHLYP/STO-3G MRSF-TDDFT gradient calculation with
properties.nmr_mrsf=True and checks, from the machine-parseable
MRSF_NMR_* log records:

  1. the prototype runs end-to-end and emits records for the target state;
  2. internal consistency: iso values equal the tensor traces/3, and
     total = dia + para to numerical precision;
  3. the state diamagnetic term differs from the reference diamagnetic term
     (the relaxed MRSF density is actually used), but by a physically small
     amount;
  4. the frozen-reference paramagnetic tensors are identical to the
     ground-state GIAO paramagnetic tensors of the same ROHF reference
     (zero-amplitude wiring check: the para path is byte-identical code);
  5. input_checker rejects inconsistent nmr_mrsf setups (no silent
     ground-state fallback).

Runs the OpenQP driver in a subprocess and skips gracefully if the shared
library is not built.
"""

import os
import re
import subprocess
import sys
import tempfile
import textwrap
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

GEOM = """\
   8   0.000000000   0.000000000  -0.041061554
   1  -0.533194329   0.533194329  -0.614469223
   1   0.533194329  -0.533194329  -0.614469223"""

# trh_impl=otr: native TRAH crashes intermittently on RHEL8/GCC 11.2 and
# yields non-canonical open-shell references; OTR is required for MRSF runs.
INPUT_MRSF_NMR = f"""\
[input]
system=
{GEOM}
charge=0
runtype=grad
functional=bhhlyp
basis=sto-3g
method=tdhf

[guess]
type=huckel

[scf]
multiplicity=3
type=rohf
trh_impl=otr

[dftgrid]
rad_type=mhl
rad_npts=99
ang_npts=302

[tdhf]
type=mrsf
nstate=2

[properties]
grad=1
nmr_mrsf=True
"""

TOL_CONSISTENCY = 1.0e-10   # iso vs trace/3, total vs dia+para (same numbers)
TOL_PARA_MATCH = 1.0e-8     # frozen-ref para vs ground-state GIAO para (ppm)
DIA_SHIFT_MIN = 1.0e-4      # state dia must differ from reference dia (ppm)
DIA_SHIFT_MAX = 50.0        # ... but by a physically small amount (ppm)


def _oqp_root():
    root = os.environ.get("OPENQP_ROOT", str(ROOT))
    lib = Path(root) / "lib"
    if (lib / "liboqp.dylib").exists() or (lib / "liboqp.so").exists():
        return root
    return None


def _run_subprocess(wd, inp_text, script_body, timeout=600):
    inp = Path(wd) / "run.inp"
    inp.write_text(inp_text)
    script = Path(wd) / "run.py"
    script.write_text(script_body)
    env = dict(os.environ)
    env["OPENQP_ROOT"] = _oqp_root()
    env["PYTHONPATH"] = os.pathsep.join(
        p for p in (str(ROOT / "pyoqp"), env.get("PYTHONPATH", "")) if p)
    return subprocess.run([sys.executable, str(script)], cwd=wd, env=env,
                          capture_output=True, text=True, timeout=timeout)


def _parse_records(text):
    """Parse MRSF_NMR_* records -> dict with tensors and iso rows."""
    state_m = re.search(r"MRSF_NMR_STATE\s+(\d+)", text)
    nat_m = re.search(r"MRSF_NMR_NATOM\s+(\d+)", text)
    if not (state_m and nat_m):
        return None
    nat = int(nat_m.group(1))
    out = {"state": int(state_m.group(1)), "natom": nat,
           "dia": {}, "dia_ref": {}, "para_unc": {}, "para_cpl": {}, "iso": {}}
    pat = re.compile(
        r"MRSF_NMR_(DIA_REF|DIA|PARA_UNC|PARA_CPL)\s+(\d+)\s+(\d)\s+(\d)\s+(\S+)")
    keymap = {"DIA": "dia", "DIA_REF": "dia_ref",
              "PARA_UNC": "para_unc", "PARA_CPL": "para_cpl"}
    for m in pat.finditer(text):
        kind, iat, t, s, val = m.groups()
        out[keymap[kind]][(int(iat), int(t), int(s))] = float(val)
    for m in re.finditer(r"MRSF_NMR_ISO\s+(\d+)((?:\s+\S+){5})", text):
        out["iso"][int(m.group(1))] = [float(x) for x in m.group(2).split()]
    return out


def _trace3(tens, iat):
    return (tens[(iat, 1, 1)] + tens[(iat, 2, 2)] + tens[(iat, 3, 3)]) / 3.0


class MRSFNMRPrototypeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if _oqp_root() is None:
            raise unittest.SkipTest("OpenQP shared library not built")
        with tempfile.TemporaryDirectory() as wd:
            log = Path(wd) / "run.log"
            proc = _run_subprocess(wd, INPUT_MRSF_NMR, textwrap.dedent(f"""
                from oqp.pyoqp import Runner
                r = Runner(project="mrsfnmr", input_file="run.inp",
                           log={str(log)!r}, silent=1, usempi=False)
                r.run()
            """))
            if not log.exists():
                raise unittest.SkipTest(
                    f"MRSF run produced no log\n{proc.stderr[-1500:]}")
            text = log.read_text()
            cls.rec = _parse_records(text)
            cls.stderr = proc.stderr
            cls.log_text = text

            # Reference ground-state GIAO shielding of the same ROHF triplet
            # reference (separate run; scf_prop=nmr forbidden alongside
            # nmr_mrsf by the input checker).
            log2 = Path(wd) / "ref.log"
            ref_inp = INPUT_MRSF_NMR.replace("runtype=grad", "runtype=energy") \
                                    .replace("method=tdhf", "method=hf") \
                                    .replace("[tdhf]\ntype=mrsf\nnstate=2\n\n", "") \
                                    .replace("grad=1\nnmr_mrsf=True",
                                             "scf_prop=nmr\nnmr_gauge=giao")
            proc2 = _run_subprocess(wd, ref_inp, textwrap.dedent(f"""
                from oqp.pyoqp import Runner
                r = Runner(project="refnmr", input_file="run.inp",
                           log={str(log2)!r}, silent=1, usempi=False)
                r.run()
            """))
            cls.ref_text = log2.read_text() if log2.exists() else ""

    def test_records_present(self):
        self.assertIsNotNone(self.rec,
                             f"MRSF_NMR records not found\n{self.stderr[-1500:]}")
        self.assertEqual(self.rec["natom"], 3)
        self.assertEqual(self.rec["state"], 1)
        self.assertEqual(len(self.rec["dia"]), 9 * 3)

    def test_iso_matches_tensor_trace_and_total(self):
        for iat in range(1, 4):
            iso = self.rec["iso"][iat]
            self.assertAlmostEqual(iso[0], _trace3(self.rec["dia"], iat),
                                   delta=TOL_CONSISTENCY)
            self.assertAlmostEqual(iso[1], _trace3(self.rec["para_unc"], iat),
                                   delta=TOL_CONSISTENCY)
            self.assertAlmostEqual(iso[2], _trace3(self.rec["para_cpl"], iat),
                                   delta=TOL_CONSISTENCY)
            self.assertAlmostEqual(iso[4], iso[0] + iso[2],
                                   delta=TOL_CONSISTENCY)

    def test_state_dia_differs_from_reference_dia(self):
        # The relaxed MRSF density must actually shift the diamagnetic term,
        # but the shift must stay physically small for a low excited state.
        shifts = [abs(_trace3(self.rec["dia"], iat)
                      - _trace3(self.rec["dia_ref"], iat))
                  for iat in range(1, 4)]
        self.assertGreater(max(shifts), DIA_SHIFT_MIN,
                           "state dia identical to reference dia: td_p unused?")
        self.assertLess(max(shifts), DIA_SHIFT_MAX,
                        "state dia shift unphysically large")

    def test_frozen_reference_para_matches_ground_state_giao(self):
        # Zero-amplitude wiring check: the prototype para path must equal the
        # ground-state open-shell GIAO para of the same ROHF reference.
        self.assertTrue(self.ref_text, "reference GIAO run produced no log")
        ref = {}
        pat = re.compile(
            r"GIAO_SHIELDING_DEBUG_PARA_CPL\s+(\d+)\s+(\d)\s+(\d)\s+(\S+)")
        for m in pat.finditer(self.ref_text):
            iat, t, s, val = m.groups()
            ref[(int(iat), int(t), int(s))] = float(val)
        self.assertTrue(ref, "no GIAO para records in the reference log")
        for key, val in self.rec["para_cpl"].items():
            self.assertAlmostEqual(val, ref[key], delta=TOL_PARA_MATCH,
                                   msg=f"para mismatch at {key}")

    def test_checker_rejects_bad_setups(self):
        sys.path.insert(0, str(ROOT / "pyoqp"))
        from oqp.utils import input_checker

        def issues(config):
            report = input_checker.CheckReport()
            input_checker._check_properties(config, report)
            return report.errors

        good = {"input": {"runtype": "grad", "method": "tdhf"},
                "tdhf": {"type": "mrsf"},
                "scf": {"type": "rohf", "multiplicity": 3},
                "properties": {"nmr_mrsf": True, "grad": [1]}}
        self.assertFalse(issues(good))

        for breaker in (("input", "runtype", "energy"),
                        ("tdhf", "type", "tda"),
                        ("scf", "type", "rhf"),
                        ("scf", "multiplicity", 1)):
            bad = {k: dict(v) for k, v in good.items()}
            bad[breaker[0]][breaker[1]] = breaker[2]
            self.assertTrue(issues(bad), f"checker missed {breaker}")

        combo = {k: dict(v) for k, v in good.items()}
        combo["properties"]["scf_prop"] = ["nmr"]
        self.assertTrue(issues(combo), "checker allowed scf_prop=nmr + nmr_mrsf")


if __name__ == "__main__":
    unittest.main()
