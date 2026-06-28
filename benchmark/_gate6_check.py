"""GATE 6 validation: external-output parsing + comparison harness."""
import os, sys
import numpy as np
from oqp.interop import parse_output, parse_pyscf_tddft, parse_oqp, compare_results, format_table

ok = True

# ---- (1) cclib parser vs committed Gaussian fixture (known reference values) ----
fixture = "benchmark/fixtures/ch2o_gaussian_td.log"
ref_exc_ev = [3.9892, 8.9000, 9.5000]          # values embedded in the fixture
parsed = parse_output(fixture, program="gaussian")
print("=== cclib (Gaussian fixture) ===")
print("  parser:", parsed.get("parser"), " scf(Ha):", parsed.get("scf_energy_ha"))
print("  excitation energies (eV):", np.round(parsed.get("excitation_energies_ev", []), 4).tolist())
print("  oscillator strengths    :", parsed.get("oscillator_strengths"))
got = parsed.get("excitation_energies_ev", [])
dmax = max(abs(a - b) for a, b in zip(got, ref_exc_ev)) if got else 1e9
print(f"  max|parsed - source| = {dmax:.2e} eV  (tol 1e-4)")
cclib_ok = dmax < 1e-4 and len(got) == 3
ok &= cclib_ok

# ---- (2) PySCF external code: real TDDFT run, native object round-trip ----
print("\n=== PySCF TDDFT (formaldehyde, real external code) ===")
from pyscf import gto, dft, tddft
mol = gto.M(atom="""C 0.0 0.0 -0.529700; O 0.0 0.0 0.677500;
                    H 0.0 0.934200 -1.124000; H 0.0 -0.934200 -1.124000""",
            basis="6-31g*", unit="Angstrom", verbose=0)
mf = dft.RKS(mol); mf.xc = "bhandhlyp"; mf.kernel()
td = tddft.TDDFT(mf); td.nstates = 3; td.kernel()
source_ev = (np.asarray(td.e) * 27.211386245988).tolist()
pp = parse_pyscf_tddft(td, mf)
rt = max(abs(a - b) for a, b in zip(pp["excitation_energies_ev"], source_ev))
print("  excitation energies (eV):", np.round(pp["excitation_energies_ev"], 4).tolist())
print(f"  native round-trip |Δ| = {rt:.2e} eV (tol 1e-4)")
ok &= rt < 1e-4

# ---- (3) comparison harness renders a table (OQP MRSF vs PySCF TDDFT) ----
print("\n=== comparison harness (OQP MRSF vs PySCF TDDFT, formaldehyde) ===")
from oqp.pyoqp import Runner
from oqp.analysis import MRSFExcitedStates
r = Runner(project="ch2o_mrsf", input_file="benchmark/inputs/ch2o_mrsf.inp",
           log="/bighome/vova/coding/oqp-misc/runs/gate6_ch2o.log", silent=1, usempi=False)
r.run()
st = MRSFExcitedStates(r.mol)
oqp_res = parse_oqp(r.mol, st)
rows, _ = compare_results(oqp_res, pp,
                          {"excitation_energies_ev": 1.0, "oscillator_strengths": 1.0},
                          ref_label="OQP-MRSF", other_label="PySCF-TDDFT")
print(format_table(rows, "OQP-MRSF", "PySCF-TDDFT"))
print("  (cross-method comparison is informational; methods differ)")
table_ok = len(rows) == 2

print("\nGATE 6:", "PASS" if (ok and table_ok) else "FAIL")
sys.exit(0 if (ok and table_ok) else 1)
