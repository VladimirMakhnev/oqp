"""GATE 5 validation: QCSchema AtomicResult + FCIDUMP -> PySCF FCI round-trip."""
import os, sys
import numpy as np
from oqp.pyoqp import Runner
from oqp.analysis import MRSFExcitedStates
from oqp.export import to_qcschema, validate_qcschema, dump_fcidump, verify_fcidump_fci

RUNDIR = "/bighome/vova/coding/oqp-misc/runs/gate5"
os.makedirs(RUNDIR, exist_ok=True)
ok = True

# ---------- QCSchema (on the formaldehyde MRSF run) ----------
inp = "benchmark/inputs/ch2o_mrsf.inp"
r = Runner(project="ch2o_mrsf", input_file=inp,
           log=os.path.join(RUNDIR, "ch2o_mrsf.log"), silent=1, usempi=False)
r.run()
st = MRSFExcitedStates(r.mol)
payload = to_qcschema(r.mol, states=st)
res = validate_qcschema(payload)
import json
with open(os.path.join(RUNDIR, "ch2o_qcschema.json"), "w") as f:
    json.dump(payload, f, indent=2)
e_round = abs(res.properties.scf_total_energy - float(r.mol.get_scf_energy()))
exc_round = max(abs(a - b) for a, b in zip(
    res.extras["oqp"]["state_energies_hartree"], st.energies.tolist()))
print("=== QCSchema ===")
print(f"  AtomicResult validated: success={res.success}  driver={res.driver}")
print(f"  scf energy round-trip |Δ| = {e_round:.2e}  (machine precision)")
print(f"  state energies round-trip |Δ| = {exc_round:.2e}")
print(f"  excitation energies (eV): {np.round(res.extras['oqp']['excitation_energies_ev'],3)}")
qok = bool(res.success) and e_round < 1e-12 and exc_round < 1e-12
ok &= qok

# ---------- FCIDUMP (H2/STO-3G) ----------
inp2 = "benchmark/inputs/h2_rhf_sto3g.inp"
r2 = Runner(project="h2_rhf_sto3g", input_file=inp2,
            log=os.path.join(RUNDIR, "h2_rhf_sto3g.log"), silent=1, usempi=False)
r2.run()
fcidump_path = os.path.join(RUNDIR, "h2_sto3g.fcidump")
meta = dump_fcidump(fcidump_path, r2.mol, source="oqp")
print("\n=== FCIDUMP (H2/STO-3G) ===")
print(f"  norb={meta['norb']} nelec={meta['nelec']} ms2={meta['ms2']}")
print(f"  OQP-vs-PySCF consistency: {meta['consistency']}")
print(f"  8-fold symmetry residual = {meta['sym_residual_8fold']:.2e}")
ver = verify_fcidump_fci(fcidump_path, r2.mol)
print(f"  E(FCIDUMP->FCI)   = {ver['e_fcidump']:.10f}")
print(f"  E(PySCF native FCI)= {ver['e_pyscf_fci']:.10f}")
print(f"  |Δ| = {ver['diff']:.2e}  (tol 1e-8)")
cons = meta["consistency"]
cons_ok = (cons["max_dS"] is None or cons["max_dS"] < 1e-6) and \
          (cons["max_dHcore"] is None or cons["max_dHcore"] < 1e-6) and \
          (cons["dEnuc"] is None or cons["dEnuc"] < 1e-8)
fok = ver["diff"] < 1e-8 and meta["sym_residual_8fold"] < 1e-10 and cons_ok
ok &= fok

print("\nGATE 5:", "PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
