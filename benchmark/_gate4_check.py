"""GATE 4 validation: cube export + grid-integral trace checks."""
import os, sys
import numpy as np
from oqp.pyoqp import Runner
from oqp.analysis import MRSFExcitedStates, AOBasis, attachment_detachment, nto_excitation
from oqp.export import CubeExporter

INP = sys.argv[1] if len(sys.argv) > 1 else "benchmark/inputs/ch2o_mrsf.inp"
TARGET = int(sys.argv[2]) if len(sys.argv) > 2 else 1
FINE = float(sys.argv[3]) if len(sys.argv) > 3 else 0.05
RUNDIR = "/bighome/vova/coding/oqp-misc/runs/gate4"
os.makedirs(RUNDIR, exist_ok=True)
proj = os.path.splitext(os.path.basename(INP))[0]
r = Runner(project=proj, input_file=INP, log=os.path.join(RUNDIR, proj + ".log"),
           silent=1, usempi=False)
r.run()
st = MRSFExcitedStates(r.mol)
ao = AOBasis(r.mol)
cub = CubeExporter(st, ao, padding=5.0, spacing=0.15)   # viz grid
print(f"=== {proj} cubes: box {cub.n}, spacing {cub.dvec[0]} (viz), fine={FINE} (trace) ===")
N = st.n_elec
checks = []

# write cubes (viz grid) and report their (coarse) integrals
homo = st.na - 1
imo = cub.mo_cube(os.path.join(RUNDIR, f"{proj}_mo_homo.cube"), homo)
isd = cub.state_density_cube(os.path.join(RUNDIR, f"{proj}_state_S{TARGET}.cube"), TARGET)
itd = cub.transition_density_cube(os.path.join(RUNDIR, f"{proj}_trans_0_{TARGET}.cube"), 0, TARGET)
ad = attachment_detachment(st, TARGET, ref=0)
ia, idet = cub.attachment_detachment_cubes(
    os.path.join(RUNDIR, f"{proj}_attach_S{TARGET}.cube"),
    os.path.join(RUNDIR, f"{proj}_detach_S{TARGET}.cube"), ad)
ntoe = nto_excitation(st, TARGET)
ihole = cub.nto_cube(os.path.join(RUNDIR, f"{proj}_nto_hole_S{TARGET}.cube"),
                     ntoe["holes_ao"][:, 0], f"NTO hole 0 S{TARGET}")
ipart = cub.nto_cube(os.path.join(RUNDIR, f"{proj}_nto_part_S{TARGET}.cube"),
                     ntoe["particles_ao"][:, 0], f"NTO particle 0 S{TARGET}")
print(f"cubes written to {RUNDIR}/")

# --- rigorous trace checks on a fine grid ---
sd = cub.integrate_density(st.state_density_ao(TARGET), spacing=FINE)
td = cub.integrate_density(st.tdm_ao(0, TARGET), spacing=FINE)
at = cub.integrate_density(ad["A_ao"], spacing=FINE)
dt = cub.integrate_density(ad["D_ao"], spacing=FINE)
mo2 = cub.integrate_orbital_sq(st.C[:, homo], spacing=FINE)
hsq = cub.integrate_orbital_sq(ntoe["holes_ao"][:, 0], spacing=FINE)
psq = cub.integrate_orbital_sq(ntoe["particles_ao"][:, 0], spacing=FINE)

print(f"\n[trace] state S{TARGET} density   int rho = {sd:.5f}   (N={N})       |d|={abs(sd-N):.2e}")
print(f"[trace] transition 0->{TARGET}    int rho = {td:.2e}  (expect 0)    |d|={abs(td):.2e}")
print(f"[trace] attachment            int rho = {at:.5f}   (n_prom={ad['n_promoted']:.5f}) |d|={abs(at-ad['n_promoted']):.2e}")
print(f"[trace] detachment            int rho = {dt:.5f}   (n_prom={ad['n_promoted']:.5f}) |d|={abs(dt-ad['n_promoted']):.2e}")
print(f"[trace] |HOMO|^2={mo2:.5f}  |NTO hole|^2={hsq:.5f}  |NTO part|^2={psq:.5f}  (expect 1)")

checks = [
    ("state_density_N", abs(sd - N) < 1e-2),
    ("transition_zero", abs(td) < 1e-2),
    ("attachment_nprom", abs(at - ad["n_promoted"]) < 1e-2),
    ("detachment_nprom", abs(dt - ad["n_promoted"]) < 1e-2),
    ("mo_norm", abs(mo2 - 1) < 1e-2),
    ("nto_hole_norm", abs(hsq - 1) < 1e-2),
    ("nto_part_norm", abs(psq - 1) < 1e-2),
]
ok = all(v for _, v in checks)
for nm, v in checks:
    print(f"   {nm}: {'PASS' if v else 'FAIL'}")
print(f"\n{proj} S{TARGET} GATE 4:", "PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
