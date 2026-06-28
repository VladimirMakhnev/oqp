"""GATE 3 validation: NTOs, attachment/detachment, descriptors. Verify-first."""
import os, sys
import numpy as np
from oqp.pyoqp import Runner
from oqp.analysis import (MRSFExcitedStates, nto_excitation, nto_transition,
                          attachment_detachment, participation_ratio,
                          tozer_lambda, fragment_ct_matrix, AOBasis, make_box_grid)

INP = sys.argv[1]
TARGET = int(sys.argv[2]) if len(sys.argv) > 2 else 1     # 0-based root index
FRAG = eval(sys.argv[3]) if len(sys.argv) > 3 else None
RUNDIR = "/bighome/vova/coding/oqp-misc/runs/gate3"
os.makedirs(RUNDIR, exist_ok=True)
proj = os.path.splitext(os.path.basename(INP))[0]
r = Runner(project=proj, input_file=INP, log=os.path.join(RUNDIR, proj + ".log"),
           silent=1, usempi=False)
r.run()
st = MRSFExcitedStates(r.mol)
ao = AOBasis(r.mol)
print(f"\n=== {proj}: nbf={st.nbf} na={st.na} nb={st.nb} nstates={st.nstates} target=S{TARGET} ===")
print("excitation energies (eV, rel S0):",
      np.round((st.energies - st.energies[0]) * 27.211386, 3))

ok = True

# --- amplitude extraction validated against the keystone difference density ---
amp_err = 0.0
for n in range(st.nstates):
    X = st.amplitude_matrix(n)
    rec = st.trans_den_from_amplitudes(X, X)
    amp_err = max(amp_err, np.max(np.abs(rec - st.diff_density_mo(n))))
# cross-state TDM from amplitudes vs exposed TDM
tdm_amp_err = 0.0
for i in range(st.nstates):
    for j in range(i + 1, st.nstates):
        rec = st.trans_den_from_amplitudes(st.amplitude_matrix(i), st.amplitude_matrix(j))
        tdm_amp_err = max(tdm_amp_err, np.max(np.abs(rec - st.tdm_mo(i, j))))
print(f"[amp]  max|trans_den(X,X) - Delta gamma| = {amp_err:.2e}  (tol 1e-9)")
print(f"[amp]  max|trans_den(Xi,Xj) - gamma^ij|  = {tdm_amp_err:.2e}  (tol 1e-9)")
ok &= amp_err < 1e-9 and tdm_amp_err < 1e-9

# --- NTO (excitation) ---
ntoe = nto_excitation(st, TARGET)
ssum = ntoe["sum_sigma2"]
sig2 = (ntoe["sigma"] ** 2).sum()
nto_consist = abs(ssum - sig2)
Xnorm = (ntoe["amplitude_matrix"] ** 2).sum()
print(f"[NTO-exc S{TARGET}] sum(sigma^2)={ssum:.8f}  ||X||^2={Xnorm:.8f}  "
      f"|sum sig^2 - ||X||^2|={abs(ssum-Xnorm):.2e}  n_sig={ntoe['n_significant']}")
ok &= abs(ssum - Xnorm) < 1e-10

# --- NTO (transition) reconstructs the dipole ---
ntot = nto_transition(st, 0, TARGET)
g_full = ntot["reconstruct_tdm"]()             # full rank
mu_recon = np.array([-np.sum((st.C @ g_full @ st.C.T) * st.R[k]) for k in range(3)])
mu_oqp = st.dip_oqp[:, 0, TARGET]
dmu = np.max(np.abs(mu_recon - mu_oqp))
print(f"[NTO-trans 0->{TARGET}] sum(sig^2)={ntot['sum_sigma2']:.6f}  "
      f"|mu_recon - mu_OQP| = {dmu:.2e} (tol 1e-8)   mu_OQP={np.round(mu_oqp,4)}")
ok &= dmu < 1e-8

# --- attachment / detachment ---
ad = attachment_detachment(st, TARGET, ref=0)
recon_delta = ad["A_mo"] - ad["D_mo"]
ad_err = np.max(np.abs(recon_delta - ad["delta_mo"]))
trAS = np.sum(ad["A_ao"] * st.S)
trDS = np.sum(ad["D_ao"] * st.S)
print(f"[A/D S{TARGET}] n_promoted={ad['n_promoted']:.6f}  Tr(A.S)={trAS:.6f}  Tr(D.S)={trDS:.6f}")
print(f"[A/D] max|A-D - Delta|={ad_err:.2e} (tol 1e-9)   |Tr(A.S)-Tr(D.S)|={abs(trAS-trDS):.2e}")
ok &= ad_err < 1e-9 and abs(trAS - ad["n_promoted"]) < 1e-6 and abs(trDS - ad["n_promoted"]) < 1e-6

# --- descriptors ---
pr = participation_ratio(ntoe["weights"])
npair = ntoe["sigma"].size
print(f"[PR S{TARGET}] PR={pr:.4f}  (1 <= PR <= {npair})")
ok &= 1.0 - 1e-9 <= pr <= npair + 1e-9

lo, n3, dvec, pts = make_box_grid(ao.coords, padding=5.0, spacing=0.12)
dV = dvec[0] * dvec[1] * dvec[2]
lam, laminfo = tozer_lambda(ao, ntoe, pts, dV)
print(f"[Lambda S{TARGET}] Tozer Lambda={lam:.4f}  (in [0,1]; high=local, low=CT)")
ok &= -1e-9 <= lam <= 1.0 + 1e-6

if FRAG is not None:
    om = fragment_ct_matrix(st, ao, TARGET, FRAG)
    print(f"[Omega S{TARGET}] total={om['total']:.6f} (||X||^2={om['amplitude_norm']:.6f})  "
          f"CT_frac={om['ct_fraction']:.3f}")
    print(f"           hole_pop={np.round(om['hole_pop'],4)}  particle_pop={np.round(om['particle_pop'],4)}")
    print("           Omega=\n", np.round(om["Omega"], 4))
    ok &= abs(om["total"] - om["amplitude_norm"]) < 1e-8

print(f"\n{proj} S{TARGET} GATE 3:", "PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
