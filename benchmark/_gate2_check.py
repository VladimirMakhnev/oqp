"""GATE 2 keystone validation: reconstruct MRSF transition dipoles from the
exposed 1-TDM and match OQP's own values; verify traces. Verify-first."""
import os, sys
import numpy as np

OPENQP_ROOT = os.environ["OPENQP_ROOT"]
INP = sys.argv[1] if len(sys.argv) > 1 else os.path.join(OPENQP_ROOT, "benchmark/inputs/ch2o_mrsf.inp")
RUNDIR = sys.argv[2] if len(sys.argv) > 2 else "/bighome/vova/coding/oqp-misc/runs/gate2"
os.makedirs(RUNDIR, exist_ok=True)
project = os.path.splitext(os.path.basename(INP))[0]
log = os.path.join(RUNDIR, project + ".log")

from oqp.pyoqp import Runner
r = Runner(project=project, input_file=INP, log=log, silent=1, usempi=False)
r.run()
mol = r.mol


def f_order(tag, shape):
    raw = np.array(mol.data[tag], copy=True).ravel(order="C")
    return raw.reshape(shape, order="F")


def unpack_lt(packed, n):
    m = np.zeros((n, n))
    rows, cols = np.tril_indices(n)
    m[rows, cols] = packed
    m[cols, rows] = packed
    return m


nbf = int(np.asarray(mol.data["OQP::VEC_MO_A"]).size ** 0.5)
na = int(np.asarray(mol.data["nelec_A"]).ravel()[0])
nb = int(np.asarray(mol.data["nelec_B"]).ravel()[0])
nstates = int(np.asarray(mol.data["OQP::td_energies"]).ravel().size)
N_elec = na + nb
print(f"nbf={nbf} na={na} nb={nb} nstates={nstates} N={N_elec}")

C = np.asarray(mol.data["OQP::VEC_MO_A"]).reshape(nbf, nbf).T          # C[ao, mo]
S = unpack_lt(np.asarray(mol.data["OQP::SM"]).ravel(), nbf)
en = np.asarray(mol.data["OQP::td_energies"]).ravel()
dip_oqp = f_order("OQP::td_trans_dipole", (3, nstates, nstates))
dip_ao_packed = f_order("OQP::td_dip_ao", (nbf * (nbf + 1) // 2, 3))
R = [unpack_lt(dip_ao_packed[:, i], nbf) for i in range(3)]            # AO dipole ints (sym)
T = f_order("OQP::td_trans_density_mo", (nbf, nbf, nstates, nstates))  # alpha-MO basis

# orthonormality of MOs: C^T S C = I  (validates SM unpacking + MO layout)
ortho_err = np.max(np.abs(C.T @ S @ C - np.eye(nbf)))
print(f"[check] max|C^T S C - I| = {ortho_err:.3e}")

# reference 1-RDM in MO basis (ROHF triplet): 2 for 1..nb, 1 for nb+1..na, else 0
gamma_ref_mo = np.zeros((nbf, nbf))
for i in range(nb):
    gamma_ref_mo[i, i] = 2.0
for i in range(nb, na):
    gamma_ref_mo[i, i] = 1.0
print(f"[check] Tr(gamma_ref_mo) = {np.trace(gamma_ref_mo):.6f} (expect {N_elec})")


def reconstruct_mu(gamma_mo, mode):
    if mode == "CDCt":
        g_ao = C @ gamma_mo @ C.T
    else:  # CtDC
        g_ao = C.T @ gamma_mo @ C
    return np.array([-np.sum(g_ao * R[i]) for i in range(3)])


print("\n  pair    |  OQP dipole (x,y,z)            |  recon(CDCt)                   max|Δ|")
worst = 0.0
checks = []
for ist in range(nstates):
    for jst in range(nstates):
        if ist >= jst:
            continue
        gamma = T[:, :, ist, jst]
        mu_o = dip_oqp[:, ist, jst]
        mu_r = reconstruct_mu(gamma, "CDCt")
        d = np.max(np.abs(mu_o - mu_r))
        worst = max(worst, d)
        tr = np.trace(gamma)
        checks.append((f"{ist+1}->{jst+1}", d, tr))
        print(f"  {ist+1}->{jst+1}   | {mu_o[0]:+.6f} {mu_o[1]:+.6f} {mu_o[2]:+.6f} | "
              f"{mu_r[0]:+.6f} {mu_r[1]:+.6f} {mu_r[2]:+.6f}   {d:.2e}   tr(gamma)={tr:+.2e}")

print(f"\n[GATE2] worst |Δμ| over all transition pairs (CDCt) = {worst:.3e}  (tol 1e-6)")

# also test the other orientation, for the record
worst_alt = 0.0
for ist in range(nstates):
    for jst in range(nstates):
        if ist >= jst:
            continue
        mu_o = dip_oqp[:, ist, jst]
        mu_r = reconstruct_mu(T[:, :, ist, jst], "CtDC")
        worst_alt = max(worst_alt, np.max(np.abs(mu_o - mu_r)))
print(f"[info] worst |Δμ| with C^T D C orientation = {worst_alt:.3e}")

# state-density trace check: Tr(gamma_n_AO S) = N for each root
print("\n[check] state densities Tr(gamma_n_AO . S):")
worst_tr = 0.0
for n in range(nstates):
    delta = T[:, :, n, n]                 # traceless difference density
    gamma_n_mo = gamma_ref_mo + delta
    g_ao = C @ gamma_n_mo @ C.T
    trval = np.sum(g_ao * S)              # Tr(g_ao S)
    worst_tr = max(worst_tr, abs(trval - N_elec))
    print(f"   root {n+1}: Tr(Delta)={np.trace(delta):+.2e}  Tr(gamma_n_AO.S)={trval:.6f} (N={N_elec})")
print(f"[GATE2] worst |Tr(gamma_n_AO.S) - N| = {worst_tr:.3e}  (tol 1e-6)")

ok = (ortho_err < 1e-8) and (worst < 1e-6) and (worst_tr < 1e-6) and all(abs(t) < 1e-8 for _, _, t in checks)
print("\nGATE 2 RESULT:", "PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
