# PROGRESS — OQP Excited-State Analysis & Export Suite

Branch: `misc-excited-analysis`. Source of truth: `oqp_misc_TASK_BRIEF.md` (re-read at
every phase / resume / before relaxing any check). On resume, read this file first
and continue from the last green gate.

## Environment (justified deviation from brief §2)
The brief specifies a conda env `oqp-misc` / Python 3.11. **This cluster has no conda.**
Per the cluster-verified `oqp-setup` skill, OQP builds with environment *modules* +
a *venv named after the working folder*. So the env is a **venv `oqp-misc`** (= folder
name, honoring "create oqp-misc, do not reuse"), Python **3.9.6**, GCC 11.2.0.

Activate everything (modules don't persist between shells — chain in one command):
```
source /etc/profile.d/modules.sh && module purge && \
module load Python/3.9.6-GCCcore-11.2.0 GCC/11.2.0 OpenMPI/4.1.1-GCC-11.2.0 CMake/3.26.3-GCCcore-11.2.0 && \
source /bighome/vova/coding/oqp-misc/oqp-misc/bin/activate && \
export OPENQP_ROOT=/bighome/vova/coding/oqp-misc/openqp
```
Build: `cmake -B build -G Ninja -DUSE_LIBINT=OFF -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_INSTALL_PREFIX=. -DENABLE_OPENMP=ON`
(NOTE: do **not** pass `-DLINALG_LIB_INT64=OFF` — broken on current main; default ILP64
links system `/usr/lib64/libopenblas64.so`.) Then `ninja -C build install`; `cd pyoqp && pip install .`.
Run single-core: `openqp file.inp --omp 1`. Known flaky: native TRAH stability check can
randomly overflow after SCF converges — workaround `trh_impl=otr` or `stability=False` in `[scf]`.

Deps installed in venv: numpy 1.26.4, scipy, cclib 1.8.1, qcelemental 0.30.1, qcengine,
pyscf 2.13.1, pytest, h5py.

## Gate status
- [x] **GATE 0** — build clean; stock MRSF example runs; test baseline recorded. (commit: pending)
- [~] **GATE 1** — inventory (see below). Essentially complete from code reading.
- [x] **GATE 2** — KEYSTONE GREEN. Formaldehyde MRSF BHHLYP/6-31G* (nbf=34, N=16, 4 states):
      `max|CᵀSC−I|=1.1e-14`; transition-dipole reconstruction (C·D·Cᵀ) worst |Δµ|=**2.3e-15**
      ≤1e-6 (Cᵀ·D·C orientation → 3.8, wrong); all TDMs traceless (~1e-17);
      `Tr(γⁿ_AO·S)=16.000000=N` (worst 7e-15). Check: `benchmark/_gate2_check.py`. (commit: pending)
- [x] **GATE 3** — GREEN on formaldehyde S1 (n→π*) and ethylene S1. `pyoqp/oqp/analysis/`:
      transition_density, nto, density_diff, descriptors, gto_grid. Amplitude extraction
      validated against the keystone to ~1e-17; NTO-trans reconstructs µ to 4e-16; A−D=Δ,
      Tr(A·S)=Tr(D·S)=n_promoted≈1; PR∈[1,Npair]; Λ=0.998 (local, low-CT ✓); Ω total=||X||².
      Self-contained GTO evaluator: analytic overlap matches OQP S to 6.7e-16, grid-S 2.6e-4
      at 0.06 bohr (≤1e-3). Check: `benchmark/_gate3_check.py`. (commit: pending)
- [x] **GATE 4** — GREEN. `pyoqp/oqp/export/cubegen.py` writes Gaussian cubes (MO, state
      density, transition density, attach/detach, NTO hole/particle). Trace checks (fine
      0.05-bohr grid): state→16.00068 (N=16, Δ=6.8e-4), transition→5e-17, attach/detach→
      n_promoted (Δ<2e-4), |orbital|²≈1, all ≤1e-2. Check: `benchmark/_gate4_check.py`. (commit: pending)
- [x] **GATE 5** — GREEN. `pyoqp/oqp/export/qcschema.py` (AtomicResult validates; SCF+state
      energies round-trip to 0.0) and `fcidump.py` (H2/STO-3G: OQP S/Hcore/E_nuc match PySCF
      to ~1e-9; 8-fold residual 5.6e-16; **FCIDUMP→PySCF-FCI vs native FCI = 1.49e-9 ≤ 1e-8**).
      Check: `benchmark/_gate5_check.py`. (commit: pending)
- [ ] GATE 6 — interop
- [ ] GATE 7 — regression + benchmark
- [ ] GATE 8 — report

## Current gate & rules in force (self-check after brief re-read)
Working GATE 2 (the keystone). Rules I must keep obeying: verify-first (numeric probe, not
code audit); red-gate stop; Python-first/surgical-Fortran (Fortran touched only to write the
1-TDM intermediates into the tagarray); respect the bridge (tagarray tags, no parallel I/O);
new code in new files; log every edit in CHANGELOG_misc.md; commit per green gate. Not drifted.

## GATE 1 — Inventory (the bridge)
The Fortran↔Python store is **tagarray**: `OQPData.__getitem__(key)` (pyoqp/oqp/molecule/oqpdata.py)
falls through to `lib.oqp_get(...)` returning a zero-copy numpy view of the named array;
`__setitem__` writes via `lib.oqp_alloc`+memmove. `mol.data['OQP::<name>']` reads any tag.
Tag-name constants are defined in `source/tagarray_driver.F90` (`OQP_<name> = OQP_prefix//"<name>"`,
prefix `OQP::`); a routine reserves with `infos%dat%reserve_data(NAME, TA_TYPE_REAL64, total, dims, comment=)`
then `tagarray_get_data(infos%dat, NAME, ptr)` and copies into `ptr`.

Key conventions (verified against pyoqp + Fortran):
- **MO coefficients** `OQP::VEC_MO_A` / `OQP::VEC_MO_B`: Fortran column-major; numpy view is
  `arr.reshape(nbf,nbf).T` → C[ao, mo] (see molecule.py `_mo_coefficients`). ROHF: alpha & beta
  share the same spatial orbitals.
- **Overlap** `OQP::SM`: packed **lower-triangle** length nbf*(nbf+1)/2 (unpack with tril_indices).
- **Response amplitudes** `OQP::td_bvec_mo`: shape (xvec_dim, nstates), occupied-index fastest;
  for sf/mrsf the spin-flip block is alpha-occ × beta-vir.
- **Excitation energies** `OQP::td_energies` (nstates).
- **MRSF transition dipoles**: printed in the .log "Summary table" and "Transition" table; now
  ALSO exposed full-precision (see GATE 2). Computed in `tdhf_mrsf_energy.F90`.

Keystone location (THE exposure point):
- `source/modules/tdhf_mrsf_energy.F90` builds `trden(nbf,nbf,nstates,nstates)` via
  `get_mrsf_transition_density` → `get_trans_den` (tdhf_mrsf_lib.F90:2583/2672), then computes
  transition dipoles with `get_transition_dipole` (tdhf_sf_lib.F90:1057).
- `trden(:,:,ist,jst)` is in the **alpha-MO basis**. Off-diagonal (ist<jst) = state-interaction
  1-TDM γ^{ist→jst}; diagonal (ist==jst) = traceless **difference** 1-RDM Δγ = γ_state − γ_ref.
  **Only occ-occ and vir-vir blocks are nonzero (occ-vir = 0)** — the hallmark of a transition
  between two *response roots* (S0 is root 1, not the reference). Only the upper triangle
  (ist≤jst) is filled; γ^{jst→ist} = (γ^{ist→jst})^T.
- Dipole (OQP's own): trden_ao = C·trden_mo·Cᵀ (alpha C), symmetrize, then
  µ_i = −Tr(γ_ao · r_ao,i), r = AO dipole integrals at center of mass.
- Reference occupations (ROHF triplet, MO basis): orbitals 1..nocb doubly occ, nocb+1..noca
  singly occ (SOMOs), where noca=nelec_A, nocb=nelec_B=noca−2 for the triplet. γ_state = γ_ref + Δγ.
- Gradient/relaxed density: `tdhf_mrsf_z_vector.F90` writes `OQP::td_mrsf_density` (7,nbf,nbf)
  (spin-pair-coupling density components for gradients) and `OQP::td_p`; the relaxed 1-RDM for
  attach/detach is reachable via the z-vector/gradient path (deferred; unrelaxed first per brief).

## GATE 2 — exposure patch (additive only; no physics changed)
New tags written at end of `tdhf_mrsf_energy`:
- `OQP::td_trans_density_mo` (nbf,nbf,nstates*nstates) — trden, pair k=ist+(jst−1)*nstates (col-major).
  Python: `flat.reshape((nbf,nbf,nstates,nstates), order='F')` then `[:,:,i-1,j-1]`.
- `OQP::td_trans_dipole` (3,nstates,nstates) — OQP's computed transition dipoles (a.u.).
- `OQP::td_dip_ao` (nbf2,3) — AO dipole integrals at center of mass (packed L-tri).
Validation plan: reconstruct µ^{1→n} = −Tr(C trden_mo Cᵀ · r_ao) and match `OQP::td_trans_dipole`
≤1e-6 for ≥2 states; Tr(γ^{1→n}_MO)≈0; CᵀSC=I; Tr(γ^n_MO)=N.

## Open questions / risks
- MO→AO transform orientation (C·D·Cᵀ vs Cᵀ·D·C) — resolved numerically by the dipole match.
- NTOs on a state-interaction TDM (no occ-vir block): design in GATE 3; will document deviation
  from the brief's closed-shell NTO wording with the physical reason.
- Overleaf git auth (GATE 8): may require token; fallback = commit locally + record push commands.

## Reference numbers (H2O BHHLYP/6-31G* MRSF, stock example, --omp 1)
SCF total = −76.0774468204. MRSF S0 (root1) = −76.3609364927. Transition dipoles (printed, a.u.):
1→2 (0.1790,0.1790,0) |µ|=0.2532 osc 0.0140 @8.947 eV; 2→3 (−1.2950,1.2950,0) |µ|=1.8315 osc 0.1303.
