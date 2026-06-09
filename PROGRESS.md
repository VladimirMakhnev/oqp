# PROGRESS.md — SSC/ZFS project state (durable memory across `claude -p` runs)

> Read `CLAUDE.md` first. This file is the single source of truth for "what is done / what is
> next." Update it at the end of every run: mark task status, record current numbers vs target,
> set the **NEXT STEP**. Status keys: ☐ todo · ◐ in-progress · ☑ done (gate passed) · ⮕ blocked.

## CURRENT STATUS

- **Branch:** `ssc-zfs` (off `SSC` @ baseline 222/225 tests passing).
- **Phase:** 1 — L1 integral work in progress. P1.1 done (route decided + closed form derived and
  numerically validated against analytic oracles). Native Fortran integral NOT started.
- **Gate cleared:** none yet. (P1.1 is a decision/prototype task; it is NOT a stage gate.)
- **NEXT STEP:** **P1.2** — implement the SS dipolar 2e integral in OpenQP via **Path A** (drive
  the ERI engine with the Hessian-of-1/r₁₂ kernel) for general angular momenta, reusing
  `rys_deriv.F90` / `grd2_rys.F90`. Then **P1.3** = the real L1 gate: FD vs OpenQP's *actual* ERI
  engine (the prototype below already proves the closed form is right; P1.3 proves the Fortran is).

---

## PHASE 0 — Scaffolding  ☑ (this setup session)
- ☑ Branch `ssc-zfs` created.
- ☑ `CLAUDE.md`, `PROGRESS.md`, `benchmarks.md`, `QUESTIONS.md` written.
- ☑ `tests/test_ssc_integrals_fd.py` skeleton (L1 placeholder, `xfail` until integrals exist).
- ☑ `docs/ssc_zfs_derivation.tex` skeleton + pushed to Overleaf.
- ☑ Reuse map recorded in `CLAUDE.md §5`.

## PHASE 1 — L1: 2e SS integral + FD validation  (gate: §7 L1)
- ☑ **P1.1** **Path A decided** (drive the ERI engine with the Hessian-of-1/r₁₂ kernel; the SS 2e
  integral is exactly that Hessian → reuse ERI machinery, no new Rys primitive needed for the
  validation). Closed form for the (s,s,s,s) quartet derived and numerically validated by a
  standalone prototype (`tests/ssc_prototype_ssss.py`), triangulated three ways. **NOT a stage
  gate** (validated vs an analytic oracle, not OpenQP's ERI engine — that is P1.3).
  - Derivation: bare-Hessian integral
    `H_kl = pref·K·[4ρ²R_kR_l F₂(T) − 2ρδ_kl F₁(T)]`, `pref=2π^{5/2}/(pq√(p+q))`, `R=P−Q`, `T=ρR²`.
  - Physical dipolar integral = **traceless part** of H: `S = H − ⅓Tr(H)·I = H + (4π/3)O·I`,
    `O=⟨ρ₁|ρ₂⟩=K(π/(p+q))^{3/2}e^{−T}`. Distributional identity:
    `∂_k∂_l(1/r) = (3r_kr_l−δ_kl r²)/r⁵ − (4π/3)δ_kl δ³(r)`; the contact term is the isotropic
    part removed by tracelessness (does not enter the ZFS D-tensor).
  - Validation numbers (run in `ssc-pyenv`): closed-form H **vs** Richardson-FD-of-ERI rel **2.3e-9**;
    closed-form H **vs** Boys-free t-quadrature rel **9.9e-14**; `Tr(S)=−3.9e-16`;
    `Tr(H)=−1.904324 = −4πO` (contact identity) ✓. All four self-checks PASS.
  - Correction logged: the t-quadrature (Gaussian transform) reproduces **H** (contact included),
    not S — same as the FD route; both equal the closed-form H. S is then traceless(H).
- ☐ **P1.2** Implement all 6 components for general angular momenta (reuse `rys_deriv.F90` /
  `grd2_rys.F90` or `comp_soc_int2_prim`).
- ☐ **P1.3** L1 FD test: compare to finite differences of the ERI engine to 6–8 sig figs;
  assert `trace = 0` to ~1e-10. **Gate L1.**

## PHASE 2 — L2: contraction + ROHF, pin `C`  (gate: §7 L2)
- ☐ **P2.1** Build the `{P_μν P_κτ − P_μκ P_ντ}` contraction as a Fock/K-like consumer
  (template: `int2_mrsf_data_t`, `tdhf_mrsf_lib.F90`). Input: ROHF `P^(α−β)` (M_S = S).
- ☐ **P2.2** Assemble the 6-component D-tensor; diagonalise → `D`, `E`, `E/D`; unit a.u.→cm⁻¹.
- ☐ **P2.3** **Pin `C` and sign NUMERICALLY on O₂ ³Σ_g⁻ @ 1.207 Å** (target `D^SS ≈ 1.44–1.6
  cm⁻¹`). Record pinned value + match in the LaTeX doc. Cross-check CH₂ ³B₁. **Gate L2.**

## PHASE 3 — L3: MRSF densities  (gate: §7 L3)  — DO NOT START before L2 passes
- ☐ **P3.1** Feed MRSF `P^(α−β)` (M_S = S via Wigner–Eckart, reuse `compute_tdm`) into the L2
  machinery. Add the UMRSF density-source flag (RO/UNO default).
- ☐ **P3.2** Reproduce the acene/radical table; target RMSD ≈ 0.0035 cm⁻¹ (RO-type). **Gate L3.**
- ☐ **P3.3** Wire Python dispatch (`runtype 'zfs'/'ssc'`), input checker, regression tests.

## OUT OF SCOPE (do not implement) — see `CLAUDE.md §1`
Z-vector / relaxed densities, response/relaxation terms, analytic gradients of D, SOC cross terms.

---

## RUNNING LOG  (newest first — one short entry per `-p` run)
- 2026-06-09 — **P1.1 done.** Decided **Path A** (SS 2e integral = Hessian of 1/r₁₂ → reuse ERI
  engine). Derived the closed form for the (s,s,s,s) quartet and wrote a standalone prototype
  `tests/ssc_prototype_ssss.py` triangulating it three independent ways: closed form (Boys F₁,F₂),
  Richardson-FD of the Coulomb ERI (independent erf-based F₀), and a Boys-free Gaussian-transform
  t-quadrature. Agreement: FD rel 2.3e-9, quad rel 9.9e-14; Tr(S)=−3.9e-16 (traceless);
  Tr(H)=−4πO contact identity confirmed. Corrected an earlier misderivation (t-quadrature gives the
  bare Hessian H, not the traceless S; S = traceless(H)). **L1 gate (P1.3) NOT cleared** — this
  validates the math vs an analytic oracle, not vs OpenQP's ERI engine. Updated derivation .tex.
  NEXT: P1.2 (native Fortran integral via Path A), then P1.3 (L1 gate vs the real ERI engine).
- 2026-06-09 — Setup session. Context gathered from `./papers/` (Sinnecker–Neese 2006 eq 9,
  Neese 2007 eq 46 — prefactor discrepancy logged in `CLAUDE.md §3`; Pokhilko–Krylov 2019 W–E
  extraction; Neese JACS 2006 mean-field). Reuse map built. Scaffolding committed. No code yet.
  NEXT: P1.1.

## OPEN ASSUMPTIONS (promote blocking ones to QUESTIONS.md)
- Final `C`/sign deferred to numerical pinning on O₂ (L2) — see `CLAUDE.md §3`. No assumption made yet.
