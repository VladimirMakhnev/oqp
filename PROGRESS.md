# PROGRESS.md — SSC/ZFS project state (durable memory across `claude -p` runs)

> Read `CLAUDE.md` first. This file is the single source of truth for "what is done / what is
> next." Update it at the end of every run: mark task status, record current numbers vs target,
> set the **NEXT STEP**. Status keys: ☐ todo · ◐ in-progress · ☑ done (gate passed) · ⮕ blocked.

## CURRENT STATUS

- **Branch:** `ssc-zfs` (off `SSC` @ baseline 222/225 tests passing).
- **Phase:** 1 — L1 integral implemented (P1.2 done). P1.3 self-test PASSES in testing; the L1
  **stage gate is NOT yet declared** — awaiting human confirmation (per session instruction).
- **Gate cleared:** **L1 ☑ (s,p,d; human-confirmed for s,p 2026-06-09; d extended same day).**
- **NOW:** starting **Phase 2 / P2.1** (the `{P_μν P_κτ − P_μκ P_ντ}` contraction + ROHF density).
  Per instruction: do NOT attempt the L2 O₂ pin yet — stop before L2.
- **NEXT STEP:** await confirmation of L1. On confirmation, mark L1 ☑ and begin **Phase 2 / P2.1**
  (the `{P_μν P_κτ − P_μκ P_ντ}` contraction + ROHF; then pin `C` on O₂ at L2). Possible follow-up
  before L2: extend the SS integral / FD self-test to **d shells** (currently s,p validated;
  spherical-harmonic d deferred).

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
- ☑ **P1.2** Native Fortran SS dipolar 2e integral implemented for general angular momenta.
  - New isolated module `source/integrals/mod_ssc_int2.F90` (does not perturb the tested SOC path):
    - `qgauss_ss` — padded Rys 2e 1D-table builder (electron 1 AND electron 2 padded +1), modelled
      on `QGaussRys2e`/`comp_soc_int2_prim` (`mod_1e_primitives.F90`); supports an electron-2 rigid
      operator displacement for the FD reference.
    - `comp_ssc_int2_prim` — assembles the 6 bare-Hessian components `H_kl` via the working identity
      `H_kl = −⟨∂_k(μν)|1/r₁₂|∂_l(κτ)⟩` (one first-derivative on each electron; reuses the SOC-style
      `soc_xyz_ij` derivative pattern as `e1d`/`e2d`/`e12d`). Physical `S = traceless(H)`.
    - `comp_eri2_prim_disp` — plain ERI with electron-2 displaced (the FD reference).
  - Key correctness fact pinned: operator displacement must enter **only** the Boys argument and
    the VRR centres, **not** the engine's `expe` Gaussian prefactor (that prefactor is a fixed-orbital
    normalisation; letting `dshift` leak into it injects a spurious `−2·ERI` second-derivative term).
- ☑ **P1.3 / Gate L1 — PASSED (human-confirmed 2026-06-09; s,p,d all validated).**
  - **d-shell extension (2026-06-09):** OpenQP uses **cartesian** Gaussians (`basis%naos = NUM_CART_BF`;
    6d, 10f) so a d shell is 6 cartesian functions and `CART_X(i,2)` gives their powers directly —
    **no spherical-harmonic transform exists in the integral engine**, so d needed only lifting the
    `am≤1` restriction (the Rys/derivative code is angular-momentum-general). Self-test now covers
    s,p,d quartets of H₂O/6-31G*: **87846/87846** comparisons agree at rel ≤ 1e-6; **worst rel diff
    over non-negligible blocks = 1.5e-8** (≈7–8 sig figs); worst overall 9.7e-8 is a vanishing-by-
    symmetry block (judged against a 1e-9 absolute floor); **worst |Tr(S)| = 2.2e-16**. The FD-block
    refactor (compute each displaced ERI block once per step) keeps it fast (~3 s).
  L1 FD self-test: `source/modules/ssc_int2_selftest.F90` (`bind(C)`; declared in
  `include/oqp.h`; driven by `tests/test_ssc_integrals_fd.py` via `oqp.ssc_int2_selftest`).
  Compares analytic `H_kl` to a **3-level Richardson FD of the engine's own ERI** for all s,p shell
  quartets of H₂O/6-31G*, and checks `Tr(S)=0`. **RESULT (run 2026-06-09, ssc-pyenv):**
  - **7776 / 7776** element comparisons agree at rel ≤ 1e-6; **worst rel diff = 1.29e-9**.
  - **worst |Tr(S)| = 1.3e-18** (traceless invariant, machine zero).
  - Independent machine-precision check: one-center `(ss|ss)` ratio `H_xx/ERI(0) = −2α/3` reproduced
    exactly vs the Python prototype, for every tested exponent.
  - Pathologically tight core primitives (exp > 100) excluded from FD (operator-displacement FD is
    roundoff-limited there; analytic path identical, covered by the prototype check). d shells
    deferred (spherical-harmonic transform). **Gate L1 — NOT declared; awaiting confirmation.**

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
- 2026-06-09 — **L1 confirmed (s,p) + extended to d; gate ☑.** Human confirmed L1 for s,p. Found
  OpenQP is fully **cartesian** (no spherical transform), so extended the SS integral + FD self-test
  to **d** shells by lifting the `am≤1` cap and refactoring the FD to compute displaced-ERI blocks
  once per step (fast). Result: **87846/87846** agree ≤1e-6, worst non-negligible **1.5e-8**,
  worst |Tr(S)| **2.2e-16**. Marked L1 ☑ (s,p,d). NEXT: P2.1 contraction + ROHF (stop before L2).
- 2026-06-09 — **P1.2 done; P1.3 self-test PASSES (gate not declared).** Implemented the native
  SS dipolar 2e integral (`source/integrals/mod_ssc_int2.F90`) via Path A and an `bind(C)` L1
  self-test (`source/modules/ssc_int2_selftest.F90`, `include/oqp.h`, `tests/test_ssc_integrals_fd.py`).
  The self-test finite-differences OpenQP's own ERI engine (3-level Richardson, operator displaced
  via rigid electron-2 shift) and compares to the analytic SS integral: **7776/7776 components agree
  to rel ≤ 1e-6, worst 1.29e-9; Tr(S)=1.3e-18**; one-center `H_xx/ERI(0)=−2α/3` matches the prototype
  exactly. Debug story: found+fixed a spurious `exp(−|P−Q|²)` term leaking the operator displacement
  into the prefactor (froze `expe` at the undisplaced geometry); core primitives excluded from FD
  (roundoff), d shells deferred. **Per instruction, L1 NOT marked passed — showing numbers, awaiting
  human confirmation.** NEXT: confirm L1 → Phase 2 (P2.1 contraction + ROHF).
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
