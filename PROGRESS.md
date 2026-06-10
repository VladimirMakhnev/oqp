# PROGRESS.md — SSC/ZFS project state (durable memory across `claude -p` runs)

> Read `CLAUDE.md` first. This file is the single source of truth for "what is done / what is
> next." Update it at the end of every run: mark task status, record current numbers vs target,
> set the **NEXT STEP**. Status keys: ☐ todo · ◐ in-progress · ☑ done (gate passed) · ⮕ blocked.

## CURRENT STATUS

- **Branch:** `ssc-zfs` (off `SSC` @ baseline 222/225 tests passing).
- **Phase:** 3 — entering L3 (MRSF densities). **L1 ☑, L2 ☑ both human-confirmed.**
- **Gate cleared:** **L1 ☑ (s,p,d).  L2 ☑ (O₂ D^SS = +1.503 cm⁻¹, C = Neese 2007 Eq. 46;
  human-confirmed 2026-06-10).** L1 gate now also enforces the absolute textbook-(ss|ss)-ERI
  normalisation check (ratio = 1.0; guards the `exp(−|P−Q|²)` regression class).
- **L2 RESULT (O₂ ³Σ_g⁻, r=1.207 Å, ROHF/6-31G*, stability=false):** axial **Dxx=Dyy=+0.1711,
  Dzz=−0.3422** a.u. (E/D=0); **D^SS = +1.5031 cm⁻¹** (target +1.44–1.6, positive ✓). **Pinned
  prefactor C = −g_e²α²/[16 S(2S−1)] = Neese 2007 Eq. 46**, matched numerically to **1.1%**
  (`C_pin/C_Eq46 = +1.011`). The Sinnecker–Neese Eq. 9 prefactor (+g_e²α²/[4S(2S−1)]) gives
  −6.01 cm⁻¹ (wrong sign, 4× too large) — it differs from the pinned C by exactly **−1/4**, the
  inter-paper discrepancy flagged in CLAUDE.md §3. **Convention resolved: use Eq. 46.**
- **NEXT STEP:** await L2 confirmation. Then (per CLAUDE.md): secondary L2 check CH₂ ³B₁
  ([CONFIRM] the reference first), then **Phase 3 / L3** (MRSF densities; Wigner–Eckart M_S=S).

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
- ☑ **P2.1 (built + structurally validated; absolute correctness deferred to L2).** Contraction
  driver `source/modules/ssc_zfs.F90`:
  - `compute_ssc_dtensor_raw` — loops all shell quartets, accumulates the contracted **cartesian**
    SS integral block over primitives (`comp_ssc_int2_prim`), makes it traceless (`S=H−⅓Tr(H)I`),
    and contracts with the ROHF spin density `P^(α−β)=DM_A−DM_B` (M_S=S, exact for a single
    determinant) in the Coulomb-like (`P_μν P_κτ`) **minus** exchange-like (`P_μκ P_ντ`) patterns,
    giving the 6 components. bfnrm absorbed by pre-scaling the density `Q_μν=P_μν·bfnrm_μ·bfnrm_ν`
    (lets the cartesian integrals be contracted directly).
  - `ssc_dtensor_selftest` (`bind(C)`, `include/oqp.h`, `tests/test_ssc_dtensor.py`).
  - **RESULT (O₂ ³Σ_g⁻ @ 1.207 Å, ROHF/6-31G*, C=1):** runs; **Tr(D)=−4.5e-15** (traceless
    invariant holds); off-diagonals ~0; raw components `Dxx=+0.521, Dyy=−0.248, Dzz=−0.273` (a.u.).
  - **FLAG for L2:** the raw tensor is **not axial about the molecular z-axis** (Dxx is the outlier,
    Dyy≈Dzz) — expected D for O₂ ³Σ_g⁻ is axial about z. This is almost certainly the **reference
    state** (plain Huckel→ROHF need not give the cylindrically-symmetric ³Σ_g⁻ π* occupation); it is
    an L2 concern (right state + `C`/unit pin), NOT a contraction bug (the 6 integral components are
    L1-validated vs FD; Tr(D)=0 holds). Resolve at L2 before trusting the magnitude.
- ☑ **P2.2** Diagonalise the 6-component D-tensor → principal values, `D`, `E`, `E/D`, a.u.→cm⁻¹.
  Symmetric-3×3 Jacobi (`jacobi3`) + ZFS ordering (`order_zfs`) in `ssc_zfs.F90`; report prints
  `D^SS`, `E^SS` (cm⁻¹) and `E/D`. O₂: principal values (+0.1711,+0.1711,−0.3422) a.u., E/D=0.
- ◐ **P2.3 / Gate L2 — PASSES in testing (O₂); not declared pending confirmation.**
  **Pinned `C = −g_e²α²/[16 S(2S−1)]` (Neese 2007, Eq. 46)**, baked into `ssc_zfs.F90`.
  O₂ ³Σ_g⁻ @ 1.207 Å, ROHF/6-31G* (stability=false): **D^SS = +1.5031 cm⁻¹** (band +1.44–1.6 ✓,
  sign + ✓), E/D=0. `C_pin/C_Eq46 = +1.011` (1.1%, within "few %"). Convention discrepancy
  (CLAUDE.md §3) **resolved**: Eq. 46, not Eq. 9 (which is off by −1/4 → −6.01 cm⁻¹).
  **Two bugs found+fixed reaching this** (see RUNNING LOG): a spurious `exp(−|P−Q|²)` ERI-prefactor
  factor (multi-centre normalisation) and an off-by-one AO index in the contraction.
  - **Secondary CH₂ ³B₁ sanity (NOT a gate — reference is [CONFIRM]):** bent triplet, ROHF/6-31G*,
    `D^SS = +0.693 cm⁻¹`, rhombic `|E/D| = 0.208` (in [0,1/3]). Positive, plausible magnitude
    (below the ~0.76 experimental *total* D, consistent with SS being dominant for a 1st-row
    carbene) and physically sensible rhombicity (in-plane a₁ + out-of-plane b₁). O₂ remains the
    binding anchor.

## PHASE 3 — L3: MRSF densities  (gate: §7 L3)  — DO NOT START before L2 passes
- ☐ **P3.1** Feed MRSF `P^(α−β)` (M_S = S via Wigner–Eckart, reuse `compute_tdm`) into the L2
  machinery. Add the UMRSF density-source flag (RO/UNO default).
- ☐ **P3.2** Reproduce the acene/radical table; target RMSD ≈ 0.0035 cm⁻¹ (RO-type). **Gate L3.**
- ☐ **P3.3** Wire Python dispatch (`runtype 'zfs'/'ssc'`), input checker, regression tests.

## OUT OF SCOPE (do not implement) — see `CLAUDE.md §1`
Z-vector / relaxed densities, response/relaxation terms, analytic gradients of D, SOC cross terms.

---

## RUNNING LOG  (newest first — one short entry per `-p` run)
- 2026-06-10 — **L2 pinned on O₂ (PASSES in testing; stopped after the pin).** Got the clean
  cylindrical ³Σ_g⁻ reference via **`scf.stability=false`** (the symmetric ROHF point is a saddle;
  OQP's stability-following otherwise escapes to a symmetry-broken non-axial state). Diagnosing the
  initial non-axiality uncovered **two real bugs**, both fixed: (1) an **off-by-one AO index** in the
  contraction (`mu=locao+i` → `locao+i-1`; `locao` is 1-based) which scrambled px/py/pz of the
  density — fixing it gave axial D and N⁴S→0; (2) a **spurious `exp(−|P−Q|²)`** factor in the Rys
  `expe` prefactor (copied from the untested, caller-less `QGaussRys2e`) that mis-normalised
  multi-centre ERIs — verified by a new textbook-ERI check (ratio 0.041→**1.0000000** after removal;
  L1 unaffected since one-centre x=0). Added P2.2 diagonalisation (`jacobi3`). **Result: O₂ D^SS =
  +1.5031 cm⁻¹** (band +1.44–1.6), E/D=0; **pinned C = Neese 2007 Eq. 46** to 1.1%; Eq. 9 off by
  −1/4 (CLAUDE.md §3 resolved). All SSC tests pass (L1 87846/87846; L2 O₂; prototype). STOP after
  the pin per instruction; await confirmation. NEXT: CH₂ ³B₁ secondary check, then L3 (MRSF).
- 2026-06-09 — **P2.1 contraction built (stopped before L2).** Implemented `source/modules/ssc_zfs.F90`
  (`compute_ssc_dtensor_raw` + `ssc_dtensor_selftest`): contracts the L1-validated SS integral with
  the ROHF spin density `DM_A−DM_B` (Coulomb − exchange), bfnrm absorbed via density pre-scaling.
  O₂ ³Σ_g⁻/6-31G* ROHF, C=1: runs, **Tr(D)=−4.5e-15** (traceless ✓), Dxx/Dyy/Dzz=+0.521/−0.248/−0.273.
  Flagged: raw tensor not axial about z → reference-state issue to fix at L2 (not a contraction bug).
  Per instruction, STOPPED before L2 (no `C` pin / no number match). NEXT: L2 (clean ³Σ_g⁻ + pin C).
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
