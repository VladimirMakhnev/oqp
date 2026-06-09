# PROGRESS.md — SSC/ZFS project state (durable memory across `claude -p` runs)

> Read `CLAUDE.md` first. This file is the single source of truth for "what is done / what is
> next." Update it at the end of every run: mark task status, record current numbers vs target,
> set the **NEXT STEP**. Status keys: ☐ todo · ◐ in-progress · ☑ done (gate passed) · ⮕ blocked.

## CURRENT STATUS

- **Branch:** `ssc-zfs` (off `SSC` @ baseline 222/225 tests passing).
- **Phase:** 0 — scaffolding complete; implementation NOT started (per session setup, deliberately).
- **Gate cleared:** none yet.
- **NEXT STEP (first `-p` task):** **P1.1** — choose integral route (Path A vs B) and stand up the
  L1 finite-difference integral unit test harness (`tests/test_ssc_integrals_fd.py`, skeleton
  already present) against the existing ERI engine. Do NOT pin `C` yet (that is L2); L1 only
  checks the *integral* values and trace = 0.

---

## PHASE 0 — Scaffolding  ☑ (this setup session)
- ☑ Branch `ssc-zfs` created.
- ☑ `CLAUDE.md`, `PROGRESS.md`, `benchmarks.md`, `QUESTIONS.md` written.
- ☑ `tests/test_ssc_integrals_fd.py` skeleton (L1 placeholder, `xfail` until integrals exist).
- ☑ `docs/ssc_zfs_derivation.tex` skeleton + pushed to Overleaf.
- ☑ Reuse map recorded in `CLAUDE.md §5`.

## PHASE 1 — L1: 2e SS integral + FD validation  (gate: §7 L1)
- ☐ **P1.1** Decide Path A (ERI engine w/ derivative kernel) vs Path B (direct Rys, model on
  `comp_soc_int2_prim`). Prototype the rank-2 dipolar kernel `(3 r_k r_l − δ_kl r²)/r⁵` for one
  shell quartet (s,s,s,s) first.
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
- 2026-06-09 — Setup session. Context gathered from `./papers/` (Sinnecker–Neese 2006 eq 9,
  Neese 2007 eq 46 — prefactor discrepancy logged in `CLAUDE.md §3`; Pokhilko–Krylov 2019 W–E
  extraction; Neese JACS 2006 mean-field). Reuse map built. Scaffolding committed. No code yet.
  NEXT: P1.1.

## OPEN ASSUMPTIONS (promote blocking ones to QUESTIONS.md)
- Final `C`/sign deferred to numerical pinning on O₂ (L2) — see `CLAUDE.md §3`. No assumption made yet.
