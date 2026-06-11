# MRSF-NMR implementation roadmap (OpenQP, branch `nmr-mrsf`)

Status: planning document, 2026-06-11.
Author: grounding session based on A. Lashkaripour's plan (Overleaf perspective
draft `main.tex` + `supporting-information.tex`, commit `15ac180`) and direct
inspection of this checkout (HEAD `0b4a91c`).

Target property (plan Eq. `state-shielding`):

    sigma^K_ab(I) = d^2 E_K / dB_a dm_{I,b}  at B=0, m=0

for a selected MRSF-TDDFT state K, with GIAO (London orbital) gauge treatment.
The plan's central proposal (plan Eq. `cp-mrsf-magnetic`) is a magnetic
analogue of the MRSF Z-vector problem,

    A_MRSF z^(B_a) = -b^(B_a),  a = x,y,z   (3 solves, nucleus-independent)

followed by per-nucleus PSO + diamagnetic contractions with the state-relaxed
density.

---

## 0. Verified code inventory (what exists at HEAD `0b4a91c`)

All file/line references verified in this checkout.

### Ground-state NMR machinery (reusable)

| Item | Location | Notes |
|---|---|---|
| CGO shielding driver | `source/modules/nmr_shielding.F90` (492 L) | closed-shell only; dia + uncoupled SOS para + coupled (exchange-only) CPHF; stores `OQP::nmr_shielding` (5 iso values/atom) |
| GIAO shielding driver | `source/modules/nmr_giao_shielding.F90` (688 L) | RHF/UHF/ROHF, HF + LDA/GGA/global hybrids; entry `nmr_giao_shielding_debug` (production binding `nmr_giao_shielding` aliases it) |
| GIAO 2e Fock derivative | `source/modules/nmr_giao_debug.F90`, `giao_h10_twoe_matrix` | Rys-only (`int2_driver%rys_only = .true.`) |
| GIAO XC derivative | `source/dftlib/dft_gridint_giao.F90`, `giao_vxc` | LDA/GGA only (no tau) |
| Magnetic 1e integrals | `source/integrals/int1.F90`: `angular_momentum_integrals` (l.309), `giao_overlap_derivative` (l.364), `giao_h10_core` (l.416), `nmr_dia_shielding` (l.474), `giao_a11part_corr` (l.546), `giao_a01gp_contract` (l.614), `pso_integrals` (l.676) | the full GIAO/PSO/dia integral set already exists; SOC L-integrals are NOT needed for NMR (NMR has its own `angular_momentum_integrals`) |
| CPHF magnetic solvers | `nmr_giao_shielding.F90`: `solve_mo1_uncoupled` (l.462), `solve_mo1_coupled` (l.488), `giao_pb_density` (l.572), `para_tensor` (l.658), `giao_para_channel` (l.589), `semicanon_orbitals` (l.627) | fixed-point CPHF; exchange-only kernel for the imaginary antisymmetric first-order density |
| Calibrated constants | `nmr_giao_shielding.F90` l.106-107: `SG=-1, SC=1, SA=-0.5, SX=-1` | "fixed vs the oracle"; must be re-derived/documented (plan Phase 1) |
| Input/dispatch | `pyoqp/oqp/library/runfunc.py::compute_scf_prop` (l.41-75); schema `pyoqp/oqp/molecule/oqpdata.py` l.190 (`nmr_gauge`, default `cgo`), l.191 (`td_prop`, bool); guards `pyoqp/oqp/utils/input_checker.py` l.923-986 | `scf_prop=nmr`, `nmr_gauge=cgo/giao`; CGO+open-shell, CAM, meta-GGA, ECP rejected |
| Output | tag `OQP::nmr_shielding` (`source/tagarray_driver.F90` l.37); JSON save `pyoqp/oqp/molecule/molecule.py` l.1240 (`reshape(-1,5)`) | isotropic-only, no tensor, no state label |
| Tests/fixtures | `tests/test_nmr_shielding.py`, `test_nmr_coupled.py` (gates 0-6), `test_nmr_giao_para_live.py`, `test_nmr_gauge_interface.py`; `tests/fixtures/nmr/{cgo,giao}_reference.json` | ground-state regression baseline |
| Examples | `examples/NMR/*.inp` (H2O RHF/PBE0; CH3 UHF/ROHF, HF and PBE0) | |

### MRSF response machinery (reusable)

| Item | Location | Notes |
|---|---|---|
| MRSF energy/amplitudes | `source/modules/tdhf_mrsf_energy.F90` | stores `OQP::td_bvec_mo` (amplitudes, all states), `OQP::td_t`, `OQP::td_energies`; requires ROHF triplet reference (`mult=3`) |
| Z-vector solver | `source/modules/tdhf_mrsf_z_vector.F90` (1740 L) | `gmres_solve` (l.103), `apply_z_operator` (l.654), CG/MINRES/GMRES/AUTO (l.1089); RHS built in `build_mrsf_zvector_rhs` (l.1556) from `mrsfcbc`/`mrsfsp`/`sfrorhs` |
| Z-vector space | `lzdim = noccb*(nsocc+nvira) + nsocc*nvira` (l.938) | ROHF orbital-rotation space, `xk` is the multiplier vector |
| Relaxed density | `build_mrsf_relaxed_density_and_w` (l.1455): `OQP::td_p` (relaxed difference density, alpha/beta, packed), `OQP::wao` (energy-weighted density W), `OQP::td_mrsf_density` (7 component densities), `OQP::td_abxc` | exactly the objects the plan wants to feed into the shielding contractions |
| Gradient contraction template | `source/modules/tdhf_mrsf_gradient.F90` | contracts td_p/wao/td_mrsf_density with derivative integrals — the pattern Gate 5 mimics for S10 |
| Excited-state property template | `source/modules/electric_moments.F90::electric_moments_excited` | relaxed state density = `DM_A + DM_B + td_p(:,1) + td_p(:,2)`; dispatched from `Gradient.gradient()` via `properties.td_prop` (`pyoqp/oqp/library/single_point.py` l.844-846) — **this is the cleanest integration-point precedent for MRSF-NMR** |
| Python<->Fortran ABI | `include/oqp.h` + cffi auto-wrap in `pyoqp/oqp/__init__.py` (l.54-58) | a new `bind(C)` subroutine + one line in `oqp.h` = callable as `oqp.<name>(mol)` |

### Critical structural observation (drives Gates 3-4)

The gradient Z-vector operator `apply_z_operator` is built for **real,
symmetric** trial densities: `sfrogen` -> AO transform -> `int2` with
`int_apb=.true.` (Coulomb + exchange) -> `symmetrize_matrix` -> `utddft_fxc`
(semilocal kernel) -> `sfrolhs`.

A magnetic perturbation is **imaginary/antisymmetric**. For such densities the
Coulomb and semilocal-XC kernel images vanish identically; the surviving
two-electron response is exact-exchange-only and antisymmetric (this is
precisely what `nmr_shielding.F90::compute_coupled_para` and
`nmr_giao_shielding.F90::solve_mo1_coupled` implement for the SCF case, with
`int_amb=.true.` instead of `int_apb`).

**Therefore plan Eq. `cp-mrsf-magnetic` cannot be implemented by literally
reusing `apply_z_operator`.** The magnetic MRSF response operator is a
different (A-B)-type, exchange-only sibling of the existing operator. What IS
reusable: the solver framework (`gmres_solve`/MINRES/CG), the orbital-rotation
vector layout (`sfrogen`/`sfrolhs` index conventions), the `int2` driver, and
the spin bookkeeping. This is the single largest deviation between the plan as
written and the code as it exists — see Open Questions Q1.

---

## Gate 0 — Reproducible ground-state baseline (plan Phase 1)

**Goal / physics.** Build this checkout, reproduce the committed CGO/GIAO
regression results and the OpenQP-vs-PySCF benchmark level of agreement
(plan Eq. `isotropic-shielding-benchmark`, Tables `tab:openqp-pyscf-nmr*`:
HF MAE 7.2e-6 ppm, PBE0 MAE 8.4e-4 ppm). Establish the "if MRSF-NMR fails
later, is it the magnetic kernels or the MRSF response?" baseline.

**Work.**
1. Build OQP on this cluster (use the established `oqp-setup` workflow:
   modules, venv, cmake/ninja, `OPENQP_ROOT`).
2. Run `tests/test_nmr_*.py` against fixtures; run all `examples/NMR/*.inp`;
   archive logs.
3. Run one CGO case with two different gauge origins and one GIAO case
   re-centered, to confirm GIAO origin-independence behavior.
4. Document (do not yet change) the calibrated constants `SG/SC/SA/SX` and the
   `pso_integrals` sign convention (`para_tensor` comment, l.679-682).

**Depends on:** nothing.
**Validation:** committed fixtures `tests/fixtures/nmr/*.json`; PySCF
cross-check optionally re-run for one molecule.
**Effort/risk:** low (days). Pure infrastructure.

## Gate 1 — Derivation: the MRSF magnetic Lagrangian (plan Phase 2)

**Goal / physics.** Produce the working equations for plan Eq.
`state-shielding` as a Lagrangian derivative. This is the *novel* content; no
published MRSF-NMR precedent exists. Deliverable: `derivation_nmr.tex` in the
Overleaf project, plus an updated code-to-theory map (extends SI Table
`tab:si-term-map`).

**Required content (minimum):**
1. The MRSF Lagrangian
   `L_K = E_K[X_K, C] + z^T (orbital stationarity) + w (normalization, S-dependence)`
   with GIAO (B-dependent) AOs, written for the ROHF triplet reference and
   `mrst=1/3` response states.
2. Mixed derivative wrt (B_a, m_Ib). Since m_I enters only the one-electron
   PSO/dia operators (linearly, no GIAO dependence on m), the derivative
   splits into:
   - diamagnetic term: relaxed state density x d2h/dB dm (existing kernels
     `nmr_dia_shielding`, `giao_a11part_corr`, `giao_a01gp_contract`);
   - paramagnetic term: (d/dB of the relaxed state density) x PSO
     (existing kernel `pso_integrals` + `para_tensor`-style contraction).
3. The B-response system. Classify exactly which first-order quantities are
   needed: orbital response U^(B), amplitude response X^(B), multiplier
   response z^(B) — or show which of them can be eliminated by the
   (2n+1)/(2n+2) rules given that L_K is stationary in X_K and C.
   **This must settle whether one Z-vector-like solve per field component
   (plan Eq. `cp-mrsf-magnetic`) is sufficient, or whether an amplitude
   response solve `(Omega - omega_K) X^(B) = -RHS` is additionally required.**
   (See Q2.)
4. The symmetry analysis: imaginary antisymmetric perturbed densities; which
   kernel pieces survive (expected: exact exchange only, as in the SCF case);
   the GIAO S10 terms entering through both h-bar^(B) and the energy-weighted
   density (mirror of `wao`*S^x in the gradient).
5. Spin structure: the orbital-Zeeman operator is spin-conserving (MS-diagonal)
   while MRSF amplitudes live in the spin-flip space; derive how the magnetic
   perturbation enters the MS=+1/-1 mixed-reference blocks (see Q4).

**Depends on:** Gate 0 only for context.
**Validation:** internal consistency checks: (a) zero-amplitude limit must
collapse to the ground-state GIAO CPHF equations as implemented; (b) sum rules
/ antisymmetry of every intermediate; (c) the derivation must reduce to the
known SF-TDDFT gradient Lagrangian structure when B -> geometric perturbation.
**Effort/risk:** high; this is the scientific bottleneck. Weeks. Everything
in Gates 3-5 is blocked on its outcome.

## Gate 2 — Standalone prototype A: state-density shielding (no new response)

**Goal / physics.** The smallest runnable artifact: contract the **existing**
relaxed MRSF state density with the **existing** shielding kernels, i.e.
compute the diamagnetic tensor of state K exactly, and the paramagnetic term
in a frozen-response approximation (ground-state-like CPHF with the same
orbitals). This is *not* MRSF-NMR (the para term is not the true state
response — must be labeled as approximation), but it:
- exercises every integration point end-to-end,
- provides the dia term that survives unchanged into the final method,
- gives the zero-amplitude regression target.

**Work.**
1. New module `source/modules/nmr_mrsf_shielding.F90`:
   - `subroutine nmr_mrsf_shielding_C(c_handle) bind(C, name="nmr_mrsf_shielding")`
   - `subroutine nmr_mrsf_shielding(infos)`:
     reads `OQP_DM_A/B` + `OQP_td_p` (pattern: `electric_moments.F90`
     l.61-73: `dens_ex = dmat_a + dmat_b + td_p(:,1) + td_p(:,2)`),
     builds the GIAO diamagnetic tensor by calling `nmr_dia_shielding`,
     `giao_a11part_corr`, `giao_a01gp_contract` with `dens_ex`,
     and (clearly-labeled) the approximate para term by reusing the
     `giao_para_channel` pipeline with semicanonicalized ROHF orbitals.
   - hard guards: abort unless `tdhf.type=mrsf`, `scf.type=rohf`, `mult=3`,
     no CAM / meta-GGA / ECP (copy the guard block from
     `nmr_giao_shielding.F90` l.180-203); abort if `OQP_td_p` missing
     (i.e. Z-vector not run) — **fail fast, no silent ground-state fallback**
     (plan Phase 3 requirement).
2. Declare in `include/oqp.h`: `void nmr_mrsf_shielding(struct oqp_handle_t *inf);`
3. Add to `source/modules/CMakeLists.txt` (or wherever module sources are
   listed — check sibling entries for `nmr_giao_shielding.F90`).
4. Python dispatch: in `pyoqp/oqp/library/single_point.py::Gradient.gradient()`
   after the Z-vector solve, next to the existing `td_prop` block (l.844-846):
   gate on a new schema key `properties.nmr_mrsf` (bool, default False) in
   `oqpdata.py`, mirrored by an `input_checker.py` rule that *requires*
   `tdhf.type=mrsf` + `runtype=grad` (or a dedicated runtype, see Q7) and
   forbids `scf_prop=nmr` in the same run.
5. New tag `OQP::nmr_shielding_mrsf` in `tagarray_driver.F90` storing, per
   atom: full 3x3 dia tensor, 3x3 para tensor, iso values, plus the state
   index and gauge as metadata; JSON save in `molecule.py` next to l.1240.
   (Per plan SI: state labels + full tensors + provenance, not 5 iso scalars.)
6. Test: `tests/test_nmr_mrsf_prototype.py` — H2O/BHHLYP/6-31G* MRSF
   (reuse `examples/MRSF-TDDFT/H2O_BHHLYP-MRSFTDDFT_GRADIENT.inp` settings,
   add the NMR flag), assert (a) it runs, (b) dia tensor is symmetric-ish and
   plausible, (c) guards fire for wrong references.

**Depends on:** Gate 0 (built code). NOT blocked on Gate 1 (the dia
contraction is already exact; the para label "approximate" is honest).
**Validation:** zero-amplitude/closed-shell limit: run a fake "MRSF with
negligible response" case and compare the dia term against the ground-state
GIAO dia output; no independent oracle exists for the approximate para term
(state explicitly: this gate validates plumbing, not physics).
**Effort/risk:** medium-low (1-2 weeks). All ingredients exist.

## Gate 3 — Magnetic MRSF response operator

**Goal / physics.** Implement the (A-B)-type, exchange-only, antisymmetric
counterpart of `apply_z_operator`, per the Gate 1 derivation: the operator
that acts on imaginary orbital-rotation vectors of the ROHF reference.

**Work.**
1. In `tdhf_mrsf_z_vector.F90` (or a new `tdhf_mrsf_magnetic_response.F90`
   module to keep the gradient path untouched):
   `subroutine apply_z_operator_magnetic(x_in, x_out, infos, basis, molGrid, int2_driver, nocca, noccb, nbf, mo_a, mo_b, mo_energy_a, fa, fb, scale_exch, dft)`
   — same signature as `apply_z_operator` (l.654) so `gmres_solve`'s
   `apply_operator` interface (l.116) accepts it. Differences inside:
   - antisymmetric density generation (an `expand`-style antisymmetric variant
     of the `sfrogen` -> `orthogonal_transform` step; do NOT call
     `symmetrize_matrix`);
   - `int2_td_data_t(..., int_apb=.false., int_amb=.true., scale_exchange=scale_exch)`
     — exchange-only, as in `solve_mo1_coupled` (l.526-527);
   - NO `utddft_fxc` call (semilocal kernel vanishes; assert this against a
     finite-difference check in the unit test rather than silently assuming —
     cheap to verify once);
   - diagonal part via `sfrolhs`-style orbital-energy differences (the
     `fa/fb` ROHF Fock blocks logic carries over).
2. Unit tests (live, like `test_nmr_giao_para_live.py` pattern):
   - operator linearity + expected (anti)symmetry on random vectors;
   - **zero-amplitude limit:** for a closed-shell-like ROHF case, the solve
     with the ground-state GIAO RHS must reproduce `solve_mo1_coupled`'s
     occ-vir response (this cross-validates against the *existing, validated*
     SCF magnetic CPHF).

**Depends on:** Gate 1 (operator definition), Gate 0.
**Validation:** as above; plus pure-functional case must reduce to the
uncoupled (orbital-energy-denominator) solve.
**Effort/risk:** medium-high. The ROHF three-block (closed/open/virtual)
rotation space with an antisymmetric perturbation is delicate; sign errors
here are the classic failure mode. 2-4 weeks.

## Gate 4 — Magnetic right-hand side b^(B) and (if required) amplitude response

**Goal / physics.** Assemble `b^(B_a)` of plan Eq. `cp-mrsf-magnetic` from the
GIAO pieces (`giao_h10_core` + `giao_h10_twoe_matrix` + `giao_vxc` + 
`giao_overlap_derivative`) contracted with the MRSF state quantities
(`td_p`, `wao`, `td_mrsf_density`, amplitudes `td_bvec_mo`), exactly as
prescribed by the Gate 1 derivation. If Gate 1 concludes an amplitude-response
solve is required, implement it here (the Davidson machinery in
`tdhf_mrsf_energy.F90` provides the Omega-matrix application to reuse for a
linear solve `(Omega - omega_K) X^(B) = rhs`).

**Work.**
1. `subroutine build_mrsf_magnetic_rhs(infos, basis, ..., b(:,3))` in the new
   module; term-by-term, each term implemented behind its own flag so it can
   be switched off in tests (mirrors how the gradient RHS is staged in
   `build_mrsf_zvector_rhs`, l.1556).
2. Per-term finite-difference oracles where possible (see Q8 for what is
   finite-differentiable): e.g. GIAO S10/h10 contractions can be checked
   against numerical d/dB of overlap/hcore matrices at complex-free level
   (the integrals themselves are real coefficient matrices of an imaginary
   perturbation — they ARE finite-differentiable through libcint/PySCF
   cross-checks, as was already done for the ground state).
3. Solve the 3 systems with `gmres_solve(apply_z_operator_magnetic, ...)`.

**Depends on:** Gates 1, 3.
**Validation:** zero-amplitude limit of the full RHS must equal the
ground-state GIAO RHS (`h1mo - s1mo*e_i` structure of `giao_para_channel`);
term-wise FD checks.
**Effort/risk:** high. This is the largest new-code surface. 4-8 weeks.

## Gate 5 — Assembly, output, integration (plan Phase 3)

**Goal / physics.** Full pipeline: MRSF energy -> Z-vector (gradient one, for
td_p) -> magnetic response solves -> PSO/dia contraction -> state-labeled
tensors. Scope per plan Phase 3: nonrelativistic, no ECP, no CAM/meta-GGA,
global hybrids (BHHLYP is the MRSF workhorse — note it IS a global hybrid, so
the exchange-only coupled response machinery applies directly), `mrst=1/3`
states only (quintets `mrst=5` deferred, Q9).

**Work.**
1. Replace the Gate 2 approximate para path with the true response:
   sigma_para^K = contraction of the B-responded relaxed density with
   `pso_integrals` per nucleus (reuse `para_tensor` shape).
2. Add the S10 x W^(B)-type and renormalization terms per Gate 1.
3. Output: extend the Gate 2 tag to {dia tensor, para tensor, total tensor,
   iso, state index K, gauge, provenance string}; human-readable table in the
   log mirroring the GIAO table (l.380-401).
4. Examples: `examples/NMR-MRSF/H2O_BHHLYP-MRSF-NMR.inp` + json.
5. Guards inventory test: every unsupported combination aborts with a
   specific message (`test_nmr_mrsf_guards.py`).

**Depends on:** Gates 2, 3, 4.
**Effort/risk:** medium. 2-3 weeks once Gate 4 lands.

## Gate 6 — Validation against independent oracles (plan Phase 4)

**Goal.** Layered validation per plan:
1. Ground-state/zero-amplitude limit -> existing GIAO results + PySCF tables
   (plan Tables `tab:openqp-pyscf-nmr`, `tab:openqp-pyscf-nmr-pbe0`).
2. Gauge-origin independence: rigid translation of the molecule must leave
   GIAO MRSF shieldings invariant to numerical tolerance.
3. Independent excited-state oracle: **the plan does not name one** (see Q8).
   Candidate strategies, in order of defensibility:
   a. finite-field cross-derivative on a *modified* axis: differentiate the
      analytic MRSF para response wrt amplitude/density perturbations
      term-wise (component oracles rather than end-to-end);
   b. cross-code excited-state NMR where any exists (e.g. CASSCF GIAO
      shieldings in Dalton for a 2-state model system; published
      excited-state shielding studies);
   c. internal consistency: S-T gap systems where the MRSF "excited" state is
      actually the closed-shell-like ground configuration of another
      reference, comparable to a ground-state calculation.
   End-to-end finite-field validation of sigma is NOT directly available in a
   real-arithmetic code (both perturbations are imaginary); this must be
   stated honestly in the paper.
4. Stress cases (plan Phase 4): near-degeneracy / spin-coupled systems.

**Depends on:** Gate 5.
**Effort/risk:** medium-high; the oracle question is partly scientific, not
just engineering.

## Gate 7 — Method paper (plan Phase 5)

Convert `derivation_nmr.tex` + benchmarks into `article_nmr.tex`. Out of
scope for coding; tracked in the Overleaf project.

---

## Open questions / gaps (surfaced, NOT silently resolved)

**Q1 — "Same A_MRSF operator" (plan Eq. `cp-mrsf-magnetic`) vs. code reality.**
The gradient Z-vector operator (`apply_z_operator`, l.654) acts on real
symmetric densities with a Coulomb + semilocal-fxc + exchange kernel
(`int_apb=.true.` + `utddft_fxc`). An imaginary magnetic perturbation requires
an antisymmetric, exchange-only operator (as the SCF-level code itself
implements in `solve_mo1_coupled`). The plan's claim that the *same* operator
is reused is not implementable as written; only the solver framework and
index conventions carry over. Gate 1 must define the correct operator;
Gate 3 implements it. *Decision needed from the theory side (Alireza):
confirm the intended operator structure.*

**Q2 — Is one Z-vector-like solve per field component sufficient?**
E_K depends on amplitudes X_K, orbitals, and the reference. For the mixed
derivative d2E/dB dm with m entering only linearly through one-electron
operators, a Handy-Schaefer-type interchange can possibly reduce the work to
3 orbital-space solves — but the amplitude response X^(B) does not obviously
disappear: the standard excited-state property route (e.g. TDDFT excited-state
NMR by analogy with excited-state gradients) needs d(X)/dB unless the
Lagrangian eliminates it via stationarity of the eigenproblem. The plan never
mentions amplitude response. Gate 1 must resolve this; if X^(B) is needed,
Gate 4 grows by an `(Omega - omega_K)` linear solve per field component
(3 more solves, possibly with singularities near state crossings).

**Q3 — GIAO dependence of the MRSF eigenproblem itself.**
With London orbitals the AO basis depends on B, so the MRSF A/B matrices and
the amplitude normalization condition acquire S10-type derivative terms. The
plan lists "GIAO derivative terms" only inside b^(B) generically. The
derivation must make explicit which GIAO terms hit the amplitude space vs the
orbital space (the gradient code's `wao`*S^x pattern is the template for the
orbital part).

**Q4 — Spin structure of the magnetic perturbation in the spin-flip manifold.**
The orbital-Zeeman/PSO operators are spin-conserving; MRSF amplitudes connect
MS=+-1 components of the mixed reference. How the magnetic response distributes
over the MRSF configuration blocks (the 7-density structure of
`OQP::td_mrsf_density`, the spin-pair couplings `spc_*` in
`build_mrsf_zvector_rhs` l.1658-1663) has no precedent. This is genuinely
novel theory; options (derive in the full MRSF ansatz vs. first derive in
plain SF-TDDFT as a stepping stone) should be decided explicitly. A plain
SF-NMR intermediate (using `tdhf_sf_z_vector.F90`) would be a publishable and
lower-risk stepping stone — plan does not discuss it.

**Q5 — Which density enters the diamagnetic term.**
Plan says "feed the relaxed MRSF density into the existing density-driven
contractions". Plausible, but the GIAO diamagnetic pieces (`a11part`, `a01gp`)
came out of a *ground-state* derivative; for a state Lagrangian there may be
additional renormalization/W-type diamagnetic terms (S10-driven). Gate 1 item.

**Q6 — Validated functional scope.**
MRSF production use is BHHLYP (global hybrid) — compatible with the
exchange-only magnetic kernel. But `tddft%HFscale` vs `dft%HFscale`
(reference vs response exchange scaling, see `build_mrsf_zvector_rhs`
l.1575-1580) raises: which c_x scales the *magnetic* exchange response?
The SCF NMR code uses `dft%HFscale`. Plan is silent. Theory decision.

**Q7 — User-facing semantics.**
Plan says only "an excited-state or MRSF-specific selector". Options:
(a) `properties.nmr_mrsf=true` piggy-backing on `runtype=grad` (minimal,
follows `td_prop` precedent, but forces a gradient run); (b) extend
`scf_prop=nmr` + `nmr_gauge` with `nmr_state=K` (cleaner for users, more
plumbing: needs Z-vector orchestration outside the gradient path);
(c) a dedicated `runtype=nmr`. Roadmap assumes (a) for the prototype and
defers the final choice. *Vladimir's call.*

**Q8 — No named numerical oracle for the excited-state result.**
Plan Phase 4 says "independent excited-state magnetic-shielding references"
without specifying any. Both perturbations are imaginary, so end-to-end
finite-field checks need complex arithmetic that OQP does not have. Term-wise
component oracles + cross-code comparisons + limits (Gate 6) are the realistic
package; the residual untested core (the MRSF-specific response coupling) must
be acknowledged. If a complex-capable reference implementation (e.g. a PySCF
toy MRSF-NMR script) is feasible, it would close this hole — large side
project, flag for discussion.

**Q9 — Quintet (`mrst=5`) and `umrsf` paths.**
The Z-vector machinery branches for quintets and a UMRSF variant exists
(`tdhf_umrsf_energy`). Prototype scope should exclude both explicitly
(guards), per plan Phase 3 minimalism — but say so in the paper.

**Q10 — GIAO entry point hygiene (plan Phase 1 ask).**
`nmr_giao_shielding.F90`'s only implementation is still named/documented
`nmr_giao_shielding_debug` with calibrated constants `SG/SC/SA/SX` "fixed vs
the oracle" and a TODO comment at l.317 that contradicts the implemented
`a01gp` term (the comment says "a01gp gauge-correction: TODO" while the term
IS included at l.332). Before building MRSF-NMR on top: rename/promote the
production entry, derive the four constants on paper, fix the stale comment,
and commit that derivation to the repo. Low effort, high trust payoff.

**Q11 — Output back-compat.**
`OQP::nmr_shielding` stores 5 isotropic scalars/atom; the MRSF result needs
full tensors + state labels. Decision: new tag (assumed here) vs extending the
old layout (breaks `molecule.py` reshape and downstream parsers).

---

## Dependency graph

```
Gate 0 (build/baseline)
  ├── Gate 1 (derivation)  ──────────────┐
  └── Gate 2 (prototype A: dia + plumbing)│
            │                             │
            ▼                             ▼
        Gate 3 (magnetic operator) ◄── needs operator definition
            │
            ▼
        Gate 4 (magnetic RHS [+ amplitude response?])
            │
            ▼
        Gate 5 (assembly/integration)
            │
            ▼
        Gate 6 (validation) ──► Gate 7 (paper)
```

Gates 0 and 2 can start immediately and in parallel with Gate 1.
