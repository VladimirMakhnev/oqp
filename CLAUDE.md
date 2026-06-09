# CLAUDE.md — Direct dipolar spin–spin contribution to the ZFS D-tensor for MRSF-TDDFT in OpenQP

**This file is a standing order for a headless (`claude -p`) instance that has NO memory of any
prior conversation.** Read it in full at the start of every run, then read `PROGRESS.md`,
`benchmarks.md`, and `QUESTIONS.md` before doing anything. Work only on branch **`ssc-zfs`**.
All Python/build/run commands go through the `ssc-pyenv` environment (see §Build/Run).

---

## 1. OBJECTIVE & SCOPE

Implement the **first-order direct dipolar spin–spin (SS) contribution to the zero-field
splitting (ZFS) D-tensor** for **MRSF-TDDFT** states in OpenQP.

- This is **first order in perturbation theory only** (the direct electron–electron magnetic
  dipole–dipole term `D^SS`). The second-order spin–orbit (SOC) contributions to D are a
  *separate, already-existing* concern (`soc_mrsf.F90`) and are **out of scope** here. `SSC` is
  the first-order companion of the existing second-order `SOC`.
- The **only physical input** is the one-particle spin-density matrix `P^(α−β)` of the state of
  interest, taken at its **M_S = S** component. Nothing else enters the working equation.
- **Explicitly OUT OF SCOPE (do not implement, do not let scope creep in):** Z-vector /
  orbital-relaxed densities, response/relaxation contributions, analytic gradients of D, SOC
  cross terms. If a task seems to require these, stop and record it in `QUESTIONS.md` — do not
  silently expand scope.

---

## 2. WORKING EQUATION (canonical form = Sinnecker–Neese)

Target form, McWeeny–Mizuno [Proc. R. Soc. London A **259**, 554 (1961)], as written in
**Sinnecker & Neese, J. Phys. Chem. A 2006, 110, 12267, Eq. 9**:

```
D^SS_kl = C · Σ_{μν} Σ_{κτ} { P^(α−β)_μν P^(α−β)_κτ − P^(α−β)_μκ P^(α−β)_ντ }
              · ⟨ μν | (3 r₁₂,k r₁₂,l − δ_kl r₁₂²) / r₁₂⁵ | κτ ⟩
```

- `k, l ∈ {x, y, z}`; the tensor is symmetric and traceless (6 independent components).
- `P^(α−β)` = one-particle spin-density matrix of the **M_S = S** state component, in the AO basis.
- `⟨μν|…|κτ⟩` is a **two-electron** integral over the rank-2 dipolar kernel (the Hessian of
  `1/r₁₂` w.r.t. the interelectronic vector); it has the same `(μν|κτ)` charge-cloud structure
  as an ERI but with a tensor operator instead of `1/r₁₂`.
- The two density products are a **Coulomb-like** term (`P_μν P_κτ`) and an **exchange-like** term
  (`P_μκ P_ντ`) — structurally identical to a Fock/K build (see Reuse Map §5).
- `C` is the prefactor/sign — **DO NOT trust the literature value blindly; see §3.**

---

## 3. CONVENTION-PINNING WARNING  (CRITICAL — read before touching `C`)

The prefactor, the overall sign, and the index grouping **differ between the source papers**, in
ways that do NOT reduce to obvious notation. Verbatim from the PDFs in `./papers/` (project root,
`/bighome/vova/coding/SSC/papers/`):

- **Sinnecker–Neese 2006, Eq. 9:**  `C = + g_e² α² / [4 S(2S−1)]`, kernel `(3 r r − δ r²)/r⁵`,
  density grouping `{P_μν P_κτ − P_μκ P_ντ}`, integral bracket `⟨μν|…|κτ⟩`.
- **Neese, J. Chem. Phys. 127, 164112 (2007), Eq. 46:**  `C = − g_e² α² / [16 S(2S−1)]`, same
  kernel sign `(3 r r − δ r²)`, same density grouping, same bracket `⟨μν|…|κτ⟩`. **This differs
  from Eq. 9 by a factor of −1/4** (sign + magnitude).
- The **same** Neese 2007, **Eq. 1**, writes the operator with the *opposite* kernel sign
  `(r² δ − 3 r r)/r⁵` and spin part `{2 ŝ_z ŝ_z − ŝ_x ŝ_x − ŝ_y ŝ_y}` — i.e. yet another
  sign/þrefactor convention layer once the spin reduction to M_S = S is done.
- The hand-off note that seeded this project also cited a `⟨μκ|…|ντ⟩` grouping for Neese 2007;
  the **printed Eq. 46 uses `⟨μν|…|κτ⟩`**. The discrepancy itself is the point.

**ORDER:** Do **not** derive `C` by algebra and do **not** adopt any single paper's number on
faith. **Pin `C` (magnitude AND sign) NUMERICALLY against the O₂ benchmark** (§6 L2). The
**canonical algebraic form** is Sinnecker–Neese Eq. 9; the **final `C`, sign, and unit factor**
are whatever reproduces the published O₂ `D^SS`. Record the pinned value, its derivation-by-
matching, and the unit convention (a.u. → cm⁻¹) in the LaTeX audit doc (§10).

---

## 4. KEY IMPLEMENTATION FACTS (from the papers)

1. **The 2e SS integral is the Hessian of `1/r₁₂`** w.r.t. the interelectronic coordinate
   (`∂²(1/r₁₂)/∂r_k ∂r_l → (3 r_k r_l − δ_kl r²)/r⁵`). It is therefore **"ERI-like"** and computable
   with the same machinery as ERIs / ERI derivatives. Two routes:
   - **Path A** — drive the existing ERI engine (`int2.F90`) with shifted angular momenta /
     derivative kernels (reuse `grd2_rys.F90`, `rys_deriv.F90`).
   - **Path B** — direct Rys quadrature for the rank-2 kernel, modelled on
     `comp_soc_int2_prim` (`source/integrals/mod_1e_primitives.F90:3649`).
   Decide A vs B during L1; whichever passes the finite-difference integral test (§6 L1) wins.
2. **Six components** of the symmetric traceless tensor are needed: `xx, yy, zz, xy, xz, yz`.
   **Trace = 0** is a free internal numerical check (the kernel `Σ_k (3 r_k r_k − r²) = 0`).
3. **Extract the M_S = S spin density via Wigner–Eckart** (Pokhilko, Epifanovsky & Krylov,
   J. Chem. Phys. 151, 034106 (2019)): the spin density of `|S, M_S⟩` is proportional to the
   single reduced one-particle (transition) density computed for *one* multiplet component
   (highest-spin component as generator; irreducible spin-tensor relations
   `[Ŝ_z, Ô^{S,M}] = M Ô^{S,M}`, `[Ŝ_±, Ô^{S,M}] = √(S(S+1)−M(M±1)) Ô^{S,M±1}`). **Reuse the
   SOC-MRSF reduced/transition densities** (`compute_tdm`, `soc_mrsf.F90:551`) — do not build a
   new density code path.
4. **Methodological risk — reference choice.** UKS/UHF **systematically overestimates** `D^SS`
   (spin polarization in the unrestricted natural orbitals), even at near-pure `⟨S²⟩`;
   **ROKS / UNO ≈ accurate** (Sinnecker–Neese; Neese 2007 §E: ROKS and UNO determinants give
   virtually identical `D^SS`). MRSF uses an **ROHF reference → favourable**. For the **UMRSF**
   branch, build in a **density-source flag** (default the RO/UNO-type density; expose UNO as an
   option) so the unrestricted overestimation can be avoided/diagnosed.

---

## 5. REUSE MAP (concrete entry points found in the tree — verify line numbers before editing)

| Need | Reuse | Location |
|------|-------|----------|
| ERI engine + "contract ERIs with a density" consumer pattern | `int2_compute_data_t` / `int2_fock_data_t` (abstract), `int2_rhf_data_t`, `int2_urohf_data_t`; driver `int2_run`/`int2_twoei` | `source/integrals/int2.F90` (`:96`, `:123`, `:129`, `:443`, `:532`) |
| MRSF density-contracted Fock-like build (closest template for the SS Coulomb/exchange contraction) | `int2_mrsf_data_t`, `int2_umrsf_data_t` (extend `int2_fock_data_t`) | `source/tdhf_mrsf_lib.F90` (`:8`, `:24`) |
| 2e SOC integral primitive (Path B template, Rys) | `comp_soc_int2_prim` (+ calling-convention doc at `:3378`) | `source/integrals/mod_1e_primitives.F90:3649` |
| Rys / Obara–Saika / HRR, incl. derivative kernels (Path A) | `int_rys.F90`, `rys.F90`, `rys_lut.F90`, `rys_deriv.F90`, `grd2_rys.F90` | `source/integrals/` |
| Reduced / transition 1-particle densities of MRSF states (M_S extraction, §4.3) | `compute_tdm`; `compute_soc_ao`, `ao2mo_soc` | `source/modules/soc_mrsf.F90` (`:551`, `:336`, `:490`) |
| MRSF state densities (boundary: relaxed density = OUT OF SCOPE, do not pull in) | `build_mrsf_relaxed_density_and_w` (Z-vector — OUT OF SCOPE marker), `grd2_mrsf_..._get_density` | `tdhf_mrsf_z_vector.F90:1455`, `tdhf_mrsf_gradient.F90:300` |
| Property-module style template (C-binding + Fortran driver + 6-component AO tensor + print) | `electric_moments` (`*_C(c_handle) bind(C, name=…)` + `subroutine electric_moments(infos)`); `compute_soc_ao`/`print_soc_ao_gamess` for the tensor pattern | `source/modules/electric_moments.F90`, `source/modules/soc_mrsf.F90:336/436` |
| Python dispatch wiring (mirror `soc`) | `runtype 'soc' → compute_soc` (`pyoqp.py:155`); `compute_soc → oqp.soc_mrsf(mol)` (`runfunc.py:119`); runtype list + `_check_soc` (`input_checker.py:15`, `:1486`); cffi cdef (`pyoqp/oqp/__init__.py`, `openqp.py`) | `pyoqp/oqp/…` |

**SSC = the first-order companion of SOC.** Wherever SOC has a piece (AO integrals → AO2MO →
state densities → assemble matrix → diagonalise/print → Python dispatch), SSC has an analogue;
prefer extending/copying these patterns over inventing new ones.

---

## 6. BENCHMARK LADDER (numerical oracle; numbers extracted from `./papers/`)

`benchmarks.md` holds the exact numbers/geometries/tolerances and is **READ-ONLY ground truth**.
Each rung isolates one error source.

- **L1 — Integral unit test (fully internal).** The 2e SS integral = Hessian of `1/r₁₂`. Verify
  each tensor component by **finite differences of the existing ERI engine** (perturb the kernel
  / nuclear-free interelectronic coordinate as appropriate) to **~6–8 significant figures**.
  Also check **trace = 0** to machine precision. No external number needed.
- **L2 — Single-determinant ROHF, where McWeeny–Mizuno is EXACT.** This is where **`C`/sign is
  pinned**. Primary pin: **O₂ ³Σ_g⁻ at r = 1.207 Å, `D^SS ≈ 1.44–1.6 cm⁻¹`** (1.44 = Vahtras
  CASSCF reference; 1.52–1.59 across BP/B3LYP/CASSCF × EPR-II/EPR-III/QZVP in Sinnecker–Neese).
  NB: the *total* O₂ `D ≈ 3.96–4.0 cm⁻¹` includes SOC — **not** our target; the **SS-only** part
  is 1.44–1.6. Secondary: **CH₂ ³B₁** (Petrenko–Neese applied MM to a single KS determinant) —
  confirm its published `D^SS` in `benchmarks.md` before using.
- **L3 — MRSF triplets.** Aromatic-triplet `D` (cm⁻¹), experimental targets: **benzene 0.1593,
  naphthalene 0.1004, anthracene 0.0702, tetracene 0.0573**; plus organic radicals/carbenes from
  Sinnecker–Neese. Target **RMSD ≈ 0.0035 cm⁻¹** for the RO-type set. **Caveat:** plain RO-DFT
  underestimates the larger polyacenes by ≈ factor 2 (benzene OK; tetracene ~0.031 calc vs 0.057
  exp) — judge MRSF against the appropriate column, and record which reference column is used.
- **L4 — Tensor invariants (free).** `trace(D) = 0`; `E/D ∈ [0, 1/3]`. Run at every level.

---

## 7. STAGE GATES (hard ordering)

```
L1  →  L2  →  L3
```

- **Do NOT** move to MRSF densities (L3) until the **ROHF level (L2) is reproduced**.
- **Do NOT** start L2 until the **integral FD test (L1) passes**.
- A gate is passed **only** when the number matches within tolerance:
  - **L1:** agreement with the ERI-engine finite difference to **machine precision relative to
    the FD step** (≈ 6–8 sig. figs.); trace = 0 to ~1e-10.
  - **L2:** reproduce the **published method value** for O₂ `D^SS` to **a few %**; sign correct.
  - **L3:** **RMSD vs the benchmark table** within the stated tolerance (≈ 0.0035 cm⁻¹, RO-type).
- A passed gate must be logged in `PROGRESS.md` with the actual number vs target.

---

## 8. AUTONOMOUS-MODE OPERATING RULES (`claude -p`, no conversational memory)

At the **start of every run**: read `CLAUDE.md`, `PROGRESS.md`, `benchmarks.md` (and `QUESTIONS.md`).

Then:
1. Pick the **next unblocked task** from `PROGRESS.md` (respect the stage gates §7).
2. Do it; **run its gate/tests**.
3. **Update `PROGRESS.md`**: what was done, current numbers vs target, what is next.
4. Make a **small commit** with a meaningful message and **push to branch `ssc-zfs`**.

Hard rules:
- **NEVER declare a gate passed without the number matching.** If a gate resists after a few
  honest attempts, write a root-cause analysis in `PROGRESS.md` and either switch to another
  unblocked subtask or record a question (§9).
- **NEVER edit `benchmarks.md` / reference values to make a test pass.** The oracle is read-only.
  Fix the code, not the target.
- **NEVER expand scope** (no Z-vector, no gradients of D — §1).
- **Be honest about uncertainty.** Do not present a plausible-but-unverified result as success.
  Distinguish "computed and matched reference" from "computed, not yet validated."

---

## 9. QUESTION / STOP PROTOCOL

Stop and ask a human **only** for **truly blocking** situations: the build is broken and cannot be
repaired, or a physical/numerical ambiguity is genuinely irreducible (e.g. two pinned `C` values
each match a different published convention and the experiment cannot discriminate).

- **Batch** questions into `QUESTIONS.md` (use its header format). Do **not** interrupt the work
  stream for every uncertainty.
- For everything else: **proceed on the best-justified assumption**, and document that assumption
  explicitly in both `PROGRESS.md` and the LaTeX audit doc (§10).

---

## 10. REPORTING / OUT-OF-BAND REVIEW (the human does NOT check the physics)

Because no human validates the physics inline, leave an auditable paper trail:

1. **LaTeX derivation/audit document** (`docs/ssc_zfs_derivation.tex`, mirrored to the Overleaf
   repo at `/bighome/vova/coding/SSC/overleaf`, pushed to `git.overleaf.com`): re-derive the
   working equation; **explicitly resolve the convention discrepancy (§3)**; record the final
   pinned `C`, sign, and unit factor and *how* they were pinned (the O₂ match); log every
   assumption. Written for an **external reviewer** to audit.
2. **RESULTS table** (value vs reference vs tolerance, pass/fail) in the repo and/or Overleaf,
   updated as gates are cleared.
3. **Regression tests** in `tests/` encoding the published numbers (L1 FD test first, then L2
   O₂, then L3 acenes).

---

## Build / Run (this server, env `ssc-pyenv`)

`ssc-pyenv` is a venv on the KNU cluster (NOT conda). Build/test on an **AVX2 node** (some nodes
are Sandy Bridge and SIGILL the module Python). The SSC branch requires **ILP64 BLAS** and the
test runner requires **Python ≥ 3.11**:

```bash
module purge
module load Python/3.11.3-GCCcore-12.3.0 GCC/12.3.0 OpenMPI/4.1.5-GCC-12.3.0 \
            CMake/3.26.3-GCCcore-12.3.0 imkl/2023.1.0
source /bighome/vova/coding/SSC/ssc-pyenv/bin/activate      # create once: python3 -m venv ssc-pyenv
pip install -U pip cffi
cmake -B build -G Ninja -DUSE_LIBINT=OFF \
  -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DCMAKE_Fortran_COMPILER=gfortran \
  -DCMAKE_INSTALL_PREFIX=. -DENABLE_OPENMP=ON \
  -DLINALG_LIB_INT64=ON -DLINALG_LIB=Intel10_64ilp
ninja -C build install
pip install ./pyoqp/.
export OPENQP_ROOT=$(pwd)
openqp some_input.inp           # run a calc
openqp --run_tests all          # full suite (keep imkl module loaded at runtime)
```

Baseline on this branch's parent: 222/225 tests pass (3 pre-existing SOC/MRSF failures, unrelated
to SSC).
