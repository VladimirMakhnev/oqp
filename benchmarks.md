# benchmarks.md — READ-ONLY ground truth for the SSC/ZFS project

> **DO NOT EDIT to make a test pass** (`CLAUDE.md §8`). These are reference numbers, geometries,
> basis/functional specs, and tolerances extracted from the papers in `./papers/`
> (`/bighome/vova/coding/SSC/papers/`). If code disagrees with a value here, **fix the code**.
> Entries marked **[CONFIRM]** must be verified against the cited source before being used as a gate.

Units: `D`, `E`, `D^SS` in cm⁻¹ unless noted. Fine-structure constant α ≈ 1/137.036 (a.u.).

---

## L1 — Integral finite-difference test (internal, no external number)
- Operator kernel (rank-2, symmetric, traceless): `T_kl(r₁₂) = (3 r₁₂,k r₁₂,l − δ_kl r₁₂²) / r₁₂⁵`.
- Identity: `T_kl = ∂²(1/r₁₂)/∂r₁₂,k ∂r₁₂,l` for `k≠l`, and the traceless diagonal combination for `k=l`
  (`Σ_k T_kk = 0` analytically).
- **Gate:** each component reproduced by finite differences of the existing ERI engine to **6–8
  significant figures**; **`Σ_k T_kk = 0` to ≤ 1e-10** for every shell quartet tested.
- Suggested FD: central differences on the ERI charge-cloud separation / kernel argument; step
  h ≈ 1e-4 a.u.; expect error ~h² (Richardson-extrapolate if needed for 8 figs).

## L2 — Single-determinant ROHF (McWeeny–Mizuno EXACT here); pins `C`/sign

### O₂  ³Σ_g⁻   (PRIMARY pin)
- **Geometry:** bond length **r(O–O) = 1.207 Å** (equilibrium; used by Vahtras et al. and Sinnecker–Neese).
- **Target SS-only `D^SS`:**
  - **1.44 cm⁻¹** — Vahtras et al. CASSCF, aug-cc-pCVTZ, CAS(10e,12o) [Sinnecker–Neese ref 16].
  - Sinnecker–Neese, single-point `D^SS` (cm⁻¹):

    | basis    | BP   | B3LYP | CASSCF |
    |----------|------|-------|--------|
    | EPR-II   | 1.52 | 1.53  | 1.55   |
    | EPR-III  | 1.57 | 1.58  | 1.57   |
    | QZVP     | 1.58 | 1.59  | 1.57   |

  - **Acceptance band for pinning: `D^SS ≈ 1.44–1.6 cm⁻¹`, sign POSITIVE.** Reproduce the value of
    the *method you run* to a few %.
- **WARNING:** the *total* O₂ `D ≈ 3.96–4.0 cm⁻¹` (Sinnecker–Neese 3.96; Neese 2007 ~4.0) **includes
  the second-order SOC term** and is **NOT** the SS-only target. Pin against the **SS-only ≈ 1.44–1.6**.
- Restricted open-shell (RODFT/ROHF) reference required; UKS/UB3LYP overestimates (see notes).

### CH₂  ³B₁   (SECONDARY single-determinant check)
- McWeeny–Mizuno applied to a single KS/HF determinant; first used by Petrenko & Neese
  [Petrenko, Petrenko, Bratus, J. Phys. Condens. Matter 14 (2002); Neese JACS 2006 ref 62].
- **`D^SS` value and geometry: [CONFIRM]** from Petrenko–Neese / Neese 2007 before using as a gate.
  (Well-known experimental total `D(CH₂ ³B₁) ≈ 0.76–0.79 cm⁻¹`; the SS-only theoretical value must be
  taken from the cited source, not assumed equal to experiment.)

## L3 — MRSF triplets

### Polyacenes — EXPERIMENTAL `D` (cm⁻¹)  [Sinnecker–Neese Table 2, refs 52,53]
| molecule    | experiment | ROBP/EPR-III (DFT) | CASSCF (mean-field) |
|-------------|-----------:|-------------------:|--------------------:|
| benzene     | **0.1593** | 0.159              | 0.146               |
| naphthalene | **0.1004** | 0.052              | 0.068               |
| anthracene  | **0.0702** | 0.042              | 0.048               |
| tetracene   | **0.0573** | 0.031              | 0.039               |
- **Caveat:** plain RO-DFT underestimates the larger acenes by ≈ factor 2 (ROHF recovers only ~½ of
  `D^SS` for these — Loboda et al.). When validating MRSF, state which reference column is compared
  and why; MRSF (multireference character) is expected to do better than single-reference RO-DFT.

### Organic radicals / carbenes / biradicals  [Sinnecker–Neese, radicals 1–15, Fig. 2/3]
- Diverse triplet set with `D < 0.10 cm⁻¹` (1,3- and 1,5-diradicals, dinitroxides, nitroxides,
  chlorophyll-a model, etc.).
- **Target accuracy:** **RMSD ≈ 0.0035 cm⁻¹** (BP/EPR-II, restricted-open-shell type) vs experiment.
  (UBP/EPR-II RMSD is much worse, 0.0772 cm⁻¹ — the RO vs U gap, see notes.)
- Carbenes/biradicals with larger `D` (radicals 16–22), `D^SS` (cm⁻¹), Sinnecker–Neese Table 3:
  16: 0.41–0.42 · 17: 0.23–0.24 · 18: 0.47–0.50 · 19: 0.52–0.54 · 20: 0.33–0.36 · 21: 0.30–0.34 ·
  22: 0.13–0.14 (BP/B3LYP × EPR-II/III spread; experimental in same table).

## L4 — Tensor invariants (free, every level)
- **`trace(D) = D_xx + D_yy + D_zz = 0`** (traceless kernel).
- **`E/D ∈ [0, 1/3]`**, with `D = D_zz − ½(D_xx + D_yy)`, `E = ½(D_xx − D_yy)`, using the
  convention `|D_zz| ≥ |D_yy| ≥ |D_xx|` (i.e. `0 ≤ E/D ≤ 1/3`).

---

## Method / convention notes (context, not gates)
- **Basis sets used in the sources:** EPR-II, EPR-III (Barone) — Sinnecker–Neese; **def2-TZVPP** —
  Neese 2007. Functionals: **BP86**, **B3LYP** (ROBP / ROB3LYP for restricted-open-shell).
- **Reference-determinant effect (critical):** UKS/UHF **overestimates** `D^SS`; **ROKS ≈ UNO**
  (spin-unrestricted natural-orbital) determinant — Neese 2007 §E. Prefer RO/UNO densities.
- **Mean-field for multiconfigurational WFs:** the 2-particle spin density factorizes into
  1-particle densities (exact for a single determinant; mean-field approximation otherwise) —
  Sinnecker–Neese eq 9 discussion, Neese JACS 2006. The 2e SS integrals themselves are computed
  **exactly** (not approximated) in ORCA — same intent here.

## Source map (`./papers/`)
- `spin-spin-contributions-...density.pdf` (+ SI `jp0643303_si_001.pdf`) = **Sinnecker & Neese,
  J. Phys. Chem. A 2006, 110, 12267** (canonical eq 9; O₂; acenes; radicals).
- `164112_1_online.pdf` = **Neese, J. Chem. Phys. 127, 164112 (2007)** (eq 46; SCF form; units; ROKS≈UNO).
- `2006-importance-...transition.pdf` (+ SI `ja061798a_si_001.pdf`) = **Neese, JACS 2006, 128, 10213**
  (mean-field; Mn(acac)₃; exact 2e SS integral code).
- `jcp-151-034106.pdf` = **Pokhilko, Epifanovsky & Krylov, J. Chem. Phys. 151, 034106 (2019)**
  (Wigner–Eckart M_S extraction from one reduced density).
