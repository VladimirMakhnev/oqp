Claude Code Task Brief — OQP Excited-State Analysis & Export Suite

Environment: oqp-misc (create it; do not reuse any existing OQP env)
Branch: create misc-excited-analysis from a fresh clone of https://github.com/Open-Quantum-Platform/openqp.git
Report target: Overleaf project 6a3407b308022f03c2a24efc (via git)

You are working autonomously and must keep working until every acceptance gate below is green or you hit a genuine, documented blocker. Commit locally after each passing gate so progress is durable and resumable.


0. Mission

Add to OQP a set of excited-state analysis and interoperability features, prove each one numerically on small benchmark molecules, and write a report to the Overleaf project. Target features:


NTOs (natural transition orbitals)
Attachment/detachment densities (unrelaxed; relaxed if the gradient relaxed density is reachable)
Transition densities (state-to-state 1-TDM, AO basis)
State-character descriptors (NTO participation ratio PR, Tozer Λ, fragment charge-transfer Ω-matrix)
Cube export (orbitals, state densities, transition density, NTO/attach/detach)
QCSchema export (AtomicResult JSON)
FCIDUMP export
Cross-check utilities (parse Gaussian/ORCA/Q-Chem/PySCF outputs and diff against OQP)


The scientific keystone is exposing the MRSF one-particle (transition) density matrix as a first-class object. Features 1–5 all derive from it. Build it first and gate everything on it (§4).


1. Operating principles (non-negotiable)


Verify-first. No feature is "done" until its numeric self-check passes (each feature below names its check). Numeric probes over code audits. Do not mark a gate green on the basis of "the code looks right."
Red-gate stop rule. If a gate's check fails, iterate on that gate. Do not proceed past a red gate. Do not weaken a tolerance to force a pass — if you change a tolerance, justify it in PROGRESS.md with the physical reason.
Python-first, surgical Fortran. Implement everything possible in new pyoqp modules operating on matrices read from the mol data store. Touch Fortran only to (a) write the 1-TDM/1-RDM intermediates into the data store, and (b) expose AO-on-grid evaluation if no Python path exists. Do not refactor unrelated code, reformat files, or alter the gradient/SOC/build logic except as strictly required.
Respect the bridge. OQP uses CFFI with a structured mol data store for Fortran↔Python exchange. Expose new quantities by writing them into that store under clearly-named tags, mirroring exactly how MO coefficients / densities / transition dipoles are already stored and read in Python. Do not invent a parallel file-based I/O path.
Reproducibility. Fix geometries, basis, functional, and convergence thresholds. The report's numbers must be regenerable by re-running one script.
Resumability. Maintain PROGRESS.md at repo root: a checklist of gates with status, the commit hash that closed each, key file paths, and any open questions. Update it at every gate. Commit per green gate.
Re-anchor on the brief. As your very first action, save this brief verbatim to oqp_misc_TASK_BRIEF.md at repo root and commit it — it is your durable source of truth if the original prompt scrolls out of context. Re-read oqp_misc_TASK_BRIEF.md in full at the top of every phase, on every resume, and any time you are about to weaken a check, relax a tolerance, or take a shortcut. After each re-read, restate to yourself (in PROGRESS.md) the current gate and the rules you must still be obeying, and confirm you have not drifted out of scope. If your working understanding ever disagrees with the file at repo root, the file wins.
No silent scope changes. Keep a running CHANGELOG_misc.md listing every file created and every existing file edited, with a one-line reason per edit.



2. Environment setup (oqp-misc)

conda create -n oqp-misc python=3.11 -y
conda activate oqp-misc
# OQP build deps are system-level (gfortran/gcc/g++, cmake, ninja, BLAS/LAPACK) —
# confirm they exist; if missing, report exact missing package, do not silently skip.
pip install numpy scipy
pip install cclib qcelemental qcengine pyscf pytest h5py

Notes:


Network-dependent steps (clone, conda, pip) must run where the cluster has egress (login node / proxy). If a proxy is required, set https_proxy/http_proxy first. If a step is blocked, report the exact command and the deny reason — do not work around it silently.
pyscf is needed as a reference for the FCIDUMP cross-check; qcelemental validates the QCSchema output; cclib does the output parsing.



3. Phase 0 — Clone, build, baseline (GATE 0)


Fresh clone, create branch misc-excited-analysis. Immediately save this brief to oqp_misc_TASK_BRIEF.md at repo root and commit it — this is the durable copy you will re-read at the top of every phase and on every resume (see §1, "Re-anchor on the brief").
Read the repo's own build docs (README, INSTALL, CMakeLists.txt, pyoqp/) and follow them. Reference build (confirm against current docs before using):


   cd openqp
   cmake -B build -G Ninja -DUSE_LIBINT=OFF \
     -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DCMAKE_Fortran_COMPILER=gfortran \
     -DCMAKE_INSTALL_PREFIX=. -DENABLE_OPENMP=ON -DLINALG_LIB_INT64=OFF
   ninja -C build install
   cd pyoqp && pip install . && cd ..

Set OPENQP_ROOT / any env the docs require.
3. Run an existing MRSF-TDDFT example from examples/ (or tests/) end-to-end. Confirm it converges and prints excitation energies + oscillator strengths.
4. Run the existing test suite and record the baseline pass/fail set. This is your regression baseline — re-run it after your changes (GATE 7).

GATE 0 passes when: clean build in oqp-misc, a stock MRSF example runs to completion, and the existing test suite baseline is recorded. Do not proceed otherwise.


4. Phase 1–2 — Inventory + the keystone (GATE 1, GATE 2)

Phase 1 — Inventory (GATE 1)

Locate and document in PROGRESS.md:


Where MRSF response amplitudes (X / Z-vectors) and MO coefficients (Cα, Cβ) are stored and how Python reads them from mol.
Where transition dipoles between states are computed (grep the MRSF modules for dipole, tdip, trans, oscillator). This is where the 1-TDM intermediate lives or is one contraction away.
Where the gradient relaxed density is formed (for relaxed attach/detach later) — the tdhf_mrsf_z_vector / gradient path.
The exact mechanism and naming convention for adding a new tagged array to the mol store and reading it via CFFI in Python.


GATE 1 passes when: you can, from a running Python session, pull Cα, Cβ, the response amplitudes, and OQP's computed transition dipole for a target state, and you have identified the minimal exposure point for the 1-TDM.

Phase 2 — Expose the MRSF 1-TDM / 1-RDM (GATE 2 — THE KEYSTONE)

Construct/expose:


γ⁰→ⁿ — the one-particle transition density matrix between the MRSF ground root S0 and excited root Sn, in MO and AO basis.
γⁿ — the excited-state 1-RDM (and γ⁰ for the ground root).


Caveats to handle explicitly (these are why standard TDDFT formulas don't transfer):


MRSF S0 is itself a response root, so S0→Sn is a genuine state-interaction TDM, not a reference→amplitude object.
Spin-flip fixes which spin blocks participate (hole in α, particle in β, or per the implemented flip convention) — build the correct blocks, don't assume closed-shell structure.
Watch amplitude normalization (the same convention trap seen in SOC-vs-GAMESS): document the convention you find and assert it.


GATE 2 — hard correctness gate (everything depends on this):


Reconstruct the transition dipole µ⁰→ⁿ = −Tr(γ⁰→ⁿ_AO · r_AO) from your exposed matrix and match OQP's directly-computed transition dipole (per component) to ≤ 1e-6 a.u. for at least two states on the benchmark molecule.
Tr(γⁿ_AO · S) = N_electrons for the state density (and Tr(γ⁰→ⁿ · S) ≈ 0, traceless).


If the reconstructed dipole does not match, the TDM construction is wrong — fix it before doing anything else. Do not build NTOs/A-D/cubes on an unverified TDM.


5. Phase 3 — Python analysis modules (GATE 3)

Create pyoqp/analysis/ with:


nto.py — SVD of the occ×virt (spin-flip) block of γ⁰→ⁿ → hole/particle NTOs + singular values, spin-resolved.
Check: Σσ² equals the amplitude norm (≈1 for a normalized single root) to ≤1e-8; reconstructing γ⁰→ⁿ from the truncated NTO pairs and re-deriving µ⁰→ⁿ reproduces the GATE 2 value.
density_diff.py — attachment/detachment from Δ = γⁿ − γ⁰; diagonalize, split by eigenvalue sign.
Check: Tr(A·S) ≈ Tr(D·S) ≈ n_promoted (≈1 for a single excitation) to ≤1e-6; A − D == Δ.
Unrelaxed first (amplitudes only). Relaxed variant only if GATE 1 confirms the gradient relaxed density is reachable — reuse it, don't re-derive.
descriptors.py — NTO participation ratio PR = (Σσ²)²/Σσ⁴; Tozer Λ (NTO overlap-weighted); fragment CT Ω-matrix from Löwdin-partitioned transition density given an atom→fragment map.
Check: PR ≥ 1 and ≤ N_pairs; Λ ∈ [0,1]; Ω row/col sums consistent with the transition normalization. Spot-check CT character against chemical expectation (e.g. n→π* of formaldehyde is local, low Λ-as-CT).


transition_density.py (the AO γ⁰→ⁿ writer/accessor) — already validated at GATE 2; expose it cleanly for the cube module.

GATE 3 passes when: all three modules' checks pass on formaldehyde S1 and at least one ethylene state.


6. Phase 4 — Cube export (GATE 4)

Create pyoqp/export/cubegen.py:


Generate a regular Cartesian box around the molecule (configurable padding/spacing; sane defaults).
Evaluate AOs on the grid. Prefer reusing OQP's existing AO-on-grid evaluator (used by the XC quadrature) via CFFI. Fallback: a small self-contained Python GTO evaluator built from the basis (shells/exponents/contraction coefficients) already in mol.
Write standard Gaussian .cube format for: a chosen MO, a state density γⁿ, the transition density γ⁰→ⁿ, and attachment/detachment/NTO densities.


Checks (do all):


If using the Python AO evaluator: numerically integrate χ_µ χ_ν over the grid and reproduce the AO overlap matrix S to ≤1e-3 (validates the evaluator before any physics).
Grid-integrate each density cube and compare to the analytic trace: state density → N (or n_promoted for A/D), transition density → ≈0. Match to ≤1e-2 (grid-resolution-limited; tighten by refining spacing).


GATE 4 passes when: cubes are produced and their grid integrals match analytic traces within tolerance.


7. Phase 5 — QCSchema + FCIDUMP export (GATE 5)

pyoqp/export/qcschema.py:


Serialize an MRSF calculation to a QCSchema AtomicResult (molecule, basis, method, driver, return_result, properties: SCF energy, excitation energies, oscillator strengths, dipole).
Check: qcelemental.models.AtomicResult(**payload) validates without error; energies in the payload equal OQP's own output to machine precision.


pyoqp/export/fcidump.py:


Dump MO 1e (h_pq) and 2e (pq|rs) integrals + E_nuc in FCIDUMP format on the high-spin reference orbitals. Use the UHF FCIDUMP layout (separate αα/ββ/αβ blocks) if dumping on the UMRSF/UHF reference; RHF-style if a closed-shell reference is used for the test.
Check (the clean cross-check): for H₂ / STO-3G, read the FCIDUMP into PySCF and run FCI; compare E_FCI to PySCF's own FCI for the same system to ≤1e-8 Hartree. Also assert 8-fold (or appropriate) permutational symmetry of the 2e integrals.


GATE 5 passes when: QCSchema validates and round-trips energies, and the FCIDUMP→PySCF-FCI energy matches reference.


8. Phase 6 — Cross-check utilities (GATE 6)

pyoqp/interop/:


parsers.py — wrap cclib to extract from Gaussian/ORCA/Q-Chem outputs: SCF energy, excitation energies, oscillator strengths, MO energies. PySCF: read native objects directly.
compare.py — a tolerance harness that diffs an OQP result against a parsed external result and prints a pass/fail table with deltas.


Check: parse a known reference output (generate a tiny Gaussian/ORCA/PySCF result for one of the benchmark molecules, or use any reference output available in the repo/test data) and confirm the extracted excitation energy matches the source value to ≤1e-4 eV. If no external code is available on the cluster, validate the parser against a committed sample output file and note the limitation.

GATE 6 passes when: the harness round-trips at least one external format and the comparison table renders.


9. Phase 7 — Regression + Benchmark / proof-of-concept (GATE 7)


Regression: re-run the existing test suite (GATE 0 baseline). Zero new failures — if any test regresses, fix it before proceeding.
Benchmark driver benchmark/run_poc.py:

Molecules: formaldehyde (clean n→π*), ethylene (π→π*), H₂/STO-3G (FCIDUMP cross-check). Geometries in Appendix A. Adapt an existing repo example for input format/keywords — substitute these geometries; do not guess the input syntax. Use a triplet (high-spin) reference for the MRSF runs; pick a standard functional (e.g. BHHLYP or a DTCAM variant the examples use) and a modest basis (6-31G* or cc-pVDZ). Fix all thresholds.
For each molecule: run MRSF, then invoke every feature and every self-check from §4–§8.
Emit benchmark/results.json (every check: name, value, tolerance, pass/fail) and benchmark/report_table.tex (a LaTeX tabular of the same, for the report to \input).



All checks must pass. A red check here means a real bug — go fix the responsible gate.


GATE 7 passes when: no regressions, and results.json shows every feature check green.


10. Phase 8 — Report to Overleaf (GATE 8)


Clone the project: git clone https://git.overleaf.com/6a3407b308022f03c2a24efc overleaf_report.

Auth: Overleaf git needs your Git authentication token (Account Settings → Git). It is expected to be available via a git credential helper or ~/.netrc. If clone/push fails on auth, do not abort the whole job — finish the report locally, stage and commit it, and write the exact git remote add / push commands plus the token requirement into PROGRESS.md for the user to run.



Inspect before writing. If the project is empty/skeleton, populate main.tex. If it already has content, add a new standalone file oqp_excited_state_analysis.tex and do not overwrite or reflow existing files (optionally add an \input{} to main only if there is an obvious, safe include point).
Report contents:

Motivation: the analysis/interop gap vs Q-Chem/ORCA/Q-Chem/PySCF and what this adds to OQP.
The keystone: the MRSF 1-TDM/1-RDM, the spin-flip / S0-as-root subtleties, and the normalization convention used.
Per feature: what was implemented, the validation method, and the numeric result — pull numbers from \input{report_table.tex}, do not hand-copy.
Benchmark setup (molecules, geometries, basis, functional, thresholds).
Limitations & next steps (e.g. relaxed A/D status, descriptor coverage, FCIDUMP reference choice, parser coverage).
Pointers: branch name, new file paths, how to regenerate (benchmark/run_poc.py).



Copy benchmark/report_table.tex (and any figures, e.g. an NTO/attach-detach cube rendering if you produce one) into the Overleaf repo. Commit with a clear message and push.


GATE 8 passes when: the report compiles-clean LaTeX with the auto-generated table is committed and pushed (or, on auth failure, committed locally with precise push instructions recorded).


11. Anticipated issues & mitigations


Build fails (missing gfortran/LAPACK/ninja, submodules). Read repo docs; report exact missing dependency; do not skip to features on a broken build (GATE 0 blocks).
1-TDM construction wrong for spin-flip. The GATE 2 dipole-reconstruction check catches this. Iterate until µ matches; suspect block selection and normalization first.
Normalization mismatch. Assert Σσ²/traces explicitly; document the convention. Known historical trap.
AO-on-grid not exposed to Python. Fallback Python GTO evaluator, validated by reproducing S (GATE 4).
FCIDUMP sign/ordering/permutation conventions. The PySCF-FCI round-trip on H₂/STO-3G is the ground truth — match the energy, not just the integral file shape.
Overleaf auth. Token required; on failure, commit locally + record commands (do not abort).
Network/proxy on cluster. Run net steps where egress exists; set proxy vars; report deny reasons rather than working around.
Long runtime / context limits. Gate-by-gate commits + PROGRESS.md make the job resumable; on resume, read PROGRESS.md first and continue from the last green gate.
Scope creep / breaking other modules. New code in new files; minimal existing-file edits logged in CHANGELOG_misc.md; GATE 7 regression run guards the rest of OQP.
Input-format guesswork. Always start from a working repo example input; never hand-author MRSF input syntax blind.



12. Definition of done (final checklist)


 oqp-misc env; clean build; stock MRSF example runs (GATE 0).
 1-TDM/1-RDM exposed; dipole reconstruction ≤1e-6; traces correct (GATE 2).
 NTOs, attach/detach, descriptors — all checks pass (GATE 3).
 Cube export — grid integrals match analytic traces (GATE 4).
 QCSchema validates + round-trips; FCIDUMP→PySCF-FCI matches on H₂/STO-3G (GATE 5).
 Cross-check parser/harness round-trips one external format (GATE 6).
 No regressions in existing test suite; results.json all green (GATE 7).
 Report committed/pushed to Overleaf with auto-generated table (GATE 8).
 PROGRESS.md, CHANGELOG_misc.md complete; final summary written.


Work through the gates in order. Stop only when all are green or you have a documented, genuine blocker — and if blocked, write a precise blocker report in PROGRESS.md (what failed, the exact command/output, what you tried, what input you need).


Appendix A — Benchmark geometries (Å)

Formaldehyde (n→π* test):

C   0.000000   0.000000  -0.529700
O   0.000000   0.000000   0.677500
H   0.000000   0.934200  -1.124000
H   0.000000  -0.934200  -1.124000

Ethylene (π→π* test):

C   0.000000   0.000000   0.668600
C   0.000000   0.000000  -0.668600
H   0.000000   0.922900   1.237200
H   0.000000  -0.922900   1.237200
H   0.000000   0.922900  -1.237200
H   0.000000  -0.922900  -1.237200

H₂ (FCIDUMP/STO-3G cross-check):

H   0.000000   0.000000   0.000000
H   0.000000   0.000000   0.741000

Use the MRSF high-spin (triplet) reference; fix functional, basis, and convergence thresholds and record them in the report.
