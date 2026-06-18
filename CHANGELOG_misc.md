# CHANGELOG (misc-excited-analysis)

Every file created or edited for this task, with a one-line reason.

## New files
- `oqp_misc_TASK_BRIEF.md` — verbatim task brief (durable source of truth).
- `PROGRESS.md` — gate checklist, inventory, env notes, resumability.
- `CHANGELOG_misc.md` — this file.

- `pyoqp/oqp/analysis/__init__.py`, `transition_density.py`, `nto.py`, `density_diff.py`,
  `descriptors.py`, `gto_grid.py` — GATE 3 analysis package (NTOs, attach/detach, PR/Tozer-Λ/
  fragment-Ω descriptors, self-contained Cartesian-GTO grid evaluator). Reads the exposed
  tags read-only via the bridge; no existing module touched. (Placed under `pyoqp/oqp/` so it
  installs/imports as `oqp.analysis`.)
- `benchmark/_gate3_check.py` — GATE 3 numeric self-check harness.

- `pyoqp/oqp/export/__init__.py`, `cubegen.py` — GATE 4 Gaussian-cube export (MO, state/
  transition densities, attach/detach, NTO) + grid-integral trace checks.
- `benchmark/_gate4_check.py` — GATE 4 self-check harness.

- `pyoqp/oqp/export/qcschema.py` — GATE 5 QCSchema AtomicResult export (+validate).
- `pyoqp/oqp/export/fcidump.py` — GATE 5 FCIDUMP export (PySCF AO engine in OQP MO basis)
  + PySCF-FCI round-trip verifier.
- `benchmark/_gate5_check.py`, `benchmark/inputs/h2_rhf_sto3g.inp` — GATE 5 self-check + input.

- `pyoqp/oqp/interop/__init__.py`, `parsers.py`, `compare.py` — GATE 6 cross-check utilities
  (cclib wrapper for Gaussian/ORCA/Q-Chem outputs; PySCF native reader; tolerance harness).
- `benchmark/_gate6_check.py`, `benchmark/fixtures/ch2o_gaussian_td.log` — GATE 6 check + fixture.

- `benchmark/run_poc.py` — GATE 7 reproducible benchmark driver (all GATE 2-6 checks over
  formaldehyde/ethylene/H2); emits `benchmark/results.json` + `benchmark/report_table.tex`.
- `benchmark/make_figure.py`, `benchmark/figures/ch2o_S1_attach_detach.png` — attach/detach figure.
- `benchmark/results.json`, `benchmark/report_table.tex` — auto-generated POC results (41/41).

## Edited existing files
- `source/tagarray_driver.F90` — added 3 tag-name constants + comments
  (`OQP_td_trans_density_mo`, `OQP_td_trans_dipole`, `OQP_td_dip_ao`) to expose the MRSF
  state-interaction densities, transition dipoles, and AO dipole integrals. No existing tag touched.
- `source/modules/tdhf_mrsf_energy.F90` — added `use int1, only: multipole_integrals`, three
  pointer/scratch declarations, and an additive write-out block right after `get_transition_dipole`
  that stores `trden`, `dip`, and the AO dipole integrals into the tagarray. No physics altered.
- `pyoqp/oqp/utils/oqp_tester.py` — replaced `ProcessPoolExecutor(max_tasks_per_child=1)`
  (Python 3.11+ only) with a non-daemonic `NoDaemonPool(maxtasksperchild=1)` so the example
  test runner works on the cluster's Python 3.9 while preserving one-calc-per-process isolation
  and allowing Hessian/NAC tests to spawn their own children. Test-harness portability only.
