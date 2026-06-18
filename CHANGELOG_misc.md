# CHANGELOG (misc-excited-analysis)

Every file created or edited for this task, with a one-line reason.

## New files
- `oqp_misc_TASK_BRIEF.md` — verbatim task brief (durable source of truth).
- `PROGRESS.md` — gate checklist, inventory, env notes, resumability.
- `CHANGELOG_misc.md` — this file.

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
