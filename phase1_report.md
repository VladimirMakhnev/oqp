# Фаза 1 — отчёт: MO reordering smoothness для UMRSF

**Статус: READY (для H₂O) — reordering гладкий, расширение не требуется.**

## Окружение
Все runtime-изменения изолированы в conda-env **`grad_env`** (python 3.11):
- native `liboqp.so` собрана (cmake/ninja, ILP64 OpenBLAS `libopenblas64_`).
- Решена проблема сборки: отсутствовал `pyconfig.h` для python3.11 (multilib stub) →
  собрали cffi-binding против python из `grad_env` (полные dev-заголовки).
- Решён конфликт MPI: `~/.local` содержал сломанный `mpi4py` → поставлен self-contained
  `openmpi`+`mpi4py` в grad_env; запуск с `PYTHONNOUSERSITE=1`.
- Доустановлены `basis_set_exchange`, `pyscf`.

**Команда запуска:**
```bash
source /opt/conda/etc/profile.d/conda.sh && conda activate grad_env
export PYTHONNOUSERSITE=1 OMP_NUM_THREADS=4
openqp <input.inp>
```

## Smoke-test: UMRSF energy H₂O (BHHLYP/6-31G*, UHF mult=3)
- UHF SCF = **−76.0022888549** Ha (13 iter).
- UMRSF state 1 (синглетное осн. сост.) = **−76.3608777063** Ha, ⟨S²⟩=0.001.
- Davidson сошёлся (err 7e-12). UMRSF-путь (`tdhf.type=umrsf`) отработал штатно.

## Smoothness probe (`tests/umrsf_grad/fd_probe.py`)
Смещение O вдоль z на ±0.005 Å, 3 точки, target state 1:

| dz (Å) | E(state 1), Ha | Jacobi segments (iter) |
|---|---|---|
| −0.005 | −76.3609906396 | seg0→4, seg1→28 |
|  0.000 | −76.3608777064 | seg0→4, seg1→28 |
| +0.005 | −76.3607206427 | seg0→4, seg1→28 |

- **PES гладкая, параболическая, без скачков.**
- **Jacobi reordering сходится идентично во всех точках** → нет discrete MO switch
  между близкими геометриями (ключевой риск из §14 .tex — не реализовался для H₂O).
- FD gradient g_z(O) = **0.01428762 Ha/Bohr**; кривизна 1.765 Ha/Å² (мала).

## Вывод и оговорки
- Для H₂O reordering **smooth → помечаю «ready», расширение (unitary continuation)
  НЕ требуется** на текущем этапе.
- Оговорка: проверено на одной малой closed-shell системе и малом окне смещения.
  Для O₂ (триплет) и ethylene@90° (conical intersection, Фаза 6) риск discrete
  switch выше — **диагностику повторить на этих системах в Фазе 6**; если там
  появятся скачки, вернуться к unitary continuation (отдельный PR).
- FD-эталонная машинерия (`fd_probe.py`) готова и переиспользуется в тестах Фазы 2
  (E.1) и Фазы 6.

## Что дальше
Фаза 2 (Q/R builders) — план в `phase2_plan.md`. Перед кодом ждут подтверждения 3
решения: (1) разделять Q/R, (2) порядок блоков rhs, (3) тест-обвязка.
