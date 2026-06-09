# Correspondence Map: UMRSF-TDDFT градиент ↔ существующий RO-MRSF код

**Фаза 0 — Reconnaissance. Deliverable. КОД НЕ ПИШЕТСЯ.**

Документ устанавливает соответствие между формулами `umrsf_gradient_theory.tex`
(далее — «.tex», ссылки вида §N / eq:label) и существующими subroutines OpenQP.
Цель — зафиксировать «золотой стандарт» RO-MRSF и точки расширения RO→UHF.

> ⚠️ Файл в задаче назван `umrsf_gradient_theory_en.tex`, фактически в репозитории —
> `umrsf_gradient_theory.tex` (русскоязычный). Содержание то же; работаю с фактическим.

---

## 0. Краткая карта файлов

| Файл | Роль |
|---|---|
| `source/tdhf_mrsf_lib.F90` (3131 стр.) | Все builders: cbc, mntoia, sp, qro*-builders, jacobi reordering, **уже есть** umrsf*-варианты для energy |
| `source/tdhf_sf_lib.F90` | SF/RO базовые операторы: sfrogen, sfrolhs, sfrorhs, sfropcal, sfromcal, sfdmat |
| `source/modules/tdhf_mrsf_z_vector.F90` (1265 стр.) | Z-vector драйвер `tdhf_mrsf_z_vector` + `apply_z_operator` (GMRES) |
| `source/modules/tdhf_mrsf_gradient.F90` (492 стр.) | Gradient драйвер `tdhf_mrsf_gradient` + Γ-сборка `grd2_mrsf_compute_data_t_get_density` |
| `source/modules/tdhf_mrsf_energy.F90` | Energy драйвер; ветвление `umrsf` уже реализовано |

---

## 1. ⚠️ КЛЮЧЕВАЯ НАХОДКА: два разных RO-пути, и `mrsfqro*` — это НЕ MRSF

В задаче `mrsfqroesum/mrsfqrorhs/mrsfqropcal/mrsfqrowcal` названы «золотым стандартом
MRSF». Фактически в `tdhf_mrsf_z_vector` они подключены к ветке **`mrst==5` —
чистый Spin-Flip TDDFT БЕЗ spin-pairing**. Настоящий MRSF-путь со spin-pairing
(`mrst==1` синглет / `mrst==3` триплет) использует ДРУГОЙ набор:

| Шаг | MRSF путь (mrst=1/3) — **с spin-pairing** | SF/«qro» путь (mrst=5) — **без spin-pairing** |
|---|---|---|
| Q/R (RHS Z) | `sfrorhs` (sf_lib:334) + `mrsfsp` (spin-pair) | `mrsfqrorhs` (mrsf_lib:2287) |
| H⁺[T] ERI builder | `mrsfcbc` (mrsf_lib:497) | `mrsfqroesum` (mrsf_lib:1538) — суммирует H⁺ в pmo |
| H[X,X] spin-pairing | `mrsfsp` (mrsf_lib:1677) | — (нет) |
| P = T+Z | `sfropcal` (sf_lib:693) | `mrsfqropcal` (mrsf_lib:2364) |
| W builder | `mrsfrowcal` (mrsf_lib:2144) | `mrsfqrowcal` (mrsf_lib:2424) |
| Z-оператор J·Z | `apply_z_operator`→`sfrogen`/`sfrolhs` | то же |
| Γ density | `grd2_mrsf_compute_data_t_get_density` (вкл. spin-pair) | `sf_2e_grad` |
| Driver | `mrsf_2e_grad` | `sf_2e_grad` |

**Вывод для проекта.** Теория .tex (§8–§13) описывает **полный MRSF со spin-pairing**
(Γ^SP, §12). Следовательно истинный референс для UMRSF-градиента — путь **mrst=1/3**:
`sfrorhs + mrsfsp + sfropcal + mrsfrowcal + mrsfcbc + Γ(get_density)`.
Routines `mrsfqro*` полезны как **шаблон явной декомпозиции Q/R/P/W** (чище записаны,
без spin-pairing), но физику spin-pairing в них надо добавить из mrst=1/3 пути.

> Это нужно учесть в названиях: задача предлагает `umrsfqrorhs/umrsfqropcal/umrsfqrowcal`.
> Содержательно они должны соответствовать **mrst=1/3** математике (sfrorhs+mrsfsp/
> mrsfrowcal), а не «голому» SF mrsfqro*. Префикс `umrsf*` оставляю; вопрос организации — ниже.

**Подтверждение из energy:** UMRSF energy уже расширил именно mrst=1/3 routines:
`mrsfcbc→umrsfcbc`, `mrsfmntoia→umrsfmntoia`, `mrsfsp→umrsfsp`, `mrsfssqu→umrsfssqu`,
`sfdmat→umrsfdmat`. Энергетический Z-vector/gradient не делался — это наш предмет.

---

## 2. Орбитальные пространства и размерности (.tex §2)

UKS, M_S=+1 референс. В коде:

| .tex | код (`tdhf_mrsf_z_vector` стр.625–630) | смысл |
|---|---|---|
| n_C (closed) | `noccb = nelec_B` | дважды занятые |
| n_O (open, SOMO) | `nsocc = nocca - noccb` | 2 SOMO (O₁,O₂) |
| n_V (virtual) | `nvira = nbf - nocca` | α-виртуальные |
| α-occ = C∪O | `nocca = nelec_A` | |
| β-virt = O∪V | `nvirb = nbf - noccb` | |

---

## 3. Z-vector: размерность и блоки (.tex §6.3, eq:Zdim)

**RO (текущий код, стр.630):**
```
lzdim = noccb*(nsocc+nvira) + nsocc*nvira
      = n_C·n_O + n_C·n_V + n_O·n_V         (3 блока: CO, CV, OV)
```
Порядок блоков в `rhs`/`z` (из `mrsfqrorhs`/`sfrolhs`):
1. **doc-socc (CO, ix):** `n_C·n_O`
2. **doc-virt (CV, ia):** `n_C·n_V`
3. **soc-virt (OV, xa):** `n_O·n_V`

**UMRSF (eq:Zdim, §6.3) — целевое:**
```
dim Z = n_O·n_V + n_C·n_O + 2·n_C·n_V        (4 активных блока)
```
Активные блоки (.tex Проверка 5):
| Блок | спин | RO-аналог |
|---|---|---|
| Z_ix_β (CO) | β | doc-socc |
| Z_xa_α (OV) | α | soc-virt |
| Z_ia_α (CV) | α | doc-virt (α-часть) |
| Z_ia_β (CV) | β | **новый** — удвоение CV |

> Разница RO→UMRSF = **+n_C·n_V** (CV-блок удваивается по спину; CO становится
> β-only, OV — α-only). ⚠️ Подводный камень из задачи: НЕ 6·n_occ·n_vir.

---

## 4. Таблица соответствия: формула .tex → subroutine

### 4.1 Unrelaxed difference density T (eq:Talpha/Tbeta, §7.1)
| .tex | код |
|---|---|
| Tα(i,j) = −Σ X X | `tij` через `dgemm` в z_vector стр.887–890 |
| Tβ(a,b) = +Σ X X | `tab` через `dgemm` стр.893–896 |

### 4.2 Операторы H⁺ и H[X,X]
| .tex | формула | RO subroutine |
|---|---|---|
| H⁺_tuσ[V] | eq:Hplus §7.2 | `mrsfcbc` (ERI A⁰) + `mntoia` + `apply_z_operator`/`sfrolhs` |
| H⁽⁰⁾, H^intra, H^inter | eq:H0/§7.4 | `mrsfsp` (mrsf_lib:1677) ← spin-pairing |
| H^σ′σ̃ (cross-spin, eq:Htilde §5.3) | — | `sfrolhs` ядро; **требует расширения на cross-spin Coulomb+f^xc_αβ** |

### 4.3 Q-builder (.tex §8, eq:Qixa…)
| .tex | RO |
|---|---|
| Q^(k)_tuσ, 18 spin-block элементов | inline: `sfrorhs`(собирает Q→R) + `mrsfqroesum`(H⁺ суммирование в pmo) |

### 4.4 R-builder = RHS Z-vector (.tex §9, eq:Rdef)
| .tex | RO |
|---|---|
| R_pqσ = Q_pqσ − Q_qpσ, 6 блоков | `sfrorhs` (sf_lib:334, mrst=1/3) / `mrsfqrorhs` (mrsf_lib:2287, mrst=5) |

### 4.5 Z-vector solve J·Z=−R (.tex §6, eq:Zvec)
| .tex | RO |
|---|---|
| генерация плотности из Z | `sfrogen` (sf_lib:463) |
| J·Z (LHS) | `apply_z_operator` (z_vector:370) → ERI + `utddft_fxc` + `sfrolhs` |
| precond | `sfromcal` (диаг) + `apply_z_precond` |
| solver | `gmres_solve` (z_vector:76) |

### 4.6 P = T+Z (.tex §11, eq:Pdef)
| .tex | RO |
|---|---|
| P_pqσ = T_pqσ + Z_pqσ | `sfropcal` (sf_lib:693) / `mrsfqropcal` (mrsf_lib:2364) |

### 4.7 W-builder (.tex §10, eq W_ix…)
| .tex | RO |
|---|---|
| W_tuσ через Q,Z,Fock | `mrsfrowcal` (mrsf_lib:2144) / `mrsfqrowcal` (mrsf_lib:2424) |

### 4.8 AO-сборка плотностей (.tex §11)
| .tex | RO |
|---|---|
| P_μνσ, D_μνσ, W_μνσ | `orthogonal_transform` (z_vector:1161–1164,1238) |
| X_μν (SF transition) | `mrsfxvec`+`iatogen`/`sfdmat`(`umrsfdmat` готов) |
| блочные X^intra/X^inter | `mrsfcbc`/`umrsfcbc` → spc(1:7) density-likes |

### 4.9 Γ two-particle density (.tex §12)
| .tex | RO (`grd2_mrsf_compute_data_t_get_density`, mrsf_gradient:300) |
|---|---|
| Γ^SF (2PD − c_HF PD − c_HF XX), eq:GammaSF | `df1`/`dq1`/`dt2`, стр.358–376 |
| Γ^SP intra αβ (eq §12.2.1) | `o21v`,`co12` блоки `db1`,`db2`, стр.379–403 |
| Γ^SP inter αα / ββ (eq §12.2.2/3) | `bco1/bco2/bo1v/bo2v` → `dc*`/`dd*`, стр.405–480 |

> ⚠️ **Новизна UMRSF:** в RO inter-блок один (`co12`,`o21v`), αα и ββ совпадают.
> В UMRSF (§12, Проверка 4) их надо **разделить** на αα и ββ (независимые MO).
> `d2(:,:,1)`/`d2(:,:,2)` в RO строятся из одного `mo_a`; в UMRSF — из `mo_a`/`mo_b`.

### 4.10 Финальный градиент Ω^ξ (.tex §13, eq:OmegaXi)
| .tex член | RO |
|---|---|
| h^ξ P (1e) | `sf_1e_grad` |
| −S^ξ W (Pulay) | `eijden`/Γ-driver, W из `wao` |
| (μν\|κλ)^ξ Γ (2e) | `mrsf_2e_grad` + Γ density |
| Ω^ξ_xc (вкл. f^xc_αβ) | `utddft_xc_gradient` (mrsf_gradient:156) |

---

## 5. RO→UHF паттерн расширения (из готового energy-кода)

UMRSF energy показывает шаблон, который повторим в градиенте:
1. Флаг `infos%tddft%umrsf`, требование `SCFTYPE=2` (UHF), `mult=3`.
2. Везде, где RO передаёт `mo_a, mo_a` → UMRSF передаёт `mo_a, mo_b` (независимые β).
3. Fock: RO усредняет (Guest-Saunders веса {2,1,1,2}); UMRSF — `fock_a`,`fock_b`
   независимо, все веса Лагранжиана = 1 (.tex §4.3).
4. Новые routines: `int2_umrsf_data_t_update` (mrsf_lib:218) уже строит UHF ERI-вклады.
5. ✅ Готово в energy: `umrsfcbc, umrsfmntoia, umrsfsp, umrsfssqu, umrsfdmat`.

---

## 6. MO reordering (Фаза 1) — статус

| subroutine | стр. | роль |
|---|---|---|
| `get_jacobi` | mrsf_lib:2708 | α/β выравнивание через max-overlap (S_MO) |
| `rotate_pair` | 2863 | jacobi rotation пары MO |
| `swap_sign_a` | 2907 | фиксация знака |
| `check_sign` | 2927 | проверка знака α/β |

Используется в UMRSF energy (`umrsf_jac`). Для градиента критична **гладкость между
геометриями** (.tex §14, «MO reordering — гладкость PES»). Фаза 1 = проверить
smoothness, при необходимости — unitary continuation. **На Фазе 0 не трогаю.**

---

## 7. Предлагаемые новые UMRSF subroutines (план, реализация — Фазы 2–5)

| Фаза | Новое | Аналог-шаблон | Формулы .tex |
|---|---|---|---|
| 2 | `umrsfqrorhs` (Q,R) | `sfrorhs`+`mrsfsp` (mrst=1/3) | §8, §9 |
| 3 | расширить `apply_z_operator`/новый `umrsf_sfrolhs` на 4 блока + cross-spin | `sfrolhs`/`sfrogen` | §5.3 eq:Htilde, §6 |
| 4 | `umrsfqropcal` (P) | `sfropcal` | §11 eq:Pdef |
| 4 | `umrsfqrowcal` (W) | `mrsfrowcal` | §10 |
| 5 | расширить Γ: split inter αα/ββ | `grd2_mrsf_..._get_density` | §12 eq:GammaSF + SP |
| 5 | `tdhf_umrsf_z_vector`, `tdhf_umrsf_gradient` драйверы | `tdhf_mrsf_z_vector`/`tdhf_mrsf_gradient` | §13 |

---

## 8. Известные подводные камни (сверка с задачей)

- ✅ Z-dim = n_O·n_V + n_C·n_O + 2·n_C·n_V (НЕ 6·n_occ·n_vir).
- ✅ W^(k)_iaα / W^(k)_xaα: член 2H_aiα[X,X]=0 удалён (.tex Проверка 2).
- ✅ Spin-pairing inter: split αα/ββ — новизна UMRSF (§12).
- ✅ f^xc_αβ cross-spin XC kernel — обязателен в Z-solver и Ω^ξ_xc.
   В коде: `utddft_fxc` (z_vector:818) и `utddft_xc_gradient` уже UHF-aware (fxa/fxb,
   dxa/dxb раздельно) — точка для cross-spin.
- ✅ Gradient_of_MRSFUHFv3.pdf НЕ использовать (в репозитории отсутствует — ок).
- ⚠️ **Новое (раздел 1):** `mrsfqro*` = SF-путь без spin-pairing; истинный референс —
  mrst=1/3 (`sfrorhs+mrsfsp+mrsfrowcal`). Нужно решение по неймингу/организации.

---

## 9. Принятые решения (после Фазы 0)

| Вопрос | Решение |
|---|---|
| Где код | **Расширять `tdhf_mrsf_lib.F90`** (рядом с umrsf energy routines) |
| Референс builders | **Полный MRSF spin-pairing, путь mrst=1/3** (`sfrorhs+mrsfsp+sfropcal+mrsfrowcal+Γ^SP`), НЕ голый SF `mrsfqro*` |
| Functional тестов | **BHHLYP** (c_HF=0.5, cross-spin f^xc_αβ тестируется) |
| Basis | 6-31G(d) |
| Запуск | **Сборка cmake + pyoqp доступны** — numerical-gradient тесты запускаю сам |

## Статус Фазы 0
✅ Reconnaissance завершён. Прочитаны: вся теория (.tex, 1324 стр.), RO Q/R/P/W builders,
`apply_z_operator`, Z-vector driver, Γ-сборка, UMRSF energy routines, jacobi reordering.
Код не написан. Решения по организации зафиксированы. Следующий шаг — Фаза 1 (проверка
smoothness MO reordering) + план Фазы 2, оба ждут одобрения перед написанием кода.
