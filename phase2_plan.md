# Фаза 2 — расширенный план: Q/R builders для UMRSF

**Референс теории:** русская версия `umrsf_gradient_theory.tex` (репо-корень).
Формулы: §8 (Q, eq:Qixa…eq:Qab), §9 (R, eq:Rdef…), §7 (T, H⁺, H[X,X]).
Все новые routines — в `source/tdhf_mrsf_lib.F90`, префикс `umrsf*`, существующие
`mrsf*`/`sf*` не трогаются. Модуль без `private` → экспорт автоматический.

> ⚠️ Фаза 2 не включает Z-solve/P/W — только Q и R на ОДНОЙ геометрии, чтобы
> проверить их против finite-difference ещё до Z-vector (Фаза 3).

---

## A. Данные на входе (уже доступны / уже строятся)

Поток в mrst=1/3 ветке (RO, `tdhf_mrsf_z_vector.F90:837–899`), который воспроизводим
для UHF. Все эти объекты УЖЕ есть либо строятся готовыми umrsf energy routines:

| Объект | Что | Источник (reuse) |
|---|---|---|
| `mo_a, mo_b` | независимые α/β MO | tagarray `OQP_VEC_MO_A/B` |
| `fa, fb` | UKS Fock в MO (α/β), `nbf×nbf` | `orthogonal_transform_sym`(fock_a/b) — для UHF брать `fock_b`, не average |
| `tij` Tα, `tab` Tβ | unrelaxed diff density, eq:Talpha/Tbeta | `dgemm` на `bvec_mo_d` (как RO 887–896) |
| `ab1_mo_a, ab1_mo_b` | H⁺[T] α/β в MO, eq:Hplus | int2 на T-плотности + `mntoia` (reuse) |
| `fmrst2` (10-comp) | spin-resolved SP density-likes | **`umrsfcbc`** (готово) |
| `hxa, hxb` | H[X,X]: H⁽⁰⁾ + spin-pairing, eq:Hfull | **`umrsfsp`** (готово, принимает ca,cb,10-comp) |

→ **Новый код Фазы 2 = только сборка Q и R из этих готовых блоков.** Вся ERI/SP-физика
переиспользуется (umrsfcbc, umrsfsp, mntoia, int2_umrsf_data_t).

---

## B. Новые subroutines — полные сигнатуры

### B.1 `umrsfqcal` — Q-builder (§8, eq:Qixa–eq:Qab)
Вычисляет Q^(k)_{tu,σ} раздельно для α и β во всех 9 блоках. Вынесен отдельно от R,
чтобы (1) сделать покомпонентный Q-тест и Check-1 тривиальными, (2) переиспользовать
Q в W-builder (Фаза 4). Соответствует §8.

```fortran
subroutine umrsfqcal(qa, qb, hpta, hptb, hxa, hxb, fa, fb, tij, tab, noca, nocb)
  real(kind=dp), intent(out), dimension(:,:) :: qa, qb      ! Q_{tu,alpha}, Q_{tu,beta}, nbf×nbf (MO)
  real(kind=dp), intent(in),  dimension(:,:) :: hpta, hptb  ! H^+[T]_alpha, _beta  (= ab1_mo_a/b)
  real(kind=dp), intent(in),  dimension(:,:) :: hxa, hxb    ! H[X,X]_alpha, _beta (из umrsfsp + H0)
  real(kind=dp), intent(in),  dimension(:,:) :: fa, fb      ! UKS Fock MO alpha/beta, nbf×nbf
  real(kind=dp), intent(in),  dimension(:,:) :: tij, tab    ! T_alpha (C∪O), T_beta (O∪V)
  integer,       intent(in)                  :: noca, nocb
end subroutine
```
Тело: для каждого блока (ix,xi,ia,ai,xa,ax,ij,xy,ab) и спина выписать формулы §8 как
naive прямую сумму: `Q = H^+[T] + 2 H[X,X] + 2 F[X,X]` с занулениями из Проверки 1.
F[X,X] (eq:FX) — локальный bilinear `-X F_β X` (α) / `+X F_α X` (β), строится inline.
Комментарий над каждым блоком: «соответствует формуле (eq:Qixa) и т.д.».

### B.2 `umrsfqrorhs` — R-builder = RHS Z-vector (§9, eq:Rdef)
Формирует 4-блочный R = Q_{pq,σ} − Q_{qp,σ}, раскладку UMRSF (раздел C ниже).
Аналог `sfrorhs`, но spin-explicit и 4 блока вместо 3.

```fortran
subroutine umrsfqrorhs(rhs, qa, qb, noca, nocb)
  real(kind=dp), intent(out), dimension(:) :: rhs          ! длина lzdim_umrsf (см. C)
  real(kind=dp), intent(in),  dimension(:,:) :: qa, qb     ! из umrsfqcal
  integer,       intent(in)                  :: noca, nocb
end subroutine
```
Тело: 4 цикла (CO_β, OV_α, CV_α, CV_β), `rhs(ij) = -(Q_{pq,σ}-Q_{qp,σ})`.

> Альтернатива: слить B.1+B.2 в один `umrsfqrorhs(rhs, hpta,…)` как у RO. Предлагаю
> РАЗДЕЛИТЬ — это и есть «naive прямая реализация» + даёт чистые Q-тесты. Жду решения.

---

## B.3 Дизайн `umrsfqcal` (финализирован) + находка по теории

**⚠️ Находка (theory vs code):** боксы §9 (напр. eq:Ria_alpha = `2(H_iaα−H_aiα)`)
**опускают** член `Hp[T]`, который golden `sfrorhs` **включает** (hpta в doc-virt блоке).
§8.3 даёт `Q_iaα = Hp_iaα[T] + 2H_iaα`, поэтому `R_iaα = Q_iaα−Q_aiα` обязан нести
`Hp_iaα[T]`. **Доверяем коду (sfrorhs), а не боксам §9.** Поскольку `umrsfqrorhs`
строит R=Q−Qᵀ из полных qa/qb, корректность обеспечивается тем, что `umrsfqcal`
кладёт Hp[T] в нужные блоки Q (а R сам даёт правильную разность).

**Контракт `umrsfqcal` (naive, full nbf×nbf MO):**
```fortran
subroutine umrsfqcal(qa, qb, hpta, hptb, hxa, hxb, fa, fb, xmo, noca, nocb)
  real(dp), intent(out) :: qa(:,:), qb(:,:)        ! Q_{tu,alpha/beta}, nbf×nbf MO
  real(dp), intent(in)  :: hpta(:,:), hptb(:,:)    ! H^+[T]_alpha/beta, FULL nbf×nbf MO
  real(dp), intent(in)  :: hxa(:,:), hxb(:,:)      ! H[X,X]_alpha/beta native scale, nbf×nbf MO
  real(dp), intent(in)  :: fa(:,:), fb(:,:)        ! UKS Fock MO alpha/beta
  real(dp), intent(in)  :: xmo(:,:)                ! U·X transition ampl., (noca, nvirb)
  integer,  intent(in)  :: noca, nocb
end subroutine
```
Сборка (§8, eq:Qixa…eq:Qab):
- `qσ(t,u) = 2·hxσ(t,u)` во всех блоках (член 2H[X,X]).
- `+ hptσ(t,u)` **только при t ≤ noca** (t∈C∪O) — маска убирает Hp[T] из блоков
  ai/ax/ab (§8.4/8.6/8.9: «Hp[T] отсутствует, т.к. первый индекс ∈ V»).
- `+ 2·F[X,X]σ(t,u)` (eq:FX), вычисляется inline; авто-зануляется вне носителя X:
  - `Fa(t,u) = −Σ_{q,s∈O∪V} xmo(t,q−nocb)·fb(q,s)·xmo(u,s−nocb)`, t,u∈[1,noca].
  - `Fb(t,u) = +Σ_{p,r∈C∪O} xmo(p,t−nocb)·fa(p,r)·xmo(r,u−nocb)`, t,u∈[nocb+1,nbf].

**Резолюция factor-ов через RO-предел (НЕ из §8 вслепую):** внутренние факторы 2.0
в `umrsfsp`/H0-dgemm драйвера → точный масштаб `hxa/hxb` подгоняется так, чтобы
тест E.3/E.4 (`½(Rα+Rβ)` == `sfrorhs` при UKS≈ROHF) сошёлся <1e-5. Мини-драйвер
строит `hxa/hxb` зеркально RO-драйверу; финальный множитель фиксируется тестом.

## C. Размерность и раскладка rhs/Z (eq:Zdim)

```
lzdim_umrsf = n_C·n_O + n_O·n_V + 2·n_C·n_V         (vs RO: n_C·n_O+n_O·n_V+n_C·n_V)
```
с `noccb=n_C, nsocc=noca-nocb=n_O, nvira=nbf-noca=n_V`. Порядок блоков (предложение —
сохранить RO-подобный, CV-блок удвоить по спину):
1. `Z_ix_β` (CO,β): `i=1..n_C, x=nocb+1..noca` → n_C·n_O
2. `Z_ia_α` (CV,α): `i=1..n_C, a=noca+1..nbf` → n_C·n_V
3. `Z_ia_β` (CV,β): `i=1..n_C, a=noca+1..nbf` → n_C·n_V
4. `Z_xa_α` (OV,α): `x=nocb+1..noca, a=noca+1..nbf` → n_O·n_V

(Финальный порядок согласуем в Фазе 3, когда J-оператор будет читать ту же раскладку.)

---

## D. Точки вставки в `tdhf_mrsf_lib.F90` (номера строк текущего файла)

| Что | Место |
|---|---|
| `umrsfqcal`, `umrsfqrorhs` | **после `end subroutine mrsfqrowcal` (стр. 2556), вставка с 2557** — рядом с RO gradient builders, перед `get_mrsf_transition_density` (2558) |

Группировка: все umrsf-gradient builders подряд после RO-gradient builders (2287–2556).
Будущие `umrsfqropcal`/`umrsfqrowcal` (Фаза 4) добавятся в тот же блок.

Драйвер-вызовы (Фаза 2 — только для теста, без полного gradient драйвера): добавлю
**временный** диагностический путь либо отдельную тестовую обвязку (см. E), не трогая
`tdhf_mrsf_z_vector`. Полная интеграция в драйвер — Фаза 5.

---

## E. Тесты Фазы 2 (deliverable, без passing — нет коммита)

Functional **BHHLYP**, basis **6-31G(d)**, UHF ref, mult=3.

### E.1 FD energy-gradient sanity (база)
- Молекула: **H₂O** (closed-shell, n_O=2 SOMO, малый).
- Скрипт `tests/umrsf_grad/fd_energy_grad.py`: ∂E/∂ξ численно (±δ, δ=1e-3 Bohr,
  central diff) через UMRSF energy на смещённых геометриях. Это эталон для будущих фаз.
- Input: `tests/umrsf_grad/h2o_umrsf_bhhlyp.inp` (energy runtype).

### E.2 Покомпонентная проверка Q_{tu,σ} (запрошено)
- Отдельно блоки **ia, ix, xa** (и спины α/β).
- Аналитический Q (из `umrsfqcal`) vs FD: Q_{tu,σ} ≈ ∂/∂(rotation) от G по соответствующей
  MO-ротации. Реализация: dump Q в JSON через тестовый драйвер; сравнить с конечно-
  разностным откликом G на pairwise MO-rotation угол θ_{tu,σ}. Критерий < 1e-5.
- Файл: `tests/umrsf_grad/q_components.py` + `q_ref_h2o.json`.

### E.3 Check 1 (Приложение A теории) — численно
- На **H₂O** при UKS-решении, близком к ROHF: проверить
  `Q̄^RO_{tu} ?= ½(Q^U_{tu,α} + Q^U_{tu,β})`.
- Берём RO Q из существующего RO-пути (mrst=1/3, `sfrorhs`-ветка) и UMRSF Q из
  `umrsfqcal`; сравнить усреднение поблочно. Критерий < 1e-5 (при совпадении UKS≈ROHF
  орбиталей; иначе документируем расхождение как контролируемое).
- Файл: `tests/umrsf_grad/check1_avg.py`.

### E.4 R = Q−Qᵀ согласованность
- Проверить, что `umrsfqrorhs(rhs, qa, qb)` даёт `rhs` совпадающий с прямым
  Q-差 поблочно; плюс что для UKS≈ROHF суммированный UMRSF R воспроизводит RO `sfrorhs`.

**Критерий всех тестов:** max|Δ| < 1e-5 (Hartree, Hartree/Bohr, или ед. Q).
Без E.1–E.4 passing — коммита Фазы 2 нет.

---

## F. Reuse vs rewrite — сводка

| RO routine | Фаза 2 | Решение |
|---|---|---|
| `umrsfcbc` (build 10-comp SP density) | вызывается | **reuse** (готово) |
| `umrsfsp` (H[X,X] spin-pairing → hxa/hxb) | вызывается | **reuse** (готово) |
| `mntoia` (H⁺[T] AO→MO) | вызывается | **reuse** |
| `int2_umrsf_data_t` (UHF ERI) | вызывается | **reuse** (готово) |
| `sfrorhs` (RO 3-блок R) | — | **rewrite** → `umrsfqrorhs` (4-блок, spin-explicit) |
| `mrsfqrorhs`/`mrsfqroesum` (SF-путь) | — | шаблон раскладки, не вызываем |
| T-сборка (dgemm на bvec_mo_d) | inline в тест-драйвере | reuse паттерн RO 887–896 |

---

## РЕШЕНИЯ (приняты, 2026-06-02)
1. ✅ **Разделять Q и R**: `umrsfqcal` (Q, оба спина) + `umrsfqrorhs` (R=Q−Qᵀ).
2. ✅ **Порядок блоков rhs/Z — α первыми, потом β:**
   ```
   1. Z_xa_α (OV_α): n_O·n_V    x=nocb+1..noca, a=noca+1..nbf
   2. Z_ia_α (CV_α): n_C·n_V    i=1..nocb,      a=noca+1..nbf
   3. Z_ix_β (CO_β): n_C·n_O    i=1..nocb,      x=nocb+1..noca
   4. Z_ia_β (CV_β): n_C·n_V    i=1..nocb,      a=noca+1..nbf
   ```
   (Эту же раскладку обязан читать J-оператор Фазы 3.)
3. ✅ **Тест-обвязка**: отдельный мини-драйвер `tdhf_umrsf_qrtest` (C-binding), дампит
   Qα/Qβ/R в tagarray/файл; production-драйвер не трогаем.

## RO-limit debugging resolution (Phase 2)

The UMRSF vs RO H[X,X] "0.189 Ha discrepancy" was a **measurement artifact, not a bug**:
`sfrorhs` declares `xhxa,xhxb` as `intent(inout)` and accumulates the Fock·T term in
place (`xhxa += 2·Fa·Tij`, `xhxb += 2·Fb·Tab`, tdhf_sf_lib.F90:364-379). The qrtest
dump captured `hxa_r` *after* `sfrorhs`, so it showed H[X,X] + Fock·T while `hxa_u`
showed plain H[X,X].

Proven correct via the mini-driver (h2o_rohf_mrsf, BHHLYP/6-31G*):
- H0 (A0 SF Fock) part:        H0A_U == H0A_R  to 0.0
- spin-pairing:                umrsfsp ≡ mrsfsp, |SPXA_U-SPXA_R| = 8.7e-19
- assembled H[X,X] (pre-sfrorhs): hxa_u == hxa_r exactly.
Also confirmed: do NOT pass an H0-filled array straight into umrsfsp/mrsfsp and rely on
intent(out) accumulation; build H[X,X] = H0(dedicated array) + sp(zeroed array) by
explicit sum.

Remaining Phase-2 work: validate umrsfqcal + umrsfqrorhs (the Q/R builders) against the
golden sfrorhs RHS. The block-norm mismatch (U.COb≈5e-4 vs R.CO≈0.15) is a Q-builder /
block-correspondence question, independent of H[X,X].
