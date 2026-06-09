# UMRSF α≠β deferred-residual review

These code paths are **not** validated by the α=β isolation tests (Phases 2–5). They are
verified against the theory §10/§12 boxes only, and numerically close at the **Phase-6
UHF (α≠β) gradient FD**. This document pairs each path with its theory box so the
derivation can be checked against the source by inspection. (File/line refer to the
current tree; `umrsf_gradient_theory.tex` §10 = "Явные формулы для W", §12 = "Γ".)

Spaces: C=closed `[1,nocb]`, O=SOMO `[nocb+1,noca]`, V=virtual `[noca+1,nbf]`.
za = α rotation density (OV_α,CV_α); zb = β (CO_β,CV_β). e_a/e_b = orbital energies.

---

## R1. CV cross-couplings in the J operator (`umrsf_sfrolhs`)

The fa-vs-fb spin label and sign of the CV↔{CO,OV} couplings. At α=β fa=fb so these are
invisible; symmetry pins the *pairing* but not the spin label.

### CV_α  (block 2), source line 2833
```fortran
do y = nocb+1, noca                 ! + fb(O,V) . z_CO   (couples CV_a -> CO_b)
  w = w + fb(y,a)*zb(i,y)
end do
```
**§10 W_iaα box (eq. CV,α):** `... + Σ_x F_xa,β Z_ix,β`.
→ `F_xa,β = fb(x,a)`, `Z_ix,β = zb(i,x)`. **Check:** code uses `fb` (β Fock) ✓, reads
`zb` (CO_β) ✓. Spin label matches box.

### CV_β  (block 4), source line 2863
```fortran
do y = nocb+1, noca                 ! - fa(O,C) . z_OV   (couples CV_b -> OV_a)
  w = w - fa(y,i)*za(y,a)
end do
```
**§10 W_iaβ box (eq. CV,β):** `... + Σ_x F_ix,α Z_xa,α`.

**CORRECTION (J vs W).** `umrsf_sfrolhs` is the **J operator**, the §10 box is **W**.
From the master eq, J·Z|_iaβ = C_iaβ[Z] − C_aiβ[Z], where C_iaβ is exactly the §10 W-box
coupling (W_iaβ = Q + C_iaβ + Hp[P]). So:
    J·Z|_iaβ (OV) = C_iaβ(OV) − C_aiβ(OV) = (+F_ix,α Z_xa,α) − C_aiβ(OV).
A literal sign flip between the J-code (−fa) and the §10 W-box (+F) is EXACTLY predicted
when C_aiβ(OV) ≠ 0 (e.g. C_aiβ(OV)=+2F ⟹ J=−F). The OV coupling is pure cross-spin
H̃^βα (Coulomb+f^xc_αβ, eq:Htilde, no exchange).

**Numerical (golden RO, α=β):** the code reproduces RO's J operator exactly:
  |JU[OVa,CVa]+JU[OVa,CVb] − JRO[OV,CV]| = 5.2e-15 ;
  |JU[CVa,OVa]+JU[CVb,OVa] − JRO[CV,OV]| = 6.9e-18.
The −fa term is INHERITED verbatim from RO sfrolhs (golden), not a free α=β choice.

**ADJUDICATED NUMERICALLY (mini-driver CAB dump):** built H̃^βα[Z_OV_α] via the shared
int2(A+B)+f^xc, full MO. Result:
  |CAB − CABᵀ| = 1.1e-15 (symmetric) ;  ratio C_aiβ(OV)/C_iaβ(OV) = 1.0000 (all 52 elts).

ratio = 1 (NOT 2, NOT 0) ⟹ the premise was wrong: the §10 +F and the code −fa are NOT the
same term. H̃^βα is symmetric, so C_iaβ = C_aiβ and the cross-spin term **CANCELS in the J
operator** (J·Z|_iaβ ⊃ C_iaβ − C_aiβ = 0). Hence:
  - §10 +F_ix,α Z_xa,α = cross-spin H̃ — REAL in W (umrsf_sfrowcal carries +fa·z_OV in
    W_iaβ, matches the box), VANISHES in the J operator (umrsf_sfrolhs).
  - code −fa·z_OV in the J operator = the SAME-SPIN generalized-Fock coupling (golden RO),
    a different term that does not cancel.
**No sign bug.** J has golden-RO same-spin −fa; W has §10 cross-spin +fa; both correct
because H̃ cancels in J but survives in W. The original §10-vs-code "contradiction" was a
W-vs-J + H̃-cancellation artifact.

**Remaining (unchanged in register):** the SPIN-BLOCK partition of the same-spin −fa·z_OV
coupling (CV_α vs CV_β). The ratio test addresses the cross-spin term only; symmetry +
RO-sum-match leave the same-spin block assignment free at α=β. Closes at Phase-6 UHF FD.

### OV_α reads z_CV as CV_β (block 1), line 2818  /  CO_β reads z_CV as CV_α (block 3), line 2851
```fortran
! OV_a:  ... - fa(j,x)*zb(j,a)          ! z_CV routed to CV_BETA (zb)  [symmetry-required]
! CO_b:  ... + fb(b,x)*za(i,b)          ! z_CV routed to CV_ALPHA (za) [symmetry-required]
```
**Pinned by J symmetry** (swapping these breaks |J−Jᵀ| by ~0.07, caught in Phase-3 test
c). The fa/fb label here is fixed by the transpose partner (CV_α/CV_β above), so it
inherits R1's status.

---

## R2. CV cross-coupling + W_ia partition in the W builder (`umrsf_sfrowcal`)

### W_ia, source lines 2992–2994
```fortran
do x = nocb+1, noca;  w = w + fa(x,i)*za(x,a);  end do   ! + fa(O,C).z_OV  -> wb (beta)
wa(i,a) = e_a(i)*za(i,a)                                  ! CV_a diagonal piece
wb(i,a) = xhxb(i,a) + hs*w + e_b(i)*zb(i,a)               ! CV_b: H[X,X]+F + z_OV coupling + diag
```
**§10 W_iaα box:** `Hp_iaα[P] + (ε_i − ε_a) Z_iaα + Σ_x F_xa,β Z_ix,β`.
**§10 W_iaβ box:** `2 H_aiβ + Hp_iaβ[P] + (ε̄_i − ε̄_a) Z_iaβ + Σ_x F_ix,α Z_xa,α`.

**Three sub-items NOT closed at α=β (the W_ia 0.203 fold residual = exactly −e_i·z_CV):**
1. **CV diagonal partition:** code puts `e_a(i)·za` in wa and `e_b(i)·zb` in wb (the
   per-spin diagonals). RO `mrsfrowcal` spin-sums one `e_i·z_CV` → fold doubles by
   `e_i·z_CV`. Which is right at α≠β: per-spin (code) is physically expected, but the
   `(ε_i − ε_a)` form means the **ε_a (virtual-energy) term is MISSING** from both wa and
   wb (RO `mrsfrowcal` only carries `e_i`, no `e_a` — inherited here).
2. **occ-virt H⁺[P]:** §10 W_iaα/β carry `Hp_iaα[P]`/`Hp_iaβ[P]` (occ-virt). RO
   `mrsfrowcal` computes H⁺[P] **occ-occ only** (ppija/ppijb are nocc×nocc); the occ-virt
   piece is **absent** in the code (inherited from RO). *Confirm whether §10 occ-virt
   Hp[P] is real or cancels.*
3. **fa-vs-fb in the z_OV coupling:** code uses `fa(x,i)` (→ wb) per §10 W_iaβ
   `Σ_x F_ix,α Z_xa,α` (α Fock). Matches box; α=β-invisible.

---

## R3. Γ inter-block αα/ββ split + weight (`grd2_umrsf_..._get_density`, Phase 5 — to write)

RO single inter block (mrsf_gradient.F90:405–480) uses one `{bco1,bco2,bo1v,bo2v}` from
`mo_a`:
```
df1 += sgnk·qfspcp3·(−dc1−dc2−dc3−dc4 + dd1+dd2+dd3+dd4)
```
**§12 eq (928):** `Γ^SP = sgn·c_HF·[ δ_σα δ_τβ Γ^intra,αβ + δ_σα δ_τα Γ^inter,αα + δ_σβ δ_τβ Γ^inter,ββ ]`.
UMRSF must evaluate the dc/dd inter terms **twice** (αα from `mo_a` densities, ββ from
`mo_b`), per §12.3 ("inter,ββ … was absent in RO … this is the new UMRSF-specific term").
The **weight** `qfspcp3` for the split (whether RO's single block = αα+ββ, or needs ½) is
**α=β-resolvable** in the Phase-5 fold test — not deferred.

---

## Summary of what each channel validates

| Item | α=β test status | Closes at |
|---|---|---|
| R1 CV fa/fb label | symmetry-pinned; **sign vs §10 open** | Phase-6 FD |
| R1 OV/CO z_CV routing | J-symmetry pinned | (pinned) |
| R2.1 CV diagonal partition + missing ε_a | fold residual = −e_i·z_CV (confirmed) | Phase-6 FD |
| R2.2 occ-virt H⁺[P] | absent (RO-inherited) | Phase-6 FD |
| R2.3 z_OV fa-label | matches §10 | (matches) |
| R3 inter αα/ββ split | structural (Phase 5) | Phase-5 fold (weight) + Phase-6 FD |

**One open theory question to resolve before Phase 6 (R1, CV_β sign):** the §10 W_iaβ box
sign (`+Σ F_ix,α Z_xa,α`) is opposite to the RO-inherited code sign (`−fa(y,i)·za`). These
disagree; at α=β both give the same RO-matching J. Need the correct α≠β sign.
