#!/usr/bin/env python
"""Phase-3 UMRSF Z-vector J-operator isolation tests (reads umrsf_qr_dump.txt).

Run test_ro_limit.py first (it builds the dump incl. dense JU/JRO). Tests:
  (b) block-order: J diagonal correlates with the (eps_v-eps_o) pattern implied by
      the umrsf_sfrogen index map (catches a mis-ordered block).
  (c) symmetry: |JU - JU^T| (orbital Hessian must be symmetric). JRO too (sanity of
      the dense-build method against the golden RO operator).
  (d') alpha=beta vs RO operator: OV_a block == RO OV, CO_b block == RO CO, and the
      CV response-sum identity pins the resp_scale (0.5) normalization.
  (f) dense solve JU.Z = -R: residual + spectrum (non-singular).
"""
import os, sys, numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))


def parse_dump(path):
    vals = {}
    with open(path) as f:
        L = f.read().splitlines()
    h = [int(x) for x in L[0].split() if x.lstrip('-').isdigit()]
    hdr = dict(nbf=h[0], nocca=h[1], noccb=h[2], lzdim_u=h[3], lzdim_ro=h[4])
    i = 1
    while i < len(L):
        t = L[i].split()
        if t and t[0] == '@MAT':
            n, a, b = t[1], int(t[2]), int(t[3])
            d = np.array([float(L[i + 1 + k]) for k in range(a * b)])
            vals[n] = d.reshape(b, a).T
            i += 1 + a * b
        elif t and t[0] == '@VEC':
            n, m = t[1], int(t[2])
            vals[n] = np.array([float(L[i + 1 + k]) for k in range(m)])
            i += 1 + m
        else:
            i += 1
    return hdr, vals


def main():
    hdr, v = parse_dump(os.path.join(BASE, 'umrsf_qr_dump.txt'))
    nbf, noca, nocb = hdr['nbf'], hdr['nocca'], hdr['noccb']
    nO, nV, nC = noca - nocb, nbf - noca, nocb
    JU, JRO = v['JU'], v['JRO']
    fa, fb = v['FA'], v['FB']
    EA = v['EA']
    rhs_u = v['RHS_U']
    C = list(range(0, nocb)); O = list(range(nocb, noca)); V = list(range(noca, nbf))
    print('dims nbf=%d C=%d O=%d V=%d  lzdim_u=%d lzdim_ro=%d'
          % (nbf, nC, nO, nV, JU.shape[0], JRO.shape[0]))

    # ---- index maps -------------------------------------------------------
    # umrsf 4-block (a outer / x|i inner), order OVa, CVa, COb, CVb
    umap = []  # (o_idx, v_idx, spin)
    for a in V:
        for x in O: umap.append((x, a, 'a'))   # OV_a
    for a in V:
        for i in C: umap.append((i, a, 'a'))   # CV_a
    for x in O:
        for i in C: umap.append((i, x, 'b'))   # CO_b
    for a in V:
        for i in C: umap.append((i, a, 'b'))   # CV_b
    OVa = slice(0, nO * nV); CVa = slice(nO * nV, nO * nV + nC * nV)
    COb = slice(nO * nV + nC * nV, nO * nV + nC * nV + nC * nO)
    CVb = slice(nO * nV + nC * nV + nC * nO, nO * nV + 2 * nC * nV + nC * nO)
    # RO 3-block order CO, CV, OV
    COr = slice(0, nC * nO); CVr = slice(nC * nO, nC * nO + nC * nV)
    OVr = slice(nC * nO + nC * nV, nC * nO + nC * nV + nO * nV)

    ok = True

    # ---- (b) block-order: diagonal vs (eps_v - eps_o) ---------------------
    ediff = np.array([(fa if s == 'a' else fb)[vv, vv] - (fa if s == 'a' else fb)[oo, oo]
                      for (oo, vv, s) in umap])
    jdiag = np.diag(JU)
    # response makes jdiag != ediff exactly, but they must correlate strongly and
    # share sign on every element if the block map is right.
    # Decisive metric: every diagonal sign must equal sign(eps_v-eps_o) under the
    # umrsf_sfrogen index map (impossible if a block is mis-ordered). corr<1 only
    # because diag(JU) = (eps_v-eps_o) + 0.5*hpz_diag (the response adds scatter).
    corr = np.corrcoef(ediff, jdiag)[0, 1]
    signmatch = np.mean(np.sign(ediff) == np.sign(jdiag))
    b_ok = (signmatch == 1.0) and corr > 0.95
    print('\n(b) block-order: sign-match=%.3f (all %d diag signs)  corr=%.4f -> %s'
          % (signmatch, len(ediff), corr, 'PASS' if b_ok else 'FAIL'))
    ok &= b_ok

    # ---- (c) symmetry (JRO sets the int2 noise floor of the dense build) --
    asu = np.max(np.abs(JU - JU.T)); asr = np.max(np.abs(JRO - JRO.T))
    c_ok = asu < max(1e-6, 10 * asr)
    print('(c) symmetry: |JU-JU^T|max=%.3e  |JRO-JRO^T|max=%.3e (noise floor) -> %s'
          % (asu, asr, 'PASS' if c_ok else 'FAIL'))
    ok &= c_ok

    # ---- (d') alpha=beta vs RO -------------------------------------------
    dOV = np.max(np.abs(JU[OVa, OVa] - JRO[OVr, OVr]))
    dCO = np.max(np.abs(JU[COb, COb] - JRO[COr, COr]))
    # CV response-sum identity: with Z_CVa=Z_CVb the perturbation density equals RO's
    # z_CV, so the shared int2 gives identical hpza/hpzb. Then
    #   sum_4(umrsf CV blocks) - JRO[CV,CV] = D + (alpha-beta cross response diagonal),
    # i.e. a PURELY DIAGONAL matrix (the off-diagonal response must cancel exactly).
    # The off-diagonal cancellation pins the 0.5 (A+B) normalization of the split CV.
    sum4 = JU[CVa, CVa] + JU[CVa, CVb] + JU[CVb, CVa] + JU[CVb, CVb]
    resid = sum4 - JRO[CVr, CVr]
    cv_offmax = np.max(np.abs(resid - np.diag(np.diag(resid))))
    # informational: diagonal carries D (eps_a-eps_i) + the physical ab cross-response
    Dcv = np.array([EA[a] - EA[i] for a in V for i in C])
    cv_diag_minus_D = np.median(np.diag(resid) - Dcv)
    d_ok = dOV < 1e-9 and dCO < 1e-9 and cv_offmax < 1e-6
    print("(d') OV_a==RO OV: %.3e | CO_b==RO CO: %.3e" % (dOV, dCO))
    print("     CV resp-sum off-diag(sum4-JRO_CV)=%.3e (pins 0.5 norm); "
          "median(diag-D)=%.3e [ab cross-resp, informational] -> %s"
          % (cv_offmax, cv_diag_minus_D, 'PASS' if d_ok else 'FAIL'))
    ok &= d_ok

    # ---- (f) dense solve + spectrum --------------------------------------
    # Z-vector eq J.Z = -R ; rhs_u stores -R already, so solve JU.Z = rhs_u.
    Z = np.linalg.solve(JU, rhs_u)
    res = np.max(np.abs(JU @ Z - rhs_u))
    ev = np.linalg.eigvals(JU)
    cond = np.linalg.cond(JU)
    f_ok = res < 1e-10
    print('(f) dense solve: residual=%.3e  cond(JU)=%.3e  Re(ev) in [%.3f,%.3f] -> %s'
          % (res, cond, ev.real.min(), ev.real.max(), 'PASS' if f_ok else 'FAIL'))
    ok &= f_ok

    print('\nPHASE-3 J-OPERATOR TESTS %s' % ('PASSED' if ok else 'FAILED'))
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
