#!/usr/bin/env python
"""Phase-4 UMRSF P=T+Z and W-builder isolation tests (reads umrsf_qr_dump.txt).

  P:  controlled-Z P == RO sfropcal P, alpha and beta (exact placement check).
  W:  spin-resolved fold W^alpha + W^beta == RO mrsfrowcal, per block. All density-
      partitionable blocks must fold exactly. The W_ia (CV) block carries the CV
      diagonal-collapse residual (RO spin-sums one e_i.z_CV; the spin-resolved scheme
      has a diagonal per spin) -- we VERIFY the residual equals exactly that known
      term, confirming no spurious error. The CV partition closes at the Phase-6 UHF
      gradient FD (deferred residual).
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
    C = list(range(nocb)); O = list(range(nocb, noca)); Vv = list(range(noca, nbf))
    EA = v['EA']; zro = v['ZRO']
    ok = True

    # ---- P = T+Z ----------------------------------------------------------
    dpa = np.max(np.abs(v['PAUC'] - v['PARO']))
    dpb = np.max(np.abs(v['PBUC'] - v['PBRO']))
    p_ok = dpa < 1e-12 and dpb < 1e-12
    print('P = T+Z (controlled Z) vs RO sfropcal:')
    print('  |PAUC-PARO|=%.3e  |PBUC-PBRO|=%.3e -> %s' % (dpa, dpb, 'PASS' if p_ok else 'FAIL'))
    ok &= p_ok

    # ---- W fold -----------------------------------------------------------
    fold = v['WA_UC'] + v['WB_UC']
    d = fold - v['WMO_RO']

    def blk(M, r, c): return M[np.ix_(r, c)]
    blocks = [('W_ij CC', C, C), ('W_xy OO', O, O), ('W_ab VV', Vv, Vv),
              ('W_ix CO', C, O), ('W_xa OV', O, Vv),
              ('W_xi OC', O, C), ('W_ai VC', Vv, C), ('W_ax VO', Vv, O)]
    print('\nW fold (W^a + W^b == RO mrsfrowcal), per block:')
    wfold_ok = True
    for nm, r, c in blocks:
        m = np.max(np.abs(blk(d, r, c)))
        st = 'OK' if m < 1e-9 else 'FAIL'
        if m >= 1e-9: wfold_ok = False
        print('  %-9s %.3e  %s' % (nm, m, st))

    # W_ia (CV): residual must equal exactly -e_i * z_CV (the CV diagonal collapse).
    # zro CV block order: a outer, i inner (RO doc-virt). Build predicted residual.
    o1 = nocb * nO  # CO size
    zCV = zro[o1:o1 + nocb * nV]                      # RO CV vector (a outer, i inner)
    pred = np.zeros((nC, nV))
    kk = 0
    for ia in range(nV):       # a outer
        for ii in range(nC):   # i inner
            pred[ii, ia] = -EA[C[ii]] * zCV[kk]
            kk += 1
    resid_ia = blk(d, C, Vv)
    cv_match = np.max(np.abs(resid_ia - pred))
    cv_ok = cv_match < 1e-9
    print('  W_ia CV   %.3e  (residual == -e_i.z_CV collapse to %.3e -> %s; partition deferred to Phase-6)'
          % (np.max(np.abs(resid_ia)), cv_match, 'CONFIRMED' if cv_ok else 'UNEXPECTED'))
    ok &= (p_ok and wfold_ok and cv_ok)

    print('\nPHASE-4 P/W TESTS %s' % ('PASSED' if ok else 'FAILED'))
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
