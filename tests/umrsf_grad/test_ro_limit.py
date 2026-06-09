#!/usr/bin/env python
"""Phase-2 UMRSF Q/R builder validation (mini-driver tdhf_umrsf_qrtest).

Validation strategy (see phase2_plan.md "RO-limit debugging resolution" and the
memory note umrsf-ro-limit-test-invalid):

  PASS criterion (asserted): with an ROHF reference (mo_a == mo_b) the UMRSF
  H[X,X] operator (A0 spin-flip Fock + spin-pairing, built by umrsfdmat/umrsfcbc/
  int2_umrsf/umrsfsp) must reduce *exactly* to the golden RO H[X,X] (mrsf path),
  for both alpha and beta.  This is the hard, genuinely-new machinery and it is
  parametrization-independent.  Tolerance 1e-10.

  REPORTED (not asserted): the umrsf 4-block Z-vector RHS (umrsfqcal+umrsfqrorhs)
  vs the golden RO sfrorhs RHS.  These do NOT correspond block-by-block: UMRSF and
  RO-MRSF parametrize the orbital response differently (independent alpha/beta
  rotations vs one spin-adapted set), so sfrorhs's CO block carries a lone
  non-antisymmetrized spin-flip term (-xhxb(j,i)) that no Q-Q^T folding reproduces.
  We print the per-block comparison for diagnostics; final numerical validation of
  Q/R is deferred to the finite-difference full-gradient test once Phases 3-6 exist.
"""
import os, sys
# IMPORTANT: import oqp before numpy (ILP64 vs LP64 BLAS clash, see project memory).
import oqp
from oqp.pyoqp import Runner
import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
TOL = 1e-10


def run_energy_and_qrtest(inp):
    project = os.path.basename(inp).replace('.inp', '')
    log = os.path.join(BASE, project + '.log')
    runner = Runner(project=project, input_file=inp, log=log, usempi=False)
    runner.run()
    oqp.tdhf_umrsf_qrtest(runner.mol)


def parse_dump(path):
    vals = {}
    with open(path) as f:
        lines = f.read().splitlines()
    hdr = [int(x) for x in lines[0].split() if x.lstrip('-').isdigit()]
    header = dict(nbf=hdr[0], nocca=hdr[1], noccb=hdr[2],
                  lzdim_u=hdr[3], lzdim_ro=hdr[4])
    i = 1
    while i < len(lines):
        t = lines[i].split()
        if t and t[0] == '@MAT':
            name, n1, n2 = t[1], int(t[2]), int(t[3])
            d = np.array([float(lines[i + 1 + k]) for k in range(n1 * n2)])
            vals[name] = d.reshape(n2, n1).T
            i += 1 + n1 * n2
        elif t and t[0] == '@VEC':
            name, m = t[1], int(t[2])
            d = np.array([float(lines[i + 1 + k]) for k in range(m)])
            vals[name] = d
            i += 1 + m
        else:
            i += 1
    return header, vals


def main():
    inp = os.path.join(BASE, 'h2o_rohf_mrsf.inp')
    run_energy_and_qrtest(inp)
    hdr, v = parse_dump(os.path.join(BASE, 'umrsf_qr_dump.txt'))
    nbf, noca, nocb = hdr['nbf'], hdr['nocca'], hdr['noccb']
    print('dims: nbf=%d nocca=%d noccb=%d' % (nbf, noca, nocb))

    # ---- ASSERTED: UMRSF H[X,X] reduces to RO golden H[X,X] ----------------
    da = float(np.max(np.abs(v['HXA_U'] - v['HXA_R'])))
    db = float(np.max(np.abs(v['HXB_U'] - v['HXB_R'])))
    print('\nH[X,X] RO-reduction (A0 SF Fock + spin-pairing):')
    print('  |HXA_U - HXA_R|max = %.3e' % da)
    print('  |HXB_U - HXB_R|max = %.3e' % db)
    ok = (da < TOL) and (db < TOL)
    print('  -> %s (tol %.0e)' % ('PASS' if ok else 'FAIL', TOL))

    # ---- REPORTED: Q/R vs sfrorhs (not a block-by-block equivalence) -------
    C = list(range(0, nocb)); O = list(range(nocb, noca)); V = list(range(noca, nbf))
    qa, qb = v['QA'], v['QB']
    nsocc, nvira = noca - nocb, nbf - noca
    ro = v['RHS_RO']; o = 0
    CO = ro[o:o + nocb * nsocc].reshape(nsocc, nocb); o += nocb * nsocc
    CV = ro[o:o + nocb * nvira].reshape(nvira, nocb); o += nocb * nvira
    OV = ro[o:o + nsocc * nvira].reshape(nvira, nsocc); o += nsocc * nvira
    Qs = qa + qb

    def fold(rows, cols):
        return np.array([[-(Qs[p, q] - Qs[q, p]) for q in cols] for p in rows])

    def rep(name, pred, ref):
        d = min(np.max(np.abs(s * pred - ref)) for s in (1, -1)) if ref.size else 0.0
        print('  %-3s -(qa+qb) fold: maxdiff=%.3e  ||pred||=%.3e ||ref||=%.3e'
              % (name, d, np.linalg.norm(pred), np.linalg.norm(ref)))

    print('\nQ/R vs golden sfrorhs (diagnostic only, see docstring):')
    rep('CV', fold(V, C), CV)
    rep('OV', fold(V, O), OV)
    rep('CO', fold(O, C), CO)

    print('\nPHASE-2 H[X,X] CHECK %s' % ('PASSED' if ok else 'FAILED'))
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
