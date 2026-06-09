#!/usr/bin/env python
"""Phase-5 UMRSF Gamma inter-block fold test (reads umrsf_qr_dump.txt).

The Gamma SF and intra parts are structurally identical to RO (fed spin-resolved
P/D/X), so only the INTER spin-pairing block changes: RO's single block ->
0.5*(inter_aa[C^a] + inter_bb[C^b]) per theory 12.3. At alpha=beta (mo_a==mo_b) the
alpha/beta inter densities coincide with RO's, so the fold pins the overall weight:
  (inter_aa + inter_bb) / RO  must be 2  ->  weight 0.5 reproduces the validated RO total.
"""
import os, sys, numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))


def parse_dump(path):
    vals = {}
    with open(path) as f:
        L = f.read().splitlines()
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
    return vals


def inter_cd(c1, c2, v1, v2, nbf, seed=0, N=4000):
    """sum |(-dc1-dc2-dc3-dc4 + dd1+dd2+dd3+dd4)| over N random AO quartets."""
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, nbf, size=(N, 4))
    tot = 0.0
    for i, j, k, l in idx:
        dc1 = c1[i,k]*v2[j,l]+c1[i,l]*v2[j,k]+c1[j,k]*v2[i,l]+c1[j,l]*v2[i,k]+c1[l,j]*v2[k,i]+c1[k,j]*v2[l,i]+c1[l,i]*v2[k,j]+c1[k,i]*v2[l,j]
        dc2 = c2[i,k]*v1[j,l]+c2[i,l]*v1[j,k]+c2[j,k]*v1[i,l]+c2[j,l]*v1[i,k]+c2[l,j]*v1[k,i]+c2[k,j]*v1[l,i]+c2[l,i]*v1[k,j]+c2[k,i]*v1[l,j]
        dc3 = v2[i,k]*c1[j,l]+v2[i,l]*c1[j,k]+v2[j,k]*c1[i,l]+v2[j,l]*c1[i,k]+v2[l,j]*c1[k,i]+v2[k,j]*c1[l,i]+v2[l,i]*c1[k,j]+v2[k,i]*c1[l,j]
        dc4 = v1[i,k]*c2[j,l]+v1[i,l]*c2[j,k]+v1[j,k]*c2[i,l]+v1[j,l]*c2[i,k]+v1[l,j]*c2[k,i]+v1[k,j]*c2[l,i]+v1[l,i]*c2[k,j]+v1[k,i]*c2[l,j]
        dd1 = c1[i,j]*v2[l,k]+c1[i,j]*v2[k,l]+c1[j,i]*v2[l,k]+c1[j,i]*v2[k,l]+c1[l,k]*v2[i,j]+c1[k,l]*v2[i,j]+c1[l,k]*v2[j,i]+c1[k,l]*v2[j,i]
        dd2 = c2[i,j]*v1[l,k]+c2[i,j]*v1[k,l]+c2[j,i]*v1[l,k]+c2[j,i]*v1[k,l]+c2[l,k]*v1[i,j]+c2[k,l]*v1[i,j]+c2[l,k]*v1[j,i]+c2[k,l]*v1[j,i]
        dd3 = v2[i,j]*c1[l,k]+v2[i,j]*c1[k,l]+v2[j,i]*c1[l,k]+v2[j,i]*c1[k,l]+v2[l,k]*c1[i,j]+v2[k,l]*c1[i,j]+v2[l,k]*c1[j,i]+v2[k,l]*c1[j,i]
        dd4 = v1[i,j]*c2[l,k]+v1[i,j]*c2[k,l]+v1[j,i]*c2[l,k]+v1[j,i]*c2[k,l]+v1[l,k]*c2[i,j]+v1[k,l]*c2[i,j]+v1[l,k]*c2[j,i]+v1[k,l]*c2[j,i]
        tot += abs(-dc1 - dc2 - dc3 - dc4 + dd1 + dd2 + dd3 + dd4)
    return tot


def main():
    v = parse_dump(os.path.join(BASE, 'umrsf_qr_dump.txt'))
    nbf = v['BCO1A'].shape[0]
    ok = True

    # alpha=beta identity of the inter densities
    print('alpha=beta inter-density identity (umrsfcbc a/b vs RO mrsfcbc):')
    idmax = 0.0
    for nm in ['BCO1', 'BCO2', 'BO1V', 'BO2V']:
        dA = np.max(np.abs(v[nm + 'A'] - v[nm + 'R']))
        dB = np.max(np.abs(v[nm + 'B'] - v[nm + 'R']))
        idmax = max(idmax, dA, dB)
        print('  %-5s |A-R|=%.2e |B-R|=%.2e' % (nm, dA, dB))
    id_ok = idmax < 1e-12
    print('  -> %s' % ('PASS' if id_ok else 'FAIL'))

    # inter contraction fold: 0.5*(aa+bb) == RO
    RO = inter_cd(v['BCO1R'], v['BCO2R'], v['BO1VR'], v['BO2VR'], nbf)
    UAA = inter_cd(v['BCO1A'], v['BCO2A'], v['BO1VA'], v['BO2VA'], nbf)
    UBB = inter_cd(v['BCO1B'], v['BCO2B'], v['BO1VB'], v['BO2VB'], nbf)
    ratio = (UAA + UBB) / RO
    weighted = 0.5 * (UAA + UBB)
    fold_err = abs(weighted - RO) / RO
    print('\ninter contraction fold (sum|inter_cd| over 4000 random AO quartets):')
    print('  RO=%.6e  aa=%.6e  bb=%.6e' % (RO, UAA, UBB))
    print('  (aa+bb)/RO = %.6f   weight = %.6f' % (ratio, RO / (UAA + UBB)))
    print('  0.5*(aa+bb) vs RO : rel err = %.3e' % fold_err)
    w_ok = abs(ratio - 2.0) < 1e-6 and fold_err < 1e-9

    ok = id_ok and w_ok
    print('\nPHASE-5 GAMMA INTER-FOLD %s (weight 0.5 set in grd2_umrsf get_density)'
          % ('PASSED' if ok else 'FAILED'))
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
