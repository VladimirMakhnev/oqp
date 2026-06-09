#!/usr/bin/env python3
"""P1.1 prototype + L1 self-check — the 2e spin-spin (SS) dipolar integral for an (s,s,s,s) quartet.

SSC/ZFS project, branch `ssc-zfs`. See CLAUDE.md (§4, §6 L1, §7) and PROGRESS.md (P1.1-P1.3).

WHAT THIS IS (and is NOT)
-------------------------
This is the standalone numerical PROTOTYPE that pins down the closed form of the rank-2 dipolar
SS integral for the simplest shell quartet, four primitive s-Gaussians. It is a SELF-CONTAINED
sanity oracle (numpy + math only; no `oqp` import). It establishes that our closed form is
algebraically correct by triangulating three INDEPENDENT computations of the same integral.

It is *NOT* the L1 stage gate (P1.3). The real L1 gate compares the native integral routine
against finite differences of OpenQP's *actual* ERI engine. That remains pending. Do not read a
pass here as an L1 pass.

THE PHYSICS
-----------
The SS dipolar operator is the rank-2 traceless kernel

    T_kl(r12) = (3 r12,k r12,l - delta_kl r12^2) / r12^5 ,   r12 = r1 - r2.

It is the *traceless part* of the Hessian of 1/r12. As a distribution the bare Hessian carries an
isotropic contact term,

    d_k d_l (1/r) = (3 r_k r_l - delta_kl r^2)/r^5  -  (4 pi / 3) delta_kl delta^3(r) .         (*)

So two related AO integrals appear, with O = overlap of the two charge clouds = <rho1|rho2>:

    H_kl   ("bare-Hessian integral") = integral of  d_k d_l (1/r12)    [trace = -4 pi O != 0]
    S_kl   ("dipolar integral")      = integral of  T_kl               [trace = 0]
    S_kl   = H_kl - (1/3) Tr(H) delta_kl = H_kl + (4 pi/3) delta_kl O   (S = traceless part of H)

S_kl is the physical SS integral entering the ZFS working equation (Sinnecker-Neese Eq. 9): the
ZFS D-tensor is traceless by construction, so the isotropic contact term does NOT contribute and
is removed. Computationally we form the bare Hessian H by any of three independent routes, then
project out its trace.

IMPORTANT (corrected vs the original prototype derivation note): the Gaussian-transform
t-quadrature reproduces the *bare Hessian* H (the contact term IS captured once the operator is
smeared against the charge clouds), exactly as the finite-difference-of-ERI route does. Both equal
the closed-form H. The identity Tr(H) = -4 pi O is verified directly.

CLOSED FORM (unnormalized primitive s-Gaussians)
------------------------------------------------
Electron 1 clouds: centers A,B exponents a,b -> p=a+b, P=(aA+bB)/p, K_AB=exp(-ab/p |A-B|^2).
Electron 2 clouds: centers C,D exponents c,d -> q=c+d, Q=(cC+dD)/q, K_CD=exp(-cd/q |C-D|^2).
rho = p q/(p+q),  R = P - Q,  T = rho |R|^2,  pref = 2 pi^{5/2} / (p q sqrt(p+q)),  K = K_AB K_CD.

    J(displacement d) = pref * K * F0(rho |R - d|^2)                 (Coulomb ERI, operator shifted)
    H_kl = d^2 J / d d_k d d_l |_{d=0} = pref*K*[ 4 rho^2 R_k R_l F2(T) - 2 rho delta_kl F1(T) ]
    O    = K * (pi/(p+q))^{3/2} * exp(-T)
    S_kl = H_kl + (4 pi/3) delta_kl O

THREE INDEPENDENT ORACLES (the triangulation)
---------------------------------------------
  (1) closed form          : H via Boys F1,F2 (convergent Kummer series); S = traceless(H).
  (2) Richardson FD of J(d): numerically differentiate the Coulomb ERI twice; uses ONLY an
                             erf-based F0 (independent Boys path). -> must reproduce H.
  (3) t-quadrature         : H via the Gaussian transform 1/r=(2/sqrt(pi)) int_0^inf e^{-t^2 r^2} dt,
                             evaluated by Gauss-Legendre. No Boys functions at all. -> must
                             reproduce H.

Then: Tr(H) == -4 pi O (contact identity), and Tr(S) == 0 by construction.
Agreement target: 6-8 significant figures (CLAUDE.md §6 L1). Trace(S) = 0 to ~1e-12.
"""

import math
import sys

import numpy as np

# --- L1 self-check tolerances (NOT the gate tolerance; this is vs an analytic oracle) ---
FD_RTOL = 1e-7        # closed-form H vs Richardson-FD of the Coulomb ERI
QUAD_RTOL = 1e-7      # closed-form S vs Boys-free t-quadrature
TRACE_ATOL = 1e-11    # |Tr S| analytic zero


# ---------------------------------------------------------------------------
# Boys functions
# ---------------------------------------------------------------------------
def boys_series(n: int, T: float, tol: float = 1e-16, kmax: int = 500) -> float:
    """F_n(T) = int_0^1 t^{2n} e^{-T t^2} dt via the all-positive Kummer series

        F_n(T) = e^{-T} * sum_{k>=0} (2n-1)!! / (2n+2k+1)!! * (2T)^k .

    No cancellation, converges for all T >= 0. Used by the CLOSED FORM (F1, F2).
    """
    # term_0 = (2n-1)!!/(2n+1)!! = 1/(2n+1)
    term = 1.0 / (2 * n + 1)
    total = term
    twoT = 2.0 * T
    k = 0
    while k < kmax:
        # term_{k+1}/term_k = (2T) / (2n + 2k + 3)
        term *= twoT / (2 * n + 2 * k + 3)
        total += term
        if term <= tol * total:
            break
        k += 1
    return math.exp(-T) * total


def boys_erf_F0(T: float) -> float:
    """F_0(T) = (1/2) sqrt(pi/T) erf(sqrt(T)), with the T->0 limit. INDEPENDENT of boys_series;
    used only by the finite-difference oracle for J(d)."""
    if T < 1e-14:
        return 1.0 - T / 3.0  # series, avoids 0/0
    sT = math.sqrt(T)
    return 0.5 * math.sqrt(math.pi / T) * math.erf(sT)


# ---------------------------------------------------------------------------
# Charge-cloud reduction
# ---------------------------------------------------------------------------
def _cloud(Z1, z1, Z2, z2):
    """Two primitive s-Gaussians -> (total exponent, product center, prefactor K)."""
    Z1 = np.asarray(Z1, float)
    Z2 = np.asarray(Z2, float)
    s = z1 + z2
    C = (z1 * Z1 + z2 * Z2) / s
    K = math.exp(-z1 * z2 / s * float(np.dot(Z1 - Z2, Z1 - Z2)))
    return s, C, K


def _reduce(A, a, B, b, C, c, D, d):
    p, P, K_AB = _cloud(A, a, B, b)
    q, Q, K_CD = _cloud(C, c, D, d)
    rho = p * q / (p + q)
    R = P - Q
    T = rho * float(np.dot(R, R))
    pref = 2.0 * math.pi ** 2.5 / (p * q * math.sqrt(p + q))
    K = K_AB * K_CD
    return dict(p=p, q=q, P=P, Q=Q, rho=rho, R=R, T=T, pref=pref, K=K)


# ---------------------------------------------------------------------------
# (1) Closed form
# ---------------------------------------------------------------------------
def closed_form_H(A, a, B, b, C, c, D, d):
    """Bare-Hessian integral H_kl (trace = -4 pi O)."""
    r = _reduce(A, a, B, b, C, c, D, d)
    rho, R, T, pref, K = r["rho"], r["R"], r["T"], r["pref"], r["K"]
    F1 = boys_series(1, T)
    F2 = boys_series(2, T)
    H = pref * K * (4.0 * rho ** 2 * np.outer(R, R) * F2 - 2.0 * rho * np.eye(3) * F1)
    return H


def cloud_overlap(A, a, B, b, C, c, D, d):
    """O = <rho1|rho2> = K (pi/(p+q))^{3/2} e^{-T}."""
    r = _reduce(A, a, B, b, C, c, D, d)
    return r["K"] * (math.pi / (r["p"] + r["q"])) ** 1.5 * math.exp(-r["T"])


def dipolar_S(A, a, B, b, C, c, D, d):
    """Physical traceless dipolar integral S_kl = H_kl + (4 pi/3) delta_kl O."""
    H = closed_form_H(A, a, B, b, C, c, D, d)
    O = cloud_overlap(A, a, B, b, C, c, D, d)
    return H + (4.0 * math.pi / 3.0) * np.eye(3) * O


# ---------------------------------------------------------------------------
# (2) Richardson finite difference of the Coulomb ERI J(d)
# ---------------------------------------------------------------------------
def _J(A, a, B, b, C, c, D, d, disp):
    """Coulomb ERI with the operator argument displaced by `disp`:
    J(disp) = pref K F0(rho |R - disp|^2). Uses the erf-based F0 (independent Boys)."""
    r = _reduce(A, a, B, b, C, c, D, d)
    Rd = r["R"] - np.asarray(disp, float)
    Td = r["rho"] * float(np.dot(Rd, Rd))
    return r["pref"] * r["K"] * boys_erf_F0(Td)


def _hessian_fd(A, a, B, b, C, c, D, d, h):
    """Central-difference Hessian of J(disp) at disp=0, step h. O(h^2)."""
    e = np.eye(3)
    Hfd = np.zeros((3, 3))
    J0 = _J(A, a, B, b, C, c, D, d, np.zeros(3))
    for k in range(3):
        Jp = _J(A, a, B, b, C, c, D, d, h * e[k])
        Jm = _J(A, a, B, b, C, c, D, d, -h * e[k])
        Hfd[k, k] = (Jp - 2.0 * J0 + Jm) / h ** 2
    for k in range(3):
        for l in range(k + 1, 3):
            Jpp = _J(A, a, B, b, C, c, D, d, h * (e[k] + e[l]))
            Jpm = _J(A, a, B, b, C, c, D, d, h * (e[k] - e[l]))
            Jmp = _J(A, a, B, b, C, c, D, d, h * (-e[k] + e[l]))
            Jmm = _J(A, a, B, b, C, c, D, d, h * (-e[k] - e[l]))
            val = (Jpp - Jpm - Jmp + Jmm) / (4.0 * h ** 2)
            Hfd[k, l] = Hfd[l, k] = val
    return Hfd


def fd_hessian_J(A, a, B, b, C, c, D, d, h=2e-2):
    """Richardson-extrapolated Hessian of J(d): combine steps h and h/2 to kill the O(h^2) error.
    Returns the bare-Hessian integral H (NOT the traceless S)."""
    Hh = _hessian_fd(A, a, B, b, C, c, D, d, h)
    Hh2 = _hessian_fd(A, a, B, b, C, c, D, d, h / 2.0)
    return (4.0 * Hh2 - Hh) / 3.0


# ---------------------------------------------------------------------------
# (3) Boys-free t-quadrature of the bare-Hessian integral H
# ---------------------------------------------------------------------------
def tquad_H(A, a, B, b, C, c, D, d, n=2000):
    """H_kl = (2/sqrt(pi)) int_0^inf [ 4 t^4 M_kl(t) - 2 t^2 delta_kl M0(t) ] dt,
    with the Gaussian-operator moments
        M0(t)  = K (pi/(p+q))^{3/2} (pi/(rho+t^2))^{3/2} exp(-rho t^2 R^2/(rho+t^2)) ,
        M_kl(t)= M0(t) [ u0_k u0_l + delta_kl/(2(rho+t^2)) ] ,  u0 = rho R/(rho+t^2).
    Gauss-Legendre on the substitution t = s/(1-s), s in [0,1). No Boys functions.
    Reproduces the *bare Hessian* H (contact term included); the traceless dipolar S is then the
    trace-projected part of H."""
    r = _reduce(A, a, B, b, C, c, D, d)
    p, q, rho, R, K = r["p"], r["q"], r["rho"], r["R"], r["K"]
    R2 = float(np.dot(R, R))
    base = K * (math.pi / (p + q)) ** 1.5

    s, w = np.polynomial.legendre.leggauss(n)        # nodes on [-1,1]
    s = 0.5 * (s + 1.0)                              # -> [0,1]
    w = 0.5 * w
    t = s / (1.0 - s)
    dtds = 1.0 / (1.0 - s) ** 2                      # dt/ds

    g = rho + t ** 2
    M0 = base * (math.pi / g) ** 1.5 * np.exp(-rho * t ** 2 * R2 / g)
    u0 = rho / g                                     # scalar; u0 vector = u0 * R
    coef = 2.0 / math.sqrt(math.pi)

    H = np.zeros((3, 3))
    eye = np.eye(3)
    for k in range(3):
        for l in range(3):
            Mkl = M0 * (u0 ** 2 * R[k] * R[l] + (eye[k, l] / (2.0 * g)))
            integrand = 4.0 * t ** 4 * Mkl - 2.0 * t ** 2 * eye[k, l] * M0
            H[k, l] = coef * float(np.sum(w * integrand * dtds))
    return H


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------
# One non-degenerate (s,s,s,s) quartet: general-position centers + distinct exponents so that
# all 6 independent tensor components are nonzero. (Bohr, a.u.)
QUARTET = dict(
    A=[0.0, 0.0, 0.0], a=1.20,
    B=[0.5, 0.1, -0.2], b=0.80,
    C=[0.3, -0.4, 0.7], c=1.50,
    D=[-0.2, 0.6, 0.25], d=0.90,
)


def _fmt(M):
    return "\n".join("  [" + "  ".join(f"{x:+.10e}" for x in row) + "]" for row in M)


def run_checks(quartet=QUARTET):
    """Run the three oracles and return a results dict (reusable by the pytest skeleton)."""
    H_cf = closed_form_H(**quartet)
    S_cf = dipolar_S(**quartet)
    O = cloud_overlap(**quartet)
    H_fd = fd_hessian_J(**quartet)
    H_q1 = tquad_H(**quartet, n=1000)
    H_q2 = tquad_H(**quartet, n=2000)

    scale_H = max(np.max(np.abs(H_cf)), 1.0)
    return dict(
        H_cf=H_cf, S_cf=S_cf, O=O, H_fd=H_fd, H_q=H_q2,
        fd_maxdiff=float(np.max(np.abs(H_cf - H_fd))),
        fd_reldiff=float(np.max(np.abs(H_cf - H_fd)) / scale_H),
        quad_maxdiff=float(np.max(np.abs(H_cf - H_q2))),
        quad_reldiff=float(np.max(np.abs(H_cf - H_q2)) / scale_H),
        quad_selfconv=float(np.max(np.abs(H_q1 - H_q2))),
        trace_S=float(np.trace(S_cf)),
        trace_H=float(np.trace(H_cf)),
        contact=4.0 * math.pi * O,   # expected -trace(H)
    )


def main():
    res = run_checks()
    print("=" * 78)
    print("SSC P1.1 prototype — (s,s,s,s) dipolar SS integral, three-way triangulation")
    print("=" * 78)
    print("\nQuartet (a.u.):")
    for key in ("A", "a", "B", "b", "C", "c", "D", "d"):
        print(f"  {key} = {QUARTET[key]}")

    print("\n[1] Closed-form bare-Hessian integral H_kl (Boys F1,F2):")
    print(_fmt(res["H_cf"]))
    print("\n[2] Richardson-FD of Coulomb ERI J(d)  (independent erf-based F0):")
    print(_fmt(res["H_fd"]))
    print("\n[3] Boys-free t-quadrature H_kl (Gauss-Legendre):")
    print(_fmt(res["H_q"]))
    print("\n[1] Physical traceless dipolar integral S_kl = H - (1/3)Tr(H) I = H + (4pi/3) O I:")
    print(_fmt(res["S_cf"]))

    print("\n" + "-" * 78)
    print("AGREEMENT (target: 6-8 significant figures, i.e. rel diff <~ 1e-7):")
    print(f"  closed-form H  vs  Richardson-FD of ERI  : max|d|={res['fd_maxdiff']:.3e}  "
          f"rel={res['fd_reldiff']:.3e}")
    print(f"  closed-form H  vs  Boys-free t-quadrature: max|d|={res['quad_maxdiff']:.3e}  "
          f"rel={res['quad_reldiff']:.3e}")
    print(f"  t-quad self-convergence (n=1000 vs 2000) : max|d|={res['quad_selfconv']:.3e}")
    print("\nTRACE / contact-term consistency:")
    print(f"  Tr(S) [traceless dipolar, target 0]      : {res['trace_S']:+.3e}")
    print(f"  Tr(H) [bare Hessian]                     : {res['trace_H']:+.6e}")
    print(f"  -4*pi*<rho1|rho2> (expected Tr(H))       : {-res['contact']:+.6e}")
    print("-" * 78)

    ok_fd = res["fd_reldiff"] <= FD_RTOL
    ok_quad = res["quad_reldiff"] <= QUAD_RTOL
    ok_trace = abs(res["trace_S"]) <= TRACE_ATOL
    ok_contact = abs(res["trace_H"] + res["contact"]) <= 1e-9 * max(res["contact"], 1.0)

    print(f"\n  FD  agreement  : {'PASS' if ok_fd else 'FAIL'} (<= {FD_RTOL:g})")
    print(f"  QUAD agreement : {'PASS' if ok_quad else 'FAIL'} (<= {QUAD_RTOL:g})")
    print(f"  Tr(S) = 0      : {'PASS' if ok_trace else 'FAIL'} (<= {TRACE_ATOL:g})")
    print(f"  Tr(H) contact  : {'PASS' if ok_contact else 'FAIL'}")
    all_ok = ok_fd and ok_quad and ok_trace and ok_contact
    print("\n  PROTOTYPE SELF-CHECK:", "PASS" if all_ok else "FAIL")
    print("\n  NOTE: this validates the CLOSED FORM against analytic oracles only.")
    print("        The L1 stage gate (P1.3, vs OpenQP's real ERI engine) is NOT cleared here.")
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
