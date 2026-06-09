#!/usr/bin/env python
"""Phase-1 smoothness probe + Phase-2 FD reference seed.

Displaces the O atom along z by +/- delta, runs UMRSF energy at each point,
parses the absolute energy of the target state, and reports:
  - finite-difference gradient g_z(O)  (Hartree/Bohr)
  - second difference (curvature) as a smoothness indicator
  - per-point Jacobi-reordering convergence (segment ... converged)

A smooth PES with consistent MO reordering => |2nd diff| stays small and the
3-point energies are monotone/parabolic with no discrete jump.
"""
import os, re, subprocess, sys, tempfile, shutil

BASE = os.path.dirname(os.path.abspath(__file__))
ANG2BOHR = 1.8897259886
TARGET_STATE = 1          # MRSF singlet ground state (row "1 ... Hartree")
DELTA_ANG = 0.005         # displacement step, Angstrom

# Base geometry (Angstrom), O first then two H.
GEOM = [
    ("O", [0.0000000000, 0.0000000000,  0.1172143154]),
    ("H", [0.0000000000, 0.7572153434, -0.4688572616]),
    ("H", [0.0000000000, -0.7572153434, -0.4688572616]),
]

INP_TEMPLATE = """[input]
system=
{geom}
charge=0
method=tdhf
runtype=energy
functional=bhhlyp
basis=6-31g*

[guess]
type=huckel
save_mol=false

[scf]
type=uhf
multiplicity=3
converger_type=diis
maxit=200
conv=1.0e-9

[dftgrid]
rad_npts=96
ang_npts=302

[tdhf]
type=umrsf
nstate=5
multiplicity=1
conv=1.0e-9
zvconv=1.0e-10
"""

SEG_RE = re.compile(r"segment\s+(\d+)\s+converged at iter\s+(\d+)")


def parse_state_energy(txt, state):
    """Absolute energy (Hartree) of `state` from the final results table.

    The table is introduced by a header line containing both 'State' and
    'Energy', followed by a 'Hartree' units line, then rows:
        <idx>  <E_Hartree>  <Exc_eV>  <Exc_rel_GS> ...
    The reference row uses idx 0. We anchor on that header to avoid matching
    SCF/iteration lines elsewhere in the log.
    """
    lines = txt.splitlines()
    hdr = None
    for i, ln in enumerate(lines):
        if "State" in ln and "Energy" in ln and "Excitation" in ln:
            hdr = i
            break
    if hdr is None:
        raise RuntimeError("results table header not found")
    row_re = re.compile(r"^\s*(\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)")
    for ln in lines[hdr + 1:hdr + 12]:
        m = row_re.match(ln)
        if m and int(m.group(1)) == state:
            return float(m.group(2))
    raise RuntimeError("state %d not found in results table" % state)


def write_geom(dz):
    lines = []
    for i, (el, xyz) in enumerate(GEOM):
        x, y, z = xyz
        if i == 0:
            z = z + dz
        lines.append(" %s  %18.10f %18.10f %18.10f" % (el, x, y, z))
    return "\n".join(lines)


def run_point(dz, tag):
    inp = INP_TEMPLATE.format(geom=write_geom(dz))
    inp_path = os.path.join(BASE, "fd_%s.inp" % tag)
    with open(inp_path, "w") as f:
        f.write(inp)
    env = dict(os.environ, PYTHONNOUSERSITE="1", OMP_NUM_THREADS="4")
    subprocess.run(["openqp", inp_path], cwd=BASE, env=env,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    log_path = os.path.join(BASE, "fd_%s.log" % tag)
    txt = open(log_path).read()
    e = parse_state_energy(txt, TARGET_STATE)
    segs = SEG_RE.findall(txt)
    return e, segs


def main():
    pts = {}
    for dz, tag in [(-DELTA_ANG, "m"), (0.0, "0"), (DELTA_ANG, "p")]:
        e, segs = run_point(dz, tag)
        pts[tag] = e
        print("dz=%+.5f A  E(state %d)=%18.10f Ha  jacobi_segments=%s"
              % (dz, TARGET_STATE, e, segs))
    d = DELTA_ANG
    g_ang = (pts["p"] - pts["m"]) / (2 * d)          # Ha/Ang
    g_bohr = g_ang / ANG2BOHR                         # Ha/Bohr
    curv = (pts["p"] - 2 * pts["0"] + pts["m"]) / d**2
    print("\nFD gradient g_z(O) = %12.8f Ha/Bohr" % g_bohr)
    print("curvature (2nd diff) = %12.6f Ha/Ang^2  (smoothness indicator)" % curv)
    print("\nSmoothness verdict: %s" %
          ("LIKELY SMOOTH (parabolic, no jump)"
           if abs(curv) < 5.0 else "CHECK: large curvature, possible reorder switch"))


if __name__ == "__main__":
    main()
