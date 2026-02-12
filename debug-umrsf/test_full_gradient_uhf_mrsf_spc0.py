#!/usr/bin/env python3
"""
Full numerical gradient verification for UHF-MRSF with SPC=0.
Tests ALL atoms and ALL coordinates (x, y, z).
Compares MRSF(SPC=0) analytical vs MRSF(SPC=0) numerical gradient.
"""

import subprocess
import os
import re

# Molecule: C2H2F2 (6 atoms)
atoms = [
    (6,  -0.0000000000,   0.6405243427,  -0.0185521420),  # C
    (6,   0.2000000000,  -0.6305243428,  -0.0165521421),  # C
    (9,  -0.1000000000,   1.2156446665,  -0.9292030778),  # F
    (9,   0.0000000000,  -1.2456446664,  -0.9192030780),  # F
    (1,  -0.1000000000,  -1.3965430797,   1.1378852200),  # H
    (1,   0.0000000000,   1.3905430798,   1.1178852199),  # H
]

h = 0.0005  # Step size in Angstrom
bohr_per_angstrom = 1.8897259886

def make_input(filename, coords, runtype='energy'):
    """Generate OQP input file for MRSF with SPC=0"""
    content = "[input]\nsystem=\n"
    for z, x, y, zc in coords:
        content += f"   {z}  {x:.10f}  {y:.10f}  {zc:.10f}\n"

    content += f"""charge=0
runtype={runtype}
basis=6-31g*
method=tdhf

[guess]
type=huckel

[scf]
multiplicity=3
type=uhf
jacobi_rotation=true

[tdhf]
type=mrsf
nstate=5
conv=1.0e-8
zvconv=1.0e-10
spc_coco=0.0
spc_ovov=0.0
spc_coov=0.0
"""
    if runtype == 'grad':
        content += """
[properties]
grad=1
"""

    with open(filename, 'w') as f:
        f.write(content)

def run_oqp(inp_file):
    subprocess.run(['openqp', inp_file], capture_output=True)
    return inp_file.replace('.inp', '.log')

def extract_energy(log_file, state=1):
    with open(log_file) as f:
        for line in f:
            if f'PyOQP state {state}' in line:
                parts = line.split()
                for p in parts:
                    if '-' in p and '.' in p:
                        try:
                            return float(p)
                        except:
                            continue
    return None

def extract_analytical_gradient(log_file, state=1):
    gradients = []
    with open(log_file) as f:
        lines = f.readlines()

    for i, line in enumerate(lines):
        if 'PyOQP electronic gradients' in line:
            for j in range(i+1, min(i+20, len(lines))):
                if f'PyOQP state {state}' in lines[j]:
                    for k in range(j+1, j+7):
                        if k >= len(lines):
                            break
                        parts = lines[k].split()
                        if len(parts) >= 4:
                            try:
                                gx = float(parts[1])
                                gy = float(parts[2])
                                gz = float(parts[3])
                                gradients.append((gx, gy, gz))
                            except (ValueError, IndexError):
                                pass
                    break
            break
    return gradients

def numerical_gradient_component(atom_idx, coord_idx, coords):
    coord_names = ['X', 'Y', 'Z']

    coords_p = [list(c) for c in coords]
    coords_p[atom_idx][coord_idx + 1] += h

    coords_m = [list(c) for c in coords]
    coords_m[atom_idx][coord_idx + 1] -= h

    inp_p = f'numgrad_spc0_a{atom_idx}_{coord_names[coord_idx]}_p.inp'
    inp_m = f'numgrad_spc0_a{atom_idx}_{coord_names[coord_idx]}_m.inp'

    make_input(inp_p, coords_p)
    make_input(inp_m, coords_m)

    log_p = run_oqp(inp_p)
    log_m = run_oqp(inp_m)

    E_p = extract_energy(log_p)
    E_m = extract_energy(log_m)

    if E_p is None or E_m is None:
        return None, None, None

    num_angstrom = (E_p - E_m) / (2 * h)
    num_bohr = num_angstrom / bohr_per_angstrom

    return num_bohr, E_p, E_m

def main():
    print("=" * 70)
    print("Full Numerical Gradient: UHF-MRSF (SPC=0, pure HF)")
    print("=" * 70)
    print(f"Step size h = {h} Angstrom")
    print(f"Molecule: C2H2F2 (6 atoms, 18 gradient components)")
    print("SPC=0: No spin-pair coupling (isolates structural bugs)")
    print()

    print("Running analytical gradient calculation...")
    grad_inp = 'numgrad_spc0_full_grad.inp'
    make_input(grad_inp, atoms, runtype='grad')
    grad_log = run_oqp(grad_inp)
    anal_grads = extract_analytical_gradient(grad_log)

    if len(anal_grads) != 6:
        print(f"ERROR: Expected 6 atoms in gradient, got {len(anal_grads)}")
        return

    print(f"Analytical gradients extracted for {len(anal_grads)} atoms")
    print()

    atom_symbols = ['C', 'C', 'F', 'F', 'H', 'H']
    coord_names = ['X', 'Y', 'Z']

    results = []
    max_error = 0.0
    max_rel_error = 0.0
    fail_count = 0

    print("Computing numerical gradients...")
    print("-" * 70)
    print(f"{'Atom':<6} {'Coord':<6} {'Numerical':>14} {'Analytical':>14} {'Abs Error':>12} {'Rel Error':>10}")
    print("-" * 70)

    for atom_idx in range(6):
        for coord_idx in range(3):
            num, E_p, E_m = numerical_gradient_component(atom_idx, coord_idx, atoms)
            anal = anal_grads[atom_idx][coord_idx]

            if num is None:
                print(f"{atom_symbols[atom_idx]}{atom_idx+1:<5} {coord_names[coord_idx]:<6} {'ERROR':>14} {anal:>14.6f}")
                continue

            abs_err = abs(num - anal)
            rel_err = abs_err / abs(anal) * 100 if anal != 0 else 0

            max_error = max(max_error, abs_err)
            max_rel_error = max(max_rel_error, rel_err)

            status = ""
            if rel_err > 1.0:
                status = " FAIL"
                fail_count += 1
            elif rel_err > 0.1:
                status = " WARN"

            print(f"{atom_symbols[atom_idx]}{atom_idx+1:<5} {coord_names[coord_idx]:<6} {num:>14.6f} {anal:>14.6f} {abs_err:>12.2e} {rel_err:>9.3f}%{status}")

    print("-" * 70)
    print()
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Max absolute error: {max_error:.2e} Hartree/Bohr")
    print(f"Max relative error: {max_rel_error:.3f}%")
    print(f"Components with >1% error: {fail_count}/18")
    print()

    if max_rel_error < 0.1:
        print("RESULT: PASS (all components < 0.1%)")
    elif max_rel_error < 1.0:
        print("RESULT: WARNING (some components 0.1-1.0%)")
    else:
        print("RESULT: FAIL (some components >= 1.0%)")

if __name__ == '__main__':
    main()
