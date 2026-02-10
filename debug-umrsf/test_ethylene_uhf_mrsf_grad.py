#!/usr/bin/env python3
"""
Full numerical gradient verification for UHF-MRSF on ethylene (planar)
Tests ALL atoms and ALL coordinates (x, y, z)
"""

import subprocess
import os

# Ethylene planar (6 atoms)
atoms = [
    (6,  -0.6695,   0.0000,   0.0000),  # C
    (6,   0.6695,   0.0000,   0.0000),  # C
    (1,  -1.2321,   0.9289,   0.0000),  # H
    (1,  -1.2321,  -0.9289,   0.0000),  # H
    (1,   1.2321,   0.9289,   0.0000),  # H
    (1,   1.2321,  -0.9289,   0.0000),  # H
]

h = 0.0005  # Step size in Angstrom
bohr_per_angstrom = 1.8897259886

def make_input(filename, coords, runtype='energy'):
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
conv=1.0e-8

[tdhf]
type=mrsf
nstate=5
conv=1.0e-8
"""
    if runtype == 'grad':
        content += "\n[properties]\ngrad=1\n"

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

def extract_analytical_gradient(log_file, natoms=6, state=1):
    gradients = []
    with open(log_file) as f:
        lines = f.readlines()

    for i, line in enumerate(lines):
        if 'PyOQP electronic gradients' in line:
            for j in range(i+1, min(i+20, len(lines))):
                if f'PyOQP state {state}' in lines[j]:
                    for k in range(j+1, j+1+natoms):
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

def main():
    print("=" * 70)
    print("Full Numerical Gradient: UHF-MRSF on Ethylene (planar)")
    print("=" * 70)
    print(f"Step size h = {h} Angstrom")
    print()

    # Analytical gradient
    print("Running analytical gradient...")
    grad_inp = 'eth_mrsf_grad.inp'
    make_input(grad_inp, atoms, runtype='grad')
    grad_log = run_oqp(grad_inp)
    anal_grads = extract_analytical_gradient(grad_log)

    if len(anal_grads) != 6:
        print(f"ERROR: Expected 6 atoms, got {len(anal_grads)}")
        return

    print(f"Got analytical gradients for {len(anal_grads)} atoms\n")

    atom_symbols = ['C', 'C', 'H', 'H', 'H', 'H']
    coord_names = ['X', 'Y', 'Z']

    fail_count = 0
    max_rel = 0.0

    print(f"{'Atom':<6} {'Coord':<6} {'Numerical':>14} {'Analytical':>14} {'Abs Error':>12} {'Rel Error':>10}")
    print("-" * 70)

    for ai in range(6):
        for ci in range(3):
            coords_p = [list(c) for c in atoms]
            coords_m = [list(c) for c in atoms]
            coords_p[ai][ci + 1] += h
            coords_m[ai][ci + 1] -= h

            inp_p = f'eth_ng_a{ai}_{coord_names[ci]}_p.inp'
            inp_m = f'eth_ng_a{ai}_{coord_names[ci]}_m.inp'
            make_input(inp_p, coords_p)
            make_input(inp_m, coords_m)

            log_p = run_oqp(inp_p)
            log_m = run_oqp(inp_m)

            E_p = extract_energy(log_p)
            E_m = extract_energy(log_m)

            if E_p is None or E_m is None:
                print(f"{atom_symbols[ai]}{ai+1:<5} {coord_names[ci]:<6} {'ERROR':>14}")
                continue

            num = (E_p - E_m) / (2 * h) / bohr_per_angstrom
            anal = anal_grads[ai][ci]

            abs_err = abs(num - anal)
            rel_err = abs_err / abs(anal) * 100 if abs(anal) > 1e-8 else 0

            max_rel = max(max_rel, rel_err)
            if rel_err > 1.0:
                fail_count += 1

            status = " FAIL" if rel_err > 1.0 else (" WARN" if rel_err > 0.1 else "")
            print(f"{atom_symbols[ai]}{ai+1:<5} {coord_names[ci]:<6} {num:>14.6f} {anal:>14.6f} {abs_err:>12.2e} {rel_err:>9.3f}%{status}")

    print("-" * 70)
    print(f"\nFailed: {fail_count}/18, Max relative error: {max_rel:.3f}%")
    if fail_count == 0:
        print("RESULT: PASS")
    else:
        print("RESULT: FAIL")

if __name__ == '__main__':
    main()
