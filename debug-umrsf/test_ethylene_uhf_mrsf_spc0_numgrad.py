#!/usr/bin/env python3
"""
Numerical gradient verification for UHF-MRSF(SPC=0) on ethylene.
Compares analytical vs numerical gradient for all 18 components.
"""

import subprocess
import os

# Planar ethylene (Angstrom)
atoms = [
    ('C', -0.6695,  0.0000,  0.0000),
    ('C',  0.6695,  0.0000,  0.0000),
    ('H', -1.2321,  0.9289,  0.0000),
    ('H', -1.2321, -0.9289,  0.0000),
    ('H',  1.2321,  0.9289,  0.0000),
    ('H',  1.2321, -0.9289,  0.0000),
]

h = 0.0005  # Angstrom
bohr_per_angstrom = 1.8897259886

def make_input(filename, coords, runtype='energy'):
    content = "[input]\nsystem=\n"
    for sym, x, y, z in coords:
        content += f" {sym}  {x:.10f}  {y:.10f}  {z:.10f}\n"
    content += f"""system2=
charge=0
basis=6-31g*
method=tdhf
runtype={runtype}

[guess]
type=huckel

[scf]
type=uhf
multiplicity=3
conv=1.0e-8

[tdhf]
type=mrsf
nstate=5
conv=1.0e-8
spc_coco=0.0
spc_ovov=0.0
spc_coov=0.0
"""
    if runtype == 'grad':
        content += "\n[properties]\ngrad=1\n"
    with open(filename, 'w') as f:
        f.write(content)

def run_oqp(inp_file):
    subprocess.run(['openqp', inp_file], capture_output=True)
    return inp_file.replace('.inp', '.log')

def extract_mrsf_state_energy(log_file, state=1):
    """Extract total energy for MRSF state from converged Davidson."""
    with open(log_file) as f:
        lines = f.readlines()
    # Find converged state energies - look for "State #  N  Energy ="
    for line in lines:
        if f'State #{state:4d}' in line and 'Energy' in line:
            # Format: State #   1  Energy =   -3.884687 eV
            parts = line.split('=')
            if len(parts) >= 2:
                ev_str = parts[1].strip().replace('eV', '').strip()
                return float(ev_str)
    return None

def extract_scf_energy(log_file):
    """Extract SCF total energy."""
    with open(log_file) as f:
        for line in f:
            if 'TOTAL energy =' in line and 'potential' not in line and 'kinetic' not in line:
                parts = line.split('=')
                return float(parts[-1].strip())
    return None

def extract_analytical_gradient(log_file):
    gradients = []
    with open(log_file) as f:
        lines = f.readlines()
    for i, line in enumerate(lines):
        if 'PyOQP electronic gradients' in line:
            for j in range(i+1, min(i+20, len(lines))):
                if 'PyOQP state' in lines[j]:
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

def main():
    print("=" * 70)
    print("Numerical Gradient: UHF-MRSF(SPC=0) on Ethylene")
    print("=" * 70)
    print(f"Step size h = {h} Angstrom")
    print()

    # Analytical gradient
    print("Running analytical gradient...")
    grad_inp = '_eth_spc0_grad.inp'
    make_input(grad_inp, atoms, runtype='grad')
    grad_log = run_oqp(grad_inp)
    anal = extract_analytical_gradient(grad_log)
    scf_e = extract_scf_energy(grad_log)
    state_ev = extract_mrsf_state_energy(grad_log, 1)
    print(f"  SCF energy: {scf_e}")
    print(f"  State 1 energy: {state_ev} eV")
    print(f"  Analytical gradient atoms: {len(anal)}")
    print()

    coord_names = ['X', 'Y', 'Z']
    fail_count = 0
    max_rel = 0

    print(f"{'Atom':<6} {'Crd':<4} {'Numerical':>14} {'Analytical':>14} {'AbsErr':>12} {'Rel%':>9}")
    print("-" * 65)

    for ai in range(6):
        for ci in range(3):
            # Skip Z for planar molecule (all zero)
            a_val = anal[ai][ci]

            coords_p = [(s, x, y, z) for s, x, y, z in atoms]
            coords_m = [(s, x, y, z) for s, x, y, z in atoms]
            lp = list(coords_p[ai])
            lm = list(coords_m[ai])
            lp[ci+1] += h
            lm[ci+1] -= h
            coords_p[ai] = tuple(lp)
            coords_m[ai] = tuple(lm)

            inp_p = f'_eth_spc0_a{ai}{coord_names[ci]}_p.inp'
            inp_m = f'_eth_spc0_a{ai}{coord_names[ci]}_m.inp'
            make_input(inp_p, coords_p)
            make_input(inp_m, coords_m)

            log_p = run_oqp(inp_p)
            log_m = run_oqp(inp_m)

            scf_p = extract_scf_energy(log_p)
            scf_m = extract_scf_energy(log_m)
            ev_p = extract_mrsf_state_energy(log_p, 1)
            ev_m = extract_mrsf_state_energy(log_m, 1)

            if ev_p is None or ev_m is None or scf_p is None or scf_m is None:
                print(f"{atoms[ai][0]}{ai+1:<5} {coord_names[ci]:<4} {'ERROR':>14} {a_val:>14.8f}")
                continue

            # Total energy = SCF + excitation(eV->Ha)
            ha_per_ev = 1.0 / 27.211386245988
            E_p = scf_p + ev_p * ha_per_ev
            E_m = scf_m + ev_m * ha_per_ev

            num_angstrom = (E_p - E_m) / (2 * h)
            num_bohr = num_angstrom / bohr_per_angstrom

            abs_err = abs(num_bohr - a_val)
            if abs(a_val) > 1e-8:
                rel = abs_err / abs(a_val) * 100
            elif abs(num_bohr) > 1e-8:
                rel = 999.9
            else:
                rel = 0.0

            max_rel = max(max_rel, rel)
            status = ""
            if rel > 1.0 and abs(a_val) > 1e-6:
                status = " FAIL"
                fail_count += 1
            elif rel > 0.5 and abs(a_val) > 1e-6:
                status = " WARN"

            sym = atoms[ai][0]
            print(f"{sym}{ai+1:<5} {coord_names[ci]:<4} {num_bohr:>14.8f} {a_val:>14.8f} {abs_err:>12.2e} {rel:>8.2f}%{status}")

    print("-" * 65)
    print(f"FAIL count: {fail_count}/18, Max relative error: {max_rel:.2f}%")

if __name__ == '__main__':
    main()
