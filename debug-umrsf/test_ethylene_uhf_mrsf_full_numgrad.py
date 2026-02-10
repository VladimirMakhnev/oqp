#!/usr/bin/env python3
"""
Numerical gradient verification for UHF-MRSF (full SPC) on ethylene.
Only X and Y components (Z=0 by symmetry).
"""

import subprocess

atoms = [
    ('C', -0.6695,  0.0000,  0.0000),
    ('C',  0.6695,  0.0000,  0.0000),
    ('H', -1.2321,  0.9289,  0.0000),
    ('H', -1.2321, -0.9289,  0.0000),
    ('H',  1.2321,  0.9289,  0.0000),
    ('H',  1.2321, -0.9289,  0.0000),
]

h = 0.0005
bohr_per_angstrom = 1.8897259886

def make_input(filename, coords, runtype='energy', spc=True):
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
"""
    if not spc:
        content += "spc_coco=0.0\nspc_ovov=0.0\nspc_coov=0.0\n"
    if runtype == 'grad':
        content += "\n[properties]\ngrad=1\n"
    with open(filename, 'w') as f:
        f.write(content)

def run_oqp(inp_file):
    subprocess.run(['openqp', inp_file], capture_output=True)
    return inp_file.replace('.inp', '.log')

def extract_scf_energy(log_file):
    with open(log_file) as f:
        for line in f:
            if 'TOTAL energy =' in line and 'potential' not in line and 'kinetic' not in line:
                return float(line.split('=')[-1].strip())
    return None

def extract_mrsf_state_energy(log_file, state=1):
    with open(log_file) as f:
        for line in f:
            if f'State #{state:4d}' in line and 'Energy' in line:
                parts = line.split('=')
                if len(parts) >= 2:
                    return float(parts[1].strip().replace('eV', '').strip())
    return None

def extract_gradient(log_file):
    grads = []
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
                                grads.append((float(parts[1]), float(parts[2]), float(parts[3])))
                            except:
                                pass
                    break
            break
    return grads

def main():
    print("=" * 70)
    print("Numerical Gradient: UHF-MRSF (full SPC) on Ethylene")
    print("=" * 70)

    # Analytical
    make_input('_eth_full_grad.inp', atoms, runtype='grad', spc=True)
    log = run_oqp('_eth_full_grad.inp')
    anal = extract_gradient(log)
    scf = extract_scf_energy(log)
    ev = extract_mrsf_state_energy(log, 1)
    print(f"SCF: {scf}, State 1: {ev} eV")
    print(f"Analytical gradient: {len(anal)} atoms")
    print()

    coord_names = ['X', 'Y']
    fail_count = 0
    max_rel = 0

    print(f"{'Atom':<6} {'Crd':<4} {'Numerical':>14} {'Analytical':>14} {'AbsErr':>12} {'Rel%':>9}")
    print("-" * 65)

    ha_per_ev = 1.0 / 27.211386245988

    for ai in range(6):
        for ci in range(2):  # Only X and Y (Z=0)
            a_val = anal[ai][ci]

            coords_p = [(s, x, y, z) for s, x, y, z in atoms]
            coords_m = [(s, x, y, z) for s, x, y, z in atoms]
            lp = list(coords_p[ai]); lp[ci+1] += h; coords_p[ai] = tuple(lp)
            lm = list(coords_m[ai]); lm[ci+1] -= h; coords_m[ai] = tuple(lm)

            inp_p = f'_eth_full_a{ai}{coord_names[ci]}_p.inp'
            inp_m = f'_eth_full_a{ai}{coord_names[ci]}_m.inp'
            make_input(inp_p, coords_p, spc=True)
            make_input(inp_m, coords_m, spc=True)

            log_p = run_oqp(inp_p)
            log_m = run_oqp(inp_m)

            scf_p = extract_scf_energy(log_p)
            scf_m = extract_scf_energy(log_m)
            ev_p = extract_mrsf_state_energy(log_p, 1)
            ev_m = extract_mrsf_state_energy(log_m, 1)

            if None in (scf_p, scf_m, ev_p, ev_m):
                print(f"{atoms[ai][0]}{ai+1:<5} {coord_names[ci]:<4} {'ERROR':>14} {a_val:>14.8f}")
                continue

            E_p = scf_p + ev_p * ha_per_ev
            E_m = scf_m + ev_m * ha_per_ev

            num = (E_p - E_m) / (2 * h) / bohr_per_angstrom

            abs_err = abs(num - a_val)
            if abs(a_val) > 1e-8:
                rel = abs_err / abs(a_val) * 100
            elif abs(num) > 1e-8:
                rel = 999.9
            else:
                rel = 0.0

            max_rel = max(max_rel, rel)
            status = ""
            if rel > 5.0 and abs(a_val) > 1e-6:
                status = " FAIL"
                fail_count += 1
            elif rel > 1.0 and abs(a_val) > 1e-6:
                status = " WARN"

            print(f"{atoms[ai][0]}{ai+1:<5} {coord_names[ci]:<4} {num:>14.8f} {a_val:>14.8f} {abs_err:>12.2e} {rel:>8.2f}%{status}")

    print("-" * 65)
    print(f"FAIL count: {fail_count}/12, Max relative error: {max_rel:.2f}%")

if __name__ == '__main__':
    main()
