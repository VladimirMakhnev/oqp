#!/usr/bin/env python3
"""Quick ROHF-MRSF gradient check on ethylene - only X components"""

import subprocess

atoms = [
    (6,  -0.6695,   0.0000,   0.0000),
    (6,   0.6695,   0.0000,   0.0000),
    (1,  -1.2321,   0.9289,   0.0000),
    (1,  -1.2321,  -0.9289,   0.0000),
    (1,   1.2321,   0.9289,   0.0000),
    (1,   1.2321,  -0.9289,   0.0000),
]

h = 0.0005
bohr = 1.8897259886

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
type=rohf
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

def run_oqp(inp):
    subprocess.run(['openqp', inp], capture_output=True)
    return inp.replace('.inp', '.log')

def extract_energy(log, state=1):
    with open(log) as f:
        for line in f:
            if f'PyOQP state {state}' in line:
                for p in line.split():
                    if '-' in p and '.' in p:
                        try: return float(p)
                        except: continue
    return None

def extract_grad(log, natoms=6, state=1):
    grads = []
    with open(log) as f:
        lines = f.readlines()
    for i, line in enumerate(lines):
        if 'PyOQP electronic gradients' in line:
            for j in range(i+1, min(i+20, len(lines))):
                if f'PyOQP state {state}' in lines[j]:
                    for k in range(j+1, j+1+natoms):
                        parts = lines[k].split()
                        if len(parts) >= 4:
                            grads.append((float(parts[1]), float(parts[2]), float(parts[3])))
                    break
            break
    return grads

def main():
    print("ROHF-MRSF gradient check on ethylene (X components only)")
    print("=" * 60)

    make_input('eth_rohf_grad.inp', atoms, 'grad')
    ag = extract_grad(run_oqp('eth_rohf_grad.inp'))

    syms = ['C', 'C', 'H', 'H', 'H', 'H']
    print(f"{'Atom':<6} {'Numerical':>14} {'Analytical':>14} {'Rel Error':>10}")
    print("-" * 50)

    for ai in range(6):
        cp = [list(c) for c in atoms]; cm = [list(c) for c in atoms]
        cp[ai][1] += h; cm[ai][1] -= h
        make_input(f'eth_ro_p{ai}.inp', cp); make_input(f'eth_ro_m{ai}.inp', cm)
        Ep = extract_energy(run_oqp(f'eth_ro_p{ai}.inp'))
        Em = extract_energy(run_oqp(f'eth_ro_m{ai}.inp'))
        num = (Ep - Em) / (2*h) / bohr
        anal = ag[ai][0]
        rel = abs(num-anal)/abs(anal)*100 if abs(anal)>1e-8 else 0
        status = " FAIL" if rel > 1 else ""
        print(f"{syms[ai]}{ai+1:<5} {num:>14.6f} {anal:>14.6f} {rel:>9.3f}%{status}")

if __name__ == '__main__':
    main()
