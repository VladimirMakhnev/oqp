#!/bin/bash
# Numerical gradient verification for UHF-MRSF (pure HF)
# Tests Atom 1 (C), Z coordinate
# NOTE: UHF-MRSF gradient may have larger errors - under development

source /opt/homebrew/Caskroom/miniforge/base/etc/profile.d/conda.sh
conda activate oqp
cd /Users/stan/Codebase/oqp_dir/openqp

h=0.0005  # Step size in Angstrom
z0=-0.0185521420  # Atom 1 C, Z coordinate

echo "=== Numerical Gradient: UHF-MRSF (pure HF) ==="
echo "Atom 1 (C), Z coordinate, h=$h Angstrom"
echo ""

for tag in p1 m1; do
    case $tag in
        p1) z=$(python3 -c "print(f'{$z0 + $h:.10f}')");;
        m1) z=$(python3 -c "print(f'{$z0 - $h:.10f}')");;
    esac

    cat > /tmp/num_uhf_mrsf_${tag}.inp << EOF
[input]
system=
   6  -0.0000000000   0.6405243427  $z
   6   0.2000000000  -0.6305243428  -0.0165521421
   9  -0.1000000000   1.2156446665  -0.9292030778
   9   0.0000000000  -1.2456446664  -0.9192030780
   1  -0.1000000000  -1.3965430797   1.1378852200
   1   0.0000000000   1.3905430798   1.1178852199
charge=0
runtype=energy
basis=6-31g*
method=tdhf

[guess]
type=huckel

[scf]
multiplicity=3
type=uhf

[tdhf]
type=mrsf
nstate=5
conv=1.0e-8
EOF
    echo "Running $tag (z=$z)..."
    openqp /tmp/num_uhf_mrsf_${tag}.inp 2>&1 > /dev/null
done

python3 << 'PY'
bohr_per_angstrom = 1.8897259886
h = 0.0005

# Extract energies
E_p1 = E_m1 = None
with open('/tmp/num_uhf_mrsf_p1.log') as f:
    for line in f:
        if 'PyOQP state 1' in line:
            parts = line.split()
            if len(parts) >= 4 and '-' in parts[3]:
                E_p1 = float(parts[3])
                break

with open('/tmp/num_uhf_mrsf_m1.log') as f:
    for line in f:
        if 'PyOQP state 1' in line:
            parts = line.split()
            if len(parts) >= 4 and '-' in parts[3]:
                E_m1 = float(parts[3])
                break

if E_p1 is None or E_m1 is None:
    print("ERROR: Could not extract energies from log files")
    exit(1)

num_angstrom = (E_p1 - E_m1) / (2*h)
num_bohr = num_angstrom / bohr_per_angstrom

# Get analytical gradient from running the gradient calculation
import subprocess
import os
os.system('''
cat > /tmp/num_uhf_mrsf_grad.inp << EOF
[input]
system=
   6  -0.0000000000   0.6405243427  -0.0185521420
   6   0.2000000000  -0.6305243428  -0.0165521421
   9  -0.1000000000   1.2156446665  -0.9292030778
   9   0.0000000000  -1.2456446664  -0.9192030780
   1  -0.1000000000  -1.3965430797   1.1378852200
   1   0.0000000000   1.3905430798   1.1178852199
charge=0
runtype=grad
basis=6-31g*
method=tdhf

[guess]
type=huckel

[scf]
multiplicity=3
type=uhf

[tdhf]
type=mrsf
nstate=5
conv=1.0e-8

[properties]
grad=1
EOF
openqp /tmp/num_uhf_mrsf_grad.inp 2>&1 > /dev/null
''')
anal = 0.0
with open('/tmp/num_uhf_mrsf_grad.log') as f:
    lines = f.readlines()
    for i, line in enumerate(lines):
        if 'PyOQP state 1' in line and i+1 < len(lines) and 'C ' in lines[i+1]:
            anal = float(lines[i+1].split()[3])
            break

print("")
print("=== UHF-MRSF (pure HF) Results ===")
print(f"E(+h) = {E_p1:.8f} Hartree")
print(f"E(-h) = {E_m1:.8f} Hartree")
print(f"")
print(f"Numerical gradient:  {num_bohr:.6f} Hartree/Bohr")
print(f"Analytical gradient: {anal:.6f} Hartree/Bohr")
print(f"Difference:          {abs(num_bohr - anal):.2e}")
if anal != 0:
    print(f"Relative error:      {100*abs(num_bohr - anal)/abs(anal):.3f}%")
    if abs(num_bohr - anal)/abs(anal) < 0.001:
        print("\nPASS: Error < 0.1%")
    elif abs(num_bohr - anal)/abs(anal) < 0.01:
        print("\nWARNING: Error 0.1-1%")
    else:
        print("\nFAIL: Error >= 1%")
PY
