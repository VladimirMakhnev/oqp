#!/bin/bash
#
# Compare ROHF-MRSF vs UHF-MRSF (both with SPC=0)
# ROHF-MRSF works correctly, UHF-MRSF has ~2.7% error
#
# Usage: ./compare_rohf_uhf_mrsf.sh
#

set -e

# Tolerance for comparison
TOL="1e-3"

echo "=============================================="
echo "  ROHF-MRSF vs UHF-MRSF (both SPC=0)"
echo "=============================================="
echo ""

# Activate conda environment
source /opt/homebrew/Caskroom/miniforge/base/etc/profile.d/conda.sh
conda activate oqp 2>/dev/null

# Create test inputs
ROHF_INPUT="test-rohf-mrsf-spc0.inp"
UHF_INPUT="test-uhf-mrsf-spc0.inp"

# Use the same molecule for both
cat > "$ROHF_INPUT" << 'EOF'
[input]
system=
 6     0.000000    -0.018552    -0.665001
 6     0.000000     0.018552     0.665001
 9    -1.140685    -0.237180    -1.256095
 9     1.140685     0.237180    -1.256095
 1    -0.010295     1.047821     1.056951
 1     0.010295    -1.047821     1.056951
charge=0
runtype=grad
basis=6-31g*
functional=
method=tdhf

[guess]
type=huckel

[scf]
multiplicity=3
type=rohf

[tdhf]
type=mrsf
nstate=5
conv=1.0e-8
zvconv=1.0e-10
spc_coco=0.0
spc_ovov=0.0
spc_coov=0.0

[properties]
grad=1
EOF

cat > "$UHF_INPUT" << 'EOF'
[input]
system=
 6     0.000000    -0.018552    -0.665001
 6     0.000000     0.018552     0.665001
 9    -1.140685    -0.237180    -1.256095
 9     1.140685     0.237180    -1.256095
 1    -0.010295     1.047821     1.056951
 1     0.010295    -1.047821     1.056951
charge=0
runtype=grad
basis=6-31g*
functional=
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
zvconv=1.0e-10
spc_coco=0.0
spc_ovov=0.0
spc_coov=0.0

[properties]
grad=1
EOF

echo "Running ROHF-MRSF(SPC=0)..."
openqp "$ROHF_INPUT" > /dev/null 2>&1
ROHF_LOG="${ROHF_INPUT%.inp}.log"

echo "Running UHF-MRSF(SPC=0)..."
openqp "$UHF_INPUT" > /dev/null 2>&1
UHF_LOG="${UHF_INPUT%.inp}.log"

echo ""
echo "Extracting and comparing debug values..."
echo ""

# Python script for extraction and comparison
python << PYTHON_EOF
import re
import sys

def extract_values(logfile, tags):
    """Extract values for given tags from logfile, preserving order"""
    values = {}
    order = []
    try:
        with open(logfile, 'r', errors='ignore') as f:
            content = f.read()

        # Extract reference (SCF) energy
        ref_match = re.search(r'Final (?:UHF|ROHF) energy is\s+([-\d.]+)', content)
        if ref_match:
            key = '[ENERGY] SCF_total'
            values[key] = float(ref_match.group(1))
            order.append(key)

        # Extract S0 energy
        s0_match = re.search(r'PyOQP state 0\s+([-\d.]+)', content)
        if s0_match:
            key = '[ENERGY] S0'
            values[key] = float(s0_match.group(1))
            order.append(key)

        # Extract S1 energy (target state)
        s1_match = re.search(r'PyOQP state 1\s+([-\d.]+)', content)
        if s1_match:
            key = '[ENERGY] S1'
            values[key] = float(s1_match.group(1))
            order.append(key)

        for tag in tags:
            # Pattern for "name norm="
            pattern_norm = rf'\[{tag}\]\s+(.*?)\s+norm\s*=\s*([-\d.E+]+)'
            for match in re.finditer(pattern_norm, content, re.IGNORECASE):
                name, val = match.groups()
                name = name.strip()
                key = f"[{tag}] {name} norm"
                if key not in values:
                    try:
                        values[key] = float(val)
                        order.append(key)
                    except:
                        pass

            # Pattern for "name=" without "norm"
            pattern_eq = rf'\[{tag}\]\s+(.*?)(?<!norm)\s*=\s*([-\d.E+]+)'
            for match in re.finditer(pattern_eq, content, re.IGNORECASE):
                name, val = match.groups()
                name = name.strip()
                if ' norm' in name.lower():
                    continue
                key = f"[{tag}] {name}"
                if key not in values:
                    try:
                        values[key] = float(val)
                        order.append(key)
                    except:
                        pass

        # Extract final gradient
        grad_match = re.search(r'Gradient \(Hartree/Bohr\).*?dE/dZ\s*\n\s*-+\s*\n(.*?)(?:\n\s*\n|\n\s*Maximum)',
                               content, re.DOTALL)
        if grad_match:
            lines = grad_match.group(1).strip().split('\n')
            for line in lines:
                parts = line.split()
                if len(parts) >= 5:
                    atom = parts[0]
                    key = f'[GRAD] atom{atom}_X'
                    values[key] = float(parts[2])
                    order.append(key)
                    key = f'[GRAD] atom{atom}_Y'
                    values[key] = float(parts[3])
                    order.append(key)
                    key = f'[GRAD] atom{atom}_Z'
                    values[key] = float(parts[4])
                    order.append(key)

    except Exception as e:
        print(f"Error reading {logfile}: {e}")
    return values, order

# Tags to extract
tags = ['ZVEC', 'MRSF_2E_GRAD', '2E_GRAD', 'MRSF_GRAD', 'UMRSF_GRD_INIT', 'SF_GRD_INIT']

rohf_log = "${ROHF_LOG}"
uhf_log = "${UHF_LOG}"

rohf_vals, rohf_order = extract_values(rohf_log, tags)
uhf_vals, uhf_order = extract_values(uhf_log, tags)

# Build ordered key list
seen = set()
all_keys = []
for k in rohf_order:
    if k not in seen:
        all_keys.append(k)
        seen.add(k)
for k in uhf_order:
    if k not in seen:
        all_keys.append(k)
        seen.add(k)

tol = float("${TOL}")
errors = 0
total = 0

print(f"{'Parameter':<55} {'ROHF-MRSF':>15} {'UHF-MRSF':>15} {'ABS Delta':>12} {'Status':>8}")
print("-" * 110)

for key in all_keys:
    rohf_val = rohf_vals.get(key, None)
    uhf_val = uhf_vals.get(key, None)

    if rohf_val is not None and uhf_val is not None:
        total += 1
        delta = abs(rohf_val - uhf_val)

        status = "OK" if delta < tol else "DIFF"

        if status == "DIFF":
            errors += 1
            color = "\033[1;33m"  # Yellow for differences
        else:
            color = "\033[0;32m"  # Green

        print(f"{key:<55} {rohf_val:>15.8e} {uhf_val:>15.8e} {delta:>12.2e} {color}{status:>8}\033[0m")
    elif rohf_val is not None:
        print(f"{key:<55} {rohf_val:>15.8e} {'N/A':>15} {'---':>12} {'ROHF':>8}")
    elif uhf_val is not None:
        print(f"{key:<55} {'N/A':>15} {uhf_val:>15.8e} {'---':>12} {'UHF':>8}")

print("-" * 110)
print(f"Total: {total} comparisons, {errors} differences (tolerance={tol})")

PYTHON_EOF

echo ""
echo "Done. Log files: $ROHF_LOG, $UHF_LOG"
