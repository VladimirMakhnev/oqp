#!/bin/bash
#
# Compare UHF-SF vs UHF-MRSF(SPC=0) debug outputs
# These methods MUST produce identical results when SPC=0
#
# Usage: ./compare_uhf_sf_mrsf.sh
#

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Tolerance for comparison (1e-3 for UHF method comparison)
TOL="1e-3"

echo "=============================================="
echo "  UHF-SF vs UHF-MRSF(SPC=0) Comparison"
echo "=============================================="
echo ""

# Activate conda environment
source /opt/homebrew/Caskroom/miniforge/base/etc/profile.d/conda.sh
conda activate oqp 2>/dev/null

# Create test inputs if they don't exist
SF_INPUT="test-uhf-sf-compare.inp"
MRSF_INPUT="test-uhf-mrsf-compare.inp"

# Use the same molecule for both
cat > "$SF_INPUT" << 'EOF'
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
type=sf
nstate=5
conv=1.0e-8
zvconv=1.0e-10

[properties]
grad=1
EOF

cat > "$MRSF_INPUT" << 'EOF'
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

echo "Running UHF-SF..."
openqp "$SF_INPUT" > /dev/null 2>&1
SF_LOG="${SF_INPUT%.inp}.log"

echo "Running UHF-MRSF(SPC=0)..."
openqp "$MRSF_INPUT" > /dev/null 2>&1
MRSF_LOG="${MRSF_INPUT%.inp}.log"

echo ""
echo "Extracting and comparing debug values..."
echo ""

# Python script for extraction and comparison
python3 << PYTHON_EOF
import re
import sys

def extract_values(logfile, tags):
    """Extract values for given tags from logfile, preserving order"""
    values = {}
    order = []  # Track order of appearance
    try:
        with open(logfile, 'r', errors='ignore') as f:
            content = f.read()

        # Extract reference (SCF) energy
        ref_match = re.search(r'Final UHF energy is\s+([-\d.]+)', content)
        if ref_match:
            key = '[ENERGY] SCF_total'
            values[key] = float(ref_match.group(1))
            order.append(key)

        # Extract ONLY S0 energy (ground state reference)
        s0_match = re.search(r'PyOQP state 0\s+([-\d.]+)', content)
        if s0_match:
            key = '[ENERGY] S0'
            values[key] = float(s0_match.group(1))
            order.append(key)

        for tag in tags:
            # Pattern for "name norm=" (e.g., [ZVEC] STEP0 mo_energy_a norm= value)
            # Must check norm= pattern BEFORE simple = pattern
            pattern_norm = rf'\[{tag}\]\s+(.*?)\s+norm\s*=\s*([-\d.E+]+)'
            for match in re.finditer(pattern_norm, content, re.IGNORECASE):
                name, val = match.groups()
                name = name.strip()
                key = f"[{tag}] {name} norm"
                if key not in values:  # First occurrence only
                    try:
                        values[key] = float(val)
                        order.append(key)
                    except:
                        pass

            # Pattern for "name=" without "norm" (e.g., [ZVEC] STEP0 target_state= value)
            pattern_eq = rf'\[{tag}\]\s+(.*?)(?<!norm)\s*=\s*([-\d.E+]+)'
            for match in re.finditer(pattern_eq, content, re.IGNORECASE):
                name, val = match.groups()
                name = name.strip()
                # Skip if this is a norm= pattern (already captured above)
                if ' norm' in name.lower():
                    continue
                key = f"[{tag}] {name}"
                if key not in values:
                    try:
                        values[key] = float(val)
                        order.append(key)
                    except:
                        pass
    except Exception as e:
        print(f"Error reading {logfile}: {e}")
    return values, order

# Tags to extract (unified format)
tags = ['ZVEC', '2E_GRAD', 'SF_2E_GRAD', 'MRSF_2E_GRAD']

sf_log = "${SF_LOG}"
mrsf_log = "${MRSF_LOG}"

sf_vals, sf_order = extract_values(sf_log, ['ZVEC', 'SF_2E_GRAD', '2E_GRAD'])
mrsf_vals, mrsf_order = extract_values(mrsf_log, ['ZVEC', 'MRSF_2E_GRAD', '2E_GRAD'])

# Normalize tag names for comparison (remove SF_/MRSF_ prefix and unify variable names)
def normalize_key(key):
    # Normalize tag prefixes
    key = key.replace('[SF_ZVEC]', '[ZVEC]')
    key = key.replace('[MRSF_ZVEC]', '[ZVEC]')
    key = key.replace('[SF_2E_GRAD]', '[2E_GRAD]')
    key = key.replace('[MRSF_2E_GRAD]', '[2E_GRAD]')
    key = key.replace('[MR2E_GRAD]', '[2E_GRAD]')
    key = key.replace('[MRZVEC]', '[ZVEC]')

    # Normalize equivalent variable names
    key = key.replace('bvec_mo(1:5,target)', 'bvec_mo(1:5)')
    key = key.replace('bvec_mo(6:10,target)', 'bvec_mo(6:10)')
    key = key.replace('bvec_mo(11:15,target)', 'bvec_mo(11:15)')
    key = key.replace('bvec_mo(target_state) norm', 'bvec_mo norm')
    key = key.replace('fmrst2(1,11) norm (=ab2 in SF)', 'H[X]_AO norm')
    key = key.replace('ab2 norm', 'H[X]_AO norm')
    key = key.replace('fmrst2(1,11,1,1:5)', 'H[X]_AO(1,1:5)')
    key = key.replace('fmrst2(1,11,5,1:5)', 'H[X]_AO(5,1:5)')
    key = key.replace('ab2(1,1:5)', 'H[X]_AO(1,1:5)')
    key = key.replace('ab2(5,1:5)', 'H[X]_AO(5,1:5)')
    key = key.replace('wrk2 norm (H[X] in MO)', 'H[X]_MO norm')
    key = key.replace('wrk3 norm (X expanded)', 'X_expanded norm')
    key = key.replace('bvec_mo_d(1:5)', 'bvec_mo(1:5)')  # SPC=0: bvec_mo_d = bvec_mo
    # Normalize INPUT parameter names (SF_ZVEC/MRSF_ZVEC already normalized to ZVEC)
    return key

sf_normalized = {normalize_key(k): v for k, v in sf_vals.items()}
mrsf_normalized = {normalize_key(k): v for k, v in mrsf_vals.items()}

# Build ordered key list based on SF order (reference), then add any MRSF-only keys
sf_ordered_normalized = [normalize_key(k) for k in sf_order]
mrsf_ordered_normalized = [normalize_key(k) for k in mrsf_order]

# Use SF order as primary, preserving appearance order
seen = set()
all_keys = []
for k in sf_ordered_normalized:
    if k not in seen:
        all_keys.append(k)
        seen.add(k)
# Add MRSF-only keys at the end
for k in mrsf_ordered_normalized:
    if k not in seen:
        all_keys.append(k)
        seen.add(k)

tol = float("${TOL}")
errors = 0
total = 0

print(f"{'Parameter':<50} {'U-SF':>15} {'U-MRSF(0)':>15} {'ABS Delta':>12} {'Status':>8}")
print("-" * 105)

for key in all_keys:
    sf_val = sf_normalized.get(key, None)
    mrsf_val = mrsf_normalized.get(key, None)

    if sf_val is not None and mrsf_val is not None:
        total += 1
        delta = abs(sf_val - mrsf_val)

        # Absolute comparison with tolerance
        status = "OK" if delta < tol else "ERROR"

        if status == "ERROR":
            errors += 1
            color = "\033[0;31m"  # Red
        else:
            color = "\033[0;32m"  # Green

        print(f"{key:<50} {sf_val:>15.8e} {mrsf_val:>15.8e} {delta:>12.2e} {color}{status:>8}\033[0m")
    elif sf_val is not None:
        print(f"{key:<50} {sf_val:>15.8e} {'N/A':>15} {'---':>12} {'MISS':>8}")
    elif mrsf_val is not None:
        print(f"{key:<50} {'N/A':>15} {mrsf_val:>15.8e} {'---':>12} {'MISS':>8}")

print("-" * 105)
print(f"Total: {total} comparisons, {errors} errors")

if errors > 0:
    print("\n\033[0;31mFAILED: Values diverge between UHF-SF and UHF-MRSF(SPC=0)\033[0m")
    sys.exit(1)
else:
    print("\n\033[0;32mPASSED: All values match within tolerance\033[0m")
    sys.exit(0)
PYTHON_EOF

echo ""
echo "Done. Log files: $SF_LOG, $MRSF_LOG"
