#!/bin/bash
#
# Test for UHF Spin-Flip TDDFT gradient (pure HF)
# Molecule: C2F2H2 (low symmetry)
# Compares: SF energies, gradients
#

# ============ CONDA SETUP ============
source /opt/homebrew/Caskroom/miniforge/base/etc/profile.d/conda.sh
conda activate oqp

# ============ CONFIG ============
INPUT="oqp-uhf-sf-grad.inp"
LOG="oqp-uhf-sf-grad.log"

# Reference values (from GAMESS, pure HF)
REF_E1="-275.4835786552"   # SF state 1 (target)

# Reference gradient state 1 (from GAMESS)
#          X              Y              Z
REF_C1="-0.024725550   0.092942530  -0.506861602"
REF_C2="-0.030086315  -0.044804482  -0.425512623"
REF_F1=" 0.038413406  -0.272832066   0.425905027"
REF_F2=" 0.058261585   0.236442185   0.348077366"
REF_H1="-0.038241032  -0.073028970   0.081610800"
REF_H2="-0.003622094   0.061280804   0.076781032"

# Tolerances
TOL_E="1e-6"
TOL_G="1e-5"

# ============ FUNCTIONS ============
check_float() {
    local calc=$1 ref=$2 tol=$3 label=$4
    local result=$(awk -v c="$calc" -v r="$ref" -v t="$tol" 'BEGIN {
        diff = c - r; if (diff < 0) diff = -diff
        if (diff < t) { printf "PASS %.2e", diff } else { printf "FAIL %.2e", diff }
    }')
    local status=$(echo "$result" | awk '{print $1}')
    local diff=$(echo "$result" | awk '{print $2}')
    printf "  %-10s %12.6f (ref: %12.6f, diff: %s) [%s]\n" "$label" "$calc" "$ref" "$diff" "$status"
    [ "$status" = "PASS" ]
}

# ============ MAIN ============
echo "============================================================"
echo "TEST: UHF Spin-Flip TDDFT Gradient (C2F2H2, pure HF)"
echo "============================================================"

if [ ! -f "$INPUT" ]; then
    echo "ERROR: Input file '$INPUT' not found!"
    exit 1
fi

echo ""
echo "Running: openqp $INPUT"
echo "----------------------------------------"
openqp "$INPUT" 2>&1 | tail -3

if [ ! -f "$LOG" ]; then
    echo "ERROR: Log file '$LOG' not found!"
    exit 1
fi

# Parse energy state 1
E1=$(grep "PyOQP state 1" "$LOG" | head -1 | awk '{print $4}')

# Parse gradient (6 atoms)
GRAD=$(grep -A 8 "dE/dX.*dE/dY.*dE/dZ" "$LOG" | grep -E "^\s+[0-9]+\s+[0-9.]+" | head -6)

FAILED=0

echo ""
echo "============================================================"
echo "ENERGY STATE 1 (tolerance: $TOL_E Hartree)"
echo "============================================================"
check_float "$E1" "$REF_E1" "$TOL_E" "State 1" || FAILED=1

echo ""
echo "============================================================"
echo "GRADIENT STATE 1 (tolerance: $TOL_G Hartree/Bohr)"
echo "============================================================"

ATOMS=("C1" "C2" "F1" "F2" "H1" "H2")
REFS=("$REF_C1" "$REF_C2" "$REF_F1" "$REF_F2" "$REF_H1" "$REF_H2")
COMPS=("X" "Y" "Z")

for i in 1 2 3 4 5 6; do
    atom=${ATOMS[$((i-1))]}
    ref_line=${REFS[$((i-1))]}
    calc_line=$(echo "$GRAD" | awk -v n=$i 'NR==n {print $3, $4, $5}')

    for j in 1 2 3; do
        comp=${COMPS[$((j-1))]}
        calc=$(echo "$calc_line" | awk -v n=$j '{print $n}')
        ref=$(echo "$ref_line" | awk -v n=$j '{print $n}')
        check_float "$calc" "$ref" "$TOL_G" "$atom($comp)" || FAILED=1
    done
done

echo ""
echo "============================================================"
if [ $FAILED -eq 0 ]; then
    echo "RESULT: ALL TESTS PASSED"
else
    echo "RESULT: SOME TESTS FAILED"
fi
echo "============================================================"
exit $FAILED
