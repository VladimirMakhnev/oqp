#!/bin/bash
#
# Test for ROKS Mixed-Reference Spin-Flip TDDFT gradient (BHHLYP)
# Molecule: C2F2H2 (low symmetry)
# Compares: MRSF energies, gradients
#

# ============ CONDA SETUP ============
source /opt/homebrew/Caskroom/miniforge/base/etc/profile.d/conda.sh
conda activate oqp

# ============ CONFIG ============
INPUT="oqp-roks-mrsf-grad.inp"
LOG="oqp-roks-mrsf-grad.log"

# Reference values (from GAMESS, BHHLYP, MRSF singlet)
REF_E1="-276.6697378954"   # MRSF state 1 (target)

# Reference gradient state 1 (from GAMESS)
#          X              Y              Z
REF_C1="-0.045206786   0.129442880  -0.481739221"
REF_C2="-0.019043606  -0.086087915  -0.403366526"
REF_F1=" 0.047717907  -0.271325647   0.415980576"
REF_F2=" 0.050602207   0.227839725   0.338518684"
REF_H1="-0.030451744  -0.056082984   0.061231461"
REF_H2="-0.003617977   0.056213941   0.069375026"

# Tolerances (DFT: slightly relaxed)
TOL_E="1e-5"
TOL_G="1e-4"

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
echo "TEST: ROKS MRSF-TDDFT Gradient (C2F2H2, BHHLYP)"
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
