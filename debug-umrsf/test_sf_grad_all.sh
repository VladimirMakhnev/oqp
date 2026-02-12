#!/bin/bash
#
# Run all SF and MRSF gradient tests
# Output: PASS if all pass, detailed report if any fail
#

TESTS=(
    "test_uhf_sf_grad.sh"
    "test_rohf_sf_grad.sh"
    "test_uks_sf_grad.sh"
    "test_roks_sf_grad.sh"
    "test_rohf_mrsf_grad.sh"
    "test_roks_mrsf_grad.sh"
    "test_rohf_mrsf_grad_gmres.sh"
)
NAMES=(
    "SF UHF (pure HF)"
    "SF ROHF (pure HF)"
    "SF UKS (BHHLYP)"
    "SF ROKS (BHHLYP)"
    "MRSF ROHF (pure HF)"
    "MRSF ROKS (BHHLYP)"
    "MRSF ROHF GMRES (pure HF)"
)

RESULTS=()
OUTPUTS=()
ALL_PASS=1
N_TESTS=${#TESTS[@]}

for ((i=0; i<N_TESTS; i++)); do
    test=${TESTS[$i]}
    output=$(./$test 2>&1)
    if [ $? -eq 0 ]; then
        RESULTS+=("PASS")
    else
        RESULTS+=("FAIL")
        ALL_PASS=0
    fi
    OUTPUTS+=("$output")
done

if [ $ALL_PASS -eq 1 ]; then
    echo "PASS"
else
    echo "============================================================"
    echo "SF/MRSF GRADIENT TESTS SUMMARY"
    echo "============================================================"
    for ((i=0; i<N_TESTS; i++)); do
        echo "  ${NAMES[$i]}: ${RESULTS[$i]}"
    done
    echo "============================================================"
    echo ""
    for ((i=0; i<N_TESTS; i++)); do
        echo "==================== ${NAMES[$i]} ===================="
        echo "${OUTPUTS[$i]}"
        echo ""
    done
    exit 1
fi
