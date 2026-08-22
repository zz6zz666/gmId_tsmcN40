#!/bin/bash
# =============================================================================
# L-level pipeline: simulate -> extract -> delete for TSMC N40 (CRN40LP)
# 1.1V core devices (nch / pch). 5 corners run in parallel per L.
#
# All machine-specific paths come from machine.env in this directory
# (no hardcoded directories or usernames in the scripts). Raw simulation
# data and final .h5 tables live under $HOME/simulation (see machine.env).
# =============================================================================

set -euo pipefail

# ----- Config ----------------------------------------------------------------
# Project dir = where this script lives (machine.env sits next to it).
SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Load machine config (env vars win over the file).
set -a
if [ -f "$SCRIPTDIR/machine.env" ]; then
    source "$SCRIPTDIR/machine.env"
fi
set +a
SCRIPTDIR="${SCRIPTDIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"

# Required settings (fail loudly if missing).
: "${RAWDIR:?machine.env must set RAWDIR (transient raw simulation output)}"
: "${MATDIR:?machine.env must set MATDIR (final .h5 output)}"
: "${PYTHON_BIN:?machine.env must set PYTHON_BIN (python>=3.10 with psf-parser)}"
: "${LOGDIR:?machine.env must set LOGDIR (pipeline logs / fail markers)}"
: "${TMPDIR:?machine.env must set TMPDIR (spectre work dirs)}"

PYTHON="$PYTHON_BIN"
FINE="--fine"
CORNERS="tt ff ss fs sf"

cleanup() {
    echo ""
    echo "[!] Interrupted. Cleaning up..."
    pkill -P $$ 2>/dev/null || true
    echo "[!] Done."
    exit 1
}
trap cleanup INT TERM

# ----- Helper: clean stale artifacts from a previous (aborted) run ----------
cleanup_stale() {
    echo "[*] Cleaning stale artifacts ..."
    mkdir -p "$LOGDIR" "$TMPDIR"
    # Leftover raw data from an aborted run
    rm -rf "$RAWDIR"/raw_tsmcN40_* 2>/dev/null || true
    # Old pipeline logs / fail markers
    rm -f "$LOGDIR"/*.log "$LOGDIR"/*.fail 2>/dev/null || true
    # Stale spectre work dirs (>= 1 day old)
    find "$TMPDIR" -maxdepth 1 -type d -name 'spec_work_*' -mtime +1 \
        -exec rm -rf {} + 2>/dev/null || true
}

# ----- Helper: run one corner simulation -------------------------------------
run_sim() {
    local corner=$1 l_start=$2 l_end=$3
    local logfile="$LOGDIR/sim_${corner}_L${l_start}.log"
    cd "$SCRIPTDIR"
    $PYTHON run_sim.py "$corner" $FINE --outdir "$RAWDIR" --L-range "$l_start" "$l_end" 2>&1 | tee "$logfile"
}

# ----- Helper: run one corner extraction -------------------------------------
run_extract() {
    local corner=$1 l_start=$2 l_end=$3
    local logfile="$LOGDIR/ext_${corner}_L${l_start}.log"
    cd "$SCRIPTDIR"
    $PYTHON extract_new.py $FINE \
        --srcdir "$RAWDIR" --outdir "$MATDIR" \
        --L-range "$l_start" "$l_end" --workers 1 "$corner" 2>&1 | tee "$logfile"
}

# ----- Helper: clean one L across all corners --------------------------------
clean_l() {
    local idx=$1
    local prefix="L$(printf '%03d' $idx)_"
    for d in "$RAWDIR"/raw_tsmcN40_*/; do
        if [ -d "$d" ]; then
            # Delete directories matching this L index
            find "$d" -maxdepth 1 -type d -name "${prefix}*" -exec rm -rf {} + 2>/dev/null || true
        fi
    done
}

# ----- Helper: check if a previous run failed --------------------------------
check_failure() {
    local stage=$1 l_idx=$2
    local failed=0
    for corner in $CORNERS; do
        if [ -f "$LOGDIR/${stage}_${corner}_L${l_idx}.fail" ]; then
            failed=1
            echo "    [ERROR] $stage failed for corner=$corner L=$l_idx"
            echo "    Log: $LOGDIR/${stage}_${corner}_L${l_idx}.log"
        fi
    done
    if [ $failed -ne 0 ]; then
        echo "*** Aborting due to failures above ***"
        exit 1
    fi
}

# ----- Main ------------------------------------------------------------------
echo "========================================"
echo "Pipeline (TSMC N40 1.1V core): 5 corners x nL, fine mode"
echo "Script:  $SCRIPTDIR"
echo "Raw:     $RAWDIR"
echo "MAT:     $MATDIR"
echo "Log:     $LOGDIR"
echo "========================================"

mkdir -p "$MATDIR" "$LOGDIR" "$TMPDIR"
cleanup_stale

cd "$SCRIPTDIR"
nL=$($PYTHON get_nl.py)

echo "[*] Total L count: $nL"
echo "[*] Corners: $CORNERS"
echo ""

t_total_start=$(date +%s)

for l in $(seq 0 $((nL - 1))); do
    l_next=$((l + 1))
    t_l_start=$(date +%s)
    echo ""
    echo "===== L=$l / $((nL - 1)) ====="

    # --- 1. Parallel simulation (5 corners) ---
    echo "[1/3] Simulating ..."
    pids=()
    for corner in $CORNERS; do
        (
            if run_sim "$corner" "$l" "$l_next"; then
                rm -f "$LOGDIR/sim_${corner}_L${l}.fail"
            else
                touch "$LOGDIR/sim_${corner}_L${l}.fail"
            fi
        ) &
        pids+=($!)
    done
    for pid in "${pids[@]}"; do wait "$pid"; done
    check_failure "sim" "$l"
    echo "      Simulation OK"

    # --- 2. Parallel extraction (5 corners) ---
    echo "[2/3] Extracting ..."
    pids=()
    for corner in $CORNERS; do
        (
            if run_extract "$corner" "$l" "$l_next"; then
                rm -f "$LOGDIR/ext_${corner}_L${l}.fail"
            else
                touch "$LOGDIR/ext_${corner}_L${l}.fail"
            fi
        ) &
        pids+=($!)
    done
    for pid in "${pids[@]}"; do wait "$pid"; done
    check_failure "ext" "$l"
    echo "      Extraction OK"

    # --- 3. Clean raw for this L ---
    echo "[3/3] Cleaning raw ..."
    clean_l "$l"
    echo "      Cleaned"

    # --- Timing ---
    t_l_end=$(date +%s)
    l_elapsed=$((t_l_end - t_l_start))
    total_elapsed=$((t_l_end - t_total_start))
    avg_l_time=$((total_elapsed / (l + 1)))
    remaining=$((nL - l - 1))
    eta=$((remaining * avg_l_time))

    printf "  [TIME] L=%d took %ds | Total %ds | Avg %ds/L | ETA %ds (~%dh%02dm)\n" \
        "$l" "$l_elapsed" "$total_elapsed" "$avg_l_time" "$eta" "$((eta / 3600))" "$((eta % 3600 / 60))"
done

t_total_end=$(date +%s)
total_time=$((t_total_end - t_total_start))
echo ""
echo "========================================"
printf "ALL DONE in %ds (~%dh%02dm%02ds)\n" "$total_time" "$((total_time / 3600))" "$((total_time % 3600 / 60))" "$((total_time % 60))"
echo "========================================"
