#!/usr/bin/env bash
set -uo pipefail

cd /home/data/t010639/projects/GBM_R9_spatial_RCTD || exit 1

LOG="logs/06_server_A2_RCTD_full_all_slices_serial_safe.log"
EXIT="logs/06_server_A2_RCTD_full_all_slices_serial_safe.exit"

rm -f "$LOG" "$EXIT"

# Limit per-process virtual memory to 50 GiB-equivalent KiB.
ulimit -v 52428800

/usr/bin/time -v Rscript scripts/06_server_A2_RCTD_full_all_slices_serial_safe.R > "$LOG" 2>&1
status=$?
echo "$status" > "$EXIT"
exit "$status"
