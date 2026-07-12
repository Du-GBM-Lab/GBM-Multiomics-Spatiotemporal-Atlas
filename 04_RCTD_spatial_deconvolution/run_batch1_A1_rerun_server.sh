#!/usr/bin/env bash
set -uo pipefail

cd /home/data/t010639/projects/GBM_R9_spatial_RCTD || exit 1
mkdir -p logs/R9_batch1_unbiased_landscape

name="40_server_batch1_A1_SVG_SPARKX.rerun4"
log="logs/R9_batch1_unbiased_landscape/${name}.log"
exit_file="logs/R9_batch1_unbiased_landscape/${name}.exit"
rm -f "$log" "$exit_file"

nohup bash -lc "/usr/bin/time -v Rscript scripts/40_server_batch1_A1_SVG_SPARKX.R > '${log}' 2>&1; echo \$? > '${exit_file}'" >/dev/null 2>&1 &
echo "${name}.pid=$!"

