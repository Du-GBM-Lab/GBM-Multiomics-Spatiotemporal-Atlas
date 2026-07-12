#!/usr/bin/env bash
set -uo pipefail

cd /home/data/t010639/projects/GBM_R9_spatial_RCTD || exit 1
mkdir -p logs/R9_batch1_unbiased_landscape

run_bg() {
  local name="$1"
  local cmd="$2"
  local log="logs/R9_batch1_unbiased_landscape/${name}.log"
  local exit_file="logs/R9_batch1_unbiased_landscape/${name}.exit"
  rm -f "$log" "$exit_file"
  nohup bash -lc "/usr/bin/time -v ${cmd} > '${log}' 2>&1; echo \$? > '${exit_file}'" >/dev/null 2>&1 &
  echo "${name}.pid=$!"
}

run_bg "40_server_batch1_A1_SVG_SPARKX.rerun2" "Rscript scripts/40_server_batch1_A1_SVG_SPARKX.R"
run_bg "42_server_batch1_A3_squidpy.rerun4" "python3 scripts/42_server_batch1_A3_squidpy.py"

