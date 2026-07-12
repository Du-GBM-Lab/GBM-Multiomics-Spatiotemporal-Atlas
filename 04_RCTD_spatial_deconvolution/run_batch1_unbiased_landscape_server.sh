#!/usr/bin/env bash
set -uo pipefail

cd /home/data/t010639/projects/GBM_R9_spatial_RCTD || exit 1

mkdir -p logs/R9_batch1_unbiased_landscape

run_one() {
  local name="$1"
  local cmd="$2"
  local log="logs/R9_batch1_unbiased_landscape/${name}.log"
  local exit_file="logs/R9_batch1_unbiased_landscape/${name}.exit"
  rm -f "$log" "$exit_file"
  echo "== [start] $name"
  /usr/bin/time -v bash -lc "$cmd" > "$log" 2>&1
  local status=$?
  echo "$status" > "$exit_file"
  echo "== [done] $name status=$status"
  return "$status"
}

# A2/A4 are R-only and should be lightweight relative to A1/A3.
run_one "41_server_batch1_A2_unbiased_domain" "Rscript scripts/41_server_batch1_A2_unbiased_domain.R"
run_one "43_server_batch1_A4_cellular_neighborhood" "Rscript scripts/43_server_batch1_A4_cellular_neighborhood.R"

# A1 and A3 are dependency-heavy. They stop early if required packages are absent.
run_one "40_server_batch1_A1_SVG_SPARKX" "Rscript scripts/40_server_batch1_A1_SVG_SPARKX.R"

if command -v python3 >/dev/null 2>&1; then
  run_one "42_server_batch1_A3_squidpy" "python3 scripts/42_server_batch1_A3_squidpy.py"
else
  echo "127" > logs/R9_batch1_unbiased_landscape/42_server_batch1_A3_squidpy.exit
  echo "python3 not found" > logs/R9_batch1_unbiased_landscape/42_server_batch1_A3_squidpy.log
fi

echo "== batch1 server wrapper complete; inspect per-task .exit and .log files."
