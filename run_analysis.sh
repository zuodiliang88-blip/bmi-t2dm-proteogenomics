#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "$0")"
export BMI_T2DM_ROOT="$PWD"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 R_DATATABLE_NUM_THREADS=4
mkdir -p results/mr/moloc_results results/finemapping/ADH1B/numcausal5 results/finemapping/ENTPD6/numcausal5
command -v Rscript >/dev/null
command -v plink >/dev/null
for script in 01_protein_bmi_mr.R 02_protein_t2dm_mr.R 03_select_candidates.R 04_coloc_bmi.R 05_coloc_t2dm.R 06_multitrait_coloc.R 07_finemap_adh1b.R 08_finemap_entpd6.R; do
  Rscript --vanilla "$script" > "results/${script%.R}.log" 2>&1
done
Rscript --vanilla 09_export_cs95.R > results/finemapping/credible_sets_95.json 2> results/09_export_cs95.log
