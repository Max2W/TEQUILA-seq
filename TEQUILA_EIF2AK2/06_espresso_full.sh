#!/bin/bash
# FULL PATH: the standard TEQUILA-seq Snakemake workflow on all reads
# (minimap2 -> samtools sort -> ESPRESSO S/C/Q), one input per flowcell per sample.
# Runs the Snakemake controller as a long, light SLURM job; it submits the
# alignment / ESPRESSO jobs itself through snakemake_profile (sbatch).
# Output: ${PIPELINE_DIR}/espresso_out/work_dir/samples_N2_R0_{abundance.esp,updated.gtf}
#
# Usage: bash 06_espresso_full.sh            write config + submit controller
#        DRYRUN=1 bash 06_espresso_full.sh   write config + snakemake -n here
set -euo pipefail

_D="${TEQ_ANALYSIS_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
source "${_D}/lib.sh"

[[ -x "${PIPELINE_ENV}/bin/snakemake" ]] || { echo "pipeline env missing: run 00_setup.sh envs" >&2; exit 1; }
[[ -L "${PIPELINE_DIR}/Snakefile" ]] || { echo "pipeline dir not set up: run 00_setup.sh pipeline" >&2; exit 1; }

make_config() {
  local fc_labels=() entry
  for entry in "${FLOWCELLS[@]}"; do fc_labels+=("${entry%%:*}"); done
  python3 "${ANALYSIS_DIR}/scripts/make_snakemake_config.py" \
    --sample-sheet "${SAMPLE_SHEET}" --fastq-dir "${FASTQ_DIR}" \
    --flowcells "$(IFS=,; echo "${fc_labels[*]}")" \
    --read-counts "${QC_DIR}/fastq_read_counts.tsv" --min-reads 1 \
    --code-dir "${CODE_DIR}" --gtf-name "${GTF_NAME}" --fasta-name "${FASTA_NAME}" \
    --threads 16 \
    --out "${PIPELINE_DIR}/snakemake_config.yaml" --index-out "${PIPELINE_DIR}/input_index.tsv"
}

# Written on the login node too (unless waiting on step 01), so problems show up at once
if [[ "${TEQ_IN_JOB:-0}" != "1" && -z "${AFTER:-}" ]]; then
  make_config
  if [[ "${DRYRUN:-0}" == "1" ]]; then
    cd "${PIPELINE_DIR}"
    "${CODE_DIR}/conda_wrapper" snakemake -n --quiet
    exit 0
  fi
fi

if submit_self eif2ak2_espresso_full 1 8G 7-00:00:00 "${_D}/06_espresso_full.sh" "$@"; then
  exit 0
fi

make_config
cd "${PIPELINE_DIR}"
"${CODE_DIR}/conda_wrapper" snakemake --profile ./snakemake_profile
