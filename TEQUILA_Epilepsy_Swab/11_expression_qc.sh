#!/bin/bash
# Full path, sample-level QC: gene counts of every read in the genome-wide BAMs of 07_qc_full.sh, the on-target rate
# on the whole probe panel (panel_genes.txt, 634 genes), PCA and sample correlation.
#   1. scripts/gene_counts.py   -> ${QC_DIR}/expression/gene_counts.tsv, panel_on_target.tsv, panel_genes_missing.txt
#   2. scripts/expression_qc.py -> ${RESULTS_DIR}/figures/expression_qc/ (panel_capture, pca, correlation,
#                                  kit_scatter; expression_qc.tsv)
# Reads every read of every sample: submitted as a job (8 CPUs). Run 10_slides.sh afterwards to add it to the deck.
# Usage: bash 11_expression_qc.sh
set -euo pipefail

_D="${TEQ_ANALYSIS_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
source "${_D}/lib.sh"

preflight() {
  compgen -G "${QC_DIR}/bam/*.bam" >/dev/null || { echo "no ${QC_DIR}/bam/*.bam: run 07_qc_full.sh first" >&2; exit 1; }
  [[ -s "${ANALYSIS_DIR}/panel_genes.txt" ]] || { echo "missing ${ANALYSIS_DIR}/panel_genes.txt" >&2; exit 1; }
}
if [[ -z "${AFTER:-}" ]]; then preflight; fi

if submit_self epi_expr_qc 8 32G 06:00:00 "${_D}/11_expression_qc.sh" "$@"; then
  exit 0
fi
preflight

activate_env "${ANALYSIS_ENV}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"
genes_csv="$(IFS=,; echo "${TARGET_GENES[*]}")"
python "${ANALYSIS_DIR}/scripts/gene_counts.py" --bam-pattern "${QC_DIR}/bam/{sample}.bam" \
  --sample-sheet "${SAMPLE_SHEET}" --gtf "${REF_DIR}/${GTF_NAME}" --panel "${ANALYSIS_DIR}/panel_genes.txt" \
  --report-genes "${genes_csv}" --threads "${THREADS}" --outdir "${QC_DIR}/expression"
python "${ANALYSIS_DIR}/scripts/expression_qc.py" --counts-dir "${QC_DIR}/expression" \
  --sample-sheet "${SAMPLE_SHEET}" --genes "${genes_csv}" --outdir "${RESULTS_DIR}/figures/expression_qc"
echo "[expression qc] done: ${RESULTS_DIR}/figures/expression_qc"
