#!/bin/bash
# Full path: what kind of RNA the swab reads come from (needs 07_qc_full.sh; markers also need 11_expression_qc.sh).
#   1. scripts/rna_composition.py  -> ${QC_DIR}/composition/ (composition.tsv, read_lengths.tsv,
#      genebody_coverage.tsv, transcript_span.tsv): every read put in one category (exonic spliced / unspliced,
#      intronic, intergenic, mitochondrial, rRNA, unmapped), read lengths per category, and coverage along canonical
#      transcripts >= 1.5 kb of genes outside the probe panel
#   2. scripts/composition_figures.py -> ${RESULTS_DIR}/figures/composition/ (composition, read_length, genebody,
#      markers; composition_figures.tsv)
# Reads every read of every sample: submitted as a job (8 CPUs). Run 10_slides.sh afterwards to add it to the deck.
# Usage: bash 12_rna_composition.sh
set -euo pipefail

_D="${TEQ_ANALYSIS_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
source "${_D}/lib.sh"

preflight() {
  compgen -G "${QC_DIR}/bam/*.bam" >/dev/null || { echo "no ${QC_DIR}/bam/*.bam: run 07_qc_full.sh first" >&2; exit 1; }
  [[ -s "${ANALYSIS_DIR}/panel_genes.txt" ]] || { echo "missing ${ANALYSIS_DIR}/panel_genes.txt" >&2; exit 1; }
  [[ -s "${REF_DIR}/${GTF_NAME}" ]] || { echo "missing ${REF_DIR}/${GTF_NAME}" >&2; exit 1; }
}
if [[ -z "${AFTER:-}" ]]; then preflight; fi

if submit_self epi_composition 8 32G 08:00:00 "${_D}/12_rna_composition.sh" "$@"; then
  exit 0
fi
preflight

activate_env "${ANALYSIS_ENV}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"
python "${ANALYSIS_DIR}/scripts/rna_composition.py" --bam-pattern "${QC_DIR}/bam/{sample}.bam" \
  --sample-sheet "${SAMPLE_SHEET}" --gtf "${REF_DIR}/${GTF_NAME}" --panel "${ANALYSIS_DIR}/panel_genes.txt" \
  --threads "${THREADS}" --outdir "${QC_DIR}/composition"
counts=()
if [[ -s "${QC_DIR}/expression/gene_counts.tsv" ]]; then
  counts=(--gene-counts "${QC_DIR}/expression/gene_counts.tsv")
else
  echo "[composition] no ${QC_DIR}/expression/gene_counts.tsv (11_expression_qc.sh): marker figure left out"
fi
python "${ANALYSIS_DIR}/scripts/composition_figures.py" --composition-dir "${QC_DIR}/composition" "${counts[@]}" \
  --expected "${ANALYSIS_DIR}/public_data/genome_expectation.tsv" --sample-sheet "${SAMPLE_SHEET}" \
  --outdir "${RESULTS_DIR}/figures/composition"
echo "[composition] done: ${RESULTS_DIR}/figures/composition"
echo "[composition] paste back: column -t ${QC_DIR}/composition/transcript_span.tsv"
