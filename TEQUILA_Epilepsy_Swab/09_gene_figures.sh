#!/bin/bash
# Figures for each target gene, from the fast-path gene BAMs (step 03):
#   * scripts/gene_figure.py: one STRIPE / Nature-style figure per gene (183 mm wide):
#       a sashimi per sample + gene model, b IGV-style reads per sample,
#       c 5'->3' coverage along the canonical transcript, d spliced / unspliced / full-length reads
#   * scripts/igv_plot.py: an IGV-style page per gene and sample over the whole gene,
#       once with all reads and once with spliced reads only
# -> ${RESULTS_DIR}/figures/<GENE>/
#
# Light work (minutes, a few GB): LOCAL=1 runs it here instead of submitting it.
# Usage: bash 09_gene_figures.sh
set -euo pipefail

_D="${TEQ_ANALYSIS_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
source "${_D}/lib.sh"

preflight() {
  for gene in "${TARGET_GENES[@]}"; do
    compgen -G "${FAST_DIR}/bam/${gene}/*.bam" >/dev/null || {
      echo "no gene BAMs in ${FAST_DIR}/bam/${gene}: run 03_bait_and_align.sh first" >&2; exit 1; }
  done
}
if [[ -z "${AFTER:-}" ]]; then preflight; fi

if submit_self epi_figures 2 8G 02:00:00 "${_D}/09_gene_figures.sh" "$@"; then
  exit 0
fi
preflight

activate_env "${ANALYSIS_ENV}"
GTF="${REF_DIR}/${GTF_NAME}"
GENOME="${REF_DIR}/${FASTA_NAME}"

for gene in "${TARGET_GENES[@]}"; do
  out="${RESULTS_DIR}/figures/${gene}"
  mkdir -p "${out}/igv"
  echo "[figures] ${gene}: gene figure"
  python "${ANALYSIS_DIR}/scripts/gene_figure.py" --gene "${gene}" --gtf "${GTF}" \
    --sample-sheet "${SAMPLE_SHEET}" --bam-pattern "${FAST_DIR}/bam/{gene}/{sample}.bam" --outdir "${out}"

  # the gene locus, and an annotation cut to it (igv_plot.py reads its GTF in full)
  read -r chrom start end < <(awk -F'\t' -v g="${gene}" \
    '$3 == "gene" && index($9, "gene_name \"" g "\";") {print $1, $4, $5; exit}' "${GTF}")
  region="${chrom}:$(( start > 500 ? start - 500 : 1 ))-$(( end + 500 ))"
  awk -F'\t' -v c="${chrom}" -v s="${start}" -v e="${end}" '$1 == c && $4 <= e + 500 && $5 >= s - 500' \
    "${GTF}" > "${out}/igv/${gene}_region.gtf"

  for sample in $(sample_ids); do
    bam="${FAST_DIR}/bam/${gene}/${sample}.bam"
    [[ -s "${bam}" ]] || continue
    echo "[figures] ${gene}: IGV view ${sample}"
    python "${ANALYSIS_DIR}/scripts/igv_plot.py" --bam "${bam}" --region "${region}" --fasta "${GENOME}" \
      --gtf "${out}/igv/${gene}_region.gtf" --outfile "${out}/igv/${gene}_${sample}_all_reads.pdf" \
      --title "${gene}  ${sample}  ${region}  (all reads)"
    python "${ANALYSIS_DIR}/scripts/igv_plot.py" --bam "${bam}" --region "${region}" --fasta "${GENOME}" \
      --gtf "${out}/igv/${gene}_region.gtf" --outfile "${out}/igv/${gene}_${sample}_spliced_reads.pdf" \
      --spliced-only --title "${gene}  ${sample}  ${region}  (spliced reads only)"
  done
done
echo "[figures] done: ${RESULTS_DIR}/figures"
