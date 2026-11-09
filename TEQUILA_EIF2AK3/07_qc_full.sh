#!/bin/bash
# QC on the full alignments from 06_espresso_full.sh (can run as soon as the
# alignment/<sample>/<i>/aligned.sorted.sam files exist; ESPRESSO need not be done).
#   * per flowcell input and per sample (flowcells merged): reads, mapped,
#     reads on EIF2AK3 / EIF2AK2, read length and quality, top genes by reads
#   * the repo's qc_analysis/plot_on_target_rates.py figure
# Outputs: ${QC_DIR}/bam/<sample>.bam (merged, indexed; also used by 05_amalgam.sh full)
#          ${QC_DIR}/per_input_qc.tsv, per_sample_qc.tsv, top_genes.tsv, *.pdf
set -euo pipefail

_D="${TEQ_ANALYSIS_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
source "${_D}/lib.sh"
preflight() {
  [[ -s "${PIPELINE_DIR}/input_index.tsv" ]] || { echo "run 06_espresso_full.sh first" >&2; exit 1; }
}
if [[ -z "${AFTER:-}" ]]; then preflight; fi

if submit_self eif2ak3_qc 8 32G 12:00:00 "${_D}/07_qc_full.sh" "$@"; then
  exit 0
fi
preflight

activate_env "${ANALYSIS_ENV}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"
GTF="${REF_DIR}/${GTF_NAME}"
mkdir -p "${QC_DIR}/bam/per_input"

# SAM -> BAM per flowcell input, then merge per sample
while IFS=$'\t' read -r sample idx fc _path; do
  [[ "${sample}" == "sample_id" ]] && continue
  sam="${PIPELINE_DIR}/alignment/${sample}/${idx}/aligned.sorted.sam"
  bam="${QC_DIR}/bam/per_input/${sample}.${fc}.bam"
  [[ -s "${sam}" ]] || { echo "missing ${sam} (alignment not finished?)" >&2; exit 1; }
  if [[ ! -s "${bam}.bai" ]]; then
    samtools view -@ "${THREADS}" -b -o "${bam}" "${sam}"
    samtools index "${bam}"
  fi
done < "${PIPELINE_DIR}/input_index.tsv"

for sample in $(sample_ids); do
  shopt -s nullglob
  parts=("${QC_DIR}/bam/per_input/${sample}".*.bam)
  out="${QC_DIR}/bam/${sample}.bam"
  (( ${#parts[@]} )) || continue
  if [[ ! -s "${out}.bai" ]]; then
    if (( ${#parts[@]} == 1 )); then cp "${parts[0]}" "${out}"; else samtools merge -f -@ "${THREADS}" "${out}" "${parts[@]}"; fi
    samtools index "${out}"
  fi
done

python "${ANALYSIS_DIR}/scripts/qc_alignments.py" \
  --input-index "${PIPELINE_DIR}/input_index.tsv" --bam-dir "${QC_DIR}/bam" \
  --sample-sheet "${SAMPLE_SHEET}" --gtf "${GTF}" --target-genes "$(IFS=,; echo "${TARGET_GENES[*]}")" \
  --threads "${THREADS}" --outdir "${QC_DIR}"

# the repo's standard on-target figure (targets = EIF2AK3 + EIF2AK2 gene bodies)
rates_in="${QC_DIR}/rates_in.tsv"
: > "${rates_in}"
while IFS=$'\t' read -r sample _no _donor _sp cell _trt _cond _bc; do
  [[ "${sample}" == "sample_id" ]] && continue
  [[ -s "${QC_DIR}/bam/${sample}.bam" ]] || continue
  printf '%s\t%s\t%s\t%s\n' "${sample}" "${QC_DIR}/bam/${sample}.bam" "${QC_DIR}/target_genes.bed" "${cell}" >> "${rates_in}"
done < "${SAMPLE_SHEET}"
python "${CODE_DIR}/qc_analysis/plot_on_target_rates.py" --mapping-file "${rates_in}" \
  --outprefix "${QC_DIR}/rates" --title "TEQUILA-seq EIF2AK3" --threads "${THREADS}"
echo "[qc] done: ${QC_DIR}"
