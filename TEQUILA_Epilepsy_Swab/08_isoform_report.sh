#!/bin/bash
# Isoform reports for every target gene and every quantification that exists:
#   ESPRESSO_fast  ${FAST_DIR}/espresso/<GENE>/      (04_espresso_targeted.sh)
#   AMALGAM_fast   ${FAST_DIR}/amalgam/<GENE>/       (05_amalgam.sh fast)
#   ESPRESSO_full  ${PIPELINE_DIR}/espresso_out/     (06_espresso_full.sh)
#   AMALGAM_full   ${PROJECT_DIR}/amalgam_full/      (05_amalgam.sh full)
# -> ${RESULTS_DIR}/<GENE>/<METHOD>/  (tables, figures, <GENE>_<METHOD>_report.md)
#    control genes (CONTROL_GENES in config.sh) go to ${RESULTS_DIR}/internal_control/<GENE>/
# -> ${RESULTS_DIR}/<GENE>/method_comparison*  (when two or more methods are done)
# -> ${RESULTS_DIR}/fragmentation/{fast,full}/  detection + fragmentation per gene and sample
# Missing runs are skipped, so this can be re-run as more of them finish.
#
# Usage: bash 08_isoform_report.sh   (submits itself; LOCAL=1 to run here)
set -euo pipefail

_D="${TEQ_ANALYSIS_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
source "${_D}/lib.sh"

if submit_self epi_report 2 16G 02:00:00 "${_D}/08_isoform_report.sh" "$@"; then
  exit 0
fi

activate_env "${ANALYSIS_ENV}"
GTF="${REF_DIR}/${GTF_NAME}"

report() {  # report <gene> <method> <matrix> <isoform gtf>
  local gene="$1" method="$2" matrix="$3" iso_gtf="$4"
  if [[ ! -s "${matrix}" || ! -s "${iso_gtf}" ]]; then
    echo "[${gene}/${method}] not available yet, skipped" >&2
    return 0
  fi
  local out="$(gene_results_dir "${gene}")/${method}"
  mkdir -p "${out}"
  python "${ANALYSIS_DIR}/scripts/isoform_report.py" \
    --matrix "${matrix}" --isoform-gtf "${iso_gtf}" --ref-gtf "${GTF}" \
    --sample-sheet "${SAMPLE_SHEET}" --gene "${gene}" --method "${method}" \
    --reference-condition "${REFERENCE_CONDITION}" \
    --find-as "${CODE_DIR}/scripts/Find_AS_Events.py" --outdir "${out}" \
    || echo "[${gene}/${method}] report failed" >&2
}

gene_results_dir() {  # results/<GENE>, or results/internal_control/<GENE> for a control
  local g
  for g in "${CONTROL_GENES[@]}"; do
    if [[ "$1" == "${g}" ]]; then echo "${RESULTS_DIR}/internal_control/$1"; return; fi
  done
  echo "${RESULTS_DIR}/$1"
}

for gene in "${TARGET_GENES[@]}"; do
  gdir="$(gene_results_dir "${gene}")"
  wd="espresso_out/work_dir"
  report "${gene}" ESPRESSO_fast "${FAST_DIR}/espresso/${gene}/${wd}/samples_N2_R0_abundance.esp" \
    "${FAST_DIR}/espresso/${gene}/${wd}/samples_N2_R0_updated.gtf"
  report "${gene}" AMALGAM_fast "${FAST_DIR}/amalgam/${gene}/amalgam_counts_matrix.tsv" \
    "${FAST_DIR}/amalgam/${gene}/amalgam.transcripts.gtf"
  report "${gene}" ESPRESSO_full "${PIPELINE_DIR}/${wd}/samples_N2_R0_abundance.esp" \
    "${PIPELINE_DIR}/${wd}/samples_N2_R0_updated.gtf"
  report "${gene}" AMALGAM_full "${PROJECT_DIR}/amalgam_full/amalgam_counts_matrix.tsv" \
    "${PROJECT_DIR}/amalgam_full/amalgam.transcripts.gtf"

  runs=()
  for method in ESPRESSO_fast AMALGAM_fast ESPRESSO_full AMALGAM_full; do
    [[ -d "${gdir}/${method}" ]] && runs+=(--run "${method}=${gdir}/${method}")
  done
  if (( ${#runs[@]} >= 4 )); then  # two or more runs (each --run is two words)
    python "${ANALYSIS_DIR}/scripts/compare_methods.py" --gene "${gene}" "${runs[@]}" \
      --out-prefix "${gdir}/method_comparison" || true
  fi
done
# detection and fragmentation (does a gene show up; how much of each transcript do reads cover)
genes_csv="$(IFS=,; echo "${TARGET_GENES[*]}")"
if compgen -G "${FAST_DIR}/bam/${TARGET_GENES[0]}/*.bam" >/dev/null; then
  python "${ANALYSIS_DIR}/scripts/fragment_report.py" --genes "${genes_csv}" --gtf "${GTF}" \
    --sample-sheet "${SAMPLE_SHEET}" --bam-pattern "${FAST_DIR}/bam/{gene}/{sample}.bam" \
    --read-counts "${QC_DIR}/fastq_read_counts.tsv" --label "fast path" \
    --outdir "${RESULTS_DIR}/fragmentation/fast" || echo "[report] fragmentation (fast) failed" >&2
fi
if compgen -G "${QC_DIR}/bam/*.bam" >/dev/null; then
  python "${ANALYSIS_DIR}/scripts/fragment_report.py" --genes "${genes_csv}" --gtf "${GTF}" \
    --sample-sheet "${SAMPLE_SHEET}" --bam-pattern "${QC_DIR}/bam/{sample}.bam" \
    --read-counts "${QC_DIR}/fastq_read_counts.tsv" --label "all reads" \
    --outdir "${RESULTS_DIR}/fragmentation/full" || echo "[report] fragmentation (full) failed" >&2
fi
echo "[report] results in ${RESULTS_DIR}"
