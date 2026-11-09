#!/bin/bash
#SBATCH --job-name=eif2ak3_fastq
# Step 1: merge the fastq_pass chunks of each sample barcode, per flowcell,
# into ${FASTQ_DIR}/<sample>/<flowcell>.fastq.gz (gzip chunks are joined, not
# recompressed) and count their reads -> ${QC_DIR}/fastq_read_counts.tsv.
# COUNT_OTHER_BARCODES=1 also counts the barcodes not in the sample sheet.
# Resumable: finished sample files are not merged again.
#
# Usage: bash 01_prepare_fastq.sh          (submits itself to SLURM)
#        LOCAL=1 bash 01_prepare_fastq.sh  (runs in the current shell)
set -euo pipefail

_D="${TEQ_ANALYSIS_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
source "${_D}/lib.sh"

if submit_self eif2ak3_fastq 8 8G 12:00:00 "${_D}/01_prepare_fastq.sh" "$@"; then
  exit 0
fi

THREADS="${SLURM_CPUS_PER_TASK:-4}"
mkdir -p "${FASTQ_DIR}" "${QC_DIR}"

# barcode -> sample_id from the sample sheet
declare -A BC2SAMPLE
while IFS=$'\t' read -r sample_id _no _donor _sp _ct _trt _cond barcode; do
  [[ "${sample_id}" == "sample_id" ]] && continue
  BC2SAMPLE["${barcode}"]="${sample_id}"
done < "${SAMPLE_SHEET}"

tasks="$(mktemp)"
for entry in "${FLOWCELLS[@]}"; do
  fc="${entry%%:*}"
  run="${entry#*:}"
  pass_dir="${RAW_DIR}/${run}/fastq_pass"
  if [[ ! -d "${pass_dir}" ]]; then
    echo "missing ${pass_dir}" >&2
    exit 1
  fi
  # the sample barcodes; with COUNT_OTHER_BARCODES=1 also the other barcode
  # folders and "unclassified" (read counts only; slow on a large run)
  for bc_dir in "${pass_dir}"/barcode* "${pass_dir}"/unclassified; do
    [[ -d "${bc_dir}" ]] || continue
    bc="$(basename "${bc_dir}")"
    sample="${BC2SAMPLE[${bc}]:--}"
    if [[ "${sample}" == "-" ]]; then
      [[ "${COUNT_OTHER_BARCODES:-0}" == "1" ]] || continue
      out="-"
    else
      out="${FASTQ_DIR}/${sample}/${fc}.fastq.gz"
    fi
    printf '%s\t%s\t%s\t%s\n' "${fc}" "${bc_dir}" "${sample}" "${out}" >> "${tasks}"
  done
  # sample barcodes that are missing on this flowcell
  for bc in "${!BC2SAMPLE[@]}"; do
    if [[ ! -d "${pass_dir}/${bc}" ]]; then
      echo "WARNING: ${fc} has no ${bc} (${BC2SAMPLE[${bc}]})" >&2
    fi
  done
done

counts="${QC_DIR}/fastq_read_counts.tsv"
printf 'flowcell\tbarcode\tsample_id\tn_files\tn_reads\tn_bases\n' > "${counts}"
# one merge job per line of ${tasks}, THREADS at a time
xargs -P "${THREADS}" -I{} \
  bash -c 'IFS=$'"'"'\t'"'"' read -r fc dir sample out <<< "$1"; bash "$2" "$fc" "$dir" "$sample" "$out"' \
  _ {} "${ANALYSIS_DIR}/scripts/merge_barcode_fastq.sh" < "${tasks}" \
  | sort -k1,1 -k2,2V >> "${counts}"
rm -f "${tasks}"

echo "Read counts written to ${counts}"
(column -t "${counts}" 2>/dev/null || cat "${counts}")
