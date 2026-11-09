#!/bin/bash
# AMALGAM isoform discovery + quantification (https://github.com/RNA-ROB/amalgam),
# following its README: StringTie per sample -> GffCompare against GENCODE ->
# Build_Transcriptome.py -> sort/bgzip/tabix -> Quantify_Transcripts.py per sample.
# (Annotate_ORF.py is skipped: ~5 h and not needed for isoform proportions.)
#
# Usage: bash 05_amalgam.sh fast   per target gene, on the baited gene BAMs (minutes)
#                                  -> ${FAST_DIR}/amalgam/<GENE>/
#        bash 05_amalgam.sh full   genome-wide, on ${QC_DIR}/bam/<sample>.bam from
#                                  07_qc_full.sh (Quantify is the slow part: hours)
#                                  -> ${PROJECT_DIR}/amalgam_full/
# Each run ends with amalgam_counts_matrix.tsv (ESPRESSO-like layout) used by 09.
set -euo pipefail

_D="${TEQ_ANALYSIS_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
source "${_D}/lib.sh"

MODE="${1:-fast}"
case "${MODE}" in
  fast) JOB_CPUS=8; JOB_MEM=32G; JOB_TIME=12:00:00 ;;
  full) JOB_CPUS=16; JOB_MEM=96G; JOB_TIME=5-00:00:00 ;;
  *) echo "usage: $0 fast|full" >&2; exit 1 ;;
esac
[[ -x "${AMALGAM_ENV}/bin/stringtie" ]] || { echo "AMALGAM env missing: run 00_setup.sh amalgam" >&2; exit 1; }
[[ -f "${AMALGAM_DIR}/scripts/Quantify_Transcripts.py" ]] || { echo "AMALGAM code missing in ${AMALGAM_DIR}" >&2; exit 1; }

if submit_self "eif2ak3_amalgam_${MODE}" "${JOB_CPUS}" "${JOB_MEM}" "${JOB_TIME}" "${_D}/05_amalgam.sh" "$@"; then
  exit 0
fi

activate_env "${AMALGAM_ENV}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"
GENOME="${REF_DIR}/${FASTA_NAME}"
GTF="${REF_DIR}/${GTF_NAME}"
echo "AMALGAM commit: $(git -C "${AMALGAM_DIR}" rev-parse HEAD)"

# run_amalgam <out dir> <reference GTF> <bam dir>
run_amalgam() {
  local out="$1" ref_gtf="$2" bam_dir="$3"
  mkdir -p "${out}/stringtie" "${out}/counts"
  local samples=() sample bam n
  for sample in $(sample_ids); do
    bam="${bam_dir}/${sample}.bam"
    if [[ ! -s "${bam}" ]]; then echo "missing ${bam}, ${sample} left out" >&2; continue; fi
    n=$(samtools view -c -F 0x904 "${bam}")
    if (( n == 0 )); then echo "${sample}: no reads, left out" >&2; continue; fi
    samples+=("${sample}")
  done
  (( ${#samples[@]} )) || { echo "no samples with reads for ${out}" >&2; return 1; }

  # 1. StringTie, per sample (README: default settings)
  for sample in "${samples[@]}"; do
    local st="${out}/stringtie/${sample}.gtf"
    [[ -s "${st}" ]] || stringtie -o "${st}" -p "${THREADS}" "${bam_dir}/${sample}.bam"
  done
  # 2. GffCompare: GENCODE first, then every sample
  echo "${ref_gtf}" > "${out}/gtf_list.txt"
  for sample in "${samples[@]}"; do echo "${out}/stringtie/${sample}.gtf" >> "${out}/gtf_list.txt"; done
  gffcompare -i "${out}/gtf_list.txt" -T -o "${out}/gffcmp"
  # 3. Build the transcriptome
  python "${AMALGAM_DIR}/scripts/Build_Transcriptome.py" -i "${out}/gffcmp" -g "${ref_gtf}" \
    -f "${GENOME}" -o "${out}/amalgam.transcripts.gtf"
  sort -k1,1V -k4,4g -k5,5g "${out}/amalgam.transcripts.gtf" | bgzip > "${out}/amalgam.transcripts.sorted.gtf.gz"
  tabix -f -p gff "${out}/amalgam.transcripts.sorted.gtf.gz"
  # 4. Quantify, per sample
  for sample in "${samples[@]}"; do
    local cnt="${out}/counts/${sample}.counts.tsv"
    [[ -s "${cnt}" ]] || python "${AMALGAM_DIR}/scripts/Quantify_Transcripts.py" \
      -i "${bam_dir}/${sample}.bam" -g "${out}/amalgam.transcripts.sorted.gtf.gz" \
      -t "${THREADS}" -o "${cnt}"
  done
  # 5. one matrix for the downstream report
  python "${ANALYSIS_DIR}/scripts/merge_amalgam_counts.py" \
    --counts-dir "${out}/counts" --sample-sheet "${SAMPLE_SHEET}" --ref-gtf "${ref_gtf}" \
    --out "${out}/amalgam_counts_matrix.tsv"
  echo "[amalgam] ${out}/amalgam_counts_matrix.tsv"
}

if [[ "${MODE}" == "fast" ]]; then
  for gene in "${TARGET_GENES[@]}"; do
    out="${FAST_DIR}/amalgam/${gene}"
    mkdir -p "${out}"
    gene_gtf="${out}/${gene}.${GTF_NAME}"
    grep -F "gene_name \"${gene}\";" "${GTF}" > "${gene_gtf}"
    run_amalgam "${out}" "${gene_gtf}" "${FAST_DIR}/bam/${gene}"
  done
else
  [[ -d "${QC_DIR}/bam" ]] || { echo "no ${QC_DIR}/bam: run 07_qc_full.sh first" >&2; exit 1; }
  run_amalgam "${PROJECT_DIR}/amalgam_full" "${GTF}" "${QC_DIR}/bam"
fi
