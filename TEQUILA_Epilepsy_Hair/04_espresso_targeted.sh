#!/bin/bash
# FAST PATH, step A3: ESPRESSO (the TEQUILA-seq Snakefile) on the baited,
# gene-locus BAMs only, with the annotation cut down to the gene. With a few
# thousand reads per sample this takes minutes instead of days.
# Output: ${FAST_DIR}/espresso/<GENE>/espresso_out/work_dir/samples_N2_R0_{abundance.esp,updated.gtf}
#
# Usage: bash 04_espresso_targeted.sh   (submits itself; LOCAL=1 to run here)
set -euo pipefail

_D="${TEQ_ANALYSIS_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
source "${_D}/lib.sh"

preflight() {
  [[ -x "${PIPELINE_ENV}/bin/snakemake" ]] || { echo "pipeline env missing: run 00_setup.sh envs" >&2; exit 1; }
  for gene in "${TARGET_GENES[@]}"; do
    [[ -d "${FAST_DIR}/bam/${gene}" ]] || { echo "no ${FAST_DIR}/bam/${gene}: run 03_bait_and_align.sh" >&2; exit 1; }
  done
}
# inputs from earlier steps may not exist yet when this job waits on them (AFTER=...)
if [[ -z "${AFTER:-}" ]]; then preflight; fi

if submit_self epi_espresso_fast 8 32G 12:00:00 "${_D}/04_espresso_targeted.sh" "$@"; then
  exit 0
fi
preflight

THREADS="${SLURM_CPUS_PER_TASK:-8}"
GTF="${REF_DIR}/${GTF_NAME}"

for gene in "${TARGET_GENES[@]}"; do
  work="${FAST_DIR}/espresso/${gene}"
  mkdir -p "${work}/references"
  cd "${work}"
  ln -sfn "${CODE_DIR}/Snakefile" Snakefile
  ln -sfn "${CODE_DIR}/scripts" scripts
  ln -sfn "${CODE_DIR}/visualization" visualization
  ln -sfn "${REF_DIR}/${FASTA_NAME}" "references/${FASTA_NAME}"
  ln -sfn "${REF_DIR}/${FASTA_NAME}.fai" "references/${FASTA_NAME}.fai"
  gene_gtf="${gene}.${GTF_NAME}"
  grep -F "gene_name \"${gene}\";" "${GTF}" > "references/${gene_gtf}"
  [[ -s "references/${gene_gtf}" ]] || { echo "${gene} not in ${GTF}" >&2; exit 1; }

  # leave out samples with no reads on the gene (ESPRESSO stops on an empty SAM)
  skip=()
  activate_env "${ANALYSIS_ENV}"
  for sample in $(sample_ids); do
    bam="${FAST_DIR}/bam/${gene}/${sample}.bam"
    if [[ ! -s "${bam}" ]] || (( $(samtools view -c "${bam}") == 0 )); then skip+=("${sample}"); fi
  done
  conda deactivate
  if (( ${#skip[@]} )); then echo "[${gene}] no reads, left out: ${skip[*]}" >&2; fi

  python3 "${ANALYSIS_DIR}/scripts/make_snakemake_config.py" \
    --sample-sheet "${SAMPLE_SHEET}" --bam-dir "${FAST_DIR}/bam/${gene}" \
    --skip-samples "$(IFS=,; echo "${skip[*]:-}")" \
    --code-dir "${CODE_DIR}" --conda-wrapper "${CONDA_WRAPPER}" --gtf-name "${gene_gtf}" --fasta-name "${FASTA_NAME}" \
    --threads "${THREADS}" --out snakemake_config.yaml --index-out input_index.tsv

  # local execution inside this job (no cluster profile): the inputs are small
  "${CONDA_WRAPPER}" snakemake --cores "${THREADS}" --latency-wait 30 --printshellcmds \
    > "snakemake.log" 2>&1 || { echo "[${gene}] ESPRESSO failed, see ${work}/snakemake.log" >&2; exit 1; }
  echo "[${gene}] done: ${work}/espresso_out/work_dir/samples_N2_R0_abundance.esp"
done
