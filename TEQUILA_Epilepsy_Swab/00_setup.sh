#!/bin/bash
# Step 0 (run once, on a login node; takes seconds).
# Code is in home; outputs go to ${PROJECT_DIR} and reference data is in ${REF_DIR}
# (both on scratch). Software (conda envs, FIRE, AMALGAM, pixi) is reused from
# ${SHARED_DIR}. This script only CHECKS those and then
#   * creates this project's folders on scratch
#   * writes ${CONDA_WRAPPER}, which runs a command in the shared pipeline env
#   * creates the Snakemake working directory ${PIPELINE_DIR}
set -euo pipefail

_TEQ_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${_TEQ_DIR}/config.sh"

missing=0
need() {  # need <path> <what>
  if [[ -e "$1" ]]; then
    echo "[setup] ok      $2: $1"
  else
    echo "[setup] MISSING $2: $1" >&2
    missing=1
  fi
}
need "${PIPELINE_ENV}/bin/snakemake" "TEQUILA-seq pipeline env"
need "${ANALYSIS_ENV}/bin/minimap2" "analysis env"
need "${AMALGAM_ENV}/bin/stringtie" "AMALGAM env"
need "${AMALGAM_DIR}/scripts/Quantify_Transcripts.py" "AMALGAM code"
need "${FIRE_DIR}/.teq_build_ok" "FIRE v0.5.0 (compiled)"
need "${REF_DIR}/${FASTA_NAME}.fai" "GRCh38 genome (indexed)"
need "${REF_DIR}/${GTF_NAME}" "GENCODE v${GENCODE_RELEASE} annotation"
if (( missing )); then
  echo "[setup] Something is missing (see MISSING above). Reference files belong in ${REF_DIR}" >&2
  echo "        (moved there from ${SHARED_DIR}/references); software in ${SHARED_DIR}." >&2
  exit 1
fi
command -v pixi >/dev/null || { echo "[setup] pixi not found (expected in ${TEQ_PIXI_DIR}/bin or on PATH)" >&2; exit 1; }

# the target genes must exist in the annotation
for gene in "${TARGET_GENES[@]}"; do
  if grep -q -F "gene_name \"${gene}\";" "${REF_DIR}/${GTF_NAME}"; then
    echo "[setup] ok      gene ${gene}"
  else
    echo "[setup] gene ${gene} (target_genes.txt) is not in GENCODE v${GENCODE_RELEASE}; check the spelling" >&2
    exit 1
  fi
done

# outputs go to scratch: check there is room (both paths need roughly 100-150 GB)
mkdir -p "${PROJECT_DIR}"
avail_gb=$(command df -P -BG "${PROJECT_DIR}" | awk 'NR == 2 {sub("G", "", $4); print $4}')
echo "[setup] ok      ${avail_gb} GB free where outputs go (${PROJECT_DIR})"
if (( avail_gb < ${MIN_FREE_GB:-200} )); then
  echo "[setup] less than ${MIN_FREE_GB:-200} GB free in ${PROJECT_DIR}; free some space first" >&2
  exit 1
fi

mkdir -p "${PIPELINE_DIR}" "${LOG_DIR}" "${FASTQ_DIR}" "${QC_DIR}" "${RESULTS_DIR}" "$(dirname "${CONDA_WRAPPER}")"

# runs its arguments inside the shared pipeline env (Snakemake and its rules)
cat > "${CONDA_WRAPPER}" <<EOF
#!/bin/bash
# written by TEQUILA_Epilepsy_Swab/00_setup.sh: run a command in ${PIPELINE_ENV}
if ! command -v conda >/dev/null 2>&1; then source "\${HOME}/.bashrc"; fi
eval "\$(conda shell.bash hook)"
conda activate "${PIPELINE_ENV}"
"\$@"
EOF
chmod +x "${CONDA_WRAPPER}"
"${CONDA_WRAPPER}" snakemake --version >/dev/null && echo "[setup] ok      ${CONDA_WRAPPER}"

# Snakemake working directory: the Snakefile uses paths relative to it
cd "${PIPELINE_DIR}"
ln -sfn "${CODE_DIR}/Snakefile" Snakefile
ln -sfn "${CODE_DIR}/scripts" scripts
ln -sfn "${CODE_DIR}/visualization" visualization
ln -sfn "${REF_DIR}" references
rm -rf snakemake_profile
cp -r "${CODE_DIR}/snakemake_profile" snakemake_profile
echo "[setup] pipeline working dir ready: ${PIPELINE_DIR}"
echo "=== SETUP COMPLETE ==="
