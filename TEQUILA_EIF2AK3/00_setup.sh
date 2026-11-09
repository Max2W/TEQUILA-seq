#!/bin/bash
# Step 0 (run once, on a login node because it downloads files):
#   * creates the conda environments (TEQUILA-seq pipeline, analysis, AMALGAM)
#   * clones AMALGAM into ${AMALGAM_DIR}
#   * downloads GRCh38 + GENCODE v44 into ${REF_DIR}
#   * creates the Snakemake working directory ${PIPELINE_DIR}
#
# Usage:  bash 00_setup.sh            # everything
#         bash 00_setup.sh envs       # only conda envs
#         bash 00_setup.sh amalgam    # only AMALGAM clone + env
#         bash 00_setup.sh fire       # only check the FIRE install and index
#         bash 00_setup.sh refs       # only references
#         bash 00_setup.sh pipeline   # only the pipeline working dir
set -euo pipefail

_TEQ_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${_TEQ_DIR}/config.sh"

# conda activate needs the shell hook
if ! command -v conda >/dev/null 2>&1; then
  source "${HOME}/.bashrc"
fi
eval "$(conda shell.bash hook)"

# Prefer mamba when available (much faster solves)
SOLVER=conda
if command -v mamba >/dev/null 2>&1; then SOLVER=mamba; fi

setup_envs() {
  if [[ ! -x "${PIPELINE_ENV}/bin/snakemake" ]]; then
    echo "[setup] creating pipeline env at ${PIPELINE_ENV}"
    ${SOLVER} create -y --prefix "${PIPELINE_ENV}" -c conda-forge -c bioconda \
      --file "${CODE_DIR}/conda_requirements.txt"
  else
    echo "[setup] pipeline env exists: ${PIPELINE_ENV}"
  fi

  if [[ ! -x "${ANALYSIS_ENV}/bin/python" ]]; then
    echo "[setup] creating analysis env at ${ANALYSIS_ENV}"
    mkdir -p "$(dirname "${ANALYSIS_ENV}")"
    ${SOLVER} create -y --prefix "${ANALYSIS_ENV}" -c conda-forge -c bioconda \
      --file "${ANALYSIS_DIR}/envs/analysis_env.txt"
  else
    echo "[setup] analysis env exists: ${ANALYSIS_ENV}"
  fi
}

setup_amalgam() {
  if [[ ! -d "${AMALGAM_DIR}/.git" ]]; then
    mkdir -p "$(dirname "${AMALGAM_DIR}")"
    git clone "${AMALGAM_REPO}" "${AMALGAM_DIR}"
  fi
  git -C "${AMALGAM_DIR}" fetch --quiet origin || true
  git -C "${AMALGAM_DIR}" checkout --quiet "${AMALGAM_COMMIT}"
  echo "[setup] AMALGAM at $(git -C "${AMALGAM_DIR}" rev-parse HEAD)"
  if [[ ! -x "${AMALGAM_ENV}/bin/stringtie" ]]; then
    echo "[setup] creating AMALGAM env at ${AMALGAM_ENV}"
    mkdir -p "$(dirname "${AMALGAM_ENV}")"
    ${SOLVER} create -y --prefix "${AMALGAM_ENV}" -c conda-forge -c bioconda \
      --file "${AMALGAM_DIR}/conda_requirements.txt"
  else
    echo "[setup] AMALGAM env exists: ${AMALGAM_ENV}"
  fi
}

check_fire() {
  if [[ -z "${FIRE_DIR}" || -z "${FIRE_INDEX}" ]]; then
    echo "[setup] FIRE_DIR / FIRE_INDEX not set in config.sh: skipping the FIRE check" \
         "(needed before 02_fire_gene_kmers.sh)"
    return 0
  fi
  [[ -f "${FIRE_INDEX}/meta.json" ]] || { echo "no FIRE index at ${FIRE_INDEX}" >&2; return 1; }
  (cd "${FIRE_DIR}" && pixi run python -m xinglab_fire version)
}

download_refs() {
  mkdir -p "${REF_DIR}"
  if [[ ! -s "${REF_DIR}/${FASTA_NAME}" ]]; then
    echo "[setup] downloading ${FASTA_URL}"
    curl -L --fail "${FASTA_URL}" -o "${REF_DIR}/${FASTA_NAME}.gz"
    gunzip "${REF_DIR}/${FASTA_NAME}.gz"
  fi
  if [[ ! -s "${REF_DIR}/${GTF_NAME}" ]]; then
    echo "[setup] downloading ${GTF_URL}"
    curl -L --fail "${GTF_URL}" -o "${REF_DIR}/${GTF_NAME}.gz"
    gunzip "${REF_DIR}/${GTF_NAME}.gz"
  fi
  if [[ ! -s "${REF_DIR}/${FASTA_NAME}.fai" ]]; then
    conda activate "${ANALYSIS_ENV}"
    samtools faidx "${REF_DIR}/${FASTA_NAME}"
    conda deactivate
  fi
  echo "[setup] references ready in ${REF_DIR}"
}

setup_pipeline_dir() {
  mkdir -p "${PIPELINE_DIR}" "${LOG_DIR}" "${FASTQ_DIR}" "${QC_DIR}" "${RESULTS_DIR}"
  cd "${PIPELINE_DIR}"
  # The Snakefile uses paths relative to its working directory, so the code is
  # linked in. The profile is copied because it writes job usage files.
  ln -sfn "${CODE_DIR}/Snakefile" Snakefile
  ln -sfn "${CODE_DIR}/scripts" scripts
  ln -sfn "${CODE_DIR}/visualization" visualization
  ln -sfn "${REF_DIR}" references
  rm -rf snakemake_profile
  cp -r "${CODE_DIR}/snakemake_profile" snakemake_profile
  echo "[setup] pipeline working dir ready: ${PIPELINE_DIR}"
}

case "${1:-all}" in
  envs) setup_envs ;;
  amalgam) setup_amalgam ;;
  fire) check_fire ;;
  refs) download_refs ;;
  pipeline) setup_pipeline_dir ;;
  all) setup_envs; setup_amalgam; download_refs; setup_pipeline_dir; check_fire
       echo "=== SETUP COMPLETE ===" ;;
  *) echo "unknown target: $1" >&2; exit 1 ;;
esac
