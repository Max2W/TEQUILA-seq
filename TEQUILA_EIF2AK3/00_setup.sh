#!/bin/bash
# Step 0 (run once, on a login node because it downloads files):
#   * creates the conda environments (TEQUILA-seq pipeline, analysis, AMALGAM)
#   * clones AMALGAM into ${AMALGAM_DIR}
#   * sets up FIRE v0.5.0 in ${FIRE_DIR}: pixi environment + compiled extensions
#     (the FIRE index itself is built later, as a SLURM job, by 02_fire_gene_kmers.sh)
#   * downloads GRCh38 + GENCODE v49 into ${REF_DIR}
#   * creates the Snakemake working directory ${PIPELINE_DIR}
# Safe to run again: finished parts are skipped, a half-made conda env is remade.
#
# Usage:  bash 00_setup.sh            # everything
#         bash 00_setup.sh envs       # only conda envs
#         bash 00_setup.sh amalgam    # only AMALGAM clone + env
#         bash 00_setup.sh fire       # only FIRE (clone, pixi install, compile)
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

# make_env <prefix> <requirements file> <check binary> [extra create args...]
# A finished env gets a marker file; an env without it (a failed earlier try) is
# removed and created again.
make_env() {
  local prefix="$1" reqs="$2" check="$3"
  shift 3
  if [[ -f "${prefix}/.teq_env_ok" && -x "${prefix}/bin/${check}" ]]; then
    echo "[setup] env exists: ${prefix}"
    return 0
  fi
  if [[ -e "${prefix}" ]]; then
    echo "[setup] removing incomplete env ${prefix}"
    rm -rf "${prefix}"
  fi
  echo "[setup] creating env ${prefix} from ${reqs}"
  mkdir -p "$(dirname "${prefix}")"
  ${SOLVER} create -y --prefix "${prefix}" "$@" --file "${reqs}"
  [[ -x "${prefix}/bin/${check}" ]] || { echo "[setup] ${prefix} has no ${check}" >&2; return 1; }
  touch "${prefix}/.teq_env_ok"
}

setup_envs() {
  # the TEQUILA-seq pipeline env (python 3.8, snakemake 5.31, ESPRESSO deps)
  if [[ -x "${PIPELINE_ENV}/bin/snakemake" && ! -f "${PIPELINE_ENV}/.teq_env_ok" ]]; then
    touch "${PIPELINE_ENV}/.teq_env_ok"   # made by an earlier version of this script
  fi
  make_env "${PIPELINE_ENV}" "${CODE_DIR}/conda_requirements.txt" snakemake -c conda-forge -c bioconda
  # conda-forge + bioconda only, so packages from 'defaults' are not mixed in
  make_env "${ANALYSIS_ENV}" "${ANALYSIS_DIR}/envs/analysis_env.txt" minimap2 \
    --override-channels -c conda-forge -c bioconda
}

setup_amalgam() {
  if [[ ! -d "${AMALGAM_DIR}/.git" ]]; then
    mkdir -p "$(dirname "${AMALGAM_DIR}")"
    git clone "${AMALGAM_REPO}" "${AMALGAM_DIR}"
  fi
  git -C "${AMALGAM_DIR}" fetch --quiet origin || true
  git -C "${AMALGAM_DIR}" checkout --quiet "${AMALGAM_COMMIT}"
  echo "[setup] AMALGAM at $(git -C "${AMALGAM_DIR}" rev-parse HEAD)"
  # AMALGAM's list pins biopython=1.78 together with python=3.11, but biopython 1.78
  # was never built for python 3.11, so that list cannot be installed. AMALGAM's
  # README asks for biopython 1.78 or newer; only that pin is relaxed, the rest is kept.
  local reqs="${PROJECT_DIR}/envs/amalgam_requirements.txt"
  mkdir -p "$(dirname "${reqs}")"
  sed 's/^biopython=1\.78$/biopython>=1.78/' "${AMALGAM_DIR}/conda_requirements.txt" > "${reqs}"
  echo "[setup] AMALGAM requirements used:"
  diff "${AMALGAM_DIR}/conda_requirements.txt" "${reqs}" || true
  make_env "${AMALGAM_ENV}" "${reqs}" stringtie \
    --override-channels -c conda-forge -c bioconda
}

setup_fire() {
  # pixi is FIRE's only supported installer. Use the one on PATH, else install
  # a private copy under tools/pixi (your shell settings are not changed).
  if ! command -v pixi >/dev/null 2>&1; then
    if [[ ! -x "${TEQ_PIXI_DIR}/bin/pixi" ]]; then
      echo "[setup] installing pixi into ${TEQ_PIXI_DIR}"
      mkdir -p "${TEQ_PIXI_DIR}"
      curl -fsSL https://pixi.sh/install.sh | PIXI_HOME="${TEQ_PIXI_DIR}" PIXI_NO_PATH_UPDATE=1 bash
    fi
    export PATH="${TEQ_PIXI_DIR}/bin:${PATH}"
  fi
  echo "[setup] pixi: $(command -v pixi) ($(pixi --version))"

  if [[ ! -d "${FIRE_DIR}/.git" ]]; then
    mkdir -p "$(dirname "${FIRE_DIR}")"
    # no password prompt here (this may run in the background); see the README if it fails
    if ! GIT_TERMINAL_PROMPT=0 git clone --branch "${FIRE_TAG}" --depth 1 "${FIRE_REPO}" "${FIRE_DIR}"; then
      echo "[setup] could not download FIRE without a password. Run this once, then run 00_setup.sh again:" >&2
      echo "  git clone --branch ${FIRE_TAG} --depth 1 git@github.com:Max2W/fire.git ${FIRE_DIR}" >&2
      return 1
    fi
  fi
  local tag
  tag="$(git -C "${FIRE_DIR}" describe --tags --exact-match 2>/dev/null || echo none)"
  [[ "${tag}" == "${FIRE_TAG}" ]] || { echo "[setup] ${FIRE_DIR} is at '${tag}', not ${FIRE_TAG}" >&2; return 1; }

  if [[ ! -f "${FIRE_DIR}/.teq_build_ok" ]]; then
    echo "[setup] FIRE: pixi install + compiling extensions (several minutes)"
    (cd "${FIRE_DIR}" && pixi install && pixi run build-ext)
    touch "${FIRE_DIR}/.teq_build_ok"
  fi
  (cd "${FIRE_DIR}" && pixi run python -m xinglab_fire version)
  [[ -f "${FIRE_PANEL}" ]] || { echo "[setup] missing ${FIRE_PANEL}" >&2; return 1; }
  echo "[setup] FIRE ${FIRE_TAG} ready in ${FIRE_DIR}; its index is built by 02_fire_gene_kmers.sh"
}

download_refs() {
  mkdir -p "${REF_DIR}"
  if [[ ! -s "${REF_DIR}/${FASTA_NAME}" ]]; then
    echo "[setup] downloading ${FASTA_URL}"
    curl -L --fail "${FASTA_URL}" -o "${REF_DIR}/${FASTA_NAME}.gz"
    gunzip -f "${REF_DIR}/${FASTA_NAME}.gz"
  fi
  if [[ ! -s "${REF_DIR}/${GTF_NAME}" ]]; then
    echo "[setup] downloading ${GTF_URL}"
    curl -L --fail "${GTF_URL}" -o "${REF_DIR}/${GTF_NAME}.gz"
    gunzip -f "${REF_DIR}/${GTF_NAME}.gz"
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
  fire) setup_fire ;;
  refs) download_refs ;;
  pipeline) setup_pipeline_dir ;;
  all) setup_envs; setup_amalgam; setup_fire; download_refs; setup_pipeline_dir
       echo "=== SETUP COMPLETE ===" ;;
  *) echo "unknown target: $1" >&2; exit 1 ;;
esac
