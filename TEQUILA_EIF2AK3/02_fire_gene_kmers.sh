#!/bin/bash
# FAST PATH, step A1 (one SLURM job):
#   1. build the FIRE v0.5.0 index from GRCh38 + GENCODE v49, if it is not there yet
#      (FIRE's measured peak is ~19 GB, wall time 1.5-2 h; skipped once built)
#   2. extract the gene-unique 31-mers of each target gene from it
#      -> ${FAST_DIR}/kmers/<GENE>/<GENE>.anchor_kmers.fa
#
# Usage: bash 02_fire_gene_kmers.sh           (submits itself to SLURM)
#        bash 02_fire_gene_kmers.sh --check   (check FIRE + references, then exit)
set -euo pipefail

_D="${TEQ_ANALYSIS_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
source "${_D}/lib.sh"

GENOME="${REF_DIR}/${FASTA_NAME}"
GTF="${REF_DIR}/${GTF_NAME}"

preflight() {
  command -v pixi >/dev/null || { echo "pixi not found: run 00_setup.sh fire" >&2; exit 1; }
  [[ -f "${FIRE_DIR}/.teq_build_ok" ]] || { echo "FIRE not set up in ${FIRE_DIR}: run 00_setup.sh fire" >&2; exit 1; }
  [[ -f "${FIRE_PANEL}" ]] || { echo "missing ${FIRE_PANEL}" >&2; exit 1; }
  [[ -s "${GENOME}" && -s "${GENOME}.fai" && -s "${GTF}" ]] || {
    echo "references missing in ${REF_DIR}: run 00_setup.sh refs" >&2; exit 1; }
  if [[ -e "${FIRE_INDEX}" && ! -f "${FIRE_INDEX}/meta.json" ]]; then
    echo "${FIRE_INDEX} exists but is not a finished FIRE index (no meta.json)." >&2
    echo "Remove it and run this step again:  rm -rf ${FIRE_INDEX}" >&2
    exit 1
  fi
}
preflight

if [[ "${1:-}" == "--check" ]]; then
  (cd "${FIRE_DIR}" && pixi run python -m xinglab_fire version)
  if [[ -f "${FIRE_INDEX}/meta.json" ]]; then
    echo "FIRE index present: ${FIRE_INDEX}"
    for gene in "${TARGET_GENES[@]}"; do
      pixi run --manifest-path "${FIRE_DIR}/pixi.toml" python \
        "${ANALYSIS_DIR}/scripts/fire_gene_kmers.py" --index "${FIRE_INDEX}" --gene "${gene}" \
        --outdir "${FAST_DIR}/kmers/${gene}" --check
    done
  else
    echo "FIRE index not built yet; the job will build it at ${FIRE_INDEX}"
  fi
  exit 0
fi

if submit_self eif2ak3_fire_kmers 16 64G 1-00:00:00 "${_D}/02_fire_gene_kmers.sh" "$@"; then
  exit 0
fi

THREADS="${SLURM_CPUS_PER_TASK:-8}"

# 1. FIRE index (FIRE builds in a temporary directory and moves it into place only
#    when it is complete and validated, so a finished index always has meta.json)
if [[ ! -f "${FIRE_INDEX}/meta.json" ]]; then
  echo "[fire] building FIRE ${FIRE_TAG} index at ${FIRE_INDEX}"
  echo "[fire] build log: ${FIRE_INDEX}.build.log"
  (cd "${FIRE_DIR}" && pixi run python -m xinglab_fire build \
    --gtf "${GTF}" --genome "${GENOME}" --panel "${FIRE_PANEL}" \
    --outdir "${FIRE_INDEX}" --n-jobs "${THREADS}")
  [[ -f "${FIRE_INDEX}/meta.json" ]] || { echo "[fire] index build did not finish" >&2; exit 1; }
else
  echo "[fire] using existing index ${FIRE_INDEX}"
fi

# 2. gene-unique k-mers
for gene in "${TARGET_GENES[@]}"; do
  out="${FAST_DIR}/kmers/${gene}"
  mkdir -p "${out}"
  pixi run --manifest-path "${FIRE_DIR}/pixi.toml" python \
    "${ANALYSIS_DIR}/scripts/fire_gene_kmers.py" \
    --index "${FIRE_INDEX}" --gene "${gene}" --outdir "${out}"
done
echo "[fire] k-mers ready in ${FAST_DIR}/kmers"
