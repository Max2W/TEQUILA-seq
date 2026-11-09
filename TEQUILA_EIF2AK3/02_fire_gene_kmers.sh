#!/bin/bash
# FAST PATH, step A1 (one SLURM job):
#   1. build the FIRE v0.5.0 index from GRCh38 + GENCODE v49, if it is not there yet.
#      FIRE_INDEX_MODE=target (default, config.sh): only the TARGET_GENES, seconds.
#      If FIRE rejects that small input, it falls back to every protein-coding
#      gene (~20 min). FIRE_INDEX_MODE=all builds the full one directly.
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
  local idx
  for idx in "${FIRE_INDEX_TARGET}" "${FIRE_INDEX_ALL}"; do
    if [[ -e "${idx}" && ! -f "${idx}/meta.json" ]]; then
      echo "${idx} exists but is not a finished FIRE index (no meta.json)." >&2
      echo "Remove it and run this step again:  rm -rf ${idx}" >&2
      exit 1
    fi
  done
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

if submit_self eif2ak3_fire_kmers 8 48G 06:00:00 "${_D}/02_fire_gene_kmers.sh" "$@"; then
  exit 0
fi

THREADS="${SLURM_CPUS_PER_TASK:-8}"

# fire_build <outdir> <gtf> <panel>
# FIRE builds in a temporary directory and moves it into place only when it is
# complete and validated, so a finished index always has meta.json.
fire_build() {
  echo "[fire] building FIRE ${FIRE_TAG} index at $1 (log: $1.build.log)"
  (cd "${FIRE_DIR}" && pixi run python -m xinglab_fire build \
    --gtf "$2" --genome "${GENOME}" --panel "$3" --outdir "$1" --n-jobs "${THREADS}")
  [[ -f "$1/meta.json" ]]
}

# 1. FIRE index
if [[ "${FIRE_INDEX_MODE}" != "all" && ! -f "${FIRE_INDEX_TARGET}/meta.json" ]]; then
  inputs="${FIRE_INDEX_TARGET}_inputs"
  mkdir -p "${inputs}"
  : > "${inputs}/genes.gtf"
  printf '%s\n' "${TARGET_GENES[@]}" > "${inputs}/panel.txt"
  for gene in "${TARGET_GENES[@]}"; do
    n=$(grep -c -F "gene_name \"${gene}\";" "${GTF}" || true)
    (( n > 0 )) || { echo "[fire] ${gene} not in ${GTF}" >&2; exit 1; }
    grep -F "gene_name \"${gene}\";" "${GTF}" >> "${inputs}/genes.gtf"
  done
  echo "[fire] target-only annotation: $(grep -c $'\ttranscript\t' "${inputs}/genes.gtf") transcripts of ${TARGET_GENES[*]}"
  if ! fire_build "${FIRE_INDEX_TARGET}" "${inputs}/genes.gtf" "${inputs}/panel.txt"; then
    echo "[fire] WARNING: the target-only build failed (see ${FIRE_INDEX_TARGET}.build.log);" \
         "falling back to the full protein-coding index" >&2
    rm -rf "${FIRE_INDEX_TARGET}"
    FIRE_INDEX="${FIRE_INDEX_ALL}"
  fi
fi
if [[ "${FIRE_INDEX}" == "${FIRE_INDEX_ALL}" && ! -f "${FIRE_INDEX_ALL}/meta.json" ]]; then
  fire_build "${FIRE_INDEX_ALL}" "${GTF}" "${FIRE_PANEL}" || { echo "[fire] index build failed" >&2; exit 1; }
fi
echo "[fire] using index ${FIRE_INDEX}"
mkdir -p "${FAST_DIR}/kmers"
echo "${FIRE_INDEX}" > "${FAST_DIR}/kmers/fire_index_used.txt"

# 2. gene-unique k-mers
for gene in "${TARGET_GENES[@]}"; do
  out="${FAST_DIR}/kmers/${gene}"
  mkdir -p "${out}"
  pixi run --manifest-path "${FIRE_DIR}/pixi.toml" python \
    "${ANALYSIS_DIR}/scripts/fire_gene_kmers.py" \
    --index "${FIRE_INDEX}" --gene "${gene}" --outdir "${out}"
done
echo "[fire] k-mers ready in ${FAST_DIR}/kmers"
