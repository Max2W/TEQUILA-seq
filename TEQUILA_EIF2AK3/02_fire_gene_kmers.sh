#!/bin/bash
# FAST PATH, step A1: extract the gene-unique 31-mers of each target gene from
# your FIRE v0.5.0 index -> ${FAST_DIR}/kmers/<GENE>/<GENE>.anchor_kmers.fa
#
# Needs FIRE_DIR and FIRE_INDEX in config.sh.
# Usage: bash 02_fire_gene_kmers.sh           (submits itself to SLURM)
#        bash 02_fire_gene_kmers.sh --check   (open the index, resolve genes, exit)
set -euo pipefail

_D="${TEQ_ANALYSIS_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
source "${_D}/lib.sh"
require_vars FIRE_DIR FIRE_INDEX

CHECK=""
if [[ "${1:-}" == "--check" ]]; then CHECK="--check"; fi

# cheap checks first, on the login node
[[ -f "${FIRE_INDEX}/meta.json" ]] || { echo "no FIRE index (meta.json) at ${FIRE_INDEX}" >&2; exit 1; }
[[ -f "${FIRE_DIR}/pixi.toml" ]] || { echo "no pixi.toml in FIRE_DIR=${FIRE_DIR}" >&2; exit 1; }
command -v pixi >/dev/null || { echo "pixi not on PATH" >&2; exit 1; }

if [[ -z "${CHECK}" ]] && submit_self eif2ak3_kmers 2 16G 02:00:00 "${_D}/02_fire_gene_kmers.sh" "$@"; then
  exit 0
fi

for gene in "${TARGET_GENES[@]}"; do
  out="${FAST_DIR}/kmers/${gene}"
  mkdir -p "${out}"
  pixi run --manifest-path "${FIRE_DIR}/pixi.toml" python \
    "${ANALYSIS_DIR}/scripts/fire_gene_kmers.py" \
    --index "${FIRE_INDEX}" --gene "${gene}" --outdir "${out}" ${CHECK}
done
