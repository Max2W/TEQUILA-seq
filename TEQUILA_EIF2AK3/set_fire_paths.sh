#!/bin/bash
# Find (or take) the FIRE v0.5.0 checkout and FIRE index, check them, and save
# them to local_settings.sh (read by config.sh, ignored by git).
#
# Usage: bash set_fire_paths.sh                       search your folders automatically
#        bash set_fire_paths.sh <fire_dir> <index_dir> use these two paths
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT="${HERE}/local_settings.sh"
# the index format of FIRE v0.5.0 (src/xinglab_fire/index/index_nova.py: INDEX_NOVA_VERSION)
WANT_INDEX_VERSION="2026.11.01"
SEARCH_ROOTS=("${HOME}" "/mnt/isilon/lin_lab_share/Feng")

is_fire_dir() { [[ -f "$1/pixi.toml" && -d "$1/src/xinglab_fire" ]]; }
is_index() { [[ -f "$1/meta.json" && -d "$1/genedb" && -d "$1/tx_fasta" && -d "$1/sshash_bundle.pcg" ]]; }
index_version() { grep -o '"version": *"[^"]*"' "$1/meta.json" | head -1 | sed 's/.*"\([^"]*\)"$/\1/'; }

if (( $# == 2 )); then
  fire_dirs=("$(cd "$1" && pwd)")
  index_dirs=("$(cd "$2" && pwd)")
else
  echo "Searching ${SEARCH_ROOTS[*]} for a FIRE checkout and FIRE indexes (can take a minute)..."
  fire_dirs=()
  index_dirs=()
  for root in "${SEARCH_ROOTS[@]}"; do
    [[ -d "${root}" ]] || continue
    while IFS= read -r f; do
      d="$(dirname "${f}")"
      if is_fire_dir "${d}"; then fire_dirs+=("${d}"); fi
    done < <(find "${root}" -maxdepth 5 -name pixi.toml -not -path '*/.pixi/*' 2>/dev/null)
    while IFS= read -r f; do
      d="$(dirname "${f}")"
      if is_index "${d}"; then index_dirs+=("${d}"); fi
    done < <(find "${root}" -maxdepth 6 -name meta.json 2>/dev/null)
  done
fi

echo
echo "FIRE checkouts found: ${#fire_dirs[@]}"
for d in "${fire_dirs[@]}"; do
  v="$(git -C "${d}" describe --tags 2>/dev/null || echo 'unknown version')"
  echo "  ${d}   (${v})"
done
echo "FIRE indexes found: ${#index_dirs[@]}"
for d in "${index_dirs[@]}"; do
  echo "  ${d}   (index format $(index_version "${d}"))"
done
echo

if (( ${#fire_dirs[@]} != 1 || ${#index_dirs[@]} != 1 )); then
  echo "Found zero or several of one of them, so nothing was saved."
  echo "Copy the right two paths from the lists above and run:"
  echo "  bash set_fire_paths.sh <FIRE checkout path> <FIRE index path>"
  exit 1
fi

fire_dir="${fire_dirs[0]}"
index_dir="${index_dirs[0]}"
is_fire_dir "${fire_dir}" || { echo "ERROR: ${fire_dir} is not a FIRE checkout (no pixi.toml / src/xinglab_fire)"; exit 1; }
is_index "${index_dir}" || { echo "ERROR: ${index_dir} is not a FIRE index (no meta.json / genedb / sshash_bundle.pcg)"; exit 1; }
got="$(index_version "${index_dir}")"
if [[ "${got}" != "${WANT_INDEX_VERSION}" ]]; then
  echo "ERROR: ${index_dir} has index format ${got}; FIRE v0.5.0 needs ${WANT_INDEX_VERSION}."
  echo "Rebuild it with FIRE v0.5.0 (fire build), or pick another index."
  exit 1
fi

# keep any other settings already in the file (e.g. SBATCH_EXTRA), replace the FIRE lines
touch "${OUT}"
grep -v -e '^FIRE_DIR=' -e '^FIRE_INDEX=' -e '^# FIRE paths written by' "${OUT}" > "${OUT}.tmp" || true
{
  cat "${OUT}.tmp"
  echo "# FIRE paths written by set_fire_paths.sh on $(date)"
  echo "FIRE_DIR=\"${fire_dir}\""
  echo "FIRE_INDEX=\"${index_dir}\""
} > "${OUT}"
rm -f "${OUT}.tmp"
echo "Saved to ${OUT}:"
cat "${OUT}"
