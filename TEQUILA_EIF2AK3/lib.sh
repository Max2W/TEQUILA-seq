#!/bin/bash
# Shared helpers for the step scripts. Sourced, not executed.

# Under sbatch the script is copied to a spool dir, so the analysis folder is
# passed in through TEQ_ANALYSIS_DIR when a step submits itself.
if [[ -n "${TEQ_ANALYSIS_DIR:-}" ]]; then
  _TEQ_DIR="${TEQ_ANALYSIS_DIR}"
else
  _TEQ_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi
source "${_TEQ_DIR}/config.sh"

# submit_self <job name> <cpus> <mem> <time> <script> [args...]
# Re-submits the calling script to SLURM unless already inside a job or
# LOCAL=1 is set (run in the current shell, e.g. inside an interactive job).
# Prints the job id. Set AFTER=<job id[:job id...]> to start only after those
# jobs finish successfully (used by run_fast_path.sh).
submit_self() {
  local name="$1" cpus="$2" mem="$3" time="$4" script="$5"
  shift 5
  if [[ "${TEQ_IN_JOB:-0}" == "1" || "${LOCAL:-0}" == "1" ]]; then
    return 1
  fi
  mkdir -p "${LOG_DIR}"
  local dep=()
  if [[ -n "${AFTER:-}" ]]; then dep=(--dependency="afterok:${AFTER}"); fi
  local job_id
  # shellcheck disable=SC2086
  job_id=$(sbatch --parsable --job-name="${name}" --cpus-per-task="${cpus}" --mem="${mem}" \
    --time="${time}" -o "${LOG_DIR}/${name}_%j.out" -e "${LOG_DIR}/${name}_%j.err" \
    "${dep[@]}" ${SBATCH_EXTRA} --export=ALL,TEQ_ANALYSIS_DIR="${_TEQ_DIR}",TEQ_IN_JOB=1 \
    "${script}" "$@")
  echo "${job_id%%;*}"
  echo "submitted ${name}: job ${job_id%%;*} (log: ${LOG_DIR}/${name}_${job_id%%;*}.out)" >&2
  return 0
}

# require_vars VAR... : stop with a clear message if any is empty
require_vars() {
  local v missing=0
  for v in "$@"; do
    if [[ -z "${!v:-}" ]]; then
      echo "ERROR: ${v} is not set (edit config.sh or export it)" >&2
      missing=1
    fi
  done
  (( missing == 0 )) || exit 1
}

# sample_ids : the sample_id column of the sample sheet
sample_ids() {
  awk -F'\t' 'NR > 1 && $1 != "" {print $1}' "${SAMPLE_SHEET}"
}

activate_env() {
  if ! command -v conda >/dev/null 2>&1; then
    source "${HOME}/.bashrc"
  fi
  eval "$(conda shell.bash hook)"
  conda activate "$1"
}
