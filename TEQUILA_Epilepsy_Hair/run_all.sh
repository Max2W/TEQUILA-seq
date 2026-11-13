#!/bin/bash
# Submit everything: step 01 once, then the fast path (FIRE k-mer bait ->
# targeted ESPRESSO + AMALGAM) and the full TEQUILA-seq path, both waiting on 01.
# Results appear in ${RESULTS_DIR} as each path finishes (fast: hours; full: days).
set -euo pipefail
_D="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${_D}/lib.sh"
require_vars FIRE_DIR FIRE_INDEX

fastq_job=$(job_id "$(bash "${_D}/01_prepare_fastq.sh")")
FASTQ_JOB="${fastq_job}" bash "${_D}/run_fast_path.sh"
FASTQ_JOB="${fastq_job}" bash "${_D}/run_full_path.sh"
