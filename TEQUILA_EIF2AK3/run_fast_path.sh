#!/bin/bash
# Submit the whole EIF2AK3 fast path as a chain of SLURM jobs, each starting
# only after the previous one succeeds:
#   01 merge fastq -> 02 FIRE k-mers -> 03 bait + align -> 04 ESPRESSO (targeted)
#                                                       -> 05 AMALGAM (targeted)
#                                                       -> 08 reports
# Run 00_setup.sh first. Skip step 01 if it already ran: SKIP_FASTQ=1
set -euo pipefail
_D="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${_D}/lib.sh"
require_vars FIRE_DIR FIRE_INDEX

after="${FASTQ_JOB:-}"   # set by run_all.sh, which submits step 01 once for both paths
if [[ -z "${after}" && "${SKIP_FASTQ:-0}" != "1" ]]; then
  after=$(job_id "$(bash "${_D}/01_prepare_fastq.sh")")
fi
# the FIRE index / k-mers do not need the fastq files, so they start right away
j_kmers=$(job_id "$(AFTER="" bash "${_D}/02_fire_gene_kmers.sh")")
j_bait=$(job_id "$(AFTER="${j_kmers}${after:+:${after}}" bash "${_D}/03_bait_and_align.sh")")
j_esp=$(job_id "$(AFTER="${j_bait}" bash "${_D}/04_espresso_targeted.sh")")
j_am=$(job_id "$(AFTER="${j_bait}" bash "${_D}/05_amalgam.sh" fast)")
j_rep=$(job_id "$(AFTER="${j_esp}:${j_am}" bash "${_D}/08_isoform_report.sh")")
echo "fast path submitted; reports job ${j_rep}. Watch: squeue -u \$USER"
