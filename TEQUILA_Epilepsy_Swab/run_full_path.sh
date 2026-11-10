#!/bin/bash
# Submit the standard (slow) TEQUILA-seq path as a chain of SLURM jobs:
#   01 merge fastq -> 06 Snakemake controller (minimap2 + ESPRESSO on all reads)
#                  -> 07 QC on the full alignments -> 05 AMALGAM genome-wide
#                  -> 08 reports (also picks up the fast-path results)
# Run 00_setup.sh first. Skip step 01 if it already ran: SKIP_FASTQ=1
# Leave out the genome-wide AMALGAM run: SKIP_AMALGAM_FULL=1
set -euo pipefail
_D="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${_D}/lib.sh"

after="${FASTQ_JOB:-}"   # set by run_all.sh, which submits step 01 once for both paths
if [[ -z "${after}" && "${SKIP_FASTQ:-0}" != "1" ]]; then
  after=$(job_id "$(bash "${_D}/01_prepare_fastq.sh")")
fi
j_full=$(job_id "$(AFTER="${after}" bash "${_D}/06_espresso_full.sh")")
j_qc=$(job_id "$(AFTER="${j_full}" bash "${_D}/07_qc_full.sh")")
last="${j_qc}"
if [[ "${SKIP_AMALGAM_FULL:-0}" != "1" ]]; then
  last=$(job_id "$(AFTER="${j_qc}" bash "${_D}/05_amalgam.sh" full)")
fi
j_rep=$(job_id "$(AFTER="${last}" bash "${_D}/08_isoform_report.sh")")
echo "full path submitted; reports job ${j_rep}. Watch: squeue -u \$USER"
