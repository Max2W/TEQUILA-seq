#!/bin/bash
# Wait until the minimap2 genome index shared with the EIF2AK3 project is
# finished (it is built by that project's bait-and-align job), then submit
# all jobs of this project with run_all.sh. Checks every 5 minutes.
set -euo pipefail
_D="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${_D}/lib.sh"
MMI="${REF_DIR}/${FASTA_NAME%.fa}.splice.k14w4.mmi"
while [[ ! -s "${MMI}" || -e "${MMI}.tmp" ]]; do
  echo "$(date '+%F %T') waiting for the shared genome index ${MMI}"
  sleep 300
done
echo "$(date '+%F %T') genome index ready; submitting jobs"
bash "${_D}/run_all.sh"
