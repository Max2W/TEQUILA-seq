#!/bin/bash
# Merge all fastq chunks of one barcode from one flowcell into a single
# fastq.gz, and print a read-count line:
#   flowcell  barcode  sample_id  n_files  n_reads  n_bases
#
# gzip chunks are joined byte for byte (a series of gzip members is a valid
# gzip file), so nothing is recompressed: the merge runs at disk speed. Only
# uncompressed chunks, if any, are compressed. Reads are counted in the same
# pass. A finished output (with its .counts file) is not redone.
#
# Usage: merge_barcode_fastq.sh <fc_label> <barcode_dir> <sample_id|-> <out.fastq.gz|->
#   sample_id "-" / out "-" : only count reads (barcodes not in the sample sheet)
set -euo pipefail

fc="$1"
bc_dir="$2"
sample="$3"
out="$4"
barcode="$(basename "${bc_dir}")"

if [[ ! -d "${bc_dir}" ]]; then
  printf '%s\t%s\t%s\t0\t0\t0\n' "${fc}" "${barcode}" "${sample}"
  exit 0
fi

shopt -s nullglob
gz_files=("${bc_dir}"/*.fastq.gz "${bc_dir}"/*.fq.gz)
plain_files=("${bc_dir}"/*.fastq "${bc_dir}"/*.fq)
n_files=$(( ${#gz_files[@]} + ${#plain_files[@]} ))
count_awk='NR%4==2 {n++; b+=length($0)} END {printf "%d\t%d\n", n, b}'

# compressed stream of all chunks of this barcode
gz_stream() {
  if (( ${#gz_files[@]} )); then cat -- "${gz_files[@]}"; fi
  if (( ${#plain_files[@]} )); then cat -- "${plain_files[@]}" | gzip -c -1; fi
}

if (( n_files == 0 )); then
  counts=$'0\t0'
elif [[ "${out}" == "-" ]]; then
  counts=$(gz_stream | gzip -dc | awk "${count_awk}")
elif [[ -s "${out}" && -s "${out}.counts" ]]; then
  counts=$(cat "${out}.counts")          # finished in an earlier run
elif [[ -s "${out}" ]]; then
  # finished by an earlier version of this script (it renamed .tmp only when
  # complete): count its reads instead of merging again
  counts=$(gzip -dc "${out}" | awk "${count_awk}")
  printf '%s\n' "${counts}" > "${out}.counts"
else
  mkdir -p "$(dirname "${out}")"
  tmp="${out}.tmp"
  counts=$(gz_stream | tee "${tmp}" | gzip -dc | awk "${count_awk}")
  printf '%s\n' "${counts}" > "${out}.counts"
  mv "${tmp}" "${out}"
fi

printf '%s\t%s\t%s\t%s\t%s\n' "${fc}" "${barcode}" "${sample}" "${n_files}" "${counts}"
