#!/bin/bash
# Merge all fastq chunks of one barcode from one flowcell into a single
# fastq.gz, and print a read-count line:
#   flowcell  barcode  sample_id  n_files  n_reads  n_bases
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

stream() {
  # Concatenated gzip members are a valid gzip file, but decompress so that
  # plain and gzipped chunks can be mixed and reads counted in one pass.
  if (( ${#gz_files[@]} )); then zcat -f -- "${gz_files[@]}"; fi
  if (( ${#plain_files[@]} )); then cat -- "${plain_files[@]}"; fi
}

count_awk='NR%4==2 {n++; b+=length($0)} END {printf "%d\t%d\n", n, b}'

if (( n_files == 0 )); then
  counts=$'0\t0'
elif [[ "${out}" == "-" ]]; then
  counts=$(stream | awk "${count_awk}")
else
  mkdir -p "$(dirname "${out}")"
  tmp="${out}.tmp"
  fifo="${out}.fifo"
  rm -f "${fifo}"
  mkfifo "${fifo}"
  gzip -c < "${fifo}" > "${tmp}" &
  gz_pid=$!
  counts=$(stream | tee "${fifo}" | awk "${count_awk}")
  wait "${gz_pid}"
  rm -f "${fifo}"
  mv "${tmp}" "${out}"
fi

printf '%s\t%s\t%s\t%s\t%s\n' "${fc}" "${barcode}" "${sample}" "${n_files}" "${counts}"
