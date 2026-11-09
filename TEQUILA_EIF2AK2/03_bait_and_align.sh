#!/bin/bash
# FAST PATH, step A2: for every sample
#   1. bait reads that carry EIF2AK3 / EIF2AK2 anchor 31-mers (all flowcells)
#   2. align only the baited reads to GRCh38 with the TEQUILA-seq minimap2 settings
#   3. cut one BAM per gene (reads overlapping the gene) for ESPRESSO / AMALGAM
# Outputs in ${FAST_DIR}:
#   bait/<sample>/<GENE>.fastq.gz, .hits.tsv, bait_summary.tsv
#   bam/all/<sample>.bam                       all baited reads, genome-wide
#   bam/<GENE>/<sample>.bam                    reads on the gene locus
#   bait_align_stats.tsv                       baited -> aligned -> on-locus, per sample and gene
#
# Usage: bash 03_bait_and_align.sh   (submits itself; LOCAL=1 to run here)
# Resumable: a sample whose gene BAMs exist is skipped (FORCE=1 to redo).
set -euo pipefail

_D="${TEQ_ANALYSIS_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
source "${_D}/lib.sh"

preflight() {
  for gene in "${TARGET_GENES[@]}"; do
    f="${FAST_DIR}/kmers/${gene}/${gene}.anchor_kmers.fa"
    [[ -s "${f}" ]] || { echo "missing ${f}: run 02_fire_gene_kmers.sh first" >&2; exit 1; }
  done
  [[ -s "${REF_DIR}/${FASTA_NAME}" && -s "${REF_DIR}/${GTF_NAME}" ]] || {
    echo "references missing in ${REF_DIR}: run 00_setup.sh refs" >&2; exit 1; }
}
# inputs from earlier steps may not exist yet when this job waits on them (AFTER=...)
if [[ -z "${AFTER:-}" ]]; then preflight; fi

if submit_self eif2ak2_bait_align 16 64G 1-00:00:00 "${_D}/03_bait_and_align.sh" "$@"; then
  exit 0
fi
preflight

activate_env "${ANALYSIS_ENV}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"
GENOME="${REF_DIR}/${FASTA_NAME}"
GTF="${REF_DIR}/${GTF_NAME}"
MMI="${REF_DIR}/${FASTA_NAME%.fa}.splice.k14w4.mmi"
JUNC_BED="${REF_DIR}/${GTF_NAME}.bed"
mkdir -p "${FAST_DIR}/bait" "${FAST_DIR}/bam/all"

# --- one-time reference prep -------------------------------------------------
if [[ ! -s "${JUNC_BED}" ]]; then
  paftools.js gff2bed "${GTF}" > "${JUNC_BED}.tmp" && mv "${JUNC_BED}.tmp" "${JUNC_BED}"
fi
if [[ ! -s "${MMI}" ]]; then
  echo "[align] building minimap2 index ${MMI}"
  minimap2 -x splice -k 14 -w 4 -t "${THREADS}" -d "${MMI}.tmp" "${GENOME}" && mv "${MMI}.tmp" "${MMI}"
fi

# gene loci from the GTF: GENE<TAB>chr:start-end
LOCI="${FAST_DIR}/target_loci.tsv"
awk -F'\t' -v genes="$(IFS=,; echo "${TARGET_GENES[*]}")" '
  BEGIN { n = split(genes, g, ","); for (i = 1; i <= n; i++) want[g[i]] = 1 }
  $3 == "gene" && match($9, /gene_name "[^"]+"/) {
    name = substr($9, RSTART + 11, RLENGTH - 12)
    if (name in want) printf "%s\t%s:%d-%d\t%s\n", name, $1, $4, $5, $7
  }' "${GTF}" > "${LOCI}"
for gene in "${TARGET_GENES[@]}"; do
  grep -q -P "^${gene}\t" "${LOCI}" || { echo "${gene} not found in ${GTF}" >&2; exit 1; }
  mkdir -p "${FAST_DIR}/bam/${gene}"
done

anchor_args=()
for gene in "${TARGET_GENES[@]}"; do
  anchor_args+=(--anchors "${gene}=${FAST_DIR}/kmers/${gene}/${gene}.anchor_kmers.fa")
done

# --- 1. bait, all samples in parallel (one core each) -------------------------
bait_one() {
  local sample="$1"
  local out="${FAST_DIR}/bait/${sample}"
  if [[ -s "${out}/bait_summary.tsv" && "${FORCE:-0}" != "1" ]]; then return 0; fi
  shopt -s nullglob
  local fq=("${FASTQ_DIR}/${sample}"/*.fastq.gz)
  if (( ${#fq[@]} == 0 )); then echo "no fastq for ${sample}" >&2; return 1; fi
  python "${ANALYSIS_DIR}/scripts/kmer_bait_reads.py" "${anchor_args[@]}" \
    --fastq "${fq[@]}" --out-dir "${out}" --min-hits "${BAIT_MIN_HITS}" \
    2> "${out}.bait.log" || { echo "bait failed for ${sample}, see ${out}.bait.log" >&2; return 1; }
}
export -f bait_one
export FAST_DIR FASTQ_DIR ANALYSIS_DIR BAIT_MIN_HITS FORCE
export ANCHOR_ARGS="${anchor_args[*]}"
# shellcheck disable=SC2016
sample_ids | xargs -P "${THREADS}" -I{} bash -c 'anchor_args=(${ANCHOR_ARGS}); bait_one "$1"' _ {}

# --- 2/3. align each sample's baited reads, then cut per-gene BAMs ------------
STATS="${FAST_DIR}/bait_align_stats.tsv"
printf 'sample_id\tgene\treads_baited\treads_aligned_primary\treads_on_locus\tfraction_on_locus\n' > "${STATS}.tmp"
for sample in $(sample_ids); do
  all_bam="${FAST_DIR}/bam/all/${sample}.bam"
  if [[ ! -s "${all_bam}.bai" || "${FORCE:-0}" == "1" ]]; then
    # union of the reads baited for any gene (a read is written once)
    for gene in "${TARGET_GENES[@]}"; do zcat "${FAST_DIR}/bait/${sample}/${gene}.fastq.gz"; done \
      | awk 'NR % 4 == 1 { id = $1; keep = !(id in seen); seen[id] = 1 } keep' \
      | minimap2 -a -x splice -ub --secondary=no --junc-bed "${JUNC_BED}" -t "${THREADS}" "${MMI}" - \
          2> "${FAST_DIR}/bam/all/${sample}.minimap2.log" \
      | samtools sort -@ 4 -o "${all_bam}.tmp" -
    mv "${all_bam}.tmp" "${all_bam}"
    samtools index "${all_bam}"
  fi
  for gene in "${TARGET_GENES[@]}"; do
    region="$(awk -F'\t' -v g="${gene}" '$1 == g {print $2; exit}' "${LOCI}")"
    gene_bam="${FAST_DIR}/bam/${gene}/${sample}.bam"
    samtools view -b -F 0x904 -o "${gene_bam}" "${all_bam}" "${region}"
    samtools index "${gene_bam}"
    baited=$(( $(zcat "${FAST_DIR}/bait/${sample}/${gene}.fastq.gz" | wc -l) / 4 ))
    # primary alignments of this gene's baited reads, and how many land on its locus
    ids="${FAST_DIR}/bait/${sample}/${gene}.ids"
    awk 'NR > 1 {print $1}' "${FAST_DIR}/bait/${sample}/${gene}.hits.tsv" > "${ids}"
    aligned=$(samtools view -c -F 0x904 -N "${ids}" "${all_bam}")
    on_locus=$(samtools view -c -F 0x904 -N "${ids}" "${all_bam}" "${region}")
    frac=$(awk -v a="${on_locus}" -v b="${baited}" 'BEGIN { printf "%.4f", (b > 0 ? a / b : 0) }')
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "${sample}" "${gene}" "${baited}" "${aligned}" "${on_locus}" "${frac}" >> "${STATS}.tmp"
  done
done
mv "${STATS}.tmp" "${STATS}"

# combined bait summary
first=1
for sample in $(sample_ids); do
  f="${FAST_DIR}/bait/${sample}/bait_summary.tsv"
  if (( first )); then awk -v s=sample_id 'NR == 1 {print s "\t" $0}' "${f}"; first=0; fi
  awk -v s="${sample}" 'NR > 1 {print s "\t" $0}' "${f}"
done > "${FAST_DIR}/bait_summary_all.tsv"

echo "[done] ${STATS}"
(column -t "${STATS}" 2>/dev/null || cat "${STATS}")
