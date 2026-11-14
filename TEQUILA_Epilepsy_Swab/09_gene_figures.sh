#!/bin/bash
# Publication figures for the target genes, from the fast-path gene BAMs (step 03):
#   1. gene read filter (scripts/gene_read_filter.py): keep only reads that belong to the gene
#      -> ${FAST_DIR}/bam_filtered/<GENE>/<sample>.bam, ${RESULTS_DIR}/figures/read_filter.tsv
#   2. detection + fragmentation on the filtered reads (scripts/fragment_report.py)
#      -> ${RESULTS_DIR}/fragmentation/filtered/
#   3. per gene -> ${RESULTS_DIR}/figures/<GENE>/
#      - scripts/gene_figure.py (gene reads): sashimi per sample, IGV-style reads, gene-body coverage, reads and
#        full-length reads -> <GENE>_figure.{pdf,png}
#      - scripts/coverage_plot.py (gene reads): overlaid coverage lines on the intron-shortened canonical
#        transcript, and IGV-style coverage tracks per sample
#      - scripts/isoform_abundance_plot.py, for each isoform quantification that exists (ESPRESSO_fast,
#        AMALGAM_fast, AMALGAM_full, ESPRESSO_full): isoform CPM stacked bars by donor with the isoform
#        structures below, and as one panel
#   4. IGV-style pages per gene and sample (scripts/igv_plot.py), gene reads, only the gene of
#      interest annotated -> ${RESULTS_DIR}/figures/<GENE>/igv/
#   5. Figures 1-2 and the paired kit comparison (scripts/pub_figures.py) -> ${RESULTS_DIR}/figures/
#
# GENE_SOURCE=full takes the gene reads from the full path instead (the genome-wide BAMs of 07_qc_full.sh, every read
# aligned): the reads overlapping the gene are copied to ${FAST_DIR}/bam/<GENE>/<sample>.bam, so steps 02-05 are not
# needed. The k-mer baiting of the fast path recovers 99.6-100% of the gene reads found this way. Isoforms then come
# from AMALGAM full (genome-wide).
# Light work (minutes, a few GB): LOCAL=1 runs it here instead of submitting it.
# Usage: bash 09_gene_figures.sh          (GENE_SET=extra GENE_SOURCE=full bash 09_gene_figures.sh)
set -euo pipefail

_D="${TEQ_ANALYSIS_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
source "${_D}/lib.sh"

GENE_SOURCE="${GENE_SOURCE:-fast}"
preflight() {
  if [[ "${GENE_SOURCE}" == "full" ]]; then
    compgen -G "${QC_DIR}/bam/*.bam" >/dev/null || { echo "no ${QC_DIR}/bam/*.bam: run 07_qc_full.sh first" >&2; exit 1; }
  else
    for gene in "${TARGET_GENES[@]}"; do
      compgen -G "${FAST_DIR}/bam/${gene}/*.bam" >/dev/null || {
        echo "no gene BAMs in ${FAST_DIR}/bam/${gene}: run 03_bait_and_align.sh first (or GENE_SOURCE=full)" >&2
        exit 1; }
    done
  fi
  [[ -s "${QC_DIR}/fastq_read_counts.tsv" ]] || { echo "missing ${QC_DIR}/fastq_read_counts.tsv (step 01)" >&2; exit 1; }
}
if [[ -z "${AFTER:-}" ]]; then preflight; fi

if submit_self epi_figures 2 8G 02:00:00 "${_D}/09_gene_figures.sh" "$@"; then
  exit 0
fi
preflight

activate_env "${ANALYSIS_ENV}"
GTF="${REF_DIR}/${GTF_NAME}"
GENOME="${REF_DIR}/${FASTA_NAME}"
FIG="${RESULTS_DIR}/figures"
FILTERED="${FAST_DIR}/bam_filtered"
mkdir -p "${FIG}" "${FILTERED}"
genes_csv="$(IFS=,; echo "${TARGET_GENES[*]}")"

# 1. gene read filter (the locus GTF is cut once per gene: the full GTF is 3 GB)
STATS="${FIG}/read_filter.tsv"
first=1
for gene in "${TARGET_GENES[@]}"; do
  read -r chrom start end < <(awk -F'\t' -v g="${gene}" \
    '$3 == "gene" && index($9, "gene_name \"" g "\";") {print $1, $4, $5; exit}' "${GTF}")
  mkdir -p "${FIG}/${gene}/igv"
  awk -F'\t' -v c="${chrom}" -v s="${start}" -v e="${end}" '$1 == c && $4 <= e + 100000 && $5 >= s - 100000' \
    "${GTF}" > "${FIG}/${gene}/${gene}_locus.gtf"
  if [[ "${GENE_SOURCE}" == "full" ]]; then   # gene reads from the genome-wide BAMs of the full path
    mkdir -p "${FAST_DIR}/bam/${gene}"
    for sample in $(sample_ids); do
      src="${QC_DIR}/bam/${sample}.bam"
      [[ -s "${src}" ]] || continue
      samtools view -b -o "${FAST_DIR}/bam/${gene}/${sample}.bam" "${src}" "${chrom}:${start}-${end}"
      samtools index "${FAST_DIR}/bam/${gene}/${sample}.bam"
    done
    echo "[figures] ${gene}: reads of ${chrom}:${start}-${end} taken from the full path" >&2
  fi
  for sample in $(sample_ids); do
    bam="${FAST_DIR}/bam/${gene}/${sample}.bam"
    [[ -s "${bam}" ]] || continue
    hdr=(); if (( first )); then hdr=(--header); first=0; fi
    python "${ANALYSIS_DIR}/scripts/gene_read_filter.py" --gene "${gene}" --gtf "${FIG}/${gene}/${gene}_locus.gtf" \
      --bam "${bam}" --out-bam "${FILTERED}/${gene}/${sample}.bam" --sample "${sample}" "${hdr[@]}"
  done
done > "${STATS}.tmp"
mv "${STATS}.tmp" "${STATS}"
echo "[figures] read filter:"; (column -t "${STATS}" 2>/dev/null || cat "${STATS}")

# 2. detection + fragmentation on the filtered reads
python "${ANALYSIS_DIR}/scripts/fragment_report.py" --genes "${genes_csv}" --gtf "${GTF}" \
  --sample-sheet "${SAMPLE_SHEET}" --bam-pattern "${FILTERED}/{gene}/{sample}.bam" \
  --read-counts "${QC_DIR}/fastq_read_counts.tsv" --label "gene-filtered reads" \
  --outdir "${RESULTS_DIR}/fragmentation/filtered"

for gene in "${TARGET_GENES[@]}"; do
  out="${FIG}/${gene}"
  # 3. sashimi gene figure, coverage figures and isoform abundance figures
  echo "[figures] ${gene}: sashimi gene figure"
  python "${ANALYSIS_DIR}/scripts/gene_figure.py" --gene "${gene}" --gtf "${GTF}" \
    --sample-sheet "${SAMPLE_SHEET}" --bam-pattern "${FILTERED}/{gene}/{sample}.bam" --outdir "${out}"
  echo "[figures] ${gene}: coverage figures"
  python "${ANALYSIS_DIR}/scripts/coverage_plot.py" --gene "${gene}" --gtf "${GTF}" \
    --sample-sheet "${SAMPLE_SHEET}" --bam-pattern "${FILTERED}/{gene}/{sample}.bam" --outdir "${out}"
  wd="espresso_out/work_dir"
  while read -r method matrix iso_gtf; do
    if [[ ! -s "${matrix}" || ! -s "${iso_gtf}" ]]; then
      echo "[figures] ${gene}: ${method} isoforms not available, skipped"
      continue
    fi
    python "${ANALYSIS_DIR}/scripts/isoform_abundance_plot.py" --gene "${gene}" --method "${method}" \
      --matrix "${matrix}" --isoform-gtf "${iso_gtf}" --ref-gtf "${GTF}" --sample-sheet "${SAMPLE_SHEET}" \
      --read-counts "${QC_DIR}/fastq_read_counts.tsv" --outdir "${out}" < /dev/null \
      || echo "[figures] ${gene}: ${method} isoform figure failed (see the message above)"
  done <<EOF_METHODS
ESPRESSO_fast ${FAST_DIR}/espresso/${gene}/${wd}/samples_N2_R0_abundance.esp ${FAST_DIR}/espresso/${gene}/${wd}/samples_N2_R0_updated.gtf
AMALGAM_fast ${FAST_DIR}/amalgam/${gene}/amalgam_counts_matrix.tsv ${FAST_DIR}/amalgam/${gene}/amalgam.transcripts.gtf
AMALGAM_full ${PROJECT_DIR}/amalgam_full/amalgam_counts_matrix.tsv ${PROJECT_DIR}/amalgam_full/amalgam.transcripts.gtf
ESPRESSO_full ${PIPELINE_DIR}/${wd}/samples_N2_R0_abundance.esp ${PIPELINE_DIR}/${wd}/samples_N2_R0_updated.gtf
EOF_METHODS

  # 4. IGV-style pages: only the gene of interest is annotated
  read -r chrom start end < <(awk -F'\t' -v g="${gene}" \
    '$3 == "gene" && index($9, "gene_name \"" g "\";") {print $1, $4, $5; exit}' "${out}/${gene}_locus.gtf")
  region="${chrom}:$(( start > 500 ? start - 500 : 1 ))-$(( end + 500 ))"
  grep -F "gene_name \"${gene}\";" "${out}/${gene}_locus.gtf" > "${out}/igv/${gene}_only.gtf"
  for sample in $(sample_ids); do
    bam="${FILTERED}/${gene}/${sample}.bam"
    [[ -s "${bam}" ]] || continue
    echo "[figures] ${gene}: IGV view ${sample}"
    for ext in pdf png; do  # png for the slide deck (step 10)
      python "${ANALYSIS_DIR}/scripts/igv_plot.py" --bam "${bam}" --region "${region}" --fasta "${GENOME}" \
        --gtf "${out}/igv/${gene}_only.gtf" --outfile "${out}/igv/${gene}_${sample}_all_reads.${ext}" \
        --title "${gene}  ${sample}  ${region}  (gene reads)"
    done
  done
done

# 5. Figures 1-2 and the paired kit comparison
python "${ANALYSIS_DIR}/scripts/pub_figures.py" --sample-sheet "${SAMPLE_SHEET}" \
  --read-counts "${QC_DIR}/fastq_read_counts.tsv" --filter-stats "${STATS}" \
  --frag-dir "${RESULTS_DIR}/fragmentation/filtered" --genes "${genes_csv}" --outdir "${FIG}"
echo "[figures] done: ${FIG}"
