#!/bin/bash
# Slide deck of the swab pilot, in the layout of the TEQUILA-seq progress decks. Run after 09_gene_figures.sh.
#   1. scripts/deck_figures.py (gene-filtered BAMs of step 09, isoforms of step 04/05):
#      reads per sample; per gene the canonical gene model, gene CPM + isoform fractions with isoform structures,
#      the same by ORF category (canonical / protein coding / PTC / incomplete), and sashimi + inclusion-level
#      figures for the best-supported splicing events -> ${RESULTS_DIR}/figures/deck/
#   2. scripts/make_slides.py: all figures of steps 09 and 10 in one .pptx, with titles and bullets computed from
#      the result tables -> ${RESULTS_DIR}/slides/TEQUILA_Epilepsy_Swab_pilot.pptx
#
# Light work (minutes): LOCAL=1 runs it here instead of submitting it.
# Usage: bash 10_slides.sh
set -euo pipefail

_D="${TEQ_ANALYSIS_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
source "${_D}/lib.sh"

preflight() {
  for gene in "${TARGET_GENES[@]}"; do
    compgen -G "${FAST_DIR}/bam_filtered/${gene}/*.bam" >/dev/null || {
      echo "no gene-filtered BAMs in ${FAST_DIR}/bam_filtered/${gene}: run 09_gene_figures.sh first" >&2; exit 1; }
  done
  [[ -s "${RESULTS_DIR}/figures/read_filter.tsv" ]] || {
    echo "missing ${RESULTS_DIR}/figures/read_filter.tsv: run 09_gene_figures.sh first" >&2; exit 1; }
}
if [[ -z "${AFTER:-}" ]]; then preflight; fi

if submit_self epi_slides 2 8G 01:00:00 "${_D}/10_slides.sh" "$@"; then
  exit 0
fi
preflight

activate_env "${ANALYSIS_ENV}"
if ! python -c "import pptx" 2>/dev/null; then
  echo "[slides] installing python-pptx into the analysis environment"
  python -m pip install --quiet "python-pptx>=0.6.21"
fi
GTF="${REF_DIR}/${GTF_NAME}"
GENOME="${REF_DIR}/${FASTA_NAME}"
FIG="${RESULTS_DIR}/figures"
genes_csv="$(IFS=,; echo "${TARGET_GENES[*]}")"

# isoform quantification per gene: the first that exists, in this order
wd="espresso_out/work_dir"
iso_args=()
for gene in "${TARGET_GENES[@]}"; do
  while read -r method matrix iso_gtf; do
    if [[ -s "${matrix}" && -s "${iso_gtf}" ]]; then
      iso_args+=(--isoforms "${gene}:${method}:${matrix}:${iso_gtf}")
      echo "[slides] ${gene}: isoforms from ${method}"
      break
    fi
  done <<EOF_METHODS
ESPRESSO_fast ${FAST_DIR}/espresso/${gene}/${wd}/samples_N2_R0_abundance.esp ${FAST_DIR}/espresso/${gene}/${wd}/samples_N2_R0_updated.gtf
AMALGAM_full ${PROJECT_DIR}/amalgam_full/amalgam_counts_matrix.tsv ${PROJECT_DIR}/amalgam_full/amalgam.transcripts.gtf
AMALGAM_fast ${FAST_DIR}/amalgam/${gene}/amalgam_counts_matrix.tsv ${FAST_DIR}/amalgam/${gene}/amalgam.transcripts.gtf
ESPRESSO_full ${PIPELINE_DIR}/${wd}/samples_N2_R0_abundance.esp ${PIPELINE_DIR}/${wd}/samples_N2_R0_updated.gtf
EOF_METHODS
done

python "${ANALYSIS_DIR}/scripts/deck_figures.py" --genes "${genes_csv}" --gtf "${GTF}" --fasta "${GENOME}" \
  --sample-sheet "${SAMPLE_SHEET}" --read-counts "${QC_DIR}/fastq_read_counts.tsv" \
  --bam-pattern "${FAST_DIR}/bam_filtered/{gene}/{sample}.bam" "${iso_args[@]}" --gene-reads "${FIG}/read_filter.tsv" \
  --outdir "${FIG}/deck"

# public data and literature figures (inputs committed under public_data/; optional rin.tsv of this pilot)
rin_args=(); rin_slide=()
RIN_TSV="${PROJECT_DIR}/rin.tsv"   # sample_id, RIN; sample data: kept on scratch, not in the repository
if [[ -s "${RIN_TSV}" ]]; then
  rin_args=(--our-rin "${RIN_TSV}"); rin_slide=(--rin "${RIN_TSV}")
fi
# GENE_SET=main: the full project deck. Any other gene set (GENE_SET=extra): a deck of that gene set's results only
MAIN=1
if [[ "${GENE_SET}" != "main" ]]; then MAIN=0; fi
deck_args=()
if (( MAIN )); then
# background figures hold published data only; this pilot's RIN is a result (rin_figures.py below)
python "${ANALYSIS_DIR}/scripts/public_data_figures.py" --public-dir "${ANALYSIS_DIR}/public_data" \
  --genes "${genes_csv}" --outdir "${FIG}/public"
# next specimen: nasal brushing vs buccal swab, from published numbers (public_data/nasal_literature.tsv)
ped=()
if [[ -s "${ANALYSIS_DIR}/public_data/nasal_pediatric.tsv" ]]; then ped=(--pediatric "${ANALYSIS_DIR}/public_data/nasal_pediatric.tsv"); fi
if [[ -s "${ANALYSIS_DIR}/public_data/nasal_gse152004_summary.tsv" ]]; then
  ped+=(--expression "${ANALYSIS_DIR}/public_data/nasal_gse152004_summary.tsv")
fi
python "${ANALYSIS_DIR}/scripts/nasal_literature_figures.py" --table "${ANALYSIS_DIR}/public_data/nasal_literature.tsv" \
  "${ped[@]}" --outdir "${FIG}/nasal"
# roadmap and accessible-tissue figures (public_data/gtex_panel_expression.tsv)
python "${ANALYSIS_DIR}/scripts/exploration_figures.py" --gtex "${ANALYSIS_DIR}/public_data/gtex_panel_expression.tsv" \
  --genes "${genes_csv}" --outdir "${FIG}/exploration"
deck_args+=(--public-figures "${FIG}/public" --nasal-figures "${FIG}/nasal" --exploration-figures "${FIG}/exploration")

# the whole probe panel in brain and accessible specimens (public_data/panel_expression_atlas.tsv), and the genes
# of extra_genes.txt: public expression, and their swab results when GENE_SET=extra has been run
mapfile -t EXTRA_GENES < <(sed 's/#.*//' "${ANALYSIS_DIR}/extra_genes.txt" 2>/dev/null | awk 'NF {print $1}')
extra_csv="$(IFS=,; echo "${EXTRA_GENES[*]}")"
python "${ANALYSIS_DIR}/scripts/panel_atlas_figures.py" --atlas "${ANALYSIS_DIR}/public_data/panel_expression_atlas.tsv" \
  --measured "${ANALYSIS_DIR}/public_data/specimen_panel_tpm.tsv" \
  --genes "${genes_csv}${extra_csv:+,${extra_csv}}" --frequency "${ANALYSIS_DIR}/public_data/dee_gene_frequency.tsv" \
  --outdir "${FIG}/atlas"
deck_args+=(--atlas-figures "${FIG}/atlas")
if [[ -n "${extra_csv}" ]]; then
  python "${ANALYSIS_DIR}/scripts/public_data_figures.py" --public-dir "${ANALYSIS_DIR}/public_data" \
    --genes "${extra_csv}" --outdir "${FIG}/public_extra"
  deck_args+=(--extra-genes "${extra_csv}" --extra-public "${FIG}/public_extra")
  if [[ -s "${PROJECT_DIR}/results_extra/figures/read_filter.tsv" ]]; then
    deck_args+=(--extra-results "${PROJECT_DIR}/results_extra")
  else
    echo "[slides] no swab results for ${extra_csv} yet (GENE_SET=extra run): public data only for these genes"
  fi
fi
else
  deck_args+=(--gene-deck)
fi

# full path (all reads aligned, 07_qc_full.sh): capture QC, top genes, and the same gene read filter on the
# genome-wide BAMs to see how many gene reads the k-mer baiting of the fast path recovered
full_args=()
if [[ -s "${QC_DIR}/per_sample_qc.tsv" ]]; then
  FULL_FILTER="${FIG}/full_read_filter.tsv"
  first=1
  for gene in "${TARGET_GENES[@]}"; do
    for sample in $(sample_ids); do
      bam="${QC_DIR}/bam/${sample}.bam"
      [[ -s "${bam}" ]] || continue
      hdr=(); if (( first )); then hdr=(--header); first=0; fi
      python "${ANALYSIS_DIR}/scripts/gene_read_filter.py" --gene "${gene}" --gtf "${FIG}/${gene}/${gene}_locus.gtf" \
        --bam "${bam}" --out-bam "${PROJECT_DIR}/full_filtered/${gene}/${sample}.bam" --sample "${sample}" "${hdr[@]}"
    done
  done > "${FULL_FILTER}.tmp"
  mv "${FULL_FILTER}.tmp" "${FULL_FILTER}"
  echo "[slides] full-path read filter:"; (column -t "${FULL_FILTER}" 2>/dev/null || cat "${FULL_FILTER}")
  python "${ANALYSIS_DIR}/scripts/full_path_figures.py" --genes "${genes_csv}" --sample-sheet "${SAMPLE_SHEET}" \
    --qc-dir "${QC_DIR}" --fast-filter "${FIG}/read_filter.tsv" --full-filter "${FULL_FILTER}" \
    --outdir "${FIG}/full_path"
  full_args=(--full-figures "${FIG}/full_path")
  if (( MAIN )) && [[ -s "${FIG}/expression_qc/expression_qc.tsv" ]]; then
    full_args+=(--expression-figures "${FIG}/expression_qc")
  else
    echo "[slides] no expression QC yet (11_expression_qc.sh): panel on-target, PCA and correlation slides left out"
  fi
else
  echo "[slides] no ${QC_DIR}/per_sample_qc.tsv (07_qc_full.sh): full-path slides left out"
fi

# RIN against read-level integrity; RNA composition (12_rna_composition.sh); sex-marker sample identity check
# (needs 11_expression_qc.sh)
extra_args=()
if [[ -s "${RIN_TSV}" && -s "${RESULTS_DIR}/fragmentation/filtered/fragmentation_summary.tsv" ]]; then
  python "${ANALYSIS_DIR}/scripts/rin_figures.py" --rin "${RIN_TSV}" \
    --fragmentation "${RESULTS_DIR}/fragmentation/filtered/fragmentation_summary.tsv" --sample-sheet "${SAMPLE_SHEET}" \
    --martorella "${ANALYSIS_DIR}/public_data/martorella2023_summary.tsv" --genes "${genes_csv}" --outdir "${FIG}/rin"
  extra_args+=(--rin-figures "${FIG}/rin")
fi
if (( MAIN )) && [[ -s "${QC_DIR}/composition/composition.tsv" ]]; then
  counts=()
  if [[ -s "${QC_DIR}/expression/gene_counts.tsv" ]]; then counts=(--gene-counts "${QC_DIR}/expression/gene_counts.tsv"); fi
  python "${ANALYSIS_DIR}/scripts/composition_figures.py" --composition-dir "${QC_DIR}/composition" "${counts[@]}" \
    --expected "${ANALYSIS_DIR}/public_data/genome_expectation.tsv" --sample-sheet "${SAMPLE_SHEET}" \
    --outdir "${FIG}/composition"
  extra_args+=(--composition-figures "${FIG}/composition")
else
  echo "[slides] no ${QC_DIR}/composition (12_rna_composition.sh): RNA composition slides left out"
fi
if (( MAIN )) && [[ -s "${QC_DIR}/expression/gene_counts.tsv" ]]; then
  python "${ANALYSIS_DIR}/scripts/sex_check.py" --counts "${QC_DIR}/expression/gene_counts.tsv" \
    --sample-sheet "${SAMPLE_SHEET}" --outdir "${FIG}/sex_check"
  extra_args+=(--sex-figures "${FIG}/sex_check")
fi

# coverage of known patient variants: the variant table holds patient information and is kept on scratch, never in
# the repository (${PROJECT_DIR}/variants/patient_variants.tsv)
var_args=()
VARIANTS="${PROJECT_DIR}/variants/patient_variants.tsv"
if (( MAIN )) && [[ -s "${VARIANTS}" ]]; then
  python "${ANALYSIS_DIR}/scripts/variant_coverage.py" --variants "${VARIANTS}" --genes "${genes_csv}" --gtf "${GTF}" \
    --fasta "${GENOME}" --sample-sheet "${SAMPLE_SHEET}" --bam-pattern "${FAST_DIR}/bam_filtered/{gene}/{sample}.bam" \
    --outdir "${FIG}/variants"
  var_args=(--variant-figures "${FIG}/variants")
else
  echo "[slides] no ${VARIANTS}: patient variant coverage slides left out"
fi

# run label for the sequencing slide: flow cell id(s) and run date(s) from the run folder names
labels=()
for fc in "${FLOWCELLS[@]}"; do
  run="${fc#*:}"
  labels+=("flow cell $(echo "${run}" | cut -d_ -f4) (run ${run%%_*})")
done
run_label="$(IFS=';'; echo "${labels[*]}" | sed 's/;/; /g')"

mkdir -p "${RESULTS_DIR}/slides"
deck="${RESULTS_DIR}/slides/TEQUILA_Epilepsy_Swab_pilot.pptx"
if (( ! MAIN )); then deck="${RESULTS_DIR}/slides/TEQUILA_Epilepsy_Swab_${GENE_SET}_genes.pptx"; fi
python "${ANALYSIS_DIR}/scripts/make_slides.py" --genes "${genes_csv}" --figures-dir "${FIG}" \
  --sample-sheet "${SAMPLE_SHEET}" --read-counts "${QC_DIR}/fastq_read_counts.tsv" --run-label "${run_label}" \
  "${deck_args[@]}" "${rin_slide[@]}" "${full_args[@]}" "${var_args[@]}" "${extra_args[@]}" --out "${deck}"
echo "[slides] done: ${deck}"
