#!/bin/bash
# Central settings for the epilepsy buccal-swab TEQUILA-seq pilot.
# Every step script sources this file. Edit here, not in the step scripts.

# ---------------------------------------------------------------------------
# Locations on the HPC
# ---------------------------------------------------------------------------
# Code lives in home; data and outputs live on scratch (home is not for data).
CODE_DIR="${CODE_DIR:-/home/wangf3/feng/projects/TEQUILA_Epilepsy_Swab/TEQUILA-seq}"
# all outputs, logs and working files of this project
PROJECT_DIR="${PROJECT_DIR:-/scr1/users/wangf3/TEQUILA_Epilepsy_Swab}"
ANALYSIS_DIR="${_TEQ_DIR:-${CODE_DIR}/TEQUILA_Epilepsy_Swab}"

# Where outputs go
PIPELINE_DIR="${PROJECT_DIR}/pipeline"      # Snakemake working dir (alignment, ESPRESSO)
FASTQ_DIR="${PROJECT_DIR}/fastq"            # per-sample merged fastq.gz
QC_DIR="${PROJECT_DIR}/qc"
RESULTS_DIR="${PROJECT_DIR}/results"
LOG_DIR="${PROJECT_DIR}/logs"

# ---------------------------------------------------------------------------
# Software installed by the earlier EIF2AK3 setup (kept in home, used read-only):
# conda environments, FIRE, AMALGAM, pixi
# ---------------------------------------------------------------------------
SHARED_DIR="${SHARED_DIR:-/home/wangf3/feng/projects/TEQUILA_EIF2AK3}"
PIPELINE_ENV="${SHARED_DIR}/TEQUILA-seq/conda_env"
ANALYSIS_ENV="${SHARED_DIR}/envs/analysis"
AMALGAM_ENV="${SHARED_DIR}/envs/amalgam"
AMALGAM_DIR="${SHARED_DIR}/tools/amalgam"
FIRE_DIR="${SHARED_DIR}/tools/fire"
TEQ_PIXI_DIR="${SHARED_DIR}/tools/pixi"
# reference data on scratch: GRCh38 + GENCODE v49, minimap2 genome index
REF_DIR="${REF_DIR:-/scr1/users/wangf3/TEQUILA_references}"
# runs a command inside PIPELINE_ENV (written by 00_setup.sh; replaces the repo's
# conda_wrapper, which only finds an env inside its own clone)
CONDA_WRAPPER="${PROJECT_DIR}/bin/pipeline_env"

# ---------------------------------------------------------------------------
# Raw data (dorado/MinKNOW fastq_pass, demultiplexed by barcode)
# ---------------------------------------------------------------------------
RAW_DIR="${RAW_DIR:-/mnt/isilon/lin_lab_share/Feng/TEQUILA_Epilepsy_Buccal_Swab_Pilot}"

# The swab pilot was sequenced on one flowcell (3A); its samples are barcodes 17-24.
# Each entry is "<short label>:<run folder inside RAW_DIR>".
FLOWCELLS=(
  "FC3A:20261104_1300_3A_PBM76834_dd62c4a3"
)

SAMPLE_SHEET="${ANALYSIS_DIR}/samples.tsv"

# ---------------------------------------------------------------------------
# Reference: GRCh38 primary assembly + GENCODE v49 (shared, see above)
# ---------------------------------------------------------------------------
GENCODE_RELEASE=49
FASTA_NAME="GRCh38.primary_assembly.genome.fa"
GTF_NAME="gencode.v${GENCODE_RELEASE}.annotation.gtf"
FASTA_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${GENCODE_RELEASE}/GRCh38.primary_assembly.genome.fa.gz"
GTF_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${GENCODE_RELEASE}/gencode.v${GENCODE_RELEASE}.annotation.gtf.gz"

# ---------------------------------------------------------------------------
# Target genes: one gene name per line in target_genes.txt (# starts a comment)
# ---------------------------------------------------------------------------
mapfile -t TARGET_GENES < <(sed 's/#.*//' "${ANALYSIS_DIR}/target_genes.txt" | awk 'NF {print $1}')
if (( ${#TARGET_GENES[@]} == 0 )); then
  echo "config.sh: no genes in ${ANALYSIS_DIR}/target_genes.txt" >&2
  exit 1
fi
# genes reported under results/internal_control/<GENE>/ instead of results/<GENE>/
CONTROL_GENES=()

# Isoform report: the swab kits are compared, paired by donor, against this kit
REFERENCE_CONDITION="Zymo"

# ---------------------------------------------------------------------------
# Fast path: FIRE k-mer bait -> targeted alignment -> ESPRESSO + AMALGAM
# ---------------------------------------------------------------------------
FIRE_REPO="${FIRE_REPO:-https://github.com/Max2W/fire.git}"
FIRE_TAG="v0.5.0"
# Which FIRE index 02_fire_gene_kmers.sh builds:
#   target (default): only the TARGET_GENES (seconds). Genome-wide specificity
#     comes from step 03, which aligns every baited read to the whole genome and
#     keeps only reads on the gene locus.
#   all: every protein-coding gene (~20 min), kept in REF_DIR.
#     Also the automatic fallback if the target-only build fails.
FIRE_INDEX_MODE="${FIRE_INDEX_MODE:-target}"
FIRE_INDEX_ALL="${REF_DIR}/fire_index_${FIRE_TAG}_gencode_v${GENCODE_RELEASE}_pcg"
FIRE_INDEX_TARGET="${PROJECT_DIR}/fire_index/fire_index_${FIRE_TAG}_gencode_v${GENCODE_RELEASE}_epilepsy_$(printf '%s\n' "${TARGET_GENES[@]}" | md5sum | cut -c1-8)"
if [[ "${FIRE_INDEX_MODE}" == "all" ]]; then FIRE_INDEX="${FIRE_INDEX_ALL}"; else FIRE_INDEX="${FIRE_INDEX_TARGET}"; fi
FIRE_PANEL="${FIRE_DIR}/config/panel_653.txt"   # used by the "all" build
if ! command -v pixi >/dev/null 2>&1 && [[ -x "${TEQ_PIXI_DIR}/bin/pixi" ]]; then
  export PATH="${TEQ_PIXI_DIR}/bin:${PATH}"
fi
# Optional site settings, e.g. SBATCH_EXTRA="--partition=defq" (git ignores this file)
if [[ -f "${ANALYSIS_DIR}/local_settings.sh" ]]; then
  source "${ANALYSIS_DIR}/local_settings.sh"
fi
# A read is baited when it carries at least this many gene-unique 31-mers
BAIT_MIN_HITS="${BAIT_MIN_HITS:-3}"
FAST_DIR="${PROJECT_DIR}/fast"

AMALGAM_REPO="https://github.com/RNA-ROB/amalgam.git"
AMALGAM_COMMIT="${AMALGAM_COMMIT:-main}"

# ---------------------------------------------------------------------------
# SLURM
# ---------------------------------------------------------------------------
# Extra sbatch flags for every job this analysis submits, e.g. "--account=xxx".
# Also passed to the Snakemake cluster jobs of 06_espresso_full.sh.
SBATCH_EXTRA="${SBATCH_EXTRA:-}"
export TEQUILA_SBATCH_EXTRA="${SBATCH_EXTRA}"
