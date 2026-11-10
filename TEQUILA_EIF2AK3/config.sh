#!/bin/bash
# Central settings for the EIF2AK3 (PERK) TEQUILA-seq analysis (Elena / EAP samples).
# Every step script sources this file. Edit here, not in the step scripts.

# ---------------------------------------------------------------------------
# Locations on the HPC
# ---------------------------------------------------------------------------
# Project folder that holds the code clone and all outputs.
# Kept separate from /home/wangf3/feng/projects/fire-figures on purpose.
PROJECT_DIR="${PROJECT_DIR:-/home/wangf3/feng/projects/TEQUILA_EIF2AK3}"

# The TEQUILA-seq code (this git repo, branch EIF2AK3) cloned inside PROJECT_DIR
CODE_DIR="${CODE_DIR:-${PROJECT_DIR}/TEQUILA-seq}"
ANALYSIS_DIR="${_TEQ_DIR:-${CODE_DIR}/TEQUILA_EIF2AK3}"

# Where outputs go
PIPELINE_DIR="${PROJECT_DIR}/pipeline"      # Snakemake working dir (alignment, ESPRESSO)
FASTQ_DIR="${PROJECT_DIR}/fastq"            # per-sample, per-flowcell merged fastq.gz
REF_DIR="${PROJECT_DIR}/references"
QC_DIR="${PROJECT_DIR}/qc"
RESULTS_DIR="${PROJECT_DIR}/results"
LOG_DIR="${PROJECT_DIR}/logs"

# Conda environments (created by 00_setup.sh)
PIPELINE_ENV="${CODE_DIR}/conda_env"
# runs a command inside PIPELINE_ENV; written by lib.sh (ensure_conda_wrapper)
CONDA_WRAPPER="${PROJECT_DIR}/bin/pipeline_env"
ANALYSIS_ENV="${PROJECT_DIR}/envs/analysis"

# ---------------------------------------------------------------------------
# Raw data (dorado/MinKNOW fastq_pass, demultiplexed by barcode)
# ---------------------------------------------------------------------------
RAW_DIR="${RAW_DIR:-/mnt/isilon/lin_lab_share/Feng/Elena_EAP/TEQUILA_Kidney_EAP}"

# Flowcell runs to use. Each entry is "<short label>:<run folder>".
# 3B and 3C are the same library sequenced twice, so both are merged per barcode
# (aligned separately, pooled by ESPRESSO / AMALGAM; per-flowcell QC is kept).
FLOWCELLS=(
  "FC3B:20261104_1300_3B_PBM77320_5544ec5f"
  "FC3C:20261104_1300_3C_PBM77349_35b0a072"
)

SAMPLE_SHEET="${ANALYSIS_DIR}/samples.tsv"

# ---------------------------------------------------------------------------
# Reference: GRCh38 primary assembly + GENCODE v49 (the release FIRE is tested with).
# One annotation for everything: FIRE index, minimap2 junctions, ESPRESSO, AMALGAM.
# ---------------------------------------------------------------------------
GENCODE_RELEASE=49
FASTA_NAME="GRCh38.primary_assembly.genome.fa"
GTF_NAME="gencode.v${GENCODE_RELEASE}.annotation.gtf"
FASTA_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${GENCODE_RELEASE}/GRCh38.primary_assembly.genome.fa.gz"
GTF_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${GENCODE_RELEASE}/gencode.v${GENCODE_RELEASE}.annotation.gtf.gz"

# ---------------------------------------------------------------------------
# Target genes
# ---------------------------------------------------------------------------
# EIF2AK3 = PERK, the gene the TEQUILA-seq probes target. EIF2AK2 (PKR) is kept
# as an off-target comparison: its reads show the background capture level.
TARGET_GENES=(EIF2AK3 EIF2AK2)
# genes in TARGET_GENES that are internal controls: processed the same way, but
# reported under results/internal_control/<GENE>/ instead of results/<GENE>/
CONTROL_GENES=(EIF2AK2)

# ---------------------------------------------------------------------------
# Fast path: FIRE k-mer bait -> targeted alignment -> ESPRESSO + AMALGAM
# ---------------------------------------------------------------------------
# FIRE v0.5.0 is downloaded and compiled by 00_setup.sh into the project folder,
# and its index is built from the GENCODE v49 reference by 02_fire_gene_kmers.sh.
FIRE_REPO="${FIRE_REPO:-https://github.com/Max2W/fire.git}"
FIRE_TAG="v0.5.0"
FIRE_DIR="${PROJECT_DIR}/tools/fire"
# Which FIRE index 02_fire_gene_kmers.sh builds:
#   target (default): only the TARGET_GENES (seconds). A k-mer is then "unique"
#     only among these genes; genome-wide specificity comes from step 03, which
#     aligns every baited read to the whole genome and keeps only reads on the
#     gene locus, so reads from paralogs (EIF2AK1/2/4) are dropped there.
#   all: every protein-coding gene (~20 min), with FIRE's panel_653.txt. Also the
#     automatic fallback if the target-only build fails.
FIRE_INDEX_MODE="${FIRE_INDEX_MODE:-target}"
FIRE_INDEX_ALL="${REF_DIR}/fire_index_${FIRE_TAG}_gencode_v${GENCODE_RELEASE}_pcg"
FIRE_INDEX_TARGET="${REF_DIR}/fire_index_${FIRE_TAG}_gencode_v${GENCODE_RELEASE}_$(IFS=_; echo "${TARGET_GENES[*]}")"
if [[ "${FIRE_INDEX_MODE}" == "all" ]]; then FIRE_INDEX="${FIRE_INDEX_ALL}"; else FIRE_INDEX="${FIRE_INDEX_TARGET}"; fi
FIRE_PANEL="${FIRE_DIR}/config/panel_653.txt"   # used by the "all" build
# pixi (FIRE's installer): the one on your PATH, else a private copy that
# 00_setup.sh installs here
TEQ_PIXI_DIR="${PROJECT_DIR}/tools/pixi"
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

# AMALGAM (https://github.com/RNA-ROB/amalgam), cloned by 00_setup.sh
AMALGAM_REPO="https://github.com/RNA-ROB/amalgam.git"
AMALGAM_COMMIT="${AMALGAM_COMMIT:-main}"   # pin a commit hash for a reproducible run
AMALGAM_DIR="${PROJECT_DIR}/tools/amalgam"
AMALGAM_ENV="${PROJECT_DIR}/envs/amalgam"

# ---------------------------------------------------------------------------
# SLURM
# ---------------------------------------------------------------------------
# Extra sbatch flags for every job this analysis submits, e.g. "--account=xxx".
# Also passed to the Snakemake cluster jobs of 06_espresso_full.sh.
SBATCH_EXTRA="${SBATCH_EXTRA:-}"
export TEQUILA_SBATCH_EXTRA="${SBATCH_EXTRA}"
