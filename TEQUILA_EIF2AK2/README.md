# TEQUILA-seq: EIF2AK2 isoforms in MDMs and iPSC neurons (EAP samples)

This folder holds the analysis of the EIF2AK2 TEQUILA-seq run
(`/mnt/isilon/lin_lab_share/Feng/Elena_EAP/TEQUILA_Kidney_EAP`, PromethION flowcells 3B and 3C,
started 2026-11-04 13:00). It runs this repository's TEQUILA-seq code and adds two things: a fast
path that works only on EIF2AK2 reads, and AMALGAM as a second isoform caller.

## Samples

`samples.tsv` holds the sample sheet. The `condition` column is what the statistics use.

| barcode | sample_id | donor | cell type | treatment | condition |
|---|---|---|---|---|---|
| NB01 | MDM_114_DMSO | EAP 1.5.114 | MDM | DMSO | Vehicle |
| NB02 | MDM_114_TG2h | EAP 1.5.114 | MDM | TG 2h | TG_2h |
| NB03 | MDM_114_LPS2h | EAP 1.5.114 | MDM | LPS 2h | LPS_2h |
| NB04-07 | MDM_130_{Ctl,Vehicle,TG2h,LPS2h} | EAP 1.5.130 | MDM | Ctl / Vehicle / TG 2h / LPS 2h | same |
| NB08-11 | MDM_136_{Ctl,Vehicle,TG2h,LPS2h} | EAP 1.5.136 | MDM | Ctl / Vehicle / TG 2h / LPS 2h | same |
| NB12-14 | iN_147_{Ctl,Vehicle,TG2h} | EAP 1.5.147 | iPSC neuron | Ctl / Vehicle / TG 2h | same |

The sheet treats DMSO in donor 114 as the vehicle. Donor 114 has no untreated control.

## Two paths

```
fastq_pass/barcodeNN (3B + 3C) --01--> fastq/<sample>/<flowcell>.fastq.gz
   |
   |-- FAST PATH (hours) -------------------------------------------------------------
   |   02  FIRE v0.5.0 index -> EIF2AK2 / EIF2AK3 anchor 31-mers (unique to the gene)
   |   03  k-mer bait every read -> align only the baited reads (minimap2, TEQUILA settings)
   |   04  ESPRESSO (this repo's Snakefile) on the gene BAMs with the gene's annotation
   |   05  AMALGAM on the gene BAMs                                (05_amalgam.sh fast)
   |
   |-- FULL PATH (days; the standard TEQUILA-seq run) ---------------------------------
   |   06  Snakemake: minimap2 + ESPRESSO S/C/Q on all reads (SLURM jobs)
   |   07  alignment QC: mapping, on-target, top captured genes, per flowcell
   |   05  AMALGAM genome-wide                                     (05_amalgam.sh full)
   |
   08  isoform reports for every finished run + cross-method comparison -> results/
```

**Why the fast path works.** FIRE's index stores every transcript of every gene, and for each 31-mer
it records whether the 31-mer occurs in only one gene. `scripts/fire_gene_kmers.py` collects the
EIF2AK2 31-mers that FIRE counts as anchors, using the same rule as FIRE's read decoder. The anchors
leave out sequence that EIF2AK2 shares with other genes, such as the EIF2AK1/3/4 kinase domain.
EIF2AK2 is not on `panel_653`, so the protein-coding (`pcg`) dictionary is used.
`scripts/kmer_bait_reads.py` then scans every read with numpy and keeps a read that has at least
`BAIT_MIN_HITS` (default 3) anchor positions on either strand. Only those reads are aligned, so
ESPRESSO and AMALGAM run on thousands of reads instead of millions. `fast/bait_align_stats.tsv` shows
how many baited reads land on the gene. A value near 1 confirms that the bait is specific.

## Running on the HPC

The project lives in its own folder, separate from `fire-figures`:

```bash
mkdir -p /home/wangf3/feng/projects/TEQUILA_EIF2AK2
```

```bash
cd /home/wangf3/feng/projects/TEQUILA_EIF2AK2 && git clone -b EIF2AK2 https://github.com/Max2W/TEQUILA-seq.git TEQUILA-seq
```

Edit `TEQUILA-seq/TEQUILA_EIF2AK2/config.sh`. Set `FIRE_DIR` (your FIRE v0.5.0 checkout, after
`pixi install && pixi run build-ext`) and `FIRE_INDEX` (an index made by `fire build` with v0.5.0).
If your cluster needs an account or partition, also set `SBATCH_EXTRA`, e.g. `--partition=defq`.
Then set up once, on a login node. This step creates three conda envs, clones AMALGAM, downloads
GRCh38 and GENCODE v44, and checks FIRE:

```bash
cd /home/wangf3/feng/projects/TEQUILA_EIF2AK2/TEQUILA-seq/TEQUILA_EIF2AK2 && bash 00_setup.sh
```

Submit everything as dependent SLURM jobs:

```bash
bash run_all.sh
```

You can also run one path at a time: `bash run_fast_path.sh`, or `SKIP_FASTQ=1 bash run_full_path.sh`.
Each numbered script can be run alone as well. It submits itself to SLURM; with `LOCAL=1` it runs in
the current shell instead, e.g. inside an interactive job. `bash 02_fire_gene_kmers.sh --check` only
opens the FIRE index and resolves the genes. `DRYRUN=1 bash 06_espresso_full.sh` writes the Snakemake
config and runs `snakemake -n`.

Folder layout on the HPC:

```
/home/wangf3/feng/projects/TEQUILA_EIF2AK2/
  TEQUILA-seq/            this repo, branch EIF2AK2 (code; conda_env/ = pipeline env)
  envs/                   analysis and AMALGAM conda envs
  tools/amalgam/          AMALGAM clone (commit printed in each AMALGAM log)
  references/             GRCh38 + GENCODE v44 (+ minimap2 .mmi, junction BED)
  fastq/                  merged fastq.gz per sample and flowcell
  fast/                   kmers/, bait/, bam/, espresso/<GENE>/, amalgam/<GENE>/
  pipeline/               full Snakemake run (alignment/, espresso_out/)
  amalgam_full/           genome-wide AMALGAM
  qc/                     fastq_read_counts.tsv, per_input_qc.tsv, top_genes.tsv, figures
  results/<GENE>/<METHOD>/   reports (METHOD = ESPRESSO_fast, AMALGAM_fast, ESPRESSO_full, AMALGAM_full)
  logs/                   one .out/.err per SLURM job
```

## What the report contains

`results/<GENE>/<METHOD>/<GENE>_<METHOD>_report.md` summarises the run. Alongside it:

- `_isoform_annotation.tsv`: one row per isoform. It says whether the isoform is annotated or novel,
  and gives its transcript type, intron chain, junctions absent from GENCODE, and closest annotated
  transcript. It also lists the alternative splicing events relative to the canonical transcript,
  from this repo's `scripts/Find_AS_Events.py`, and the mean proportion per cell type and condition.
- `_isoform_counts.tsv`, `_isoform_proportions.tsv`, `_gene_reads.tsv`.
- `_contrasts.tsv`: each condition against Vehicle within a cell type, paired by donor. It gives the
  change in isoform proportion per donor and the mean change. It also gives a paired t-test across
  donors (MDM TG and LPS: n = 3; MDM Ctl: n = 2, so no t-test) and a Cochran-Mantel-Haenszel test
  on read counts with donor as the stratum. Both tests come with Benjamini-Hochberg FDR. The neurons
  have one donor, so their contrasts are counts-only. Use them as preliminary data.
- Figures: `_structure_heatmap` (exon structures beside the proportion of each isoform per sample),
  `_stacked_proportions`, `_isoform_by_donor` (one line per donor across conditions), and
  `_gene_reads`.
- `results/<GENE>/method_comparison*`: isoforms matched by intron chain across ESPRESSO and AMALGAM
  and across the fast and full paths, with correlations of the per-sample proportions.

Samples with fewer than 20 gene reads are shown, but they are left out of the statistics.

## Things to check first

1. **EIF2AK2 or EIF2AK3.** The request says the probes are for EIF2AK3 (PERK), but the branch and
   folder are named EIF2AK2 (PKR). Both genes are analysed. `qc/per_sample_qc.tsv`,
   `qc/top_genes.tsv` and `fast/bait_summary_all.tsv` show which gene the panel actually captured.
2. **Flowcells 3B and 3C.** Both have barcode01-14 and are merged per barcode by default. The run
   folder is named `TEQUILA_Kidney_EAP`, and both flowcells also carry barcodes 15-24. If 3B and 3C
   are different libraries, `qc/fastq_read_counts.tsv` and the per-flowcell rows of
   `qc/per_input_qc.tsv` will show it. In that case, remove the wrong line from `FLOWCELLS` in
   `config.sh` and re-run.
3. **AMALGAM's StringTie step** uses default settings, as AMALGAM's README specifies. AMALGAM is under
   active development, so set `AMALGAM_COMMIT` in `config.sh` to pin the commit you report.

## Software

- TEQUILA-seq pipeline env: `conda_requirements.txt` at the repo root (snakemake 5.31, minimap2 2.17,
  ESPRESSO alpha1.2.2 from `ESPRESSO_alpha1.2.2/`).
- Analysis env: `envs/analysis_env.txt` (python 3.11, pandas, scipy, matplotlib, pysam, networkx,
  samtools 1.20, minimap2 2.28).
- AMALGAM: its own `conda_requirements.txt` (StringTie 3.0.0, GffCompare 0.12.6, htslib).
- FIRE v0.5.0: your pixi environment. It is used only to read the index in step 02.

The scripts were tested in a sandbox without the data or the cluster: the read bait was checked
against a brute-force k-mer count, the report and AMALGAM merge on synthetic EIF2AK2 data, the
FIRE k-mer extraction against a mock of the index API, and the job chaining with a stand-in
`sbatch`. Nothing has yet been run on Respublica, against your real FIRE index, or on the reads.
