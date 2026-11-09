# TEQUILA-seq: EIF2AK3 (PERK) isoforms in MDMs and iPSC neurons (EAP samples)

This folder holds the analysis of the EIF2AK3 (PERK) TEQUILA-seq run
(`/mnt/isilon/lin_lab_share/Feng/Elena_EAP/TEQUILA_Kidney_EAP`, PromethION flowcells 3B and 3C,
started 2026-11-04 13:00). The probes target EIF2AK3. EIF2AK2 (PKR) is also analysed, as an
off-target comparison gene. The code is on branch `EIF2AK3` of this repository.

The analysis runs this repository's TEQUILA-seq code and adds two things: a fast path that works
only on the target gene's reads, and AMALGAM as a second isoform caller.

3B and 3C are the same library sequenced twice. Each barcode's reads from both flowcells are
aligned separately and pooled for quantification. QC is still reported per flowcell.

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
   |   02  FIRE v0.5.0 index -> EIF2AK3 / EIF2AK2 anchor 31-mers (unique to the gene)
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
EIF2AK3 31-mers that FIRE counts as anchors, using the same rule as FIRE's read decoder. The anchors
leave out sequence that EIF2AK3 shares with other genes, such as the EIF2AK1/2/4 kinase domain.
Neither EIF2AK3 nor EIF2AK2 is on `panel_653`, so the protein-coding (`pcg`) dictionary is used.
`scripts/kmer_bait_reads.py` then scans every read with numpy and keeps a read that has at least
`BAIT_MIN_HITS` (default 3) anchor positions on either strand. Only those reads are aligned, so
ESPRESSO and AMALGAM run on thousands of reads instead of millions. `fast/bait_align_stats.tsv` shows
how many baited reads land on the gene. A value near 1 confirms that the bait is specific.

## Running on the HPC, step by step

Paste each block into your terminal on the HPC login node, one block at a time, in this order.
Every block starts with `cd`, so it also works in a new terminal window.

### Step 1. Download the code (once)

```bash
mkdir -p /home/wangf3/feng/projects/TEQUILA_EIF2AK3
cd /home/wangf3/feng/projects/TEQUILA_EIF2AK3
git clone -b EIF2AK3 https://github.com/Max2W/TEQUILA-seq.git TEQUILA-seq
```

If git asks for a username and password, press Ctrl+C and use this block instead. It uses your
SSH key for GitHub:

```bash
cd /home/wangf3/feng/projects/TEQUILA_EIF2AK3
git clone -b EIF2AK3 git@github.com:Max2W/TEQUILA-seq.git TEQUILA-seq
```

### Step 2. Tell the analysis where FIRE is (once)

This searches your home folder and `/mnt/isilon/lin_lab_share/Feng` for your FIRE v0.5.0 checkout
and FIRE index, checks them, and saves the two paths:

```bash
cd /home/wangf3/feng/projects/TEQUILA_EIF2AK3/TEQUILA-seq/TEQUILA_EIF2AK3
bash set_fire_paths.sh
```

It finishes with `Saved to ...`. If it lists several checkouts or indexes instead, run it again
with the two paths you want, copied from its list:

```bash
cd /home/wangf3/feng/projects/TEQUILA_EIF2AK3/TEQUILA-seq/TEQUILA_EIF2AK3
bash set_fire_paths.sh /PATH/TO/FIRE/CHECKOUT /PATH/TO/FIRE/INDEX
```

### Step 3. Install software and download the genome (once, up to about an hour)

This creates three conda environments, downloads AMALGAM, GRCh38 and GENCODE v44, and checks FIRE.
It runs in the background, so a dropped connection does not stop it:

```bash
cd /home/wangf3/feng/projects/TEQUILA_EIF2AK3/TEQUILA-seq/TEQUILA_EIF2AK3
nohup bash 00_setup.sh > setup.log 2>&1 &
```

To watch it (press Ctrl+C to stop watching; setup keeps running):

```bash
tail -f /home/wangf3/feng/projects/TEQUILA_EIF2AK3/TEQUILA-seq/TEQUILA_EIF2AK3/setup.log
```

It is finished when this prints `=== SETUP COMPLETE ===`:

```bash
tail -3 /home/wangf3/feng/projects/TEQUILA_EIF2AK3/TEQUILA-seq/TEQUILA_EIF2AK3/setup.log
```

If the last lines show an error instead, send them to Claude.

### Step 4. Start the analysis

This submits all jobs to SLURM. Each job waits for the ones it needs.

```bash
cd /home/wangf3/feng/projects/TEQUILA_EIF2AK3/TEQUILA-seq/TEQUILA_EIF2AK3
bash run_all.sh
```

If `sbatch` complains about a missing account or partition, run this with your partition name
(here `defq`) and then run Step 4 again:

```bash
cd /home/wangf3/feng/projects/TEQUILA_EIF2AK3/TEQUILA-seq/TEQUILA_EIF2AK3
echo 'SBATCH_EXTRA="--partition=defq"' >> local_settings.sh
```

### Step 5. Check on the jobs

Your jobs and their state (PD = waiting, R = running; a job that is gone has finished):

```bash
squeue -u $USER
```

The five newest log files, and the end of the newest one:

```bash
cd /home/wangf3/feng/projects/TEQUILA_EIF2AK3/logs
ls -t | head -5
tail -20 "$(ls -t | head -1)"
```

A job that failed has an error at the end of its `.err` file. Send that file's last lines to Claude.

### Step 6. Look at the results

The fast path finishes in a few hours; the full path takes days. QC first:

```bash
cd /home/wangf3/feng/projects/TEQUILA_EIF2AK3
column -t fast/bait_align_stats.tsv
column -t qc/fastq_read_counts.tsv
```

The EIF2AK3 report from the fast path:

```bash
cd /home/wangf3/feng/projects/TEQUILA_EIF2AK3
cat results/EIF2AK3/ESPRESSO_fast/EIF2AK3_ESPRESSO_fast_report.md
cat results/EIF2AK3/AMALGAM_fast/EIF2AK3_AMALGAM_fast_report.md
```

All figures (PDF/PNG) and tables:

```bash
ls /home/wangf3/feng/projects/TEQUILA_EIF2AK3/results/EIF2AK3/*/
```

### Getting code updates later

```bash
cd /home/wangf3/feng/projects/TEQUILA_EIF2AK3/TEQUILA-seq
git pull
```

### Running single steps (optional)

Each numbered script can also be run alone with `bash <script>`; it submits itself to SLURM.
`bash 02_fire_gene_kmers.sh --check` only opens the FIRE index. `DRYRUN=1 bash 06_espresso_full.sh`
writes the Snakemake config and runs `snakemake -n`. `bash 08_isoform_report.sh` remakes the
reports from whatever has finished.

### Folder layout on the HPC

```
/home/wangf3/feng/projects/TEQUILA_EIF2AK3/
  TEQUILA-seq/            this repo, branch EIF2AK3 (code; conda_env/ = pipeline env)
    TEQUILA_EIF2AK3/      these scripts, plus local_settings.sh (your FIRE paths) and setup.log
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

1. **Capture.** `qc/per_sample_qc.tsv` and `qc/top_genes.tsv` should show EIF2AK3 as the top
   captured gene. EIF2AK2 reads show the background level.
2. **Flowcell agreement.** 3B and 3C are the same library, so the per-flowcell rows of
   `qc/per_input_qc.tsv` should show similar on-target rates for each barcode.
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
against a brute-force k-mer count, the report and AMALGAM merge on synthetic data, the
FIRE k-mer extraction against a mock of the index API, and the job chaining with a stand-in
`sbatch`. Nothing has yet been run on Respublica, against your real FIRE index, or on the reads.
