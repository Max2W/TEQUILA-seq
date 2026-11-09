# TEQUILA-seq: epilepsy genes in buccal swab RNA (pilot)

Lab-member donors A-D each gave buccal swabs, collected with two swab / RNA-extraction kits, and the
RNA was captured with TEQUILA-seq probes for epilepsy genes, including STXBP1 and SYNGAP1. Swab RNA
is degraded (low RIN), but RIN does not say how much usable RNA remains for these genes. This
analysis asks:

1. **Detection.** Are the target genes detected in each sample, and at what depth?
2. **Fragmentation.** How much of each transcript do the reads cover? Are any reads full length?
3. **Isoforms.** Which transcripts are recovered, and do the two kits differ (paired by donor)?

Data: `/mnt/isilon/lin_lab_share/Feng/TEQUILA_Epilepsy_Buccal_Swab_Pilot`, flowcell 3A
(`20261104_1300_3A_PBM76834_dd62c4a3`), barcodes 17-24 (SQK-NBD114.24).

## Samples

| barcode | sample_id | donor | swab kit | RNA extraction | condition |
|---|---|---|---|---|---|
| NB17-20 | Zymo_A ... Zymo_D | A-D | Zymo R1107-1 | Zymo Quick-RNA Microprep | Zymo |
| NB21-24 | Genotek_A ... Genotek_D | A-D | Genotek OR-100 | Qiagen miRNeasy | Genotek |

## Genes

The probe panel is two genes, STXBP1 and SYNGAP1, listed in `target_genes.txt`. Every step reads
that file.

## Shared with the EIF2AK3 project

This is a separate project with its own data. It only shares software and reference files with the
EIF2AK3 project: the conda environments, FIRE v0.5.0, AMALGAM, pixi, GRCh38 + GENCODE v49 and the
minimap2 index already installed in `/home/wangf3/feng/projects/TEQUILA_EIF2AK3`. Nothing
is installed twice, and this project never writes into that folder. The only exception is the
shared minimap2 index: it is built once there, and a lock stops two jobs from building it at the
same time. Its own outputs go to `/home/wangf3/feng/projects/TEQUILA_Epilepsy_Swab`.

## The steps

The steps are the same as for EIF2AK3:

- **01** merges the barcode files.
- **Fast path:**
  - 02 builds a FIRE index of the target genes and extracts their gene-unique 31-mers;
  - 03 baits the reads that carry them and aligns those reads to the whole genome, keeping reads
    on each gene;
  - 04 / 05 run ESPRESSO and AMALGAM on those reads.
- **Full path:** 06 runs the standard TEQUILA-seq ESPRESSO run on all reads, 07 does QC, and 05 runs
  AMALGAM genome-wide.
- **08** writes the reports.

The new piece for this project is `scripts/fragment_report.py`, which 08 runs. For each gene and
sample it reports:

- reads on the gene, and reads per million of all reads sequenced;
- median read length, and **transcript coverage**: the fraction of the canonical (MANE)
  transcript between a read's 5'-most and 3'-most aligned positions;
- reads that reach the last exon (3') and the first exon (5'), and **full-length reads**, which
  reach both;
- a **gene-body coverage curve** (5' -> 3'). RNA that is fragmented and captured from the 3' end
  shows coverage falling toward the 5' end.

The isoform reports compare Genotek with Zymo, paired by donor (n = 4).

## Running on the HPC, step by step

Paste each block on the HPC login node, one at a time, in this order.

### Step 1. Download the code (once)

```bash
mkdir -p /home/wangf3/feng/projects/TEQUILA_Epilepsy_Swab
cd /home/wangf3/feng/projects/TEQUILA_Epilepsy_Swab
git clone -b Epilepsy_Swab https://github.com/Max2W/TEQUILA-seq.git TEQUILA-seq
```

### Step 2. Set up (seconds)

This checks the shared tools and genes and creates this project's folders. It should end with
`=== SETUP COMPLETE ===`:

```bash
cd /home/wangf3/feng/projects/TEQUILA_Epilepsy_Swab/TEQUILA-seq/TEQUILA_Epilepsy_Swab
bash 00_setup.sh
```

### Step 3. Start the analysis

The genome index is shared with the EIF2AK3 project. This block waits in the background until that
index is ready, then submits all jobs; its messages go to `start.log`:

```bash
cd /home/wangf3/feng/projects/TEQUILA_Epilepsy_Swab/TEQUILA-seq/TEQUILA_Epilepsy_Swab
nohup bash wait_and_run_all.sh > start.log 2>&1 &
```

To see whether the jobs were submitted:

```bash
cat /home/wangf3/feng/projects/TEQUILA_Epilepsy_Swab/TEQUILA-seq/TEQUILA_Epilepsy_Swab/start.log
```

### Step 4. Check on the jobs

```bash
squeue -u $USER -o "%.10i %.28j %.3t %.10M %R" | grep epi_
```

```bash
sacct -u $USER -S today --format=JobID%10,JobName%28,State%14,Elapsed,MaxRSS | grep epi_ | grep -v "\.ba\|\.ex"
```

### Step 5. Look at the results

Detection and fragmentation (fast path; ready a few hours after the start):

```bash
cd /home/wangf3/feng/projects/TEQUILA_Epilepsy_Swab
column -t qc/fastq_read_counts.tsv
column -t fast/bait_align_stats.tsv
cat results/fragmentation/fast/fragmentation_report.md
```

Isoform reports, one per gene:

```bash
cd /home/wangf3/feng/projects/TEQUILA_Epilepsy_Swab
cat results/STXBP1/ESPRESSO_fast/STXBP1_ESPRESSO_fast_report.md
cat results/SYNGAP1/ESPRESSO_fast/SYNGAP1_ESPRESSO_fast_report.md
```

All figures and tables:

```bash
ls /home/wangf3/feng/projects/TEQUILA_Epilepsy_Swab/results/*/ /home/wangf3/feng/projects/TEQUILA_Epilepsy_Swab/results/fragmentation/*/
```

### Getting code updates later

```bash
cd /home/wangf3/feng/projects/TEQUILA_Epilepsy_Swab/TEQUILA-seq
git pull
```

## Folder layout on the HPC

```
/home/wangf3/feng/projects/TEQUILA_Epilepsy_Swab/
  TEQUILA-seq/          this repo, branch Epilepsy_Swab (code)
    TEQUILA_Epilepsy_Swab/   these scripts, samples.tsv, target_genes.txt
  bin/pipeline_env      runs a command in the shared TEQUILA-seq conda env
  fire_index/           FIRE index of the target genes
  fastq/                merged fastq.gz per sample
  fast/                 kmers/, bait/, bam/<GENE>/, espresso/<GENE>/, amalgam/<GENE>/
  pipeline/             full Snakemake run
  amalgam_full/         genome-wide AMALGAM
  qc/                   read counts, per-sample QC, top captured genes
  results/
    fragmentation/fast/ and full/   detection + fragmentation report and figures
    <GENE>/<METHOD>/                isoform reports (ESPRESSO / AMALGAM, fast / full)
  logs/                 one .out/.err per SLURM job (job names start with epi_)
```

## Notes

- **Full length** is judged against the canonical transcript. A read of another isoform that skips
  the canonical first or last exon is not counted as full length; the isoform reports cover those.
- The scripts were tested in a sandbox on synthetic reads (coverage, full-length counts and the
  kit contrast came out as planted) and with a stand-in `sbatch`. They have not run on the real data.
