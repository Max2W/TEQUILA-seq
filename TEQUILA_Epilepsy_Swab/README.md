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

## Where things are on the HPC

- **Code** (this repository) is in home: `/home/wangf3/feng/projects/TEQUILA_Epilepsy_Swab/TEQUILA-seq`.
- **Outputs, logs and working files** are on scratch: `/scr1/users/wangf3/TEQUILA_Epilepsy_Swab`.
- **Reference data** is on scratch: `/scr1/users/wangf3/TEQUILA_references`. It holds GRCh38 +
  GENCODE v49 and the minimap2 genome index.
- **Software** stays in home, in `/home/wangf3/feng/projects/TEQUILA_EIF2AK3`: the conda
  environments, FIRE v0.5.0, AMALGAM and pixi from the earlier EIF2AK3 setup.

Nothing is written to home except the code. Both paths together need roughly 100-150 GB on
scratch; `00_setup.sh` stops if less than 200 GB is free there. The full path deletes each
unsorted alignment as soon as its sorted copy exists.

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

### Step 1. Get the code (once)

```bash
mkdir -p /home/wangf3/feng/projects/TEQUILA_Epilepsy_Swab
cd /home/wangf3/feng/projects/TEQUILA_Epilepsy_Swab
git clone -b Epilepsy_Swab https://github.com/Max2W/TEQUILA-seq.git TEQUILA-seq
```

If it is already there, update it instead:

```bash
cd /home/wangf3/feng/projects/TEQUILA_Epilepsy_Swab/TEQUILA-seq
git pull
```

### Step 2. Set up (seconds)

It should end with `=== SETUP COMPLETE ===`:

```bash
cd /home/wangf3/feng/projects/TEQUILA_Epilepsy_Swab/TEQUILA-seq/TEQUILA_Epilepsy_Swab
bash 00_setup.sh
```

### Step 3. Start both paths

```bash
cd /home/wangf3/feng/projects/TEQUILA_Epilepsy_Swab/TEQUILA-seq/TEQUILA_Epilepsy_Swab
bash run_all.sh
```

### Step 4. Check on the jobs

```bash
squeue -u $USER -o "%.10i %.28j %.3t %.10M %R" | grep -E "epi_|snakejob"
```

```bash
sacct -u $USER -S today --format=JobID%10,JobName%28,State%14,ExitCode,Elapsed | grep epi_ | grep -v "\.ba\|\.ex"
```

If a job failed, the end of its error log names the failing line:

```bash
cd /scr1/users/wangf3/TEQUILA_Epilepsy_Swab/logs
for f in $(ls -t *.err | head -3); do echo "=== $f"; tail -15 "$f"; done
```

### Step 5. Results

```bash
cd /scr1/users/wangf3/TEQUILA_Epilepsy_Swab
column -t qc/fastq_read_counts.tsv
column -t fast/bait_align_stats.tsv
cat results/fragmentation/fast/fragmentation_report.md
cat results/STXBP1/ESPRESSO_fast/STXBP1_ESPRESSO_fast_report.md
cat results/SYNGAP1/ESPRESSO_fast/SYNGAP1_ESPRESSO_fast_report.md
```

### Step 6. Gene read filter and publication figures

`09_gene_figures.sh` runs these steps in order, using the fast-path gene BAMs:
1. **Gene read filter** (`scripts/gene_read_filter.py`). A read is kept for a gene only if it meets all of these:
   - at least 50% of its aligned bases lie on the gene's exons;
   - it has more bases on this gene's exons than on any other gene's;
   - if spliced, its transcript strand (minimap2 `ts`) is the gene's strand.

   The kept reads go to `fast/bam_filtered/`, and the numbers per sample to `results/figures/read_filter.tsv`.
2. **Detection and fragmentation** on the filtered reads, written to `results/fragmentation/filtered/`.
3. **Coverage and isoform figures** per gene, in `results/figures/<GENE>/`:
   - `scripts/coverage_plot.py` (filtered reads):
     - `<GENE>_coverage_lines`: one coverage line per sample, all reads (top) and spliced reads (bottom), on the
       canonical transcript with introns shortened; Zymo in blues, Genotek in reds, one shade per donor;
     - `<GENE>_coverage_tracks` and `_tracks_spliced`: IGV-style, one track per sample with its own `[0 - max]`
       scale, one colour per donor, genomic coordinates.
   - `scripts/isoform_abundance_plot.py`, for each isoform quantification that exists (ESPRESSO_fast,
     AMALGAM_fast, AMALGAM_full, ESPRESSO_full). CPM = isoform reads per million reads sequenced; the 4 most
     abundant isoforms are drawn and the rest summed as "Other":
     - `<GENE>_<METHOD>_isoform_abundance_by_donor`: one panel per donor, one bar per kit, isoform structures below;
     - `<GENE>_<METHOD>_isoform_abundance`: all samples in one panel.

   Every drawn number is also written to a table (`*_isoform_abundance.tsv`, `<GENE>_coverage.tsv.gz`).
4. **IGV-style pages** (`scripts/igv_plot.py`) per gene and sample, with filtered reads and only the gene of
   interest annotated, in `results/figures/<GENE>/igv/`.
5. **Figures 1 and 2** (`scripts/pub_figures.py`) and `kit_comparison.tsv`, which compares Genotek with Zymo,
   paired by donor.

The work is light, so it runs here on the login node in a few minutes:

```bash
cd /home/wangf3/feng/projects/TEQUILA_Epilepsy_Swab/TEQUILA-seq/TEQUILA_Epilepsy_Swab
LOCAL=1 bash 09_gene_figures.sh
```

### Step 7. Slide deck

`10_slides.sh` runs after step 6. It makes the figures of the TEQUILA-seq progress-deck style
(`scripts/deck_figures.py`), then puts them and the step 6 figures into one PowerPoint file
(`scripts/make_slides.py`). The new figures, in `results/figures/deck/`, are:

- **reads per sample**: sequencing reads per sample, sorted;
- **gene model**: the canonical (MANE) transcript of each gene;
- **isoform fraction**: gene abundance (CPM) above isoform fractions (%), with the isoform structures;
- **isoform category**: the same, grouped by open reading frame:
  - canonical;
  - protein coding;
  - premature termination codon (stop codon more than 50 nt before the last exon-exon junction);
  - incomplete ORF;
- **splicing events**: sashimi plots and inclusion (or splice-site usage) levels per sample. Events are cassette
  exons and alternative 3'/5' splice sites, found from the gene-filtered reads.

As in the progress decks:
- isoforms are shown individually when, in at least one sample, they reach a fraction of 5% or more of at least
  20 gene reads and have at least 3 reads;
- inclusion levels need at least 20 reads in the denominator.

Every number in a slide title or bullet is computed from the result tables when the deck is built. Isoforms come
from the first available of ESPRESSO_fast, AMALGAM_full, AMALGAM_fast or ESPRESSO_full. The first run installs
`python-pptx` into the analysis environment.

```bash
cd /home/wangf3/feng/projects/TEQUILA_Epilepsy_Swab/TEQUILA-seq/TEQUILA_Epilepsy_Swab
LOCAL=1 bash 10_slides.sh
```

The deck is `results/slides/TEQUILA_Epilepsy_Swab_pilot.pptx`.

The deck also has public-data slides from `scripts/public_data_figures.py`. Their inputs are committed in
`public_data/`, so no internet is needed:
- GTEx v10 and Human Protein Atlas expression of both genes;
- swab RNA quality from Martorella et al., BMC Genomics 2023: RIN, library QC failures and OMIM gene capture.

To add this pilot's RIN values to that figure and slide, put a `rin.tsv` (columns `sample_id`, `RIN`) in this folder
and run step 10 again.

### Folder layout on scratch

```
/scr1/users/wangf3/TEQUILA_Epilepsy_Swab/
  bin/pipeline_env      runs a command in the TEQUILA-seq conda env
  fire_index/           FIRE index of STXBP1 + SYNGAP1
  fastq/                merged fastq.gz per sample (needed by the full path)
  fast/                 kmers/, bait/, bam/<GENE>/, espresso/<GENE>/, amalgam/<GENE>/
  pipeline/             full Snakemake run (alignment/, espresso_out/)
  amalgam_full/         genome-wide AMALGAM
  qc/                   read counts, per-sample QC, top captured genes
  results/
    fragmentation/fast/ and full/   detection + fragmentation report and figures
    <GENE>/<METHOD>/                isoform reports (ESPRESSO / AMALGAM, fast / full)
    figures/                        figure1/2 (pdf, png, data), kit_comparison.tsv, read_filter.tsv
    figures/deck/                   progress-deck style figures and tables (step 10)
    slides/                         TEQUILA_Epilepsy_Swab_pilot.pptx (step 10)
    figures/<GENE>/                 coverage and isoform figures (pdf, png, data tables) and igv/ pages per sample
  logs/                 one .out/.err per SLURM job (job names start with epi_)
/scr1/users/wangf3/TEQUILA_references/   GRCh38, GENCODE v49, minimap2 index
```

## Notes

- **Full length** is judged against the canonical transcript. A read of another isoform that skips
  the canonical first or last exon is not counted as full length; the isoform reports cover those.
- The scripts were tested in a sandbox on synthetic reads (coverage, full-length counts and the
  kit contrast came out as planted) and with a stand-in `sbatch`. They have not run on the real data.
