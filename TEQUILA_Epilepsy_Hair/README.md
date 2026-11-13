# TEQUILA-seq: epilepsy genes in plucked hair-follicle RNA (pilot)

The lab-member donors A-D of the buccal-swab pilot (branch `Epilepsy_Swab`) give plucked scalp hair follicles. Each
donor's follicles are split between two storage conditions: frozen immediately (-80 C), or kept in RNAlater at room
temperature for 2-3 days, which is what a mailed clinical sample would go through. A third tube per donor is kept for
DNA. RNA is captured with the same TEQUILA-seq probe panel and analysed exactly like the swabs. Questions:

1. **Integrity.** Is hair-follicle RNA intact (RIN, full-length reads), and does RNAlater storage at room
   temperature preserve it?
2. **Detection and coverage.** Are the target genes detected, and how much of each transcript, and how many known
   patient variants, do the reads cover?
3. **Hair vs swab.** How does each donor's hair compare with the same donor's swabs?

Data: run folder and barcodes are set in `config.sh` and `samples.tsv` once the run exists (they read `SET_ME`
until then, and step 01 stops with a message).

## Samples

| sample_id | donor | collection and storage | condition |
|---|---|---|---|
| Frozen_A ... Frozen_D | A-D | plucked scalp hair follicles, frozen immediately (-80 C) | Frozen |
| RNAlater_A ... RNAlater_D | A-D | plucked scalp hair follicles, RNAlater at room temperature 2-3 days | RNAlater |

## Collection and extraction (per donor)

- **Collection.**
  - Pluck 25-30 scalp hairs from the back of the head (occipital) with tweezers, firmly at the base.
  - Keep only hairs with a visible whitish root sheath (growth-phase, anagen, follicles). Cut off the shaft about
    5 mm above the root.
  - Make three tubes of about 10 follicles each: 1 Frozen (RNA), 2 RNAlater (RNA), 3 DNA.
  - Record per tube: number of follicles with a sheath, hair type (fine or coarse), and minutes from plucking to
    freezing or RNAlater.
- **Storage.**
  - Tube 1 goes straight onto dry ice, then -80 C.
  - Tube 2 goes into RNAlater, stays at room temperature for 2-3 days, then -80 C or straight to extraction.
- **RNA.**
  - Homogenise with beads, then use Zymo Quick-RNA Microprep with the on-column DNase I step; TRI Reagent with Zymo
    Direct-zol MicroPrep is the published alternative.
  - Measure RIN on a high-sensitivity chip and record the yield.
- **DNA.**
  - Tube 3 with a DNA kit for small inputs (for example Qiagen DNeasy or QIAamp DNA Micro). Pool the follicles:
    single hairs give little DNA.
- **What to expect** (Martorella et al., BMC Genomics 2023): 10 follicles gave a median of 1,444 ng RNA with a mean
  RIN of 8.9. Coarse hair keeps its follicles better than fine hair.

Put the RIN values in `/scr1/users/wangf3/TEQUILA_Epilepsy_Hair/rin.tsv` (columns `sample_id`, `RIN`; sample
data stay on scratch). Step 10 compares each donor's hair with the same donor's swabs (`scripts/compare_specimens.py`,
read-only on the swab project).

## Genes

The probe panel is the 634-gene Epilepsy & NDD panel (`panel_genes.txt`, from Epilepsy_NDD634.xlsx). The genes of
interest are STXBP1 and SYNGAP1, listed in `target_genes.txt`; every step reads that file. The panel list is
used for the on-target rate (step 11).

## Where things are on the HPC

- **Code** (this repository) is in home: `/home/wangf3/feng/projects/TEQUILA_Epilepsy_Hair/TEQUILA-seq`.
- **Outputs, logs and working files** are on scratch: `/scr1/users/wangf3/TEQUILA_Epilepsy_Hair`.
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
mkdir -p /home/wangf3/feng/projects/TEQUILA_Epilepsy_Hair
cd /home/wangf3/feng/projects/TEQUILA_Epilepsy_Hair
git clone -b Epilepsy_Hair https://github.com/Max2W/TEQUILA-seq.git TEQUILA-seq
```

If it is already there, update it instead:

```bash
cd /home/wangf3/feng/projects/TEQUILA_Epilepsy_Hair/TEQUILA-seq
git pull
```

### Step 2. Set up (seconds)

It should end with `=== SETUP COMPLETE ===`:

```bash
cd /home/wangf3/feng/projects/TEQUILA_Epilepsy_Hair/TEQUILA-seq/TEQUILA_Epilepsy_Hair
bash 00_setup.sh
```

### Step 3. Start both paths

```bash
cd /home/wangf3/feng/projects/TEQUILA_Epilepsy_Hair/TEQUILA-seq/TEQUILA_Epilepsy_Hair
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
cd /scr1/users/wangf3/TEQUILA_Epilepsy_Hair/logs
for f in $(ls -t *.err | head -3); do echo "=== $f"; tail -15 "$f"; done
```

### Step 5. Results

```bash
cd /scr1/users/wangf3/TEQUILA_Epilepsy_Hair
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
cd /home/wangf3/feng/projects/TEQUILA_Epilepsy_Hair/TEQUILA-seq/TEQUILA_Epilepsy_Hair
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
cd /home/wangf3/feng/projects/TEQUILA_Epilepsy_Hair/TEQUILA-seq/TEQUILA_Epilepsy_Hair
LOCAL=1 bash 10_slides.sh
```

The deck is `results/slides/TEQUILA_Epilepsy_Hair_pilot.pptx`.

The deck also has public-data slides from `scripts/public_data_figures.py`. Their inputs are committed in
`public_data/`, so no internet is needed:
- GTEx v10 and Human Protein Atlas expression of both genes;
- noninvasive-specimen RNA quality from Martorella et al., BMC Genomics 2023: RIN, library QC failures and OMIM gene capture.

When the full path has run (`07_qc_full.sh`), step 10 also applies the gene read filter to the genome-wide BAMs
(`results/figures/full_read_filter.tsv`). It adds capture QC slides from `scripts/full_path_figures.py`:
- mapping rate, on-target rate and read length;
- the genes collecting the most reads;
- how many gene reads the k-mer baiting of the fast path recovered.

With `/scr1/users/wangf3/TEQUILA_Epilepsy_Hair/rin.tsv` (columns `sample_id`, `RIN`; sample data, kept on scratch),
step 10 adds:
- the pilot's RIN to the public-data RIN figure and slide;
- RIN against the read-level integrity of the target genes (`scripts/rin_figures.py`).

When step 11 has run, step 10 also adds a sex-marker sample identity check (`scripts/sex_check.py`; XIST against
Y-linked genes; the two samples of a donor must agree). Its output stays with the results on scratch.

### Patient variant coverage (step 10)

If `/scr1/users/wangf3/TEQUILA_Epilepsy_Hair/variants/patient_variants.tsv` exists, step 10 runs
`scripts/variant_coverage.py`, which counts how many gene reads per sample cover each known patient variant of the target
genes. The file needs columns gene, nucleotide (c.), protein and pathogenicity, separated by tabs or `|`.
- c. positions are read on the MANE Select transcript, and the reference base is checked against the genome.
- For splice-site variants, reads splicing across the variant's intron also count.

The output, in `results/figures/variants/`, is a table, a figure per gene and a summary. The variant table holds patient
information: keep it on scratch, never in this repository.

### Step 8. Expression QC of the full run (panel on-target rate, PCA, correlation)

`11_expression_qc.sh` reads every read of the genome-wide BAMs (`qc/bam/`, step 07). It runs as a job:
1. `scripts/gene_counts.py` counts reads per gene (assigned by alignment midpoint) and the on-target rate on the
   634 panel genes. It writes `qc/expression/gene_counts.tsv` and `panel_on_target.tsv`, and lists any panel symbols
   missing from the GTF in `panel_genes_missing.txt`.
2. `scripts/expression_qc.py` draws the panel on-target rate and where the reads go, PCA (all genes and panel
   genes), sample correlation, and panel genes Zymo vs Genotek per donor, in `results/figures/expression_qc/`.

Run step 10 again afterwards to add these slides to the deck.

```bash
cd /home/wangf3/feng/projects/TEQUILA_Epilepsy_Hair/TEQUILA-seq/TEQUILA_Epilepsy_Hair
bash 11_expression_qc.sh
```

### Folder layout on scratch

```
/scr1/users/wangf3/TEQUILA_Epilepsy_Hair/
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
    slides/                         TEQUILA_Epilepsy_Hair_pilot.pptx (step 10)
    figures/<GENE>/                 coverage and isoform figures (pdf, png, data tables) and igv/ pages per sample
  logs/                 one .out/.err per SLURM job (job names start with epi_)
/scr1/users/wangf3/TEQUILA_references/   GRCh38, GENCODE v49, minimap2 index
```

## Notes

- **Full length** is judged against the canonical transcript. A read of another isoform that skips
  the canonical first or last exon is not counted as full length; the isoform reports cover those.
- The scripts were tested in a sandbox on synthetic reads (coverage, full-length counts and the
  kit contrast came out as planted) and with a stand-in `sbatch`. They have not run on the real data.
