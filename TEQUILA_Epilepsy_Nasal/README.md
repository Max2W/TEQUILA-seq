# TEQUILA-seq: epilepsy genes in nasal swab RNA (pilot)

**Status: planned.** Samples are not collected yet. The run folder and barcodes read `SET_ME` in `config.sh` and
`samples.tsv`, and step 01 stops with a message until they are set.

The buccal-swab pilot (branch `Epilepsy_Swab`) detected STXBP1 and SYNGAP1 in every swab, but the RNA was too
fragmented for full-length isoforms. The literature points to the nose as the better accessible tissue:
- **Integrity:** RNA from a brush on the inferior turbinate is intact (median RIN 8.9 in Lai 2015), because it
  comes from living airway epithelium, not from shed squamous cells.
- **Expression:** public nasal brushings of 695 children and young people express both genes in every sample
  (GEO GSE152004; median 37 CPM for STXBP1 and 71 CPM for SYNGAP1). Figures: `scripts/nasal_literature_figures.py`,
  with numbers in `public_data/`.

This pilot uses the swab kits already in the lab: the Zymo R1107 kit contains a flocked swab, and Zymo lists the
nose as a sample site. Questions:

1. **Integrity.** Is nasal swab RNA more intact than the same donor's buccal swab RNA (RIN, transcript coverage,
   full-length reads)?
2. **Detection.** Are the target genes and the panel detected, at what depth, and how many patient variants are
   covered?
3. **Stabiliser.** Does the Genotek stabiliser do as well as Zymo DNA/RNA Shield for nasal swabs?

## Samples

Same donors A-D as the buccal pilot, one session each, one nostril per sample:

| sample_id | donor | swab and site | stabiliser | RNA extraction | condition |
|---|---|---|---|---|---|
| Zymo_A ... Zymo_D | A-D | Zymo R1107 flocked swab, inferior turbinate, left | Zymo DNA/RNA Shield (1 mL, kit tube) | Zymo Quick-RNA Microprep | Zymo |
| Genotek_A ... Genotek_D | A-D | Zymo R1107 flocked swab, inferior turbinate, right | Genotek OR-100 tube (sponge not used) | Qiagen miRNeasy | Genotek |

Sample IDs match the buccal pilot (`<condition>_<donor>`). Step 10 therefore draws each donor's nasal and buccal
samples side by side (`scripts/compare_specimens.py`, read only on `/scr1/users/wangf3/TEQUILA_Epilepsy_Swab`).

**Check the Genotek arm on the bench first,** with one spare kit:
- The OR-100 sponge is for the mouth and is not for the nose.
- Break a Zymo swab at its breakpoint into the Genotek tube, cut the sponge off the Genotek cap, and close it.
- If the swab does not fit or the cap does not seal, use a second Zymo kit for the right nostril instead, and set
  `condition` to `Zymo2`.

## Collection (adult lab donors)

Each step names its source. Sources are in `public_data/nasal_literature.tsv` and `nasal_pediatric.tsv`.

1. **Do not sample** after a cold or other respiratory infection in the last 2 weeks (Chu 2016 postponed sampling
   by about 2 weeks). Also do not sample anyone with nosebleeds, a bleeding disorder or anticoagulants, or nasal
   surgery or polyps.
2. **Prepare.** The donor blows their nose gently (Brewington 2018). They sit with the head slightly back, under
   good light.
3. **Left nostril, Zymo tube.**
   - Insert the flocked swab along the floor of the nostril, parallel to the palate, towards the inferior
     turbinate (the lower fold on the side wall). Stop at gentle resistance and do not force it.
   - Press the tip against the turbinate and rotate it for 3-5 s, moving it slightly in and out 4-5 times (Chu
     2016; Brewington 2018).
4. **Store at once.** Put the swab straight into its tube, break it at the breakpoint, close the cap and invert
   5 times.
5. **Right nostril, Genotek tube (or the second Zymo kit).** Take a new swab the same way. Do not sample the same
   nostril twice (Ulm 2016).
6. **Record:**
   - time of collection and side;
   - discomfort from 0 to 10 (the scale of Lai 2015, where turbinate brushing scored a median of 3.5);
   - blood on the swab, and any nosebleed;
   - watch the donor for 5 minutes (Chu 2016).
7. **Storage.** Keep the Zymo tubes at room temperature (Zymo states RNA stability for at least 1 month). Keep the
   Genotek tubes at room temperature and not at 4 C (Genotek). Extract both on the same day, within a week.

## RNA extraction and QC

- **Zymo tube:**
  - Use the whole 1 mL of DNA/RNA Shield. The nasal cell pellet is small and the Shield dilutes it.
  - Optionally digest with Proteinase K in the Shield first. Then add 1 volume of RNA Lysis Buffer and load the
    Quick-RNA Microprep column in several passes.
  - Do the on-column DNase I step, and elute in about 15 uL.
- **Genotek tube:** use the same Qiagen miRNeasy procedure as for the buccal Genotek swabs, so that only the
  sampling site changes.
- **QC:**
  - RIN on a high-sensitivity chip, and yield by Qubit RNA HS.
  - Write the RIN values to `/scr1/users/wangf3/TEQUILA_Epilepsy_Nasal/rin.tsv` (columns `sample_id`, `RIN`).
    Sample data stay on scratch, not in this repository.
- **Library:** use the same TEQUILA-seq panel, protocol and input as the buccal pilot, with 8 barcodes on one flow
  cell.

## Bench add-on: is buccal RNA degraded inside cells or outside them? (no sequencing)

Sequencing cannot tell cell-free RNA from RNA inside cells. Take one extra Zymo buccal swab per donor and split it:
1. Vortex the tube, press the swab against the wall and remove it.
2. Spin at low speed to pellet the cells (for example 500 x g, 5 min). Move the supernatant to a new tube.
3. Clear the supernatant at high speed (for example 16,000 x g, 10 min).
4. Extract the cell pellet (resuspended in Shield) and the cleared supernatant separately, and measure RIN and yield
   of each.

If most of the RNA, and the low RIN, sit in the pellet, the RNA is degraded inside the shed cells. If the
supernatant holds a large share, extracellular RNA matters. The spin speeds are a suggestion, not a published
protocol.

## Decision rules (set before the data)

The nasal swab is taken forward to children when, in at least 3 of 4 donors, it beats the same donor's buccal swab
on all of these:
- STXBP1 median transcript coverage per read and full-length STXBP1 reads per million are higher (step 08 /
  `compare_specimens.py`);
- RIN is at least 2 units higher;
- STXBP1 and SYNGAP1 gene reads per million are not more than 2-fold lower;
- discomfort is 3 or less, and no nosebleed needs treatment.

If integrity does not improve, the next test is a cytology brush on the inferior turbinate (Lai 2015: RIN 8.9 vs 1.1
for a swab at the front of the nose). It needs brushes bought.

## Genes

The probe panel is the 634-gene Epilepsy & NDD panel (`panel_genes.txt`, from Epilepsy_NDD634.xlsx). The genes of
interest are STXBP1 and SYNGAP1, listed in `target_genes.txt`; every step reads that file. The panel list is
used for the on-target rate (step 11).

## Where things are on the HPC

- **Code** (this repository) is in home: `/home/wangf3/feng/projects/TEQUILA_Epilepsy_Nasal/TEQUILA-seq`.
- **Outputs, logs and working files** are on scratch: `/scr1/users/wangf3/TEQUILA_Epilepsy_Nasal`.
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
mkdir -p /home/wangf3/feng/projects/TEQUILA_Epilepsy_Nasal
cd /home/wangf3/feng/projects/TEQUILA_Epilepsy_Nasal
git clone -b Epilepsy_Nasal https://github.com/Max2W/TEQUILA-seq.git TEQUILA-seq
```

If it is already there, update it instead:

```bash
cd /home/wangf3/feng/projects/TEQUILA_Epilepsy_Nasal/TEQUILA-seq
git pull
```

### Step 2. Set up (seconds)

It should end with `=== SETUP COMPLETE ===`:

```bash
cd /home/wangf3/feng/projects/TEQUILA_Epilepsy_Nasal/TEQUILA-seq/TEQUILA_Epilepsy_Nasal
bash 00_setup.sh
```

### Step 3. Start both paths

```bash
cd /home/wangf3/feng/projects/TEQUILA_Epilepsy_Nasal/TEQUILA-seq/TEQUILA_Epilepsy_Nasal
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
cd /scr1/users/wangf3/TEQUILA_Epilepsy_Nasal/logs
for f in $(ls -t *.err | head -3); do echo "=== $f"; tail -15 "$f"; done
```

### Step 5. Results

```bash
cd /scr1/users/wangf3/TEQUILA_Epilepsy_Nasal
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
cd /home/wangf3/feng/projects/TEQUILA_Epilepsy_Nasal/TEQUILA-seq/TEQUILA_Epilepsy_Nasal
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
cd /home/wangf3/feng/projects/TEQUILA_Epilepsy_Nasal/TEQUILA-seq/TEQUILA_Epilepsy_Nasal
LOCAL=1 bash 10_slides.sh
```

The deck is `results/slides/TEQUILA_Epilepsy_Nasal_pilot.pptx`.

The deck also has public-data slides from `scripts/public_data_figures.py`. Their inputs are committed in
`public_data/`, so no internet is needed:
- GTEx v10 and Human Protein Atlas expression of both genes;
- swab RNA quality from Martorella et al., BMC Genomics 2023: RIN, library QC failures and OMIM gene capture.

When the full path has run (`07_qc_full.sh`), step 10 also applies the gene read filter to the genome-wide BAMs
(`results/figures/full_read_filter.tsv`). It adds capture QC slides from `scripts/full_path_figures.py`:
- mapping rate, on-target rate and read length;
- the genes collecting the most reads;
- how many gene reads the k-mer baiting of the fast path recovered.

With `/scr1/users/wangf3/TEQUILA_Epilepsy_Nasal/rin.tsv` (columns `sample_id`, `RIN`; sample data, kept on scratch),
step 10 adds:
- the pilot's RIN to the public-data RIN figure and slide;
- RIN against the read-level integrity of the target genes (`scripts/rin_figures.py`).

When step 11 has run, step 10 also adds a sex-marker sample identity check (`scripts/sex_check.py`; XIST against
Y-linked genes; the two swabs of a donor must agree). Its output stays with the results on scratch.

### Patient variant coverage (step 10)

If `/scr1/users/wangf3/TEQUILA_Epilepsy_Nasal/variants/patient_variants.tsv` exists, step 10 runs
`scripts/variant_coverage.py`, which counts how many gene reads per swab cover each known patient variant of the target
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
cd /home/wangf3/feng/projects/TEQUILA_Epilepsy_Nasal/TEQUILA-seq/TEQUILA_Epilepsy_Nasal
bash 11_expression_qc.sh
```

### Step 9. RNA composition: what kind of RNA the swab reads come from

The question: is the low-RIN swab RNA mature mRNA from cells, broken into pieces, or something else? `12_rna_composition.sh`
reads every read of the genome-wide BAMs (`qc/bam/`, step 07) and runs as a job (8 CPUs):

1. `scripts/rna_composition.py` puts each read in one category: exonic spliced (processed mRNA), exonic unspliced,
   intronic (pre-mRNA or DNA), intergenic, mitochondrial, rRNA, or unmapped. Each read is also marked panel or
   off-panel; off-panel reads were not selected by the probes, so they are the closest view of the input RNA. It
   also records read lengths per category and the coverage of reads along canonical transcripts of at least 1.5 kb
   outside the panel (5' to 3'). The tables are written to `qc/composition/`.
2. `scripts/composition_figures.py` draws the categories, read lengths, transcript coverage, and cell-type marker
   genes (squamous epithelium, salivary gland, immune cells, airway epithelium, olfactory neurons) in
   `results/figures/composition/`. The marker figure needs `qc/expression/gene_counts.tsv` from step 8.

Reads cannot tell RNA inside cells from cell-free RNA. That needs a wet-lab test: spin the swab eluate, then extract
the cell pellet and the supernatant separately and measure the RIN and yield of each.

The same step runs unchanged on a nasal-brushing pilot, so the two specimens can be compared category by category.
Run step 10 again afterwards to add these slides to the deck.

```bash
cd /home/wangf3/feng/projects/TEQUILA_Epilepsy_Nasal/TEQUILA-seq/TEQUILA_Epilepsy_Nasal
bash 12_rna_composition.sh
```

When it finishes, paste back the output of:

```bash
column -t /scr1/users/wangf3/TEQUILA_Epilepsy_Nasal/qc/composition/transcript_span.tsv
```

### Folder layout on scratch

```
/scr1/users/wangf3/TEQUILA_Epilepsy_Nasal/
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
    slides/                         TEQUILA_Epilepsy_Nasal_pilot.pptx (step 10)
    figures/<GENE>/                 coverage and isoform figures (pdf, png, data tables) and igv/ pages per sample
  logs/                 one .out/.err per SLURM job (job names start with epi_)
/scr1/users/wangf3/TEQUILA_references/   GRCh38, GENCODE v49, minimap2 index
```

## Notes

- **Full length** is judged against the canonical transcript. A read of another isoform that skips
  the canonical first or last exon is not counted as full length; the isoform reports cover those.
- The scripts were tested in a sandbox on synthetic reads (coverage, full-length counts and the
  kit contrast came out as planted) and with a stand-in `sbatch`. They have not run on the real data.
