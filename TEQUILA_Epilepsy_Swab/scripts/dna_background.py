#!/usr/bin/env python3
"""How many of a gene's reads could be genomic DNA? An upper bound per gene and swab.

The off-panel reads of the swabs look largely like genomic DNA (12_rna_composition.sh). Within each target gene:
  * DNA density: reads lying wholly inside the gene's introns (no exonic base, no splice), per intronic base, from
    the genome-wide BAM of the full path. Intronic reads are DNA or pre-mRNA, so this is an upper bound for DNA.
  * DNA reads expected among the gene reads: density x the gene's exonic length (merged exons of all transcripts),
    i.e. fragments of DNA that fall on exons.
  * gene reads (gene-filtered BAM of step 09): spliced reads (with an intron) are certainly RNA; unspliced reads
    are RNA or DNA.
Outputs: dna_background.tsv (gene, sample, gene reads, spliced, unspliced, intronic reads, intronic bp, exonic bp,
expected DNA reads, upper-bound DNA share of gene reads) and dna_background.{pdf,png}.
"""
import argparse
import csv
import os
import sys
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figure_style import ATTR  # noqa: E402

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
    "font.size": 9, "axes.labelsize": 9, "xtick.labelsize": 8, "ytick.labelsize": 8, "axes.linewidth": 0.8,
    "pdf.fonttype": 42, "savefig.facecolor": "white"})


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--genes", required=True)
    p.add_argument("--gtf-pattern", required=True, help="per-gene locus GTF with {gene} (step 09 <GENE>_locus.gtf)")
    p.add_argument("--full-bam-pattern", required=True, help="genome-wide BAM with {sample}")
    p.add_argument("--gene-bam-pattern", required=True, help="gene-filtered BAM with {gene} and {sample}")
    p.add_argument("--sample-sheet", required=True)
    p.add_argument("--outdir", required=True)
    return p.parse_args()


def merge(iv):
    out = []
    for s, e in sorted(iv):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


def gene_model(gtf, gene):
    exons, span = [], None
    with open(gtf) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or dict(ATTR.findall(f[8])).get("gene_name") != gene:
                continue
            if f[2] == "gene":
                span = (f[0], int(f[3]) - 1, int(f[4]))
            elif f[2] == "exon":
                exons.append((int(f[3]) - 1, int(f[4])))
    return span, merge(exons)


def overlaps_any(blocks, ivs):
    return any(s < ie and is_ < e for s, e in blocks for is_, ie in ivs)


def main():
    a = parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    with open(a.sample_sheet) as fh:
        samples = [r["sample_id"] for r in csv.DictReader(fh, delimiter="\t")]
    rows = []
    for g in a.genes.split(","):
        gtf = a.gtf_pattern.format(gene=g)
        if not os.path.exists(gtf):
            continue
        span, exons = gene_model(gtf, g)
        if not span:
            continue
        chrom, gs, ge = span
        exon_bp = sum(e - s for s, e in exons)
        intron_bp = (ge - gs) - exon_bp
        for s in samples:
            fb, gb = a.full_bam_pattern.format(sample=s), a.gene_bam_pattern.format(gene=g, sample=s)
            if not (os.path.exists(fb) and os.path.exists(gb)):
                continue
            n_intr = 0
            with pysam.AlignmentFile(fb) as bam:
                for r in bam.fetch(chrom, gs, ge):
                    if r.is_secondary or r.is_supplementary or r.is_unmapped:
                        continue
                    if r.reference_start < gs or r.reference_end > ge:
                        continue
                    if any(op == 3 for op, _ in (r.cigartuples or ())):
                        continue
                    if not overlaps_any(r.get_blocks(), exons):
                        n_intr += 1
            spl = uns = 0
            with pysam.AlignmentFile(gb) as bam:
                for r in bam.fetch(until_eof=True):
                    if r.is_secondary or r.is_supplementary or r.is_unmapped:
                        continue
                    if any(op == 3 and ln >= 30 for op, ln in (r.cigartuples or ())):
                        spl += 1
                    else:
                        uns += 1
            dens = n_intr / intron_bp if intron_bp > 0 else 0
            exp = dens * exon_bp
            tot = spl + uns
            rows.append((g, s, tot, spl, uns, n_intr, intron_bp, exon_bp, exp,
                         min(1.0, exp / tot) if tot else np.nan))
    with open(os.path.join(a.outdir, "dna_background.tsv"), "w") as fh:
        fh.write("gene\tsample_id\tgene_reads\tspliced\tunspliced\tintronic_reads\tintronic_bp\texonic_bp\t"
                 "expected_dna_reads\tmax_dna_share\n")
        for r in rows:
            fh.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.1f}\t{:.3f}\n".format(*r))
    genes = list(dict.fromkeys(r[0] for r in rows))
    if not genes:
        return
    fig, axes = plt.subplots(1, len(genes), figsize=(3.6 * len(genes) + 1.0, 3.8), squeeze=False)
    fig.subplots_adjust(left=0.08, right=0.98, bottom=0.3, top=0.86, wspace=0.35)
    for ax, g in zip(axes[0], genes):
        rr = [r for r in rows if r[0] == g]
        x = np.arange(len(rr))
        ax.bar(x, [r[3] for r in rr], color="#1B7837", width=0.75, label="spliced (RNA)")
        ax.bar(x, [r[4] for r in rr], bottom=[r[3] for r in rr], color="#A6DBA0", width=0.75,
               label="unspliced (RNA or DNA)")
        ax.scatter(x, [r[8] for r in rr], marker="_", s=180, color="#B2182B", lw=2, zorder=3,
                   label="DNA expected (upper bound)")
        for xi, r in zip(x, rr):
            ax.text(xi, (r[3] + r[4]) * 1.02, "{:.0%}".format(r[9]) if not np.isnan(r[9]) else "", ha="center",
                    fontsize=6.5, color="#B2182B")
        ax.set_xticks(x)
        ax.set_xticklabels([r[1].replace("_", " ") for r in rr], rotation=40, ha="right", rotation_mode="anchor")
        ax.set_ylabel("Gene reads")
        ax.set_title(g, loc="left", fontsize=10, fontweight="bold", fontstyle="italic")
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    axes[0][-1].legend(frameon=False, fontsize=7.5, loc="upper left", bbox_to_anchor=(1.0, 1.0))
    fig.text(0.08, 0.01, "Red line and %: DNA reads expected on the exons from the density of intronic reads "
                         "(intronic reads are DNA or pre-mRNA, so an upper bound).", fontsize=7.5, color="#555555")
    fig.savefig(os.path.join(a.outdir, "dna_background.pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(a.outdir, "dna_background.png"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("DNA background -> {}".format(a.outdir))
    for g in genes:
        sh = [r[9] for r in rows if r[0] == g and not np.isnan(r[9])]
        print("  {}: upper-bound DNA share of gene reads {:.0%}-{:.0%}".format(g, min(sh), max(sh)))


if __name__ == "__main__":
    main()
