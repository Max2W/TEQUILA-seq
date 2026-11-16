#!/usr/bin/env python3
"""Fragmentation in one figure: for each gene, the fraction of the gene's reads that cover each position of the
canonical (MANE) transcript, 5' to 3', one line per sample (colour: kit; line style: donor).

Input: gene_body_coverage.tsv and fragmentation_summary.tsv of scripts/fragment_report.py (step 09).
Output: genebody.{pdf,png} and genebody.tsv (the drawn values, plus reads and full-length reads per gene and sample).
"""
import argparse
import csv
import os
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
    "font.size": 9, "axes.labelsize": 10, "xtick.labelsize": 9, "ytick.labelsize": 9, "axes.linewidth": 0.8,
    "pdf.fonttype": 42, "savefig.facecolor": "white"})
KIT_COLORS = ["#3188BD", "#EC7B02", "#57B45B", "#C9A0DC"]
STYLES = ["-", "--", ":", "-."]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--frag-dir", required=True, help="results/fragmentation/filtered")
    p.add_argument("--sample-sheet", required=True)
    p.add_argument("--genes", required=True)
    p.add_argument("--min-reads", type=int, default=20)
    p.add_argument("--outdir", required=True)
    return p.parse_args()


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main():
    a = parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    sheet = read_tsv(a.sample_sheet)
    kits = list(dict.fromkeys(r["condition"] for r in sheet))
    donors = list(dict.fromkeys(r["donor"] for r in sheet))
    meta = {r["sample_id"]: r for r in sheet}
    prof = defaultdict(lambda: np.zeros(100))
    for r in read_tsv(os.path.join(a.frag_dir, "gene_body_coverage.tsv")):
        prof[(r["gene"], r["sample_id"])][int(r["bin_pct_5to3"])] = float(r["fraction_of_reads_covering"])
    summ = {(r["gene"], r["sample_id"]): r for r in read_tsv(os.path.join(a.frag_dir, "fragmentation_summary.tsv"))}
    genes = [g for g in a.genes.split(",") if any(k[0] == g for k in prof)]
    nc = min(len(genes), 2) if len(genes) <= 2 else 3
    nr = int(np.ceil(len(genes) / nc))
    fig, axes = plt.subplots(nr, nc, figsize=(5.4 * nc + 1.6, 3.6 * nr), squeeze=False)
    fig.subplots_adjust(left=0.07, right=0.84, bottom=0.14 / nr + 0.06, top=0.9, wspace=0.25, hspace=0.45)
    out = []
    x = np.arange(100) + 0.5
    for ax, g in zip(axes.flat, genes):
        fl, nreads = 0, 0
        for s, m in meta.items():
            row = summ.get((g, s))
            n = int(row["reads_on_canonical_exons"]) if row else 0
            if (g, s) not in prof or n < a.min_reads:
                continue
            fl += int(row["full_length_reads"])
            nreads += n
            ax.plot(x, prof[(g, s)], color=KIT_COLORS[kits.index(m["condition"]) % 4],
                    ls=STYLES[donors.index(m["donor"]) % 4], lw=1.3)
            out += [(g, s, i, "{:.4f}".format(v)) for i, v in enumerate(prof[(g, s)])]
        ax.set_xlim(0, 100)
        ax.set_ylim(0, 1.0)
        ax.set_xlabel("Position along the canonical transcript (% of length, 5' to 3')")
        ax.set_ylabel("Fraction of reads covering")
        ax.set_title("{}  ({:,} reads, {:,} full length)".format(g, nreads, fl), loc="left", fontsize=10,
                     fontstyle="italic")
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    for ax in list(axes.flat)[len(genes):]:
        ax.axis("off")
    handles = [plt.Line2D([], [], color=KIT_COLORS[i % 4], lw=2, label=k) for i, k in enumerate(kits)]
    handles += [plt.Line2D([], [], color="#555555", ls=STYLES[i % 4], lw=1.3, label="Donor {}".format(d))
                for i, d in enumerate(donors)]
    fig.legend(handles=handles, loc="center left", bbox_to_anchor=(0.85, 0.5), frameon=False, fontsize=9)
    fig.text(0.07, 0.01, "Swabs with >= {} reads on the canonical exons. A full-length read reaches the first and the "
                         "last canonical exon.".format(a.min_reads), fontsize=8, color="#555555")
    fig.savefig(os.path.join(a.outdir, "genebody.pdf"))
    fig.savefig(os.path.join(a.outdir, "genebody.png"), dpi=300)
    plt.close(fig)
    with open(os.path.join(a.outdir, "genebody.tsv"), "w") as fh:
        fh.write("gene\tsample_id\tbin_pct_5to3\tfraction_of_reads_covering\n")
        for r in out:
            fh.write("\t".join(map(str, r)) + "\n")
    print("gene-body coverage figure -> {}".format(a.outdir))


if __name__ == "__main__":
    main()
