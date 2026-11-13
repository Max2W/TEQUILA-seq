#!/usr/bin/env python3
"""Figures from the full path (all reads aligned to the genome; 07_qc_full.sh), in the progress-deck style.

  capture_qc        per sample: (a) mapping rate, (b) reads on each target gene as % of all reads (on-target rate),
                    (c) median read length
  top_genes         genes collecting the most mapped reads, all samples pooled; target genes highlighted
  bait_recovery     gene-filtered reads found by the fast path (k-mer baiting, step 03) against those found by
                    aligning every read (full path), per gene and sample; the same for spliced reads
Inputs: qc/per_sample_qc.tsv and qc/top_genes.tsv (07_qc_full.sh), the gene read filter tables of the fast path
(results/figures/read_filter.tsv) and of the full path (full_read_filter.tsv, same filter on qc/bam/<sample>.bam).
Every drawn number goes to full_path_figures.tsv.
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
    "xtick.major.size": 3, "ytick.major.size": 3, "pdf.fonttype": 42, "savefig.facecolor": "white"})
SLATE = "#708090"
GENE_COLORS = [SLATE, "#E2C485", "#3A8DC4", "#57B45B"]
KIT_COLORS = ["#3188BD", "#EC7B02"]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--genes", required=True)
    p.add_argument("--sample-sheet", required=True)
    p.add_argument("--qc-dir", required=True)
    p.add_argument("--fast-filter", required=True, help="results/figures/read_filter.tsv")
    p.add_argument("--full-filter", help="the same table for the full-path BAMs")
    p.add_argument("--top-n", type=int, default=15)
    p.add_argument("--outdir", required=True)
    return p.parse_args()


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def classic(ax):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def labels(ax, names):
    ax.set_xticks(np.arange(len(names)))
    ax.set_xticklabels([n.replace("_", " ") for n in names], rotation=-30, ha="left", rotation_mode="anchor")


def save(fig, path):
    fig.savefig(path + ".pdf")
    fig.savefig(path + ".png", dpi=300)
    plt.close(fig)


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    genes = args.genes.split(",")
    sheet = read_tsv(args.sample_sheet)
    samples = [r["sample_id"] for r in sheet]
    kits = list(dict.fromkeys(r["condition"] for r in sheet))
    kit_of = {r["sample_id"]: r["condition"] for r in sheet}
    rows = []

    # ---- capture QC
    qc = {r["sample_id"]: r for r in read_tsv(os.path.join(args.qc_dir, "per_sample_qc.tsv"))}
    order = [s for s in samples if s in qc]
    fig, axes = plt.subplots(1, 3, figsize=(12.3, 3.6), gridspec_kw={"width_ratios": [1, 1.3, 1]})
    fig.subplots_adjust(left=0.06, right=0.99, top=0.88, bottom=0.25, wspace=0.35)
    x = np.arange(len(order))
    colors = [KIT_COLORS[kits.index(kit_of[s]) % 2] for s in order]
    mr = [float(qc[s]["mapping_rate"]) * 100 for s in order]
    axes[0].bar(x, mr, width=0.8, color=colors)
    axes[0].set_ylim(0, 100)
    axes[0].set_ylabel("Reads mapped to the genome (%)")
    w = 0.8 / len(genes)
    for k, g in enumerate(genes):
        v = [float(qc[s].get("{}_rate".format(g), 0)) * 100 for s in order]
        axes[1].bar(x - 0.4 + w * (k + 0.5), v, width=w, color=GENE_COLORS[k % 4], label=g)
        rows += [("on_target_percent", g, s, "{:.4f}".format(val)) for s, val in zip(order, v)]
    leg = axes[1].legend(frameon=False, fontsize=9, loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=len(genes))
    for t in leg.get_texts():
        t.set_fontstyle("italic")
    axes[1].set_ylabel("Reads on the gene (% of all reads)")
    ml = [float(qc[s]["median_read_length"]) for s in order]
    axes[2].bar(x, ml, width=0.8, color=colors)
    axes[2].set_ylabel("Median read length (nt)")
    for i, (ax, letter) in enumerate(zip(axes, "abc")):
        labels(ax, order)
        classic(ax)
        ax.set_title(letter, loc="left", fontweight="bold")
    handles = [plt.Rectangle((0, 0), 1, 1, color=KIT_COLORS[i % 2]) for i in range(len(kits))]
    axes[2].legend(handles, kits, frameon=False, fontsize=9, loc="upper right", bbox_to_anchor=(1.0, 1.15), ncol=2)
    save(fig, os.path.join(args.outdir, "capture_qc"))
    rows += [("mapping_rate_percent", "", s, "{:.2f}".format(v)) for s, v in zip(order, mr)]
    rows += [("median_read_length", "", s, "{:.0f}".format(v)) for s, v in zip(order, ml)]

    # ---- top genes, pooled over samples
    top = defaultdict(int)
    for r in read_tsv(os.path.join(args.qc_dir, "top_genes.tsv")):
        if r["flowcell"] == "all" and r["sample_id"] in qc:
            top[r["gene"]] += int(r["reads"])
    mapped = sum(int(qc[s]["mapped"]) for s in order)
    best = sorted(top.items(), key=lambda kv: -kv[1])[:args.top_n]
    fig, ax = plt.subplots(figsize=(7.0, 4.6))
    fig.subplots_adjust(left=0.22, right=0.97, top=0.97, bottom=0.12)
    y = np.arange(len(best))[::-1]
    ax.barh(y, [n / mapped * 100 for _, n in best], color=[GENE_COLORS[1] if g in genes else SLATE for g, _ in best],
            height=0.75)
    ax.set_yticks(y)
    ax.set_yticklabels([g for g, _ in best], fontstyle="italic")
    ax.set_xlabel("Reads assigned to the gene (% of mapped reads, all swabs)")
    classic(ax)
    save(fig, os.path.join(args.outdir, "top_genes"))
    rows += [("top_gene_percent_of_mapped", g, "all", "{:.4f}".format(n / mapped * 100)) for g, n in best]

    # ---- fast-path recovery
    if args.full_filter and os.path.exists(args.full_filter):
        fast = {(r["gene"], r["sample_id"]): r for r in read_tsv(args.fast_filter)}
        full = {(r["gene"], r["sample_id"]): r for r in read_tsv(args.full_filter)}
        fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.0))
        fig.subplots_adjust(left=0.09, right=0.98, top=0.9, bottom=0.14, wspace=0.3)
        for ax, col, what, letter in ((axes[0], "kept", "gene reads", "a"),
                                      (axes[1], "spliced_kept", "spliced gene reads", "b")):
            lo, hi = 1, 1
            for k, g in enumerate(genes):
                pts = [(int(full[(g, s)][col]), int(fast[(g, s)][col])) for s in samples
                       if (g, s) in full and (g, s) in fast]
                if not pts:
                    continue
                a, b = np.array(pts, float).T
                ax.scatter(np.maximum(a, 0.8), np.maximum(b, 0.8), s=28, color=GENE_COLORS[k % 4], label=g,
                           edgecolor="black", linewidth=0.4, zorder=3)
                hi = max(hi, a.max(), b.max())
                rows += [("{}_full_vs_fast".format(col), g, s, "{}/{}".format(int(fa), int(fb)))
                         for s, (fa, fb) in zip([s for s in samples if (g, s) in full and (g, s) in fast], pts)]
            ax.plot([0.8, hi * 1.5], [0.8, hi * 1.5], color="#999999", linewidth=0.8, linestyle="--", zorder=1)
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.set_xlim(0.8, hi * 1.5)
            ax.set_ylim(0.8, hi * 1.5)
            ax.set_xlabel("Full path: all reads aligned ({})".format(what))
            ax.set_ylabel("Fast path: k-mer baited ({})".format(what))
            ax.set_title(letter, loc="left", fontweight="bold")
            leg = ax.legend(frameon=False, fontsize=9, loc="upper left")
            for t in leg.get_texts():
                t.set_fontstyle("italic")
            classic(ax)
        save(fig, os.path.join(args.outdir, "bait_recovery"))

    with open(os.path.join(args.outdir, "full_path_figures.tsv"), "w") as fh:
        fh.write("measure\tgene\tsample_id\tvalue\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")
    print("full-path figures -> {}".format(args.outdir))


if __name__ == "__main__":
    main()
