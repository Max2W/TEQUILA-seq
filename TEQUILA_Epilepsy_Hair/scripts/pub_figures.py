#!/usr/bin/env python3
"""Publication Figures 1 and 2 for the hair-follicle TEQUILA-seq pilot (STRIPE / Nature style, 183 mm wide).

Figure 1  Sequencing and detection
  a  reads sequenced per sample (bars) and mean read length (number above)
  b  read filtering per gene: reads at the locus kept for the gene vs removed, per sample
  c  detection: gene reads per million sequenced reads, Zymo vs Genotek, one line per donor
Figure 2  RNA integrity at the target genes
  a  full-length reads (reach the first and the last canonical exon) per sample
  b  read length of the gene's reads per sample

Every number drawn is written to figure1_data.tsv / figure2_data.tsv; kit_comparison.tsv holds the paired
(by donor) Genotek vs Zymo comparison of each measure: per-donor values, log2 ratio, mean and a paired t-test on the
log2 values (n = donors; exploratory).

Inputs: the sample sheet, qc/fastq_read_counts.tsv, the filter table from gene_read_filter.py and the outputs of
fragment_report.py run on the filtered BAMs.
"""
import argparse
import csv
import gzip
import math
import os
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.ticker  # noqa: E402
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle
from scipy import stats

KIT_COLORS = ["#3188BD", "#EC7B02", "#90AFB1", "#C2A5CF"]
REMOVED = "#BEBEBE"
NATURE_W_IN = 183 / 25.4
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Nimbus Sans", "Liberation Sans", "Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 6.5, "axes.linewidth": 0.5, "axes.edgecolor": "black", "xtick.major.width": 0.5,
    "ytick.major.width": 0.5, "xtick.major.size": 2, "ytick.major.size": 2, "axes.titlesize": 7,
    "axes.labelsize": 6.5, "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 5.5,
    "pdf.fonttype": 42, "savefig.facecolor": "white",
})


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sample-sheet", required=True)
    p.add_argument("--read-counts", required=True)
    p.add_argument("--filter-stats", required=True)
    p.add_argument("--frag-dir", required=True, help="fragment_report.py output on the filtered BAMs")
    p.add_argument("--genes", required=True, help="comma separated")
    p.add_argument("--outdir", required=True)
    return p.parse_args()


def theme(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def letter(ax, s, dx=-0.22, dy=1.04):
    ax.text(dx, dy, s, transform=ax.transAxes, fontsize=8, fontweight="bold", va="bottom", ha="left")


def read_tsv(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def paired_rows(values, samples, meta, kits):
    """values: sample -> number; returns per donor (ref, other) for the first two kits."""
    by = defaultdict(dict)
    for s in samples:
        if s in values:
            by[meta[s]["donor"]][meta[s]["condition"]] = values[s]
    return {d: (v.get(kits[0]), v.get(kits[1])) for d, v in by.items()}


def compare(measure, gene, pairs, kits, log=True):
    rows = []
    a, b = [], []
    for donor, (x, y) in sorted(pairs.items()):
        if x is None or y is None:
            continue
        rows.append((donor, x, y))
        a.append(x)
        b.append(y)
    out = {"measure": measure, "gene": gene, "n_donors": len(rows),
           "per_donor": ";".join("{}:{:.4g}/{:.4g}".format(d, x, y) for d, x, y in rows)}
    if len(rows) >= 2:
        if log:
            la = np.log2(np.array(a) + 1e-9)
            lb = np.log2(np.array(b) + 1e-9)
            diff = lb - la
            out["mean_log2_ratio_{}_vs_{}".format(kits[1], kits[0])] = round(float(diff.mean()), 3)
        else:
            diff = np.array(b) - np.array(a)
            out["mean_difference_{}_minus_{}".format(kits[1], kits[0])] = round(float(diff.mean()), 4)
        out["donors_higher_in_{}".format(kits[1])] = int((diff > 0).sum())
        out["paired_t_p"] = round(float(stats.ttest_1samp(diff, 0.0).pvalue), 4) if len(rows) >= 3 else ""
    return out


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    genes = [g for g in args.genes.split(",") if g]
    sheet = read_tsv(args.sample_sheet)
    samples = [r["sample_id"] for r in sheet]
    meta = {r["sample_id"]: r for r in sheet}
    kits = list(dict.fromkeys(r["condition"] for r in sheet))
    donors = list(dict.fromkeys(r["donor"] for r in sheet))
    color = {k: KIT_COLORS[i % 4] for i, k in enumerate(kits)}

    reads_total, bases_total = defaultdict(int), defaultdict(int)
    for r in read_tsv(args.read_counts):
        if r["sample_id"] in meta:
            reads_total[r["sample_id"]] += int(r["n_reads"])
            bases_total[r["sample_id"]] += int(r["n_bases"])
    filt = {(r["gene"], r["sample_id"]): r for r in read_tsv(args.filter_stats)}
    summ = {(r["gene"], r["sample_id"]): r for r in read_tsv(os.path.join(args.frag_dir, "fragmentation_summary.tsv"))}
    prof = defaultdict(list)
    for r in read_tsv(os.path.join(args.frag_dir, "gene_body_coverage.tsv")):
        prof[(r["gene"], r["sample_id"])].append(float(r["fraction_of_reads_covering"]))
    per_read = defaultdict(list)
    for r in read_tsv(os.path.join(args.frag_dir, "read_tx_coverage.tsv.gz")):
        per_read[(r["gene"], r["sample_id"])].append((int(r["read_length"]), float(r["tx_coverage"]),
                                                       r["tx_5prime_pos"] != "NA"))
    comparisons = []
    data1, data2 = [], []

    # ------------------------------------------------------------------ Figure 1
    W = NATURE_W_IN
    fig = plt.figure(figsize=(W, 2.5 + 0.0))
    gs = fig.add_gridspec(1, 2 + len(genes), width_ratios=[1.3, 1.6] + [0.62] * len(genes), wspace=0.6,
                          left=0.07, right=0.98, bottom=0.27, top=0.88)
    ax = fig.add_subplot(gs[0, 0])
    x = np.arange(len(samples))
    vals = [reads_total[s] / 1e6 for s in samples]
    ax.bar(x, vals, color=[color[meta[s]["condition"]] for s in samples], width=0.7, lw=0)
    for i, s in enumerate(samples):
        mean_len = bases_total[s] / reads_total[s] if reads_total[s] else 0
        ax.text(i, vals[i], "{:.0f}".format(mean_len), ha="center", va="bottom", fontsize=4.8, rotation=90)
        data1.append(("a", s, "reads_sequenced", reads_total[s]))
        data1.append(("a", s, "mean_read_length_nt", round(mean_len, 1)))
    ax.set_xticks(x)
    ax.set_xticklabels(samples, rotation=60, ha="right")
    ax.set_ylabel("reads sequenced (millions)")
    ax.set_ylim(0, max(vals) * 1.25)
    ax.set_title("number above bar:\nmean read length (nt)", fontsize=5, loc="left")
    theme(ax)
    letter(ax, "a")

    ax = fig.add_subplot(gs[0, 1])
    width = 0.38
    for gi, gene in enumerate(genes):
        kept = np.array([int(filt.get((gene, s), {}).get("kept", 0)) for s in samples])
        rin = np.array([int(filt.get((gene, s), {}).get("reads_in", 0)) for s in samples])
        xs = x + (gi - (len(genes) - 1) / 2) * width
        ax.bar(xs, kept, width=width * 0.9, color=[color[meta[s]["condition"]] for s in samples], lw=0,
               alpha=1.0 if gi == 0 else 0.55)
        ax.bar(xs, rin - kept, bottom=kept, width=width * 0.9, color=REMOVED, lw=0)
        for i, s in enumerate(samples):
            data1.append(("b", s, "{}_reads_at_locus".format(gene), int(rin[i])))
            data1.append(("b", s, "{}_reads_kept".format(gene), int(kept[i])))
    ax.set_xticks(x)
    ax.set_xticklabels(samples, rotation=60, ha="right")
    ax.set_ylabel("reads at the gene locus")
    ax.set_title(("left bar {}, right bar {}\n".format(*genes[:2]) if len(genes) >= 2 else "")
                 + "grey = removed by the gene read filter", fontsize=5, loc="left")
    theme(ax)
    letter(ax, "b")

    # c: gene reads per million
    panels = [("kept", "reads per million sequenced", "c", "reads_per_million")]
    for pi, (col, ylabel, panel, measure) in enumerate(panels):
        for gi, gene in enumerate(genes):
            ax = fig.add_subplot(gs[0, 2 + pi * len(genes) + gi])
            rpm = {s: 1e6 * int(filt.get((gene, s), {}).get(col, 0)) / reads_total[s]
                   for s in samples if reads_total[s]}
            pairs = paired_rows(rpm, samples, meta, kits)
            for d, (a, b) in pairs.items():
                if a is not None and b is not None:
                    ax.plot([0, 1], [a, b], color="#7F7F7F", lw=0.6, zorder=1)
                    ax.text(1.08, b, d, fontsize=5, va="center")
            for s, v in rpm.items():
                ax.scatter(kits.index(meta[s]["condition"]), v, s=14, color=color[meta[s]["condition"]], zorder=2,
                           lw=0)
                data1.append((panel, s, "{}_{}".format(gene, measure), round(v, 3)))
            positive = [v for v in rpm.values() if v > 0]
            if positive:
                ax.set_yscale("log")
                ax.yaxis.set_major_locator(matplotlib.ticker.LogLocator(base=10, subs=(1, 2, 5)))
                ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: "{:g}".format(v)))
                ax.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
                ax.set_ylim(min(positive) / 1.6, max(positive) * 1.6)
            ax.set_xticks([0, 1])
            ax.set_xticklabels(kits[:2], rotation=45, ha="right")
            ax.set_xlim(-0.4, 1.4)
            ax.set_title(gene, fontstyle="italic")
            if gi == 0:
                ax.set_ylabel(ylabel)
                letter(ax, panel, dx=-0.75)
            theme(ax)
            comparisons.append(compare(measure, gene, pairs, kits))
    fig.savefig(os.path.join(args.outdir, "figure1.pdf"))
    fig.savefig(os.path.join(args.outdir, "figure1.png"), dpi=300)
    plt.close(fig)

    # ------------------------------------------------------------------ Figure 2
    fig = plt.figure(figsize=(W, 2.9))
    gs = fig.add_gridspec(1, 2, wspace=0.28, left=0.07, right=0.98, bottom=0.27, top=0.86)
    width = 0.38
    x = np.arange(len(samples))
    # a: full-length reads per sample
    ax = fig.add_subplot(gs[0, 0])
    for gi, gene in enumerate(genes):
        vals = np.array([int(summ.get((gene, s), {}).get("full_length_reads", 0)) for s in samples])
        xs = x + (gi - (len(genes) - 1) / 2) * width
        ax.bar(xs, vals, width=width * 0.9, color=[color[meta[s]["condition"]] for s in samples], lw=0,
               alpha=1.0 if gi == 0 else 0.45)
        for i, s in enumerate(samples):
            ax.text(xs[i], vals[i], str(vals[i]), ha="center", va="bottom", fontsize=5)
            data2.append(("a", gene, s, "full_length_reads", int(vals[i])))
    ax.set_xticks(x)
    ax.set_xticklabels(samples, rotation=45, ha="right")
    ax.set_ylabel("full-length reads")
    if len(genes) >= 2:
        ax.legend(handles=[Rectangle((0, 0), 1, 1, color="#555555"), Rectangle((0, 0), 1, 1, color="#555555",
                                                                                alpha=0.45)],
                  labels=genes[:2], frameon=False, loc="upper right", fontsize=6, prop={"style": "italic",
                                                                                         "size": 6})
    theme(ax)
    letter(ax, "a", dx=-0.1)
    # b: read length of the gene reads
    ax = fig.add_subplot(gs[0, 1])
    pos, data, cols = [], [], []
    for i, s in enumerate(samples):
        for gi, gene in enumerate(genes):
            lens = [ln for (ln, _, _) in per_read.get((gene, s), [])]
            if not lens:
                continue
            pos.append(i + (gi - (len(genes) - 1) / 2) * width)
            data.append(lens)
            cols.append(color[meta[s]["condition"]])
            data2.append(("b", gene, s, "median_read_length_nt", float(np.median(lens))))
    if data:
        bp = ax.boxplot(data, positions=pos, widths=width * 0.8, showfliers=False, patch_artist=True,
                        medianprops=dict(color="black", lw=0.6), whiskerprops=dict(lw=0.5),
                        capprops=dict(lw=0.5), boxprops=dict(lw=0.4))
        for k, (patch, c) in enumerate(zip(bp["boxes"], cols)):
            patch.set_facecolor(c)
            patch.set_alpha(0.85 if k % len(genes) == 0 else 0.45)
    ax.set_xticks(x)
    ax.set_xticklabels(samples, rotation=45, ha="right")
    ax.set_xlim(-0.6, len(samples) - 0.4)
    ax.set_ylabel("read length (nt)")
    if len(genes) >= 2:
        ax.set_title("left box {}, right box {}".format(*genes[:2]), fontsize=6, loc="left")
    theme(ax)
    letter(ax, "b", dx=-0.1)
    fig.savefig(os.path.join(args.outdir, "figure2.pdf"))
    fig.savefig(os.path.join(args.outdir, "figure2.png"), dpi=300)
    plt.close(fig)

    for gene in genes:
        fl = {s: int(summ.get((gene, s), {}).get("full_length_reads", 0)) for s in samples}
        comparisons.append(compare("full_length_reads", gene, paired_rows(fl, samples, meta, kits), kits,
                                   log=False))
        med = {s: float(summ[(gene, s)]["median_tx_coverage"]) for s in samples if (gene, s) in summ}
        comparisons.append(compare("median_transcript_coverage", gene, paired_rows(med, samples, meta, kits), kits,
                                   log=False))

    with open(os.path.join(args.outdir, "figure1_data.tsv"), "w") as fh:
        fh.write("panel\tsample_id\tmeasure\tvalue\n")
        fh.writelines("\t".join(map(str, r)) + "\n" for r in data1)
    with open(os.path.join(args.outdir, "figure2_data.tsv"), "w") as fh:
        fh.write("panel\tgene\tgroup\tmeasure\tvalue\n")
        fh.writelines("\t".join(map(str, r)) + "\n" for r in data2)
    keys = []
    for c in comparisons:
        for k in c:
            if k not in keys:
                keys.append(k)
    with open(os.path.join(args.outdir, "kit_comparison.tsv"), "w") as fh:
        fh.write("\t".join(keys) + "\n")
        for c in comparisons:
            fh.write("\t".join(str(c.get(k, "")) for k in keys) + "\n")
    print("figures 1-2 -> {}".format(args.outdir))


if __name__ == "__main__":
    main()
