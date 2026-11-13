#!/usr/bin/env python3
"""Sample-level QC from the gene counts of the full path (scripts/gene_counts.py), in the progress-deck style.

  panel_capture     (a) on-target rate: reads overlapping any probe-panel gene, % of all reads; (b) where the reads
                    go: panel genes, other genes, outside annotated genes, unmapped (each read counted once, by its
                    alignment midpoint)
  pca               PCA of log2(CPM + 1), (a) all genes and (b) panel genes; genes with CPM >= 1 in >= 2 samples
  correlation       Pearson correlation between samples of log2(CPM + 1), all genes and panel genes
  kit_scatter       per donor, panel genes in Zymo against Genotek, log2(CPM + 1); target genes labelled
CPM uses the reads assigned to genes (intergenic and unmapped reads left out). Tables: expression_qc.tsv.
"""
import argparse
import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
    "font.size": 9, "axes.labelsize": 10, "xtick.labelsize": 9, "ytick.labelsize": 9, "axes.linewidth": 0.8,
    "xtick.major.size": 3, "ytick.major.size": 3, "pdf.fonttype": 42, "savefig.facecolor": "white"})
SLATE = "#708090"
KIT_COLORS = ["#3188BD", "#EC7B02"]
HEAT = LinearSegmentedColormap.from_list("deck", ["#F3E79B", "#FAB27B", "#EB7F86", "#CE6693", "#A059A0", "#5C53A5"])


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--counts-dir", required=True, help="output of scripts/gene_counts.py")
    p.add_argument("--sample-sheet", required=True)
    p.add_argument("--genes", required=True, help="target genes to label")
    p.add_argument("--outdir", required=True)
    return p.parse_args()


def classic(ax):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def save(fig, path):
    fig.savefig(path + ".pdf")
    fig.savefig(path + ".png", dpi=300)
    plt.close(fig)


def label(s):
    return s.replace("_", " ")


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    genes = args.genes.split(",")
    with open(args.sample_sheet) as fh:
        sheet = list(csv.DictReader(fh, delimiter="\t"))
    kits = list(dict.fromkeys(r["condition"] for r in sheet))
    kit_of = {r["sample_id"]: r["condition"] for r in sheet}
    donor_of = {r["sample_id"]: r["donor"] for r in sheet}
    tab = pd.read_csv(os.path.join(args.counts_dir, "gene_counts.tsv"), sep="\t", index_col=0)
    samples = [r["sample_id"] for r in sheet if r["sample_id"] in tab.columns]
    on = pd.read_csv(os.path.join(args.counts_dir, "panel_on_target.tsv"), sep="\t", index_col=0).loc[samples]
    rows = []

    # ---- capture
    counts = tab[samples]
    in_panel = tab["in_panel"].astype(str) == "True"
    genic = counts.drop(index=[g for g in ("intergenic", "unmapped") if g in counts.index])
    comp = pd.DataFrame({
        "Panel genes": counts[in_panel].sum(),
        "Other genes": genic[~in_panel.reindex(genic.index).fillna(False)].sum(),
        "Outside annotated genes": counts.loc["intergenic"] if "intergenic" in counts.index else 0,
        "Unmapped": counts.loc["unmapped"] if "unmapped" in counts.index else 0}).T
    comp_pct = comp / comp.sum() * 100
    fig, axes = plt.subplots(1, 2, figsize=(12.3, 3.8), gridspec_kw={"width_ratios": [1, 1.25]})
    fig.subplots_adjust(left=0.06, right=0.8, top=0.88, bottom=0.25, wspace=0.3)
    x = np.arange(len(samples))
    rate = on["panel_rate"].astype(float) * 100
    axes[0].bar(x, rate, width=0.8, color=[KIT_COLORS[kits.index(kit_of[s]) % 2] for s in samples])
    for i, v in enumerate(rate):
        axes[0].text(i, v, "{:.1f}".format(v), ha="center", va="bottom", fontsize=8)
    axes[0].set_ylabel("On-target rate (% of all reads\non {} panel genes)".format(int(in_panel.sum())))
    handles = [plt.Rectangle((0, 0), 1, 1, color=KIT_COLORS[i % 2]) for i in range(len(kits))]
    axes[0].legend(handles, kits, frameon=False, fontsize=9, loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2)
    colors = {"Panel genes": "#3A8DC4", "Other genes": SLATE, "Outside annotated genes": "#BEBEBE",
              "Unmapped": "#E6E6E6"}
    bottom = np.zeros(len(samples))
    for part in comp_pct.index:
        v = comp_pct.loc[part, samples].values
        axes[1].bar(x, v, width=0.8, bottom=bottom, color=colors[part], label=part)
        bottom += v
    axes[1].set_ylim(0, 100)
    axes[1].set_yticks([0, 25, 50, 75, 100])
    axes[1].set_ylabel("Reads (%)")
    axes[1].legend(frameon=False, fontsize=9, loc="center left", bbox_to_anchor=(1.0, 0.5))
    for ax, letter in zip(axes, "ab"):
        ax.set_xticks(x)
        ax.set_xticklabels([label(s) for s in samples], rotation=-30, ha="left", rotation_mode="anchor")
        ax.set_title(letter, loc="left", fontweight="bold")
        classic(ax)
    save(fig, os.path.join(args.outdir, "panel_capture"))
    for s in samples:
        rows.append(("on_target_percent", s, "{:.3f}".format(rate[s])))
        for part in comp_pct.index:
            rows.append(("reads_percent_" + part.replace(" ", "_"), s, "{:.3f}".format(comp_pct.loc[part, s])))

    # ---- expression matrices
    expr = genic.loc[:, samples]
    cpm = expr / expr.sum() * 1e6
    keep = (cpm >= 1).sum(axis=1) >= 2
    logc = np.log2(cpm[keep] + 1)
    panel_log = logc[in_panel.reindex(logc.index).fillna(False)]
    sets = [("All genes", logc), ("Panel genes", panel_log)]
    rows.append(("genes_used", "all", str(len(logc))))
    rows.append(("genes_used", "panel", str(len(panel_log))))

    # PCA
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2))
    fig.subplots_adjust(left=0.08, right=0.83, top=0.9, bottom=0.14, wspace=0.32)
    for ax, (name, m), letter in zip(axes, sets, "ab"):
        X = m.T.values - m.T.values.mean(axis=0)
        u, sv, _ = np.linalg.svd(X, full_matrices=False)
        pcs = u * sv
        var = sv ** 2 / (sv ** 2).sum() * 100
        for i, s in enumerate(samples):
            ax.scatter(pcs[i, 0], pcs[i, 1], s=70, color=KIT_COLORS[kits.index(kit_of[s]) % 2],
                       edgecolor="black", linewidth=0.5, zorder=3)
            ax.annotate(donor_of[s], (pcs[i, 0], pcs[i, 1]), xytext=(6, 4), textcoords="offset points", fontsize=9)
            rows.append(("PC1_PC2_" + name.split()[0].lower(), s, "{:.3f},{:.3f}".format(pcs[i, 0], pcs[i, 1])))
        ax.set_xlabel("PC1 ({:.1f}%)".format(var[0]))
        ax.set_ylabel("PC2 ({:.1f}%)".format(var[1]))
        ax.set_title("{}  {} ({:,} genes)".format(letter, name, len(m)), loc="left", fontweight="bold", fontsize=10)
        classic(ax)
    axes[1].legend([plt.Line2D([], [], marker="o", linestyle="", markersize=8, markerfacecolor=KIT_COLORS[i % 2],
                                markeredgecolor="black") for i in range(len(kits))], kits, frameon=False,
                   loc="center left", bbox_to_anchor=(1.02, 0.5), title="Kit (letter = donor)", title_fontsize=9)
    save(fig, os.path.join(args.outdir, "pca"))

    # correlation heatmaps
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.6))
    fig.subplots_adjust(left=0.09, right=0.9, top=0.9, bottom=0.2, wspace=0.45)
    for ax, (name, m), letter in zip(axes, sets, "ab"):
        r = np.corrcoef(m.T.values)
        im = ax.imshow(r, cmap=HEAT, vmin=min(0.5, np.floor(r.min() * 20) / 20), vmax=1)
        for i in range(len(samples)):
            for j in range(len(samples)):
                ax.text(j, i, "{:.2f}".format(r[i, j]), ha="center", va="center", fontsize=7,
                        color="white" if r[i, j] > (r.min() + 1) / 2 + 0.25 * (1 - r.min()) else "black")
                if i < j:
                    rows.append(("pearson_" + name.split()[0].lower(), "{}|{}".format(samples[i], samples[j]),
                                 "{:.4f}".format(r[i, j])))
        ax.set_xticks(range(len(samples)))
        ax.set_yticks(range(len(samples)))
        ax.set_xticklabels([label(s) for s in samples], rotation=-35, ha="left", rotation_mode="anchor")
        ax.set_yticklabels([label(s) for s in samples])
        ax.set_title("{}  {} ({:,} genes)".format(letter, name, len(m)), loc="left", fontweight="bold", fontsize=10)
        cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
        cb.set_label("Pearson r", fontsize=9)
    save(fig, os.path.join(args.outdir, "correlation"))

    # per donor: Zymo vs Genotek, panel genes
    donors = list(dict.fromkeys(donor_of[s] for s in samples))
    by = {(donor_of[s], kit_of[s]): s for s in samples}
    pairs = [(d, by[(d, kits[0])], by[(d, kits[1])]) for d in donors
             if len(kits) == 2 and (d, kits[0]) in by and (d, kits[1]) in by]
    if pairs:
        fig, axes = plt.subplots(1, len(pairs), figsize=(3.0 * len(pairs), 3.3), squeeze=False)
        fig.subplots_adjust(left=0.07, right=0.99, top=0.88, bottom=0.17, wspace=0.35)
        hi = float(panel_log.values.max()) + 0.5
        for ax, (d, a, b) in zip(axes[0], pairs):
            xa, yb = panel_log[a], panel_log[b]
            ax.scatter(xa, yb, s=6, color=SLATE, alpha=0.6, linewidth=0)
            for g in genes:
                if g in panel_log.index:
                    ax.scatter(xa[g], yb[g], s=30, color="#C00000", zorder=3)
                    ax.annotate(g, (xa[g], yb[g]), xytext=(4, -10), textcoords="offset points", fontsize=8,
                                fontstyle="italic", color="#C00000")
            ax.plot([0, hi], [0, hi], color="#999999", linewidth=0.8, linestyle="--")
            r = np.corrcoef(xa, yb)[0, 1]
            ax.set_xlim(0, hi)
            ax.set_ylim(0, hi)
            ax.set_title("Donor {}  r = {:.2f}".format(d, r), fontsize=10)
            ax.set_xlabel("{}, log2(CPM + 1)".format(kits[0]))
            ax.set_ylabel("{}, log2(CPM + 1)".format(kits[1]))
            classic(ax)
            rows.append(("kit_pearson_panel", d, "{:.4f}".format(r)))
        save(fig, os.path.join(args.outdir, "kit_scatter"))

    with open(os.path.join(args.outdir, "expression_qc.tsv"), "w") as fh:
        fh.write("measure\tsample_or_pair\tvalue\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")
    print("expression QC -> {}".format(args.outdir))


if __name__ == "__main__":
    main()
