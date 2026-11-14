#!/usr/bin/env python3
"""Does the cell mix of a swab explain its RNA integrity? (exploratory, n = number of swabs)

Cell-type marker levels per swab come from the composition figures (composition_figures.tsv, measure
marker_group_cpm: counts per million off-panel gene reads, mean of each marker group). Integrity per swab and gene
comes from the fragmentation summary (median transcript coverage of a read, % of reads covering >= 50%), plus RIN when
given. Each panel plots one marker group against one integrity measure, with Spearman rho.
Outputs: cellmix_integrity.{pdf,png}, cellmix_integrity.tsv.
"""
import argparse
import csv
import os
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import spearmanr

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
    "font.size": 9, "axes.labelsize": 9, "xtick.labelsize": 8, "ytick.labelsize": 8, "axes.linewidth": 0.8,
    "pdf.fonttype": 42, "savefig.facecolor": "white"})
KIT_COLORS = ["#3188BD", "#EC7B02"]
GROUPS = ["Squamous epithelium", "Salivary gland", "Immune cells"]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--composition-figures", required=True, help="composition_figures.tsv")
    p.add_argument("--fragmentation", required=True, help="fragmentation_summary.tsv")
    p.add_argument("--rin", help="TSV with sample_id, RIN (optional)")
    p.add_argument("--sample-sheet", required=True)
    p.add_argument("--genes", required=True)
    p.add_argument("--outdir", required=True)
    return p.parse_args()


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"))


def main():
    a = parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    sheet = read_tsv(a.sample_sheet)
    kits = list(dict.fromkeys(r["condition"] for r in sheet))
    kit_of = {r["sample_id"]: r["condition"] for r in sheet}
    donor_of = {r["sample_id"]: r["donor"] for r in sheet}
    mk = defaultdict(dict)
    for r in read_tsv(a.composition_figures):
        if r["measure"] == "marker_group_cpm":
            mk[r["key"]][r["sample_id"]] = float(r["value"])
    groups = [g for g in GROUPS if g in mk]
    if not groups:
        raise SystemExit("no marker_group_cpm in {} (needs gene_counts.tsv for the composition figures)".format(
            a.composition_figures))
    ys = []
    frag = read_tsv(a.fragmentation)
    for g in a.genes.split(","):
        v = {r["sample_id"]: float(r["median_tx_coverage"]) * 100 for r in frag
             if r["gene"] == g and int(float(r.get("reads_on_canonical_exons") or 20)) >= 20}
        if len(v) >= 4:
            ys.append(("{} median transcript coverage (%)".format(g), v))
    if a.rin and os.path.exists(a.rin):
        ys.append(("RIN", {r["sample_id"]: float(r["RIN"]) for r in read_tsv(a.rin) if r.get("RIN")}))
    fig, axes = plt.subplots(len(ys), len(groups), figsize=(3.3 * len(groups), 2.6 * len(ys) + 0.4), squeeze=False)
    fig.subplots_adjust(left=0.1, right=0.98, bottom=0.12, top=0.93, wspace=0.35, hspace=0.6)
    out = []
    for i, (ylab, yv) in enumerate(ys):
        for j, grp in enumerate(groups):
            ax = axes[i, j]
            ss = [s for s in yv if s in mk[grp]]
            x = np.array([mk[grp][s] for s in ss])
            y = np.array([yv[s] for s in ss])
            for s, xi, yi in zip(ss, x, y):
                ax.scatter(xi, yi, s=30, color=KIT_COLORS[kits.index(kit_of[s]) % 2], edgecolor="black", lw=0.4)
                ax.annotate(donor_of.get(s, ""), (xi, yi), xytext=(4, 0), textcoords="offset points", fontsize=7,
                            va="center")
            rho, pv = spearmanr(x, y) if len(ss) >= 4 else (np.nan, np.nan)
            ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: "{:,.0f}".format(v)))
            ax.tick_params(axis="x", labelsize=7)
            ax.set_title("rho {:.2f} (p {:.2f}, n = {})".format(rho, pv, len(ss)), fontsize=8, loc="left")
            if i == len(ys) - 1:
                ax.set_xlabel("{} markers (CPM)".format(grp))
            if j == 0:
                ax.set_ylabel(ylab, fontsize=8)
            for side in ("top", "right"):
                ax.spines[side].set_visible(False)
            out.append((ylab, grp, len(ss), round(float(rho), 3) if not np.isnan(rho) else "",
                        round(float(pv), 3) if not np.isnan(pv) else ""))
    for k, kit in enumerate(kits):
        axes[0, -1].scatter([], [], s=30, color=KIT_COLORS[k % 2], edgecolor="black", lw=0.4, label=kit)
    h, l = axes[0, -1].get_legend_handles_labels()
    fig.legend(h, l, frameon=False, fontsize=8, loc="lower right", bbox_to_anchor=(0.98, 0.0), ncol=len(kits))
    fig.text(0.1, 0.01, "Letters: donors. Exploratory: Spearman correlation on {} swabs.".format(len(sheet)),
             fontsize=7.5, color="#555555")
    fig.savefig(os.path.join(a.outdir, "cellmix_integrity.pdf"))
    fig.savefig(os.path.join(a.outdir, "cellmix_integrity.png"), dpi=300)
    plt.close(fig)
    with open(os.path.join(a.outdir, "cellmix_integrity.tsv"), "w") as fh:
        fh.write("integrity_measure\tmarker_group\tn\tspearman_rho\tp_value\n")
        for r in out:
            fh.write("\t".join(map(str, r)) + "\n")
    print("cell mix vs integrity -> {}".format(a.outdir))
    for r in out:
        print("  ", r)


if __name__ == "__main__":
    main()
