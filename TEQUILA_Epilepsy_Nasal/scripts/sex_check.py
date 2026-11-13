#!/usr/bin/env python3
"""Sample identity check from sex-chromosome gene expression (XIST against Y-linked genes).

Two swabs of the same donor must agree. The genes used (XIST; RPS4Y1, DDX3Y, KDM5D, UTY, EIF1AY, USP9Y, ZFY) are the
standard expression-based sex check. Per swab: XIST and summed Y-gene expression in CPM, and a call (female: XIST
high, Y genes near zero; male: the reverse; unclear otherwise). Output stays with the results (not the repository):
sex_check.tsv, sex_check.{pdf,png}.
Input: gene_counts.tsv from scripts/gene_counts.py.
"""
import argparse
import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
    "font.size": 9, "axes.labelsize": 10, "axes.linewidth": 0.8, "pdf.fonttype": 42, "savefig.facecolor": "white"})
Y_GENES = ["RPS4Y1", "DDX3Y", "KDM5D", "UTY", "EIF1AY", "USP9Y", "ZFY"]
KIT_COLORS = ["#3188BD", "#EC7B02"]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--counts", required=True, help="gene_counts.tsv")
    p.add_argument("--sample-sheet", required=True)
    p.add_argument("--outdir", required=True)
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    with open(args.sample_sheet) as fh:
        sheet = list(csv.DictReader(fh, delimiter="\t"))
    kits = list(dict.fromkeys(r["condition"] for r in sheet))
    tab = pd.read_csv(args.counts, sep="\t", index_col=0)
    samples = [r["sample_id"] for r in sheet if r["sample_id"] in tab.columns]
    genic = tab.drop(index=[g for g in ("intergenic", "unmapped") if g in tab.index])[samples]
    cpm = genic / genic.sum() * 1e6
    xist = cpm.loc["XIST"] if "XIST" in cpm.index else pd.Series(0.0, index=samples)
    ygenes = [g for g in Y_GENES if g in cpm.index]
    ysum = cpm.loc[ygenes].sum() if ygenes else pd.Series(0.0, index=samples)
    calls = {}
    for s in samples:
        x, y = xist[s], ysum[s]
        calls[s] = "female" if x >= 5 * max(y, 1) else ("male" if y >= 5 * max(x, 1) else "unclear")
    meta = {r["sample_id"]: r for r in sheet}
    agree = {}
    for d in dict.fromkeys(meta[s]["donor"] for s in samples):
        cs = {calls[s] for s in samples if meta[s]["donor"] == d}
        agree[d] = len(cs) == 1  # same call in every swab of the donor ("unclear" in both also agrees)
    with open(os.path.join(args.outdir, "sex_check.tsv"), "w") as fh:
        fh.write("sample_id\tdonor\tkit\tXIST_cpm\tY_genes_cpm\tcall\tdonor_swabs_agree\n")
        for s in samples:
            fh.write("{}\t{}\t{}\t{:.2f}\t{:.2f}\t{}\t{}\n".format(s, meta[s]["donor"], meta[s]["condition"], xist[s],
                                                                   ysum[s], calls[s], agree[meta[s]["donor"]]))
    fig, ax = plt.subplots(figsize=(5.0, 4.2))
    fig.subplots_adjust(left=0.16, right=0.97, top=0.95, bottom=0.14)
    for s in samples:
        k = meta[s]["condition"]
        ax.scatter(xist[s] + 0.1, ysum[s] + 0.1, s=60, color=KIT_COLORS[kits.index(k) % 2], edgecolor="black",
                   linewidth=0.5, zorder=3)
        ax.annotate(meta[s]["donor"], (xist[s] + 0.1, ysum[s] + 0.1), xytext=(5, 4), textcoords="offset points")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("XIST (CPM + 0.1)")
    ax.set_ylabel("Y-linked genes, summed (CPM + 0.1)")
    ax.legend([plt.Line2D([], [], marker="o", ls="", color=KIT_COLORS[i % 2], markeredgecolor="black")
               for i in range(len(kits))], kits, frameon=False, fontsize=8, title="letter = donor",
              title_fontsize=8)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    fig.savefig(os.path.join(args.outdir, "sex_check.pdf"))
    fig.savefig(os.path.join(args.outdir, "sex_check.png"), dpi=300)
    plt.close(fig)
    unclear = sorted({meta[s]["donor"] for s in samples if calls[s] == "unclear"})
    print("sex check: the swabs agree for {} of {} donors; unclear calls: {} -> {}".format(
        sum(agree.values()), len(agree), ", ".join(unclear) or "none", args.outdir))


if __name__ == "__main__":
    main()
