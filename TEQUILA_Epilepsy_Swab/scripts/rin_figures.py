#!/usr/bin/env python3
"""RNA integrity number (RIN) of the swabs against the read-level integrity of the target genes.

  rin_integrity   (a) RIN per swab, coloured by kit, with the mean RIN of buccal swabs reported by Martorella et al.
                  2023 (dashed); (b) RIN against the median transcript coverage of a gene read; (c) RIN against the share
                  of gene reads covering at least half the canonical transcript. Letters mark donors; Spearman rho is
                  exploratory (n = swabs).
Inputs: rin.tsv (sample_id, RIN), results/fragmentation/filtered/fragmentation_summary.tsv, the sample sheet and
public_data/martorella2023_summary.tsv. Every drawn number goes to rin_integrity.tsv.
"""
import argparse
import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
    "font.size": 9, "axes.labelsize": 10, "axes.linewidth": 0.8, "pdf.fonttype": 42, "savefig.facecolor": "white"})
KIT_COLORS = ["#3188BD", "#EC7B02"]
MARKERS = ["o", "s", "^", "D"]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--rin", required=True)
    p.add_argument("--fragmentation", required=True, help="fragmentation_summary.tsv (gene reads)")
    p.add_argument("--sample-sheet", required=True)
    p.add_argument("--martorella", required=True, help="public_data/martorella2023_summary.tsv")
    p.add_argument("--genes", required=True)
    p.add_argument("--outdir", required=True)
    return p.parse_args()


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"))


def classic(ax):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    genes = args.genes.split(",")
    sheet = read_tsv(args.sample_sheet)
    kits = list(dict.fromkeys(r["condition"] for r in sheet))
    meta = {r["sample_id"]: r for r in sheet}
    rin = {r["sample_id"]: float(r["RIN"]) for r in read_tsv(args.rin) if r.get("RIN")}
    samples = [r["sample_id"] for k in kits for r in sheet if r["condition"] == k and r["sample_id"] in rin]
    frag = {(r["gene"], r["sample_id"]): r for r in read_tsv(args.fragmentation)}
    buccal = [float(r["value"]) for r in read_tsv(args.martorella)
              if r["measure"] == "mean_RIN" and r["sample_type"].startswith("Buccal")]
    color = {k: KIT_COLORS[i % 2] for i, k in enumerate(kits)}
    rows = []

    fig, axes = plt.subplots(1, 3, figsize=(12.3, 3.7), gridspec_kw={"width_ratios": [1.2, 1, 1]})
    fig.subplots_adjust(left=0.06, right=0.98, top=0.86, bottom=0.25, wspace=0.35)
    ax = axes[0]
    x = np.arange(len(samples))
    ax.bar(x, [rin[s] for s in samples], width=0.75, color=[color[meta[s]["condition"]] for s in samples])
    for i, s in enumerate(samples):
        ax.text(i, rin[s] + 0.1, "{:.1f}".format(rin[s]), ha="center", va="bottom", fontsize=8)
        rows.append(("RIN", "", s, rin[s]))
    handles = [plt.Rectangle((0, 0), 1, 1, color=color[k]) for k in kits]
    names = list(kits)
    if buccal:
        ax.axhline(buccal[0], color="#555555", lw=0.9, ls="--")
        handles.append(plt.Line2D([], [], color="#555555", lw=0.9, ls="--"))
        names.append("Buccal swabs, Martorella 2023 (mean {:.1f})".format(buccal[0]))
    ax.set_ylim(0, 10)
    ax.set_ylabel("RNA integrity number (RIN)")
    ax.set_xticks(x)
    ax.set_xticklabels([s.replace("_", " ") for s in samples], rotation=-30, ha="left", rotation_mode="anchor")
    ax.legend(handles, names, frameon=False, fontsize=8, loc="upper left")
    ax.set_title("a", loc="left", fontweight="bold")
    classic(ax)
    for ax, (col, ylabel, scale), letter in zip(
            axes[1:], (("median_tx_coverage", "Median transcript coverage\nof a gene read (%)", 100),
                       ("pct_reads_cover_ge50pct", "Gene reads covering >= 50%\nof the transcript (%)", 1)), "bc"):
        for gi, g in enumerate(genes):
            pts = [(rin[s], float(frag[(g, s)][col]) * scale, s) for s in samples if (g, s) in frag]
            if not pts:
                continue
            for r_, v, s in pts:
                ax.scatter(r_, v, s=36, marker=MARKERS[gi % 4], color=color[meta[s]["condition"]],
                           edgecolor="black", linewidth=0.4, zorder=3)
                ax.annotate(meta[s]["donor"], (r_, v), xytext=(4, 3), textcoords="offset points", fontsize=7)
                rows.append((col, g, s, round(v, 3)))
            rho, p = stats.spearmanr([a for a, _, _ in pts], [b for _, b, _ in pts])
            rows.append(("spearman_rho_" + col, g, "all", "{:.3f} (p = {:.3g})".format(rho, p)))
            ax.scatter([], [], marker=MARKERS[gi % 4], color="white", edgecolor="black",
                       label="{}  rho = {:.2f}".format(g, rho))
        ax.set_xlabel("RIN")
        ax.set_ylabel(ylabel)
        leg = ax.legend(frameon=False, fontsize=8, loc="upper left")
        for t in leg.get_texts():
            t.set_fontstyle("italic")
        ax.set_title(letter, loc="left", fontweight="bold")
        classic(ax)
    fig.savefig(os.path.join(args.outdir, "rin_integrity.pdf"))
    fig.savefig(os.path.join(args.outdir, "rin_integrity.png"), dpi=300)
    plt.close(fig)
    with open(os.path.join(args.outdir, "rin_integrity.tsv"), "w") as fh:
        fh.write("measure\tgene\tsample_id\tvalue\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")
    print("RIN figures -> {}".format(args.outdir))


if __name__ == "__main__":
    main()
