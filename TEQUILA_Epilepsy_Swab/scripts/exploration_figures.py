#!/usr/bin/env python3
"""Figures for the exploration part of the deck: where the project has been and where it is going.

  roadmap          pilot 1 (buccal swabs) -> accessible tissues considered (status of each) -> pilot 2 (nasal)
  panel_by_tissue  GTEx v8: (a) share of the probe-panel genes with median TPM >= 10 and (b) median TPM of the target
                   genes, in brain and in the tissue closest to each accessible specimen (public_data/
                   gtex_panel_expression.tsv)
Table: exploration_figures.tsv (values drawn).
"""
import argparse
import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Patch

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
    "font.size": 9, "axes.labelsize": 10, "xtick.labelsize": 9, "ytick.labelsize": 9, "axes.linewidth": 0.8,
    "pdf.fonttype": 42, "savefig.facecolor": "white"})
STATUS = {"tested": ("Tested (pilot 1)", "#708090"), "dropped": ("Considered, dropped", "#D6604D"),
          "later": ("Kept for later", "#F4A261"), "chosen": ("Chosen (pilot 2)", "#1B7837"),
          "reference": ("Reference", "#2B2B2B")}
# GTEx tissue -> (label, specimen it stands in for, status)
TISSUES = [
    ("Brain - Cortex", "Brain cortex", "target tissue", "reference"),
    ("Lung", "Lung", "airway; closest to nasal", "chosen"),
    ("Esophagus - Mucosa", "Esophagus mucosa", "squamous; closest to buccal", "tested"),
    ("Minor Salivary Gland", "Minor salivary gland", "closest to saliva", "dropped"),
    ("Skin - Not Sun Exposed (Suprapubic)", "Skin", "closest to hair follicle", "dropped"),
    ("Cells - Cultured fibroblasts", "Cultured fibroblasts", "from a skin biopsy", "later"),
    ("Cells - EBV-transformed lymphocytes", "EBV lymphocytes", "cell line from blood", "later"),
    ("Whole Blood", "Whole blood", "blood draw", "later"),
]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--gtex", required=True, help="public_data/gtex_panel_expression.tsv")
    p.add_argument("--genes", required=True)
    p.add_argument("--outdir", required=True)
    return p.parse_args()


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"))


def save(fig, path):
    fig.savefig(path + ".pdf")
    fig.savefig(path + ".png", dpi=300)
    plt.close(fig)


def box(ax, x, y, w, h, text, color, title=None, fs=8.5, tc="black"):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.012,rounding_size=0.02", facecolor=color,
                                edgecolor="none", alpha=0.95))
    if title:
        ax.text(x + w / 2, y + h - 0.035, title, ha="center", va="top", fontsize=fs + 1, fontweight="bold", color=tc)
        ax.text(x + w / 2, y + h - 0.105, text, ha="center", va="top", fontsize=fs, color=tc, linespacing=1.3)
    else:
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs, color=tc, linespacing=1.3)


def arrow(ax, a, b, color="#555555"):
    ax.add_patch(FancyArrowPatch(a, b, arrowstyle="-|>", mutation_scale=14, lw=1.4, color=color))


def roadmap(path):
    fig = plt.figure(figsize=(12.6, 5.6))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    # row 1: pilot 1
    box(ax, 0.02, 0.70, 0.27, 0.25, "Two kits (Zymo, Genotek),\n4 lab donors, 634-gene panel,\nlong reads",
        "#DCE6F2", title="Pilot 1: buccal swabs")
    box(ax, 0.36, 0.70, 0.27, 0.25, "Both genes detected in every swab\nRNA degraded; almost no full-length\n"
        "reads; off-panel reads mostly\nintronic / intergenic", "#DCE6F2", title="What we learned")
    box(ax, 0.70, 0.70, 0.28, 0.25, "Which accessible tissue gives\nintact RNA that expresses the\n"
        "genes, and can be collected\nfrom young children?", "#FFF4D6", title="New question")
    arrow(ax, (0.295, 0.825), (0.355, 0.825))
    arrow(ax, (0.635, 0.825), (0.695, 0.825))
    # row 2: options
    opts = [("Hair follicles", "RIN 8.9 published,\nbut 25-30 plucked hairs;\nnot for under-3s", "dropped"),
            ("Saliva", "RIN 2.6; 43% of\nlibraries passed QC", "dropped"),
            ("Blood", "30% of panel genes\n>= 10 TPM; STXBP1\nnear absent", "later"),
            ("Skin biopsy", "Fibroblasts / iNeurons;\ninvasive, weeks of\nculture", "later"),
            ("Urine cells", "Cell pellets: 70% of\nOMIM genes captured;\nresearch backup", "later"),
            ("Nasal epithelium", "Brush RIN 8.9; both\ngenes in every public\nnasal brushing", "chosen")]
    w, gap, x = 0.135, 0.028, 0.02
    for name, txt, st in opts:
        box(ax, x, 0.33, w, 0.26, txt, STATUS[st][1], title=name, fs=8, tc="white")
        x += w + gap
    arrow(ax, (0.9, 0.695), (0.9, 0.60))
    ax.text(0.5, 0.62, "Accessible tissues considered", ha="center", fontsize=10, fontweight="bold", color="#333333")
    # row 3: pilot 2
    box(ax, 0.40, 0.03, 0.30, 0.21, "Same donors; nasal flocked swab vs\ncytology brush at the inferior\n"
        "turbinate; rules set before the data", "#E5F4E4", title="Pilot 2 (planned): nasal")
    box(ax, 0.75, 0.03, 0.23, 0.21, "Flocked swab at every age;\nbrush from 5 years, or during\nan anaesthetic",
        "#F0EFEC", title="Then: children")
    arrow(ax, (0.9, 0.325), (0.62, 0.245), color=STATUS["chosen"][1])
    arrow(ax, (0.705, 0.135), (0.745, 0.135))
    ax.legend(handles=[Patch(color=STATUS[k][1], label=STATUS[k][0]) for k in ("tested", "dropped", "later",
                                                                               "chosen")],
              loc="lower left", bbox_to_anchor=(0.01, 0.03), frameon=False, fontsize=8.5, ncol=1)
    save(fig, path)


def panel_by_tissue(args, rows_out):
    tab = {r["tissue"]: r for r in read_tsv(args.gtex)}
    genes = args.genes.split(",")
    items = [(lab, proxy, st, tab[t]) for t, lab, proxy, st in TISSUES if t in tab]
    fig, axes = plt.subplots(1, 2, figsize=(11.6, 3.9), gridspec_kw={"width_ratios": [1, 1.25]})
    fig.subplots_adjust(left=0.2, right=0.98, bottom=0.2, top=0.88, wspace=0.08)
    y = list(range(len(items)))[::-1]
    ax = axes[0]
    for yi, (lab, proxy, st, r) in zip(y, items):
        v = float(r["percent_ge10_tpm"])
        ax.barh(yi, v, color=STATUS[st][1], height=0.68)
        ax.text(v + 1, yi, "{:.0f}%".format(v), va="center", fontsize=8.5)
        rows_out.append(("panel_percent_ge10_tpm", lab, v))
    ax.set_yticks(y)
    ax.set_yticklabels(["{}\n({})".format(lab, proxy) for lab, proxy, _, _ in items], fontsize=8.5)
    ax.set_xlim(0, 80)
    ax.set_xlabel("% of panel genes with median TPM >= 10")
    ax.set_title("a  Probe-panel genes expressed", loc="left", fontsize=10, fontweight="bold")
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax = axes[1]
    h = 0.8 / len(genes)
    for k, g in enumerate(genes):
        for yi, (lab, proxy, st, r) in zip(y, items):
            v = float(r["{}_tpm".format(g)])
            yy = yi + 0.4 - h * (k + 0.5)
            ax.barh(yy, v, height=h * 0.9, color=STATUS[st][1], alpha=1 if k == 0 else 0.55)
            ax.text(v * 1.08 + 0.05, yy, "{:.1f}".format(v) if v < 10 else "{:.0f}".format(v), va="center",
                    fontsize=7.5)
            rows_out.append(("{}_tpm".format(g), lab, v))
    ax.set_yticks(y)
    ax.set_yticklabels([])
    ax.set_xscale("symlog", linthresh=1)
    ax.set_xlim(0, 600)
    ax.set_xlabel("Median TPM (symmetric log scale)")
    ax.set_title("b  {} (dark) and {} (light)".format(*genes) if len(genes) == 2 else "b  Target genes",
                 loc="left", fontsize=10, fontweight="bold")
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    fig.legend(handles=[Patch(color=STATUS[k][1], label=STATUS[k][0]) for k in STATUS], loc="lower center",
               ncol=5, frameon=False, fontsize=8, bbox_to_anchor=(0.58, 0.0))
    save(fig, os.path.join(args.outdir, "panel_by_tissue"))


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    roadmap(os.path.join(args.outdir, "roadmap"))
    rows = []
    panel_by_tissue(args, rows)
    with open(os.path.join(args.outdir, "exploration_figures.tsv"), "w") as fh:
        fh.write("measure\ttissue\tvalue\n")
        for r in rows:
            fh.write("{}\t{}\t{}\n".format(*r))
    print("exploration figures -> {}".format(args.outdir))


if __name__ == "__main__":
    main()
