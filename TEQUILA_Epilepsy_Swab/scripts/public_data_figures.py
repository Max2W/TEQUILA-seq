#!/usr/bin/env python3
"""Figures from public data and literature, for context on the swab pilot (drawn in the progress-deck style).

Inputs are committed under public_data/ (no internet needed):
  raw/gtex_v10_samples_<GENE>.json   GTEx Portal API v2, expression/geneExpression, dataset gtex_v10 (per-sample TPM)
  raw/gtex_v10_median_<GENE>.json    GTEx Portal API v2, expression/medianGeneExpression, dataset gtex_v10
  raw/hpa_tissue_<GENE>.json         Human Protein Atlas search_download API, tissue RNA nTPM
  martorella2023_summary.tsv         numbers transcribed from Martorella et al., BMC Genomics 2023 (full text in raw/)

Figures (<outdir>/):
  public_martorella2023   a mean RIN by noninvasive sample type (with this pilot's swabs when --our-rin is given),
                          b samples failing pre-sequencing QC with the Illumina kit, c OMIM genes captured (median %)
  public_gtex_<GENE>      expression (TPM) across GTEx v10 tissues, boxes sorted by median; accessible tissues in colour
  public_accessible       median TPM of each gene in brain and in accessible tissues (GTEx v10)
Tables: public_expression.tsv (GTEx medians and HPA nTPM shown on the slides).
"""
import argparse
import csv
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
    "font.size": 9, "axes.labelsize": 10, "xtick.labelsize": 9, "ytick.labelsize": 9, "axes.linewidth": 0.8,
    "xtick.major.size": 3, "ytick.major.size": 3, "pdf.fonttype": 42, "savefig.facecolor": "white"})
SLATE = "#708090"
ACCESSIBLE = {  # GTEx tissue id -> label, colour
    "Minor_Salivary_Gland": ("Minor Salivary Gland", "#F6CF2C"),
    "Esophagus_Mucosa": ("Esophagus - Mucosa", "#F4A261"),
    "Skin_Sun_Exposed_Lower_leg": ("Skin - Sun Exposed", "#9CCB3B"),
    "Skin_Not_Sun_Exposed_Suprapubic": ("Skin - Not Sun Exposed", "#B8DDA8"),
    "Cells_Cultured_fibroblasts": ("Cultured fibroblasts", "#57B45B"),
    "Cells_EBV-transformed_lymphocytes": ("EBV lymphocytes", "#E8957B"),
    "Whole_Blood": ("Whole Blood", "#C0392B"),
}
BAR_TISSUES = [("Brain_Cerebellum", "Brain - Cerebellum"), ("Brain_Cortex", "Brain - Cortex"),
               ("Minor_Salivary_Gland", "Minor Salivary Gland"), ("Esophagus_Mucosa", "Esophagus - Mucosa"),
               ("Skin_Sun_Exposed_Lower_leg", "Skin - Sun Exposed"), ("Cells_Cultured_fibroblasts", "Fibroblasts"),
               ("Cells_EBV-transformed_lymphocytes", "EBV lymphocytes"), ("Whole_Blood", "Whole Blood")]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--public-dir", required=True, help="public_data/ of this analysis folder")
    p.add_argument("--genes", required=True, help="comma-separated")
    p.add_argument("--our-rin", help="optional TSV with columns sample_id, RIN (this pilot's swab RNA)")
    p.add_argument("--outdir", required=True)
    return p.parse_args()


def classic(ax):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def nice(tissue_id):
    short = {"Cells_EBV-transformed_lymphocytes": "EBV lymphocytes", "Cells_Cultured_fibroblasts": "Cultured fibroblasts"}
    return short.get(tissue_id) or tissue_id.replace("_", " ").replace("Brain ", "Brain - ")


def save(fig, path):
    fig.savefig(path + ".pdf")
    fig.savefig(path + ".png", dpi=300)
    plt.close(fig)


def martorella(args):
    rows = list(csv.DictReader((l for l in open(os.path.join(args.public_dir, "martorella2023_summary.tsv"))
                                if not l.startswith("#")), delimiter="\t"))
    rin = [(r["sample_type"], float(r["value"])) for r in rows if r["measure"] == "mean_RIN"]
    qc = [(r["sample_type"], int(r["n_fail"]), int(r["n_total"])) for r in rows if r["measure"] == "illumina_preseq_qc"]
    omim = [(r["sample_type"], float(r["value"])) for r in rows if r["measure"] == "omim_capture_percent"]
    ours = []
    if args.our_rin and os.path.exists(args.our_rin):
        ours = [float(r["RIN"]) for r in csv.DictReader(open(args.our_rin), delimiter="\t") if r.get("RIN")]

    fig, axes = plt.subplots(1, 3, figsize=(12.3, 3.6), gridspec_kw={"width_ratios": [1.2, 1, 1]})
    fig.subplots_adjust(left=0.06, right=0.99, top=0.88, bottom=0.3, wspace=0.42)
    ax = axes[0]
    names = [n for n, _ in rin] + (["This pilot (buccal swabs)"] if ours else [])
    vals = [v for _, v in rin] + ([float(np.mean(ours))] if ours else [])
    cols = ["#BEBEBE" if "Hair" in n else SLATE for n in names[:len(rin)]] + (["#3A8DC4"] if ours else [])
    ax.bar(range(len(vals)), vals, width=0.7, color=cols)
    if ours:
        ax.scatter(np.full(len(ours), len(vals) - 1) + np.random.default_rng(1).uniform(-0.18, 0.18, len(ours)),
                   ours, s=12, color="black", zorder=3)
    for i, v in enumerate(vals):
        ax.text(i, v + 0.2, "{:.1f}".format(v), ha="center", va="bottom", fontsize=8)
    ax.axhline(7, color="#C00000", linewidth=0.8, linestyle="--")
    ax.text(len(vals) - 0.5, 7.15, "RIN 7", color="#C00000", fontsize=8, ha="right", va="bottom")
    ax.set_ylim(0, 10.5)
    ax.set_ylabel("RNA integrity number (mean)")
    ax.set_xticks(range(len(vals)))
    ax.set_xticklabels(names, rotation=-30, ha="left", rotation_mode="anchor")
    ax.set_title("a", loc="left", fontweight="bold")
    classic(ax)
    ax = axes[1]
    x = np.arange(len(qc))
    fails = [f / n * 100 for _, f, n in qc]
    ax.bar(x, fails, width=0.7, color=SLATE, label="failed")
    ax.bar(x, [100 - f for f in fails], bottom=fails, width=0.7, color="#D9D9D9", label="passed")
    for i, (_, f, n) in enumerate(qc):
        ax.text(i, fails[i] / 2, "{}/{}".format(f, n), ha="center", va="center", fontsize=8, color="white")
    ax.set_ylim(0, 100)
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_ylabel("Samples failing pre-sequencing\nQC, Illumina kit (%)")
    ax.set_xticks(x)
    ax.set_xticklabels([n for n, _, _ in qc], rotation=-30, ha="left", rotation_mode="anchor")
    ax.legend(frameon=False, fontsize=8, loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2)
    ax.set_title("b", loc="left", fontweight="bold", x=-0.25)
    classic(ax)
    ax = axes[2]
    ax.bar(range(len(omim)), [v for _, v in omim], width=0.7,
           color=[SLATE if "Buccal" in n else "#BEBEBE" for n, _ in omim])
    for i, (_, v) in enumerate(omim):
        ax.text(i, v + 1.5, "{:.0f}%".format(v), ha="center", va="bottom", fontsize=8)
    ax.set_ylim(0, 100)
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_ylabel("OMIM genes captured\n(median %, > 0.1 TPM)")
    ax.set_xticks(range(len(omim)))
    ax.set_xticklabels([n for n, _ in omim], rotation=-30, ha="left", rotation_mode="anchor")
    ax.set_title("c", loc="left", fontweight="bold")
    classic(ax)
    fig.text(0.06, 0.02, "Redrawn from values reported in Martorella et al., BMC Genomics 2023 (19 donors; "
             "short-read RNA-seq)" + ("; this pilot: per-sample RIN, mean bar" if ours else ""),
             fontsize=7.5, color="#444444")
    save(fig, os.path.join(args.outdir, "public_martorella2023"))
    return {"rin": dict(rin), "qc": {n: (f, t) for n, f, t in qc}, "omim": dict(omim), "ours": ours}


def gtex_box(args, gene):
    d = json.load(open(os.path.join(args.public_dir, "raw", "gtex_v10_samples_{}.json".format(gene))))["data"]
    tissues = sorted(d, key=lambda t: -np.median(t["data"]))
    fig, ax = plt.subplots(figsize=(12.3, 4.6))
    fig.subplots_adjust(left=0.06, right=0.955, top=0.95, bottom=0.42)
    data = [t["data"] for t in tissues]
    bp = ax.boxplot(data, widths=0.7, showfliers=False, patch_artist=True,
                    medianprops=dict(color="black", linewidth=1.0), whiskerprops=dict(linewidth=0.6),
                    capprops=dict(linewidth=0.6), boxprops=dict(linewidth=0.6))
    for patch, t in zip(bp["boxes"], tissues):
        patch.set_facecolor(ACCESSIBLE.get(t["tissueSiteDetailId"], (None, "white"))[1])
    ax.set_xticks(range(1, len(tissues) + 1))
    ax.set_xticklabels([nice(t["tissueSiteDetailId"]) for t in tissues], rotation=-35, ha="left",
                       rotation_mode="anchor", fontsize=7)
    for lbl, t in zip(ax.get_xticklabels(), tissues):
        if t["tissueSiteDetailId"] in ACCESSIBLE:
            lbl.set_fontweight("bold")
            lbl.set_color(ACCESSIBLE[t["tissueSiteDetailId"]][1] if t["tissueSiteDetailId"] != "Minor_Salivary_Gland"
                          else "#B8860B")
    # accessible tissues can sit far below brain (STXBP1): a full-height band in their colour keeps them visible
    for i, t in enumerate(tissues, 1):
        if t["tissueSiteDetailId"] in ACCESSIBLE:
            ax.axvspan(i - 0.5, i + 0.5, color=ACCESSIBLE[t["tissueSiteDetailId"]][1], alpha=0.25, lw=0, zorder=0)
    ax.set_ylabel("{} expression (TPM)".format(gene))
    ax.set_xlim(0.3, len(tissues) + 0.7)
    classic(ax)
    handles = [Patch(facecolor=c, edgecolor="black", linewidth=0.6, label=l) for l, c in ACCESSIBLE.values()]
    ax.legend(handles=handles, frameon=False, fontsize=8, loc="upper center", bbox_to_anchor=(0.42, 1.0),
              title="Accessible tissues", title_fontsize=8, ncol=2)
    save(fig, os.path.join(args.outdir, "public_gtex_{}".format(gene)))
    return {t["tissueSiteDetailId"]: (float(np.median(t["data"])), len(t["data"])) for t in tissues}


def accessible_bar(args, genes, medians):
    fig, ax = plt.subplots(figsize=(7.5, 3.6))
    fig.subplots_adjust(left=0.11, right=0.98, top=0.93, bottom=0.33)
    x = np.arange(len(BAR_TISSUES))
    w = 0.8 / len(genes)
    colors = [SLATE, "#E2C485", "#3A8DC4"]
    for k, g in enumerate(genes):
        v = [medians[g][t][0] for t, _ in BAR_TISSUES]
        ax.bar(x - 0.4 + w * (k + 0.5), v, width=w, color=colors[k % 3], label=g)
        for i, val in enumerate(v):
            ax.text(x[i] - 0.4 + w * (k + 0.5), val * 1.08 + 0.05, "{:.1f}".format(val), ha="center",
                    va="bottom", fontsize=6.5, rotation=90)
    ax.set_yscale("log")
    ax.set_ylim(0.1, 2000)
    ax.set_ylabel("Median expression (TPM)")
    ax.set_xticks(x)
    ax.set_xticklabels([l for _, l in BAR_TISSUES], rotation=-30, ha="left", rotation_mode="anchor")
    leg = ax.legend(frameon=False, fontsize=9, loc="upper right")
    for t in leg.get_texts():
        t.set_fontstyle("italic")
    classic(ax)
    fig.text(0.11, 0.01, "GTEx v10, median TPM per tissue", fontsize=7.5, color="#444444")
    save(fig, os.path.join(args.outdir, "public_accessible"))


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    genes = args.genes.split(",")
    martorella(args)
    medians = {g: gtex_box(args, g) for g in genes}
    accessible_bar(args, genes, medians)
    with open(os.path.join(args.outdir, "public_expression.tsv"), "w") as fh:
        fh.write("source\ttissue\tgene\tvalue\tunit\tn_samples\n")
        for g in genes:
            for t, (m, n) in sorted(medians[g].items()):
                fh.write("GTEx v10\t{}\t{}\t{:.3f}\tmedian TPM\t{}\n".format(nice(t), g, m, n))
            path = os.path.join(args.public_dir, "raw", "hpa_tissue_{}.json".format(g))
            if os.path.exists(path):
                for rec in json.load(open(path)):
                    if rec.get("Gene") != g:
                        continue
                    for k, v in rec.items():
                        if k.startswith("Tissue RNA - "):
                            fh.write("Human Protein Atlas\t{}\t{}\t{}\tnTPM\t\n".format(
                                k[len("Tissue RNA - "):].replace(" [nTPM]", ""), g, v))
    print("public data figures -> {}".format(args.outdir))


if __name__ == "__main__":
    main()
