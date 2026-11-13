#!/usr/bin/env python3
"""Figures of the RNA composition analysis (scripts/rna_composition.py), in the progress-deck style.

  composition      where the reads come from, per sample: (a) all reads, (b) off-panel reads (not selected by the
                   probes, the closest view of the input RNA)
  read_length      read length per category (median and interquartile range per sample)
  genebody         coverage along canonical transcripts >= 1.5 kb (off-panel genes), 5' -> 3', per sample
  markers          cell-type marker genes, counts per million off-panel gene reads (from 11_expression_qc.sh
                   gene_counts.tsv): squamous epithelium, salivary gland, immune cells, airway epithelium, and OMP
                   (olfactory neurons). Airway and olfactory markers are the comparison for a nasal sample.
Table: composition_figures.tsv (every value drawn, plus the summary measures used in slide titles).
"""
import argparse
import csv
import os
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
    "font.size": 9, "axes.labelsize": 10, "xtick.labelsize": 9, "ytick.labelsize": 9, "axes.linewidth": 0.8,
    "xtick.major.size": 3, "ytick.major.size": 3, "pdf.fonttype": 42, "savefig.facecolor": "white"})
KIT_COLORS = ["#3188BD", "#EC7B02"]
CAT_COLORS = {"exonic, spliced": "#1B7837", "exonic, unspliced": "#A6DBA0", "intronic": "#9970AB",
              "intergenic": "#C2A5CF", "mitochondrial": "#E08214", "rRNA": "#FDB863", "unmapped": "#BDBDBD"}
CATEGORIES = list(CAT_COLORS)
HEAT = LinearSegmentedColormap.from_list("deck", ["#FFFFFF", "#F3E79B", "#FAB27B", "#EB7F86", "#CE6693", "#A059A0",
                                                  "#5C53A5"])
MARKERS = [
    ("Squamous epithelium", ["KRT13", "KRT4", "KRT6A", "SPRR3", "CRNN", "TGM3"]),
    ("Salivary gland", ["STATH", "HTN3", "PRB1", "MUC7", "SMR3B"]),
    ("Immune cells", ["PTPRC", "FCGR3B", "CSF3R", "CXCR2", "LYZ", "CD74"]),
    ("Airway epithelium", ["FOXJ1", "CAPS", "TPPP3", "SCGB1A1", "BPIFA1", "MUC5AC"]),
    ("Olfactory neurons", ["OMP"]),
]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--composition-dir", required=True, help="output of scripts/rna_composition.py")
    p.add_argument("--gene-counts", help="gene_counts.tsv of scripts/gene_counts.py (marker genes)")
    p.add_argument("--sample-sheet", required=True)
    p.add_argument("--outdir", required=True)
    return p.parse_args()


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


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
    sheet = read_tsv(args.sample_sheet)
    cd = args.composition_dir
    comp = defaultdict(dict)
    for r in read_tsv(os.path.join(cd, "composition.tsv")):
        comp[(r["sample_id"], r["scope"])][r["category"]] = (int(r["reads"]), float(r["percent"]))
    present = {s for s, _ in comp}
    sheet = [r for r in sheet if r["sample_id"] in present]
    kits = list(dict.fromkeys(r["condition"] for r in sheet))
    samples = [r["sample_id"] for k in kits for r in sheet if r["condition"] == k]
    kit_of = {r["sample_id"]: r["condition"] for r in sheet}
    kcol = {k: KIT_COLORS[i % 2] for i, k in enumerate(kits)}
    rows = []

    # ---- composition: stacked bars, all reads and off-panel reads
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.1), sharey=True)
    fig.subplots_adjust(left=0.07, right=0.80, bottom=0.22, top=0.86, wspace=0.12)
    for ax, scope, title in zip(axes, ("all", "off_panel"), ("a  All reads", "b  Off-panel reads")):
        bottom = np.zeros(len(samples))
        for c in CATEGORIES:
            v = np.array([comp[(s, scope)].get(c, (0, 0))[1] for s in samples])
            ax.bar(range(len(samples)), v, bottom=bottom, color=CAT_COLORS[c], edgecolor="white", lw=0.4, width=0.8,
                   label=c)
            bottom += v
            rows += [("percent_" + scope, c, s, "{:.2f}".format(x)) for s, x in zip(samples, v)]
        ax.set_xticks(range(len(samples)))
        ax.set_xticklabels([label(s) for s in samples], rotation=30, ha="right", rotation_mode="anchor")
        ax.set_ylim(0, 100)
        ax.set_title(title, loc="left", fontsize=10, fontweight="bold", pad=16)
        classic(ax)
        n = [sum(x[0] for x in comp[(s, scope)].values()) for s in samples]
        for i, x in enumerate(n):
            ax.text(i, 101, "{:.2g}M".format(x / 1e6) if x >= 1e5 else "{:,}".format(x), ha="center", va="bottom",
                    fontsize=6.5, color="#555555")
    axes[0].set_ylabel("% of reads")
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h[::-1], l[::-1], loc="center left", bbox_to_anchor=(0.81, 0.55), frameon=False, fontsize=8.5)
    save(fig, os.path.join(args.outdir, "composition"))

    # ---- read length per category
    lens = defaultdict(dict)
    for r in read_tsv(os.path.join(cd, "read_lengths.tsv")):
        if r["median"] != "NA" and int(r["reads"]) >= 20:
            lens[r["category"]][r["sample_id"]] = (int(r["q25"]), int(r["median"]), int(r["q75"]))
    cats = [c for c in CATEGORIES if lens.get(c)]
    fig, ax = plt.subplots(figsize=(8.6, 3.6))
    fig.subplots_adjust(left=0.09, right=0.98, bottom=0.2, top=0.86)
    w = 0.8 / max(len(samples), 1)
    for ci, c in enumerate(cats):
        for si, s in enumerate(samples):
            if s not in lens[c]:
                continue
            q1, m, q3 = lens[c][s]
            x = ci - 0.4 + w * (si + 0.5)
            ax.plot([x, x], [q1, q3], color=kcol[kit_of[s]], lw=1.2)
            ax.scatter(x, m, s=14, color=kcol[kit_of[s]], edgecolor="black", lw=0.3, zorder=3)
            rows.append(("read_length_median", c, s, m))
    ax.set_xticks(range(len(cats)))
    ax.set_xticklabels(cats, rotation=20, ha="right", rotation_mode="anchor")
    ax.set_ylabel("Read length (nt)")
    ax.set_ylim(bottom=0)
    for k in kits:
        ax.scatter([], [], s=14, color=kcol[k], edgecolor="black", lw=0.3, label=k)
    ax.legend(frameon=False, ncol=len(kits), loc="lower center", bbox_to_anchor=(0.5, 1.0), fontsize=8.5)
    ax.text(1.0, 1.02, "points: median; lines: interquartile range; one per sample", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=7.5, color="#555555")
    classic(ax)
    save(fig, os.path.join(args.outdir, "read_length"))

    # ---- coverage along transcripts (off-panel genes)
    cov = defaultdict(lambda: np.zeros(100))
    for r in read_tsv(os.path.join(cd, "genebody_coverage.tsv")):
        if r["scope"] == "off_panel":
            cov[r["sample_id"]][int(r["bin"])] = float(r["reads"])
    span = {r["sample_id"]: r for r in read_tsv(os.path.join(cd, "transcript_span.tsv")) if r["scope"] == "off_panel"}
    fig, ax = plt.subplots(figsize=(6.2, 3.8))
    fig.subplots_adjust(left=0.12, right=0.75, bottom=0.15, top=0.92)
    styles = ["-", "--", ":", "-."]
    donors = list(dict.fromkeys(r["donor"] for r in sheet))
    donor_of = {r["sample_id"]: r["donor"] for r in sheet}
    for s in samples:
        c = cov[s]
        if c.max() <= 0:
            continue
        ax.plot(np.arange(100) + 0.5, c / c.max(), color=kcol[kit_of[s]], lw=1.2,
                ls=styles[donors.index(donor_of[s]) % 4], label=label(s))
        rows += [("genebody_normalised", str(i), s, "{:.4f}".format(v)) for i, v in enumerate(c / c.max())]
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 1.05)
    ax.set_xlabel("Position along transcript (% of length, 5' to 3')")
    ax.set_ylabel("Read coverage (scaled to max)")
    ax.legend(frameon=False, fontsize=8, loc="center left", bbox_to_anchor=(1.01, 0.5))
    classic(ax)
    save(fig, os.path.join(args.outdir, "genebody"))
    for s, r in span.items():
        rows += [("tx_reads_off_panel", "", s, r["reads"]), ("median_tx_fraction_off_panel", "", s,
                                                             r["median_tx_fraction"]),
                 ("percent_reads_ge50_off_panel", "", s, r["percent_reads_ge50"]),
                 ("ratio_3p_5p_off_panel", "", s, r["ratio_3p_5p"])]

    # ---- cell-type markers
    if args.gene_counts and os.path.exists(args.gene_counts):
        tab = {}
        with open(args.gene_counts) as fh:
            rd = csv.DictReader(fh, delimiter="\t")
            for r in rd:
                tab[r["gene"]] = r
        denom = {s: sum(int(r[s]) for g, r in tab.items() if g not in ("intergenic", "unmapped")
                        and r["in_panel"] != "True" and s in r) for s in samples}
        groups, genes_shown, missing = [], [], []
        for grp, gs in MARKERS:
            keep = [g for g in gs if g in tab and tab[g]["in_panel"] != "True"]
            missing += [g for g in gs if g not in keep]
            if keep:
                groups.append((grp, len(genes_shown), len(genes_shown) + len(keep)))
                genes_shown += keep
        m = np.array([[int(tab[g][s]) / denom[s] * 1e6 if denom[s] else 0 for s in samples] for g in genes_shown])
        height = 0.26 * len(genes_shown) + 1.3
        fig, ax = plt.subplots(figsize=(7.6, height))
        fig.subplots_adjust(left=0.32, right=0.86, bottom=1.0 / height, top=1 - 0.15 / height)
        vals = np.log10(m + 1)
        im = ax.imshow(vals, aspect="auto", cmap=HEAT, vmin=0, vmax=max(vals.max(), 1))
        for i in range(len(genes_shown)):
            for j in range(len(samples)):
                v = m[i, j]
                ax.text(j, i, "{:.0f}".format(v) if v >= 1 else ("<1" if v > 0 else "0"), ha="center",
                        va="center", fontsize=6.5, color="white" if vals[i, j] > 0.6 * vals.max() else "black")
                rows.append(("marker_cpm", genes_shown[i], samples[j], "{:.2f}".format(v)))
        ax.set_yticks(range(len(genes_shown)))
        ax.set_yticklabels(genes_shown, fontstyle="italic", fontsize=8)
        ax.set_xticks(range(len(samples)))
        ax.set_xticklabels([label(s) for s in samples], rotation=30, ha="right", rotation_mode="anchor")
        for grp, a, b in groups:
            ax.axhline(b - 0.5, color="black", lw=0.8)
            ax.annotate(grp, xy=(-0.27, 1 - (a + b) / 2 / len(genes_shown)), xycoords="axes fraction", ha="right",
                        va="center", fontsize=8.5, fontweight="bold")
            rows.append(("marker_group_mean_cpm", grp, "all", "{:.2f}".format(m[a:b].mean())))
            for j, s in enumerate(samples):
                rows.append(("marker_group_cpm", grp, s, "{:.2f}".format(m[a:b, j].mean())))
        cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
        cb.set_label("log10(CPM + 1)", fontsize=8)
        save(fig, os.path.join(args.outdir, "markers"))
        if missing:
            print("markers left out (not in the counts, or panel genes): " + ", ".join(missing))

    # summary measures (off-panel reads, per sample) for slide titles
    for s in samples:
        p = {c: comp[(s, "off_panel")].get(c, (0, 0))[1] for c in CATEGORIES}
        rows += [("summary_exonic_spliced", "", s, "{:.2f}".format(p["exonic, spliced"])),
                 ("summary_intronic_intergenic", "", s, "{:.2f}".format(p["intronic"] + p["intergenic"])),
                 ("summary_mito_rrna", "", s, "{:.2f}".format(p["mitochondrial"] + p["rRNA"])),
                 ("summary_unmapped", "", s, "{:.2f}".format(p["unmapped"]))]
    with open(os.path.join(args.outdir, "composition_figures.tsv"), "w") as fh:
        fh.write("measure\tkey\tsample_id\tvalue\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")
    print("composition figures -> {}".format(args.outdir))


if __name__ == "__main__":
    main()
