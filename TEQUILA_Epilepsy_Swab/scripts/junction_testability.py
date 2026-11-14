#!/usr/bin/env python3
"""Can a variant be tested with RNA from an accessible tissue? Read support of the canonical splice junctions in
brain and accessible tissues (GTEx v8), and, when a variant table is given, of the junctions each patient variant
sits on.

Public part (no patient data): for each gene, every intron of the MANE Select transcript (--gtf) is looked up in the
GTEx v8 median junction read counts (public_data/raw/gtex_v8_junctions_<GENE>.tsv, from scripts/
fetch_gtex_isoforms.py). Outputs: junction_testability.tsv (gene x tissue: junctions, % with median >= 20 and >= 5
reads), junctions_<GENE>.tsv (median reads per junction and tissue), junction_testability.{pdf,png}.

Patient part (--variant-coverage: variant_coverage.tsv of scripts/variant_coverage.py, kept on scratch): each variant
is tied to the junctions that show it:
  exonic variant           the two junctions flanking its exon (the higher count: reads that include the exon)
  intronic / splice-site   the junction of the intron that holds it (reads of the splicing the variant may change)
A variant is testable in a tissue when that junction has a median of >= 20 reads (5-19: low; < 5: not testable),
the threshold used for the swab reads. The swab column is the median, over the swabs, of the reads at the variant in
this pilot (variant_coverage.tsv). Outputs: variant_testability.tsv and variant_testability.{pdf,png}. These hold
patient variants: keep them out of the repository.
GTEx medians come from short-read libraries of about 50-80 million reads; counts scale with depth.
"""
import argparse
import csv
import os
import sys
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.patches import Patch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figure_style import canonical_model  # noqa: E402

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
    "font.size": 9, "axes.labelsize": 10, "xtick.labelsize": 9, "ytick.labelsize": 9, "axes.linewidth": 0.8,
    "pdf.fonttype": 42, "savefig.facecolor": "white"})
TISSUES = [("Brain_Cortex", "Brain cortex"), ("Brain_Cerebellum", "Brain cerebellum"), ("Lung", "Lung (airway)"),
           ("Cells_Cultured_fibroblasts", "Fibroblasts"), ("Minor_Salivary_Gland", "Salivary gland"),
           ("Skin_Not_Sun_Exposed_Suprapubic", "Skin"), ("Esophagus_Mucosa", "Esophagus mucosa"),
           ("Cells_EBV-transformed_lymphocytes", "EBV lymphocytes"), ("Whole_Blood", "Whole blood")]
BINS = [0, 5, 20, 1e9]
BIN_COLORS = ["#F0F0F0", "#FDB863", "#1B7837"]
BIN_LABELS = ["< 5 reads", "5-19", ">= 20 (testable)"]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--gtf", required=True)
    p.add_argument("--genes", required=True)
    p.add_argument("--junctions-dir", required=True, help="public_data/raw")
    p.add_argument("--variant-coverage", help="variant_coverage.tsv (patient variants; optional)")
    p.add_argument("--outdir", required=True)
    return p.parse_args()


def save(fig, path):
    fig.savefig(path + ".pdf")
    fig.savefig(path + ".png", dpi=300)
    plt.close(fig)


def read_junctions(path):
    jx = defaultdict(dict)
    with open(path) as fh:
        for r in csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"):
            c, s, e = r["junction"].split("_")
            jx[(c, int(s), int(e))][r["tissue"]] = float(r["median_reads"])
    return jx


def binned(v):
    return 0 if v is None or np.isnan(v) else int(np.digitize(v, BINS[1:-1]))


def main():
    a = parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    genes = [g for g in a.genes.split(",") if os.path.exists(
        os.path.join(a.junctions_dir, "gtex_v8_junctions_{}.tsv".format(g)))]
    models, introns, med = {}, {}, {}
    summary = []
    for g in genes:
        m = canonical_model(a.gtf, g)
        models[g] = m
        ex = m["exons"]
        intr = [(m["chrom"], ex[i][1] + 1, ex[i + 1][0] - 1) for i in range(len(ex) - 1)]
        if m["strand"] == "-":
            intr = intr[::-1]          # 5' -> 3' along the transcript
        introns[g] = intr
        jx = read_junctions(os.path.join(a.junctions_dir, "gtex_v8_junctions_{}.tsv".format(g)))
        M = np.array([[jx.get(j, {}).get(t, 0.0) for t, _ in TISSUES] for j in intr]) if intr else np.zeros((0, 9))
        med[g] = M
        with open(os.path.join(a.outdir, "junctions_{}.tsv".format(g)), "w") as fh:
            fh.write("intron\tjunction\t" + "\t".join(t for t, _ in TISSUES) + "\n")
            for k, (j, row) in enumerate(zip(intr, M)):
                fh.write("{}\t{}_{}_{}\t{}\n".format(k + 1, *j, "\t".join("{:g}".format(v) for v in row)))
        for ti, (t, _) in enumerate(TISSUES):
            col = M[:, ti] if len(M) else np.array([])
            summary.append((g, t, len(col), round(float((col >= 20).mean() * 100), 1) if len(col) else "",
                            round(float((col >= 5).mean() * 100), 1) if len(col) else "",
                            round(float(np.median(col)), 1) if len(col) else ""))
    with open(os.path.join(a.outdir, "junction_testability.tsv"), "w") as fh:
        fh.write("gene\ttissue\tmane_junctions\tpercent_ge20\tpercent_ge5\tmedian_reads\n")
        for r in summary:
            fh.write("\t".join(map(str, r)) + "\n")

    # ---- public figure: % of MANE junctions with >= 20 reads; per-junction tiles for the first two genes
    P = np.array([[next(r[3] for r in summary if r[0] == g and r[1] == t) or 0 for t, _ in TISSUES]
                  for g in genes], dtype=float)
    two = genes[:2]
    fig = plt.figure(figsize=(13.0, 4.8))
    ax = fig.add_axes([0.07, 0.34, 0.3, 0.56])
    cmap = ListedColormap(["#F0F0F0", "#D9F0D3", "#A6DBA0", "#5AAE61", "#1B7837"])
    im = ax.imshow(P, aspect="auto", cmap=cmap, vmin=0, vmax=100)
    for i in range(P.shape[0]):
        for j in range(P.shape[1]):
            ax.text(j, i, "{:.0f}".format(P[i, j]), ha="center", va="center", fontsize=7.5,
                    color="white" if P[i, j] > 70 else "black")
    ax.set_xticks(range(len(TISSUES)))
    ax.set_xticklabels([l for _, l in TISSUES], rotation=45, ha="right", rotation_mode="anchor", fontsize=7.5)
    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels(genes, fontstyle="italic")
    ax.set_title("a  % of MANE junctions with >= 20 reads", loc="left", fontsize=10, fontweight="bold")
    cm2 = ListedColormap(BIN_COLORS)
    norm = BoundaryNorm([0, 1, 2, 3], 3)
    for k, g in enumerate(two):
        ax = fig.add_axes([0.45 + k * 0.28, 0.34, 0.24, 0.56])
        M = med[g]
        B = np.vectorize(binned)(M.T) if len(M) else np.zeros((len(TISSUES), 0))
        ax.imshow(B + 0.5, aspect="auto", cmap=cm2, norm=norm)
        ax.set_yticks(range(len(TISSUES)))
        ax.set_yticklabels([l for _, l in TISSUES] if k == 0 else [], fontsize=7.5)
        ax.set_xticks(range(M.shape[0]))
        ax.set_xticklabels([str(i + 1) for i in range(M.shape[0])], fontsize=6.5)
        ax.set_xlabel("Intron of the MANE transcript (5' to 3')", fontsize=8)
        ax.set_title("{}  {}".format("bc"[k], g), loc="left", fontsize=10, fontweight="bold", fontstyle="italic")
    fig.legend(handles=[Patch(color=c, label=l) for c, l in zip(BIN_COLORS, BIN_LABELS)], loc="lower right",
               bbox_to_anchor=(0.97, 0.01), ncol=3, frameon=False, fontsize=8, title="median junction reads (GTEx)",
               title_fontsize=8)
    fig.text(0.07, 0.01, "GTEx v8 median junction reads per tissue (short-read libraries, ~50-80 M reads).",
             fontsize=7.5, color="#555555")
    save(fig, os.path.join(a.outdir, "junction_testability"))
    print("junction testability -> {}".format(a.outdir))

    # ---- patient variants
    if not a.variant_coverage or not os.path.exists(a.variant_coverage):
        return
    with open(a.variant_coverage) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    swab_cols = [c for c in (rows[0] if rows else {}) if c.endswith("_reads") and not c.endswith("junction_reads")
                 and c not in ("all_swabs_reads", "all_samples_reads")]
    res = []
    for r in rows:
        g = r["gene"]
        if g not in models:
            continue
        m = models[g]
        pos = int(r["start"])
        ex = m["exons"]
        intr_g = [(m["chrom"], ex[i][1] + 1, ex[i + 1][0] - 1) for i in range(len(ex) - 1)]
        use = []
        if r["location"].startswith("intronic"):
            use = [j for j in intr_g if j[1] <= pos <= j[2]]
            if not use:   # splice-site base just outside the intron (e.g. c.N+1 written on the exon edge)
                use = [j for j in intr_g if abs(j[1] - pos) <= 3 or abs(j[2] - pos) <= 3]
        else:
            k = next((i for i, (s, e) in enumerate(ex) if s <= pos <= e), None)
            if k is not None:
                use = [intr_g[i] for i in (k - 1, k) if 0 <= i < len(intr_g)]
        jx = read_junctions(os.path.join(a.junctions_dir, "gtex_v8_junctions_{}.tsv".format(g)))
        depth = {t: max([jx.get(j, {}).get(t, 0.0) for j in use] or [np.nan]) for t, _ in TISSUES}
        swab = float(np.median([float(r[c]) for c in swab_cols])) if swab_cols else np.nan
        res.append((g, r["variant"], r["type"], r["location"], pos, swab, depth))
    with open(os.path.join(a.outdir, "variant_testability.tsv"), "w") as fh:
        fh.write("gene\tvariant\ttype\tlocation\tstart\tswab_median_reads\t" + "\t".join(t for t, _ in TISSUES) + "\n")
        for g, v, ty, loc, pos, swab, d in res:
            fh.write("{}\t{}\t{}\t{}\t{}\t{:g}\t{}\n".format(g, v, ty, loc, pos, swab,
                                                         "\t".join("{:g}".format(d[t]) for t, _ in TISSUES)))
    vg = [g for g in genes if any(x[0] == g for x in res)]
    if not vg:
        return
    cols = [("swab", "Swabs, this pilot")] + TISSUES
    fig, axes = plt.subplots(1, len(vg) + 1, figsize=(4.2 * len(vg) + 4.6, 6.4),
                             gridspec_kw={"width_ratios": [1] * len(vg) + [1.25]})
    fig.subplots_adjust(left=0.06, right=0.98, bottom=0.27, top=0.9, wspace=0.35)
    summ = []
    for ax, g in zip(axes, vg):
        vv = sorted([x for x in res if x[0] == g], key=lambda x: x[4] if models[g]["strand"] == "+" else -x[4])
        B = np.array([[binned(x[5])] + [binned(x[6][t]) if not np.isnan(x[6][t]) else 0 for t, _ in TISSUES]
                      for x in vv])
        ax.imshow(B + 0.5, aspect="auto", cmap=cm2, norm=norm, interpolation="nearest")
        ax.axvline(0.5, color="black", lw=1)
        ax.set_xticks(range(len(cols)))
        ax.set_xticklabels([l for _, l in cols], rotation=45, ha="right", rotation_mode="anchor", fontsize=7.5)
        ax.set_yticks([])
        ax.set_ylabel("{} patient variants, 5' to 3'".format(len(vv)), fontsize=8)
        ax.set_title("{}  {}".format("abc"[vg.index(g)], g), loc="left", fontsize=10, fontweight="bold",
                     fontstyle="italic")
        for j in range(len(cols)):
            summ.append((g, cols[j][0], (B[:, j] == 2).mean() * 100, len(vv)))
    ax = axes[-1]
    w = 0.8 / len(vg)
    for k, g in enumerate(vg):
        v = [s[2] for s in summ if s[0] == g]
        ax.bar(np.arange(len(cols)) - 0.4 + w * (k + 0.5), v, width=w, color=["#2166AC", "#B2182B"][k % 2],
               label=g)
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels([l for _, l in cols], rotation=45, ha="right", rotation_mode="anchor", fontsize=7.5)
    ax.set_ylim(0, 105)
    ax.set_ylabel("% of variants testable (>= 20 reads)")
    ax.set_title("{}  Testable variants".format("abc"[len(vg)]), loc="left", fontsize=10, fontweight="bold")
    leg = ax.legend(frameon=False, fontsize=8)
    for t in leg.get_texts():
        t.set_fontstyle("italic")
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    fig.legend(handles=[Patch(color=c, label=l) for c, l in zip(BIN_COLORS, BIN_LABELS)], loc="lower left",
               bbox_to_anchor=(0.06, 0.0), ncol=3, frameon=False, fontsize=8)
    fig.text(0.98, 0.01, "Swabs: median reads at the variant across the swabs of this pilot. Tissues: GTEx v8 median "
                         "reads of the junction(s) the variant sits on.", ha="right", fontsize=7.5, color="#555555")
    save(fig, os.path.join(a.outdir, "variant_testability"))
    with open(os.path.join(a.outdir, "variant_testability_summary.tsv"), "w") as fh:
        fh.write("gene\tspecimen\tpercent_testable\tvariants\n")
        for s in summ:
            fh.write("{}\t{}\t{:.2f}\t{}\n".format(*s))
    for s in summ:
        print("  {} {}: {:.0f}% of {} variants testable".format(*s))


if __name__ == "__main__":
    main()
