#!/usr/bin/env python3
"""Whole-panel view: how the probe-panel genes are expressed in brain and in accessible specimens (public data,
public_data/panel_expression_atlas.tsv from scripts/panel_atlas.py).

  atlas_coverage     (a) % of panel genes with median TPM >= 1 and >= 10 per specimen; (b) of the panel genes expressed
                     in brain cortex (>= 10 TPM), % that also reach 10 TPM in each specimen
  atlas_heatmap      every panel gene (rows) x specimen (columns), log10(TPM + 1); genes grouped by how many accessible
                     specimens express them; highlighted genes labelled
  atlas_vs_brain     each accessible specimen against brain cortex, one point per panel gene; Spearman rho
  atlas_reach        for brain-expressed panel genes: in how many accessible specimens each reaches 10 TPM, and how
                     many are covered by the best one and two specimens together
  targets_specimens  the highlighted genes x specimens, median TPM printed in each cell
Table: panel_atlas_figures.tsv (values drawn and summary numbers for the slides).
"""
import argparse
import csv
import itertools
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import spearmanr

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
    "font.size": 9, "axes.labelsize": 10, "xtick.labelsize": 9, "ytick.labelsize": 9, "axes.linewidth": 0.8,
    "pdf.fonttype": 42, "savefig.facecolor": "white"})
HEAT = LinearSegmentedColormap.from_list("deck", ["#FFFFFF", "#F3E79B", "#FAB27B", "#EB7F86", "#CE6693", "#A059A0",
                                                  "#5C53A5"])
# column -> (short label, what it is, colour); measured specimens first, then GTEx tissues standing in for others
SPEC = {
    "Brain cortex": ("Brain cortex", "GTEx", "#2B2B2B"),
    "Brain cerebellum": ("Brain cerebellum", "GTEx", "#555555"),
    "Nasal brushing": ("Nasal brushing", "GSE152004", "#1B7837"),
    "Nasal brushing, children": ("Nasal brushing, children", "GSE118761", "#5AAE61"),
    "Buccal brushing": ("Buccal brushing", "GSE330022", "#708090"),
    "Plucked hair follicles": ("Plucked hair", "GSE224010", "#D6604D"),
    "Whole blood (PAXgene)": ("Whole blood", "GSE112087", "#C0392B"),
    "Urine-derived cells (cultured)": ("Urine cells, cultured", "GSE235813", "#F4A261"),
    "Cultured fibroblasts": ("Fibroblasts", "GTEx; skin biopsy", "#F6B26B"),
    "EBV lymphocytes": ("EBV lymphocytes", "GTEx; from blood", "#F6C28B"),
    "Minor salivary gland": ("Salivary gland", "GTEx; for saliva", "#E8957B"),
    "Lung": ("Lung", "GTEx; airway", "#A6DBA0"),
    "Esophagus mucosa": ("Esophagus mucosa", "GTEx; squamous", "#B0BEC5"),
    "Skin": ("Skin", "GTEx", "#F2C4B6"),
    "Whole blood": ("Whole blood, GTEx", "GTEx", "#E6A0A0"),
}
MEASURED = ["Nasal brushing", "Nasal brushing, children", "Buccal brushing", "Plucked hair follicles",
            "Whole blood (PAXgene)", "Urine-derived cells (cultured)"]
# the accessible specimens compared: measured ones, plus GTEx for those without an open dataset
ACCESSIBLE = MEASURED + ["Cultured fibroblasts", "EBV lymphocytes", "Minor salivary gland"]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--atlas", required=True, help="public_data/panel_expression_atlas.tsv")
    p.add_argument("--measured", help="public_data/specimen_panel_tpm.tsv (scripts/specimen_expression.py)")
    p.add_argument("--genes", help="order of the highlighted genes (default: as in the table)")
    p.add_argument("--frequency", help="public_data/dee_gene_frequency.tsv: how often each gene explains epilepsy")
    p.add_argument("--outdir", required=True)
    return p.parse_args()


def classic(ax):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def save(fig, path):
    fig.savefig(path + ".pdf")
    fig.savefig(path + ".png", dpi=300)
    plt.close(fig)


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    with open(args.atlas) as fh:
        rows = [r for r in csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t")]
    if args.measured and os.path.exists(args.measured):
        with open(args.measured) as fh:
            meas = {r["gene"]: r for r in csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t")}
        for r in rows:   # a gene missing from a dataset's matrix had no reads there: 0
            m = meas.get(r["gene"], {})
            for c in MEASURED[1:]:
                r[c] = m.get(c, "0") if m.get(c, "NA") != "NA" else "0"
    cols = [c for c in SPEC if c in rows[0]]
    rows = [r for r in rows if all(r[c] != "NA" for c in cols)]
    n_meas = sum(1 for c in cols if c in MEASURED)
    genes = [r["gene"] for r in rows]
    targets = [r["gene"] for r in rows if r["target"] == "True"]
    if args.genes:
        targets = [g for g in args.genes.split(",") if g in targets] + [g for g in targets
                                                                         if g not in args.genes.split(",")]
    X = np.array([[float(r[c]) for c in cols] for r in rows])
    ci = {c: i for i, c in enumerate(cols)}
    out = [("panel_genes", "", len(genes))]

    # ---- coverage
    brain = X[:, ci["Brain cortex"]] >= 10
    fig, axes = plt.subplots(1, 2, figsize=(12.4, 4.3), gridspec_kw={"width_ratios": [1.5, 1]})
    fig.subplots_adjust(left=0.06, right=0.98, bottom=0.34, top=0.88, wspace=0.2)
    ax = axes[0]
    x = np.arange(len(cols))
    for k, (thr, alpha, lab) in enumerate(((1, 0.45, ">= 1 TPM"), (10, 1.0, ">= 10 TPM"))):
        v = (X >= thr).mean(axis=0) * 100
        ax.bar(x + (k - 0.5) * 0.38, v, width=0.38, color=[SPEC[c][2] for c in cols], alpha=alpha,
               edgecolor="black", lw=0.3, label=lab)
        for xi, vi in zip(x, v):
            ax.text(xi + (k - 0.5) * 0.38, vi + 1, "{:.0f}".format(vi), ha="center", fontsize=6.5)
        out += [("percent_ge{}".format(thr), c, round(vv, 1)) for c, vv in zip(cols, v)]
    ax.set_xticks(x)
    ax.set_xticklabels([SPEC[c][0] + (" *" if c not in MEASURED else "") for c in cols], rotation=40, ha="right",
                       rotation_mode="anchor", fontsize=8)
    ax.set_ylim(0, 105)
    ax.set_ylabel("% of panel genes")
    ax.set_title("a  {} panel genes: light >= 1 TPM, dark >= 10 TPM".format(len(genes)), loc="left", fontsize=10,
                 fontweight="bold")
    classic(ax)
    ax = axes[1]
    acc = [c for c in ACCESSIBLE if c in ci]
    v = [(X[brain, ci[c]] >= 10).mean() * 100 for c in acc]
    order = np.argsort(v)[::-1]
    ax.bar(range(len(acc)), [v[i] for i in order], color=[SPEC[acc[i]][2] for i in order], edgecolor="black",
           lw=0.3)
    for k, i in enumerate(order):
        ax.text(k, v[i] + 1, "{:.0f}%".format(v[i]), ha="center", fontsize=8)
    ax.set_xticks(range(len(acc)))
    ax.set_xticklabels([SPEC[acc[i]][0] + (" *" if acc[i] not in MEASURED else "") for i in order], rotation=40,
                       ha="right", rotation_mode="anchor", fontsize=8)
    ax.set_ylim(0, 105)
    ax.set_ylabel("% also >= 10 TPM")
    ax.set_title("b  The {} panel genes expressed in brain cortex".format(int(brain.sum())), loc="left", fontsize=10,
                 fontweight="bold")
    classic(ax)
    out += [("brain_expressed_genes", "", int(brain.sum()))]
    out += [("brain_expressed_recovered_percent_all", c, round((X[brain, ci[c]] >= 10).mean() * 100, 1))
            for c in cols if not c.startswith("Brain")]
    out += [("brain_expressed_recovered_percent", acc[i], round(v[i], 1)) for i in order]
    fig.text(0.07, 0.01, "* GTEx tissue (brain, or standing in for a specimen without an open dataset); others: "
                         "measured in the specimen (GEO)", fontsize=7.5, color="#555555")
    save(fig, os.path.join(args.outdir, "atlas_coverage"))

    # ---- reach: how many accessible specimens express each brain-expressed gene; best one and two specimens
    expr = X[:, [ci[c] for c in acc]] >= 10
    n_spec = expr[brain].sum(axis=1)
    best1 = max(acc, key=lambda c: expr[brain, acc.index(c)].mean())
    best2 = max(itertools.combinations(acc, 2),
                key=lambda p: (expr[brain][:, [acc.index(p[0]), acc.index(p[1])]].any(axis=1)).mean())
    cov1 = expr[brain, acc.index(best1)].mean() * 100
    cov2 = expr[brain][:, [acc.index(best2[0]), acc.index(best2[1])]].any(axis=1).mean() * 100
    covall = expr[brain].any(axis=1).mean() * 100
    none = [g for g, b, e in zip(genes, brain, expr.any(axis=1)) if b and not e]
    out += [("best_single_specimen", best1, round(cov1, 1)), ("best_pair", " + ".join(best2), round(cov2, 1)),
            ("any_accessible", "", round(covall, 1)), ("brain_only_genes", "", len(none))]
    fig, axes = plt.subplots(1, 2, figsize=(11.6, 3.6), gridspec_kw={"width_ratios": [1.2, 1]})
    fig.subplots_adjust(left=0.07, right=0.98, bottom=0.2, top=0.86, wspace=0.3)
    ax = axes[0]
    counts = [(n_spec == k).sum() for k in range(len(acc) + 1)]
    ax.bar(range(len(acc) + 1), counts, color=["#C0392B"] + ["#708090"] * len(acc), edgecolor="black", lw=0.3)
    for k, cnt in enumerate(counts):
        ax.text(k, cnt + 2, str(cnt), ha="center", fontsize=8)
    ax.set_xticks(range(len(acc) + 1))
    ax.set_xlabel("Accessible specimens where the gene reaches 10 TPM")
    ax.set_ylabel("Brain-expressed panel genes")
    ax.set_title("a  {} genes reach 10 TPM in no accessible specimen".format(len(none)), loc="left", fontsize=10,
                 fontweight="bold")
    classic(ax)
    ax = axes[1]
    labels = ["Best single:\n{}".format(SPEC[best1][0]), "Best pair:\n{} +\n{}".format(SPEC[best2[0]][0],
                                                                                       SPEC[best2[1]][0]),
              "Any of the {}\nspecimens".format(len(acc))]
    vals = [cov1, cov2, covall]
    ax.bar(range(3), vals, color=["#1B7837", "#4D9A57", "#9E9E9E"], edgecolor="black", lw=0.3, width=0.6)
    for k, vv in enumerate(vals):
        ax.text(k, vv + 1.5, "{:.0f}%".format(vv), ha="center", fontsize=9)
    ax.set_xticks(range(3))
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylim(0, 105)
    ax.set_ylabel("% of brain-expressed panel genes")
    ax.set_title("b  Covered at >= 10 TPM", loc="left", fontsize=10, fontweight="bold")
    classic(ax)
    save(fig, os.path.join(args.outdir, "atlas_reach"))
    with open(os.path.join(args.outdir, "brain_only_genes.txt"), "w") as fh:
        fh.write("\n".join(none) + "\n")

    # ---- heatmap
    acc_n = expr.sum(axis=1)
    order = sorted(range(len(genes)), key=lambda i: (-acc_n[i], -X[i, ci["Brain cortex"]]))
    L = np.log10(X[order] + 1)
    fig = plt.figure(figsize=(8.4, 9.6))
    ax = fig.add_axes([0.08, 0.25, 0.66, 0.72])
    im = ax.imshow(L, aspect="auto", cmap=HEAT, vmin=0, vmax=3, interpolation="nearest")
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels(["{} ({})".format(SPEC[c][0], SPEC[c][1]) for c in cols], rotation=50, ha="right",
                       rotation_mode="anchor", fontsize=7.5)
    ax.axvline(1.5 + n_meas, color="black", lw=1, ls="--")
    ax.set_yticks([])
    ax.set_ylabel("{} panel genes, grouped by the number of accessible specimens >= 10 TPM".format(len(genes)))
    ax.axvline(1.5, color="black", lw=1)
    # group boundaries
    srt = acc_n[order]
    for k in range(1, len(srt)):
        if srt[k] != srt[k - 1]:
            ax.axhline(k - 0.5, color="black", lw=0.4)
    pos = {genes[i]: k for k, i in enumerate(order)}
    tlist = sorted(targets, key=lambda g: pos[g])
    for j, g in enumerate(tlist):
        ytxt = (j + 0.5) / len(tlist) * len(genes)
        ax.annotate(g, xy=(len(cols) - 0.5, pos[g]), xytext=(len(cols) + 0.9, ytxt), fontsize=8, fontstyle="italic",
                    va="center", annotation_clip=False,
                    arrowprops=dict(arrowstyle="-", color="#2166AC", lw=0.7))
    cax = fig.add_axes([0.08, 0.045, 0.35, 0.013])
    cb = fig.colorbar(im, cax=cax, orientation="horizontal")
    cb.set_label("log10(median TPM + 1)", fontsize=8)
    save(fig, os.path.join(args.outdir, "atlas_heatmap"))

    # ---- each accessible specimen vs brain cortex
    acc2 = [c for c in ACCESSIBLE if c in ci]
    nc = 5
    nr = int(np.ceil(len(acc2) / nc))
    fig, axes = plt.subplots(nr, nc, figsize=(13.0, 2.9 * nr), squeeze=False)
    fig.subplots_adjust(left=0.06, right=0.99, bottom=0.1, top=0.93, wspace=0.28, hspace=0.45)
    bx = np.log10(X[:, ci["Brain cortex"]] + 1)
    for ax, c in zip(axes.flat, acc2):
        y = np.log10(X[:, ci[c]] + 1)
        rho = spearmanr(bx, y).correlation
        ax.scatter(bx, y, s=5, color=SPEC[c][2], alpha=0.5, lw=0)
        low = 0
        for g in sorted(targets, key=lambda g: bx[genes.index(g)]):
            i = genes.index(g)
            ax.scatter(bx[i], y[i], s=16, color="#2166AC", edgecolor="black", lw=0.4, zorder=3)
            ty = y[i]
            if y[i] < 0.3:   # stack the labels of genes near zero
                ty = 0.25 + 0.22 * low
                low += 1
                ax.plot([bx[i], bx[i] + 0.04], [y[i], ty], color="#2166AC", lw=0.4)
            ax.text(bx[i] + 0.05, ty, g, fontsize=6, fontstyle="italic", va="center")
        ax.plot([0, 3.5], [0, 3.5], color="#999999", lw=0.7, ls="--")
        ax.axhline(np.log10(11), color="#BBBBBB", lw=0.6)
        ax.set_xlim(0, 3.6)
        ax.set_ylim(0, 3.6)
        ax.set_title("{}{} (rho {:.2f})".format(SPEC[c][0], "" if c in MEASURED else " *", rho), fontsize=9,
                     loc="left")
        ax.set_xlabel("Brain cortex, log10(TPM + 1)", fontsize=8)
        ax.set_ylabel("log10(TPM + 1)", fontsize=8)
        classic(ax)
        out.append(("spearman_vs_brain_cortex", c, round(rho, 3)))
    for ax in list(axes.flat)[len(acc2):]:
        ax.axis("off")
    save(fig, os.path.join(args.outdir, "atlas_vs_brain"))

    # ---- highlighted genes x specimens
    if targets:
        T = np.array([X[genes.index(g)] for g in targets])
        h = 0.42 * len(targets) + 1.9
        fig, ax = plt.subplots(figsize=(12.4, h))
        fig.subplots_adjust(left=0.08, right=0.94, bottom=1.45 / h, top=1 - 0.1 / h)
        im = ax.imshow(np.log10(T + 1), aspect="auto", cmap=HEAT, vmin=0, vmax=2.6)
        for i in range(len(targets)):
            for j in range(len(cols)):
                v = T[i, j]
                ax.text(j, i, "{:.0f}".format(v) if v >= 10 else "{:.1f}".format(v) if v >= 0.1 else
                        "<0.1" if v > 0 else "0",
                        ha="center", va="center", fontsize=8, color="white" if np.log10(v + 1) > 1.8 else "black")
                out.append(("target_tpm", "{}|{}".format(targets[i], cols[j]), round(v, 2)))
        ax.set_yticks(range(len(targets)))
        ax.set_yticklabels(targets, fontstyle="italic")
        ax.set_xticks(range(len(cols)))
        ax.axvline(1.5 + n_meas, color="black", lw=1, ls="--")
        ax.set_xticklabels(["{}\n({})".format(SPEC[c][0], SPEC[c][1]) for c in cols], rotation=40, ha="right",
                           rotation_mode="anchor", fontsize=8)
        ax.axvline(1.5, color="black", lw=1)
        cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.01)
        cb.set_label("log10(TPM + 1)", fontsize=8)
        save(fig, os.path.join(args.outdir, "targets_specimens"))

    # ---- why these genes: published frequency
    if args.frequency and os.path.exists(args.frequency):
        with open(args.frequency) as fh:
            fq = [r for r in csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t")]
        fig, axes = plt.subplots(1, 2, figsize=(11.6, 3.9), gridspec_kw={"width_ratios": [1.6, 1]})
        fig.subplots_adjust(left=0.07, right=0.98, bottom=0.24, top=0.86, wspace=0.25)
        for ax, m, ylab, title in ((axes[0], "truty_posmd", "Molecular diagnoses",
                                    "a  Diagnoses per gene, 9,769 people tested (Truty 2019)"),
                                   (axes[1], "symonds_incidence", "Per 100,000 live births",
                                    "b  Incidence, onset < 3 years (Symonds 2019)")):
            rr = [r for r in fq if r["measure"] == m]
            v = [float(r["value"]) for r in rr]
            col = ["#2166AC" if r["gene"] in targets else "#BDBDBD" for r in rr]
            ax.bar(range(len(rr)), v, color=col, edgecolor="black", lw=0.3)
            for k, vv in enumerate(v):
                ax.text(k, vv * 1.01, "{:g}".format(vv), ha="center", va="bottom", fontsize=7)
            ax.set_xticks(range(len(rr)))
            ax.set_xticklabels([r["gene"] for r in rr], rotation=45, ha="right", rotation_mode="anchor",
                               fontstyle="italic")
            ax.set_ylabel(ylab)
            ax.set_title(title, loc="left", fontsize=10, fontweight="bold")
            classic(ax)
            out += [("frequency_" + m, r["gene"], r["value"]) for r in rr]
        fig.text(0.07, 0.01, "Blue: the seven genes analysed. PRRT2 and DEPDC5 mostly cause self-limited or focal "
                             "epilepsy; PCDH19 affects females and mosaic males.", fontsize=7.5, color="#555555")
        save(fig, os.path.join(args.outdir, "gene_frequency"))

    with open(os.path.join(args.outdir, "panel_atlas_figures.tsv"), "w") as fh:
        fh.write("measure\tkey\tvalue\n")
        for r in out:
            fh.write("{}\t{}\t{}\n".format(*r))
    print("panel atlas figures -> {}".format(args.outdir))
    for r in out:
        if r[0] in ("brain_expressed_genes", "best_single_specimen", "best_pair", "any_accessible",
                    "brain_only_genes") or r[0].startswith("brain_expressed_recovered") or r[0].startswith("spearman"):
            print("  ", r)


if __name__ == "__main__":
    main()
