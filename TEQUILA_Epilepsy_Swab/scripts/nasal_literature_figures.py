#!/usr/bin/env python3
"""Literature figures for the next specimen, nasal epithelial brushing, against buccal swabs (progress-deck style).

Input: public_data/nasal_literature.tsv (numbers transcribed from the full texts, with study and location of each
value); optional --our-rin (this pilot's swab RIN, kept on scratch) and --pediatric (public_data/nasal_pediatric.tsv).
Figures (<outdir>/):
  nasal_rin         RNA integrity (RIN) by sampling site and device
  nasal_yield       total RNA per sample (ng, log scale)
  nasal_cells       what is collected: (a) % respiratory vs squamous cells by nasal site and device (Lai 2015),
                    (b) epithelial vs immune cells in nasopharyngeal swabs (Ziegler 2021), (c) viability of buccal
                    epithelial cells (Rudney 2006)
  nasal_success     samples with enough RNA / passing library QC
  nasal_pediatric   (with --pediatric) success by age group
Table: nasal_literature_figures.tsv (every value drawn).
"""
import argparse
import csv
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
SITES = {"mouth": ("Mouth", "#EC7B02"), "nose_front": ("Front of the nose", "#9ECAE1"),
         "nose_deep": ("Inside the nose (turbinate, nasopharynx)", "#2166AC"), "scalp": ("Scalp hair", "#9E9E9E"),
         "pilot": ("This pilot", "#000000")}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--table", required=True, help="public_data/nasal_literature.tsv")
    p.add_argument("--pediatric", help="public_data/nasal_pediatric.tsv (optional)")
    p.add_argument("--expression", help="public_data/nasal_gse152004_summary.tsv (optional)")
    p.add_argument("--our-rin", help="TSV with sample_id, RIN of this pilot's buccal swabs (optional)")
    p.add_argument("--outdir", required=True)
    return p.parse_args()


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"))


def classic(ax):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def save(fig, path):
    fig.savefig(path + ".pdf")
    fig.savefig(path + ".png", dpi=300)
    plt.close(fig)


def num(x):
    return float(x) if x not in ("", None) else None


def site_legend(fig, sites, note=None):
    """site colours at the bottom of the figure (fixed 1.0 in band), optional footnote under them."""
    h = fig.get_size_inches()[1]
    fig.legend(handles=[Patch(color=SITES[s][1], label=SITES[s][0]) for s in sites], frameon=False, ncol=3,
               fontsize=7.5, loc="lower left", bbox_to_anchor=(0.01, 0.30 / h))
    if note:
        fig.text(0.01, 0.06 / h, note, fontsize=7, color="#555555")


def dot_chart(rows, xlabel, path, log=False, extra=None, ref=None):
    """horizontal dots with intervals, sorted by value; rows are dicts of the table."""
    items = [(r["label"], r["site"], float(r["value"]), num(r["lo"]), num(r["hi"]), r["interval"],
              "{}, n = {}".format(r["study"], r["n"])) for r in rows]
    if extra:
        items.append(extra)
    items.sort(key=lambda t: t[2])
    h = 0.42 * len(items) + 1.9
    fig, ax = plt.subplots(figsize=(8.2, h))
    fig.subplots_adjust(left=0.30, right=0.72, bottom=1.45 / h, top=1 - 0.25 / h)
    for i, (lab, site, v, lo, hi, kind, src) in enumerate(items):
        col = SITES[site][1]
        if lo is not None and hi is not None:
            ax.plot([lo, hi], [i, i], color=col, lw=2.2, solid_capstyle="butt", alpha=0.7)
        ax.scatter(v, i, s=46, color=col, edgecolor="black", lw=0.5, zorder=3)
        ax.annotate(src, xy=(1.01, i), xycoords=("axes fraction", "data"), va="center", fontsize=7.5,
                    color="#555555")
    ax.set_yticks(range(len(items)))
    ax.set_yticklabels([t[0] for t in items])
    ax.set_ylim(-0.6, len(items) - 0.4)
    ax.set_xlabel(xlabel)
    if log:
        ax.set_xscale("log")
    if ref is not None:
        ax.axvline(ref[0], color="#888888", lw=0.8, ls="--", zorder=0)
        ax.text(ref[0], len(items) - 0.45, ref[1], ha="center", va="bottom", fontsize=7.5, color="#555555")
    classic(ax)
    return fig, ax, items


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    tab = read_tsv(args.table)
    out = []

    # ---- RIN
    rows = [r for r in tab if r["measure"] == "RIN"]
    extra = None
    if args.our_rin and os.path.exists(args.our_rin):
        rin = [float(r["RIN"]) for r in read_tsv(args.our_rin) if r.get("RIN")]
        if rin:
            extra = ("Buccal swabs, this pilot", "pilot", float(np.mean(rin)), min(rin), max(rin), "range",
                     "this pilot, n = {}".format(len(rin)))
    fig, ax, items = dot_chart(rows, "RIN (RNA integrity number; 10 = intact)", None, extra=extra,
                               ref=(7, "RIN 7"))
    ax.set_xlim(0, 10)
    sites = [s for s in SITES if any(t[1] == s for t in items)]
    site_legend(fig, sites, "Dots: median or mean as reported; bars: IQR, mean +/- SD or range "
                            "(public_data/nasal_literature.tsv)")
    save(fig, os.path.join(args.outdir, "nasal_rin"))
    out += [("RIN", t[0], t[2], t[3], t[4], t[6]) for t in items]

    # ---- yield
    rows = [r for r in tab if r["measure"] == "yield_ng"]
    fig, ax, items = dot_chart(rows, "Total RNA per sample (ng, log scale)", None, log=True)
    sites = [s for s in SITES if any(t[1] == s for t in items)]
    site_legend(fig, sites, "Dots: median or mean as reported; bars: IQR or mean +/- SD; Chu 2016: mean, SD 437 ng")
    save(fig, os.path.join(args.outdir, "nasal_yield"))
    out += [("yield_ng", t[0], t[2], t[3], t[4], t[6]) for t in items]

    # ---- cells
    fig, axes = plt.subplots(1, 3, figsize=(11.0, 3.9), gridspec_kw={"width_ratios": [2.6, 1.3, 1.3]})
    fig.subplots_adjust(left=0.06, right=0.97, bottom=0.36, top=0.84, wspace=0.75)
    ax = axes[0]
    labels = ["Anterior nares swab", "Anterior nares brush", "Inferior turbinate brush"]
    get = lambda m, lab: next((float(r["value"]) for r in tab if r["measure"] == m and r["label"] == lab), None)
    resp = [get("pct_respiratory", l) for l in labels]
    squa = [get("pct_squamous", l) for l in labels]
    for i, (rp, sq) in enumerate(zip(resp, squa)):
        sq = sq or 0.0
        ax.bar(i, rp, color="#2166AC", width=0.65)
        ax.bar(i, sq, bottom=rp, color="#EC7B02", width=0.65)
        other = max(0.0, 100 - rp - sq)
        if other > 0.5:
            ax.bar(i, other, bottom=rp + sq, color="#D9D9D9", width=0.65)
        ax.text(i, rp / 2, "{:.0f}%".format(rp), ha="center", va="center", color="white", fontsize=8.5)
        if sq >= 5:
            ax.text(i, rp + sq / 2, "{:.0f}%".format(sq), ha="center", va="center", fontsize=8.5)
        out += [("pct_respiratory", labels[i], rp, "", "", "Lai 2015"), ("pct_squamous", labels[i], sq, "", "",
                                                                         "Lai 2015")]
    ax.set_xticks(range(3))
    ax.set_xticklabels(labels, rotation=20, ha="right", rotation_mode="anchor")
    ax.set_ylim(0, 100)
    ax.set_ylabel("% of cells (median)")
    ax.set_title("a  Nasal site and device (Lai 2015, n = 12)", loc="left", fontsize=9, fontweight="bold")
    ax.legend(handles=[Patch(color="#2166AC", label="respiratory epithelial"), Patch(color="#EC7B02",
                                                                                  label="squamous"),
                       Patch(color="#D9D9D9", label="other / not reported")], frameon=False, fontsize=7.5,
              loc="upper center", bbox_to_anchor=(0.5, -0.33), ncol=3)
    classic(ax)
    ax = axes[1]
    ep, im = get("cells_epithelial", "Nasopharyngeal swab (scRNA-seq)"), get("cells_immune",
                                                                              "Nasopharyngeal swab (scRNA-seq)")
    tot = ep + im
    ax.bar(0, ep / tot * 100, color="#2166AC", width=0.6)
    ax.bar(0, im / tot * 100, bottom=ep / tot * 100, color="#8C6BB1", width=0.6)
    ax.text(0, ep / tot * 50, "epithelial\n{:.0f}%".format(ep / tot * 100), ha="center", va="center",
            color="white", fontsize=8)
    ax.text(0.42, ep / tot * 100 + im / tot * 50, "immune {:.0f}%".format(im / tot * 100), ha="left",
            va="center", fontsize=8)
    ax.set_xlim(-0.6, 1.6)
    ax.set_xticks([0])
    ax.set_xticklabels(["Nasopharyngeal\nswab"])
    ax.set_ylim(0, 100)
    ax.set_ylabel("% of cells sequenced")
    ax.set_title("b  scRNA-seq\n(Ziegler 2021, {:,} cells)".format(int(tot)), loc="left", fontsize=9,
                 fontweight="bold")
    classic(ax)
    out += [("pct_epithelial", "Nasopharyngeal swab", ep / tot * 100, "", "", "Ziegler 2021")]
    ax = axes[2]
    parts = [("viable", get("pct_viable", "Buccal epithelial cells"), "#1B7837"),
             ("apoptotic", get("pct_apoptotic", "Buccal epithelial cells"), "#FDB863"),
             ("dead", get("pct_dead", "Buccal epithelial cells"), "#B2182B")]
    b = 0
    for name, v, col in parts:
        ax.bar(0, v, bottom=b, color=col, width=0.6)
        ax.text(0.42, b + v / 2, "{} {:.0f}%".format(name, v), ha="left", va="center", fontsize=8)
        b += v
        out.append(("pct_" + name, "Buccal epithelial cells", v, "", "", "Rudney 2006"))
    ax.set_xlim(-0.6, 1.6)
    ax.set_xticks([0])
    ax.set_xticklabels(["Buccal\nepithelial cells"])
    ax.set_ylim(0, 100)
    ax.set_ylabel("% of cells")
    ax.set_title("c  Buccal cells\n(Rudney 2006, n = 36)", loc="left", fontsize=9, fontweight="bold")
    classic(ax)
    save(fig, os.path.join(args.outdir, "nasal_cells"))

    # ---- success
    rows = []
    for r in tab:
        if r["measure"] == "pass_n":
            k, n = float(r["value"]), float(r["lo"])
            rows.append((r["label"], r["site"], k / n * 100, "{:.0f}/{:.0f}".format(k, n), r["study"]))
        elif r["measure"] == "pass_pct":
            rows.append((r["label"], r["site"], float(r["value"]), "{:.0f}%".format(float(r["value"])),
                         r["study"]))
    rows.sort(key=lambda t: t[2])
    h = 0.42 * len(rows) + 1.9
    fig, ax = plt.subplots(figsize=(8.2, h))
    fig.subplots_adjust(left=0.30, right=0.80, bottom=1.45 / h, top=1 - 0.15 / h)
    for i, (lab, site, pct, txt, study) in enumerate(rows):
        ax.barh(i, pct, color=SITES[site][1], height=0.65)
        ax.text(pct + 1, i, txt, va="center", fontsize=8)
        ax.annotate(study, xy=(1.03, i), xycoords=("axes fraction", "data"), va="center", fontsize=7.5,
                    color="#555555")
        out.append(("percent_passing", lab, round(pct, 1), "", "", study))
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([t[0] for t in rows])
    ax.set_xlim(0, 112)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_xlabel("% of samples with enough RNA or passing library QC")
    sites = [s for s in SITES if any(t[1] == s for t in rows)]
    site_legend(fig, sites, "Lai 2015: enough RNA for a microarray; the others: libraries passing QC or adequate "
                            "RNA for sequencing")
    classic(ax)
    save(fig, os.path.join(args.outdir, "nasal_success"))

    # ---- pediatric (optional table: age_group, measure, value, n, study, location)
    if args.pediatric and os.path.exists(args.pediatric):
        ped = read_tsv(args.pediatric)
        outcome = {"tolerated": ("tolerated", "#9ECAE1"), "culture": ("cells obtained / grew in culture", "#4292C6"),
                   "rna": ("enough cells for RNA", "#08519C"), "library_qc": ("RNA library passed QC", "#1B7837")}
        h = 0.5 * len(ped) + 1.5
        fig, ax = plt.subplots(figsize=(8.2, h))
        fig.subplots_adjust(left=0.36, right=0.74, bottom=1.15 / h, top=1 - 0.15 / h)
        for i, r in enumerate(ped):
            v = float(r["value"])
            ax.barh(i, v, color=outcome.get(r["measure"], ("", "#708090"))[1], height=0.65)
            ax.text(v + 1, i, "{:.0f}%  ({})".format(v, r["n"]), va="center", fontsize=8)
            ax.annotate(r["study"], xy=(1.18, i), xycoords=("axes fraction", "data"), va="center", fontsize=7.5,
                        color="#555555")
            out.append(("pediatric_" + r["measure"], "{} | {}".format(r["age_group"], r["label"]), v, "", "",
                        r["study"]))
        ax.set_yticks(range(len(ped)))
        ax.set_yticklabels(["{}: {}".format(r["age_group"], r["label"]) for r in ped], fontsize=8)
        ax.set_xlim(0, 118)
        ax.set_xticks([0, 25, 50, 75, 100])
        ax.set_xlabel("% of participants")
        ax.invert_yaxis()
        used = list(dict.fromkeys(r["measure"] for r in ped))
        fig.legend(handles=[Patch(color=outcome[m][1], label=outcome[m][0]) for m in used if m in outcome],
                   frameon=False, ncol=4, fontsize=7.5, loc="lower left", bbox_to_anchor=(0.01, 0.1 / h))
        classic(ax)
        save(fig, os.path.join(args.outdir, "nasal_pediatric"))

    # ---- expression of the target genes in nasal brushings (GSE152004)
    if args.expression and os.path.exists(args.expression):
        ex = {r["gene"]: r for r in read_tsv(args.expression) if r["median"] != "NA"}
        targets = [g for g, r in ex.items() if "target" in r["set"]]
        show = targets + [g for g in ("FOXJ1", "SCGB1A1", "KRT5", "KRT13", "PTPRC", "OMP") if g in ex]
        fig, axes = plt.subplots(1, 2, figsize=(10.6, 3.9), gridspec_kw={"width_ratios": [1.5, 1]})
        fig.subplots_adjust(left=0.07, right=0.98, bottom=0.2, top=0.86, wspace=0.3)
        ax = axes[0]
        for i, g in enumerate(show):
            r = ex[g]
            q = [float(r[k]) + 0.01 for k in ("q10", "q25", "median", "q75", "q90")]
            col = "#2166AC" if g in targets else "#9E9E9E"
            ax.plot([i, i], [q[0], q[4]], color=col, lw=1)
            ax.add_patch(plt.Rectangle((i - 0.3, q[1]), 0.6, q[3] - q[1], facecolor=col, alpha=0.35,
                                       edgecolor=col))
            ax.plot([i - 0.3, i + 0.3], [q[2], q[2]], color=col, lw=2)
            ax.text(i, q[4] * 1.4, "{:.0f}".format(float(r["median"])) if float(r["median"]) >= 1 else
                    "{:.1f}".format(float(r["median"])), ha="center", fontsize=7.5)
            out.append(("nasal_cpm_median", g, float(r["median"]), r["q25"], r["q75"], "GSE152004"))
        ax.set_yscale("log")
        ax.set_xticks(range(len(show)))
        ax.set_xticklabels(show, fontstyle="italic", rotation=30, ha="right", rotation_mode="anchor")
        ax.set_xlim(-0.6, len(show) - 0.4)
        ax.set_ylabel("Counts per million")
        ax.set_title("a  Target genes and cell markers", loc="left", fontsize=9, fontweight="bold")
        classic(ax)
        ax = axes[1]
        pm = np.array([float(r["median"]) for g, r in ex.items() if "panel" in r["set"]])
        bins = np.logspace(-2, 4, 37)
        ax.hist(np.clip(pm, 0.011, None), bins=bins, color="#708090")
        ax.set_xscale("log")
        ax.axvline(10, color="#888888", ls="--", lw=0.8)
        share = (pm >= 10).mean() * 100
        top = ax.get_ylim()[1]
        ax.set_ylim(0, top * 1.25)
        ax.text(0.99, 0.50, "{:.0f}% of panel genes\n>= 10 CPM".format(share), transform=ax.transAxes,
                fontsize=8, va="center", ha="right")
        ax.text(0.03, 0.40, "<= 0.01\n(incl. not\ndetected)", transform=ax.transAxes, fontsize=7, va="bottom",
                color="#555555")
        for k, g in enumerate(targets):
            v = float(ex[g]["median"])
            ax.annotate(g, xy=(v, top * 1.02), xytext=(v * (0.12 if k % 2 == 0 else 8), top * 1.12), ha="center",
                        fontsize=8, fontstyle="italic", color="#2166AC",
                        arrowprops=dict(arrowstyle="->", color="#2166AC", lw=0.8))
        ax.set_xlabel("Median CPM across samples")
        ax.set_ylabel("Panel genes")
        ax.set_title("b  {} probe-panel genes".format(len(pm)), loc="left", fontsize=9, fontweight="bold")
        classic(ax)
        out.append(("nasal_panel_percent_ge10cpm", "panel", round(share, 1), "", "", "GSE152004"))
        fig.text(0.07, 0.015, "Nasal airway brushings, n = 695 (GEO GSE152004; Sajuthi 2020); boxes: IQR, line: "
                              "median, whiskers: 10th-90th percentile", fontsize=7, color="#555555")
        save(fig, os.path.join(args.outdir, "nasal_expression"))

    with open(os.path.join(args.outdir, "nasal_literature_figures.tsv"), "w") as fh:
        fh.write("measure\tlabel\tvalue\tlo\thi\tsource\n")
        for r in out:
            fh.write("\t".join("" if x is None else str(x) for x in r) + "\n")
    print("nasal literature figures -> {}".format(args.outdir))


if __name__ == "__main__":
    main()
