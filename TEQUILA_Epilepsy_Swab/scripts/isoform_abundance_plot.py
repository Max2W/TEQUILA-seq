#!/usr/bin/env python3
"""Isoform abundance (CPM) of one gene: stacked bars of its most abundant isoforms plus "Other".

Two figures, from the same numbers:
  <prefix>_isoform_abundance_by_donor   one panel per donor (grey strip on top), one bar per kit; a horizontal
                                        legend; below it the structures of the plotted isoforms, the canonical
                                        (MANE) transcript first in white when it is not among them
  <prefix>_isoform_abundance            one panel, one bar per sample (kits in turn, donors within), legend on
                                        the right
CPM = isoform reads / all reads sequenced for the sample (qc/fastq_read_counts.tsv) x 10^6.
The isoforms shown are the --top with the highest mean CPM over the samples; all other isoforms of the gene are
summed as "Other". Stack order, top to bottom, follows the legend. Labels are the transcript IDs of the matrix
(GENCODE ENST IDs, or ESPRESSO / AMALGAM IDs for novel isoforms).

Inputs: the ESPRESSO (samples_N2_R0_abundance.esp + samples_N2_R0_updated.gtf) or AMALGAM (counts matrix +
transcripts GTF) output, the GENCODE GTF and the sample sheet. Writes <prefix>_isoform_abundance.tsv with every
drawn number.
"""
import argparse
import csv
import os
import sys
from collections import defaultdict

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figure_style import IntronShrink, canonical_model, ordered, plt, read_sheet, save  # noqa: E402
from isoform_report import overlaps, read_gtf_exons, strip_version  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

# by-donor figure: pink, yellow, periwinkle, green (+ more for --top > 4); Other grey
PALETTE_DONOR = ["#F0527B", "#FED976", "#7A87C6", "#A8CC3A", "#66C2A5", "#E78AC3", "#A6761D"]
OTHER_DONOR = "#808080"
# single-panel figure: Okabe-Ito; Other light grey
PALETTE_SINGLE = ["#56B4E9", "#E69F00", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"]
OTHER_SINGLE = "#999999"


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--gene", required=True)
    p.add_argument("--matrix", required=True, help="transcript_ID, transcript_name, gene_ID, one column per sample")
    p.add_argument("--isoform-gtf", required=True, help="structures of the isoforms in the matrix")
    p.add_argument("--ref-gtf", required=True, help="GENCODE GTF")
    p.add_argument("--sample-sheet", required=True)
    p.add_argument("--read-counts", required=True, help="fastq_read_counts.tsv (sample_id, n_reads)")
    p.add_argument("--method", default="ESPRESSO", help="label for the file names and the data table")
    p.add_argument("--top", type=int, default=4, help="isoforms drawn; the rest are 'Other'")
    p.add_argument("--outdir", required=True)
    return p.parse_args()


def gene_isoforms(args):
    """counts (isoform x sample) of the gene's isoforms, their structures, and the canonical model."""
    ref = read_gtf_exons(args.ref_gtf, keep=lambda f, a: a.get("gene_name") == args.gene)
    if not ref:
        sys.exit("{} not in {}".format(args.gene, args.ref_gtf))
    any_tx = next(iter(ref.values()))
    chrom, strand = any_tx["chrom"], any_tx["strand"]
    gene_ids = {strip_version(r["attrs"].get("gene_id", "")) for r in ref.values()}
    gene_exons = sorted({e for r in ref.values() for e in r["exons"]})

    m = pd.read_csv(args.matrix, sep="\t", dtype={0: str, 1: str, 2: str})
    m.columns = [c.lstrip("#") for c in m.columns]
    id_col, gene_col = m.columns[0], m.columns[2]
    structures = read_gtf_exons(args.isoform_gtf)

    def belongs(row):
        rec = structures.get(row[id_col])
        if rec is None:
            return False
        if {strip_version(g) for g in str(row[gene_col]).split(",")} & gene_ids:
            return True
        return rec["chrom"] == chrom and rec["strand"] == strand and overlaps(rec["exons"], gene_exons)

    sub = m[m.apply(belongs, axis=1)].drop_duplicates(subset=id_col).set_index(id_col)
    if sub.empty:
        sys.exit("no isoform of {} in {}".format(args.gene, args.matrix))
    return sub, {t: structures[t] for t in sub.index}, ref


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    prefix = os.path.join(args.outdir, "{}_{}".format(args.gene, args.method))
    sheet = read_sheet(args.sample_sheet)
    conditions = ordered(r["condition"] for r in sheet)
    donors = ordered(r["donor"] for r in sheet)
    by_key = {(r["donor"], r["condition"]): r["sample_id"] for r in sheet}
    samples = [by_key[(d, c)] for c in conditions for d in donors if (d, c) in by_key]

    total = defaultdict(int)
    with open(args.read_counts) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            total[row["sample_id"]] += int(row["n_reads"])
    missing = [s for s in samples if not total[s]]
    if missing:
        sys.exit("no read count for {} in {}".format(", ".join(missing), args.read_counts))

    sub, structures, ref = gene_isoforms(args)
    counts = pd.DataFrame({s: sub[s].astype(float) if s in sub.columns else 0.0 for s in samples}, index=sub.index)
    cpm = counts.div(pd.Series({s: total[s] for s in samples})) * 1e6
    top = cpm.mean(axis=1).sort_values(ascending=False).index[:args.top].tolist()
    other = [t for t in cpm.index if t not in top]
    shown = top + (["Other"] if other else [])
    plot = cpm.loc[top].copy()
    if other:
        plot.loc["Other"] = cpm.loc[other].sum(axis=0)

    rows = []
    for t in shown:
        for s in samples:
            n = counts.loc[other, s].sum() if t == "Other" else counts.loc[t, s]
            rows.append([args.gene, args.method, t, len(other) if t == "Other" else 1, s, round(float(n), 2),
                         total[s], round(float(plot.loc[t, s]), 4)])
    with open(prefix + "_isoform_abundance.tsv", "w") as fh:
        fh.write("gene\tmethod\tisoform\tn_isoforms\tsample_id\treads\treads_sequenced\tcpm\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")

    canon = canonical_model(args.ref_gtf, args.gene)
    canon_id = canon["tx_id"]
    tx_rows = [] if any(strip_version(t) == strip_version(canon_id) for t in top) else [(canon_id, canon["exons"], None)]
    tx_rows += [(t, structures[t]["exons"], i) for i, t in enumerate(top)]

    sample_meta = {r["sample_id"]: r for r in sheet}
    by_donor_figure(args, prefix, plot, shown, donors, conditions, by_key, tx_rows)
    single_figure(args, prefix, plot, shown, samples, sample_meta)
    print("{} {}: {} isoforms, {} shown -> {}_isoform_abundance*.pdf/.png".format(
        args.gene, args.method, len(cpm), len(top), prefix))


def nice_ylim(top_value):
    top_value = top_value if top_value > 0 else 1.0
    return -0.05 * top_value, 1.05 * top_value


def stack(ax, x, values_by_isoform, shown, colors, width=0.9):
    """bars stacked bottom -> top in the reverse of the legend order, so the legend reads top -> bottom."""
    bottom = np.zeros(len(x))
    for t in reversed(shown):
        v = np.asarray(values_by_isoform[t], float)
        ax.bar(x, v, width=width, bottom=bottom, color=colors[t], edgecolor="black", linewidth=0.7, zorder=2)
        bottom += v
    return bottom.max() if len(bottom) else 0.0


def text_width_in(fig, text, **kw):
    t = fig.text(0, 0, text, **kw)
    w = t.get_window_extent(renderer=fig.canvas.get_renderer()).width / fig.dpi
    t.remove()
    return w


def by_donor_figure(args, prefix, plot, shown, donors, conditions, by_key, tx_rows):
    colors = {t: PALETTE_DONOR[i % len(PALETTE_DONOR)] for i, t in enumerate(s for s in shown if s != "Other")}
    colors["Other"] = OTHER_DONOR
    W = 6.7
    left, right, gap = 0.62, 0.08, 0.10
    strip_h, panel_h, xlab_h = 0.28, 2.0, 0.75
    row_h, line_h, sq = 0.42, 0.32, 0.2
    fig = plt.figure(figsize=(W, 4))
    # legend items ("<GENE> isoforms" title, then square + label), wrapped into rows that fit the width
    gene_w = text_width_in(fig, args.gene, fontsize=11, style="italic")
    items = [("title", gene_w + text_width_in(fig, " isoforms", fontsize=11))] + \
            [(t, sq + 0.06 + text_width_in(fig, t, fontsize=11)) for t in shown]
    rows, cur, cur_w = [], [], 0.0
    for name, w in items:
        if cur and cur_w + 0.18 + w > W - 0.2:
            rows.append((cur, cur_w))
            cur, cur_w = [], 0.0
        cur_w += (0.18 if cur else 0.0) + w
        cur.append((name, w))
    rows.append((cur, cur_w))
    legend_h = 0.2 + line_h * len(rows)
    H = 0.08 + strip_h + panel_h + xlab_h + legend_h + row_h * len(tx_rows) + 0.15
    fig.set_size_inches(W, H)
    pw = (W - left - right - gap * (len(donors) - 1)) / len(donors)
    ymax = max(plot[[by_key[(d, c)] for d in donors for c in conditions if (d, c) in by_key]].sum(axis=0).max(), 0)
    ylim = nice_ylim(ymax)
    panel_bottom = H - 0.08 - strip_h - panel_h
    first = None
    for i, d in enumerate(donors):
        x0 = left + i * (pw + gap)
        ax = fig.add_axes([x0 / W, panel_bottom / H, pw / W, panel_h / H], sharey=first)
        first = first or ax
        conds = [c for c in conditions if (d, c) in by_key]
        ids = [by_key[(d, c)] for c in conds]
        stack(ax, np.arange(len(ids)), {t: plot.loc[t, ids].values for t in shown}, shown, colors)
        ax.set_xlim(-0.6, len(ids) - 0.4)
        ax.set_ylim(*ylim)
        ax.set_xticks(np.arange(len(ids)))
        ax.set_xticklabels(conds, rotation=45, ha="right", rotation_mode="anchor")
        if i:
            ax.tick_params(axis="y", labelleft=False, left=False)
        else:
            ax.set_ylabel("Isoform abundance (CPM)")
        strip = fig.add_axes([x0 / W, (panel_bottom + panel_h) / H, pw / W, strip_h / H])
        strip.set_xticks([]), strip.set_yticks([])
        strip.set_facecolor("#D9D9D9")
        strip.text(0.5, 0.5, "Donor {}".format(d), ha="center", va="center", fontsize=11)

    # legend row
    for k, (row, row_w) in enumerate(rows):
        ly = panel_bottom - xlab_h - 0.1 - line_h * (k + 0.5)
        x = max(0.05, (W - row_w) / 2)
        for name, w in row:
            if name == "title":
                fig.text(x / W, ly / H, args.gene, fontsize=11, style="italic", va="center")
                fig.text((x + gene_w) / W, ly / H, " isoforms", fontsize=11, va="center")
            else:
                fig.patches.append(Rectangle((x / W, (ly - sq / 2) / H), sq / W, sq / H, transform=fig.transFigure,
                                             facecolor=colors[name], edgecolor="black", linewidth=0.7))
                fig.text((x + sq + 0.06) / W, ly / H, name, fontsize=11, va="center")
            x += w + 0.18

    # structures
    label_w = max(text_width_in(fig, t, fontsize=11) for t, _, _ in tx_rows) + 0.15
    sx0 = max(left, label_w + 0.1)
    shrink = IntronShrink([e for _, ex, _ in tx_rows for e in ex])
    top_y = panel_bottom - xlab_h - legend_h
    ax = fig.add_axes([sx0 / W, 0.12 / H, (W - sx0 - right) / W, (top_y - 0.12) / H])
    ax.set_xlim(-0.01, 1.01)
    ax.set_ylim(len(tx_rows) - 0.5, -0.5)
    ax.axis("off")
    order = [t for t in shown if t != "Other"]
    for r, (tid, exons, idx) in enumerate(tx_rows):
        fill = "white" if idx is None else colors[order[idx]]
        xs = [(shrink(s), shrink(e + 1)) for s, e in exons]
        ax.plot([xs[0][0], xs[-1][1]], [r, r], color="black", linewidth=1.0, zorder=1, solid_capstyle="butt")
        for a, b in xs:
            ax.add_patch(Rectangle((a, r - 0.2), max(b - a, 0.002), 0.4, facecolor=fill, edgecolor="black",
                                   linewidth=0.7, zorder=2))
        ax.annotate(tid, xy=(0, r), xycoords=("axes fraction", "data"), xytext=(-8, 0),
                    textcoords="offset points", ha="right", va="center", fontsize=11, annotation_clip=False)
    save(fig, prefix + "_isoform_abundance_by_donor")


def single_figure(args, prefix, plot, shown, samples, sample_meta):
    colors = {t: PALETTE_SINGLE[i % len(PALETTE_SINGLE)] for i, t in enumerate(s for s in shown if s != "Other")}
    colors["Other"] = OTHER_SINGLE
    labels = ["{} (donor {})".format(sample_meta[s]["condition"], sample_meta[s]["donor"]) for s in samples]
    fig = plt.figure(figsize=(6.7, 4.4))
    tmp_w = max(text_width_in(fig, t, fontsize=13) for t in shown) + 0.45
    legend_w = max(tmp_w, text_width_in(fig, args.gene + " isoforms", fontsize=14)) + 0.2
    xlab_h = max(text_width_in(fig, l, fontsize=13) for l in labels) + 0.25
    left, right, top = 0.75, legend_w + 0.25, 0.1
    pw, ph = 6.7 - left - right, 4.4 - top - xlab_h
    if ph < 1.5:
        fig.set_size_inches(6.7, 4.4 + 1.5 - ph)
        ph = 1.5
    H = fig.get_size_inches()[1]
    ax = fig.add_axes([left / 6.7, xlab_h / H, pw / 6.7, ph / H])
    x = np.arange(len(samples))
    ymax = stack(ax, x, {t: plot.loc[t, samples].values for t in shown}, shown, colors)
    ax.set_xlim(-0.6, len(samples) - 0.4)
    ax.set_ylim(*nice_ylim(ymax))
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=90, fontsize=13)
    ax.tick_params(axis="y", labelsize=13)
    ax.set_ylabel("Isoform abundance (CPM)", fontsize=15)
    handles = [Rectangle((0, 0), 1, 1, facecolor=colors[t], edgecolor="black", linewidth=0.7) for t in shown]
    leg = ax.legend(handles, shown, loc="center left", bbox_to_anchor=(1.03, 0.6), frameon=False, fontsize=13,
                    handlelength=1.1, handleheight=1.1, title="$\\it{" + args.gene + "}$ isoforms",
                    title_fontsize=14, alignment="left")
    leg._legend_box.align = "left"
    save(fig, prefix + "_isoform_abundance")


if __name__ == "__main__":
    main()
