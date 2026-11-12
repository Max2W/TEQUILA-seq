#!/usr/bin/env python3
"""Read coverage of one gene in every sample, from the gene-filtered BAMs.

Figures
  <GENE>_coverage_lines        one line per sample, overlaid: all reads (top) and spliced reads only (bottom), on the
                               canonical (MANE) transcript with its introns shortened; kits in blue / red shades, one
                               shade per donor; the gene model below (UTR white, CDS grey)
  <GENE>_coverage_tracks       IGV-style: one coverage track per sample, each scaled to its own maximum ("[0 - max]"),
                               one colour per donor, genomic coordinates; the gene model below (CDS thick, UTR thin)
  <GENE>_coverage_tracks_spliced   the same for spliced reads only
Coverage counts aligned bases only (no deletions, no introns). Every drawn value is written to
<GENE>_coverage.tsv.gz (bin, genomic start/end, sample, all-read and spliced-read coverage, for both x-axes).
"""
import argparse
import gzip
import os
import sys

import numpy as np
import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figure_style import (CONDITION_SHADES, DONOR_COLORS, IntronShrink, canonical_model, ordered, plt,  # noqa: E402
                          read_sheet, save)
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

W = 6.7


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--gene", required=True)
    p.add_argument("--gtf", required=True, help="GENCODE GTF")
    p.add_argument("--sample-sheet", required=True)
    p.add_argument("--bam-pattern", required=True, help="path with {gene} and {sample}")
    p.add_argument("--outdir", required=True)
    p.add_argument("--bins", type=int, default=1500, help="points along the x-axis")
    return p.parse_args()


def coverage(bam_path, chrom, lo, hi):
    """per-base coverage over lo..hi (1-based, inclusive): all primary reads, and spliced reads only."""
    n = hi - lo + 1
    d_all, d_spl = np.zeros(n + 1, np.int64), np.zeros(n + 1, np.int64)
    with pysam.AlignmentFile(bam_path) as bam:
        if chrom not in bam.references:
            return np.zeros(n), np.zeros(n)
        for r in bam.fetch(chrom, lo - 1, hi):
            if r.is_unmapped or r.is_secondary or r.is_supplementary:
                continue
            spliced = any(op == 3 for op, _ in r.cigartuples or ())
            for s, e in r.get_blocks():  # 0-based half-open
                a, b = max(s + 1, lo) - lo, min(e, hi) - lo + 1
                if a >= b:
                    continue
                d_all[a] += 1
                d_all[b] -= 1
                if spliced:
                    d_spl[a] += 1
                    d_spl[b] -= 1
    return np.cumsum(d_all)[:n].astype(float), np.cumsum(d_spl)[:n].astype(float)


def bin_mean(x, values, bins):
    """mean of values over the positions whose x falls in each of `bins` equal bins of [0, 1]; empty bins filled."""
    idx = np.minimum((x * bins).astype(int), bins - 1)
    n = np.bincount(idx, minlength=bins)
    s = np.bincount(idx, weights=values, minlength=bins)
    out = np.full(bins, np.nan)
    out[n > 0] = s[n > 0] / n[n > 0]
    ok = ~np.isnan(out)
    return np.interp(np.arange(bins), np.flatnonzero(ok), out[ok]) if ok.any() else np.zeros(bins)


def segments(model):
    """the canonical exons split into UTR / CDS pieces: (start, end, is_cds), 1-based inclusive."""
    out = []
    for s, e in model["exons"]:
        cuts = sorted({s, e + 1} | {c for cs, ce in model["cds"] for c in (cs, ce + 1) if s < c <= e})
        for a, b in zip(cuts, cuts[1:]):
            out.append((a, b - 1, any(cs <= a and b - 1 <= ce for cs, ce in model["cds"])))
    return out


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    sheet = read_sheet(args.sample_sheet)
    conditions = ordered(r["condition"] for r in sheet)
    donors = ordered(r["donor"] for r in sheet)
    by_key = {(r["donor"], r["condition"]): r["sample_id"] for r in sheet}
    model = canonical_model(args.gtf, args.gene)
    chrom = model["chrom"]
    lo, hi = model["exons"][0][0], model["exons"][-1][1]
    pad = int(0.02 * (hi - lo))
    glo, ghi = lo - pad, hi + pad
    pos = np.arange(glo, ghi + 1)
    shrink = IntronShrink(model["exons"])
    in_tx = (pos >= lo) & (pos <= hi)
    x_shr = shrink(pos[in_tx] + 0.5)
    x_gen = (pos - glo + 0.5) / (ghi - glo + 1)

    cov = {}
    for r in sheet:
        bam = args.bam_pattern.format(gene=args.gene, sample=r["sample_id"])
        if not os.path.exists(bam):
            print("  {}: no BAM {}, left out".format(r["sample_id"], bam), file=sys.stderr)
            continue
        a, s = coverage(bam, chrom, glo, ghi)
        cov[r["sample_id"]] = {
            "lines": (bin_mean(x_shr, a[in_tx], args.bins), bin_mean(x_shr, s[in_tx], args.bins)),
            "tracks": (bin_mean(x_gen, a, args.bins), bin_mean(x_gen, s, args.bins)),
        }
    if not cov:
        sys.exit("no BAM found for {}".format(args.gene))

    edges_shr = np.interp(np.linspace(0, 1, args.bins + 1), shrink.new, shrink.old)
    edges_gen = glo + np.linspace(0, 1, args.bins + 1) * (ghi - glo + 1)
    with gzip.open(os.path.join(args.outdir, "{}_coverage.tsv.gz".format(args.gene)), "wt") as fh:
        fh.write("gene\tx_axis\tbin\tchrom\tstart\tend\tsample_id\tcoverage_all\tcoverage_spliced\n")
        for sid, c in cov.items():
            for axis, edges in (("lines", edges_shr), ("tracks", edges_gen)):
                a, s = c[axis]
                for i in range(args.bins):
                    fh.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.3f}\t{:.3f}\n".format(
                        args.gene, "intron_shortened" if axis == "lines" else "genomic", i, chrom,
                        int(edges[i]), int(edges[i + 1]) - 1, sid, a[i], s[i]))

    lines_figure(args, model, shrink, cov, conditions, donors, by_key)
    for k, suffix in ((0, ""), (1, "_spliced")):
        tracks_figure(args, model, cov, conditions, donors, by_key, glo, ghi, k, suffix)
    print("{}: {} samples -> {}/{}_coverage_*.pdf/.png".format(args.gene, len(cov), args.outdir, args.gene))


def lines_figure(args, model, shrink, cov, conditions, donors, by_key):
    H = 5.3
    fig = plt.figure(figsize=(W, H))
    left, right = 0.95, 0.15
    panel_h, model_h = 1.7, 0.45
    xs = (np.arange(args.bins) + 0.5) / args.bins
    order = []
    for ci, c in enumerate(conditions):
        for di, d in enumerate(donors):
            if by_key.get((d, c)) in cov:
                order.append((by_key[(d, c)], c, d, CONDITION_SHADES[ci % 2][di % 4]))
    for k, (label, y0) in enumerate((("Read coverage (all)", H - 0.75 - panel_h),
                                     ("Read coverage (spliced)", H - 0.75 - 2 * panel_h - 0.25))):
        ax = fig.add_axes([left / W, y0 / H, (W - left - right) / W, panel_h / H])
        for sid, c, d, color in order:
            ax.plot(xs, cov[sid]["lines"][k], color=color, linewidth=1.1, alpha=0.85, solid_joinstyle="round")
        ymax = max(max(cov[sid]["lines"][k].max() for sid, *_ in order), 1)
        ax.set_ylim(-0.03 * ymax, 1.05 * ymax)
        ax.set_xlim(-0.01, 1.01)
        for side in ("top", "right", "bottom"):
            ax.spines[side].set_visible(False)
        ax.spines["left"].set_position(("outward", 6))
        ax.spines["left"].set_linewidth(1.2)
        ax.set_xticks([])
        ax.set_ylabel(label, fontsize=12)
        ax.yaxis.set_label_coords(-0.105, 0.5)
    # legend: one row per kit, bold kit name then its donors
    for ci, c in enumerate(conditions):
        y = (H - 0.2 - ci * 0.24) / H
        x = (W - 0.95 - 1.15 * len(donors)) / 2 + 0.1
        fig.text(x / W, y, c + ":", fontsize=11, fontweight="bold", va="center")
        x += 0.95
        for di, d in enumerate(donors):
            if by_key.get((d, c)) not in cov:
                continue
            fig.lines.append(Line2D([x / W, (x + 0.3) / W], [y, y], transform=fig.transFigure,
                                    color=CONDITION_SHADES[ci % 2][di % 4], linewidth=1.6, alpha=0.85))
            fig.text((x + 0.38) / W, y, "donor {}".format(d), fontsize=11, va="center")
            x += 1.15
    # gene model
    ax = fig.add_axes([left / W, 0.4 / H, (W - left - right) / W, model_h / H])
    draw_model(ax, model, shrink, style="lines")
    fig.text((left + (W - left - right) / 2) / W, 0.12 / H, args.gene, fontsize=13, style="italic",
             ha="center", va="center")
    save(fig, os.path.join(args.outdir, "{}_coverage_lines".format(args.gene)))


def draw_model(ax, model, to_x, style):
    ax.set_xlim(-0.01, 1.01)
    ax.set_ylim(-0.5, 0.5)
    ax.axis("off")
    xs = [to_x(s) for s, _ in model["exons"]] + [to_x(e + 1) for _, e in model["exons"]]
    ax.plot([min(xs), max(xs)], [0, 0], color="black", linewidth=1.2 if style == "lines" else 1.0, zorder=1,
            solid_capstyle="butt")
    for s, e, is_cds in segments(model):
        a, b = to_x(s), to_x(e + 1)
        if style == "lines":
            ax.add_patch(Rectangle((a, -0.42), max(b - a, 0.0015), 0.84, facecolor="#BFBFBF" if is_cds else "white",
                                   edgecolor="black", linewidth=1.2, zorder=2))
        else:
            h = 0.5 if is_cds else 0.25
            ax.add_patch(Rectangle((a, -h / 2), max(b - a, 0.001), h, facecolor="black", edgecolor="black",
                                   linewidth=0.4, zorder=2))


def tracks_figure(args, model, cov, conditions, donors, by_key, glo, ghi, k, suffix):
    order = []
    for di, d in enumerate(donors):
        for c in conditions:
            if by_key.get((d, c)) in cov:
                order.append((by_key[(d, c)], "Donor {}, {}".format(d, c), DONOR_COLORS[di % len(DONOR_COLORS)]))
    track_h, model_h = 0.40, 0.55
    H = 0.1 + track_h * len(order) + model_h + 0.35
    fig = plt.figure(figsize=(W, H))
    left, right = 1.45, 0.05
    edges = np.linspace(0, 1, args.bins + 1)
    for i, (sid, label, color) in enumerate(order):
        y0 = H - 0.1 - (i + 1) * track_h
        ax = fig.add_axes([left / W, y0 / H, (W - left - right) / W, track_h / H])
        v = cov[sid]["tracks"][k]
        vmax = max(v.max(), 1e-9)
        ax.fill_between(edges, np.append(v, v[-1]), step="post", color=color, linewidth=0)
        ax.axhline(0, color="#7F7F7F", linewidth=0.8)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, vmax * 1.12)
        ax.axis("off")
        ax.text(0.005, 0.98, "[0 - {}]".format(int(round(v.max()))), transform=ax.transAxes, fontsize=6,
                va="top", ha="left")
        fig.text((left - 0.08) / W, (y0 + track_h * 0.45) / H, label, fontsize=11, ha="right", va="center")
    ax = fig.add_axes([left / W, 0.38 / H, (W - left - right) / W, model_h / H])
    draw_model(ax, model, lambda g: (np.asarray(g, float) - glo) / (ghi - glo + 1), style="tracks")
    ax.set_xlim(0, 1)
    fig.text((left + (W - left - right) / 2) / W, 0.15 / H, args.gene, fontsize=13, style="italic",
             ha="center", va="center")
    save(fig, os.path.join(args.outdir, "{}_coverage_tracks{}".format(args.gene, suffix)))


if __name__ == "__main__":
    main()
