#!/usr/bin/env python3
"""One multi-panel figure per gene for the swab pilot, in the STRIPE / Nature style.

Layout and drawing rules follow the shared stripe-style figure kit (plot_stripe_style.py, Xing lab, after Wang et
al., Sci. Adv. 2026): 183 mm wide (Nature double column), 7-pt sans-serif text, bold lower-case panel letters, the
STRIPE palette, STRIPE's intron shrinking (features < 50 nt widened to 50; features > 300 nt drawn as
300 + sqrt(length)), UTRs as white boxes and CDS as grey boxes, and genomic orientation kept left to right.

Panels
  a  sashimi plot per sample: coverage plus junction arcs labelled with their read counts (arcs below the per-sample
     floor are hidden; the floor is printed), above the canonical (MANE) gene model
  b  IGV-style reads: up to --reads-per-sample reads per sample (random, seeded), drawn as aligned blocks with thin
     intron lines; reads that reach both the first and the last canonical exon (full length) are drawn darker
  c  gene-body coverage along the canonical transcript, 5' -> 3': fraction of the sample's reads covering each 1%
  d  reads per sample on canonical exons; the number above each bar is the count of
     full-length reads

Outputs: <outdir>/<GENE>_figure.pdf, .png (300 dpi) and <GENE>_figure_data.tsv (every number drawn in c and d, and the
junctions shown in a).
"""
import argparse
import csv
import math
import os
import random
import re
import sys
from collections import Counter, defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pysam
from matplotlib.patches import PathPatch, Rectangle
from matplotlib.path import Path as MPath

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fragment_report import TxMap  # noqa: E402

# STRIPE palette (plot_stripe_style.py); one colour per kit
KIT_COLORS = ["#3188BD", "#EC7B02", "#90AFB1", "#C2A5CF"]
INTRON_GREY = "#9E9E9E"
CDS_FILL = "#BFBFBF"
NATURE_W_IN = 183 / 25.4
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Nimbus Sans", "Liberation Sans", "Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 6.5, "axes.linewidth": 0.5, "axes.edgecolor": "black", "axes.labelcolor": "black",
    "xtick.color": "black", "ytick.color": "black", "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.major.size": 2, "ytick.major.size": 2, "axes.titlesize": 7, "pdf.fonttype": 42,
    "savefig.facecolor": "white",
})
ATTR = re.compile(r'(\S+) "([^"]*)"')


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--gene", required=True)
    p.add_argument("--gtf", required=True, help="GENCODE GTF")
    p.add_argument("--sample-sheet", required=True, help="sample_id, donor, condition columns (condition = kit)")
    p.add_argument("--bam-pattern", required=True, help="path with {gene} and {sample}")
    p.add_argument("--outdir", required=True)
    p.add_argument("--reads-per-sample", type=int, default=25, help="reads drawn per sample in panel b")
    p.add_argument("--min-junction-reads", type=int, default=3,
                   help="floor for arcs in panel a (also at least 5%% of the sample's largest junction)")
    p.add_argument("--seed", type=int, default=1)
    return p.parse_args()


# ---------------------------------------------------------------- annotation

def gene_model(gtf, gene):
    """gene span and the canonical transcript (MANE Select, else Ensembl_canonical): exons and CDS, 1-based."""
    span, tx_tags, exons, cds = None, {}, defaultdict(list), defaultdict(list)
    needle = 'gene_name "{}";'.format(gene)
    with open(gtf) as fh:
        for line in fh:
            if needle not in line or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            a = dict(ATTR.findall(f[8]))
            if f[2] == "gene":
                span = (f[0], int(f[3]), int(f[4]), f[6])
            elif f[2] == "transcript":
                tx_tags[a["transcript_id"]] = (a.get("transcript_name", ""), "MANE_Select" in f[8],
                                               "Ensembl_canonical" in f[8])
            elif f[2] == "exon":
                exons[a["transcript_id"]].append((int(f[3]), int(f[4])))
            elif f[2] == "CDS":
                cds[a["transcript_id"]].append((int(f[3]), int(f[4])))
    if span is None:
        sys.exit("{} not in {}".format(gene, gtf))
    ranked = sorted(tx_tags, key=lambda t: (not tx_tags[t][1], not tx_tags[t][2]))
    if not ranked or not (tx_tags[ranked[0]][1] or tx_tags[ranked[0]][2]):
        sys.exit("{}: no MANE Select / Ensembl_canonical transcript".format(gene))
    tid = ranked[0]
    chrom, lo, hi, strand = span
    return {"chrom": chrom, "lo": lo, "hi": hi, "strand": strand, "tx_id": tid, "tx_name": tx_tags[tid][0],
            "exons": sorted(exons[tid]), "cds": sorted(cds[tid])}


# ---------------------------------------------------------------- STRIPE coordinate shrinking

def rescale_feature(x):
    return 50 if 0 < x < 50 else (300 + round(math.sqrt(x)) if x > 300 else x)


class Shrink:
    """Piecewise-linear genomic -> [0, 1] map with breakpoints at feature boundaries (STRIPE)."""

    def __init__(self, breakpoints):
        bp = sorted(set(breakpoints))
        new = [0.0]
        for a, b in zip(bp, bp[1:]):
            new.append(new[-1] + rescale_feature(b - a))
        self.old, self.new = np.array(bp, float), np.array(new, float) / new[-1]

    def __call__(self, g):
        return np.interp(g, self.old, self.new)


# ---------------------------------------------------------------- reads

def blocks_of(read):
    """aligned blocks (1-based inclusive) and introns (first, last intron base) of a read."""
    blocks, introns, pos = [], [], read.reference_start
    for op, ln in read.cigartuples or ():
        if op in (0, 7, 8):
            blocks.append((pos + 1, pos + ln))
            pos += ln
        elif op == 2:
            if blocks:
                blocks[-1] = (blocks[-1][0], pos + ln)
            pos += ln
        elif op == 3:
            introns.append((pos + 1, pos + ln))
            pos += ln
    return blocks, introns


def sample_reads(bam_path, m, tmap):
    """per read: blocks, introns, transcript 5'/3' positions (or None) and whether it is full length."""
    out = []
    n_exons = len(m["exons"])
    with pysam.AlignmentFile(bam_path) as bam:
        if m["chrom"] not in bam.references:
            return out
        for r in bam.fetch(m["chrom"], m["lo"] - 1, m["hi"]):
            if r.is_unmapped or r.is_secondary or r.is_supplementary:
                continue
            blocks, introns = blocks_of(r)
            pos = []
            for bs, be in blocks:
                for s, e in tmap.exons:
                    lo, hi = max(bs, s), min(be, e)
                    if lo <= hi:
                        pos += [tmap.tx_pos(lo), tmap.tx_pos(hi)]
            ranks = tmap.first_last_exon_hit(blocks)
            out.append({"blocks": blocks, "introns": introns, "t5": min(pos) if pos else None,
                        "t3": max(pos) if pos else None,
                        "full": (0 in ranks) and (n_exons - 1 in ranks)})
    return out


def coverage(bam_path, chrom, lo, hi):
    with pysam.AlignmentFile(bam_path) as bam:
        if chrom not in bam.references:
            return np.zeros(hi - lo + 1, dtype=np.int64)
        return np.sum(np.array(bam.count_coverage(chrom, lo - 1, hi, quality_threshold=0, read_callback="nofilter"),
                               dtype=np.int64), axis=0)


# ---------------------------------------------------------------- drawing

def theme_classic(ax, x=True, y=True):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_visible(x)
    ax.spines["left"].set_visible(y)
    if not x:
        ax.set_xticks([])
    if not y:
        ax.set_yticks([])


def panel_letter(ax, letter, dx=-0.075, dy=1.0):
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=8, fontweight="bold", va="bottom", ha="right")


def draw_sashimi(ax, cov, jc, lo, sh, color, label, floor):
    xs = sh(np.arange(lo, lo + len(cov)))
    top = max(int(cov.max()), 1)
    ax.fill_between(xs, cov, 0, color=color, lw=0, step="mid")
    shown = sorted((k, n) for k, n in jc.items() if n >= floor)
    nmax = max([n for _, n in shown] or [1])
    for i, ((s, e), n) in enumerate(shown):
        x0, x1 = sh(s - 1), sh(e + 1)
        y0 = cov[min(max(s - 1 - lo, 0), len(cov) - 1)]
        y1 = cov[min(max(e + 1 - lo, 0), len(cov) - 1)]
        up = i % 2 == 0
        h = top * (0.28 + 0.12 * (i % 3))
        apex = max(y0, y1) + h if up else -h
        path = MPath([(x0, y0 if up else 0), ((x0 + x1) / 2, apex * 1.6 - (max(y0, y1) * 0.6 if up else 0)),
                      (x1, y1 if up else 0)], [MPath.MOVETO, MPath.CURVE3, MPath.CURVE3])
        ax.add_patch(PathPatch(path, facecolor="none", edgecolor=color, lw=0.4 + 1.6 * min(1.0, n / nmax), zorder=2))
        ax.text((x0 + x1) / 2, apex * 0.98, "{:,}".format(n), ha="center", va="bottom" if up else "top",
                fontsize=5, zorder=4)
    ax.set_xlim(-0.002, 1.002)
    ax.set_ylim(-top * (0.75 if len(shown) > 1 else 0.08), top * 1.75)
    ax.set_yticks([0, top])
    ax.set_yticklabels(["0", "{:,}".format(top)])
    ax.tick_params(axis="y", labelsize=5, pad=1)
    theme_classic(ax, x=False)
    ax.spines["left"].set_bounds(0, top)
    ax.text(0.998, 0.98, "{}  (arcs >= {})".format(label, floor), transform=ax.transAxes, ha="right", va="top",
            fontsize=5.5)


def draw_gene_model(ax, m, sh, lo, hi):
    ex = [(max(s, lo), min(e, hi)) for s, e in m["exons"] if e >= lo and s <= hi]
    ax.plot([sh(lo), sh(hi)], [0, 0], color="black", lw=0.6, zorder=1)
    for s, e in ex:
        ax.add_patch(Rectangle((sh(s), -0.25), sh(e + 1) - sh(s), 0.5, facecolor="white", edgecolor="black",
                               lw=0.5, zorder=2))
    for s, e in m["cds"]:
        s, e = max(s, lo), min(e, hi)
        if s <= e:
            ax.add_patch(Rectangle((sh(s), -0.25), sh(e + 1) - sh(s), 0.5, facecolor=CDS_FILL, edgecolor="black",
                                   lw=0.5, zorder=3))
    ends = ("5'", "3'") if m["strand"] == "+" else ("3'", "5'")
    ax.text(-0.004, 0, ends[0], ha="right", va="center", fontsize=6, transform=ax.get_yaxis_transform())
    ax.text(1.004, 0, ends[1], ha="left", va="center", fontsize=6, transform=ax.get_yaxis_transform())
    ax.text(0.5, -0.55, "{} {} ({}), {}:{:,}-{:,} ({} strand); introns shrunk".format(
        m["gene"], m["tx_name"], m["tx_id"], m["chrom"], lo, hi, m["strand"]),
        ha="center", va="top", fontsize=5.5, fontstyle="italic", transform=ax.transData)
    ax.set_xlim(-0.002, 1.002)
    ax.set_ylim(-0.9, 0.5)
    ax.axis("off")


def draw_reads(ax, groups, sh, colors, m):
    """groups: list of (sample, kit, reads); one row per read, groups separated by a gap."""
    y, ticks = 0, []
    for sample, kit, reads in groups:
        top = y
        for r in reads:
            dark = r["full"]
            for (a, b) in r["introns"]:
                ax.plot([sh(a), sh(b)], [y + 0.5, y + 0.5], color=INTRON_GREY, lw=0.25, zorder=1)
            for (a, b) in r["blocks"]:
                ax.add_patch(Rectangle((sh(a), y + 0.1), max(sh(b + 1) - sh(a), 0.0008), 0.8, lw=0,
                                       facecolor=colors[kit], alpha=1.0 if dark else 0.55, zorder=2))
            y += 1
        ticks.append(((top + y) / 2, sample, len(reads)))
        y += 2
    for s, e in m["exons"]:
        ax.axvspan(sh(s), sh(e + 1), color="#000000", alpha=0.035, lw=0, zorder=0)
    ax.set_xlim(-0.002, 1.002)
    ax.set_ylim(y, -1)
    ax.set_yticks([t for t, _, _ in ticks])
    ax.set_yticklabels(["{} (n={})".format(s, n) for _, s, n in ticks], fontsize=5)
    ax.tick_params(axis="y", length=0, pad=1)
    theme_classic(ax, x=False)
    ax.spines["left"].set_visible(False)


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    rng = random.Random(args.seed)
    with open(args.sample_sheet) as fh:
        sheet = list(csv.DictReader(fh, delimiter="\t"))
    kits = list(dict.fromkeys(r["condition"] for r in sheet))
    colors = {k: KIT_COLORS[i % len(KIT_COLORS)] for i, k in enumerate(kits)}

    m = gene_model(args.gtf, args.gene)
    m["gene"] = args.gene
    tmap = TxMap(m["exons"], m["strand"])
    lo, hi = m["lo"] - 100, m["hi"] + 100

    per = []
    breakpoints = [lo, hi + 1] + [x for s, e in m["exons"] for x in (s, e + 1)]
    for row in sheet:
        sample = row["sample_id"]
        bam = args.bam_pattern.format(gene=args.gene, sample=sample)
        if not os.path.exists(bam):
            print("NOTE: no BAM {}".format(bam), file=sys.stderr)
            continue
        reads = sample_reads(bam, m, tmap)
        jc = Counter(j for r in reads for j in r["introns"] if lo <= j[0] and j[1] <= hi)
        floor = max(args.min_junction_reads, int(0.05 * max(jc.values() or [0])))
        for (s, e), n in jc.items():
            if n >= floor:
                breakpoints += [s, e + 1]
        per.append({"sample": sample, "kit": row["condition"], "donor": row.get("donor", ""), "reads": reads,
                    "jc": jc, "floor": floor, "cov": coverage(bam, m["chrom"], lo, hi)})
    if not per:
        sys.exit("no BAMs found for {}".format(args.gene))
    sh = Shrink(breakpoints)

    # summaries for panels c and d
    nb = 100
    for p in per:
        on = [r for r in p["reads"] if r["t5"] is not None]
        prof = np.zeros(nb)
        for r in on:
            b0, b1 = int(r["t5"] * nb / tmap.length), int(r["t3"] * nb / tmap.length)
            prof[b0:b1 + 1] += 1
        p["profile"] = prof / len(on) if on else prof
        p["n_on"] = len(on)
        p["n_spliced"] = sum(1 for r in on if r["introns"])
        p["n_full"] = sum(1 for r in on if r["full"])

    # ---------------------------------------------------------------- layout
    n = len(per)
    track_h, model_h = 0.36, 0.42
    reads_h = 0.012 * sum(min(len(p["reads"]), args.reads_per_sample) + 2 for p in per) + 0.2
    bottom_h = 1.55
    gap = 0.32
    H = 0.25 + n * track_h + model_h + gap + reads_h + gap + bottom_h + 0.55
    W = NATURE_W_IN
    fig = plt.figure(figsize=(W, H))
    left, right = 1.05 / W, 1 - 0.15 / W
    y = H - 0.25

    def add(h, x0=left, x1=right):
        nonlocal y
        ax = fig.add_axes([x0, (y - h) / H, x1 - x0, h / H])
        y -= h
        return ax

    first = None
    for p in per:
        ax = add(track_h)
        first = first or ax
        draw_sashimi(ax, p["cov"], p["jc"], lo, sh, colors[p["kit"]], p["sample"], p["floor"])
    panel_letter(first, "a", dx=-0.06)
    first.text(-0.035, 0.5, "coverage", transform=first.transAxes, rotation=90, ha="center", va="center",
               fontsize=6)
    draw_gene_model(add(model_h), m, sh, lo, hi)
    y -= gap

    groups = []
    for p in per:
        reads = list(p["reads"])
        if len(reads) > args.reads_per_sample:
            reads = rng.sample(reads, args.reads_per_sample)
        reads.sort(key=lambda r: (r["blocks"][0][0] if r["blocks"] else 0))
        groups.append((p["sample"], p["kit"], reads))
    ax_b = add(reads_h)
    draw_reads(ax_b, groups, sh, colors, m)
    panel_letter(ax_b, "b", dx=-0.06)
    ax_b.text(0.998, -0.01, "reads (up to {} per sample, random); darker = full length (first + last exon)".format(
        args.reads_per_sample), transform=ax_b.transAxes, ha="right", va="top", fontsize=5.5)
    y -= gap

    # c and d side by side
    y_cd = y
    ax_c = fig.add_axes([left, (y_cd - bottom_h) / H, 0.40, bottom_h / H])
    ax_d = fig.add_axes([left + 0.50, (y_cd - bottom_h) / H, right - left - 0.50, bottom_h / H])
    styles = ["-", "--", ":", "-."]
    donors = list(dict.fromkeys(p["donor"] for p in per))
    xs = np.arange(nb) + 0.5
    for p in per:
        ax_c.plot(xs, p["profile"], color=colors[p["kit"]], lw=0.8, ls=styles[donors.index(p["donor"]) % 4],
                  label="{} (n={:,})".format(p["sample"], p["n_on"]))
    ax_c.set_xlim(0, 100)
    ax_c.set_ylim(0, 1.02)
    ax_c.set_xlabel("position along {} (% of {:,} nt, 5' -> 3')".format(m["tx_name"], tmap.length), fontsize=6)
    ax_c.set_ylabel("fraction of reads covering", fontsize=6)
    ax_c.legend(fontsize=4.8, frameon=False, ncol=2, loc="upper left", handlelength=1.8)
    theme_classic(ax_c)
    panel_letter(ax_c, "c", dx=-0.16)

    xb = np.arange(n)
    tot = np.array([p["n_on"] for p in per])
    ax_d.bar(xb, tot, color=[colors[p["kit"]] for p in per], width=0.7, lw=0)
    for i, p in enumerate(per):
        ax_d.text(i, tot[i], "{:,}".format(p["n_full"]), ha="center", va="bottom", fontsize=5)
    ax_d.set_xticks(xb)
    ax_d.set_xticklabels([p["sample"] for p in per], rotation=45, ha="right", fontsize=5.5)
    ax_d.set_ylabel("reads on canonical exons", fontsize=6)
    ax_d.set_ylim(0, max(1, tot.max()) * 1.18)
    ax_d.set_title("number above bar = full-length reads", fontsize=5, loc="right", pad=2)
    theme_classic(ax_d)
    panel_letter(ax_d, "d", dx=-0.16)

    base = os.path.join(args.outdir, "{}_figure".format(args.gene))
    fig.savefig(base + ".pdf")
    fig.savefig(base + ".png", dpi=300)
    plt.close(fig)

    with open(base + "_data.tsv", "w") as fh:
        fh.write("section\tsample_id\tkit\tkey\tvalue\n")
        for p in per:
            for key in ("n_on", "n_spliced", "n_full", "floor"):
                fh.write("summary\t{}\t{}\t{}\t{}\n".format(p["sample"], p["kit"], key, p[key]))
            for b, v in enumerate(p["profile"]):
                fh.write("profile\t{}\t{}\tbin_{}\t{:.4f}\n".format(p["sample"], p["kit"], b, v))
            for (s, e), c in sorted(p["jc"].items()):
                if c >= p["floor"]:
                    fh.write("junction\t{}\t{}\t{}_{}_{}\t{}\n".format(p["sample"], p["kit"], m["chrom"], s, e, c))
    print("{}: {} samples, {:.0f} x {:.0f} mm -> {}.pdf/.png".format(args.gene, n, W * 25.4, H * 25.4, base))


if __name__ == "__main__":
    main()
