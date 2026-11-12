#!/usr/bin/env python3
"""Figures in the style of the TEQUILA-seq haploinsufficient-gene progress deck, for the swab pilot.

  reads_per_sample          sequencing reads per sample (millions), sorted, slate bars
  <GENE>_gene_model         the canonical (MANE) transcript: exons, CDS, length, protein length
  <GENE>_isoform_fraction   gene abundance (CPM) above isoform fractions (%), samples sorted by the canonical
                            isoform's fraction; the structures of the isoforms shown on the right (CDS coloured,
                            UTR white)
  <GENE>_isoform_category   the same, with isoforms grouped by ORF: canonical / protein coding / premature
                            termination codon (PTC) / incomplete ORF; samples sorted by the PTC fraction
  <GENE>_event<k>           one splicing event: sashimi plots of the event region for the samples with the most
                            event reads, and the inclusion (or splice-site usage) level per sample, sorted

Rules (as in the deck)
  * isoforms shown individually: isoform fraction >= 5 % (gene reads >= 20 in that sample) and >= 3 reads, in at least
    one sample; the rest are "Other". The canonical isoform is always shown.
  * ORF: translated from the canonical start codon; PTC when the stop codon ends > 50 nt upstream of the last
    exon-exon junction; "incomplete ORF" when the isoform lacks the canonical start codon or has no stop codon.
  * splicing events, relative to the canonical transcript, from the junctions of the gene-filtered reads:
      cassette exon      an exon (annotated, in the isoform GTF, or supported by a junction pair) inside a canonical
                         intron; inclusion level = reads with an inclusion junction / (those + reads with the
                         skipping junction)
      alternative 3'/5'  a junction sharing one end with a canonical intron, its other end within 300 nt of the
      splice site        canonical site; usage = reads with it / (those + reads with the canonical junction)
    levels need >= 20 reads in the denominator. Events are ranked by the pooled reads of the alternative form; the top
    --max-events per gene with >= 20 pooled reads in the denominator and >= 3 pooled alternative reads are drawn.
Every number drawn is written to a table (<GENE>_isoform_fraction.tsv, <GENE>_isoforms.tsv,
<GENE>_splicing_events.tsv, reads_per_sample.tsv).
"""
import argparse
import csv
import os
import sys
import types
from collections import Counter, defaultdict

import numpy as np
import pandas as pd
import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from coverage_plot import bin_mean, coverage  # noqa: E402
from figure_style import IntronShrink, canonical_model, ordered, plt, read_sheet, save  # noqa: E402
from isoform_abundance_plot import gene_isoforms  # noqa: E402
from isoform_report import read_gtf_exons, strip_version  # noqa: E402
from matplotlib.patches import FancyArrow, PathPatch, Rectangle  # noqa: E402
from matplotlib.path import Path as MPath  # noqa: E402
from matplotlib.ticker import MaxNLocator  # noqa: E402

plt.rcParams.update({"font.size": 9, "axes.labelsize": 10, "xtick.labelsize": 9, "ytick.labelsize": 9,
                     "axes.linewidth": 0.8, "xtick.major.size": 3, "ytick.major.size": 3})
SLATE = "#708090"
ISO_COLORS = ["#3A8DC4", "#E2C485", "#8FD3C9", "#57B45B", "#2A9D8F", "#3D6B5B", "#C9A0DC", "#E8957B",
              "#9DB6D3", "#B5B867"]
OTHER = "#BEBEBE"
CATEGORIES = ["canonical", "Protein coding", "Premature termination codon", "Incomplete ORF"]
CAT_COLORS = {"canonical": "#3A8DC4", "Protein coding": "#66C2A5", "Premature termination codon": "#FFC0CB",
              "Incomplete ORF": "#D9D9D9"}
KIT_SASHIMI = ["#F6CF2C", "#9CCB3B", "#B8DDA8", "#F9E27D"]
MIN_DENOM, MIN_FRAC, MIN_READS = 20, 0.05, 3


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--genes", required=True, help="comma-separated")
    p.add_argument("--gtf", required=True, help="GENCODE GTF")
    p.add_argument("--fasta", required=True, help="genome FASTA (indexed)")
    p.add_argument("--sample-sheet", required=True)
    p.add_argument("--read-counts", required=True, help="fastq_read_counts.tsv")
    p.add_argument("--bam-pattern", required=True, help="gene-filtered BAMs, path with {gene} and {sample}")
    p.add_argument("--isoforms", action="append", default=[],
                   help="GENE:METHOD:matrix:isoform_gtf (repeat per gene); the isoform figures use it")
    p.add_argument("--max-events", type=int, default=3)
    p.add_argument("--outdir", required=True)
    return p.parse_args()


def label(sample_id):
    return sample_id.replace("_", " ")


def classic(ax):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def rotated_labels(ax, labels):
    ax.set_xticks(np.arange(len(labels)))
    ax.set_xticklabels(labels, rotation=-30, ha="left", rotation_mode="anchor")


# ------------------------------------------------------------------------------------------- reads per sample

def reads_figure(args, sheet, total):
    samples = sorted((r["sample_id"] for r in sheet), key=lambda s: -total[s])
    with open(os.path.join(args.outdir, "reads_per_sample.tsv"), "w") as fh:
        fh.write("sample_id\treads\n")
        for s in samples:
            fh.write("{}\t{}\n".format(s, total[s]))
    fig, ax = plt.subplots(figsize=(5.2, 3.1))
    fig.subplots_adjust(left=0.16, right=0.97, top=0.95, bottom=0.27)
    ax.bar(np.arange(len(samples)), [total[s] / 1e6 for s in samples], width=0.82, color=SLATE)
    classic(ax)
    rotated_labels(ax, [label(s) for s in samples])
    ax.set_xlim(-0.6, len(samples) - 0.4)
    ax.set_ylabel("Number of sequencing reads\n(per million)")
    save(fig, os.path.join(args.outdir, "reads_per_sample"))


# ------------------------------------------------------------------------------------------- ORF and categories

COMP = str.maketrans("ACGTN", "TGCAN")
STOPS = {"TAA", "TAG", "TGA"}


def tx_sequence(fasta, chrom, exons, strand):
    seq = "".join(fasta.fetch(chrom, s - 1, e).upper() for s, e in exons)
    return seq.translate(COMP)[::-1] if strand == "-" else seq


def tx_to_genomic(exons, strand, t0, t1):
    """transcript interval [t0, t1) -> genomic intervals (1-based inclusive)."""
    order = exons if strand == "+" else exons[::-1]
    out, pos = [], 0
    for s, e in order:
        n = e - s + 1
        a, b = max(t0, pos), min(t1, pos + n)
        if a < b:
            if strand == "+":
                out.append((s + a - pos, s + b - pos - 1))
            else:
                out.append((e - (b - pos) + 1, e - (a - pos)))
        pos += n
    return sorted(out)


def genomic_to_tx(exons, strand, g):
    order = exons if strand == "+" else exons[::-1]
    pos = 0
    for s, e in order:
        if s <= g <= e:
            return pos + (g - s if strand == "+" else e - g)
        pos += e - s + 1
    return None


def orf(fasta, chrom, exons, strand, start_codon):
    """(category without 'canonical', CDS genomic intervals) from the canonical start codon."""
    t0 = genomic_to_tx(exons, strand, start_codon)
    if t0 is None:
        return "Incomplete ORF", []
    seq = tx_sequence(fasta, chrom, exons, strand)
    stop = None
    for i in range(t0, len(seq) - 2, 3):
        if seq[i:i + 3] in STOPS:
            stop = i + 3
            break
    if stop is None:
        return "Incomplete ORF", tx_to_genomic(exons, strand, t0, len(seq) - (len(seq) - t0) % 3)
    lengths = [e - s + 1 for s, e in (exons if strand == "+" else exons[::-1])]
    last_junction = sum(lengths[:-1])
    cat = "Premature termination codon" if len(exons) > 1 and last_junction - stop > 50 else "Protein coding"
    return cat, tx_to_genomic(exons, strand, t0, stop)


# ------------------------------------------------------------------------------------------- isoform figures

def draw_structures(ax, rows, shrink, label_side, fontsize):
    """rows: (label, exons, cds, color). label_side 'above' (labels over each transcript) or 'left'."""
    n = len(rows)
    ax.set_xlim(-0.01, 1.01)
    ax.set_ylim(n - 0.4, -0.7)
    ax.axis("off")
    h = 0.32 if label_side == "above" else 0.5
    for r, (name, exons, cds, color) in enumerate(rows):
        xs = [shrink(exons[0][0]), shrink(exons[-1][1] + 1)]
        ax.plot(xs, [r, r], color="black", linewidth=0.6, zorder=1, solid_capstyle="butt")
        for s, e in exons:
            cuts = sorted({s, e + 1} | {c for cs, ce in cds for c in (cs, ce + 1) if s < c <= e})
            for a, b in zip(cuts, cuts[1:]):
                in_cds = any(cs <= a and b - 1 <= ce for cs, ce in cds)
                x0, x1 = shrink(a), shrink(b)
                ax.add_patch(Rectangle((x0, r - h / 2), max(x1 - x0, 0.0015), h, linewidth=0.5,
                                       facecolor=color if in_cds else "white", edgecolor="black", zorder=2))
        if label_side == "above":
            ax.text(0.0, r - h / 2 - 0.06, name, fontsize=fontsize, ha="left", va="bottom")
        else:
            ax.annotate(name, xy=(0, r), xycoords=("axes fraction", "data"), xytext=(-6, 0),
                        textcoords="offset points", ha="right", va="center", fontsize=fontsize,
                        annotation_clip=False)


def abundance_and_fraction(fig, box, samples, cpm, frac_rows, colors, ok, sort_note=None):
    """top: gene CPM bars; bottom: stacked fractions (frac_rows: name -> per-sample fraction, legend order)."""
    x0, y0, w, h = box
    ax1 = fig.add_axes([x0, y0 + h * 0.62, w, h * 0.33])
    ax2 = fig.add_axes([x0, y0, w, h * 0.55])
    x = np.arange(len(samples))
    ax1.bar(x, [cpm[s] for s in samples], width=0.8, color=SLATE)
    ax1.set_ylabel("Gene abundance (CPM)")
    ax1.set_xticks(x)
    ax1.set_xticklabels([])
    ax1.set_xlim(-0.6, len(samples) - 0.4)
    classic(ax1)
    bottom = np.zeros(len(samples))
    for name in reversed(list(frac_rows)):
        v = np.array([frac_rows[name][s] * 100 if ok[s] else 0.0 for s in samples])
        ax2.bar(x, v, width=0.8, bottom=bottom, color=colors[name], linewidth=0)
        bottom += v
    for i, s in enumerate(samples):
        if not ok[s]:
            ax2.text(i, 2, "<20 reads", rotation=90, ha="center", va="bottom", fontsize=7, color="#555555")
    ax2.set_ylim(-3, 103)
    ax2.set_yticks([0, 25, 50, 75, 100])
    ax2.set_ylabel("Isoform fraction (%)")
    ax2.set_xlim(-0.6, len(samples) - 0.4)
    rotated_labels(ax2, [label(s) for s in samples])
    classic(ax2)
    return ax2


def isoform_figures(args, gene, spec, sheet, total, canon, fasta, out):
    method, matrix, iso_gtf = spec
    ns = types.SimpleNamespace(gene=gene, ref_gtf=args.gtf, matrix=matrix, isoform_gtf=iso_gtf)
    sub, structures, ref = gene_isoforms(ns)
    samples_all = [r["sample_id"] for r in sheet]
    counts = pd.DataFrame({s: sub[s].astype(float) if s in sub.columns else 0.0 for s in samples_all},
                          index=sub.index)
    gene_total = counts.sum(axis=0)
    ok = {s: gene_total[s] >= MIN_DENOM for s in samples_all}
    frac = counts.div(gene_total.replace(0, np.nan), axis=1).fillna(0.0)
    cpm = {s: gene_total[s] / total[s] * 1e6 for s in samples_all}
    canon_ids = [t for t in counts.index if strip_version(t) == strip_version(canon["tx_id"])]
    canon_iso = canon_ids[0] if canon_ids else None

    # ORF / category for every isoform
    ref_cds = defaultdict(list)
    with open(args.gtf) as fh:
        needle = 'gene_name "{}";'.format(gene)
        for line in fh:
            if needle in line and "\tCDS\t" in line:
                f = line.split("\t")
                tid = f[8].split('transcript_id "')[1].split('"')[0]
                ref_cds[strip_version(tid)].append((int(f[3]), int(f[4])))
    start_codon = min(s for s, _ in canon["cds"]) if canon["strand"] == "+" else max(e for _, e in canon["cds"])
    info = {}
    for t in counts.index:
        ex = structures[t]["exons"]
        cat, cds = orf(fasta, canon["chrom"], ex, canon["strand"], start_codon)
        if strip_version(t) in ref_cds:
            cds = sorted(ref_cds[strip_version(t)])
        if t == canon_iso:
            cat = "canonical"
        info[t] = (cat, cds)

    passing = [t for t in counts.index
               if any(ok[s] and frac.loc[t, s] >= MIN_FRAC and counts.loc[t, s] >= MIN_READS for s in samples_all)]
    shown = ([canon_iso] if canon_iso else []) + sorted(
        [t for t in passing if t != canon_iso], key=lambda t: -frac.loc[t, [s for s in samples_all if ok[s]]].mean())
    other = [t for t in counts.index if t not in shown]

    with open(os.path.join(out, "{}_isoforms.tsv".format(gene)), "w") as fh:
        fh.write("gene\tmethod\tisoform\tcategory\tshown\tn_exons\tcds\t" + "\t".join(
            "{}_reads\t{}_fraction".format(s, s) for s in samples_all) + "\n")
        for t in counts.index:
            cat, cds = info[t]
            fh.write("\t".join([gene, method, t, cat, str(t in shown), str(len(structures[t]["exons"])),
                                ";".join("{}-{}".format(a, b) for a, b in cds) or "NA"] + [
                "{:.2f}\t{}".format(counts.loc[t, s], "{:.4f}".format(frac.loc[t, s]) if ok[s] else "NA")
                for s in samples_all]) + "\n")

    def label_of(t):
        return t + (" (canonical)" if t == canon_iso else "")

    # ---- per-isoform figure (slide 10 style)
    rows = {t: frac.loc[t] for t in shown}
    if other:
        rows["Other"] = frac.loc[other].sum(axis=0)
    colors = {t: ISO_COLORS[i % len(ISO_COLORS)] for i, t in enumerate(shown)}
    colors["Other"] = OTHER
    key = canon_iso if canon_iso else (shown[0] if shown else None)
    oks = [s for s in samples_all if ok[s]]
    order = sorted(oks, key=lambda s: frac.loc[key, s] if key else 0) + [s for s in samples_all if not ok[s]]
    W, H = 12.3, 5.4
    fig = plt.figure(figsize=(W, H))
    abundance_and_fraction(fig, (0.75 / W, 1.75 / H, 4.9 / W, 3.5 / H), order, cpm, rows, colors, ok)
    handles = [Rectangle((0, 0), 1, 1, facecolor=colors[t], edgecolor="none") for t in rows]
    names = [label_of(t) if t != "Other" else "Other" for t in rows]
    fig.legend(handles, names, loc="upper left", bbox_to_anchor=(1.1 / W, 0.95 / H), ncol=2, frameon=False,
               fontsize=8, handlelength=1.6, handleheight=1.1, columnspacing=1.5)
    srows = []
    if canon_iso is None:
        srows.append((canon["tx_id"] + " (canonical)", canon["exons"], canon["cds"], "white"))
    srows += [(label_of(t), structures[t]["exons"], info[t][1], colors[t]) for t in shown]
    shrink = IntronShrink([e for _, ex, _, _ in srows for e in ex], intron_fraction=0.5)
    rh = min(0.55, 4.9 / max(len(srows), 1))
    ax = fig.add_axes([6.6 / W, (H - 0.15 - rh * len(srows)) / H, 5.5 / W, rh * len(srows) / H])
    draw_structures(ax, srows, shrink, "above", 8)
    save(fig, os.path.join(out, "{}_isoform_fraction".format(gene)))

    with open(os.path.join(out, "{}_isoform_fraction.tsv".format(gene)), "w") as fh:
        fh.write("gene\tmethod\tsample_id\tgene_reads\tgene_cpm\tisoform\tfraction\n")
        for s in order:
            for t in rows:
                fh.write("{}\t{}\t{}\t{:.2f}\t{:.3f}\t{}\t{}\n".format(
                    gene, method, s, gene_total[s], cpm[s], t, "{:.4f}".format(rows[t][s]) if ok[s] else "NA"))

    # ---- category figure (slide 15 style)
    cats = [c for c in CATEGORIES if any(info[t][0] == c for t in counts.index)]
    crow = {c: frac.loc[[t for t in counts.index if info[t][0] == c]].sum(axis=0) for c in cats}
    ccolors = dict(CAT_COLORS)
    order_c = sorted(oks, key=lambda s: crow.get("Premature termination codon", pd.Series(0, index=samples_all))[s]) \
        + [s for s in samples_all if not ok[s]]
    fig = plt.figure(figsize=(W, H))
    abundance_and_fraction(fig, (0.75 / W, 1.75 / H, 4.9 / W, 3.5 / H), order_c, cpm, crow, ccolors, ok)
    names = [(canon_iso or canon["tx_id"]) + " (canonical)" if c == "canonical" else c for c in cats]
    fig.legend([Rectangle((0, 0), 1, 1, facecolor=ccolors[c], edgecolor="none") for c in cats], names,
               loc="upper left", bbox_to_anchor=(1.1 / W, 0.95 / H), ncol=1, frameon=False, fontsize=8,
               handlelength=1.6, handleheight=1.1)
    shown_c = sorted(shown, key=lambda t: (CATEGORIES.index(info[t][0]), -frac.loc[t, oks].mean() if oks else 0))
    srows = [(t + (" (canonical)" if t == canon_iso else ""), structures[t]["exons"], info[t][1],
              ccolors[info[t][0]]) for t in shown_c]
    if srows:
        shrink = IntronShrink([e for _, ex, _, _ in srows for e in ex], intron_fraction=0.5)
        rh = min(0.3, 5.0 / len(srows))
        ax = fig.add_axes([8.2 / W, (H - 0.2 - rh * len(srows)) / H, 3.9 / W, rh * len(srows) / H])
        draw_structures(ax, srows, shrink, "left", 7)
    save(fig, os.path.join(out, "{}_isoform_category".format(gene)))
    return {"method": method, "n_isoforms": len(counts), "shown": shown, "canonical": canon_iso,
            "samples_ok": oks, "categories": {c: [t for t in counts.index if info[t][0] == c] for c in cats}}


# ------------------------------------------------------------------------------------------- gene model

def gene_model_figure(gene, canon, out):
    W, H = 6.0, 1.3
    fig = plt.figure(figsize=(W, H))
    ax = fig.add_axes([0.03, 0.35, 0.94, 0.4])
    shrink = IntronShrink(canon["exons"], intron_fraction=0.4)
    draw_structures(ax, [("", canon["exons"], canon["cds"], "#3A8DC4")], shrink, "above", 8)
    ax.set_ylim(0.5, -0.5)
    tx_len = sum(e - s + 1 for s, e in canon["exons"])
    aa = sum(e - s + 1 for s, e in canon["cds"]) // 3 - 1
    fig.text(0.03, 0.88, "{} {} ({}), {}:{:,}-{:,} ({} strand)".format(
        gene, canon["tx_id"], "MANE Select", canon["chrom"], canon["exons"][0][0], canon["exons"][-1][1],
        canon["strand"]), fontsize=8, style="normal")
    fig.text(0.03, 0.08, "{} exons, {:,} nt transcript, {:,} aa protein; introns shortened".format(
        len(canon["exons"]), tx_len, aa), fontsize=8, color="#444444")
    save(fig, os.path.join(out, "{}_gene_model".format(gene)))
    return {"n_exons": len(canon["exons"]), "tx_len": tx_len, "aa": aa}


# ------------------------------------------------------------------------------------------- splicing events

def read_junctions(bam_path, chrom, lo, hi):
    """one set of introns (1-based start, end) per spliced primary read."""
    reads = []
    with pysam.AlignmentFile(bam_path) as bam:
        if chrom not in bam.references:
            return reads
        for r in bam.fetch(chrom, lo - 1, hi):
            if r.is_unmapped or r.is_secondary or r.is_supplementary:
                continue
            pos, js = r.reference_start, set()
            for op, ln in r.cigartuples or ():
                if op == 3 and ln > 0:
                    js.add((pos + 1, pos + ln))
                if op in (0, 2, 3, 7, 8):
                    pos += ln
            if js:
                reads.append(js)
    return reads


def find_events(canon, ann_exons, gencode_exons, ann_introns, pooled):
    exons = canon["exons"]
    introns = [(exons[i][1] + 1, exons[i + 1][0] - 1) for i in range(len(exons) - 1)]
    events = []
    for k, (s, e) in enumerate(introns):
        up, down = exons[k], exons[k + 1]
        # cassette exons: annotated exons inside the intron, or a supported junction pair
        cands = {x for x in ann_exons if s < x[0] and x[1] < e}
        starts = {b for (a, b), n in pooled.items() if a == s and b < e and n >= 2}
        ends = {a for (a, b), n in pooled.items() if b == e and a > s and n >= 2}
        cands |= {(j1 + 1, j2 - 1) for j1 in starts for j2 in ends if j1 + 1 <= j2 - 1}
        for xs, xe in sorted(cands):
            events.append({"type": "Cassette exon", "intron": (s, e), "up": up, "down": down, "alt": (xs, xe),
                           "inc": [(s, xs - 1), (xe + 1, e)], "skip": [(s, e)],
                           "annotated": (xs, xe) in gencode_exons, "site": (xs + xe) // 2})
        # alternative splice sites: one end shared, the other within 300 nt of the canonical site
        for (a, b), n in pooled.items():
            if n < 2 or (a, b) == (s, e):
                continue
            if a == s and b != e and abs(b - e) <= 300:
                side = "3'" if canon["strand"] == "+" else "5'"
                site = b + 1
            elif b == e and a != s and abs(a - s) <= 300:
                side = "5'" if canon["strand"] == "+" else "3'"
                site = a - 1
            else:
                continue
            events.append({"type": "Alternative {} splice site".format(side), "intron": (s, e), "up": up,
                           "down": down, "alt": (a, b), "inc": [(a, b)], "skip": [(s, e)],
                           "annotated": (a, b) in ann_introns, "site": site})
    return events


def event_levels(ev, reads_by_sample):
    out = {}
    inc, skip = set(ev["inc"]), set(ev["skip"])
    for s, reads in reads_by_sample.items():
        i = sum(1 for js in reads if js & inc)
        k = sum(1 for js in reads if js & skip and not js & inc)
        out[s] = (i, k)
    return out


def event_id(gene, ev, chrom):
    a, b = ev["alt"]
    kind = "CE" if ev["type"] == "Cassette exon" else ("A3SS" if "3'" in ev["type"] else "A5SS")
    return "{}_{}_{}:{}-{}".format(gene, kind, chrom, a, b)


def sashimi_and_bar(gene, ev, levels, sheet, cov_by_sample, junc_by_sample, chrom, strand, out_path, kit_colors):
    samples_all = [r["sample_id"] for r in sheet]
    meta = {r["sample_id"]: r for r in sheet}
    up = (max(ev["up"][0], ev["up"][1] - 249), ev["up"][1])
    down = (ev["down"][0], min(ev["down"][1], ev["down"][0] + 249))
    seg = [up, down]
    if ev["type"] == "Cassette exon":
        seg.append(ev["alt"])
    else:
        a, b = ev["alt"]
        s, e = ev["intron"]
        if a == s and b < e:
            seg.append((b + 1, e))
        elif b == e and a > s:
            seg.append((s, a - 1))
    seg = sorted(seg)
    lo, hi = up[0], down[1]
    shrink = IntronShrink(seg, intron_fraction=0.45)
    # panels: the samples with the most event reads, at most two per kit
    by_kit = defaultdict(list)
    for s in sorted(samples_all, key=lambda s: -sum(levels[s])):
        by_kit[meta[s]["condition"]].append(s)
    conds = ordered(r["condition"] for r in sheet)
    panels = []
    for c in conds:
        good = [s for s in by_kit[c] if sum(levels[s]) >= MIN_DENOM][:2]
        panels += good or [s for s in by_kit[c] if sum(levels[s]) > 0][:1]
    psi = {s: (levels[s][0] / sum(levels[s]) * 100 if sum(levels[s]) >= MIN_DENOM else None) for s in samples_all}

    W, H = 12.3, 5.2
    fig = plt.figure(figsize=(W, H))
    ph = (H - 1.25) / max(len(panels), 1)
    pos = np.arange(lo, hi + 1)
    xg = shrink(pos + 0.5)
    for i, s in enumerate(panels):
        color = kit_colors[conds.index(meta[s]["condition"]) % len(kit_colors)]
        ax = fig.add_axes([0.85 / W, (H - 0.15 - (i + 1) * ph + 0.12) / H, 4.6 / W, (ph - 0.2) / H])
        cv = cov_by_sample[s][lo - cov_by_sample[s + "_lo"]: hi - cov_by_sample[s + "_lo"] + 1]
        ys = bin_mean(xg, cv, 600)
        edges = np.linspace(0, 1, 601)
        ax.fill_between(edges, np.append(ys, ys[-1]), step="post", color=color, linewidth=0)
        ymax = max(ys.max(), 1)
        js = [(j, n) for j, n in junc_by_sample[s].items() if j in set(ev["inc"]) | set(ev["skip"])]
        nmax = max([n for _, n in js] + [1])
        for k, ((a, b), n) in enumerate(sorted(js)):
            x1, x2 = float(shrink(a)), float(shrink(b + 1))
            above = k % 2 == 0
            y1 = ys[max(min(int(x1 * 600) - 1, 599), 0)] if above else 0
            y2 = ys[min(int(x2 * 600), 599)] if above else 0
            peak = (ymax * 1.25) if above else -ymax * 0.35
            path = MPath([(x1, y1), (x1, peak), (x2, peak), (x2, y2)],
                         [MPath.MOVETO, MPath.CURVE4, MPath.CURVE4, MPath.CURVE4])
            ax.add_patch(PathPatch(path, facecolor="none", edgecolor=color, linewidth=0.6 + 2.4 * n / nmax))
            ax.text((x1 + x2) / 2, peak * 0.76 if above else peak * 0.78, str(n), ha="center", va="center",
                    fontsize=7, bbox=dict(boxstyle="round,pad=0.15", facecolor="white", edgecolor="none"))
        ax.set_xlim(-0.01, 1.01)
        ax.set_ylim(-ymax * 0.45, ymax * 1.3)
        ticks = [t for t in MaxNLocator(nbins=3, steps=[1, 2, 2.5, 5, 10]).tick_values(0, ymax) if 0 <= t <= ymax]
        ax.set_yticks(ticks if len(ticks) > 1 else [0, ymax])
        for side in ("top", "right", "bottom"):
            ax.spines[side].set_visible(False)
        ax.spines["left"].set_bounds(0, ymax)
        ax.set_xticks([])
        ax.text(1.0, 1.0, label(s), transform=ax.transAxes, ha="right", va="top", fontsize=9)
    # model: canonical exons in black, red arrow at the event
    axm = fig.add_axes([0.85 / W, 0.35 / H, 4.6 / W, 0.55 / H])
    axm.set_xlim(-0.01, 1.01)
    axm.set_ylim(-1, 1)
    axm.axis("off")
    axm.plot([shrink(up[0]), shrink(down[1] + 1)], [0.45, 0.45], color="black", linewidth=0.8)
    for s_, e_ in (up, down):
        axm.add_patch(Rectangle((shrink(s_), 0.2), shrink(e_ + 1) - shrink(s_), 0.5, color="black"))
    xa = float(shrink(ev["site"]))
    axm.add_patch(FancyArrow(xa, -1.0, 0, 0.85, width=0.012, head_width=0.035, head_length=0.4,
                             length_includes_head=True, color="#C00000", transform=axm.transData))
    # bar: level per sample, sorted
    have = sorted([s for s in samples_all if psi[s] is not None], key=lambda s: -psi[s])
    low = [s for s in samples_all if psi[s] is None]
    order = have + low
    ax = fig.add_axes([6.5 / W, 1.4 / H, 5.5 / W, 2.6 / H])
    x = np.arange(len(order))
    ax.bar(x[:len(have)], [psi[s] for s in have], width=0.8, color=SLATE)
    for i, s in enumerate(order):
        n = sum(levels[s])
        if psi[s] is None:
            ax.text(i, 2, "n = {}".format(n), rotation=90, ha="center", va="bottom", fontsize=7, color="#555555")
        else:
            ax.text(i, psi[s] + 2, "n = {}".format(n), ha="center", va="bottom", fontsize=6.5, color="#555555")
    ax.set_ylim(-3, 105)
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_xlim(-0.6, len(order) - 0.4)
    rotated_labels(ax, [label(s) for s in order])
    classic(ax)
    is_ce = ev["type"] == "Cassette exon"
    ax.set_ylabel("Exon inclusion level (%)" if is_ce else "{} usage (%)".format(ev["type"]))
    what = "{} {}".format(gene, "cassette exon" if is_ce else ev["type"].lower())
    fig.text(9.25 / W, 4.75 / H, "{} levels for {} were calculated directly\nfrom the long-read alignments "
             "(n = reads in the denominator)*".format("Inclusion" if is_ce else "Usage", what),
             ha="center", va="center", fontsize=9,
             bbox=dict(boxstyle="square,pad=0.4", facecolor="white", edgecolor="black", linewidth=0.8))
    save(fig, out_path)
    return psi, panels


# ------------------------------------------------------------------------------------------- main

def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    sheet = read_sheet(args.sample_sheet)
    samples_all = [r["sample_id"] for r in sheet]
    total = defaultdict(int)
    with open(args.read_counts) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            total[row["sample_id"]] += int(row["n_reads"])
    missing = [s for s in samples_all if not total[s]]
    if missing:
        sys.exit("no read count for {} in {}".format(", ".join(missing), args.read_counts))
    reads_figure(args, sheet, total)
    fasta = pysam.FastaFile(args.fasta)
    specs = {}
    for item in args.isoforms:
        g, method, matrix, gtf = item.split(":", 3) if item.count(":") == 3 else (None, None, None, None)
        if g is None:
            sys.exit("--isoforms takes GENE:METHOD:matrix:isoform_gtf, got {}".format(item))
        specs.setdefault(g, (method, matrix, gtf))
    conds = ordered(r["condition"] for r in sheet)

    for gene in args.genes.split(","):
        canon = canonical_model(args.gtf, gene)
        chrom, strand = canon["chrom"], canon["strand"]
        gene_model_figure(gene, canon, args.outdir)
        if gene in specs:
            r = isoform_figures(args, gene, specs[gene], sheet, total, canon, fasta, args.outdir)
            print("{}: isoform figures from {} ({} isoforms, {} shown, {} samples with >= 20 reads)".format(
                gene, r["method"], r["n_isoforms"], len(r["shown"]), len(r["samples_ok"])))
        else:
            print("{}: no isoform quantification given, isoform figures skipped".format(gene))

        # splicing events
        ref = read_gtf_exons(args.gtf, keep=lambda f, a: a.get("gene_name") == gene)
        gencode_exons = {e for t in ref.values() for e in t["exons"]}
        ann_exons = set(gencode_exons)
        ann_introns = {(t["exons"][i][1] + 1, t["exons"][i + 1][0] - 1) for t in ref.values()
                       for i in range(len(t["exons"]) - 1)}
        if gene in specs:
            iso_all = read_gtf_exons(specs[gene][2])
            ann_exons |= {e for t in iso_all.values() if t["chrom"] == chrom and t["strand"] == strand
                          for e in t["exons"]}
        lo, hi = canon["exons"][0][0], canon["exons"][-1][1]
        reads_by_sample, junc_by_sample, cov_by_sample = {}, {}, {}
        for s in samples_all:
            bam = args.bam_pattern.format(gene=gene, sample=s)
            reads_by_sample[s] = read_junctions(bam, chrom, lo, hi) if os.path.exists(bam) else []
            junc_by_sample[s] = Counter(j for js in reads_by_sample[s] for j in js)
        pooled = sum(junc_by_sample.values(), Counter())
        events = find_events(canon, ann_exons, gencode_exons, ann_introns, pooled)
        rows, ranked = [], []
        for ev in events:
            lv = event_levels(ev, reads_by_sample)
            alt_pool = sum(v[0] for v in lv.values())
            den_pool = sum(sum(v) for v in lv.values())
            ev["id"] = event_id(gene, ev, chrom)
            for s in samples_all:
                i, k = lv[s]
                rows.append([gene, ev["id"], ev["type"], chrom, ev["alt"][0], ev["alt"][1], ev["annotated"], s,
                             i, k, i + k, "{:.2f}".format(i / (i + k) * 100) if i + k >= MIN_DENOM else "NA"])
            if den_pool >= MIN_DENOM and alt_pool >= MIN_READS:
                ranked.append((alt_pool, den_pool, ev, lv))
        ranked.sort(key=lambda r: (-r[0], -r[1]))
        with open(os.path.join(args.outdir, "{}_splicing_events.tsv".format(gene)), "w") as fh:
            fh.write("gene\tevent_id\ttype\tchrom\tstart\tend\tin_gencode\tsample_id\talt_reads\tcanonical_reads\t"
                     "denominator\tlevel_percent\n")
            for r in rows:
                fh.write("\t".join(map(str, r)) + "\n")
        with open(os.path.join(args.outdir, "{}_events_drawn.tsv".format(gene)), "w") as fh:
            fh.write("rank\tevent_id\ttype\tin_gencode\tpooled_alt_reads\tpooled_denominator\tfigure\n")
            for k, (alt_pool, den_pool, ev, lv) in enumerate(ranked[:args.max_events], 1):
                for s in samples_all:
                    if s not in cov_by_sample:
                        bam = args.bam_pattern.format(gene=gene, sample=s)
                        a, _ = coverage(bam, chrom, lo, hi) if os.path.exists(bam) else (np.zeros(hi - lo + 1), 0)
                        cov_by_sample[s], cov_by_sample[s + "_lo"] = a, lo
                name = "{}_event{}".format(gene, k)
                sashimi_and_bar(gene, ev, lv, sheet, cov_by_sample, junc_by_sample, chrom, strand,
                                os.path.join(args.outdir, name), KIT_SASHIMI)
                fh.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(k, ev["id"], ev["type"], ev["annotated"], alt_pool,
                                                               den_pool, name))
        print("{}: {} splicing events found, {} with >= {} pooled reads, {} drawn".format(
            gene, len(events), len(ranked), MIN_DENOM, min(len(ranked), args.max_events)))


if __name__ == "__main__":
    main()
