"""Shared drawing rules for the sample figures (isoform_abundance_plot.py, coverage_plot.py).

The style follows the isoform and coverage figures of the TEQUILA-seq / SMG1i (NMD inhibition) studies: ggplot-like
panels with black outlines, Arial-type text, transcripts drawn on a common axis whose introns are shortened.
"""
import re
import sys
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ATTR = re.compile(r'(\S+) "([^"]*)"')

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "Nimbus Sans", "DejaVu Sans"],
    "font.size": 11, "axes.linewidth": 0.8, "axes.edgecolor": "black", "axes.labelcolor": "black",
    "xtick.color": "black", "ytick.color": "black", "xtick.major.width": 0.8, "ytick.major.width": 0.8,
    "xtick.major.size": 4, "ytick.major.size": 4, "xtick.labelsize": 11, "ytick.labelsize": 11,
    "axes.labelsize": 13, "pdf.fonttype": 42, "svg.fonttype": "none", "savefig.facecolor": "white",
})

# kits: reference condition first (blues), the other second (reds), one shade per donor (dark -> light)
CONDITION_SHADES = [
    ["#2166AC", "#4393C3", "#74ADD1", "#9ECAE1"],
    ["#B2182B", "#D6604D", "#F4A582", "#F7B9A0"],
]
# one colour per donor (IGV-style tracks)
DONOR_COLORS = ["#BE1E5A", "#EFA06C", "#8DB040", "#4A88B4", "#9B5BA5", "#D4B42C", "#DE85AE"]


def save(fig, path_no_ext):
    fig.savefig(path_no_ext + ".pdf")
    fig.savefig(path_no_ext + ".png", dpi=300)
    plt.close(fig)


def read_sheet(path):
    import csv
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def ordered(values):
    return list(dict.fromkeys(values))


def canonical_model(gtf, gene):
    """the gene's MANE Select transcript (else Ensembl_canonical): id, chrom, strand, exons and CDS (1-based)."""
    needle = 'gene_name "{}";'.format(gene)
    tags, exons, cds, where = {}, defaultdict(list), defaultdict(list), {}
    with open(gtf) as fh:
        for line in fh:
            if needle not in line or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[2] not in ("transcript", "exon", "CDS"):
                continue
            tid = dict(ATTR.findall(f[8]))["transcript_id"]
            if f[2] == "transcript":
                tags[tid] = ("MANE_Select" in f[8], "Ensembl_canonical" in f[8])
                where[tid] = (f[0], f[6])
            elif f[2] == "exon":
                exons[tid].append((int(f[3]), int(f[4])))
            else:
                cds[tid].append((int(f[3]), int(f[4])))
    ranked = sorted((t for t in tags if any(tags[t])), key=lambda t: (not tags[t][0], not tags[t][1]))
    if not ranked:
        sys.exit("{}: no MANE Select / Ensembl_canonical transcript in {}".format(gene, gtf))
    tid = ranked[0]
    return {"tx_id": tid, "chrom": where[tid][0], "strand": where[tid][1],
            "exons": sorted(exons[tid]), "cds": sorted(cds[tid])}


def merge(intervals):
    out = []
    for s, e in sorted(intervals):
        if out and s <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [tuple(x) for x in out]


class IntronShrink:
    """Genomic position -> [0, 1]. Exonic stretches keep their length; each gap between exons is drawn with a width
    proportional to the square root of its length, all gaps together taking --intron-fraction of the axis (never
    wider than the gap itself)."""

    def __init__(self, exons, intron_fraction=0.4):
        ex = merge(exons)
        old, new = [ex[0][0]], [0.0]
        exon_total = sum(e - s + 1 for s, e in ex)
        gaps = [ex[i + 1][0] - ex[i][1] - 1 for i in range(len(ex) - 1)]
        weight = sum(np.sqrt(g) for g in gaps if g > 0)
        scale = (intron_fraction / (1 - intron_fraction) * exon_total / weight) if weight else 0.0
        for i, (s, e) in enumerate(ex):
            old.append(e + 1)
            new.append(new[-1] + (e - s + 1))
            if i < len(gaps) and gaps[i] > 0:
                old.append(ex[i + 1][0])
                new.append(new[-1] + min(gaps[i], scale * np.sqrt(gaps[i])))
        self.old, self.new = np.array(old, float), np.array(new, float) / new[-1]
        self.lo, self.hi = ex[0][0], ex[-1][1]

    def __call__(self, g):
        return np.interp(g, self.old, self.new)
