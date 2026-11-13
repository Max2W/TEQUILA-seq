#!/usr/bin/env python3
"""Compare two specimen types collected from the same donors (here: plucked hair follicles vs buccal swabs).

Each project directory is a finished run of this pipeline. Sample IDs are <condition>_<donor> in both, so donors
are matched by the part after the first underscore. Measures (each drawn when present in both projects):
  RIN                             <dir>/rin.tsv
  panel on-target rate (%)        <dir>/qc/expression/panel_on_target.tsv
  gene reads per million          <dir>/results/figures/read_filter.tsv and <dir>/qc/fastq_read_counts.tsv
  median transcript coverage (%)  <dir>/results/fragmentation/filtered/fragmentation_summary.tsv
  full-length reads               same file
  patient variants with >= 20 reads   <dir>/results/figures/variants/variant_coverage.tsv
Output: specimen_comparison.{pdf,png} (one panel per measure; conditions on the x-axis, specimen by specimen; a
line per donor) and specimen_comparison.tsv (every value drawn).
"""
import argparse
import csv
import os
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
    "font.size": 8.5, "axes.labelsize": 9, "axes.linewidth": 0.8, "pdf.fonttype": 42, "savefig.facecolor": "white"})
SPECIMEN_COLORS = [["#2A9D8F", "#8FD3C9"], ["#3188BD", "#EC7B02"]]   # first project: teals; second: blue / orange
DONOR_COLORS = ["#555555", "#888888", "#AAAAAA", "#333333"]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--project", action="append", required=True, help="LABEL=directory (give two)")
    p.add_argument("--genes", required=True)
    p.add_argument("--outdir", required=True)
    return p.parse_args()


def read_tsv(path):
    if not os.path.exists(path):
        return []
    with open(path) as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"))


def split(sample_id):
    cond, _, donor = sample_id.partition("_")
    return cond, donor


def measures(d, genes):
    """measure name -> {sample_id: value}."""
    out = defaultdict(dict)
    for r in read_tsv(os.path.join(d, "rin.tsv")):
        if r.get("RIN"):
            out["RIN"][r["sample_id"]] = float(r["RIN"])
    for r in read_tsv(os.path.join(d, "qc", "expression", "panel_on_target.tsv")):
        out["Panel on-target rate (%)"][r["sample_id"]] = float(r["panel_rate"]) * 100
    total = defaultdict(int)
    for r in read_tsv(os.path.join(d, "qc", "fastq_read_counts.tsv")):
        total[r["sample_id"]] += int(r["n_reads"])
    for r in read_tsv(os.path.join(d, "results", "figures", "read_filter.tsv")):
        if r["gene"] in genes and total[r["sample_id"]]:
            out["{} reads per million".format(r["gene"])][r["sample_id"]] = int(r["kept"]) / total[r["sample_id"]] * 1e6
    for r in read_tsv(os.path.join(d, "results", "fragmentation", "filtered", "fragmentation_summary.tsv")):
        if r["gene"] in genes:
            out["{} median transcript coverage (%)".format(r["gene"])][r["sample_id"]] = \
                float(r["median_tx_coverage"]) * 100
            out["{} full-length reads".format(r["gene"])][r["sample_id"]] = int(r["full_length_reads"])
    rows = read_tsv(os.path.join(d, "results", "figures", "variants", "variant_coverage.tsv"))
    if rows:
        samples = [c[:-len("_reads")] for c in rows[0] if c.endswith("_reads") and not c.endswith("junction_reads")
                   and c not in ("all_swabs_reads", "all_samples_reads")]
        for g in genes:
            vs = [r for r in rows if r["gene"] == g]
            if vs:
                for s in samples:
                    out["{} patient variants with >= 20 reads".format(g)][s] = sum(int(r[s + "_reads"]) >= 20
                                                                                   for r in vs)
    return out


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    genes = args.genes.split(",")
    projects = []
    for item in args.project:
        label, _, d = item.partition("=")
        projects.append((label, d, measures(d, genes)))
    names = [m for m in projects[0][2] if all(m in p[2] for p in projects)]
    if not names:
        raise SystemExit("no measure is available in every project yet")
    # x positions: conditions of project 1, then of project 2
    groups = []
    for pi, (label, _, ms) in enumerate(projects):
        conds = list(dict.fromkeys(split(s)[0] for m in ms.values() for s in m))
        groups += [(pi, label, c) for c in conds]
    ncol = min(4, len(names))
    nrow = (len(names) + ncol - 1) // ncol
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.1 * ncol, 2.9 * nrow), squeeze=False)
    fig.subplots_adjust(left=0.07, right=0.98, top=0.92, bottom=0.12, wspace=0.42, hspace=0.75)
    rows = []
    for ax, name in zip(axes.flat, names):
        by_donor = defaultdict(dict)
        for gi, (pi, label, cond) in enumerate(groups):
            for s, v in projects[pi][2][name].items():
                c, donor = split(s)
                if c == cond:
                    by_donor[donor][gi] = v
                    rows.append((name, label, cond, donor, s, v))
                    ax.scatter(gi, v, s=26, color=SPECIMEN_COLORS[pi % 2][list(dict.fromkeys(
                        g[2] for g in groups if g[0] == pi)).index(cond) % 2], edgecolor="black", lw=0.4, zorder=3)
        for k, (donor, pts) in enumerate(sorted(by_donor.items())):
            xs = sorted(pts)
            ax.plot(xs, [pts[x] for x in xs], color=DONOR_COLORS[k % 4], lw=0.6, zorder=1)
            ax.text(xs[-1] + 0.12, pts[xs[-1]], donor, fontsize=7, va="center")
        ax.set_xticks(range(len(groups)))
        ax.set_xticklabels(["{}\n{}".format(label, cond) for _, label, cond in groups], fontsize=7)
        ax.set_xlim(-0.5, len(groups) - 0.2)
        ax.axvline(len([g for g in groups if g[0] == 0]) - 0.5, color="#CCCCCC", lw=0.8, ls="--")
        ax.set_title(name, fontsize=8.5)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    for ax in list(axes.flat)[len(names):]:
        ax.axis("off")
    fig.suptitle("{} vs {} in the same donors (lines join a donor's samples)".format(projects[0][0], projects[1][0]),
                 fontsize=10)
    fig.savefig(os.path.join(args.outdir, "specimen_comparison.pdf"))
    fig.savefig(os.path.join(args.outdir, "specimen_comparison.png"), dpi=250)
    plt.close(fig)
    with open(os.path.join(args.outdir, "specimen_comparison.tsv"), "w") as fh:
        fh.write("measure\tspecimen\tcondition\tdonor\tsample_id\tvalue\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")
    print("specimen comparison ({} measures) -> {}".format(len(names), args.outdir))


if __name__ == "__main__":
    main()
