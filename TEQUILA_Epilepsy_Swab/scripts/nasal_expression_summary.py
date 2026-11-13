#!/usr/bin/env python3
"""Summarise gene expression in nasal airway epithelial brushings from a public RNA-seq cohort, for the next-specimen
slides: GEO GSE152004 (GALA II; Sajuthi SP et al., Nat Commun 2020, doi:10.1038/s41467-020-18781-2), inferior
turbinate brushings of 695 Puerto Rican children and young people (441 with asthma, 254 healthy controls), Illumina
whole-transcriptome RNA-seq, 250 ng RNA per library.

Reads GSE152004_695_raw_counts.txt.gz (genes x samples; download it with --download), computes counts per million
(CPM) per sample, and writes per-gene quantiles across the 695 samples for the probe-panel genes, the target genes
and cell-type markers. Output (committed under public_data/): nasal_gse152004_summary.tsv. The raw file is not kept.
"""
import argparse
import gzip
import os
import urllib.request

import numpy as np
import pandas as pd

URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE152nnn/GSE152004/suppl/GSE152004_695_raw_counts.txt.gz"
MARKERS = ["OMP", "FOXJ1", "CAPS", "SCGB1A1", "BPIFA1", "MUC5AC", "KRT5", "KRT13", "KRT4", "SPRR3", "PTPRC", "ACTB",
           "GAPDH"]
ALIASES = {"SPATA5": "AFG2A"}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--counts", required=True, help="GSE152004_695_raw_counts.txt.gz")
    p.add_argument("--download", action="store_true", help="fetch --counts from GEO first")
    p.add_argument("--panel", required=True)
    p.add_argument("--genes", required=True, help="target genes, comma-separated")
    p.add_argument("--out", required=True)
    return p.parse_args()


def main():
    a = parse_args()
    if a.download and not os.path.exists(a.counts):
        urllib.request.urlretrieve(URL, a.counts)
    with gzip.open(a.counts, "rt") as fh:
        d = pd.read_csv(fh, sep="\t", index_col=0)
    cpm = d.div(d.sum(axis=0), axis=1) * 1e6
    q = cpm.quantile([0.1, 0.25, 0.5, 0.75, 0.9], axis=1).T
    rank = q[0.5].rank(pct=True) * 100
    panel = [ALIASES.get(l.strip(), l.strip()) for l in open(a.panel) if l.strip() and not l.startswith("#")]
    targets = a.genes.split(",")
    with open(a.out, "w") as fh:
        fh.write("# GEO GSE152004 (Sajuthi et al., Nat Commun 2020): nasal airway epithelial brushings, n = {} samples, "
                 "median library {:.1f} M reads. CPM per sample; quantiles across samples. Made by "
                 "scripts/nasal_expression_summary.py.\n".format(d.shape[1], d.sum(axis=0).median() / 1e6))
        fh.write("gene\tset\tq10\tq25\tmedian\tq75\tq90\tpercentile_rank\tpercent_samples_detected\n")
        for g in list(dict.fromkeys(targets + MARKERS + panel)):
            sets = [s for s, members in (("target", targets), ("marker", MARKERS), ("panel", panel)) if g in members]
            if g not in q.index:
                fh.write("{}\t{}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n".format(g, ",".join(sets)))
                continue
            fh.write("{}\t{}\t{}\t{:.1f}\t{:.1f}\n".format(
                g, ",".join(sets), "\t".join("{:.2f}".format(v) for v in q.loc[g]), rank[g],
                (d.loc[g] > 0).mean() * 100))
    pm = q.loc[[g for g in panel if g in q.index], 0.5]
    print("panel genes in the data: {}/{}; median CPM >= 10: {:.0f}%".format(len(pm), len(panel),
                                                                           (pm >= 10).mean() * 100))
    for g in targets:
        if g in q.index:
            print("{}: median CPM {:.1f} (IQR {:.1f}-{:.1f}), percentile {:.0f}".format(g, q.loc[g, 0.5],
                                                                                      q.loc[g, 0.25], q.loc[g, 0.75],
                                                                                      rank[g]))


if __name__ == "__main__":
    main()
