#!/usr/bin/env python3
"""Expression of every probe-panel gene in brain and in accessible specimens, from public data, in one table.

Sources:
  GTEx v8 gene median TPM (GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz): brain and the GTEx tissue
    closest to each accessible specimen; each tissue rescaled so that its protein-coding genes sum to 1e6
  GEO GSE152004 (Sajuthi et al., Nat Commun 2020): raw counts of 695 nasal airway epithelial brushings, turned into TPM
    with the length of each gene's merged exons in the GENCODE annotation (--gtf), scaled over protein-coding genes
    (as in scripts/specimen_expression.py); median across samples
Output (committed under public_data/): panel_expression_atlas.tsv, one row per panel gene with the median TPM in
each column, plus whether the gene is a highlighted target.
"""
import argparse
import gzip
from collections import defaultdict

import numpy as np
import pandas as pd

COLUMNS = [  # output column -> GTEx tissue
    ("Brain cortex", "Brain - Cortex"), ("Brain cerebellum", "Brain - Cerebellum"), ("Lung", "Lung"),
    ("Esophagus mucosa", "Esophagus - Mucosa"), ("Minor salivary gland", "Minor Salivary Gland"),
    ("Skin", "Skin - Not Sun Exposed (Suprapubic)"), ("Cultured fibroblasts", "Cells - Cultured fibroblasts"),
    ("EBV lymphocytes", "Cells - EBV-transformed lymphocytes"), ("Whole blood", "Whole Blood")]
# current HGNC symbol -> older symbols used by GTEx v8 or GSE152004
OLD = {"AFG2A": ["SPATA5"], "AARS1": ["AARS"], "EPRS1": ["EPRS"], "QARS1": ["QARS"], "H1-4": ["HIST1H1E"],
       "H3-3A": ["H3F3A"], "NEXMIF": ["KIAA2022"], "LINS1": ["LINS"], "PLPBP": ["PROSC"], "KMT5B": ["SUV420H1"],
       "NSD2": ["WHSC1"], "SPATA5": ["AFG2A"]}


def lookup(table, gene):
    for g in [gene] + OLD.get(gene, []):
        if g in table:
            return g
    return None


CODING = set()


def exon_lengths(gtf):
    iv = defaultdict(list)
    with gzip.open(gtf, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if f[2] != "exon":
                continue
            name = f[8].split('gene_name "', 1)[1].split('"', 1)[0]
            iv[(name, f[0])].append((int(f[3]) - 1, int(f[4])))
            if 'gene_type "protein_coding"' in f[8]:
                CODING.add(name)
    length = {}
    for (name, _), v in iv.items():
        v.sort()
        tot, cs, ce = 0, -1, -1
        for s, e in v:
            if s > ce:
                tot += max(0, ce - cs)
                cs, ce = s, e
            else:
                ce = max(ce, e)
        tot += max(0, ce - cs)
        length[name] = max(length.get(name, 0), tot)
    return length


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--gtex-gct", required=True)
    p.add_argument("--nasal-counts", required=True)
    p.add_argument("--gtf", required=True)
    p.add_argument("--panel", required=True)
    p.add_argument("--genes", required=True, help="highlighted genes, comma-separated")
    p.add_argument("--out", required=True)
    a = p.parse_args()
    panel = [l.strip() for l in open(a.panel) if l.strip() and not l.startswith("#")]
    length = exon_lengths(a.gtf)
    gtex = pd.read_csv(a.gtex_gct, sep="\t", skiprows=2)
    gtex = gtex.groupby("Description").max(numeric_only=True)
    # rescale each tissue's median TPM over protein-coding genes, as for the nasal and specimen data
    gc = gtex.index.isin(CODING)
    gtex = gtex * (1e6 / gtex[gc].sum(axis=0))
    counts = pd.read_csv(a.nasal_counts, sep="\t", index_col=0)
    lens = pd.Series({g: length.get(g, length.get(next((n for n, olds in OLD.items() if g in olds), ""), np.nan))
                      for g in counts.index})
    coding = pd.Series([g in CODING or any(g in olds for n, olds in OLD.items() if n in CODING)
                        for g in counts.index], index=counts.index)
    keep = lens.notna() & (lens > 0) & coding
    rpk = counts[keep].div(lens[keep] / 1000, axis=0)
    tpm = rpk.div(rpk.sum(axis=0), axis=1) * 1e6
    nasal_med, nasal_det = tpm.median(axis=1), (counts > 0).mean(axis=1) * 100
    print("nasal TPM from {} of {} genes with a GENCODE exon length".format(int(keep.sum()), len(counts)))
    targets = a.genes.split(",")
    with open(a.out, "w") as fh:
        fh.write("# Median TPM per panel gene, scaled over protein-coding genes. GTEx v8 gene median TPM; nasal brushings: GEO GSE152004 (n = {}), TPM "
                 "from raw counts and GENCODE merged-exon lengths, over protein-coding genes. NA: gene not found. Made by "
                 "scripts/panel_atlas.py.\n".format(counts.shape[1]))
        fh.write("gene\ttarget\t" + "\t".join(c for c, _ in COLUMNS) + "\tNasal brushing\tnasal_percent_detected\n")
        for g in panel:
            gg = lookup(gtex.index, g)
            vals = ["{:.2f}".format(gtex.loc[gg, t]) if gg is not None else "NA" for _, t in COLUMNS]
            ng = lookup(nasal_med.index, g)
            nv = "{:.2f}".format(nasal_med[ng]) if ng is not None and ng in nasal_med.index else "NA"
            nd = "{:.1f}".format(nasal_det[ng]) if ng is not None else "NA"
            fh.write("{}\t{}\t{}\t{}\t{}\n".format(g, g in targets, "\t".join(vals), nv, nd))
    print("written {}".format(a.out))


if __name__ == "__main__":
    main()
