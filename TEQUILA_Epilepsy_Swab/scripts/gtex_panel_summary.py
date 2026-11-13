#!/usr/bin/env python3
"""Share of the probe-panel genes expressed in each GTEx v8 tissue, for the accessible-tissue slides.

Input: GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz (median TPM per tissue; GTEx Portal
downloads, https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/), the panel list and the target genes.
Output (committed under public_data/): gtex_panel_expression.tsv with, per tissue, the panel genes found, the share
with median TPM >= 1 and >= 10, and the median TPM of each target gene.
"""
import argparse
import gzip

import numpy as np

ALIASES = {"AFG2A": "SPATA5", "AARS1": "AARS", "EPRS1": "EPRS", "QARS1": "QARS", "H1-4": "HIST1H1E",
           "H3-3A": "H3F3A", "NEXMIF": "KIAA2022", "LINS1": "LINS", "PLPBP": "PROSC", "KMT5B": "SUV420H1",
           "NSD2": "WHSC1"}


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--gct", required=True)
    p.add_argument("--panel", required=True)
    p.add_argument("--genes", required=True)
    p.add_argument("--out", required=True)
    a = p.parse_args()
    with gzip.open(a.gct, "rt") as fh:
        fh.readline(), fh.readline()
        tissues = fh.readline().rstrip("\n").split("\t")[2:]
        tpm = {}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            v = np.array(f[2:], dtype=float)
            if f[1] not in tpm or v.sum() > tpm[f[1]].sum():   # duplicated symbols: keep the most expressed
                tpm[f[1]] = v
    panel = [l.strip() for l in open(a.panel) if l.strip() and not l.startswith("#")]
    found = [ALIASES.get(g, g) if g not in tpm else g for g in panel]
    found = [g for g in found if g in tpm]
    m = np.array([tpm[g] for g in found])
    targets = a.genes.split(",")
    with open(a.out, "w") as fh:
        fh.write("# GTEx v8 median TPM per tissue (gene_median_tpm.gct). {} of {} panel genes found (older symbols "
                 "mapped). Made by scripts/gtex_panel_summary.py.\n".format(len(found), len(panel)))
        fh.write("tissue\tpanel_genes\tpercent_ge1_tpm\tpercent_ge10_tpm\t" +
                 "\t".join("{}_tpm".format(g) for g in targets) + "\n")
        for j, t in enumerate(tissues):
            fh.write("{}\t{}\t{:.1f}\t{:.1f}\t{}\n".format(t, len(found), (m[:, j] >= 1).mean() * 100,
                                                          (m[:, j] >= 10).mean() * 100,
                                                          "\t".join("{:.1f}".format(tpm[g][j]) for g in targets)))
    print("{} panel genes found; written {}".format(len(found), a.out))


if __name__ == "__main__":
    main()
