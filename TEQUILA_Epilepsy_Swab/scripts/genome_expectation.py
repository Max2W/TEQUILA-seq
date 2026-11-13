#!/usr/bin/env python3
"""Read categories expected if the reads were random fragments of genomic DNA, as a reference for the RNA
composition figure (scripts/composition_figures.py).

Base fractions of the ungapped primary chromosomes (chr1-22, X, Y) of GRCh38 that are exonic (exons of any GENCODE
gene, merged), intronic (inside a gene body but not exonic) or intergenic. DNA gives no spliced reads.
Inputs: GENCODE annotation GTF (gz), UCSC hg38 chromInfo.txt.gz and gap.txt.gz
(https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/). Output (committed under public_data/):
genome_expectation.tsv.
"""
import argparse
import gzip
from collections import defaultdict

PRIMARY = ["chr{}".format(c) for c in list(range(1, 23)) + ["X", "Y"]]


def covered(iv):
    iv.sort()
    tot, cs, ce = 0, -1, -1
    for s, e in iv:
        if s > ce:
            tot += max(0, ce - cs)
            cs, ce = s, e
        else:
            ce = max(ce, e)
    return tot + max(0, ce - cs)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--gtf", required=True)
    p.add_argument("--chrom-info", required=True)
    p.add_argument("--gap", required=True)
    p.add_argument("--label", default="GENCODE v49, GRCh38")
    p.add_argument("--out", required=True)
    a = p.parse_args()
    size, gap = {}, defaultdict(int)
    for line in gzip.open(a.chrom_info, "rt"):
        c, n = line.split("\t")[:2]
        if c in PRIMARY:
            size[c] = int(n)
    for line in gzip.open(a.gap, "rt"):
        f = line.split("\t")
        if f[1] in PRIMARY:
            gap[f[1]] += int(f[3]) - int(f[2])
    exons, genes = defaultdict(list), defaultdict(list)
    for line in gzip.open(a.gtf, "rt"):
        if line.startswith("#"):
            continue
        f = line.split("\t", 5)
        if f[0] in size and f[2] in ("exon", "gene"):
            (exons if f[2] == "exon" else genes)[f[0]].append((int(f[3]) - 1, int(f[4])))
    e = sum(covered(exons[c]) for c in size)
    g = sum(covered(genes[c]) for c in size)
    u = sum(size[c] - gap[c] for c in size)
    with open(a.out, "w") as fh:
        fh.write("# Expected read categories for random genomic DNA: base fractions of the ungapped primary "
                 "chromosomes ({:.3f} Gb), {}. Made by scripts/genome_expectation.py.\n".format(u / 1e9, a.label))
        fh.write("category\tpercent\n")
        for cat, v in (("exonic, spliced", 0), ("exonic, unspliced", e / u * 100), ("intronic", (g - e) / u * 100),
                       ("intergenic", (u - g) / u * 100)):
            fh.write("{}\t{:.2f}\n".format(cat, v))
    print("exonic {:.1f}%, intronic {:.1f}%, intergenic {:.1f}%".format(e / u * 100, (g - e) / u * 100,
                                                                       (u - g) / u * 100))


if __name__ == "__main__":
    main()
