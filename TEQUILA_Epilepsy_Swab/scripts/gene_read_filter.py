#!/usr/bin/env python3
"""Keep only the reads that belong to one gene.

The gene BAMs of step 03 hold every read overlapping the gene's locus, which can include reads of neighbouring or
overlapping genes and reads that run mostly through introns. A read is kept for the gene when all of these hold:

  1. at least --min-exonic-fraction of its aligned bases fall on the gene's annotated exons (any transcript);
  2. it has more aligned bases on this gene's exons than on the exons of any other single gene;
  3. if it is spliced and minimap2 recorded the transcript strand (ts tag), that strand is the gene's strand
     (transcript strand = ts tag, flipped when the read aligned to the reverse strand).

Only primary, mapped alignments are considered. Outputs:
  <out_bam>            the kept reads, sorted and indexed
  one line on stdout   gene, sample, reads_in, kept, removed by each rule, spliced reads kept
"""
import argparse
import os
import re
import sys
from collections import defaultdict

import pysam

ATTR = re.compile(r'(\S+) "([^"]*)"')


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--gene", required=True)
    p.add_argument("--gtf", required=True, help="GENCODE GTF (or a part of it covering the locus)")
    p.add_argument("--bam", required=True)
    p.add_argument("--out-bam", required=True)
    p.add_argument("--sample", default="")
    p.add_argument("--min-exonic-fraction", type=float, default=0.5)
    p.add_argument("--flank", type=int, default=100000, help="other genes this close to the locus are considered")
    p.add_argument("--header", action="store_true", help="print the column names first")
    return p.parse_args()


def merge(intervals):
    out = []
    for s, e in sorted(intervals):
        if out and s <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [(s, e) for s, e in out]


def load_exons(gtf, gene, flank):
    """exons (merged, 1-based) of the gene and of every other gene near it; gene strand and span."""
    needle = 'gene_name "{}";'.format(gene)
    span = None
    with open(gtf) as fh:
        for line in fh:
            if needle in line and "\tgene\t" in line:
                f = line.split("\t")
                span = (f[0], int(f[3]), int(f[4]), f[6])
                break
    if span is None:
        sys.exit("{} not in {}".format(gene, gtf))
    chrom, lo, hi, strand = span
    exons = defaultdict(list)
    with open(gtf) as fh:
        for line in fh:
            if not line.startswith(chrom + "\t"):
                continue
            f = line.split("\t", 9)
            if f[2] != "exon":
                continue
            s, e = int(f[3]), int(f[4])
            if e < lo - flank or s > hi + flank:
                continue
            name = dict(ATTR.findall(f[8])).get("gene_name") or dict(ATTR.findall(f[8])).get("gene_id")
            exons[name].append((s, e))
    target = merge(exons.pop(gene, []))
    others = {g: merge(v) for g, v in exons.items()}
    return span, target, others


def overlap_bases(blocks, intervals):
    n = 0
    for bs, be in blocks:
        for s, e in intervals:
            if e < bs:
                continue
            if s > be:
                break
            n += min(be, e) - max(bs, s) + 1
    return n


def aligned_blocks(read):
    blocks, pos, spliced = [], read.reference_start, False
    for op, ln in read.cigartuples or ():
        if op in (0, 7, 8):
            blocks.append((pos + 1, pos + ln))
            pos += ln
        elif op == 2:
            pos += ln
        elif op == 3:
            spliced = True
            pos += ln
    return blocks, spliced


def main():
    args = parse_args()
    (chrom, lo, hi, strand), target, others = load_exons(args.gtf, args.gene, args.flank)
    counts = dict(reads_in=0, kept=0, low_exonic=0, other_gene=0, antisense_spliced=0, spliced_kept=0)
    os.makedirs(os.path.dirname(os.path.abspath(args.out_bam)), exist_ok=True)
    tmp = args.out_bam + ".unsorted.bam"
    with pysam.AlignmentFile(args.bam) as bam, pysam.AlignmentFile(tmp, "wb", template=bam) as out:
        if chrom in bam.references:
            for r in bam.fetch(chrom, lo - 1, hi):
                if r.is_unmapped or r.is_secondary or r.is_supplementary:
                    continue
                counts["reads_in"] += 1
                blocks, spliced = aligned_blocks(r)
                aligned = sum(e - s + 1 for s, e in blocks) or 1
                on_target = overlap_bases(blocks, target)
                if on_target < args.min_exonic_fraction * aligned:
                    counts["low_exonic"] += 1
                    continue
                if any(overlap_bases(blocks, iv) > on_target for iv in others.values()):
                    counts["other_gene"] += 1
                    continue
                if spliced and r.has_tag("ts"):
                    ts = r.get_tag("ts")
                    tx_strand = ts if not r.is_reverse else ("-" if ts == "+" else "+")
                    if tx_strand != strand:
                        counts["antisense_spliced"] += 1
                        continue
                counts["kept"] += 1
                counts["spliced_kept"] += int(spliced)
                out.write(r)
    pysam.sort("-o", args.out_bam, tmp)
    pysam.index(args.out_bam)
    os.remove(tmp)
    cols = ["gene", "sample_id", "reads_in", "kept", "low_exonic", "other_gene", "antisense_spliced", "spliced_kept"]
    if args.header:
        print("\t".join(cols))
    print("\t".join([args.gene, args.sample] + [str(counts[c]) for c in cols[2:]]))


if __name__ == "__main__":
    main()
