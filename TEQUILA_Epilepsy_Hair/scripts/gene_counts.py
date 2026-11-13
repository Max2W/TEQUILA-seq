#!/usr/bin/env python3
"""Gene-level read counts and the panel on-target rate from the genome-wide BAMs of the full path (07_qc_full.sh).

For every sample BAM (primary alignments only):
  * each mapped read is assigned to the gene containing its alignment midpoint (protein-coding genes first, then the
    shortest gene), or to "intergenic"
  * a read is on target when its aligned span overlaps the gene body of any probe-panel gene
Outputs (--outdir):
  gene_counts.tsv       gene, gene_type, in_panel, one read-count column per sample (plus intergenic / unmapped rows)
  panel_on_target.tsv   sample_id, reads, mapped, unmapped, panel_reads, panel_rate, reads on each --report-genes gene
  panel_genes_missing.txt  panel symbols not found as gene names in the GTF
"""
import argparse
import csv
import os
from collections import Counter, defaultdict
from multiprocessing import Pool

import pysam

BIN = 100_000
# panel symbols renamed by HGNC since the panel list was made (old -> current GENCODE gene_name)
ALIASES = {"SPATA5": "AFG2A"}
_BINS = _PANEL = _REPORT = None


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bam-pattern", required=True, help="path with {sample}")
    p.add_argument("--sample-sheet", required=True)
    p.add_argument("--gtf", required=True)
    p.add_argument("--panel", required=True, help="panel gene symbols, one per line (# comments allowed)")
    p.add_argument("--report-genes", default="", help="comma-separated genes reported in panel_on_target.tsv")
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--outdir", required=True)
    return p.parse_args()


def load_genes(gtf):
    genes = []
    with open(gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 9)
            if len(f) < 9 or f[2] != "gene":
                continue
            a = f[8]
            name = a.split('gene_name "')[1].split('"')[0] if 'gene_name "' in a else a.split('gene_id "')[1].split('"')[0]
            gtype = a.split('gene_type "')[1].split('"')[0] if 'gene_type "' in a else "NA"
            genes.append((f[0], int(f[3]) - 1, int(f[4]), name, gtype))
    return genes


def index(intervals):
    bins = defaultdict(list)
    for chrom, start, end, name, gtype in intervals:
        for b in range(start // BIN, end // BIN + 1):
            bins[(chrom, b)].append((start, end, name, gtype == "protein_coding"))
    return bins


def _init(genes, panel, report):
    global _BINS, _PANEL, _REPORT
    _BINS = index(genes)
    _PANEL = index([g for g in genes if g[3] in panel])
    _REPORT = set(report)


def gene_at(chrom, pos):
    hits = [g for g in _BINS.get((chrom, pos // BIN), ()) if g[0] <= pos < g[1]]
    if not hits:
        return "intergenic"
    hits.sort(key=lambda g: (not g[3], g[1] - g[0]))
    return hits[0][2]


def panel_hit(chrom, start, end):
    """names of panel genes whose body overlaps [start, end)."""
    out = set()
    for b in range(start // BIN, (end - 1) // BIN + 1):
        for s, e, name, _ in _PANEL.get((chrom, b), ()):
            if s < end and start < e:
                out.add(name)
    return out


def count(job):
    sample, path = job
    counts, report = Counter(), Counter()
    reads = mapped = on_panel = 0
    with pysam.AlignmentFile(path) as bam:
        for r in bam.fetch(until_eof=True):
            if r.is_secondary or r.is_supplementary:
                continue
            reads += 1
            if r.is_unmapped:
                continue
            mapped += 1
            s, e = r.reference_start, r.reference_end
            counts[gene_at(r.reference_name, (s + e) // 2)] += 1
            hits = panel_hit(r.reference_name, s, e)
            if hits:
                on_panel += 1
                for g in hits & _REPORT:
                    report[g] += 1
    return sample, reads, mapped, on_panel, counts, report


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    with open(args.sample_sheet) as fh:
        samples = [r["sample_id"] for r in csv.DictReader(fh, delimiter="\t")]
    panel = [ALIASES.get(l.strip(), l.strip()) for l in open(args.panel) if l.strip() and not l.startswith("#")]
    report = [g for g in args.report_genes.split(",") if g]
    genes = load_genes(args.gtf)
    names = {g[3] for g in genes}
    missing = [g for g in panel if g not in names]
    with open(os.path.join(args.outdir, "panel_genes_missing.txt"), "w") as fh:
        fh.write("\n".join(missing) + ("\n" if missing else ""))
    print("panel: {} genes, {} found in the GTF, {} missing (panel_genes_missing.txt)".format(
        len(panel), len(panel) - len(missing), len(missing)))
    jobs = [(s, args.bam_pattern.format(sample=s)) for s in samples if os.path.exists(args.bam_pattern.format(sample=s))]
    with Pool(min(args.threads, len(jobs)), initializer=_init, initargs=(genes, set(panel), report)) as pool:
        results = {r[0]: r for r in pool.map(count, jobs, chunksize=1)}
    order = [s for s in samples if s in results]

    gtype = {}
    for _c, _s, _e, name, t in genes:
        gtype.setdefault(name, t)
    allg = sorted(set().union(*(results[s][4] for s in order)))
    with open(os.path.join(args.outdir, "gene_counts.tsv"), "w") as fh:
        fh.write("gene\tgene_type\tin_panel\t" + "\t".join(order) + "\n")
        for g in allg:
            fh.write("{}\t{}\t{}\t{}\n".format(g, gtype.get(g, "NA"), g in set(panel),
                                               "\t".join(str(results[s][4].get(g, 0)) for s in order)))
        fh.write("unmapped\tNA\tFalse\t" + "\t".join(str(results[s][1] - results[s][2]) for s in order) + "\n")
    with open(os.path.join(args.outdir, "panel_on_target.tsv"), "w") as fh:
        fh.write("sample_id\treads\tmapped\tunmapped\tintergenic\tpanel_reads\tpanel_rate\t" +
                 "\t".join("{}_reads".format(g) for g in report) + "\n")
        for s in order:
            _, reads, mapped, on_panel, counts, rep = results[s]
            fh.write("{}\t{}\t{}\t{}\t{}\t{}\t{:.6f}\t{}\n".format(
                s, reads, mapped, reads - mapped, counts.get("intergenic", 0), on_panel,
                on_panel / reads if reads else 0, "\t".join(str(rep.get(g, 0)) for g in report)))
    print("gene counts and panel on-target -> {}".format(args.outdir))


if __name__ == "__main__":
    main()
