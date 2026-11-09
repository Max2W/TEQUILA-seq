#!/usr/bin/env python3
"""Alignment QC for the epilepsy buccal-swab TEQUILA-seq run.

For every flowcell input BAM and every merged per-sample BAM:
  reads (primary records), mapped, mapping rate, reads overlapping each target gene and the
  on-target rate, median read length, median read quality (mean per-read error, Phred), and the
  genes that collect the most reads (a read is assigned to the gene containing its alignment
  midpoint). The top-gene table shows directly what the probe panel captured, and the per-flowcell
  rows show whether flowcells 3B and 3C carry the same library.

Outputs (in --outdir): per_input_qc.tsv, per_sample_qc.tsv, top_genes.tsv, target_genes.bed,
                       qc_on_target.pdf, qc_read_length.pdf
"""
import argparse
import csv
import os
import re
from collections import Counter, defaultdict
from multiprocessing import Pool

import numpy as np
import pysam

BIN = 100_000
MAX_SAMPLED = 200_000


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--input-index', required=True, help='sample_id, input_index, flowcell, path')
    parser.add_argument('--bam-dir', required=True, help='has <sample>.bam and per_input/<sample>.<fc>.bam')
    parser.add_argument('--sample-sheet', required=True)
    parser.add_argument('--gtf', required=True)
    parser.add_argument('--target-genes', required=True, help='comma separated gene names')
    parser.add_argument('--top-n', type=int, default=20)
    parser.add_argument('--threads', type=int, default=4)
    parser.add_argument('--outdir', required=True)
    return parser.parse_args()


def load_genes(gtf):
    pat_name = re.compile(r'gene_name "([^"]+)"')
    pat_type = re.compile(r'gene_type "([^"]+)"')
    genes = list()
    with open(gtf) as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            f = line.split('\t', 9)
            if len(f) < 9 or f[2] != 'gene':
                continue
            name = pat_name.search(f[8])
            gtype = pat_type.search(f[8])
            genes.append((f[0], int(f[3]) - 1, int(f[4]), f[6], name.group(1) if name else 'NA',
                          gtype.group(1) if gtype else 'NA'))
    return genes


_BINS = None
_TARGETS = None


def _init(genes, targets):
    global _BINS, _TARGETS
    bins = defaultdict(list)
    for chrom, start, end, _strand, name, gtype in genes:
        for b in range(start // BIN, end // BIN + 1):
            bins[(chrom, b)].append((start, end, name, gtype == 'protein_coding'))
    _BINS = bins
    _TARGETS = targets


def gene_at(chrom, pos):
    hits = [g for g in _BINS.get((chrom, pos // BIN), ()) if g[0] <= pos < g[1]]
    if not hits:
        return 'intergenic'
    hits.sort(key=lambda g: (not g[3], g[1] - g[0]))  # protein-coding first, then shortest
    return hits[0][2]


def qc_bam(job):
    label, sample, flowcell, path = job
    n_reads = n_mapped = 0
    gene_counts = Counter()
    target_ids = {name: set() for name, _region in _TARGETS}
    lengths, quals = list(), list()
    with pysam.AlignmentFile(path, 'rb') as bam:
        for rec in bam.fetch(until_eof=True):
            if rec.is_secondary or rec.is_supplementary:
                continue
            n_reads += 1
            if len(lengths) < MAX_SAMPLED:
                lengths.append(rec.query_length)
                q = rec.query_qualities
                if q is not None and len(q):
                    err = np.power(10.0, -np.frombuffer(bytes(q), dtype=np.uint8) / 10.0).mean()
                    quals.append(-10.0 * np.log10(err))
            if rec.is_unmapped:
                continue
            n_mapped += 1
            mid = (rec.reference_start + rec.reference_end) // 2
            gene_counts[gene_at(rec.reference_name, mid)] += 1
        for name, (chrom, start, end) in _TARGETS:
            for rec in bam.fetch(chrom, start, end):
                if not (rec.is_secondary or rec.is_supplementary or rec.is_unmapped):
                    target_ids[name].add(rec.query_name)
    return {
        'label': label, 'sample_id': sample, 'flowcell': flowcell, 'reads': n_reads, 'mapped': n_mapped,
        'targets': {name: len(ids) for name, ids in target_ids.items()},
        'median_length': float(np.median(lengths)) if lengths else 0.0,
        'median_q': float(np.median(quals)) if quals else 0.0,
        'lengths': lengths, 'top': gene_counts.most_common(),
    }


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    genes = load_genes(args.gtf)
    wanted = [g for g in args.target_genes.split(',') if g]
    targets = list()
    with open(os.path.join(args.outdir, 'target_genes.bed'), 'w') as bed:
        for name in wanted:
            rec = [g for g in genes if g[4] == name]
            if not rec:
                raise SystemExit('{} not in {}'.format(name, args.gtf))
            chrom, start, end, strand = rec[0][:4]
            targets.append((name, (chrom, start, end)))
            bed.write('{}\t{}\t{}\t{}\t0\t{}\n'.format(chrom, start, end, name, strand))

    with open(args.sample_sheet) as handle:
        sheet = {row['sample_id']: row for row in csv.DictReader(handle, delimiter='\t')}
    jobs = list()
    with open(args.input_index) as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            path = os.path.join(args.bam_dir, 'per_input', '{}.{}.bam'.format(row['sample_id'], row['flowcell']))
            if os.path.exists(path):
                jobs.append(('input', row['sample_id'], row['flowcell'], path))
    for sample in sheet:
        path = os.path.join(args.bam_dir, '{}.bam'.format(sample))
        if os.path.exists(path):
            jobs.append(('sample', sample, 'all', path))

    with Pool(args.threads, initializer=_init, initargs=(genes, targets)) as pool:
        results = pool.map(qc_bam, jobs, chunksize=1)

    cols = ['sample_id', 'flowcell', 'cell_type', 'condition', 'reads', 'mapped', 'mapping_rate']
    for name in wanted:
        cols += ['{}_reads'.format(name), '{}_rate'.format(name)]
    cols += ['median_read_length', 'median_read_q']
    for level in ('input', 'sample'):
        out = os.path.join(args.outdir, 'per_{}_qc.tsv'.format(level))
        with open(out, 'w') as handle:
            handle.write('\t'.join(cols) + '\n')
            for r in results:
                if r['label'] != level:
                    continue
                meta = sheet.get(r['sample_id'], {})
                vals = [r['sample_id'], r['flowcell'], meta.get('cell_type', ''), meta.get('condition', ''),
                        r['reads'], r['mapped'], '{:.4f}'.format(r['mapped'] / r['reads'] if r['reads'] else 0)]
                for name in wanted:
                    n = r['targets'][name]
                    vals += [n, '{:.4f}'.format(n / r['reads'] if r['reads'] else 0)]
                vals += ['{:.0f}'.format(r['median_length']), '{:.1f}'.format(r['median_q'])]
                handle.write('\t'.join(str(v) for v in vals) + '\n')

    with open(os.path.join(args.outdir, 'top_genes.tsv'), 'w') as handle:
        handle.write('sample_id\tflowcell\trank\tgene\treads\tfraction_of_mapped\n')
        for r in results:
            for rank, (gene, n) in enumerate(r['top'][:args.top_n], 1):
                handle.write('{}\t{}\t{}\t{}\t{}\t{:.4f}\n'.format(
                    r['sample_id'], r['flowcell'], rank, gene, n, n / r['mapped'] if r['mapped'] else 0))

    plot(results, wanted, sheet, args.outdir)


def plot(results, wanted, sheet, outdir):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    matplotlib.rcParams['pdf.fonttype'] = 42

    per_sample = [r for r in results if r['label'] == 'sample']
    order = [s for s in sheet if any(r['sample_id'] == s for r in per_sample)]
    per_sample.sort(key=lambda r: order.index(r['sample_id']))
    if not per_sample:
        return
    x = np.arange(len(per_sample))
    width = 0.8 / max(len(wanted), 1)
    fig, ax = plt.subplots(figsize=(max(6, 0.55 * len(x) + 2), 4))
    for i, name in enumerate(wanted):
        rates = [100.0 * r['targets'][name] / r['reads'] if r['reads'] else 0 for r in per_sample]
        ax.bar(x + (i - (len(wanted) - 1) / 2) * width, rates, width, label=name)
    ax.set_xticks(x)
    ax.set_xticklabels([r['sample_id'] for r in per_sample], rotation=60, ha='right', fontsize=8)
    ax.set_ylabel('% of reads on gene')
    ax.legend(frameon=False)
    ax.spines[['top', 'right']].set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, 'qc_on_target.pdf'))
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(max(6, 0.55 * len(x) + 2), 4))
    ax.boxplot([np.asarray(r['lengths']) for r in per_sample], showfliers=False)
    ax.set_xticks(x + 1)
    ax.set_xticklabels([r['sample_id'] for r in per_sample], rotation=60, ha='right', fontsize=8)
    ax.set_ylabel('read length (nt)')
    ax.spines[['top', 'right']].set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, 'qc_read_length.pdf'))
    plt.close(fig)


if __name__ == '__main__':
    main()
