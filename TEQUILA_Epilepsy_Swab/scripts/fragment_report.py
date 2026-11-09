#!/usr/bin/env python3
"""Detection and fragmentation report for target genes (degraded RNA, e.g. buccal swabs).

For every target gene and sample it answers two questions:
  1. Is the gene detected?  reads on the gene, and reads per million of all reads sequenced
  2. How intact are its molecules?  each read is placed on the gene's canonical transcript
     (Ensembl_canonical / MANE Select) and summarised by
       - transcript coverage: the fraction of the canonical transcript between the read's 5'-most
         and 3'-most aligned transcript positions
       - full length: the read reaches both the first and the last exon of the canonical transcript
       - 3' anchored / 5' anchored: the read reaches the last / first exon
     plus a gene-body coverage profile (5' -> 3', in 100 bins). Fragmented RNA captured from the 3'
     end shows coverage that falls off toward the 5' end.

Reads are taken from --bam-pattern, with {gene} and {sample} filled in (gene BAMs of the fast path,
or genome-wide BAMs); only primary, mapped alignments overlapping the gene are used.

Outputs (in --outdir)
  fragmentation_summary.tsv        one row per gene and sample
  read_tx_coverage.tsv.gz          one row per read: gene, sample, length, 5'/3' transcript positions
  gene_body_coverage.tsv           coverage profile per gene, sample and bin
  <GENE>_fragmentation.pdf/.png    coverage profile, coverage-fraction curves, full-length counts
  overview.pdf/.png                reads per million and % full length, all genes x samples
  fragmentation_report.md          the tables in readable form
"""
import argparse
import csv
import gzip
import os
import re
import sys
from collections import defaultdict

import numpy as np
import pysam

ATTR = re.compile(r'(\S+) "([^"]*)"')
N_BINS = 100


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--genes', required=True, help='comma separated gene names')
    parser.add_argument('--gtf', required=True, help='GENCODE GTF')
    parser.add_argument('--sample-sheet', required=True)
    parser.add_argument('--bam-pattern', required=True,
                        help='path with {gene} and/or {sample}, e.g. fast/bam/{gene}/{sample}.bam')
    parser.add_argument('--read-counts', help='fastq_read_counts.tsv (total reads per sample, for RPM)')
    parser.add_argument('--label', default='', help='added to titles, e.g. "fast path"')
    parser.add_argument('--outdir', required=True)
    return parser.parse_args()


def canonical_models(gtf, genes):
    """gene -> dict(chrom, strand, start, end, tx_id, tx_name, exons [(s, e) 1-based, genomic order])."""
    wanted = set(genes)
    gene_span, tx_info, exons = dict(), dict(), defaultdict(list)
    with open(gtf) as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9:
                continue
            attrs = dict(ATTR.findall(f[8]))
            name = attrs.get('gene_name')
            if name not in wanted:
                continue
            if f[2] == 'gene':
                gene_span[name] = (f[0], int(f[3]), int(f[4]), f[6])
            elif f[2] == 'transcript' and 'tag "Ensembl_canonical"' in f[8]:
                tx_info.setdefault(name, []).append((attrs['transcript_id'], attrs.get('transcript_name', ''),
                                                     'tag "MANE_Select"' in f[8]))
            elif f[2] == 'exon':
                exons[attrs['transcript_id']].append((int(f[3]), int(f[4])))
    models = dict()
    for gene in genes:
        if gene not in gene_span or gene not in tx_info:
            sys.exit('{}: no gene record or no Ensembl_canonical transcript in {}'.format(gene, args_gtf(gtf)))
        cands = sorted(tx_info[gene], key=lambda t: not t[2])   # MANE first
        tx_id, tx_name, _ = cands[0]
        chrom, start, end, strand = gene_span[gene]
        models[gene] = {'chrom': chrom, 'strand': strand, 'start': start, 'end': end,
                        'tx_id': tx_id, 'tx_name': tx_name, 'exons': sorted(exons[tx_id])}
    return models


def args_gtf(path):
    return os.path.basename(path)


class TxMap:
    """Genomic <-> transcript coordinates of one transcript, transcript positions counted 5' -> 3'."""

    def __init__(self, exons, strand):
        self.exons = exons                          # genomic order, 1-based inclusive
        self.strand = strand
        self.lengths = [e - s + 1 for s, e in exons]
        self.length = sum(self.lengths)
        offs = np.cumsum([0] + self.lengths[:-1])
        self.offsets = offs                         # offset of each exon's genomic-left base (plus order)
        # exon index (in transcript 5'->3' order) of each exon
        n = len(exons)
        self.tx_rank = list(range(n)) if strand == '+' else list(range(n - 1, -1, -1))

    def tx_pos(self, gpos):
        """1-based genomic position -> 0-based transcript position (5'->3'), or None if not exonic."""
        for i, (s, e) in enumerate(self.exons):
            if s <= gpos <= e:
                p = self.offsets[i] + (gpos - s)
                return int(p) if self.strand == '+' else int(self.length - 1 - p)
        return None

    def first_last_exon_hit(self, blocks):
        """Which canonical exons (transcript order) the aligned blocks touch."""
        hit = set()
        for bs, be in blocks:                       # 1-based inclusive genomic
            for i, (s, e) in enumerate(self.exons):
                if bs <= e and be >= s:
                    hit.add(self.tx_rank[i])
        return hit


def read_records(bam_path, model):
    """Yield (read_length, tx_first, tx_last, exon_ranks_hit) for primary reads on the gene."""
    tmap = model['map']
    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        if model['chrom'] not in bam.references:
            return
        for rec in bam.fetch(model['chrom'], model['start'] - 1, model['end']):
            if rec.is_unmapped or rec.is_secondary or rec.is_supplementary:
                continue
            blocks = [(s + 1, e) for s, e in rec.get_blocks()]   # 1-based inclusive
            positions = list()
            for bs, be in blocks:
                for s, e in tmap.exons:
                    lo, hi = max(bs, s), min(be, e)
                    if lo <= hi:
                        positions.append(tmap.tx_pos(lo))
                        positions.append(tmap.tx_pos(hi))
            length = rec.query_length or rec.infer_read_length() or 0
            if not positions:
                yield length, None, None, set()
                continue
            yield length, min(positions), max(positions), tmap.first_last_exon_hit(blocks)


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    genes = [g for g in args.genes.split(',') if g]
    with open(args.sample_sheet) as handle:
        sheet = list(csv.DictReader(handle, delimiter='\t'))
    samples = [r['sample_id'] for r in sheet]
    meta = {r['sample_id']: r for r in sheet}

    total_reads = defaultdict(int)
    if args.read_counts and os.path.exists(args.read_counts):
        with open(args.read_counts) as handle:
            for row in csv.DictReader(handle, delimiter='\t'):
                if row['sample_id'] in meta:
                    total_reads[row['sample_id']] += int(row['n_reads'])

    models = canonical_models(args.gtf, genes)
    for m in models.values():
        m['map'] = TxMap(m['exons'], m['strand'])

    summary, profiles, ecdf = list(), dict(), dict()
    per_read = gzip.open(os.path.join(args.outdir, 'read_tx_coverage.tsv.gz'), 'wt')
    per_read.write('gene\tsample_id\tread_length\ttx_5prime_pos\ttx_3prime_pos\ttx_coverage\t'
                   'first_exon\tlast_exon\n')
    for gene in genes:
        m = models[gene]
        tmap = m['map']
        n_exons = len(m['exons'])
        for sample in samples:
            bam = args.bam_pattern.format(gene=gene, sample=sample)
            prof = np.zeros(N_BINS)
            cov_fracs, lengths = list(), list()
            n_reads = n_canon = n_full = n_first = n_last = 0
            if os.path.exists(bam):
                for length, t5, t3, ranks in read_records(bam, m):
                    n_reads += 1
                    lengths.append(length)
                    if t5 is None:
                        per_read.write('{}\t{}\t{}\tNA\tNA\t0\tFalse\tFalse\n'.format(gene, sample, length))
                        continue
                    n_canon += 1
                    frac = (t3 - t5 + 1) / tmap.length
                    cov_fracs.append(frac)
                    first, last = 0 in ranks, (n_exons - 1) in ranks
                    n_first += first
                    n_last += last
                    n_full += first and last
                    b0 = int(t5 * N_BINS / tmap.length)
                    b1 = int(t3 * N_BINS / tmap.length)
                    prof[b0:b1 + 1] += 1
                    per_read.write('{}\t{}\t{}\t{}\t{}\t{:.4f}\t{}\t{}\n'.format(
                        gene, sample, length, t5, t3, frac, first, last))
            else:
                print('NOTE: no BAM {}'.format(bam), file=sys.stderr)
            tot = total_reads.get(sample, 0)
            cf = np.asarray(cov_fracs)
            summary.append({
                'gene': gene, 'sample_id': sample, 'donor': meta[sample].get('donor', ''),
                'condition': meta[sample].get('condition', ''),
                'canonical_transcript': '{} ({})'.format(m['tx_name'], m['tx_id']),
                'transcript_length': tmap.length, 'n_exons': n_exons,
                'total_reads_sequenced': tot,
                'reads_on_gene': n_reads,
                'reads_per_million': round(1e6 * n_reads / tot, 3) if tot else '',
                'reads_on_canonical_exons': n_canon,
                'median_read_length': int(np.median(lengths)) if lengths else 0,
                'median_tx_coverage': round(float(np.median(cf)), 3) if len(cf) else 0,
                'pct_reads_cover_ge50pct': round(100 * float((cf >= 0.5).mean()), 1) if len(cf) else 0,
                'pct_reads_cover_ge80pct': round(100 * float((cf >= 0.8).mean()), 1) if len(cf) else 0,
                'reads_reach_last_exon_3prime': n_last,
                'reads_reach_first_exon_5prime': n_first,
                'full_length_reads': n_full,
                'pct_full_length': round(100 * n_full / n_canon, 1) if n_canon else 0,
            })
            profiles[(gene, sample)] = prof / n_canon if n_canon else prof
            ecdf[(gene, sample)] = np.sort(cf)
    per_read.close()

    cols = list(summary[0].keys())
    with open(os.path.join(args.outdir, 'fragmentation_summary.tsv'), 'w') as handle:
        handle.write('\t'.join(cols) + '\n')
        for row in summary:
            handle.write('\t'.join(str(row[c]) for c in cols) + '\n')
    with open(os.path.join(args.outdir, 'gene_body_coverage.tsv'), 'w') as handle:
        handle.write('gene\tsample_id\tbin_pct_5to3\tfraction_of_reads_covering\n')
        for (gene, sample), prof in profiles.items():
            for b, v in enumerate(prof):
                handle.write('{}\t{}\t{}\t{:.4f}\n'.format(gene, sample, b, v))

    plot(args, genes, samples, meta, summary, profiles, ecdf, models)
    write_md(args, genes, summary)
    print('fragmentation report written to {}'.format(args.outdir))


def plot(args, genes, samples, meta, summary, profiles, ecdf, models):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['font.size'] = 8
    conds = list(dict.fromkeys(meta[s].get('condition', '') for s in samples))
    donors = list(dict.fromkeys(meta[s].get('donor', '') for s in samples))
    palette = plt.get_cmap('tab10').colors
    color = {c: palette[i % 10] for i, c in enumerate(conds)}
    styles = ['-', '--', ':', '-.']
    style = {d: styles[i % len(styles)] for i, d in enumerate(donors)}
    rows = {(r['gene'], r['sample_id']): r for r in summary}

    for gene in genes:
        m = models[gene]
        fig, axes = plt.subplots(1, 3, figsize=(13, 3.6))
        ax = axes[0]
        x = (np.arange(N_BINS) + 0.5)
        for s in samples:
            if rows[(gene, s)]['reads_on_canonical_exons'] == 0:
                continue
            ax.plot(x, profiles[(gene, s)], color=color[meta[s]['condition']], ls=style[meta[s]['donor']],
                    lw=1, label=s)
        ax.set_xlabel("position along {} (% of transcript, 5' -> 3')".format(m['tx_name']))
        ax.set_ylabel('fraction of reads covering')
        ax.set_ylim(0, 1.02)
        ax.set_title('{}: gene-body coverage{}'.format(gene, ' ({})'.format(args.label) if args.label else ''))
        ax.legend(fontsize=6, frameon=False)
        ax = axes[1]
        for s in samples:
            v = ecdf[(gene, s)]
            if len(v) == 0:
                continue
            ax.plot(v, 1 - np.arange(len(v)) / len(v), color=color[meta[s]['condition']],
                    ls=style[meta[s]['donor']], lw=1)
        ax.set_xlabel('transcript coverage of a read (fraction of {} nt)'.format(m['map'].length))
        ax.set_ylabel('fraction of reads with at least this coverage')
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1.02)
        ax = axes[2]
        xs = np.arange(len(samples))
        tot = [rows[(gene, s)]['reads_on_canonical_exons'] for s in samples]
        full = [rows[(gene, s)]['full_length_reads'] for s in samples]
        ax.bar(xs, tot, color='0.8', label='reads on gene')
        ax.bar(xs, full, color=[color[meta[s]['condition']] for s in samples], label='full length')
        for i, s in enumerate(samples):
            ax.text(i, max(tot[i], 1), '{}%'.format(rows[(gene, s)]['pct_full_length']), ha='center',
                    va='bottom', fontsize=6)
        ax.set_yscale('symlog', linthresh=10)
        ax.set_xticks(xs)
        ax.set_xticklabels(samples, rotation=60, ha='right')
        ax.set_ylabel('reads')
        ax.set_title('full-length reads (first + last exon)')
        ax.legend(fontsize=6, frameon=False)
        for a in axes:
            a.spines[['top', 'right']].set_visible(False)
        fig.tight_layout()
        for ext in ('pdf', 'png'):
            fig.savefig(os.path.join(args.outdir, '{}_fragmentation.{}'.format(gene, ext)), dpi=200)
        plt.close(fig)

    # overview heatmaps
    fig, axes = plt.subplots(1, 2, figsize=(4 + 0.9 * len(samples), 1.2 + 0.35 * len(genes)))
    for ax, key, title, fmt in ((axes[0], 'reads_per_million', 'reads per million sequenced', '{:.1f}'),
                                (axes[1], 'pct_full_length', '% full-length reads', '{:.0f}')):
        mat = np.array([[float(rows[(g, s)][key] or 0) for s in samples] for g in genes])
        im = ax.imshow(np.log10(mat + 1) if key == 'reads_per_million' else mat, aspect='auto', cmap='viridis')
        for (i, j), v in np.ndenumerate(mat):
            ax.text(j, i, fmt.format(v), ha='center', va='center', fontsize=6, color='white')
        ax.set_xticks(range(len(samples)))
        ax.set_xticklabels(samples, rotation=60, ha='right')
        ax.set_yticks(range(len(genes)))
        ax.set_yticklabels(genes)
        ax.set_title(title + (' (log scale colour)' if key == 'reads_per_million' else ''))
    fig.tight_layout()
    for ext in ('pdf', 'png'):
        fig.savefig(os.path.join(args.outdir, 'overview.{}'.format(ext)), dpi=200)
    plt.close(fig)


def write_md(args, genes, summary):
    lines = ['# Detection and fragmentation of target genes{}'.format(
        ' ({})'.format(args.label) if args.label else ''), '',
        'Full length = the read reaches both the first and the last exon of the canonical transcript. '
        'Transcript coverage = fraction of the canonical transcript between the read\'s 5\'-most and '
        '3\'-most aligned positions.', '']
    for gene in genes:
        rows = [r for r in summary if r['gene'] == gene]
        lines += ['## {} ({}; {} nt, {} exons)'.format(gene, rows[0]['canonical_transcript'],
                                                       rows[0]['transcript_length'], rows[0]['n_exons']), '',
                  '| sample | reads on gene | per million | median read length | median coverage | '
                  '>=80% covered | reach 3\' exon | reach 5\' exon | full length |',
                  '|---|---|---|---|---|---|---|---|---|']
        for r in rows:
            lines.append('| {} | {} | {} | {} | {:.2f} | {}% | {} | {} | {} ({}%) |'.format(
                r['sample_id'], r['reads_on_gene'], r['reads_per_million'], r['median_read_length'],
                r['median_tx_coverage'], r['pct_reads_cover_ge80pct'], r['reads_reach_last_exon_3prime'],
                r['reads_reach_first_exon_5prime'], r['full_length_reads'], r['pct_full_length']))
        lines += ['', '![{0}]({0}_fragmentation.png)'.format(gene), '']
    with open(os.path.join(args.outdir, 'fragmentation_report.md'), 'w') as handle:
        handle.write('\n'.join(lines) + '\n')


if __name__ == '__main__':
    main()
