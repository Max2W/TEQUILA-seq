#!/usr/bin/env python3
"""Extract the 31-mers that identify one gene from a FIRE v0.5.0 index.

Run inside FIRE's pixi environment, e.g.
    pixi run --manifest-path $FIRE_DIR/pixi.toml python fire_gene_kmers.py \
        --index $FIRE_INDEX --gene EIF2AK3 --outdir kmers/EIF2AK3

How the k-mers are chosen
-------------------------
* Every transcript of the gene is taken from the index's own transcript FASTA (`tx_fasta/`), using the
  gene table (`genedb/`) to list them, so the k-mers match exactly what FIRE screens with.
* Each distinct 31-mer is looked up in the index's protein-coding k-mer dictionary (`sshash_bundle.pcg`).
  A k-mer is kept as an ANCHOR only when the index marks it globally unique (it occurs in this one gene
  across the full reference, coding and non-coding) and its gene set is exactly this gene. These are
  the k-mers FIRE itself treats as
  anchors, so a read that carries them comes from this gene and not from a paralog (for EIF2AK3:
  EIF2AK1/2/4 share kinase-domain sequence). This is the same rule as FIRE's read decoder
  (`decode_anchors_and_colour_classes`): the string's global-uniqueness bit, and its one gene.

Outputs (in --outdir)
---------------------
<GENE>.anchor_kmers.fa   one record per anchor k-mer (input for the read bait step)
<GENE>.kmers.tsv         every distinct k-mer of the gene with its status and the genes that share it
<GENE>.transcripts.tsv   per transcript: length, number of k-mers, number and fraction of anchors
<GENE>.summary.json      counts, the index path and what the index was built from
"""
import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--index', required=True, help='FIRE v0.5.0 index directory (has meta.json)')
    parser.add_argument('--gene', required=True, help='protein-coding gene symbol, e.g. EIF2AK3')
    parser.add_argument('--scope', default='pcg', help='k-mer dictionary scope (default: pcg)')
    parser.add_argument('--outdir', required=True)
    parser.add_argument('--check', action='store_true',
                        help='only open the index and resolve the gene, then exit')
    return parser.parse_args()


def main():
    args = parse_args()
    try:
        from xinglab_fire.index.index_nova import FireIndexNova
    except ImportError as e:
        sys.exit('cannot import FIRE ({}). Run this through `pixi run --manifest-path '
                 '<fire>/pixi.toml python ...` after `pixi run build-ext`.'.format(e))

    index = FireIndexNova(Path(args.index), verify_deep=False)
    gene_db = index.gene_db
    k = index.k
    if args.scope not in index.meta.scopes:
        sys.exit('scope {!r} not in this index ({})'.format(args.scope, index.meta.scopes))
    if args.gene not in gene_db.d_gene_symbol_to_n:
        sys.exit('{} is not a protein-coding gene in this index'.format(args.gene))
    gene_n = gene_db.d_gene_symbol_to_n[args.gene]
    locus = gene_db.genes[gene_n]
    transcripts = gene_db.transcripts_of(args.gene)
    contig = gene_db.contig_names[locus.contig_n]
    print('{}: gene_n={} ids={} {}:{}-{}({}) {} transcripts, k={}'.format(
        args.gene, gene_n, ','.join(locus.gene_ids), contig, locus.start_0b + 1, locus.end_0b,
        '+' if locus.strand else '-', len(transcripts), k))
    if args.check:
        return

    bundle = index.sshash_bundle[args.scope]
    gene_symbols = index.gene_symbols

    # distinct k-mers of all transcripts, with the transcripts that hold each
    kmer_to_tx = defaultdict(set)
    tx_kmers = dict()
    tx_len = dict()
    for tx in transcripts:
        seq = index.tx_fasta.fetch(tx.transcript_id).upper()
        tx_len[tx.transcript_id] = len(seq)
        kmers = set()
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i + k]
            if 'N' in kmer:
                continue
            kmers.add(kmer)
            kmer_to_tx[kmer].add(tx.transcript_id)
        tx_kmers[tx.transcript_id] = kmers

    kmer_list = sorted(kmer_to_tx)
    if not kmer_list:
        sys.exit('no k-mers for {}'.format(args.gene))
    kmer_id, _orientation, string_n = bundle.sshash_dict.lookup_string(
        ''.join(kmer_list).encode('ascii'), len(kmer_list))
    kmer_id = np.asarray(kmer_id)
    string_n = np.asarray(string_n)

    rows = list()
    anchors = list()
    status_counts = defaultdict(int)
    for i, kmer in enumerate(kmer_list):
        if kmer_id[i] < 0:
            status, others = 'absent', ''
        else:
            s = int(string_n[i])
            genes = [int(g) for g in bundle.genes_for_string(s)]
            others = ','.join(gene_symbols[g] for g in genes if g != gene_n)
            if bundle.is_globally_unique(s) and genes == [gene_n]:
                status = 'anchor'
            elif genes == [gene_n]:
                # only this gene among protein-coding genes, but also in a non-coding gene
                status = 'unique_in_scope_only'
            else:
                status = 'shared'
        status_counts[status] += 1
        if status == 'anchor':
            anchors.append(kmer)
        rows.append((kmer, status, len(kmer_to_tx[kmer]), others))

    anchor_set = set(anchors)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    prefix = outdir / args.gene
    with open('{}.anchor_kmers.fa'.format(prefix), 'w') as handle:
        for i, kmer in enumerate(anchors):
            handle.write('>{}_anchor_{}\n{}\n'.format(args.gene, i, kmer))
    with open('{}.kmers.tsv'.format(prefix), 'w') as handle:
        handle.write('kmer\tstatus\tn_transcripts\tshared_with\n')
        for row in rows:
            handle.write('{}\t{}\t{}\t{}\n'.format(*row))
    with open('{}.transcripts.tsv'.format(prefix), 'w') as handle:
        handle.write('transcript_id\tcanonical\tlength\tn_kmers\tn_anchor_kmers\tanchor_fraction\n')
        for tx in transcripts:
            kmers = tx_kmers[tx.transcript_id]
            n_anchor = len(kmers & anchor_set)
            handle.write('{}\t{}\t{}\t{}\t{}\t{:.3f}\n'.format(
                tx.transcript_id, tx.ensembl_canonical, tx_len[tx.transcript_id], len(kmers), n_anchor,
                n_anchor / len(kmers) if kmers else 0.0))

    summary = {
        'gene': args.gene,
        'gene_ids': list(locus.gene_ids),
        'locus': '{}:{}-{}'.format(contig, locus.start_0b + 1, locus.end_0b),
        'strand': '+' if locus.strand else '-',
        'k': k,
        'scope': args.scope,
        'n_transcripts': len(transcripts),
        'n_distinct_kmers': len(kmer_list),
        'status_counts': dict(status_counts),
        'n_anchor_kmers': len(anchors),
        'index': str(Path(args.index).resolve()),
        'index_inputs': dict(index.meta.inputs),
    }
    with open('{}.summary.json'.format(prefix), 'w') as handle:
        json.dump(summary, handle, indent=2, default=str)
    print(json.dumps({key: summary[key] for key in
                      ('gene', 'locus', 'n_transcripts', 'n_distinct_kmers', 'status_counts')}, indent=2))
    if not anchors:
        sys.exit('no anchor k-mers found for {}; the bait step would select nothing'.format(args.gene))


if __name__ == '__main__':
    main()
