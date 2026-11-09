#!/usr/bin/env python3
"""Pull out the reads that carry a gene's anchor 31-mers (one or several genes in a single pass).

Reads are scanned in batches with numpy: every 31-mer of every read is 2-bit encoded, made canonical
(the smaller of the k-mer and its reverse complement) and looked up in the sorted anchor set, so both
read orientations are found. A read is kept for a gene when it has at least --min-hits anchor k-mer
positions of that gene. Nanopore errors break many 31-mers, so a true read carries far fewer hits than
its length, but --min-hits 3 is already well above what a read from an unrelated gene produces
(anchors are unique to their gene in the FIRE index).

Inputs : --anchors GENE=anchor_kmers.fa (repeatable; from fire_gene_kmers.py), fastq(.gz) files
Outputs: <out_dir>/<GENE>.fastq.gz     selected reads, unchanged
         <out_dir>/<GENE>.hits.tsv     read_id, read_length, anchor_hits, source file
         <out_dir>/bait_summary.tsv    reads scanned / selected per gene and input file, hit histogram
"""
import argparse
import os
import subprocess
import sys

import numpy as np

K = 31
_CODE = np.full(256, 4, dtype=np.uint8)
for _base, _value in zip(b'ACGTacgt', (0, 1, 2, 3, 0, 1, 2, 3)):
    _CODE[_base] = _value


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--anchors', required=True, action='append',
                        help='GENE=path/to/GENE.anchor_kmers.fa; repeat for several genes')
    parser.add_argument('--fastq', required=True, nargs='+', help='input fastq(.gz) files')
    parser.add_argument('--out-dir', required=True)
    parser.add_argument('--min-hits', type=int, default=3)
    parser.add_argument('--batch-bases', type=int, default=20_000_000,
                        help='bases per numpy batch (memory ~ 60 bytes per base)')
    return parser.parse_args()


def canonical_kmers(codes):
    """Canonical 2-bit codes of every k-mer window of `codes` (uint8, 4 = invalid base or separator).

    Returns (canonical uint64 array of length len(codes)-K+1, boolean validity array)."""
    n_windows = len(codes) - K + 1
    if n_windows <= 0:
        return np.zeros(0, dtype=np.uint64), np.zeros(0, dtype=bool)
    invalid = (codes == 4)
    cum = np.concatenate(([0], np.cumsum(invalid, dtype=np.int64)))
    valid = (cum[K:] - cum[:-K]) == 0
    c = np.where(invalid, 0, codes).astype(np.uint64)
    r = np.uint64(3) - c
    # Build k-mers by doubling (1 -> 2 -> 4 -> 8 -> 16 bases), then join 16+8+4+2+1 = 31.
    # f_n[j] encodes bases j..j+n-1 with the first base most significant; r_n[j] encodes the
    # reverse complement of the same bases.
    f = {1: c}
    rv = {1: r}
    for n in (2, 4, 8, 16):
        h = n // 2
        a, b = f[h], rv[h]
        f[n] = (a[:-h] << np.uint64(2 * h)) | a[h:]
        rv[n] = b[:-h] | (b[h:] << np.uint64(2 * h))
    fwd = np.zeros(n_windows, dtype=np.uint64)
    rev = np.zeros(n_windows, dtype=np.uint64)
    offset = 0
    for n in (16, 8, 4, 2, 1):
        remaining = K - offset - n
        fwd |= f[n][offset:offset + n_windows] << np.uint64(2 * remaining)
        rev |= rv[n][offset:offset + n_windows] << np.uint64(2 * offset)
        offset += n
    return np.minimum(fwd, rev), valid


def encode_strings(seqs):
    codes = list()
    for seq in seqs:
        codes.append(_CODE[np.frombuffer(seq, dtype=np.uint8)])
        codes.append(np.array([4], dtype=np.uint8))
    return np.concatenate(codes) if codes else np.zeros(0, dtype=np.uint8)


def load_anchor_fasta(path):
    seqs = list()
    with open(path, 'rb') as handle:
        for line in handle:
            line = line.strip()
            if line and not line.startswith(b'>'):
                if len(line) != K:
                    sys.exit('anchor {!r} in {} is not {} bases'.format(line, path, K))
                seqs.append(line)
    if not seqs:
        sys.exit('no anchors in {}'.format(path))
    canon, valid = canonical_kmers(encode_strings(seqs))
    return np.unique(canon[valid])


def load_anchors(specs):
    """Returns gene names, sorted anchor codes and the gene index of each anchor."""
    genes, codes, owners = list(), list(), list()
    for spec in specs:
        if '=' not in spec:
            sys.exit('--anchors must be GENE=path, got {!r}'.format(spec))
        gene, path = spec.split('=', 1)
        anchors = load_anchor_fasta(path)
        print('{}: {} anchor k-mers from {}'.format(gene, len(anchors), path), file=sys.stderr)
        codes.append(anchors)
        owners.append(np.full(len(anchors), len(genes), dtype=np.int64))
        genes.append(gene)
    codes = np.concatenate(codes)
    owners = np.concatenate(owners)
    order = np.argsort(codes, kind='stable')
    codes, owners = codes[order], owners[order]
    # a k-mer claimed by two genes identifies neither: drop it
    dup = np.zeros(len(codes), dtype=bool)
    same = codes[1:] == codes[:-1]
    dup[1:] |= same
    dup[:-1] |= same
    if dup.any():
        print('dropping {} k-mers shared between genes'.format(int(dup.sum())), file=sys.stderr)
    return genes, codes[~dup], owners[~dup]


def read_fastq(path):
    # decompress in a separate process: much faster than Python's gzip module
    if path.endswith('.gz'):
        proc = subprocess.Popen(['gzip', '-dc', path], stdout=subprocess.PIPE, bufsize=1 << 20)
        handle = proc.stdout
    else:
        proc = None
        handle = open(path, 'rb', buffering=1 << 20)
    with handle:
        while True:
            header = handle.readline()
            if not header:
                break
            seq = handle.readline().rstrip(b'\r\n')
            plus = handle.readline()
            qual = handle.readline().rstrip(b'\r\n')
            if not plus.startswith(b'+') or len(seq) != len(qual):
                sys.exit('malformed fastq record in {}: {!r}'.format(path, header[:80]))
            yield header.rstrip(b'\r\n'), seq, qual
    if proc is not None and proc.wait() != 0:
        sys.exit('gzip failed on {}'.format(path))


def count_hits(seqs, anchors, owners, n_genes):
    """Anchor hits per read and gene: int array of shape (len(seqs), n_genes)."""
    codes = encode_strings(seqs)
    canon, valid = canonical_kmers(codes)
    pos = np.searchsorted(anchors, canon)
    pos[pos == len(anchors)] = 0
    hit = valid & (anchors[pos] == canon)
    # window start -> read number (each read is followed by one separator)
    lengths = np.fromiter((len(s) + 1 for s in seqs), dtype=np.int64, count=len(seqs))
    read_of_pos = np.repeat(np.arange(len(seqs)), lengths)[:len(canon)]
    flat = read_of_pos[hit] * n_genes + owners[pos[hit]]
    return np.bincount(flat, minlength=len(seqs) * n_genes).reshape(len(seqs), n_genes)


HIST_EDGES = [0, 1, 2, 3, 5, 10, 20, 50, 100, 200, 500, 10 ** 9]


def hist_labels():
    labels = list()
    for lo, hi in zip(HIST_EDGES[:-1], HIST_EDGES[1:]):
        labels.append('hits_{}'.format(lo) if hi - 1 == lo else 'hits_{}-{}'.format(lo, hi - 1))
    labels[-1] = 'hits_{}+'.format(HIST_EDGES[-2])
    return labels


def main():
    args = parse_args()
    genes, anchors, owners = load_anchors(args.anchors)
    n_genes = len(genes)
    os.makedirs(args.out_dir, exist_ok=True)

    out_fq, out_hits = list(), list()
    for gene in genes:
        prefix = os.path.join(args.out_dir, gene)
        # gzip in a separate process so compression does not slow the scan
        proc = subprocess.Popen(['gzip', '-c', '-4'], stdin=subprocess.PIPE,
                                stdout=open(prefix + '.fastq.gz.tmp', 'wb'))
        out_fq.append(proc)
        handle = open(prefix + '.hits.tsv.tmp', 'w')
        handle.write('read_id\tread_length\tanchor_hits\tsource\n')
        out_hits.append(handle)

    summary = list()
    for path in args.fastq:
        source = os.path.basename(path)
        n_scanned = 0
        n_selected = np.zeros(n_genes, dtype=np.int64)
        bases_selected = np.zeros(n_genes, dtype=np.int64)
        hist = np.zeros((n_genes, len(HIST_EDGES) - 1), dtype=np.int64)
        batch = list()
        batch_bases = 0

        def flush():
            nonlocal n_scanned
            if not batch:
                return
            hits = count_hits([rec[1] for rec in batch], anchors, owners, n_genes)
            n_scanned += len(batch)
            for g in range(n_genes):
                hist[g] += np.histogram(hits[:, g], bins=HIST_EDGES)[0]
                for i in np.nonzero(hits[:, g] >= args.min_hits)[0]:
                    header, seq, qual = batch[i]
                    n_selected[g] += 1
                    bases_selected[g] += len(seq)
                    out_fq[g].stdin.write(header + b'\n' + seq + b'\n+\n' + qual + b'\n')
                    read_id = header[1:].split()[0].decode()
                    out_hits[g].write('{}\t{}\t{}\t{}\n'.format(read_id, len(seq), int(hits[i, g]), source))
            batch.clear()

        for record in read_fastq(path):
            batch.append(record)
            batch_bases += len(record[1])
            if batch_bases >= args.batch_bases:
                flush()
                batch_bases = 0
        flush()
        for g, gene in enumerate(genes):
            summary.append((gene, source, n_scanned, int(n_selected[g]), int(bases_selected[g]), hist[g]))
            print('{} {}: {} reads scanned, {} selected'.format(gene, source, n_scanned, n_selected[g]),
                  file=sys.stderr)

    for g, gene in enumerate(genes):
        prefix = os.path.join(args.out_dir, gene)
        out_fq[g].stdin.close()
        if out_fq[g].wait() != 0:
            sys.exit('gzip failed for {}'.format(gene))
        out_hits[g].close()
        os.replace(prefix + '.fastq.gz.tmp', prefix + '.fastq.gz')
        os.replace(prefix + '.hits.tsv.tmp', prefix + '.hits.tsv')

    with open(os.path.join(args.out_dir, 'bait_summary.tsv'), 'w') as handle:
        handle.write('gene\tsource\treads_scanned\treads_selected\tfraction_selected\tbases_selected\t'
                     'min_hits\t' + '\t'.join(hist_labels()) + '\n')
        for gene, source, n_scanned, n_sel, bases, hist in summary:
            handle.write('{}\t{}\t{}\t{}\t{:.6f}\t{}\t{}\t{}\n'.format(
                gene, source, n_scanned, n_sel, n_sel / n_scanned if n_scanned else 0.0, bases,
                args.min_hits, '\t'.join(str(x) for x in hist)))


if __name__ == '__main__':
    main()
