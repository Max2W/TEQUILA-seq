#!/usr/bin/env python3
"""Merge AMALGAM per-sample Quantify_Transcripts.py outputs (gene_id, transcript_id, count) into one
matrix laid out like ESPRESSO's abundance.esp: transcript_ID, transcript_name, gene_ID, <samples...>.
Sample columns follow the sample sheet order; a sample without a counts file is left out."""
import argparse
import csv
import os
import re
import sys
from collections import OrderedDict


def transcript_names(gtf):
    names = dict()
    pat_id = re.compile(r'transcript_id "([^"]+)"')
    pat_name = re.compile(r'transcript_name "([^"]+)"')
    with open(gtf) as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9 or fields[2] != 'transcript':
                continue
            m_id, m_name = pat_id.search(fields[8]), pat_name.search(fields[8])
            if m_id and m_name:
                names[m_id.group(1)] = m_name.group(1)
    return names


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--counts-dir', required=True, help='<sample>.counts.tsv files')
    parser.add_argument('--sample-sheet', required=True)
    parser.add_argument('--ref-gtf', required=True, help='GENCODE GTF, for transcript names')
    parser.add_argument('--out', required=True)
    args = parser.parse_args()

    with open(args.sample_sheet) as handle:
        samples = [row['sample_id'] for row in csv.DictReader(handle, delimiter='\t')]
    names = transcript_names(args.ref_gtf)

    matrix = OrderedDict()
    used = list()
    for sample in samples:
        path = os.path.join(args.counts_dir, '{}.counts.tsv'.format(sample))
        if not os.path.exists(path):
            continue
        used.append(sample)
        with open(path) as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            missing = {'gene_id', 'transcript_id', 'count'} - set(reader.fieldnames or [])
            if missing:
                sys.exit('{} lacks columns {}'.format(path, sorted(missing)))
            for row in reader:
                key = (row['transcript_id'], row['gene_id'])
                matrix.setdefault(key, dict())[sample] = float(row['count'])
    if not used:
        sys.exit('no counts files in {}'.format(args.counts_dir))

    with open(args.out, 'w') as handle:
        handle.write('\t'.join(['transcript_ID', 'transcript_name', 'gene_ID'] + used) + '\n')
        for (tx, gene), counts in matrix.items():
            handle.write('\t'.join([tx, names.get(tx, 'NA'), gene] +
                                   ['{:.2f}'.format(counts.get(s, 0.0)) for s in used]) + '\n')
    print('wrote {} ({} transcripts x {} samples)'.format(args.out, len(matrix), len(used)))


if __name__ == '__main__':
    main()
