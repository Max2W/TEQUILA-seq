#!/usr/bin/env python3
"""Compare isoform proportions of one gene between quantification runs (e.g. ESPRESSO fast / full,
AMALGAM fast / full). Isoforms are matched by intron chain, since novel isoform IDs differ between tools
(mono-exonic isoforms are matched only by ID).

Input : --run LABEL=<dir with <GENE>_<LABEL>_isoform_annotation.tsv and _isoform_proportions.tsv> (repeat)
Output: <out_prefix>.tsv (one row per intron chain, mean proportion per run)
        <out_prefix>_pairs.tsv (per pair of runs: Pearson r over per-sample proportions, shared chains)
        <out_prefix>.pdf (scatter of per-sample proportions for each pair of runs)
"""
import argparse
import itertools
import os
import sys

import numpy as np
import pandas as pd


def load(gene, label, path):
    ann_path = os.path.join(path, '{}_{}_isoform_annotation.tsv'.format(gene, label))
    prop_path = os.path.join(path, '{}_{}_isoform_proportions.tsv'.format(gene, label))
    if not (os.path.exists(ann_path) and os.path.exists(prop_path)):
        print('skip {}: no report in {}'.format(label, path), file=sys.stderr)
        return None
    ann = pd.read_csv(ann_path, sep='\t', index_col=0)
    props = pd.read_csv(prop_path, sep='\t', index_col=0).drop(columns=['label'])
    key = ann['intron_chain'].where(ann['intron_chain'] != 'monoexonic', 'mono:' + ann.index.astype(str))
    props.index = key.reindex(props.index)
    props = props.groupby(level=0).sum(min_count=1)
    names = ann.set_index(key)['label'].groupby(level=0).first()
    return props, names


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--gene', required=True)
    parser.add_argument('--run', action='append', required=True)
    parser.add_argument('--out-prefix', required=True)
    args = parser.parse_args()

    runs = dict()
    for spec in args.run:
        label, path = spec.split('=', 1)
        loaded = load(args.gene, label, path)
        if loaded:
            runs[label] = loaded
    if len(runs) < 2:
        sys.exit('need at least two runs with reports to compare')

    chains = sorted(set().union(*(p.index for p, _ in runs.values())))
    table = pd.DataFrame(index=chains)
    for label, (props, names) in runs.items():
        table['label_' + label] = names.reindex(chains)
        table['mean_prop_' + label] = props.reindex(chains).mean(axis=1)
    table.index.name = 'intron_chain'
    table.sort_values([c for c in table.columns if c.startswith('mean_prop_')][0], ascending=False) \
        .to_csv(args.out_prefix + '.tsv', sep='\t', float_format='%.4f')

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    matplotlib.rcParams['pdf.fonttype'] = 42
    pairs = list(itertools.combinations(runs, 2))
    fig, axes = plt.subplots(1, len(pairs), figsize=(3.4 * len(pairs), 3.3), squeeze=False)
    rows = list()
    for ax, (a, b) in zip(axes[0], pairs):
        pa, pb = runs[a][0], runs[b][0]
        samples = [s for s in pa.columns if s in pb.columns]
        shared = sorted(set(pa.index) & set(pb.index))
        x = pa.reindex(chains)[samples].fillna(0).to_numpy().ravel()
        y = pb.reindex(chains)[samples].fillna(0).to_numpy().ravel()
        ok = ~(np.isnan(x) | np.isnan(y))
        r = float(np.corrcoef(x[ok], y[ok])[0, 1]) if ok.sum() > 2 else np.nan
        rows.append({'run_a': a, 'run_b': b, 'pearson_r': r, 'chains_a': len(pa), 'chains_b': len(pb),
                     'chains_shared': len(shared)})
        ax.scatter(x, y, s=6, alpha=0.6)
        ax.plot([0, 1], [0, 1], color='0.6', lw=0.8)
        ax.set_xlabel(a)
        ax.set_ylabel(b)
        ax.set_title('{}: r = {:.2f}'.format(args.gene, r), fontsize=8)
        ax.spines[['top', 'right']].set_visible(False)
    pd.DataFrame(rows).to_csv(args.out_prefix + '_pairs.tsv', sep='\t', index=False, float_format='%.4f')
    fig.tight_layout()
    fig.savefig(args.out_prefix + '.pdf')
    print(pd.DataFrame(rows).to_string(index=False))


if __name__ == '__main__':
    main()
