#!/usr/bin/env python3
"""Isoform report for one gene (EIF2AK3 by default) from ESPRESSO or AMALGAM output.

Inputs
  --matrix       counts matrix: transcript_ID, transcript_name, gene_ID, <one column per sample>
                 (ESPRESSO samples_N2_R0_abundance.esp, or amalgam_counts_matrix.tsv)
  --isoform-gtf  isoform structures (ESPRESSO samples_N2_R0_updated.gtf, or AMALGAM transcripts GTF)
  --ref-gtf      GENCODE GTF (gene, transcript names/types/tags, annotated intron chains)
  --sample-sheet samples.tsv (sample_id, donor, cell_type, condition, ...)

What it does
  * keeps the isoforms of the gene: same strand and overlapping one of the gene's annotated exons
    (works for ESPRESSO and AMALGAM IDs, including novel isoforms)
  * annotates each isoform: annotated vs novel, intron chain, junctions absent from the annotation,
    closest annotated transcript, and (with --find-as) the alternative splicing events relative to the
    canonical transcript using the repo's scripts/Find_AS_Events.py
  * counts -> within-gene isoform proportions per sample
  * contrasts against Vehicle within each cell type, paired by donor:
      - delta proportion per donor, mean delta, paired t-test across donors (MDMs: n = 3 for TG / LPS,
        n = 2 for Ctl; not computed for n < 3)
      - Cochran-Mantel-Haenszel test on read counts with donor as stratum (one stratum for the
        iPSC neurons). Reads are not independent molecules' worth of evidence across donors, so treat
        CMH p-values as a screen, not as final statistics. BH-FDR within each contrast.
  * figures: structures + proportion heatmap, stacked proportions, per-isoform donor plots,
    gene read depth
"""
import argparse
import csv
import os
import re
import subprocess
import sys
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy import stats

ATTR = re.compile(r'(\S+) "([^"]*)"')
CONDITION_ORDER = ['Ctl', 'Vehicle', 'TG_2h', 'LPS_2h']
CELL_ORDER = ['MDM', 'iPSC_neuron']


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--matrix', required=True)
    parser.add_argument('--isoform-gtf', required=True)
    parser.add_argument('--ref-gtf', required=True)
    parser.add_argument('--sample-sheet', required=True)
    parser.add_argument('--gene', default='EIF2AK3')
    parser.add_argument('--method', default='ESPRESSO', help='label used in titles and file names')
    parser.add_argument('--outdir', required=True)
    parser.add_argument('--find-as', help='path to scripts/Find_AS_Events.py (optional)')
    parser.add_argument('--min-gene-reads', type=float, default=20,
                        help='samples with fewer gene reads are shown but excluded from statistics')
    parser.add_argument('--min-prop', type=float, default=0.05,
                        help='an isoform is plotted when it reaches this proportion in any sample')
    parser.add_argument('--max-plot', type=int, default=12, help='at most this many isoforms in figures')
    parser.add_argument('--min-test-prop', type=float, default=0.01,
                        help='isoforms with a mean proportion below this are not tested')
    return parser.parse_args()


# ----------------------------------------------------------------------------------------------------
# GTF parsing
# ----------------------------------------------------------------------------------------------------
def read_gtf_exons(path, keep=None):
    """transcript_id -> dict(chrom, strand, exons=[(start1, end1)...] sorted, attrs of transcript line)."""
    tx = dict()
    with open(path) as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9 or f[2] not in ('exon', 'transcript'):
                continue
            attrs = dict(ATTR.findall(f[8]))
            tid = attrs.get('transcript_id')
            if tid is None or (keep is not None and not keep(f, attrs)):
                continue
            rec = tx.setdefault(tid, {'chrom': f[0], 'strand': f[6], 'exons': [], 'attrs': {}, 'tags': set()})
            if f[2] == 'exon':
                rec['exons'].append((int(f[3]), int(f[4])))
            else:
                rec['attrs'] = attrs
                rec['tags'] = set(re.findall(r'tag "([^"]+)"', f[8]))
    for rec in tx.values():
        rec['exons'].sort()
    return tx


def gene_lines(ref_gtf, gene):
    needle = 'gene_name "{}";'.format(gene)
    with open(ref_gtf) as handle:
        return [line for line in handle if not line.startswith('#') and needle in line]


def introns(exons):
    return tuple((exons[i][1] + 1, exons[i + 1][0] - 1) for i in range(len(exons) - 1))


def strip_version(x):
    return x.split('.')[0]


# ----------------------------------------------------------------------------------------------------
# isoform selection and annotation
# ----------------------------------------------------------------------------------------------------
def overlaps(exons_a, exons_b):
    i = j = 0
    while i < len(exons_a) and j < len(exons_b):
        if exons_a[i][1] < exons_b[j][0]:
            i += 1
        elif exons_b[j][1] < exons_a[i][0]:
            j += 1
        else:
            return True
    return False


def annotate(isoforms, ref_tx, canonical_id):
    ann_chains = {tid: introns(r['exons']) for tid, r in ref_tx.items()}
    ann_introns = set(i for chain in ann_chains.values() for i in chain)
    chain_to_ref = defaultdict(list)
    for tid, chain in ann_chains.items():
        chain_to_ref[chain].append(tid)
    rows = list()
    for tid, rec in isoforms.items():
        chain = introns(rec['exons'])
        ref_ids_v = {strip_version(t): t for t in ref_tx}
        is_annotated = tid in ref_tx or strip_version(tid) in ref_ids_v
        ref_id = tid if tid in ref_tx else ref_ids_v.get(strip_version(tid))
        same_chain = [t for t in chain_to_ref.get(chain, []) if len(chain) > 0]
        best, best_j = 'NA', 0.0
        for rid, rchain in ann_chains.items():
            a, b = set(chain), set(rchain)
            if not a and not b:
                continue
            j = len(a & b) / len(a | b) if (a | b) else 0.0
            if j > best_j:
                best, best_j = rid, j
        novel_junctions = [i for i in chain if i not in ann_introns]
        attrs = ref_tx[ref_id]['attrs'] if ref_id else {}
        rows.append({
            'isoform_id': tid,
            'status': 'annotated' if is_annotated else 'novel',
            'transcript_name': attrs.get('transcript_name', 'NA'),
            'transcript_type': attrs.get('transcript_type', 'NA'),
            'tags': ','.join(sorted(ref_tx[ref_id]['tags'])) if ref_id else '',
            'is_canonical': bool(ref_id) and ref_id == canonical_id,
            'chrom': rec['chrom'], 'strand': rec['strand'],
            'start': rec['exons'][0][0], 'end': rec['exons'][-1][1],
            'n_exons': len(rec['exons']),
            'length': sum(e - s + 1 for s, e in rec['exons']),
            'intron_chain': ';'.join('{}-{}'.format(*i) for i in chain) or 'monoexonic',
            'same_intron_chain_as': ','.join(t for t in same_chain if t != ref_id) or '',
            'closest_annotated': best, 'closest_jaccard': round(best_j, 3),
            'n_novel_junctions': len(novel_junctions),
            'novel_junctions': ';'.join('{}-{}'.format(*i) for i in novel_junctions),
        })
    return pd.DataFrame(rows).set_index('isoform_id')


def label_isoforms(ann):
    """Annotated: transcript name (* = canonical). Novel: novel_<rank by mean proportion>."""
    labels = dict()
    novel_n = 0
    for tid, row in ann.iterrows():
        if row['status'] == 'annotated':
            name = row['transcript_name'] if row['transcript_name'] != 'NA' else tid
            labels[tid] = name + (' *' if row['is_canonical'] else '')
        else:
            novel_n += 1
            labels[tid] = 'novel_{}'.format(novel_n)
    return labels


def run_find_as(script, isoforms, ann, gene_id, ref_gene_lines, outdir):
    in_gtf = os.path.join(outdir, 'isoforms_for_AS.gtf')
    ref_gtf = os.path.join(outdir, 'reference_gene.gtf')
    with open(in_gtf, 'w') as handle:
        for tid, rec in isoforms.items():
            attr = 'gene_id "{}"; transcript_id "{}";'.format(gene_id, tid)
            handle.write('\t'.join([rec['chrom'], 'isoform', 'transcript', str(rec['exons'][0][0]),
                                    str(rec['exons'][-1][1]), '.', rec['strand'], '.', attr]) + '\n')
            for s, e in rec['exons']:
                handle.write('\t'.join([rec['chrom'], 'isoform', 'exon', str(s), str(e), '.', rec['strand'],
                                        '.', attr]) + '\n')
    with open(ref_gtf, 'w') as handle:
        handle.writelines(ref_gene_lines)
    out = os.path.join(outdir, 'as_events.tsv')
    result = subprocess.run([sys.executable, script, '-i', in_gtf, '-r', ref_gtf, '-o', out],
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    if result.returncode != 0 or not os.path.exists(out):
        print('WARNING: Find_AS_Events.py failed:\n' + result.stdout[-2000:], file=sys.stderr)
        return None
    events = pd.read_csv(out, sep='\t')
    summary = events.groupby('txID')['event_type'].apply(
        lambda x: ';'.join('{}x{}'.format(n, t) if n > 1 else t for t, n in x.value_counts().items()))
    return summary


# ----------------------------------------------------------------------------------------------------
# statistics
# ----------------------------------------------------------------------------------------------------
def cmh(tables):
    """Cochran-Mantel-Haenszel for 2x2 tables [[a, b], [c, d]]; returns (MH odds ratio, chi2, p)."""
    sum_a = sum_e = sum_v = num = den = 0.0
    for (a, b), (c, d) in tables:
        n = a + b + c + d
        if n < 2:
            continue
        sum_a += a
        sum_e += (a + b) * (a + c) / n
        sum_v += (a + b) * (c + d) * (a + c) * (b + d) / (n * n * (n - 1))
        num += a * d / n
        den += b * c / n
    if sum_v <= 0:
        return np.nan, np.nan, np.nan
    chi2 = (abs(sum_a - sum_e) - 0.5) ** 2 / sum_v
    odds = num / den if den > 0 else np.inf
    return odds, chi2, float(stats.chi2.sf(chi2, 1))


def bh(pvalues):
    p = np.asarray(pvalues, dtype=float)
    out = np.full(len(p), np.nan)
    ok = ~np.isnan(p)
    if ok.sum() == 0:
        return out
    pv = p[ok]
    order = np.argsort(pv)
    ranked = pv[order] * len(pv) / (np.arange(len(pv)) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    res = np.empty(len(pv))
    res[order] = np.minimum(ranked, 1.0)
    out[ok] = res
    return out


def contrasts(counts, props, gene_total, sheet, tested, min_gene_reads):
    rows = list()
    for cell in sheet['cell_type'].unique():
        cs = sheet[sheet['cell_type'] == cell]
        for cond in [c for c in CONDITION_ORDER if c != 'Vehicle'] + \
                [c for c in cs['condition'].unique() if c not in CONDITION_ORDER]:
            pairs = list()
            for donor, ds in cs.groupby('donor'):
                treated = ds[ds['condition'] == cond]['sample_id'].tolist()
                vehicle = ds[ds['condition'] == 'Vehicle']['sample_id'].tolist()
                if treated and vehicle and gene_total[treated[0]] >= min_gene_reads \
                        and gene_total[vehicle[0]] >= min_gene_reads:
                    pairs.append((donor, treated[0], vehicle[0]))
            if not pairs:
                continue
            block = list()
            for iso in tested:
                deltas = [props.at[iso, t] - props.at[iso, v] for _, t, v in pairs]
                tables = [((counts.at[iso, t], gene_total[t] - counts.at[iso, t]),
                           (counts.at[iso, v], gene_total[v] - counts.at[iso, v])) for _, t, v in pairs]
                odds, chi2, p_cmh = cmh(tables)
                p_t = float(stats.ttest_1samp(deltas, 0.0).pvalue) if len(deltas) >= 3 else np.nan
                block.append({
                    'cell_type': cell, 'contrast': '{}_vs_Vehicle'.format(cond), 'isoform_id': iso,
                    'n_donors': len(pairs), 'donors': ','.join(d for d, _, _ in pairs),
                    'mean_prop_treated': np.mean([props.at[iso, t] for _, t, _ in pairs]),
                    'mean_prop_vehicle': np.mean([props.at[iso, v] for _, _, v in pairs]),
                    'mean_delta_prop': float(np.mean(deltas)),
                    'per_donor_delta': ','.join('{:+.3f}'.format(x) for x in deltas),
                    'same_direction_all_donors': bool(len(deltas) > 1 and
                                                      (all(x > 0 for x in deltas) or all(x < 0 for x in deltas))),
                    'paired_t_p': p_t, 'cmh_odds_ratio': odds, 'cmh_chi2': chi2, 'cmh_p': p_cmh,
                })
            frame = pd.DataFrame(block)
            frame['paired_t_fdr'] = bh(frame['paired_t_p'])
            frame['cmh_fdr'] = bh(frame['cmh_p'])
            rows.append(frame)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


# ----------------------------------------------------------------------------------------------------
# figures
# ----------------------------------------------------------------------------------------------------
def compressed_axis(isoforms, ids, intron_width=0.04):
    """Map genomic positions to a plot axis where every exonic base counts and each intron gap is fixed."""
    blocks = sorted(e for tid in ids for e in isoforms[tid]['exons'])
    merged = list()
    for s, e in blocks:
        if merged and s <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    exonic = sum(e - s + 1 for s, e in merged)
    gap = intron_width * exonic
    starts, x = list(), 0.0
    for s, e in merged:
        starts.append(x)
        x += (e - s + 1) + gap
    total = x - gap

    def to_x(pos):
        for (s, e), x0 in zip(merged, starts):
            if pos <= e:
                return x0 + max(pos - s, 0)
        return total
    return to_x, total


def plot_all(args, isoforms, plot_ids, labels, props, counts, gene_total, sheet, contrast_df, strand):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['font.size'] = 8

    samples = sheet['sample_id'].tolist()
    prefix = os.path.join(args.outdir, '{}_{}'.format(args.gene, args.method))
    palette = plt.get_cmap('tab20').colors
    shown = [i for i in plot_ids if i != 'Other']
    colors = {iso: palette[k % len(palette)] for k, iso in enumerate(shown)}
    colors['Other'] = (0.75, 0.75, 0.75)

    # 1. structures + heatmap
    to_x, total = compressed_axis(isoforms, shown)
    height = 0.35 * len(shown) + 1.6
    fig, (ax_s, ax_h) = plt.subplots(1, 2, figsize=(7 + 0.38 * len(samples), height),
                                     gridspec_kw={'width_ratios': [1.3, 0.38 * len(samples) / 2.2]})
    for row, iso in enumerate(shown):
        y = len(shown) - 1 - row
        exons = isoforms[iso]['exons']
        ax_s.plot([to_x(exons[0][0]), to_x(exons[-1][1])], [y, y], color='0.4', lw=0.8, zorder=1)
        for s, e in exons:
            ax_s.add_patch(Rectangle((to_x(s), y - 0.3), max(to_x(e) - to_x(s), total * 0.002), 0.6,
                                     color=colors[iso], zorder=2))
    ax_s.set_yticks(range(len(shown)))
    ax_s.set_yticklabels([labels[i] for i in reversed(shown)])
    ax_s.set_xlim(-total * 0.01, total * 1.01)
    ax_s.set_ylim(-0.7, len(shown) - 0.3)
    ax_s.set_xticks([])
    ax_s.set_title('{} isoforms ({} strand; introns shrunk; 5\'{})'.format(
        args.gene, strand, ' on the right' if strand == '-' else ' on the left'), fontsize=8)
    for side in ('top', 'right', 'bottom'):
        ax_s.spines[side].set_visible(False)
    if strand == '-':
        ax_s.invert_xaxis()
    mat = props.loc[shown, samples].to_numpy(dtype=float)
    im = ax_h.imshow(mat, aspect='auto', cmap='viridis', vmin=0, vmax=max(0.01, np.nanmax(mat)))
    ax_h.set_xticks(range(len(samples)))
    ax_h.set_xticklabels(['{} (n={:.0f})'.format(s, gene_total[s]) for s in samples], rotation=70,
                         ha='right')
    ax_h.set_yticks([])
    for (r, c), v in np.ndenumerate(mat):
        if not np.isnan(v):
            ax_h.text(c, r, '{:.0f}'.format(100 * v), ha='center', va='center', fontsize=6,
                      color='white' if v < 0.6 * np.nanmax(mat) else 'black')
    fig.colorbar(im, ax=ax_h, fraction=0.04, pad=0.02, label='isoform proportion')
    ax_h.set_title('% of {} reads ({})'.format(args.gene, args.method), fontsize=8)
    fig.tight_layout()
    fig.savefig(prefix + '_structure_heatmap.pdf')
    fig.savefig(prefix + '_structure_heatmap.png', dpi=200)
    plt.close(fig)

    # 2. stacked proportions
    fig, ax = plt.subplots(figsize=(1.5 + 0.45 * len(samples), 4))
    bottom = np.zeros(len(samples))
    for iso in plot_ids:
        vals = np.nan_to_num(props.loc[iso, samples].to_numpy(dtype=float))
        ax.bar(range(len(samples)), vals, bottom=bottom, color=colors[iso],
               label=labels.get(iso, iso), width=0.8)
        bottom += vals
    ax.set_xticks(range(len(samples)))
    ax.set_xticklabels(samples, rotation=70, ha='right')
    ax.set_ylabel('isoform proportion')
    ax.set_ylim(0, 1)
    ax.legend(bbox_to_anchor=(1.01, 1), loc='upper left', frameon=False, fontsize=6)
    ax.spines[['top', 'right']].set_visible(False)
    fig.tight_layout()
    fig.savefig(prefix + '_stacked_proportions.pdf')
    plt.close(fig)

    # 3. per-isoform donor plots
    n = len(shown)
    ncol = min(4, n)
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.0 * ncol, 2.4 * nrow), squeeze=False)
    for k, iso in enumerate(shown):
        ax = axes[k // ncol][k % ncol]
        offset = 0
        ticks, ticklabels = list(), list()
        for cell in [c for c in CELL_ORDER if c in set(sheet['cell_type'])] + \
                [c for c in sheet['cell_type'].unique() if c not in CELL_ORDER]:
            cs = sheet[sheet['cell_type'] == cell]
            conds = [c for c in CONDITION_ORDER if c in set(cs['condition'])]
            for donor, ds in cs.groupby('donor'):
                xs, ys = list(), list()
                for j, cond in enumerate(conds):
                    sids = ds[ds['condition'] == cond]['sample_id'].tolist()
                    if sids and gene_total[sids[0]] >= args.min_gene_reads:
                        xs.append(offset + j)
                        ys.append(props.at[iso, sids[0]])
                ax.plot(xs, ys, marker='o', ms=3, lw=0.8, label=donor)
            ticks += [offset + j for j in range(len(conds))]
            ticklabels += ['{}\n{}'.format(c, 'MDM' if cell == 'MDM' else 'iN') for c in conds]
            offset += len(conds) + 0.8
        ax.set_xticks(ticks)
        ax.set_xticklabels(ticklabels, fontsize=6)
        ax.set_title(labels[iso], fontsize=7)
        ax.set_ylabel('proportion')
        ax.spines[['top', 'right']].set_visible(False)
    for k in range(n, nrow * ncol):
        axes[k // ncol][k % ncol].axis('off')
    handles, names = axes[0][0].get_legend_handles_labels()
    fig.legend(handles, names, loc='upper right', fontsize=6, frameon=False)
    fig.tight_layout()
    fig.savefig(prefix + '_isoform_by_donor.pdf')
    plt.close(fig)

    # 4. gene depth
    fig, ax = plt.subplots(figsize=(1.5 + 0.45 * len(samples), 3))
    ax.bar(range(len(samples)), [gene_total[s] for s in samples], color='0.4')
    ax.axhline(args.min_gene_reads, color='red', lw=0.8, ls='--')
    ax.set_yscale('log')
    ax.set_xticks(range(len(samples)))
    ax.set_xticklabels(samples, rotation=70, ha='right')
    ax.set_ylabel('{} reads ({})'.format(args.gene, args.method))
    ax.spines[['top', 'right']].set_visible(False)
    fig.tight_layout()
    fig.savefig(prefix + '_gene_reads.pdf')
    plt.close(fig)


# ----------------------------------------------------------------------------------------------------
def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    prefix = os.path.join(args.outdir, '{}_{}'.format(args.gene, args.method))

    sheet = pd.read_csv(args.sample_sheet, sep='\t', dtype=str)
    ref_lines = gene_lines(args.ref_gtf, args.gene)
    if not ref_lines:
        sys.exit('{} not found in {}'.format(args.gene, args.ref_gtf))
    gene_ids = sorted({dict(ATTR.findall(l.split('\t')[8]))['gene_id'] for l in ref_lines})
    ref_tx = read_gtf_exons(args.ref_gtf, keep=lambda f, a: a.get('gene_name') == args.gene)
    strand = next(iter(ref_tx.values()))['strand']
    chrom = next(iter(ref_tx.values()))['chrom']
    gene_exons = sorted({e for r in ref_tx.values() for e in r['exons']})
    canonical = [t for t, r in ref_tx.items() if 'Ensembl_canonical' in r['tags']]
    canonical_id = canonical[0] if canonical else None

    matrix = pd.read_csv(args.matrix, sep='\t', dtype={0: str, 1: str, 2: str})
    matrix.columns = [c.lstrip('#') for c in matrix.columns]
    id_col, name_col, gene_col = matrix.columns[:3]
    present = [s for s in sheet['sample_id'] if s in matrix.columns]
    missing = [s for s in sheet['sample_id'] if s not in matrix.columns]
    if not present:
        sys.exit('no sample of the sample sheet is a column of {}'.format(args.matrix))
    if missing:
        print('NOTE: not in matrix (counted as 0 reads): {}'.format(', '.join(missing)), file=sys.stderr)

    all_tx = read_gtf_exons(args.isoform_gtf)
    gene_id_core = {strip_version(g) for g in gene_ids}

    def belongs(row):
        tid = row[id_col]
        ids = {strip_version(g) for g in str(row[gene_col]).split(',')}
        rec = all_tx.get(tid)
        if rec is None:
            return False
        if ids & gene_id_core:
            return True
        return rec['chrom'] == chrom and rec['strand'] == strand and overlaps(rec['exons'], gene_exons)
    keep = matrix.apply(belongs, axis=1)
    sub = matrix[keep].copy()
    if sub.empty:
        sys.exit('no isoform of {} in {}'.format(args.gene, args.matrix))
    sub = sub.drop_duplicates(subset=id_col).set_index(id_col)
    counts = sub[present].astype(float)
    for s in missing:
        counts[s] = 0.0
    counts = counts[sheet['sample_id'].tolist()]
    isoforms = {tid: all_tx[tid] for tid in counts.index}

    ann = annotate(isoforms, ref_tx, canonical_id)
    if args.find_as:
        as_summary = run_find_as(args.find_as, isoforms, ann, gene_ids[0], ref_lines, args.outdir)
        if as_summary is not None:
            ann['as_events_vs_canonical'] = ann.index.map(lambda t: as_summary.get(t, ''))

    gene_total = counts.sum(axis=0)
    props = counts.div(gene_total.replace(0, np.nan), axis=1)
    ok_samples = [s for s in counts.columns if gene_total[s] >= args.min_gene_reads]
    ann['total_reads'] = counts.sum(axis=1)
    ann['mean_prop'] = props[ok_samples].mean(axis=1) if ok_samples else np.nan
    ann['max_prop'] = props[ok_samples].max(axis=1) if ok_samples else np.nan
    ann['n_samples_detected'] = (counts[ok_samples] >= 2).sum(axis=1) if ok_samples else 0
    ann = ann.sort_values(['mean_prop', 'total_reads'], ascending=False)
    labels = label_isoforms(ann)
    ann.insert(0, 'label', [labels[t] for t in ann.index])

    # condition means
    long = props.T.join(sheet.set_index('sample_id')[['cell_type', 'condition']])
    long = long[long.index.isin(ok_samples)]
    cond_mean = long.groupby(['cell_type', 'condition']).mean().T
    cond_mean.columns = ['{}|{}'.format(c, d) for c, d in cond_mean.columns]

    # isoforms to plot
    big = ann[(ann['max_prop'] >= args.min_prop)].index.tolist()[:args.max_plot]
    if not big:
        big = ann.index.tolist()[:min(args.max_plot, len(ann))]
    plot_ids = big + (['Other'] if len(big) < len(ann) else [])
    props_plot = props.loc[big].copy()
    if 'Other' in plot_ids:
        props_plot.loc['Other'] = props.drop(index=big).sum(axis=0)
    labels['Other'] = 'other isoforms ({})'.format(len(ann) - len(big))

    tested = ann[ann['mean_prop'] >= args.min_test_prop].index.tolist()
    contrast_df = contrasts(counts, props, gene_total, sheet, tested, args.min_gene_reads)
    if not contrast_df.empty:
        contrast_df.insert(3, 'label', contrast_df['isoform_id'].map(labels))
        contrast_df = contrast_df.sort_values(['cell_type', 'contrast', 'cmh_p'])

    # tables
    out_counts = counts.copy()
    out_counts.insert(0, 'label', [labels[t] for t in out_counts.index])
    out_counts.to_csv(prefix + '_isoform_counts.tsv', sep='\t', index_label='isoform_id', float_format='%.2f')
    out_props = props.copy()
    out_props.insert(0, 'label', [labels[t] for t in out_props.index])
    out_props.to_csv(prefix + '_isoform_proportions.tsv', sep='\t', index_label='isoform_id', float_format='%.4f')
    ann.join(cond_mean).to_csv(prefix + '_isoform_annotation.tsv', sep='\t', index_label='isoform_id',
                               float_format='%.4f')
    depth = sheet.set_index('sample_id')[['donor', 'cell_type', 'condition']].copy()
    depth['gene_reads'] = gene_total
    depth['used_in_statistics'] = depth.index.isin(ok_samples)
    depth.to_csv(prefix + '_gene_reads.tsv', sep='\t', float_format='%.2f')
    if not contrast_df.empty:
        contrast_df.to_csv(prefix + '_contrasts.tsv', sep='\t', index=False, float_format='%.4g')

    plot_all(args, isoforms, plot_ids, labels, props_plot, counts, gene_total, sheet, contrast_df, strand)
    write_markdown(args, prefix, ann, depth, contrast_df, labels, gene_ids, canonical_id)
    print('report written to {}'.format(args.outdir))


def write_markdown(args, prefix, ann, depth, contrast_df, labels, gene_ids, canonical_id):
    lines = ['# {} isoforms ({})'.format(args.gene, args.method), '',
             'Gene: {}; canonical transcript: {} (marked * in labels).'.format(', '.join(gene_ids), canonical_id),
             '', '## Reads on {} per sample'.format(args.gene), '',
             '| sample | donor | cell type | condition | reads | in statistics |', '|---|---|---|---|---|---|']
    for sid, row in depth.iterrows():
        lines.append('| {} | {} | {} | {} | {:.0f} | {} |'.format(sid, row['donor'], row['cell_type'],
                                                                   row['condition'], row['gene_reads'],
                                                                   'yes' if row['used_in_statistics'] else 'no'))
    lines += ['', '## Isoforms (mean proportion over samples with >= {:.0f} reads)'.format(args.min_gene_reads),
              '', '| isoform | status | type | exons | mean % | max % | novel junctions | AS events vs canonical |',
              '|---|---|---|---|---|---|---|---|']
    for tid, row in ann.head(25).iterrows():
        lines.append('| {} | {} | {} | {} | {:.1f} | {:.1f} | {} | {} |'.format(
            row['label'], row['status'], row['transcript_type'], row['n_exons'], 100 * (row['mean_prop'] or 0),
            100 * (row['max_prop'] or 0), row['n_novel_junctions'], row.get('as_events_vs_canonical', '')))
    if not contrast_df.empty:
        lines += ['', '## Condition effects (vs Vehicle, paired by donor; exploratory)', '',
                  '| cell type | contrast | isoform | donors | mean delta % | per-donor delta | paired t p | CMH p | CMH FDR |',
                  '|---|---|---|---|---|---|---|---|---|']
        top = contrast_df[contrast_df['cmh_fdr'] < 0.1]
        if top.empty:
            top = contrast_df.groupby(['cell_type', 'contrast']).head(3)
        for _, r in top.iterrows():
            lines.append('| {} | {} | {} | {} | {:+.1f} | {} | {} | {:.2g} | {:.2g} |'.format(
                r['cell_type'], r['contrast'], r['label'], r['n_donors'], 100 * r['mean_delta_prop'],
                r['per_donor_delta'], 'n/a' if np.isnan(r['paired_t_p']) else '{:.2g}'.format(r['paired_t_p']),
                r['cmh_p'], r['cmh_fdr']))
    lines += ['', 'Figures: `{0}_structure_heatmap.pdf`, `{0}_stacked_proportions.pdf`, '
                  '`{0}_isoform_by_donor.pdf`, `{0}_gene_reads.pdf`.'.format(os.path.basename(prefix))]
    with open(prefix + '_report.md', 'w') as handle:
        handle.write('\n'.join(lines) + '\n')


if __name__ == '__main__':
    main()
