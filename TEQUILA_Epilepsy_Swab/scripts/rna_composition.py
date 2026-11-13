#!/usr/bin/env python3
"""What kind of RNA the swab reads come from: read categories, read lengths and coverage along transcripts, from the
genome-wide BAMs of the full path (07_qc_full.sh).

Question: is the low-RIN swab RNA mature mRNA from cells that is broken into pieces, or something else (pre-mRNA,
mitochondrial or ribosomal RNA, DNA, reads that do not map)? Reads cannot tell RNA inside cells from RNA outside
cells (cell-free RNA); that needs a wet-lab test (see the README).

Every primary read of a sample is put in one category:
  unmapped              not aligned to GRCh38
  mitochondrial         aligned to chrM
  rRNA                  overlaps a GENCODE rRNA / rRNA_pseudogene gene (the 45S rDNA repeats are not in the primary
                        assembly, so most 45S rRNA reads end up unmapped)
  exonic, spliced       >= 50% of aligned bases on exons of any gene, and at least one intron (N >= 30 nt) in the
                        alignment: processed mRNA
  exonic, unspliced     >= 50% of aligned bases on exons, no intron: short mRNA pieces, single-exon genes, or DNA
  intronic              alignment midpoint inside a gene, < 50% exonic: pre-mRNA or genomic DNA
  intergenic            outside annotated genes: genomic DNA, unannotated transcription
Each read is also marked panel (aligned span overlaps a probe-panel gene) or off-panel. Off-panel reads were not
selected by the probes, so they are the best view of the RNA that went into the library.

Coverage along transcripts: exonic reads of protein-coding genes whose canonical transcript (MANE Select, else
Ensembl_canonical) is >= --min-tx-length nt are placed on that transcript (the transcript with most overlapping exon
bases); the span from the first to the last covered transcript base is added to 100 bins from 5' to 3'.

Outputs (--outdir):
  composition.tsv        sample_id, scope (all / off_panel), category, reads, percent
  read_lengths.tsv       sample_id, category, read length quantiles (10, 25, 50, 75, 90%) and mean
  genebody_coverage.tsv  sample_id, scope, bin (0-99, 5' -> 3'), reads covering the bin
  transcript_span.tsv    sample_id, scope, reads, median fraction of the transcript a read covers, % covering >= 50%,
                         3'/5' coverage ratio (mean of the last 20 bins / first 20 bins)
"""
import argparse
import bisect
import csv
import os
from collections import Counter, defaultdict
from multiprocessing import Pool

import numpy as np
import pysam

BIN = 100_000
NBINS = 100
CATEGORIES = ["exonic, spliced", "exonic, unspliced", "intronic", "intergenic", "mitochondrial", "rRNA", "unmapped"]
LEN_EDGES = np.concatenate([np.arange(0, 2000, 10), np.arange(2000, 20001, 100)])
_G = None


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bam-pattern", required=True, help="path with {sample}")
    p.add_argument("--sample-sheet", required=True)
    p.add_argument("--gtf", required=True)
    p.add_argument("--panel", required=True, help="probe-panel gene symbols, one per line")
    p.add_argument("--min-tx-length", type=int, default=1500)
    p.add_argument("--max-reads", type=int, default=0, help="stop after this many primary reads per sample (0 = all)")
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--outdir", required=True)
    return p.parse_args()


def attr(a, key):
    k = key + ' "'
    return a.split(k, 1)[1].split('"', 1)[0] if k in a else None


def merge(iv):
    out = []
    for s, e in sorted(iv):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


def bin_index(items):
    """items: (chrom, start, end, payload) -> {(chrom, bin): [(start, end, payload)]}"""
    idx = defaultdict(list)
    for chrom, s, e, x in items:
        for b in range(s // BIN, (e - 1) // BIN + 1):
            idx[(chrom, b)].append((s, e, x))
    return idx


def load(gtf, panel, min_tx_length):
    genes, exons, rrna = [], defaultdict(list), []
    tx_info, tx_exons, tx_tags = {}, defaultdict(list), {}
    with open(gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            chrom, feat, s, e, strand, a = f[0], f[2], int(f[3]) - 1, int(f[4]), f[6], f[8]
            if feat == "gene":
                name, gtype = attr(a, "gene_name") or attr(a, "gene_id"), attr(a, "gene_type")
                genes.append((chrom, s, e, name))
                if gtype in ("rRNA", "rRNA_pseudogene"):
                    rrna.append((chrom, s, e, name))
            elif feat == "exon":
                exons[chrom].append((s, e))
                tid = attr(a, "transcript_id")
                tx_exons[tid].append((s, e))
                if tid not in tx_info:
                    tx_info[tid] = (chrom, strand, attr(a, "gene_name"), attr(a, "gene_type"))
                    tx_tags[tid] = ("MANE_Select" in a, "Ensembl_canonical" in a)
    merged = {}
    for chrom, iv in exons.items():
        m = merge(iv)
        merged[chrom] = ([x[0] for x in m], [x[1] for x in m])
    # one canonical transcript per protein-coding gene (MANE Select first)
    best = {}
    for tid, (chrom, strand, gname, gtype) in tx_info.items():
        mane, canon = tx_tags[tid]
        if gtype != "protein_coding" or not (mane or canon):
            continue
        if gname not in best or (mane and not tx_tags[best[gname]][0]):
            best[gname] = tid
    canon_items, tx_model = [], {}
    for gname, tid in best.items():
        chrom, strand, _, _ = tx_info[tid]
        ex = sorted(tx_exons[tid])
        length = sum(e - s for s, e in ex)
        if length < min_tx_length:
            continue
        # cumulative transcript offset of each exon, 5' -> 3'
        order = ex if strand == "+" else ex[::-1]
        off, cum = {}, 0
        for s, e in order:
            off[(s, e)] = cum
            cum += e - s
        tx_model[tid] = (strand, length, off, gname in panel)
        for s, e in ex:
            canon_items.append((chrom, s, e, (tid, s, e)))
    return {"genes": bin_index(genes), "rrna": bin_index(rrna), "exons": merged, "canon": bin_index(canon_items),
            "tx": tx_model, "panel": bin_index([g for g in genes if g[3] in panel])}


def overlaps(idx, chrom, s, e):
    seen = set()
    for b in range(s // BIN, (e - 1) // BIN + 1):
        for gs, ge, x in idx.get((chrom, b), ()):
            if gs < e and s < ge and (gs, ge, x) not in seen:
                seen.add((gs, ge, x))
                yield gs, ge, x


def exonic_bases(chrom, blocks):
    if chrom not in _G["exons"]:
        return 0
    starts, ends = _G["exons"][chrom]
    n = 0
    for s, e in blocks:
        i = bisect.bisect_right(starts, s) - 1
        i = max(i, 0)
        while i < len(starts) and starts[i] < e:
            n += max(0, min(e, ends[i]) - max(s, starts[i]))
            i += 1
    return n


def in_gene(chrom, pos):
    return any(True for _ in overlaps(_G["genes"], chrom, pos, pos + 1))


def tx_span(chrom, blocks):
    """(tid, first, last) transcript coordinates covered by the read on its best canonical transcript, or None."""
    per_tx = defaultdict(list)
    for s, e in blocks:
        for gs, ge, (tid, es, ee) in overlaps(_G["canon"], chrom, s, e):
            a, b = max(s, es), min(e, ee)
            if b > a:
                per_tx[tid].append((a, b, es, ee))
    if not per_tx:
        return None
    tid = max(per_tx, key=lambda t: sum(b - a for a, b, _, _ in per_tx[t]))
    strand, length, off, _ = _G["tx"][tid]
    pos = []
    for a, b, es, ee in per_tx[tid]:
        o = off[(es, ee)]
        if strand == "+":
            pos += [o + a - es, o + b - es]
        else:
            pos += [o + ee - b, o + ee - a]
    return tid, min(pos), max(pos)


def _init(g):
    global _G
    _G = g


def run(job):
    sample, path, max_reads = job
    cats = {"all": Counter(), "off_panel": Counter()}
    lengths = {c: np.zeros(len(LEN_EDGES), dtype=np.int64) for c in CATEGORIES}
    cov = {"all": np.zeros(NBINS), "off_panel": np.zeros(NBINS)}
    frac = {"all": [], "off_panel": []}
    n = 0
    with pysam.AlignmentFile(path) as bam:
        for r in bam.fetch(until_eof=True):
            if r.is_secondary or r.is_supplementary:
                continue
            n += 1
            if max_reads and n > max_reads:
                break
            qlen = r.query_length or r.infer_read_length() or 0
            panel = False
            if r.is_unmapped or r.reference_end is None:
                cat = "unmapped"
            else:
                chrom, s, e = r.reference_name, r.reference_start, r.reference_end
                panel = any(True for _ in overlaps(_G["panel"], chrom, s, e))
                if chrom in ("chrM", "MT", "M"):
                    cat = "mitochondrial"
                elif any(True for _ in overlaps(_G["rrna"], chrom, s, e)):
                    cat = "rRNA"
                else:
                    blocks = r.get_blocks()
                    aligned = sum(b - a for a, b in blocks) or 1
                    if exonic_bases(chrom, blocks) / aligned >= 0.5:
                        spliced = any(op == 3 and ln >= 30 for op, ln in (r.cigartuples or ()))
                        cat = "exonic, spliced" if spliced else "exonic, unspliced"
                        span = tx_span(chrom, blocks)
                        if span:
                            tid, a, b = span
                            length = _G["tx"][tid][1]
                            lo, hi = int(a / length * NBINS), int(np.ceil(b / length * NBINS))
                            scopes = ["all"] if _G["tx"][tid][3] else ["all", "off_panel"]
                            for sc in scopes:
                                cov[sc][lo:max(hi, lo + 1)] += 1
                                frac[sc].append((b - a) / length)
                    elif in_gene(chrom, (s + e) // 2):
                        cat = "intronic"
                    else:
                        cat = "intergenic"
            cats["all"][cat] += 1
            if not panel:
                cats["off_panel"][cat] += 1
            lengths[cat][min(np.searchsorted(LEN_EDGES, qlen, side="right") - 1, len(LEN_EDGES) - 1)] += 1
    return sample, cats, lengths, cov, frac


def quantiles(hist, qs):
    total = hist.sum()
    if not total:
        return ["NA"] * len(qs)
    c = np.cumsum(hist)
    return [int(LEN_EDGES[np.searchsorted(c, q * total)]) for q in qs]


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    with open(args.sample_sheet) as fh:
        samples = [r["sample_id"] for r in csv.DictReader(fh, delimiter="\t")]
    panel = {l.strip() for l in open(args.panel) if l.strip() and not l.startswith("#")} | {"AFG2A"}
    g = load(args.gtf, panel, args.min_tx_length)
    print("{} canonical protein-coding transcripts >= {} nt for transcript coverage".format(len(g["tx"]),
                                                                                         args.min_tx_length))
    jobs = [(s, args.bam_pattern.format(sample=s), args.max_reads) for s in samples
            if os.path.exists(args.bam_pattern.format(sample=s))]
    if not jobs:
        raise SystemExit("no BAM found for {}".format(args.bam_pattern))
    with Pool(min(args.threads, len(jobs)), initializer=_init, initargs=(g,)) as pool:
        res = pool.map(run, jobs, chunksize=1)

    out = lambda name: open(os.path.join(args.outdir, name), "w")
    with out("composition.tsv") as fh:
        fh.write("sample_id\tscope\tcategory\treads\tpercent\n")
        for sample, cats, _, _, _ in res:
            for sc in ("all", "off_panel"):
                tot = sum(cats[sc].values())
                for c in CATEGORIES:
                    fh.write("{}\t{}\t{}\t{}\t{:.3f}\n".format(sample, sc, c, cats[sc][c],
                                                                cats[sc][c] / tot * 100 if tot else 0))
    with out("read_lengths.tsv") as fh:
        fh.write("sample_id\tcategory\treads\tq10\tq25\tmedian\tq75\tq90\tmean\n")
        mids = np.append((LEN_EDGES[:-1] + LEN_EDGES[1:]) / 2, LEN_EDGES[-1])
        for sample, _, lengths, _, _ in res:
            for c in CATEGORIES:
                h = lengths[c]
                mean = "{:.0f}".format((h * mids).sum() / h.sum()) if h.sum() else "NA"
                fh.write("{}\t{}\t{}\t{}\t{}\n".format(sample, c, int(h.sum()),
                                                       "\t".join(map(str, quantiles(h, [.1, .25, .5, .75, .9]))), mean))
    with out("genebody_coverage.tsv") as fh:
        fh.write("sample_id\tscope\tbin\treads\n")
        for sample, _, _, cov, _ in res:
            for sc in ("all", "off_panel"):
                for i, v in enumerate(cov[sc]):
                    fh.write("{}\t{}\t{}\t{:.0f}\n".format(sample, sc, i, v))
    with out("transcript_span.tsv") as fh:
        fh.write("sample_id\tscope\treads\tmedian_tx_fraction\tpercent_reads_ge50\tratio_3p_5p\n")
        for sample, _, _, cov, frac in res:
            for sc in ("all", "off_panel"):
                f = np.array(frac[sc])
                c = cov[sc]
                ratio = c[80:].mean() / c[:20].mean() if c[:20].mean() > 0 else float("nan")
                fh.write("{}\t{}\t{}\t{}\t{}\t{:.2f}\n".format(
                    sample, sc, len(f), "{:.3f}".format(np.median(f)) if len(f) else "NA",
                    "{:.1f}".format((f >= 0.5).mean() * 100) if len(f) else "NA", ratio))
    print("RNA composition -> {}".format(args.outdir))
    for sample, cats, _, _, _ in res:
        tot = sum(cats["off_panel"].values()) or 1
        print("  {}: off-panel reads {:,}; ".format(sample, tot) + ", ".join(
            "{} {:.1f}%".format(c, cats["off_panel"][c] / tot * 100) for c in CATEGORIES))


if __name__ == "__main__":
    main()
