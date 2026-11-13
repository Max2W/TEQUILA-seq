#!/usr/bin/env python3
import argparse
import gzip
import random
import re
import sys
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pysam
from matplotlib.patches import Rectangle

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Liberation Sans", "Arial", "Helvetica", "DejaVu Sans"]
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["svg.fonttype"] = "none"

BASE_COLORS = {"A": "#00a600", "C": "#0000ff", "G": "#d17105", "T": "#ff0000", "N": "#888888"}
READ_COLOR = "#bdbdbd"
DEL_COLOR = "#000000"
INS_COLOR = "#8a2be2"
GAP_COLOR = "#9e9e9e"
COV_COLOR = "#a6a6a6"
TX_COLOR = "#1f3a93"
SKIP_FLAGS = 0x4 | 0x100 | 0x200 | 0x400
AF_THRESHOLD = 0.2
AF_MIN_DEPTH = 10
FONT = 7.0
PAGE_W, PAGE_H = 11.0, 8.5
DPI = 300


def parse_args():
    p = argparse.ArgumentParser(
        description="Draw an IGV-style view of a BAM region: a coverage track (positions where more than 20%% of "
                    "reads differ from the reference are coloured by base), the alignments (mismatches coloured "
                    "by base, deletions as black lines, insertions as purple ticks, introns as thin grey lines), "
                    "the reference sequence for regions of 150 bp or less, and optionally transcripts from a GTF. "
                    "Unmapped, secondary, QC-fail and duplicate reads are skipped. The page is 11 x 8.5 in "
                    "(US letter, landscape); the output format follows the --outfile extension (.pdf, .png, .svg).")
    p.add_argument("--bam", required=True, help="Indexed BAM or CRAM (CRAM needs --fasta).")
    p.add_argument("--region", required=True,
                   help="chrom:start-end, 1-based and inclusive; commas are allowed (chr1:1,000,000-1,002,000).")
    p.add_argument("--max-alignments", type=int, default=500,
                   help="Alignments drawn; a larger region is randomly downsampled to this many (default 500). "
                        "Coverage always uses every read.")
    p.add_argument("--outfile", required=True)
    p.add_argument("--fasta", required=True, help="Reference FASTA (indexed with samtools faidx).")
    p.add_argument("--gtf", default=None,
                   help="Annotation GTF (optionally .gz). Draws one transcript per gene overlapping the region (the "
                        "Ensembl_canonical one, or the longest).")
    p.add_argument("--min-indel", type=int, default=10,
                   help="Draw deletions and insertions at least this long (default 10); shorter ones are drawn as "
                        "plain read, which hides the small-indel noise of long reads.")
    p.add_argument("--spliced-only", action="store_true", help="Only use spliced reads (coverage and alignments).")
    p.add_argument("--title", default=None, help="Plot title (default: the BAM name and the region).")
    p.add_argument("--seed", type=int, default=1, help="Random seed for downsampling.")
    return p.parse_args()


def parse_region(text):
    m = re.fullmatch(r"\s*([^:\s]+):([\d,]+)-([\d,]+)\s*", text)
    if not m:
        sys.exit(f"ERROR: --region must look like chr1:1000-2000, got {text!r}")
    start, end = int(m.group(2).replace(",", "")), int(m.group(3).replace(",", ""))
    if end < start:
        sys.exit("ERROR: region end is before its start.")
    return m.group(1), start - 1, end


def resolve_contig(bam, chrom):
    names = set(bam.references)
    for c in (chrom, chrom[3:] if chrom.startswith("chr") else "chr" + chrom):
        if c in names:
            return c
    sys.exit(f"ERROR: contig {chrom!r} is not in the BAM header.")


def keep_read(spliced_only):
    def keep(r):
        if r.flag & SKIP_FLAGS:
            return False
        return not spliced_only or any(op == 3 for op, _ in (r.cigartuples or ()))
    return keep


def collect_reads(bam, chrom, start, end, max_reads, seed, keep):
    rng = random.Random(seed)
    chosen, n = [], 0
    for r in bam.fetch(chrom, start, end):
        if not keep(r):
            continue
        n += 1
        if len(chosen) < max_reads:
            chosen.append(r)
        else:
            k = rng.randrange(n)
            if k < max_reads:
                chosen[k] = r
    chosen.sort(key=lambda r: (r.reference_start, r.reference_end or 0))
    return chosen, n


def fetch_reference(fasta, chrom, start, end):
    names = set(fasta.references)
    c = chrom if chrom in names else (chrom[3:] if chrom.startswith("chr") else "chr" + chrom)
    if c not in names:
        print(f"WARNING: {chrom} is not in the FASTA; mismatches are not shown.", file=sys.stderr)
        return ""
    seq = fasta.fetch(c, start, end).upper()
    return seq + "N" * (end - start - len(seq))


_BASE_INDEX = np.full(256, -1, dtype=np.int64)
for _i, _b in enumerate("ACGT"):
    _BASE_INDEX[ord(_b)] = _i


def fmt_pos(x, _=None):
    return f"{int(x):,}"


def draw_coverage(ax, acgt, ref_seq, start, end):
    depth = acgt.sum(axis=0)
    L = end - start
    if L > 4000:
        edges = np.linspace(0, L, 2001).astype(int)
        mx = np.maximum.reduceat(depth, edges[:-1])
        ax.fill_between(np.repeat(start + edges, 2)[1:-1], np.repeat(mx, 2), color=COV_COLOR, linewidth=0,
                        step=None, rasterized=True)
    else:
        ax.bar(np.arange(start, end), depth, width=1.0, align="edge", color=COV_COLOR, linewidth=0)
    top = max(int(depth.max()) if depth.size else 0, 1)
    if ref_seq:
        ref_idx = _BASE_INDEX[np.frombuffer(ref_seq.encode("latin-1", "replace"), dtype=np.uint8)]
        refcount = np.where(ref_idx >= 0, acgt[np.clip(ref_idx, 0, 3), np.arange(L)], 0)
        with np.errstate(divide="ignore", invalid="ignore"):
            af = np.where(depth > 0, (depth - refcount) / depth, 0)
        w = max(1.0, L / 500.0)
        for i in np.where((af > AF_THRESHOLD) & (depth >= AF_MIN_DEPTH))[0]:
            bottom = 0
            for b, base in enumerate("ACGT"):
                c = acgt[b, i]
                if c:
                    ax.add_patch(Rectangle((start + i + 0.5 - w / 2, bottom), w, c, color=BASE_COLORS[base],
                                           linewidth=0))
                    bottom += c
    ax.set_xlim(start, end)
    ax.set_ylim(0, top)
    ax.set_yticks([])
    ax.text(0.003, 0.97, f"[0-{top:,}]", transform=ax.transAxes, va="top", ha="left", fontsize=FONT - 1)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_linewidth(0.4)
    ax.xaxis.tick_top()
    ax.xaxis.set_major_locator(plt.MaxNLocator(5, integer=True))
    ax.xaxis.set_major_formatter(plt.FuncFormatter(fmt_pos))
    ax.tick_params(axis="x", labelsize=FONT - 1, length=2, pad=1)
    ax.set_ylabel("coverage", fontsize=FONT, rotation=0, ha="right", va="center", labelpad=4, fontweight="bold")


def rgba(c):
    return tuple(int(round(v * 255)) for v in matplotlib.colors.to_rgb(c)) + (255,)


def pack_rows(reads, start, end):
    gap = max(1, int((end - start) * 0.003))
    row_ends, placed = [], []
    for r in reads:
        rs, re_ = r.reference_start, r.reference_end or r.reference_start
        row = next((i for i, e in enumerate(row_ends) if e + gap <= rs), None)
        if row is None:
            row = len(row_ends)
            row_ends.append(re_)
        else:
            row_ends[row] = re_
        placed.append((r, row))
    return placed, max(len(row_ends), 8)


def draw_alignments(ax, placed, nrows, n_total, ref_seq, start, end, width_in, height_in, dpi, min_indel):
    L = max(end - start, 1)
    W = max(min(int(round(width_in * dpi * 2)), 8000), 1)
    H = max(min(int(round(height_in * dpi * 2)), 8000), 1)
    img = np.full((H, W, 4), 255, dtype=np.uint8)
    row_px = H / nrows
    padpx = int(row_px * 0.1) if row_px >= 10 else 0
    ppb = W / L
    read_c, gap_c, del_c, ins_c = rgba(READ_COLOR), rgba(GAP_COLOR), rgba(DEL_COLOR), rgba(INS_COLOR)
    base_c = {b: rgba(c) for b, c in BASE_COLORS.items()}
    mm_w = max(1, int(round(ppb)), int(W / 1500))
    line_h = max(1, int(round(row_px * 0.15)))

    def cols(a, b):
        c0 = int((a - start) * ppb)
        return max(c0, 0), min(max(int(round((b - start) * ppb)), c0 + 1), W)

    for r, row in placed:
        r0 = int(round(row * row_px)) + padpx
        r1 = max(int(round((row + 1) * row_px)) - padpx, r0 + 1)
        mid = (r0 + r1) // 2
        l0, l1 = max(mid - line_h // 2, r0), min(mid - line_h // 2 + line_h, r1)
        seq = r.query_sequence
        pos, qpos = r.reference_start, 0
        for op, ln in r.cigartuples or []:
            if op in (0, 7, 8):
                a, b = max(pos, start), min(pos + ln, end)
                if a < b:
                    c0, c1 = cols(a, b)
                    img[r0:r1, c0:c1] = read_c
                    if seq and ref_seq:
                        q = np.frombuffer(seq[qpos + (a - pos): qpos + (b - pos)].upper().encode(), dtype="S1")
                        rf = np.frombuffer(ref_seq[a - start: b - start].encode(), dtype="S1")
                        for k in np.nonzero(q != rf)[0]:
                            c = int((a + k + 0.5 - start) * ppb)
                            c0m = max(c - mm_w // 2, 0)
                            img[r0:r1, c0m:min(c0m + mm_w, W)] = base_c.get(q[k].decode(), base_c["N"])
                pos += ln
                qpos += ln
            elif op == 1:
                if start <= pos <= end and ln >= min_indel:
                    c = min(int((pos - start) * ppb), W - 1)
                    img[r0:r1, max(c - 1, 0):c + 1] = ins_c
                qpos += ln
            elif op == 2:
                a, b = max(pos, start), min(pos + ln, end)
                if a < b:
                    c0, c1 = cols(a, b)
                    if ln >= min_indel:
                        img[l0:l1, c0:c1] = del_c
                    else:
                        img[r0:r1, c0:c1] = read_c
                pos += ln
            elif op == 3:
                a, b = max(pos, start), min(pos + ln, end)
                if a < b:
                    c0, c1 = cols(a, b)
                    img[l0:l1, c0:c1] = gap_c
                pos += ln
            elif op == 4:
                qpos += ln
    ax.imshow(img, extent=(start, end, nrows, 0), aspect="auto", interpolation="none")
    ax.set_xlim(start, end)
    ax.set_ylim(nrows, 0)
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_ylabel("alignments", fontsize=FONT, rotation=0, ha="right", va="center", labelpad=4, fontweight="bold")
    if n_total > len(placed):
        ax.text(0.995, 0.99, f"downsampled {n_total:,} → {len(placed):,}", transform=ax.transAxes, ha="right",
                va="top", fontsize=FONT - 1.5, color="#333",
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="#bbb", lw=0.4, alpha=0.9))


def draw_reference(ax, ref_seq, start, end):
    for i, b in enumerate(ref_seq):
        ax.text(start + i + 0.5, 0.5, b, ha="center", va="center", fontsize=FONT - 0.5,
                color=BASE_COLORS.get(b, "#555"), fontweight="bold", family="monospace")
    ax.set_xlim(start, end)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.text(-0.003, 0.5, "reference", transform=ax.transAxes, ha="right", va="center", fontsize=FONT,
            fontweight="bold")


_ATTR = {k: re.compile(k + r' "([^"]+)"') for k in ("gene_name", "gene_id", "transcript_id")}


def load_transcripts(gtf, chrom, start, end):
    alt = chrom[3:] if chrom.startswith("chr") else "chr" + chrom
    txs = defaultdict(lambda: {"exons": [], "cds": [], "canonical": False})
    opener = gzip.open if gtf.endswith(".gz") else open
    with opener(gtf, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[0] not in (chrom, alt) or parts[2] not in ("exon", "CDS"):
                continue
            s, e = int(parts[3]) - 1, int(parts[4])
            t = _ATTR["transcript_id"].search(parts[8])
            if not t:
                continue
            g = _ATTR["gene_name"].search(parts[8]) or _ATTR["gene_id"].search(parts[8])
            rec = txs[t.group(1)]
            rec["gene"], rec["strand"], rec["id"] = g.group(1) if g else t.group(1), parts[6], t.group(1)
            if 'tag "Ensembl_canonical"' in parts[8]:
                rec["canonical"] = True
            rec["exons" if parts[2] == "exon" else "cds"].append((s, e))
    by_gene = defaultdict(list)
    for r in txs.values():
        if not r["exons"]:
            continue
        r["exons"].sort()
        r["cds"].sort()
        if r["exons"][-1][1] <= start or r["exons"][0][0] >= end:
            continue
        by_gene[r["gene"]].append(r)
    out = []
    for gene in sorted(by_gene, key=lambda g: min(r["exons"][0][0] for r in by_gene[g])):
        rs = by_gene[gene]
        canon = [r for r in rs if r["canonical"]]
        out.append(canon[0] if canon else max(rs, key=lambda r: sum(e - s for s, e in r["exons"])))
    return out


def draw_transcript(ax, tx, start, end):
    ax.set_xlim(start, end)
    ax.set_ylim(0, 1)
    ax.axis("off")
    exons, cds = tx["exons"], tx["cds"]
    a, b = max(exons[0][0], start), min(exons[-1][1], end)
    if a < b:
        ax.plot([a, b], [0.5, 0.5], color=TX_COLOR, linewidth=0.6, zorder=1)
        marker = ">" if tx["strand"] == "+" else "<"
        n = max(3, int(round(ax.get_position().width * ax.figure.get_figwidth() * 1.3)))
        for x in np.linspace(a, b, n + 2)[1:-1]:
            inside = any(s <= x < e for s, e in exons)
            ax.plot(x, 0.5, marker=marker, markersize=2.2, color="white" if inside else TX_COLOR,
                    zorder=5 if inside else 2)
    for s, e in exons:
        if e > start and s < end:
            ax.add_patch(Rectangle((s, 0.32), e - s, 0.36, color=TX_COLOR, linewidth=0, zorder=3))
    for s, e in cds:
        if e > start and s < end:
            ax.add_patch(Rectangle((s, 0.15), e - s, 0.70, color=TX_COLOR, linewidth=0, zorder=3))
    label = tx["gene"] if tx["gene"] == tx["id"] else f"{tx['gene']}\n{tx['id']}"
    ax.text(-0.003, 0.5, label, transform=ax.transAxes, ha="right", va="center", fontsize=FONT - 1,
            linespacing=1.0)


def main():
    args = parse_args()
    fasta = pysam.FastaFile(args.fasta)
    mode = "rc" if args.bam.endswith(".cram") else "rb"
    bam = pysam.AlignmentFile(args.bam, mode, reference_filename=args.fasta)
    chrom, start, end = parse_region(args.region)
    contig = resolve_contig(bam, chrom)
    end = min(end, bam.get_reference_length(contig))
    if end <= start:
        sys.exit("ERROR: the region lies beyond the end of the contig.")

    keep = keep_read(args.spliced_only)
    acgt = np.array(bam.count_coverage(contig, start, end, quality_threshold=0, read_callback=keep), dtype=np.int64)
    reads, n_total = collect_reads(bam, contig, start, end, args.max_alignments, args.seed, keep)
    placed, nrows = pack_rows(reads, start, end)
    ref_seq = fetch_reference(fasta, contig, start, end)
    show_ref = bool(ref_seq) and end - start <= 150
    txs = load_transcripts(args.gtf, contig, start, end) if args.gtf else []

    margin_t, margin_b, left, right = 0.4, 0.4, 1.25, 0.4
    title_h, tick_h, cov_h, ref_h, gap = 0.3, 0.2, 1.1, (0.18 if show_ref else 0.0), 0.05
    tx_h = min(0.3, 1.5 / max(len(txs), 1)) if txs else 0.0
    aln_h = PAGE_H - margin_t - margin_b - title_h - tick_h - cov_h - ref_h - tx_h * len(txs) - 3 * gap
    width = PAGE_W - left - right
    fig = plt.figure(figsize=(PAGE_W, PAGE_H))

    def axes(top, h):
        return fig.add_axes([left / PAGE_W, (top - h) / PAGE_H, width / PAGE_W, h / PAGE_H])

    title = args.title if args.title is not None else \
        f"{args.bam.split('/')[-1]}   {contig}:{start + 1:,}-{end:,}  ({end - start:,} bp)"
    fig.text(left / PAGE_W, 1 - margin_t / PAGE_H, title, fontsize=FONT + 3, fontweight="bold", va="top")
    y = PAGE_H - margin_t - title_h - tick_h
    draw_coverage(axes(y, cov_h), acgt, ref_seq, start, end)
    y -= cov_h + gap
    draw_alignments(axes(y, aln_h), placed, nrows, n_total, ref_seq, start, end, width, aln_h, DPI, args.min_indel)
    y -= aln_h + gap
    if show_ref:
        draw_reference(axes(y, ref_h), ref_seq, start, end)
        y -= ref_h + gap
    for tx in txs:
        draw_transcript(axes(y, tx_h), tx, start, end)
        y -= tx_h
    fig.savefig(args.outfile, dpi=DPI)
    plt.close(fig)
    bam.close()
    print(f"{contig}:{start + 1}-{end}: {n_total:,} read(s), {len(reads):,} drawn in {nrows} row(s), "
          f"{len(txs)} transcript(s) -> {args.outfile}")


if __name__ == "__main__":
    main()